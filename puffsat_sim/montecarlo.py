"""Monte Carlo / open-loop dispersion harness — the JVM-side run loop (ADR 0002).

``run_ensemble`` samples per-run inputs (:mod:`puffsat_sim.dispersion`), builds the
perturbed full-force run, applies the injection Δv to the apogee deployment state,
propagates to the 200 km interception crossing, and records the miss (in the
nominal-crossing RTN frame), the time-of-arrival error, and the osculating perigee.

The Stage-1 capstone (design doc §13) is this harness with ``control=None``; Rung D
supplies a controller through the same hook (§14.1).  Per-run replay (§14.2):
``replay_inputs(master_seed, spec, run_index)`` reconstructs any run's draws.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.frames import FramesFactory
from org.orekit.orbits import CartesianOrbit, KeplerianOrbit
from org.orekit.propagation.events import AltitudeDetector
from org.orekit.propagation.events.handlers import StopOnDecreasing
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions, TimeStampedPVCoordinates

from org.hipparchus.geometry.euclidean.threed import Vector3D

from puffsat_sim import mission, presets
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.constants import EARTH_RADIUS_M
from puffsat_sim.control import ControlPlan, Controller, Target
from puffsat_sim.dispersion import (
    Basis,
    DispersionSpec,
    RunInputs,
    Vec3,
    replay_inputs,
    rtn_basis,
    rtn_components,
    rtn_to_cartesian,
    summarize,
)
from puffsat_sim.forces import (
    AtmosphericDrag,
    Geopotential,
    Relativity,
    SolarRadiation,
    ThirdBody,
)
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import build_propagator, build_propagator_from_orbit
from puffsat_sim.records import EnsembleResult, RunRecord
from puffsat_sim.sink import append_record, plan_resume, read_records

# Caps the adaptive integrator step on the terminal descent so a smooth low-drag arc
# cannot overstep the 200 km altitude event below the surface ("point is inside
# ellipsoid").  Interim fix for the §6.2 fragility, pending regime-switched
# propagation; nominal and perturbed runs share it so the miss stays common-mode.
_TERMINAL_MAX_STEP_S: float = 30.0


def physics_from_inputs(
    inputs: RunInputs, geopotential_degree: int = 8, geopotential_order: int = 8
) -> PhysicsConfig:
    """Full-force truth config carrying this run's drawn coefficients and space weather."""
    return PhysicsConfig(
        (
            Geopotential(degree=geopotential_degree, order=geopotential_order),
            ThirdBody(),
            SolarRadiation(cr_area_over_mass=inputs.cr_area_over_mass),
            AtmosphericDrag(
                cd_area_over_mass=inputs.cd_area_over_mass, f10p7=inputs.f10p7, ap=inputs.ap
            ),
            Relativity(),
        )
    )


def _to_absolute_date(dt: datetime) -> Any:
    utc = TimeScalesFactory.getUTC()
    return AbsoluteDate(dt.year, dt.month, dt.day, dt.hour, dt.minute, float(dt.second), utc)


def _earth() -> Any:
    return OneAxisEllipsoid(
        Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
        Constants.WGS84_EARTH_FLATTENING,
        FramesFactory.getITRF(IERSConventions.IERS_2010, True),
    )


@dataclass(frozen=True)
class _Crossing:
    position_m: Vec3
    velocity_m_s: Vec3
    toa_s: float
    perigee_alt_m: float


def _vec3(v: Any) -> Vec3:
    return (float(v.getX()), float(v.getY()), float(v.getZ()))


def _propagate_to_interception(propagator: Any, epoch: Any, period: float, earth: Any) -> _Crossing:
    """Stop the descending arc at the 200 km crossing and read off the state + perigee."""
    propagator.addEventDetector(
        AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(
            StopOnDecreasing()  # type: ignore[no-untyped-call]
        )
    )
    state = propagator.propagate(epoch.shiftedBy(period))
    pv = state.getPVCoordinates()
    orbit = KeplerianOrbit(state.getOrbit())
    perigee_alt = float(orbit.getA()) * (1.0 - float(orbit.getE())) - EARTH_RADIUS_M
    return _Crossing(
        position_m=_vec3(pv.getPosition()),
        velocity_m_s=_vec3(pv.getVelocity()),
        toa_s=float(state.getDate().durationFrom(epoch)),
        perigee_alt_m=perigee_alt,
    )


def _apogee_state(orbital_config: OrbitalConfig) -> tuple[Any, Any, Any]:
    """Nominal deployment state at apogee (epoch, mean anomaly π): (date, position, velocity)."""
    state = build_propagator(orbital_config, presets.two_body()).getInitialState()
    pv = state.getPVCoordinates()
    return state.getDate(), pv.getPosition(), pv.getVelocity()


@dataclass(frozen=True)
class _RunContext:
    """Per-ensemble constants shared by every run (built once before the loop)."""

    apo_date: Any
    apo_pos: Any
    apo_vel: Any
    apo_basis: Basis
    frame: Any
    mu: float
    epoch: Any
    period: float
    earth: Any
    nominal: _Crossing
    nominal_basis: Basis
    target: Target


def _run_record(ctx: _RunContext, inputs: RunInputs, control: Controller | None) -> RunRecord:
    """Propagate one run: apply injection, solve+execute the control plan, record the miss.

    ``predict`` (the corrector's onboard model) and ``execute`` (truth) are the same
    full-force physics at Rung A (ADR 0003), so a converged plan lands the recorded
    crossing on the nominal aim to machine precision.  The injection Δv is baked into
    the closure, so the corrector solves for the *correction* alone, starting from zero.
    """
    physics = physics_from_inputs(inputs)
    injection_dv_eme = rtn_to_cartesian(inputs.dv_rtn_m_s, ctx.apo_basis)

    def make_crossing(correction_rtn: Vec3) -> _Crossing:
        corr_eme = rtn_to_cartesian(correction_rtn, ctx.apo_basis)
        vel = ctx.apo_vel.add(
            Vector3D(injection_dv_eme[0], injection_dv_eme[1], injection_dv_eme[2])
        ).add(Vector3D(corr_eme[0], corr_eme[1], corr_eme[2]))
        orbit = CartesianOrbit(
            TimeStampedPVCoordinates(ctx.apo_date, ctx.apo_pos, vel), ctx.frame, ctx.mu
        )
        prop = build_propagator_from_orbit(orbit, physics, _TERMINAL_MAX_STEP_S)
        return _propagate_to_interception(prop, ctx.epoch, ctx.period, ctx.earth)

    if control is None:
        plan = ControlPlan(actions=(), converged=True, iterations=0)
    else:
        plan = control(lambda c: make_crossing(c).position_m, ctx.target)

    # Execute the commanded plan against truth.  A1's single action is at the apogee
    # node (elapsed_s=0), so it folds into the initial velocity; downstream multi-node
    # execution (ImpulseManeuver events) is an A2 addition.
    applied_rtn: Vec3 = plan.actions[0].dv_rtn_m_s if plan.actions else (0.0, 0.0, 0.0)
    crossing = make_crossing(applied_rtn)

    miss_vec: Vec3 = (
        crossing.position_m[0] - ctx.nominal.position_m[0],
        crossing.position_m[1] - ctx.nominal.position_m[1],
        crossing.position_m[2] - ctx.nominal.position_m[2],
    )
    return RunRecord(
        inputs=inputs,
        miss_rtn_m=rtn_components(miss_vec, ctx.nominal_basis),
        toa_miss_s=crossing.toa_s - ctx.nominal.toa_s,
        perigee_alt_m=crossing.perigee_alt_m,
        crossing_position_m=crossing.position_m,
        crossing_velocity_m_s=crossing.velocity_m_s,
        control_log=plan.actions,
        total_dv_m_s=plan.total_dv_m_s,
        converged=plan.converged,
        iterations=plan.iterations,
    )


def run_ensemble(
    spec: DispersionSpec,
    n: int,
    master_seed: int,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    control: Controller | None = None,
    sink_path: Path | None = None,
) -> EnsembleResult:
    """Run a dispersion ensemble and aggregate it.

    ``control`` is the §14.1 hook (ADR 0003): ``None`` is the open-loop capstone; a
    ``Controller`` (e.g. the Rung A1 ``solve_apogee_correction``) closes the loop, each
    run solving and executing its plan against the same full-force truth.

    ``sink_path`` enables run-granular checkpoint/resume: completed records stream to a
    JSONL sink keyed by ``run_index``; on restart only the missing indices are run
    (the present ones must match this ``master_seed``/``spec``).  The caller must resume
    with the same ``control`` — inputs are control-independent, so it is not auto-checked.
    """
    earth = _earth()
    epoch = _to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)
    frame = FramesFactory.getEME2000()
    mu: float = Constants.WGS84_EARTH_MU

    # Nominal (unperturbed) crossing — the reference the miss is measured against and
    # the corrector's target.
    nominal_prop = build_propagator(orbital_config, presets.full_force(), _TERMINAL_MAX_STEP_S)
    nominal = _propagate_to_interception(nominal_prop, epoch, period, earth)
    nominal_basis: Basis = rtn_basis(nominal.position_m, nominal.velocity_m_s)

    apo_date, apo_pos, apo_vel = _apogee_state(orbital_config)
    apo_basis: Basis = rtn_basis(_vec3(apo_pos), _vec3(apo_vel))

    ctx = _RunContext(
        apo_date=apo_date,
        apo_pos=apo_pos,
        apo_vel=apo_vel,
        apo_basis=apo_basis,
        frame=frame,
        mu=mu,
        epoch=epoch,
        period=period,
        earth=earth,
        nominal=nominal,
        nominal_basis=nominal_basis,
        target=Target(nominal.position_m),
    )

    reuse: dict[int, RunRecord] = {}
    todo = list(range(n))
    if sink_path is not None:
        reuse, todo = plan_resume(read_records(sink_path), master_seed, spec, n)
    records_by_index: dict[int, RunRecord] = dict(reuse)
    for i in todo:
        record = _run_record(ctx, replay_inputs(master_seed, spec, i), control)
        if sink_path is not None:
            append_record(sink_path, record)
        records_by_index[i] = record
    records = [records_by_index[i] for i in range(n)]

    stats = summarize(
        np.array([r.miss_rtn_m for r in records], dtype=np.float64).reshape(n, 3),
        np.array([r.toa_miss_s for r in records], dtype=np.float64),
        np.array([r.perigee_alt_m for r in records], dtype=np.float64),
        np.array([r.total_dv_m_s for r in records], dtype=np.float64),
        np.array([r.converged for r in records], dtype=np.bool_),
    )
    return EnsembleResult(
        master_seed=master_seed,
        nominal_perigee_alt_m=nominal.perigee_alt_m,
        nominal_toa_s=nominal.toa_s,
        records=tuple(records),
        stats=stats,
    )


def format_summary(result: EnsembleResult) -> str:
    """Human-readable one-screen summary of an ensemble (the capstone / Rung A report)."""
    s = result.stats
    mean = s.miss_rtn_mean_m
    std = s.miss_rtn_std_m
    controlled = any(r.control_log for r in result.records)
    title = "Closed-loop dispersion ensemble" if controlled else "Open-loop dispersion capstone"
    lines = [
        f"{title} — N={s.n}, master_seed={result.master_seed}",
        f"  Nominal: perigee {result.nominal_perigee_alt_m / 1e3:.1f} km,"
        f" coast {result.nominal_toa_s / 3600:.2f} h",
        "  Interception miss vs nominal, RTN frame [m] (T = dr_p/dv_a lever):",
        f"    bias R/T/N = {mean[0]:+.1f} / {mean[1]:+.1f} / {mean[2]:+.1f}",
        f"    std  R/T/N = {std[0]:.1f} / {std[1]:.1f} / {std[2]:.1f}",
        f"  Time-of-arrival miss: {s.toa_miss_mean_s:+.2f} ± {s.toa_miss_std_s:.2f} s",
        f"  Perigee (diagnostic, low=good): {s.perigee_alt_mean_m / 1e3:.1f}"
        f" ± {s.perigee_alt_std_m / 1e3:.2f} km"
        f" [min {s.perigee_alt_min_m / 1e3:.1f}, max {s.perigee_alt_max_m / 1e3:.1f}]",
    ]
    if controlled:
        lines.append(
            f"  Correction Δv [m/s]: mean {s.total_dv_mean_m_s:.4f},"
            f" max {s.total_dv_max_m_s:.4f} (std {s.total_dv_std_m_s:.4f})"
        )
        lines.append(f"  Corrector converged: {s.converged_fraction * 100:.0f}% of {s.n} runs")
    return "\n".join(lines)


def main() -> None:
    # Smoke-sized ensemble (design doc §10.4: N=50 is a smoke test, not the result;
    # the resolved-tail controllability result needs N=10³–10⁴, a longer job).
    result = run_ensemble(DispersionSpec(), n=50, master_seed=20260608)
    print(format_summary(result))


if __name__ == "__main__":
    main()
