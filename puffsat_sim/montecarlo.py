"""Monte Carlo / open-loop dispersion harness — the JVM-side run loop (ADR 0002).

``run_ensemble`` samples per-run inputs (:mod:`puffsat_sim.dispersion`), builds the
perturbed full-force run, applies the injection Δv to the apogee deployment state,
propagates to the 200 km interception crossing, and records the miss (in the
nominal-crossing RTN frame), the time-of-arrival error, and the osculating perigee.

The Stage-1 capstone (design doc §13) is this harness with ``control=None``; Rung D
supplies a controller through the same hook (§14.1).  Per-run replay (§14.2):
``replay_inputs(master_seed, spec, run_index)`` reconstructs any run's draws.

The harness surface is public (ADR 0017): the per-slice run modules in
:mod:`puffsat_sim.runs` drive :func:`run_record` over :func:`build_context` with
their own input generators and a :class:`RunVariant`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.attitudes import FrameAlignedProvider
from org.orekit.forces.maneuvers import ConstantThrustManeuver
from org.orekit.frames import FramesFactory
from org.orekit.orbits import CartesianOrbit
from org.orekit.utils import Constants, TimeStampedPVCoordinates

from org.hipparchus.geometry.euclidean.threed import Vector3D

from puffsat_sim import mission, presets
from puffsat_sim.actuator import Actuator, plan_burn
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.control import ControlPlan, Controller, Target, passes_toa_gate
from puffsat_sim.descent import (
    COAST_MAX_STEP_S,
    Crossing,
    apogee_state,
    descend,
    earth_model,
    to_absolute_date,
    vec3,
)
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
from puffsat_sim.navigation import Vec6
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import build_propagator
from puffsat_sim.records import EnsembleResult, RunRecord
from puffsat_sim.sink import append_record, guard_manifest, plan_resume, read_records

# Bump when a model change (physics, actuator realization, truth path) makes prior sink
# records incomparable, so a resume onto a stale sink fails loud (guard_manifest) rather than
# mixing incompatible records.  It is a manual proxy for a source fingerprint, not automatic.
MODEL_VERSION: str = "1"

# Finite-burn execution (B1, ADR 0008): the propagator runs at a fictitious 1 kg so the
# lumped Cd·(A/m) / Cr·(A/m) scale drag/SRP correctly, so the burn thrust is scaled to that
# mass (F·m_p/m_wet) — reproducing the real a=F/m and burn duration of the 25 kg / 400 mN
# actuator — and fires at a sentinel Isp so the executed arc is constant-mass (Isp-free
# trajectory).  Real propellant is the pure Tsiolkovsky transform at the actuator's Isp
# (puffsat_sim.actuator), per ADR 0004 decision 2.  Real-mass depletion stays unmodelled:
# B3a/C3a falsified the "large anti-drag burn" premise (~1.5 g, ADR 0014 findings).
PROPAGATOR_MASS_KG: float = 1.0
BURN_ISP_SENTINEL_S: float = 1.0e12


def scale_thrust_for_propagator_mass(thrust_n: float, wet_mass_kg: float) -> float:
    """Scale a commanded thrust from the real wet mass onto ``PROPAGATOR_MASS_KG`` (B1).

    Shared by every finite-burn call site (the A1 corrector's burn, the C3a feedforward,
    the C3b terminal loop) so the F·m_p/m_wet scaling is written once.
    """
    return thrust_n * PROPAGATOR_MASS_KG / wet_mass_kg


# Nominal truth inputs (ADR 0011 decision 5): zero injection + nominal coefficients
# (perfect model), so x_true is the nominal apogee state and the only predict-vs-execute
# divergence is whatever the RunVariant injects.
_NOMINAL_CD_AREA_OVER_MASS: float = 0.04
_NOMINAL_CR_AREA_OVER_MASS: float = 0.02
_NOMINAL_F10P7: float = 150.0
_NOMINAL_AP: float = 15.0


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


def nominal_inputs(run_index: int) -> RunInputs:
    """Zero-injection, nominal-coefficient truth inputs (ADR 0011 decision 5)."""
    return RunInputs(
        run_index=run_index,
        dv_rtn_m_s=(0.0, 0.0, 0.0),
        cd_area_over_mass=_NOMINAL_CD_AREA_OVER_MASS,
        cr_area_over_mass=_NOMINAL_CR_AREA_OVER_MASS,
        f10p7=_NOMINAL_F10P7,
        ap=_NOMINAL_AP,
    )


@dataclass(frozen=True)
class RunContext:
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
    nominal: Crossing
    nominal_basis: Basis
    target: Target


@dataclass(frozen=True)
class RunVariant:
    """How one run's plan is solved and executed — the rung knobs bundled (ADR 0017).

    Each field is one rung's default-off divergence between predict and execute:

    ``control`` — the §14.1 hook (A1, ADR 0003): ``None`` is open loop; a ``Controller``
    solves the plan against predict and the harness executes it against truth.

    ``toa_window_s`` — the A3 spurious-far-root gate (ADR 0007 decision 3iii): a
    converged plan whose crossing falls outside ±window of the nominal ToA is
    recorded non-converged.

    ``actuator`` — B1 (ADR 0008): execute the commanded Δv as a finite burn while
    predict stays impulsive; the residual miss is the actuator-realism erosion.

    ``nav_offset_rtn6`` — C0 (ADR 0011): predict-side apogee-RTN nav-error offset
    (position R/T/N then velocity R/T/N); the corrector plans from ``x_true + offset``
    while execute stays on truth, so the residual is the sensitivity Φ times the
    nav error.
    """

    control: Controller | None = None
    toa_window_s: float | None = None
    actuator: Actuator | None = None
    nav_offset_rtn6: Vec6 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def _finite_burn_maneuver(actuator: Actuator, correction_rtn: Vec3, ctx: RunContext) -> Any:
    """The ConstantThrustManeuver that executes the corrector's Δv as a finite burn (B1).

    Thrust is scaled to the propagator's fictitious 1 kg (F·m_p/m_wet) so a=F/m and the burn
    duration match the real 25 kg / 400 mN actuator; the sentinel Isp keeps the executed arc
    constant-mass.  Direction is the inertial Δv, held fixed by a frame-aligned attitude.
    """
    burn = plan_burn(actuator, correction_rtn)
    corr_eme = rtn_to_cartesian(correction_rtn, ctx.apo_basis)
    norm = (corr_eme[0] ** 2 + corr_eme[1] ** 2 + corr_eme[2] ** 2) ** 0.5  # nonzero: caller guards
    direction = Vector3D(corr_eme[0] / norm, corr_eme[1] / norm, corr_eme[2] / norm)
    thrust_eff = scale_thrust_for_propagator_mass(actuator.max_thrust_n, actuator.wet_mass_kg)
    return ConstantThrustManeuver(
        ctx.apo_date,
        burn.duration_s,
        thrust_eff,
        BURN_ISP_SENTINEL_S,
        FrameAlignedProvider(ctx.frame),
        direction,
    )


def run_record(ctx: RunContext, inputs: RunInputs, variant: RunVariant) -> RunRecord:
    """Propagate one run: apply injection, solve+execute the control plan, record the miss.

    ``predict`` (the corrector's onboard model) and ``execute`` (truth) are the same
    full-force physics at Rung A (ADR 0003), so a converged plan lands the recorded
    crossing on the nominal aim to machine precision.  The injection Δv is baked into
    the closure, so the corrector solves for the *correction* alone, starting from zero.
    The per-rung predict/execute divergences are the :class:`RunVariant` knobs.
    """
    physics = physics_from_inputs(inputs)
    injection_dv_eme = rtn_to_cartesian(inputs.dv_rtn_m_s, ctx.apo_basis)
    injection = Vector3D(injection_dv_eme[0], injection_dv_eme[1], injection_dv_eme[2])

    def make_crossing(
        correction_rtn: Vec3, apo_offset_rtn6: Vec6 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ) -> Crossing:
        pos_off = rtn_to_cartesian(
            (apo_offset_rtn6[0], apo_offset_rtn6[1], apo_offset_rtn6[2]), ctx.apo_basis
        )
        vel_off = rtn_to_cartesian(
            (apo_offset_rtn6[3], apo_offset_rtn6[4], apo_offset_rtn6[5]), ctx.apo_basis
        )
        corr_eme = rtn_to_cartesian(correction_rtn, ctx.apo_basis)
        position = ctx.apo_pos.add(Vector3D(pos_off[0], pos_off[1], pos_off[2]))
        vel = (
            ctx.apo_vel.add(injection)
            .add(Vector3D(corr_eme[0], corr_eme[1], corr_eme[2]))
            .add(Vector3D(vel_off[0], vel_off[1], vel_off[2]))
        )
        orbit = CartesianOrbit(
            TimeStampedPVCoordinates(ctx.apo_date, position, vel), ctx.frame, ctx.mu
        )
        return descend(orbit, physics, ctx.epoch, ctx.period, ctx.earth)

    if variant.control is None:
        plan = ControlPlan(actions=(), converged=True, iterations=0)
    else:
        plan = variant.control(
            lambda c: make_crossing(c, variant.nav_offset_rtn6).position_m, ctx.target
        )

    # Execute the commanded plan against truth.  A1's single action is at the apogee
    # node (elapsed_s=0), so it folds into the initial velocity; downstream multi-node
    # execution (ImpulseManeuver events) is an A2 addition.
    applied_rtn: Vec3 = plan.actions[0].dv_rtn_m_s if plan.actions else (0.0, 0.0, 0.0)
    if variant.actuator is not None and plan.actions and plan.actions[0].dv_mag_m_s > 0.0:
        # B1: fire the correction as a finite burn (injection alone folds into the velocity).
        orbit = CartesianOrbit(
            TimeStampedPVCoordinates(ctx.apo_date, ctx.apo_pos, ctx.apo_vel.add(injection)),
            ctx.frame,
            ctx.mu,
        )
        crossing = descend(
            orbit,
            physics,
            ctx.epoch,
            ctx.period,
            ctx.earth,
            maneuver=_finite_burn_maneuver(variant.actuator, applied_rtn, ctx),
        )
    else:
        crossing = make_crossing(applied_rtn)

    miss_vec: Vec3 = (
        crossing.position_m[0] - ctx.nominal.position_m[0],
        crossing.position_m[1] - ctx.nominal.position_m[1],
        crossing.position_m[2] - ctx.nominal.position_m[2],
    )
    toa_miss_s = crossing.toa_s - ctx.nominal.toa_s
    return RunRecord(
        inputs=inputs,
        miss_rtn_m=rtn_components(miss_vec, ctx.nominal_basis),
        toa_miss_s=toa_miss_s,
        perigee_alt_m=crossing.perigee_alt_m,
        crossing_position_m=crossing.position_m,
        crossing_velocity_m_s=crossing.velocity_m_s,
        control_log=plan.actions,
        total_dv_m_s=plan.total_dv_m_s,
        converged=passes_toa_gate(plan.converged, toa_miss_s, variant.toa_window_s),
        iterations=plan.iterations,
    )


def build_context(orbital_config: OrbitalConfig) -> RunContext:
    """Build the per-run constants shared by every run (the nominal crossing + apogee frame).

    Factored out of ``run_ensemble`` so the deterministic graders (``runs.sweep``,
    ``runs.navigation``) reuse the exact same nominal-crossing setup and
    :class:`RunContext` (ADR 0007 decision 5).
    """
    earth = earth_model()
    epoch = to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)
    frame = FramesFactory.getEME2000()
    mu: float = Constants.WGS84_EARTH_MU

    # Nominal (unperturbed) crossing — the reference the miss is measured against and
    # the corrector's target.
    nominal_orbit = (
        build_propagator(orbital_config, presets.full_force(), COAST_MAX_STEP_S)
        .getInitialState()
        .getOrbit()
    )
    nominal = descend(nominal_orbit, presets.full_force(), epoch, period, earth)
    nominal_basis: Basis = rtn_basis(nominal.position_m, nominal.velocity_m_s)

    apo_date, apo_pos, apo_vel = apogee_state(orbital_config)
    apo_basis: Basis = rtn_basis(vec3(apo_pos), vec3(apo_vel))

    return RunContext(
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


def _experiment_manifest(
    master_seed: int,
    spec: DispersionSpec,
    orbital_config: OrbitalConfig,
    variant: RunVariant,
) -> dict[str, str]:
    """The sink's experiment identity (see ``guard_manifest``).

    ``n`` is deliberately absent — extending an ensemble to a larger ``n`` on the same
    experiment is a valid resume.  The controller id is module.qualname, so it does not
    distinguish two controllers that differ only by closed-over parameters.
    """
    control = variant.control
    if control is None:
        control_id = "open-loop"
    else:
        module = getattr(control, "__module__", "?")
        qualname = getattr(control, "__qualname__", repr(control))
        control_id = f"{module}.{qualname}"
    return {
        "model_version": MODEL_VERSION,
        "master_seed": str(master_seed),
        "spec": repr(spec),
        "orbital_config": repr(orbital_config),
        "control": control_id,
        "actuator": "impulsive" if variant.actuator is None else repr(variant.actuator),
    }


def run_ensemble(
    spec: DispersionSpec,
    n: int,
    master_seed: int,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    control: Controller | None = None,
    sink_path: Path | None = None,
    actuator: Actuator | None = None,
) -> EnsembleResult:
    """Run a dispersion ensemble and aggregate it.

    ``control`` is the §14.1 hook (ADR 0003): ``None`` is the open-loop capstone; a
    ``Controller`` (e.g. the Rung A1 ``solve_apogee_correction``) closes the loop, each
    run solving and executing its plan against the same full-force truth.

    ``actuator`` (B1, ADR 0008) executes each commanded Δv as a finite mass-depleting burn
    instead of an impulse (``predict`` stays impulsive); the residual miss is then the
    actuator-realism erosion.  Requires a ``control``; ``None`` keeps impulsive execution.

    ``sink_path`` enables run-granular checkpoint/resume: completed records stream to a
    JSONL sink keyed by ``run_index``; on restart only the missing indices are run.  A
    resume is admitted only if the sink's replayed inputs match this ``master_seed``/``spec``
    (``plan_resume``) *and* its sidecar manifest matches this controller / actuator / orbital
    config / ``MODEL_VERSION`` (``guard_manifest``) — so a stale sink cannot be silently
    extended under a different experiment.
    """
    ctx = build_context(orbital_config)
    variant = RunVariant(control=control, actuator=actuator)

    reuse: dict[int, RunRecord] = {}
    todo = list(range(n))
    if sink_path is not None:
        guard_manifest(sink_path, _experiment_manifest(master_seed, spec, orbital_config, variant))
        reuse, todo = plan_resume(read_records(sink_path), master_seed, spec, n)
    records_by_index: dict[int, RunRecord] = dict(reuse)
    for i in todo:
        record = run_record(ctx, replay_inputs(master_seed, spec, i), variant)
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
        nominal_perigee_alt_m=ctx.nominal.perigee_alt_m,
        nominal_toa_s=ctx.nominal.toa_s,
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
