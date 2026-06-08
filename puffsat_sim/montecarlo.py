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

from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
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
from puffsat_sim.dispersion import (
    Basis,
    DispersionSpec,
    EnsembleStats,
    RunInputs,
    Vec3,
    rtn_basis,
    rtn_components,
    rtn_to_cartesian,
    sample_run_inputs,
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

# Optional controller hook (§14.1).  Absent for the open-loop capstone; Rung D fills
# it.  Left untyped beyond a callable since the controller interface is a Rung-D design.
Controller = Callable[..., Any]

# Caps the adaptive integrator step on the terminal descent so a smooth low-drag arc
# cannot overstep the 200 km altitude event below the surface ("point is inside
# ellipsoid").  Interim fix for the §6.2 fragility, pending regime-switched
# propagation; nominal and perturbed runs share it so the miss stays common-mode.
_TERMINAL_MAX_STEP_S: float = 30.0


@dataclass(frozen=True)
class RunRecord:
    """One run's outcome: its inputs, the RTN miss, ToA error, and perigee."""

    inputs: RunInputs
    miss_rtn_m: Vec3
    toa_miss_s: float
    perigee_alt_m: float
    crossing_position_m: Vec3
    crossing_velocity_m_s: Vec3


@dataclass(frozen=True)
class EnsembleResult:
    """An ensemble's per-run records, aggregate statistics, and the nominal reference."""

    master_seed: int
    nominal_perigee_alt_m: float
    nominal_toa_s: float
    records: tuple[RunRecord, ...]
    stats: EnsembleStats


def replay_inputs(master_seed: int, spec: DispersionSpec, run_index: int) -> RunInputs:
    """Reconstruct a single run's draws standalone (§14.2), without the ensemble."""
    rng = np.random.default_rng(_child_seed(master_seed, run_index))
    return sample_run_inputs(rng, spec, run_index)


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


def _child_seed(master_seed: int, run_index: int) -> Any:
    # spawn_key=(i,) reproduces SeedSequence(master_seed).spawn(n)[i] standalone.
    return np.random.SeedSequence(entropy=master_seed, spawn_key=(run_index,))


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


def run_ensemble(
    spec: DispersionSpec,
    n: int,
    master_seed: int,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    control: Controller | None = None,
) -> EnsembleResult:
    """Run an open-loop dispersion ensemble and aggregate it.

    ``control`` is the §14.1 hook; it must be None here (the open-loop capstone).
    """
    if control is not None:
        raise NotImplementedError(
            "Closed-loop control is a Rung-D addition; capstone is open-loop."
        )

    earth = _earth()
    epoch = _to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)
    frame = FramesFactory.getEME2000()
    mu: float = Constants.WGS84_EARTH_MU

    # Nominal (unperturbed) crossing — the reference the miss is measured against.
    nominal_prop = build_propagator(orbital_config, presets.full_force(), _TERMINAL_MAX_STEP_S)
    nominal = _propagate_to_interception(nominal_prop, epoch, period, earth)
    nominal_basis: Basis = rtn_basis(nominal.position_m, nominal.velocity_m_s)

    # Apogee deployment state and its RTN basis (for injecting the Δv).
    apo_date, apo_pos, apo_vel = _apogee_state(orbital_config)
    apo_basis: Basis = rtn_basis(_vec3(apo_pos), _vec3(apo_vel))

    records: list[RunRecord] = []
    for i in range(n):
        rng = np.random.default_rng(_child_seed(master_seed, i))
        inputs = sample_run_inputs(rng, spec, i)

        dv_eme = rtn_to_cartesian(inputs.dv_rtn_m_s, apo_basis)
        perturbed_vel = apo_vel.add(Vector3D(dv_eme[0], dv_eme[1], dv_eme[2]))
        orbit = CartesianOrbit(
            TimeStampedPVCoordinates(apo_date, apo_pos, perturbed_vel), frame, mu
        )
        prop = build_propagator_from_orbit(orbit, physics_from_inputs(inputs), _TERMINAL_MAX_STEP_S)
        crossing = _propagate_to_interception(prop, epoch, period, earth)

        miss_vec: Vec3 = (
            crossing.position_m[0] - nominal.position_m[0],
            crossing.position_m[1] - nominal.position_m[1],
            crossing.position_m[2] - nominal.position_m[2],
        )
        records.append(
            RunRecord(
                inputs=inputs,
                miss_rtn_m=rtn_components(miss_vec, nominal_basis),
                toa_miss_s=crossing.toa_s - nominal.toa_s,
                perigee_alt_m=crossing.perigee_alt_m,
                crossing_position_m=crossing.position_m,
                crossing_velocity_m_s=crossing.velocity_m_s,
            )
        )

    stats = summarize(
        np.array([r.miss_rtn_m for r in records], dtype=np.float64).reshape(n, 3),
        np.array([r.toa_miss_s for r in records], dtype=np.float64),
        np.array([r.perigee_alt_m for r in records], dtype=np.float64),
    )
    return EnsembleResult(
        master_seed=master_seed,
        nominal_perigee_alt_m=nominal.perigee_alt_m,
        nominal_toa_s=nominal.toa_s,
        records=tuple(records),
        stats=stats,
    )


def format_summary(result: EnsembleResult) -> str:
    """Human-readable one-screen summary of an ensemble (the capstone report)."""
    s = result.stats
    mean = s.miss_rtn_mean_m
    std = s.miss_rtn_std_m
    return "\n".join(
        (
            f"Open-loop dispersion capstone — N={s.n}, master_seed={result.master_seed}",
            f"  Nominal: perigee {result.nominal_perigee_alt_m / 1e3:.1f} km,"
            f" coast {result.nominal_toa_s / 3600:.2f} h",
            "  Interception miss vs nominal, RTN frame [m] (T = dr_p/dv_a lever):",
            f"    bias R/T/N = {mean[0]:+.1f} / {mean[1]:+.1f} / {mean[2]:+.1f}",
            f"    std  R/T/N = {std[0]:.1f} / {std[1]:.1f} / {std[2]:.1f}",
            f"  Time-of-arrival miss: {s.toa_miss_mean_s:+.2f} ± {s.toa_miss_std_s:.2f} s",
            f"  Perigee (diagnostic, low=good): {s.perigee_alt_mean_m / 1e3:.1f}"
            f" ± {s.perigee_alt_std_m / 1e3:.2f} km"
            f" [min {s.perigee_alt_min_m / 1e3:.1f}, max {s.perigee_alt_max_m / 1e3:.1f}]",
        )
    )


def main() -> None:
    # Smoke-sized ensemble (design doc §10.4: N=50 is a smoke test, not the result;
    # the resolved-tail controllability result needs N=10³–10⁴, a longer job).
    result = run_ensemble(DispersionSpec(), n=50, master_seed=20260608)
    print(format_summary(result))


if __name__ == "__main__":
    main()
