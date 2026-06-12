"""C1 sensing feasibility — the JVM run for :mod:`puffsat_sim.nav_feasibility` (ADR 0012)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.frames import FramesFactory
from org.orekit.orbits import CartesianOrbit
from org.orekit.utils import Constants, TimeStampedPVCoordinates

from puffsat_sim import mission
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.control import report_controller
from puffsat_sim.descent import COAST_MAX_STEP_S, apogee_state
from puffsat_sim.montecarlo import nominal_inputs, physics_from_inputs
from puffsat_sim.nav_feasibility import (
    NavFeasibilityResult,
    NavFeasibilitySpec,
    NavValidationOutcome,
    format_nav_feasibility,
    format_nav_validation,
    sweep_nav_feasibility,
    validate_cell,
)
from puffsat_sim.navigation import NavSweepSpec, summarize_nav_requirement
from puffsat_sim.propagator import build_propagator_from_orbit
from puffsat_sim.runs.navigation import run_nav_sweep


def truth_arc_to_apogee(
    cadence_hz: float,
    arc_duration_s: float,
    physics: PhysicsConfig,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> NDArray[np.float64]:
    """Full-force truth states at each measurement epoch of the C1 coast arc (ADR 0012).

    The arc *ends at the apogee deployment state*: the propagator runs backward
    ``arc_duration_s`` from apogee, then forward through the cadence epochs — so
    the final row is the apogee correction node where Σ is judged.  Rows are
    EME2000 ``[x y z vx vy vz]``, the arc start plus one row per epoch — the
    ``truth_states`` contract of :func:`puffsat_sim.nav_feasibility.validate_cell`.
    """
    apo_date, apo_pos, apo_vel = apogee_state(orbital_config)
    frame = FramesFactory.getEME2000()
    orbit = CartesianOrbit(
        TimeStampedPVCoordinates(apo_date, apo_pos, apo_vel),
        frame,
        float(Constants.WGS84_EARTH_MU),
    )
    propagator = build_propagator_from_orbit(orbit, physics, COAST_MAX_STEP_S)

    dt_s = 1.0 / cadence_hz
    n_epochs = max(1, round(arc_duration_s * cadence_hz))
    states = np.empty((n_epochs + 1, 6), dtype=np.float64)
    for k in range(n_epochs + 1):
        date = apo_date.shiftedBy(float(k * dt_s - n_epochs * dt_s))
        pv = propagator.propagate(date).getPVCoordinates()
        p, v = pv.getPosition(), pv.getVelocity()
        states[k] = (p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ())
    return states


def run_nav_feasibility(
    spec: NavFeasibilitySpec,
    phi: NDArray[np.float64],
    seed: int,
    validation_indices: tuple[int, ...] | None = None,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> tuple[NavFeasibilityResult, tuple[NavValidationOutcome, ...]]:
    """The C1 JVM seam: pure LinCov sweep + seeded UKF truth runs judged by NEES.

    Layer 1 (the envelope) is entirely pure — only the layer-2 validation touches
    the JVM, generating the full-force truth arc the filter's two-body+J2 model
    must stay honest against (the model gap q absorbs, ADR 0012 decisions 6/7).
    ``phi`` is C0's measured sensitivity; ``validation_indices`` defaults to the
    nominal cell plus the range-only cell (the two headline claims).  Truth arcs
    are cached per cadence — they depend on sampling, not on the cell's noises.
    """
    result = sweep_nav_feasibility(spec, phi)
    if validation_indices is None:
        range_only = tuple(
            o.cell.cell_index for o in result.outcomes if o.cell.doppler_sigma_m_s is None
        )
        validation_indices = (0,) + range_only

    physics = physics_from_inputs(nominal_inputs(0))
    arcs: dict[float, NDArray[np.float64]] = {}
    validations = []
    for index in validation_indices:
        cell = result.outcomes[index].cell
        if cell.cadence_hz not in arcs:
            arcs[cell.cadence_hz] = truth_arc_to_apogee(
                cell.cadence_hz, spec.arc_duration_s, physics, orbital_config
            )
        validations.append(validate_cell(cell, spec, arcs[cell.cadence_hz], seed))
    return result, tuple(validations)


def nav_feasibility_report(
    spec: NavFeasibilitySpec | None = None,
    seed: int = 20260610,
    nav_spec: NavSweepSpec | None = None,
) -> str:
    """Run the full C1 report: measured Φ → LinCov envelope → NEES validation (ADR 0012).

    Φ is re-derived live from a minimal C0 sweep (``points_per_sign=1`` — C0
    measured the residual linear across the full swept range, so the smallest
    ± pair per axis already determines every column).
    """
    nav_result = run_nav_sweep(nav_spec or NavSweepSpec(points_per_sign=1), report_controller)
    requirement = summarize_nav_requirement(nav_result)
    result, validations = run_nav_feasibility(spec or NavFeasibilitySpec(), requirement.phi, seed)
    return format_nav_feasibility(result) + "\n" + format_nav_validation(validations)
