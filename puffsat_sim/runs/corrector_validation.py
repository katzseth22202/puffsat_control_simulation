"""D1.x corrector-in-loop validation — the JVM run for :mod:`puffsat_sim.corrector_validation`.

Runs the real midcourse corrector (the C0 path: :func:`montecarlo.run_record` with
``control=report_controller`` and a predict-side nav offset) over a brute-force batch of
**combined** nav draws sampled from the C1 apogee Σ, and measures two things D1.1 only approximated
(ADR 0018 decision 6):

* the interception **crossing miss** (does the real corrector reproduce the linear ``Φ·δ`` that
  D1.1's sampled entry assumes, with all six axes perturbed at once?), and
* the actual **800 km hand-off lateral displacement** of the corrected trajectory (vs the crossing
  budget D1.1 fed the terminal loop as a conservative proxy).

Nav leg only (the C2a Cr-prior mismatch leg and the Φ-Jacobian quasi-Newton speedup are separate
follow-ons).  Φ comes from a minimal C0 sweep and Σ from the C1 nominal cell — both reused, not
rebuilt; the pure reduction lives in :mod:`puffsat_sim.corrector_validation`.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.orbits import CartesianOrbit
from org.orekit.utils import TimeStampedPVCoordinates

from org.hipparchus.geometry.euclidean.threed import Vector3D

from puffsat_sim import mission
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.control import report_controller
from puffsat_sim.corrector_validation import (
    CorrectorValidationFinding,
    format_corrector_validation,
    summarize_corrector_validation,
)
from puffsat_sim.descent import COAST_MAX_STEP_S, coast_to_handoff, vec3
from puffsat_sim.dispersion import Basis, Vec3, rtn_basis, rtn_components, rtn_to_cartesian
from puffsat_sim.montecarlo import (
    RunContext,
    RunVariant,
    build_context,
    nominal_inputs,
    physics_from_inputs,
    run_record,
)
from puffsat_sim.nav_feasibility import NavFeasibilitySpec, evaluate_cell, nav_feasibility_cells
from puffsat_sim.navigation import NavSweepSpec, sample_nav_error, summarize_nav_requirement
from puffsat_sim.propagator import build_propagator_from_orbit
from puffsat_sim.records import RunRecord
from puffsat_sim.runs.navigation import run_nav_sweep
from puffsat_sim.train import ENTRY_LATERAL_PERUNIT_M

VALIDATION_MASTER_SEED: int = 20260614
SMOKE_N: int = 24


def _c1_nominal_sigma6(phi: NDArray[np.float64]) -> NDArray[np.float64]:
    """The C1 nominal-cell apogee-RTN nav covariance, as a diagonal Σ from its per-axis σ.

    LinCov is pure given Φ, so this needs no extra JVM.  The full 6×6 Σ is computed inside
    ``evaluate_cell`` but not exposed; the per-draw linearity check is correlation-insensitive, so
    the diagonal at the honest per-axis magnitudes is the right sampling distribution here.
    """
    spec = NavFeasibilitySpec()
    outcome = evaluate_cell(nav_feasibility_cells(spec)[0], spec, phi)
    variances = [s**2 for s in outcome.pos_sigma_rtn_m] + [s**2 for s in outcome.vel_sigma_rtn_m_s]
    return np.diag(np.array(variances, dtype=np.float64))


def _coast_to_handoff_offset(
    ctx: RunContext,
    correction_rtn: Vec3,
    physics: PhysicsConfig,
    nominal_handoff_pos: Vec3,
    nominal_handoff_basis: Basis,
) -> tuple[float, float]:
    """The corrected execute trajectory's ⊥v (T, N) lateral offset at the 800 km hand-off.

    Builds the same corrected apogee state ``run_record`` executes (apogee velocity + the solved
    correction; nominal injection is zero) and coasts it to the hand-off, returning the lateral
    displacement vs the nominal hand-off state — the funnel-entry the terminal loop would actually
    see, against D1.1's crossing-budget proxy.
    """
    corr_eme = rtn_to_cartesian(correction_rtn, ctx.apo_basis)
    vel = ctx.apo_vel.add(Vector3D(corr_eme[0], corr_eme[1], corr_eme[2]))
    orbit = CartesianOrbit(
        TimeStampedPVCoordinates(ctx.apo_date, ctx.apo_pos, vel), ctx.frame, ctx.mu
    )
    coast = build_propagator_from_orbit(orbit, physics, COAST_MAX_STEP_S)
    pos = vec3(
        coast_to_handoff(coast, ctx.epoch, ctx.period, ctx.earth).getPVCoordinates().getPosition()
    )
    offset: Vec3 = (
        pos[0] - nominal_handoff_pos[0],
        pos[1] - nominal_handoff_pos[1],
        pos[2] - nominal_handoff_pos[2],
    )
    _, t, n = rtn_components(offset, nominal_handoff_basis)
    return (t, n)


def _nominal_handoff(ctx: RunContext, physics: PhysicsConfig) -> tuple[Vec3, Basis]:
    """Coast the unperturbed apogee state to the 800 km hand-off; its position + RTN basis."""
    orbit = CartesianOrbit(
        TimeStampedPVCoordinates(ctx.apo_date, ctx.apo_pos, ctx.apo_vel), ctx.frame, ctx.mu
    )
    coast = build_propagator_from_orbit(orbit, physics, COAST_MAX_STEP_S)
    pv = coast_to_handoff(coast, ctx.epoch, ctx.period, ctx.earth).getPVCoordinates()
    pos, vel = vec3(pv.getPosition()), vec3(pv.getVelocity())
    return pos, rtn_basis(pos, vel)


def _correction_rtn(record: RunRecord) -> Vec3:
    return record.control_log[0].dv_rtn_m_s if record.control_log else (0.0, 0.0, 0.0)


def run_corrector_validation(
    n: int = SMOKE_N,
    master_seed: int = VALIDATION_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> CorrectorValidationFinding:
    """Fly the real corrector over ``n`` combined C1-Σ nav draws and reduce vs the linear Φ/Σ map.

    Each draw runs the corrector from the perturbed estimate and executes against truth (nominal
    coefficients, zero injection — only the nav offset diverges predict from execute, C0's regime),
    recording the crossing miss + the corrected hand-off offset.
    """
    ctx = build_context(orbital_config)
    phi = summarize_nav_requirement(
        run_nav_sweep(NavSweepSpec(points_per_sign=1), report_controller)
    ).phi
    sigma6 = _c1_nominal_sigma6(phi)
    physics = physics_from_inputs(nominal_inputs(0))
    nominal_handoff_pos, nominal_handoff_basis = _nominal_handoff(ctx, physics)

    draws: list[tuple[float, ...]] = []
    crossing_misses: list[Vec3] = []
    handoff_laterals: list[tuple[float, float]] = []
    for i in range(n):
        delta = sample_nav_error(np.random.default_rng((master_seed, i)), sigma6)
        record = run_record(
            ctx, nominal_inputs(i), RunVariant(control=report_controller, nav_offset_rtn6=delta)
        )
        draws.append(delta)
        crossing_misses.append(record.miss_rtn_m)
        handoff_laterals.append(
            _coast_to_handoff_offset(
                ctx, _correction_rtn(record), physics, nominal_handoff_pos, nominal_handoff_basis
            )
        )

    return summarize_corrector_validation(
        np.array(draws, dtype=np.float64),
        np.array(crossing_misses, dtype=np.float64),
        np.array(handoff_laterals, dtype=np.float64),
        phi,
        sigma6,
        ENTRY_LATERAL_PERUNIT_M,
    )


def corrector_validation_report(
    n: int = SMOKE_N,
    master_seed: int = VALIDATION_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Run the brute-force corrector-in-loop validation and format the D1.x finding."""
    return format_corrector_validation(run_corrector_validation(n, master_seed, orbital_config))


if __name__ == "__main__":
    print(corrector_validation_report())
