"""Pure truth-model validation core — a Rung-D pre-gate (ADR 0018, no JVM in this module).

Rung D's verdict is only as trustworthy as the Orekit truth it rides on, and the coast is
~99 % of the trajectory — a truth-model bug (wrong frame, μ, J2 sign, a leaking integrator)
would hide there. This gate is the confirmation (non-blocking) that it does not, in two tiers:

* **Tier 1 — integrator health, on the conservative dynamics.** On a *numerical* two-body coast
  (the analytic Kepler route bypassed), the specific orbital energy ``v²/2 − μ/r`` and the
  angular-momentum magnitude ``|r×v|`` are constants of motion, so any drift over the arc is pure
  numerical leakage. Plus **tolerance-halving**: re-fly the same coast at a tenth the integrator
  relative tolerance — if the trajectory barely moves, truncation error sits far below the
  dispersion scale Rung D measures.

* **Tier 2 — an independent cross-check, on the dominant perturbed dynamics.** Re-fly the J2 coast
  with an *independent* pure-Python RK4 Cowell (`estimation.two_body_j2_flow`, the C1 onboard
  model — a separate implementation that shares only the pinned constants) and compare
  trajectories. Agreement confirms Orekit's force assembly / frame / integrator *setup* is right
  in the dynamics that dominate the coast; the non-conservative forces (drag/SRP/third-body) are
  validated separately by the Rung-A force-signature tests.

This module owns the pure invariants, the independent propagation, and the comparison reductions;
the live Orekit coasts are flown by the JVM glue in :mod:`puffsat_sim.runs.truth_validation`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.constants import WGS84_MU
from puffsat_sim.estimation import two_body_j2_flow

# Verdict thresholds (fractional / relative to the orbit scale), set against the measured
# drift on the reference coast — generous headroom over the floor the integrator actually hits.
CONSERVATION_TOL_FRAC: float = 1e-7
CONVERGENCE_TOL_FRAC: float = 1e-7
CROSSCHECK_TOL_FRAC: float = 1e-4


def specific_energy_j_per_kg(
    states: NDArray[np.float64], mu: float = WGS84_MU
) -> NDArray[np.float64]:
    """Two-body specific orbital energy ``v²/2 − μ/r`` per sample (a (N,6) state history)."""
    s = np.asarray(states, dtype=np.float64).reshape(-1, 6)
    r_mag = np.linalg.norm(s[:, :3], axis=1)
    v_sq = np.sum(s[:, 3:] ** 2, axis=1)
    return np.asarray(0.5 * v_sq - mu / r_mag, dtype=np.float64)


def angular_momentum_magnitude(states: NDArray[np.float64]) -> NDArray[np.float64]:
    """Specific angular-momentum magnitude ``|r×v|`` per sample."""
    s = np.asarray(states, dtype=np.float64).reshape(-1, 6)
    return np.asarray(np.linalg.norm(np.cross(s[:, :3], s[:, 3:]), axis=1), dtype=np.float64)


def max_fractional_drift(values: NDArray[np.float64]) -> float:
    """The largest deviation from the first sample, relative to it: ``max|v − v₀| / |v₀|``."""
    v = np.asarray(values, dtype=np.float64)
    return float(np.max(np.abs(v - v[0])) / abs(v[0]))


@dataclass(frozen=True)
class ConservationDrift:
    """Fractional drift of the two-body constants of motion over a coast arc."""

    energy_frac: float
    ang_mom_frac: float

    @property
    def worst(self) -> float:
        return max(self.energy_frac, self.ang_mom_frac)


def conservation_drift(states: NDArray[np.float64], mu: float = WGS84_MU) -> ConservationDrift:
    """Energy and angular-momentum fractional drift over a numerical two-body coast."""
    return ConservationDrift(
        energy_frac=max_fractional_drift(specific_energy_j_per_kg(states, mu)),
        ang_mom_frac=max_fractional_drift(angular_momentum_magnitude(states)),
    )


def independent_coast(
    state0: NDArray[np.float64], times_s: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Propagate ``state0`` through the sample epochs with the independent RK4 Cowell.

    Chained segment-by-segment from the single initial state (never re-seeded from Orekit),
    so the result is a genuinely independent trajectory to compare against the truth samples.
    """
    t = np.asarray(times_s, dtype=np.float64)
    states = [np.asarray(state0, dtype=np.float64)]
    for i in range(1, t.size):
        states.append(two_body_j2_flow(states[-1], float(t[i] - t[i - 1])))
    return np.asarray(states, dtype=np.float64)


def max_position_divergence_m(
    states_a: NDArray[np.float64], states_b: NDArray[np.float64]
) -> float:
    """Largest position separation between two state histories sampled at common epochs."""
    a = np.asarray(states_a, dtype=np.float64).reshape(-1, 6)
    b = np.asarray(states_b, dtype=np.float64).reshape(-1, 6)
    return float(np.max(np.linalg.norm(a[:, :3] - b[:, :3], axis=1)))


@dataclass(frozen=True)
class TruthValidationFinding:
    """The truth-model validation verdict: integrator health (Tier 1) + cross-check (Tier 2)."""

    conservation: ConservationDrift
    convergence_divergence_m: float
    crosscheck_divergence_m: float
    orbit_scale_m: float
    span_s: float
    n_samples: int
    conservation_tol_frac: float = CONSERVATION_TOL_FRAC
    convergence_tol_frac: float = CONVERGENCE_TOL_FRAC
    crosscheck_tol_frac: float = CROSSCHECK_TOL_FRAC

    @property
    def convergence_frac(self) -> float:
        return self.convergence_divergence_m / self.orbit_scale_m

    @property
    def crosscheck_frac(self) -> float:
        return self.crosscheck_divergence_m / self.orbit_scale_m

    @property
    def conservation_ok(self) -> bool:
        return self.conservation.worst <= self.conservation_tol_frac

    @property
    def convergence_ok(self) -> bool:
        return self.convergence_frac <= self.convergence_tol_frac

    @property
    def crosscheck_ok(self) -> bool:
        return self.crosscheck_frac <= self.crosscheck_tol_frac

    @property
    def validated(self) -> bool:
        return self.conservation_ok and self.convergence_ok and self.crosscheck_ok


def format_truth_validation(finding: TruthValidationFinding) -> str:
    """One-screen truth-validation report — Tier 1 health + Tier 2 cross-check (non-blocking)."""

    def verdict(ok: bool) -> str:
        return "OK" if ok else "DRIFT"

    lines = [
        "Truth-model validation — coast integrator health + independent cross-check"
        " (ADR 0018; non-blocking)",
        f"  Arc: {finding.span_s / 3600.0:.2f} h coast, {finding.n_samples} samples,"
        f" orbit scale {finding.orbit_scale_m / 1e3:.0f} km.",
        "  Tier 1 — numerical two-body coast (energy & angular momentum are constants of motion):",
        f"    energy drift {finding.conservation.energy_frac:.2e},"
        f" |h| drift {finding.conservation.ang_mom_frac:.2e}"
        f" vs {finding.conservation_tol_frac:.0e} — {verdict(finding.conservation_ok)}.",
        f"    tolerance-halving (×0.1 rel_tol): {finding.convergence_divergence_m:.3e} m"
        f" = {finding.convergence_frac:.2e} of the orbit scale"
        f" vs {finding.convergence_tol_frac:.0e} — {verdict(finding.convergence_ok)}.",
        "  Tier 2 — independent RK4 Cowell (two-body + J2) vs the Orekit J2 coast:",
        f"    max position divergence {finding.crosscheck_divergence_m:.3e} m"
        f" = {finding.crosscheck_frac:.2e} of the orbit scale"
        f" vs {finding.crosscheck_tol_frac:.0e} — {verdict(finding.crosscheck_ok)}.",
        f"  → {'VALIDATED' if finding.validated else 'NOT VALIDATED'}:"
        " the coast truth model conserves, converges, and matches an independent propagator.",
    ]
    return "\n".join(lines)
