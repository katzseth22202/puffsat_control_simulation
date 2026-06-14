"""Pure corrector-in-loop validation core — a Rung-D / D1.x sub-slice (ADR 0018, no JVM).

D1.1 did not run the midcourse corrector per unit: it *sampled* each PuffSat's hand-off lateral
entry offset from the linear C0/C1/C2a budget (``train.ENTRY_LATERAL_PERUNIT_M`` 141 m nav,
``ENTRY_LATERAL_SHARED_M`` 149 m Cr-prior) and flew only the C3b terminal loop.  That bundled two
corrector-side approximations the ADR 0018 D1.1 finding flagged as D1.x refinements:

* **Linear superposition / unbiasedness.** C0 (:mod:`puffsat_sim.navigation`) proved the corrector
  residual is linear in nav error (``miss = Φ·δ``), but only by sweeping *one axis at a time*.
  D1.1's sampled entry assumes that superposes across all six axes combined at realistic C1
  magnitudes — never tested with the real corrector on **combined stochastic draws**.
* **Crossing budget applied at the hand-off.** The 141 m is the lateral *crossing* (200 km)
  sensitivity, fed to the terminal loop as the *800 km hand-off* displacement — a conservative
  proxy, the true hand-off-lateral magnitude a D1.x refinement.

This module is the pure *reduction* for the brute-force validation batch (ADR 0018 decision 6,
"validated by a brute-force batch (A)"): it takes the JVM batch's per-draw nav errors, measured
interception misses, and measured hand-off offsets, and reduces them against the linear ``Φ``/``Σ``
prediction the sampled entry relies on.  The scope is the **nav leg** (the Cr-prior mismatch leg and
the Φ-Jacobian quasi-Newton speedup are separate follow-ons).  Like :mod:`puffsat_sim.navigation`
the math is pure (the JVM glue lives in :mod:`puffsat_sim.runs.corrector_validation`).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.navigation import induced_miss_covariance

# The interception miss is in the nominal-crossing RTN frame; its radial component is pinned ~0 by
# the 200 km altitude-event crossing, so "lateral" is the (Transverse, Normal) pair — C0's
# convention (navigation._LATERAL_AXES / nav_feasibility's √(induced[1,1]+induced[2,2])).
_LATERAL_AXES: tuple[int, int] = (1, 2)

# Irreducible model-tolerance floor on how close the measured crossing-miss σ must sit to the linear
# ΦΣΦᵀ prediction (Φ from a minimal sweep, the diagonal-Σ sampling approximation).
SIGMA_CONSISTENCY_TOL: float = 0.15
# An RMS estimate from N samples has relative standard error ≈ 1/√(2N); the consistency band is the
# looser of the model floor and this 3σ sampling band, so a smoke-sized batch is judged honestly
# (a flat tolerance would false-fail small N on pure sampling noise).
_SIGMA_SAMPLING_K: float = 3.0


def _lateral_rms_m(lateral: NDArray[np.float64]) -> float:
    """RMS lateral miss magnitude ``√⟨|d_⊥|²⟩`` for an ``(N, 2)`` ⊥v miss set.

    Matches C0/C1's ``√(σ_T² + σ_N²)`` convention (``nav_feasibility.lateral_miss_1sigma_m``), so a
    measured RMS lateral lines up directly with ``√tr(ΦΣΦᵀ)_TN`` and ``ENTRY_LATERAL_PERUNIT_M``.
    """
    return math.sqrt(float(np.mean(np.sum(lateral**2, axis=1)))) if lateral.size else 0.0


@dataclass(frozen=True)
class CorrectorValidationFinding:
    """The D1.x verdict: does the real corrector reproduce D1.1's linear sampled entry?"""

    n_samples: int
    linearity_residual_rel: float
    bias_m: tuple[float, float, float]
    measured_crossing_lateral_sigma_m: float
    predicted_crossing_lateral_sigma_m: float
    measured_handoff_lateral_sigma_m: float
    crossing_proxy_m: float
    rel_tol: float

    @property
    def bias_lateral_m(self) -> float:
        """The lateral (⊥v) component of the mean interception miss — should vanish in N."""
        return float(math.hypot(self.bias_m[_LATERAL_AXES[0]], self.bias_m[_LATERAL_AXES[1]]))

    @property
    def linear(self) -> bool:
        """The real corrector's per-draw miss tracks ``Φ·δ`` — combined-draw superposition holds."""
        return self.linearity_residual_rel < self.rel_tol

    @property
    def unbiased(self) -> bool:
        """The mean lateral miss sits within the sample-mean noise of zero.

        A *systematic* bias is caught more sharply by :attr:`linear` (a constant offset makes
        ``miss ≠ Φ·δ`` on every draw, inflating the residual); this is the belt-and-suspenders read
        on the sample mean, at a generous 4σ/√N so a clean batch never false-fails on tail noise.
        """
        standard_error = 4.0 * self.measured_crossing_lateral_sigma_m / math.sqrt(self.n_samples)
        return self.bias_lateral_m <= standard_error

    @property
    def crossing_sigma_consistent(self) -> bool:
        """The measured crossing-miss σ matches the linear ΦΣΦᵀ prediction (the 141 m end check).

        Judged against the looser of the model floor and the 3σ RMS-sampling band ``3/√(2N)``, so a
        smoke-sized batch is not failed for sampling noise the population prediction cannot have.
        """
        if self.predicted_crossing_lateral_sigma_m == 0.0:
            return self.measured_crossing_lateral_sigma_m == 0.0
        band = max(SIGMA_CONSISTENCY_TOL, _SIGMA_SAMPLING_K / math.sqrt(2 * self.n_samples))
        rel = abs(self.measured_crossing_lateral_sigma_m - self.predicted_crossing_lateral_sigma_m)
        return rel / self.predicted_crossing_lateral_sigma_m < band

    @property
    def conservatism_factor(self) -> float:
        """How much larger the fed crossing-proxy entry is than the true hand-off displacement.

        ``> 1`` means feeding the fully-developed crossing miss (141 m) as the hand-off offset
        over-stresses the terminal loop relative to the smaller actual hand-off displacement — so
        D1.1's proxy is conservative (its verdict is pessimistic, the safe direction).
        """
        if self.measured_handoff_lateral_sigma_m == 0.0:
            return math.inf
        return self.measured_crossing_lateral_sigma_m / self.measured_handoff_lateral_sigma_m

    @property
    def handoff_conservative(self) -> bool:
        return self.measured_handoff_lateral_sigma_m <= self.crossing_proxy_m

    @property
    def validated(self) -> bool:
        """D1.1's sampled-entry shortcut is confirmed: linear, unbiased, consistent, safe."""
        return (
            self.linear
            and self.unbiased
            and self.crossing_sigma_consistent
            and self.handoff_conservative
        )


def summarize_corrector_validation(
    nav_draws: NDArray[np.float64],
    crossing_misses_m: NDArray[np.float64],
    handoff_laterals_m: NDArray[np.float64],
    phi: NDArray[np.float64],
    sigma6: NDArray[np.float64],
    crossing_proxy_m: float,
    rel_tol: float = 0.05,
) -> CorrectorValidationFinding:
    """Reduce a corrector-in-loop batch against the linear Φ/Σ prediction (ADR 0018 decision 6).

    ``nav_draws`` is ``(N, 6)`` apogee-RTN nav errors; ``crossing_misses_m`` the ``(N, 3)`` measured
    interception misses (nominal-crossing RTN); ``handoff_laterals_m`` the ``(N, 2)`` measured ⊥v
    hand-off offsets.  ``phi`` is C0's 3×6 sensitivity and ``sigma6`` the apogee-RTN Σ the draws
    came from; ``crossing_proxy_m`` is D1.1's fed entry (``ENTRY_LATERAL_PERUNIT_M``).
    """
    draws = np.asarray(nav_draws, dtype=np.float64).reshape(-1, 6)
    measured = np.asarray(crossing_misses_m, dtype=np.float64).reshape(-1, 3)
    handoff = np.asarray(handoff_laterals_m, dtype=np.float64).reshape(-1, 2)
    n = measured.shape[0]

    predicted = draws @ np.asarray(phi, dtype=np.float64).T
    residual_rms = math.sqrt(float(np.mean(np.sum((measured - predicted) ** 2, axis=1))))
    measured_rms = math.sqrt(float(np.mean(np.sum(measured**2, axis=1))))
    linearity_residual_rel = residual_rms / measured_rms if measured_rms > 0.0 else 0.0

    bias = measured.mean(axis=0)
    induced = induced_miss_covariance(phi, sigma6)
    predicted_lateral = math.sqrt(
        float(
            induced[_LATERAL_AXES[0], _LATERAL_AXES[0]]
            + induced[_LATERAL_AXES[1], _LATERAL_AXES[1]]
        )
    )

    return CorrectorValidationFinding(
        n_samples=n,
        linearity_residual_rel=linearity_residual_rel,
        bias_m=(float(bias[0]), float(bias[1]), float(bias[2])),
        measured_crossing_lateral_sigma_m=_lateral_rms_m(measured[:, list(_LATERAL_AXES)]),
        predicted_crossing_lateral_sigma_m=predicted_lateral,
        measured_handoff_lateral_sigma_m=_lateral_rms_m(handoff),
        crossing_proxy_m=crossing_proxy_m,
        rel_tol=rel_tol,
    )


def format_corrector_validation(finding: CorrectorValidationFinding) -> str:
    """One-screen D1.x corrector-in-loop validation report."""
    verdict = "VALIDATED" if finding.validated else "NOT validated"
    lines = [
        "Rung D / D1.x — corrector-in-loop validation (nav leg; ADR 0018 decision 6)",
        f"  Real corrector over N={finding.n_samples} combined nav draws (C1 Σ); the C3b loop's"
        " sampled entry stands in for this.",
        f"  Linearity: per-draw |miss − Φ·δ| = {finding.linearity_residual_rel * 100:.2f}%"
        f" of |miss| (< {finding.rel_tol * 100:g}%"
        f" ⇒ {'linear' if finding.linear else 'NONLINEAR'});"
        f" lateral bias {finding.bias_lateral_m:.2f} m"
        f" — {'unbiased' if finding.unbiased else 'BIASED'}.",
        f"  Crossing miss σ: measured {finding.measured_crossing_lateral_sigma_m:.1f} m vs"
        f" ΦΣΦᵀ {finding.predicted_crossing_lateral_sigma_m:.1f} m vs D1.1 proxy"
        f" {finding.crossing_proxy_m:.0f} m"
        f" — {'consistent' if finding.crossing_sigma_consistent else 'INCONSISTENT'}"
        " (the per-unit entry magnitude, end-to-end).",
        f"  Hand-off displacement σ {finding.measured_handoff_lateral_sigma_m:.1f} m vs the"
        f" {finding.crossing_proxy_m:.0f} m crossing proxy → {finding.conservatism_factor:.2f}×"
        f" conservative — feeding the crossing miss at the hand-off"
        f" {'over-stresses' if finding.handoff_conservative else 'UNDER-stresses'} the loop.",
        f"  Verdict: D1.1's Φ-composed sampled entry is {verdict}"
        " (linear, unbiased, σ-consistent, conservative).",
    ]
    return "\n".join(lines)
