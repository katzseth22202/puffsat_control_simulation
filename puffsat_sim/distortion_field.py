"""Pure distortion-field study — the terminal cross-track hedge (no JVM).

The σ_θ tracker budget (:mod:`puffsat_sim.tracker_budget`) treats the focal-plane distortion as
a **scalar** 3 µrad floor — the load-bearing bench-calibration residual that the whole terminal
cross-track-to-target accuracy (σ_θ·R vs the 5 m plate) rides on.  D1 is already *feasible past
the noise knee* (the tail is entry/authority-limited, ADR 0021), so this module is **not** closing
a capture gap: it is **defense-in-depth on that 3 µrad number being optimistic**, plus extra
capture margin (the hedge program, ``todos/improve_terminal_crosstrack.md``).

The single idea: model distortion as a spatial *field* (an RMS amplitude and a correlation
structure) rather than a scalar floor.  Two parameters of that field answer two long-standing
questions — the same field read two ways — and the deliverable is a *sensitivity curve* plus a
**bench-characterization requirement**, honest about the fact that there is no bench data (this is
pure simulation: the spectrum is an explicit swept input, not a measured fact).

* **Readout A — per-detector differential gain (the "Gaia trick").**  Measure the beacon's bearing
  *relative to reference stars in the same FOV/exposure* (inertial directions known to µas) instead
  of absolutely on the focal plane.  The distortion common to beacon and nearby star cancels;
  only the *gradient* over the target-to-star separation survives.  This is genuine **uncredited**
  headroom — today the differential is invoked only to cancel bus *smear* and to justify
  cross-detector independence, never to cut the distortion floor itself.  The payoff is
  **spectrum-contingent and sign-flips**: a *smooth* (long correlation length) field leaves a
  gradient ≪ the floor → a several-× win; a *rough* (pixel-scale) field makes differencing two
  uncorrelated errors √2 *worse*.  The break-even correlation length is the headline, and it
  *outputs a bench-characterization requirement* (measure the distortion correlation length; bank
  the gain only if it exceeds the break-even).

* **Readout B — array common-mode tolerance (audit of ADR 0019).**  The fusion gate
  (:mod:`puffsat_sim.tracker_fusion`) banks the *entire* 3 µrad distortion in the **independent**
  bucket — zero assumed common-mode distortion — so the 1.62 µrad 5-array grade is the *optimistic*
  end.  The audit varies the detector-to-detector distortion correlation ρ and asks how much
  correlated distortion the array tolerates before its grade exceeds the D1.1 capture grade.  (The
  mechanism that *earns* the √N credit is Readout A: each detector differencing against its **own**
  stars leaves only its own optics-specific gradient, which is what makes the residual independent.)

Only the field's second-order statistics enter every σ here, so the correlation function — not a
sampled realization — fixes the result exactly; the module is analytic and math-only (a literal
random-field draw is a test-side cross-check of these formulas, not a production dependency).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from puffsat_sim.tracker_budget import (
    STAR_DENSITY_PER_SR,
    TrackerHardware,
    gyro_bridge_sigma_theta_rad,
)
from puffsat_sim.tracker_fusion import (
    ARRAY_N_DETECTORS,
    D1_CAPTURE_GRADE_SIGMA_THETA_RAD,
    SMEAR_COMMON_SIGMA_RAD,
)

_REFERENCE_HW = TrackerHardware()

# The absolute bench-calibration residual the σ_θ budget treats as a scalar floor — here the RMS
# amplitude of the distortion modeled as a spatial field.
DISTORTION_RMS_RAD: float = _REFERENCE_HW.distortion_floor_rad

# The per-detector independent term (besides distortion) that still averages with √N: the gyro
# inter-frame bridge.  The photon term (~0.01 µrad) is negligible and omitted from the audit.
DETECTOR_GYRO_SIGMA_RAD: float = gyro_bridge_sigma_theta_rad(
    _REFERENCE_HW.gyro_arw_rad_sqrt_s, _REFERENCE_HW.frame_interval_s
)

# Readout-A sweep: distortion correlation lengths from pixel-scale (rough) to several FOV (smooth).
DEFAULT_CORRELATION_LENGTHS_RAD: tuple[float, ...] = (
    0.3e-3,
    1.0e-3,
    2.0e-3,
    3.0e-3,
    5.0e-3,
    10.0e-3,
    30.0e-3,
)

# Readout-B sweep: the detector-to-detector distortion correlation fraction, 0 (the ADR 0019
# optimistic split) to 1 (fully common-mode — no √N gain on the distortion term).
DEFAULT_CORRELATION_FRACTIONS: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)


def nearest_star_separation_rad(star_density_per_sr: float = STAR_DENSITY_PER_SR) -> float:
    """Mean nearest-neighbor reference-star separation: ``1/(2√n)`` for a 2-D Poisson field.

    The differential measurement references the beacon against its *nearest* reference star, so this
    separation is the lag over which the distortion gradient acts.  A multi-star plate fit only
    shortens the effective lag, so the nearest-star reading is the conservative one.
    """
    return 1.0 / (2.0 * math.sqrt(star_density_per_sr))


def gaussian_correlation(separation_rad: float, correlation_length_rad: float) -> float:
    """Distortion spatial correlation ``ρ(Δθ) = exp(−Δθ²/2L²)`` at lag ``separation_rad``.

    A Gaussian autocorrelation: a *smooth* (low-spatial-frequency) field has a long correlation
    length ``L``, a *pixel-scale* (high-frequency) field a short one.  The differential variance
    depends on the field only through this second-order statistic, so the shape is the model.
    """
    return math.exp(-(separation_rad**2) / (2.0 * correlation_length_rad**2))


def differential_residual_rad(
    distortion_rms_rad: float, separation_rad: float, correlation_length_rad: float
) -> float:
    """Differential bearing residual ``σ_d·√(2(1−ρ))`` — beacon minus nearest reference star.

    ``Var(d_beacon − d_star) = 2σ_d²(1 − ρ(Δθ))``: the distortion *common* to beacon and star
    cancels (the Gaia trick), leaving only what de-correlates over their separation.  ρ→1 (smooth)
    drives it to zero; ρ→0 (rough) leaves ``√2·σ_d`` — *worse* than the absolute floor, because
    differencing two uncorrelated errors adds their variances.
    """
    rho = gaussian_correlation(separation_rad, correlation_length_rad)
    return distortion_rms_rad * math.sqrt(2.0 * (1.0 - rho))


def breakeven_correlation_length_rad(separation_rad: float) -> float:
    """Correlation length where the differential just equals the absolute floor (``ρ = 0.5``).

    ``σ_diff = σ_d`` ⇔ ``2(1−ρ) = 1`` ⇔ ``ρ = 0.5``; for the Gaussian ρ that is
    ``L = Δθ/√(2 ln 2)``.  Above this length the differential *cuts* the floor; below it the
    differential *backfires* — this length is the bench-characterization requirement.
    """
    return separation_rad / math.sqrt(2.0 * math.log(2.0))


@dataclass(frozen=True)
class DifferentialPoint:
    """One point on the Readout-A sensitivity curve: a correlation length and its diff gain."""

    correlation_length_rad: float
    separation_rad: float
    distortion_rms_rad: float

    @property
    def correlation(self) -> float:
        return gaussian_correlation(self.separation_rad, self.correlation_length_rad)

    @property
    def residual_rad(self) -> float:
        return differential_residual_rad(
            self.distortion_rms_rad, self.separation_rad, self.correlation_length_rad
        )

    @property
    def gain(self) -> float:
        """Absolute floor ÷ differential residual: >1 the Gaia trick wins, <1 it backfires."""
        return self.distortion_rms_rad / self.residual_rad

    @property
    def improves(self) -> bool:
        return self.residual_rad < self.distortion_rms_rad


def differential_curve(
    correlation_lengths_rad: Sequence[float] = DEFAULT_CORRELATION_LENGTHS_RAD,
    *,
    separation_rad: float | None = None,
    distortion_rms_rad: float = DISTORTION_RMS_RAD,
) -> tuple[DifferentialPoint, ...]:
    """The Readout-A residual-vs-correlation-length curve (separation defaults to nearest-star)."""
    sep = separation_rad if separation_rad is not None else nearest_star_separation_rad()
    return tuple(
        DifferentialPoint(length, sep, distortion_rms_rad) for length in correlation_lengths_rad
    )


@dataclass(frozen=True)
class DetectorBudget:
    """The per-detector σ_θ split for the array audit — only distortion carries a correlation knob.

    The smear residual is the rigid-body **common** floor (never averages across detectors); the
    gyro bridge is **independent** per detector (averages with √N); the distortion is the term whose
    detector-to-detector correlation the audit varies.  At zero distortion correlation this
    reproduces the ADR 0019 independent/common split (the optimistic 1.62 µrad 5-array grade).
    """

    sigma_distortion_rad: float = DISTORTION_RMS_RAD
    sigma_indep_other_rad: float = DETECTOR_GYRO_SIGMA_RAD
    sigma_common_rad: float = SMEAR_COMMON_SIGMA_RAD


def array_sigma_theta_rad(
    budget: DetectorBudget, n_detectors: int, distortion_correlation: float
) -> float:
    """Array σ_θ when the distortion is correlated detector-to-detector by ``ρ`` (equicorrelated).

    The mean of ``N`` equicorrelated errors has variance ``σ²·(1+(N−1)ρ)/N``: ρ=0 → ``σ²/N`` (full
    √N gain), ρ=1 → ``σ²`` (no gain — the array collapses to a single detector on this term).  The
    smear floor never averages; the gyro term averages fully.  So
    ``σ_array² = σ_smear² + σ_other²/N + σ_dist²·(1+(N−1)ρ)/N``.
    """
    n = float(n_detectors)
    dist_var = budget.sigma_distortion_rad**2 * (1.0 + (n - 1.0) * distortion_correlation) / n
    return math.sqrt(budget.sigma_common_rad**2 + budget.sigma_indep_other_rad**2 / n + dist_var)


def common_mode_tolerance(
    budget: DetectorBudget, n_detectors: int, capture_grade_sigma_theta_rad: float
) -> float:
    """The distortion correlation fraction ρ at which the array grade reaches the capture grade.

    ``σ_array²(ρ)`` is linear in ρ (= ``C + Bρ``), so ``ρ_tol = (grade² − C)/B``.  A value ≥ 1 means
    the array stays capture-grade across the *entire* correlation range — the verdict is robust to
    the zero-common-mode assumption, because the worst case (ρ=1) collapses to ~the single-detector
    grade, which is itself at the capture threshold.
    """
    n = float(n_detectors)
    c = (
        budget.sigma_common_rad**2
        + budget.sigma_indep_other_rad**2 / n
        + budget.sigma_distortion_rad**2 / n
    )
    b = budget.sigma_distortion_rad**2 * (n - 1.0) / n
    return (capture_grade_sigma_theta_rad**2 - c) / b


@dataclass(frozen=True)
class ArrayTolerancePoint:
    """One point on the Readout-B curve: a distortion correlation fraction and the array grade."""

    distortion_correlation: float
    n_detectors: int
    budget: DetectorBudget
    capture_grade_sigma_theta_rad: float = D1_CAPTURE_GRADE_SIGMA_THETA_RAD

    @property
    def sigma_theta_rad(self) -> float:
        return array_sigma_theta_rad(self.budget, self.n_detectors, self.distortion_correlation)

    @property
    def meets_capture_grade(self) -> bool:
        return self.sigma_theta_rad <= self.capture_grade_sigma_theta_rad


def tolerance_curve(
    correlation_fractions: Sequence[float] = DEFAULT_CORRELATION_FRACTIONS,
    *,
    n_detectors: int = ARRAY_N_DETECTORS,
    budget: DetectorBudget | None = None,
    capture_grade_sigma_theta_rad: float = D1_CAPTURE_GRADE_SIGMA_THETA_RAD,
) -> tuple[ArrayTolerancePoint, ...]:
    """The Readout-B array-grade-vs-distortion-correlation curve."""
    b = budget if budget is not None else DetectorBudget()
    return tuple(
        ArrayTolerancePoint(rho, n_detectors, b, capture_grade_sigma_theta_rad)
        for rho in correlation_fractions
    )


@dataclass(frozen=True)
class DistortionFieldFinding:
    """The distortion-field hedge: Readout A (differential gain) + Readout B (array tolerance)."""

    differential_points: tuple[DifferentialPoint, ...]
    separation_rad: float
    breakeven_length_rad: float
    distortion_rms_rad: float
    tolerance_points: tuple[ArrayTolerancePoint, ...]
    budget: DetectorBudget
    n_detectors: int
    capture_grade_sigma_theta_rad: float = D1_CAPTURE_GRADE_SIGMA_THETA_RAD

    @property
    def best_case_independent_sigma_rad(self) -> float:
        """The optimistic array grade — zero common-mode distortion (the ADR 0019 banked value)."""
        return array_sigma_theta_rad(self.budget, self.n_detectors, 0.0)

    @property
    def worst_case_correlated_sigma_rad(self) -> float:
        """The pessimistic grade — fully common-mode distortion (no √N gain on the floor term)."""
        return array_sigma_theta_rad(self.budget, self.n_detectors, 1.0)

    @property
    def common_mode_tolerance(self) -> float:
        return common_mode_tolerance(
            self.budget, self.n_detectors, self.capture_grade_sigma_theta_rad
        )

    @property
    def common_mode_robust(self) -> bool:
        """Capture-grade across the full [0, 1] distortion-correlation range (tolerance ≥ 1)."""
        return self.common_mode_tolerance >= 1.0

    @property
    def best_case_margin(self) -> float:
        return self.capture_grade_sigma_theta_rad / self.best_case_independent_sigma_rad

    @property
    def worst_case_margin(self) -> float:
        return self.capture_grade_sigma_theta_rad / self.worst_case_correlated_sigma_rad


def distortion_field_finding(
    *,
    correlation_lengths_rad: Sequence[float] = DEFAULT_CORRELATION_LENGTHS_RAD,
    correlation_fractions: Sequence[float] = DEFAULT_CORRELATION_FRACTIONS,
    separation_rad: float | None = None,
    distortion_rms_rad: float = DISTORTION_RMS_RAD,
    n_detectors: int = ARRAY_N_DETECTORS,
    budget: DetectorBudget | None = None,
    capture_grade_sigma_theta_rad: float = D1_CAPTURE_GRADE_SIGMA_THETA_RAD,
) -> DistortionFieldFinding:
    """Assemble the distortion-field hedge finding (the pure runner)."""
    sep = separation_rad if separation_rad is not None else nearest_star_separation_rad()
    b = budget if budget is not None else DetectorBudget()
    return DistortionFieldFinding(
        differential_points=differential_curve(
            correlation_lengths_rad, separation_rad=sep, distortion_rms_rad=distortion_rms_rad
        ),
        separation_rad=sep,
        breakeven_length_rad=breakeven_correlation_length_rad(sep),
        distortion_rms_rad=distortion_rms_rad,
        tolerance_points=tolerance_curve(
            correlation_fractions,
            n_detectors=n_detectors,
            budget=b,
            capture_grade_sigma_theta_rad=capture_grade_sigma_theta_rad,
        ),
        budget=b,
        n_detectors=n_detectors,
        capture_grade_sigma_theta_rad=capture_grade_sigma_theta_rad,
    )


def format_distortion_field(finding: DistortionFieldFinding) -> str:
    """One-screen distortion-field hedge report — both readouts and the bench requirement."""
    floor_urad = finding.distortion_rms_rad * 1e6
    sep_mrad = finding.separation_rad * 1e3
    breakeven_mrad = finding.breakeven_length_rad * 1e3

    lines = [
        "Distortion-field study — terminal cross-track hedge (hardening the 3 µrad floor)",
        "  Hedge, not a gap: D1 is feasible past the noise knee (ADR 0021) — this hardens the",
        "  load-bearing bench-calibration number and adds capture margin. No bench data: the",
        "  distortion spectrum is a swept input, so the deliverable is a sensitivity curve.",
        f"  Distortion RMS (the σ_θ budget's scalar floor): {floor_urad:.2f} µrad;"
        f" nearest reference-star separation Δθ {sep_mrad:.2f} mrad.",
        "",
        "  Readout A — differential ('Gaia trick') gain vs distortion correlation length L:",
    ]
    for p in finding.differential_points:
        mark = "win " if p.improves else "WORSE"
        lines.append(
            f"    L {p.correlation_length_rad * 1e3:6.2f} mrad → ρ {p.correlation:.3f},"
            f" residual {p.residual_rad * 1e6:.2f} µrad, gain {p.gain:.2f}× [{mark}]"
        )
    lines.append(
        f"    Break-even L {breakeven_mrad:.2f} mrad (ρ = 0.5): a SMOOTH field (L above this) cuts"
        " the floor;"
    )
    lines.append(
        "    a ROUGH field (L below) makes the differential √2 worse → BENCH REQUIREMENT: measure"
    )
    lines.append(
        f"    the focal-plane distortion correlation length; bank the gain only if L ≳"
        f" {breakeven_mrad:.2f} mrad."
    )
    lines.append("")
    lines.append(
        "  Readout B — array common-mode tolerance"
        f" ({finding.n_detectors}-detector array, audit of ADR 0019's zero-common-mode split):"
    )
    for tp in finding.tolerance_points:
        verdict = "ok" if tp.meets_capture_grade else "OVER"
        lines.append(
            f"    distortion correlation ρ {tp.distortion_correlation:.2f} →"
            f" array σ_θ {tp.sigma_theta_rad * 1e6:.2f} µrad [{verdict}]"
        )
    cap_urad = finding.capture_grade_sigma_theta_rad * 1e6
    best_urad = finding.best_case_independent_sigma_rad * 1e6
    worst_urad = finding.worst_case_correlated_sigma_rad * 1e6
    lines.append(
        f"    Banked (ρ=0) {best_urad:.2f} µrad ({finding.best_case_margin:.1f}× inside the"
        f" {cap_urad:.1f} µrad capture grade) is the OPTIMISTIC end;"
    )
    lines.append(
        f"    fully common-mode (ρ=1) {worst_urad:.2f} µrad"
        f" ({finding.worst_case_margin:.2f}× inside) is the worst case."
    )
    if finding.common_mode_robust:
        lines.append(
            "    → ROBUST: capture-grade across the entire correlation range (worst case collapses"
            " to ~the single-detector grade);"
        )
        lines.append(
            "    the √N margin is best-case only, so a physically independent platform (the"
            " co-flyer) is what restores margin under pessimistic distortion correlation."
        )
    else:
        lines.append(
            f"    → tolerance: capture-grade holds up to distortion correlation"
            f" ρ ≤ {finding.common_mode_tolerance:.2f}."
        )
    return "\n".join(lines)
