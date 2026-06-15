"""Unit tests for the pure distortion-field study — the terminal cross-track hedge."""

from __future__ import annotations

import math

import numpy as np

from puffsat_sim.distortion_field import (
    DEFAULT_CORRELATION_FRACTIONS,
    DEFAULT_CORRELATION_LENGTHS_RAD,
    DISTORTION_RMS_RAD,
    ArrayTolerancePoint,
    DetectorBudget,
    DifferentialPoint,
    DistortionFieldFinding,
    array_sigma_theta_rad,
    breakeven_correlation_length_rad,
    common_mode_tolerance,
    differential_curve,
    differential_residual_rad,
    distortion_field_finding,
    format_distortion_field,
    gaussian_correlation,
    nearest_star_separation_rad,
    tolerance_curve,
)
from puffsat_sim.tracker_fusion import (
    D1_CAPTURE_GRADE_SIGMA_THETA_RAD,
    target_array,
    tracker_fusion_finding,
)


def test_nearest_star_separation_is_half_inverse_sqrt_density() -> None:
    assert math.isclose(nearest_star_separation_rad(2.8e4), 1.0 / (2.0 * math.sqrt(2.8e4)))
    # Denser star fields → closer reference stars → smaller separation.
    assert nearest_star_separation_rad(1e5) < nearest_star_separation_rad(1e4)


def test_gaussian_correlation_is_one_at_zero_lag_and_decays_with_separation() -> None:
    assert math.isclose(gaussian_correlation(0.0, 3e-3), 1.0)
    assert gaussian_correlation(3e-3, 3e-3) < 1.0
    # A longer correlation length holds the field more correlated at a fixed lag.
    assert gaussian_correlation(3e-3, 10e-3) > gaussian_correlation(3e-3, 1e-3)


def test_differential_is_sqrt2_worse_for_a_rough_field_and_vanishes_for_a_smooth_one() -> None:
    sep = 3e-3
    # L ≪ Δθ (rough, uncorrelated): residual → √2·σ_d — worse than the absolute floor.
    rough = differential_residual_rad(DISTORTION_RMS_RAD, sep, 1e-5)
    assert math.isclose(rough, DISTORTION_RMS_RAD * math.sqrt(2.0), rel_tol=1e-6)
    # L ≫ Δθ (smooth, correlated): residual → 0 — the common distortion cancels.
    smooth = differential_residual_rad(DISTORTION_RMS_RAD, sep, 1.0)
    assert smooth < 0.05 * DISTORTION_RMS_RAD


def test_breakeven_length_makes_the_differential_equal_the_floor() -> None:
    sep = 3e-3
    length = breakeven_correlation_length_rad(sep)
    # At the break-even length ρ = 0.5 and the differential exactly equals the absolute floor.
    assert math.isclose(gaussian_correlation(sep, length), 0.5, rel_tol=1e-9)
    residual = differential_residual_rad(DISTORTION_RMS_RAD, sep, length)
    assert math.isclose(residual, DISTORTION_RMS_RAD)
    # A longer (smoother) field wins; a shorter (rougher) field backfires.
    assert differential_residual_rad(DISTORTION_RMS_RAD, sep, 2.0 * length) < DISTORTION_RMS_RAD
    assert differential_residual_rad(DISTORTION_RMS_RAD, sep, 0.5 * length) > DISTORTION_RMS_RAD


def test_differential_point_gain_and_improves_are_consistent() -> None:
    sep = nearest_star_separation_rad()
    win = DifferentialPoint(30e-3, sep, DISTORTION_RMS_RAD)
    lose = DifferentialPoint(0.3e-3, sep, DISTORTION_RMS_RAD)
    assert win.gain > 1.0 and win.improves
    assert lose.gain < 1.0 and not lose.improves
    assert math.isclose(win.gain, win.distortion_rms_rad / win.residual_rad)


def test_differential_curve_improves_monotonically_with_correlation_length() -> None:
    points = differential_curve()
    assert len(points) == len(DEFAULT_CORRELATION_LENGTHS_RAD)
    residuals = [p.residual_rad for p in points]
    # Smoother fields (longer L) always leave a smaller differential residual.
    assert residuals == sorted(residuals, reverse=True)


def test_array_at_zero_correlation_reproduces_the_committed_fusion_grade() -> None:
    # Readout B at ρ=0 must reproduce ADR 0019's optimistic 5-array grade (the audit pins, not
    # moves, the committed number).
    sigma = array_sigma_theta_rad(DetectorBudget(), 5, 0.0)
    fusion = tracker_fusion_finding([target_array(5)]).effective_sigma_theta_rad
    assert math.isclose(sigma, fusion, rel_tol=2e-2)
    assert math.isclose(sigma * 1e6, 1.62, rel_tol=1e-2)


def test_fully_correlated_distortion_collapses_the_array_toward_a_single_detector() -> None:
    one = array_sigma_theta_rad(DetectorBudget(), 1, 0.0)
    array_indep = array_sigma_theta_rad(DetectorBudget(), 5, 0.0)
    array_common = array_sigma_theta_rad(DetectorBudget(), 5, 1.0)
    # Independent distortion → √N gain; fully common-mode → ~no gain (back to a single detector).
    assert array_indep < array_common
    assert array_common < one
    assert math.isclose(array_common, one, rel_tol=5e-2)


def test_array_sigma_increases_monotonically_with_distortion_correlation() -> None:
    fractions = (0.0, 0.25, 0.5, 0.75, 1.0)
    sigmas = [array_sigma_theta_rad(DetectorBudget(), 5, rho) for rho in fractions]
    assert sigmas == sorted(sigmas)


def test_common_mode_tolerance_solves_the_grade_crossing() -> None:
    # A budget whose distortion floor crosses the capture grade somewhere in ρ ∈ (0, 1): array-grade
    # at ρ=0 (full √N) stays inside, at ρ=1 (no gain) falls outside.
    budget = DetectorBudget(sigma_distortion_rad=5e-6)
    rho = common_mode_tolerance(budget, 5, D1_CAPTURE_GRADE_SIGMA_THETA_RAD)
    assert 0.0 < rho < 1.0
    # At the tolerance correlation, the array grade exactly equals the capture grade.
    assert math.isclose(
        array_sigma_theta_rad(budget, 5, rho), D1_CAPTURE_GRADE_SIGMA_THETA_RAD, rel_tol=1e-9
    )


def test_default_budget_is_robust_across_the_full_correlation_range() -> None:
    # The headline audit result: the worst case (fully common-mode) collapses to ~the
    # single-detector grade, which is itself at the capture threshold → tolerance ≥ 1.
    finding = distortion_field_finding()
    assert finding.common_mode_robust
    assert finding.common_mode_tolerance >= 1.0
    assert finding.worst_case_correlated_sigma_rad <= finding.capture_grade_sigma_theta_rad
    assert all(p.meets_capture_grade for p in finding.tolerance_points)


def test_finding_margins_bracket_the_optimistic_and_pessimistic_ends() -> None:
    finding = distortion_field_finding()
    # The banked √N margin (~2×) is the optimistic end; full common-mode squeezes it to ~1×.
    assert math.isclose(finding.best_case_margin, 2.0, rel_tol=0.1)
    assert 1.0 <= finding.worst_case_margin < 1.1
    assert finding.best_case_margin > finding.worst_case_margin


def test_tolerance_curve_default_counts_and_capture_flags() -> None:
    points = tolerance_curve()
    assert len(points) == len(DEFAULT_CORRELATION_FRACTIONS)
    assert all(isinstance(p, ArrayTolerancePoint) for p in points)
    assert all(p.meets_capture_grade for p in points)


def test_finding_is_a_frozen_value_type() -> None:
    finding = distortion_field_finding()
    assert isinstance(finding, DistortionFieldFinding)
    assert isinstance(finding.budget, DetectorBudget)
    assert finding.distortion_rms_rad == DISTORTION_RMS_RAD


def test_format_reports_both_readouts_breakeven_and_robust_verdict() -> None:
    text = format_distortion_field(distortion_field_finding())
    assert "Readout A" in text
    assert "Readout B" in text
    assert "Break-even" in text
    assert "BENCH REQUIREMENT" in text
    assert "ROBUST" in text
    assert "µrad" in text


def test_random_field_realization_matches_the_analytic_differential() -> None:
    # The module is analytic because only the field's second-order statistics enter; validate that
    # a literal Gaussian random-field draw with correlation length L reproduces the σ_diff formula.
    rng = np.random.default_rng(20260615)
    sigma_d = DISTORTION_RMS_RAD
    length = 4e-3
    sep = nearest_star_separation_rad()

    # A small line of focal-plane points; the beacon sits at 0, its nearest star at the separation.
    coords = np.array([0.0, sep, 2.0 * sep, 3.0 * sep])
    lag = coords[:, None] - coords[None, :]
    cov = sigma_d**2 * np.exp(-(lag**2) / (2.0 * length**2))
    chol = np.linalg.cholesky(cov + 1e-18 * np.eye(coords.size))

    draws = chol @ rng.standard_normal((coords.size, 200_000))
    differential = draws[0, :] - draws[1, :]  # beacon minus nearest star

    analytic = differential_residual_rad(sigma_d, sep, length)
    assert math.isclose(float(differential.std()), analytic, rel_tol=2e-2)
