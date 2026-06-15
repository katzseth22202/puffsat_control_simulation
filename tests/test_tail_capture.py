"""Unit tests for the pure IS tail P(capture) core (ADR 0018 D1.x, decision 6).

The flown C3b loop is replaced by an analytic surrogate with a *closed-form* tail probability:
the arrival lateral is ``a·entry + noise``, an isotropic 2-D Gaussian of per-axis scale
``s = √((a·σ_axis)² + σ_n²)``, so ``P(escape at r) = exp(−r²/2s²)`` (Rayleigh).  That lets us
assert the IS estimator hits the true tail probability, the weights average to one, and the
brute-force-vs-IS agreement reduction fires — all without a JVM.
"""

from __future__ import annotations

import math

import numpy as np

from puffsat_sim.tail_capture import (
    TailCaptureFinding,
    TailCaptureSpec,
    TailEstimate,
    escaped_flags,
    format_tail_capture,
    sample_entry_is,
    sample_entry_nominal,
    summarize_tail_capture,
    weighted_tail_probability,
)

# Surrogate so the nominal arrival scale is ~1.35 m (the D1.1 σ at the achievable grade): the entry
# leg contributes scale 1.0 per axis (a·σ_axis = 1), the fresh noise 0.9.
_SIGMA_N: float = 0.9


def _surrogate_arrival(
    entries: np.ndarray, rng: np.random.Generator, spec: TailCaptureSpec
) -> tuple[np.ndarray, np.ndarray]:
    a = 1.0 / spec.sigma_axis_m
    noise = rng.normal(0.0, _SIGMA_N, size=entries.shape)
    arrival = a * entries + noise
    lateral_norm = np.hypot(arrival[:, 0], arrival[:, 1])
    toa = np.zeros(entries.shape[0], dtype=np.float64)
    return lateral_norm, toa


def _analytic_escape(spec: TailCaptureSpec, radius_m: float) -> float:
    s_sq = 1.0**2 + _SIGMA_N**2
    return math.exp(-(radius_m**2) / (2.0 * s_sq))


def _run_surrogate(
    spec: TailCaptureSpec, is_n: int, bf_n: int, seed: int = 1234
) -> TailCaptureFinding:
    rng = np.random.default_rng(seed)
    is_entries, is_weights = sample_entry_is(rng, spec, is_n)
    is_lat, is_toa = _surrogate_arrival(is_entries, rng, spec)
    bf_entries = sample_entry_nominal(rng, spec, bf_n)
    bf_lat, bf_toa = _surrogate_arrival(bf_entries, rng, spec)
    return summarize_tail_capture(is_weights, is_lat, is_toa, bf_lat, bf_toa, spec)


def test_is_estimate_matches_analytic_tail() -> None:
    spec = TailCaptureSpec()
    finding = _run_surrogate(spec, is_n=8000, bf_n=4000)
    analytic = _analytic_escape(spec, spec.plate_radius_m)
    # The IS plate estimator hits the closed-form tail within a few of its own σ and a factor 2.
    assert abs(finding.is_plate.probability - analytic) < 4.0 * finding.is_plate.std_error
    assert 0.5 < finding.is_plate.probability / analytic < 2.0
    assert 0.8 < finding.is_plate.weight_sum_ratio < 1.2
    # The surrogate is exactly Gaussian, so the flown tail is not heavier than its extrapolation.
    assert 0.5 < finding.tail_excess_factor < 2.0


def test_brute_force_and_is_agree_at_validation_radius() -> None:
    spec = TailCaptureSpec()
    finding = _run_surrogate(spec, is_n=8000, bf_n=4000)
    analytic_rval = _analytic_escape(spec, spec.r_val_m)
    assert abs(finding.bf_rval.probability - analytic_rval) < 4.0 * finding.bf_rval.std_error
    assert finding.validated


def test_corrupted_weights_fail_validation() -> None:
    # A biased reweighting (2× the true likelihood ratio) must be caught by the r_val agreement.
    spec = TailCaptureSpec()
    rng = np.random.default_rng(99)
    is_entries, is_weights = sample_entry_is(rng, spec, 8000)
    is_lat, is_toa = _surrogate_arrival(is_entries, rng, spec)
    bf_entries = sample_entry_nominal(rng, spec, 4000)
    bf_lat, bf_toa = _surrogate_arrival(bf_entries, rng, spec)
    finding = summarize_tail_capture(2.0 * is_weights, is_lat, is_toa, bf_lat, bf_toa, spec)
    assert not finding.validated


def test_weighted_estimator_reduces_to_plain_fraction() -> None:
    escaped = np.array([True] * 10 + [False] * 990)
    est = weighted_tail_probability(np.ones(1000), escaped)
    assert math.isclose(est.probability, 0.01, rel_tol=1e-9)
    assert math.isclose(est.ess, 1000.0, rel_tol=1e-9)
    assert math.isclose(est.weight_sum_ratio, 1.0, rel_tol=1e-9)
    # Wald σ ≈ √(p(1−p)/n) for the unit-weight case.
    assert math.isclose(est.std_error, math.sqrt(0.01 * 0.99 / 1000), rel_tol=0.02)


def test_ess_drops_with_weight_disparity() -> None:
    equal = weighted_tail_probability(np.ones(100), np.ones(100, dtype=bool))
    skewed_w = np.concatenate([np.full(1, 100.0), np.ones(99)])
    skewed = weighted_tail_probability(skewed_w, np.ones(100, dtype=bool))
    assert math.isclose(equal.ess, 100.0, rel_tol=1e-9)
    assert skewed.ess < 60.0


def test_zero_escapes_reads_as_under_resolved() -> None:
    est = weighted_tail_probability(np.ones(500), np.zeros(500, dtype=bool))
    assert est.probability == 0.0
    assert math.isinf(est.relative_error)


def test_escaped_flags_catch_toa_window() -> None:
    lateral = np.array([1.0, 1.0])
    toa = np.array([0.0, 0.05])  # second arrival is inside the plate but late
    flags = escaped_flags(lateral, toa, radius_m=5.0, toa_limit_s=0.010)
    assert flags.tolist() == [False, True]


def _estimate(p: float, half: float, n: int = 8000) -> TailEstimate:
    return TailEstimate(
        probability=p,
        std_error=half / 1.96,
        ci_low=max(0.0, p - half),
        ci_high=p + half,
        ess=float(n),
        weight_sum_ratio=1.0,
        n=n,
    )


def test_capture_ci_and_criterion() -> None:
    spec = TailCaptureSpec(target_capture=0.99)
    tight = _estimate(0.001, 0.0004)
    loose = _estimate(0.05, 0.02)
    meets = TailCaptureFinding(
        spec=spec, is_plate=tight, is_rval=tight, bf_plate=tight, bf_rval=tight, scatter_sigma_m=1.3
    )
    below = TailCaptureFinding(
        spec=spec, is_plate=loose, is_rval=loose, bf_plate=loose, bf_rval=loose, scatter_sigma_m=1.7
    )
    assert meets.p_capture > 0.99
    assert meets.meets_criterion
    assert not below.meets_criterion
    # The headline P(capture) comes from the (robust) brute-force plate estimate, not the IS one.
    assert math.isclose(meets.p_escape, meets.bf_plate.probability)


def test_tail_excess_factor_flags_heavier_than_gaussian() -> None:
    # An empirical escape of 1.6% against a σ=1.57 m Gaussian (~0.6%) reads ~2.5× heavier.
    spec = TailCaptureSpec()
    bf = _estimate(0.016, 0.006, n=400)
    finding = TailCaptureFinding(
        spec=spec, is_plate=bf, is_rval=bf, bf_plate=bf, bf_rval=bf, scatter_sigma_m=1.57
    )
    assert finding.gaussian_escape < finding.p_escape
    assert finding.tail_excess_factor > 1.5


def test_format_runs() -> None:
    spec = TailCaptureSpec()
    text = format_tail_capture(_run_surrogate(spec, is_n=4000, bf_n=2000))
    assert "importance-sampling tail" in text
    assert "P(capture)" in text
    assert "Validation @ r_val" in text
