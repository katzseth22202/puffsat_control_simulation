"""Unit tests for the pure corrector-in-loop validation reduction (ADR 0018 D1.x)."""

from __future__ import annotations

import numpy as np

from puffsat_sim.corrector_validation import (
    CorrectorValidationFinding,
    format_corrector_validation,
    summarize_corrector_validation,
)

# A toy 3×6 sensitivity: transverse velocity (col 4) dominates the lateral (T,N) rows, the C0
# structure in miniature.  ‖Φ_TN col4‖ ≈ 2.0e5 m per m/s, so a 0.66 mm/s draw → ~141 m lateral.
_PHI = np.zeros((3, 6), dtype=np.float64)
_PHI[1, 4] = 1.5e5  # T-miss per T-velocity
_PHI[2, 4] = 1.0e5  # N-miss per T-velocity (so ‖Φ_TN col4‖ ≈ 1.8e5)
_PHI[0, 0] = 0.6  # a weak radial-pos → radial-miss term (pinned axis)

# A diagonal apogee-RTN Σ dominated by transverse velocity (the C1 binding axis).
_SIGMA = np.diag([5.0**2, 5.0**2, 5.0**2, 1e-4**2, 6.6e-4**2, 1e-4**2]).astype(np.float64)
_PROXY_M = 141.0


def _draws(n: int, seed: int = 7) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.multivariate_normal(np.zeros(6), _SIGMA, size=n)


def test_perfectly_linear_corrector_reads_as_validated() -> None:
    draws = _draws(200)
    crossing = draws @ _PHI.T  # measured == Φ·δ exactly
    handoff = 0.4 * crossing[:, 1:3]  # hand-off displacement ~2.5× smaller than the crossing miss
    finding = summarize_corrector_validation(draws, crossing, handoff, _PHI, _SIGMA, _PROXY_M)
    assert finding.linear
    assert finding.unbiased
    assert finding.crossing_sigma_consistent
    assert finding.handoff_conservative
    assert finding.validated
    assert finding.linearity_residual_rel < 1e-9


def test_measured_crossing_sigma_tracks_phi_sigma_phi() -> None:
    draws = _draws(500)
    crossing = draws @ _PHI.T
    finding = summarize_corrector_validation(
        draws, crossing, 0.5 * crossing[:, 1:3], _PHI, _SIGMA, _PROXY_M
    )
    # The measured RMS lateral and the ΦΣΦᵀ prediction must agree (the end-to-end number check).
    assert finding.predicted_crossing_lateral_sigma_m > 100.0
    assert np.isclose(
        finding.measured_crossing_lateral_sigma_m,
        finding.predicted_crossing_lateral_sigma_m,
        rtol=0.1,
    )


def test_sigma_consistency_band_is_sample_size_aware() -> None:
    # A 30% σ deviation is sampling noise at tiny N (band ≈ 3/√(2N)) but a real gap at large N.
    draws_small = _draws(8)
    crossing_small = 1.3 * (draws_small @ _PHI.T)  # measured σ 30% above the ΦΣΦᵀ prediction
    small = summarize_corrector_validation(
        draws_small, crossing_small, 0.5 * crossing_small[:, 1:3], _PHI, _SIGMA, _PROXY_M
    )
    assert small.crossing_sigma_consistent  # 30% < 3/√16 ≈ 75% band

    draws_big = _draws(400)
    crossing_big = 1.3 * (draws_big @ _PHI.T)
    big = summarize_corrector_validation(
        draws_big, crossing_big, 0.5 * crossing_big[:, 1:3], _PHI, _SIGMA, _PROXY_M
    )
    assert not big.crossing_sigma_consistent  # 30% > max(15%, 3/√800 ≈ 11%) band


def test_nonlinear_residual_flips_the_linear_read() -> None:
    draws = _draws(200)
    rng = np.random.default_rng(11)
    # A 20%-of-signal scatter that Φ·δ cannot explain — the corrector residual is not linear.
    crossing = draws @ _PHI.T
    crossing = crossing + rng.normal(0.0, 0.2 * np.std(crossing), size=crossing.shape)
    finding = summarize_corrector_validation(
        draws, crossing, 0.5 * crossing[:, 1:3], _PHI, _SIGMA, _PROXY_M
    )
    assert finding.linearity_residual_rel > 0.05
    assert not finding.linear
    assert not finding.validated


def test_a_real_bias_is_caught() -> None:
    draws = _draws(200)
    crossing = draws @ _PHI.T
    crossing = crossing + np.array([0.0, 80.0, 80.0])  # a fixed lateral offset on every miss
    finding = summarize_corrector_validation(
        draws, crossing, 0.5 * crossing[:, 1:3], _PHI, _SIGMA, _PROXY_M
    )
    assert finding.bias_lateral_m > 50.0
    assert not finding.unbiased
    assert not finding.validated


def test_handoff_offset_above_the_proxy_is_not_conservative() -> None:
    draws = _draws(200)
    crossing = draws @ _PHI.T
    handoff = 1.5 * crossing[:, 1:3]  # actual hand-off bigger than the fed crossing proxy
    finding = summarize_corrector_validation(draws, crossing, handoff, _PHI, _SIGMA, _PROXY_M)
    assert finding.measured_handoff_lateral_sigma_m > _PROXY_M
    assert not finding.handoff_conservative
    assert finding.conservatism_factor < 1.0


def test_conservatism_factor_is_the_crossing_over_handoff_ratio() -> None:
    draws = _draws(200)
    crossing = draws @ _PHI.T
    handoff = 0.25 * crossing[:, 1:3]  # hand-off 4× smaller than the crossing miss
    finding = summarize_corrector_validation(draws, crossing, handoff, _PHI, _SIGMA, _PROXY_M)
    assert np.isclose(finding.conservatism_factor, 4.0, rtol=1e-6)


def test_finding_is_a_frozen_value_type() -> None:
    draws = _draws(8)
    crossing = draws @ _PHI.T
    finding = summarize_corrector_validation(
        draws, crossing, 0.5 * crossing[:, 1:3], _PHI, _SIGMA, _PROXY_M
    )
    assert isinstance(finding, CorrectorValidationFinding)
    assert finding.n_samples == 8


def test_format_reports_linearity_sigma_and_conservatism() -> None:
    draws = _draws(64)
    crossing = draws @ _PHI.T
    text = format_corrector_validation(
        summarize_corrector_validation(
            draws, crossing, 0.4 * crossing[:, 1:3], _PHI, _SIGMA, _PROXY_M
        )
    )
    assert "Linearity" in text
    assert "conservative" in text.lower()
    assert "VALIDATED" in text or "validated" in text.lower()
