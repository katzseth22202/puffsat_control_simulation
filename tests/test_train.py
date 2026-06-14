"""Unit tests for the pure train-mode dispersion core (D1.0/D1.1; ADR 0016/0018)."""

from __future__ import annotations

import math

import numpy as np

from puffsat_sim.guidance import CAPTURE_SIGMA_MAX_M, PLATE_RADIUS_M, PlateMiss
from puffsat_sim.terminal import FeedforwardPlan
from puffsat_sim.train import (
    TrainCaptureStats,
    TrainDispersionSpec,
    TrainEnsembleFinding,
    format_train_capture,
    format_train_ensemble,
    replay_train_entry_offset,
    replay_train_unit,
    sample_train,
    sample_train_entry_offsets,
    summarize_train_capture,
    summarize_train_ensemble,
)


def _shared_only_spec(n_units: int = 8) -> TrainDispersionSpec:
    """All per-unit σ zero: every unit of a train is identical (only shared draws act)."""
    return TrainDispersionSpec(
        n_units=n_units,
        sigma_cd_bias_frac=0.2,
        sigma_cr_bias_frac=0.2,
        sigma_f10p7_frac=0.15,
        sigma_ap_frac=0.5,
        sigma_dv_systematic_transverse_m_s=0.1,
        sigma_cd_spread_frac=0.0,
        sigma_cr_spread_frac=0.0,
        sigma_dv_scatter_radial_m_s=0.0,
        sigma_dv_scatter_transverse_m_s=0.0,
        sigma_dv_scatter_normal_m_s=0.0,
    )


def _per_unit_only_spec(n_units: int = 8) -> TrainDispersionSpec:
    """All shared σ zero: f10p7/ap stay nominal, coefficients vary per unit, dv is scatter only."""
    return TrainDispersionSpec(
        n_units=n_units,
        sigma_cd_bias_frac=0.0,
        sigma_cr_bias_frac=0.0,
        sigma_f10p7_frac=0.0,
        sigma_ap_frac=0.0,
        sigma_dv_systematic_radial_m_s=0.0,
        sigma_dv_systematic_transverse_m_s=0.0,
        sigma_dv_systematic_normal_m_s=0.0,
        sigma_cd_spread_frac=0.2,
        sigma_cr_spread_frac=0.2,
        sigma_dv_scatter_transverse_m_s=0.1,
    )


def test_per_unit_zero_makes_a_train_internally_identical() -> None:
    units = sample_train(20260613, _shared_only_spec(), train_index=0)
    first = units[0]
    for u in units[1:]:
        assert u.cd_area_over_mass == first.cd_area_over_mass
        assert u.cr_area_over_mass == first.cr_area_over_mass
        assert u.f10p7 == first.f10p7
        assert u.ap == first.ap
        assert u.dv_rtn_m_s == first.dv_rtn_m_s


def test_shared_zero_keeps_space_weather_nominal_and_varies_coefficients() -> None:
    spec = _per_unit_only_spec()
    units = sample_train(20260613, spec, train_index=0)
    for u in units:
        assert u.f10p7 == spec.f10p7
        assert u.ap == spec.ap
    # Per-unit spread makes the coefficients differ across units.
    assert len({u.cd_area_over_mass for u in units}) == len(units)


def test_space_weather_is_shared_across_a_train_even_with_per_unit_draws() -> None:
    spec = TrainDispersionSpec(n_units=6, sigma_f10p7_frac=0.15, sigma_cd_spread_frac=0.2)
    units = sample_train(7, spec, train_index=3)
    assert len({u.f10p7 for u in units}) == 1
    assert len({u.ap for u in units}) == 1


def test_injection_is_systematic_shared_plus_per_unit_scatter() -> None:
    # Scatter zero: every unit carries exactly the shared deployer systematic.
    sys_only = TrainDispersionSpec(n_units=5, sigma_dv_systematic_transverse_m_s=0.2)
    units = sample_train(11, sys_only, train_index=0)
    assert len({u.dv_rtn_m_s for u in units}) == 1
    # Systematic zero: injection is per-unit scatter, distinct and mean ~0.
    scatter_only = TrainDispersionSpec(n_units=200, sigma_dv_scatter_transverse_m_s=0.2)
    s_units = sample_train(11, scatter_only, train_index=0)
    t_components = [u.dv_rtn_m_s[1] for u in s_units]
    assert len(set(t_components)) == len(t_components)
    assert abs(float(np.mean(t_components))) < 0.05


def test_run_index_is_the_flat_train_unit_index() -> None:
    spec = TrainDispersionSpec(n_units=4)
    units = sample_train(1, spec, train_index=3)
    assert [u.run_index for u in units] == [12, 13, 14, 15]


def test_replay_reconstructs_a_unit_standalone() -> None:
    spec = TrainDispersionSpec(
        n_units=5, sigma_cd_bias_frac=0.2, sigma_cd_spread_frac=0.2, sigma_f10p7_frac=0.15
    )
    units = sample_train(42, spec, train_index=2)
    for j, u in enumerate(units):
        assert replay_train_unit(42, spec, 2, j) == u


def test_sampling_is_deterministic() -> None:
    spec = TrainDispersionSpec(n_units=4, sigma_cd_bias_frac=0.2, sigma_cd_spread_frac=0.2)
    assert sample_train(99, spec, 0) == sample_train(99, spec, 0)


def test_marginal_log_variance_is_bias_plus_spread() -> None:
    # ln(cd/cd_nom) = ln(bias) + ln(spread); both N(0, s) and independent, so the marginal
    # log-variance over many units across many trains is s_bias² + s_spread².
    spec = TrainDispersionSpec(n_units=20, sigma_cd_bias_frac=0.2, sigma_cd_spread_frac=0.2)
    logs = [
        math.log(u.cd_area_over_mass / spec.cd_area_over_mass)
        for t in range(120)
        for u in sample_train(2024, spec, t)
    ]
    s_sq = math.log(1.0 + 0.2**2)
    assert math.isclose(float(np.var(logs)), 2.0 * s_sq, rel_tol=0.15)


def _miss(x: float, y: float, toa: float = 0.0) -> PlateMiss:
    return PlateMiss(lateral_m=(x, y), toa_error_s=toa)


def test_centroid_is_the_mean_lateral_miss() -> None:
    spec = TrainDispersionSpec(n_units=2)
    stats = summarize_train_capture([_miss(10.0, 0.0), _miss(20.0, 4.0)], spec)
    assert math.isclose(stats.centroid_m[0], 15.0)
    assert math.isclose(stats.centroid_m[1], 2.0)
    assert math.isclose(stats.centroid_drift_m, math.hypot(15.0, 2.0))


def test_common_mode_shift_is_absorbed_by_the_retarget_not_the_plate() -> None:
    # A 1.5 km common-mode lateral shift + tight per-unit scatter: the plane retargets the
    # centroid (drift ≤ 2 km), and only the scatter faces the plate.
    spec = TrainDispersionSpec(n_units=4, centroid_retarget_m=2000.0)
    misses = [
        _miss(1500.0 + dx, dy) for dx, dy in [(0.5, 0.0), (-0.5, 0.4), (0.0, -0.4), (0.0, 0.0)]
    ]
    stats = summarize_train_capture(misses, spec)
    assert math.isclose(stats.centroid_drift_m, 1500.0, abs_tol=1.0)
    assert stats.retarget_ok
    assert stats.capture_about_centroid == 1.0
    assert stats.capture_absolute == 0.0  # 1.5 km ≫ the 5 m plate


def test_retarget_fails_beyond_capability() -> None:
    spec = TrainDispersionSpec(n_units=1, centroid_retarget_m=2000.0)
    stats = summarize_train_capture([_miss(2500.0, 0.0)], spec)
    assert not stats.retarget_ok


def test_scatter_sigma_is_per_axis_and_checked_against_the_capture_criterion() -> None:
    # Isotropic scatter about a zero centroid: E[|d|²] = 2σ², so σ = √(mean|d|²/2).
    spec = TrainDispersionSpec(n_units=4)
    misses = [_miss(1.0, 0.0), _miss(-1.0, 0.0), _miss(0.0, 1.0), _miss(0.0, -1.0)]
    stats = summarize_train_capture(misses, spec)
    assert math.isclose(stats.centroid_drift_m, 0.0, abs_tol=1e-9)
    assert math.isclose(stats.scatter_sigma_m, math.sqrt(0.5))
    assert stats.scatter_sigma_ok  # √0.5 ≈ 0.71 < 1.65


def test_scatter_sigma_ok_trips_above_the_criterion() -> None:
    spec = TrainDispersionSpec(n_units=2)
    big = CAPTURE_SIGMA_MAX_M * 4.0
    stats = summarize_train_capture([_miss(big, 0.0), _miss(-big, 0.0)], spec)
    assert not stats.scatter_sigma_ok


def test_toa_splits_into_centroid_drift_and_scatter() -> None:
    spec = TrainDispersionSpec(n_units=2)
    stats = summarize_train_capture([_miss(0.0, 0.0, 0.02), _miss(0.0, 0.0, 0.04)], spec)
    assert math.isclose(stats.toa_centroid_drift_s, 0.03)
    assert math.isclose(stats.toa_scatter_rms_s, 0.01)


def test_format_reports_centroid_scatter_and_capture() -> None:
    spec = TrainDispersionSpec(n_units=4)
    stats = summarize_train_capture([_miss(1.0, 0.0), _miss(-1.0, 0.0), _miss(0.0, 1.0)], spec)
    text = format_train_capture(stats)
    assert "centroid" in text.lower()
    assert "scatter" in text.lower()
    assert "capture" in text.lower()
    assert str(PLATE_RADIUS_M) in text or "5" in text


def test_capture_stats_is_a_frozen_value_type() -> None:
    spec = TrainDispersionSpec(n_units=1)
    stats = summarize_train_capture([_miss(0.0, 0.0)], spec)
    assert isinstance(stats, TrainCaptureStats)


# --- D1.1: hand-off entry offsets + ensemble finding ---


def test_entry_offsets_share_the_shared_component_when_per_unit_is_zero() -> None:
    spec = TrainDispersionSpec(
        n_units=6, sigma_entry_lateral_shared_m=149.0, sigma_entry_lateral_perunit_m=0.0
    )
    offsets = sample_train_entry_offsets(20260613, spec, train_index=0)
    assert len(offsets) == 6
    for off in offsets[1:]:
        assert off == offsets[0]


def test_entry_offsets_differ_per_unit_and_replay_standalone() -> None:
    spec = TrainDispersionSpec(
        n_units=5, sigma_entry_lateral_shared_m=149.0, sigma_entry_lateral_perunit_m=141.0
    )
    offsets = sample_train_entry_offsets(7, spec, train_index=2)
    assert len({off for off in offsets}) == len(offsets)
    for j, off in enumerate(offsets):
        assert replay_train_entry_offset(7, spec, 2, j) == off


def test_entry_offsets_are_independent_of_the_coefficient_draws() -> None:
    # The entry-offset seed tree is masked off the coefficient/injection tree, so a train's
    # shared entry offset is not the same byte-stream as its shared coefficient bias.
    spec = TrainDispersionSpec(n_units=3, sigma_entry_lateral_perunit_m=0.0)
    inputs = sample_train(123, spec, train_index=0)
    offsets = sample_train_entry_offsets(123, spec, train_index=0)
    # Shared coefficient factor (cd/cd_nom) and the shared entry x must not coincide.
    assert inputs[0].cd_area_over_mass / spec.cd_area_over_mass != offsets[0][0]


def test_entry_offset_variance_splits_into_shared_and_per_unit() -> None:
    # σ is the 2-D lateral magnitude, sampled isotropically (per-axis σ/√2), so the per-axis var
    # of (shared + per-unit) = (σ_shared² + σ_perunit²)/2; the per-train centroid (mean over units)
    # keeps only the shared half as n_units grows.
    spec = TrainDispersionSpec(
        n_units=40, sigma_entry_lateral_shared_m=150.0, sigma_entry_lateral_perunit_m=141.0
    )
    xs = [off[0] for t in range(120) for off in sample_train_entry_offsets(99, spec, t)]
    centroids_x = [
        float(np.mean([off[0] for off in sample_train_entry_offsets(99, spec, t)]))
        for t in range(120)
    ]
    assert math.isclose(float(np.var(xs)), (150.0**2 + 141.0**2) / 2.0, rel_tol=0.15)
    assert math.isclose(float(np.var(centroids_x)), 150.0**2 / 2.0, rel_tol=0.2)


def _plan(dv_m_s: float) -> FeedforwardPlan:
    return FeedforwardPlan(
        commands=(), dv_m_s=dv_m_s, mass_kg=25.0, saturated=False, peak_slew_rate_deg_s=0.0
    )


def test_ensemble_finding_assembles_capture_propellant_and_perigee() -> None:
    spec = TrainDispersionSpec(n_units=3, midcourse_dv_m_s=2.19)
    misses = [_miss(0.5, 0.0), _miss(-0.5, 0.3), _miss(0.0, -0.3)]
    plans = [_plan(1.0), _plan(1.5), _plan(1.2)]
    perigees = [52_000.0, 48_000.0, 50_000.0]
    finding = summarize_train_ensemble(misses, plans, perigees, spec)
    assert isinstance(finding, TrainEnsembleFinding)
    assert finding.capture.n_units == 3
    assert math.isclose(finding.terminal_dv_max_m_s, 1.5)
    # Worst-case total Δv = midcourse + max terminal; under 2% at the conservative Isp anchor.
    assert math.isclose(finding.total_dv_worst_m_s, 2.19 + 1.5)
    assert finding.within_budget
    assert finding.perigee_min_m == 48_000.0
    assert finding.perigee_max_m == 52_000.0


def test_ensemble_format_reports_capture_propellant_and_perigee() -> None:
    spec = TrainDispersionSpec(n_units=2)
    finding = summarize_train_ensemble(
        [_miss(1.0, 0.0), _miss(-1.0, 0.0)], [_plan(1.0), _plan(1.0)], [50_000.0, 50_000.0], spec
    )
    text = format_train_ensemble(finding)
    assert "capture" in text.lower()
    assert "propellant" in text.lower()
    assert "perigee" in text.lower()
