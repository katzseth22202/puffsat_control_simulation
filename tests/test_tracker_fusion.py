"""Unit tests for the pure multi-tracker fusion gate (ADR 0019)."""

from __future__ import annotations

import math

from puffsat_sim.tracker_fusion import (
    COFLYER_RANGE_M,
    D1_CAPTURE_GRADE_SIGMA_THETA_RAD,
    DETECTOR_INDEP_SIGMA_RAD,
    SMEAR_COMMON_SIGMA_RAD,
    TARGET_RANGE_M,
    Tracker,
    TrackerFusionFinding,
    angular_sigma_theta_rad,
    coflyer,
    effective_sigma_theta_rad,
    format_tracker_fusion,
    fuse_lateral_sigma_m,
    lateral_sigma_m,
    target_array,
    tracker_fusion_finding,
)


def test_single_detector_reproduces_the_sigma_theta_gate() -> None:
    # One detector at the per-detector split reproduces the σ_θ budget gate's 3.2 µrad RSS.
    t = Tracker(range_m=TARGET_RANGE_M, n_detectors=1)
    assert math.isclose(
        angular_sigma_theta_rad(t),
        math.hypot(SMEAR_COMMON_SIGMA_RAD, DETECTOR_INDEP_SIGMA_RAD),
    )
    assert math.isclose(angular_sigma_theta_rad(t) * 1e6, 3.2, rel_tol=5e-2)


def test_n_detectors_average_the_independent_term_toward_the_common_floor() -> None:
    one = angular_sigma_theta_rad(Tracker(range_m=TARGET_RANGE_M, n_detectors=1))
    five = angular_sigma_theta_rad(Tracker(range_m=TARGET_RANGE_M, n_detectors=5))
    many = angular_sigma_theta_rad(Tracker(range_m=TARGET_RANGE_M, n_detectors=10_000))
    assert five < one
    # The independent term averages as 1/√N; the common floor does not.
    assert math.isclose(
        five, math.sqrt(SMEAR_COMMON_SIGMA_RAD**2 + DETECTOR_INDEP_SIGMA_RAD**2 / 5)
    )
    assert math.isclose(many, SMEAR_COMMON_SIGMA_RAD, rel_tol=1e-3)


def test_lateral_sigma_is_angular_times_range_rss_relgeom() -> None:
    t = Tracker(range_m=500e3, n_detectors=3, rel_geom_sigma_m=2.0)
    expected = math.hypot(angular_sigma_theta_rad(t) * 500e3, 2.0)
    assert math.isclose(lateral_sigma_m(t), expected)


def test_fuse_is_inverse_variance() -> None:
    a = Tracker(range_m=1.0, n_detectors=1, sigma_indep_rad=3.0, sigma_common_rad=0.0)
    b = Tracker(range_m=1.0, n_detectors=1, sigma_indep_rad=4.0, sigma_common_rad=0.0)
    # lateral σ = 3 and 4 (range 1 m) → fused 1/√(1/9+1/16) = 2.4.
    assert math.isclose(fuse_lateral_sigma_m([a, b]), 2.4)


def test_effective_sigma_theta_is_lateral_over_design_range() -> None:
    assert math.isclose(effective_sigma_theta_rad(5.206, TARGET_RANGE_M), 5.206 / TARGET_RANGE_M)


def test_target_array_alone_recovers_capture_grade_without_the_coflyer() -> None:
    # 5× 10 µrad detectors on the target → effective σ_θ inside the D1.1 capture-grade.
    finding = tracker_fusion_finding([target_array(5)])
    assert finding.effective_sigma_theta_rad < D1_CAPTURE_GRADE_SIGMA_THETA_RAD
    assert finding.meets_capture_grade
    assert math.isclose(finding.effective_sigma_theta_rad * 1e6, 1.6, rel_tol=0.1)


def test_close_coflyer_beats_a_far_single_target_detector() -> None:
    # A single detector 5× closer (500 vs 2603 km) lands less lateral error even with the
    # GNSS relative-geometry floor — the range advantage is the strong lever.
    far = lateral_sigma_m(Tracker(range_m=TARGET_RANGE_M, n_detectors=1))
    close = lateral_sigma_m(Tracker(range_m=COFLYER_RANGE_M, n_detectors=1, rel_geom_sigma_m=2.0))
    assert close < far


def test_fused_architecture_is_strongly_capture_grade() -> None:
    finding = tracker_fusion_finding([target_array(5), coflyer(3)])
    assert finding.meets_capture_grade
    # Fusing the close co-flyer with the target array beats either alone.
    assert (
        finding.effective_sigma_theta_rad
        < tracker_fusion_finding([target_array(5)]).effective_sigma_theta_rad
    )
    assert math.isclose(finding.effective_sigma_theta_rad * 1e6, 0.76, rel_tol=0.15)


def test_a_worse_per_detector_floor_fails_single_but_a_close_coflyer_recovers() -> None:
    # If the per-detector distortion floor is a full 10 µrad, one target detector fails D1.1...
    bad = Tracker(range_m=TARGET_RANGE_M, n_detectors=1, sigma_indep_rad=10e-6)
    assert not tracker_fusion_finding([bad]).meets_capture_grade
    # ...but the same crude detector on the co-flyer at 500 km is capture-grade (range lever).
    bad_close = Tracker(
        range_m=COFLYER_RANGE_M, n_detectors=1, sigma_indep_rad=10e-6, rel_geom_sigma_m=2.0
    )
    assert tracker_fusion_finding([bad_close]).meets_capture_grade


def test_finding_is_a_frozen_value_type_with_a_homing_floor() -> None:
    finding = tracker_fusion_finding([target_array(5), coflyer(3)])
    assert isinstance(finding, TrackerFusionFinding)
    assert finding.homing_floor_m > 0.0


def test_format_reports_effective_grade_and_verdict() -> None:
    text = format_tracker_fusion(tracker_fusion_finding([target_array(5), coflyer(3)]))
    assert "effective" in text.lower()
    assert "µrad" in text or "urad" in text
    assert "capture" in text.lower()
