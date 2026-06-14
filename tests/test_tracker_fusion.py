"""Unit tests for the pure multi-tracker fusion gate (ADR 0019)."""

from __future__ import annotations

import math

from puffsat_sim.tracker_fusion import (
    ARRAY_N_DETECTORS,
    COFLYER_N_DETECTORS,
    COFLYER_RANGE_M,
    D1_CAPTURE_GRADE_SIGMA_THETA_RAD,
    DEFAULT_RANGING_SIGMA_M,
    DETECTOR_INDEP_SIGMA_RAD,
    GPS_CEILING_M,
    SMEAR_COMMON_SIGMA_RAD,
    TARGET_RANGE_M,
    CoflyerPhasing,
    Tracker,
    TrackerFusionFinding,
    angular_sigma_theta_rad,
    array_with_coflyer,
    coflyer,
    effective_sigma_theta_rad,
    format_coflyer_phasing,
    format_tracker_fusion,
    fuse_lateral_sigma_m,
    fused_tracker_grade,
    lateral_sigma_m,
    phasing_verdict,
    single_target_detector,
    target_array,
    target_array_only,
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


def test_architectures_pin_the_adr0019_detector_counts() -> None:
    assert single_target_detector() == (target_array(1),)
    (array,) = target_array_only()
    assert array.n_detectors == ARRAY_N_DETECTORS and array.range_m == TARGET_RANGE_M
    target, rocket = array_with_coflyer()
    assert target.n_detectors == ARRAY_N_DETECTORS
    assert rocket.n_detectors == COFLYER_N_DETECTORS and rocket.range_m == COFLYER_RANGE_M


def test_fused_tracker_grade_carries_the_effective_sigma_theta_and_ranging_sigma() -> None:
    # The re-key (ADR 0019 dec 4): the loop's σ_θ is the architecture's fused effective grade.
    trackers = array_with_coflyer()
    grade = fused_tracker_grade(trackers)
    assert grade.sigma_range_m == DEFAULT_RANGING_SIGMA_M
    assert grade.sigma_theta_rad is not None
    assert math.isclose(
        grade.sigma_theta_rad, tracker_fusion_finding(trackers).effective_sigma_theta_rad
    )


def test_fused_grade_recovers_capture_where_the_single_tracker_ceiling_fails() -> None:
    # The 10 µrad single-tracker ceiling that D1.1 failed at...
    single_ceiling = fused_tracker_grade(
        [Tracker(range_m=TARGET_RANGE_M, n_detectors=1, sigma_indep_rad=10e-6)]
    )
    assert single_ceiling.sigma_theta_rad is not None
    assert single_ceiling.sigma_theta_rad > D1_CAPTURE_GRADE_SIGMA_THETA_RAD
    # ...is recovered to capture-grade by the fused target array + co-flyer.
    fused = fused_tracker_grade(array_with_coflyer())
    assert fused.sigma_theta_rad is not None
    assert fused.sigma_theta_rad <= D1_CAPTURE_GRADE_SIGMA_THETA_RAD


def test_fused_tracker_grade_respects_an_overridden_ranging_sigma() -> None:
    grade = fused_tracker_grade(target_array_only(), sigma_range_m=0.2)
    assert grade.sigma_range_m == 0.2


def test_phasing_verdict_reduces_window_samples_to_max_range_and_alt_band() -> None:
    finding = phasing_verdict(
        window_ranges_m=[120e3, 300e3, 220e3],
        window_rocket_alts_m=[800e3, 400e3, 160e3],
        window_alt_hi_m=800e3,
        window_alt_lo_m=200e3,
    )
    assert math.isclose(finding.max_range_m, 300e3)
    assert math.isclose(finding.max_rocket_alt_m, 800e3)
    assert math.isclose(finding.min_rocket_alt_m, 160e3)


def test_phasing_is_feasible_when_close_and_in_the_gps_volume() -> None:
    # Within 500 km of the centroid and well below the GPS constellation ceiling.
    finding = phasing_verdict([120e3, 300e3], [800e3, 160e3], 800e3, 200e3)
    assert finding.range_ok
    assert finding.gps_ok
    assert finding.feasible


def test_phasing_fails_when_drifting_beyond_the_angle_useful_range() -> None:
    finding = phasing_verdict([400e3, 900e3], [600e3, 200e3], 800e3, 200e3)
    assert not finding.range_ok
    assert not finding.feasible


def test_phasing_fails_when_above_the_gps_ceiling() -> None:
    finding = phasing_verdict([100e3], [GPS_CEILING_M + 1e3], 800e3, 200e3)
    assert finding.range_ok
    assert not finding.gps_ok
    assert not finding.feasible


def test_format_coflyer_phasing_reports_both_legs_and_verdict() -> None:
    text = format_coflyer_phasing(phasing_verdict([120e3], [400e3], 800e3, 200e3))
    assert "phasing" in text.lower()
    assert "GPS" in text
    assert "range" in text.lower()


def test_coflyer_phasing_is_a_frozen_value_type() -> None:
    finding = phasing_verdict([1.0], [2.0], 800e3, 200e3)
    assert isinstance(finding, CoflyerPhasing)
    assert finding.angle_useful_range_m == COFLYER_RANGE_M
    assert finding.gps_ceiling_m == GPS_CEILING_M
