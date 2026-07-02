"""Unit tests for the pure surveyor-anchored centering budget — the Tier-2 ~10 cm sizing module."""

from __future__ import annotations

import math
from dataclasses import replace

from puffsat_sim.centering_budget import (
    COMMITTED_PLATE_RADIUS_M,
    DEFAULT_LINK_RANGE_M,
    DIVERSE_CAMERA_CORRELATION,
    IDENTICAL_CAMERA_CORRELATION,
    R_OVER_SIGMA,
    STRETCH_PLATE_RADIUS_M,
    STROBE_RATE_HZ,
    TARGET_PLATE_RADIUS_M,
    CameraModel,
    CenteringSpec,
    QSwitchedBeacon,
    arrival_sigma_m,
    camera_reads,
    camera_sigma_theta_rad,
    centering_finding,
    format_centering_budget,
    hoop_sweep,
    per_unit_scatter_m,
    plate_radius_for_sigma,
    sigma_for_plate_radius,
)
from puffsat_sim.guidance import CAPTURE_SIGMA_MAX_M, PLATE_RADIUS_M
from puffsat_sim.tracker_budget import INTERCEPT_SPEED_M_S, LED_BEAM_HALFANGLE_RAD


def test_r_over_sigma_matches_the_committed_capture_criterion() -> None:
    # The 5 m plate ↔ σ ≤ 1.65 m criterion (2-D Rayleigh, 99 %) is the ratio that sizes any plate.
    assert math.isclose(R_OVER_SIGMA, PLATE_RADIUS_M / CAPTURE_SIGMA_MAX_M)
    assert math.isclose(R_OVER_SIGMA, 3.0, abs_tol=0.05)


def test_plate_radius_and_sigma_are_inverses() -> None:
    assert math.isclose(plate_radius_for_sigma(0.03), R_OVER_SIGMA * 0.03)
    assert math.isclose(sigma_for_plate_radius(plate_radius_for_sigma(0.02)), 0.02)


def test_arrival_sigma_is_the_rss_of_the_two_legs() -> None:
    assert math.isclose(arrival_sigma_m(0.01, 0.03), math.hypot(0.01, 0.03))
    # A hoop-dominated point and a scatter-dominated point both reduce to the larger leg.
    assert arrival_sigma_m(0.01, 0.0) == 0.01
    assert arrival_sigma_m(0.0, 0.03) == 0.03


def test_default_link_range_is_the_committed_strobe_cadence() -> None:
    # σ_δ = σ_θ·v/f: the committed band is "cm at 2–4 Hz"; the conservative default is 2 Hz.
    assert math.isclose(DEFAULT_LINK_RANGE_M, INTERCEPT_SPEED_M_S / STROBE_RATE_HZ)
    assert 2.0 <= STROBE_RATE_HZ <= 4.0
    # A 2 Hz cadence puts the nearest follower ~5.4 km back.
    assert 4.0e3 < DEFAULT_LINK_RANGE_M < 6.0e3


def test_qswitched_beacon_collects_pulse_energy_not_average_power() -> None:
    beacon = QSwitchedBeacon()
    # The gated measurement integrates one pulse: energy = avg / rep, duration = energy / peak.
    assert math.isclose(beacon.pulse_energy_j, beacon.avg_power_w / beacon.rep_rate_hz)
    assert math.isclose(beacon.pulse_duration_s, beacon.pulse_energy_j / beacon.peak_power_w)
    # A ~100 kW peak, few-hundred-mW average beacon is a ns-class, tiny-duty strobe.
    assert 1.0e-9 < beacon.pulse_duration_s < 1.0e-6
    assert beacon.duty_cycle < 1.0e-4


def test_qswitched_beacon_makes_the_link_distortion_limited() -> None:
    finding = centering_finding()
    # The whole thesis: the bright pulse pushes the photon term under the calibrated distortion.
    assert finding.photon_negligible
    assert finding.photon_sigma_theta_rad < finding.spec.camera.distortion_floor_rad
    # So the single-camera grade is the distortion floor to within a few percent (photon is small).
    single = replace(CameraModel(), n_cameras=1)
    floor = single.distortion_floor_rad
    assert math.isclose(camera_sigma_theta_rad(single, DEFAULT_LINK_RANGE_M), floor, rel_tol=0.05)


def test_a_naive_wide_cone_cw_beacon_would_be_photon_limited() -> None:
    # The counterfactual that motivates the design: a naive *wide-cone* 1 W CW beacon on the
    # gram-scale optic at the km-class link is photon-limited (over the distortion floor). Either
    # lever — the bright Q-switched pulse or the directional beam — clears it; the naive case has
    # neither.
    naive = replace(
        CameraModel(),
        beam_half_angle_rad=LED_BEAM_HALFANGLE_RAD,
        beacon=QSwitchedBeacon(peak_power_w=1.0, avg_power_w=1.0e-3, rep_rate_hz=1.0),
    )
    finding = centering_finding(CenteringSpec(camera=naive))
    assert not finding.photon_negligible
    assert finding.photon_sigma_theta_rad > finding.spec.camera.distortion_floor_rad


def test_scatter_is_sigma_theta_times_range_and_grows_with_range() -> None:
    camera = CameraModel()
    near = per_unit_scatter_m(camera, 1.0e3)
    far = per_unit_scatter_m(camera, DEFAULT_LINK_RANGE_M)
    assert math.isclose(near, camera_sigma_theta_rad(camera, 1.0e3) * 1.0e3)
    # Distortion-dominated at the design link, so scatter grows ~linearly with range.
    assert far > near


def test_identical_cameras_do_not_reduce_the_distortion_bias_but_diverse_ones_do() -> None:
    single = replace(CameraModel(), n_cameras=1)
    identical = replace(
        CameraModel(), n_cameras=3, distortion_correlation=IDENTICAL_CAMERA_CORRELATION
    )
    diverse = replace(CameraModel(), n_cameras=3, distortion_correlation=DIVERSE_CAMERA_CORRELATION)

    sigma_single = camera_sigma_theta_rad(single, DEFAULT_LINK_RANGE_M)
    sigma_identical = camera_sigma_theta_rad(identical, DEFAULT_LINK_RANGE_M)
    sigma_diverse = camera_sigma_theta_rad(diverse, DEFAULT_LINK_RANGE_M)

    # Voting three identical copies barely moves σ_θ (the shared distortion bias does not average).
    assert math.isclose(sigma_identical, sigma_single, rel_tol=0.05)
    # Physically diverse copies average the distortion with √N — a real reduction.
    assert sigma_diverse < 0.7 * sigma_single


def test_nominal_finding_reaches_the_ten_cm_target_and_is_far_below_five_metres() -> None:
    finding = centering_finding()
    # The committed claim is ~10 cm robust; the nominal point should land at or under it.
    assert finding.meets_target
    assert finding.plate_radius_m <= TARGET_PLATE_RADIUS_M + 1e-9
    # An order of magnitude smaller than the Tier-1 5 m plate (the "much stronger argument").
    assert finding.improvement_over_committed > 40.0
    assert math.isclose(
        finding.improvement_over_committed, COMMITTED_PLATE_RADIUS_M / finding.plate_radius_m
    )


def test_stretch_needs_tighter_metrology_than_nominal() -> None:
    finding = centering_finding()
    # The 5 cm stretch is contingent (σ_hoop ≤ 1 cm + calibrated camera): nominal need not reach it.
    assert STRETCH_PLATE_RADIUS_M < TARGET_PLATE_RADIUS_M
    if not finding.meets_stretch:
        assert finding.plate_radius_m > STRETCH_PLATE_RADIUS_M


def test_a_faster_train_is_a_closer_sharper_link() -> None:
    # σ_δ = σ_θ·v/f: a 4 Hz cadence halves the link range vs 2 Hz, so scatter (and plate) shrink.
    slow = centering_finding(CenteringSpec(link_range_m=INTERCEPT_SPEED_M_S / 2.0))
    fast = centering_finding(CenteringSpec(link_range_m=INTERCEPT_SPEED_M_S / 4.0))
    assert fast.spec.strobe_rate_hz > slow.spec.strobe_rate_hz
    assert fast.plate_radius_m < slow.plate_radius_m


def test_hoop_sweep_plate_grows_monotonically_with_hoop_precision() -> None:
    spec = CenteringSpec()
    points = hoop_sweep(spec, (0.002, 0.01, 0.05))
    radii = [p.plate_radius_m for p in points]
    assert radii == sorted(radii)
    # A tighter hoop always yields a smaller (or equal) plate.
    assert points[0].plate_radius_m < points[-1].plate_radius_m


def test_max_hoop_for_target_is_the_break_even_and_is_self_consistent() -> None:
    finding = centering_finding()
    if not finding.scatter_limits_target:
        # A point exactly at the break-even hoop must sit right on the 10 cm target.
        breakeven = finding.max_hoop_for_target_m
        plate = plate_radius_for_sigma(arrival_sigma_m(breakeven, finding.spec.scatter_m))
        assert math.isclose(plate, TARGET_PLATE_RADIUS_M, rel_tol=1e-6)


def test_scatter_limited_when_camera_alone_overruns_the_target() -> None:
    # A far link with a coarse (uncalibrated) camera pushes scatter past the 10 cm budget alone.
    coarse = replace(CameraModel(), distortion_floor_rad=30.0e-6, n_cameras=1)
    spec = CenteringSpec(sigma_hoop_m=0.001, camera=coarse, link_range_m=20.0e3)
    finding = centering_finding(spec)
    assert finding.scatter_limits_target
    assert finding.max_hoop_for_target_m == 0.0
    assert not finding.meets_target


def test_reference_points_span_meeting_and_missing_the_target() -> None:
    finding = centering_finding()
    labels = {p.label: p for p in finding.reference_points}
    # The calibrated-tracker classes reach the target; the uncalibrated micro-optic should not.
    calibrated = labels["mm laser rangefinder + calibrated tracker"]
    uncalibrated = labels["docking lidar + uncalibrated micro-optic"]
    assert calibrated.plate_radius_m < uncalibrated.plate_radius_m
    assert calibrated.meets_target


def test_camera_reads_cover_single_identical_and_diverse() -> None:
    reads = camera_reads(CenteringSpec())
    assert len(reads) == 3
    assert reads[0].n_cameras == 1
    # Identical voting ≈ single; diverse voting strictly better.
    assert math.isclose(reads[1].sigma_theta_rad, reads[0].sigma_theta_rad, rel_tol=0.05)
    assert reads[2].sigma_theta_rad < reads[1].sigma_theta_rad


def test_format_runs_and_mentions_the_key_reads() -> None:
    report = format_centering_budget(centering_finding())
    assert "Surveyor-anchored centering budget" in report
    assert "Q-switched beacon" in report
    assert "plate" in report
    assert "voting" in report.lower()
    assert "Reference hardware" in report
