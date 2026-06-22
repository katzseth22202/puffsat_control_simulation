"""Unit tests for the pure σ_θ tracker-budget gate (ADR 0018; requirement ADR 0015)."""

from __future__ import annotations

import dataclasses
import math

from puffsat_sim.constants import PLANCK_J_S, SPEED_OF_LIGHT_M_S
from puffsat_sim.guidance import CAPTURE_SIGMA_MAX_M, homing_floor_m
from puffsat_sim.tracker_budget import (
    LASER_REFERENCE_BEAM_HALFANGLE_RAD,
    REQUIRED_SIGMA_THETA_RAD,
    TARGET_SIGMA_THETA_RAD,
    BudgetTerm,
    PeakPowerSpec,
    TrackerHardware,
    beacon_intensity_w_sr,
    detector_pixels_across,
    diffraction_spot_rad,
    error_budget,
    format_peak_power_sweep,
    format_tracker_budget,
    gyro_bridge_sigma_theta_rad,
    peak_power_sweep,
    photon_energy_j,
    photon_sigma_theta_rad,
    photon_snr,
    photons_for_snr,
    power_for_centroid_floor,
    reference_star_fov_halfangle_rad,
    required_fov_halfangle_rad,
    rss_sigma_theta_rad,
    signal_photons,
    smear_sigma_theta_rad,
    tracker_budget_finding,
)


def test_photon_energy_is_hc_over_lambda() -> None:
    assert math.isclose(photon_energy_j(1.064e-6), PLANCK_J_S * SPEED_OF_LIGHT_M_S / 1.064e-6)
    # Nd:YAG photon ≈ 1.87e-19 J.
    assert math.isclose(photon_energy_j(1.064e-6), 1.868e-19, rel_tol=1e-3)


def test_beacon_intensity_is_power_over_beam_solid_angle() -> None:
    assert math.isclose(beacon_intensity_w_sr(1.0, 2.0e-3), 1.0 / (math.pi * 2.0e-3**2))


def test_diffraction_spot_is_lambda_over_d() -> None:
    assert math.isclose(diffraction_spot_rad(0.05, 1.064e-6), 1.064e-6 / 0.05)


def test_signal_photons_fall_off_as_inverse_square_range() -> None:
    hw = TrackerHardware()
    near = signal_photons(hw, 100_000.0)
    far = signal_photons(hw, 200_000.0)
    assert math.isclose(near / far, 4.0, rel_tol=1e-9)


def test_photon_snr_is_root_n_when_signal_dominates() -> None:
    hw = TrackerHardware()
    n_sig = signal_photons(hw, 300_000.0)
    # Signal swamps the 100-photon noise floor, so SNR ≈ √N.
    assert math.isclose(photon_snr(hw, 300_000.0), math.sqrt(n_sig), rel_tol=1e-4)


def test_photon_term_is_negligible_with_an_active_beacon() -> None:
    # The 1 W beacon drives the photon-limited term far below the requirement.
    assert photon_sigma_theta_rad(TrackerHardware(), 300_000.0) < 0.1e-6


def test_smear_is_streak_over_root_twelve() -> None:
    assert math.isclose(smear_sigma_theta_rad(3.0e-3, 1.0e-3), 3.0e-3 * 1.0e-3 / math.sqrt(12.0))


def test_gyro_bridge_is_arw_times_root_interval() -> None:
    assert math.isclose(gyro_bridge_sigma_theta_rad(5.8e-7, 4.0), 5.8e-7 * 2.0)


def test_error_budget_has_four_terms_and_rss_combines_them() -> None:
    terms = error_budget(TrackerHardware(), 300_000.0)
    assert len(terms) == 4
    expected = math.sqrt(sum(t.sigma_theta_rad**2 for t in terms))
    assert math.isclose(rss_sigma_theta_rad(terms), expected)


def test_rss_dominated_by_the_distortion_floor() -> None:
    terms = (BudgetTerm("a", 3.0e-6), BudgetTerm("b", 1.0e-8))
    # A negligible second term leaves the RSS at the dominant one.
    assert math.isclose(rss_sigma_theta_rad(terms), 3.0e-6, rel_tol=1e-4)


def test_required_fov_is_n_sigma_lateral_over_range() -> None:
    assert math.isclose(required_fov_halfangle_rad(141.0, 300_000.0, 3.0), 3.0 * 141.0 / 300_000.0)


def test_reference_star_fov_inverts_the_density() -> None:
    theta = reference_star_fov_halfangle_rad(3, 2.8e4)
    # The cone π·θ² holds exactly n_stars at the given density.
    assert math.isclose(math.pi * theta**2 * 2.8e4, 3.0, rel_tol=1e-9)


def test_detector_pixels_rounds_up() -> None:
    assert detector_pixels_across(5.0e-3, 10.0e-6) == 1000
    assert detector_pixels_across(5.0001e-3, 10.0e-6) == 1001


def test_default_hardware_passes_the_gate_and_meets_the_target() -> None:
    f = tracker_budget_finding()
    assert f.meets_requirement
    assert f.meets_target
    assert f.capture_floor_met
    assert f.requirement_margin > 1.0
    # The conservative default lands ~3.2 µrad, distortion-floor-limited.
    assert math.isclose(f.achieved_sigma_theta_rad, 3.18e-6, rel_tol=1e-2)


def test_driving_term_is_the_calibration_floor() -> None:
    f = tracker_budget_finding()
    assert "distortion floor" in f.driving_term.label


def test_homing_floor_uses_the_achieved_grade() -> None:
    f = tracker_budget_finding()
    assert math.isclose(
        f.homing_floor_m,
        homing_floor_m(f.achieved_sigma_theta_rad, f.speed_m_s, f.a_max_m_s2),
    )
    # The bare requirement reproduces ADR 0015's thin-margin 1.46 m reference.
    assert math.isclose(f.required_grade_floor_m, 1.45, rel_tol=2e-2)
    assert f.homing_floor_m < CAPTURE_SIGMA_MAX_M


def test_requirement_and_target_constants_match_adr_0015() -> None:
    assert REQUIRED_SIGMA_THETA_RAD == 10e-6
    assert TARGET_SIGMA_THETA_RAD == 5e-6


def test_binding_fov_is_the_larger_requirement_and_beam_covers_acquisition() -> None:
    f = tracker_budget_finding()
    assert f.fov_halfangle_rad == max(
        f.acquisition_fov_halfangle_rad, f.reference_star_fov_halfangle_rad
    )
    # Reference-star availability is the binding FOV, not acquisition.
    assert f.fov_halfangle_rad == f.reference_star_fov_halfangle_rad
    assert f.beam_covers_acquisition


def test_a_coarse_distortion_floor_fails_the_blocking_gate() -> None:
    hw = dataclasses.replace(TrackerHardware(), distortion_floor_rad=12.0e-6)
    f = tracker_budget_finding(hw)
    assert not f.meets_requirement
    assert not f.meets_target


def test_a_mid_distortion_floor_meets_requirement_but_not_target() -> None:
    hw = dataclasses.replace(TrackerHardware(), distortion_floor_rad=6.0e-6)
    f = tracker_budget_finding(hw)
    assert f.meets_requirement
    assert not f.meets_target


def test_format_reports_verdict_floor_and_acquisition() -> None:
    text = format_tracker_budget(tracker_budget_finding())
    assert "tracker budget" in text
    assert "PASS" in text
    assert "Homing floor" in text
    assert "Acquisition" in text
    assert "reference stars" in text


def test_format_marks_a_failing_gate() -> None:
    hw = dataclasses.replace(TrackerHardware(), distortion_floor_rad=12.0e-6)
    text = format_tracker_budget(tracker_budget_finding(hw))
    assert "FAIL" in text
    assert "blocks Rung D" in text


# --- Peak-power sizing sweep -------------------------------------------------------------


def test_photons_for_snr_inverts_photon_snr() -> None:
    hw = TrackerHardware()
    n = photons_for_snr(70.0, hw.noise_photons)
    powered = dataclasses.replace(hw, beacon_power_w=1.0)
    # Reconstruct the SNR from a power that yields exactly n photons at some range.
    scale = n / signal_photons(powered, 300_000.0)
    hw_scaled = dataclasses.replace(hw, beacon_power_w=scale)
    assert math.isclose(photon_snr(hw_scaled, 300_000.0), 70.0, rel_tol=1e-9)


def test_power_for_centroid_floor_round_trips_to_the_target() -> None:
    hw = TrackerHardware()
    target = 0.3e-6
    power = power_for_centroid_floor(hw, 300_000.0, target)
    powered = dataclasses.replace(hw, beacon_power_w=power)
    assert math.isclose(photon_sigma_theta_rad(powered, 300_000.0), target, rel_tol=1e-9)


def test_required_power_scales_as_range_squared() -> None:
    hw = TrackerHardware()
    near = power_for_centroid_floor(hw, 100_000.0, 0.3e-6)
    far = power_for_centroid_floor(hw, 200_000.0, 0.3e-6)
    assert math.isclose(far / near, 4.0, rel_tol=1e-9)


def test_required_power_scales_as_beam_solid_angle() -> None:
    narrow = dataclasses.replace(TrackerHardware(), beam_half_angle_rad=2.0e-3)
    wide = dataclasses.replace(TrackerHardware(), beam_half_angle_rad=4.0e-3)
    p_narrow = power_for_centroid_floor(narrow, 300_000.0, 0.3e-6)
    p_wide = power_for_centroid_floor(wide, 300_000.0, 0.3e-6)
    assert math.isclose(p_wide / p_narrow, 4.0, rel_tol=1e-9)


def test_a_bigger_aperture_needs_less_power() -> None:
    small = dataclasses.replace(TrackerHardware(), aperture_m=0.05)
    big = dataclasses.replace(TrackerHardware(), aperture_m=0.10)
    p_small = power_for_centroid_floor(small, 2_603_000.0, 0.3e-6)
    p_big = power_for_centroid_floor(big, 2_603_000.0, 0.3e-6)
    # Roughly 1/D⁴ in the shot-noise-dominated regime: doubling the dish cuts power ~16×.
    assert p_big < p_small
    assert 12.0 < p_small / p_big < 18.0


def test_sweep_grid_size_and_monotonicity() -> None:
    finding = peak_power_sweep()
    spec = finding.spec
    assert len(finding.cells) == len(spec.ranges_m) * len(spec.apertures_m)
    # Every photon term lands on the target (good-centroiding) line by construction.
    for cell in finding.cells:
        assert math.isclose(
            cell.photon_sigma_theta_rad, finding.target_sigma_theta_rad, rel_tol=1e-9
        )
    # Power rises with range (fixed aperture) and falls with aperture (fixed range).
    for aperture in spec.apertures_m:
        col = [c for c in finding.cells if c.aperture_m == aperture]
        powers = [c.peak_power_w for c in col]
        assert powers == sorted(powers)
    for range_m in spec.ranges_m:
        row = [c for c in finding.cells if c.range_m == range_m]
        powers = [c.peak_power_w for c in row]
        assert powers == sorted(powers, reverse=True)


def test_close_swarm_is_milliwatts_and_far_wide_led_is_kilowatts() -> None:
    finding = peak_power_sweep()
    by_key = {(c.range_m, c.aperture_m): c for c in finding.cells}
    # Surveyor↔PuffSat at 1 km, 5 cm dish: well under a watt.
    assert by_key[(1.0e3, 0.05)].peak_power_w < 1.0e-2
    # The far target link with a wide LED cone and a small dish is the only kW-class point.
    assert by_key[(2.603e6, 0.02)].peak_power_w > 1.0e3
    assert finding.min_peak_power_w < 1.0e-3
    assert finding.max_peak_power_w > 1.0e3


def test_narrow_beam_factor_is_the_solid_angle_ratio() -> None:
    finding = peak_power_sweep()
    expected = (finding.spec.beam_half_angle_rad / LASER_REFERENCE_BEAM_HALFANGLE_RAD) ** 2
    assert math.isclose(finding.narrow_beam_factor, expected)


def test_custom_spec_shrinks_the_grid() -> None:
    spec = PeakPowerSpec(ranges_m=(1.0e3,), apertures_m=(0.05,))
    finding = peak_power_sweep(spec)
    assert len(finding.cells) == 1


def test_format_peak_power_reports_table_and_reference() -> None:
    text = format_peak_power_sweep(peak_power_sweep())
    assert "peak-power sizing" in text
    assert "range \\ aperture" in text
    assert "narrow the beam" in text
    assert "LED-vs-laser" in text
