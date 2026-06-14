"""Unit tests for the pure apogee nav-constellation sizing (ADR 0020, Lever 3)."""

from __future__ import annotations

import math

import numpy as np

from puffsat_sim.apogee_nav import (
    APOGEE_RADIUS_M,
    DEFAULT_SWEEP_COUNTS,
    KA_FREQ_HZ,
    NORMAL_AXIS,
    PASSIVE_RECEIVER,
    REQUIRED_TVEL_SIGMA_M_S,
    TARGET_TVEL_SIGMA_M_S,
    TRANSVERSE_AXIS,
    ApogeeNavFinding,
    GdopSweepFinding,
    apogee_nav_finding,
    apogee_position,
    axis_observable,
    carrier_phase_velocity_sigma_m_s,
    carrier_to_noise_density_dbhz,
    constellation_los_units,
    downlink_cn0_dbhz,
    format_apogee_nav,
    format_gdop_sweep,
    free_space_path_loss_db,
    gdop_sweep,
    gdop_sweep_finding,
    line_of_sight_unit_vectors,
    min_members_for_target,
    parabolic_gain_dbi,
    transverse_velocity_sigma_m_s,
    velocity_dop,
)


def test_free_space_path_loss_matches_the_standard_formula() -> None:
    # FSPL(dB) = 92.45 + 20log10(f_GHz) + 20log10(R_km).
    fspl = free_space_path_loss_db(150_000e3, 30e9)
    expected = 92.45 + 20.0 * math.log10(30.0) + 20.0 * math.log10(150_000.0)
    assert math.isclose(fspl, expected, abs_tol=0.05)


def test_parabolic_gain_of_a_1m_ka_dish_is_about_48_dbi() -> None:
    assert math.isclose(parabolic_gain_dbi(1.0, 30e9, 0.6), 47.7, abs_tol=0.5)


def test_cn0_is_eirp_minus_path_loss_plus_gain_minus_noise() -> None:
    cn0 = carrier_to_noise_density_dbhz(
        eirp_dbw=57.7, range_m=150_000e3, freq_hz=30e9, rx_gain_dbi=0.0, system_temp_k=400.0
    )
    # ~35 dB-Hz at the apogee range with a modest dish (the downlink the PuffSat receives).
    assert 30.0 < cn0 < 40.0


def test_default_downlink_closes_at_the_apogee_range() -> None:
    assert 30.0 < downlink_cn0_dbhz() < 40.0
    # Farther members (larger range) only lose link — monotone in range.
    assert downlink_cn0_dbhz(2 * APOGEE_RADIUS_M) < downlink_cn0_dbhz(APOGEE_RADIUS_M)


def test_carrier_phase_velocity_sigma_is_sub_mm_s_and_improves_with_cn0_and_time() -> None:
    sigma = carrier_phase_velocity_sigma_m_s(35.0, KA_FREQ_HZ, 1.0)
    assert (
        sigma < TARGET_TVEL_SIGMA_M_S
    )  # already under the requirement on a single 1 s measurement
    # T^-3/2 in σ_v and 10^(-dB/20): more C/N0 and longer integration both shrink it.
    assert carrier_phase_velocity_sigma_m_s(45.0, KA_FREQ_HZ, 1.0) < sigma
    assert carrier_phase_velocity_sigma_m_s(35.0, KA_FREQ_HZ, 4.0) < sigma


def test_a_higher_carrier_gives_finer_velocity_at_equal_cn0() -> None:
    # σ_v = λ·σ_f, so Ka (30 GHz) beats L-band (1.5 GHz) at the same C/N0.
    ka = carrier_phase_velocity_sigma_m_s(35.0, 30e9, 1.0)
    el = carrier_phase_velocity_sigma_m_s(35.0, 1.5e9, 1.0)
    assert ka < el


def test_occulted_members_behind_the_earth_are_dropped() -> None:
    observer = apogee_position()
    # One member just past Earth on the anti-observer side (LOS through Earth) + one to the side.
    behind = np.array([-APOGEE_RADIUS_M, 0.0, 0.0])
    beside = np.array([0.0, APOGEE_RADIUS_M, 0.0])
    los = line_of_sight_unit_vectors(observer, np.array([behind, beside]))
    assert len(los) == 1  # the occulted one is removed


def test_a_colocated_member_is_dropped_without_dividing_by_zero() -> None:
    observer = apogee_position()
    los = line_of_sight_unit_vectors(
        observer, np.array([observer, np.array([0.0, APOGEE_RADIUS_M, 0.0])])
    )
    assert len(los) == 1


def test_velocity_dop_falls_as_members_are_added() -> None:
    sigma_radial = 1.0
    few = transverse_velocity_sigma_m_s(constellation_los_units(4, "shell"), sigma_radial)
    many = transverse_velocity_sigma_m_s(constellation_los_units(20, "shell"), sigma_radial)
    assert many < few  # more members → tighter (lower DOP)


def test_coplanar_ring_leaves_the_normal_axis_unobservable_but_covers_transverse() -> None:
    ring = constellation_los_units(8, "ring")
    assert axis_observable(ring, TRANSVERSE_AXIS)  # the binding in-plane axis is observable
    assert not axis_observable(ring, NORMAL_AXIS)  # the out-of-plane axis is not
    # The transverse DOP is still finite for the ring.
    assert math.isfinite(velocity_dop(ring, TRANSVERSE_AXIS))


def test_shell_observes_all_three_axes() -> None:
    shell = constellation_los_units(12, "shell")
    assert axis_observable(shell, TRANSVERSE_AXIS)
    assert axis_observable(shell, NORMAL_AXIS)


def test_minimum_members_is_small_because_the_link_gives_ample_radial_precision() -> None:
    sigma_radial = carrier_phase_velocity_sigma_m_s(downlink_cn0_dbhz())
    n = min_members_for_target(sigma_radial, TARGET_TVEL_SIGMA_M_S, "shell")
    assert n is not None
    # Even a minimal 3-D-observable shell meets the C1 target — the match-not-beat headline.
    assert 3 <= n <= 6


def test_passive_receiver_is_a_negligible_mass_fraction_and_the_asic_is_not_the_driver() -> None:
    finding = apogee_nav_finding()
    assert finding.puffsat_mass_g < 100.0
    assert finding.mass_fraction < 0.01  # well under 1% of the 25 kg bus
    assert not finding.crypto_asic_is_mass_driver  # the crypto ASIC is sub-gram, not the driver


def test_finding_meets_the_c1_matching_target_with_margin() -> None:
    finding = apogee_nav_finding()
    assert finding.meets_target
    assert finding.meets_requirement
    assert finding.shell_transverse_velocity_sigma_m_s <= TARGET_TVEL_SIGMA_M_S
    assert finding.target_margin > 1.0
    # Match-not-beat: the requirement is looser than the target we design to.
    assert REQUIRED_TVEL_SIGMA_M_S > TARGET_TVEL_SIGMA_M_S


def test_finding_is_a_frozen_value_type() -> None:
    finding = apogee_nav_finding()
    assert isinstance(finding, ApogeeNavFinding)
    assert finding.components == PASSIVE_RECEIVER


def test_format_reports_link_velocity_gdop_and_mass() -> None:
    text = format_apogee_nav(apogee_nav_finding())
    assert "C/N0" in text
    assert "transverse" in text.lower()
    assert "mass driver" in text.lower()
    assert "match-not-beat" in text.lower() or "match" in text.lower()


# --- ring-vs-shell GDOP / min-member sweep -------------------------------------------------------
def test_sweep_has_one_point_per_swept_member_count() -> None:
    sigma_radial = carrier_phase_velocity_sigma_m_s(downlink_cn0_dbhz())
    points = gdop_sweep(sigma_radial)
    assert tuple(p.n_members for p in points) == DEFAULT_SWEEP_COUNTS


def test_sweep_transverse_sigma_falls_monotonically_with_members() -> None:
    sigma_radial = carrier_phase_velocity_sigma_m_s(downlink_cn0_dbhz())
    points = gdop_sweep(sigma_radial)
    shell = [p.shell_transverse_sigma_m_s for p in points]
    assert shell == sorted(shell, reverse=True)  # more members → tighter, no reversals


def test_sweep_marks_underdetermined_ring_at_n3_as_unobservable() -> None:
    # The coplanar ring at N=3 occults antipodally to 2 usable LOS — below the 3-LOS solve floor.
    sigma_radial = carrier_phase_velocity_sigma_m_s(downlink_cn0_dbhz())
    n3 = gdop_sweep(sigma_radial, (3,))[0]
    assert n3.ring_los_count == 2
    assert not n3.ring_transverse_observable
    assert math.isinf(n3.ring_transverse_sigma_m_s)
    # The shell at N=3 is well-posed (no antipodal member sits on the apse line).
    assert n3.shell_transverse_observable
    assert n3.shell_normal_observable


def test_sweep_ring_never_observes_the_orbit_normal_axis() -> None:
    sigma_radial = carrier_phase_velocity_sigma_m_s(downlink_cn0_dbhz())
    assert all(not p.ring_normal_observable for p in gdop_sweep(sigma_radial))


def test_finding_minimum_members_is_small_for_both_geometries() -> None:
    finding = gdop_sweep_finding()
    assert finding.min_members_shell == 3
    assert finding.min_members_ring == 4  # one more than shell: the N=3 antipodal occultation
    assert finding.shell_meets_target_at_min


def test_finding_ring_is_tighter_than_shell_on_the_binding_transverse_axis() -> None:
    # A coplanar ring pins the binding in-plane axis with every member; the shell spends members on
    # the (weak) orbit-normal axis — so the ring is tighter per member on transverse.
    finding = gdop_sweep_finding()
    assert finding.ring_transverse_advantage > 1.0


def test_finding_dop_follows_inverse_sqrt_n_diminishing_returns() -> None:
    finding = gdop_sweep_finding()
    assert finding.inverse_sqrt_n_scaling  # TDOP·√N flat on the asymptote
    assert 1.5 < finding.shell_tdop_sqrt_n < 2.5  # the ~2 DOP constant


def test_sweep_finding_is_a_frozen_value_type() -> None:
    finding = gdop_sweep_finding()
    assert isinstance(finding, GdopSweepFinding)


def test_format_sweep_reports_the_table_and_both_reads() -> None:
    text = format_gdop_sweep(gdop_sweep_finding())
    assert "ring vs shell" in text.lower()
    assert "minimum members" in text.lower()
    assert "1/√N" in text or "tighter" in text.lower()
    assert "match-not-beat" in text.lower()
