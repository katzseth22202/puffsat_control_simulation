"""Unit tests for the pure torque-margin confirmation gate (ADR 0018)."""

from __future__ import annotations

import dataclasses
import math

from puffsat_sim.constants import EARTH_RADIUS_M
from puffsat_sim.mission import APOGEE_ALT_M, PERIGEE_ALT_M
from puffsat_sim.orbital_math import keplerian_elements, perigee_speed
from puffsat_sim.torque_margin import (
    MEASURED_DEMAND_SLEW_DEG_S,
    RAIL_SLEW_DEG_S,
    TorqueMarginFinding,
    aero_disturbance_torque_n_m,
    angular_accel_deg_s2,
    couple_torque_n_m,
    format_torque_margin,
    inertia_from_gyradius,
    perigee_los_rate_rad_s,
    reach_rate_time_s,
    torque_margin_finding,
)


def test_perigee_los_rate_is_vp_over_rp() -> None:
    rate = perigee_los_rate_rad_s(PERIGEE_ALT_M, APOGEE_ALT_M)
    a, _ = keplerian_elements(PERIGEE_ALT_M, APOGEE_ALT_M)
    expected = perigee_speed(a, PERIGEE_ALT_M) / (EARTH_RADIUS_M + PERIGEE_ALT_M)
    assert math.isclose(rate, expected)
    # The design-doc anchor: ~0.1 °/s near perigee.
    assert math.isclose(math.degrees(rate), 0.097, rel_tol=2e-2)


def test_inertia_from_gyradius_is_m_k_squared() -> None:
    assert math.isclose(inertia_from_gyradius(25.0, 0.4), 25.0 * 0.16)


def test_couple_torque_is_two_f_l() -> None:
    assert math.isclose(couple_torque_n_m(0.05, 0.5), 2.0 * 0.05 * 0.5)


def test_angular_accel_is_torque_over_inertia_in_degrees() -> None:
    assert math.isclose(angular_accel_deg_s2(0.05, 5.0), math.degrees(0.01))


def test_reach_rate_time_is_rate_over_accel() -> None:
    assert math.isclose(reach_rate_time_s(1.0, 0.5), 2.0)


def test_aero_disturbance_torque_is_force_times_arm() -> None:
    assert math.isclose(aero_disturbance_torque_n_m(0.0167, 0.15), 0.0167 * 0.15)


def test_default_gate_confirms_with_margin() -> None:
    f = torque_margin_finding()
    assert f.confirmed
    assert f.demand_within_rail
    assert f.reaches_demand_promptly
    assert f.holds_against_drag
    # ~10× rail margin over the perigee rate, matching design doc §13 line 344.
    assert math.isclose(f.rate_margin, 10.3, rel_tol=5e-2)


def test_measured_demand_carries_more_margin_than_the_perigee_estimate() -> None:
    f = torque_margin_finding()
    assert f.measured_rate_margin > f.rate_margin
    assert math.isclose(f.measured_rate_margin, RAIL_SLEW_DEG_S / MEASURED_DEMAND_SLEW_DEG_S)


def test_reach_times_order_demand_before_rail() -> None:
    f = torque_margin_finding()
    assert f.reach_demand_time_s < f.reach_rail_time_s
    assert f.reach_demand_time_s <= f.control_period_s


def test_disturbance_margin_is_control_over_disturbance() -> None:
    f = torque_margin_finding()
    assert math.isclose(f.disturbance_margin, f.control_torque_n_m / f.disturbance_torque_n_m)


def test_breakeven_inertia_reaches_demand_in_exactly_one_period() -> None:
    f = torque_margin_finding()
    at_breakeven = dataclasses.replace(f, inertia_kg_m2=f.breakeven_inertia_kg_m2)
    assert math.isclose(at_breakeven.reach_demand_time_s, f.control_period_s, rel_tol=1e-9)
    # The default inertia sits well under the break-even.
    assert f.inertia_kg_m2 < f.breakeven_inertia_kg_m2


def test_breakeven_control_torque_closes_the_disturbance_margin() -> None:
    f = torque_margin_finding()
    at_breakeven = dataclasses.replace(f, control_torque_n_m=f.breakeven_control_torque_n_m)
    assert math.isclose(at_breakeven.disturbance_margin, 1.0, rel_tol=1e-9)


def test_a_heavy_low_torque_case_fails_the_demand_agility() -> None:
    # A large inertia with a weak wheel cannot reach the demand rate within a period.
    f = TorqueMarginFinding(
        demand_slew_deg_s=0.097,
        measured_demand_slew_deg_s=MEASURED_DEMAND_SLEW_DEG_S,
        rail_slew_deg_s=RAIL_SLEW_DEG_S,
        inertia_kg_m2=200.0,
        control_torque_n_m=1.0e-3,
        disturbance_torque_n_m=2.5e-3,
    )
    assert not f.reaches_demand_promptly
    assert not f.holds_against_drag
    assert not f.confirmed


def test_format_reports_demand_actuator_and_verdict() -> None:
    text = format_torque_margin(torque_margin_finding())
    assert "Torque-margin confirmation" in text
    assert "CONFIRMED" in text
    assert "Break-even" in text
    assert "rail" in text
