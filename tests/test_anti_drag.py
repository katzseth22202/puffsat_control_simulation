"""Tests for the pure B3a anti-drag profile reductions — synthetic series, no JVM."""

import math

import pytest

from puffsat_sim.anti_drag import summarize_anti_drag


def test_anti_drag_dv_integrates_acceleration_magnitude_over_time() -> None:
    """Anti-drag Δv is ∫|a_drag| dt — the velocity the burn must add back to cancel drag."""
    times = [0.0, 50.0, 100.0]
    accel = [(-0.01, 0.0, 0.0)] * 3  # constant 0.01 m/s² for 100 s → 1.0 m/s
    profile = summarize_anti_drag(times, accel, mass_kg=25.0)
    assert profile.anti_drag_dv_m_s == pytest.approx(0.01 * 100.0)


def test_peak_thrust_is_max_drag_force_the_burn_must_match() -> None:
    """Peak required thrust is the largest drag force (max|a_drag|·mass) over the descent."""
    times = [0.0, 1.0, 2.0]
    accel = [(-0.002, 0.0, 0.0), (-0.01, 0.0, 0.0), (-0.004, 0.0, 0.0)]
    profile = summarize_anti_drag(times, accel, mass_kg=25.0)
    assert profile.peak_thrust_n == pytest.approx(0.01 * 25.0)


def test_peak_slew_rate_is_the_fastest_sweep_of_the_anti_drag_direction() -> None:
    """Peak slew rate is the max angular rate of the thrust direction (how fast it must turn)."""
    times = [0.0, 1.0, 2.0]
    # constant |a|, direction rotating 0.5° per second in the xy-plane → 0.5 °/s throughout.
    accel = [
        (0.01 * math.cos(math.radians(d)), 0.01 * math.sin(math.radians(d)), 0.0)
        for d in (0.0, 0.5, 1.0)
    ]
    profile = summarize_anti_drag(times, accel, mass_kg=25.0)
    assert profile.peak_slew_rate_deg_s == pytest.approx(0.5, rel=1e-6)


def test_duration_is_the_sampled_descent_span() -> None:
    times = [10.0, 30.0, 70.0]
    accel = [(-0.01, 0.0, 0.0)] * 3
    profile = summarize_anti_drag(times, accel, mass_kg=25.0)
    assert profile.duration_s == pytest.approx(60.0)


def test_near_zero_drag_samples_do_not_register_spurious_slew() -> None:
    """Where drag is negligible the direction is meaningless; it must not count as fast slew."""
    times = [0.0, 1.0, 2.0]
    accel = [
        (1e-12, 1e-12, 0.0),  # negligible drag, ~45° direction — must be ignored
        (0.01, 0.0, 0.0),
        (0.01 * math.cos(math.radians(0.5)), 0.01 * math.sin(math.radians(0.5)), 0.0),
    ]
    profile = summarize_anti_drag(times, accel, mass_kg=25.0)
    assert profile.peak_slew_rate_deg_s == pytest.approx(0.5, rel=1e-6)
