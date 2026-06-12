"""Tests for the pure C3a terminal feedforward planning — synthetic series, no JVM."""

import math

import pytest

from puffsat_sim.terminal import (
    TerminalFeedforwardFinding,
    format_terminal_feedforward,
    plan_feedforward,
)


def test_zoh_holds_the_sample_at_or_last_before_each_control_tick() -> None:
    """With a finer sample grid, each command holds the profile value at its own tick."""
    times = [0.0, 0.4, 1.0, 2.0]
    accel = [(-0.001, 0.0, 0.0), (-0.002, 0.0, 0.0), (-0.003, 0.0, 0.0), (-0.004, 0.0, 0.0)]
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert [cmd.thrust_n for cmd in plan.commands] == [25.0 * 0.001, 25.0 * 0.003]


def test_plan_dv_is_the_commanded_impulse_over_wet_mass() -> None:
    """The plan's Δv is what the ZOH commands actually deliver: Σ (F/m)·dt, after the cap."""
    times = [0.0, 1.0, 2.0, 3.0]
    accel = [(-0.01, 0.0, 0.0)] * 4
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert plan.dv_m_s == pytest.approx(0.01 * 3.0)


def test_commands_are_clipped_at_the_actuator_thrust_limit() -> None:
    """Where the drag force exceeds the 400 mN actuator, the command saturates at the cap."""
    times = [0.0, 1.0]
    accel = [(-1.0, 0.0, 0.0)] * 2  # 1 m/s² × 25 kg = 25 N demanded, far over the cap
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert plan.commands[0].thrust_n == 0.4


def test_plan_reports_saturation_when_any_command_hits_the_cap() -> None:
    """Saturation means the feedforward cannot fully cancel drag — the report must flag it."""
    times = [0.0, 1.0, 2.0]
    weak = [(-0.001, 0.0, 0.0)] * 3
    strong = [(-1.0, 0.0, 0.0)] * 3
    assert not plan_feedforward(times, weak, mass_kg=25.0, control_period_s=1.0).saturated
    assert plan_feedforward(times, strong, mass_kg=25.0, control_period_s=1.0).saturated


def test_peak_slew_rate_is_the_fastest_turn_between_consecutive_commands() -> None:
    """The commanded thrust direction must not turn faster than the 1 °/s actuator loop."""
    times = [0.0, 1.0, 2.0, 3.0]
    accel = [
        (-0.01 * math.cos(math.radians(d)), -0.01 * math.sin(math.radians(d)), 0.0)
        for d in (0.0, 0.5, 1.0, 1.5)
    ]
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert plan.peak_slew_rate_deg_s == pytest.approx(0.5, rel=1e-6)


def test_near_zero_commands_do_not_register_spurious_slew() -> None:
    """Where drag is negligible the commanded direction is noise; it must not gate the slew."""
    times = [0.0, 1.0, 2.0, 3.0]
    accel = [
        (1e-12, 1e-12, 0.0),  # negligible, ~45° off — must be ignored
        (-0.01, 0.0, 0.0),
        (-0.01 * math.cos(math.radians(0.5)), -0.01 * math.sin(math.radians(0.5)), 0.0),
        (-0.01 * math.cos(math.radians(1.0)), -0.01 * math.sin(math.radians(1.0)), 0.0),
    ]
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert plan.peak_slew_rate_deg_s == pytest.approx(0.5, rel=1e-6)


def test_trailing_partial_step_is_covered_by_a_short_final_command() -> None:
    """The actuator holds the last tick's command until the crossing; the sub-period tail
    (where drag peaks on a descent) must not be silently dropped from the plan."""
    times = [0.0, 1.0, 2.0, 2.5]
    accel = [(-0.01, 0.0, 0.0)] * 4
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert [cmd.start_s for cmd in plan.commands] == [0.0, 1.0, 2.0]
    assert plan.commands[-1].duration_s == pytest.approx(0.5)
    assert plan.dv_m_s == pytest.approx(0.01 * 2.5)


def test_zero_drag_step_commands_zero_thrust_with_a_finite_direction() -> None:
    """An exactly-zero drag sample means engine off for that step — never a NaN direction."""
    times = [0.0, 1.0, 2.0]
    accel = [(0.0, 0.0, 0.0), (-0.01, 0.0, 0.0), (-0.01, 0.0, 0.0)]
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert plan.commands[0].thrust_n == 0.0
    assert all(math.isfinite(c) for c in plan.commands[0].direction)


def test_constant_drag_profile_becomes_one_zoh_command_per_control_step() -> None:
    """A constant drag profile turns into per-step commands opposing drag at m·|a|."""
    times = [0.0, 1.0, 2.0, 3.0]
    accel = [(-0.01, 0.0, 0.0)] * 4  # drag decelerating along -x
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert len(plan.commands) == 3  # three full steps cover the 3 s span
    for k, cmd in enumerate(plan.commands):
        assert cmd.start_s == float(k)
        assert cmd.duration_s == 1.0
        assert cmd.thrust_n == 25.0 * 0.01
        assert cmd.direction == (1.0, 0.0, 0.0)  # thrust opposes the drag acceleration


def _finding(
    plan_accel: float = 0.0005,
    drag_displacement_m: float = 2.0,
    executed_residual_m: float = 0.5,
    slew_deg_per_step: float = 0.0,
) -> TerminalFeedforwardFinding:
    times = [float(t) for t in range(5)]
    accel = [
        (
            -plan_accel * math.cos(math.radians(slew_deg_per_step * k)),
            -plan_accel * math.sin(math.radians(slew_deg_per_step * k)),
            0.0,
        )
        for k in range(5)
    ]
    return TerminalFeedforwardFinding(
        plan=plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0),
        equivalence_pin_m=0.003,
        equivalence_pin_toa_s=1e-5,
        drag_displacement_m=drag_displacement_m,
        executed_residual_m=executed_residual_m,
        executed_residual_toa_s=2e-5,
    )


def test_format_passes_both_actuator_gates_for_a_weak_smooth_plan() -> None:
    """Under-cap thrust and zero slew clear the ADR 0004 gates."""
    out = format_terminal_feedforward(_finding())
    assert out.count("PASS") == 2
    assert "FAIL" not in out


def test_format_fails_the_thrust_gate_when_the_plan_saturates() -> None:
    """A saturated plan means the actuator cannot match drag — the thrust gate must FAIL."""
    out = format_terminal_feedforward(_finding(plan_accel=1.0))
    assert "FAIL" in out
    assert "saturated" in out


def test_format_fails_the_slew_gate_when_direction_turns_faster_than_the_loop() -> None:
    """A commanded direction sweeping 2 °/s exceeds the ~1 °/s actuator loop."""
    out = format_terminal_feedforward(_finding(slew_deg_per_step=2.0))
    assert "FAIL" in out


def test_format_reports_the_rejection_ratio_and_the_pin() -> None:
    """The headline numbers: equivalence pin, displacement vs residual, rejection factor."""
    out = format_terminal_feedforward(_finding(drag_displacement_m=2.0, executed_residual_m=0.5))
    assert "4.0×" in out
    assert "0.003" in out  # the fixed-vs-adaptive pin distance


def test_format_reports_propellant_within_the_2_percent_budget() -> None:
    """The tiny anti-drag Δv sits far under the <2 % line at every Isp anchor."""
    out = format_terminal_feedforward(_finding())
    assert "within 2%" in out


def test_plan_peak_thrust_is_the_largest_command() -> None:
    """The gate reads the largest commanded thrust, post-cap."""
    times = [0.0, 1.0, 2.0]
    accel = [(-0.001, 0.0, 0.0), (-0.002, 0.0, 0.0), (-0.002, 0.0, 0.0)]
    plan = plan_feedforward(times, accel, mass_kg=25.0, control_period_s=1.0)
    assert plan.peak_thrust_n == pytest.approx(25.0 * 0.002)


def test_executed_plan_summarizes_a_command_history() -> None:
    """The C3b loop assembles its executed ZOH history into the same plan value type."""
    from puffsat_sim.terminal import ThrustCommand, executed_plan

    commands = (
        ThrustCommand(start_s=0.0, duration_s=1.0, thrust_n=0.25, direction=(1.0, 0.0, 0.0)),
        ThrustCommand(start_s=1.0, duration_s=1.0, thrust_n=0.4, direction=(0.0, 1.0, 0.0)),
    )
    plan = executed_plan(commands, mass_kg=25.0, saturated=True)

    assert plan.dv_m_s == pytest.approx((0.25 + 0.4) / 25.0)
    assert plan.peak_thrust_n == 0.4
    assert plan.saturated
    assert plan.peak_slew_rate_deg_s == pytest.approx(90.0)
