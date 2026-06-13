"""Unit tests for the pure C4 control-loop latency core (ADR 0014, §16.8)."""

from __future__ import annotations

import math

from puffsat_sim.latency import (
    DEGRADATION_FACTOR,
    MIDCOURSE_LOOP,
    SPEED_OF_LIGHT_M_S,
    TERMINAL_LOOP,
    ControlLoop,
    DeadTimeBuffer,
    LatencyFinding,
    LatencySource,
    TauSweepPoint,
    comms_one_way_s,
    fly_terminal_loop,
    format_latency,
    latency_finding,
    tau_sweep,
)


def test_comms_one_way_is_light_time() -> None:
    assert comms_one_way_s(SPEED_OF_LIGHT_M_S) == 1.0
    # The §16.8 anchor: 2000 km one-way ≈ 6.7 ms.
    assert math.isclose(comms_one_way_s(2_000_000.0), 6.671e-3, rel_tol=1e-3)


def test_loop_tau_sums_its_sources() -> None:
    loop = ControlLoop(
        "x", 1.0, (LatencySource("a", 1.0e-3), LatencySource("b", 5.0e-3), LatencySource("c", 1e-3))
    )
    assert math.isclose(loop.tau_s, 7.0e-3, rel_tol=1e-12)


def test_phase_loss_is_omega_c_tau_in_degrees() -> None:
    loop = ControlLoop("x", 1.0, (LatencySource("a", 7.333564e-3),))
    assert math.isclose(
        loop.phase_margin_loss_deg, math.degrees(2.0 * math.pi * 1.0 * 7.333564e-3), rel_tol=1e-12
    )
    # The measured terminal budget at 1 Hz lands at single-digit degrees (~2.64°).
    assert TERMINAL_LOOP.phase_margin_loss_deg is not None
    assert math.isclose(TERMINAL_LOOP.phase_margin_loss_deg, 2.64, abs_tol=0.05)


def test_discrete_loop_has_no_phase_to_erode() -> None:
    # The midcourse replan is impulsive — no continuous loop, so no phase-margin loss.
    assert MIDCOURSE_LOOP.bandwidth_hz is None
    assert MIDCOURSE_LOOP.phase_margin_loss_deg is None
    # But it still carries the larger dead-time (comms round-trip dominated).
    assert MIDCOURSE_LOOP.tau_s > TERMINAL_LOOP.tau_s


def test_dead_time_buffer_zero_ticks_is_identity() -> None:
    buffer = DeadTimeBuffer(0)
    assert buffer.step((1.0, 2.0, 3.0)) == (1.0, 2.0, 3.0)
    assert buffer.step((4.0, 5.0, 6.0)) == (4.0, 5.0, 6.0)


def test_dead_time_buffer_holds_then_lags_by_delay() -> None:
    buffer = DeadTimeBuffer(2)
    # The line is not yet full — the loop acts on the first value it ever saw.
    assert buffer.step((1.0, 0.0, 0.0)) == (1.0, 0.0, 0.0)
    assert buffer.step((2.0, 0.0, 0.0)) == (1.0, 0.0, 0.0)
    assert buffer.step((3.0, 0.0, 0.0)) == (1.0, 0.0, 0.0)
    # Now full: the acted value lags the input by exactly 2 steps.
    assert buffer.step((4.0, 0.0, 0.0)) == (2.0, 0.0, 0.0)
    assert buffer.step((5.0, 0.0, 0.0)) == (3.0, 0.0, 0.0)


def test_noiseless_loop_with_no_offset_barely_moves() -> None:
    # No funnel-entry error and no noise → no command, so the miss stays at the origin.
    assert fly_terminal_loop(0.0, 0.0, rng=None) == 0.0
    assert fly_terminal_loop(0.0, 20.0, rng=None) == 0.0


def test_dead_time_degrades_the_homing_miss() -> None:
    # A large dead-time (many ticks stale) wrecks the miss; the budget-scale delay does not.
    baseline = fly_terminal_loop(400.0, 0.0, rng=None)
    assert fly_terminal_loop(400.0, 20.0, rng=None) > 5.0 * baseline


def test_sub_tick_latency_reproduces_the_zero_delay_miss() -> None:
    # At a 1 s control period the ms-class budget rounds to 0 ticks — the identity buffer —
    # so the miss is byte-identical to the zero-delay run (the budget is structurally invisible).
    points = tau_sweep((0.0, 0.0073, 0.4), entry_m=400.0, dt_s=1.0)
    assert all(p.delay_ticks == 0 for p in points)
    assert points[0].lateral_m == points[1].lateral_m == points[2].lateral_m


def test_sweep_delay_ticks_round_to_the_cadence() -> None:
    points = tau_sweep((1.0, 2.0, 5.0), entry_m=400.0, dt_s=1.0)
    assert [p.delay_ticks for p in points] == [1, 2, 5]


def _finding() -> LatencyFinding:
    # A synthetic sweep: flat through 1 s, doubling (degrading) at 2 s.
    sweep = (
        TauSweepPoint(0.0, 0, 1.34),
        TauSweepPoint(0.0073, 0, 1.34),
        TauSweepPoint(1.0, 1, 1.01),
        TauSweepPoint(2.0, 2, 3.38),
        TauSweepPoint(5.0, 5, 10.64),
    )
    return LatencyFinding(
        loops=(TERMINAL_LOOP, MIDCOURSE_LOOP),
        sweep=sweep,
        sweep_cadence_hz=1.0,
        sweep_entry_m=400.0,
    )


def test_finding_reads_the_budget_and_baseline() -> None:
    f = _finding()
    assert f.terminal_budget_tau_s == TERMINAL_LOOP.tau_s
    assert f.control_period_s == 1.0
    assert f.baseline_lateral_m == 1.34  # the zero-delay point


def test_tolerated_and_breakdown_bracket_the_degradation() -> None:
    f = _finding()
    # 1.01 m (τ=1 s) is within 2× of 1.34; 3.38 m (τ=2 s) is not.
    assert f.tolerated_latency_s == 1.0
    assert f.breakdown_latency_s == 2.0
    # ~137× the ms-class budget is absorbed before the miss doubles.
    assert math.isclose(f.budget_margin, 1.0 / TERMINAL_LOOP.tau_s, rel_tol=1e-12)
    assert f.budget_margin > 100.0


def test_breakdown_is_none_when_no_swept_tau_degrades() -> None:
    flat = LatencyFinding(
        loops=(TERMINAL_LOOP, MIDCOURSE_LOOP),
        sweep=(TauSweepPoint(0.0, 0, 1.0), TauSweepPoint(1.0, 1, 1.5)),
        sweep_cadence_hz=1.0,
        sweep_entry_m=400.0,
    )
    assert flat.breakdown_latency_s is None
    assert flat.tolerated_latency_s == 1.0


def test_degradation_factor_is_the_doubling_threshold() -> None:
    assert DEGRADATION_FACTOR == 2.0


def test_latency_finding_runner_carries_baseline_and_budget() -> None:
    f = latency_finding(extra_latencies_s=(1.0, 2.0))
    latencies = [p.latency_s for p in f.sweep]
    assert 0.0 in latencies  # the baseline is always present
    assert any(math.isclose(t, TERMINAL_LOOP.tau_s, rel_tol=1e-9) for t in latencies)
    assert latencies == sorted(latencies)


def test_format_reports_budget_phase_and_rung_d_deferral() -> None:
    text = format_latency(_finding())
    assert "C4" in text
    assert "ω_c·τ" in text
    assert "discrete replan" in text  # the midcourse loop's framing
    assert "sub-tick" in text
    assert "flat" in text
    assert "DEGRADES" in text  # the τ=2 s point crosses the doubling threshold
    assert "Rung D" in text  # the deferred combined-stress note
