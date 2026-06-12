"""Integration test for the C3a executed terminal feedforward (live JVM, one descent set)."""

from __future__ import annotations

import pytest

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.runs.terminal import run_terminal_feedforward
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

pytestmark = pytest.mark.integration


def test_executed_feedforward_cancels_drag_within_the_measured_envelope() -> None:
    """The C3a slice end-to-end (ADR 0014 decisions 4/6), bounds set from the measured run.

    Measured 2026-06-11: pin 5.5 mm, drag displacement 0.085 m (far under ADR 0014's ~1–2 m
    estimate — drag concentrates in the final seconds, where it cannot integrate into
    position), executed residual 2 mm (rejection ~45×), Δv 0.0145 m/s ≈ B3a's 0.015 m/s,
    peak thrust 16 mN, peak slew 0.05 °/s.  Bounds are loosened ~an order so the test pins
    the physics, not the platform.
    """
    finding = run_terminal_feedforward()

    # Equivalence pin first (decision 4): the fixed-step Cowell terminal must reproduce the
    # proven adaptive-30 s descent before any burn is trusted.
    assert finding.equivalence_pin_m < 0.05
    assert abs(finding.equivalence_pin_toa_s) < 1e-3

    # The disease and the cure: uncompensated drag displaces the crossing at cm scale, and
    # the open-loop ZOH burn must reject the bulk of it.
    assert 0.01 < finding.drag_displacement_m < 1.0
    assert finding.executed_residual_m < finding.drag_displacement_m / 5.0
    assert finding.executed_residual_m < 0.02

    # The plan realizes B3a's profile: same Δv scale, one command per control-clock second.
    plan = finding.plan
    assert 0.005 < plan.dv_m_s < 0.05
    assert len(plan.commands) > 100

    # ADR 0004 actuator gates, now on the *executed* burn (B3a only measured the demand).
    assert not plan.saturated
    assert plan.peak_thrust_n < 0.1  # measured ~16 mN, cap 400 mN
    assert plan.peak_slew_rate_deg_s < 1.0
