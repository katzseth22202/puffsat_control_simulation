"""C3a terminal feedforward — pure ZOH burn planning over a sampled drag profile (no JVM).

ADR 0014 decision 6: C3a executes B3a's anti-drag feedforward as a real zero-order-hold
burn on the control clock.  The planning is pure: a sampled drag-acceleration history
(produced JVM-side by :mod:`puffsat_sim.montecarlo`) is held over each control step as a
thrust command opposing the drag.  The Orekit maneuver segments that execute the commands
live on the JVM side.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from puffsat_sim.anti_drag import PEAK_SLEW_LIMIT_DEG_S, PEAK_THRUST_LIMIT_N
from puffsat_sim.dispersion import Vec3
from puffsat_sim.propellant import propellant_curve


@dataclass(frozen=True)
class ThrustCommand:
    """One zero-order-hold burn segment: thrust held fixed over a control step."""

    start_s: float
    duration_s: float
    thrust_n: float
    direction: Vec3


@dataclass(frozen=True)
class FeedforwardPlan:
    """The ZOH realization of the anti-drag feedforward over the instrumented descent."""

    commands: tuple[ThrustCommand, ...]
    dv_m_s: float
    mass_kg: float
    saturated: bool
    peak_slew_rate_deg_s: float

    @property
    def peak_thrust_n(self) -> float:
        """The largest commanded thrust, post-cap (the ADR 0004 gate reads this)."""
        return max((cmd.thrust_n for cmd in self.commands), default=0.0)


def plan_feedforward(
    times_s: Sequence[float],
    drag_accel_m_s2: Sequence[Vec3],
    mass_kg: float,
    control_period_s: float,
    max_thrust_n: float = PEAK_THRUST_LIMIT_N,
) -> FeedforwardPlan:
    """Hold the sampled drag profile over each control step as an opposing thrust command.

    Each command starts on a control-clock tick and carries the thrust ``mass·|a_drag|``
    anti-parallel to the drag sampled at (or last before) the tick — the zero-order hold
    of the B3a profile — saturated at the actuator's ``max_thrust_n``.  The actuator holds
    the final tick's command until the end of the span, so the last command may be shorter
    than a control period (drag peaks there on a descent; dropping the tail would
    under-deliver the feedforward where it matters most).
    """
    accel = np.asarray(drag_accel_m_s2, dtype=np.float64).reshape(-1, 3)
    times = np.asarray(times_s, dtype=np.float64)

    commands: list[ThrustCommand] = []
    saturated = False
    start = float(times[0])
    span = float(times[-1] - times[0])
    # The -1e-9 keeps float noise in an exact multiple from spawning a ~zero-length step.
    steps = math.ceil(span / control_period_s - 1e-9)
    direction: Vec3 = (1.0, 0.0, 0.0)
    for k in range(steps):
        tick = start + k * control_period_s
        sample = int(np.searchsorted(times, tick, side="right")) - 1
        mag = float(np.linalg.norm(accel[sample]))
        if mag > 0.0:
            direction = (
                float(-accel[sample][0] / mag),
                float(-accel[sample][1] / mag),
                float(-accel[sample][2] / mag),
            )
        saturated = saturated or mass_kg * mag > max_thrust_n
        commands.append(
            ThrustCommand(
                start_s=tick,
                duration_s=min(control_period_s, start + span - tick),
                thrust_n=min(mass_kg * mag, max_thrust_n),
                direction=direction,
            )
        )
    dv = sum(cmd.thrust_n / mass_kg * cmd.duration_s for cmd in commands)
    return FeedforwardPlan(
        commands=tuple(commands),
        dv_m_s=dv,
        mass_kg=mass_kg,
        saturated=saturated,
        peak_slew_rate_deg_s=_peak_slew_rate_deg_s(commands),
    )


# Below this fraction of the plan's peak thrust the commanded direction is numerical
# noise (negligible drag), so direction changes there must not gate the slew rate.
_NEGLIGIBLE_THRUST_FRACTION: float = 1e-6


def _peak_slew_rate_deg_s(commands: Sequence[ThrustCommand]) -> float:
    floor = _NEGLIGIBLE_THRUST_FRACTION * max((cmd.thrust_n for cmd in commands), default=0.0)
    peak_rad_s = 0.0
    for prev, cmd in zip(commands, commands[1:], strict=False):
        if prev.thrust_n <= floor or cmd.thrust_n <= floor:
            continue
        cos_angle = (
            prev.direction[0] * cmd.direction[0]
            + prev.direction[1] * cmd.direction[1]
            + prev.direction[2] * cmd.direction[2]
        )
        angle = math.acos(min(1.0, max(-1.0, cos_angle)))
        peak_rad_s = max(peak_rad_s, angle / (cmd.start_s - prev.start_s))
    return math.degrees(peak_rad_s)


@dataclass(frozen=True)
class TerminalFeedforwardFinding:
    """The C3a executed-feedforward measurement set (pure container for the JVM numbers).

    Distances are 3D crossing-position separations at each trajectory's own 200 km
    event; "drag-free" is the same hand-off state descended with the drag perturbation
    removed, so ``executed_residual_m`` → 0 means the ZOH burn cancelled drag exactly.
    """

    plan: FeedforwardPlan
    equivalence_pin_m: float
    equivalence_pin_toa_s: float
    drag_displacement_m: float
    executed_residual_m: float
    executed_residual_toa_s: float


def format_terminal_feedforward(finding: TerminalFeedforwardFinding) -> str:
    """One-screen C3a report: pin, displacement vs executed residual, ADR 0004 gates, propellant."""
    plan = finding.plan
    span = (
        plan.commands[-1].start_s + plan.commands[-1].duration_s - plan.commands[0].start_s
        if plan.commands
        else 0.0
    )
    rejection = (
        finding.drag_displacement_m / finding.executed_residual_m
        if finding.executed_residual_m > 0.0
        else math.inf
    )
    thrust_verdict = "FAIL (saturated)" if plan.saturated else "PASS"
    slew_verdict = "PASS" if plan.peak_slew_rate_deg_s <= PEAK_SLEW_LIMIT_DEG_S else "FAIL"
    points = propellant_curve(plan.dv_m_s)
    curve = ", ".join(f"{p.fraction * 100:.4f}% @Isp{p.isp_s:.0f}" for p in points)
    budget = (
        "all within 2%"
        if all(p.within_budget for p in points)
        else "OVER 2% at Isp " + "/".join(f"{p.isp_s:.0f}" for p in points if not p.within_budget)
    )
    return "\n".join(
        [
            "C3a terminal feedforward — executed ZOH anti-drag burn (ADR 0014)",
            f"  Equivalence pin (fixed-step Cowell vs adaptive-30 s, unburned):"
            f" {finding.equivalence_pin_m:.4f} m, ToA {finding.equivalence_pin_toa_s:+.6f} s",
            f"  Drag displacement at crossing (unburned vs drag-free):"
            f" {finding.drag_displacement_m:.3f} m",
            f"  Executed residual (burned vs drag-free): {finding.executed_residual_m:.3f} m"
            f" → rejection {rejection:.1f}×",
            f"  ToA residual vs drag-free: {finding.executed_residual_toa_s:+.6f} s",
            f"  Plan: {len(plan.commands)} commands over {span:.1f} s, Δv {plan.dv_m_s:.6f} m/s",
            f"  Gates (ADR 0004): peak thrust {plan.peak_thrust_n * 1e3:.2f} mN"
            f" vs {PEAK_THRUST_LIMIT_N * 1e3:.0f} mN — {thrust_verdict}",
            f"                    peak slew {plan.peak_slew_rate_deg_s:.3f} °/s"
            f" vs {PEAK_SLEW_LIMIT_DEG_S:.1f} °/s — {slew_verdict}",
            f"  Propellant (fraction of wet mass): {curve} — {budget}",
        ]
    )
