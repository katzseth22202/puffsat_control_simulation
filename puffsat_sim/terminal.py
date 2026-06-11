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

from puffsat_sim.anti_drag import PEAK_THRUST_LIMIT_N
from puffsat_sim.dispersion import Vec3


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


def plan_feedforward(
    times_s: Sequence[float],
    drag_accel_m_s2: Sequence[Vec3],
    mass_kg: float,
    control_period_s: float,
    max_thrust_n: float = PEAK_THRUST_LIMIT_N,
) -> FeedforwardPlan:
    """Hold the sampled drag profile over each control step as an opposing thrust command.

    Each command starts on a control-clock tick, lasts one control period, and carries the
    thrust ``mass·|a_drag|`` anti-parallel to the drag sampled at (or last before) the tick
    — the zero-order hold of the B3a profile — saturated at the actuator's ``max_thrust_n``.
    """
    accel = np.asarray(drag_accel_m_s2, dtype=np.float64).reshape(-1, 3)
    times = np.asarray(times_s, dtype=np.float64)

    commands: list[ThrustCommand] = []
    saturated = False
    start = float(times[0])
    span = float(times[-1] - times[0])
    steps = int(span / control_period_s)
    for k in range(steps):
        tick = start + k * control_period_s
        sample = int(np.searchsorted(times, tick, side="right")) - 1
        mag = float(np.linalg.norm(accel[sample]))
        direction: Vec3 = (
            float(-accel[sample][0] / mag),
            float(-accel[sample][1] / mag),
            float(-accel[sample][2] / mag),
        )
        saturated = saturated or mass_kg * mag > max_thrust_n
        commands.append(
            ThrustCommand(
                start_s=tick,
                duration_s=control_period_s,
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
