"""B3a anti-drag profile — pure reductions over a sampled descent (no JVM).

The §13/ADR 0008 B3 deliverable is a *feedforward cost baseline*: instrument the
known-drag descent 600 → 200 km and report what an anti-drag burn would have to
deliver — the Δv to cancel drag, the peak thrust it demands of the actuator, and how
fast its direction sweeps — then check those against the ADR 0004 actuator limits
(400 mN, ~1°/s) and the paper's GOCE-ANFO estimate (`sec:estimate_cold_gas`, 374 g /
400 mN).  The numbers come from sampling the truth descent on the JVM side
(:mod:`puffsat_sim.montecarlo`); the reductions here are pure so they unit-test
against synthetic series without booting Orekit.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from puffsat_sim.dispersion import Vec3

# Actuator limits the B3a feasibility check is read against (ADR 0004 / paper §13).
PEAK_THRUST_LIMIT_N: float = 0.4  # 400 mN cold-gas max
PEAK_SLEW_LIMIT_DEG_S: float = 1.0  # direction-loop rate limit


@dataclass(frozen=True)
class AntiDragProfile:
    """What an anti-drag burn must deliver over the instrumented descent."""

    anti_drag_dv_m_s: float
    peak_thrust_n: float
    peak_slew_rate_deg_s: float
    duration_s: float


def summarize_anti_drag(
    times_s: Sequence[float],
    drag_accel_m_s2: Sequence[Vec3],
    mass_kg: float,
    slew_floor_frac: float = 0.05,
) -> AntiDragProfile:
    """Reduce a sampled drag-acceleration history to the anti-drag burn requirement.

    ``anti_drag_dv`` is the trapezoidal ∫|a_drag| dt — the speed the burn must add back.
    ``peak_thrust`` is the largest drag force (max|a_drag|·mass) the actuator must match.
    ``peak_slew_rate`` is the fastest the thrust direction (anti-parallel to drag) must
    turn; it is measured only where the drag exceeds ``slew_floor_frac`` of its peak, since
    the direction is ill-defined where drag is negligible (the burn is effectively off).
    """
    accel = np.asarray(drag_accel_m_s2, dtype=np.float64).reshape(-1, 3)
    times = np.asarray(times_s, dtype=np.float64)
    if times.size == 0:
        return AntiDragProfile(0.0, 0.0, 0.0, 0.0)

    mag = np.linalg.norm(accel, axis=1)
    anti_drag_dv = float(np.trapezoid(mag, times))
    peak_thrust = float(np.max(mag)) * mass_kg
    duration = float(times[-1] - times[0])

    floor = slew_floor_frac * float(np.max(mag))
    peak_slew_rad_s = 0.0
    for i in range(len(mag) - 1):
        dt = float(times[i + 1] - times[i])
        if dt <= 0.0 or mag[i] < floor or mag[i + 1] < floor:
            continue
        u0 = accel[i] / mag[i]
        u1 = accel[i + 1] / mag[i + 1]
        angle = float(np.arccos(np.clip(np.dot(u0, u1), -1.0, 1.0)))
        peak_slew_rad_s = max(peak_slew_rad_s, angle / dt)

    return AntiDragProfile(
        anti_drag_dv_m_s=anti_drag_dv,
        peak_thrust_n=peak_thrust,
        peak_slew_rate_deg_s=math.degrees(peak_slew_rad_s),
        duration_s=duration,
    )
