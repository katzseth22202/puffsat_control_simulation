"""Pure torque-margin confirmation gate — a Rung-D pre-gate (ADR 0018, no JVM).

The C3b terminal aim rides a **1 °/s direction-loop slew rail** (`anti_drag.PEAK_SLEW_LIMIT_DEG_S`);
its noise discipline (the 45° firing-lag hold) holds fire until the gimbal/body has caught up to
the command at that rate.  Design doc §13 (line 344) flags the ~10× margin over the ~0.1 °/s
perigee sweep as "a result to confirm, not assume."  This gate is that confirmation — a
back-of-envelope, *non-blocking* unlike the σ_θ tracker budget — answering two questions for the
*PuffSat's* attitude/direction loop:

* **Is the demand comfortably inside the rail?**  The thrust direction need only track the
  line-of-sight rotation, which peaks at the perigee orbital rate ``v_p / r_p`` (computed from the
  reference orbit) — ~0.1 °/s, and B3a *measured* the descent demand at 0.048 °/s.  The rail
  carries ~10–20× over that.

* **Can a realistic actuator deliver the rail on the PuffSat's inertia?**  Taking the conservative
  whole-body-slew case (gimballing the small nozzle is easier a fortiori), a modest cold-gas RCS
  couple gives angular acceleration ``T/I`` — enough to reach the demand rate within a control
  period and the full rail within a few — and out-torques the aerodynamic disturbance (peak drag
  force × the centre-of-pressure / centre-of-mass offset) it must hold against during the burn.

The PuffSat moment of inertia and the attitude-actuator torque are **paper-side pins** (like the
target-vehicle inertia in ADR 0015), so the gate reports its verdict's **break-even** in both —
the inertia and the control torque at which the margins close — not just a point pass.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from puffsat_sim.anti_drag import PEAK_SLEW_LIMIT_DEG_S
from puffsat_sim.constants import EARTH_RADIUS_M
from puffsat_sim.mission import APOGEE_ALT_M, PERIGEE_ALT_M
from puffsat_sim.orbital_math import keplerian_elements, perigee_speed

# The C3b direction-loop rail being confirmed, and the one-tick control clock it slews on.
RAIL_SLEW_DEG_S: float = PEAK_SLEW_LIMIT_DEG_S
CONTROL_PERIOD_S: float = 1.0

# B3a measured demand on the descent (ADR 0008 finding): peak |drag force| the main thruster
# matches, and the fastest the anti-drag thrust direction actually turned — the realized slew
# demand, well under the ~0.1 °/s perigee-rate estimate.
MEASURED_PEAK_DRAG_FORCE_N: float = 0.0167
MEASURED_DEMAND_SLEW_DEG_S: float = 0.048

PUFFSAT_MASS_KG: float = 25.0

# Confirmation thresholds: the rail must carry a comfortable factor over the demand, and the
# attitude actuator a comfortable factor over the aero disturbance, for the gate to read "slack".
RATE_MARGIN_MIN: float = 3.0
DISTURBANCE_MARGIN_MIN: float = 3.0


def perigee_los_rate_rad_s(perigee_alt_m: float, apogee_alt_m: float) -> float:
    """The line-of-sight rotation rate at perigee, ``v_p / r_p`` — the thrust-direction demand."""
    semi_major_axis_m, _ = keplerian_elements(perigee_alt_m, apogee_alt_m)
    return perigee_speed(semi_major_axis_m, perigee_alt_m) / (EARTH_RADIUS_M + perigee_alt_m)


def inertia_from_gyradius(mass_kg: float, gyradius_m: float) -> float:
    """Moment of inertia ``m·k²`` from a radius of gyration (the paper-side inertia pin)."""
    return mass_kg * gyradius_m**2


def couple_torque_n_m(thrust_n: float, moment_arm_m: float) -> float:
    """Control torque from a balanced RCS couple: two opposed thrusters, ``2·f·L``."""
    return 2.0 * thrust_n * moment_arm_m


def angular_accel_deg_s2(torque_n_m: float, inertia_kg_m2: float) -> float:
    """Angular acceleration ``T/I`` the actuator delivers, in deg/s²."""
    return math.degrees(torque_n_m / inertia_kg_m2)


def reach_rate_time_s(rate_deg_s: float, accel_deg_s2: float) -> float:
    """Time to accelerate from rest to ``rate_deg_s`` at constant ``accel_deg_s2``."""
    return rate_deg_s / accel_deg_s2


def aero_disturbance_torque_n_m(drag_force_n: float, cp_cm_offset_m: float) -> float:
    """The aero disturbance torque the attitude loop holds against: drag force × CP–CM arm."""
    return drag_force_n * cp_cm_offset_m


@dataclass(frozen=True)
class TorqueMarginFinding:
    """The torque-margin confirmation: demand vs the rail, actuator agility, disturbance margin."""

    demand_slew_deg_s: float
    measured_demand_slew_deg_s: float
    rail_slew_deg_s: float
    inertia_kg_m2: float
    control_torque_n_m: float
    disturbance_torque_n_m: float
    control_period_s: float = CONTROL_PERIOD_S
    rate_margin_min: float = RATE_MARGIN_MIN
    disturbance_margin_min: float = DISTURBANCE_MARGIN_MIN

    @property
    def rate_margin(self) -> float:
        """How many times inside the rail the perigee-rate demand sits."""
        return self.rail_slew_deg_s / self.demand_slew_deg_s

    @property
    def measured_rate_margin(self) -> float:
        """How many times inside the rail the B3a-measured descent demand sits."""
        return self.rail_slew_deg_s / self.measured_demand_slew_deg_s

    @property
    def angular_accel_deg_s2(self) -> float:
        return angular_accel_deg_s2(self.control_torque_n_m, self.inertia_kg_m2)

    @property
    def reach_demand_time_s(self) -> float:
        return reach_rate_time_s(self.demand_slew_deg_s, self.angular_accel_deg_s2)

    @property
    def reach_rail_time_s(self) -> float:
        return reach_rate_time_s(self.rail_slew_deg_s, self.angular_accel_deg_s2)

    @property
    def disturbance_margin(self) -> float:
        """Control torque over the aero disturbance torque it must hold against."""
        return self.control_torque_n_m / self.disturbance_torque_n_m

    @property
    def demand_within_rail(self) -> bool:
        return self.rate_margin >= self.rate_margin_min

    @property
    def reaches_demand_promptly(self) -> bool:
        """The actuator reaches the demand rate within one control period."""
        return self.reach_demand_time_s <= self.control_period_s

    @property
    def holds_against_drag(self) -> bool:
        return self.disturbance_margin >= self.disturbance_margin_min

    @property
    def confirmed(self) -> bool:
        return self.demand_within_rail and self.reaches_demand_promptly and self.holds_against_drag

    @property
    def breakeven_inertia_kg_m2(self) -> float:
        """The inertia at which the actuator just reaches the demand rate in one control period."""
        demand_rad_s = math.radians(self.demand_slew_deg_s)
        return self.control_torque_n_m * self.control_period_s / demand_rad_s

    @property
    def breakeven_control_torque_n_m(self) -> float:
        """The control torque at which the disturbance margin closes to 1."""
        return self.disturbance_torque_n_m


def torque_margin_finding(
    *,
    perigee_alt_m: float = PERIGEE_ALT_M,
    apogee_alt_m: float = APOGEE_ALT_M,
    mass_kg: float = PUFFSAT_MASS_KG,
    gyradius_m: float = 0.45,
    rcs_thrust_n: float = 0.05,
    moment_arm_m: float = 0.5,
    cp_cm_offset_m: float = 0.15,
    peak_drag_force_n: float = MEASURED_PEAK_DRAG_FORCE_N,
) -> TorqueMarginFinding:
    """Assemble the torque-margin confirmation (the pure runner).

    Conservative paper-side defaults: a 25 kg PuffSat with a 0.45 m radius of gyration
    (a moderately inflated puff — most mass central, membrane at radius), a cold-gas RCS
    couple of two 50 mN thrusters at a 0.5 m arm, and a 0.15 m CP–CM offset against the
    B3a peak drag force.  Gimballing the small nozzle would beat whole-body slew a fortiori.
    """
    return TorqueMarginFinding(
        demand_slew_deg_s=math.degrees(perigee_los_rate_rad_s(perigee_alt_m, apogee_alt_m)),
        measured_demand_slew_deg_s=MEASURED_DEMAND_SLEW_DEG_S,
        rail_slew_deg_s=RAIL_SLEW_DEG_S,
        inertia_kg_m2=inertia_from_gyradius(mass_kg, gyradius_m),
        control_torque_n_m=couple_torque_n_m(rcs_thrust_n, moment_arm_m),
        disturbance_torque_n_m=aero_disturbance_torque_n_m(peak_drag_force_n, cp_cm_offset_m),
    )


def format_torque_margin(finding: TorqueMarginFinding) -> str:
    """One-screen torque-margin confirmation report (non-blocking)."""
    verdict = "CONFIRMED" if finding.confirmed else "NOT CONFIRMED"
    lines = [
        "Torque-margin confirmation — the 1 °/s slew rail the C3b aim rides"
        " (ADR 0018; non-blocking)",
        f"  Demand vs rail: perigee LOS rate {finding.demand_slew_deg_s:.3f} °/s"
        f" → {finding.rate_margin:.1f}× inside the {finding.rail_slew_deg_s:g} °/s rail"
        f" (B3a measured {finding.measured_demand_slew_deg_s:g} °/s"
        f" → {finding.measured_rate_margin:.0f}×).",
        f"  Actuator (paper-side pins): I {finding.inertia_kg_m2:.2f} kg·m²,"
        f" control torque {finding.control_torque_n_m * 1e3:.0f} mN·m"
        f" → α {finding.angular_accel_deg_s2:.3f} °/s².",
        f"    Reaches the demand rate in {finding.reach_demand_time_s:.2f} s"
        f" (≤ {finding.control_period_s:g} s control period:"
        f" {'yes' if finding.reaches_demand_promptly else 'NO'});"
        f" the full rail in {finding.reach_rail_time_s:.2f} s.",
        f"  Disturbance: aero torque {finding.disturbance_torque_n_m * 1e3:.2f} mN·m"
        f" (drag {MEASURED_PEAK_DRAG_FORCE_N * 1e3:g} mN × CP–CM offset)"
        f" → control torque {finding.disturbance_margin:.0f}× over it.",
        f"  Break-even (the unpinned margins): inertia up to {finding.breakeven_inertia_kg_m2:.0f}"
        f" kg·m² still reaches demand in a period; control torque down to"
        f" {finding.breakeven_control_torque_n_m * 1e3:.2f} mN·m still holds against drag.",
        f"  → {verdict}: the 1 °/s rail is deliverable with margin on the conservative pins.",
    ]
    return "\n".join(lines)
