"""B1 finite-burn actuator — pure Tsiolkovsky burn kinematics (no JVM).

Maps a commanded impulsive Δv (the unchanged A1 corrector's output, ADR 0003) to a
finite, mass-depleting burn: the ADR 0004 proportional cold-gas thruster (400 mN,
Isp sweep, Tsiolkovsky depletion).  Pure so the kinematics are unit-tested without a
JVM; the Orekit ``ConstantThrustManeuver`` that executes the burn lives on the JVM
side (:mod:`puffsat_sim.montecarlo`).

The burn the propagator executes is *constant-mass* and therefore Isp-free — its
duration is Δv·m/F (ADR 0008 mass convention: scale thrust to the propagator's
fictitious 1 kg, keep the lumped drag/SRP coefficients).  Isp enters only the
reported propellant (the ADR 0004 decision-2 transform: Δv computed once, propellant
swept across Isp at B2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from puffsat_sim.constants import STANDARD_GRAVITY_M_S2


@dataclass(frozen=True)
class Actuator:
    """ADR 0004 cold-gas actuator: max thrust, Isp, and PuffSat wet mass."""

    isp_s: float
    max_thrust_n: float = 0.4  # 400 mN (ADR 0004)
    wet_mass_kg: float = 25.0  # PuffSat mass (paper; the <2%-of-25 kg claim)


@dataclass(frozen=True)
class FiniteBurn:
    """The finite realization of a commanded Δv: firing time and propellant spent."""

    duration_s: float
    propellant_kg: float


def plan_burn(actuator: Actuator, dv_rtn: tuple[float, float, float]) -> FiniteBurn:
    """Map a commanded Δv (RTN) to its finite-burn realization.

    Constant-mass duration Δv·m/F: the executed burn delivers the commanded Δv at the
    full thrust, so its duration is Isp-free (Isp only scales the propellant ledger).
    Propellant is the Tsiolkovsky mass at the actuator's Isp (the ADR 0004 ledger).
    """
    dv = (dv_rtn[0] ** 2 + dv_rtn[1] ** 2 + dv_rtn[2] ** 2) ** 0.5
    duration_s = dv * actuator.wet_mass_kg / actuator.max_thrust_n
    propellant_kg = actuator.wet_mass_kg * (
        1.0 - math.exp(-dv / (actuator.isp_s * STANDARD_GRAVITY_M_S2))
    )
    return FiniteBurn(duration_s=duration_s, propellant_kg=propellant_kg)
