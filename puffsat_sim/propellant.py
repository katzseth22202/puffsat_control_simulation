"""B2 propellant ledger — the fraction-vs-Isp curve against the <2% claim (pure, no JVM).

Rung B's headline deliverable (§13 / ADR 0004, 0008): aggregate the mission Δv (B1 apogee
corrections + the B3a anti-drag burn) and report the propellant fraction it costs at each
Isp anchor, against the paper's "<2% of 25 kg" line.  The model is the ADR 0004 *linear*
relation ``fraction = Δv/(Isp·g₀)`` — the exact inverse of the A3 ``budget_dv_m_s`` — so the
ledger and the controllability budget share one 2% line.  (The dry mass cancels, so this is
mass-independent; for the exact propellant *mass* of a single burn use ``actuator.plan_burn``,
which differs only at the <0.5% Tsiolkovsky correction these small Δv never reach.)
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from puffsat_sim.constants import STANDARD_GRAVITY_M_S2

# The ADR 0004 propellant-budget anchors — defined alongside the A3 controllability budget
# (`sweep.budget_dv_m_s`) so the ledger and that budget share one Isp sweep and one 2% line.
from puffsat_sim.sweep import ISP_ANCHORS_S, PROPELLANT_BUDGET_FRACTION


@dataclass(frozen=True)
class PropellantPoint:
    """One Isp anchor on the ledger: the fraction it costs and whether it clears the budget."""

    isp_s: float
    fraction: float
    within_budget: bool


def propellant_fraction(dv_m_s: float, isp_s: float) -> float:
    """Propellant as a fraction of wet mass to deliver ``dv_m_s`` at ``isp_s`` (ADR 0004 linear)."""
    return dv_m_s / (isp_s * STANDARD_GRAVITY_M_S2)


def propellant_curve(
    dv_m_s: float,
    isp_anchors_s: Sequence[float] = ISP_ANCHORS_S,
    budget_fraction: float = PROPELLANT_BUDGET_FRACTION,
) -> tuple[PropellantPoint, ...]:
    """The fraction-vs-Isp curve for a total Δv, each anchor flagged against the <2% line."""
    return tuple(
        PropellantPoint(
            isp_s=isp,
            fraction=propellant_fraction(dv_m_s, isp),
            within_budget=propellant_fraction(dv_m_s, isp) <= budget_fraction,
        )
        for isp in isp_anchors_s
    )
