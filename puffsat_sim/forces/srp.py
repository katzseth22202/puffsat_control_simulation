"""Solar radiation pressure perturbation — spec + analytic acceleration (pure)."""

from __future__ import annotations

from dataclasses import dataclass

from puffsat_sim.constants import SRP_P0_PA, SUN_MEAN_DISTANCE_M


@dataclass(frozen=True)
class SolarRadiation:
    """Cannonball solar radiation pressure.

    cr_area_over_mass is the lumped Cr·(A/m) [m²/kg]; the physical reflectivity
    Cr is folded in here so the Orekit model uses an effective coefficient of 1.
    """

    cr_area_over_mass: float


def srp_acceleration(
    cr_area_over_mass: float,
    sun_distance_m: float = SUN_MEAN_DISTANCE_M,
) -> float:
    """SRP acceleration magnitude [m/s²] at a given distance from the Sun.

    a_srp = P₀ · (d₀/r)² · (Cr·A/m)  where P₀ = 4.56×10⁻⁶ Pa at 1 AU.
    """
    return SRP_P0_PA * (SUN_MEAN_DISTANCE_M / sun_distance_m) ** 2 * cr_area_over_mass
