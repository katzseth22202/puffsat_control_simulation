"""Geopotential perturbation — spec + analytic J2 secular rates (pure)."""

from __future__ import annotations

import math
from dataclasses import dataclass

from puffsat_sim.constants import EARTH_RADIUS_M, J2, WGS84_MU


@dataclass(frozen=True)
class Geopotential:
    """Spherical-harmonic Earth gravity field beyond the point-mass monopole.

    degree is the maximum harmonic degree.  order is the maximum tesseral order
    and must not exceed degree; order=0 selects a zonal-only field, so
    degree=2, order=0 is pure J2 (no C21/S21/C22/S22 sectorial terms).
    """

    degree: int
    order: int = 0

    def __post_init__(self) -> None:
        if self.order > self.degree:
            raise ValueError(
                f"geopotential order ({self.order}) cannot exceed degree ({self.degree})."
            )


def j2_nodal_regression_rate(
    semi_major_axis_m: float, eccentricity: float, inclination_rad: float
) -> float:
    """First-order secular RAAN drift rate due to J2 [rad/s].

    dΩ/dt = -3/2 · n · J2 · (Rₑ/p)² · cos i
    """
    n = math.sqrt(WGS84_MU / semi_major_axis_m**3)
    p = semi_major_axis_m * (1.0 - eccentricity**2)
    return -1.5 * n * J2 * (EARTH_RADIUS_M / p) ** 2 * math.cos(inclination_rad)


def j2_apsidal_precession_rate(
    semi_major_axis_m: float, eccentricity: float, inclination_rad: float
) -> float:
    """First-order secular argument-of-perigee drift rate due to J2 [rad/s].

    dω/dt = 3/4 · n · J2 · (Rₑ/p)² · (5 cos²i − 1)
    """
    n = math.sqrt(WGS84_MU / semi_major_axis_m**3)
    p = semi_major_axis_m * (1.0 - eccentricity**2)
    return 0.75 * n * J2 * (EARTH_RADIUS_M / p) ** 2 * (5.0 * math.cos(inclination_rad) ** 2 - 1.0)
