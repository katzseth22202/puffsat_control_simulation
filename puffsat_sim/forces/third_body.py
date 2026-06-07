"""Third-body perturbation — spec + analytic tidal-strength ratios (pure)."""
from __future__ import annotations

from dataclasses import dataclass

from puffsat_sim.constants import (
    EARTH_RADIUS_M,
    MOON_MEAN_DISTANCE_M,
    MOON_MU,
    SUN_MEAN_DISTANCE_M,
    SUN_MU,
    WGS84_MU,
)


@dataclass(frozen=True)
class ThirdBody:
    """Gravitational pull of the Sun and/or Moon on the PuffSat."""

    sun: bool = True
    moon: bool = True


def tidal_acceleration_ratio(
    apogee_alt_m: float, body_mu_m3_s2: float, body_distance_m: float
) -> float:
    """Ratio of third-body tidal acceleration to Earth monopole gravity at apogee.

    Uses the Hill (tidal) approximation: a_tidal ≈ 2·μ_body·r_apogee / d_body³.
    Ratio = a_tidal / a_earth = 2·μ_body·r_apogee³ / (d_body³·μ_earth).

    For the reference orbit (apogee 150 000 km) the Moon gives ~0.17% and
    the Sun gives ~0.08%, consistent with the design doc "~0.1%" benchmark.
    """
    r_a = EARTH_RADIUS_M + apogee_alt_m
    a_tidal = 2.0 * body_mu_m3_s2 * r_a / body_distance_m**3
    a_earth = WGS84_MU / r_a**2
    return a_tidal / a_earth


def lunar_tidal_ratio(apogee_alt_m: float) -> float:
    """Tidal acceleration ratio for the Moon at mean distance."""
    return tidal_acceleration_ratio(apogee_alt_m, MOON_MU, MOON_MEAN_DISTANCE_M)


def solar_tidal_ratio(apogee_alt_m: float) -> float:
    """Tidal acceleration ratio for the Sun at 1 AU."""
    return tidal_acceleration_ratio(apogee_alt_m, SUN_MU, SUN_MEAN_DISTANCE_M)
