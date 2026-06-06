"""Pure-Python orbital mechanics helpers — no JVM dependency.

Used by truth_model.py and tested independently of Orekit.
"""
from __future__ import annotations

import math
from typing import Final

_EARTH_RADIUS_M: Final[float] = 6_378_137.0   # WGS84 equatorial radius [m]
_WGS84_MU: Final[float] = 3.986_004_418e14    # Earth gravitational parameter [m³/s²]


def keplerian_elements(perigee_alt_m: float, apogee_alt_m: float) -> tuple[float, float]:
    """Return (semi-major axis [m], eccentricity) from altitude above the surface."""
    r_p = _EARTH_RADIUS_M + perigee_alt_m
    r_a = _EARTH_RADIUS_M + apogee_alt_m
    a = (r_p + r_a) / 2.0
    e = (r_a - r_p) / (r_a + r_p)
    return a, e


def keplerian_period(semi_major_axis_m: float) -> float:
    """Return orbital period [s] from semi-major axis [m]."""
    return 2.0 * math.pi * math.sqrt(semi_major_axis_m**3 / _WGS84_MU)


def perigee_speed(semi_major_axis_m: float, perigee_alt_m: float) -> float:
    """Return speed at perigee [m/s] via the vis-viva equation."""
    r_p = _EARTH_RADIUS_M + perigee_alt_m
    return math.sqrt(_WGS84_MU * (2.0 / r_p - 1.0 / semi_major_axis_m))
