"""Foundational two-body orbital mechanics helpers — pure, force-agnostic.

Per-force analytic signatures (J2 rates, tidal ratios, SRP/drag) live with their
perturbation in :mod:`puffsat_sim.forces`; the great-circle orbital-plane builder
lives in :mod:`puffsat_sim.orbital_plane`.
"""
from __future__ import annotations

import math

from puffsat_sim.constants import EARTH_RADIUS_M, WGS84_MU


def keplerian_elements(perigee_alt_m: float, apogee_alt_m: float) -> tuple[float, float]:
    """Return (semi-major axis [m], eccentricity) from altitude above the surface."""
    r_p = EARTH_RADIUS_M + perigee_alt_m
    r_a = EARTH_RADIUS_M + apogee_alt_m
    a = (r_p + r_a) / 2.0
    e = (r_a - r_p) / (r_a + r_p)
    return a, e


def keplerian_period(semi_major_axis_m: float) -> float:
    """Return orbital period [s] from semi-major axis [m]."""
    return 2.0 * math.pi * math.sqrt(semi_major_axis_m**3 / WGS84_MU)


def wrap_to_pi(angle_rad: float) -> float:
    """Wrap an angle [rad] to (-π, π].

    Used to take the true signed difference of two angles that may straddle the
    0/2π branch cut — e.g. a small secular drift that Orekit reports as the raw
    difference ~2π − ε instead of ~−ε.
    """
    wrapped = (angle_rad + math.pi) % (2.0 * math.pi) - math.pi
    # (a + π) % 2π lands in [0, 2π); the subtraction yields [-π, π).  Map the
    # -π endpoint to +π so the interval is the conventional (-π, π].
    return math.pi if wrapped == -math.pi else wrapped


def perigee_speed(semi_major_axis_m: float, perigee_alt_m: float) -> float:
    """Return speed at perigee [m/s] via the vis-viva equation."""
    r_p = EARTH_RADIUS_M + perigee_alt_m
    return math.sqrt(WGS84_MU * (2.0 / r_p - 1.0 / semi_major_axis_m))
