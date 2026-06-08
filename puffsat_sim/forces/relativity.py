"""Relativistic (Schwarzschild) perturbation — spec + analytic apsidal advance (pure).

The Schwarzschild post-Newtonian correction to Earth's point-mass gravity is
conservative and tiny (~1e-9 of the monopole acceleration), but it accumulates as
a *prograde apsidal precession*.  On this high-eccentricity orbit that precession
reaches ~cm per pass at perigee — negligible at the orbit-level (km-to-m)
controllability scale, but right at the deferred ~5 cm terminal-centering budget,
which is why it is carried in the truth model.

There is no separate "special relativity" force: the velocity-dependent terms are
part of the same post-Newtonian (Schwarzschild) expansion that Orekit's
``Relativity`` model applies, wired in :mod:`puffsat_sim.forces.build`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from puffsat_sim.constants import SPEED_OF_LIGHT_M_S, WGS84_MU


@dataclass(frozen=True)
class Relativity:
    """Schwarzschild post-Newtonian correction to Earth's point-mass gravity.

    Parameter-free: its only inputs are the central body's μ and the speed of
    light, both physical constants, so the spec carries no per-run parameters.
    """


def schwarzschild_apsidal_advance_per_orbit(semi_major_axis_m: float, eccentricity: float) -> float:
    """Prograde apsidal (perigee) advance per orbit due to relativity [rad].

    Δϖ = 6π·GM / (c²·a·(1 − e²)) — the same expression as Mercury's perihelion
    advance, evaluated for an Earth-orbiting body.
    """
    return (
        6.0
        * math.pi
        * WGS84_MU
        / (SPEED_OF_LIGHT_M_S**2 * semi_major_axis_m * (1.0 - eccentricity**2))
    )


def schwarzschild_perigee_advance_per_orbit_m(
    semi_major_axis_m: float, eccentricity: float
) -> float:
    """Along-track displacement at perigee from one orbit's apsidal advance [m].

    The apse line rotates by Δϖ; a point at perigee radius r_p = a(1 − e) is
    carried ≈ r_p·Δϖ across-track.  For the reference orbit this is ~4 cm/orbit.
    """
    r_p = semi_major_axis_m * (1.0 - eccentricity)
    return r_p * schwarzschild_apsidal_advance_per_orbit(semi_major_axis_m, eccentricity)
