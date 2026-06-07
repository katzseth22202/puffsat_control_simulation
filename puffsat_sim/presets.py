"""Force-model presets — named PhysicsConfig bundles.

Named by which perturbations are active, not by the design-doc build-ladder step
(which collides with the control-stage Rung A–D ladder).  The cumulative bundles
mirror the "validate one force at a time" story: each adds one force to the last.
"""

from __future__ import annotations

from puffsat_sim.config import PhysicsConfig
from puffsat_sim.forces import AtmosphericDrag, Geopotential, SolarRadiation, ThirdBody

_DEFAULT_CR_AREA_OVER_MASS = 0.02  # Cr·(A/m) [m²/kg]
_DEFAULT_CD_AREA_OVER_MASS = 0.04  # Cd·(A/m) [m²/kg]


def two_body() -> PhysicsConfig:
    """Point-mass / Keplerian — no perturbations."""
    return PhysicsConfig()


def j2() -> PhysicsConfig:
    """Zonal J2 geopotential only — degree 2, order 0."""
    return PhysicsConfig((Geopotential(degree=2),))


def j2_third_body() -> PhysicsConfig:
    """J2 + third-body Sun and Moon."""
    return PhysicsConfig((Geopotential(degree=2), ThirdBody()))


def j2_third_body_srp(
    cr_area_over_mass: float = _DEFAULT_CR_AREA_OVER_MASS,
) -> PhysicsConfig:
    """J2 + third-body + SRP cannonball."""
    return PhysicsConfig(
        (
            Geopotential(degree=2),
            ThirdBody(),
            SolarRadiation(cr_area_over_mass=cr_area_over_mass),
        )
    )


def full_force(
    cr_area_over_mass: float = _DEFAULT_CR_AREA_OVER_MASS,
    cd_area_over_mass: float = _DEFAULT_CD_AREA_OVER_MASS,
) -> PhysicsConfig:
    """Full force model: J2 + third-body + SRP + NRLMSISE-00 atmospheric drag."""
    return PhysicsConfig(
        (
            Geopotential(degree=2),
            ThirdBody(),
            SolarRadiation(cr_area_over_mass=cr_area_over_mass),
            AtmosphericDrag(cd_area_over_mass=cd_area_over_mass),
        )
    )
