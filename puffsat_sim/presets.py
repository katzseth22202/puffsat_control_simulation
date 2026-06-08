"""Force-model presets — named PhysicsConfig bundles.

Named by which perturbations are active, not by the design-doc build-ladder step
(which collides with the control-stage Rung A–D ladder).  The cumulative bundles
mirror the "validate one force at a time" story: each adds one force to the last.
"""

from __future__ import annotations

from puffsat_sim.config import PhysicsConfig
from puffsat_sim.forces import (
    AtmosphericDrag,
    Geopotential,
    Relativity,
    SolarRadiation,
    ThirdBody,
)

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


_DEFAULT_GEOPOTENTIAL_DEGREE = 8  # truth-model field degree/order (non-J2 gravity)
_DEFAULT_GEOPOTENTIAL_ORDER = 8


def full_force(
    cr_area_over_mass: float = _DEFAULT_CR_AREA_OVER_MASS,
    cd_area_over_mass: float = _DEFAULT_CD_AREA_OVER_MASS,
    geopotential_degree: int = _DEFAULT_GEOPOTENTIAL_DEGREE,
    geopotential_order: int = _DEFAULT_GEOPOTENTIAL_ORDER,
) -> PhysicsConfig:
    """Full truth model: 8×8 geopotential + third-body + SRP + drag + relativity.

    The geopotential defaults to degree/order 8 (not just J2): harmonics fall off
    as (Rₑ/r)^ℓ, so beyond J2 they are dead at the 150 000 km apogee and only bite
    in the last few thousand km of descent, where the low zonals reach ~metre scale
    at the 200 km interception.  The field degree is a tunable parameter — raise it
    for the regime-switched terminal phase / deferred 5 cm work; the coast does not
    need it.

    Relativity (Schwarzschild) is conservative and negligible at the orbit-level
    (km-to-m) scale, but is carried here because it reaches ~cm per pass on this
    high-e orbit — the deferred 5 cm terminal-centering budget.  Truth-only: it is
    below the estimator's noise floor, so it is not an onboard/filter force.
    """
    return PhysicsConfig(
        (
            Geopotential(degree=geopotential_degree, order=geopotential_order),
            ThirdBody(),
            SolarRadiation(cr_area_over_mass=cr_area_over_mass),
            AtmosphericDrag(cd_area_over_mass=cd_area_over_mass),
            Relativity(),
        )
    )
