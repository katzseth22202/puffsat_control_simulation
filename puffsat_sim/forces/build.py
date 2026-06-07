"""Force-model construction — the JVM side of the perturbation seam.

``to_force_models`` turns one pure :class:`Perturbation` spec into the Orekit
``ForceModel``(s) that implement it.  The shared frames and bodies the models act
in are bundled in :class:`Environment`, built once per propagator so the WGS84
ellipsoid and Sun are not re-derived for every force.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, assert_never

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
from org.orekit.forces.drag import DragForce, IsotropicDrag
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel, ThirdBodyAttraction
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient, SolarRadiationPressure
from org.orekit.frames import FramesFactory
from org.orekit.models.earth.atmosphere import NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.utils import Constants, IERSConventions

from puffsat_sim.forces import (
    AtmosphericDrag,
    Geopotential,
    Perturbation,
    SolarRadiation,
    ThirdBody,
)


@dataclass(frozen=True)
class Environment:
    """Frames and bodies the force models act in, built once per propagator."""

    itrf: Any
    earth_ellipsoid: Any
    sun: Any
    moon: Any

    @classmethod
    def build(cls) -> Environment:
        itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        earth_ellipsoid = OneAxisEllipsoid(
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.WGS84_EARTH_FLATTENING,
            itrf,
        )
        return cls(
            itrf=itrf,
            earth_ellipsoid=earth_ellipsoid,
            sun=CelestialBodyFactory.getSun(),
            moon=CelestialBodyFactory.getMoon(),
        )


def to_force_models(perturbation: Perturbation, env: Environment) -> list[Any]:
    """Return the Orekit ForceModel(s) implementing one perturbation.

    The mass=1.0 kg convention in the propagator means the lumped Cr·(A/m) and
    Cd·(A/m) are passed as effective cross-sections with Orekit coefficient 1.
    """
    match perturbation:
        case Geopotential(degree=degree, order=order):
            provider = GravityFieldFactory.getNormalizedProvider(degree, order)
            return [HolmesFeatherstoneAttractionModel(env.itrf, provider)]
        case ThirdBody(sun=sun, moon=moon):
            models: list[Any] = []
            if sun:
                models.append(ThirdBodyAttraction(env.sun))
            if moon:
                models.append(ThirdBodyAttraction(env.moon))
            return models
        case SolarRadiation(cr_area_over_mass=cr_area_over_mass):
            srp_model = IsotropicRadiationSingleCoefficient(cr_area_over_mass, 1.0)
            return [SolarRadiationPressure(env.sun, env.earth_ellipsoid, srp_model)]
        case AtmosphericDrag(cd_area_over_mass=cd_area_over_mass):
            # NRLMSISE-00 driven by real/predicted CSSI space weather (1957–2096).
            # f10p7 / ap are reserved for the Monte Carlo per-run drag bias (Rung D).
            space_weather = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")
            atmosphere = NRLMSISE00(space_weather, env.sun, env.earth_ellipsoid)
            drag_model = IsotropicDrag(cd_area_over_mass, 1.0)
            return [DragForce(atmosphere, drag_model)]
        case _:
            assert_never(perturbation)
