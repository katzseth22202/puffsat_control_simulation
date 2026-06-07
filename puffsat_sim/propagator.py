"""Orekit propagator factory — JVM boundary module.

Importing this module starts the JVM and loads Orekit data.  All
org.orekit.* imports elsewhere in the package must come after this module
has been imported.

build_propagator() is the sole public entry point.  Adding a new force model
(Rung 2b, 2c, 2d) means adding one branch in _build_numerical_propagator();
the rest of the stack is unchanged.
"""
from __future__ import annotations

from typing import Any, Final

import orekit_jpype

_VM = orekit_jpype.initVM(
    vmargs="--enable-native-access=ALL-UNNAMED"
)  # must precede any org.orekit import

from orekit_jpype.pyhelpers import setup_orekit_curdir

setup_orekit_curdir()  # loads orekit-data.zip from the current working directory

from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
from org.orekit.forces.drag import DragForce, IsotropicDrag
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel, ThirdBodyAttraction
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient, SolarRadiationPressure
from org.orekit.models.earth.atmosphere import NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit, OrbitType, PositionAngleType
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.analytical import KeplerianPropagator
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions

from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.orbital_math import keplerian_elements

# ---------------------------------------------------------------------------
# Numerical integrator settings — coast phase (Rung 2).
# Tight relative tolerance avoids integration error masquerading as physical
# dispersion in the sensitivity analysis (design doc §11.2).
# ---------------------------------------------------------------------------
_MIN_STEP_S: Final[float] = 1.0      # 1 s — handles fast perigee pass
_MAX_STEP_S: Final[float] = 600.0    # 10 min — safe for the slow apogee region
_ABS_TOL_M: Final[float] = 1e-3     # 1 mm absolute position / velocity tolerance
_REL_TOL: Final[float] = 1e-10      # relative tolerance


def build_propagator(orbital_config: OrbitalConfig, physics_config: PhysicsConfig) -> Any:
    """Return a configured Orekit propagator for the given run.

    Keplerian (PhysicsConfig.is_keplerian=True): analytical KeplerianPropagator.
    Otherwise: NumericalPropagator with force models selected by PhysicsConfig.

    The propagator is initialised at orbital_config.epoch with mean anomaly
    orbital_config.mean_anomaly_at_epoch_rad (default π = apogee).

    orbital_config.epoch must be a whole-second UTC datetime.
    """
    utc = TimeScalesFactory.getUTC()
    frame = FramesFactory.getEME2000()
    mu: float = Constants.WGS84_EARTH_MU

    a, e = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)

    ep = orbital_config.epoch
    if ep.tzinfo is None:
        raise ValueError("OrbitalConfig.epoch must be timezone-aware UTC.")
    epoch = AbsoluteDate(ep.year, ep.month, ep.day, ep.hour, ep.minute, float(ep.second), utc)

    orbit = KeplerianOrbit(
        a,
        e,
        orbital_config.inclination_rad,
        orbital_config.raan_rad,
        orbital_config.arg_of_perigee_rad,
        orbital_config.mean_anomaly_at_epoch_rad,
        PositionAngleType.MEAN,
        frame,
        epoch,
        mu,
    )

    if physics_config.is_keplerian:
        return KeplerianPropagator(orbit)

    return _build_numerical_propagator(orbit, physics_config)


def _build_numerical_propagator(orbit: Any, physics_config: PhysicsConfig) -> Any:
    """Build a NumericalPropagator and attach the requested force models."""
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

    integrator = DormandPrince853Integrator(_MIN_STEP_S, _MAX_STEP_S, _ABS_TOL_M, _REL_TOL)
    propagator = NumericalPropagator(integrator)
    propagator.setOrbitType(OrbitType.KEPLERIAN)
    propagator.setPositionAngleType(PositionAngleType.MEAN)
    # mass=1.0 kg so that cr_area_over_mass and cd_area_over_mass are used directly
    # as effective cross-sections — the physical Cr/Cd are already folded in.
    propagator.setInitialState(SpacecraftState(orbit, 1.0))

    if physics_config.geopotential_degree > 0:
        provider = GravityFieldFactory.getNormalizedProvider(
            physics_config.geopotential_degree,
            physics_config.geopotential_degree,
        )
        propagator.addForceModel(HolmesFeatherstoneAttractionModel(itrf, provider))

    if physics_config.third_body:
        propagator.addForceModel(ThirdBodyAttraction(CelestialBodyFactory.getSun()))
        propagator.addForceModel(ThirdBodyAttraction(CelestialBodyFactory.getMoon()))

    if physics_config.srp_cr_area_over_mass is not None:
        earth_ellipsoid = OneAxisEllipsoid(
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.WGS84_EARTH_FLATTENING,
            itrf,
        )
        sun = CelestialBodyFactory.getSun()
        # absorptionCoeff=1.0 → effective Cr_orekit = 1; the physical Cr is already
        # encoded in srp_cr_area_over_mass (= Cr · A/m).  With SpacecraftState
        # mass=1 kg, Orekit computes a = P · 1 · srp_cr_area_over_mass / 1 = correct.
        srp_model = IsotropicRadiationSingleCoefficient(
            physics_config.srp_cr_area_over_mass, 1.0
        )
        propagator.addForceModel(SolarRadiationPressure(sun, earth_ellipsoid, srp_model))

    if physics_config.drag_cd_area_over_mass is not None:
        # NRLMSISE-00 driven by real/predicted CSSI space weather (covers 1957–2096).
        # PhysicsConfig.f10p7 / .ap are reserved for the Monte Carlo per-run
        # multiplicative bias applied on top of this nominal model (Rung D).
        earth_drag = OneAxisEllipsoid(
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.WGS84_EARTH_FLATTENING,
            itrf,
        )
        sun = CelestialBodyFactory.getSun()
        sw = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")
        atmosphere = NRLMSISE00(sw, sun, earth_drag)
        # Same mass=1.0 convention as SRP: cd_area_over_mass = Cd·(A/m) with Cd=1.
        drag_model = IsotropicDrag(physics_config.drag_cd_area_over_mass, 1.0)
        propagator.addForceModel(DragForce(atmosphere, drag_model))

    return propagator
