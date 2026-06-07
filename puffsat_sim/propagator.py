"""Orekit propagator factory.

build_propagator() is the sole public entry point.  Which force models attach is
determined by the PhysicsConfig's perturbations; each is turned into Orekit force
models on the JVM side by :mod:`puffsat_sim.forces.build`.
"""
from __future__ import annotations

from typing import Any, Final

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit, OrbitType, PositionAngleType
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.analytical import KeplerianPropagator
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants

from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.forces.build import Environment, to_force_models
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

    No perturbations (PhysicsConfig.is_keplerian): analytical KeplerianPropagator.
    Otherwise: NumericalPropagator with force models built from the perturbations.

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
    """Build a NumericalPropagator and attach the perturbations' force models."""
    integrator = DormandPrince853Integrator(_MIN_STEP_S, _MAX_STEP_S, _ABS_TOL_M, _REL_TOL)
    propagator = NumericalPropagator(integrator)
    propagator.setOrbitType(OrbitType.KEPLERIAN)
    propagator.setPositionAngleType(PositionAngleType.MEAN)
    # mass=1.0 kg so that the lumped Cr·(A/m) and Cd·(A/m) are used directly as
    # effective cross-sections — the physical Cr/Cd are already folded into them.
    propagator.setInitialState(SpacecraftState(orbit, 1.0))

    env = Environment.build()
    for perturbation in physics_config.perturbations:
        for force_model in to_force_models(perturbation, env):
            propagator.addForceModel(force_model)

    return propagator
