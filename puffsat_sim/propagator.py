"""Orekit propagator factory — JVM boundary module.

Importing this module starts the JVM and loads Orekit data.  All
org.orekit.* imports elsewhere in the package must come after this module
has been imported.

build_propagator() is the sole public entry point for Rung 1.  Force models
(Rung 2+) are added here when PhysicsConfig gains non-Keplerian settings.
"""
from __future__ import annotations

from typing import Any

import orekit_jpype

_VM = orekit_jpype.initVM(
    vmargs="--enable-native-access=ALL-UNNAMED"
)  # must precede any org.orekit import

from orekit_jpype.pyhelpers import setup_orekit_curdir

setup_orekit_curdir()  # loads orekit-data.zip from the current working directory

from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.propagation.analytical import KeplerianPropagator
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants

from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.orbital_math import keplerian_elements


def build_propagator(orbital_config: OrbitalConfig, physics_config: PhysicsConfig) -> Any:
    """Return a configured Orekit propagator for the given run.

    Rung 1 (Keplerian): returns KeplerianPropagator — analytical, no force models.
    Rung 2+ (NumericalPropagator with force models): not yet implemented.

    The propagator epoch is set to orbital_config.epoch and the satellite is
    placed at mean anomaly = orbital_config.mean_anomaly_at_epoch_rad
    (default π = apogee, the nominal deployment point).

    orbital_config.epoch must be a whole-second UTC datetime (microseconds ignored).
    """
    if not physics_config.is_keplerian:
        raise NotImplementedError(
            "NumericalPropagator with force models is not yet implemented (Rung 2+)."
        )

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

    return KeplerianPropagator(orbit)
