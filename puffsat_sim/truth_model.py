"""Rung A truth model: Keplerian propagation of the PuffSat reference orbit.

Verifies the Orekit / JVM bridge and the reference orbit parameters before
perturbation force models are added (Rung A of the design doc build ladder).

Run with:
    make run
or:
    python -m puffsat_sim.truth_model

The orbit used matches the near-term architecture from the paper:
  - periapsis 50 km (orbit periapsis; PuffSat burns up here after impact)
  - interception at 200 km during descent, before periapsis
  - apogee  ~150 000 km altitude (recommended deployment apogee from design doc)
  - eccentricity ~0.921, period ~2.68 days
  - perigee speed ~10.91 km/s
"""
from __future__ import annotations

import math
from typing import Final

# orekit_jpype is the conda-forge JPype-based Orekit wrapper.
# orekit_jpype and org.orekit.* have no type stubs; mypy suppresses all errors
# for these modules via the overrides in pyproject.toml.  Our own code is fully typed.
import orekit_jpype

_VM = orekit_jpype.initVM(
    vmargs="--enable-native-access=ALL-UNNAMED"
)  # start the JVM — must precede any org.orekit import

from orekit_jpype.pyhelpers import setup_orekit_curdir

# setup_orekit_curdir looks for orekit-data.zip in the current working directory.
# Run `make data` from the project root once to download it.
setup_orekit_curdir()

from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.propagation.analytical import KeplerianPropagator
from org.orekit.propagation.events import AltitudeDetector
from org.orekit.propagation.events.handlers import StopOnDecreasing
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions

from puffsat_sim.orbital_math import keplerian_elements, perigee_speed

# ---------------------------------------------------------------------------
# Orbital parameters — near-term PuffSat architecture (paper §2 / design doc §3)
# ---------------------------------------------------------------------------

# Orbit periapsis: 50 km — intentionally below the Kármán line so the PuffSat
# burns up after impact (debris disposal).  Interception occurs at 200 km during
# descent, before this periapsis is reached (design doc §3).
_PERIGEE_ALT_M: Final[float] = 50_000.0

# Deployment apogee: design doc recommends ~150 000 km altitude as a balance
# between perigee-speed (~10.8 km/s, nearly independent of apogee) and
# controllability (solar tidal force <0.1% of Earth gravity at apogee).
_APOGEE_ALT_M: Final[float] = 150_000_000.0

# Interception altitude: 200 km during descent, before the 50 km periapsis.
# This is the control target for the simulation (design doc §1, §6.3, §10.3).
_INTERCEPTION_ALT_M: Final[float] = 200_000.0


def propagate_one_period() -> None:
    """Propagate the PuffSat reference orbit for one Keplerian period."""
    utc = TimeScalesFactory.getUTC()
    frame = FramesFactory.getEME2000()
    mu: float = Constants.WGS84_EARTH_MU

    a, e = keplerian_elements(_PERIGEE_ALT_M, _APOGEE_ALT_M)

    # Start at apogee (mean anomaly = π) so the clock ticks from deployment
    # toward interception — more representative of the mission timeline.
    epoch = AbsoluteDate(2026, 1, 1, 0, 0, 0.0, utc)
    orbit = KeplerianOrbit(
        a,
        e,
        math.radians(28.5),   # inclination: mid-latitude launch site
        math.radians(0.0),    # RAAN
        math.radians(0.0),    # argument of perigee
        math.radians(180.0),  # mean anomaly = π → start at apogee
        PositionAngleType.MEAN,
        frame,
        epoch,
        mu,
    )

    period: float = orbit.getKeplerianPeriod()
    v_perigee: float = perigee_speed(a, _PERIGEE_ALT_M)

    propagator = KeplerianPropagator(orbit)
    final_state = propagator.propagate(epoch.shiftedBy(period))
    pos = final_state.getPVCoordinates().getPosition()
    vel = final_state.getPVCoordinates().getVelocity()

    # Residual after one Keplerian period (should be ~floating-point zero)
    r_0 = orbit.getPVCoordinates().getPosition()
    v_0 = orbit.getPVCoordinates().getVelocity()
    dr = math.sqrt(
        (pos.getX() - r_0.getX()) ** 2
        + (pos.getY() - r_0.getY()) ** 2
        + (pos.getZ() - r_0.getZ()) ** 2
    )
    dv = math.sqrt(
        (vel.getX() - v_0.getX()) ** 2
        + (vel.getY() - v_0.getY()) ** 2
        + (vel.getZ() - v_0.getZ()) ** 2
    )

    print("PuffSat Control Simulation — Rung A: Keplerian reference orbit")
    print("  Orekit / JVM : OK")
    print()
    print("  Reference orbit (near-term architecture):")
    print(
        f"    Orbit periapsis  : {_PERIGEE_ALT_M / 1e3:.0f} km"
        "  (burns up here; interception at 200 km during descent)"
    )
    print(f"    Apogee altitude  : {_APOGEE_ALT_M / 1e6:.0f} × 10³ km  (deployment)")
    print(f"    Semi-major axis  : {a / 1e3:.1f} km")
    print(f"    Eccentricity     : {e:.6f}")
    print("    Inclination      : 28.5°")
    print(f"    Orbital period   : {period:.1f} s  ({period / 86400:.2f} days)")
    print(f"    Perigee speed    : {v_perigee / 1e3:.3f} km/s")
    print()
    print("  One-period propagation residual (Keplerian → should be ~0):")
    print(f"    |Δr| = {dr:.3e} m")
    print(f"    |Δv| = {dv:.3e} m/s")


def propagate_to_interception() -> None:
    """Propagate from apogee, stopping at the 200 km descent crossing (interception).

    Uses AltitudeDetector with StopOnDecreasing: the g-function is
    (altitude − 200 km), which decreases through zero as the PuffSat descends.
    Starting at apogee, the first zero-crossing is the descending one, so no
    additional filtering is needed.
    """
    utc = TimeScalesFactory.getUTC()
    frame = FramesFactory.getEME2000()
    mu: float = Constants.WGS84_EARTH_MU

    a, e = keplerian_elements(_PERIGEE_ALT_M, _APOGEE_ALT_M)
    epoch = AbsoluteDate(2026, 1, 1, 0, 0, 0.0, utc)
    orbit = KeplerianOrbit(
        a,
        e,
        math.radians(28.5),
        math.radians(0.0),
        math.radians(0.0),
        math.radians(180.0),  # mean anomaly = π → start at apogee
        PositionAngleType.MEAN,
        frame,
        epoch,
        mu,
    )

    # WGS84 ellipsoid in the Earth-fixed frame — used by AltitudeDetector to
    # compute geodetic altitude above the surface.
    earth = OneAxisEllipsoid(
        Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
        Constants.WGS84_EARTH_FLATTENING,
        FramesFactory.getITRF(IERSConventions.IERS_2010, True),
    )
    detector = AltitudeDetector(_INTERCEPTION_ALT_M, earth).withHandler(
        StopOnDecreasing()  # type: ignore[no-untyped-call]
    )

    period: float = orbit.getKeplerianPeriod()
    propagator = KeplerianPropagator(orbit)
    propagator.addEventDetector(detector)

    # Upper bound is one full period; event fires ~halfway through (apo → peri descent).
    final_state = propagator.propagate(epoch.shiftedBy(period))

    elapsed: float = final_state.getDate().durationFrom(epoch)
    pos = final_state.getPVCoordinates().getPosition()
    vel = final_state.getPVCoordinates().getVelocity()
    v: float = math.sqrt(
        vel.getX() ** 2 + vel.getY() ** 2 + vel.getZ() ** 2
    )
    geodetic = earth.transform(pos, frame, final_state.getDate())
    alt_km: float = geodetic.getAltitude() / 1e3

    print("  Propagation to interception (200 km descent crossing):")
    print(f"    Coast time from apogee : {elapsed / 3600:.3f} h  ({elapsed / 86400:.3f} days)")
    print(f"    Altitude at stop       : {alt_km:.3f} km  (event target: 200 km)")
    print(f"    Speed at interception  : {v / 1e3:.3f} km/s")


def main() -> None:
    propagate_one_period()
    print()
    propagate_to_interception()


if __name__ == "__main__":
    main()
