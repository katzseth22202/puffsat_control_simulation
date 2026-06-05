"""Hello Orekit: verify the Python/JVM bridge with a PuffSat reference orbit.

Run with:
    make run
or:
    python -m puffsat_sim.hello_orekit

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

from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.propagation.analytical import KeplerianPropagator
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants

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

    print("PuffSat Control Simulation — Orekit / JVM bridge hello world")
    print("  Python/JVM : OK")
    print()
    print("  Reference orbit (near-term architecture):")
    print(f"    Orbit periapsis  : {_PERIGEE_ALT_M / 1e3:.0f} km  (burns up here; interception at 200 km during descent)")
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


def main() -> None:
    propagate_one_period()


if __name__ == "__main__":
    main()
