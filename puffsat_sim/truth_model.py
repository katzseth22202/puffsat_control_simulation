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
from datetime import UTC, datetime
from typing import Final

from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.orbital_math import (
    j2_apsidal_precession_rate,
    j2_nodal_regression_rate,
    keplerian_elements,
    keplerian_period,
    orbital_config_from_cities,
    perigee_speed,
)

# Importing propagator starts the JVM and loads Orekit data.
# All org.orekit.* imports must follow this line.
from puffsat_sim.propagator import build_propagator  # noqa: E402

from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit
from org.orekit.propagation.events import AltitudeDetector
from org.orekit.propagation.events.handlers import StopOnDecreasing
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions

# ---------------------------------------------------------------------------
# Mission constants — fixed across all runs and Monte Carlo draws
# ---------------------------------------------------------------------------

_PERIGEE_ALT_M: Final[float] = 50_000.0        # orbit periapsis [m]; debris disposal below Kármán
_APOGEE_ALT_M: Final[float] = 150_000_000.0    # deployment apogee [m]
_INTERCEPTION_ALT_M: Final[float] = 200_000.0  # control target: 200 km descent crossing

# ---------------------------------------------------------------------------
# Nominal orbital plane — defined by a great circle through two surface points.
#
# The specific cities are arbitrary: we want a realistic mid-to-high inclination
# (~70°) for the perturbation study.  Ground tracks are not modelled and the
# simulation is not sensitive to which locations are used — only the resulting
# inclination and RAAN matter.  The epoch sets the RAAN via GMST.
# ---------------------------------------------------------------------------

_TOKYO: Final[tuple[float, float]] = (35.6762, 139.6503)    # (lat°N, lon°E)
_NEW_YORK: Final[tuple[float, float]] = (40.7128, -74.0060)  # (lat°N, lon°W)
_EPOCH: Final[datetime] = datetime(2026, 6, 2, 0, 0, 0, tzinfo=UTC)

_NOMINAL_CONFIG: Final[OrbitalConfig] = orbital_config_from_cities(
    *_TOKYO,
    *_NEW_YORK,
    epoch=_EPOCH,
    perigee_alt_m=_PERIGEE_ALT_M,
    apogee_alt_m=_APOGEE_ALT_M,
)


def _to_absolute_date(dt: datetime) -> AbsoluteDate:
    utc = TimeScalesFactory.getUTC()
    return AbsoluteDate(dt.year, dt.month, dt.day, dt.hour, dt.minute, float(dt.second), utc)


def propagate_one_period(orbital_config: OrbitalConfig, physics_config: PhysicsConfig) -> None:
    """Propagate the PuffSat reference orbit for one Keplerian period."""
    a, e = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    v_perigee = perigee_speed(a, orbital_config.perigee_alt_m)

    epoch = _to_absolute_date(orbital_config.epoch)
    propagator = build_propagator(orbital_config, physics_config)
    initial_pv = propagator.getInitialState().getPVCoordinates()
    r_0 = initial_pv.getPosition()
    v_0 = initial_pv.getVelocity()

    final_state = propagator.propagate(epoch.shiftedBy(period))
    pos = final_state.getPVCoordinates().getPosition()
    vel = final_state.getPVCoordinates().getVelocity()

    # Residual after one Keplerian period (should be ~floating-point zero)
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
        f"    Orbit periapsis  : {orbital_config.perigee_alt_m / 1e3:.0f} km"
        "  (burns up here; interception at 200 km during descent)"
    )
    print(f"    Apogee altitude  : {orbital_config.apogee_alt_m / 1e6:.0f} × 10³ km  (deployment)")
    print(f"    Semi-major axis  : {a / 1e3:.1f} km")
    print(f"    Eccentricity     : {e:.6f}")
    print(f"    Inclination      : {math.degrees(orbital_config.inclination_rad):.1f}°")
    print(f"    Orbital period   : {period:.1f} s  ({period / 86400:.2f} days)")
    print(f"    Perigee speed    : {v_perigee / 1e3:.3f} km/s")
    print()
    print("  One-period propagation residual (Keplerian → should be ~0):")
    print(f"    |Δr| = {dr:.3e} m")
    print(f"    |Δv| = {dv:.3e} m/s")


def propagate_to_interception(
    orbital_config: OrbitalConfig, physics_config: PhysicsConfig
) -> None:
    """Propagate from apogee, stopping at the 200 km descent crossing (interception).

    Uses AltitudeDetector with StopOnDecreasing: the g-function is
    (altitude − 200 km), which decreases through zero as the PuffSat descends.
    Starting at apogee, the first zero-crossing is the descending one, so no
    additional filtering is needed.
    """
    frame = FramesFactory.getEME2000()
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)

    epoch = _to_absolute_date(orbital_config.epoch)
    propagator = build_propagator(orbital_config, physics_config)

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
    propagator.addEventDetector(detector)

    # Upper bound is one full period; event fires ~halfway through (apo → peri descent).
    final_state = propagator.propagate(epoch.shiftedBy(period))

    elapsed: float = final_state.getDate().durationFrom(epoch)
    pos = final_state.getPVCoordinates().getPosition()
    vel = final_state.getPVCoordinates().getVelocity()
    v: float = math.sqrt(vel.getX() ** 2 + vel.getY() ** 2 + vel.getZ() ** 2)
    geodetic = earth.transform(pos, frame, final_state.getDate())
    alt_km: float = geodetic.getAltitude() / 1e3

    print("  Propagation to interception (200 km descent crossing):")
    print(f"    Coast time from apogee : {elapsed / 3600:.3f} h  ({elapsed / 86400:.3f} days)")
    print(f"    Altitude at stop       : {alt_km:.3f} km  (event target: 200 km)")
    print(f"    Speed at interception  : {v / 1e3:.3f} km/s")


def validate_j2_signatures(orbital_config: OrbitalConfig) -> None:
    """Propagate one period with J2 and compare ΔRAAN / Δω against analytic predictions.

    Checks that the NumericalPropagator + HolmesFeatherstoneAttractionModel (degree 2)
    produces nodal regression and apsidal precession rates consistent with the
    first-order J2 secular formulas in orbital_math.  Agreement within ~1% confirms
    the force model, integrator tolerances, and rate formulas are mutually consistent.
    """
    a, e = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    i = orbital_config.inclination_rad

    epoch = _to_absolute_date(orbital_config.epoch)
    propagator = build_propagator(orbital_config, PhysicsConfig.rung_2a())

    initial_state = propagator.getInitialState()
    initial_orbit = KeplerianOrbit(initial_state.getOrbit())
    raan_0: float = initial_orbit.getRightAscensionOfAscendingNode()
    omega_0: float = initial_orbit.getPerigeeArgument()

    final_state = propagator.propagate(epoch.shiftedBy(period))
    final_orbit = KeplerianOrbit(final_state.getOrbit())
    raan_f: float = final_orbit.getRightAscensionOfAscendingNode()
    omega_f: float = final_orbit.getPerigeeArgument()

    d_raan_num = math.degrees(raan_f - raan_0)
    d_omega_num = math.degrees(omega_f - omega_0)

    rate_raan = j2_nodal_regression_rate(a, e, i)
    rate_omega = j2_apsidal_precession_rate(a, e, i)
    d_raan_analytic = math.degrees(rate_raan * period)
    d_omega_analytic = math.degrees(rate_omega * period)

    pct_raan = abs(d_raan_num - d_raan_analytic) / abs(d_raan_analytic) * 100.0
    pct_omega = abs(d_omega_num - d_omega_analytic) / abs(d_omega_analytic) * 100.0

    print("  J2 signature validation (Rung 2a):")
    print(
        f"    ΔRAAN  numeric: {d_raan_num:+.4f}°  "
        f"analytic: {d_raan_analytic:+.4f}°  agree: {100.0 - pct_raan:.2f}%"
    )
    print(
        f"    Δω     numeric: {d_omega_num:+.4f}°  "
        f"analytic: {d_omega_analytic:+.4f}°  agree: {100.0 - pct_omega:.2f}%"
    )


def main() -> None:
    propagate_one_period(_NOMINAL_CONFIG, PhysicsConfig.rung_keplerian())
    print()
    propagate_to_interception(_NOMINAL_CONFIG, PhysicsConfig.rung_keplerian())
    print()
    validate_j2_signatures(_NOMINAL_CONFIG)


if __name__ == "__main__":
    main()
