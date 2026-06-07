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
    _EARTH_RADIUS_M,
    drag_deceleration,
    j2_apsidal_precession_rate,
    j2_nodal_regression_rate,
    keplerian_elements,
    keplerian_period,
    lunar_tidal_ratio,
    orbital_config_from_cities,
    perigee_speed,
    solar_tidal_ratio,
    srp_acceleration,
    wrap_to_pi,
)

# Importing propagator starts the JVM and loads Orekit data.
# All org.orekit.* imports must follow this line.
from puffsat_sim.propagator import build_propagator  # noqa: E402

from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit
from org.orekit.propagation.events import AltitudeDetector, EclipseDetector, EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent, StopOnDecreasing
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions

# ---------------------------------------------------------------------------
# Mission constants — fixed across all runs and Monte Carlo draws
# ---------------------------------------------------------------------------

_PERIGEE_ALT_M: Final[float] = 50_000.0        # orbit periapsis [m]; debris disposal below Kármán
_APOGEE_ALT_M: Final[float] = 150_000_000.0    # deployment apogee [m]
_INTERCEPTION_ALT_M: Final[float] = 200_000.0  # control target: 200 km descent crossing
_SUN_RADIUS_M: Final[float] = 6.957e8          # solar radius [m] — used by EclipseDetector

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

    # Wrap to (-π, π]: the osculating angles straddle the 0/2π branch cut, so a
    # raw subtraction reports a small retrograde drift as ~+2π instead of ~−ε.
    d_raan_num = math.degrees(wrap_to_pi(raan_f - raan_0))
    d_omega_num = math.degrees(wrap_to_pi(omega_f - omega_0))

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


def validate_third_body_signatures(orbital_config: OrbitalConfig) -> None:
    """Verify third-body Sun + Moon perturbations (Rung 2b).

    Confirms that ThirdBodyAttraction is active by comparing one-period
    propagation with J2-only vs J2+Sun+Moon, then reports the analytic tidal
    acceleration ratios at apogee (Hill approximation, design doc §2 benchmark).
    """
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    epoch = _to_absolute_date(orbital_config.epoch)

    def _final_pos(physics_config: PhysicsConfig) -> tuple[float, float, float]:
        prop = build_propagator(orbital_config, physics_config)
        state = prop.propagate(epoch.shiftedBy(period))
        pv = state.getPVCoordinates().getPosition()
        return float(pv.getX()), float(pv.getY()), float(pv.getZ())

    pos_j2 = _final_pos(PhysicsConfig.rung_2a())
    pos_j2_tb = _final_pos(PhysicsConfig.rung_2b())

    dr = math.sqrt(sum((a - b) ** 2 for a, b in zip(pos_j2, pos_j2_tb, strict=True)))

    moon_pct = lunar_tidal_ratio(orbital_config.apogee_alt_m) * 100.0
    sun_pct = solar_tidal_ratio(orbital_config.apogee_alt_m) * 100.0

    print("  Third-body signature validation (Rung 2b):")
    print(f"    Tidal ratio at apogee  — Moon: {moon_pct:.3f}%  Sun: {sun_pct:.3f}%")
    print(f"    One-period position drift (J2 → J2+Sun+Moon): {dr / 1e3:.3f} km")


def validate_srp_signatures(orbital_config: OrbitalConfig) -> None:
    """Verify SRP cannonball force model and shadow/eclipse detection (Rung 2c).

    Confirms SRP is active by comparing one-period propagation with J2+third_body
    vs J2+third_body+SRP, then counts eclipse entry/exit events using EclipseDetector
    to verify the cylindrical shadow model is wired.  Prints the analytic SRP
    acceleration at 1 AU as the order-of-magnitude benchmark.
    """
    physics_2c = PhysicsConfig.rung_2c()
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    epoch = _to_absolute_date(orbital_config.epoch)

    def _final_pos(physics_config: PhysicsConfig) -> tuple[float, float, float]:
        prop = build_propagator(orbital_config, physics_config)
        state = prop.propagate(epoch.shiftedBy(period))
        pv = state.getPVCoordinates().getPosition()
        return float(pv.getX()), float(pv.getY()), float(pv.getZ())

    pos_2b = _final_pos(PhysicsConfig.rung_2b())
    pos_2c = _final_pos(physics_2c)
    dr = math.sqrt(
        sum((a - b) ** 2 for a, b in zip(pos_2b, pos_2c, strict=True))
    )

    # Eclipse detection: count shadow entry/exit crossings during one period.
    sun = CelestialBodyFactory.getSun()
    earth = OneAxisEllipsoid(
        Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
        Constants.WGS84_EARTH_FLATTENING,
        FramesFactory.getITRF(IERSConventions.IERS_2010, True),
    )
    eclipse_prop = build_propagator(orbital_config, physics_2c)
    logger: EventsLogger = EventsLogger()  # type: ignore[no-untyped-call]
    # EclipseDetector defaults to StopOnIncreasing, which would halt the arc at
    # the first umbra exit; ContinueOnEvent lets it log every crossing across the
    # whole period so the count is meaningful.
    eclipse_prop.addEventDetector(
        logger.monitorDetector(
            EclipseDetector(sun, _SUN_RADIUS_M, earth).withHandler(
                ContinueOnEvent()  # type: ignore[no-untyped-call]
            )
        )
    )
    eclipse_prop.propagate(epoch.shiftedBy(period))
    n_eclipse = len(logger.getLoggedEvents())

    cr_am = physics_2c.srp_cr_area_over_mass
    assert cr_am is not None
    accel = srp_acceleration(cr_am)

    print("  SRP signature validation (Rung 2c):")
    print(f"    Analytic SRP accel (1 AU, Cr·A/m={cr_am}): {accel:.3e} m/s²")
    print(f"    One-period divergence (J2+3body → +SRP): {dr / 1e3:.3f} km")
    print(f"    Eclipse crossings detected in one period: {n_eclipse}")


def validate_drag_signatures(orbital_config: OrbitalConfig) -> None:
    """Verify NRLMSISE-00 atmospheric drag (Rung 2d).

    Propagates from apogee to the 200 km descent crossing with SRP-only (rung_2c)
    and full-force (rung_2d).  The difference in orbital specific energy at that
    altitude confirms drag is removing energy during the terminal descent; the
    analytic deceleration at key altitudes gives the expected order of magnitude.

    Design doc §4: drag "bites below ~300-400 km" — validated by comparing
    energy at 200 km (interception altitude) with and without drag.
    """
    physics_2d = PhysicsConfig.rung_2d()
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    mu: float = float(Constants.WGS84_EARTH_MU)

    def _energy_and_speed_at_200km(physics_config: PhysicsConfig) -> tuple[float, float]:
        earth = OneAxisEllipsoid(
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.WGS84_EARTH_FLATTENING,
            FramesFactory.getITRF(IERSConventions.IERS_2010, True),
        )
        prop = build_propagator(orbital_config, physics_config)
        prop.addEventDetector(
            AltitudeDetector(_INTERCEPTION_ALT_M, earth).withHandler(
                StopOnDecreasing()  # type: ignore[no-untyped-call]
            )
        )
        epoch = _to_absolute_date(orbital_config.epoch)
        state = prop.propagate(epoch.shiftedBy(period))
        pv = state.getPVCoordinates()
        pos = pv.getPosition()
        vel = pv.getVelocity()
        r = math.sqrt(pos.getX() ** 2 + pos.getY() ** 2 + pos.getZ() ** 2)
        v_sq = vel.getX() ** 2 + vel.getY() ** 2 + vel.getZ() ** 2
        energy = 0.5 * v_sq - mu / r
        return energy, math.sqrt(v_sq)

    # Both arcs stop at the same 200 km altitude (the event fixes position and time),
    # so the meaningful drag signatures are the loss of specific orbital energy and the
    # resulting speed reduction at the crossing — not the position, which barely moves.
    energy_2c, speed_2c = _energy_and_speed_at_200km(PhysicsConfig.rung_2c())
    energy_2d, speed_2d = _energy_and_speed_at_200km(physics_2d)
    d_energy = energy_2d - energy_2c
    d_speed = speed_2d - speed_2c

    cd_am = physics_2d.drag_cd_area_over_mass
    assert cd_am is not None
    v_200 = math.sqrt(mu * (2.0 / (_EARTH_RADIUS_M + 200_000.0) - 1.0 / a))
    v_300 = math.sqrt(mu * (2.0 / (_EARTH_RADIUS_M + 300_000.0) - 1.0 / a))
    a_drag_200 = drag_deceleration(cd_am, v_200, 200_000.0)
    a_drag_300 = drag_deceleration(cd_am, v_300, 300_000.0)

    print("  Drag signature validation (Rung 2d, NRLMSISE-00):")
    print(f"    Analytic drag decel @ 200 km: {a_drag_200:.3e} m/s²")
    print(f"    Analytic drag decel @ 300 km: {a_drag_300:.3e} m/s²")
    print(f"    Orbital energy at 200 km — no drag: {energy_2c:.1f} J/kg")
    print(f"    Orbital energy at 200 km — with drag: {energy_2d:.1f} J/kg")
    print(f"    ΔE from drag: {d_energy:.1f} J/kg  (negative = drag removed energy)")
    print(f"    Speed reduction at 200 km from drag: {d_speed * 100.0:+.3f} cm/s")


def main() -> None:
    propagate_one_period(_NOMINAL_CONFIG, PhysicsConfig.rung_keplerian())
    print()
    propagate_to_interception(_NOMINAL_CONFIG, PhysicsConfig.rung_keplerian())
    print()
    validate_j2_signatures(_NOMINAL_CONFIG)
    print()
    validate_third_body_signatures(_NOMINAL_CONFIG)
    print()
    validate_srp_signatures(_NOMINAL_CONFIG)
    print()
    validate_drag_signatures(_NOMINAL_CONFIG)


if __name__ == "__main__":
    main()
