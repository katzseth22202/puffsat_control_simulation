"""Integration tests: exercise build_propagator against analytic physics.

These require a live Orekit JVM and orekit-data.zip in the working directory;
the module skips itself cleanly when either is unavailable.  They assert the
quantities that puffsat_sim.truth_model only prints, so the force models are
actually verified rather than eyeballed.

Run with:
    make test-integration
or:
    pytest -m integration tests/integration
"""
from __future__ import annotations

import math
from datetime import UTC, datetime
from typing import Any

import pytest

from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.orbital_math import (
    j2_apsidal_precession_rate,
    j2_nodal_regression_rate,
    keplerian_elements,
    keplerian_period,
    orbital_config_from_cities,
    wrap_to_pi,
)

try:
    # Importing propagator boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.propagator import build_propagator

    from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
    from org.orekit.frames import FramesFactory
    from org.orekit.orbits import KeplerianOrbit
    from org.orekit.propagation.events import (
        AltitudeDetector,
        EclipseDetector,
        EventsLogger,
    )
    from org.orekit.propagation.events.handlers import ContinueOnEvent, StopOnDecreasing
    from org.orekit.time import AbsoluteDate, TimeScalesFactory
    from org.orekit.utils import Constants, IERSConventions
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

pytestmark = pytest.mark.integration

_EPOCH = datetime(2026, 6, 2, 0, 0, 0, tzinfo=UTC)
_PERIGEE_ALT_M = 50_000.0
_APOGEE_ALT_M = 150_000_000.0
_INTERCEPTION_ALT_M = 200_000.0
_SUN_RADIUS_M = 6.957e8


@pytest.fixture(scope="module")
def config() -> OrbitalConfig:
    return orbital_config_from_cities(
        35.6762,
        139.6503,
        40.7128,
        -74.0060,
        epoch=_EPOCH,
        perigee_alt_m=_PERIGEE_ALT_M,
        apogee_alt_m=_APOGEE_ALT_M,
    )


def _abs_date(dt: datetime) -> Any:
    utc = TimeScalesFactory.getUTC()
    return AbsoluteDate(dt.year, dt.month, dt.day, dt.hour, dt.minute, float(dt.second), utc)


def _period(config: OrbitalConfig) -> float:
    a, _ = keplerian_elements(config.perigee_alt_m, config.apogee_alt_m)
    return keplerian_period(a)


def _earth() -> Any:
    return OneAxisEllipsoid(
        Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
        Constants.WGS84_EARTH_FLATTENING,
        FramesFactory.getITRF(IERSConventions.IERS_2010, True),
    )


def _final_position(config: OrbitalConfig, physics: PhysicsConfig) -> tuple[float, float, float]:
    prop = build_propagator(config, physics)
    state = prop.propagate(_abs_date(config.epoch).shiftedBy(_period(config)))
    p = state.getPVCoordinates().getPosition()
    return float(p.getX()), float(p.getY()), float(p.getZ())


@pytest.fixture(scope="module")
def pos_j2(config: OrbitalConfig) -> tuple[float, float, float]:
    return _final_position(config, PhysicsConfig.rung_2a())


@pytest.fixture(scope="module")
def pos_j2_third_body(config: OrbitalConfig) -> tuple[float, float, float]:
    return _final_position(config, PhysicsConfig.rung_2b())


@pytest.fixture(scope="module")
def pos_srp(config: OrbitalConfig) -> tuple[float, float, float]:
    return _final_position(config, PhysicsConfig.rung_2c())


def test_keplerian_one_period_round_trip(config: OrbitalConfig) -> None:
    """A pure-Keplerian orbit must return to its start after exactly one period."""
    prop = build_propagator(config, PhysicsConfig.rung_keplerian())
    pv0 = prop.getInitialState().getPVCoordinates()
    state = prop.propagate(_abs_date(config.epoch).shiftedBy(_period(config)))
    pv = state.getPVCoordinates()
    assert float(pv0.getPosition().distance(pv.getPosition())) < 1.0  # actual ~1e-8 m
    assert float(pv0.getVelocity().distance(pv.getVelocity())) < 1e-3


def test_interception_stops_at_200km(config: OrbitalConfig) -> None:
    """AltitudeDetector must stop the descending arc at the 200 km crossing."""
    earth = _earth()
    prop = build_propagator(config, PhysicsConfig.rung_keplerian())
    prop.addEventDetector(
        AltitudeDetector(_INTERCEPTION_ALT_M, earth).withHandler(StopOnDecreasing())
    )
    epoch = _abs_date(config.epoch)
    period = _period(config)
    state = prop.propagate(epoch.shiftedBy(period))

    elapsed = float(state.getDate().durationFrom(epoch))
    geodetic = earth.transform(
        state.getPVCoordinates().getPosition(), FramesFactory.getEME2000(), state.getDate()
    )
    assert geodetic.getAltitude() == pytest.approx(_INTERCEPTION_ALT_M, abs=10.0)
    # Event fired mid-orbit (apogee -> descent), well before a full period elapsed.
    assert 0.0 < elapsed < period


def test_j2_secular_rates_match_analytic(config: OrbitalConfig) -> None:
    """Numerical J2 nodal/apsidal drift over one period must match the analytic rates.

    This also guards the wrap_to_pi fix: the raw osculating ω difference straddles
    the 0/2π branch cut, so the test only passes if the difference is wrapped.
    """
    a, e = keplerian_elements(config.perigee_alt_m, config.apogee_alt_m)
    period = _period(config)
    i = config.inclination_rad

    prop = build_propagator(config, PhysicsConfig.rung_2a())
    orbit_0 = KeplerianOrbit(prop.getInitialState().getOrbit())
    orbit_f = KeplerianOrbit(prop.propagate(_abs_date(config.epoch).shiftedBy(period)).getOrbit())

    d_raan = wrap_to_pi(
        orbit_f.getRightAscensionOfAscendingNode() - orbit_0.getRightAscensionOfAscendingNode()
    )
    d_omega = wrap_to_pi(orbit_f.getPerigeeArgument() - orbit_0.getPerigeeArgument())

    assert d_raan == pytest.approx(j2_nodal_regression_rate(a, e, i) * period, rel=0.05)
    assert d_omega == pytest.approx(j2_apsidal_precession_rate(a, e, i) * period, rel=0.10)


def test_third_body_perturbs_orbit(
    pos_j2: tuple[float, float, float],
    pos_j2_third_body: tuple[float, float, float],
) -> None:
    """Adding Sun + Moon must move the one-period endpoint by a measurable amount."""
    # Actual drift at 150 000 km apogee is ~100 km; require >1 km to confirm it is wired.
    assert math.dist(pos_j2, pos_j2_third_body) > 1_000.0


def test_srp_perturbs_orbit(
    pos_j2_third_body: tuple[float, float, float],
    pos_srp: tuple[float, float, float],
) -> None:
    """Enabling SRP must shift the one-period endpoint (cannonball model active)."""
    # Actual divergence is ~400 m over one period; require >10 m.
    assert math.dist(pos_j2_third_body, pos_srp) > 10.0


def test_drag_removes_orbital_energy(config: OrbitalConfig) -> None:
    """Specific orbital energy at the 200 km crossing must be lower with drag on."""
    mu = float(Constants.WGS84_EARTH_MU)
    earth = _earth()

    def energy_at_interception(physics: PhysicsConfig) -> float:
        prop = build_propagator(config, physics)
        prop.addEventDetector(
            AltitudeDetector(_INTERCEPTION_ALT_M, earth).withHandler(StopOnDecreasing())
        )
        state = prop.propagate(_abs_date(config.epoch).shiftedBy(_period(config)))
        pv = state.getPVCoordinates()
        r = float(pv.getPosition().getNorm())
        v = float(pv.getVelocity().getNorm())
        return 0.5 * v * v - mu / r

    energy_no_drag = energy_at_interception(PhysicsConfig.rung_2c())
    energy_with_drag = energy_at_interception(PhysicsConfig.rung_2d())
    assert energy_with_drag < energy_no_drag


def test_eclipse_detector_runs_full_period(config: OrbitalConfig) -> None:
    """ContinueOnEvent must let the eclipse logger span the whole period.

    With EclipseDetector's default StopOnIncreasing handler the arc would halt at
    the first umbra exit; this asserts the arc completes and that crossings pair
    up (entries and exits) over a closed period.
    """
    sun = CelestialBodyFactory.getSun()
    earth = _earth()
    prop = build_propagator(config, PhysicsConfig.rung_2c())
    logger = EventsLogger()
    prop.addEventDetector(
        logger.monitorDetector(
            EclipseDetector(sun, _SUN_RADIUS_M, earth).withHandler(ContinueOnEvent())
        )
    )
    epoch = _abs_date(config.epoch)
    period = _period(config)
    final = prop.propagate(epoch.shiftedBy(period))

    assert float(final.getDate().durationFrom(epoch)) == pytest.approx(period, abs=1.0)
    assert len(logger.getLoggedEvents()) % 2 == 0
