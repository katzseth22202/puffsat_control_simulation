"""Integration tests for the truth-path kernel (live JVM, ADR 0017).

descent.py is the regime-switched descent seam every JVM harness consumes; these tests
cover its public surface directly instead of only incidentally through montecarlo/runs.
"""

from __future__ import annotations

import datetime as dt
import math

import pytest

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.descent import (
        HANDOFF_ALT_M,
        TERMINAL_MAX_STEP_S,
        apogee_state,
        coast_to_altitude,
        coast_to_handoff,
        descend,
        earth_model,
        propagate_to_interception,
        to_absolute_date,
        vec3,
    )
    from puffsat_sim.propagator import build_propagator, build_propagator_from_orbit
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

from org.hipparchus.geometry.euclidean.threed import Vector3D  # noqa: E402
from org.orekit.frames import FramesFactory  # noqa: E402
from org.orekit.utils import Constants  # noqa: E402

from puffsat_sim import mission, presets  # noqa: E402
from puffsat_sim.constants import EARTH_RADIUS_M, WGS84_MU  # noqa: E402
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period  # noqa: E402

pytestmark = pytest.mark.integration

_ALT_TOL_M = 1.0  # AltitudeDetector converges to sub-mm in practice; 1 m is a generous gate


def _altitude_m(position: object, date: object, earth: object) -> float:
    frame = FramesFactory.getEME2000()
    return float(earth.transform(position, frame, date).getAltitude())


def _vector3d(v: tuple[float, float, float]) -> Vector3D:
    return Vector3D(v[0], v[1], v[2])


def _nominal_period() -> float:
    semi_major, _ = keplerian_elements(mission.PERIGEE_ALT_M, mission.APOGEE_ALT_M)
    return keplerian_period(semi_major)


def test_to_absolute_date_round_trip() -> None:
    """Whole-second offsets from the Python side land at the same offset on the Orekit side."""
    base = to_absolute_date(mission.EPOCH)
    same = to_absolute_date(mission.EPOCH)
    later = to_absolute_date(mission.EPOCH + dt.timedelta(seconds=3661))

    assert same.durationFrom(base) == 0.0
    assert later.durationFrom(base) == pytest.approx(3661.0)


def test_earth_model_is_wgs84_in_itrf() -> None:
    """The altitude reference every event detector shares is the standard WGS84 ellipsoid."""
    earth = earth_model()

    assert earth.getEquatorialRadius() == pytest.approx(Constants.WGS84_EARTH_EQUATORIAL_RADIUS)
    assert earth.getFlattening() == pytest.approx(Constants.WGS84_EARTH_FLATTENING)


def test_apogee_state_is_at_apogee() -> None:
    """Mean anomaly π puts the deployment state at the design apogee radius, moving purely
    tangentially (zero radial velocity component)."""
    _, pos, vel = apogee_state(mission.NOMINAL_CONFIG)
    p = vec3(pos)
    v = vec3(vel)

    r = math.dist((0.0, 0.0, 0.0), p)
    speed = math.dist((0.0, 0.0, 0.0), v)
    semi_major, _ = keplerian_elements(mission.PERIGEE_ALT_M, mission.APOGEE_ALT_M)
    r_a = EARTH_RADIUS_M + mission.APOGEE_ALT_M
    expected_speed = math.sqrt(WGS84_MU * (2.0 / r_a - 1.0 / semi_major))

    assert r == pytest.approx(r_a, abs=1.0)
    assert speed == pytest.approx(expected_speed, rel=1e-9)
    radial_component = sum(p[i] * v[i] for i in range(3)) / r
    assert abs(radial_component) < 1e-6


def test_coast_to_altitude_stops_at_an_arbitrary_altitude() -> None:
    """coast_to_altitude generalizes past the 800 km hand-off (C3c uses it at trim-node alts)."""
    earth = earth_model()
    epoch = to_absolute_date(mission.EPOCH)
    period = _nominal_period()
    orbit = build_propagator(mission.NOMINAL_CONFIG, presets.j2()).getInitialState().getOrbit()
    coast = build_propagator_from_orbit(orbit, presets.j2(), 600.0)

    target_alt_m = 30_000_000.0
    state = coast_to_altitude(coast, epoch, period, earth, target_alt_m)

    assert _altitude_m(state.getPVCoordinates().getPosition(), state.getDate(), earth) == (
        pytest.approx(target_alt_m, abs=_ALT_TOL_M)
    )


def test_coast_to_handoff_stops_at_800km() -> None:
    """The §6.2 hand-off specialization stops exactly at HANDOFF_ALT_M."""
    earth = earth_model()
    epoch = to_absolute_date(mission.EPOCH)
    period = _nominal_period()
    orbit = build_propagator(mission.NOMINAL_CONFIG, presets.j2()).getInitialState().getOrbit()
    coast = build_propagator_from_orbit(orbit, presets.j2(), 600.0)

    state = coast_to_handoff(coast, epoch, period, earth)

    assert _altitude_m(state.getPVCoordinates().getPosition(), state.getDate(), earth) == (
        pytest.approx(HANDOFF_ALT_M, abs=_ALT_TOL_M)
    )


def test_propagate_to_interception_stops_at_200km() -> None:
    """The descending-altitude event fires at the 200 km control target, not the first crossing
    below it, and reports a perigee below the interception altitude."""
    earth = earth_model()
    epoch = to_absolute_date(mission.EPOCH)
    period = _nominal_period()
    orbit = build_propagator(mission.NOMINAL_CONFIG, presets.j2()).getInitialState().getOrbit()
    handoff = coast_to_handoff(
        build_propagator_from_orbit(orbit, presets.j2(), 600.0), epoch, period, earth
    )
    terminal = build_propagator_from_orbit(handoff.getOrbit(), presets.j2(), TERMINAL_MAX_STEP_S)

    crossing = propagate_to_interception(terminal, epoch, period, earth)

    crossing_alt = _altitude_m(
        _vector3d(crossing.position_m), epoch.shiftedBy(crossing.toa_s), earth
    )
    assert crossing_alt == pytest.approx(mission.INTERCEPTION_ALT_M, abs=_ALT_TOL_M)
    assert crossing.toa_s > 0.0
    assert 0.0 < crossing.perigee_alt_m < mission.INTERCEPTION_ALT_M


def test_descend_two_body_reproduces_the_design_perigee() -> None:
    """With no perturbations active, the regime switch introduces no error: the descent
    reaches exactly the design 200 km crossing and the design 50 km perigee."""
    earth = earth_model()
    epoch = to_absolute_date(mission.EPOCH)
    period = _nominal_period()
    orbit = (
        build_propagator(mission.NOMINAL_CONFIG, presets.two_body()).getInitialState().getOrbit()
    )

    crossing = descend(orbit, presets.two_body(), epoch, period, earth)

    crossing_alt = _altitude_m(
        _vector3d(crossing.position_m), epoch.shiftedBy(crossing.toa_s), earth
    )
    assert crossing_alt == pytest.approx(mission.INTERCEPTION_ALT_M, abs=_ALT_TOL_M)
    assert crossing.perigee_alt_m == pytest.approx(mission.PERIGEE_ALT_M, abs=1e-3)


def test_descend_matches_continuous_fine_step_propagation() -> None:
    """The coast(600s)->800km->terminal(30s) regime switch is a performance trick, not a
    physics change: it must reproduce a single continuous propagation held at the terminal's
    tight step for the whole arc, to tight numerical agreement (ADR 0017 / §6.2 B0)."""
    earth = earth_model()
    epoch = to_absolute_date(mission.EPOCH)
    period = _nominal_period()
    orbit = build_propagator(mission.NOMINAL_CONFIG, presets.j2()).getInitialState().getOrbit()

    regime_switched = descend(orbit, presets.j2(), epoch, period, earth)

    continuous_prop = build_propagator_from_orbit(orbit, presets.j2(), TERMINAL_MAX_STEP_S)
    continuous = propagate_to_interception(continuous_prop, epoch, period, earth)

    assert math.dist(regime_switched.position_m, continuous.position_m) < 1.0
    assert math.dist(regime_switched.velocity_m_s, continuous.velocity_m_s) < 1e-3
    assert regime_switched.toa_s == pytest.approx(continuous.toa_s, abs=1e-3)
    assert regime_switched.perigee_alt_m == pytest.approx(continuous.perigee_alt_m, abs=1.0)
