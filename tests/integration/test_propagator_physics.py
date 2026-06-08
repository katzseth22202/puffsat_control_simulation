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
from datetime import datetime
from typing import Any

import pytest

from puffsat_sim import mission, presets
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.constants import SUN_RADIUS_M
from puffsat_sim.forces import AtmosphericDrag, Geopotential, Relativity
from puffsat_sim.forces.geopotential import (
    j2_apsidal_precession_rate,
    j2_nodal_regression_rate,
)
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period, wrap_to_pi

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


@pytest.fixture(scope="module")
def config() -> OrbitalConfig:
    return mission.NOMINAL_CONFIG


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
    return _final_position(config, presets.j2())


@pytest.fixture(scope="module")
def pos_j2_third_body(config: OrbitalConfig) -> tuple[float, float, float]:
    return _final_position(config, presets.j2_third_body())


@pytest.fixture(scope="module")
def pos_srp(config: OrbitalConfig) -> tuple[float, float, float]:
    return _final_position(config, presets.j2_third_body_srp())


@pytest.fixture(scope="module")
def pos_j2_relativity(config: OrbitalConfig) -> tuple[float, float, float]:
    return _final_position(config, PhysicsConfig((Geopotential(degree=2), Relativity())))


@pytest.fixture(scope="module")
def pos_geopotential_8x8(config: OrbitalConfig) -> tuple[float, float, float]:
    return _final_position(config, PhysicsConfig((Geopotential(degree=8, order=8),)))


def test_keplerian_one_period_round_trip(config: OrbitalConfig) -> None:
    """A pure-Keplerian orbit must return to its start after exactly one period."""
    prop = build_propagator(config, presets.two_body())
    pv0 = prop.getInitialState().getPVCoordinates()
    state = prop.propagate(_abs_date(config.epoch).shiftedBy(_period(config)))
    pv = state.getPVCoordinates()
    assert float(pv0.getPosition().distance(pv.getPosition())) < 1.0  # actual ~1e-8 m
    assert float(pv0.getVelocity().distance(pv.getVelocity())) < 1e-3


def test_interception_stops_at_200km(config: OrbitalConfig) -> None:
    """AltitudeDetector must stop the descending arc at the 200 km crossing."""
    earth = _earth()
    prop = build_propagator(config, presets.two_body())
    prop.addEventDetector(
        AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(StopOnDecreasing())
    )
    epoch = _abs_date(config.epoch)
    period = _period(config)
    state = prop.propagate(epoch.shiftedBy(period))

    elapsed = float(state.getDate().durationFrom(epoch))
    geodetic = earth.transform(
        state.getPVCoordinates().getPosition(), FramesFactory.getEME2000(), state.getDate()
    )
    assert geodetic.getAltitude() == pytest.approx(mission.INTERCEPTION_ALT_M, abs=10.0)
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

    prop = build_propagator(config, presets.j2())
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


def test_higher_degree_gravity_perturbs_orbit(
    pos_j2: tuple[float, float, float],
    pos_geopotential_8x8: tuple[float, float, float],
) -> None:
    """The 8×8 truth field must move the one-period endpoint by ~km vs pure J2.

    Beyond J2 the harmonics bite only near perigee, but over a pass they reach
    ~km scale at orbit level (measured ~1.8 km), so the truth model carries an
    8×8 field rather than J2 alone.
    """
    assert 500.0 < math.dist(pos_j2, pos_geopotential_8x8) < 5_000.0


def test_relativity_perturbs_orbit(
    pos_j2: tuple[float, float, float],
    pos_j2_relativity: tuple[float, float, float],
) -> None:
    """Schwarzschild relativity must shift the one-period endpoint by ~1 m.

    Differencing two runs that share J2 and integrator settings makes the
    integration error common-mode, isolating the relativistic apsidal-advance
    signal (measured ~1.04 m at the apogee endpoint).  Negligible at orbit level
    but carried in the truth model for the deferred 5 cm terminal budget.
    """
    assert 0.1 < math.dist(pos_j2, pos_j2_relativity) < 10.0


def test_drag_removes_orbital_energy(config: OrbitalConfig) -> None:
    """Specific orbital energy at the 200 km crossing must be lower with drag on.

    Drag is isolated by comparing full_force against itself with only the drag
    force removed (all conservative forces + SRP held identical), so the energy
    difference is the work drag removes — not a geopotential/relativity mismatch,
    which would confound a comparison against the lower-fidelity j2_third_body_srp.
    """
    mu = float(Constants.WGS84_EARTH_MU)
    earth = _earth()

    def energy_at_interception(physics: PhysicsConfig) -> float:
        prop = build_propagator(config, physics)
        prop.addEventDetector(
            AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(StopOnDecreasing())
        )
        state = prop.propagate(_abs_date(config.epoch).shiftedBy(_period(config)))
        pv = state.getPVCoordinates()
        r = float(pv.getPosition().getNorm())
        v = float(pv.getVelocity().getNorm())
        return 0.5 * v * v - mu / r

    full = presets.full_force()
    no_drag = PhysicsConfig(
        tuple(p for p in full.perturbations if not isinstance(p, AtmosphericDrag))
    )
    assert energy_at_interception(full) < energy_at_interception(no_drag)


def test_f10p7_drives_drag_density(config: OrbitalConfig) -> None:
    """Higher F10.7/Ap must inflate thermospheric density → more drag energy loss.

    Confirms the spec's f10p7/ap reach NRLMSISE-00 through the constant space-weather
    provider: an active-sun draw removes more orbital energy by the 200 km crossing
    than a quiet-sun draw, all else equal.  Built on full_force (swapping only the
    drag spec): its rich near-perigee dynamics keep the adaptive integrator stable
    through the descent, where a sparse J2+drag config can overshoot (design doc §6.2
    regime-switching).
    """
    mu = float(Constants.WGS84_EARTH_MU)
    earth = _earth()

    def energy_at_interception(f10p7: float, ap: float) -> float:
        physics = PhysicsConfig(
            tuple(
                AtmosphericDrag(cd_area_over_mass=0.04, f10p7=f10p7, ap=ap)
                if isinstance(p, AtmosphericDrag)
                else p
                for p in presets.full_force().perturbations
            )
        )
        prop = build_propagator(config, physics)
        prop.addEventDetector(
            AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(StopOnDecreasing())
        )
        state = prop.propagate(_abs_date(config.epoch).shiftedBy(_period(config)))
        pv = state.getPVCoordinates()
        r = float(pv.getPosition().getNorm())
        v = float(pv.getVelocity().getNorm())
        return 0.5 * v * v - mu / r

    assert energy_at_interception(250.0, 50.0) < energy_at_interception(70.0, 5.0)


def test_eclipse_detector_runs_full_period(config: OrbitalConfig) -> None:
    """ContinueOnEvent must let the eclipse logger span the whole period.

    With EclipseDetector's default StopOnIncreasing handler the arc would halt at
    the first umbra exit; this asserts the arc completes and that crossings pair
    up (entries and exits) over a closed period.
    """
    sun = CelestialBodyFactory.getSun()
    earth = _earth()
    prop = build_propagator(config, presets.j2_third_body_srp())
    logger = EventsLogger()
    prop.addEventDetector(
        logger.monitorDetector(
            EclipseDetector(sun, SUN_RADIUS_M, earth).withHandler(ContinueOnEvent())
        )
    )
    epoch = _abs_date(config.epoch)
    period = _period(config)
    final = prop.propagate(epoch.shiftedBy(period))

    assert float(final.getDate().durationFrom(epoch)) == pytest.approx(period, abs=1.0)
    assert len(logger.getLoggedEvents()) % 2 == 0
