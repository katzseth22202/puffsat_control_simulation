"""Truth-path kernel — the regime-switched descent every JVM harness consumes (ADR 0017).

The §6.2 regime switch (B0, ADR 0008): the adaptive integrator oversteps the 200 km
event below the surface ("point is inside ellipsoid") whenever drag is too weak to
force a small step there — the low-drag dispersion tail AND, decisively, the orbits
the corrector probes (large re-phasing Δv → wildly varying perigee).  A single global
cap cannot be both fast in the long coast and safe in the stiff terminal phase, so
hand off at an altitude event: coast on the big adaptive step, then descend the last
leg on a tight cap.  The terminal cap matches the old 30 s global cap (proven safe),
while the coast runs at 600 s.  Nominal and perturbed runs share this path so the
interception miss stays common-mode.  The *fixed-step* Cowell terminal phase exists
for the executed C3a burn (:mod:`puffsat_sim.runs.terminal`, ADR 0014); the paths
here stay adaptive.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Any

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit
from org.orekit.propagation.events import AltitudeDetector
from org.orekit.propagation.events.handlers import StopOnDecreasing
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions

from puffsat_sim import mission, presets
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.constants import EARTH_RADIUS_M
from puffsat_sim.dispersion import Vec3
from puffsat_sim.propagator import build_propagator, build_propagator_from_orbit

COAST_MAX_STEP_S: float = 600.0
HANDOFF_ALT_M: float = 800_000.0  # §6.3 drag-on guard band
TERMINAL_MAX_STEP_S: float = 30.0


def to_absolute_date(dt: datetime) -> Any:
    """The Orekit UTC date for a whole-second timezone-aware datetime."""
    utc = TimeScalesFactory.getUTC()
    return AbsoluteDate(dt.year, dt.month, dt.day, dt.hour, dt.minute, float(dt.second), utc)


def earth_model() -> Any:
    """The WGS84 ellipsoid in ITRF — the altitude reference for every event detector."""
    return OneAxisEllipsoid(
        Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
        Constants.WGS84_EARTH_FLATTENING,
        FramesFactory.getITRF(IERSConventions.IERS_2010, True),
    )


@dataclass(frozen=True)
class Crossing:
    """The state read off at the 200 km interception crossing."""

    position_m: Vec3
    velocity_m_s: Vec3
    toa_s: float
    perigee_alt_m: float


def vec3(v: Any) -> Vec3:
    """A pure tuple from an Orekit Vector3D."""
    return (float(v.getX()), float(v.getY()), float(v.getZ()))


def propagate_to_interception(propagator: Any, epoch: Any, period: float, earth: Any) -> Crossing:
    """Stop the descending arc at the 200 km crossing and read off the state + perigee."""
    propagator.addEventDetector(
        AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(
            StopOnDecreasing()  # type: ignore[no-untyped-call]
        )
    )
    end = epoch.shiftedBy(period)
    state = propagator.propagate(end)
    # A stop-on-decreasing altitude event ends the arc *before* one period; reaching `end`
    # means the arc never crossed 200 km, so `state` is not an interception — fail loud
    # rather than report a full-period state as a crossing.
    if state.getDate().durationFrom(end) >= -1.0:
        raise ValueError(
            f"descending arc never crossed {mission.INTERCEPTION_ALT_M / 1e3:.0f} km "
            "within one period — no interception to read"
        )
    pv = state.getPVCoordinates()
    orbit = KeplerianOrbit(state.getOrbit())
    perigee_alt = float(orbit.getA()) * (1.0 - float(orbit.getE())) - EARTH_RADIUS_M
    return Crossing(
        position_m=vec3(pv.getPosition()),
        velocity_m_s=vec3(pv.getVelocity()),
        toa_s=float(state.getDate().durationFrom(epoch)),
        perigee_alt_m=perigee_alt,
    )


def coast_to_altitude(
    coast_prop: Any, epoch: Any, period: float, earth: Any, altitude_m: float
) -> Any:
    """Run the smooth coast on the big adaptive step, stopping at a descending altitude event.

    The §6.2 coast guard generalized to an arbitrary altitude — the C3c authority/trim
    sweep (ADR 0014 decision 5/6) stops the coast at each burn-start and trim-node
    altitude; :func:`coast_to_handoff` is the 800 km hand-off specialization.
    """
    coast_prop.addEventDetector(
        AltitudeDetector(altitude_m, earth).withHandler(
            StopOnDecreasing()  # type: ignore[no-untyped-call]
        )
    )
    return coast_prop.propagate(epoch.shiftedBy(period))


def coast_to_handoff(coast_prop: Any, epoch: Any, period: float, earth: Any) -> Any:
    """Run the smooth coast on the big adaptive step, stopping at the 800 km hand-off (§6.2)."""
    return coast_to_altitude(coast_prop, epoch, period, earth, HANDOFF_ALT_M)


def descend(
    orbit: Any,
    physics: PhysicsConfig,
    epoch: Any,
    period: float,
    earth: Any,
    maneuver: Any = None,
) -> Crossing:
    """Regime-switched descent to the 200 km crossing: coast (600 s) → 800 km → terminal (30 s).

    ``maneuver`` (B1) is an optional finite burn attached to the coast leg — the apogee
    correction fires entirely above the hand-off, so the terminal leg is unaffected.
    """
    coast = build_propagator_from_orbit(orbit, physics, COAST_MAX_STEP_S)
    if maneuver is not None:
        coast.addForceModel(maneuver)
    handoff_state = coast_to_handoff(coast, epoch, period, earth)
    terminal = build_propagator_from_orbit(handoff_state.getOrbit(), physics, TERMINAL_MAX_STEP_S)
    return propagate_to_interception(terminal, epoch, period, earth)


def apogee_state(orbital_config: OrbitalConfig) -> tuple[Any, Any, Any]:
    """Nominal deployment state at apogee (epoch, mean anomaly π): (date, position, velocity)."""
    state = build_propagator(orbital_config, presets.two_body()).getInitialState()
    pv = state.getPVCoordinates()
    return state.getDate(), pv.getPosition(), pv.getVelocity()
