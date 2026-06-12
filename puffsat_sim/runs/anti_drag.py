"""B3a anti-drag instrument — the JVM run for :mod:`puffsat_sim.anti_drag` (ADR 0008)."""

from __future__ import annotations

from typing import Any

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.propagation.events import AltitudeDetector
from org.orekit.propagation.events.handlers import StopOnDecreasing

from puffsat_sim import mission, presets
from puffsat_sim.anti_drag import AntiDragProfile, summarize_anti_drag
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.descent import (
    COAST_MAX_STEP_S,
    TERMINAL_MAX_STEP_S,
    coast_to_handoff,
    earth_model,
    to_absolute_date,
)
from puffsat_sim.dispersion import Vec3
from puffsat_sim.forces import AtmosphericDrag
from puffsat_sim.forces.build import Environment, to_force_models
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import build_propagator, build_propagator_from_orbit

INSTRUMENT_BELOW_ALT_M: float = 600_000.0  # §13 B3 window: the anti-drag burn runs 600 → 200 km
PUFFSAT_WET_MASS_KG: float = 25.0  # real mass for force/thrust; matches Actuator.wet_mass_kg


def drag_force_model(physics: PhysicsConfig, env: Any) -> Any:
    """The Orekit DragForce for this run's atmospheric-drag perturbation (to evaluate a_drag)."""
    for perturbation in physics.perturbations:
        if isinstance(perturbation, AtmosphericDrag):
            return to_force_models(perturbation, env)[0]
    raise ValueError("physics config carries no AtmosphericDrag perturbation")


def sample_drag_window(
    ephemeris: Any,
    start: Any,
    end: Any,
    physics: PhysicsConfig,
    earth: Any,
    epoch: Any,
    sample_dt_s: float,
) -> tuple[list[float], list[Vec3]]:
    """Sample the truth drag acceleration through the 600 → 200 km window (§13 B3).

    Times are seconds from ``epoch`` (the apogee deployment date), so they line up with
    the dates the executed maneuver segments are anchored to.
    """
    drag_force = drag_force_model(physics, Environment.build())
    params = drag_force.getParameters()
    span = float(end.durationFrom(start))

    times: list[float] = []
    accels: list[Vec3] = []
    steps = int(span / sample_dt_s)
    for k in range(steps + 1):
        date = start.shiftedBy(min(float(k) * sample_dt_s, span))
        state = ephemeris.propagate(date)
        position = state.getPVCoordinates().getPosition()
        altitude = float(earth.transform(position, state.getFrame(), date).getAltitude())
        if altitude > INSTRUMENT_BELOW_ALT_M:
            continue
        accel = drag_force.acceleration(state, params)
        times.append(float(date.durationFrom(epoch)))
        accels.append((float(accel.getX()), float(accel.getY()), float(accel.getZ())))
    return times, accels


def instrument_anti_drag(
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    mass_kg: float = PUFFSAT_WET_MASS_KG,
    sample_dt_s: float = 1.0,
) -> AntiDragProfile:
    """Instrument the nominal known-drag descent and reduce it to the anti-drag requirement (B3a).

    Feedforward cost baseline (ADR 0008/0009, perfect knowledge): descend the nominal trajectory,
    sample the truth drag acceleration through the 600 → 200 km window, and report what an
    anti-drag burn must deliver (Δv, peak thrust, peak direction-slew) — measured, not executed
    (the executed burn is C3a's :func:`puffsat_sim.runs.terminal.run_terminal_feedforward`; the
    closed loop is C3b).  Drag is evaluated at the propagator's 1 kg, which yields the real
    a_drag directly (the lumped Cd·(A/m) is the real coefficient, ADR 0009); peak thrust then
    scales by the real ``mass_kg``.
    """
    physics = presets.full_force()
    earth = earth_model()
    epoch = to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)

    apogee_orbit = (
        build_propagator(orbital_config, physics, COAST_MAX_STEP_S).getInitialState().getOrbit()
    )
    coast = build_propagator_from_orbit(apogee_orbit, physics, COAST_MAX_STEP_S)
    handoff_state = coast_to_handoff(coast, epoch, period, earth)

    terminal = build_propagator_from_orbit(handoff_state.getOrbit(), physics, TERMINAL_MAX_STEP_S)
    generator = terminal.getEphemerisGenerator()
    terminal.addEventDetector(
        AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(
            StopOnDecreasing()  # type: ignore[no-untyped-call]
        )
    )
    end_state = terminal.propagate(epoch.shiftedBy(period))
    ephemeris = generator.getGeneratedEphemeris()

    times, accels = sample_drag_window(
        ephemeris, handoff_state.getDate(), end_state.getDate(), physics, earth, epoch, sample_dt_s
    )
    return summarize_anti_drag(times, accels, mass_kg)
