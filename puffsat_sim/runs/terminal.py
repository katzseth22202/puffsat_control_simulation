"""C3a executed terminal feedforward — the JVM run for :mod:`puffsat_sim.terminal` (ADR 0014)."""

from __future__ import annotations

import math

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.attitudes import FrameAlignedProvider
from org.orekit.forces.maneuvers import ConstantThrustManeuver
from org.orekit.frames import FramesFactory

from org.hipparchus.geometry.euclidean.threed import Vector3D

from puffsat_sim import mission, presets
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.descent import (
    COAST_MAX_STEP_S,
    TERMINAL_MAX_STEP_S,
    coast_to_handoff,
    earth_model,
    propagate_to_interception,
    to_absolute_date,
)
from puffsat_sim.forces import AtmosphericDrag
from puffsat_sim.montecarlo import BURN_ISP_SENTINEL_S, scale_thrust_for_propagator_mass
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import (
    build_fixed_step_propagator_from_orbit,
    build_propagator,
    build_propagator_from_orbit,
)
from puffsat_sim.runs.anti_drag import PUFFSAT_WET_MASS_KG, sample_drag_window
from puffsat_sim.terminal import (
    TerminalFeedforwardFinding,
    format_terminal_feedforward,
    plan_feedforward,
)

# C3a fixed-step terminal (ADR 0014 decision 4): the integrator step equals the control
# period, so every zero-order-hold command boundary lands exactly on the integrator grid.
TERMINAL_FIXED_STEP_S: float = 1.0


def physics_without_drag(physics: PhysicsConfig) -> PhysicsConfig:
    """The same force model with the drag perturbation removed (the C3a drag-free reference)."""
    return PhysicsConfig(
        tuple(p for p in physics.perturbations if not isinstance(p, AtmosphericDrag))
    )


def run_terminal_feedforward(
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    control_period_s: float = TERMINAL_FIXED_STEP_S,
    step_s: float = TERMINAL_FIXED_STEP_S,
    sample_dt_s: float = 1.0,
) -> TerminalFeedforwardFinding:
    """Execute B3a's anti-drag feedforward as a real ZOH burn on the fixed-step terminal (C3a).

    ADR 0014 decisions 4/6, in order: (1) descend the unburned terminal leg from the 800 km
    hand-off on both the proven adaptive-30 s config and the fixed-step Cowell — their
    crossing separation is the equivalence pin, measured before any burn; (2) sample the
    truth drag profile on the unburned fixed-step descent and ZOH-plan it
    (:func:`puffsat_sim.terminal.plan_feedforward`, real units: 25 kg, 400 mN cap);
    (3) re-descend the same hand-off state with the commands attached as finite maneuver
    segments, thrust scaled to the 1 kg propagator mass at the sentinel Isp as in B1.
    The drag-free descent of the same hand-off state is the reference: the unburned
    distance to it is the drag displacement (the disease), the burned distance is the
    executed residual (what the open-loop feedforward leaves for C3b's feedback).
    """
    physics = presets.full_force()
    earth = earth_model()
    epoch = to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)
    frame = FramesFactory.getEME2000()

    apogee_orbit = (
        build_propagator(orbital_config, physics, COAST_MAX_STEP_S).getInitialState().getOrbit()
    )
    coast = build_propagator_from_orbit(apogee_orbit, physics, COAST_MAX_STEP_S)
    handoff = coast_to_handoff(coast, epoch, period, earth)

    adaptive = propagate_to_interception(
        build_propagator_from_orbit(handoff.getOrbit(), physics, TERMINAL_MAX_STEP_S),
        epoch,
        period,
        earth,
    )

    fixed_prop = build_fixed_step_propagator_from_orbit(handoff.getOrbit(), physics, step_s)
    generator = fixed_prop.getEphemerisGenerator()
    unburned = propagate_to_interception(fixed_prop, epoch, period, earth)
    ephemeris = generator.getGeneratedEphemeris()

    dragfree = propagate_to_interception(
        build_fixed_step_propagator_from_orbit(
            handoff.getOrbit(), physics_without_drag(physics), step_s
        ),
        epoch,
        period,
        earth,
    )

    times, accels = sample_drag_window(
        ephemeris,
        handoff.getDate(),
        epoch.shiftedBy(unburned.toa_s),
        physics,
        earth,
        epoch,
        sample_dt_s,
    )
    plan = plan_feedforward(
        times, accels, mass_kg=PUFFSAT_WET_MASS_KG, control_period_s=control_period_s
    )

    burned_prop = build_fixed_step_propagator_from_orbit(handoff.getOrbit(), physics, step_s)
    for cmd in plan.commands:
        if cmd.thrust_n <= 0.0:
            continue
        burned_prop.addForceModel(
            ConstantThrustManeuver(
                epoch.shiftedBy(cmd.start_s),
                cmd.duration_s,
                scale_thrust_for_propagator_mass(cmd.thrust_n, PUFFSAT_WET_MASS_KG),
                BURN_ISP_SENTINEL_S,
                FrameAlignedProvider(frame),
                Vector3D(cmd.direction[0], cmd.direction[1], cmd.direction[2]),
            )
        )
    burned = propagate_to_interception(burned_prop, epoch, period, earth)

    return TerminalFeedforwardFinding(
        plan=plan,
        equivalence_pin_m=math.dist(unburned.position_m, adaptive.position_m),
        equivalence_pin_toa_s=unburned.toa_s - adaptive.toa_s,
        drag_displacement_m=math.dist(unburned.position_m, dragfree.position_m),
        executed_residual_m=math.dist(burned.position_m, dragfree.position_m),
        executed_residual_toa_s=burned.toa_s - dragfree.toa_s,
    )


def terminal_feedforward_report(
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    control_period_s: float = TERMINAL_FIXED_STEP_S,
) -> str:
    """Run the C3a executed-feedforward measurement and format the one-screen report."""
    return format_terminal_feedforward(
        run_terminal_feedforward(orbital_config, control_period_s=control_period_s)
    )
