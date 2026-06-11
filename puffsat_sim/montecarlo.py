"""Monte Carlo / open-loop dispersion harness — the JVM-side run loop (ADR 0002).

``run_ensemble`` samples per-run inputs (:mod:`puffsat_sim.dispersion`), builds the
perturbed full-force run, applies the injection Δv to the apogee deployment state,
propagates to the 200 km interception crossing, and records the miss (in the
nominal-crossing RTN frame), the time-of-arrival error, and the osculating perigee.

The Stage-1 capstone (design doc §13) is this harness with ``control=None``; Rung D
supplies a controller through the same hook (§14.1).  Per-run replay (§14.2):
``replay_inputs(master_seed, spec, run_index)`` reconstructs any run's draws.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.attitudes import FrameAlignedProvider
from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.forces.maneuvers import ConstantThrustManeuver
from org.orekit.frames import FramesFactory
from org.orekit.orbits import CartesianOrbit, KeplerianOrbit
from org.orekit.propagation.events import AltitudeDetector
from org.orekit.propagation.events.handlers import StopOnDecreasing
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions, TimeStampedPVCoordinates

from org.hipparchus.geometry.euclidean.threed import Vector3D

from puffsat_sim import mission, presets
from puffsat_sim.actuator import Actuator, plan_burn
from puffsat_sim.anti_drag import AntiDragProfile, summarize_anti_drag
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.constants import EARTH_RADIUS_M
from puffsat_sim.control import (
    ControlPlan,
    Controller,
    PredictFn,
    Target,
    passes_toa_gate,
    solve_apogee_correction,
)
from puffsat_sim.dispersion import (
    Basis,
    DispersionSpec,
    RunInputs,
    Vec3,
    replay_inputs,
    rtn_basis,
    rtn_components,
    rtn_to_cartesian,
    summarize,
)
from puffsat_sim.coeff_requirement import (
    format_coeff_requirement,
    summarize_coeff_requirement,
)
from puffsat_sim.forces import (
    AtmosphericDrag,
    Geopotential,
    Relativity,
    SolarRadiation,
    ThirdBody,
)
from puffsat_sim.forces.build import Environment, to_force_models
from puffsat_sim.nav_feasibility import (
    NavFeasibilityResult,
    NavFeasibilitySpec,
    NavValidationOutcome,
    format_nav_feasibility,
    format_nav_validation,
    sweep_nav_feasibility,
    validate_cell,
)
from puffsat_sim.navigation import (
    NavSweepResult,
    NavSweepSpec,
    Vec6,
    format_nav_requirement,
    nav_grid_offsets,
    summarize_nav_requirement,
)
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import (
    build_fixed_step_propagator_from_orbit,
    build_propagator,
    build_propagator_from_orbit,
)
from puffsat_sim.terminal import (
    TerminalFeedforwardFinding,
    format_terminal_feedforward,
    plan_feedforward,
)
from puffsat_sim.records import EnsembleResult, RunRecord
from puffsat_sim.sink import append_record, plan_resume, read_records
from puffsat_sim.sweep import SweepResult, SweepSpec, grid_inputs

# Regime-switched descent (B0, ADR 0008 / design §6.2): the adaptive integrator oversteps
# the 200 km event below the surface ("point is inside ellipsoid") whenever drag is too
# weak to force a small step there — the low-drag dispersion tail AND, decisively, the
# orbits the corrector probes (large re-phasing Δv → wildly varying perigee). A single
# global cap cannot be both fast in the long coast and safe in the stiff terminal phase,
# so hand off at an altitude event: coast on the big adaptive step, then descend the last
# leg on a tight cap.  The terminal cap matches the old 30 s global cap (proven safe),
# while the coast runs at 600 s — recovering the coast tax.  Nominal and perturbed runs
# share this path so the interception miss stays common-mode.  The *fixed-step* Cowell
# terminal phase exists for the executed C3a burn (run_terminal_feedforward, ADR 0014);
# the dispersion path here stays adaptive.
_COAST_MAX_STEP_S: float = 600.0
_HANDOFF_ALT_M: float = 800_000.0  # §6.3 drag-on guard band
_TERMINAL_MAX_STEP_S: float = 30.0

# Finite-burn execution (B1, ADR 0008): the propagator runs at a fictitious 1 kg so the
# lumped Cd·(A/m) / Cr·(A/m) scale drag/SRP correctly, so the burn thrust is scaled to that
# mass (F·m_p/m_wet) — reproducing the real a=F/m and burn duration of the 25 kg / 400 mN
# actuator — and fires at a sentinel Isp so the executed arc is constant-mass (Isp-free
# trajectory).  Real propellant is the pure Tsiolkovsky transform at the actuator's Isp
# (puffsat_sim.actuator), per ADR 0004 decision 2.  Real-mass depletion coupled to descent
# drag is deferred to B3 (its large anti-drag burn is the first consumer).
_PROPAGATOR_MASS_KG: float = 1.0
_BURN_ISP_SENTINEL_S: float = 1.0e12


def physics_from_inputs(
    inputs: RunInputs, geopotential_degree: int = 8, geopotential_order: int = 8
) -> PhysicsConfig:
    """Full-force truth config carrying this run's drawn coefficients and space weather."""
    return PhysicsConfig(
        (
            Geopotential(degree=geopotential_degree, order=geopotential_order),
            ThirdBody(),
            SolarRadiation(cr_area_over_mass=inputs.cr_area_over_mass),
            AtmosphericDrag(
                cd_area_over_mass=inputs.cd_area_over_mass, f10p7=inputs.f10p7, ap=inputs.ap
            ),
            Relativity(),
        )
    )


def _to_absolute_date(dt: datetime) -> Any:
    utc = TimeScalesFactory.getUTC()
    return AbsoluteDate(dt.year, dt.month, dt.day, dt.hour, dt.minute, float(dt.second), utc)


def _earth() -> Any:
    return OneAxisEllipsoid(
        Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
        Constants.WGS84_EARTH_FLATTENING,
        FramesFactory.getITRF(IERSConventions.IERS_2010, True),
    )


@dataclass(frozen=True)
class _Crossing:
    position_m: Vec3
    velocity_m_s: Vec3
    toa_s: float
    perigee_alt_m: float


def _vec3(v: Any) -> Vec3:
    return (float(v.getX()), float(v.getY()), float(v.getZ()))


def _propagate_to_interception(propagator: Any, epoch: Any, period: float, earth: Any) -> _Crossing:
    """Stop the descending arc at the 200 km crossing and read off the state + perigee."""
    propagator.addEventDetector(
        AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(
            StopOnDecreasing()  # type: ignore[no-untyped-call]
        )
    )
    state = propagator.propagate(epoch.shiftedBy(period))
    pv = state.getPVCoordinates()
    orbit = KeplerianOrbit(state.getOrbit())
    perigee_alt = float(orbit.getA()) * (1.0 - float(orbit.getE())) - EARTH_RADIUS_M
    return _Crossing(
        position_m=_vec3(pv.getPosition()),
        velocity_m_s=_vec3(pv.getVelocity()),
        toa_s=float(state.getDate().durationFrom(epoch)),
        perigee_alt_m=perigee_alt,
    )


def _coast_to_handoff(coast_prop: Any, epoch: Any, period: float, earth: Any) -> Any:
    """Run the smooth coast on the big adaptive step, stopping at the 800 km hand-off (§6.2)."""
    coast_prop.addEventDetector(
        AltitudeDetector(_HANDOFF_ALT_M, earth).withHandler(
            StopOnDecreasing()  # type: ignore[no-untyped-call]
        )
    )
    return coast_prop.propagate(epoch.shiftedBy(period))


def _descend(
    orbit: Any,
    physics: PhysicsConfig,
    epoch: Any,
    period: float,
    earth: Any,
    maneuver: Any = None,
) -> _Crossing:
    """Regime-switched descent to the 200 km crossing: coast (600 s) → 800 km → terminal (30 s).

    ``maneuver`` (B1) is an optional finite burn attached to the coast leg — the apogee
    correction fires entirely above the hand-off, so the terminal leg is unaffected.
    """
    coast = build_propagator_from_orbit(orbit, physics, _COAST_MAX_STEP_S)
    if maneuver is not None:
        coast.addForceModel(maneuver)
    handoff_state = _coast_to_handoff(coast, epoch, period, earth)
    terminal = build_propagator_from_orbit(handoff_state.getOrbit(), physics, _TERMINAL_MAX_STEP_S)
    return _propagate_to_interception(terminal, epoch, period, earth)


def _finite_burn_maneuver(actuator: Actuator, correction_rtn: Vec3, ctx: _RunContext) -> Any:
    """The ConstantThrustManeuver that executes the corrector's Δv as a finite burn (B1).

    Thrust is scaled to the propagator's fictitious 1 kg (F·m_p/m_wet) so a=F/m and the burn
    duration match the real 25 kg / 400 mN actuator; the sentinel Isp keeps the executed arc
    constant-mass.  Direction is the inertial Δv, held fixed by a frame-aligned attitude.
    """
    burn = plan_burn(actuator, correction_rtn)
    corr_eme = rtn_to_cartesian(correction_rtn, ctx.apo_basis)
    norm = (corr_eme[0] ** 2 + corr_eme[1] ** 2 + corr_eme[2] ** 2) ** 0.5  # nonzero: caller guards
    direction = Vector3D(corr_eme[0] / norm, corr_eme[1] / norm, corr_eme[2] / norm)
    thrust_eff = _PROPAGATOR_MASS_KG * actuator.max_thrust_n / actuator.wet_mass_kg
    return ConstantThrustManeuver(
        ctx.apo_date,
        burn.duration_s,
        thrust_eff,
        _BURN_ISP_SENTINEL_S,
        FrameAlignedProvider(ctx.frame),
        direction,
    )


_INSTRUMENT_BELOW_ALT_M: float = 600_000.0  # §13 B3 window: the anti-drag burn runs 600 → 200 km
_PUFFSAT_WET_MASS_KG: float = 25.0  # real mass for force/thrust; matches Actuator.wet_mass_kg


def _drag_force_model(physics: PhysicsConfig, env: Any) -> Any:
    """The Orekit DragForce for this run's atmospheric-drag perturbation (to evaluate a_drag)."""
    for perturbation in physics.perturbations:
        if isinstance(perturbation, AtmosphericDrag):
            return to_force_models(perturbation, env)[0]
    raise ValueError("physics config carries no AtmosphericDrag perturbation")


def instrument_anti_drag(
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    mass_kg: float = _PUFFSAT_WET_MASS_KG,
    sample_dt_s: float = 1.0,
) -> AntiDragProfile:
    """Instrument the nominal known-drag descent and reduce it to the anti-drag requirement (B3a).

    Feedforward cost baseline (ADR 0008/0009, perfect knowledge): descend the nominal trajectory,
    sample the truth drag acceleration through the 600 → 200 km window, and report what an
    anti-drag burn must deliver (Δv, peak thrust, peak direction-slew) — measured, not executed
    (the executed burn is C3a's :func:`run_terminal_feedforward`; the closed loop is C3b).  Drag
    is evaluated at the propagator's 1 kg, which yields the real a_drag directly (the lumped
    Cd·(A/m) is the real coefficient, ADR 0009); peak thrust then scales by the real ``mass_kg``.
    """
    physics = presets.full_force()
    earth = _earth()
    epoch = _to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)

    apogee_orbit = (
        build_propagator(orbital_config, physics, _COAST_MAX_STEP_S).getInitialState().getOrbit()
    )
    coast = build_propagator_from_orbit(apogee_orbit, physics, _COAST_MAX_STEP_S)
    handoff_state = _coast_to_handoff(coast, epoch, period, earth)

    terminal = build_propagator_from_orbit(handoff_state.getOrbit(), physics, _TERMINAL_MAX_STEP_S)
    generator = terminal.getEphemerisGenerator()
    terminal.addEventDetector(
        AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(
            StopOnDecreasing()  # type: ignore[no-untyped-call]
        )
    )
    end_state = terminal.propagate(epoch.shiftedBy(period))
    ephemeris = generator.getGeneratedEphemeris()

    times, accels = _sample_drag_window(
        ephemeris, handoff_state.getDate(), end_state.getDate(), physics, earth, epoch, sample_dt_s
    )
    return summarize_anti_drag(times, accels, mass_kg)


def _sample_drag_window(
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
    drag_force = _drag_force_model(physics, Environment.build())
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
        if altitude > _INSTRUMENT_BELOW_ALT_M:
            continue
        accel = drag_force.acceleration(state, params)
        times.append(float(date.durationFrom(epoch)))
        accels.append((float(accel.getX()), float(accel.getY()), float(accel.getZ())))
    return times, accels


# C3a fixed-step terminal (ADR 0014 decision 4): the integrator step equals the control
# period, so every zero-order-hold command boundary lands exactly on the integrator grid.
_TERMINAL_FIXED_STEP_S: float = 1.0


def _physics_without_drag(physics: PhysicsConfig) -> PhysicsConfig:
    """The same force model with the drag perturbation removed (the C3a drag-free reference)."""
    return PhysicsConfig(
        tuple(p for p in physics.perturbations if not isinstance(p, AtmosphericDrag))
    )


def run_terminal_feedforward(
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    control_period_s: float = _TERMINAL_FIXED_STEP_S,
    step_s: float = _TERMINAL_FIXED_STEP_S,
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
    earth = _earth()
    epoch = _to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)
    frame = FramesFactory.getEME2000()

    apogee_orbit = (
        build_propagator(orbital_config, physics, _COAST_MAX_STEP_S).getInitialState().getOrbit()
    )
    coast = build_propagator_from_orbit(apogee_orbit, physics, _COAST_MAX_STEP_S)
    handoff = _coast_to_handoff(coast, epoch, period, earth)

    adaptive = _propagate_to_interception(
        build_propagator_from_orbit(handoff.getOrbit(), physics, _TERMINAL_MAX_STEP_S),
        epoch,
        period,
        earth,
    )

    fixed_prop = build_fixed_step_propagator_from_orbit(handoff.getOrbit(), physics, step_s)
    generator = fixed_prop.getEphemerisGenerator()
    unburned = _propagate_to_interception(fixed_prop, epoch, period, earth)
    ephemeris = generator.getGeneratedEphemeris()

    dragfree = _propagate_to_interception(
        build_fixed_step_propagator_from_orbit(
            handoff.getOrbit(), _physics_without_drag(physics), step_s
        ),
        epoch,
        period,
        earth,
    )

    times, accels = _sample_drag_window(
        ephemeris,
        handoff.getDate(),
        epoch.shiftedBy(unburned.toa_s),
        physics,
        earth,
        epoch,
        sample_dt_s,
    )
    plan = plan_feedforward(
        times, accels, mass_kg=_PUFFSAT_WET_MASS_KG, control_period_s=control_period_s
    )

    burned_prop = build_fixed_step_propagator_from_orbit(handoff.getOrbit(), physics, step_s)
    for cmd in plan.commands:
        if cmd.thrust_n <= 0.0:
            continue
        burned_prop.addForceModel(
            ConstantThrustManeuver(
                epoch.shiftedBy(cmd.start_s),
                cmd.duration_s,
                cmd.thrust_n * _PROPAGATOR_MASS_KG / _PUFFSAT_WET_MASS_KG,
                _BURN_ISP_SENTINEL_S,
                FrameAlignedProvider(frame),
                Vector3D(cmd.direction[0], cmd.direction[1], cmd.direction[2]),
            )
        )
    burned = _propagate_to_interception(burned_prop, epoch, period, earth)

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
    control_period_s: float = _TERMINAL_FIXED_STEP_S,
) -> str:
    """Run the C3a executed-feedforward measurement and format the one-screen report."""
    return format_terminal_feedforward(
        run_terminal_feedforward(orbital_config, control_period_s=control_period_s)
    )


def _apogee_state(orbital_config: OrbitalConfig) -> tuple[Any, Any, Any]:
    """Nominal deployment state at apogee (epoch, mean anomaly π): (date, position, velocity)."""
    state = build_propagator(orbital_config, presets.two_body()).getInitialState()
    pv = state.getPVCoordinates()
    return state.getDate(), pv.getPosition(), pv.getVelocity()


@dataclass(frozen=True)
class _RunContext:
    """Per-ensemble constants shared by every run (built once before the loop)."""

    apo_date: Any
    apo_pos: Any
    apo_vel: Any
    apo_basis: Basis
    frame: Any
    mu: float
    epoch: Any
    period: float
    earth: Any
    nominal: _Crossing
    nominal_basis: Basis
    target: Target


def _run_record(
    ctx: _RunContext,
    inputs: RunInputs,
    control: Controller | None,
    toa_window_s: float | None = None,
    actuator: Actuator | None = None,
    nav_offset_rtn6: Vec6 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
) -> RunRecord:
    """Propagate one run: apply injection, solve+execute the control plan, record the miss.

    ``predict`` (the corrector's onboard model) and ``execute`` (truth) are the same
    full-force physics at Rung A (ADR 0003), so a converged plan lands the recorded
    crossing on the nominal aim to machine precision.  The injection Δv is baked into
    the closure, so the corrector solves for the *correction* alone, starting from zero.

    ``actuator`` (B1, ADR 0008) makes ``execute`` a finite burn while ``predict`` stays
    impulsive: the corrector solves the impulsive commanded Δv, truth fires it as a finite
    maneuver, and the residual miss is the actuator-realism erosion.  ``None`` keeps the
    Rung-A impulsive execution (commanded == applied).

    ``nav_offset_rtn6`` (C0, ADR 0011) is a predict-side apogee-RTN navigation-error offset
    (position R/T/N then velocity R/T/N): the corrector plans from ``x_true + offset`` while
    execute stays on truth, so the residual miss is the apogee→crossing sensitivity Φ times
    the nav error.  Zero-default leaves A/B untouched.

    ``toa_window_s`` (default off) is the A3 spurious-far-root gate: a converged plan whose
    crossing falls outside ±window of the nominal ToA is recorded non-converged (ADR 0007
    decision 3iii).  ``run_ensemble`` / the capstone leave it ``None``.
    """
    physics = physics_from_inputs(inputs)
    injection_dv_eme = rtn_to_cartesian(inputs.dv_rtn_m_s, ctx.apo_basis)
    injection = Vector3D(injection_dv_eme[0], injection_dv_eme[1], injection_dv_eme[2])

    def make_crossing(
        correction_rtn: Vec3, apo_offset_rtn6: Vec6 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ) -> _Crossing:
        pos_off = rtn_to_cartesian(
            (apo_offset_rtn6[0], apo_offset_rtn6[1], apo_offset_rtn6[2]), ctx.apo_basis
        )
        vel_off = rtn_to_cartesian(
            (apo_offset_rtn6[3], apo_offset_rtn6[4], apo_offset_rtn6[5]), ctx.apo_basis
        )
        corr_eme = rtn_to_cartesian(correction_rtn, ctx.apo_basis)
        position = ctx.apo_pos.add(Vector3D(pos_off[0], pos_off[1], pos_off[2]))
        vel = (
            ctx.apo_vel.add(injection)
            .add(Vector3D(corr_eme[0], corr_eme[1], corr_eme[2]))
            .add(Vector3D(vel_off[0], vel_off[1], vel_off[2]))
        )
        orbit = CartesianOrbit(
            TimeStampedPVCoordinates(ctx.apo_date, position, vel), ctx.frame, ctx.mu
        )
        return _descend(orbit, physics, ctx.epoch, ctx.period, ctx.earth)

    if control is None:
        plan = ControlPlan(actions=(), converged=True, iterations=0)
    else:
        plan = control(lambda c: make_crossing(c, nav_offset_rtn6).position_m, ctx.target)

    # Execute the commanded plan against truth.  A1's single action is at the apogee
    # node (elapsed_s=0), so it folds into the initial velocity; downstream multi-node
    # execution (ImpulseManeuver events) is an A2 addition.
    applied_rtn: Vec3 = plan.actions[0].dv_rtn_m_s if plan.actions else (0.0, 0.0, 0.0)
    if actuator is not None and plan.actions and plan.actions[0].dv_mag_m_s > 0.0:
        # B1: fire the correction as a finite burn (injection alone folds into the velocity).
        orbit = CartesianOrbit(
            TimeStampedPVCoordinates(ctx.apo_date, ctx.apo_pos, ctx.apo_vel.add(injection)),
            ctx.frame,
            ctx.mu,
        )
        crossing = _descend(
            orbit,
            physics,
            ctx.epoch,
            ctx.period,
            ctx.earth,
            maneuver=_finite_burn_maneuver(actuator, applied_rtn, ctx),
        )
    else:
        crossing = make_crossing(applied_rtn)

    miss_vec: Vec3 = (
        crossing.position_m[0] - ctx.nominal.position_m[0],
        crossing.position_m[1] - ctx.nominal.position_m[1],
        crossing.position_m[2] - ctx.nominal.position_m[2],
    )
    toa_miss_s = crossing.toa_s - ctx.nominal.toa_s
    return RunRecord(
        inputs=inputs,
        miss_rtn_m=rtn_components(miss_vec, ctx.nominal_basis),
        toa_miss_s=toa_miss_s,
        perigee_alt_m=crossing.perigee_alt_m,
        crossing_position_m=crossing.position_m,
        crossing_velocity_m_s=crossing.velocity_m_s,
        control_log=plan.actions,
        total_dv_m_s=plan.total_dv_m_s,
        converged=passes_toa_gate(plan.converged, toa_miss_s, toa_window_s),
        iterations=plan.iterations,
    )


def _build_context(orbital_config: OrbitalConfig) -> _RunContext:
    """Build the per-run constants shared by every run (the nominal crossing + apogee frame).

    Factored out of ``run_ensemble`` so the deterministic ``run_sweep`` reuses the exact same
    nominal-crossing setup and ``_RunContext`` (ADR 0007 decision 5).
    """
    earth = _earth()
    epoch = _to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)
    frame = FramesFactory.getEME2000()
    mu: float = Constants.WGS84_EARTH_MU

    # Nominal (unperturbed) crossing — the reference the miss is measured against and
    # the corrector's target.
    nominal_orbit = (
        build_propagator(orbital_config, presets.full_force(), _COAST_MAX_STEP_S)
        .getInitialState()
        .getOrbit()
    )
    nominal = _descend(nominal_orbit, presets.full_force(), epoch, period, earth)
    nominal_basis: Basis = rtn_basis(nominal.position_m, nominal.velocity_m_s)

    apo_date, apo_pos, apo_vel = _apogee_state(orbital_config)
    apo_basis: Basis = rtn_basis(_vec3(apo_pos), _vec3(apo_vel))

    return _RunContext(
        apo_date=apo_date,
        apo_pos=apo_pos,
        apo_vel=apo_vel,
        apo_basis=apo_basis,
        frame=frame,
        mu=mu,
        epoch=epoch,
        period=period,
        earth=earth,
        nominal=nominal,
        nominal_basis=nominal_basis,
        target=Target(nominal.position_m),
    )


def run_ensemble(
    spec: DispersionSpec,
    n: int,
    master_seed: int,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    control: Controller | None = None,
    sink_path: Path | None = None,
    actuator: Actuator | None = None,
) -> EnsembleResult:
    """Run a dispersion ensemble and aggregate it.

    ``control`` is the §14.1 hook (ADR 0003): ``None`` is the open-loop capstone; a
    ``Controller`` (e.g. the Rung A1 ``solve_apogee_correction``) closes the loop, each
    run solving and executing its plan against the same full-force truth.

    ``actuator`` (B1, ADR 0008) executes each commanded Δv as a finite mass-depleting burn
    instead of an impulse (``predict`` stays impulsive); the residual miss is then the
    actuator-realism erosion.  Requires a ``control``; ``None`` keeps impulsive execution.

    ``sink_path`` enables run-granular checkpoint/resume: completed records stream to a
    JSONL sink keyed by ``run_index``; on restart only the missing indices are run
    (the present ones must match this ``master_seed``/``spec``).  The caller must resume
    with the same ``control`` — inputs are control-independent, so it is not auto-checked.
    """
    ctx = _build_context(orbital_config)

    reuse: dict[int, RunRecord] = {}
    todo = list(range(n))
    if sink_path is not None:
        reuse, todo = plan_resume(read_records(sink_path), master_seed, spec, n)
    records_by_index: dict[int, RunRecord] = dict(reuse)
    for i in todo:
        record = _run_record(ctx, replay_inputs(master_seed, spec, i), control, actuator=actuator)
        if sink_path is not None:
            append_record(sink_path, record)
        records_by_index[i] = record
    records = [records_by_index[i] for i in range(n)]

    stats = summarize(
        np.array([r.miss_rtn_m for r in records], dtype=np.float64).reshape(n, 3),
        np.array([r.toa_miss_s for r in records], dtype=np.float64),
        np.array([r.perigee_alt_m for r in records], dtype=np.float64),
        np.array([r.total_dv_m_s for r in records], dtype=np.float64),
        np.array([r.converged for r in records], dtype=np.bool_),
    )
    return EnsembleResult(
        master_seed=master_seed,
        nominal_perigee_alt_m=ctx.nominal.perigee_alt_m,
        nominal_toa_s=ctx.nominal.toa_s,
        records=tuple(records),
        stats=stats,
    )


def run_sweep(
    spec: SweepSpec,
    control: Controller,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    toa_window_s: float | None = None,
) -> SweepResult:
    """Run the deterministic A3 controllability grid against a fixed targeter (ADR 0007).

    Same physics path as :func:`run_ensemble` (shared :func:`_build_context` /
    :func:`_run_record` / nominal crossing), but the inputs are the deterministic
    :func:`~puffsat_sim.sweep.grid_inputs` (zero injection, swept Cd/Cr) rather than
    stochastic draws, and ``control`` is required — A3 maps the *required Δv*, so there is no
    open-loop variant.  ``toa_window_s`` arms the spurious-far-root gate (decision 3iii); the
    caller sizes it off the capstone's open-loop ToA dispersion.

    ``SweepResult.nominal`` is a dedicated factor-(1,1) reference run (zero coefficient error,
    zero injection) so the perigee/ToA overlays have a baseline even when the grid does not
    land a point exactly on nominal.
    """
    ctx = _build_context(orbital_config)
    records = tuple(_run_record(ctx, inputs, control, toa_window_s) for inputs in grid_inputs(spec))
    nominal_inputs = RunInputs(
        run_index=-1,
        dv_rtn_m_s=(0.0, 0.0, 0.0),
        cd_area_over_mass=spec.cd_area_over_mass,
        cr_area_over_mass=spec.cr_area_over_mass,
        f10p7=spec.f10p7,
        ap=spec.ap,
    )
    nominal_record = _run_record(ctx, nominal_inputs, control, toa_window_s)
    return SweepResult(spec=spec, records=records, nominal=nominal_record)


# C0 nominal truth (ADR 0011 decision 5): zero injection + nominal coefficients (perfect model),
# so x_true is the nominal apogee state and the only predict-vs-execute divergence is the nav error.
_NOMINAL_CD_AREA_OVER_MASS: float = 0.04
_NOMINAL_CR_AREA_OVER_MASS: float = 0.02
_NOMINAL_F10P7: float = 150.0
_NOMINAL_AP: float = 15.0


def _nav_nominal_inputs(run_index: int) -> RunInputs:
    """Zero-injection, nominal-coefficient truth for a C0 cell (ADR 0011 decision 5)."""
    return RunInputs(
        run_index=run_index,
        dv_rtn_m_s=(0.0, 0.0, 0.0),
        cd_area_over_mass=_NOMINAL_CD_AREA_OVER_MASS,
        cr_area_over_mass=_NOMINAL_CR_AREA_OVER_MASS,
        f10p7=_NOMINAL_F10P7,
        ap=_NOMINAL_AP,
    )


def run_nav_sweep(
    spec: NavSweepSpec,
    control: Controller,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> NavSweepResult:
    """Run the deterministic C0 navigation-error sensitivity sweep (ADR 0011).

    Same physics path as :func:`run_sweep` / :func:`run_ensemble` (shared
    :func:`_build_context` / :func:`_run_record` / nominal crossing), but the perturbation is a
    **predict-side** apogee-RTN nav-error offset, with zero injection and nominal coefficients
    (perfect model) — so the only predict-vs-execute divergence is the nav error.  Each cell
    runs the corrector from the perturbed estimate and executes against truth; the recorded
    ``miss_rtn_m`` is the residual ``−Φδ`` (uncontrollable at the apogee node) and
    ``total_dv_m_s`` the phantom correction the corrector burned chasing it.  The zero cell is
    the on-target reference (residual ~0).  ``control`` is required — the corrector is C0's
    subject, not an option.
    """
    ctx = _build_context(orbital_config)
    cells = nav_grid_offsets(spec)
    records = tuple(
        _run_record(
            ctx,
            _nav_nominal_inputs(cell.cell_index),
            control,
            nav_offset_rtn6=cell.offset_rtn6,
        )
        for cell in cells
    )
    return NavSweepResult(spec=spec, cells=cells, records=records)


def _c0_controller(predict: PredictFn, target: Target) -> ControlPlan:
    """C0 corrector: tol tight enough that the predict-null floor sits well below the residual −Φδ.

    LM damping + budget-scale step cap as in the A3 sweep (ADR 0007 decision 3); ``tol_m=0.01``
    so the corrector's 1-cm predict-null floor is far under even the few-metre smallest residuals,
    keeping the assembled Φ clean.
    """
    return solve_apogee_correction(
        predict, target, tol_m=0.01, lm=True, max_step_m_s=50.0, max_iter=15
    )


def nav_requirement_report(
    spec: NavSweepSpec | None = None,
    catch_radii_m: Sequence[float] = (5_000.0, 1_000.0, 100.0),
) -> str:
    """Run the C0 nav-error sweep and reduce it to the navigation-requirement report (ADR 0011)."""
    result = run_nav_sweep(spec or NavSweepSpec(), control=_c0_controller)
    req = summarize_nav_requirement(result)
    phi_lines = ["  Φ (3×6 apogee→crossing sensitivity), rows R/T/N, cols R/T/N-pos R/T/N-vel:"]
    phi_lines += ["    " + "  ".join(f"{v:+.3e}" for v in row) for row in req.phi]
    return format_nav_requirement(req, catch_radii_m) + "\n" + "\n".join(phi_lines)


def _truth_arc_to_apogee(
    cadence_hz: float,
    arc_duration_s: float,
    physics: PhysicsConfig,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> NDArray[np.float64]:
    """Full-force truth states at each measurement epoch of the C1 coast arc (ADR 0012).

    The arc *ends at the apogee deployment state*: the propagator runs backward
    ``arc_duration_s`` from apogee, then forward through the cadence epochs — so
    the final row is the apogee correction node where Σ is judged.  Rows are
    EME2000 ``[x y z vx vy vz]``, the arc start plus one row per epoch — the
    ``truth_states`` contract of :func:`puffsat_sim.nav_feasibility.validate_cell`.
    """
    apo_date, apo_pos, apo_vel = _apogee_state(orbital_config)
    frame = FramesFactory.getEME2000()
    orbit = CartesianOrbit(
        TimeStampedPVCoordinates(apo_date, apo_pos, apo_vel),
        frame,
        float(Constants.WGS84_EARTH_MU),
    )
    propagator = build_propagator_from_orbit(orbit, physics, _COAST_MAX_STEP_S)

    dt_s = 1.0 / cadence_hz
    n_epochs = max(1, round(arc_duration_s * cadence_hz))
    states = np.empty((n_epochs + 1, 6), dtype=np.float64)
    for k in range(n_epochs + 1):
        date = apo_date.shiftedBy(float(k * dt_s - n_epochs * dt_s))
        pv = propagator.propagate(date).getPVCoordinates()
        p, v = pv.getPosition(), pv.getVelocity()
        states[k] = (p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ())
    return states


def run_nav_feasibility(
    spec: NavFeasibilitySpec,
    phi: NDArray[np.float64],
    seed: int,
    validation_indices: tuple[int, ...] | None = None,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> tuple[NavFeasibilityResult, tuple[NavValidationOutcome, ...]]:
    """The C1 JVM seam: pure LinCov sweep + seeded UKF truth runs judged by NEES.

    Layer 1 (the envelope) is entirely pure — only the layer-2 validation touches
    the JVM, generating the full-force truth arc the filter's two-body+J2 model
    must stay honest against (the model gap q absorbs, ADR 0012 decisions 6/7).
    ``phi`` is C0's measured sensitivity; ``validation_indices`` defaults to the
    nominal cell plus the range-only cell (the two headline claims).  Truth arcs
    are cached per cadence — they depend on sampling, not on the cell's noises.
    """
    result = sweep_nav_feasibility(spec, phi)
    if validation_indices is None:
        range_only = tuple(
            o.cell.cell_index for o in result.outcomes if o.cell.doppler_sigma_m_s is None
        )
        validation_indices = (0,) + range_only

    physics = physics_from_inputs(_nav_nominal_inputs(0))
    arcs: dict[float, NDArray[np.float64]] = {}
    validations = []
    for index in validation_indices:
        cell = result.outcomes[index].cell
        if cell.cadence_hz not in arcs:
            arcs[cell.cadence_hz] = _truth_arc_to_apogee(
                cell.cadence_hz, spec.arc_duration_s, physics, orbital_config
            )
        validations.append(validate_cell(cell, spec, arcs[cell.cadence_hz], seed))
    return result, tuple(validations)


def nav_feasibility_report(
    spec: NavFeasibilitySpec | None = None,
    seed: int = 20260610,
    nav_spec: NavSweepSpec | None = None,
) -> str:
    """Run the full C1 report: measured Φ → LinCov envelope → NEES validation (ADR 0012).

    Φ is re-derived live from a minimal C0 sweep (``points_per_sign=1`` — C0
    measured the residual linear across the full swept range, so the smallest
    ± pair per axis already determines every column).
    """
    nav_result = run_nav_sweep(nav_spec or NavSweepSpec(points_per_sign=1), _c0_controller)
    requirement = summarize_nav_requirement(nav_result)
    result, validations = run_nav_feasibility(spec or NavFeasibilitySpec(), requirement.phi, seed)
    return format_nav_feasibility(result) + "\n" + format_nav_validation(validations)


def coeff_requirement_report(
    catch_radius_m: float = 5_000.0,
    prior_sigma_factor: float = 0.2,
    cut_points: int = 3,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Run the C2a coefficient-knowledge requirement (ADR 0013): measured cuts → tolerance vs prior.

    Thin glue over already-covered harnesses (the B2 no-integration-test precedent):
    Φ from a minimal C0 sweep, ``∂Δv/∂c`` from two 1D A3 cuts under the same LM-damped
    corrector, the verdict from the pure :mod:`puffsat_sim.coeff_requirement` chain.
    The coast for the analytic SRP cross-check is the apogee→perigee half period.
    """
    nav_result = run_nav_sweep(NavSweepSpec(points_per_sign=1), _c0_controller, orbital_config)
    phi = summarize_nav_requirement(nav_result).phi
    cd_cut = run_sweep(SweepSpec(cd_points=cut_points, cr_points=1), _c0_controller, orbital_config)
    cr_cut = run_sweep(SweepSpec(cd_points=1, cr_points=cut_points), _c0_controller, orbital_config)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    requirement = summarize_coeff_requirement(
        phi,
        cd_cut=cd_cut,
        cr_cut=cr_cut,
        catch_radius_m=catch_radius_m,
        prior_sigma_factor=prior_sigma_factor,
        coast_duration_s=keplerian_period(semi_major) / 2.0,
    )
    return format_coeff_requirement(requirement)


def format_summary(result: EnsembleResult) -> str:
    """Human-readable one-screen summary of an ensemble (the capstone / Rung A report)."""
    s = result.stats
    mean = s.miss_rtn_mean_m
    std = s.miss_rtn_std_m
    controlled = any(r.control_log for r in result.records)
    title = "Closed-loop dispersion ensemble" if controlled else "Open-loop dispersion capstone"
    lines = [
        f"{title} — N={s.n}, master_seed={result.master_seed}",
        f"  Nominal: perigee {result.nominal_perigee_alt_m / 1e3:.1f} km,"
        f" coast {result.nominal_toa_s / 3600:.2f} h",
        "  Interception miss vs nominal, RTN frame [m] (T = dr_p/dv_a lever):",
        f"    bias R/T/N = {mean[0]:+.1f} / {mean[1]:+.1f} / {mean[2]:+.1f}",
        f"    std  R/T/N = {std[0]:.1f} / {std[1]:.1f} / {std[2]:.1f}",
        f"  Time-of-arrival miss: {s.toa_miss_mean_s:+.2f} ± {s.toa_miss_std_s:.2f} s",
        f"  Perigee (diagnostic, low=good): {s.perigee_alt_mean_m / 1e3:.1f}"
        f" ± {s.perigee_alt_std_m / 1e3:.2f} km"
        f" [min {s.perigee_alt_min_m / 1e3:.1f}, max {s.perigee_alt_max_m / 1e3:.1f}]",
    ]
    if controlled:
        lines.append(
            f"  Correction Δv [m/s]: mean {s.total_dv_mean_m_s:.4f},"
            f" max {s.total_dv_max_m_s:.4f} (std {s.total_dv_std_m_s:.4f})"
        )
        lines.append(f"  Corrector converged: {s.converged_fraction * 100:.0f}% of {s.n} runs")
    return "\n".join(lines)


def main() -> None:
    # Smoke-sized ensemble (design doc §10.4: N=50 is a smoke test, not the result;
    # the resolved-tail controllability result needs N=10³–10⁴, a longer job).
    result = run_ensemble(DispersionSpec(), n=50, master_seed=20260608)
    print(format_summary(result))


if __name__ == "__main__":
    main()
