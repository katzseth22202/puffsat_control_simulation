"""C3b closed ZEM terminal loop — the JVM run for :mod:`puffsat_sim.guidance` (ADR 0014/0015).

The loop closes tick by tick on the C3a fixed-step terminal machinery: from the 800 km
hand-off, each control step reads the truth state, injects the tracker-grade position
noise (σ_rel(R) = σ_θ·R), predicts the zero-effort miss with the onboard two-body+J2
model, and holds the capped ZEM + drag-feedforward command over the step as a finite
maneuver segment.  The aim point is the drag-free nominal crossing (the planned target;
consistent with the executed-feedforward nominal to the C3a 2 mm residual).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.attitudes import FrameAlignedProvider
from org.orekit.forces.maneuvers import ConstantThrustManeuver
from org.orekit.frames import FramesFactory
from org.orekit.orbits import CartesianOrbit
from org.orekit.utils import Constants, TimeStampedPVCoordinates

from org.hipparchus.geometry.euclidean.threed import Vector3D

from puffsat_sim import mission, presets
from puffsat_sim.anti_drag import PEAK_THRUST_LIMIT_N
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.descent import (
    COAST_MAX_STEP_S,
    coast_to_handoff,
    earth_model,
    propagate_to_interception,
    to_absolute_date,
    vec3,
)
from puffsat_sim.dispersion import Vec3, rtn_basis, rtn_to_cartesian
from puffsat_sim.guidance import (
    ZEM_GAIN,
    GuidanceCell,
    GuidanceRun,
    GuidanceSweepSpec,
    NavNoiseProcess,
    TerminalGuidanceFinding,
    TrackerGrade,
    format_terminal_guidance,
    plate_frame_miss,
    predicted_zem,
    terminal_tick,
)
from puffsat_sim.montecarlo import (
    BURN_ISP_SENTINEL_S,
    PROPAGATOR_MASS_KG,
    nominal_inputs,
    physics_from_inputs,
)
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import (
    build_fixed_step_propagator_from_orbit,
    build_propagator,
    build_propagator_from_orbit,
)
from puffsat_sim.runs.anti_drag import PUFFSAT_WET_MASS_KG, sample_drag_window
from puffsat_sim.runs.terminal import TERMINAL_FIXED_STEP_S, physics_without_drag
from puffsat_sim.terminal import FeedforwardPlan, ThrustCommand, executed_plan, plan_feedforward


@dataclass(frozen=True)
class GuidanceContext:
    """Per-sweep constants: the nominal hand-off, the aim point, and the drag profile."""

    frame: Any
    mu: float
    earth: Any
    epoch: Any
    period: float
    physics: PhysicsConfig
    handoff: Any
    target_position_m: Vec3
    target_toa_s: float
    target_speed_m_s: float
    ff_times_s: tuple[float, ...]
    ff_accels_m_s2: tuple[Vec3, ...]


def build_guidance_context(
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> GuidanceContext:
    """Coast to the 800 km hand-off once and fix the aim point + nominal drag profile."""
    physics = presets.full_force()
    earth = earth_model()
    epoch = to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)

    apogee_orbit = (
        build_propagator(orbital_config, physics, COAST_MAX_STEP_S).getInitialState().getOrbit()
    )
    coast = build_propagator_from_orbit(apogee_orbit, physics, COAST_MAX_STEP_S)
    handoff = coast_to_handoff(coast, epoch, period, earth)

    target = propagate_to_interception(
        build_fixed_step_propagator_from_orbit(
            handoff.getOrbit(), physics_without_drag(physics), TERMINAL_FIXED_STEP_S
        ),
        epoch,
        period,
        earth,
    )

    fixed = build_fixed_step_propagator_from_orbit(
        handoff.getOrbit(), physics, TERMINAL_FIXED_STEP_S
    )
    generator = fixed.getEphemerisGenerator()
    unburned = propagate_to_interception(fixed, epoch, period, earth)
    times, accels = sample_drag_window(
        generator.getGeneratedEphemeris(),
        handoff.getDate(),
        epoch.shiftedBy(unburned.toa_s),
        physics,
        earth,
        epoch,
        TERMINAL_FIXED_STEP_S,
    )

    return GuidanceContext(
        frame=FramesFactory.getEME2000(),
        mu=float(Constants.WGS84_EARTH_MU),
        earth=earth,
        epoch=epoch,
        period=period,
        physics=physics,
        handoff=handoff,
        target_position_m=target.position_m,
        target_toa_s=target.toa_s,
        target_speed_m_s=math.hypot(*target.velocity_m_s),
        ff_times_s=tuple(times),
        ff_accels_m_s2=tuple(accels),
    )


def _feedforward_accel(plan: FeedforwardPlan, starts_s: np.ndarray, t_s: float) -> Vec3:
    """The anti-drag feedforward acceleration commanded at time ``t_s`` (zero outside)."""
    idx = int(np.searchsorted(starts_s, t_s, side="right")) - 1
    if idx < 0:
        return (0.0, 0.0, 0.0)
    cmd = plan.commands[idx]
    if t_s >= cmd.start_s + cmd.duration_s:
        return (0.0, 0.0, 0.0)
    accel = cmd.thrust_n / plan.mass_kg
    return (accel * cmd.direction[0], accel * cmd.direction[1], accel * cmd.direction[2])


def run_guidance(
    ctx: GuidanceContext,
    entry_offset_m: float = 0.0,
    grade: TrackerGrade | None = None,
    control_period_s: float = 1.0,
    rng: np.random.Generator | None = None,
    truth_physics: PhysicsConfig | None = None,
    gain: float = ZEM_GAIN,
) -> GuidanceRun:
    """Fly one closed-loop terminal descent and read the arrival in the plate frame.

    ``entry_offset_m`` displaces the hand-off state along the RTN normal axis (pure
    lateral, exactly ⊥ v) — the funnel-entry error the loop must null.  ``grade`` is the
    injected nav noise (``None`` = perfect knowledge); ``truth_physics`` diverges the
    truth drag from the nominal-planned feedforward (the dispersed-drag cells).
    """
    physics = truth_physics if truth_physics is not None else ctx.physics
    step = min(TERMINAL_FIXED_STEP_S, control_period_s)

    pv = ctx.handoff.getPVCoordinates()
    basis = rtn_basis(vec3(pv.getPosition()), vec3(pv.getVelocity()))
    off = rtn_to_cartesian((0.0, 0.0, entry_offset_m), basis)
    orbit = CartesianOrbit(
        TimeStampedPVCoordinates(
            ctx.handoff.getDate(),
            pv.getPosition().add(Vector3D(off[0], off[1], off[2])),
            pv.getVelocity(),
        ),
        ctx.frame,
        ctx.mu,
    )
    prop = build_fixed_step_propagator_from_orbit(orbit, physics, step)

    ff_plan = plan_feedforward(
        ctx.ff_times_s,
        ctx.ff_accels_m_s2,
        mass_kg=PUFFSAT_WET_MASS_KG,
        control_period_s=control_period_s,
    )
    ff_starts = np.array([cmd.start_s for cmd in ff_plan.commands], dtype=np.float64)
    attitude = FrameAlignedProvider(ctx.frame)

    noise = NavNoiseProcess(grade, rng) if grade is not None and rng is not None else None
    commands: list[ThrustCommand] = []
    ticks = 0
    saturated_ticks = 0
    attitude_dir: Vec3 | None = None
    state = prop.getInitialState()
    t = float(state.getDate().durationFrom(ctx.epoch))
    while ctx.target_toa_s - t > control_period_s:
        t_go = ctx.target_toa_s - t
        state_pv = state.getPVCoordinates()
        pos = vec3(state_pv.getPosition())
        vel = vec3(state_pv.getVelocity())
        measured = pos
        sigma_m = 0.0
        if noise is not None and grade is not None:
            los = (
                ctx.target_position_m[0] - pos[0],
                ctx.target_position_m[1] - pos[1],
                ctx.target_position_m[2] - pos[2],
            )
            err = noise.sample(los, dt_s=control_period_s)
            measured = (pos[0] + err[0], pos[1] + err[1], pos[2] + err[2])
            sigma_m = (
                grade.sigma_theta_rad * math.hypot(*los)
                if grade.sigma_theta_rad is not None
                else grade.sigma_range_m
            )
        onboard = np.array([*measured, *vel], dtype=np.float64)
        zem = predicted_zem(onboard, ctx.target_position_m, t_go)
        tick = terminal_tick(
            zem,
            knowledge_sigma_m=sigma_m,
            t_go_s=t_go,
            feedforward_m_s2=_feedforward_accel(ff_plan, ff_starts, t),
            attitude_dir=attitude_dir,
            control_period_s=control_period_s,
            mass_kg=PUFFSAT_WET_MASS_KG,
            gain=gain,
        )
        attitude_dir = tick.attitude_dir
        ticks += 1
        saturated_ticks += int(tick.saturated)
        if tick.fire and attitude_dir is not None:
            prop.addForceModel(
                ConstantThrustManeuver(
                    state.getDate(),
                    control_period_s,
                    tick.thrust_n * PROPAGATOR_MASS_KG / PUFFSAT_WET_MASS_KG,
                    BURN_ISP_SENTINEL_S,
                    attitude,
                    Vector3D(attitude_dir[0], attitude_dir[1], attitude_dir[2]),
                )
            )
            commands.append(
                ThrustCommand(
                    start_s=t,
                    duration_s=control_period_s,
                    thrust_n=tick.thrust_n,
                    direction=attitude_dir,
                )
            )
        state = prop.propagate(state.getDate().shiftedBy(control_period_s))
        t = float(state.getDate().durationFrom(ctx.epoch))

    crossing = propagate_to_interception(prop, ctx.epoch, ctx.period, ctx.earth)
    return GuidanceRun(
        miss=plate_frame_miss(
            crossing.position_m,
            crossing.velocity_m_s,
            crossing.toa_s,
            ctx.target_position_m,
            ctx.target_toa_s,
        ),
        plan=executed_plan(tuple(commands), PUFFSAT_WET_MASS_KG, saturated=saturated_ticks > 0),
        saturated_fraction=saturated_ticks / ticks if ticks else 0.0,
    )


def _noise_cell(
    ctx: GuidanceContext,
    spec: GuidanceSweepSpec,
    label: str,
    grade: TrackerGrade,
    cell_index: int,
    n_runs: int,
    control_period_s: float,
    axis_value: float | None,
) -> GuidanceCell:
    runs = tuple(
        run_guidance(
            ctx,
            grade=grade,
            control_period_s=control_period_s,
            rng=np.random.default_rng((spec.master_seed, cell_index, i)),
            gain=spec.gain,
        )
        for i in range(n_runs)
    )
    return GuidanceCell(label=label, runs=runs, axis_value=axis_value)


def run_terminal_guidance(
    spec: GuidanceSweepSpec | None = None,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> TerminalGuidanceFinding:
    """Run the C3b one-axis-at-a-time sweep over one shared context (ADR 0014 decision 6)."""
    spec = spec if spec is not None else GuidanceSweepSpec()
    ctx = build_guidance_context(orbital_config)

    entry_cells = tuple(
        GuidanceCell(
            label=f"{offset:.0f} m",
            runs=(
                run_guidance(
                    ctx,
                    entry_offset_m=offset,
                    control_period_s=spec.control_period_s,
                    gain=spec.gain,
                ),
            ),
            axis_value=offset,
        )
        for offset in spec.entry_offsets_m
    )

    grade_cells = tuple(
        _noise_cell(
            ctx,
            spec,
            label=f"σ_θ = {sigma * 1e6:.0f} µrad",
            grade=TrackerGrade(sigma_theta_rad=sigma, sigma_range_m=spec.sigma_range_m),
            cell_index=i,
            n_runs=spec.n_noise_runs,
            control_period_s=spec.control_period_s,
            axis_value=sigma,
        )
        for i, sigma in enumerate(spec.sigma_thetas_rad)
    )
    if spec.constant_sigma_m is not None:
        grade_cells += (
            _noise_cell(
                ctx,
                spec,
                label=f"σ_rel = {spec.constant_sigma_m:.0f} m const",
                grade=TrackerGrade(sigma_theta_rad=None, sigma_range_m=spec.constant_sigma_m),
                cell_index=len(spec.sigma_thetas_rad),
                n_runs=spec.n_noise_runs,
                control_period_s=spec.control_period_s,
                axis_value=None,
            ),
        )

    nominal_grade = TrackerGrade(
        sigma_theta_rad=spec.nominal_sigma_theta_rad, sigma_range_m=spec.sigma_range_m
    )
    cadence_cells = tuple(
        _noise_cell(
            ctx,
            spec,
            label=f"{hz:g} Hz",
            grade=nominal_grade,
            cell_index=100 + i,
            n_runs=spec.n_cadence_runs,
            control_period_s=1.0 / hz,
            axis_value=hz,
        )
        for i, hz in enumerate(spec.cadences_hz)
    )

    drag_settings: list[tuple[str, float, float, float]] = [
        (f"Cd ×{factor:g}", factor, 0.0, 0.0) for factor in spec.drag_factors
    ]
    if spec.storm_f10p7_ap is not None:
        drag_settings.append(("storm F10.7/Ap", 1.0, *spec.storm_f10p7_ap))
    drag_cells = []
    base = nominal_inputs(0)
    for label, factor, f10p7, ap in drag_settings:
        inputs = replace(
            base,
            cd_area_over_mass=base.cd_area_over_mass * factor,
            f10p7=f10p7 or base.f10p7,
            ap=ap or base.ap,
        )
        drag_cells.append(
            GuidanceCell(
                label=label,
                runs=(
                    run_guidance(
                        ctx,
                        control_period_s=spec.control_period_s,
                        truth_physics=physics_from_inputs(inputs),
                        gain=spec.gain,
                    ),
                ),
                axis_value=factor,
            )
        )

    return TerminalGuidanceFinding(
        entry_cells=entry_cells,
        grade_cells=grade_cells,
        cadence_cells=cadence_cells,
        drag_cells=tuple(drag_cells),
        cadence_hz=1.0 / spec.control_period_s,
        gain=spec.gain,
        a_max_m_s2=PEAK_THRUST_LIMIT_N / PUFFSAT_WET_MASS_KG,
        speed_m_s=ctx.target_speed_m_s,
    )


def terminal_guidance_report(
    spec: GuidanceSweepSpec | None = None,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Run the C3b sweep and format the one-screen report."""
    return format_terminal_guidance(run_terminal_guidance(spec, orbital_config))
