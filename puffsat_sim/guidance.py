"""Pure C3b terminal-guidance core (ADR 0014/0015, no JVM).

The dumb fixed terminal law of ADR 0014 decision 4: at each control step, predict the
crossing miss under no-further-thrust, command ``a = k·ZEM/t_go²`` capped at the
actuator, hold over the step (ZOH).  The closed loop that executes the commands lives
JVM-side in :mod:`puffsat_sim.runs.guidance`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.anti_drag import PEAK_SLEW_LIMIT_DEG_S, PEAK_THRUST_LIMIT_N, THRUST_FLOOR_N
from puffsat_sim.dispersion import Vec3
from puffsat_sim.estimation import two_body_j2_flow
from puffsat_sim.propellant import propellant_curve
from puffsat_sim.terminal import FeedforwardPlan

ZEM_GAIN: float = 3.0  # minimum-energy gain for the double integrator

# Plate-capture criterion (ADR 0015 decision 1): ≥99 % per-PuffSat capture on a 5 m
# plate ↔ σ_lateral ≤ 1.65 m (2D Rayleigh, R/σ = 3.03); ToA from the pulse cadence.
PLATE_RADIUS_M: float = 5.0
CAPTURE_SIGMA_MAX_M: float = 1.65
TOA_LIMIT_S: float = 0.010


@dataclass(frozen=True)
class TrackerGrade:
    """One terminal nav-noise grade for the σ_rel sweep (ADR 0015 decision 4).

    The angle grade is the target-side astrometric tracker: lateral position noise
    σ_θ·R on the two axes perpendicular to the line of sight, two-way-ranging noise
    along it.  ``sigma_theta_rad=None`` is the range-independent continuity grade
    (ADR 0014's code-differential point): isotropic ``sigma_range_m`` on all axes.
    """

    sigma_theta_rad: float | None
    sigma_range_m: float


def _perp_basis(unit: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """A deterministic orthonormal pair perpendicular to a unit vector."""
    seed = (0.0, 0.0, 1.0) if abs(unit[2]) < 0.9 else (0.0, 1.0, 0.0)
    e1 = np.cross(unit, seed)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(unit, e1)


def position_noise(grade: TrackerGrade, rng: np.random.Generator, los_m: Vec3) -> Vec3:
    """One draw of relative-position measurement noise for a line of sight to the target."""
    if grade.sigma_theta_rad is None:
        draw = rng.normal(0.0, grade.sigma_range_m, 3)
        return (float(draw[0]), float(draw[1]), float(draw[2]))
    los = np.asarray(los_m, dtype=np.float64)
    range_m = float(np.linalg.norm(los))
    unit = los / range_m
    e1, e2 = _perp_basis(unit)
    sigma_lat = grade.sigma_theta_rad * range_m
    noise = (
        rng.normal(0.0, grade.sigma_range_m) * unit
        + rng.normal(0.0, sigma_lat) * e1
        + rng.normal(0.0, sigma_lat) * e2
    )
    return (float(noise[0]), float(noise[1]), float(noise[2]))


# Tracker knowledge-error correlation time: the astrometric track is gyro-bridged and
# smoothed over ~tens of 1 Hz frames (ADR 0015 decision 3), so the injected state error
# is a Gauss–Markov process riding the σ_θ·R envelope — NOT white per-tick measurement
# noise, which would let the loop thrash on raw draws no real filter would pass through.
NAV_CORRELATION_TIME_S: float = 10.0


class NavNoiseProcess:
    """Correlated nav knowledge error for one run: unit Gauss–Markov × the grade envelope.

    Each call advances the 3-axis unit process by ``dt_s`` and scales it onto the
    line-of-sight frame: σ_θ·R on the two lateral axes, ranging σ along the LOS
    (isotropic σ for the constant continuity grade).
    """

    def __init__(
        self,
        grade: TrackerGrade,
        rng: np.random.Generator,
        tau_s: float = NAV_CORRELATION_TIME_S,
    ) -> None:
        self._grade = grade
        self._rng = rng
        self._tau_s = tau_s
        self._g: NDArray[np.float64] = rng.standard_normal(3)

    def sample(self, los_m: Vec3, dt_s: float) -> Vec3:
        rho = float(np.exp(-dt_s / self._tau_s))
        self._g = rho * self._g + np.sqrt(1.0 - rho**2) * self._rng.standard_normal(3)
        if self._grade.sigma_theta_rad is None:
            noise = self._grade.sigma_range_m * self._g
            return (float(noise[0]), float(noise[1]), float(noise[2]))
        los = np.asarray(los_m, dtype=np.float64)
        range_m = float(np.linalg.norm(los))
        unit = los / range_m
        e1, e2 = _perp_basis(unit)
        sigma_lat = self._grade.sigma_theta_rad * range_m
        noise = (
            self._grade.sigma_range_m * self._g[0] * unit
            + sigma_lat * self._g[1] * e1
            + sigma_lat * self._g[2] * e2
        )
        return (float(noise[0]), float(noise[1]), float(noise[2]))


def predicted_zem(state_eme: NDArray[np.float64], target_position_m: Vec3, t_go_s: float) -> Vec3:
    """The zero-effort miss toward the target under the onboard no-further-thrust model.

    The onboard predictor is two-body + J2 (the C1 flight-filter dynamics, ADR 0012):
    drag is omitted because the feedforward cancels it (C3a measured the uncompensated
    displacement at 8.5 cm), and re-prediction at every control step suppresses the
    remaining model bias (ADR 0014).
    """
    predicted = two_body_j2_flow(state_eme, t_go_s)[:3]
    return (
        float(target_position_m[0] - predicted[0]),
        float(target_position_m[1] - predicted[1]),
        float(target_position_m[2] - predicted[2]),
    )


def homing_floor_m(sigma_theta_rad: float, speed_m_s: float, a_max_m_s2: float) -> float:
    """The angle-noise homing floor ``σ_miss ≈ 2σ_θ²v²/a_max`` (ADR 0015).

    Knowledge improves as σ_θ·R while authority dies as ½·a_max·(R/v)²; the crossover
    range sets the smallest miss the loop can act on.
    """
    return 2.0 * sigma_theta_rad**2 * speed_m_s**2 / a_max_m_s2


def zem_acceleration(zem_m: Vec3, t_go_s: float, gain: float = ZEM_GAIN) -> Vec3:
    """The ZEM-law commanded acceleration ``k·ZEM/t_go²`` toward the target (uncapped).

    ``zem_m`` is the zero-effort miss *toward* the target (target minus predicted
    no-further-thrust position), so the command accelerates along it.  The actuator
    cap is applied by the caller on the *combined* command (ZEM + drag feedforward
    share the single thruster).
    """
    factor = gain / t_go_s**2
    return (factor * zem_m[0], factor * zem_m[1], factor * zem_m[2])


# Noise discipline for the ZEM command (measured on the double-integrator harness,
# C3b findings).  Closing the raw law on σ_θ·R-envelope knowledge wrecks the loop:
# phantom estimates at large R swing the commanded direction faster than the 1 °/s
# gimbal, so full-magnitude thrust fires along stale directions and rectifies
# zero-mean noise into ~150 m RMS of real error at 10 µrad (vs the 1.45 m closed-form
# floor).  Three declared constants fix what can be fixed:
#  - GATE_SIGMAS: stay silent unless the estimate is significant (|ZEM| > n·σ).  The
#    gate decides *whether* to act, never *how much* — acting on the excess
#    (soft-thresholding) parks the real error at the decaying floor n·σ_θ·v·t_go,
#    whose tracking saturates the actuator below t_go ≈ 2n·σ_θ·v/a_max and strands
#    n²× the floor (measured 12.4 m at n = 3, 10 µrad).
#  - TRACK_WINDOW_S: below this t_go the gate is held open — a phantom firing costs
#    a_max·t_go·dt of real error, cheap late, and the endgame must track the estimate
#    continuously or the chip-down stalls one gate-width above the floor.
#  - MAX_FIRING_LAG_DEG: hold fire while the gimbal hasn't caught up to the command;
#    burning mid-slew is the rectifier itself.
# The measured result of the package: lateral RMS ≈ 3× the closed-form homing floor
# (the residual price of correlated knowledge error on a slew-limited single engine).
GATE_SIGMAS: float = 3.0
TRACK_WINDOW_S: float = 35.0
MAX_FIRING_LAG_DEG: float = 45.0


def significant_zem(zem_est_m: Vec3, knowledge_sigma_m: float, t_go_s: float) -> Vec3:
    """The estimated ZEM if it is worth acting on, else zero (see the constants above)."""
    if t_go_s <= TRACK_WINDOW_S:
        return zem_est_m
    mag = (zem_est_m[0] ** 2 + zem_est_m[1] ** 2 + zem_est_m[2] ** 2) ** 0.5
    if mag <= GATE_SIGMAS * knowledge_sigma_m:
        return (0.0, 0.0, 0.0)
    return zem_est_m


@dataclass(frozen=True)
class TickCommand:
    """One control tick realized on the actuator: throttle, gimbal state, fire decision."""

    thrust_n: float
    attitude_dir: Vec3 | None
    fire: bool
    saturated: bool


def terminal_tick(
    zem_est_m: Vec3,
    knowledge_sigma_m: float,
    t_go_s: float,
    feedforward_m_s2: Vec3,
    attitude_dir: Vec3 | None,
    control_period_s: float,
    mass_kg: float,
    gain: float = ZEM_GAIN,
) -> TickCommand:
    """One tick of the C3b terminal law, from ZEM estimate to actuator command.

    Gate the estimate (:func:`significant_zem`), add the anti-drag feedforward (the
    single thruster serves both), cap on the actuator, slew the gimbal toward the
    command at the ADR 0004 rate, and fire only when the demand clears the
    proportional floor *and* the gimbal has caught up to within
    ``MAX_FIRING_LAG_DEG`` — burning mid-slew converts knowledge noise into real
    error.  Pure, so the JVM loop and the unit harness share one law.
    """
    zem = significant_zem(zem_est_m, knowledge_sigma_m, t_go_s)
    a_zem = zem_acceleration(zem, t_go_s, gain)
    thrust_n, commanded_dir, saturated = thrust_command(
        (
            a_zem[0] + feedforward_m_s2[0],
            a_zem[1] + feedforward_m_s2[1],
            a_zem[2] + feedforward_m_s2[2],
        ),
        mass_kg=mass_kg,
        max_thrust_n=PEAK_THRUST_LIMIT_N,
    )
    if thrust_n <= 0.0:
        return TickCommand(0.0, attitude_dir, False, False)
    attitude = (
        commanded_dir
        if attitude_dir is None
        else slew_limited_direction(
            attitude_dir, commanded_dir, PEAK_SLEW_LIMIT_DEG_S * control_period_s
        )
    )
    lag_cos = (
        attitude[0] * commanded_dir[0]
        + attitude[1] * commanded_dir[1]
        + attitude[2] * commanded_dir[2]
    )
    fire = thrust_n >= THRUST_FLOOR_N and lag_cos >= np.cos(np.radians(MAX_FIRING_LAG_DEG))
    return TickCommand(thrust_n, attitude, fire, saturated)


@dataclass(frozen=True)
class PlateMiss:
    """One arrival in the plate frame: 2D lateral miss ⊥ v_rel + time-of-arrival error.

    At closest approach the miss vector is perpendicular to the relative velocity by
    definition (ADR 0015); the along-track axis is *time*, not position.  Negative
    ``toa_error_s`` means early.
    """

    lateral_m: tuple[float, float]
    toa_error_s: float

    @property
    def lateral_norm_m(self) -> float:
        return float(np.hypot(self.lateral_m[0], self.lateral_m[1]))


def plate_frame_miss(
    position_m: Vec3,
    velocity_m_s: Vec3,
    toa_s: float,
    target_position_m: Vec3,
    target_toa_s: float,
) -> PlateMiss:
    """Decompose a crossing state vs the target point into the plate frame.

    The crossing is read at the 200 km altitude event, not at closest approach to the
    target, so the along-velocity component of the position offset is converted to the
    coasted closest-approach time: ``dt = −d·v/|v|²`` leaves ``d + v·dt ⊥ v``.
    """
    d = np.asarray(position_m, dtype=np.float64) - np.asarray(target_position_m, dtype=np.float64)
    v = np.asarray(velocity_m_s, dtype=np.float64)
    speed_sq = float(v @ v)
    dt = -float(d @ v) / speed_sq
    lateral = d + v * dt

    unit = v / np.sqrt(speed_sq)
    e1, e2 = _perp_basis(unit)
    return PlateMiss(
        lateral_m=(float(lateral @ e1), float(lateral @ e2)),
        toa_error_s=(toa_s + dt) - target_toa_s,
    )


def capture_fraction(
    misses: Sequence[PlateMiss],
    plate_radius_m: float = PLATE_RADIUS_M,
    toa_limit_s: float = TOA_LIMIT_S,
) -> float:
    """The fraction of arrivals captured: lateral within the plate AND ToA within the window."""
    captured = sum(
        1
        for miss in misses
        if miss.lateral_norm_m <= plate_radius_m and abs(miss.toa_error_s) <= toa_limit_s
    )
    return captured / len(misses)


def capture_curve(
    misses: Sequence[PlateMiss], radii_m: Sequence[float]
) -> tuple[tuple[float, float], ...]:
    """Capture fraction vs plate radius (the ADR 0015 economics curve; Rung D's P(capture))."""
    return tuple((radius, capture_fraction(misses, plate_radius_m=radius)) for radius in radii_m)


def slew_limited_direction(previous: Vec3, commanded: Vec3, max_angle_deg: float) -> Vec3:
    """Rotate the executed thrust direction toward the command, capped at the slew budget.

    The ADR 0004 direction loop tracks at ≤1 °/s; commands turning faster than the
    budget since the last firing are executed lagging, and the lag's effect lands in
    the measured miss instead of an unphysical command history.
    """
    prev = np.asarray(previous, dtype=np.float64)
    cmd = np.asarray(commanded, dtype=np.float64)
    cos_angle = float(np.clip(prev @ cmd, -1.0, 1.0))
    angle = float(np.arccos(cos_angle))
    max_angle = np.radians(max_angle_deg)
    if angle <= max_angle:
        return commanded
    axis = np.cross(prev, cmd)
    norm = float(np.linalg.norm(axis))
    if norm < 1e-12:
        axis = np.cross(prev, _perp_basis(prev)[0])
        norm = float(np.linalg.norm(axis))
    axis /= norm
    rotated = prev * np.cos(max_angle) + np.cross(axis, prev) * np.sin(max_angle)
    return (float(rotated[0]), float(rotated[1]), float(rotated[2]))


@dataclass(frozen=True)
class GuidanceRun:
    """One closed-loop terminal descent: the plate-frame arrival + the executed burn."""

    miss: PlateMiss
    plan: FeedforwardPlan
    saturated_fraction: float
    perigee_alt_m: float = 0.0  # the crossing perigee (Rung-D diagnostic, low = good)


@dataclass(frozen=True)
class GuidanceCell:
    """One sweep cell: repeated seeded runs at a single setting, with the summary reads.

    ``axis_value`` is the swept quantity in its own unit (entry offset [m], σ_θ [rad],
    cadence [Hz], drag factor); ``None`` for off-axis cells like the constant grade.
    """

    label: str
    runs: tuple[GuidanceRun, ...]
    axis_value: float | None = None

    @property
    def rms_lateral_m(self) -> float:
        return float(np.sqrt(np.mean([run.miss.lateral_norm_m**2 for run in self.runs])))

    @property
    def max_lateral_m(self) -> float:
        return max(run.miss.lateral_norm_m for run in self.runs)

    @property
    def rms_toa_s(self) -> float:
        return float(np.sqrt(np.mean([run.miss.toa_error_s**2 for run in self.runs])))

    @property
    def capture(self) -> float:
        return capture_fraction(tuple(run.miss for run in self.runs))

    @property
    def max_dv_m_s(self) -> float:
        return max(run.plan.dv_m_s for run in self.runs)


def measured_catch_radius_m(entry_cells: Sequence[GuidanceCell]) -> float | None:
    """The largest swept entry offset whose worst residual is still capture-grade.

    This is the measured knee of the residual-vs-entry curve (ADR 0014 decision 1):
    the radius the thrust-limited terminal burn actually buys, judged at the ADR 0015
    σ_lateral ≤ 1.65 m capture requirement.  ``None`` means no swept offset qualified.
    """
    qualifying = [
        cell.axis_value
        for cell in entry_cells
        if cell.axis_value is not None and cell.max_lateral_m <= CAPTURE_SIGMA_MAX_M
    ]
    return max(qualifying) if qualifying else None


@dataclass(frozen=True)
class GuidanceSweepSpec:
    """The C3b one-axis-at-a-time sweep grid (ADR 0014 decision 6, ADR 0015 decision 4).

    Entry offsets are noiseless single runs (the residual-vs-entry curve whose knee is
    the measured catch radius); tracker grades are seeded repeats at zero offset (the
    nav floor); cadence cells repeat at the nominal grade; drag cells are noiseless
    truth-model mismatches against the nominal-planned feedforward.
    """

    entry_offsets_m: tuple[float, ...] = (0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0)
    # Tracker grades are architectures, not error bars (ADR 0015): 2/5 µrad are aspirational
    # and design-target target-side optical astrometry, 10 µrad the optical requirement ceiling
    # (the validated nominal), 50 µrad ground-station instrumentation radar — kept as the
    # measured refutation of ground terminal tracking ("funnel-entry class, not terminal class").
    sigma_thetas_rad: tuple[float, ...] = (2e-6, 5e-6, 10e-6, 50e-6)
    constant_sigma_m: float | None = 1.0
    nominal_sigma_theta_rad: float = 10e-6
    sigma_range_m: float = 1.0
    cadences_hz: tuple[float, ...] = (0.1, 10.0)
    drag_factors: tuple[float, ...] = (0.5, 2.0)
    storm_f10p7_ap: tuple[float, float] | None = (250.0, 100.0)
    n_noise_runs: int = 8
    n_cadence_runs: int = 4
    master_seed: int = 20260612
    control_period_s: float = 1.0
    gain: float = ZEM_GAIN


@dataclass(frozen=True)
class TerminalGuidanceFinding:
    """The C3b sweep measurement set (pure container for the JVM numbers)."""

    entry_cells: tuple[GuidanceCell, ...]
    grade_cells: tuple[GuidanceCell, ...]
    cadence_cells: tuple[GuidanceCell, ...]
    drag_cells: tuple[GuidanceCell, ...]
    cadence_hz: float
    gain: float
    a_max_m_s2: float
    speed_m_s: float


def _cell_line(cell: GuidanceCell, floor_m: float | None = None) -> str:
    floor = f", floor {floor_m:.2f} m" if floor_m is not None else ""
    sat = max(run.saturated_fraction for run in cell.runs)
    return (
        f"    {cell.label:>16} → RMS lateral {cell.rms_lateral_m:.3f} m"
        f" (max {cell.max_lateral_m:.3f}), RMS ToA {cell.rms_toa_s * 1e3:.2f} ms,"
        f" capture {cell.capture * 100:.0f}%{floor},"
        f" Δv {cell.max_dv_m_s:.3f} m/s [sat {sat * 100:.0f}%]"
    )


# The executed gimbal history is rail-limited by construction (slew_limited_direction
# caps every step at the budget), but recomputing the step angle through acos(dot)
# carries ~1e-12 of round-off — riding the rail must not read as a gate violation.
_SLEW_GATE_ROUNDOFF: float = 1e-9


def format_terminal_guidance(finding: TerminalGuidanceFinding) -> str:
    """One-screen C3b report: entry curve + catch radius, grade floors, cadence/drag, gates."""
    worst = max(
        (run for cell in finding.entry_cells + finding.grade_cells for run in cell.runs),
        key=lambda run: run.plan.dv_m_s,
    )
    peak_thrust = max(
        run.plan.peak_thrust_n
        for cell in finding.entry_cells + finding.grade_cells
        for run in cell.runs
    )
    peak_slew = max(
        run.plan.peak_slew_rate_deg_s
        for cell in finding.entry_cells + finding.grade_cells
        for run in cell.runs
    )
    radius = measured_catch_radius_m(finding.entry_cells)
    radius_text = f"{radius:.0f} m" if radius is not None else "below the smallest swept offset"
    points = propellant_curve(worst.plan.dv_m_s)
    curve = ", ".join(f"{p.fraction * 100:.4f}% @Isp{p.isp_s:.0f}" for p in points)
    thrust_verdict = "PASS" if peak_thrust <= PEAK_THRUST_LIMIT_N else "FAIL"
    slew_verdict = (
        "PASS" if peak_slew <= PEAK_SLEW_LIMIT_DEG_S * (1.0 + _SLEW_GATE_ROUNDOFF) else "FAIL"
    )

    lines = [
        "C3b terminal guidance — closed ZEM loop (ADR 0014/0015)",
        f"  Law: a = {finding.gain:.0f}·ZEM/t_go² (ZOH @ {finding.cadence_hz:g} Hz),"
        f" a_max {finding.a_max_m_s2:.3f} m/s², v_rel {finding.speed_m_s:.0f} m/s",
        f"  Entry-offset curve (σ ≤ {CAPTURE_SIGMA_MAX_M} m is capture-grade):",
    ]
    lines.extend(_cell_line(cell) for cell in finding.entry_cells)
    lines.append(f"    → measured catch radius: {radius_text}")
    lines.append("  Tracker grades (σ_rel(R) = σ_θ·R; closed-form floor 2σ_θ²v²/a_max):")
    lines.extend(
        _cell_line(
            cell,
            floor_m=(
                homing_floor_m(cell.axis_value, finding.speed_m_s, finding.a_max_m_s2)
                if cell.axis_value is not None
                else None
            ),
        )
        for cell in finding.grade_cells
    )
    lines.append("  Control cadence:")
    lines.extend(_cell_line(cell) for cell in finding.cadence_cells)
    lines.append("  Dispersed drag (truth ≠ feedforward):")
    lines.extend(_cell_line(cell) for cell in finding.drag_cells)
    lines.append(
        f"  Gates (ADR 0004): peak thrust {peak_thrust * 1e3:.1f} mN"
        f" vs {PEAK_THRUST_LIMIT_N * 1e3:.0f} mN — {thrust_verdict};"
        f" peak slew {peak_slew:.3f} °/s vs {PEAK_SLEW_LIMIT_DEG_S:.1f} °/s — {slew_verdict}"
    )
    lines.append(f"  Propellant (worst run Δv {worst.plan.dv_m_s:.4f} m/s): {curve}")
    return "\n".join(lines)


def thrust_command(
    accel_m_s2: Vec3, mass_kg: float, max_thrust_n: float, floor_n: float = 0.0
) -> tuple[float, Vec3, bool]:
    """Realize a commanded acceleration on the single actuator: ``(thrust, direction, saturated)``.

    The cap binds on the *combined* command (ZEM + drag feedforward share one thruster),
    so callers sum accelerations before calling.  Demands below ``floor_n`` (the ADR 0004
    ~5 mN proportional floor) are engine-off — firing them would be direction chatter the
    real thruster cannot deliver.  A zero command keeps a placeholder direction (the C3a
    zero-drag convention).
    """
    mag = (accel_m_s2[0] ** 2 + accel_m_s2[1] ** 2 + accel_m_s2[2] ** 2) ** 0.5
    if mag == 0.0 or mass_kg * mag < floor_n:
        return 0.0, (1.0, 0.0, 0.0), False
    direction = (accel_m_s2[0] / mag, accel_m_s2[1] / mag, accel_m_s2[2] / mag)
    saturated = mass_kg * mag > max_thrust_n
    return min(mass_kg * mag, max_thrust_n), direction, saturated
