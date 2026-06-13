"""Pure C4 control-loop latency core (ADR 0014, design doc §16.8, no JVM).

Dead-time in the loop (measure → coordinator compute → uplink → actuate) erodes stability
margin, most in the fast terminal phase.  §16.8's thesis and this module's two complementary
checks:

* **Deliverable — the analytic phase budget.**  What decides whether dead-time bites is the
  loop **bandwidth**, not the sample rate: dead-time costs phase margin as ``ω_c·τ``, and the
  drag-rejection bandwidth is ~1 Hz even though the inner loop *samples* at 100 Hz, so the
  tens-of-ms terminal budget costs single-digit degrees.  And the two latency sources live on
  **different loops**: the large comms round-trip rides the slow midcourse loop (discrete
  replan — no continuous phase loop to destabilize), while the fast terminal loop carries only
  the small onboard sensor/compute/valve delay.  This — the per-loop budget plus the ``ω_c·τ``
  check — is the back-of-envelope that settles §16.8.

* **Confirmation — a noiseless τ-sweep on the terminal homing loop.**  The dead-time mechanism
  is a buffer holding the loop's nav fix stale for τ (equivalent loop ``e^{-sτ}``); the sweep
  flies the C3b ZEM law on the double integrator (dead-time is a loop-transfer effect, so the
  Orekit physics add nothing the stability question needs).  It is run **noiseless on purpose**:
  the terminal nav floor under tracker noise is already C3b's measured finding (~1.07 m at the
  10 µrad grade on the full Orekit loop), so superimposing it here would only obscure the pure
  dead-time signal — and the *combined* entry-offset × tracker-noise stress at the dispersion
  tail is a Rung-D full-Monte-Carlo question by construction (C3b measured catch radius and nav
  floor on **separate** axes).  The sweep therefore reads **relative degradation** vs the
  zero-delay run, not absolute capture: it confirms the homing miss is flat through dead-time
  far past the budget, degrading only at τ on the order of a whole control period.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass

import numpy as np

from puffsat_sim.anti_drag import PEAK_THRUST_LIMIT_N
from puffsat_sim.dispersion import Vec3
from puffsat_sim.guidance import (
    ZEM_GAIN,
    NavNoiseProcess,
    TrackerGrade,
    terminal_tick,
)

SPEED_OF_LIGHT_M_S: float = 299_792_458.0


def comms_one_way_s(range_m: float) -> float:
    """One-way light-time over a link of length ``range_m`` (the irreducible comms dead-time)."""
    return range_m / SPEED_OF_LIGHT_M_S


@dataclass(frozen=True)
class LatencySource:
    """One named contributor to a loop's dead-time (seconds)."""

    label: str
    tau_s: float


@dataclass(frozen=True)
class ControlLoop:
    """A control loop's dead-time budget and the phase it erodes.

    ``bandwidth_hz`` is the closed-loop crossover (``ω_c = 2π·f``); ``None`` marks a
    *discrete* loop (the midcourse replan fires impulsive corrections, not a continuous
    feedback loop, so there is no phase margin to erode — its latency only shifts when a
    correction is computed against an hours-scale maneuver timeline).
    """

    name: str
    bandwidth_hz: float | None
    sources: tuple[LatencySource, ...]

    @property
    def tau_s(self) -> float:
        return sum(s.tau_s for s in self.sources)

    @property
    def phase_margin_loss_deg(self) -> float | None:
        """Phase eroded by the dead-time: ``ω_c·τ`` in degrees (``None`` for a discrete loop)."""
        if self.bandwidth_hz is None:
            return None
        return math.degrees(2.0 * math.pi * self.bandwidth_hz * self.tau_s)


# The measured back-of-envelope budget (§16.8).  The drag-rejection bandwidth is ~1 Hz
# (drag varies over the ~3-min descent) — the inner loop samples the *magnitude* at 100 Hz
# but the bandwidth that sets stability is 1 Hz.
DRAG_REJECTION_BANDWIDTH_HZ: float = 1.0
# Terminal closing baseline (ADR 0014/0015: ≤100 km) and a representative coordinator-node
# range (design doc §16.8: distance/c ≈ 6.7 ms one-way at 2000 km).
TERMINAL_CROSSLINK_RANGE_M: float = 100_000.0
COORDINATOR_NODE_RANGE_M: float = 2_000_000.0
# A healthy continuous loop carries 30–60° of phase margin; the budget's erosion is read
# against the conservative end.
NOMINAL_PHASE_MARGIN_DEG: float = 30.0

# Inner loop: the C3b terminal tracking loop — target-side astrometric tracker → crosslink →
# onboard ZEM → valve, all fast and local (the crosslink is the ≤100 km closing hop).
TERMINAL_LOOP: ControlLoop = ControlLoop(
    name="terminal inner loop (drag rejection / ZEM aim)",
    bandwidth_hz=DRAG_REJECTION_BANDWIDTH_HZ,
    sources=(
        LatencySource("target tracker exposure", 1.0e-3),
        LatencySource("crosslink one-way (≤100 km)", comms_one_way_s(TERMINAL_CROSSLINK_RANGE_M)),
        LatencySource("onboard compute", 1.0e-3),
        LatencySource("valve open (ms-class)", 5.0e-3),
    ),
)
# Outer loop: the midcourse replan — coordinator-node ranging → off-board UKF/corrector →
# uplink → impulsive burn.  Larger τ, but discrete: no phase margin to erode.
MIDCOURSE_LOOP: ControlLoop = ControlLoop(
    name="midcourse outer loop (discrete replan)",
    bandwidth_hz=None,
    sources=(
        LatencySource(
            "node comms round-trip (2000 km)", 2.0 * comms_one_way_s(COORDINATOR_NODE_RANGE_M)
        ),
        LatencySource("coordinator compute (UKF + corrector)", 50.0e-3),
        LatencySource("uplink one-way (2000 km)", comms_one_way_s(COORDINATOR_NODE_RANGE_M)),
    ),
)


class DeadTimeBuffer:
    """A pure FIFO delay line: the value the loop acts on is ``delay_ticks`` steps stale.

    The lumped per-loop dead-time (§16.8): the actuation at step k reflects nav information
    from step ``k − delay_ticks`` (measure → compute → uplink → actuate).  Until the line
    fills it holds the first value — at hand-off the loop knows its midcourse-delivered state
    and the τ-delayed terminal stream has not yet arrived.  ``delay_ticks = 0`` is identity.
    """

    def __init__(self, delay_ticks: int) -> None:
        self._line: deque[Vec3] = deque(maxlen=delay_ticks + 1)

    def step(self, value: Vec3) -> Vec3:
        self._line.append(value)
        return self._line[0]


# The measured C3b terminal geometry (ADR 0014 C3b findings): interception speed, the
# 400 mN / 25 kg actuator, the ~247 s descent from the 800 km hand-off.
TERMINAL_SPEED_M_S: float = 10_780.0
TERMINAL_MASS_KG: float = 25.0
TERMINAL_A_MAX_M_S2: float = PEAK_THRUST_LIMIT_N / TERMINAL_MASS_KG
TERMINAL_SPAN_S: float = 247.0
# The nominal tracker grade (ADR 0015): sets the significance gate's threshold so the
# noiseless sweep keeps the C3b loop's gating dynamics, with no stochastic noise injected.
NOMINAL_SIGMA_THETA_RAD: float = 10e-6


def fly_terminal_loop(
    entry_m: float,
    latency_s: float,
    *,
    rng: np.random.Generator | None = None,
    sigma_theta_rad: float = NOMINAL_SIGMA_THETA_RAD,
    dt_s: float = 1.0,
    span_s: float = TERMINAL_SPAN_S,
    speed_m_s: float = TERMINAL_SPEED_M_S,
    a_max_m_s2: float = TERMINAL_A_MAX_M_S2,
    mass_kg: float = TERMINAL_MASS_KG,
    gain: float = ZEM_GAIN,
) -> float:
    """Fly the C3b ZEM loop on the double integrator with a dead-time, return the lateral miss.

    Mirrors the C3b loop (`terminal_tick` over a ZOH double integrator) and inserts a
    :class:`DeadTimeBuffer` on the measured position — the loop acts on a τ-stale nav fix.
    ``sigma_theta_rad`` sets the significance-gate threshold (the gate dynamics are part of
    the loop); ``rng=None`` is the noiseless run (no stochastic nav error injected, the
    default for the τ-sweep), while a generator injects the Gauss–Markov tracker noise for
    the deferred Rung-D combined-stress case.
    """
    noise = NavNoiseProcess(TrackerGrade(sigma_theta_rad, 1.0), rng) if rng is not None else None
    buffer = DeadTimeBuffer(round(latency_s / dt_s))
    x = np.array([0.0, entry_m, 0.0])
    v = np.zeros(3)
    attitude: Vec3 | None = None
    for k in range(int(span_s / dt_s)):
        t_go = span_s - k * dt_s
        range_m = speed_m_s * t_go
        err = (
            np.asarray(noise.sample((range_m, 0.0, 0.0), dt_s))
            if noise is not None
            else np.zeros(3)
        )
        measured = x + err
        delayed = np.asarray(buffer.step((measured[0], measured[1], measured[2])))
        zem_est = -(delayed + v * t_go)
        tick = terminal_tick(
            (zem_est[0], zem_est[1], zem_est[2]),
            knowledge_sigma_m=sigma_theta_rad * range_m,
            t_go_s=t_go,
            feedforward_m_s2=(0.0, 0.0, 0.0),
            attitude_dir=attitude,
            control_period_s=dt_s,
            mass_kg=mass_kg,
            gain=gain,
        )
        attitude = tick.attitude_dir
        a = tick.thrust_n / mass_kg * np.asarray(attitude) if tick.fire else np.zeros(3)
        x = x + v * dt_s + 0.5 * a * dt_s**2
        v = v + a * dt_s
    return float(np.hypot(x[1], x[2]))


# The homing miss has "materially" grown — the dead-time has started to bite — once it
# doubles vs the zero-delay run; below that the sweep reads as flat.
DEGRADATION_FACTOR: float = 2.0


@dataclass(frozen=True)
class TauSweepPoint:
    """One swept dead-time: the noiseless lateral miss and its depth in control ticks."""

    latency_s: float
    delay_ticks: int
    lateral_m: float


def tau_sweep(
    latencies_s: tuple[float, ...],
    *,
    entry_m: float = 400.0,
    dt_s: float = 1.0,
    span_s: float = TERMINAL_SPAN_S,
) -> tuple[TauSweepPoint, ...]:
    """Sweep the dead-time on the noiseless C3b loop; one deterministic lateral miss per τ.

    Run at the C3b-validated cadence so each point is a clean dead-time perturbation of that
    operating loop.  Noiseless, so the miss is a deterministic function of τ — the sweep
    isolates the dead-time mechanism from the (separately characterized) tracker-noise floor.
    Sub-tick latencies (``round(τ/dt) == 0``) are the identity buffer and reproduce the
    zero-delay miss exactly — the structural statement that the ms-class budget is invisible
    at this cadence.
    """
    return tuple(
        TauSweepPoint(
            latency_s=latency,
            delay_ticks=round(latency / dt_s),
            lateral_m=fly_terminal_loop(entry_m, latency, rng=None, dt_s=dt_s, span_s=span_s),
        )
        for latency in latencies_s
    )


@dataclass(frozen=True)
class LatencyFinding:
    """The C4 measurement set: the per-loop budget + phase check, and the noiseless τ-sweep."""

    loops: tuple[ControlLoop, ...]
    sweep: tuple[TauSweepPoint, ...]
    sweep_cadence_hz: float
    sweep_entry_m: float
    degradation_factor: float = DEGRADATION_FACTOR

    @property
    def control_period_s(self) -> float:
        return 1.0 / self.sweep_cadence_hz

    @property
    def terminal_budget_tau_s(self) -> float:
        """The fast (terminal) loop's dead-time budget — the operating point the sweep dwarfs."""
        return self.loops[0].tau_s

    @property
    def baseline_lateral_m(self) -> float:
        """The zero-delay homing miss the sweep reads degradation against."""
        return min(self.sweep, key=lambda p: p.latency_s).lateral_m

    @property
    def tolerated_latency_s(self) -> float:
        """The largest swept τ whose miss is still within the degradation factor of baseline."""
        threshold = self.degradation_factor * self.baseline_lateral_m
        return max(p.latency_s for p in self.sweep if p.lateral_m <= threshold)

    @property
    def breakdown_latency_s(self) -> float | None:
        """The smallest swept τ whose miss exceeds the degradation factor (``None`` if none do)."""
        threshold = self.degradation_factor * self.baseline_lateral_m
        failed = [p.latency_s for p in self.sweep if p.lateral_m > threshold]
        return min(failed) if failed else None

    @property
    def budget_margin(self) -> float:
        """How many dead-time budgets the homing loop absorbs before the miss doubles."""
        return self.tolerated_latency_s / self.terminal_budget_tau_s


def latency_finding(
    extra_latencies_s: tuple[float, ...] = (0.5, 1.0, 2.0, 5.0, 10.0, 20.0),
    *,
    entry_m: float = 400.0,
    cadence_hz: float = DRAG_REJECTION_BANDWIDTH_HZ,
    loops: tuple[ControlLoop, ...] = (TERMINAL_LOOP, MIDCOURSE_LOOP),
) -> LatencyFinding:
    """Assemble the C4 finding: the per-loop budget + the noiseless τ-sweep (the pure runner).

    The sweep grid always carries the zero-delay baseline and the terminal budget itself (so
    the report can show the budget landing sub-tick) alongside the τ ≫ budget probes.
    """
    budget = loops[0].tau_s
    latencies = tuple(sorted({0.0, budget, *extra_latencies_s}))
    return LatencyFinding(
        loops=loops,
        sweep=tau_sweep(latencies, entry_m=entry_m, dt_s=1.0 / cadence_hz),
        sweep_cadence_hz=cadence_hz,
        sweep_entry_m=entry_m,
    )


def _loop_lines(loop: ControlLoop) -> list[str]:
    lines = [f"  {loop.name}: τ = {loop.tau_s * 1e3:.1f} ms"]
    loss = loop.phase_margin_loss_deg
    if loss is None:
        lines[0] += " (discrete replan — no phase loop to erode)"
    else:
        lines[0] += f", bandwidth {loop.bandwidth_hz:g} Hz → phase loss ω_c·τ = {loss:.2f}°"
    lines += [f"      {s.label}: {s.tau_s * 1e3:.2f} ms" for s in loop.sources]
    return lines


def format_latency(finding: LatencyFinding) -> str:
    """One-screen C4 report: the per-loop budget + phase check, then the noiseless τ-sweep."""
    period = finding.control_period_s
    budget = finding.terminal_budget_tau_s
    baseline = finding.baseline_lateral_m
    phase_loss = finding.loops[0].phase_margin_loss_deg
    midcourse_tau_ms = finding.loops[1].tau_s * 1e3

    lines = ["C4 control-loop latency — dead-time budget + phase margin + τ-sweep (§16.8)"]
    lines.append("  Per-loop dead-time budget and phase erosion:")
    for loop in finding.loops:
        lines += _loop_lines(loop)

    lines.append(
        f"  Terminal-homing τ-sweep — noiseless (isolates dead-time from the C3b nav floor),"
        f" entry {finding.sweep_entry_m:.0f} m at {finding.sweep_cadence_hz:g} Hz:"
    )
    for p in finding.sweep:
        ratio = p.lateral_m / baseline
        verdict = "flat" if ratio <= finding.degradation_factor else "DEGRADES"
        budget_tag = "  ← the budget" if math.isclose(p.latency_s, budget, rel_tol=1e-9) else ""
        lines.append(
            f"    τ = {p.latency_s:6.3f} s ({p.delay_ticks:2d} tk):"
            f" lateral {p.lateral_m:7.3f} m ({ratio:5.1f}× baseline) — {verdict}{budget_tag}"
        )

    tolerated = finding.tolerated_latency_s
    breakdown = finding.breakdown_latency_s
    breakdown_text = f"{breakdown:g} s" if breakdown is not None else "no swept τ"
    lines.append(
        f"  → Terminal budget {budget * 1e3:.1f} ms = {budget / period * 100:.2f}%"
        f" of the {period:g} s control period — sub-tick, structurally zero."
    )
    if phase_loss is not None:
        lines.append(
            f"    Phase erosion ω_c·τ = {phase_loss:.2f}° on the"
            f" {finding.loops[0].bandwidth_hz:g} Hz drag-rejection loop — negligible vs a"
            f" {NOMINAL_PHASE_MARGIN_DEG:.0f}–60° margin; the {midcourse_tau_ms:.0f} ms midcourse"
            " latency rides a discrete replan (no phase loop to erode)."
        )
    lines.append(
        f"    Noiseless homing miss flat through τ ≈ {tolerated:g} s"
        f" (~{tolerated / period:.0f} control period ≈ {finding.budget_margin:.0f}× the budget);"
        f" degrades at {breakdown_text}."
    )
    lines.append(
        "    Deferred to Rung D: the combined entry-offset × tracker-noise stress at the"
        " dispersion tail (measured on separate axes at C3b)."
    )
    return "\n".join(lines)
