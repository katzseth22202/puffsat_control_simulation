"""Pure multi-tracker fusion gate — the early-error reduction of ADR 0019 (no JVM).

D1.1 found the combined entry×noise stress needs an **effective** terminal-nav grade of ~3 µrad,
and the σ_θ budget gate (ADR 0018, :mod:`puffsat_sim.tracker_budget`) showed a single conservative
detector achieves 3.2 µrad — resting on one 3 µrad bench-calibratable distortion floor. This gate
quantifies how a *multi-tracker* architecture of redundant detectors of that same bench-grade
(3.2 µrad) class recovers capture-grade nav *with margin* — hedging the floor proving optimistic,
not rescuing 10 µrad hardware (five true 10 µrad detectors fuse only to ~4.6 µrad and still fail;
the 10 µrad class needs N ≈ 12 to go marginal) — by two independent levers (ADR 0019):

* **Averaging (√N):** N independent detectors on one platform reduce the *independent* part of the
  per-detector σ_θ as ``1/√N`` while the *common-mode* floor (correlated distortion + beacon-shape
  asymmetry; here the smear residual) does not average.
* **Range (σ_θ·R):** a *closer* platform (the co-flying launch rocket at ~500 km vs the target's
  2603 km) has proportionally less lateral error, plus a relative-geometry floor (the GNSS-pinned
  rocket→target vector) added in quadrature.

The trackers fuse by inverse-variance into one **effective σ_θ** at the target design range — the
grade the C3b/D1.1 loop sees, read against the D1.1-validated capture-grade.  The phasing
feasibility of the co-flyer lever is a separate JVM check (:mod:`puffsat_sim.runs.coflyer`).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from puffsat_sim.anti_drag import PEAK_THRUST_LIMIT_N
from puffsat_sim.guidance import TrackerGrade, homing_floor_m

# Per-detector σ_θ split (reproducing the σ_θ budget gate's 3.2 µrad RSS, ADR 0018): the
# independent part (distortion 3.0 ⊕ gyro 0.58 ⊕ photon 0.01 µrad — separately calibrated, own
# star astrometry → averages with √N) and the common-mode floor (the smear residual — a rigid-body
# rate shared across detectors; √N cannot cross it).  Beacon-shape asymmetry would add to the
# floor; the design rule (a clean point beacon) keeps it at ~0, so it is not defaulted in.
DETECTOR_INDEP_SIGMA_RAD: float = 3.06e-6
SMEAR_COMMON_SIGMA_RAD: float = 0.87e-6

# The geometry pins from D1.1 / ADR 0019.
TARGET_RANGE_M: float = 2603e3  # the hand-off→target range (D1.1 measured)
COFLYER_RANGE_M: float = 500e3  # the co-flying launch rocket's design range to the train centroid
COFLYER_RELGEOM_SIGMA_M: float = 2.0  # GNSS-pinned rocket→target lateral (low-altitude terminal)

# The ADR 0019 architecture sizing + the C3b along-LOS ranging σ the fused grade carries unchanged
# (fusion only re-keys the lateral σ_θ; ranging σ matches GuidanceSweepSpec / TrainDispersionSpec).
ARRAY_N_DETECTORS: int = 5  # the target X-pattern (decision 1: 5 independent detectors)
COFLYER_N_DETECTORS: int = 3  # the co-flying rocket's detector count (Stage 1's 0.76 µrad point)
DEFAULT_RANGING_SIGMA_M: float = 1.0

# The co-flyer phasing gate (Lever 2; the JVM check is :mod:`puffsat_sim.runs.coflyer`): through
# the PuffSat terminal window the rocket must stay close enough that its angle is the strong lever
# (≤ the design range, so the σ_θ·R advantage the fusion credits actually holds) AND low enough to
# sit in the unlocked-spaceborne-GNSS volume (below the GPS constellation) — which is what pins its
# rocket→target vector independently of the long baseline.
GPS_CEILING_M: float = 20_200e3  # GPS constellation altitude — unlocked spaceborne receiver volume

# The D1.1-validated capture-grade: the sweep cleared the σ ≤ 1.65 m criterion at 3.2 µrad (and
# failed at 10 µrad), so an effective σ_θ at or under this reads as capture-grade.
D1_CAPTURE_GRADE_SIGMA_THETA_RAD: float = 3.2e-6

# Homing-floor diagnostic inputs (the 400 mN / 25 kg actuator; ~interception closing speed).
PUFFSAT_MASS_KG: float = 25.0
A_MAX_M_S2: float = PEAK_THRUST_LIMIT_N / PUFFSAT_MASS_KG
CLOSING_SPEED_M_S: float = 10_780.0


@dataclass(frozen=True)
class Tracker:
    """One sensing platform: N independent detectors at a range, with a relative-geometry floor."""

    range_m: float
    n_detectors: int
    sigma_indep_rad: float = DETECTOR_INDEP_SIGMA_RAD
    sigma_common_rad: float = SMEAR_COMMON_SIGMA_RAD
    rel_geom_sigma_m: float = 0.0


def angular_sigma_theta_rad(tracker: Tracker) -> float:
    """Per-platform σ_θ ``√(σ_common² + σ_indep²/N)`` — only the independent part averages."""
    return math.sqrt(tracker.sigma_common_rad**2 + tracker.sigma_indep_rad**2 / tracker.n_detectors)


def lateral_sigma_m(tracker: Tracker) -> float:
    """One platform's lateral position σ at the target: ``σ_θ·R ⊕ relative-geometry floor``."""
    return math.hypot(angular_sigma_theta_rad(tracker) * tracker.range_m, tracker.rel_geom_sigma_m)


def fuse_lateral_sigma_m(trackers: Sequence[Tracker]) -> float:
    """Inverse-variance fusion of the platforms' lateral σ (independent measurements)."""
    inv_var = sum(1.0 / lateral_sigma_m(t) ** 2 for t in trackers)
    return 1.0 / math.sqrt(inv_var)


def effective_sigma_theta_rad(fused_lateral_m: float, design_range_m: float) -> float:
    """The single-tracker grade equivalent at ``design_range_m`` — read against the D1.1 sweep."""
    return fused_lateral_m / design_range_m


def target_array(n_detectors: int) -> Tracker:
    """N independent bench-grade (3.2 µrad-class) detectors on the target (rel-geom 0)."""
    return Tracker(range_m=TARGET_RANGE_M, n_detectors=n_detectors)


def coflyer(n_detectors: int) -> Tracker:
    """The co-flying launch rocket: a close platform with the GNSS-pinned rel-geometry floor."""
    return Tracker(
        range_m=COFLYER_RANGE_M, n_detectors=n_detectors, rel_geom_sigma_m=COFLYER_RELGEOM_SIGMA_M
    )


def single_target_detector() -> tuple[Tracker, ...]:
    """The baseline σ_θ-gate architecture: one target detector (the D1.1 single-tracker grade)."""
    return (target_array(1),)


def target_array_only() -> tuple[Tracker, ...]:
    """ADR 0019 Lever 1 alone: the N-detector target array (no co-flyer)."""
    return (target_array(ARRAY_N_DETECTORS),)


def array_with_coflyer() -> tuple[Tracker, ...]:
    """The full ADR 0019 architecture: the target array fused with the close co-flying rocket."""
    return (target_array(ARRAY_N_DETECTORS), coflyer(COFLYER_N_DETECTORS))


def fused_tracker_grade(
    trackers: Sequence[Tracker],
    design_range_m: float = TARGET_RANGE_M,
    sigma_range_m: float = DEFAULT_RANGING_SIGMA_M,
) -> TrackerGrade:
    """Re-key a multi-tracker architecture into the C3b/D1.1 loop's :class:`TrackerGrade`.

    ADR 0019 decision 4: the terminal-nav grade is re-read as *per-detector*; the system grade
    the ZEM loop sees is the **fused effective σ_θ** at the target design range.  The C3b noise
    model already consumes a scalar ``sigma_theta_rad`` (its σ_θ·R lateral draw), so fusion needs
    no new noise code — only this grade construction.  Ranging σ is unchanged by fusion.
    """
    fused = fuse_lateral_sigma_m(trackers)
    return TrackerGrade(
        sigma_theta_rad=effective_sigma_theta_rad(fused, design_range_m),
        sigma_range_m=sigma_range_m,
    )


@dataclass(frozen=True)
class TrackerFusionFinding:
    """The fused effective terminal-nav grade and its capture-grade verdict (ADR 0019)."""

    trackers: tuple[Tracker, ...]
    design_range_m: float
    fused_lateral_sigma_m: float
    effective_sigma_theta_rad: float
    capture_grade_sigma_theta_rad: float = D1_CAPTURE_GRADE_SIGMA_THETA_RAD

    @property
    def meets_capture_grade(self) -> bool:
        return self.effective_sigma_theta_rad <= self.capture_grade_sigma_theta_rad

    @property
    def margin(self) -> float:
        """How many times inside the D1.1 capture-grade the effective σ_θ sits."""
        return self.capture_grade_sigma_theta_rad / self.effective_sigma_theta_rad

    @property
    def homing_floor_m(self) -> float:
        """The zero-entry homing floor at the effective grade (diagnostic; D1.1's combined stress
        is worse, so the verdict reads against the effective σ_θ, not this floor)."""
        return homing_floor_m(self.effective_sigma_theta_rad, CLOSING_SPEED_M_S, A_MAX_M_S2)


def tracker_fusion_finding(
    trackers: Sequence[Tracker], design_range_m: float = TARGET_RANGE_M
) -> TrackerFusionFinding:
    """Fuse a set of trackers into the effective grade + capture verdict (the pure runner)."""
    fused = fuse_lateral_sigma_m(trackers)
    return TrackerFusionFinding(
        trackers=tuple(trackers),
        design_range_m=design_range_m,
        fused_lateral_sigma_m=fused,
        effective_sigma_theta_rad=effective_sigma_theta_rad(fused, design_range_m),
    )


def format_tracker_fusion(finding: TrackerFusionFinding) -> str:
    """One-screen multi-tracker fusion report (ADR 0019)."""
    verdict = "CAPTURE-GRADE" if finding.meets_capture_grade else "NOT capture-grade"
    lines = [
        "Multi-tracker fusion — effective terminal-nav grade (ADR 0019)",
    ]
    for t in finding.trackers:
        lines.append(
            f"  platform: {t.n_detectors}× detector @ {t.range_m / 1e3:.0f} km"
            f" → σ_θ {angular_sigma_theta_rad(t) * 1e6:.2f} µrad,"
            f" lateral {lateral_sigma_m(t):.2f} m"
            + (f" (rel-geom {t.rel_geom_sigma_m:g} m)" if t.rel_geom_sigma_m else "")
        )
    lines += [
        f"  Fused lateral σ {finding.fused_lateral_sigma_m:.2f} m"
        f" → effective σ_θ {finding.effective_sigma_theta_rad * 1e6:.2f} µrad"
        f" at the {finding.design_range_m / 1e3:.0f} km design range.",
        f"  vs the D1.1 capture-grade {finding.capture_grade_sigma_theta_rad * 1e6:g} µrad:"
        f" {finding.margin:.1f}× inside — homing floor {finding.homing_floor_m:.2f} m.",
        f"  → {verdict}: redundant bench-grade detectors fused vs the ~3 µrad D1.1 requirement.",
    ]
    return "\n".join(lines)


@dataclass(frozen=True)
class CoflyerPhasing:
    """The co-flying rocket's phasing feasibility over the PuffSat 800→200 km window (ADR 0019).

    The Lever-2 gate: a constant-semi-major-axis maneuver (raise perigee, lower apogee) keeps the
    rocket phase-locked to the descending train (same period → no secular drift); the question is
    whether, through the terminal window, it stays close (the σ_θ·R lever holds) and low (in the
    GNSS volume that pins its rocket→target vector).  Fed by :mod:`puffsat_sim.runs.coflyer`.
    """

    max_range_m: float  # peak rocket↔centroid separation through the window
    max_rocket_alt_m: float  # peak rocket altitude in the window (vs the GPS ceiling)
    min_rocket_alt_m: float  # lowest rocket altitude — is it still aloft at interception?
    window_alt_hi_m: float
    window_alt_lo_m: float
    angle_useful_range_m: float = COFLYER_RANGE_M
    gps_ceiling_m: float = GPS_CEILING_M

    @property
    def range_ok(self) -> bool:
        """Stays within the design range the fusion credits the co-flyer's σ_θ·R advantage at."""
        return self.max_range_m <= self.angle_useful_range_m

    @property
    def gps_ok(self) -> bool:
        """Stays inside the unlocked-spaceborne-GNSS volume (below the GPS constellation)."""
        return self.max_rocket_alt_m <= self.gps_ceiling_m

    @property
    def feasible(self) -> bool:
        return self.range_ok and self.gps_ok


def phasing_verdict(
    window_ranges_m: Sequence[float],
    window_rocket_alts_m: Sequence[float],
    window_alt_hi_m: float,
    window_alt_lo_m: float,
) -> CoflyerPhasing:
    """Reduce the sampled rocket↔centroid range and rocket altitude over the window to a verdict."""
    return CoflyerPhasing(
        max_range_m=max(window_ranges_m),
        max_rocket_alt_m=max(window_rocket_alts_m),
        min_rocket_alt_m=min(window_rocket_alts_m),
        window_alt_hi_m=window_alt_hi_m,
        window_alt_lo_m=window_alt_lo_m,
    )


def format_coflyer_phasing(finding: CoflyerPhasing) -> str:
    """One-screen co-flyer phasing-feasibility report (ADR 0019, the Lever-2 gate)."""
    verdict = "PHASING-FEASIBLE" if finding.feasible else "NOT phasing-feasible"
    range_mark = "✓" if finding.range_ok else "✗"
    gps_mark = "✓" if finding.gps_ok else "✗"
    return "\n".join(
        [
            "Co-flying rocket phasing feasibility (ADR 0019, Lever 2)",
            f"  PuffSat window {finding.window_alt_hi_m / 1e3:.0f}"
            f"→{finding.window_alt_lo_m / 1e3:.0f} km altitude.",
            f"  {range_mark} max rocket↔centroid range {finding.max_range_m / 1e3:.0f} km"
            f" vs the {finding.angle_useful_range_m / 1e3:.0f} km angle-useful design range.",
            f"  {gps_mark} rocket altitude {finding.min_rocket_alt_m / 1e3:.0f}–"
            f"{finding.max_rocket_alt_m / 1e3:.0f} km"
            f" vs the {finding.gps_ceiling_m / 1e3:.0f} km GPS ceiling.",
            f"  → {verdict}: the co-flyer lever {'holds' if finding.feasible else 'is dropped'}.",
        ]
    )
