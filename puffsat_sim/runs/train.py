"""Rung D / D1.1 closed-loop train ensemble — the JVM run for :mod:`puffsat_sim.train` (ADR 0018).

Strings the C-rung pieces into one closed-loop run per PuffSat and reduces a **train** of them
into the headline **P(capture) about the centroid**.  Per the user-confirmed D1.1 architecture
(Φ-composed entry + flown terminal): each unit's hand-off lateral entry offset is *sampled* from
the characterized C0/C1/C2a budget (the midcourse residual is provably linear in nav/coefficient
error — C0's Φ), and the C3b ZEM terminal loop is *flown* (:func:`run_guidance`) through the
catch-radius cliff, the significance-gate rectification, and the σ_θ·R tracker noise — the binding
nonlinearity.  The corrector-in-every-run brute-force validation, MCC-2 scheduling, and the
node-count Σ sweep are later D1.x sub-slices (ADR 0018 decisions 4/6).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace

import numpy as np

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from puffsat_sim import mission
from puffsat_sim.config import OrbitalConfig
from puffsat_sim.guidance import PlateMiss, TrackerGrade
from puffsat_sim.montecarlo import physics_from_inputs
from puffsat_sim.runs.guidance import GuidanceContext, build_guidance_context, run_guidance
from puffsat_sim.tracker_fusion import (
    Tracker,
    array_with_coflyer,
    fused_tracker_grade,
    single_target_detector,
    target_array_only,
    tracker_fusion_finding,
)
from puffsat_sim.train import (
    PooledCaptureFinding,
    TrainDispersionSpec,
    TrainEnsembleFinding,
    format_pooled_capture,
    format_train_ensemble,
    sample_train,
    sample_train_entry_offsets,
    summarize_pooled_capture,
    summarize_train_ensemble,
)

D1_MASTER_SEED: int = 20260613
SMOKE_N_UNITS: int = 16
SMOKE_N_TRAINS: int = 2
# The ADR 0019 single-tracker ceiling the re-run reproduces as the FAIL reference (D1.1 finding).
LEGACY_SINGLE_TRACKER_SIGMA_THETA_RAD: float = 10e-6


def run_train_dispersion(
    spec: TrainDispersionSpec,
    train_index: int,
    master_seed: int = D1_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    ctx: GuidanceContext | None = None,
    trackers: Sequence[Tracker] | None = None,
) -> TrainEnsembleFinding:
    """Fly one train's units through the closed terminal loop and reduce them (D1.1).

    Each unit carries its train-sampled coefficients/space-weather into the *truth* drag and its
    Φ-composed lateral entry offset into the funnel; the loop injects the σ_θ tracker noise.  Pass
    a shared ``ctx`` to fly several trains over one hand-off/aim-point build.  ``trackers`` re-keys
    the noise grade to a fused multi-tracker architecture (ADR 0019); ``None`` flies the spec's
    single-tracker ``tracker_sigma_theta_rad``.
    """
    ctx = ctx if ctx is not None else build_guidance_context(orbital_config)
    grade = (
        fused_tracker_grade(trackers, sigma_range_m=spec.tracker_sigma_range_m)
        if trackers is not None
        else TrackerGrade(
            sigma_theta_rad=spec.tracker_sigma_theta_rad, sigma_range_m=spec.tracker_sigma_range_m
        )
    )
    inputs = sample_train(master_seed, spec, train_index)
    entries = sample_train_entry_offsets(master_seed, spec, train_index)

    runs = [
        run_guidance(
            ctx,
            entry_offset_lateral_m=entries[j],
            grade=grade,
            control_period_s=spec.control_period_s,
            rng=np.random.default_rng((master_seed, train_index, j)),
            truth_physics=physics_from_inputs(inputs[j]),
        )
        for j in range(spec.n_units)
    ]
    return summarize_train_ensemble(
        [r.miss for r in runs],
        [r.plan for r in runs],
        [r.perigee_alt_m for r in runs],
        spec,
    )


def train_dispersion_report(
    n_trains: int = SMOKE_N_TRAINS,
    spec: TrainDispersionSpec | None = None,
    master_seed: int = D1_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Fly ``n_trains`` independent trains over one shared context and format each D1.1 finding."""
    spec = spec if spec is not None else TrainDispersionSpec(n_units=SMOKE_N_UNITS)
    ctx = build_guidance_context(orbital_config)
    return "\n\n".join(
        f"Train {t} (seed {master_seed}):\n"
        + format_train_ensemble(run_train_dispersion(spec, t, master_seed, ctx=ctx))
        for t in range(n_trains)
    )


def _rerun_block(
    label: str,
    finding: TrainEnsembleFinding,
    effective_sigma_theta_rad: float,
    margin: float | None,
) -> str:
    grade = f"effective σ_θ {effective_sigma_theta_rad * 1e6:.2f} µrad"
    grade += f" ({margin:.1f}× inside capture-grade)" if margin is not None else ""
    return f"[{label}] {grade}\n" + format_train_ensemble(finding)


def fused_train_rerun_report(
    n_units: int = SMOKE_N_UNITS,
    spec: TrainDispersionSpec | None = None,
    master_seed: int = D1_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Re-run D1.1 at the fused multi-tracker grade (ADR 0019 decision 4).

    Flies one train per architecture over one shared context: the legacy single-tracker 10 µrad
    ceiling (the D1.1 FAIL reference), then the fused target array and array+co-flyer — showing the
    multi-tracker architecture recovers the ~3 µrad effective grade D1.1 needs.
    """
    spec = spec if spec is not None else TrainDispersionSpec(n_units=n_units)
    ctx = build_guidance_context(orbital_config)

    legacy_spec = replace(spec, tracker_sigma_theta_rad=LEGACY_SINGLE_TRACKER_SIGMA_THETA_RAD)
    blocks = [
        _rerun_block(
            "legacy single-tracker 10 µrad ceiling",
            run_train_dispersion(legacy_spec, 0, master_seed, ctx=ctx),
            LEGACY_SINGLE_TRACKER_SIGMA_THETA_RAD,
            None,
        )
    ]
    architectures = [
        ("single target detector (σ_θ gate)", single_target_detector()),
        ("target 5-array (Lever 1)", target_array_only()),
        ("target 5-array + co-flyer (Levers 1+2)", array_with_coflyer()),
    ]
    for train_index, (label, trackers) in enumerate(architectures, start=1):
        fusion = tracker_fusion_finding(trackers)
        blocks.append(
            _rerun_block(
                label,
                run_train_dispersion(spec, train_index, master_seed, ctx=ctx, trackers=trackers),
                fusion.effective_sigma_theta_rad,
                fusion.margin,
            )
        )
    header = "ADR 0019 — D1.1 re-run at the fused terminal-nav grade (C3b noise re-keyed)"
    return header + "\n\n" + "\n\n".join(blocks)


def run_pooled_train_capture(
    n_trains: int = SMOKE_N_TRAINS,
    spec: TrainDispersionSpec | None = None,
    master_seed: int = D1_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    ctx: GuidanceContext | None = None,
    trackers: Sequence[Tracker] | None = None,
) -> PooledCaptureFinding:
    """Fly ``n_trains`` independent trains and pool every unit into one counted capture figure.

    The escape tail is driven by the *shared* (per-train) entry leg, so a single train's units
    cannot resolve it however many there are — this pools across trains (ADR 0018, 2026-07-22).
    Pass ``replace(spec, sigma_entry_lateral_shared_m=0.0)`` to fly the **retarget-credited**
    convention, in which the ±2 km centroid retarget absorbs that shared leg (ADR 0016).
    """
    spec = spec if spec is not None else TrainDispersionSpec(n_units=SMOKE_N_UNITS)
    ctx = ctx if ctx is not None else build_guidance_context(orbital_config)
    grade = (
        fused_tracker_grade(trackers, sigma_range_m=spec.tracker_sigma_range_m)
        if trackers is not None
        else TrackerGrade(
            sigma_theta_rad=spec.tracker_sigma_theta_rad, sigma_range_m=spec.tracker_sigma_range_m
        )
    )
    misses: list[PlateMiss] = []
    entries: list[tuple[float, float]] = []
    for train in range(1, n_trains + 1):
        inputs = sample_train(master_seed, spec, train)
        train_entries = sample_train_entry_offsets(master_seed, spec, train)
        entries.extend(train_entries)
        misses.extend(
            run_guidance(
                ctx,
                entry_offset_lateral_m=train_entries[j],
                grade=grade,
                control_period_s=spec.control_period_s,
                rng=np.random.default_rng((master_seed, train, j)),
                truth_physics=physics_from_inputs(inputs[j]),
            ).miss
            for j in range(spec.n_units)
        )
    return summarize_pooled_capture(misses, entries, n_trains)


def pooled_train_capture_report(
    n_trains: int = SMOKE_N_TRAINS,
    n_units: int = SMOKE_N_UNITS,
    master_seed: int = D1_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Both entry conventions at the committed 5-array grade, side by side (ADR 0018)."""
    as_flown = TrainDispersionSpec(n_units=n_units)
    retargeted = replace(as_flown, sigma_entry_lateral_shared_m=0.0)
    ctx = build_guidance_context(orbital_config)
    trackers = target_array_only()
    blocks = [
        f"[{label}]\n"
        + format_pooled_capture(
            run_pooled_train_capture(n_trains, spec, master_seed, ctx=ctx, trackers=trackers)
        )
        for label, spec in (
            ("as flown — funnel nulls shared ⊕ per-unit entry (conservative)", as_flown),
            ("retarget credited — funnel nulls per-unit entry only (ADR 0016 spec)", retargeted),
        )
    ]
    header = "ADR 0018 — pooled per-unit capture at the 5-array fused grade, both entry conventions"
    return header + "\n\n" + "\n\n".join(blocks)


if __name__ == "__main__":
    print(fused_train_rerun_report())
