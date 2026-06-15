"""Rung D / D1.x IS tail P(capture) — the JVM run for :mod:`puffsat_sim.tail_capture`.

Resolves the per-unit capture-failure tail that D1.1's 16-unit empirical could not (ADR 0018
decision 6).  Flies the same C3b ZEM terminal loop D1.1 uses (:func:`run_guidance`) but draws each
unit's hand-off entry offset from the variance-inflated **importance-sampling** proposal, recording
the plate-frame arrival; a parallel **brute-force** batch (nominal entries) cross-checks the
reweighting at a shallow radius.  The tracker noise is sampled fresh per trajectory (its own
``rng``), so it stays out of the importance weight; truth physics are held at nominal (the
feedforward then fully rejects drag), isolating the binding entry×noise tail driver.  The pure
reduction (the weighted estimator, CI/ESS, and the brute-force agreement check) lives in
:mod:`puffsat_sim.tail_capture`.
"""

from __future__ import annotations

import numpy as np

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from puffsat_sim import mission
from puffsat_sim.config import OrbitalConfig
from puffsat_sim.guidance import TrackerGrade
from puffsat_sim.runs.guidance import GuidanceContext, build_guidance_context, run_guidance
from puffsat_sim.tail_capture import (
    TailCaptureFinding,
    TailCaptureSpec,
    format_tail_capture,
    sample_entry_is,
    sample_entry_nominal,
    summarize_tail_capture,
)
from puffsat_sim.tracker_fusion import D1_CAPTURE_GRADE_SIGMA_THETA_RAD

TAIL_MASTER_SEED: int = 20260615
# Smoke sizes: the IS batch resolves the deep plate tail, the brute-force batch sees the shallow
# r_val tail for the unbiasedness cross-check. Scale both up for a publication-grade estimate.
SMOKE_IS_N: int = 400
SMOKE_BF_N: int = 400


def _fly_batch(
    ctx: GuidanceContext,
    entries: np.ndarray,
    grade: TrackerGrade,
    control_period_s: float,
    master_seed: int,
    batch_id: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Fly a batch of entry offsets through the terminal loop → (lateral norms, ToA errors)."""
    lateral = np.empty(entries.shape[0], dtype=np.float64)
    toa = np.empty(entries.shape[0], dtype=np.float64)
    for j, entry in enumerate(entries):
        run = run_guidance(
            ctx,
            entry_offset_lateral_m=(float(entry[0]), float(entry[1])),
            grade=grade,
            control_period_s=control_period_s,
            rng=np.random.default_rng((master_seed, batch_id, j)),
        )
        lateral[j] = run.miss.lateral_norm_m
        toa[j] = run.miss.toa_error_s
    return lateral, toa


def run_tail_capture(
    spec: TailCaptureSpec | None = None,
    is_n: int = SMOKE_IS_N,
    bf_n: int = SMOKE_BF_N,
    master_seed: int = TAIL_MASTER_SEED,
    grade: TrackerGrade | None = None,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    ctx: GuidanceContext | None = None,
) -> TailCaptureFinding:
    """Fly the IS + brute-force batches and reduce them to the tail-resolved P(capture) (D1.x).

    ``grade`` defaults to the achievable ~3.2 µrad capture grade (the D1.1 conditional); the
    entry draws come from :mod:`puffsat_sim.tail_capture` so the proposal/weights stay pure.
    """
    spec = spec if spec is not None else TailCaptureSpec()
    ctx = ctx if ctx is not None else build_guidance_context(orbital_config)
    grade = grade if grade is not None else TrackerGrade(D1_CAPTURE_GRADE_SIGMA_THETA_RAD, 1.0)

    is_entries, is_weights = sample_entry_is(np.random.default_rng((master_seed, 0)), spec, is_n)
    bf_entries = sample_entry_nominal(np.random.default_rng((master_seed, 1)), spec, bf_n)

    is_lat, is_toa = _fly_batch(ctx, is_entries, grade, spec.control_period_s, master_seed, 0)
    bf_lat, bf_toa = _fly_batch(ctx, bf_entries, grade, spec.control_period_s, master_seed, 1)

    return summarize_tail_capture(is_weights, is_lat, is_toa, bf_lat, bf_toa, spec)


def tail_capture_report(
    spec: TailCaptureSpec | None = None,
    is_n: int = SMOKE_IS_N,
    bf_n: int = SMOKE_BF_N,
    master_seed: int = TAIL_MASTER_SEED,
    grade: TrackerGrade | None = None,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Fly the tail batches at the achievable grade and format the D1.x finding."""
    finding = run_tail_capture(spec, is_n, bf_n, master_seed, grade, orbital_config=orbital_config)
    sigma_theta = grade.sigma_theta_rad if grade is not None else D1_CAPTURE_GRADE_SIGMA_THETA_RAD
    grade_str = f"{sigma_theta * 1e6:.1f} µrad" if sigma_theta is not None else "range-only"
    header = f"Terminal-nav grade: σ_θ {grade_str} (the D1.1 achievable capture grade)"
    return header + "\n" + format_tail_capture(finding)


if __name__ == "__main__":
    print(tail_capture_report())
