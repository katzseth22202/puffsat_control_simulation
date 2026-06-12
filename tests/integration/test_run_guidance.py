"""Integration test for the C3b closed ZEM terminal loop (live JVM, two descents)."""

from __future__ import annotations

import numpy as np
import pytest

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.runs.guidance import build_guidance_context, run_guidance
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

from puffsat_sim.guidance import CAPTURE_SIGMA_MAX_M, TrackerGrade

pytestmark = pytest.mark.integration


def test_closed_zem_loop_nulls_an_entry_offset_and_homes_through_tracker_noise() -> None:
    """The C3b slice end-to-end (ADR 0014/0015), bounds set from the measured sweep.

    Measured 2026-06-12: a noiseless 400 m lateral entry offset lands 0.019 m from the
    aim point (catch radius 500 m, cliff at 600 m); the nominal 10 µrad tracker grade
    lands RMS 1.07 m / max 2.7 m over 8 seeds, capture-grade vs the σ ≤ 1.65 m plate
    requirement.  Bounds are loosened ~an order so the test pins the physics, not the
    platform.
    """
    ctx = build_guidance_context()

    offset = run_guidance(ctx, entry_offset_m=400.0)
    assert offset.miss.lateral_norm_m < CAPTURE_SIGMA_MAX_M
    assert abs(offset.miss.toa_error_s) < 1e-3
    assert offset.plan.dv_m_s < 4.0
    assert offset.saturated_fraction < 1.0

    noisy = run_guidance(
        ctx,
        grade=TrackerGrade(sigma_theta_rad=10e-6, sigma_range_m=1.0),
        rng=np.random.default_rng(20260612),
    )
    assert noisy.miss.lateral_norm_m < 5.0
    assert abs(noisy.miss.toa_error_s) < 1e-3
    assert noisy.plan.peak_thrust_n <= 0.4
