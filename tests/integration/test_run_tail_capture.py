"""Integration test for the Rung D / D1.x IS tail P(capture) slice (live JVM, tiny batches)."""

from __future__ import annotations

import math

import pytest

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.runs.tail_capture import run_tail_capture
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

from puffsat_sim.tail_capture import TailCaptureSpec

pytestmark = pytest.mark.integration


def test_tail_capture_flies_both_batches_and_reduces_end_to_end() -> None:
    """The IS + brute-force batches fly through the C3b loop and reduce to a well-formed finding.

    Tiny N pins the chain mechanics only — the tail-resolved verdict (ESS, CI tightness) is the
    measured finding from a full-size run, not something a smoke batch can assert.
    """
    finding = run_tail_capture(TailCaptureSpec(), is_n=24, bf_n=24)

    # The weighted estimator and its diagnostics are well-formed.
    assert 0.0 <= finding.p_escape <= 1.0
    assert 0.0 <= finding.p_capture <= 1.0
    assert finding.is_plate.ess > 0.0
    assert 0.5 < finding.is_plate.weight_sum_ratio < 1.5  # the likelihood ratio averages to ~1
    assert isinstance(finding.validated, bool)

    # The per-unit P(capture) interval brackets the point estimate, and the σ read is physical.
    lo, hi = finding.p_capture_ci
    assert lo <= finding.p_capture <= hi
    assert math.isfinite(finding.scatter_sigma_m)
    assert finding.scatter_sigma_m > 0.0
