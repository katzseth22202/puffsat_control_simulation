"""Integration test for the Rung D / D1.1 closed-loop train ensemble (live JVM, a tiny train)."""

from __future__ import annotations

import math

import pytest

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.runs.train import run_train_dispersion
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

from puffsat_sim.train import TrainDispersionSpec

pytestmark = pytest.mark.integration


def test_d1_1_train_ensemble_flies_and_reduces_end_to_end() -> None:
    """The D1.1 slice end-to-end: a 3-unit train homes through the funnel, reduced about its
    centroid.  Flown at the *achievable* 3 µrad grade (σ_θ budget gate) so capture is reliably
    high — this pins the chain mechanics; the grade-sensitivity verdict is the measured finding.
    """
    finding = run_train_dispersion(
        TrainDispersionSpec(n_units=3, tracker_sigma_theta_rad=3e-6), train_index=0
    )

    # Capture reduction is well-formed and capture-grade about the centroid.
    assert finding.capture.n_units == 3
    assert math.isfinite(finding.capture.centroid_drift_m)
    assert 0.0 <= finding.capture.capture_about_centroid <= 1.0
    assert finding.capture.capture_about_centroid >= 0.5

    # Propellant ledger and the deorbit-diagnostic perigee are physical.
    assert finding.terminal_dv_max_m_s > 0.0
    assert finding.within_budget  # worst-unit mission Δv under 2% at Isp 50
    assert finding.perigee_max_m < 100_000.0  # below the Kármán line — deorbit-good
