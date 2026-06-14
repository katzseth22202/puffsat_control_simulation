"""Integration test for the ADR 0019 Lever-2 co-flyer phasing run (live JVM)."""

from __future__ import annotations

import pytest

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.runs.coflyer import coflyer_config, run_coflyer_phasing
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

from puffsat_sim import mission

pytestmark = pytest.mark.integration


def test_constant_a_maneuver_keeps_the_rocket_phase_locked() -> None:
    # +perigee / −apogee by the same amount holds the semi-major axis (period) exactly.
    rocket = coflyer_config(mission.NOMINAL_CONFIG)
    nominal_sum = mission.NOMINAL_CONFIG.perigee_alt_m + mission.NOMINAL_CONFIG.apogee_alt_m
    assert rocket.perigee_alt_m + rocket.apogee_alt_m == pytest.approx(nominal_sum)
    assert rocket.mean_anomaly_at_epoch_rad == mission.NOMINAL_CONFIG.mean_anomaly_at_epoch_rad


def test_coflyer_stays_close_and_in_the_gps_volume_through_the_window() -> None:
    """The phase-locked rocket holds angle-useful range and the GNSS volume across the window."""
    finding = run_coflyer_phasing(n_samples=15)

    assert finding.feasible
    assert finding.range_ok  # stays within the angle-useful design range
    assert finding.gps_ok  # stays inside the GPS constellation volume
    # It is genuinely close (the σ_θ·R lever) and genuinely aloft (above the 200 km crossing).
    assert finding.max_range_m < finding.angle_useful_range_m
    assert finding.min_rocket_alt_m > mission.INTERCEPTION_ALT_M
