"""Unit tests for the pure truth-model validation core (ADR 0018)."""

from __future__ import annotations

import math

import numpy as np

from puffsat_sim.constants import WGS84_MU
from puffsat_sim.estimation import two_body_j2_flow
from puffsat_sim.truth_validation import (
    ConservationDrift,
    TruthValidationFinding,
    angular_momentum_magnitude,
    conservation_drift,
    independent_coast,
    max_fractional_drift,
    max_position_divergence_m,
    specific_energy_j_per_kg,
)


def _circular_state(radius_m: float = 7.0e6) -> np.ndarray:
    """A circular-orbit EME2000 state in the x–y plane."""
    speed = math.sqrt(WGS84_MU / radius_m)
    return np.array([radius_m, 0.0, 0.0, 0.0, speed, 0.0])


def test_specific_energy_matches_vis_viva_for_a_circular_orbit() -> None:
    r = 7.0e6
    energy = specific_energy_j_per_kg(_circular_state(r).reshape(1, 6))
    # Circular orbit: ε = −μ / 2r.
    assert math.isclose(float(energy[0]), -WGS84_MU / (2.0 * r), rel_tol=1e-12)


def test_angular_momentum_is_r_cross_v_magnitude() -> None:
    r = 7.0e6
    speed = math.sqrt(WGS84_MU / r)
    h = angular_momentum_magnitude(_circular_state(r).reshape(1, 6))
    assert math.isclose(float(h[0]), r * speed, rel_tol=1e-12)


def test_max_fractional_drift_is_peak_deviation_over_first() -> None:
    values = np.array([100.0, 101.0, 99.5, 102.0])
    assert math.isclose(max_fractional_drift(values), 2.0 / 100.0)


def test_two_body_coast_conserves_energy_and_angular_momentum() -> None:
    # A pure two-body history (j2=0) holds both constants of motion to RK4 precision.
    s0 = _circular_state()
    states = np.array([two_body_j2_flow(s0, t, j2=0.0) for t in np.linspace(0.0, 6000.0, 12)])
    drift = conservation_drift(states)
    assert drift.energy_frac < 1e-9
    assert drift.ang_mom_frac < 1e-9


def test_independent_coast_chains_the_flow_from_the_initial_state() -> None:
    s0 = _circular_state()
    times = np.array([0.0, 300.0, 900.0])
    coast = independent_coast(s0, times)
    assert coast.shape == (3, 6)
    np.testing.assert_array_equal(coast[0], s0)
    # The second sample is the flow over the first segment.
    np.testing.assert_allclose(coast[1], two_body_j2_flow(s0, 300.0))


def test_position_divergence_is_the_max_separation() -> None:
    a = np.array([[0.0, 0.0, 0.0, 0, 0, 0], [10.0, 0.0, 0.0, 0, 0, 0]])
    b = np.array([[0.0, 0.0, 0.0, 0, 0, 0], [13.0, 4.0, 0.0, 0, 0, 0]])
    assert math.isclose(max_position_divergence_m(a, b), 5.0)


def _passing_finding(**overrides: float) -> TruthValidationFinding:
    base: dict[str, float] = dict(
        convergence_divergence_m=0.0,
        crosscheck_divergence_m=15.0,
        orbit_scale_m=1.56e8,
        span_s=115_000.0,
        n_samples=25,
    )
    base.update(overrides)
    return TruthValidationFinding(
        conservation=ConservationDrift(energy_frac=5e-15, ang_mom_frac=6e-16),
        convergence_divergence_m=base["convergence_divergence_m"],
        crosscheck_divergence_m=base["crosscheck_divergence_m"],
        orbit_scale_m=base["orbit_scale_m"],
        span_s=base["span_s"],
        n_samples=int(base["n_samples"]),
    )


def test_default_thresholds_validate_the_reference_coast() -> None:
    f = _passing_finding()
    assert f.conservation_ok
    assert f.convergence_ok
    assert f.crosscheck_ok
    assert f.validated


def test_a_leaking_integrator_fails_conservation() -> None:
    f = TruthValidationFinding(
        conservation=ConservationDrift(energy_frac=1e-5, ang_mom_frac=6e-16),
        convergence_divergence_m=0.0,
        crosscheck_divergence_m=15.0,
        orbit_scale_m=1.56e8,
        span_s=115_000.0,
        n_samples=25,
    )
    assert not f.conservation_ok
    assert not f.validated


def test_a_gross_crosscheck_divergence_fails_tier_two() -> None:
    # A megametre divergence (e.g. a frame/μ bug) blows the cross-check threshold.
    f = _passing_finding(crosscheck_divergence_m=1.0e6)
    assert not f.crosscheck_ok
    assert not f.validated


def test_crosscheck_and_convergence_fracs_normalize_by_orbit_scale() -> None:
    f = _passing_finding(crosscheck_divergence_m=156.0, orbit_scale_m=1.56e8)
    assert math.isclose(f.crosscheck_frac, 1e-6)


def test_format_reports_both_tiers_and_the_verdict() -> None:
    from puffsat_sim.truth_validation import format_truth_validation

    text = format_truth_validation(_passing_finding())
    assert "Tier 1" in text
    assert "Tier 2" in text
    assert "VALIDATED" in text
    assert "independent" in text
