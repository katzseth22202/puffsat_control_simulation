"""Unit tests for the pure C3c terminal-authority / tail-correction core (ADR 0014)."""

from __future__ import annotations

import math

from puffsat_sim.authority import (
    AuthorityPoint,
    TailAuthorityFinding,
    TailAuthoritySpec,
    authority_point,
    cheapest_trim,
    dv_per_km_m_s,
    format_tail_authority,
    funnel_growth_dv_m_s,
    handoff_alt_to_cover_tail_m,
    lateral_lever_m_per_m_s,
    saturation_dv_m_s,
    thrust_limited_radius_m,
    trim_dv_to_cover_tail_m_s,
    trim_point,
)
from puffsat_sim.coeff_requirement import BudgetEntry, rss_lateral

# The ADR 0004 actuator at the 25 kg wet mass: a_max = 0.4 N / 25 kg.
A_MAX = 0.016

# The full upstream lateral 1σ budget (ADR 0013): C1 nav + B1 erosion + Cr coefficient
# prior → RSS ≈ 224 m, the budget ADR 0014 decision 5 reads the 3σ tail against.
_BUDGET = (
    BudgetEntry("nav-induced lateral (C1)", 141.0),
    BudgetEntry("finite-burn erosion (B1)", 89.0),
    BudgetEntry("coefficient prior (Cr)", 149.0),
)


def test_thrust_limited_radius_is_half_a_t_squared() -> None:
    # C3b's anchor: ~246 s of descent from the 800 km hand-off → ~484 m isotropic funnel.
    assert thrust_limited_radius_m(A_MAX, 246.0) == 0.5 * A_MAX * 246.0**2
    assert math.isclose(thrust_limited_radius_m(A_MAX, 246.0), 484.128, rel_tol=1e-4)


def test_saturation_dv_is_the_a_t_ceiling() -> None:
    # The Δv ceiling a_max·t is C3b's residual cliff (~3.94 m/s at 246 s).
    assert saturation_dv_m_s(A_MAX, 246.0) == A_MAX * 246.0
    assert math.isclose(saturation_dv_m_s(A_MAX, 246.0), 3.936, rel_tol=1e-4)


def test_funnel_growth_dv_is_the_saturation_ceiling_at_the_grow_altitude() -> None:
    # Nulling a 671 m entry at the grown-funnel edge costs √(2·a·entry) = a·t where ½at² = entry.
    assert math.isclose(
        funnel_growth_dv_m_s(A_MAX, 671.0), math.sqrt(2.0 * A_MAX * 671.0), rel_tol=1e-12
    )
    t = math.sqrt(2.0 * 671.0 / A_MAX)
    assert math.isclose(
        funnel_growth_dv_m_s(A_MAX, 671.0), saturation_dv_m_s(A_MAX, t), rel_tol=1e-12
    )


def test_radius_grows_quadratically_with_descent_time() -> None:
    # Doubling the descent time (start the burn higher) quadruples the funnel.
    assert math.isclose(
        thrust_limited_radius_m(A_MAX, 492.0),
        4.0 * thrust_limited_radius_m(A_MAX, 246.0),
        rel_tol=1e-12,
    )


def test_lateral_lever_strips_the_along_velocity_component() -> None:
    # v along x̂; a ±node-Δv that shifts the crossing (2 m along x, 200 m along y) per 1 m/s
    # central step must report only the 100 m/(m/s) lateral (⊥ v) lever.
    lever = lateral_lever_m_per_m_s(
        crossing_plus_m=(1.0, 100.0, 0.0),
        crossing_minus_m=(-1.0, -100.0, 0.0),
        dv_m_s=1.0,
        crossing_velocity_m_s=(7000.0, 0.0, 0.0),
    )
    assert math.isclose(lever, 100.0, rel_tol=1e-12)


def test_dv_per_km_is_the_inverse_lever() -> None:
    assert math.isclose(dv_per_km_m_s(100.0), 10.0, rel_tol=1e-12)  # 1000 m / 100 (m/(m/s))
    assert math.isclose(dv_per_km_m_s(1.0e4), 0.1, rel_tol=1e-12)
    assert dv_per_km_m_s(0.0) == math.inf  # a dead node has no lateral authority


def test_authority_point_packages_radius_and_ceiling() -> None:
    pt = authority_point(handoff_alt_m=800_000.0, t_descent_s=246.0, a_max_m_s2=A_MAX)
    assert pt.handoff_alt_m == 800_000.0
    assert pt.t_descent_s == 246.0
    assert pt.radius_m == thrust_limited_radius_m(A_MAX, 246.0)
    assert pt.dv_ceiling_m_s == saturation_dv_m_s(A_MAX, 246.0)


def test_trim_point_carries_lever_and_cost() -> None:
    pt = trim_point(node_alt_m=30_000_000.0, lateral_lever_m_per_m_s=1.0e4)
    assert pt.node_alt_m == 30_000_000.0
    assert pt.lateral_lever_m_per_m_s == 1.0e4
    assert math.isclose(pt.dv_per_km_m_s, 0.1, rel_tol=1e-12)


def _finding() -> TailAuthorityFinding:
    # A synthetic finding: the funnel covers 480 m at the 800 km hand-off and grows with
    # altitude; the trim is cheap high and dead low.
    authority = (
        AuthorityPoint(600_000.0, 180.0, thrust_limited_radius_m(A_MAX, 180.0), A_MAX * 180.0),
        AuthorityPoint(800_000.0, 246.0, thrust_limited_radius_m(A_MAX, 246.0), A_MAX * 246.0),
        AuthorityPoint(1_500_000.0, 360.0, thrust_limited_radius_m(A_MAX, 360.0), A_MAX * 360.0),
    )
    trims = (
        trim_point(30_000_000.0, 1.0e4),
        trim_point(3_000_000.0, 1.0e2),
        trim_point(1_000_000.0, 5.0),
    )
    return TailAuthorityFinding(
        authority_points=authority,
        trim_points=trims,
        a_max_m_s2=A_MAX,
        measured_radius_m=500.0,
        measured_radius_alt_m=800_000.0,
        budget_entries=_BUDGET,
        tail_sigma=3.0,
    )


def test_tail_is_three_sigma_of_the_budget_rss() -> None:
    f = _finding()
    expected_tail = 3.0 * rss_lateral(_BUDGET)  # ≈ 671 m
    assert math.isclose(f.tail_m, expected_tail, rel_tol=1e-12)
    assert math.isclose(f.budget_rss_1sigma_m, rss_lateral(_BUDGET), rel_tol=1e-12)
    # The 500 m measured funnel leaves the rest uncovered.
    assert math.isclose(f.uncovered_tail_m, expected_tail - 500.0, rel_tol=1e-9)


def test_handoff_altitude_to_cover_tail_interpolates_the_authority_curve() -> None:
    f = _finding()
    alt = handoff_alt_to_cover_tail_m(f)
    assert alt is not None
    # the ~671 m tail sits between the 800 km (~484 m) and 1500 km (~1037 m) funnels, so the
    # covering burn-start interpolates strictly between those two swept altitudes
    assert 800_000.0 < alt < 1_500_000.0


def test_handoff_altitude_is_none_when_no_swept_funnel_reaches_the_tail() -> None:
    f = _finding()
    huge = TailAuthorityFinding(
        authority_points=f.authority_points,
        trim_points=f.trim_points,
        a_max_m_s2=A_MAX,
        measured_radius_m=500.0,
        measured_radius_alt_m=800_000.0,
        budget_entries=(BudgetEntry("huge", 5_000.0),),  # 3σ = 15 km, past the largest funnel
        tail_sigma=3.0,
    )
    assert handoff_alt_to_cover_tail_m(huge) is None


def test_cheapest_trim_is_the_highest_authority_node() -> None:
    f = _finding()
    cheap = cheapest_trim(f)
    assert cheap.node_alt_m == 30_000_000.0  # smallest dv_per_km
    # Killing the 672 m tail from the cheap node is a sliver of Δv.
    dv = trim_dv_to_cover_tail_m_s(f, cheap)
    assert math.isclose(dv, cheap.dv_per_km_m_s * f.tail_m / 1000.0, rel_tol=1e-12)
    assert dv < 0.1


def test_format_reports_both_curves_and_the_crossover() -> None:
    text = format_tail_authority(_finding())
    assert "C3c" in text
    assert "Authority curve" in text
    assert "MCC-2 cost curve" in text
    assert "Verdict" in text
    assert "MCC-2 vindicated" in text
    # the budget breakdown, the C3b anchor, and the cost contrast are surfaced
    assert "coefficient prior" in text
    assert "500 m" in text
    assert "cheaper" in text


def test_spec_defaults_sweep_handoff_and_node_altitudes() -> None:
    spec = TailAuthoritySpec()
    assert spec.handoff_altitudes_m[0] < spec.handoff_altitudes_m[-1]
    # nodes run high → low so the cost curve shows authority dying as the node drops
    assert spec.node_altitudes_m[0] > spec.node_altitudes_m[-1]
    assert spec.trim_dv_m_s > 0.0
    assert spec.tail_sigma == 3.0
