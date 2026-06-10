"""Tests for the pure C1 nav-feasibility sweep harness (ADR 0012) — cell grid,
node-constellation geometry, apogee reference, RTN covariance rotation, the
per-cell LinCov evaluation, and the envelope report (no JVM)."""

import numpy as np
import pytest

from puffsat_sim.constants import EARTH_RADIUS_M, WGS84_MU
from puffsat_sim.mission import APOGEE_ALT_M, PERIGEE_ALT_M
from puffsat_sim.nav_feasibility import (
    NavFeasibilityCell,
    NavFeasibilitySpec,
    apogee_state,
    covariance_to_rtn,
    evaluate_cell,
    format_nav_feasibility,
    nav_feasibility_cells,
    node_directions_rtn,
    sweep_nav_feasibility,
)
from puffsat_sim.orbital_math import keplerian_elements


class TestNavFeasibilityCells:
    def test_nominal_cell_plus_one_axis_at_a_time_variations(self) -> None:
        spec = NavFeasibilitySpec()
        cells = nav_feasibility_cells(spec)

        nominal = cells[0]
        assert nominal.axis == "nominal"
        assert nominal.range_sigma_m == spec.nominal_range_sigma_m
        assert nominal.doppler_sigma_m_s == spec.nominal_doppler_sigma_m_s
        assert nominal.cadence_hz == spec.nominal_cadence_hz

        swept_axes = {c.axis for c in cells[1:]}
        assert swept_axes == {
            "range_sigma_m",
            "doppler_sigma_m_s",
            "cadence_hz",
            "cone_half_angle_rad",
            "n_nodes",
            "q_accel_m_s2",
        }
        for cell in cells[1:]:
            off_axis_knobs = {
                "range_sigma_m": cell.range_sigma_m == spec.nominal_range_sigma_m,
                "doppler_sigma_m_s": cell.doppler_sigma_m_s == spec.nominal_doppler_sigma_m_s,
                "cadence_hz": cell.cadence_hz == spec.nominal_cadence_hz,
                "cone_half_angle_rad": cell.cone_half_angle_rad == spec.nominal_cone_half_angle_rad,
                "n_nodes": cell.n_nodes == spec.nominal_n_nodes,
                "q_accel_m_s2": cell.q_accel_m_s2 == spec.nominal_q_accel_m_s2,
            }
            assert all(held for axis, held in off_axis_knobs.items() if axis != cell.axis)
            assert not off_axis_knobs[cell.axis]

    def test_range_only_appears_as_a_none_doppler_cell(self) -> None:
        cells = nav_feasibility_cells(NavFeasibilitySpec())

        assert any(c.doppler_sigma_m_s is None for c in cells)

    def test_cell_indices_are_sequential(self) -> None:
        cells = nav_feasibility_cells(NavFeasibilitySpec())

        assert [c.cell_index for c in cells] == list(range(len(cells)))


class TestNodeDirections:
    def test_directions_sit_on_the_cone_at_equal_azimuth_spacing(self) -> None:
        half_angle = np.deg2rad(45.0)
        directions = node_directions_rtn(4, half_angle)

        assert len(directions) == 4
        for d in directions:
            assert np.linalg.norm(d) == pytest.approx(1.0)
            assert d[0] == pytest.approx(np.cos(half_angle))
        azimuths = sorted(np.arctan2(d[2], d[1]) % (2 * np.pi) for d in directions)
        gaps = np.diff(azimuths + [azimuths[0] + 2 * np.pi])
        np.testing.assert_allclose(gaps, np.pi / 2, rtol=1e-12)

    def test_zero_half_angle_collapses_to_a_collinear_constellation(self) -> None:
        directions = node_directions_rtn(4, 0.0)

        for d in directions:
            np.testing.assert_allclose(d, [1.0, 0.0, 0.0], rtol=0, atol=1e-15)


class TestApogeeState:
    def test_apogee_state_matches_the_reference_orbit(self) -> None:
        x = apogee_state()

        r_a = EARTH_RADIUS_M + APOGEE_ALT_M
        a, _ = keplerian_elements(PERIGEE_ALT_M, APOGEE_ALT_M)
        assert np.linalg.norm(x[:3]) == pytest.approx(r_a)
        assert x[:3] @ x[3:6] == pytest.approx(0.0)
        vis_viva = np.sqrt(WGS84_MU * (2.0 / r_a - 1.0 / a))
        assert np.linalg.norm(x[3:6]) == pytest.approx(vis_viva)


class TestCovarianceToRtn:
    def test_pure_radial_position_variance_lands_on_the_r_pos_axis(self) -> None:
        x = apogee_state()
        cov = np.zeros((6, 6))
        cov[0, 0] = 9.0

        rtn = covariance_to_rtn(cov, x)

        assert rtn[0, 0] == pytest.approx(9.0)
        assert np.sum(np.abs(rtn)) == pytest.approx(9.0)

    def test_rotation_preserves_the_trace(self) -> None:
        rng = np.random.default_rng(7)
        m = rng.normal(size=(6, 6))
        cov = m @ m.T
        inclined = np.array([1.5e8, 2e7, 1e7, -50.0, 400.0, 120.0])

        rtn = covariance_to_rtn(cov, inclined)

        assert np.trace(rtn) == pytest.approx(np.trace(cov))


_SHORT_ARC_SPEC = NavFeasibilitySpec(arc_duration_s=2_000.0)

_T_VEL_SENSITIVITY = 2.15e5


def _t_vel_phi() -> np.ndarray:
    """A synthetic Φ with only the C0-dominant entry: lateral miss per T-vel error."""
    phi = np.zeros((3, 6))
    phi[1, 4] = _T_VEL_SENSITIVITY
    return phi


def _cell_with(**overrides: object) -> NavFeasibilityCell:
    spec = _SHORT_ARC_SPEC
    knobs: dict[str, object] = {
        "cell_index": 0,
        "axis": "nominal",
        "range_sigma_m": spec.nominal_range_sigma_m,
        "doppler_sigma_m_s": spec.nominal_doppler_sigma_m_s,
        "cadence_hz": spec.nominal_cadence_hz,
        "cone_half_angle_rad": spec.nominal_cone_half_angle_rad,
        "n_nodes": spec.nominal_n_nodes,
        "q_accel_m_s2": spec.nominal_q_accel_m_s2,
    }
    knobs.update(overrides)
    return NavFeasibilityCell(**knobs)  # type: ignore[arg-type]


class TestEvaluateCell:
    def test_nominal_cell_pins_velocity_and_threads_phi_to_the_lateral_miss(self) -> None:
        outcome = evaluate_cell(_cell_with(), _SHORT_ARC_SPEC, _t_vel_phi())

        t_vel_sigma = outcome.vel_sigma_rtn_m_s[1]
        assert t_vel_sigma < _SHORT_ARC_SPEC.prior_vel_sigma_m_s / 10.0
        assert outcome.lateral_miss_1sigma_m == pytest.approx(
            _T_VEL_SENSITIVITY * t_vel_sigma, rel=1e-6
        )
        assert outcome.meets_catch_radius == (
            outcome.lateral_miss_1sigma_m < _SHORT_ARC_SPEC.catch_radius_m
        )

    def test_dropping_doppler_costs_transverse_velocity_knowledge(self) -> None:
        with_doppler = evaluate_cell(_cell_with(), _SHORT_ARC_SPEC, _t_vel_phi())
        range_only = evaluate_cell(
            _cell_with(axis="doppler_sigma_m_s", doppler_sigma_m_s=None),
            _SHORT_ARC_SPEC,
            _t_vel_phi(),
        )

        assert range_only.vel_sigma_rtn_m_s[1] > with_doppler.vel_sigma_rtn_m_s[1]


_TINY_SWEEP_SPEC = NavFeasibilitySpec(
    arc_duration_s=1_000.0,
    range_sigma_values_m=(100.0,),
    doppler_sigma_values_m_s=(None,),
    cadence_values_hz=(0.003,),
    cone_half_angle_values_rad=(0.0,),
    n_nodes_values=(3,),
    q_accel_values_m_s2=(5e-7,),
)


class TestSweepNavFeasibility:
    def test_outcomes_align_with_the_cell_grid(self) -> None:
        result = sweep_nav_feasibility(_TINY_SWEEP_SPEC, _t_vel_phi())

        cells = nav_feasibility_cells(_TINY_SWEEP_SPEC)
        assert len(result.outcomes) == len(cells)
        assert tuple(o.cell for o in result.outcomes) == cells

    def test_report_carries_the_verdict_per_cell(self) -> None:
        result = sweep_nav_feasibility(_TINY_SWEEP_SPEC, _t_vel_phi())

        report = format_nav_feasibility(result)

        assert "C1 navigation feasibility" in report
        assert "nominal" in report
        assert "range-only" in report
        for axis in ("cadence_hz", "cone_half_angle_rad", "n_nodes", "q_accel_m_s2"):
            assert axis in report
