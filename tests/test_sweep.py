"""Tests for the pure A3 sweep core — grid generation, reshape, and overlays (no JVM)."""

import math

import numpy as np
import pytest

from puffsat_sim.constants import STANDARD_GRAVITY_M_S2
from puffsat_sim.dispersion import RunInputs
from puffsat_sim.records import RunRecord
from puffsat_sim.sweep import (
    Controllability,
    SweepGrid,
    SweepSpec,
    budget_dv_m_s,
    classify_controllability,
    grid_inputs,
    sigma_equivalent,
    to_grid,
)


def _record(
    run_index: int, cd: float, cr: float, *, dv: float, converged: bool, perigee: float = 50_000.0
) -> RunRecord:
    """A minimal RunRecord carrying just the fields the grid/overlays read."""
    return RunRecord(
        inputs=RunInputs(run_index, (0.0, 0.0, 0.0), cd, cr, 150.0, 15.0),
        miss_rtn_m=(0.0, 0.0, 0.0),
        toa_miss_s=0.0,
        perigee_alt_m=perigee,
        crossing_position_m=(0.0, 0.0, 0.0),
        crossing_velocity_m_s=(0.0, 0.0, 0.0),
        control_log=(),
        total_dv_m_s=dv,
        converged=converged,
        iterations=1,
    )


class TestGridInputs:
    def test_1d_cd_cut_sweeps_cd_holds_cr_zeroes_injection(self) -> None:
        spec = SweepSpec(
            cd_area_over_mass=0.04,
            cr_area_over_mass=0.02,
            cd_factor_range=(0.5, 2.0),
            cd_points=5,
            cr_points=1,
        )
        inputs = grid_inputs(spec)

        assert len(inputs) == 5
        expected_cd = 0.04 * np.geomspace(0.5, 2.0, 5)  # log-spaced ⇒ linear in σ-equivalent
        for ri, cd in zip(inputs, expected_cd, strict=True):
            assert ri.cd_area_over_mass == pytest.approx(cd)
            assert ri.cr_area_over_mass == pytest.approx(0.02)  # held at nominal
            assert ri.dv_rtn_m_s == (0.0, 0.0, 0.0)  # zero injection

    def test_2d_grid_is_full_cd_by_cr_product(self) -> None:
        spec = SweepSpec(
            cd_factor_range=(0.5, 2.0), cr_factor_range=(0.8, 1.25), cd_points=3, cr_points=2
        )
        inputs = grid_inputs(spec)
        assert len(inputs) == 6  # 3 × 2
        assert len({ri.cd_area_over_mass for ri in inputs}) == 3  # cd swept
        assert len({ri.cr_area_over_mass for ri in inputs}) == 2  # cr swept

    def test_single_point_axis_sits_at_nominal(self) -> None:
        spec = SweepSpec(cd_points=4, cr_points=1)
        inputs = grid_inputs(spec)
        assert all(ri.cr_area_over_mass == pytest.approx(spec.cr_area_over_mass) for ri in inputs)

    def test_log_symmetric_range_with_odd_points_straddles_nominal(self) -> None:
        spec = SweepSpec(
            cd_area_over_mass=0.04, cd_factor_range=(0.5, 2.0), cd_points=5, cr_points=1
        )
        cds = [ri.cd_area_over_mass for ri in grid_inputs(spec)]
        assert any(c == pytest.approx(0.04) for c in cds)  # nominal is a grid point
        assert min(cds) < 0.04 < max(cds)  # straddles below and above

    def test_run_index_is_sequential(self) -> None:
        inputs = grid_inputs(SweepSpec(cd_points=3, cr_points=3))
        assert [ri.run_index for ri in inputs] == list(range(9))

    def test_space_weather_held_nominal(self) -> None:
        spec = SweepSpec(cd_points=2, cr_points=2)
        inputs = grid_inputs(spec)
        assert all(ri.f10p7 == pytest.approx(spec.f10p7) for ri in inputs)
        assert all(ri.ap == pytest.approx(spec.ap) for ri in inputs)


class TestToGrid:
    def test_reshapes_row_major_cd_outer_cr_inner(self) -> None:
        # total_dv = run_index, so the grid cell at [cd_i, cr_j] must equal cd_i·cr_points + cr_j.
        spec = SweepSpec(cd_points=2, cr_points=3)
        records = tuple(
            _record(
                ri.run_index,
                ri.cd_area_over_mass,
                ri.cr_area_over_mass,
                dv=float(ri.run_index),
                converged=True,
            )
            for ri in grid_inputs(spec)
        )
        grid = to_grid(records, spec)

        assert grid.required_dv_m_s.shape == (2, 3)
        np.testing.assert_array_equal(grid.required_dv_m_s, [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])

    def test_sorts_out_of_order_records_by_run_index(self) -> None:
        spec = SweepSpec(cd_points=2, cr_points=3)
        ordered = [
            _record(
                ri.run_index,
                ri.cd_area_over_mass,
                ri.cr_area_over_mass,
                dv=float(ri.run_index),
                converged=True,
            )
            for ri in grid_inputs(spec)
        ]
        shuffled = tuple([ordered[4], ordered[0], ordered[5], ordered[2], ordered[3], ordered[1]])
        grid = to_grid(shuffled, spec)
        np.testing.assert_array_equal(grid.required_dv_m_s, [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])

    def test_axes_are_factor_scaled_nominal_coefficients(self) -> None:
        spec = SweepSpec(cd_area_over_mass=0.04, cr_area_over_mass=0.02, cd_points=3, cr_points=2)
        records = tuple(
            _record(
                ri.run_index, ri.cd_area_over_mass, ri.cr_area_over_mass, dv=0.0, converged=True
            )
            for ri in grid_inputs(spec)
        )
        grid = to_grid(records, spec)
        np.testing.assert_allclose(grid.cd_area_over_mass, 0.04 * np.geomspace(0.5, 2.0, 3))
        np.testing.assert_allclose(grid.cr_area_over_mass, 0.02 * np.geomspace(0.5, 2.0, 2))

    def test_converged_and_perigee_ride_along(self) -> None:
        spec = SweepSpec(cd_points=1, cr_points=2)
        records = (
            _record(0, 0.04, 0.01, dv=1.0, converged=True, perigee=48_000.0),
            _record(1, 0.04, 0.04, dv=99.0, converged=False, perigee=120_000.0),
        )
        grid = to_grid(records, spec)
        np.testing.assert_array_equal(grid.converged, [[True, False]])
        np.testing.assert_array_equal(grid.perigee_alt_m, [[48_000.0, 120_000.0]])

    def test_wrong_record_count_raises(self) -> None:
        spec = SweepSpec(cd_points=2, cr_points=2)  # needs 4
        records = (_record(0, 0.04, 0.02, dv=0.0, converged=True),)
        with pytest.raises(ValueError, match="expected 4 records"):
            to_grid(records, spec)


class TestSigmaEquivalent:
    def test_nominal_factor_is_zero_sigma(self) -> None:
        assert sigma_equivalent(np.array([1.0]), cv=0.20)[0] == pytest.approx(0.0)

    def test_one_s_factor_is_one_sigma(self) -> None:
        cv = 0.20
        s = math.sqrt(math.log(1.0 + cv * cv))
        assert sigma_equivalent(np.array([math.exp(s)]), cv)[0] == pytest.approx(1.0)

    def test_reciprocal_factors_are_opposite_sigma(self) -> None:
        k = sigma_equivalent(np.array([0.5, 2.0]), cv=0.20)
        assert k[0] == pytest.approx(-k[1])  # 0.5 = 1/2 ⇒ ln symmetric


class TestBudget:
    def test_linear_isp_relation(self) -> None:
        assert budget_dv_m_s(50.0) == pytest.approx(0.02 * 50.0 * STANDARD_GRAVITY_M_S2)
        assert budget_dv_m_s(200.0) == pytest.approx(0.02 * 200.0 * STANDARD_GRAVITY_M_S2)

    def test_anchors_bracket_the_paper_budget(self) -> None:
        # 50 s ≈ 9.8 m/s (paper claim fails ~3×), 200 s ≈ 39 m/s (ADR 0007 ~32–40 budget).
        assert budget_dv_m_s(50.0) == pytest.approx(9.807, abs=0.01)
        assert budget_dv_m_s(200.0) == pytest.approx(39.227, abs=0.01)


class TestClassifyControllability:
    def _grid(self, dv: list[float], converged: list[bool]) -> SweepGrid:
        return SweepGrid(
            cd_area_over_mass=np.array([0.04]),
            cr_area_over_mass=np.asarray([0.02] * len(dv)),
            cd_factors=np.array([1.0]),
            cr_factors=np.ones(len(dv)),
            required_dv_m_s=np.array([dv]),
            converged=np.array([converged], dtype=np.bool_),
            perigee_alt_m=np.zeros((1, len(dv))),
        )

    def test_three_regions_at_one_isp(self) -> None:
        # At Isp=70, budget ≈ 13.73 m/s. dv=5 under, dv=15 over; the dv=50 point is
        # uncontrollable (no root) and must NOT read as merely over-budget.
        grid = self._grid(dv=[5.0, 15.0, 50.0], converged=[True, True, False])
        labels = classify_controllability(grid, isp_s=70.0)
        np.testing.assert_array_equal(
            labels,
            [
                [
                    Controllability.CONTROLLABLE,
                    Controllability.OVER_BUDGET,
                    Controllability.UNCONTROLLABLE,
                ]
            ],
        )

    def test_higher_isp_pulls_points_back_under_budget(self) -> None:
        # The same 30 m/s point is over-budget at 70 s but controllable at 200 s (≈39 m/s).
        grid = self._grid(dv=[30.0], converged=[True])
        assert classify_controllability(grid, isp_s=70.0)[0, 0] == Controllability.OVER_BUDGET
        assert classify_controllability(grid, isp_s=200.0)[0, 0] == Controllability.CONTROLLABLE
