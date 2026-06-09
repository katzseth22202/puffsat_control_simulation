"""Tests for the pure A3 sweep core — deterministic grid generation (no JVM)."""

import numpy as np
import pytest

from puffsat_sim.sweep import SweepSpec, grid_inputs


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
