"""Tests for the pure C0 navigation module (ADR 0011) — nav-error sweep grid, Φ
assembly, covariance/tolerance, and the seeded Gaussian sampler (no JVM)."""

import numpy as np
import pytest

from puffsat_sim.dispersion import RunInputs, Vec3
from puffsat_sim.navigation import (
    NavSweepResult,
    NavSweepSpec,
    assemble_sensitivity,
    axis_tolerance,
    induced_miss_covariance,
    linearity_range,
    nav_grid_offsets,
    sample_nav_error,
    summarize_nav_requirement,
)
from puffsat_sim.records import RunRecord


def _nav_record(run_index: int, miss_rtn_m: Vec3, dv: float, converged: bool = True) -> RunRecord:
    """A minimal RunRecord carrying just the fields the nav summary reads."""
    return RunRecord(
        inputs=RunInputs(run_index, (0.0, 0.0, 0.0), 0.04, 0.02, 150.0, 15.0),
        miss_rtn_m=miss_rtn_m,
        toa_miss_s=0.0,
        perigee_alt_m=50_000.0,
        crossing_position_m=(0.0, 0.0, 0.0),
        crossing_velocity_m_s=(0.0, 0.0, 0.0),
        control_log=(),
        total_dv_m_s=dv,
        converged=converged,
        iterations=1,
    )


class TestNavGridOffsets:
    def test_each_cell_perturbs_one_component_plus_a_single_zero_cell(self) -> None:
        cells = nav_grid_offsets(NavSweepSpec(points_per_sign=2))

        zero_cells = [c for c in cells if all(v == 0.0 for v in c.offset_rtn6)]
        assert len(zero_cells) == 1  # exactly one nominal/zero cell

        for c in cells:
            nonzero = [v for v in c.offset_rtn6 if v != 0.0]
            assert len(nonzero) <= 1  # one-component-at-a-time

    def test_axes_use_position_then_velocity_ranges_both_signs(self) -> None:
        spec = NavSweepSpec(pos_range_m=(0.1, 100.0), vel_range_m_s=(1e-3, 1.0), points_per_sign=3)
        cells = nav_grid_offsets(spec)

        assert len(cells) == 6 * 2 * 3 + 1  # 6 axes × both signs × points + one zero cell

        pos_mags = np.geomspace(0.1, 100.0, 3)
        vel_mags = np.geomspace(1e-3, 1.0, 3)
        for axis in (0, 1, 2):  # R/T/N position
            mags = sorted(c.magnitude for c in cells if c.axis == axis)
            assert mags == pytest.approx(sorted([*pos_mags, *(-pos_mags)]))
        for axis in (3, 4, 5):  # R/T/N velocity
            mags = sorted(c.magnitude for c in cells if c.axis == axis)
            assert mags == pytest.approx(sorted([*vel_mags, *(-vel_mags)]))

        for c in cells:
            if c.axis >= 0:
                assert c.offset_rtn6[c.axis] == pytest.approx(c.magnitude)  # offset on its axis


class TestAssembleSensitivity:
    def test_recovers_known_phi_from_linear_responses(self) -> None:
        # For a linear truth map miss = Φ·offset, the central-difference slope at zero
        # recovers Φ exactly — the empirical statement of the ADR 0011 cancellation.
        rng = np.random.default_rng(0)
        phi_true = rng.normal(size=(3, 6))
        cells = nav_grid_offsets(NavSweepSpec(points_per_sign=3))
        misses = np.array([phi_true @ np.asarray(c.offset_rtn6, dtype=np.float64) for c in cells])

        phi = assemble_sensitivity(cells, misses)

        np.testing.assert_allclose(phi, phi_true, rtol=1e-9, atol=1e-9)


class TestInducedMissCovariance:
    def test_is_phi_sigma_phi_transpose_and_symmetric(self) -> None:
        rng = np.random.default_rng(1)
        phi = rng.normal(size=(3, 6))
        root = rng.normal(size=(6, 6))
        sigma = root @ root.T  # any SPD apogee-RTN nav covariance

        cov = induced_miss_covariance(phi, sigma)

        assert cov.shape == (3, 3)
        np.testing.assert_allclose(cov, phi @ sigma @ phi.T)
        np.testing.assert_allclose(cov, cov.T)


class TestAxisTolerance:
    def test_is_catch_radius_over_lateral_column_norm_radial_excluded(self) -> None:
        phi = np.zeros((3, 6))
        phi[1, 4] = 3.0  # transverse-velocity axis: lateral (T,N) column = (3, 4) → norm 5
        phi[2, 4] = 4.0
        phi[0, 4] = 100.0  # radial sensitivity must NOT count (R is pinned at the crossing)

        tol = axis_tolerance(phi, catch_radius_m=10.0)

        assert tol[4] == pytest.approx(10.0 / 5.0)

    def test_axis_with_no_lateral_sensitivity_is_unconstrained(self) -> None:
        phi = np.zeros((3, 6))
        phi[0, 2] = 7.0  # purely radial response on axis 2 → no lateral miss to bound
        assert np.isinf(axis_tolerance(phi, catch_radius_m=5.0)[2])


class TestLinearityRange:
    def test_perfectly_linear_response_spans_the_full_range(self) -> None:
        slope = np.array([2.0, -1.0, 0.5])
        mags = np.array([-100.0, -10.0, -1.0, 1.0, 10.0, 100.0])
        misses = np.outer(mags, slope)  # miss = slope · magnitude, exactly linear

        assert linearity_range(mags, misses, rel_tol=1e-6) == pytest.approx(100.0)

    def test_truncates_where_response_bends_beyond_tolerance(self) -> None:
        slope = np.array([2.0, 0.0, 0.0])
        mags = np.array([-100.0, -10.0, -1.0, 1.0, 10.0, 100.0])
        bend = np.outer(mags**2, np.array([0.0, 0.0, 0.01]))  # quadratic, bites only at large |m|
        misses = np.outer(mags, slope) + bend
        # relative deviation = 0.01·m² / (2|m|) = 0.005·|m|: within 1% only out to |m| = 1
        assert linearity_range(mags, misses, rel_tol=0.01) == pytest.approx(1.0)


class TestSampleNavError:
    def test_same_seed_reproduces_the_draw(self) -> None:
        sigma = np.diag([1.0, 2.0, 3.0, 0.1, 0.2, 0.3])
        first = sample_nav_error(np.random.default_rng(42), sigma)
        second = sample_nav_error(np.random.default_rng(42), sigma)
        assert first == second

    def test_sample_covariance_converges_to_sigma(self) -> None:
        # A correlated Σ (off-diagonal position block) so a diagonal-only bug would fail.
        sigma = np.array(
            [
                [4.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                [1.0, 9.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1e-2, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 4e-2, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 9e-2],
            ]
        )
        rng = np.random.default_rng(7)
        draws = np.array([sample_nav_error(rng, sigma) for _ in range(30000)])
        np.testing.assert_allclose(np.cov(draws, rowvar=False), sigma, atol=0.2)


class TestSummarizeNavRequirement:
    def _linear_result(self, phi_true: np.ndarray, dv: float) -> NavSweepResult:
        spec = NavSweepSpec(points_per_sign=3)
        cells = nav_grid_offsets(spec)
        records = tuple(
            _nav_record(
                c.cell_index,
                tuple(phi_true @ np.asarray(c.offset_rtn6, dtype=np.float64)),  # type: ignore[arg-type]
                dv=0.0 if c.axis < 0 else dv,
            )
            for c in cells
        )
        return NavSweepResult(spec=spec, cells=cells, records=records)

    def test_recovers_phi_lateral_sensitivity_full_linearity_and_phantom_dv(self) -> None:
        rng = np.random.default_rng(3)
        phi_true = rng.normal(size=(3, 6))
        req = summarize_nav_requirement(self._linear_result(phi_true, dv=0.2), rel_tol=1e-6)

        np.testing.assert_allclose(req.phi, phi_true, atol=1e-9)
        np.testing.assert_allclose(
            req.axis_lateral_sensitivity, np.linalg.norm(phi_true[1:3, :], axis=0)
        )
        # perfectly linear → linear out to the full swept magnitude (1e4 m pos, 1.0 m/s vel)
        np.testing.assert_allclose(
            req.axis_linearity_range, [1e4, 1e4, 1e4, 1.0, 1.0, 1.0], rtol=1e-9
        )
        assert req.phantom_dv_max_m_s == pytest.approx(0.2)
        assert req.converged_fraction == pytest.approx(1.0)
