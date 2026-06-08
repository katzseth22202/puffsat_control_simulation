"""Tests for the pure dispersion core — sampling, RTN geometry, statistics (no JVM)."""

import math

import numpy as np
import pytest

from puffsat_sim.dispersion import (
    DispersionSpec,
    rtn_basis,
    rtn_components,
    rtn_to_cartesian,
    sample_run_inputs,
    summarize,
)


class TestSampleRunInputs:
    def test_same_seed_reproduces_inputs(self) -> None:
        spec = DispersionSpec()
        a = sample_run_inputs(np.random.default_rng(123), spec, 7)
        b = sample_run_inputs(np.random.default_rng(123), spec, 7)
        assert a == b

    def test_different_seed_differs(self) -> None:
        spec = DispersionSpec()
        a = sample_run_inputs(np.random.default_rng(1), spec, 0)
        b = sample_run_inputs(np.random.default_rng(2), spec, 0)
        assert a != b

    def test_run_index_recorded(self) -> None:
        ri = sample_run_inputs(np.random.default_rng(0), DispersionSpec(), 42)
        assert ri.run_index == 42

    def test_zero_sigma_returns_nominal(self) -> None:
        spec = DispersionSpec(
            sigma_dv_radial_m_s=0.0,
            sigma_dv_transverse_m_s=0.0,
            sigma_dv_normal_m_s=0.0,
            sigma_cd_frac=0.0,
            sigma_cr_frac=0.0,
            sigma_f10p7_frac=0.0,
            sigma_ap_frac=0.0,
        )
        ri = sample_run_inputs(np.random.default_rng(99), spec, 0)
        assert ri.dv_rtn_m_s == (0.0, 0.0, 0.0)
        assert ri.cd_area_over_mass == pytest.approx(spec.cd_area_over_mass)
        assert ri.cr_area_over_mass == pytest.approx(spec.cr_area_over_mass)
        assert ri.f10p7 == pytest.approx(spec.f10p7)
        assert ri.ap == pytest.approx(spec.ap)

    def test_lognormal_draws_are_positive(self) -> None:
        spec = DispersionSpec()
        rng = np.random.default_rng(5)
        for _ in range(200):
            ri = sample_run_inputs(rng, spec, 0)
            assert ri.cd_area_over_mass > 0.0
            assert ri.f10p7 > 0.0
            assert ri.ap > 0.0

    def test_lognormal_median_is_nominal(self) -> None:
        # median-nominal convention: sample median of F10.7 ≈ nominal (150).
        spec = DispersionSpec()
        rng = np.random.default_rng(2024)
        draws = [sample_run_inputs(rng, spec, 0).f10p7 for _ in range(4000)]
        assert float(np.median(draws)) == pytest.approx(spec.f10p7, rel=0.03)


class TestRtnBasis:
    def test_orthonormal(self) -> None:
        r_hat, t_hat, n_hat = rtn_basis((7.0e6, 1.0e6, -2.0e6), (1.0e3, 7.0e3, 5.0e2))
        vecs = [np.asarray(v) for v in (r_hat, t_hat, n_hat)]
        for v in vecs:
            assert float(np.linalg.norm(v)) == pytest.approx(1.0, abs=1e-12)
        for i in range(3):
            for j in range(i + 1, 3):
                assert float(vecs[i] @ vecs[j]) == pytest.approx(0.0, abs=1e-12)

    def test_transverse_equals_velocity_at_apogee(self) -> None:
        # velocity ⟂ position (apsis): T must align with the velocity direction.
        _, t_hat, _ = rtn_basis((7.0e6, 0.0, 0.0), (0.0, 6.0e3, 0.0))
        assert t_hat == pytest.approx((0.0, 1.0, 0.0), abs=1e-12)

    def test_radial_points_outward(self) -> None:
        r_hat, _, _ = rtn_basis((7.0e6, 0.0, 0.0), (0.0, 6.0e3, 0.0))
        assert r_hat == pytest.approx((1.0, 0.0, 0.0), abs=1e-12)


class TestRtnProjection:
    def test_components_then_back_round_trips(self) -> None:
        basis = rtn_basis((7.0e6, 1.0e6, -2.0e6), (1.0e3, 7.0e3, 5.0e2))
        original = (123.0, -45.0, 6.0)
        comps = rtn_components(original, basis)
        recovered = rtn_to_cartesian(comps, basis)
        assert recovered == pytest.approx(original, abs=1e-9)

    def test_components_along_transverse_only(self) -> None:
        basis = rtn_basis((7.0e6, 0.0, 0.0), (0.0, 6.0e3, 0.0))  # T = +y
        c_r, c_t, c_n = rtn_components((0.0, 50.0, 0.0), basis)
        assert (c_r, c_t, c_n) == pytest.approx((0.0, 50.0, 0.0), abs=1e-9)


class TestSummarize:
    def test_two_runs_mean_and_extents(self) -> None:
        miss = np.array([[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]])
        toa = np.array([10.0, -10.0])
        perigee = np.array([48_000.0, 52_000.0])
        stats = summarize(miss, toa, perigee)
        assert stats.n == 2
        assert stats.miss_rtn_mean_m == pytest.approx((2.0, 2.0, 2.0))
        assert stats.toa_miss_mean_s == pytest.approx(0.0)
        assert stats.perigee_alt_min_m == pytest.approx(48_000.0)
        assert stats.perigee_alt_max_m == pytest.approx(52_000.0)
        # Sample std (ddof=1) of [1,3] is sqrt(2).
        assert stats.miss_rtn_std_m[0] == pytest.approx(math.sqrt(2.0))

    def test_single_run_has_zero_spread(self) -> None:
        stats = summarize(np.array([[5.0, -3.0, 1.0]]), np.array([2.0]), np.array([50_000.0]))
        assert stats.n == 1
        assert stats.miss_rtn_std_m == (0.0, 0.0, 0.0)
        assert stats.toa_miss_std_s == 0.0
        assert stats.perigee_alt_min_m == pytest.approx(50_000.0)
