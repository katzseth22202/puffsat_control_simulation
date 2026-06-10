"""Unit tests for the pure C2a coefficient-knowledge requirement core (ADR 0013).

The chain under test: a 1D A3 cut's per-cell Δv vectors → the gradient ∂Δv/∂c at
nominal → lateral miss sensitivity through C0's Φ → tolerance vs the ground prior,
plus the analytic SRP cross-check and the RSS error-budget ledger.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from puffsat_sim.coeff_requirement import (
    MEASURED_BUDGET,
    BudgetEntry,
    analytic_srp_dv,
    coefficient_sensitivity,
    coefficient_tolerance,
    cut_dv_vectors,
    dv_gradient,
    format_coeff_requirement,
    rss_lateral,
    summarize_coeff_requirement,
)
from puffsat_sim.constants import SRP_P0_PA
from puffsat_sim.control import ControlAction
from puffsat_sim.dispersion import RunInputs, Vec3
from puffsat_sim.records import RunRecord
from puffsat_sim.sweep import SweepResult, SweepSpec


def _record(run_index: int, dv_rtn: Vec3, *, converged: bool = True) -> RunRecord:
    mag = float(np.linalg.norm(dv_rtn))
    return RunRecord(
        inputs=RunInputs(
            run_index=run_index,
            dv_rtn_m_s=(0.0, 0.0, 0.0),
            cd_area_over_mass=0.04,
            cr_area_over_mass=0.02,
            f10p7=150.0,
            ap=15.0,
        ),
        miss_rtn_m=(0.0, 0.0, 0.0),
        toa_miss_s=0.0,
        perigee_alt_m=50_000.0,
        crossing_position_m=(0.0, 0.0, 0.0),
        crossing_velocity_m_s=(0.0, 0.0, 0.0),
        control_log=(ControlAction("apogee", 0.0, dv_rtn, mag),),
        total_dv_m_s=mag,
        converged=converged,
        iterations=3,
    )


def _cr_cut(dv_by_factor: dict[float, Vec3], *, converged: bool = True) -> SweepResult:
    """A 1D Cr-cut SweepResult whose factors are geomspace(0.5, 2, n) = the dict keys."""
    spec = SweepSpec(cd_points=1, cr_points=len(dv_by_factor))
    ordered = [dv_by_factor[f] for f in sorted(dv_by_factor)]
    records = tuple(_record(i, dv, converged=converged or i != 1) for i, dv in enumerate(ordered))
    return SweepResult(spec=spec, records=records, nominal=_record(-1, (0.0, 0.0, 0.0)))


def _linear_dv(factor: float, gradient: Vec3) -> Vec3:
    return (
        gradient[0] * (factor - 1.0),
        gradient[1] * (factor - 1.0),
        gradient[2] * (factor - 1.0),
    )


def _dominant_phi(t_vel_sensitivity: float = 2.15e5) -> np.ndarray:
    phi = np.zeros((3, 6))
    phi[1, 4] = t_vel_sensitivity  # lateral T-row, apogee T-vel column (C0 finding)
    return phi


class TestCutDvVectors:
    def test_extracts_factors_and_vectors_in_factor_order(self) -> None:
        g: Vec3 = (1e-3, 1e-2, 0.0)
        cut = _cr_cut({f: _linear_dv(f, g) for f in (0.5, 1.0, 2.0)})

        factors, dvs = cut_dv_vectors(cut)

        assert np.allclose(factors, [0.5, 1.0, 2.0])
        assert dvs.shape == (3, 3)
        assert np.allclose(dvs[0], _linear_dv(0.5, g))
        assert np.allclose(dvs[2], _linear_dv(2.0, g))

    def test_sums_multiple_control_actions(self) -> None:
        cut = _cr_cut({f: _linear_dv(f, (0.0, 1e-2, 0.0)) for f in (0.5, 1.0, 2.0)})
        doubled = tuple(
            RunRecord(
                **{
                    **rec.__dict__,
                    "control_log": rec.control_log + rec.control_log,
                }
            )
            for rec in cut.records
        )
        cut = SweepResult(spec=cut.spec, records=doubled, nominal=cut.nominal)

        _, dvs = cut_dv_vectors(cut)

        assert np.allclose(dvs[2], np.asarray(_linear_dv(2.0, (0.0, 1e-2, 0.0))) * 2.0)

    def test_rejects_a_2d_grid(self) -> None:
        spec = SweepSpec(cd_points=3, cr_points=3)
        grid = SweepResult(
            spec=spec,
            records=tuple(_record(i, (0.0, 0.0, 0.0)) for i in range(9)),
            nominal=_record(-1, (0.0, 0.0, 0.0)),
        )
        with pytest.raises(ValueError, match="1D cut"):
            cut_dv_vectors(grid)

    def test_rejects_non_converged_points(self) -> None:
        cut = _cr_cut(
            {f: _linear_dv(f, (0.0, 1e-2, 0.0)) for f in (0.5, 1.0, 2.0)}, converged=False
        )
        with pytest.raises(ValueError, match="non-converged"):
            cut_dv_vectors(cut)


class TestDvGradient:
    def test_recovers_a_linear_law_exactly(self) -> None:
        g = np.array([1e-3, 1e-2, -2e-3])
        factors = np.array([0.5, 1.0, 2.0])
        dvs = np.outer(factors - 1.0, g)

        assert np.allclose(dv_gradient(factors, dvs), g)

    def test_uses_the_bracketing_pair_nearest_nominal(self) -> None:
        factors = np.array([0.5, 0.8, 1.25, 2.0])
        dvs = np.zeros((4, 3))
        dvs[1] = (0.0, -1.0, 0.0)  # f = 0.8
        dvs[2] = (0.0, 1.0, 0.0)  # f = 1.25

        gradient = dv_gradient(factors, dvs)

        assert gradient[1] == pytest.approx(2.0 / 0.45)

    def test_requires_points_straddling_nominal(self) -> None:
        with pytest.raises(ValueError, match="straddle"):
            dv_gradient(np.array([1.1, 1.5, 2.0]), np.zeros((3, 3)))


class TestSensitivityAndTolerance:
    def test_dominant_entry_sensitivity(self) -> None:
        sens = coefficient_sensitivity(_dominant_phi(), np.array([0.0, 1e-2, 0.0]))
        assert sens == pytest.approx(2.15e3)

    def test_only_the_lateral_velocity_block_counts(self) -> None:
        phi = np.zeros((3, 6))
        phi[0, 4] = 1e6  # radial row — pinned by the altitude event
        phi[1, 1] = 1e6  # position column — a Δv error has no position component
        assert coefficient_sensitivity(phi, np.array([1e-2, 1e-2, 1e-2])) == 0.0

    def test_tolerance_is_radius_over_sensitivity(self) -> None:
        tol = coefficient_tolerance(_dominant_phi(), np.array([0.0, 1e-2, 0.0]), 5_000.0)
        assert tol == pytest.approx(5_000.0 / 2.15e3)

    def test_zero_gradient_is_unconstrained(self) -> None:
        tol = coefficient_tolerance(_dominant_phi(), np.zeros(3), 5_000.0)
        assert tol == math.inf


class TestAnalyticSrp:
    def test_impulse_is_pressure_times_coefficient_times_coast(self) -> None:
        dv = analytic_srp_dv(0.02, 115_700.0)
        assert dv == pytest.approx(0.02 * SRP_P0_PA * 115_700.0)
        assert dv == pytest.approx(1.055e-2, rel=1e-3)


class TestErrorBudget:
    def test_rss_is_pythagorean(self) -> None:
        entries = (BudgetEntry("a", 3.0), BudgetEntry("b", 4.0))
        assert rss_lateral(entries) == pytest.approx(5.0)

    def test_measured_budget_carries_the_known_findings(self) -> None:
        labels = " ".join(e.label for e in MEASURED_BUDGET)
        assert "C1" in labels and "B1" in labels


class TestSummarize:
    def _requirement(self, prior: float = 0.2):  # noqa: ANN202 - helper
        g_cr: Vec3 = (0.0, 1e-2, 0.0)
        cr_cut = _cr_cut({f: _linear_dv(f, g_cr) for f in (0.5, 1.0, 2.0)})
        cd_spec = SweepSpec(cd_points=3, cr_points=1)
        cd_cut = SweepResult(
            spec=cd_spec,
            records=tuple(_record(i, (0.0, 0.0, 0.0)) for i in range(3)),
            nominal=_record(-1, (0.0, 0.0, 0.0)),
        )
        return summarize_coeff_requirement(
            _dominant_phi(),
            cd_cut=cd_cut,
            cr_cut=cr_cut,
            catch_radius_m=5_000.0,
            prior_sigma_factor=prior,
            coast_duration_s=115_700.0,
        )

    def test_cr_axis_tolerance_and_coverage(self) -> None:
        req = self._requirement()

        assert req.cr.lateral_sensitivity_m == pytest.approx(2.15e3)
        assert req.cr.tolerance_factor == pytest.approx(5_000.0 / 2.15e3)
        assert req.cr.prior_lateral_miss_m == pytest.approx(0.2 * 2.15e3)
        assert req.cr.covered_by_prior
        assert req.cr.margin == pytest.approx((5_000.0 / 2.15e3) / 0.2)

    def test_cd_axis_is_unconstrained(self) -> None:
        req = self._requirement()

        assert req.cd.tolerance_factor == math.inf
        assert req.cd.covered_by_prior
        assert req.cd.prior_lateral_miss_m == 0.0

    def test_a_tight_radius_flags_not_covered(self) -> None:
        g_cr: Vec3 = (0.0, 1e-2, 0.0)
        cr_cut = _cr_cut({f: _linear_dv(f, g_cr) for f in (0.5, 1.0, 2.0)})
        cd_cut = SweepResult(
            spec=SweepSpec(cd_points=3, cr_points=1),
            records=tuple(_record(i, (0.0, 0.0, 0.0)) for i in range(3)),
            nominal=_record(-1, (0.0, 0.0, 0.0)),
        )
        req = summarize_coeff_requirement(
            _dominant_phi(),
            cd_cut=cd_cut,
            cr_cut=cr_cut,
            catch_radius_m=100.0,
            prior_sigma_factor=0.2,
        )

        assert not req.cr.covered_by_prior

    def test_analytic_cross_check_ratio(self) -> None:
        req = self._requirement()

        assert req.analytic_srp_dv_m_s == pytest.approx(0.02 * SRP_P0_PA * 115_700.0)
        assert req.measured_over_analytic == pytest.approx(1e-2 / (0.02 * SRP_P0_PA * 115_700.0))


class TestFormat:
    def test_report_carries_verdicts_and_budget(self) -> None:
        g_cr: Vec3 = (0.0, 1e-2, 0.0)
        cr_cut = _cr_cut({f: _linear_dv(f, g_cr) for f in (0.5, 1.0, 2.0)})
        cd_cut = SweepResult(
            spec=SweepSpec(cd_points=3, cr_points=1),
            records=tuple(_record(i, (0.0, 0.0, 0.0)) for i in range(3)),
            nominal=_record(-1, (0.0, 0.0, 0.0)),
        )
        req = summarize_coeff_requirement(
            _dominant_phi(),
            cd_cut=cd_cut,
            cr_cut=cr_cut,
            catch_radius_m=5_000.0,
            prior_sigma_factor=0.2,
            coast_duration_s=115_700.0,
        )

        report = format_coeff_requirement(req)

        assert "Cr·(A/m)" in report and "Cd·(A/m)" in report
        assert "COVERED" in report
        assert "UNCONSTRAINED" in report
        assert "RSS" in report
        assert "C1" in report and "B1" in report
        assert "analytic" in report
