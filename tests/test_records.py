"""Tests for the pure result value types (no JVM).

Living under ``make test`` (which boots no JVM and has no orekit-data) is itself the
proof of the ADR 0003 relocation: ``RunRecord`` / ``EnsembleResult`` import and
construct without Orekit, so the forthcoming resume sink can serialize them purely.
"""

import dataclasses

import pytest

from puffsat_sim.control import ControlAction
from puffsat_sim.dispersion import EnsembleStats, RunInputs
from puffsat_sim.records import EnsembleResult, RunRecord


def _run_record(converged: bool = True) -> RunRecord:
    return RunRecord(
        inputs=RunInputs(0, (0.0, 0.1, 0.0), 0.04, 0.02, 150.0, 15.0),
        miss_rtn_m=(1.0, 2.0, 3.0),
        toa_miss_s=0.5,
        perigee_alt_m=50_000.0,
        crossing_position_m=(7.0e6, 0.0, 0.0),
        crossing_velocity_m_s=(0.0, 7.5e3, 0.0),
        control_log=(ControlAction("apogee", 0.0, (0.0, 0.1, 0.0), 0.1),),
        total_dv_m_s=0.1,
        converged=converged,
        iterations=2,
    )


def _stats() -> EnsembleStats:
    z3 = (0.0, 0.0, 0.0)
    return EnsembleStats(
        n=1,
        miss_rtn_mean_m=z3,
        miss_rtn_std_m=z3,
        miss_rtn_cov_m2=(z3, z3, z3),
        toa_miss_mean_s=0.0,
        toa_miss_std_s=0.0,
        perigee_alt_mean_m=50_000.0,
        perigee_alt_std_m=0.0,
        perigee_alt_min_m=50_000.0,
        perigee_alt_max_m=50_000.0,
        total_dv_mean_m_s=0.1,
        total_dv_std_m_s=0.0,
        total_dv_max_m_s=0.1,
        converged_fraction=1.0,
    )


class TestRunRecord:
    def test_constructs_without_jvm(self) -> None:
        r = _run_record()
        assert r.total_dv_m_s == pytest.approx(0.1)
        assert r.control_log[0].node_label == "apogee"

    def test_is_frozen(self) -> None:
        r = _run_record()
        with pytest.raises(dataclasses.FrozenInstanceError):
            r.converged = False  # type: ignore[misc]

    def test_is_deeply_immutable(self) -> None:
        # Tuples (not lists) throughout, so it crosses threads/processes safely.
        r = _run_record()
        assert isinstance(r.miss_rtn_m, tuple)
        assert isinstance(r.control_log, tuple)
        assert isinstance(r.control_log[0], ControlAction)


class TestEnsembleResult:
    def test_holds_records_and_stats(self) -> None:
        result = EnsembleResult(
            master_seed=1,
            nominal_perigee_alt_m=50_000.0,
            nominal_toa_s=100.0,
            records=(_run_record(), _run_record(converged=False)),
            stats=_stats(),
        )
        assert result.stats.n == 1
        assert len(result.records) == 2
        assert isinstance(result.records, tuple)
