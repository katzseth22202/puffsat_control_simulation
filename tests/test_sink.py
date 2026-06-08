"""Tests for the pure JSONL resume sink — serialization and resume planning (no JVM)."""

import json
from pathlib import Path

import pytest

from puffsat_sim.control import ControlAction
from puffsat_sim.dispersion import DispersionSpec, RunInputs, replay_inputs
from puffsat_sim.records import RunRecord
from puffsat_sim.sink import (
    append_record,
    plan_resume,
    read_records,
    record_from_dict,
    record_to_dict,
)


def _record(inputs: RunInputs) -> RunRecord:
    return RunRecord(
        inputs=inputs,
        miss_rtn_m=(1.0, -2.0, 3.5),
        toa_miss_s=0.5,
        perigee_alt_m=50_000.0,
        crossing_position_m=(7.0e6, 1.0, -2.0),
        crossing_velocity_m_s=(0.0, 7.5e3, 0.0),
        control_log=(ControlAction("apogee", 0.0, (0.01, -0.02, 0.03), 0.037),),
        total_dv_m_s=0.037,
        converged=True,
        iterations=2,
    )


class TestSerialization:
    def test_dict_round_trip(self) -> None:
        r = _record(replay_inputs(7, DispersionSpec(), 0))
        assert record_from_dict(record_to_dict(r)) == r

    def test_dict_is_json_serializable(self) -> None:
        d = record_to_dict(_record(replay_inputs(7, DispersionSpec(), 0)))
        assert json.loads(json.dumps(d)) == d

    def test_empty_control_log_round_trips(self) -> None:
        # The open-loop capstone logs an empty plan.
        r = RunRecord(
            inputs=replay_inputs(7, DispersionSpec(), 0),
            miss_rtn_m=(0.0, 0.0, 0.0),
            toa_miss_s=0.0,
            perigee_alt_m=50_000.0,
            crossing_position_m=(7.0e6, 0.0, 0.0),
            crossing_velocity_m_s=(0.0, 7.5e3, 0.0),
            control_log=(),
            total_dv_m_s=0.0,
            converged=True,
            iterations=0,
        )
        assert record_from_dict(record_to_dict(r)) == r


class TestJsonlFile:
    def test_append_then_read_round_trips(self, tmp_path: Path) -> None:
        path = tmp_path / "sink.jsonl"
        recs = [_record(replay_inputs(7, DispersionSpec(), i)) for i in range(3)]
        for r in recs:
            append_record(path, r)
        assert read_records(path) == recs

    def test_missing_file_reads_empty(self, tmp_path: Path) -> None:
        assert read_records(tmp_path / "absent.jsonl") == []

    def test_torn_final_line_is_tolerated(self, tmp_path: Path) -> None:
        path = tmp_path / "sink.jsonl"
        r = _record(replay_inputs(7, DispersionSpec(), 0))
        append_record(path, r)
        with path.open("a", encoding="utf-8") as f:
            f.write('{"inputs": {"run_index": 1, ')  # half-written append, no newline
        assert read_records(path) == [r]

    def test_corrupt_nonfinal_line_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "sink.jsonl"
        good = json.dumps(record_to_dict(_record(replay_inputs(7, DispersionSpec(), 0))))
        path.write_text(f"not json\n{good}\n", encoding="utf-8")
        with pytest.raises(ValueError, match="corrupt"):
            read_records(path)


class TestPlanResume:
    def test_reuses_matching_and_lists_complement(self) -> None:
        seed, spec, n = 42, DispersionSpec(), 4
        existing = [_record(replay_inputs(seed, spec, i)) for i in (0, 2)]
        reuse, todo = plan_resume(existing, seed, spec, n)
        assert set(reuse) == {0, 2}
        assert todo == [1, 3]

    def test_rejects_seed_or_spec_mismatch(self) -> None:
        seed, spec, n = 42, DispersionSpec(), 4
        stale = _record(replay_inputs(seed + 1, spec, 0))  # index 0 for a different seed
        with pytest.raises(ValueError, match="does not match"):
            plan_resume([stale], seed, spec, n)

    def test_ignores_indices_at_or_above_n(self) -> None:
        seed, spec = 42, DispersionSpec()
        existing = [_record(replay_inputs(seed, spec, 5))]
        reuse, todo = plan_resume(existing, seed, spec, n=2)
        assert reuse == {}
        assert todo == [0, 1]
