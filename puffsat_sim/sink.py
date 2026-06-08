"""JSONL resume sink for the Monte Carlo harness — pure (no JVM).

Completed ``RunRecord``s stream to a newline-delimited JSON file keyed by
``run_index`` (ADR 0003).  Recovery is run-granular: on restart, run only the indices
missing from the sink and re-summarize — leaning on per-run seed reproducibility
(``replay_inputs``), never within-run integrator snapshots.  Serialization, read-back,
and resume planning are all pure, so they unit-test without booting Orekit; only
``run_ensemble``'s append-as-you-go loop touches the JVM.

Process-sharded parallel writers are deferred (ADR 0003); the per-index keying is
chosen so that becomes additive (shards are concatenated JSONL).
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from puffsat_sim.control import ControlAction
from puffsat_sim.dispersion import DispersionSpec, RunInputs, Vec3, replay_inputs
from puffsat_sim.records import RunRecord


def _vec3(seq: Any) -> Vec3:
    return (float(seq[0]), float(seq[1]), float(seq[2]))


def _action_to_dict(action: ControlAction) -> dict[str, Any]:
    return {
        "node_label": action.node_label,
        "elapsed_s": action.elapsed_s,
        "dv_rtn_m_s": list(action.dv_rtn_m_s),
        "dv_mag_m_s": action.dv_mag_m_s,
    }


def _action_from_dict(d: Any) -> ControlAction:
    return ControlAction(
        node_label=d["node_label"],
        elapsed_s=d["elapsed_s"],
        dv_rtn_m_s=_vec3(d["dv_rtn_m_s"]),
        dv_mag_m_s=d["dv_mag_m_s"],
    )


def record_to_dict(record: RunRecord) -> dict[str, Any]:
    """A JSON-serializable dict for one RunRecord (the persisted schema)."""
    i = record.inputs
    return {
        "inputs": {
            "run_index": i.run_index,
            "dv_rtn_m_s": list(i.dv_rtn_m_s),
            "cd_area_over_mass": i.cd_area_over_mass,
            "cr_area_over_mass": i.cr_area_over_mass,
            "f10p7": i.f10p7,
            "ap": i.ap,
        },
        "miss_rtn_m": list(record.miss_rtn_m),
        "toa_miss_s": record.toa_miss_s,
        "perigee_alt_m": record.perigee_alt_m,
        "crossing_position_m": list(record.crossing_position_m),
        "crossing_velocity_m_s": list(record.crossing_velocity_m_s),
        "control_log": [_action_to_dict(a) for a in record.control_log],
        "total_dv_m_s": record.total_dv_m_s,
        "converged": record.converged,
        "iterations": record.iterations,
    }


def record_from_dict(d: Any) -> RunRecord:
    """Reconstruct a RunRecord from :func:`record_to_dict`'s output."""
    i = d["inputs"]
    return RunRecord(
        inputs=RunInputs(
            run_index=i["run_index"],
            dv_rtn_m_s=_vec3(i["dv_rtn_m_s"]),
            cd_area_over_mass=i["cd_area_over_mass"],
            cr_area_over_mass=i["cr_area_over_mass"],
            f10p7=i["f10p7"],
            ap=i["ap"],
        ),
        miss_rtn_m=_vec3(d["miss_rtn_m"]),
        toa_miss_s=d["toa_miss_s"],
        perigee_alt_m=d["perigee_alt_m"],
        crossing_position_m=_vec3(d["crossing_position_m"]),
        crossing_velocity_m_s=_vec3(d["crossing_velocity_m_s"]),
        control_log=tuple(_action_from_dict(a) for a in d["control_log"]),
        total_dv_m_s=d["total_dv_m_s"],
        converged=d["converged"],
        iterations=d["iterations"],
    )


def append_record(path: Path, record: RunRecord) -> None:
    """Append one record as a JSON line (the run-granular checkpoint)."""
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record_to_dict(record)) + "\n")


def read_records(path: Path) -> list[RunRecord]:
    """Read all records from a sink, tolerating a torn final line (interrupted append)."""
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    records: list[RunRecord] = []
    for n, raw in enumerate(lines):
        line = raw.strip()
        if not line:
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            if n == len(lines) - 1:
                break  # a half-written final line from a killed append — just re-run it
            raise ValueError(f"corrupt sink record at line {n + 1} of {path}") from exc
        records.append(record_from_dict(payload))
    return records


def plan_resume(
    existing: Iterable[RunRecord], master_seed: int, spec: DispersionSpec, n: int
) -> tuple[dict[int, RunRecord], list[int]]:
    """Split an ensemble into reusable records and the indices still to run.

    A sink record is reused only if its inputs match what ``(master_seed, spec)``
    replays for its index, so resuming with a different seed or spec is caught rather
    than silently mixed; records with index >= n are ignored.  (Controller consistency
    is the caller's responsibility — inputs are control-independent.)
    """
    reuse: dict[int, RunRecord] = {}
    for record in existing:
        idx = record.inputs.run_index
        if idx >= n:
            continue
        if record.inputs != replay_inputs(master_seed, spec, idx):
            raise ValueError(
                f"sink record {idx} does not match this master_seed/spec — "
                "use a fresh sink path or matching parameters"
            )
        reuse[idx] = record
    todo = [i for i in range(n) if i not in reuse]
    return reuse, todo
