"""Per-run and per-ensemble result value types — pure (no JVM).

``RunRecord`` and ``EnsembleResult`` are *produced* by the JVM-side run loop
(:mod:`puffsat_sim.montecarlo`) but carry no Orekit objects, so they live here in a
pure module: this lets the forthcoming resume sink (ADR 0003) serialize and replay
them — and be unit-tested — without booting the JVM.  Dependency layering stays
acyclic: ``dispersion`` (base) ← ``control`` ← ``records``.
"""

from __future__ import annotations

from dataclasses import dataclass

from puffsat_sim.control import ControlAction
from puffsat_sim.dispersion import EnsembleStats, RunInputs, Vec3


@dataclass(frozen=True)
class RunRecord:
    """One run's outcome: inputs, the RTN miss, ToA error, perigee, and the control plan.

    For the open-loop capstone (``control=None``) ``control_log`` is empty,
    ``total_dv_m_s`` is 0, and ``converged`` is True — a superset of the open-loop
    record.  Stays deeply immutable (tuples, frozen ``ControlAction``) so it logs to a
    line-oriented sink and crosses processes safely (ADR 0003).
    """

    inputs: RunInputs
    miss_rtn_m: Vec3
    toa_miss_s: float
    perigee_alt_m: float
    crossing_position_m: Vec3
    crossing_velocity_m_s: Vec3
    control_log: tuple[ControlAction, ...]
    total_dv_m_s: float
    converged: bool
    iterations: int


@dataclass(frozen=True)
class EnsembleResult:
    """An ensemble's per-run records, aggregate statistics, and the nominal reference."""

    master_seed: int
    nominal_perigee_alt_m: float
    nominal_toa_s: float
    records: tuple[RunRecord, ...]
    stats: EnsembleStats
