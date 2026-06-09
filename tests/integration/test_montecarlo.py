"""Integration tests for the open-loop dispersion harness (live JVM, small N)."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pytest

from puffsat_sim.actuator import Actuator
from puffsat_sim.anti_drag import PEAK_SLEW_LIMIT_DEG_S, PEAK_THRUST_LIMIT_N
from puffsat_sim.control import solve_apogee_correction
from puffsat_sim.dispersion import DispersionSpec
from puffsat_sim.sink import read_records, record_to_dict

try:
    # Importing montecarlo boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.montecarlo import instrument_anti_drag, replay_inputs, run_ensemble
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

pytestmark = pytest.mark.integration


def _miss_magnitude(miss_rtn_m: tuple[float, float, float]) -> float:
    return sum(c * c for c in miss_rtn_m) ** 0.5


def test_ensemble_smoke() -> None:
    """A small ensemble runs end-to-end, aggregates, and replays."""
    spec = DispersionSpec()
    result = run_ensemble(spec, n=4, master_seed=20260608)

    assert result.stats.n == 4
    assert len(result.records) == 4
    # Nominal osculating perigee at the crossing is the low debris-disposal orbit.
    assert 30_000.0 < result.nominal_perigee_alt_m < 90_000.0
    for r in result.records:
        # Debris-disposal safety: every dispersed perigee stays low enough to deorbit.
        assert 0.0 < r.perigee_alt_m < 120_000.0
        # Bounded miss (catches an overshoot/garbage state).
        assert _miss_magnitude(r.miss_rtn_m) < 1.0e6
    # Replay (§14.2): a standalone reconstruction matches the run's recorded inputs.
    assert replay_inputs(result.master_seed, spec, 0) == result.records[0].inputs


def test_zero_dispersion_returns_to_nominal() -> None:
    """With every σ = 0 the single run reproduces the nominal crossing (miss ≈ 0)."""
    spec = DispersionSpec(
        sigma_dv_radial_m_s=0.0,
        sigma_dv_transverse_m_s=0.0,
        sigma_dv_normal_m_s=0.0,
        sigma_cd_frac=0.0,
        sigma_cr_frac=0.0,
        sigma_f10p7_frac=0.0,
        sigma_ap_frac=0.0,
    )
    result = run_ensemble(spec, n=1, master_seed=1)
    (record,) = result.records
    assert _miss_magnitude(record.miss_rtn_m) < 1.0
    assert record.perigee_alt_m == pytest.approx(result.nominal_perigee_alt_m, abs=1.0)


def test_closed_loop_nulls_correctable_runs_and_records_authority() -> None:
    """Rung A1: correctable runs null to ≈0 at physical Δv; an uncorrectable tail run is
    recorded as non-converged (the authority boundary), never a spurious huge Δv.

    Run 1 of this seed is a 2.4σ radial-injection draw whose ~28 km along-track miss has
    no sub-budget single-apogee solution (ADR 0003); the corrector must reject it, not
    emit the ~88 m/s re-phasing root.
    """
    result = run_ensemble(
        DispersionSpec(), n=2, master_seed=20260608, control=solve_apogee_correction
    )

    assert 0.0 < result.stats.converged_fraction <= 1.0  # the loop closes on ≥1 run
    converged = [r for r in result.records if r.converged]
    assert converged
    for r in converged:
        assert _miss_magnitude(r.miss_rtn_m) < 2.0  # nulled onto the nominal aim
        assert 0.0 < r.total_dv_m_s < 5.0  # a real local correction, not a far re-phase
    for r in result.records:
        # Ledger is internally consistent: one apogee action, magnitude = |Δv| (commanded
        # == applied at Rung A), whether or not the run converged.
        (action,) = r.control_log
        assert action.node_label == "apogee"
        assert action.dv_mag_m_s == pytest.approx(_miss_magnitude(action.dv_rtn_m_s), rel=1e-9)
        assert r.total_dv_m_s == pytest.approx(action.dv_mag_m_s)


def test_resume_reuses_sink_records_and_runs_only_the_complement(tmp_path: Path) -> None:
    """A sink-backed ensemble resumes: present indices are reused, not recomputed."""
    spec = DispersionSpec()
    seed = 20260608
    sink = tmp_path / "ensemble.jsonl"

    # First pass writes index 0.  Tamper a *non-input* field with a sentinel so that a
    # reused record is distinguishable from a recomputed one (the replay guard only
    # checks inputs, so the tampered record still validates on resume).
    run_ensemble(spec, n=1, master_seed=seed, sink_path=sink)
    rec0 = read_records(sink)[0]
    sentinel = dataclasses.replace(rec0, perigee_alt_m=-12345.0)
    sink.write_text(json.dumps(record_to_dict(sentinel)) + "\n", encoding="utf-8")

    # Resume at n=2: index 0 comes back from the sink (sentinel intact), index 1 is run.
    result = run_ensemble(spec, n=2, master_seed=seed, sink_path=sink)
    assert len(result.records) == 2
    assert result.records[0].perigee_alt_m == -12345.0  # reused, not recomputed
    assert result.records[1].perigee_alt_m > 0.0  # index 1 actually propagated
    assert {r.inputs.run_index for r in read_records(sink)} == {0, 1}


def test_finite_burn_reproduces_the_impulsive_null_within_a_small_erosion() -> None:
    """B1 tracer (§13 / ADR 0008): executing the A1 corrector's Δv as a finite apogee burn
    lands the interception within a small, measured erosion of the impulsive null.

    The corrector nulls its *impulsive* prediction onto the nominal aim, so the residual
    miss of the *finite*-executed run is precisely the predict≠execute erosion.  At apogee
    (a 2.7-day orbit, a ~tens-of-seconds burn) it is expected small; this records it, not
    asserts it away.
    """
    spec = DispersionSpec()
    seed = 20260608
    impulsive = run_ensemble(spec, n=1, master_seed=seed, control=solve_apogee_correction)
    finite = run_ensemble(
        spec,
        n=1,
        master_seed=seed,
        control=solve_apogee_correction,
        actuator=Actuator(isp_s=50.0),
    )

    assert impulsive.records[0].converged and finite.records[0].converged
    # Same commanded plan either way (predict is impulsive in both); only execution differs.
    assert finite.records[0].total_dv_m_s == pytest.approx(impulsive.records[0].total_dv_m_s)
    impulsive_miss = _miss_magnitude(impulsive.records[0].miss_rtn_m)
    finite_miss = _miss_magnitude(finite.records[0].miss_rtn_m)
    assert impulsive_miss < 2.0  # corrector nulls the impulsive prediction
    # Finite execution erodes the null by a real, bounded amount (this run: ~89 m, almost
    # all along-track — the burn centroid lands ~68 s past the apogee node).  Small vs the
    # km-scale open-loop dispersion, far above the 0.7 m impulsive residual: a measured
    # finding, not noise (ADR 0008 — a finite-aware targeter that centers the burn is deferred).
    assert impulsive_miss < finite_miss < 200.0


def test_anti_drag_feedforward_clears_the_actuator_limits_with_margin() -> None:
    """B3a (§13 / ADR 0009): the known-drag 600→200 km feedforward anti-drag requirement
    clears the 400 mN / 1°/s actuator limits with large margin under the *conservative*
    cannonball coefficient — confirming the paper's <2% propellant claim a fortiori.

    Measured (nominal descent): peak thrust ~17 mN (~24× under 400 mN), peak slew ~0.05°/s
    (~20× under 1°/s), anti-drag Δv ~0.015 m/s. The paper's 374 g / 400 mN is a deliberately
    stacked-pessimistic upper bound; the physical NRLMSISE requirement is far below it.
    """
    profile = instrument_anti_drag()

    assert profile.duration_s > 0.0
    assert 0.0 < profile.anti_drag_dv_m_s < 0.5  # tiny vs the paper's ~0.7 m/s baseline Δv
    assert profile.peak_thrust_n < PEAK_THRUST_LIMIT_N  # ~17 mN, well under 400 mN
    assert profile.peak_slew_rate_deg_s < PEAK_SLEW_LIMIT_DEG_S  # ~0.05°/s, well under 1°/s


def test_closed_loop_beats_open_loop_for_the_same_run() -> None:
    """The same dispersed run misses by km open-loop but ≈0 closed-loop (same seed)."""
    spec = DispersionSpec()
    seed = 20260608
    open_loop = run_ensemble(spec, n=1, master_seed=seed)
    closed_loop = run_ensemble(spec, n=1, master_seed=seed, control=solve_apogee_correction)

    # Same seed → identical injection / coefficients, so this isolates the corrector.
    assert open_loop.records[0].inputs == closed_loop.records[0].inputs
    assert _miss_magnitude(open_loop.records[0].miss_rtn_m) > 100.0
    assert _miss_magnitude(closed_loop.records[0].miss_rtn_m) < 2.0
