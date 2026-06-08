"""Integration tests for the open-loop dispersion harness (live JVM, small N)."""

from __future__ import annotations

import pytest

from puffsat_sim.control import solve_apogee_correction
from puffsat_sim.dispersion import DispersionSpec

try:
    # Importing montecarlo boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.montecarlo import replay_inputs, run_ensemble
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
