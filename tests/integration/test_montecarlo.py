"""Integration tests for the open-loop dispersion harness (live JVM, small N)."""

from __future__ import annotations

import pytest

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


def test_control_hook_rejected() -> None:
    """The §14.1 control hook must refuse a non-None controller (Rung-D only)."""
    with pytest.raises(NotImplementedError):
        run_ensemble(DispersionSpec(), n=1, master_seed=0, control=lambda *a: None)
