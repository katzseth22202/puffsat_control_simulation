"""Integration test for the A3 controllability sweep (live JVM, small grid)."""

from __future__ import annotations

import pytest

from puffsat_sim.control import ControlPlan, PredictFn, Target, solve_apogee_correction
from puffsat_sim.sweep import (
    Controllability,
    SweepSpec,
    classify_controllability,
    grid_inputs,
    to_grid,
)

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.runs.sweep import run_sweep
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

pytestmark = pytest.mark.integration


def _miss_magnitude(miss_rtn_m: tuple[float, float, float]) -> float:
    return sum(c * c for c in miss_rtn_m) ** 0.5


def _a3_controller(predict: PredictFn, target: Target) -> ControlPlan:
    # A3 solver config (ADR 0007 decision 3): LM damping regularizes the near-singular
    # altitude-event direction; step cap raised to budget scale; max_iter headroom passed
    # here (not the global default) so A1/A2 keep their committed Newton config.
    return solve_apogee_correction(predict, target, lm=True, max_step_m_s=50.0, max_iter=15)


def test_sweep_maps_required_dv_against_coefficient_error() -> None:
    """A 3×3 perfect-model grid produces a sane controllability map (ADR 0007).

    Re-tuned corrector finding: a transverse apogee Δv has strong authority over the
    along-track crossing, so coefficient dispersion is cheaply correctable — the map comes
    out *controllable everywhere* at Δv far under even the 50 s-Isp budget, and ~flat in Cd
    (drag is negligible at this apogee) with the small cost driven by Cr (SRP over the coast).
    """
    spec = SweepSpec(cd_points=3, cr_points=3)
    result = run_sweep(spec, control=_a3_controller)

    # One record per grid point, in grid order, carrying the swept inputs and zero injection.
    assert len(result.records) == spec.cd_points * spec.cr_points
    for rec, inp in zip(result.records, grid_inputs(spec), strict=True):
        assert rec.inputs.cd_area_over_mass == pytest.approx(inp.cd_area_over_mass)
        assert rec.inputs.cr_area_over_mass == pytest.approx(inp.cr_area_over_mass)
        assert rec.inputs.dv_rtn_m_s == (0.0, 0.0, 0.0)

    # The dedicated factor-(1,1) reference: perfect model, no error → nulls to ~0 Δv.
    assert result.nominal.converged
    assert result.nominal.total_dv_m_s < 0.5
    assert _miss_magnitude(result.nominal.miss_rtn_m) < 2.0

    grid = to_grid(result.records, spec)
    assert grid.required_dv_m_s.shape == (3, 3)

    # Every point is controllable and nulls onto the aim — the re-tune fixed the false
    # non-convergence from the over-damped config.
    assert grid.converged.all()
    for rec in result.records:
        assert _miss_magnitude(rec.miss_rtn_m) < 2.0
        # Debris-disposal safety: every swept perigee stays low enough to deorbit.
        assert 0.0 < rec.perigee_alt_m < 130_000.0

    # The on-nominal centre cell is the cheapest; coefficient error costs a little Δv, but
    # the whole map sits far under even the conservative 50 s-Isp budget (~9.8 m/s).
    assert grid.required_dv_m_s[1, 1] < 0.1
    assert grid.required_dv_m_s.max() > grid.required_dv_m_s[1, 1]
    assert grid.required_dv_m_s.max() < 1.0

    # The map is Cr-driven and ~flat in Cd: spread along the Cd axis is far below the Cr axis.
    assert grid.required_dv_m_s.std(axis=0).max() < grid.required_dv_m_s.std(axis=1).max()

    # Post-processing budget overlay: controllable everywhere even at the conservative anchor.
    labels = classify_controllability(grid, isp_s=50.0)
    assert (labels == Controllability.CONTROLLABLE).all()
