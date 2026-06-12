"""Integration test for the C0 navigation-error sensitivity sweep (live JVM, ADR 0011)."""

from __future__ import annotations

import numpy as np
import pytest

from puffsat_sim.control import ControlPlan, PredictFn, Target, solve_apogee_correction
from puffsat_sim.navigation import NavSweepSpec, assemble_sensitivity, axis_tolerance

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.runs.navigation import run_nav_sweep
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

pytestmark = pytest.mark.integration

_TRANSVERSE_VELOCITY_AXIS = 4  # apogee-RTN: 0-2 position R/T/N, 3-5 velocity R/T/N


def _miss_magnitude(miss_rtn_m: tuple[float, float, float]) -> float:
    return sum(c * c for c in miss_rtn_m) ** 0.5


def _c0_controller(predict: PredictFn, target: Target) -> ControlPlan:
    # Same robust config as the A3 sweep (ADR 0007 decision 3): LM damping for the
    # near-singular altitude-event direction, budget-scale step cap, iteration headroom.
    return solve_apogee_correction(predict, target, lm=True, max_step_m_s=50.0, max_iter=15)


def test_nav_sweep_transmits_nav_error_to_a_residual_miss() -> None:
    """A small predict-side nav error produces an uncontrollable residual miss (ADR 0011).

    The zero cell is on target (truth == nominal), while a perturbed apogee state leaves a
    residual the corrector cannot null — the seam carries the nav error through to a measurable
    interception miss, and the assembled Φ / tolerance are finite and physical.
    """
    spec = NavSweepSpec(pos_range_m=(100.0, 100.0), vel_range_m_s=(0.01, 0.01), points_per_sign=1)
    result = run_nav_sweep(spec, control=_c0_controller)

    assert len(result.records) == 6 * 2 * 1 + 1  # 6 axes × both signs × 1 point + one zero cell
    assert result.records[0].inputs.dv_rtn_m_s == (0.0, 0.0, 0.0)  # zero injection (isolation)

    # The zero cell is the on-target reference: x_true == nominal, so the corrector nulls to ~0.
    zero_cell = next(r for c, r in zip(result.cells, result.records, strict=True) if c.axis == -1)
    assert zero_cell.converged
    assert _miss_magnitude(zero_cell.miss_rtn_m) < 2.0
    assert zero_cell.total_dv_m_s < 0.5

    # A transverse-velocity nav error (the dr_p/dv_a lever) transmits to a residual far above
    # the on-target reference — the corrector burns a phantom Δv yet cannot remove the miss.
    tvel = [r for c, r in zip(result.cells, result.records, strict=True) if c.axis == 4]
    assert tvel and all(_miss_magnitude(r.miss_rtn_m) > 5.0 for r in tvel)
    assert all(r.total_dv_m_s > 0.0 for r in tvel)  # phantom correction (ADR 0011 decision 5)

    # Φ assembly + tolerance are finite and the lever axis carries real lateral sensitivity.
    misses = np.array([r.miss_rtn_m for r in result.records], dtype=np.float64)
    phi = assemble_sensitivity(result.cells, misses)
    assert np.all(np.isfinite(phi))
    tol = axis_tolerance(phi, catch_radius_m=5_000.0)
    assert 0.0 < tol[_TRANSVERSE_VELOCITY_AXIS] < np.inf
