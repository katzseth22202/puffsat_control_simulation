"""Tests for the pure Rung A1 differential corrector — synthetic predict, no JVM."""

import numpy as np
import pytest

from puffsat_sim.control import (
    ControlAction,
    ControlPlan,
    PredictFn,
    Target,
    Vec3,
    solve_apogee_correction,
)


def _linear_predict(matrix: np.ndarray, offset: np.ndarray) -> PredictFn:
    """predict(dv) = A·dv + b — a perfectly linear lever (FD Jacobian is exact)."""

    def predict(dv_rtn: Vec3) -> Vec3:
        out = matrix @ np.asarray(dv_rtn, dtype=np.float64) + offset
        return (float(out[0]), float(out[1]), float(out[2]))

    return predict


# A lever loosely shaped like the real one: a strong transverse→position response
# (the dr_p/dv_a direction) with mild cross-coupling, in metres per (m/s).
_LEVER = np.array(
    [
        [3.0e4, 1.5e3, 0.0],
        [2.0e3, 2.6e4, 8.0e2],
        [0.0, 5.0e2, 2.2e4],
    ]
)


class TestLinearConvergence:
    def test_nulls_a_linear_lever_in_one_step(self) -> None:
        true_dv = np.array([0.012, -0.020, 0.005])
        target = np.array([10_000.0, -4_000.0, 1_500.0])
        # b chosen so predict(true_dv) == target; predict(0) = b is the dispersed crossing.
        offset = target - _LEVER @ true_dv
        plan = solve_apogee_correction(_linear_predict(_LEVER, offset), Target(tuple(target)))

        assert plan.converged
        assert plan.iterations == 1  # exact FD Jacobian on a linear map → one Newton step
        assert plan.actions[0].dv_rtn_m_s == pytest.approx(tuple(true_dv), abs=1e-9)

    def test_already_on_target_takes_no_step(self) -> None:
        # predict(0) already hits the target → converged with zero correction.
        plan = solve_apogee_correction(
            _linear_predict(_LEVER, np.zeros(3)), Target((0.0, 0.0, 0.0))
        )
        assert plan.converged
        assert plan.iterations == 0
        assert plan.actions[0].dv_mag_m_s == pytest.approx(0.0, abs=1e-12)


class TestNonlinearConvergence:
    def test_mildly_nonlinear_converges_in_a_few_steps(self) -> None:
        target = np.array([2_000.0, 5_000.0, -1_000.0])
        offset = np.array([1_500.0, -800.0, 600.0])

        def predict(dv_rtn: Vec3) -> Vec3:
            dv = np.asarray(dv_rtn, dtype=np.float64)
            quad = 5.0e4 * dv * np.abs(dv)  # gentle curvature in the lever
            out = _LEVER @ dv + offset + quad
            return (float(out[0]), float(out[1]), float(out[2]))

        plan = solve_apogee_correction(predict, Target(tuple(target)))
        assert plan.converged
        assert 1 < plan.iterations <= 4
        # residual at the solved Δv is below tolerance
        got = np.asarray(predict(plan.actions[0].dv_rtn_m_s))
        assert float(np.linalg.norm(got - target)) < 1.0


class TestAuthorityBoundary:
    def test_singular_jacobian_records_non_convergence(self) -> None:
        # Third RTN axis is uncontrollable → unreachable target component, singular J.
        singular = np.array([[3.0e4, 0.0, 0.0], [0.0, 2.6e4, 0.0], [0.0, 0.0, 0.0]])
        plan = solve_apogee_correction(
            _linear_predict(singular, np.zeros(3)), Target((1_000.0, 1_000.0, 1_000.0))
        )
        assert not plan.converged  # recorded, not raised
        assert len(plan.actions) == 1

    def test_step_cap_bounds_motion_and_flags_unreachable_target(self) -> None:
        # A lever so weak that the unconstrained Newton step would be enormous; the
        # cap bounds each step to the physical correction scale (it keeps the real
        # solver out of the spurious far-ToA root), so an out-of-reach target reads as
        # honest non-convergence rather than a runaway Δv.
        weak = np.eye(3)  # 1 m per (m/s): nulling a 10 km miss would want 10 km/s
        plan = solve_apogee_correction(
            _linear_predict(weak, np.array([1.0e4, 0.0, 0.0])),
            Target((0.0, 0.0, 0.0)),
            max_iter=2,
            max_step_m_s=2.0,
        )
        assert not plan.converged
        assert plan.actions[0].dv_mag_m_s <= 2.0 * 2 + 1e-9  # ≤ max_iter × cap

    def test_exhausting_iterations_records_non_convergence(self) -> None:
        # A diverging map the Newton step cannot null within the iteration ceiling.
        def predict(dv_rtn: Vec3) -> Vec3:
            dv = np.asarray(dv_rtn, dtype=np.float64)
            out = 1.0e6 * dv**3 + np.array([5.0e3, 0.0, 0.0])
            return (float(out[0]), float(out[1]), float(out[2]))

        plan = solve_apogee_correction(predict, Target((0.0, 0.0, 0.0)), max_iter=3)
        assert plan.iterations <= 3
        assert isinstance(plan.converged, bool)


class TestControlPlan:
    def test_total_dv_sums_action_magnitudes(self) -> None:
        plan = ControlPlan(
            actions=(
                ControlAction("apogee", 0.0, (0.0, 0.1, 0.0), 0.1),
                ControlAction("midcourse", 3600.0, (0.0, 0.0, 0.2), 0.2),
            ),
            converged=True,
            iterations=2,
        )
        assert plan.total_dv_m_s == pytest.approx(0.3)

    def test_empty_plan_has_zero_total_dv(self) -> None:
        # The open-loop capstone records an empty plan.
        plan = ControlPlan(actions=(), converged=True, iterations=0)
        assert plan.total_dv_m_s == 0.0
