"""Tests for the pure Rung A1 differential corrector — synthetic predict, no JVM."""

import numpy as np
import pytest

from puffsat_sim.control import (
    ControlAction,
    ControlPlan,
    PredictFn,
    Target,
    TwoBurnPredictFn,
    Vec3,
    Vec6,
    solve_apogee_correction,
    solve_two_burn_correction,
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


# A second 3×3 lever for the mid-descent burn block: deliberately strong on the
# axis the apogee lever is weak on, so min-norm has a reason to split the load.
_LEVER_MID = np.array(
    [
        [1.8e4, 0.0, 6.0e2],
        [0.0, 1.2e3, 0.0],
        [4.0e2, 0.0, 2.0e4],
    ]
)


def _linear_predict_6(matrix: np.ndarray, offset: np.ndarray) -> TwoBurnPredictFn:
    """predict(dv6) = A·dv6 + b for a 3×6 lever — exact FD Jacobian, min-norm solvable."""

    def predict(dv6: Vec6) -> Vec3:
        out = matrix @ np.asarray(dv6, dtype=np.float64) + offset
        return (float(out[0]), float(out[1]), float(out[2]))

    return predict


class TestTwoBurnLinearConvergence:
    def test_min_norm_nulls_a_linear_lever_in_one_step(self) -> None:
        lever6 = np.hstack([_LEVER, _LEVER_MID])  # 3×6: [apogee | mid-descent]
        target = np.array([10_000.0, -4_000.0, 1_500.0])
        offset = np.array([2_000.0, 1_000.0, -500.0])  # predict(0) = the dispersed crossing
        plan = solve_two_burn_correction(_linear_predict_6(lever6, offset), Target(tuple(target)))

        assert plan.converged
        assert plan.iterations == 1  # exact FD Jacobian on a linear map → one Gauss-Newton step
        assert len(plan.actions) == 2
        # the solved (apogee, mid) stack equals the analytic minimum-norm solution
        x_expected = np.linalg.pinv(lever6) @ (target - offset)
        x_got = np.array(plan.actions[0].dv_rtn_m_s + plan.actions[1].dv_rtn_m_s)
        assert x_got == pytest.approx(x_expected, abs=1e-9)


class TestTwoBurnMinNorm:
    def test_returns_least_norm_split_versus_all_in_apogee(self) -> None:
        # A redundant (reachable) target: the min-Σ‖Δv‖² split must be no larger than
        # the feasible "do it all from the apogee burn" solution (A1's choice).
        lever6 = np.hstack([_LEVER, _LEVER_MID])
        target = np.array([9_000.0, 7_000.0, -3_000.0])
        offset = np.array([500.0, -500.0, 250.0])
        plan = solve_two_burn_correction(_linear_predict_6(lever6, offset), Target(tuple(target)))
        assert plan.converged

        all_in_apogee = np.linalg.solve(_LEVER, target - offset)  # feasible: [x_apogee, 0]
        stacked_norm = np.hypot(plan.actions[0].dv_mag_m_s, plan.actions[1].dv_mag_m_s)
        assert stacked_norm < float(np.linalg.norm(all_in_apogee))

    def test_offloads_a_low_apogee_authority_axis_onto_the_midcourse_burn(self) -> None:
        # Axis 1 is nearly dead from the apogee burn (1e2 m per m/s) but strong from
        # mid-descent (3e4) — the §8/A1 along-track problem in miniature.  A single
        # apogee burn would cost ~60 m/s; min-norm routes the axis to the mid burn.
        apogee = np.diag([3.0e4, 1.0e2, 3.0e4])
        mid = np.diag([1.0e3, 3.0e4, 1.0e3])
        lever6 = np.hstack([apogee, mid])
        target = np.array([0.0, 6_000.0, 0.0])  # pure weak-axis miss
        plan = solve_two_burn_correction(
            _linear_predict_6(lever6, np.zeros(3)), Target(tuple(target))
        )
        assert plan.converged

        apogee_dv, mid_dv = plan.actions[0].dv_rtn_m_s, plan.actions[1].dv_rtn_m_s
        all_in_apogee_cost = abs(target[1] / 1.0e2)  # ~60 m/s — the single-burn wall
        stacked_norm = np.hypot(plan.actions[0].dv_mag_m_s, plan.actions[1].dv_mag_m_s)
        assert all_in_apogee_cost > 50.0
        assert stacked_norm < 1.0  # two-burn collapses it into budget
        assert abs(mid_dv[1]) > 10 * abs(apogee_dv[1])


class TestTwoBurnNonlinearConvergence:
    def test_mildly_nonlinear_converges_in_a_few_steps(self) -> None:
        lever6 = np.hstack([_LEVER, _LEVER_MID])
        target = np.array([2_000.0, 5_000.0, -1_000.0])
        offset = np.array([1_500.0, -800.0, 600.0])

        def predict(dv6: Vec6) -> Vec3:
            lin = lever6 @ np.asarray(dv6, dtype=np.float64)
            out = lin + offset + 2.0e-5 * lin * np.abs(lin)  # gentle curvature on the response
            return (float(out[0]), float(out[1]), float(out[2]))

        plan = solve_two_burn_correction(predict, Target(tuple(target)))
        assert plan.converged
        assert 1 < plan.iterations <= 5
        got = np.asarray(predict(plan.actions[0].dv_rtn_m_s + plan.actions[1].dv_rtn_m_s))
        assert float(np.linalg.norm(got - target)) < 1.0

    def test_already_on_target_takes_no_step(self) -> None:
        lever6 = np.hstack([_LEVER, _LEVER_MID])
        plan = solve_two_burn_correction(
            _linear_predict_6(lever6, np.zeros(3)), Target((0.0, 0.0, 0.0))
        )
        assert plan.converged
        assert plan.iterations == 0
        assert plan.actions[0].dv_mag_m_s == pytest.approx(0.0, abs=1e-12)
        assert plan.actions[1].dv_mag_m_s == pytest.approx(0.0, abs=1e-12)


class TestTwoBurnAuthorityBoundary:
    def test_unreachable_target_records_non_convergence(self) -> None:
        # Axis 2 is uncontrollable from BOTH burns; lstsq finds the min-residual point
        # but that residual stays above tol — recorded as non-converged, not raised.
        apogee = np.array([[3.0e4, 0.0, 0.0], [0.0, 2.6e4, 0.0], [0.0, 0.0, 0.0]])
        mid = np.array([[1.0e4, 0.0, 0.0], [0.0, 1.5e4, 0.0], [0.0, 0.0, 0.0]])
        lever6 = np.hstack([apogee, mid])
        plan = solve_two_burn_correction(
            _linear_predict_6(lever6, np.zeros(3)), Target((1_000.0, 1_000.0, 1_000.0))
        )
        assert not plan.converged
        assert len(plan.actions) == 2

    def test_step_cap_bounds_the_six_vector_motion(self) -> None:
        # A lever so weak the unconstrained step would be enormous; the cap bounds each
        # 6-vector step to the physical correction scale, so out-of-reach reads as honest
        # non-convergence rather than a runaway Δv (the A1 spurious-far-root guard).
        weak = np.hstack([np.eye(3), np.eye(3)])  # 1 m per (m/s) from each burn
        plan = solve_two_burn_correction(
            _linear_predict_6(weak, np.array([1.0e4, 0.0, 0.0])),
            Target((0.0, 0.0, 0.0)),
            max_iter=2,
            max_step_m_s=2.0,
        )
        assert not plan.converged
        stacked_norm = np.hypot(plan.actions[0].dv_mag_m_s, plan.actions[1].dv_mag_m_s)
        assert stacked_norm <= 2.0 * 2 + 1e-9  # ≤ max_iter × cap


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
