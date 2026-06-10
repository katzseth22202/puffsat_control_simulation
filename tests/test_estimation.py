"""Tests for the pure C1 estimation core (ADR 0012) — sigma points, unscented
transform, UKF predict/update, filter dynamics, measurement models, LinCov
recursion, and NEES consistency bounds (no JVM)."""

import numpy as np
import pytest

from puffsat_sim.constants import WGS84_MU
from puffsat_sim.estimation import (
    FilterState,
    LincovEpoch,
    MeasurementModel,
    NodeState,
    UnscentedSpec,
    average_nees_bounds,
    gnss_position_fix,
    nees,
    run_lincov,
    los_velocity_to_node,
    merwe_sigma_points,
    range_to_node,
    process_noise_white_accel,
    two_body_j2_flow,
    ukf_predict,
    ukf_update,
    unscented_transform,
)
from puffsat_sim.forces.geopotential import j2_nodal_regression_rate
from puffsat_sim.orbital_math import keplerian_period


def _circular_state(radius_m: float) -> np.ndarray:
    speed = np.sqrt(WGS84_MU / radius_m)
    return np.array([radius_m, 0.0, 0.0, 0.0, speed, 0.0])


class TestUnscentedTransform:
    def test_sigma_points_round_trip_mean_and_covariance(self) -> None:
        x = np.array([1.0, -2.0, 0.5])
        cov = np.array([[4.0, 1.0, 0.0], [1.0, 3.0, 0.5], [0.0, 0.5, 2.0]])

        points = merwe_sigma_points(x, cov, UnscentedSpec())
        mean, recovered = unscented_transform(points)

        np.testing.assert_allclose(mean, x, rtol=0, atol=1e-9)
        np.testing.assert_allclose(recovered, cov, rtol=1e-9, atol=1e-9)


class TestUkfPredict:
    def test_linear_flow_matches_analytic_propagation(self) -> None:
        a = np.array([[1.0, 0.5], [0.0, 1.0]])
        q = np.diag([0.01, 0.02])
        state = FilterState(x=np.array([1.0, 2.0]), cov=np.array([[2.0, 0.3], [0.3, 1.0]]))

        predicted = ukf_predict(state, lambda x: a @ x, q, UnscentedSpec())

        np.testing.assert_allclose(predicted.x, a @ state.x, rtol=0, atol=1e-12)
        np.testing.assert_allclose(predicted.cov, a @ state.cov @ a.T + q, rtol=1e-12, atol=1e-12)


class TestUkfUpdate:
    def test_linear_measurement_matches_analytic_kalman_update(self) -> None:
        h = np.array([[1.0, 0.0]])
        r = np.array([[0.25]])
        state = FilterState(x=np.array([1.0, 2.0]), cov=np.array([[2.0, 0.3], [0.3, 1.0]]))
        z = np.array([1.7])

        updated = ukf_update(state, z, lambda x: h @ x, r, UnscentedSpec())

        gain = state.cov @ h.T @ np.linalg.inv(h @ state.cov @ h.T + r)
        expected_x = state.x + gain @ (z - h @ state.x)
        expected_cov = (np.eye(2) - gain @ h) @ state.cov
        np.testing.assert_allclose(updated.x, expected_x, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(updated.cov, expected_cov, rtol=1e-12, atol=1e-12)


class TestTwoBodyJ2Flow:
    def test_two_body_circular_orbit_closes_after_one_period(self) -> None:
        radius = 7_000_000.0
        x0 = _circular_state(radius)
        period = keplerian_period(radius)

        x1 = two_body_j2_flow(x0, period, j2=0.0)

        np.testing.assert_allclose(x1[:3], x0[:3], rtol=0, atol=1.0)
        np.testing.assert_allclose(x1[3:], x0[3:], rtol=0, atol=1e-3)

    def test_inclined_orbit_node_regresses_at_analytic_j2_rate(self) -> None:
        radius = 7_000_000.0
        inclination = np.deg2rad(60.0)
        speed = np.sqrt(WGS84_MU / radius)
        x0 = np.array(
            [radius, 0.0, 0.0, 0.0, speed * np.cos(inclination), speed * np.sin(inclination)]
        )
        n_orbits = 20
        elapsed = n_orbits * keplerian_period(radius)

        x1 = two_body_j2_flow(x0, elapsed)

        h = np.cross(x1[:3], x1[3:])
        node = np.cross([0.0, 0.0, 1.0], h)
        measured_raan = np.arctan2(node[1], node[0])
        expected_raan = j2_nodal_regression_rate(radius, 0.0, inclination) * elapsed
        np.testing.assert_allclose(measured_raan, expected_raan, rtol=0.02)

    def test_components_beyond_six_pass_through_unchanged(self) -> None:
        x0 = np.concatenate([_circular_state(7_000_000.0), [0.04, 0.02]])

        x1 = two_body_j2_flow(x0, 600.0)

        assert x1.shape == (8,)
        np.testing.assert_array_equal(x1[6:], x0[6:])


class TestProcessNoise:
    def test_white_acceleration_q_has_analytic_blocks(self) -> None:
        q_accel = 1e-7
        dt = 100.0

        q = process_noise_white_accel(q_accel, dt)

        var = q_accel**2
        np.testing.assert_allclose(np.diag(q)[:3], np.full(3, var * dt**3 / 3.0))
        np.testing.assert_allclose(np.diag(q)[3:], np.full(3, var * dt))
        np.testing.assert_allclose(np.diag(q[:3, 3:]), np.full(3, var * dt**2 / 2.0))
        np.testing.assert_allclose(q, q.T)

    def test_extra_states_get_zero_process_noise(self) -> None:
        q = process_noise_white_accel(1e-7, 100.0, n_states=8)

        assert q.shape == (8, 8)
        np.testing.assert_array_equal(q[6:, :], np.zeros((2, 8)))


class TestMeasurementModels:
    def test_range_and_los_velocity_follow_the_geometry(self) -> None:
        node = NodeState(position_m=np.array([1e6, 0.0, 0.0]), velocity_m_s=np.zeros(3))
        x = np.array([4e6, 4e6, 0.0, 0.0, 25.0, 0.0])

        assert range_to_node(x, node) == pytest.approx(5e6)
        assert los_velocity_to_node(x, node) == pytest.approx(0.8 * 25.0)

    def test_los_velocity_is_blind_to_transverse_motion(self) -> None:
        node = NodeState(position_m=np.zeros(3), velocity_m_s=np.zeros(3))
        x = np.array([1e6, 0.0, 0.0, 0.0, 100.0, 0.0])

        assert los_velocity_to_node(x, node) == pytest.approx(0.0)

    def test_gnss_position_fix_returns_the_position_components(self) -> None:
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

        np.testing.assert_array_equal(gnss_position_fix(x), x[:3])


class TestNees:
    def test_nees_is_the_squared_mahalanobis_distance(self) -> None:
        cov = np.diag([4.0, 1.0])
        error = np.array([2.0, 3.0])

        assert nees(error, cov) == pytest.approx(1.0 + 9.0)

    def test_average_nees_bounds_match_the_chi_square_table(self) -> None:
        lo, hi = average_nees_bounds(dim=2, n_samples=50)
        assert lo == pytest.approx(74.2219 / 50, rel=5e-3)
        assert hi == pytest.approx(129.561 / 50, rel=5e-3)

        lo6, hi6 = average_nees_bounds(dim=6, n_samples=1)
        assert lo6 == pytest.approx(1.2373, rel=3e-2)
        assert hi6 == pytest.approx(14.4494, rel=3e-2)


class TestRunLincov:
    def test_linear_system_matches_analytic_kalman_recursion_with_pinned_mean(self) -> None:
        a = np.array([[1.0, 1.0], [0.0, 1.0]])
        q = np.diag([1e-4, 1e-4])
        h = np.array([[1.0, 0.0]])
        r = np.array([[0.5]])
        initial = FilterState(x=np.array([0.0, 1.0]), cov=np.diag([4.0, 4.0]))
        epochs = tuple(
            LincovEpoch(
                dt_s=1.0,
                process_noise=q,
                measurements=(MeasurementModel(h=lambda x: h @ x, noise_cov=r),),
            )
            for _ in range(5)
        )

        states = run_lincov(initial, lambda x, dt: a @ x, epochs, UnscentedSpec())

        x_ref, cov = initial.x, initial.cov
        for _ in range(5):
            x_ref = a @ x_ref
            cov = a @ cov @ a.T + q
            gain = cov @ h.T @ np.linalg.inv(h @ cov @ h.T + r)
            cov = (np.eye(2) - gain @ h) @ cov
        assert len(states) == 5
        np.testing.assert_allclose(states[-1].x, x_ref, rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(states[-1].cov, cov, rtol=1e-9, atol=1e-12)


def _lincov_with_node_offsets(offsets_m: list[np.ndarray]) -> FilterState:
    """LinCov over a slow high arc with range+Doppler from nodes at given offsets."""
    radius = 1.5e8
    x0 = _circular_state(radius)
    initial = FilterState(x=x0, cov=np.diag([1e4**2] * 3 + [1.0**2] * 3))
    dt = 33.0
    spec = UnscentedSpec()

    def epoch_at(k: int) -> LincovEpoch:
        x_ref = two_body_j2_flow(x0, (k + 1) * dt)
        measurements: list[MeasurementModel] = []
        for offset in offsets_m:
            node = NodeState(position_m=x_ref[:3] + offset, velocity_m_s=x_ref[3:6])
            measurements.append(
                MeasurementModel(
                    h=lambda x, n=node: np.array([range_to_node(x, n)]),
                    noise_cov=np.array([[1.0**2]]),
                )
            )
            measurements.append(
                MeasurementModel(
                    h=lambda x, n=node: np.array([los_velocity_to_node(x, n)]),
                    noise_cov=np.array([[1e-3**2]]),
                )
            )
        return LincovEpoch(
            dt_s=dt,
            process_noise=process_noise_white_accel(5e-8, dt),
            measurements=tuple(measurements),
        )

    epochs = tuple(epoch_at(k) for k in range(30))
    return run_lincov(initial, lambda x, dt_s: two_body_j2_flow(x, dt_s), epochs, spec)[-1]


class TestApogeeObservability:
    def test_los_diverse_nodes_pin_velocity_below_the_c0_threshold(self) -> None:
        spread = [
            np.array([3e6, 0.0, 0.0]),
            np.array([0.0, 3e6, 0.0]),
            np.array([0.0, 0.0, 3e6]),
            np.array([-2e6, -2e6, -2e6]),
        ]

        final = _lincov_with_node_offsets(spread)

        velocity_sigma = np.sqrt(np.diag(final.cov)[3:6])
        assert velocity_sigma.max() < 2.3e-2

    def test_collinear_nodes_leave_perpendicular_velocity_unobserved(self) -> None:
        collinear = [np.array([d, 0.0, 0.0]) for d in (3e6, 4e6, 5e6, 6e6)]
        spread = [
            np.array([3e6, 0.0, 0.0]),
            np.array([0.0, 3e6, 0.0]),
            np.array([0.0, 0.0, 3e6]),
            np.array([-2e6, -2e6, -2e6]),
        ]

        blind = _lincov_with_node_offsets(collinear)
        diverse = _lincov_with_node_offsets(spread)

        blind_perp = np.sqrt(np.diag(blind.cov)[4:6])
        diverse_perp = np.sqrt(np.diag(diverse.cov)[4:6])
        assert blind_perp.min() > 5.0 * diverse_perp.max()
