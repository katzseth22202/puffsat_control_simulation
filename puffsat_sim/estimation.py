"""Pure C1 estimation core — UKF, measurement models, LinCov, NEES (no JVM).

C1 (design doc §13, ADR 0012) asks whether a realistic measurement suite —
two-way range + carrier Doppler from co-flying coordinator nodes, GNSS near
perigee — can know the apogee state well enough for the C0 threshold
(``Φ Σ Φᵀ`` inside the terminal catch radius).  This module is the pure side:
the owned, fully-typed UKF (sigma points, unscented transform, predict/update),
the onboard two-body+J2 filter dynamics, the measurement models ``h(x)`` with
known-ephemeris node states as inputs, the LinCov covariance recursion (the
sweep engine — measurement-value-independent), and the NEES consistency bounds
the seeded truth runs are judged by.  The JVM ``run_nav_feasibility`` (in
:mod:`puffsat_sim.montecarlo`) consumes it.

State convention: ``x[:3]`` position [m], ``x[3:6]`` velocity [m/s] in EME2000;
components beyond 6 (C2's lumped coefficients) pass through the dynamics
unchanged.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.constants import EARTH_RADIUS_M, J2, WGS84_MU

StateFlow = Callable[[NDArray[np.float64]], NDArray[np.float64]]
TimedFlow = Callable[[NDArray[np.float64], float], NDArray[np.float64]]
MeasurementFn = Callable[[NDArray[np.float64]], NDArray[np.float64]]


@dataclass(frozen=True)
class FilterState:
    """A filter estimate: state mean ``x`` and error covariance ``cov``."""

    x: NDArray[np.float64]
    cov: NDArray[np.float64]


@dataclass(frozen=True)
class UnscentedSpec:
    """Merwe scaled sigma-point parameters (ADR 0012 decision 5).

    Defaults use ``alpha=1`` (full-spread points), not the canonical ``1e-3``:
    the tiny-alpha convention makes the center weight ~``-1e6`` and recovers the
    covariance only through large-weight cancellation — needless conditioning
    risk at orbital scales.  ``beta=2`` is Gaussian-optimal; ``kappa=0`` keeps
    the spread ``√n``-scaled for any state dimension.
    """

    alpha: float = 1.0
    beta: float = 2.0
    kappa: float = 0.0


@dataclass(frozen=True)
class SigmaPoints:
    """A sigma-point set: ``points[i]`` rows with mean/covariance weights."""

    points: NDArray[np.float64]
    w_mean: NDArray[np.float64]
    w_cov: NDArray[np.float64]


def merwe_sigma_points(
    x: NDArray[np.float64], cov: NDArray[np.float64], spec: UnscentedSpec
) -> SigmaPoints:
    """The 2n+1 Merwe scaled sigma points of a Gaussian ``(x, cov)``."""
    mean = np.asarray(x, dtype=np.float64)
    n = mean.shape[0]
    lam = spec.alpha**2 * (n + spec.kappa) - n
    spread = np.linalg.cholesky((n + lam) * np.asarray(cov, dtype=np.float64))

    points = np.empty((2 * n + 1, n), dtype=np.float64)
    points[0] = mean
    points[1 : n + 1] = mean + spread.T
    points[n + 1 :] = mean - spread.T

    w_mean = np.full(2 * n + 1, 1.0 / (2.0 * (n + lam)), dtype=np.float64)
    w_cov = w_mean.copy()
    w_mean[0] = lam / (n + lam)
    w_cov[0] = lam / (n + lam) + (1.0 - spec.alpha**2 + spec.beta)
    return SigmaPoints(points=points, w_mean=w_mean, w_cov=w_cov)


def unscented_transform(
    sigma: SigmaPoints,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Weighted mean and covariance of a sigma-point set."""
    mean = sigma.w_mean @ sigma.points
    centered = sigma.points - mean
    cov = (sigma.w_cov * centered.T) @ centered
    return np.asarray(mean, dtype=np.float64), np.asarray(cov, dtype=np.float64)


def ukf_predict(
    state: FilterState,
    flow: StateFlow,
    process_noise: NDArray[np.float64],
    spec: UnscentedSpec,
) -> FilterState:
    """UKF time update: push the sigma points through ``flow``, add process noise.

    ``flow`` maps a state vector across the prediction interval (dt is the
    caller's, baked into the closure — for C1 a ``two_body_j2_flow`` partial);
    ``process_noise`` is the discrete Q for that same interval.
    """
    sigma = merwe_sigma_points(state.x, state.cov, spec)
    propagated = np.apply_along_axis(flow, 1, sigma.points)
    mean, cov = unscented_transform(
        SigmaPoints(points=propagated, w_mean=sigma.w_mean, w_cov=sigma.w_cov)
    )
    return FilterState(x=mean, cov=cov + np.asarray(process_noise, dtype=np.float64))


def two_body_j2_acceleration(
    position_m: NDArray[np.float64], j2: float = J2
) -> NDArray[np.float64]:
    """Central-body + J2 acceleration [m/s²] at an EME2000 position (the onboard model).

    The cheap flight-filter dynamics of ADR 0012 decision 6 — everything the
    truth model adds beyond this (third body, SRP, drag) is the model gap the
    process noise Q must absorb at C1 (and the estimated coefficients at C2).
    """
    r = float(np.linalg.norm(position_m))
    central = -WGS84_MU / r**3 * position_m

    z_over_r_sq = (position_m[2] / r) ** 2
    factor = -1.5 * j2 * WGS84_MU * EARTH_RADIUS_M**2 / r**5
    oblateness = factor * position_m * (1.0 - 5.0 * z_over_r_sq)
    oblateness[2] = factor * position_m[2] * (3.0 - 5.0 * z_over_r_sq)
    return np.asarray(central + oblateness, dtype=np.float64)


def _state_derivative(x: NDArray[np.float64], j2: float) -> NDArray[np.float64]:
    dx = np.zeros_like(x)
    dx[:3] = x[3:6]
    dx[3:6] = two_body_j2_acceleration(x[:3], j2)
    return dx


def two_body_j2_flow(
    x: NDArray[np.float64],
    dt_s: float,
    *,
    max_step_s: float = 10.0,
    j2: float = J2,
) -> NDArray[np.float64]:
    """Propagate a state across ``dt_s`` under two-body + J2 with RK4 substeps."""
    state = np.asarray(x, dtype=np.float64).copy()
    steps = max(1, int(np.ceil(abs(dt_s) / max_step_s)))
    h = dt_s / steps
    for _ in range(steps):
        k1 = _state_derivative(state, j2)
        k2 = _state_derivative(state + 0.5 * h * k1, j2)
        k3 = _state_derivative(state + 0.5 * h * k2, j2)
        k4 = _state_derivative(state + h * k3, j2)
        state = state + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return state


def nees(error: NDArray[np.float64], cov: NDArray[np.float64]) -> float:
    """Normalized estimation error squared ``eᵀ P⁻¹ e`` — the filter-consistency statistic.

    For a consistent filter the error is N(0, P), so NEES is χ²(dim)-distributed;
    the seeded truth runs of ADR 0012 decision 7 are judged by its average
    against :func:`average_nees_bounds`.
    """
    e = np.asarray(error, dtype=np.float64)
    return float(e @ np.linalg.solve(np.asarray(cov, dtype=np.float64), e))


def _normal_quantile(p: float) -> float:
    """Standard-normal quantile via bisection on the erf-based CDF (no scipy)."""
    lo, hi = -10.0, 10.0
    for _ in range(80):
        mid = (lo + hi) / 2.0
        if 0.5 * (1.0 + math.erf(mid / math.sqrt(2.0))) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def _chi_square_quantile(p: float, dof: int) -> float:
    """χ² quantile via the Wilson–Hilferty cube approximation.

    Within ~0.1% for dof ≳ 30 (the NEES-gate regime, dof = runs × dim); the
    worst case at dof ~ 6 is ~3% in the lower tail — immaterial for a
    consistency gate and documented in the tests' tolerances.
    """
    z = _normal_quantile(p)
    spread = 2.0 / (9.0 * dof)
    return dof * (1.0 - spread + z * math.sqrt(spread)) ** 3


def average_nees_bounds(dim: int, n_samples: int, confidence: float = 0.95) -> tuple[float, float]:
    """The two-sided acceptance interval for the *average* NEES of ``n_samples`` draws.

    ``n_samples · NEES_avg ~ χ²(n_samples · dim)`` for a consistent filter; an
    average outside these bounds means the claimed covariance is fiction (too
    small → optimistic filter, the dangerous direction; too large → pessimistic).
    """
    alpha = 1.0 - confidence
    dof = n_samples * dim
    return (
        _chi_square_quantile(alpha / 2.0, dof) / n_samples,
        _chi_square_quantile(1.0 - alpha / 2.0, dof) / n_samples,
    )


@dataclass(frozen=True)
class NodeState:
    """A coordinator node's known ephemeris at a measurement epoch (ADR 0012 decision 1).

    An *input* to the measurement models, never a filter state — the node's own
    nav error is carried as an inflation term in the measurement noise R.
    """

    position_m: NDArray[np.float64]
    velocity_m_s: NDArray[np.float64]


def range_to_node(x: NDArray[np.float64], node: NodeState) -> float:
    """Geometric range [m] from the PuffSat state to a node (two-way range observable)."""
    return float(np.linalg.norm(x[:3] - node.position_m))


def los_velocity_to_node(x: NDArray[np.float64], node: NodeState) -> float:
    """Range-rate [m/s] along the node line of sight (two-way carrier-Doppler observable).

    Positive = receding.  The projection is what makes node LOS *diversity* the
    velocity-observability question of ADR 0012 decision 2 — a single LOS is
    blind to the perpendicular velocity components.
    """
    offset = x[:3] - node.position_m
    los = offset / np.linalg.norm(offset)
    return float((x[3:6] - node.velocity_m_s) @ los)


def gnss_position_fix(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """The GNSS observable: the position components, directly (ADR 0012 decision 4)."""
    return np.asarray(x[:3], dtype=np.float64)


def process_noise_white_accel(
    q_accel_m_s2: float, dt_s: float, n_states: int = 6
) -> NDArray[np.float64]:
    """Discrete Q for a white unmodeled acceleration of 1σ ``q_accel_m_s2`` over ``dt_s``.

    The classic continuous-white-noise-acceleration discretization (per axis:
    ``[[dt³/3, dt²/2], [dt²/2, dt]]·q²``) — the explicitly swept knob of ADR 0012
    decision 6 that absorbs the truth−filter model gap.  States beyond the 6
    position/velocity components (C2's coefficients) get zero process noise here;
    their random-walk q is C2's to add.
    """
    var = q_accel_m_s2**2
    q = np.zeros((n_states, n_states), dtype=np.float64)
    eye3 = np.eye(3, dtype=np.float64)
    q[:3, :3] = var * dt_s**3 / 3.0 * eye3
    q[:3, 3:6] = var * dt_s**2 / 2.0 * eye3
    q[3:6, :3] = var * dt_s**2 / 2.0 * eye3
    q[3:6, 3:6] = var * dt_s * eye3
    return q


def ukf_update(
    state: FilterState,
    z: NDArray[np.float64],
    h: MeasurementFn,
    noise_cov: NDArray[np.float64],
    spec: UnscentedSpec,
) -> FilterState:
    """UKF measurement update: fuse measurement ``z`` with model ``h`` and noise ``R``.

    The posterior covariance is symmetrized (``(P + Pᵀ)/2``) — the Kalman algebra
    is symmetric in exact arithmetic, and re-imposing it each update is the cheap
    half of the conditioning defenses ADR 0012 decision 5 reserves (square-root
    form stays the escalation if validation runs ever break Cholesky).
    """
    sigma = merwe_sigma_points(state.x, state.cov, spec)
    predicted = np.apply_along_axis(h, 1, sigma.points)

    z_mean = sigma.w_mean @ predicted
    z_centered = predicted - z_mean
    innovation_cov = (sigma.w_cov * z_centered.T) @ z_centered + np.asarray(
        noise_cov, dtype=np.float64
    )
    x_centered = sigma.points - state.x
    cross_cov = (sigma.w_cov * x_centered.T) @ z_centered

    gain = np.linalg.solve(innovation_cov, cross_cov.T).T
    x_post = state.x + gain @ (np.asarray(z, dtype=np.float64) - z_mean)
    cov_post = state.cov - gain @ innovation_cov @ gain.T
    return FilterState(x=x_post, cov=(cov_post + cov_post.T) / 2.0)


@dataclass(frozen=True)
class MeasurementModel:
    """A measurement type at an epoch: the observable ``h(x)`` and its noise R."""

    h: MeasurementFn
    noise_cov: NDArray[np.float64]


@dataclass(frozen=True)
class LincovEpoch:
    """One LinCov step: propagate ``dt_s`` (with that interval's Q), then update."""

    dt_s: float
    process_noise: NDArray[np.float64]
    measurements: tuple[MeasurementModel, ...]


def _flow_over(flow: TimedFlow, dt_s: float) -> StateFlow:
    def step(x: NDArray[np.float64]) -> NDArray[np.float64]:
        return flow(x, dt_s)

    return step


def run_lincov(
    initial: FilterState,
    flow: TimedFlow,
    epochs: tuple[LincovEpoch, ...],
    spec: UnscentedSpec,
) -> tuple[FilterState, ...]:
    """The LinCov covariance recursion along a reference trajectory (ADR 0012 decision 7).

    The Kalman covariance recursion does not depend on measurement *values* —
    each update is fed ``z = h(x̂)`` (zero innovation), so the mean stays pinned
    to the reference while Σ evolves exactly as the filter's would.  This is the
    C1 sweep engine: one deterministic pass per sweep cell, no truth propagation,
    reusing the same UKF machinery the seeded validation runs exercise.
    """
    states: list[FilterState] = []
    state = initial
    for epoch in epochs:
        state = ukf_predict(state, _flow_over(flow, epoch.dt_s), epoch.process_noise, spec)
        for measurement in epoch.measurements:
            state = ukf_update(
                state, measurement.h(state.x), measurement.h, measurement.noise_cov, spec
            )
        states.append(state)
    return tuple(states)
