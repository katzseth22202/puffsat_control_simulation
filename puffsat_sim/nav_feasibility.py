"""Pure C1 nav-feasibility sweep harness — the LinCov requirement envelope (no JVM).

C1 (design doc §13, ADR 0012 decision 7) measures which sensing architectures meet
C0's navigation requirement: for each sweep cell — (range σ, Doppler σ with a
range-only point, cadence, node geometry, process-noise q) varied one axis at a
time around a nominal, the A3/C0 grid discipline — run the LinCov covariance
recursion (:func:`puffsat_sim.estimation.run_lincov`) along the coast arc ending
at the apogee correction node, rotate Σ into apogee-RTN, and threshold the
``Φ Σ Φᵀ`` lateral interception miss against the terminal catch radius through
C0's measured Φ (:func:`puffsat_sim.navigation.induced_miss_covariance`).

Everything here is deterministic and JVM-free: the reference trajectory is the
*filter's own* two-body+J2 dynamics (the onboard model — ADR 0012 decision 6),
and node ephemerides are rigid RTN-frame offsets co-flying with the reference
(known-ephemeris beacons, decision 1).  The seeded truth runs that *validate*
this envelope (NEES, decision 7 layer 2) live behind the JVM seam in
:mod:`puffsat_sim.montecarlo`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.constants import EARTH_RADIUS_M, WGS84_MU
from puffsat_sim.dispersion import rtn_basis, rtn_components
from puffsat_sim.estimation import (
    FilterState,
    LincovEpoch,
    MeasurementModel,
    NodeState,
    UnscentedSpec,
    average_nees_bounds,
    los_velocity_to_node,
    nees,
    process_noise_white_accel,
    range_to_node,
    run_lincov,
    two_body_j2_flow,
    ukf_predict,
    ukf_update,
)
from puffsat_sim.mission import APOGEE_ALT_M, PERIGEE_ALT_M
from puffsat_sim.navigation import induced_miss_covariance
from puffsat_sim.orbital_math import keplerian_elements

_DEFAULT_UNSCENTED = UnscentedSpec()


@dataclass(frozen=True)
class NavFeasibilityCell:
    """One sweep cell: a fully-resolved sensing architecture to evaluate.

    ``axis`` names the knob this cell varies from the nominal (``"nominal"`` for
    the center cell); ``doppler_sigma_m_s=None`` is the range-only architecture
    (ADR 0012 decision 3 — "is Doppler load-bearing?" is measured, not assumed).
    """

    cell_index: int
    axis: str
    range_sigma_m: float
    doppler_sigma_m_s: float | None
    cadence_hz: float
    cone_half_angle_rad: float
    n_nodes: int
    q_accel_m_s2: float


@dataclass(frozen=True)
class NavFeasibilitySpec:
    """The C1 one-axis-at-a-time sweep around a nominal architecture (ADR 0012).

    Nominal values sit at the grilled leans: ~1 m two-way range, ~1 mm/s carrier
    Doppler (OCXO-bound link budget, decision 3), 0.03 Hz coast cadence (§16.4),
    4 nodes at 45° LOS half-angle on ~3000 km baselines.  Sweep brackets run each
    knob from comfortably-better to clearly-degraded so the envelope edge is
    inside the swept range.

    Nominal q is the **NEES-validated** value (measured 2026-06-10, ADR 0012
    findings): the truth−filter gap at apogee is dominated by the third-body
    tidal acceleration (~3e-5 m/s², Moon+Sun — three orders above SRP), and the
    seeded truth runs test honest only at q ≈ 1e-4.  The original SRP-scale
    guess (5e-8) stays as a swept point to document the trap: its envelope is
    fiction (claimed 6 µm/s vs actual 0.09 m/s — above the C0 requirement), and
    3e-5 marks the measured consistency crossover.
    """

    nominal_range_sigma_m: float = 1.0
    nominal_doppler_sigma_m_s: float | None = 1e-3
    nominal_cadence_hz: float = 0.03
    nominal_cone_half_angle_rad: float = math.radians(45.0)
    nominal_n_nodes: int = 4
    nominal_q_accel_m_s2: float = 1e-4

    range_sigma_values_m: tuple[float, ...] = (0.01, 100.0)
    doppler_sigma_values_m_s: tuple[float | None, ...] = (1e-4, 1e-2, None)
    cadence_values_hz: tuple[float, ...] = (0.003, 0.3)
    cone_half_angle_values_rad: tuple[float, ...] = (
        math.radians(5.0),
        math.radians(15.0),
        math.radians(90.0),
    )
    n_nodes_values: tuple[int, ...] = (3, 6)
    q_accel_values_m_s2: tuple[float, ...] = (5e-8, 3e-5)

    baseline_m: float = 3e6
    arc_duration_s: float = 43_200.0
    prior_pos_sigma_m: float = 1e4
    prior_vel_sigma_m_s: float = 1.0
    catch_radius_m: float = 5e3
    flow_max_step_s: float = 60.0


def nav_feasibility_cells(spec: NavFeasibilitySpec) -> tuple[NavFeasibilityCell, ...]:
    """Enumerate the nominal cell plus the one-axis-at-a-time variations."""

    def cell(index: int, axis: str, **overrides: object) -> NavFeasibilityCell:
        knobs: dict[str, object] = {
            "range_sigma_m": spec.nominal_range_sigma_m,
            "doppler_sigma_m_s": spec.nominal_doppler_sigma_m_s,
            "cadence_hz": spec.nominal_cadence_hz,
            "cone_half_angle_rad": spec.nominal_cone_half_angle_rad,
            "n_nodes": spec.nominal_n_nodes,
            "q_accel_m_s2": spec.nominal_q_accel_m_s2,
        }
        knobs.update(overrides)
        return NavFeasibilityCell(cell_index=index, axis=axis, **knobs)  # type: ignore[arg-type]

    cells = [cell(0, "nominal")]
    swept: tuple[tuple[str, tuple[object, ...]], ...] = (
        ("range_sigma_m", spec.range_sigma_values_m),
        ("doppler_sigma_m_s", spec.doppler_sigma_values_m_s),
        ("cadence_hz", spec.cadence_values_hz),
        ("cone_half_angle_rad", spec.cone_half_angle_values_rad),
        ("n_nodes", spec.n_nodes_values),
        ("q_accel_m_s2", spec.q_accel_values_m_s2),
    )
    index = 1
    for axis, values in swept:
        for value in values:
            cells.append(cell(index, axis, **{axis: value}))
            index += 1
    return tuple(cells)


def apogee_state() -> NDArray[np.float64]:
    """The reference-orbit apogee state (EME2000 Cartesian, equatorial convention).

    Position along +x̂, velocity along +ŷ (so RTN at this state is axis-aligned);
    speed from vis-viva on the mission orbit.  The C1 sweep's Σ is evaluated
    here — the apogee correction node where C0's Φ applies.
    """
    r_a = EARTH_RADIUS_M + APOGEE_ALT_M
    a, _ = keplerian_elements(PERIGEE_ALT_M, APOGEE_ALT_M)
    v_a = math.sqrt(WGS84_MU * (2.0 / r_a - 1.0 / a))
    return np.array([r_a, 0.0, 0.0, 0.0, v_a, 0.0], dtype=np.float64)


def covariance_to_rtn(cov6: NDArray[np.float64], state: NDArray[np.float64]) -> NDArray[np.float64]:
    """Rotate a 6×6 EME2000 state covariance into the RTN frame of ``state``.

    The bridge from the filter's Cartesian Σ to C0's apogee-RTN Φ basis: both
    position and velocity blocks rotate by the same RTN basis (an instantaneous
    frame change, no transport terms — the same convention C0's nav offsets use).
    """
    basis = np.array(
        rtn_basis(
            (float(state[0]), float(state[1]), float(state[2])),
            (float(state[3]), float(state[4]), float(state[5])),
        ),
        dtype=np.float64,
    )
    rotation = np.zeros((6, 6), dtype=np.float64)
    rotation[:3, :3] = basis
    rotation[3:, 3:] = basis
    return np.asarray(rotation @ np.asarray(cov6, dtype=np.float64) @ rotation.T, dtype=np.float64)


@dataclass(frozen=True)
class NavFeasibilityOutcome:
    """One evaluated cell: the achieved apogee knowledge and its C0 verdict.

    ``vel_sigma_rtn_m_s[1]`` (transverse velocity) is the C0-binding axis;
    ``lateral_miss_1sigma_m`` is the RMS lateral interception miss ``√tr(ΦΣΦᵀ)_TN``
    the achieved Σ induces through C0's measured Φ.
    """

    cell: NavFeasibilityCell
    pos_sigma_rtn_m: tuple[float, float, float]
    vel_sigma_rtn_m_s: tuple[float, float, float]
    lateral_miss_1sigma_m: float
    meets_catch_radius: bool


def _node_states_at(
    x_ref: NDArray[np.float64], cell: NavFeasibilityCell, spec: NavFeasibilitySpec
) -> tuple[NodeState, ...]:
    basis = np.array(
        rtn_basis(
            (float(x_ref[0]), float(x_ref[1]), float(x_ref[2])),
            (float(x_ref[3]), float(x_ref[4]), float(x_ref[5])),
        ),
        dtype=np.float64,
    )
    nodes = []
    for direction in node_directions_rtn(cell.n_nodes, cell.cone_half_angle_rad):
        offset = spec.baseline_m * (direction @ basis)
        nodes.append(NodeState(position_m=x_ref[:3] + offset, velocity_m_s=x_ref[3:6]))
    return tuple(nodes)


def _epoch_measurements(
    nodes: tuple[NodeState, ...], cell: NavFeasibilityCell
) -> tuple[MeasurementModel, ...]:
    measurements: list[MeasurementModel] = []
    for node in nodes:
        measurements.append(
            MeasurementModel(
                h=lambda x, n=node: np.array([range_to_node(x, n)]),  # type: ignore[misc]
                noise_cov=np.array([[cell.range_sigma_m**2]]),
            )
        )
        if cell.doppler_sigma_m_s is not None:
            measurements.append(
                MeasurementModel(
                    h=lambda x, n=node: np.array([los_velocity_to_node(x, n)]),  # type: ignore[misc]
                    noise_cov=np.array([[cell.doppler_sigma_m_s**2]]),
                )
            )
    return tuple(measurements)


def evaluate_cell(
    cell: NavFeasibilityCell,
    spec: NavFeasibilitySpec,
    phi: NDArray[np.float64],
    unscented: UnscentedSpec = _DEFAULT_UNSCENTED,
) -> NavFeasibilityOutcome:
    """LinCov one cell along the coast arc ending at apogee; threshold via C0's Φ.

    The reference trajectory is the filter's own two-body+J2 flow (the onboard
    model — what LinCov linearizes about); nodes co-fly as rigid RTN offsets at
    the cell's geometry.  Σ at the final epoch (the apogee correction node) is
    rotated to apogee-RTN and pushed through ``Φ Σ Φᵀ``; the verdict compares the
    RMS lateral (T–N) miss against the catch radius (C0's threshold convention —
    the radial component is pinned by the altitude-event crossing).
    """
    dt_s = 1.0 / cell.cadence_hz
    n_epochs = max(1, round(spec.arc_duration_s * cell.cadence_hz))

    def flow(x: NDArray[np.float64], dt: float) -> NDArray[np.float64]:
        return two_body_j2_flow(x, dt, max_step_s=spec.flow_max_step_s)

    x_start = flow(apogee_state(), -spec.arc_duration_s)
    q = process_noise_white_accel(cell.q_accel_m_s2, dt_s)

    epochs: list[LincovEpoch] = []
    x_ref = x_start
    for _ in range(n_epochs):
        x_ref = flow(x_ref, dt_s)
        epochs.append(
            LincovEpoch(
                dt_s=dt_s,
                process_noise=q,
                measurements=_epoch_measurements(_node_states_at(x_ref, cell, spec), cell),
            )
        )

    prior = np.diag([spec.prior_pos_sigma_m**2] * 3 + [spec.prior_vel_sigma_m_s**2] * 3).astype(
        np.float64
    )
    final = run_lincov(FilterState(x=x_start, cov=prior), flow, tuple(epochs), unscented)[-1]

    sigma_rtn = covariance_to_rtn(final.cov, final.x)
    sigmas = np.sqrt(np.diag(sigma_rtn))
    induced = induced_miss_covariance(np.asarray(phi, dtype=np.float64), sigma_rtn)
    lateral_miss = float(math.sqrt(induced[1, 1] + induced[2, 2]))
    return NavFeasibilityOutcome(
        cell=cell,
        pos_sigma_rtn_m=(float(sigmas[0]), float(sigmas[1]), float(sigmas[2])),
        vel_sigma_rtn_m_s=(float(sigmas[3]), float(sigmas[4]), float(sigmas[5])),
        lateral_miss_1sigma_m=lateral_miss,
        meets_catch_radius=lateral_miss < spec.catch_radius_m,
    )


@dataclass(frozen=True)
class NavFeasibilityResult:
    """A completed C1 feasibility sweep: the spec and per-cell outcomes, grid-aligned."""

    spec: NavFeasibilitySpec
    outcomes: tuple[NavFeasibilityOutcome, ...]


def sweep_nav_feasibility(
    spec: NavFeasibilitySpec,
    phi: NDArray[np.float64],
    unscented: UnscentedSpec = _DEFAULT_UNSCENTED,
) -> NavFeasibilityResult:
    """Evaluate every cell of the C1 sweep — the measured feasibility envelope.

    ``phi`` is C0's measured apogee→crossing sensitivity (ADR 0011); the JVM
    report runner re-derives it via ``run_nav_sweep`` and hands it in here.
    """
    cells = nav_feasibility_cells(spec)
    return NavFeasibilityResult(
        spec=spec,
        outcomes=tuple(evaluate_cell(cell, spec, phi, unscented) for cell in cells),
    )


def _cell_value_label(cell: NavFeasibilityCell) -> str:
    if cell.axis == "nominal":
        return "nominal"
    if cell.axis == "doppler_sigma_m_s" and cell.doppler_sigma_m_s is None:
        return f"{cell.axis}=range-only"
    value = getattr(cell, cell.axis)
    return f"{cell.axis}={value:.4g}"


def format_nav_feasibility(result: NavFeasibilityResult) -> str:
    """Human-readable C1 envelope report — per-cell knowledge and catch-radius verdict."""
    spec = result.spec
    lines = [
        "C1 navigation feasibility — LinCov envelope at the apogee node "
        f"(catch radius {spec.catch_radius_m / 1e3:g} km, arc {spec.arc_duration_s / 3600:g} h)",
        "  cell | T-vel σ [m/s] | vel σ RTN [m/s] | pos σ RTN [m] | lateral 1σ miss [m] | verdict",
    ]
    for outcome in result.outcomes:
        vel = " ".join(f"{s:.3g}" for s in outcome.vel_sigma_rtn_m_s)
        pos = " ".join(f"{s:.3g}" for s in outcome.pos_sigma_rtn_m)
        verdict = "MEETS" if outcome.meets_catch_radius else "fails"
        lines.append(
            f"    {_cell_value_label(outcome.cell)}: {outcome.vel_sigma_rtn_m_s[1]:.3g} | "
            f"{vel} | {pos} | {outcome.lateral_miss_1sigma_m:.3g} | {verdict}"
        )
    return "\n".join(lines)


def node_directions_rtn(
    n_nodes: int, cone_half_angle_rad: float
) -> tuple[NDArray[np.float64], ...]:
    """Unit LOS directions (RTN components) on a cone of given half-angle around R̂.

    The one-parameter LOS-diversity family of ADR 0012 decision 2: equal azimuth
    spacing in the T–N plane, half-angle 0 = collinear (blind to transverse
    velocity), 90° = full in-plane spread.  What the sweep varies is exactly the
    GDOP-style geometry quality the derived node requirement is stated in.
    """
    directions: list[NDArray[np.float64]] = []
    for k in range(n_nodes):
        azimuth = 2.0 * math.pi * k / n_nodes
        directions.append(
            np.array(
                [
                    math.cos(cone_half_angle_rad),
                    math.sin(cone_half_angle_rad) * math.cos(azimuth),
                    math.sin(cone_half_angle_rad) * math.sin(azimuth),
                ],
                dtype=np.float64,
            )
        )
    return tuple(directions)


@dataclass(frozen=True)
class NavValidationOutcome:
    """One seeded UKF truth run judged by NEES (ADR 0012 decision 7, layer 2).

    ``consistent`` is the two-sided chi-square gate on the time-averaged NEES —
    the claim "the filter's actual error is statistically inside its claimed Σ".
    Failing *high* is the dangerous direction (optimistic filter: the LinCov
    envelope would be fiction and q needs retuning); failing low is pessimism.
    ``claimed_t_vel_sigma_m_s`` vs ``actual_t_vel_error_m_s`` is the C0-binding
    axis diagnostic at the final (apogee) epoch.
    """

    cell: NavFeasibilityCell
    n_epochs: int
    average_nees: float
    nees_bounds: tuple[float, float]
    consistent: bool
    claimed_t_vel_sigma_m_s: float
    actual_t_vel_error_m_s: float


def validate_cell(
    cell: NavFeasibilityCell,
    spec: NavFeasibilitySpec,
    truth_states: NDArray[np.float64],
    seed: int,
    unscented: UnscentedSpec = _DEFAULT_UNSCENTED,
) -> NavValidationOutcome:
    """Run the real UKF against a truth arc with seeded measurements; judge by NEES.

    ``truth_states`` is the ``(n_epochs+1, 6)`` truth at the arc start plus each
    measurement epoch, sampled at the cell's cadence (the JVM seam supplies an
    Orekit full-force arc; tests supply synthetic ones).  Nodes ride the truth as
    rigid RTN offsets and the filter knows their ephemerides exactly (decision 1
    — node error lives in R, not in the geometry).  The filter is initialized at
    truth + a seeded prior draw, so prior, measurement noise, and (if the truth
    has none beyond the filter model) dynamics are all consistently modeled and
    the time-averaged NEES should sit inside :func:`average_nees_bounds` — the
    epochs are treated as independent draws, the standard (approximate) NEES
    convention.
    """
    truth = np.asarray(truth_states, dtype=np.float64)
    n_epochs = truth.shape[0] - 1
    dt_s = 1.0 / cell.cadence_hz
    rng = np.random.default_rng(seed)

    def flow(x: NDArray[np.float64]) -> NDArray[np.float64]:
        return two_body_j2_flow(x, dt_s, max_step_s=spec.flow_max_step_s)

    prior = np.diag([spec.prior_pos_sigma_m**2] * 3 + [spec.prior_vel_sigma_m_s**2] * 3).astype(
        np.float64
    )
    state = FilterState(x=truth[0] + rng.multivariate_normal(np.zeros(6), prior), cov=prior)
    q = process_noise_white_accel(cell.q_accel_m_s2, dt_s)

    nees_values = np.empty(n_epochs, dtype=np.float64)
    for k in range(1, n_epochs + 1):
        state = ukf_predict(state, flow, q, unscented)
        for model in _epoch_measurements(_node_states_at(truth[k], cell, spec), cell):
            noise = rng.normal(0.0, math.sqrt(float(model.noise_cov[0, 0])), size=1)
            state = ukf_update(
                state, model.h(truth[k]) + noise, model.h, model.noise_cov, unscented
            )
        nees_values[k - 1] = nees(state.x - truth[k], state.cov)

    average = float(nees_values.mean())
    bounds = average_nees_bounds(dim=6, n_samples=n_epochs)
    basis = rtn_basis(
        (float(truth[-1, 0]), float(truth[-1, 1]), float(truth[-1, 2])),
        (float(truth[-1, 3]), float(truth[-1, 4]), float(truth[-1, 5])),
    )
    error_vel = state.x[3:6] - truth[-1, 3:6]
    actual_t_vel = rtn_components(
        (float(error_vel[0]), float(error_vel[1]), float(error_vel[2])), basis
    )[1]
    sigma_rtn = covariance_to_rtn(state.cov, truth[-1])
    return NavValidationOutcome(
        cell=cell,
        n_epochs=n_epochs,
        average_nees=average,
        nees_bounds=bounds,
        consistent=bounds[0] <= average <= bounds[1],
        claimed_t_vel_sigma_m_s=float(math.sqrt(sigma_rtn[4, 4])),
        actual_t_vel_error_m_s=abs(actual_t_vel),
    )


def format_nav_validation(outcomes: tuple[NavValidationOutcome, ...]) -> str:
    """Human-readable NEES validation report — the layer-2 verdict per validated cell.

    "OPTIMISTIC" (average NEES above the upper bound) is the dangerous direction:
    the filter's claimed Σ — and hence the LinCov envelope — would be fiction
    until q is retuned.  "pessimistic" (below the lower bound) wastes margin but
    does not invalidate a MEETS verdict.
    """
    lines = ["NEES validation — seeded UKF truth runs vs claimed covariance (ADR 0012)"]
    for outcome in outcomes:
        lo, hi = outcome.nees_bounds
        if outcome.consistent:
            verdict = "consistent"
        elif outcome.average_nees > hi:
            verdict = "OPTIMISTIC (claimed Σ too small — retune q)"
        else:
            verdict = "pessimistic (claimed Σ too large)"
        lines.append(
            f"    {_cell_value_label(outcome.cell)}: avg NEES {outcome.average_nees:.2f} "
            f"vs [{lo:.2f}, {hi:.2f}] over {outcome.n_epochs} epochs → {verdict}; "
            f"T-vel σ claimed {outcome.claimed_t_vel_sigma_m_s:.3g} m/s, "
            f"actual |err| {outcome.actual_t_vel_error_m_s:.3g} m/s"
        )
    return "\n".join(lines)
