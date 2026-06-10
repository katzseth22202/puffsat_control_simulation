"""Pure C0 navigation-requirement core — nav-error sweep grid (no JVM).

C0 (design doc §13, ADR 0011) asks how accurately the PuffSat must know its apogee
state for the interception to survive.  Because the apogee corrector *cancels* a
navigation error to first order (it applies the same correction in predict and
execute), the residual interception miss is the apogee→crossing sensitivity Φ (a 3×6
STM) times the nav error — so C0 is a deterministic per-component **sensitivity** sweep,
not a sampled ensemble.  This module is the pure grid side; the JVM ``run_nav_sweep``
(in :mod:`puffsat_sim.montecarlo`) consumes it.

The nav error is a 6-vector in the **apogee-RTN** frame — position (R/T/N) then
velocity (R/T/N), the same frame the injection is sampled in (T is the ``dr_p/dv_a``
lever).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.records import RunRecord

Vec6 = tuple[float, float, float, float, float, float]

_N_AXES = 6  # 0-2 = position R/T/N, 3-5 = velocity R/T/N
_FIRST_VELOCITY_AXIS = 3


@dataclass(frozen=True)
class NavSweepSpec:
    """A deterministic 6-axis apogee-RTN nav-error sweep (ADR 0011).

    Each of the 6 components (R/T/N position, R/T/N velocity) is perturbed one at a
    time over a log-spaced magnitude range, both signs.  Defaults bracket what
    coordinator-node ranging can plausibly deliver through to gross failure (ADR 0011
    decision 7): position ~cm floor → 10 km catch-radius scale; velocity ~0.1 mm/s
    Doppler floor → 1 m/s.
    """

    pos_range_m: tuple[float, float] = (1e-1, 1e4)
    vel_range_m_s: tuple[float, float] = (1e-4, 1e0)
    points_per_sign: int = 5


@dataclass(frozen=True)
class NavOffsetCell:
    """One sweep cell: a single-component apogee-RTN nav-error offset.

    ``axis`` is the perturbed component (0-5; ``-1`` for the zero/nominal cell);
    ``magnitude`` is the signed offset along it; ``offset_rtn6`` is the full 6-vector
    (exactly one nonzero, or all-zero for the nominal cell).
    """

    cell_index: int
    axis: int
    magnitude: float
    offset_rtn6: Vec6


@dataclass(frozen=True)
class NavSweepResult:
    """A completed C0 nav-error sensitivity sweep (ADR 0011): the spec, cells, and per-cell records.

    Pure value type — produced by the JVM ``run_nav_sweep`` but JVM-free here, so Φ assembly
    and tolerance post-processing stay unit-testable.  ``records[k]`` is the corrector-in-loop
    outcome for ``cells[k]``: its ``miss_rtn_m`` is the residual ``−Φδ`` and ``total_dv_m_s`` the
    phantom correction the corrector burned chasing the (unobserved) nav error.
    """

    spec: NavSweepSpec
    cells: tuple[NavOffsetCell, ...]
    records: tuple[RunRecord, ...]


def _axis_magnitudes(value_range: tuple[float, float], points: int) -> list[float]:
    """Log-spaced magnitudes across ``value_range`` (one per sweep point, per sign)."""
    lo, hi = value_range
    return [float(m) for m in np.geomspace(lo, hi, points)]


def _offset_with(axis: int, magnitude: float) -> Vec6:
    components = [0.0] * _N_AXES
    components[axis] = magnitude
    return (
        components[0],
        components[1],
        components[2],
        components[3],
        components[4],
        components[5],
    )


def nav_grid_offsets(spec: NavSweepSpec) -> tuple[NavOffsetCell, ...]:
    """Enumerate the one-component-at-a-time nav-error cells plus a single zero cell.

    The zero cell (``axis=-1``, ``magnitude=0``) is the nominal reference whose residual
    is ~0; every other cell perturbs exactly one apogee-RTN component by a signed,
    log-spaced magnitude (position axes use ``pos_range_m``, velocity axes
    ``vel_range_m_s``), so the residual response per axis builds that column of Φ.
    """
    cells: list[NavOffsetCell] = [NavOffsetCell(0, -1, 0.0, _offset_with(0, 0.0))]
    index = 1
    for axis in range(_N_AXES):
        value_range = spec.vel_range_m_s if axis >= _FIRST_VELOCITY_AXIS else spec.pos_range_m
        for mag in _axis_magnitudes(value_range, spec.points_per_sign):
            for signed in (mag, -mag):
                cells.append(NavOffsetCell(index, axis, signed, _offset_with(axis, signed)))
                index += 1
    return tuple(cells)


def _central_slope(
    magnitudes: NDArray[np.float64], misses_m: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Central-difference slope at zero from the smallest |magnitude| ± pair (one Φ column).

    The single definition of "the sensitivity at zero" shared by :func:`assemble_sensitivity`
    and :func:`linearity_range`, so Φ and its linearity check agree by construction:
    ``(miss(+m) − miss(−m)) / (m₊ − m₋)``.
    """
    by_magnitude = list(np.argsort(np.abs(magnitudes)))
    pos = next(i for i in by_magnitude if magnitudes[i] > 0.0)
    neg = next(i for i in by_magnitude if magnitudes[i] < 0.0)
    return np.asarray(
        (misses_m[pos] - misses_m[neg]) / (magnitudes[pos] - magnitudes[neg]), dtype=np.float64
    )


def assemble_sensitivity(
    cells: Sequence[NavOffsetCell], misses_m: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Assemble the 3×6 apogee→crossing sensitivity Φ from per-cell residual misses.

    Each column is the small-offset slope at zero along that apogee-RTN axis (ADR 0011
    decision 1).  ``misses_m`` is ``(N, 3)`` in the nominal-crossing RTN frame, aligned to
    ``cells``.
    """
    misses = np.asarray(misses_m, dtype=np.float64)
    magnitudes = np.array([c.magnitude for c in cells], dtype=np.float64)
    phi = np.empty((3, _N_AXES), dtype=np.float64)
    for axis in range(_N_AXES):
        on_axis = [i for i, c in enumerate(cells) if c.axis == axis]
        phi[:, axis] = _central_slope(magnitudes[on_axis], misses[on_axis])
    return phi


# The interception miss is in the nominal-crossing RTN frame; its radial (R) component is
# pinned ~0 by the 200 km altitude-event crossing, so the catch-radius threshold is on the
# lateral (Transverse, Normal) miss (ADR 0011 decision 4).
_LATERAL_AXES: tuple[int, int] = (1, 2)


def induced_miss_covariance(
    phi: NDArray[np.float64], sigma6: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Interception-miss covariance ``Φ Σ Φᵀ`` from an apogee-RTN nav covariance Σ (ADR 0011).

    Because the residual is linear in the nav error, Φ is the complete requirement: any Σ
    (diagonal or geometry-structured) maps to its 3×3 induced miss covariance here, with no
    new propagation — the C1 UKF's achieved covariance is checked against the threshold via
    this one matrix product.
    """
    phi_arr = np.asarray(phi, dtype=np.float64)
    sigma = np.asarray(sigma6, dtype=np.float64)
    return np.asarray(phi_arr @ sigma @ phi_arr.T, dtype=np.float64)


def axis_tolerance(phi: NDArray[np.float64], catch_radius_m: float) -> NDArray[np.float64]:
    """Per-axis nav tolerance: the 1σ apogee-RTN error each axis can carry within the catch radius.

    A single-axis nav error of magnitude σ induces a lateral interception miss
    ``‖Φ_lateral[:, axis]‖ · σ``; holding that within ``catch_radius_m`` (the §9 terminal
    authority, ADR 0011 decision 4) gives ``tol[axis] = catch_radius / ‖Φ_lateral[:, axis]‖``.
    Only the lateral (T, N) rows count — the radial miss is pinned ~0 by the altitude-event
    crossing.  An axis with no lateral sensitivity is unconstrained (``inf``).
    """
    phi_arr = np.asarray(phi, dtype=np.float64)
    lateral_col_norms = np.linalg.norm(phi_arr[_LATERAL_AXES, :], axis=0)
    with np.errstate(divide="ignore"):
        return np.asarray(catch_radius_m / lateral_col_norms, dtype=np.float64)


def linearity_range(
    magnitudes: NDArray[np.float64], misses_m: NDArray[np.float64], rel_tol: float
) -> float:
    """Largest |magnitude| over which a single-axis response stays linear (ADR 0011).

    Validates the ``−Φδ`` approximation: take the slope from the smallest |magnitude| pair,
    then walk outward in increasing |magnitude| and return the extent of the contiguous
    region where ``‖miss − slope·m‖ ≤ rel_tol · ‖slope·m‖``.  The first cell that departs
    (in increasing |magnitude|) caps the linear region — beyond it the corrector's phantom
    correction has driven the executed arc nonlinear.  ``magnitudes`` is the signed 1D axis
    sweep; ``misses_m`` the aligned ``(N, 3)`` crossing-frame misses.
    """
    mags = np.asarray(magnitudes, dtype=np.float64)
    miss = np.asarray(misses_m, dtype=np.float64)
    slope = _central_slope(mags, miss)

    extent = 0.0
    for i in list(np.argsort(np.abs(mags))):
        predicted = slope * mags[i]
        predicted_norm = float(np.linalg.norm(predicted))
        if predicted_norm == 0.0:
            continue
        residual = float(np.linalg.norm(miss[i] - predicted))
        if residual <= rel_tol * predicted_norm:
            extent = max(extent, abs(float(mags[i])))
        else:
            break
    return extent


def sample_nav_error(rng: np.random.Generator, sigma6: NDArray[np.float64]) -> Vec6:
    """One seeded Gaussian apogee-RTN nav-error draw (mean 0, covariance Σ).

    The only place randomness enters C0 — used for the single representative-Σ cross-check
    that ``Φ Σ Φᵀ`` predicts the sample interception-miss covariance (ADR 0011 decision 2).
    """
    draw = rng.multivariate_normal(np.zeros(_N_AXES), np.asarray(sigma6, dtype=np.float64))
    return (
        float(draw[0]),
        float(draw[1]),
        float(draw[2]),
        float(draw[3]),
        float(draw[4]),
        float(draw[5]),
    )


@dataclass(frozen=True)
class NavRequirement:
    """The C0 navigation-accuracy requirement reduced from a sweep (ADR 0011).

    ``axis_lateral_sensitivity`` is ``‖Φ_lateral[:, j]‖`` per axis — the dominance metric
    (lateral interception-miss per unit nav error); ``axis_linearity_range`` the extent over
    which ``−Φδ`` holds; ``phantom_dv_*`` the Δv the corrector burned chasing the unobserved
    error (ADR 0011 decision 5).
    """

    phi: NDArray[np.float64]
    axis_lateral_sensitivity: NDArray[np.float64]
    axis_linearity_range: NDArray[np.float64]
    phantom_dv_max_m_s: float
    phantom_dv_mean_m_s: float
    converged_fraction: float


def summarize_nav_requirement(result: NavSweepResult, rel_tol: float = 0.1) -> NavRequirement:
    """Reduce a nav-error sweep to the C0 requirement (ADR 0011).

    Φ, per-axis dominance and linearity range, and the phantom-Δv the corrector burned.
    """
    misses = np.array([r.miss_rtn_m for r in result.records], dtype=np.float64)
    phi = assemble_sensitivity(result.cells, misses)
    lateral = np.linalg.norm(phi[_LATERAL_AXES, :], axis=0)

    magnitudes = np.array([c.magnitude for c in result.cells], dtype=np.float64)
    axes = np.array([c.axis for c in result.cells])
    linear = np.empty(_N_AXES, dtype=np.float64)
    for axis in range(_N_AXES):
        on_axis = axes == axis
        linear[axis] = linearity_range(magnitudes[on_axis], misses[on_axis], rel_tol)

    perturbed = [
        rec for cell, rec in zip(result.cells, result.records, strict=True) if cell.axis >= 0
    ]
    phantom = np.array([rec.total_dv_m_s for rec in perturbed], dtype=np.float64)
    converged = np.array([rec.converged for rec in perturbed], dtype=np.bool_)
    return NavRequirement(
        phi=phi,
        axis_lateral_sensitivity=np.asarray(lateral, dtype=np.float64),
        axis_linearity_range=linear,
        phantom_dv_max_m_s=float(phantom.max()),
        phantom_dv_mean_m_s=float(phantom.mean()),
        converged_fraction=float(converged.mean()),
    )


_AXIS_NAMES: tuple[str, ...] = ("R-pos", "T-pos", "N-pos", "R-vel", "T-vel", "N-vel")
_AXIS_UNITS: tuple[str, ...] = ("m", "m", "m", "m/s", "m/s", "m/s")


def format_nav_requirement(req: NavRequirement, catch_radii_m: Sequence[float]) -> str:
    """Human-readable C0 requirement report — dominance ranking + per-axis tolerance table."""
    ranking = np.argsort(req.axis_lateral_sensitivity)[::-1]
    tol_by_radius = [(r, axis_tolerance(req.phi, r)) for r in catch_radii_m]

    lines = [
        "C0 navigation requirement — apogee-RTN error → 200 km lateral interception miss",
        f"  Corrector converged: {req.converged_fraction * 100:.0f}% of perturbed cells",
        f"  Phantom Δv (chasing the unobserved nav error): mean {req.phantom_dv_mean_m_s:.4f},"
        f" max {req.phantom_dv_max_m_s:.4f} m/s",
        "  Dominance (lateral miss per unit error ‖Φ_TN‖, largest first): "
        + ", ".join(f"{_AXIS_NAMES[j]} {req.axis_lateral_sensitivity[j]:.3g}" for j in ranking),
        "  Per-axis: sensitivity [m miss/unit] | linear-to [unit] | "
        + " | ".join(f"tol@{r / 1e3:g}km" for r in catch_radii_m),
    ]
    for axis in range(_N_AXES):
        tols = " | ".join(f"{tol[axis]:.3g}" for _, tol in tol_by_radius)
        lines.append(
            f"    {_AXIS_NAMES[axis]} [{_AXIS_UNITS[axis]}]: "
            f"{req.axis_lateral_sensitivity[axis]:.3g} | "
            f"{req.axis_linearity_range[axis]:.3g} | {tols}"
        )
    return "\n".join(lines)
