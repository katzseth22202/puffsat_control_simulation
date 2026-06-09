"""Pure A3 controllability-sweep core — deterministic grid generation (no JVM).

A3 (design doc §13, ADR 0007) holds the targeter fixed and sweeps ``Cd·(A/m)`` and
``Cr·(A/m)`` across a *deterministic* factor grid (not the stochastic
:class:`~puffsat_sim.dispersion.DispersionSpec` draws), injection zeroed, to map
required Δv against coefficient error.  This module is the pure grid side; the JVM
``run_sweep`` (in :mod:`puffsat_sim.montecarlo`) consumes it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import IntEnum
from typing import Final

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.constants import STANDARD_GRAVITY_M_S2
from puffsat_sim.dispersion import RunInputs
from puffsat_sim.records import RunRecord

# ADR 0004: Isp is a *reported sweep*, not a fixed value, and the propellant budget is
# 2% of the 25 kg PuffSat mass (paper claim).  The dry mass cancels in the Δv budget
# (`budget_dv = prop_frac · Isp · g₀`), so only the fraction and Isp anchors are needed.
ISP_ANCHORS_S: Final[tuple[float, float, float]] = (50.0, 70.0, 200.0)
PROPELLANT_BUDGET_FRACTION: Final[float] = 0.02


@dataclass(frozen=True)
class SweepSpec:
    """A deterministic ``Cd``/``Cr`` controllability grid (ADR 0007).

    Coefficients are swept as multiplicative factors on the nominal, log-spaced so
    the axis is linear in σ-equivalent.  ``*_points == 1`` holds that axis at nominal
    (factor 1.0), so 1D cuts and the full 2D grid are the same type.  ``f10p7``/``ap``
    stay nominal — A3 sweeps only the two lumped coefficients.
    """

    cd_area_over_mass: float = 0.04
    cr_area_over_mass: float = 0.02
    cd_factor_range: tuple[float, float] = (0.5, 2.0)
    cr_factor_range: tuple[float, float] = (0.5, 2.0)
    cd_points: int = 9
    cr_points: int = 9
    f10p7: float = 150.0
    ap: float = 15.0


def _axis_factors(factor_range: tuple[float, float], points: int) -> list[float]:
    """Log-spaced factors across ``factor_range``; a single point is nominal (1.0)."""
    if points == 1:
        return [1.0]
    lo, hi = factor_range
    return [float(f) for f in np.geomspace(lo, hi, points)]


def grid_inputs(spec: SweepSpec) -> tuple[RunInputs, ...]:
    """Enumerate the deterministic grid as zero-injection ``RunInputs`` (ADR 0007).

    Cd is the outer axis and Cr the inner, so ``run_index = cd_i·cr_points + cr_i`` —
    the row-major ordering :func:`to_grid` relies on to reshape records back to 2D.
    """
    cd_factors = _axis_factors(spec.cd_factor_range, spec.cd_points)
    cr_factors = _axis_factors(spec.cr_factor_range, spec.cr_points)
    inputs: list[RunInputs] = []
    index = 0
    for cd_f in cd_factors:
        for cr_f in cr_factors:
            inputs.append(
                RunInputs(
                    run_index=index,
                    dv_rtn_m_s=(0.0, 0.0, 0.0),
                    cd_area_over_mass=spec.cd_area_over_mass * cd_f,
                    cr_area_over_mass=spec.cr_area_over_mass * cr_f,
                    f10p7=spec.f10p7,
                    ap=spec.ap,
                )
            )
            index += 1
    return tuple(inputs)


@dataclass(frozen=True)
class SweepResult:
    """A completed controllability sweep (ADR 0007): the spec, per-point records, nominal.

    ``records`` are the grid points (one per :func:`grid_inputs` entry); ``nominal`` is the
    factor-(1,1) reference run (perfect aim, ~0 required Δv) the perigee/ToA overlays read
    against.  Pure value type — produced by the JVM ``run_sweep`` but JVM-free here.
    """

    spec: SweepSpec
    records: tuple[RunRecord, ...]
    nominal: RunRecord


@dataclass(frozen=True)
class SweepGrid:
    """The sweep reshaped onto its 2D ``Cd × Cr`` axes — the controllability-map arrays.

    Axes are row-major ``(cd_points, cr_points)``: ``[i, j]`` is Cd factor ``i``, Cr factor
    ``j``.  ``required_dv_m_s`` is the mapped cost (ADR 0007 decision 2); ``converged`` marks
    the uncontrollable region; ``perigee_alt_m`` is the debris-disposal overlay.
    """

    cd_area_over_mass: NDArray[np.float64]  # (cd_points,) swept Cd·(A/m) values
    cr_area_over_mass: NDArray[np.float64]  # (cr_points,) swept Cr·(A/m) values
    cd_factors: NDArray[np.float64]  # (cd_points,) multiplicative factors on nominal
    cr_factors: NDArray[np.float64]  # (cr_points,)
    required_dv_m_s: NDArray[np.float64]  # (cd_points, cr_points)
    converged: NDArray[np.bool_]  # (cd_points, cr_points)
    perigee_alt_m: NDArray[np.float64]  # (cd_points, cr_points)


def to_grid(records: tuple[RunRecord, ...], spec: SweepSpec) -> SweepGrid:
    """Reshape flat per-point records into the 2D ``Cd × Cr`` controllability map.

    Records are ordered by ``inputs.run_index`` (so an out-of-order / resumed set still
    lands on the right cell) and reshaped row-major to ``(cd_points, cr_points)`` — the
    same layout :func:`grid_inputs` emits.  Raises ``ValueError`` if the record count does
    not match the grid (a partial sweep is not a map).
    """
    expected = spec.cd_points * spec.cr_points
    if len(records) != expected:
        raise ValueError(f"expected {expected} records for the grid, got {len(records)}")
    ordered = sorted(records, key=lambda r: r.inputs.run_index)
    shape = (spec.cd_points, spec.cr_points)
    cd_factors = np.asarray(_axis_factors(spec.cd_factor_range, spec.cd_points), dtype=np.float64)
    cr_factors = np.asarray(_axis_factors(spec.cr_factor_range, spec.cr_points), dtype=np.float64)
    return SweepGrid(
        cd_area_over_mass=spec.cd_area_over_mass * cd_factors,
        cr_area_over_mass=spec.cr_area_over_mass * cr_factors,
        cd_factors=cd_factors,
        cr_factors=cr_factors,
        required_dv_m_s=np.array([r.total_dv_m_s for r in ordered]).reshape(shape),
        converged=np.array([r.converged for r in ordered], dtype=np.bool_).reshape(shape),
        perigee_alt_m=np.array([r.perigee_alt_m for r in ordered]).reshape(shape),
    )


def sigma_equivalent(factor: NDArray[np.float64], cv: float) -> NDArray[np.float64]:
    """Map a multiplicative factor to its σ in the Rung-D log-normal (ADR 0007 decision 4).

    Inverts ``factor = exp(k·s)`` with ``s = √(ln(1+cv²))`` — the same dispersion the
    capstone samples (``dispersion._lognormal_factor``) — so a grid axis reads in σ of the
    coefficient error: ``k = ln(factor)/s``.  Factor 1 maps to 0; ``factor`` and ``1/factor``
    are ±k.
    """
    s = math.sqrt(math.log(1.0 + cv * cv))
    return np.asarray(np.log(factor) / s, dtype=np.float64)


def budget_dv_m_s(isp_s: float, prop_frac: float = PROPELLANT_BUDGET_FRACTION) -> float:
    """Δv a propellant fraction buys at this Isp (ADR 0004: ``prop_frac ≈ Δv/(Isp·g₀)``)."""
    return prop_frac * isp_s * STANDARD_GRAVITY_M_S2


class Controllability(IntEnum):
    """A grid point's mission-meaning region (ADR 0007 decision 2)."""

    CONTROLLABLE = 0  # a solution exists within the Δv budget
    OVER_BUDGET = 1  # a solution exists but costs more propellant than allowed ("buy Isp/mass")
    UNCONTROLLABLE = 2  # no valid solution at all ("physically unreachable")


def classify_controllability(
    grid: SweepGrid, isp_s: float, prop_frac: float = PROPELLANT_BUDGET_FRACTION
) -> NDArray[np.int_]:
    """Label each grid cell controllable / over-budget / uncontrollable at one Isp anchor.

    The budget is post-processing (ADR 0007 decision 2): the solver returns raw required Δv,
    and the Isp/fraction draw the contour.  Non-convergence (no valid root) takes precedence
    over the budget test — an uncontrollable point is never merely "over budget".
    """
    budget = budget_dv_m_s(isp_s, prop_frac)
    labels = np.full(grid.required_dv_m_s.shape, int(Controllability.CONTROLLABLE), dtype=int)
    labels[grid.converged & (grid.required_dv_m_s > budget)] = int(Controllability.OVER_BUDGET)
    labels[~grid.converged] = int(Controllability.UNCONTROLLABLE)
    return labels
