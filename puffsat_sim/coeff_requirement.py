"""Pure C2a coefficient-knowledge requirement core (no JVM, ADR 0013).

C2a asks how well the lumped ``Cd·(A/m)`` / ``Cr·(A/m)`` must be *known at the apogee
burn* for the interception to survive.  Because deploy ≈ apogee ≈ burn (B1) and
along-track authority is apogee-bound (ADR 0006), the knowledge the corrector acts on
is the **ground prior** — so the requirement verdict is tolerance vs prior, not filter
convergence.  The chain: the corrector burns ``Δv(ĉ)`` believing ĉ while truth flies c,
so the residual interception miss ≈ ``Φ_lat,vel · (∂Δv/∂c) · δc``, with ``∂Δv/∂c``
measured from a 1D A3 cut's per-cell Δv vectors (``RunRecord.control_log``) and Φ from
C0 (ADR 0011).  This module is the pure side; the JVM ``coeff_requirement_report``
(in :mod:`puffsat_sim.montecarlo`) feeds it the live cuts.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.constants import SRP_P0_PA
from puffsat_sim.sweep import SweepResult, axis_factors

# The interception miss is judged on the lateral (T, N) crossing components — the radial
# is pinned ~0 by the 200 km altitude-event crossing (ADR 0011 decision 4) — and a
# coefficient-induced Δv error enters through Φ's velocity columns only.
_LATERAL_AXES: Final[list[int]] = [1, 2]
_VELOCITY_COLS: Final[slice] = slice(3, 6)


def cut_dv_vectors(
    result: SweepResult,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Extract (factors, per-cell RTN Δv vectors) from a 1D A3 cut (ADR 0013 decision 1).

    Exactly one of the spec's axes must be swept (a cut, not a grid).  Records are
    ordered by ``run_index`` (the ``grid_inputs`` row-major convention, which for a 1D
    cut is factor order); each cell's Δv is the vector sum of its commanded actions.
    A non-converged point has no trustworthy Δv, so it is an error, not a NaN.
    """
    spec = result.spec
    if (spec.cd_points > 1) == (spec.cr_points > 1):
        raise ValueError("expected a 1D cut: exactly one of cd_points/cr_points > 1")
    swept_cd = spec.cd_points > 1
    factor_range = spec.cd_factor_range if swept_cd else spec.cr_factor_range
    points = spec.cd_points if swept_cd else spec.cr_points
    factors = np.asarray(axis_factors(factor_range, points), dtype=np.float64)

    ordered = sorted(result.records, key=lambda r: r.inputs.run_index)
    if len(ordered) != points:
        raise ValueError(f"expected {points} cut records, got {len(ordered)}")
    for factor, record in zip(factors, ordered, strict=True):
        if not record.converged:
            raise ValueError(f"non-converged cut point at factor {factor:g}")

    dvs = np.array(
        [
            np.sum([action.dv_rtn_m_s for action in record.control_log], axis=0)
            if record.control_log
            else np.zeros(3)
            for record in ordered
        ],
        dtype=np.float64,
    )
    return factors, dvs


def dv_gradient(factors: NDArray[np.float64], dv_rtn: NDArray[np.float64]) -> NDArray[np.float64]:
    """``∂Δv/∂c`` at nominal: central difference over the bracketing pair nearest factor 1.

    The same "slope at zero from the nearest ± pair" philosophy as C0's
    ``navigation._central_slope``, in factor space: the cut must straddle the nominal
    (factor 1.0), and the slope is taken between the largest factor below it and the
    smallest above it.  Returns the RTN Δv gradient per unit factor.
    """
    f = np.asarray(factors, dtype=np.float64)
    dv = np.asarray(dv_rtn, dtype=np.float64)
    below = f[f < 1.0]
    above = f[f > 1.0]
    if below.size == 0 or above.size == 0:
        raise ValueError("cut factors must straddle nominal (factor 1.0)")
    i_lo = int(np.flatnonzero(f == below.max())[0])
    i_hi = int(np.flatnonzero(f == above.min())[0])
    return np.asarray((dv[i_hi] - dv[i_lo]) / (f[i_hi] - f[i_lo]), dtype=np.float64)


def coefficient_sensitivity(phi: NDArray[np.float64], gradient: NDArray[np.float64]) -> float:
    """Lateral interception miss per unit coefficient-factor error: ``‖Φ_lat,vel · G‖``."""
    phi_arr = np.asarray(phi, dtype=np.float64)
    g = np.asarray(gradient, dtype=np.float64)
    return float(np.linalg.norm(phi_arr[_LATERAL_AXES, _VELOCITY_COLS] @ g))


def coefficient_tolerance(
    phi: NDArray[np.float64], gradient: NDArray[np.float64], catch_radius_m: float
) -> float:
    """The 1σ coefficient-factor error the catch radius can carry (``inf`` if insensitive)."""
    sensitivity = coefficient_sensitivity(phi, gradient)
    return math.inf if sensitivity == 0.0 else catch_radius_m / sensitivity


def analytic_srp_dv(
    cr_area_over_mass: float,
    coast_duration_s: float,
    srp_pressure_pa: float = SRP_P0_PA,
) -> float:
    """Analytic SRP Δv impulse over the coast per unit Cr factor: ``Cr·(A/m)·P₀·t``.

    The order-of-magnitude cross-check on the measured ``‖∂Δv/∂Cr‖`` (ADR 0013): the
    full SRP acceleration at the nominal coefficient, integrated over the coast, is the
    Δv the corrector must re-aim per 100% coefficient error (an upper bound — eclipse
    and geometry projection trim it).
    """
    return cr_area_over_mass * srp_pressure_pa * coast_duration_s


@dataclass(frozen=True)
class BudgetEntry:
    """One lateral-miss contribution in the catch-radius error budget (ADR 0013 decision 5)."""

    label: str
    lateral_miss_m: float


# Measured lateral-miss contributions banked by earlier rungs: C1's nav-induced lateral
# at the NEES-honest q (ADR 0012 findings) and B1's finite-burn erosion (ADR 0008).
MEASURED_BUDGET: Final[tuple[BudgetEntry, ...]] = (
    BudgetEntry("nav-induced lateral (C1, honest q)", 141.0),
    BudgetEntry("finite-burn erosion (B1)", 89.0),
)


def rss_lateral(entries: Sequence[BudgetEntry]) -> float:
    """Root-sum-square of the budget's lateral contributions (independent error sources)."""
    return math.sqrt(sum(e.lateral_miss_m**2 for e in entries))


@dataclass(frozen=True)
class AxisRequirement:
    """One coefficient axis reduced to its knowledge requirement (ADR 0013 decision 1).

    ``tolerance_factor`` and ``prior_sigma_factor``-derived fields are in A3's
    multiplicative-factor units (1.0 = 100% of nominal); ``margin`` is
    tolerance / prior (``inf`` when the axis is unconstrained).
    """

    axis: str
    dv_gradient_m_s: tuple[float, float, float]
    lateral_sensitivity_m: float
    tolerance_factor: float
    prior_lateral_miss_m: float
    covered_by_prior: bool
    margin: float


@dataclass(frozen=True)
class CoeffRequirement:
    """The C2a verdict: both axes, the prior, and the analytic SRP cross-check."""

    cd: AxisRequirement
    cr: AxisRequirement
    catch_radius_m: float
    prior_sigma_factor: float
    analytic_srp_dv_m_s: float | None
    measured_over_analytic: float | None


def _axis_requirement(
    axis: str,
    phi: NDArray[np.float64],
    cut: SweepResult,
    catch_radius_m: float,
    prior_sigma_factor: float,
) -> AxisRequirement:
    factors, dvs = cut_dv_vectors(cut)
    gradient = dv_gradient(factors, dvs)
    sensitivity = coefficient_sensitivity(phi, gradient)
    tolerance = coefficient_tolerance(phi, gradient, catch_radius_m)
    margin = math.inf if tolerance == math.inf else tolerance / prior_sigma_factor
    return AxisRequirement(
        axis=axis,
        dv_gradient_m_s=(float(gradient[0]), float(gradient[1]), float(gradient[2])),
        lateral_sensitivity_m=sensitivity,
        tolerance_factor=tolerance,
        prior_lateral_miss_m=prior_sigma_factor * sensitivity,
        covered_by_prior=prior_sigma_factor <= tolerance,
        margin=margin,
    )


def summarize_coeff_requirement(
    phi: NDArray[np.float64],
    *,
    cd_cut: SweepResult,
    cr_cut: SweepResult,
    catch_radius_m: float = 5_000.0,
    prior_sigma_factor: float = 0.2,
    coast_duration_s: float | None = None,
) -> CoeffRequirement:
    """Reduce the two 1D cuts + C0's Φ to the C2a requirement verdict (ADR 0013).

    ``prior_sigma_factor`` defaults to 0.2 — conservative for a manufactured balloon
    area (few %) plus material ``Cr`` (~10–20%).  When ``coast_duration_s`` is given,
    the measured ``‖∂Δv/∂Cr‖`` is cross-checked against the analytic SRP impulse.
    """
    cd = _axis_requirement("Cd", phi, cd_cut, catch_radius_m, prior_sigma_factor)
    cr = _axis_requirement("Cr", phi, cr_cut, catch_radius_m, prior_sigma_factor)

    analytic: float | None = None
    ratio: float | None = None
    if coast_duration_s is not None:
        analytic = analytic_srp_dv(cr_cut.spec.cr_area_over_mass, coast_duration_s)
        ratio = float(np.linalg.norm(cr.dv_gradient_m_s)) / analytic
    return CoeffRequirement(
        cd=cd,
        cr=cr,
        catch_radius_m=catch_radius_m,
        prior_sigma_factor=prior_sigma_factor,
        analytic_srp_dv_m_s=analytic,
        measured_over_analytic=ratio,
    )


def _axis_line(req: AxisRequirement, prior_sigma_factor: float) -> str:
    gradient_norm = float(np.linalg.norm(req.dv_gradient_m_s))
    if req.tolerance_factor == math.inf:
        verdict = "UNCONSTRAINED (no lateral sensitivity at apogee)"
    elif req.covered_by_prior:
        verdict = f"COVERED by prior ~{req.margin:.1f}x"
    else:
        verdict = f"NOT COVERED (prior {prior_sigma_factor:g} > tolerance)"
    return (
        f"    {req.axis}·(A/m): |dDv/dc| {gradient_norm:.3g} m/s per 1.0 factor | "
        f"lateral {req.lateral_sensitivity_m:.3g} m per factor | "
        f"tolerance {req.tolerance_factor:.3g} | "
        f"prior miss {req.prior_lateral_miss_m:.3g} m | {verdict}"
    )


def format_coeff_requirement(
    req: CoeffRequirement, measured_budget: Sequence[BudgetEntry] = MEASURED_BUDGET
) -> str:
    """Human-readable C2a report — per-axis verdicts, SRP cross-check, RSS error budget."""
    lines = [
        "C2a coefficient-knowledge requirement — prior vs tolerance at the apogee burn (ADR 0013)",
        f"  catch radius {req.catch_radius_m:g} m | prior sigma "
        f"{req.prior_sigma_factor:g} (factor units)",
        _axis_line(req.cr, req.prior_sigma_factor),
        _axis_line(req.cd, req.prior_sigma_factor),
    ]
    if req.analytic_srp_dv_m_s is not None and req.measured_over_analytic is not None:
        lines.append(
            f"  analytic SRP impulse {req.analytic_srp_dv_m_s:.3g} m/s per factor | "
            f"measured/analytic {req.measured_over_analytic:.2f}"
        )
    ledger = (
        *measured_budget,
        BudgetEntry("coefficient prior (Cr)", req.cr.prior_lateral_miss_m),
        BudgetEntry("coefficient prior (Cd)", req.cd.prior_lateral_miss_m),
    )
    total = rss_lateral(ledger)
    lines.append(f"  error budget (lateral 1-sigma, RSS vs catch radius {req.catch_radius_m:g} m):")
    lines += [f"    {e.label}: {e.lateral_miss_m:.3g} m" for e in ledger]
    headroom = math.inf if total == 0.0 else req.catch_radius_m / total
    lines.append(f"    RSS {total:.3g} m -> headroom ~{headroom:.1f}x")
    return "\n".join(lines)
