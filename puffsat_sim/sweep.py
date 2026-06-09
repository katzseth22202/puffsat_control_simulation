"""Pure A3 controllability-sweep core — deterministic grid generation (no JVM).

A3 (design doc §13, ADR 0007) holds the targeter fixed and sweeps ``Cd·(A/m)`` and
``Cr·(A/m)`` across a *deterministic* factor grid (not the stochastic
:class:`~puffsat_sim.dispersion.DispersionSpec` draws), injection zeroed, to map
required Δv against coefficient error.  This module is the pure grid side; the JVM
``run_sweep`` (in :mod:`puffsat_sim.montecarlo`) consumes it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from puffsat_sim.dispersion import RunInputs


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
    """Enumerate the deterministic grid as zero-injection ``RunInputs`` (ADR 0007)."""
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
