"""Atmospheric drag perturbation — spec + analytic density/deceleration (pure).

The piecewise-exponential atmosphere here is an order-of-magnitude *analytic*
estimate (good to ~factor-of-2); the truth propagation uses Orekit's full
NRLMSISE-00 model, wired in :mod:`puffsat_sim.forces.build`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True)
class AtmosphericDrag:
    """Cannonball atmospheric drag.

    cd_area_over_mass is the lumped Cd·(A/m) [m²/kg]; Cd is folded in so the
    Orekit model uses an effective coefficient of 1.  f10p7 / ap are the solar
    flux and geomagnetic indices for the Monte Carlo per-run drag bias (Rung D).
    """

    cd_area_over_mass: float
    f10p7: float = 150.0  # solar flux index
    ap: float = 15.0  # geomagnetic index


# (base_altitude_m, density_kg_m3, scale_height_m) — calibrated to NRLMSISE-00
# at moderate solar activity (F10.7≈150, Ap≈15).  Good to ~factor-of-2 for
# order-of-magnitude checks; not a substitute for the full model in build.py.
_STD_ATM_LAYERS: Final[tuple[tuple[float, float, float], ...]] = (
    (0, 1.225, 8_500),
    (25_000, 3.9e-2, 6_700),
    (50_000, 1.0e-3, 7_200),
    (80_000, 1.0e-5, 9_000),
    (100_000, 5.6e-7, 15_500),
    (150_000, 2.2e-9, 22_000),
    (200_000, 2.5e-10, 29_000),
    (300_000, 8.0e-12, 37_000),
    (500_000, 5.5e-13, 60_000),
    (700_000, 3.6e-14, 73_000),
)


def std_atm_density(altitude_m: float) -> float:
    """Piecewise-exponential atmospheric density [kg/m³] at the given altitude.

    Calibrated to NRLMSISE-00 at moderate solar activity (F10.7≈150, Ap≈15).
    Accurate to roughly a factor of 2 below ~800 km; use for sanity checks only.
    """
    for i in range(len(_STD_ATM_LAYERS) - 1):
        if altitude_m < _STD_ATM_LAYERS[i + 1][0]:
            h0, rho0, scale_h = _STD_ATM_LAYERS[i]
            return rho0 * math.exp(-(altitude_m - h0) / scale_h)
    h0, rho0, scale_h = _STD_ATM_LAYERS[-1]
    return rho0 * math.exp(-(altitude_m - h0) / scale_h)


def drag_deceleration(cd_area_over_mass: float, speed_m_s: float, altitude_m: float) -> float:
    """Drag deceleration magnitude [m/s²] using piecewise-exponential density.

    a_drag = ½ · ρ(h) · v² · (Cd·A/m)
    """
    return 0.5 * std_atm_density(altitude_m) * speed_m_s**2 * cd_area_over_mass
