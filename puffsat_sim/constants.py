"""Physical constants — single source of truth for the pure-Python analytic side.

The Earth values match Orekit's ``org.orekit.utils.Constants`` (WGS84); an
integration test pins them together so the analytic predictions and the Orekit
truth propagation cannot silently drift apart.

Bulky model data tied to one function (e.g. the piecewise atmosphere table)
stays with that function rather than living here.
"""

from __future__ import annotations

from typing import Final

EARTH_RADIUS_M: Final[float] = 6_378_137.0  # WGS84 equatorial radius [m]
WGS84_MU: Final[float] = 3.986_004_418e14  # Earth gravitational parameter [m³/s²]
J2: Final[float] = 1.08262668e-3  # EGM2008 zonal harmonic J2
SPEED_OF_LIGHT_M_S: Final[float] = 299_792_458.0  # exact; for the relativistic correction
STANDARD_GRAVITY_M_S2: Final[float] = 9.80665  # g₀ in the Isp→Δv rocket relation (exact, CGPM)
PLANCK_J_S: Final[float] = 6.62607015e-34  # exact (SI 2019); photon energy in the tracker budget

# Solar radiation pressure at 1 AU [Pa = N/m²] — used for analytic SRP estimates.
SRP_P0_PA: Final[float] = 4.56e-6

# Third-body gravitational parameters and mean distances.  Distances are mean
# values; instantaneous geometry varies but these suffice for order-of-magnitude
# tidal ratio checks.
MOON_MU: Final[float] = 4.9048695e12  # lunar gravitational parameter [m³/s²]
SUN_MU: Final[float] = 1.32712440018e20  # solar gravitational parameter [m³/s²]
MOON_MEAN_DISTANCE_M: Final[float] = 3.84400e8  # mean Earth–Moon distance [m]
SUN_MEAN_DISTANCE_M: Final[float] = 1.495978707e11  # 1 AU [m]
SUN_RADIUS_M: Final[float] = 6.957e8  # solar radius [m] — used by EclipseDetector
