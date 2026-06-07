"""Reference mission scenario — the nominal PuffSat orbit under study.

Pure Python.  These fixed mission parameters are defined here once and shared by
the truth-model runner and its tests, so the reference orbit lives in a single
place.  Monte Carlo injection-state dispersion (Rung D) samples around
NOMINAL_CONFIG; the altitudes below are not varied.
"""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Final

from puffsat_sim.config import OrbitalConfig
from puffsat_sim.orbital_plane import orbital_config_from_cities

PERIGEE_ALT_M: Final[float] = 50_000.0  # orbit periapsis [m]; debris disposal below Kármán
APOGEE_ALT_M: Final[float] = 150_000_000.0  # deployment apogee [m]
INTERCEPTION_ALT_M: Final[float] = 200_000.0  # control target: 200 km descent crossing

# Nominal orbital plane — a great circle through two surface points.  The cities
# are arbitrary: they only fix a realistic mid-to-high inclination (~70°) for the
# perturbation study.  Ground tracks are not modelled and the simulation is not
# sensitive to which locations are used — only the resulting inclination and RAAN
# matter.  The epoch sets the RAAN via GMST.
_TOKYO: Final[tuple[float, float]] = (35.6762, 139.6503)  # (lat°N, lon°E)
_NEW_YORK: Final[tuple[float, float]] = (40.7128, -74.0060)  # (lat°N, lon°W)
EPOCH: Final[datetime] = datetime(2026, 6, 2, 0, 0, 0, tzinfo=UTC)

NOMINAL_CONFIG: Final[OrbitalConfig] = orbital_config_from_cities(
    *_TOKYO,
    *_NEW_YORK,
    epoch=EPOCH,
    perigee_alt_m=PERIGEE_ALT_M,
    apogee_alt_m=APOGEE_ALT_M,
)
