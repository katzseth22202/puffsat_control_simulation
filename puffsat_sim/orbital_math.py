"""Pure-Python orbital mechanics helpers — no JVM dependency.

Used by truth_model.py and tested independently of Orekit.
"""
from __future__ import annotations

import math
from datetime import UTC, datetime
from typing import Final

from puffsat_sim.config import OrbitalConfig

_EARTH_RADIUS_M: Final[float] = 6_378_137.0   # WGS84 equatorial radius [m]
_WGS84_MU: Final[float] = 3.986_004_418e14    # Earth gravitational parameter [m³/s²]
_J2: Final[float] = 1.08262668e-3             # EGM2008 zonal harmonic J2

# J2000.0 reference epoch (2000-01-01 12:00:00 UTC) for GMST calculation.
_J2000_UTC: Final[datetime] = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)


def keplerian_elements(perigee_alt_m: float, apogee_alt_m: float) -> tuple[float, float]:
    """Return (semi-major axis [m], eccentricity) from altitude above the surface."""
    r_p = _EARTH_RADIUS_M + perigee_alt_m
    r_a = _EARTH_RADIUS_M + apogee_alt_m
    a = (r_p + r_a) / 2.0
    e = (r_a - r_p) / (r_a + r_p)
    return a, e


def keplerian_period(semi_major_axis_m: float) -> float:
    """Return orbital period [s] from semi-major axis [m]."""
    return 2.0 * math.pi * math.sqrt(semi_major_axis_m**3 / _WGS84_MU)


def perigee_speed(semi_major_axis_m: float, perigee_alt_m: float) -> float:
    """Return speed at perigee [m/s] via the vis-viva equation."""
    r_p = _EARTH_RADIUS_M + perigee_alt_m
    return math.sqrt(_WGS84_MU * (2.0 / r_p - 1.0 / semi_major_axis_m))


# ---------------------------------------------------------------------------
# J2 secular perturbation rates (first-order, zonal only)
# ---------------------------------------------------------------------------

def j2_nodal_regression_rate(
    semi_major_axis_m: float, eccentricity: float, inclination_rad: float
) -> float:
    """First-order secular RAAN drift rate due to J2 [rad/s].

    dΩ/dt = -3/2 · n · J2 · (Rₑ/p)² · cos i
    """
    n = math.sqrt(_WGS84_MU / semi_major_axis_m**3)
    p = semi_major_axis_m * (1.0 - eccentricity**2)
    return -1.5 * n * _J2 * (_EARTH_RADIUS_M / p) ** 2 * math.cos(inclination_rad)


def j2_apsidal_precession_rate(
    semi_major_axis_m: float, eccentricity: float, inclination_rad: float
) -> float:
    """First-order secular argument-of-perigee drift rate due to J2 [rad/s].

    dω/dt = 3/4 · n · J2 · (Rₑ/p)² · (5 cos²i − 1)
    """
    n = math.sqrt(_WGS84_MU / semi_major_axis_m**3)
    p = semi_major_axis_m * (1.0 - eccentricity**2)
    return 0.75 * n * _J2 * (_EARTH_RADIUS_M / p) ** 2 * (
        5.0 * math.cos(inclination_rad) ** 2 - 1.0
    )


# ---------------------------------------------------------------------------
# Great-circle orbital plane helpers
# ---------------------------------------------------------------------------

def _unit_vec(lat_deg: float, lon_deg: float) -> tuple[float, float, float]:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    return math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)


def _cross3(
    a: tuple[float, float, float], b: tuple[float, float, float]
) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm3(v: tuple[float, float, float]) -> float:
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def _gmst_rad(epoch: datetime) -> float:
    """Greenwich Mean Sidereal Time at epoch [rad], IAU 1982 approximation.

    Accurate to ~0.1° for dates near J2000; sufficient for setting a
    nominal RAAN from a geographic great-circle plane.
    """
    if epoch.tzinfo is None:
        raise ValueError("epoch must be timezone-aware UTC")
    dt_days = (epoch - _J2000_UTC).total_seconds() / 86400.0
    gmst_deg = (280.46061837 + 360.98564736629 * dt_days) % 360.0
    return math.radians(gmst_deg)


def orbital_config_from_cities(
    lat_a_deg: float,
    lon_a_deg: float,
    lat_b_deg: float,
    lon_b_deg: float,
    epoch: datetime,
    perigee_alt_m: float,
    apogee_alt_m: float,
    arg_of_perigee_rad: float = 0.0,
    mean_anomaly_at_epoch_rad: float = math.pi,
) -> OrbitalConfig:
    """Build an OrbitalConfig whose orbital plane contains two surface points.

    The two geographic coordinates define a great circle (the intersection of
    a plane through Earth's centre with the surface).  That plane becomes the
    orbital plane: inclination and RAAN are derived from it.

    The cities are chosen purely to specify a realistic orbital plane; they are
    otherwise arbitrary.  Ground tracks are not modelled and the simulation is
    not sensitive to which specific locations are used — only the resulting
    inclination and RAAN matter.

    City order determines orbit direction: the satellite travels from A toward B
    (prograde if inclination < 90°).

    epoch must be a timezone-aware UTC datetime.  RAAN is computed by rotating
    the ECEF ascending-node longitude into the EME2000 inertial frame using the
    Greenwich Mean Sidereal Time (IAU 1982) at epoch.
    """
    r_a = _unit_vec(lat_a_deg, lon_a_deg)
    r_b = _unit_vec(lat_b_deg, lon_b_deg)

    pole_raw = _cross3(r_a, r_b)
    mag = _norm3(pole_raw)
    if mag < 1e-12:
        raise ValueError(
            "Cities must not be coincident or antipodal — orbital plane is undefined."
        )
    pole = (pole_raw[0] / mag, pole_raw[1] / mag, pole_raw[2] / mag)

    inclination_rad = math.acos(max(-1.0, min(1.0, pole[2])))

    # Ascending node in ECEF: Ẑ × h_hat (standard definition).
    asc_raw = _cross3((0.0, 0.0, 1.0), pole)
    asc_mag = _norm3(asc_raw)
    if asc_mag < 1e-12:
        # Equatorial orbit: ascending node undefined; RAAN conventionally 0.
        raan_rad = 0.0
    else:
        asc = (asc_raw[0] / asc_mag, asc_raw[1] / asc_mag, asc_raw[2] / asc_mag)
        lon_asc_rad = math.atan2(asc[1], asc[0])
        # Rotate from ECEF longitude to ECI right ascension via GMST.
        raan_rad = (lon_asc_rad + _gmst_rad(epoch)) % (2.0 * math.pi)

    return OrbitalConfig(
        perigee_alt_m=perigee_alt_m,
        apogee_alt_m=apogee_alt_m,
        inclination_rad=inclination_rad,
        raan_rad=raan_rad,
        arg_of_perigee_rad=arg_of_perigee_rad,
        mean_anomaly_at_epoch_rad=mean_anomaly_at_epoch_rad,
        epoch=epoch,
    )
