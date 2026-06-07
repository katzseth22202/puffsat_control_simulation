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

# Solar radiation pressure at 1 AU [Pa = N/m²] — used for analytic SRP estimates.
_SRP_P0_PA: Final[float] = 4.56e-6

# Third-body constants for analytic perturbation estimates.
# Distances are mean values; instantaneous geometry varies but these are
# sufficient for order-of-magnitude tidal ratio checks.
_MOON_MU: Final[float] = 4.9048695e12         # lunar gravitational parameter [m³/s²]
_SUN_MU: Final[float] = 1.32712440018e20      # solar gravitational parameter [m³/s²]
_MOON_MEAN_DISTANCE_M: Final[float] = 3.84400e8   # mean Earth–Moon distance [m]
_SUN_MEAN_DISTANCE_M: Final[float] = 1.495978707e11  # 1 AU [m]

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


def wrap_to_pi(angle_rad: float) -> float:
    """Wrap an angle [rad] to (-π, π].

    Used to take the true signed difference of two angles that may straddle the
    0/2π branch cut — e.g. a small secular drift that Orekit reports as the raw
    difference ~2π − ε instead of ~−ε.
    """
    wrapped = (angle_rad + math.pi) % (2.0 * math.pi) - math.pi
    # (a + π) % 2π lands in [0, 2π); the subtraction yields [-π, π).  Map the
    # -π endpoint to +π so the interval is the conventional (-π, π].
    return math.pi if wrapped == -math.pi else wrapped


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
# Third-body tidal perturbation strength (analytic, Hill approximation)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Atmospheric drag (analytic, piecewise-exponential density estimate)
# ---------------------------------------------------------------------------

# (base_altitude_m, density_kg_m3, scale_height_m) — calibrated to NRLMSISE-00
# at moderate solar activity (F10.7≈150, Ap≈15).  Good to ~factor-of-2 for
# order-of-magnitude checks; not a substitute for the full model in propagator.py.
_STD_ATM_LAYERS: Final[tuple[tuple[float, float, float], ...]] = (
    (0,        1.225,    8_500),
    (25_000,   3.9e-2,   6_700),
    (50_000,   1.0e-3,   7_200),
    (80_000,   1.0e-5,   9_000),
    (100_000,  5.6e-7,  15_500),
    (150_000,  2.2e-9,  22_000),
    (200_000,  2.5e-10, 29_000),
    (300_000,  8.0e-12, 37_000),
    (500_000,  5.5e-13, 60_000),
    (700_000,  3.6e-14, 73_000),
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


def drag_deceleration(
    cd_area_over_mass: float, speed_m_s: float, altitude_m: float
) -> float:
    """Drag deceleration magnitude [m/s²] using piecewise-exponential density.

    a_drag = ½ · ρ(h) · v² · (Cd·A/m)
    """
    return 0.5 * std_atm_density(altitude_m) * speed_m_s**2 * cd_area_over_mass


def srp_acceleration(
    cr_area_over_mass: float,
    sun_distance_m: float = _SUN_MEAN_DISTANCE_M,
) -> float:
    """SRP acceleration magnitude [m/s²] at a given distance from the Sun.

    a_srp = P₀ · (d₀/r)² · (Cr·A/m)  where P₀ = 4.56×10⁻⁶ Pa at 1 AU.
    """
    return _SRP_P0_PA * (_SUN_MEAN_DISTANCE_M / sun_distance_m) ** 2 * cr_area_over_mass


def tidal_acceleration_ratio(
    apogee_alt_m: float, body_mu_m3_s2: float, body_distance_m: float
) -> float:
    """Ratio of third-body tidal acceleration to Earth monopole gravity at apogee.

    Uses the Hill (tidal) approximation: a_tidal ≈ 2·μ_body·r_apogee / d_body³.
    Ratio = a_tidal / a_earth = 2·μ_body·r_apogee³ / (d_body³·μ_earth).

    For the reference orbit (apogee 150 000 km) the Moon gives ~0.17% and
    the Sun gives ~0.08%, consistent with the design doc "~0.1%" benchmark.
    """
    r_a = _EARTH_RADIUS_M + apogee_alt_m
    a_tidal = 2.0 * body_mu_m3_s2 * r_a / body_distance_m**3
    a_earth = _WGS84_MU / r_a**2
    return a_tidal / a_earth


def lunar_tidal_ratio(apogee_alt_m: float) -> float:
    """Tidal acceleration ratio for the Moon at mean distance."""
    return tidal_acceleration_ratio(apogee_alt_m, _MOON_MU, _MOON_MEAN_DISTANCE_M)


def solar_tidal_ratio(apogee_alt_m: float) -> float:
    """Tidal acceleration ratio for the Sun at 1 AU."""
    return tidal_acceleration_ratio(apogee_alt_m, _SUN_MU, _SUN_MEAN_DISTANCE_M)


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
