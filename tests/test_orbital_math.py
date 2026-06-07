"""Tests for pure-Python orbital mechanics helpers."""
import math
from datetime import UTC, datetime

import pytest

from puffsat_sim.config import OrbitalConfig
from puffsat_sim.orbital_math import (
    drag_deceleration,
    j2_apsidal_precession_rate,
    j2_nodal_regression_rate,
    keplerian_elements,
    keplerian_period,
    lunar_tidal_ratio,
    orbital_config_from_cities,
    perigee_speed,
    solar_tidal_ratio,
    srp_acceleration,
    std_atm_density,
    tidal_acceleration_ratio,
)

# Reference orbit: 50 km orbit periapsis (debris disposal), 150 000 km apogee
# Interception occurs at 200 km during descent (design doc §3).
_PERIGEE_ALT_M = 50_000.0
_APOGEE_ALT_M = 150_000_000.0
_EARTH_RADIUS_M = 6_378_137.0


class TestKeplerianElements:
    def test_circular_orbit_zero_eccentricity(self) -> None:
        a, e = keplerian_elements(400_000.0, 400_000.0)
        assert e == pytest.approx(0.0, abs=1e-12)
        assert a == pytest.approx(_EARTH_RADIUS_M + 400_000.0)

    def test_reference_orbit_semi_major_axis(self) -> None:
        a, _ = keplerian_elements(_PERIGEE_ALT_M, _APOGEE_ALT_M)
        r_p = _EARTH_RADIUS_M + _PERIGEE_ALT_M
        r_a = _EARTH_RADIUS_M + _APOGEE_ALT_M
        assert a == pytest.approx((r_p + r_a) / 2.0)

    def test_reference_orbit_eccentricity(self) -> None:
        _, e = keplerian_elements(_PERIGEE_ALT_M, _APOGEE_ALT_M)
        # Design doc §3: e ≈ 0.921 for the near-term architecture (50 km periapsis)
        assert e == pytest.approx(0.921033, rel=1e-4)

    def test_eccentricity_bounded(self) -> None:
        _, e = keplerian_elements(_PERIGEE_ALT_M, _APOGEE_ALT_M)
        assert 0.0 <= e < 1.0


class TestKeplerianPeriod:
    def test_reference_orbit_period(self) -> None:
        a, _ = keplerian_elements(_PERIGEE_ALT_M, _APOGEE_ALT_M)
        period = keplerian_period(a)
        # Design doc §3: ~2.68 days (50 km periapsis)
        assert period == pytest.approx(231138.7, rel=1e-4)

    def test_iss_altitude_period(self) -> None:
        # ISS at ~420 km: period ~92 min = 5520 s
        a, _ = keplerian_elements(420_000.0, 420_000.0)
        period = keplerian_period(a)
        assert period == pytest.approx(5554.0, rel=0.01)

    def test_period_scales_with_altitude(self) -> None:
        a_low, _ = keplerian_elements(400_000.0, 400_000.0)
        a_high, _ = keplerian_elements(800_000.0, 800_000.0)
        assert keplerian_period(a_low) < keplerian_period(a_high)


class TestOrbitalConfigFromCities:
    _EPOCH = datetime(2026, 6, 2, 0, 0, 0, tzinfo=UTC)
    _PERIGEE = 50_000.0
    _APOGEE = 150_000_000.0

    def _cfg(self, lat_a: float, lon_a: float, lat_b: float, lon_b: float) -> OrbitalConfig:
        return orbital_config_from_cities(
            lat_a, lon_a, lat_b, lon_b,
            epoch=self._EPOCH,
            perigee_alt_m=self._PERIGEE,
            apogee_alt_m=self._APOGEE,
        )

    def test_tokyo_nyc_inclination_approx_70_deg(self) -> None:
        cfg = self._cfg(35.6762, 139.6503, 40.7128, -74.0060)
        assert math.degrees(cfg.inclination_rad) == pytest.approx(70.0, abs=0.5)

    def test_equatorial_cities_inclination_zero(self) -> None:
        cfg = self._cfg(0.0, 0.0, 0.0, 90.0)
        assert math.degrees(cfg.inclination_rad) == pytest.approx(0.0, abs=1e-9)

    def test_pole_equator_inclination_90_deg(self) -> None:
        cfg = self._cfg(0.0, 0.0, 90.0, 0.0)
        assert math.degrees(cfg.inclination_rad) == pytest.approx(90.0, abs=1e-9)

    def test_periapsis_apoapsis_preserved(self) -> None:
        cfg = self._cfg(35.6762, 139.6503, 40.7128, -74.0060)
        assert isinstance(cfg, OrbitalConfig)
        assert cfg.perigee_alt_m == self._PERIGEE
        assert cfg.apogee_alt_m == self._APOGEE
        assert cfg.mean_anomaly_at_epoch_rad == pytest.approx(math.pi)
        assert cfg.epoch == self._EPOCH

    def test_raan_in_range(self) -> None:
        cfg = self._cfg(35.6762, 139.6503, 40.7128, -74.0060)
        assert 0.0 <= cfg.raan_rad < 2 * math.pi

    def test_coincident_cities_raises(self) -> None:
        with pytest.raises(ValueError, match="coincident or antipodal"):
            orbital_config_from_cities(
                35.0, 139.0, 35.0, 139.0,
                epoch=self._EPOCH,
                perigee_alt_m=self._PERIGEE,
                apogee_alt_m=self._APOGEE,
            )

    def test_naive_epoch_raises(self) -> None:
        with pytest.raises(ValueError, match="timezone-aware"):
            orbital_config_from_cities(
                35.0, 139.0, 40.0, -74.0,
                epoch=datetime(2026, 6, 2),  # naive
                perigee_alt_m=self._PERIGEE,
                apogee_alt_m=self._APOGEE,
            )


class TestJ2Rates:
    """Verify J2 secular rate formulas against known values for the reference orbit.

    Reference orbit: 50 km periapsis, 150 000 km apogee, ~70° inclination.
    At i ≈ 70°, 5cos²i − 1 ≈ 5·(0.342)² − 1 ≈ −0.415, so apsidal precession
    is retrograde.  Nodal regression is always retrograde for prograde orbits.

    Values derived from the J2 analytic formulas with the same constants used
    in orbital_math.py (_J2, _WGS84_MU, _EARTH_RADIUS_M).
    """

    _A, _E = keplerian_elements(_PERIGEE_ALT_M, _APOGEE_ALT_M)
    _PERIOD = keplerian_period(_A)
    _I_70 = math.radians(70.0)
    _I_63_4 = math.radians(63.435)  # critical inclination: apsidal precession ≈ 0

    def test_nodal_regression_retrograde_prograde(self) -> None:
        rate = j2_nodal_regression_rate(self._A, self._E, self._I_70)
        assert rate < 0.0, "Nodal regression must be retrograde (negative) for prograde orbit"

    def test_apsidal_precession_retrograde_at_70_deg(self) -> None:
        rate = j2_apsidal_precession_rate(self._A, self._E, self._I_70)
        assert rate < 0.0, "Apsidal precession must be retrograde at i=70° (5cos²i−1 < 0)"

    def test_apsidal_precession_zero_at_critical_inclination(self) -> None:
        rate = j2_apsidal_precession_rate(self._A, self._E, self._I_63_4)
        assert abs(rate) < 1e-12, "Apsidal precession must vanish at critical inclination ~63.4°"

    def test_apsidal_precession_prograde_at_low_inclination(self) -> None:
        rate = j2_apsidal_precession_rate(self._A, self._E, math.radians(28.5))
        assert rate > 0.0, "Apsidal precession must be prograde at i=28.5° (5cos²i−1 > 0)"

    def test_nodal_drift_per_period_order_of_magnitude(self) -> None:
        # Reference orbit i≈70°: expect ΔRAAN ≈ −0.04° to −0.07° per period
        rate = j2_nodal_regression_rate(self._A, self._E, self._I_70)
        d_raan_deg = math.degrees(rate * self._PERIOD)
        assert -0.10 < d_raan_deg < -0.02

    def test_apsidal_drift_per_period_order_of_magnitude(self) -> None:
        # Reference orbit i≈70°: expect Δω ≈ −0.02° to −0.05° per period
        rate = j2_apsidal_precession_rate(self._A, self._E, self._I_70)
        d_omega_deg = math.degrees(rate * self._PERIOD)
        assert -0.06 < d_omega_deg < -0.01

    def test_equatorial_orbit_max_nodal_regression(self) -> None:
        # i=0: nodal regression is most negative; cos(0)=1, maximum magnitude
        rate_eq = j2_nodal_regression_rate(self._A, self._E, 0.0)
        rate_70 = j2_nodal_regression_rate(self._A, self._E, self._I_70)
        assert rate_eq < rate_70, "Equatorial orbit must have faster nodal regression than 70°"

    def test_polar_orbit_zero_nodal_regression(self) -> None:
        rate = j2_nodal_regression_rate(self._A, self._E, math.radians(90.0))
        assert abs(rate) < 1e-16, "Polar orbit (i=90°) must have zero nodal regression"


class TestDragDeceleration:
    """Verify piecewise-exponential atmosphere and drag deceleration helper.

    Calibrated to NRLMSISE-00 at moderate solar activity (F10.7≈150, Ap≈15).
    Design doc §4: drag "bites below ~300-400 km."
    """

    _CD_AM = 0.04    # rung_2d default Cd·(A/m) [m²/kg]
    _V_10KMS = 10_000.0

    def test_density_decreases_with_altitude(self) -> None:
        assert std_atm_density(200_000) > std_atm_density(300_000)
        assert std_atm_density(300_000) > std_atm_density(500_000)

    def test_density_at_200km_matches_nrlmsise(self) -> None:
        # NRLMSISE-00 at 200 km, F10.7=150, Ap=15 ≈ 2.5e-10 kg/m³
        rho = std_atm_density(200_000)
        assert 1e-10 < rho < 5e-10

    def test_density_at_surface_approx_1kg_m3(self) -> None:
        assert std_atm_density(0) == pytest.approx(1.225, rel=0.01)

    def test_drag_zero_cd_am(self) -> None:
        assert drag_deceleration(0.0, self._V_10KMS, 200_000) == pytest.approx(0.0)

    def test_drag_scales_with_speed_squared(self) -> None:
        a1 = drag_deceleration(self._CD_AM, 10_000.0, 200_000)
        a2 = drag_deceleration(self._CD_AM, 20_000.0, 200_000)
        assert a2 == pytest.approx(4.0 * a1, rel=1e-9)

    def test_drag_scales_with_cd_am(self) -> None:
        a1 = drag_deceleration(0.04, self._V_10KMS, 200_000)
        a2 = drag_deceleration(0.08, self._V_10KMS, 200_000)
        assert a2 == pytest.approx(2.0 * a1, rel=1e-9)

    def test_drag_at_200km_order_of_magnitude(self) -> None:
        # a_drag = 0.5 * 2.5e-10 * (10000)^2 * 0.04 ≈ 5e-4 m/s²
        a = drag_deceleration(self._CD_AM, self._V_10KMS, 200_000)
        assert 1e-4 < a < 5e-3

    def test_drag_at_300km_much_smaller_than_200km(self) -> None:
        a_200 = drag_deceleration(self._CD_AM, self._V_10KMS, 200_000)
        a_300 = drag_deceleration(self._CD_AM, self._V_10KMS, 300_000)
        # Density ratio ~30×, so drag ratio should be ~30×
        assert a_200 / a_300 > 10

    def test_drag_at_1000km_negligible(self) -> None:
        # Above ~800 km drag is negligible (<1e-8 m/s²) for any reasonable Cd·A/m
        a = drag_deceleration(self._CD_AM, self._V_10KMS, 1_000_000)
        assert a < 1e-8


class TestSrpAcceleration:
    """Verify SRP acceleration helper against the design doc benchmark (Rung 2c).

    At 1 AU with Cr·A/m = 0.02 m²/kg (PhysicsConfig.rung_2c default):
        a = 4.56e-6 Pa × 0.02 m²/kg ≈ 9.12e-8 m/s²
    This is ~9.3× larger than the Moon's tidal acceleration at apogee, making
    SRP the dominant non-gravitational force in the coast phase.
    """

    _CR_AM = 0.02         # rung_2c default [m²/kg]
    _ONE_AU = 1.495978707e11

    def test_default_at_1au(self) -> None:
        a = srp_acceleration(self._CR_AM)
        assert a == pytest.approx(4.56e-6 * self._CR_AM, rel=1e-9)

    def test_zero_cr_am(self) -> None:
        assert srp_acceleration(0.0) == pytest.approx(0.0)

    def test_scales_inverse_square_with_distance(self) -> None:
        a1 = srp_acceleration(self._CR_AM, self._ONE_AU)
        a2 = srp_acceleration(self._CR_AM, 2 * self._ONE_AU)
        assert a1 == pytest.approx(4.0 * a2, rel=1e-9)

    def test_order_of_magnitude_at_1au(self) -> None:
        # Expected ~9e-8 m/s² for Cr·A/m = 0.02
        a = srp_acceleration(self._CR_AM)
        assert 5e-8 < a < 2e-7

    def test_increases_closer_to_sun(self) -> None:
        a_near = srp_acceleration(self._CR_AM, 0.5 * self._ONE_AU)
        a_far = srp_acceleration(self._CR_AM, 2.0 * self._ONE_AU)
        assert a_near > a_far


class TestTidalAccelerationRatio:
    """Verify the Hill-approximation tidal ratio helper against design doc benchmarks.

    Design doc §2: at 150 000 km apogee, solar-tidal acceleration is ~0.1% of
    Earth's gravity.  Moon and Sun combined should be a few tenths of a percent.
    """

    _MOON_MU = 4.9048695e12
    _SUN_MU = 1.32712440018e20
    _MOON_DIST = 3.84400e8
    _SUN_DIST = 1.495978707e11

    def test_moon_ratio_at_reference_apogee_order_of_magnitude(self) -> None:
        # Moon tidal ratio at 150 000 km apogee should be ~0.1–0.3%
        ratio = lunar_tidal_ratio(_APOGEE_ALT_M)
        assert 0.001 < ratio < 0.003

    def test_sun_ratio_at_reference_apogee_order_of_magnitude(self) -> None:
        # Sun tidal ratio at 150 000 km apogee should be ~0.05–0.15%
        ratio = solar_tidal_ratio(_APOGEE_ALT_M)
        assert 0.0005 < ratio < 0.002

    def test_moon_dominates_sun_at_150000_km(self) -> None:
        assert lunar_tidal_ratio(_APOGEE_ALT_M) > solar_tidal_ratio(_APOGEE_ALT_M)

    def test_ratio_increases_with_apogee_altitude(self) -> None:
        low = tidal_acceleration_ratio(100_000_000.0, self._MOON_MU, self._MOON_DIST)
        high = tidal_acceleration_ratio(200_000_000.0, self._MOON_MU, self._MOON_DIST)
        assert high > low

    def test_ratio_positive(self) -> None:
        assert tidal_acceleration_ratio(_APOGEE_ALT_M, self._MOON_MU, self._MOON_DIST) > 0.0

    def test_low_altitude_ratio_negligible(self) -> None:
        # At 400 km LEO apogee the Moon's tidal effect is completely negligible (<1e-6)
        ratio = tidal_acceleration_ratio(400_000.0, self._MOON_MU, self._MOON_DIST)
        assert ratio < 1e-6


class TestPerigeeSpeed:
    def test_reference_orbit_perigee_speed(self) -> None:
        a, _ = keplerian_elements(_PERIGEE_ALT_M, _APOGEE_ALT_M)
        v = perigee_speed(a, _PERIGEE_ALT_M)
        # Design doc §3: ~10.91 km/s at 50 km periapsis
        assert v == pytest.approx(10914.2, rel=1e-3)

    def test_circular_orbit_speed(self) -> None:
        # Circular orbit: vis-viva gives v = sqrt(mu/r)
        alt = 500_000.0
        a, _ = keplerian_elements(alt, alt)
        v = perigee_speed(a, alt)
        r = _EARTH_RADIUS_M + alt
        mu = 3.986_004_418e14
        assert v == pytest.approx(math.sqrt(mu / r), rel=1e-9)
