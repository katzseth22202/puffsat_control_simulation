"""Tests for pure-Python orbital mechanics helpers."""
import math
from datetime import UTC, datetime

import pytest

from puffsat_sim.config import OrbitalConfig
from puffsat_sim.orbital_math import (
    keplerian_elements,
    keplerian_period,
    orbital_config_from_cities,
    perigee_speed,
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
