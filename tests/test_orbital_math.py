"""Tests for pure-Python orbital mechanics helpers."""
import math

import pytest

from puffsat_sim.orbital_math import keplerian_elements, keplerian_period, perigee_speed

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
