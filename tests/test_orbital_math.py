"""Tests for foundational two-body orbital mechanics helpers."""

import math

import pytest

from puffsat_sim.mission import APOGEE_ALT_M, PERIGEE_ALT_M
from puffsat_sim.orbital_math import (
    keplerian_elements,
    keplerian_period,
    perigee_speed,
    wrap_to_pi,
)

# Independent literal (not imported from constants) so these tests verify
# keplerian_elements rather than re-using the same radius the code uses.
_EARTH_RADIUS_M = 6_378_137.0


class TestKeplerianElements:
    def test_circular_orbit_zero_eccentricity(self) -> None:
        a, e = keplerian_elements(400_000.0, 400_000.0)
        assert e == pytest.approx(0.0, abs=1e-12)
        assert a == pytest.approx(_EARTH_RADIUS_M + 400_000.0)

    def test_reference_orbit_semi_major_axis(self) -> None:
        a, _ = keplerian_elements(PERIGEE_ALT_M, APOGEE_ALT_M)
        r_p = _EARTH_RADIUS_M + PERIGEE_ALT_M
        r_a = _EARTH_RADIUS_M + APOGEE_ALT_M
        assert a == pytest.approx((r_p + r_a) / 2.0)

    def test_reference_orbit_eccentricity(self) -> None:
        _, e = keplerian_elements(PERIGEE_ALT_M, APOGEE_ALT_M)
        # Design doc §3: e ≈ 0.921 for the near-term architecture (50 km periapsis)
        assert e == pytest.approx(0.921033, rel=1e-4)

    def test_eccentricity_bounded(self) -> None:
        _, e = keplerian_elements(PERIGEE_ALT_M, APOGEE_ALT_M)
        assert 0.0 <= e < 1.0


class TestKeplerianPeriod:
    def test_reference_orbit_period(self) -> None:
        a, _ = keplerian_elements(PERIGEE_ALT_M, APOGEE_ALT_M)
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


class TestWrapToPi:
    """Verify angle wrapping used by the J2 secular-drift validation.

    Orekit reports osculating RAAN/ω on a 0/2π branch, so a small retrograde
    drift surfaces as a raw difference near +2π.  wrap_to_pi must recover the
    true signed value.
    """

    def test_zero_unchanged(self) -> None:
        assert wrap_to_pi(0.0) == pytest.approx(0.0)

    def test_small_drift_reported_near_2pi(self) -> None:
        # A −0.03 rad drift surfacing as +2π−0.03 must wrap back to −0.03.
        raw = 2.0 * math.pi - 0.03
        assert wrap_to_pi(raw) == pytest.approx(-0.03, abs=1e-12)

    def test_small_positive_unchanged(self) -> None:
        assert wrap_to_pi(0.05) == pytest.approx(0.05, abs=1e-12)

    def test_pi_maps_to_plus_pi(self) -> None:
        # Endpoint convention: interval is (-π, π], so −π folds to +π.
        assert wrap_to_pi(math.pi) == pytest.approx(math.pi)
        assert wrap_to_pi(-math.pi) == pytest.approx(math.pi)

    def test_result_in_canonical_interval(self) -> None:
        for k in range(-4, 5):
            for base in (0.1, 1.0, 2.5, 3.1):
                w = wrap_to_pi(base + 2.0 * math.pi * k)
                assert -math.pi < w <= math.pi + 1e-12

    def test_idempotent_on_wrapped_value(self) -> None:
        for raw in (0.3, -0.3, 3.0, -3.0, 6.0, -6.0):
            once = wrap_to_pi(raw)
            assert wrap_to_pi(once) == pytest.approx(once, abs=1e-12)


class TestPerigeeSpeed:
    def test_reference_orbit_perigee_speed(self) -> None:
        a, _ = keplerian_elements(PERIGEE_ALT_M, APOGEE_ALT_M)
        v = perigee_speed(a, PERIGEE_ALT_M)
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
