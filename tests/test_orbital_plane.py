"""Tests for the great-circle orbital-plane builder."""
import math
from datetime import UTC, datetime

import pytest

from puffsat_sim.config import OrbitalConfig
from puffsat_sim.orbital_plane import orbital_config_from_cities


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
