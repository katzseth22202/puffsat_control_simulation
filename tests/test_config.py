"""Tests for OrbitalConfig and PhysicsConfig dataclasses."""
import dataclasses
import math
from datetime import UTC, datetime

import pytest

from puffsat_sim.config import OrbitalConfig, PhysicsConfig

_EPOCH = datetime(2026, 6, 2, tzinfo=UTC)


class TestPhysicsConfigIsKeplerian:
    def test_default_is_keplerian(self) -> None:
        assert PhysicsConfig().is_keplerian

    def test_rung_keplerian_is_keplerian(self) -> None:
        assert PhysicsConfig.rung_keplerian().is_keplerian

    def test_nonzero_geopotential_not_keplerian(self) -> None:
        assert not PhysicsConfig(geopotential_degree=2).is_keplerian

    def test_third_body_not_keplerian(self) -> None:
        assert not PhysicsConfig(third_body=True).is_keplerian

    def test_srp_enabled_not_keplerian(self) -> None:
        assert not PhysicsConfig(srp_cr_area_over_mass=0.02).is_keplerian

    def test_drag_enabled_not_keplerian(self) -> None:
        assert not PhysicsConfig(drag_cd_area_over_mass=0.04).is_keplerian


class TestPhysicsConfigClassMethods:
    def test_rung_2a_geopotential_only(self) -> None:
        cfg = PhysicsConfig.rung_2a()
        assert cfg.geopotential_degree == 2
        assert not cfg.third_body
        assert cfg.srp_cr_area_over_mass is None
        assert cfg.drag_cd_area_over_mass is None

    def test_rung_2b_adds_third_body(self) -> None:
        cfg = PhysicsConfig.rung_2b()
        assert cfg.geopotential_degree == 2
        assert cfg.third_body

    def test_rung_2c_custom_cr(self) -> None:
        cfg = PhysicsConfig.rung_2c(cr_area_over_mass=0.05)
        assert cfg.srp_cr_area_over_mass == pytest.approx(0.05)
        assert cfg.third_body
        assert cfg.drag_cd_area_over_mass is None

    def test_rung_2d_full_model(self) -> None:
        cfg = PhysicsConfig.rung_2d(cr_area_over_mass=0.05, cd_area_over_mass=0.08)
        assert cfg.srp_cr_area_over_mass == pytest.approx(0.05)
        assert cfg.drag_cd_area_over_mass == pytest.approx(0.08)
        assert cfg.third_body
        assert cfg.geopotential_degree == 2

    def test_rung_2d_not_keplerian(self) -> None:
        assert not PhysicsConfig.rung_2d().is_keplerian


class TestOrbitalConfigFrozen:
    def _make(self) -> OrbitalConfig:
        return OrbitalConfig(
            perigee_alt_m=50_000.0,
            apogee_alt_m=150_000_000.0,
            inclination_rad=math.radians(70.0),
            raan_rad=0.0,
            arg_of_perigee_rad=0.0,
            mean_anomaly_at_epoch_rad=math.pi,
            epoch=_EPOCH,
        )

    def test_fields_set_correctly(self) -> None:
        cfg = self._make()
        assert cfg.perigee_alt_m == pytest.approx(50_000.0)
        assert cfg.apogee_alt_m == pytest.approx(150_000_000.0)
        assert cfg.inclination_rad == pytest.approx(math.radians(70.0))
        assert cfg.mean_anomaly_at_epoch_rad == pytest.approx(math.pi)
        assert cfg.epoch == _EPOCH

    def test_immutable(self) -> None:
        cfg = self._make()
        with pytest.raises(dataclasses.FrozenInstanceError):
            cfg.inclination_rad = 1.0  # type: ignore[misc]
