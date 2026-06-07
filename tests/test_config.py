"""Tests for OrbitalConfig and PhysicsConfig dataclasses."""
import dataclasses
import math
from datetime import UTC, datetime

import pytest

from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.forces import Geopotential, SolarRadiation, ThirdBody

_EPOCH = datetime(2026, 6, 2, tzinfo=UTC)


class TestPhysicsConfigIsKeplerian:
    def test_default_is_keplerian(self) -> None:
        assert PhysicsConfig().is_keplerian

    def test_empty_perturbations_is_keplerian(self) -> None:
        assert PhysicsConfig(()).is_keplerian

    def test_geopotential_not_keplerian(self) -> None:
        assert not PhysicsConfig((Geopotential(degree=2),)).is_keplerian

    def test_third_body_not_keplerian(self) -> None:
        assert not PhysicsConfig((ThirdBody(),)).is_keplerian

    def test_srp_not_keplerian(self) -> None:
        assert not PhysicsConfig((SolarRadiation(cr_area_over_mass=0.02),)).is_keplerian


class TestPhysicsConfigPerturbations:
    def test_holds_perturbations_in_order(self) -> None:
        cfg = PhysicsConfig((Geopotential(degree=2), ThirdBody()))
        assert cfg.perturbations == (Geopotential(degree=2), ThirdBody())

    def test_immutable(self) -> None:
        cfg = PhysicsConfig((Geopotential(degree=2),))
        with pytest.raises(dataclasses.FrozenInstanceError):
            cfg.perturbations = ()  # type: ignore[misc]


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
