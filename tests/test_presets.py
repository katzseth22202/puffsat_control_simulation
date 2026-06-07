"""Tests for the force-model presets — content, not build-ladder numbers."""

import pytest

from puffsat_sim import presets
from puffsat_sim.forces import AtmosphericDrag, Geopotential, SolarRadiation, ThirdBody


class TestPresets:
    def test_two_body_is_keplerian(self) -> None:
        cfg = presets.two_body()
        assert cfg.is_keplerian
        assert cfg.perturbations == ()

    def test_j2_zonal_only(self) -> None:
        (geo,) = presets.j2().perturbations
        assert isinstance(geo, Geopotential)
        assert geo.degree == 2
        assert geo.order == 0  # zonal J2, no tesseral terms

    def test_j2_third_body_adds_sun_and_moon(self) -> None:
        cfg = presets.j2_third_body()
        assert {type(p) for p in cfg.perturbations} == {Geopotential, ThirdBody}

    def test_j2_third_body_srp_custom_cr(self) -> None:
        cfg = presets.j2_third_body_srp(cr_area_over_mass=0.05)
        srp = next(p for p in cfg.perturbations if isinstance(p, SolarRadiation))
        assert srp.cr_area_over_mass == pytest.approx(0.05)
        assert not any(isinstance(p, AtmosphericDrag) for p in cfg.perturbations)

    def test_full_force_custom_coefficients(self) -> None:
        cfg = presets.full_force(cr_area_over_mass=0.05, cd_area_over_mass=0.08)
        srp = next(p for p in cfg.perturbations if isinstance(p, SolarRadiation))
        drag = next(p for p in cfg.perturbations if isinstance(p, AtmosphericDrag))
        assert srp.cr_area_over_mass == pytest.approx(0.05)
        assert drag.cd_area_over_mass == pytest.approx(0.08)

    def test_full_force_has_all_four_forces(self) -> None:
        cfg = presets.full_force()
        assert {type(p) for p in cfg.perturbations} == {
            Geopotential,
            ThirdBody,
            SolarRadiation,
            AtmosphericDrag,
        }

    def test_full_force_not_keplerian(self) -> None:
        assert not presets.full_force().is_keplerian
