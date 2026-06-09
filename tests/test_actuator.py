"""Tests for the pure B1 actuator — Tsiolkovsky burn kinematics, no JVM."""

import math

import pytest

from puffsat_sim.actuator import Actuator, plan_burn
from puffsat_sim.constants import STANDARD_GRAVITY_M_S2


def test_burn_duration_is_dv_times_mass_over_thrust() -> None:
    """A commanded Δv maps to a constant-mass burn of duration Δv·m/F (Isp-free)."""
    actuator = Actuator(isp_s=50.0, max_thrust_n=0.4, wet_mass_kg=25.0)
    burn = plan_burn(actuator, (0.5, 0.0, 0.0))
    assert burn.duration_s == pytest.approx(0.5 * 25.0 / 0.4)


def test_propellant_is_tsiolkovsky_at_the_spec_isp() -> None:
    """Propellant is the Tsiolkovsky mass m·(1 − e^{−Δv/Isp·g₀}) at the actuator's Isp."""
    actuator = Actuator(isp_s=50.0, max_thrust_n=0.4, wet_mass_kg=25.0)
    burn = plan_burn(actuator, (0.0, 0.6, 0.0))
    dv = 0.6
    expected = 25.0 * (1.0 - math.exp(-dv / (50.0 * STANDARD_GRAVITY_M_S2)))
    assert burn.propellant_kg == pytest.approx(expected)


def test_small_dv_propellant_matches_the_linear_isp_model() -> None:
    """In the small-Δv limit the Tsiolkovsky mass collapses to the ADR 0004 ledger Δv/(Isp·g₀)."""
    actuator = Actuator(isp_s=200.0, max_thrust_n=0.4, wet_mass_kg=25.0)
    dv = 0.05
    burn = plan_burn(actuator, (dv, 0.0, 0.0))
    linear = 25.0 * dv / (200.0 * STANDARD_GRAVITY_M_S2)
    assert burn.propellant_kg == pytest.approx(linear, rel=1e-4)


def test_zero_dv_is_no_burn() -> None:
    """A null command costs no time and no propellant (the open-loop / converged-to-zero case)."""
    burn = plan_burn(Actuator(isp_s=50.0), (0.0, 0.0, 0.0))
    assert burn.duration_s == 0.0
    assert burn.propellant_kg == 0.0
