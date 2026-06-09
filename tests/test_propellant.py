"""Tests for the pure B2 propellant ledger — Δv → fraction-vs-Isp, no JVM."""

import pytest

from puffsat_sim.constants import STANDARD_GRAVITY_M_S2
from puffsat_sim.propellant import propellant_curve, propellant_fraction
from puffsat_sim.sweep import ISP_ANCHORS_S, budget_dv_m_s


def test_propellant_fraction_is_the_inverse_of_budget_dv() -> None:
    """Fraction = Δv/(Isp·g₀) (ADR 0004 linear): the exact inverse of the A3 budget Δv."""
    isp = 50.0
    # the Δv that the 2% budget buys at 50 s must read back as exactly 2%.
    assert propellant_fraction(budget_dv_m_s(isp), isp) == pytest.approx(0.02)
    assert propellant_fraction(1.0, 50.0) == pytest.approx(1.0 / (50.0 * STANDARD_GRAVITY_M_S2))


def test_propellant_curve_flags_each_isp_anchor_against_the_2pct_line() -> None:
    """The curve reports the fraction at each Isp anchor and whether it clears the 2% budget."""
    points = propellant_curve(15.0)  # over budget at 50/70 s, under at 200 s
    by_isp = {p.isp_s: p for p in points}
    assert set(by_isp) == set(ISP_ANCHORS_S)
    assert by_isp[50.0].fraction == pytest.approx(15.0 / (50.0 * STANDARD_GRAVITY_M_S2))
    assert not by_isp[50.0].within_budget
    assert not by_isp[70.0].within_budget
    assert by_isp[200.0].within_budget
