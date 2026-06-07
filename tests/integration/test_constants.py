"""Pin the pure-Python physical constants to Orekit's WGS84 values.

The analytic signatures use puffsat_sim.constants; the truth propagation uses
Orekit's org.orekit.utils.Constants.  If the two drift apart, analytic
cross-checks would pass or fail spuriously.  This test fails loudly instead.
"""

from __future__ import annotations

import pytest

from puffsat_sim import constants

try:
    import puffsat_sim.jvm  # noqa: F401  boots the JVM

    from org.orekit.utils import Constants
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

pytestmark = pytest.mark.integration


def test_earth_radius_matches_orekit() -> None:
    assert constants.EARTH_RADIUS_M == pytest.approx(float(Constants.WGS84_EARTH_EQUATORIAL_RADIUS))


def test_earth_mu_matches_orekit() -> None:
    assert constants.WGS84_MU == pytest.approx(float(Constants.WGS84_EARTH_MU))
