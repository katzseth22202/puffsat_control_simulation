"""Constant-F10.7/Ap space-weather provider for NRLMSISE-00 (JVM boundary).

Implements Orekit's ``NRLMSISE00InputParameters`` in Python via JPype so the truth
model's drag density is driven by the ``AtmosphericDrag`` spec's ``f10p7`` / ``ap``
— a per-run Monte Carlo draw (design doc §16.7) — instead of calendar-tied CSSI
data.  This is the one place a Java interface is implemented in Python: JPype's
``@JOverride`` surfaces as an untyped decorator and the override methods must keep
Orekit's camelCase Java names, so this module alone is exempted from
``untyped-decorator`` (mypy) and ``N802`` (ruff) in ``pyproject.toml``; everything
else stays strict.
"""

from __future__ import annotations

from typing import Any

from jpype import JImplements, JOverride  # type: ignore[attr-defined]  # absent from JPype stubs

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.models.earth.atmosphere import NRLMSISE00InputParameters
from org.orekit.time import AbsoluteDate


@JImplements(NRLMSISE00InputParameters)
class _ConstantSpaceWeather:
    """Returns the same F10.7 (daily = 81-day average) and Ap at every date."""

    def __init__(self, f10p7: float, ap: float) -> None:
        self._f10p7 = f10p7
        self._ap = ap

    @JOverride
    def getMinDate(self) -> Any:
        return AbsoluteDate.PAST_INFINITY

    @JOverride
    def getMaxDate(self) -> Any:
        return AbsoluteDate.FUTURE_INFINITY

    @JOverride
    def getDailyFlux(self, date: Any) -> float:
        return self._f10p7

    @JOverride
    def getAverageFlux(self, date: Any) -> float:
        return self._f10p7

    @JOverride
    def getAp(self, date: Any) -> Any:
        # NRLMSISE-00 expects 7 Ap values (daily + 3-hourly history); flat here.
        return [self._ap] * 7


def constant_space_weather(f10p7: float, ap: float) -> Any:
    """An NRLMSISE00InputParameters with constant flux/Ap.

    Returned as Any: the JPype proxy is not visible to mypy as the Java interface,
    so callers (forces/build.py) stay clear of a spurious call-overload error.
    """
    return _ConstantSpaceWeather(f10p7, ap)
