"""Integration test for the C3c authority + MCC-2 cost sweep (live JVM)."""

from __future__ import annotations

import pytest

try:
    # Importing any JVM-side module boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.runs.authority import run_tail_authority
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

from puffsat_sim.authority import (
    TailAuthoritySpec,
    cheapest_trim,
    trim_dv_to_cover_tail_m_s,
)

pytestmark = pytest.mark.integration


def test_authority_and_trim_curves_have_the_expected_shape() -> None:
    """The C3c slice end-to-end (ADR 0014 decision 5/6), bounds set from the measured sweep.

    A reduced grid (two burn-start altitudes, a high and a low trim node) pins the physics,
    not the platform: the 800 km funnel reproduces C3b's ~480 m (½·a_max·t² at ~246 s),
    raising the burn-start grows it quadratically, and the impulsive trim's lateral
    authority collapses as the node drops toward the §16.6 dead zone.
    """
    spec = TailAuthoritySpec(
        handoff_altitudes_m=(800e3, 1500e3),
        node_altitudes_m=(30_000e3, 1000e3),
        trim_dv_m_s=0.1,
    )
    finding = run_tail_authority(spec)

    # Authority anchor: the descent time from the 800 km hand-off and the funnel it buys
    # reproduce C3b's measured terminal phase.
    a800, a1500 = finding.authority_points
    assert a800.handoff_alt_m == 800e3
    assert 200.0 < a800.t_descent_s < 320.0
    assert 350.0 < a800.radius_m < 700.0
    # A higher burn-start leaves more descent time → a bigger funnel (quadratic).
    assert a1500.t_descent_s > a800.t_descent_s
    assert a1500.radius_m > a800.radius_m

    # The impulsive trim is far cheaper per km from a high node than a low one.
    high, low = finding.trim_points
    assert high.node_alt_m == 30_000e3
    assert high.dv_per_km_m_s < low.dv_per_km_m_s

    # Killing the ~672 m upstream tail from the cheapest node is a sliver of Δv.
    assert trim_dv_to_cover_tail_m_s(finding, cheapest_trim(finding)) < 1.0
