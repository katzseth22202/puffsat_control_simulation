"""Integration test for the C1 nav-feasibility seam (live JVM, ADR 0012).

Exercises the truth-arc generation and the seeded UKF validation against
full-force Orekit truth on a deliberately tiny arc; the physics conclusions
(q sizing, envelope verdicts) belong to the findings run, not this test.
"""

from __future__ import annotations

import numpy as np
import pytest

from puffsat_sim.nav_feasibility import NavFeasibilitySpec, apogee_state

try:
    # Importing montecarlo boots the JVM and loads orekit-data.zip from the cwd.
    from puffsat_sim.montecarlo import _truth_arc_to_apogee, run_nav_feasibility
except Exception as exc:  # pragma: no cover - environment guard
    pytest.skip(f"Orekit unavailable: {exc}", allow_module_level=True)

from puffsat_sim import presets
from puffsat_sim.config import PhysicsConfig

pytestmark = pytest.mark.integration

_TINY_SPEC = NavFeasibilitySpec(
    arc_duration_s=600.0,
    range_sigma_values_m=(),
    doppler_sigma_values_m_s=(None,),
    cadence_values_hz=(),
    cone_half_angle_values_rad=(),
    n_nodes_values=(),
    q_accel_values_m_s2=(),
)


def _phi() -> np.ndarray:
    phi = np.zeros((3, 6))
    phi[1, 4] = 2.15e5
    return phi


def test_truth_arc_ends_at_the_apogee_state() -> None:
    physics: PhysicsConfig = presets.full_force()
    arc = _truth_arc_to_apogee(0.03, 600.0, physics)

    assert arc.shape == (19, 6)
    apogee = apogee_state()
    # The Orekit truth flies the mission orientation (28.5° inclination) while the
    # pure apogee_state() uses the equatorial convention — the validation is judged
    # in RTN, so only the rotation-invariant facts must agree: the radius and speed
    # at the apogee node, and zero radial velocity there.
    assert np.linalg.norm(arc[-1, :3]) == pytest.approx(np.linalg.norm(apogee[:3]), rel=1e-9)
    assert np.linalg.norm(arc[-1, 3:]) == pytest.approx(np.linalg.norm(apogee[3:]), rel=1e-6)
    radial_speed = arc[-1, :3] @ arc[-1, 3:] / np.linalg.norm(arc[-1, :3])
    assert abs(radial_speed) < 0.5


def test_run_nav_feasibility_returns_envelope_and_validation() -> None:
    result, validations = run_nav_feasibility(_TINY_SPEC, _phi(), seed=20260610)

    assert len(result.outcomes) == 2  # nominal + the range-only cell
    assert len(validations) == 2
    for outcome in validations:
        assert outcome.n_epochs == 18
        assert np.isfinite(outcome.average_nees)
        assert outcome.claimed_t_vel_sigma_m_s < _TINY_SPEC.prior_vel_sigma_m_s
