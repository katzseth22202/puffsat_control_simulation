"""Truth-model validation — the JVM run for :mod:`puffsat_sim.truth_validation` (ADR 0018).

Flies the reference apogee→800 km descent coast three ways and feeds the pure checks: a
numerical two-body coast (Tier 1 integrator health + tolerance-halving) and a numerical J2 coast
cross-checked against the independent pure-Python RK4 Cowell (Tier 2).
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.frames import FramesFactory

from puffsat_sim import mission, presets
from puffsat_sim.descent import COAST_MAX_STEP_S, coast_to_handoff, earth_model, to_absolute_date
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import _REL_TOL, _build_numerical_propagator, build_propagator
from puffsat_sim.truth_validation import (
    TruthValidationFinding,
    conservation_drift,
    format_truth_validation,
    independent_coast,
    max_position_divergence_m,
)

N_SAMPLES: int = 25
CONVERGENCE_TOL_FACTOR: float = 0.1  # tolerance-halving: re-fly at a tenth the relative tolerance


def coast_states(propagator: Any, sample_dates: list[Any], frame: Any) -> NDArray[np.float64]:
    """Sample a coast's EME2000 state vectors at the given dates via its ephemeris."""
    generator = propagator.getEphemerisGenerator()
    propagator.propagate(sample_dates[-1])
    ephemeris = generator.getGeneratedEphemeris()
    rows: list[list[float]] = []
    for date in sample_dates:
        pv = ephemeris.propagate(date).getPVCoordinates(frame)
        p, v = pv.getPosition(), pv.getVelocity()
        rows.append([p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()])
    return np.asarray(rows, dtype=np.float64)


def run_truth_validation(
    orbital_config: Any = mission.NOMINAL_CONFIG, n_samples: int = N_SAMPLES
) -> TruthValidationFinding:
    """Fly the reference coast and assemble the truth-validation finding (Tier 1 + Tier 2)."""
    frame = FramesFactory.getEME2000()
    earth = earth_model()
    epoch = to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)
    orbit = build_propagator(orbital_config, presets.two_body()).getInitialState().getOrbit()

    handoff_date = coast_to_handoff(
        _build_numerical_propagator(orbit, presets.j2(), COAST_MAX_STEP_S), epoch, period, earth
    ).getDate()
    span_s = float(handoff_date.durationFrom(epoch))
    sample_times_s = np.linspace(0.0, span_s, n_samples)
    sample_dates = [epoch.shiftedBy(float(t)) for t in sample_times_s]

    states_2b = coast_states(
        _build_numerical_propagator(orbit, presets.two_body(), COAST_MAX_STEP_S),
        sample_dates,
        frame,
    )
    states_2b_fine = coast_states(
        _build_numerical_propagator(
            orbit, presets.two_body(), COAST_MAX_STEP_S, _REL_TOL * CONVERGENCE_TOL_FACTOR
        ),
        sample_dates,
        frame,
    )
    states_j2 = coast_states(
        _build_numerical_propagator(orbit, presets.j2(), COAST_MAX_STEP_S), sample_dates, frame
    )
    independent = independent_coast(states_j2[0], sample_times_s)

    return TruthValidationFinding(
        conservation=conservation_drift(states_2b),
        convergence_divergence_m=max_position_divergence_m(states_2b, states_2b_fine),
        crosscheck_divergence_m=max_position_divergence_m(states_j2, independent),
        orbit_scale_m=float(np.linalg.norm(states_2b[0, :3])),
        span_s=span_s,
        n_samples=n_samples,
    )


def truth_validation_report() -> None:
    print(format_truth_validation(run_truth_validation()))


if __name__ == "__main__":
    truth_validation_report()
