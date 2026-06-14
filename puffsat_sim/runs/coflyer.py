"""Co-flyer phasing feasibility — the JVM run for the ADR 0019 Lever-2 gate.

Lever 2 reuses the already-on-orbit launch rocket as a *close* terminal tracker.  A small apogee
maneuver — raise perigee, lower apogee by the same amount — keeps the **semi-major axis constant**,
so the rocket holds the PuffSat train's period and stays **phase-locked** (no secular drift).  This
run flies that modified rocket orbit alongside the nominal descent and measures, through the PuffSat
800→200 km terminal window, the rocket↔train-centroid range and the rocket's own altitude — feeding
the pure verdict :func:`puffsat_sim.tracker_fusion.phasing_verdict`:

* **range** — does the rocket stay within the ~500 km design range the fusion credits its σ_θ·R
  advantage at?
* **GNSS volume** — does the rocket stay below the GPS constellation, where an unlocked spaceborne
  receiver pins its rocket→target vector independently of the long-baseline angle?

The rocket orbit is *constructed directly* at the constant-a / lower-e elements (co-phased in mean
anomaly with the train): the phasing question is geometric, so the maneuver that realizes those
elements is a separate propulsion detail.  Physics is J2 — the geometry of a ~32 h coast, not the
near-perigee drag (a feasibility-of-phasing check, consistent with the truth-validation coast).
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.frames import FramesFactory

from puffsat_sim import mission, presets
from puffsat_sim.config import OrbitalConfig
from puffsat_sim.descent import (
    COAST_MAX_STEP_S,
    HANDOFF_ALT_M,
    coast_to_altitude,
    earth_model,
    to_absolute_date,
)
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import build_propagator
from puffsat_sim.tracker_fusion import CoflyerPhasing, format_coflyer_phasing, phasing_verdict

COFLYER_PERIGEE_RAISE_M: float = 100_000.0  # raise the rocket's perigee (survives the window)
COFLYER_APOGEE_LOWER_M: float = 100_000.0  # lower apogee the same → constant a → phase-locked
WINDOW_ALT_HI_M: float = HANDOFF_ALT_M  # 800 km — the terminal-window top (hand-off)
WINDOW_ALT_LO_M: float = mission.INTERCEPTION_ALT_M  # 200 km — the interception crossing
N_SAMPLES: int = 25


def coflyer_config(orbital_config: OrbitalConfig) -> OrbitalConfig:
    """The launch rocket's modified orbit: +perigee / −apogee (constant a), co-phased in M."""
    return replace(
        orbital_config,
        perigee_alt_m=orbital_config.perigee_alt_m + COFLYER_PERIGEE_RAISE_M,
        apogee_alt_m=orbital_config.apogee_alt_m - COFLYER_APOGEE_LOWER_M,
    )


def _altitude_m(earth: Any, position: Any, frame: Any, date: Any) -> float:
    """Ellipsoid altitude of a position sampled in ``frame`` at ``date`` (the window reference)."""
    return float(earth.transform(position, frame, date).getAltitude())


def _sample_positions(propagator: Any, frame: Any, final_date: Any, dates: list[Any]) -> list[Any]:
    """Sample a propagator's positions in ``frame`` at ``dates`` via its generated ephemeris."""
    generator = propagator.getEphemerisGenerator()
    propagator.propagate(final_date)
    ephemeris = generator.getGeneratedEphemeris()
    return [ephemeris.propagate(d).getPVCoordinates(frame).getPosition() for d in dates]


def run_coflyer_phasing(
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG, n_samples: int = N_SAMPLES
) -> CoflyerPhasing:
    """Fly the phase-locked rocket alongside the descent and reduce the window to a verdict."""
    frame = FramesFactory.getEME2000()
    earth = earth_model()
    epoch = to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)

    hi_date = coast_to_altitude(
        build_propagator(orbital_config, presets.j2(), COAST_MAX_STEP_S),
        epoch,
        period,
        earth,
        WINDOW_ALT_HI_M,
    ).getDate()
    lo_date = coast_to_altitude(
        build_propagator(orbital_config, presets.j2(), COAST_MAX_STEP_S),
        epoch,
        period,
        earth,
        WINDOW_ALT_LO_M,
    ).getDate()

    hi_s = float(hi_date.durationFrom(epoch))
    lo_s = float(lo_date.durationFrom(epoch))
    window_dates = [epoch.shiftedBy(float(t)) for t in np.linspace(hi_s, lo_s, n_samples)]

    puffsat_pos = _sample_positions(
        build_propagator(orbital_config, presets.j2(), COAST_MAX_STEP_S),
        frame,
        lo_date,
        window_dates,
    )
    rocket_pos = _sample_positions(
        build_propagator(coflyer_config(orbital_config), presets.j2(), COAST_MAX_STEP_S),
        frame,
        lo_date,
        window_dates,
    )

    ranges_m = [
        float(r.subtract(p).getNorm()) for p, r in zip(puffsat_pos, rocket_pos, strict=True)
    ]
    rocket_alts_m = [
        _altitude_m(earth, r, frame, d) for r, d in zip(rocket_pos, window_dates, strict=True)
    ]

    return phasing_verdict(ranges_m, rocket_alts_m, WINDOW_ALT_HI_M, WINDOW_ALT_LO_M)


def coflyer_phasing_report() -> None:
    print(format_coflyer_phasing(run_coflyer_phasing()))


if __name__ == "__main__":
    coflyer_phasing_report()
