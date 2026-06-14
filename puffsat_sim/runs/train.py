"""Rung D / D1.1 closed-loop train ensemble — the JVM run for :mod:`puffsat_sim.train` (ADR 0018).

Strings the C-rung pieces into one closed-loop run per PuffSat and reduces a **train** of them
into the headline **P(capture) about the centroid**.  Per the user-confirmed D1.1 architecture
(Φ-composed entry + flown terminal): each unit's hand-off lateral entry offset is *sampled* from
the characterized C0/C1/C2a budget (the midcourse residual is provably linear in nav/coefficient
error — C0's Φ), and the C3b ZEM terminal loop is *flown* (:func:`run_guidance`) through the
catch-radius cliff, the significance-gate rectification, and the σ_θ·R tracker noise — the binding
nonlinearity.  The corrector-in-every-run brute-force validation, MCC-2 scheduling, and the
node-count Σ sweep are later D1.x sub-slices (ADR 0018 decisions 4/6).
"""

from __future__ import annotations

import numpy as np

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from puffsat_sim import mission
from puffsat_sim.config import OrbitalConfig
from puffsat_sim.guidance import TrackerGrade
from puffsat_sim.montecarlo import physics_from_inputs
from puffsat_sim.runs.guidance import GuidanceContext, build_guidance_context, run_guidance
from puffsat_sim.train import (
    TrainDispersionSpec,
    TrainEnsembleFinding,
    format_train_ensemble,
    sample_train,
    sample_train_entry_offsets,
    summarize_train_ensemble,
)

D1_MASTER_SEED: int = 20260613
SMOKE_N_UNITS: int = 16
SMOKE_N_TRAINS: int = 2


def run_train_dispersion(
    spec: TrainDispersionSpec,
    train_index: int,
    master_seed: int = D1_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    ctx: GuidanceContext | None = None,
) -> TrainEnsembleFinding:
    """Fly one train's units through the closed terminal loop and reduce them (D1.1).

    Each unit carries its train-sampled coefficients/space-weather into the *truth* drag and its
    Φ-composed lateral entry offset into the funnel; the loop injects the σ_θ tracker noise.  Pass
    a shared ``ctx`` to fly several trains over one hand-off/aim-point build.
    """
    ctx = ctx if ctx is not None else build_guidance_context(orbital_config)
    grade = TrackerGrade(
        sigma_theta_rad=spec.tracker_sigma_theta_rad, sigma_range_m=spec.tracker_sigma_range_m
    )
    inputs = sample_train(master_seed, spec, train_index)
    entries = sample_train_entry_offsets(master_seed, spec, train_index)

    runs = [
        run_guidance(
            ctx,
            entry_offset_lateral_m=entries[j],
            grade=grade,
            control_period_s=spec.control_period_s,
            rng=np.random.default_rng((master_seed, train_index, j)),
            truth_physics=physics_from_inputs(inputs[j]),
        )
        for j in range(spec.n_units)
    ]
    return summarize_train_ensemble(
        [r.miss for r in runs],
        [r.plan for r in runs],
        [r.perigee_alt_m for r in runs],
        spec,
    )


def train_dispersion_report(
    n_trains: int = SMOKE_N_TRAINS,
    spec: TrainDispersionSpec | None = None,
    master_seed: int = D1_MASTER_SEED,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Fly ``n_trains`` independent trains over one shared context and format each D1.1 finding."""
    spec = spec if spec is not None else TrainDispersionSpec(n_units=SMOKE_N_UNITS)
    ctx = build_guidance_context(orbital_config)
    return "\n\n".join(
        f"Train {t} (seed {master_seed}):\n"
        + format_train_ensemble(run_train_dispersion(spec, t, master_seed, ctx=ctx))
        for t in range(n_trains)
    )


if __name__ == "__main__":
    print(train_dispersion_report())
