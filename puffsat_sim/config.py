"""Mission configuration dataclasses — pure Python, no JVM dependency."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

from puffsat_sim.forces import Perturbation


@dataclass(frozen=True)
class OrbitalConfig:
    """Keplerian orbital elements for one simulation run.

    perigee_alt_m and apogee_alt_m are fixed mission parameters not varied in
    the Monte Carlo.  inclination_rad, raan_rad, arg_of_perigee_rad,
    mean_anomaly_at_epoch_rad, and epoch represent injection state and are
    sampled across Monte Carlo runs.

    epoch must be a timezone-aware UTC datetime.
    mean_anomaly_at_epoch_rad=π places the satellite at apogee (deployment point).
    """

    perigee_alt_m: float
    apogee_alt_m: float
    inclination_rad: float
    raan_rad: float
    arg_of_perigee_rad: float
    mean_anomaly_at_epoch_rad: float
    epoch: datetime  # UTC, must be timezone-aware


@dataclass(frozen=True)
class PhysicsConfig:
    """Force model selection for one simulation run — a set of active Perturbations.

    An empty tuple selects pure Keplerian (point-mass / two-body) propagation.
    Each Perturbation carries only the parameters its force needs; per-run Monte
    Carlo draws (e.g. log-normal Cd·(A/m) / Cr·(A/m)) are applied by constructing
    the relevant spec with the sampled value and held fixed for the run, because a
    constant coefficient error integrates coherently over the coast and drives
    perigee dispersion.

    Use :mod:`puffsat_sim.presets` for the named, content-described bundles.
    """

    perturbations: tuple[Perturbation, ...] = ()

    @property
    def is_keplerian(self) -> bool:
        """True when no perturbations are active (pure two-body propagation)."""
        return not self.perturbations
