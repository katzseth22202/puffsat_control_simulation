"""Mission configuration dataclasses — pure Python, no JVM dependency."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime


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
    """Force model selection and per-run parameters for one simulation run.

    geopotential_degree=0 with all other fields at defaults selects pure
    Keplerian (point-mass) propagation.

    srp_cr_area_over_mass and drag_cd_area_over_mass are None when the
    corresponding force model is inactive; a float (m²/kg) enables it.
    In the Monte Carlo both are drawn per-trajectory from a multiplicative
    log-normal (~10–30% 1-sigma) and held fixed for the run, because a
    constant coefficient error integrates coherently over the coast and
    drives perigee dispersion.
    """

    geopotential_degree: int = 0
    third_body: bool = False
    srp_cr_area_over_mass: float | None = None   # Cr·(A/m) [m²/kg]; None = SRP off
    drag_cd_area_over_mass: float | None = None  # Cd·(A/m) [m²/kg]; None = drag off
    f10p7: float = 150.0  # solar flux index for atmospheric drag model
    ap: float = 15.0      # geomagnetic index for atmospheric drag model

    @property
    def is_keplerian(self) -> bool:
        """True when no perturbations are active (pure two-body propagation)."""
        return (
            self.geopotential_degree == 0
            and not self.third_body
            and self.srp_cr_area_over_mass is None
            and self.drag_cd_area_over_mass is None
        )

    @classmethod
    def rung_keplerian(cls) -> PhysicsConfig:
        """Point-mass / Keplerian — no perturbations (Rung 1 baseline)."""
        return cls()

    @classmethod
    def rung_2a(cls) -> PhysicsConfig:
        """J2 geopotential only (Rung 2a)."""
        return cls(geopotential_degree=2)

    @classmethod
    def rung_2b(cls) -> PhysicsConfig:
        """J2 + third-body Sun and Moon (Rung 2b)."""
        return cls(geopotential_degree=2, third_body=True)

    @classmethod
    def rung_2c(cls, cr_area_over_mass: float = 0.02) -> PhysicsConfig:
        """J2 + third-body + SRP cannonball (Rung 2c).

        cr_area_over_mass: Cr·(A/m) [m²/kg].
        """
        return cls(
            geopotential_degree=2,
            third_body=True,
            srp_cr_area_over_mass=cr_area_over_mass,
        )

    @classmethod
    def rung_2d(
        cls,
        cr_area_over_mass: float = 0.02,
        cd_area_over_mass: float = 0.04,
    ) -> PhysicsConfig:
        """Full force model: J2 + third-body + SRP + drag (Rung 2d).

        cr_area_over_mass: Cr·(A/m) [m²/kg].
        cd_area_over_mass: Cd·(A/m) [m²/kg].
        """
        return cls(
            geopotential_degree=2,
            third_body=True,
            srp_cr_area_over_mass=cr_area_over_mass,
            drag_cd_area_over_mass=cd_area_over_mass,
        )
