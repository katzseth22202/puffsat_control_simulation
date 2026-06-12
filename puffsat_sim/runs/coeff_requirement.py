"""C2a coefficient requirement — the JVM run for :mod:`puffsat_sim.coeff_requirement` (ADR 0013)."""

from __future__ import annotations

from puffsat_sim import mission
from puffsat_sim.coeff_requirement import (
    format_coeff_requirement,
    summarize_coeff_requirement,
)
from puffsat_sim.config import OrbitalConfig
from puffsat_sim.control import report_controller
from puffsat_sim.navigation import NavSweepSpec, summarize_nav_requirement
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.runs.navigation import run_nav_sweep
from puffsat_sim.runs.sweep import run_sweep
from puffsat_sim.sweep import SweepSpec


def coeff_requirement_report(
    catch_radius_m: float = 5_000.0,
    prior_sigma_factor: float = 0.2,
    cut_points: int = 3,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Run the C2a coefficient-knowledge requirement (ADR 0013): measured cuts → tolerance vs prior.

    Thin glue over already-covered harnesses (the B2 no-integration-test precedent):
    Φ from a minimal C0 sweep, ``∂Δv/∂c`` from two 1D A3 cuts under the same LM-damped
    corrector, the verdict from the pure :mod:`puffsat_sim.coeff_requirement` chain.
    The coast for the analytic SRP cross-check is the apogee→perigee half period.
    """
    nav_result = run_nav_sweep(NavSweepSpec(points_per_sign=1), report_controller, orbital_config)
    phi = summarize_nav_requirement(nav_result).phi
    cd_cut = run_sweep(
        SweepSpec(cd_points=cut_points, cr_points=1), report_controller, orbital_config
    )
    cr_cut = run_sweep(
        SweepSpec(cd_points=1, cr_points=cut_points), report_controller, orbital_config
    )
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    requirement = summarize_coeff_requirement(
        phi,
        cd_cut=cd_cut,
        cr_cut=cr_cut,
        catch_radius_m=catch_radius_m,
        prior_sigma_factor=prior_sigma_factor,
        coast_duration_s=keplerian_period(semi_major) / 2.0,
    )
    return format_coeff_requirement(requirement)
