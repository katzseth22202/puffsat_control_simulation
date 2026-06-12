"""A3 controllability sweep — the JVM run for :mod:`puffsat_sim.sweep` (ADR 0007)."""

from __future__ import annotations

from puffsat_sim import mission
from puffsat_sim.config import OrbitalConfig
from puffsat_sim.control import Controller
from puffsat_sim.dispersion import RunInputs
from puffsat_sim.montecarlo import RunVariant, build_context, run_record
from puffsat_sim.sweep import SweepResult, SweepSpec, grid_inputs


def run_sweep(
    spec: SweepSpec,
    control: Controller,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
    toa_window_s: float | None = None,
) -> SweepResult:
    """Run the deterministic A3 controllability grid against a fixed targeter (ADR 0007).

    Same physics path as :func:`puffsat_sim.montecarlo.run_ensemble` (shared
    :func:`~puffsat_sim.montecarlo.build_context` / :func:`~puffsat_sim.montecarlo.run_record`
    / nominal crossing), but the inputs are the deterministic
    :func:`~puffsat_sim.sweep.grid_inputs` (zero injection, swept Cd/Cr) rather than
    stochastic draws, and ``control`` is required — A3 maps the *required Δv*, so there is no
    open-loop variant.  ``toa_window_s`` arms the spurious-far-root gate (decision 3iii); the
    caller sizes it off the capstone's open-loop ToA dispersion.

    ``SweepResult.nominal`` is a dedicated factor-(1,1) reference run (zero coefficient error,
    zero injection) so the perigee/ToA overlays have a baseline even when the grid does not
    land a point exactly on nominal.
    """
    ctx = build_context(orbital_config)
    variant = RunVariant(control=control, toa_window_s=toa_window_s)
    records = tuple(run_record(ctx, inputs, variant) for inputs in grid_inputs(spec))
    nominal_run_inputs = RunInputs(
        run_index=-1,
        dv_rtn_m_s=(0.0, 0.0, 0.0),
        cd_area_over_mass=spec.cd_area_over_mass,
        cr_area_over_mass=spec.cr_area_over_mass,
        f10p7=spec.f10p7,
        ap=spec.ap,
    )
    nominal_record = run_record(ctx, nominal_run_inputs, variant)
    return SweepResult(spec=spec, records=records, nominal=nominal_record)
