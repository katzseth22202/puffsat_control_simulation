"""C0 navigation-requirement sweep — the JVM run for :mod:`puffsat_sim.navigation` (ADR 0011)."""

from __future__ import annotations

from collections.abc import Sequence

from puffsat_sim import mission
from puffsat_sim.config import OrbitalConfig
from puffsat_sim.control import Controller, report_controller
from puffsat_sim.montecarlo import RunVariant, build_context, nominal_inputs, run_record
from puffsat_sim.navigation import (
    NavSweepResult,
    NavSweepSpec,
    format_nav_requirement,
    nav_grid_offsets,
    summarize_nav_requirement,
)


def run_nav_sweep(
    spec: NavSweepSpec,
    control: Controller,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> NavSweepResult:
    """Run the deterministic C0 navigation-error sensitivity sweep (ADR 0011).

    Same physics path as :func:`puffsat_sim.montecarlo.run_ensemble` (shared
    :func:`~puffsat_sim.montecarlo.build_context` / :func:`~puffsat_sim.montecarlo.run_record`
    / nominal crossing), but the perturbation is a **predict-side** apogee-RTN nav-error
    offset, with zero injection and nominal coefficients (perfect model) — so the only
    predict-vs-execute divergence is the nav error.  Each cell runs the corrector from the
    perturbed estimate and executes against truth; the recorded ``miss_rtn_m`` is the
    residual ``−Φδ`` (uncontrollable at the apogee node) and ``total_dv_m_s`` the phantom
    correction the corrector burned chasing it.  The zero cell is the on-target reference
    (residual ~0).  ``control`` is required — the corrector is C0's subject, not an option.
    """
    ctx = build_context(orbital_config)
    cells = nav_grid_offsets(spec)
    records = tuple(
        run_record(
            ctx,
            nominal_inputs(cell.cell_index),
            RunVariant(control=control, nav_offset_rtn6=cell.offset_rtn6),
        )
        for cell in cells
    )
    return NavSweepResult(spec=spec, cells=cells, records=records)


def nav_requirement_report(
    spec: NavSweepSpec | None = None,
    catch_radii_m: Sequence[float] = (5_000.0, 1_000.0, 100.0),
) -> str:
    """Run the C0 nav-error sweep and reduce it to the navigation-requirement report (ADR 0011)."""
    result = run_nav_sweep(spec or NavSweepSpec(), control=report_controller)
    req = summarize_nav_requirement(result)
    phi_lines = ["  Φ (3×6 apogee→crossing sensitivity), rows R/T/N, cols R/T/N-pos R/T/N-vel:"]
    phi_lines += ["    " + "  ".join(f"{v:+.3e}" for v in row) for row in req.phi]
    return format_nav_requirement(req, catch_radii_m) + "\n" + "\n".join(phi_lines)
