"""C3c terminal authority + MCC-2 tail correction — the JVM run for :mod:`puffsat_sim.authority`.

Two measurements over the proven §6.2 descent (ADR 0014 decision 5/6):

* the **authority curve** — for each burn-start altitude, the descent time left to the
  200 km crossing (``full descent − coast-to-altitude``, both on the same nominal path),
  fed through the C3b-validated ½·a_max·t² funnel model; and
* the **MCC-2 cost curve** — at each high node, an out-of-plane ±Δv central difference
  descended to the crossing gives the lateral lever (the kept corrector's Jacobian column
  for an impulsive node burn), and so the Δv-per-km of the tail-correcting trim.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import puffsat_sim.jvm  # noqa: F401  boots the JVM before any org.orekit import

from org.orekit.frames import FramesFactory
from org.orekit.orbits import CartesianOrbit
from org.orekit.utils import Constants, TimeStampedPVCoordinates

from org.hipparchus.geometry.euclidean.threed import Vector3D

from puffsat_sim import mission, presets
from puffsat_sim.anti_drag import PEAK_THRUST_LIMIT_N
from puffsat_sim.authority import (
    AuthorityPoint,
    TailAuthorityFinding,
    TailAuthoritySpec,
    TrimPoint,
    authority_point,
    format_tail_authority,
    lateral_lever_m_per_m_s,
    trim_point,
)
from puffsat_sim.coeff_requirement import MEASURED_BUDGET, BudgetEntry
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.descent import (
    COAST_MAX_STEP_S,
    HANDOFF_ALT_M,
    TERMINAL_MAX_STEP_S,
    Crossing,
    coast_to_altitude,
    descend,
    earth_model,
    propagate_to_interception,
    to_absolute_date,
    vec3,
)
from puffsat_sim.dispersion import Vec3, rtn_basis, rtn_to_cartesian
from puffsat_sim.orbital_math import keplerian_elements, keplerian_period
from puffsat_sim.propagator import build_propagator, build_propagator_from_orbit
from puffsat_sim.runs.anti_drag import PUFFSAT_WET_MASS_KG

# The C3b measured funnel at the 800 km hand-off (ADR 0014 C3b findings): the
# authority-curve anchor the ½·a_max·t² model is read against.
C3B_MEASURED_RADIUS_M: float = 500.0

# The C2a Cr coefficient-prior lateral miss banked from ADR 0013 (0.2 factor prior ×
# 745 m/factor measured sensitivity); MEASURED_BUDGET carries only the nav + erosion
# contributions, so the full upstream tail the funnel must catch adds this third leg —
# RSS(141, 89, 149) = 224 m, the budget ADR 0014 decision 5 reads the 3σ tail against.
C2A_CR_PRIOR_LATERAL_M: float = 149.0
UPSTREAM_BUDGET: tuple[BudgetEntry, ...] = (
    *MEASURED_BUDGET,
    BudgetEntry("coefficient prior (Cr, C2a)", C2A_CR_PRIOR_LATERAL_M),
)


@dataclass(frozen=True)
class AuthorityContext:
    """Per-sweep constants: the proven nominal descent and the actuator/frame handles."""

    physics: PhysicsConfig
    earth: Any
    epoch: Any
    period: float
    frame: Any
    mu: float
    a_max_m_s2: float
    apogee_orbit: Any
    full_descent_toa_s: float
    crossing_velocity_m_s: Vec3


def build_authority_context(
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> AuthorityContext:
    """Fix the nominal descent (apogee → 200 km) once; its time/velocity anchor both curves."""
    physics = presets.full_force()
    earth = earth_model()
    epoch = to_absolute_date(orbital_config.epoch)
    semi_major, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(semi_major)
    apogee_orbit = (
        build_propagator(orbital_config, physics, COAST_MAX_STEP_S).getInitialState().getOrbit()
    )
    nominal = descend(apogee_orbit, physics, epoch, period, earth)
    return AuthorityContext(
        physics=physics,
        earth=earth,
        epoch=epoch,
        period=period,
        frame=FramesFactory.getEME2000(),
        mu=float(Constants.WGS84_EARTH_MU),
        a_max_m_s2=PEAK_THRUST_LIMIT_N / PUFFSAT_WET_MASS_KG,
        apogee_orbit=apogee_orbit,
        full_descent_toa_s=nominal.toa_s,
        crossing_velocity_m_s=nominal.velocity_m_s,
    )


def _coast_from_apogee(ctx: AuthorityContext, altitude_m: float) -> Any:
    """Coast the proven big-step arc from apogee to a descending altitude event."""
    coast = build_propagator_from_orbit(ctx.apogee_orbit, ctx.physics, COAST_MAX_STEP_S)
    return coast_to_altitude(coast, ctx.epoch, ctx.period, ctx.earth, altitude_m)


def _descend_from_node(ctx: AuthorityContext, orbit: Any, node_alt_m: float) -> Crossing:
    """Descend a (possibly perturbed) node orbit to the 200 km crossing via the regime switch.

    Above the 800 km hand-off the proven coast→800 km→terminal kernel applies; at/below it
    the state is already in the stiff region, so descend on the terminal cap straight down.
    """
    if node_alt_m > HANDOFF_ALT_M:
        return descend(orbit, ctx.physics, ctx.epoch, ctx.period, ctx.earth)
    return propagate_to_interception(
        build_propagator_from_orbit(orbit, ctx.physics, TERMINAL_MAX_STEP_S),
        ctx.epoch,
        ctx.period,
        ctx.earth,
    )


def measure_authority_point(ctx: AuthorityContext, handoff_alt_m: float) -> AuthorityPoint:
    """The funnel a terminal burn started at ``handoff_alt_m`` buys (½·a_max·t_descent²)."""
    handoff = _coast_from_apogee(ctx, handoff_alt_m)
    t_handoff = float(handoff.getDate().durationFrom(ctx.epoch))
    return authority_point(handoff_alt_m, ctx.full_descent_toa_s - t_handoff, ctx.a_max_m_s2)


def measure_trim_point(ctx: AuthorityContext, node_alt_m: float, dv_m_s: float) -> TrimPoint:
    """The lateral lever of an impulsive out-of-plane trim at ``node_alt_m`` → its Δv-per-km."""
    node = _coast_from_apogee(ctx, node_alt_m)
    pv = node.getPVCoordinates()
    basis = rtn_basis(vec3(pv.getPosition()), vec3(pv.getVelocity()))
    crossings: list[Crossing] = []
    for sign in (1.0, -1.0):
        dv_vec = rtn_to_cartesian((0.0, 0.0, sign * dv_m_s), basis)  # out-of-plane (RTN normal)
        velocity = pv.getVelocity().add(Vector3D(dv_vec[0], dv_vec[1], dv_vec[2]))
        orbit = CartesianOrbit(
            TimeStampedPVCoordinates(node.getDate(), pv.getPosition(), velocity),
            ctx.frame,
            ctx.mu,
        )
        crossings.append(_descend_from_node(ctx, orbit, node_alt_m))
    lever = lateral_lever_m_per_m_s(
        crossings[0].position_m, crossings[1].position_m, dv_m_s, ctx.crossing_velocity_m_s
    )
    return trim_point(node_alt_m, lever)


def run_tail_authority(
    spec: TailAuthoritySpec | None = None,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> TailAuthorityFinding:
    """Run the C3c authority + MCC-2 cost sweep over one shared nominal descent (ADR 0014)."""
    spec = spec if spec is not None else TailAuthoritySpec()
    ctx = build_authority_context(orbital_config)
    return TailAuthorityFinding(
        authority_points=tuple(measure_authority_point(ctx, h) for h in spec.handoff_altitudes_m),
        trim_points=tuple(
            measure_trim_point(ctx, n, spec.trim_dv_m_s) for n in spec.node_altitudes_m
        ),
        a_max_m_s2=ctx.a_max_m_s2,
        measured_radius_m=C3B_MEASURED_RADIUS_M,
        measured_radius_alt_m=HANDOFF_ALT_M,
        budget_entries=UPSTREAM_BUDGET,
        tail_sigma=spec.tail_sigma,
    )


def tail_authority_report(
    spec: TailAuthoritySpec | None = None,
    orbital_config: OrbitalConfig = mission.NOMINAL_CONFIG,
) -> str:
    """Run the C3c sweep and format the one-screen report."""
    return format_tail_authority(run_tail_authority(spec, orbital_config))
