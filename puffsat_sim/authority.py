"""Pure C3c terminal-authority / tail-correction core (ADR 0014 decision 5/6, no JVM).

C3b measured the terminal funnel: a thrust-limited catch radius of ~500 m at the 800 km
hand-off, set by the actuator (½·a_max·t² of isotropic authority over the descent time).
But the upstream 1σ error budget (ADR 0013: RSS 224 m) has a ~3σ tail (~672 m) that
exceeds that funnel, so a sliver of Rung-D runs would saturate it.  C3c quantifies the
two levers the architecture chose between (ADR 0014 decision 5):

* the **authority curve** — how the funnel grows if the aim burn starts higher
  (``½·a_max·t_descent²`` vs burn-start altitude), and
* the **MCC-2 cost curve** — Δv-per-km of an impulsive high-node trim that pre-shrinks
  the entry error, measured as the out-of-plane node-Δv → crossing-lateral lever vs node
  altitude (the kept corrector's Jacobian physics; A2's table located the authority at
  ~km scale: cheap high, dead low).

The funnel model is validated against C3b's measured 500 m at the 800 km point; the JVM
glue in :mod:`puffsat_sim.runs.authority` supplies the measured descent times and levers.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from puffsat_sim.coeff_requirement import BudgetEntry, rss_lateral
from puffsat_sim.dispersion import Vec3
from puffsat_sim.propellant import propellant_curve


def thrust_limited_radius_m(a_max_m_s2: float, t_descent_s: float) -> float:
    """The isotropic ½·a_max·t² lateral catch radius — the C3b-validated funnel model.

    The terminal burn has ``t_descent`` to null a lateral entry error at the actuator's
    ``a_max``; the largest it can erase is ½·a_max·t² (C3b bracketed this at 480 m vs the
    measured 500 m from the 800 km hand-off).
    """
    return 0.5 * a_max_m_s2 * t_descent_s**2


def saturation_dv_m_s(a_max_m_s2: float, t_descent_s: float) -> float:
    """The Δv ceiling ``a_max·t`` — the residual cliff C3b saw at ~3.94 m/s (800 km)."""
    return a_max_m_s2 * t_descent_s


def funnel_growth_dv_m_s(a_max_m_s2: float, entry_m: float) -> float:
    """Δv to null a lateral entry of ``entry_m`` at the grown-funnel edge: ``√(2·a_max·entry)``.

    Catching an entry the size of the upstream tail needs the funnel grown until ½·a_max·t²
    just reaches it (``t = √(2·entry/a_max)``), and nulling it then costs the full ceiling
    ``a_max·t`` — the saturation-edge price.  This is the propellant a *raise-the-burn*
    tail fix pays, the foil for the cheap high-node trim.
    """
    return math.sqrt(2.0 * a_max_m_s2 * entry_m)


def lateral_lever_m_per_m_s(
    crossing_plus_m: Vec3,
    crossing_minus_m: Vec3,
    dv_m_s: float,
    crossing_velocity_m_s: Vec3,
) -> float:
    """Crossing lateral (⊥ v_rel) shift per unit node Δv, from a ± central difference.

    A small out-of-plane Δv at the node moves the 200 km crossing; the component along
    the relative velocity is time-of-arrival (handled separately, the plate frame), so
    only the ⊥-v magnitude is the lateral correction lever.  This is the column of the
    kept corrector's Jacobian for an impulsive node burn, isolated to the lateral plane.
    """
    delta = np.asarray(crossing_plus_m, dtype=np.float64) - np.asarray(
        crossing_minus_m, dtype=np.float64
    )
    lever = delta / (2.0 * dv_m_s)
    v = np.asarray(crossing_velocity_m_s, dtype=np.float64)
    v_hat = v / np.linalg.norm(v)
    lateral = lever - (lever @ v_hat) * v_hat
    return float(np.linalg.norm(lateral))


def dv_per_km_m_s(lateral_lever_m_per_m_s: float) -> float:
    """Δv to move the crossing 1 km laterally from a node: ``1000 / lever`` (∞ if dead)."""
    if lateral_lever_m_per_m_s <= 0.0:
        return math.inf
    return 1000.0 / lateral_lever_m_per_m_s


@dataclass(frozen=True)
class AuthorityPoint:
    """One burn-start altitude on the authority curve: its descent time, funnel, ceiling."""

    handoff_alt_m: float
    t_descent_s: float
    radius_m: float
    dv_ceiling_m_s: float


def authority_point(handoff_alt_m: float, t_descent_s: float, a_max_m_s2: float) -> AuthorityPoint:
    """Build an authority-curve point from a measured descent time and the actuator a_max."""
    return AuthorityPoint(
        handoff_alt_m=handoff_alt_m,
        t_descent_s=t_descent_s,
        radius_m=thrust_limited_radius_m(a_max_m_s2, t_descent_s),
        dv_ceiling_m_s=saturation_dv_m_s(a_max_m_s2, t_descent_s),
    )


@dataclass(frozen=True)
class TrimPoint:
    """One node on the MCC-2 cost curve: its lateral lever and the Δv-per-km it implies."""

    node_alt_m: float
    lateral_lever_m_per_m_s: float
    dv_per_km_m_s: float


def trim_point(node_alt_m: float, lateral_lever_m_per_m_s: float) -> TrimPoint:
    """Build a trim-cost point from a measured lateral lever at a node altitude."""
    return TrimPoint(
        node_alt_m=node_alt_m,
        lateral_lever_m_per_m_s=lateral_lever_m_per_m_s,
        dv_per_km_m_s=dv_per_km_m_s(lateral_lever_m_per_m_s),
    )


@dataclass(frozen=True)
class TailAuthoritySpec:
    """The C3c measurement grid (ADR 0014 decision 5/6).

    ``handoff_altitudes_m`` are the burn-start altitudes for the authority curve, all at or
    above the 800 km hand-off (low → high) — coasting below it is the B0-unsafe stiff
    region, and 800 km is the chosen floor (ADR 0014 decision 1), so the curve measures
    the *raise-the-burn* direction.  ``node_altitudes_m`` are the impulsive-trim nodes for
    the cost curve (high → low, so the report reads authority dying as the node drops
    toward the §16.6 dead zone).  ``trim_dv_m_s`` is the out-of-plane central-difference
    impulse (small enough to stay in C0's linear regime, large enough to clear the ~cm
    integrator floor even at a near-dead low node).
    """

    handoff_altitudes_m: tuple[float, ...] = (800e3, 1000e3, 1500e3, 2000e3, 3000e3)
    node_altitudes_m: tuple[float, ...] = (100_000e3, 30_000e3, 10_000e3, 3000e3, 1000e3)
    trim_dv_m_s: float = 0.1
    tail_sigma: float = 3.0


@dataclass(frozen=True)
class TailAuthorityFinding:
    """The C3c measurement set: both curves plus the budget tail they are read against.

    ``measured_radius_m`` is C3b's measured funnel at ``measured_radius_alt_m`` (the
    authority-curve anchor); ``budget_entries`` is the full upstream lateral 1σ ledger
    (the C2a RSS 224 m: C1 nav + B1 erosion + the Cr coefficient prior) and ``tail_sigma``
    the multiple defining the entry tail the funnel must catch.
    """

    authority_points: tuple[AuthorityPoint, ...]
    trim_points: tuple[TrimPoint, ...]
    a_max_m_s2: float
    measured_radius_m: float
    measured_radius_alt_m: float
    budget_entries: tuple[BudgetEntry, ...]
    tail_sigma: float

    @property
    def budget_rss_1sigma_m(self) -> float:
        """The upstream lateral 1σ RSS (ADR 0013's 224 m)."""
        return rss_lateral(self.budget_entries)

    @property
    def tail_m(self) -> float:
        """The upstream error tail the funnel must cover: ``tail_sigma × budget RSS``."""
        return self.tail_sigma * self.budget_rss_1sigma_m

    @property
    def uncovered_tail_m(self) -> float:
        """How far the tail overruns the measured funnel (0 if the funnel already covers it)."""
        return max(0.0, self.tail_m - self.measured_radius_m)


def handoff_alt_to_cover_tail_m(finding: TailAuthorityFinding) -> float | None:
    """The burn-start altitude whose funnel just covers the tail, interpolated on the curve.

    ``None`` if no swept altitude's funnel reaches the tail (covering it by funnel growth
    alone is off the measured range — the authority-curve verdict that pushes the tail
    onto the impulsive trim instead).
    """
    tail = finding.tail_m
    points = sorted(finding.authority_points, key=lambda p: p.handoff_alt_m)
    for lo, hi in zip(points, points[1:], strict=False):
        if lo.radius_m <= tail <= hi.radius_m:
            frac = (tail - lo.radius_m) / (hi.radius_m - lo.radius_m)
            return lo.handoff_alt_m + frac * (hi.handoff_alt_m - lo.handoff_alt_m)
    return None


def cheapest_trim(finding: TailAuthorityFinding) -> TrimPoint:
    """The highest-authority trim node (smallest Δv-per-km) on the cost curve."""
    return min(finding.trim_points, key=lambda t: t.dv_per_km_m_s)


def trim_dv_to_cover_tail_m_s(finding: TailAuthorityFinding, trim: TrimPoint) -> float:
    """Δv for an impulsive trim at ``trim`` to remove the full tail (``dv/km × tail``)."""
    return trim.dv_per_km_m_s * finding.tail_m / 1000.0


def _km(alt_m: float) -> str:
    return f"{alt_m / 1e3:,.0f} km"


def format_tail_authority(finding: TailAuthorityFinding) -> str:
    """One-screen C3c report: the authority curve, the MCC-2 cost curve, the crossover verdict."""
    tail = finding.tail_m
    cheap = cheapest_trim(finding)
    cover_dv = trim_dv_to_cover_tail_m_s(finding, cheap)
    cover_points = propellant_curve(cover_dv)
    cover_curve = ", ".join(f"{p.fraction * 100:.4f}% @Isp{p.isp_s:.0f}" for p in cover_points)
    grow_alt = handoff_alt_to_cover_tail_m(finding)
    grow_text = _km(grow_alt) if grow_alt is not None else "above the largest swept burn-start"
    grow_dv = funnel_growth_dv_m_s(finding.a_max_m_s2, tail)
    ratio = grow_dv / cover_dv if cover_dv > 0.0 else math.inf

    budget = ", ".join(f"{e.label} {e.lateral_miss_m:.0f} m" for e in finding.budget_entries)
    lines = [
        "C3c terminal authority + MCC-2 tail correction (ADR 0014 decision 5/6)",
        f"  Funnel model ½·a_max·t² (a_max {finding.a_max_m_s2:.3f} m/s²); C3b anchor"
        f" {finding.measured_radius_m:.0f} m at {_km(finding.measured_radius_alt_m)}",
        f"  Upstream tail: {finding.tail_sigma:.0f}σ of the C2a budget RSS"
        f" {finding.budget_rss_1sigma_m:.0f} m = {tail:.0f} m"
        f" (800 km funnel leaves {finding.uncovered_tail_m:.0f} m uncovered)",
        f"    budget (lateral 1σ): {budget}",
        "  Authority curve — funnel vs burn-start altitude:",
    ]
    for p in sorted(finding.authority_points, key=lambda q: q.handoff_alt_m):
        covers = "covers tail" if p.radius_m >= tail else ""
        lines.append(
            f"    {_km(p.handoff_alt_m):>10}: t_descent {p.t_descent_s:6.1f} s →"
            f" funnel {p.radius_m:7.1f} m, ceiling {p.dv_ceiling_m_s:5.2f} m/s {covers}"
        )
    lines.append(
        f"    → funnel growth reaches the {tail:.0f} m tail at burn-start {grow_text},"
        f" but at the saturation edge (~{grow_dv:.2f} m/s, the ceiling)"
    )
    lines.append("  MCC-2 cost curve — impulsive out-of-plane trim vs node altitude:")
    for t in sorted(finding.trim_points, key=lambda q: q.node_alt_m, reverse=True):
        lines.append(
            f"    {_km(t.node_alt_m):>12}: lever {t.lateral_lever_m_per_m_s:10.1f} m/(m/s) →"
            f" {t.dv_per_km_m_s:9.4f} m/s per km"
        )
    lines.append(
        f"    → killing the {tail:.0f} m tail from {_km(cheap.node_alt_m)} costs"
        f" {cover_dv:.4f} m/s: {cover_curve}"
    )
    lines.append(
        f"  Verdict: the high-node impulsive trim covers the tail ~{ratio:.0f}× cheaper than"
        f" the funnel-growth edge ({cover_dv:.2f} vs {grow_dv:.2f} m/s) — MCC-2 vindicated"
        " (ADR 0014 decision 5); trim scheduling defers to Rung D"
    )
    return "\n".join(lines)
