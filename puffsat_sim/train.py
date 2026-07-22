"""Pure train-mode dispersion core — the first Rung-D / D1 sub-slice (ADR 0016/0018, no JVM).

A **train** is the batch of PuffSats deployed together (CONTEXT: *Train*). The single-PuffSat
sampler (:mod:`puffsat_sim.dispersion`) draws every axis independently per run — valid for
per-unit dispersion, silent on the train-relative split a fleet claim needs. ADR 0016 splits the
draws:

* **Shared (one per train, common-mode):** coefficient model *bias* (the systematic offset of the
  cannonball prior — ADR 0009), the F10.7/Ap *drivers* (one atmosphere for the whole train), and
  the deployer *systematic* injection. The §16.7 multiplicative-density gap collapses here: the
  *common* density component is the shared space-weather driver.
* **Per-unit (one per PuffSat):** coefficient unit *spread*, injection *scatter*. The per-unit
  density error ≈ the per-unit ``Cd·(A/m)`` spread to first order (drag ∝ ρ·Cd·A/m), so it needs
  no separate axis.

The verdict splits accordingly (ADR 0016, CONTEXT: *Centroid retarget*): the train **centroid
drift** is absorbed by the plane's ±2 km centroid retarget, while the **scatter about the
centroid** is what each PuffSat's terminal burn must fit inside the plate (ADR 0015). This module
owns the pure sampling (composing the same :class:`~puffsat_sim.dispersion.RunInputs` the JVM
``run_record`` already consumes, so D1.1 wires in with no record change) and the train-relative
reduction (reusing the :mod:`puffsat_sim.guidance` plate-capture machinery).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.dispersion import RunInputs, Vec3, lognormal_factor
from puffsat_sim.guidance import (
    CAPTURE_SIGMA_MAX_M,
    PLATE_RADIUS_M,
    TOA_LIMIT_S,
    PlateMiss,
    capture_fraction,
)
from puffsat_sim.propellant import PropellantPoint, propellant_curve
from puffsat_sim.terminal import FeedforwardPlan

Vec2 = tuple[float, float]

# The plane's launch-window retarget capability (CONTEXT: Centroid retarget; ADR 0006/0016).
CENTROID_RETARGET_M: float = 2000.0

# Masks the D1.1 hand-off entry-offset seed tree off the coefficient/injection tree, so a
# train's shared entry offset is an independent draw from its shared coefficient bias.
_ENTRY_SEED_MASK: int = 0x5EC0FFEE

# Characterized C0/C1/C2a budget legs the hand-off entry offset is sampled from (swept in D1.x),
# as 2-D lateral RMS *magnitudes*: the per-unit is C1's validated nav lateral (141 m, in-plane /
# nav-dominated); the shared is C2a's coefficient-bias lateral (0.2 Cr prior × 745 m/factor ≈
# 149 m), the common-mode centroid.  Sampled as an isotropic 2-D Gaussian, so each axis draws at
# σ/√2 (E[|d|²] = 2·(σ/√2)² = σ²) to match the characterized magnitude.
ENTRY_LATERAL_PERUNIT_M: float = 141.0
ENTRY_LATERAL_SHARED_M: float = 149.0
_PER_AXIS_FROM_MAGNITUDE: float = 1.0 / math.sqrt(2.0)
# The B2 representative mission midcourse Δv (correction 2.17 + anti-drag 0.015 m/s); the terminal
# aim Δv (per-unit, flown) stacks on top for the <2 % propellant check.
MIDCOURSE_DV_M_S: float = 2.19


@dataclass(frozen=True)
class TrainDispersionSpec:
    """Train-mode dispersion: each :class:`DispersionSpec` σ split into shared + per-unit.

    Coefficient and space-weather σ are fractional (coefficient of variation) for the
    median-nominal multiplicative log-normals; injection σ are per-RTN-axis [m/s].  The
    shared coefficient bias and per-unit spread *compose multiplicatively*, so the marginal
    per-unit log-variance is the sum (``s_bias² + s_spread²``).  ``centroid_retarget_m`` and
    ``n_units`` are reduction pins (the plane capability and the train size), not draws.
    """

    n_units: int = 1
    cd_area_over_mass: float = 0.04
    cr_area_over_mass: float = 0.02
    f10p7: float = 150.0
    ap: float = 15.0
    # Shared (per-train, common-mode).
    sigma_cd_bias_frac: float = 0.0
    sigma_cr_bias_frac: float = 0.0
    sigma_f10p7_frac: float = 0.0
    sigma_ap_frac: float = 0.0
    sigma_dv_systematic_radial_m_s: float = 0.0
    sigma_dv_systematic_transverse_m_s: float = 0.0
    sigma_dv_systematic_normal_m_s: float = 0.0
    # Per-unit (per-PuffSat).
    sigma_cd_spread_frac: float = 0.0
    sigma_cr_spread_frac: float = 0.0
    sigma_dv_scatter_radial_m_s: float = 0.0
    sigma_dv_scatter_transverse_m_s: float = 0.0
    sigma_dv_scatter_normal_m_s: float = 0.0
    # D1.1 hand-off entry offset (the Φ-composed midcourse residual the terminal funnel flies
    # from) and the flown terminal nav grade; defaults are the characterized C0/C1/C2a budget legs.
    sigma_entry_lateral_shared_m: float = ENTRY_LATERAL_SHARED_M
    sigma_entry_lateral_perunit_m: float = ENTRY_LATERAL_PERUNIT_M
    tracker_sigma_theta_rad: float = 10e-6
    tracker_sigma_range_m: float = 1.0
    control_period_s: float = 1.0
    midcourse_dv_m_s: float = MIDCOURSE_DV_M_S
    # Reduction pin.
    centroid_retarget_m: float = CENTROID_RETARGET_M


@dataclass(frozen=True)
class _SharedDraw:
    """One train's common-mode draw (drawn once, shared by every unit)."""

    bias_cd: float
    bias_cr: float
    f10p7: float
    ap: float
    sys_dv_rtn: Vec3


def _shared_rng(master_seed: int, train_index: int) -> np.random.Generator:
    # spawn_key=(train_index,) — arity 1 keeps the shared stream distinct from the
    # per-unit (train_index, unit_index) streams, so a unit replays standalone (§14.2).
    return np.random.default_rng(
        np.random.SeedSequence(entropy=master_seed, spawn_key=(train_index,))
    )


def _unit_rng(master_seed: int, train_index: int, unit_index: int) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence(entropy=master_seed, spawn_key=(train_index, unit_index))
    )


def _sample_shared(master_seed: int, spec: TrainDispersionSpec, train_index: int) -> _SharedDraw:
    """Draw one train's common-mode part in a fixed order (so a seed reproduces it)."""
    rng = _shared_rng(master_seed, train_index)
    bias_cd = lognormal_factor(rng, spec.sigma_cd_bias_frac)
    bias_cr = lognormal_factor(rng, spec.sigma_cr_bias_frac)
    f10p7 = spec.f10p7 * lognormal_factor(rng, spec.sigma_f10p7_frac)
    ap = spec.ap * lognormal_factor(rng, spec.sigma_ap_frac)
    sys_dv: Vec3 = (
        float(rng.normal(0.0, spec.sigma_dv_systematic_radial_m_s)),
        float(rng.normal(0.0, spec.sigma_dv_systematic_transverse_m_s)),
        float(rng.normal(0.0, spec.sigma_dv_systematic_normal_m_s)),
    )
    return _SharedDraw(bias_cd, bias_cr, f10p7, ap, sys_dv)


def _compose_unit(
    master_seed: int,
    spec: TrainDispersionSpec,
    train_index: int,
    unit_index: int,
    shared: _SharedDraw,
) -> RunInputs:
    """Combine one train's shared draw with unit ``unit_index``'s per-unit draw."""
    rng = _unit_rng(master_seed, train_index, unit_index)
    spread_cd = lognormal_factor(rng, spec.sigma_cd_spread_frac)
    spread_cr = lognormal_factor(rng, spec.sigma_cr_spread_frac)
    scatter: Vec3 = (
        float(rng.normal(0.0, spec.sigma_dv_scatter_radial_m_s)),
        float(rng.normal(0.0, spec.sigma_dv_scatter_transverse_m_s)),
        float(rng.normal(0.0, spec.sigma_dv_scatter_normal_m_s)),
    )
    dv_rtn: Vec3 = (
        shared.sys_dv_rtn[0] + scatter[0],
        shared.sys_dv_rtn[1] + scatter[1],
        shared.sys_dv_rtn[2] + scatter[2],
    )
    return RunInputs(
        run_index=train_index * spec.n_units + unit_index,
        dv_rtn_m_s=dv_rtn,
        cd_area_over_mass=spec.cd_area_over_mass * shared.bias_cd * spread_cd,
        cr_area_over_mass=spec.cr_area_over_mass * shared.bias_cr * spread_cr,
        f10p7=shared.f10p7,
        ap=shared.ap,
    )


def sample_train(
    master_seed: int, spec: TrainDispersionSpec, train_index: int
) -> tuple[RunInputs, ...]:
    """Draw one train: one shared common-mode draw, then ``n_units`` per-unit draws."""
    shared = _sample_shared(master_seed, spec, train_index)
    return tuple(
        _compose_unit(master_seed, spec, train_index, j, shared) for j in range(spec.n_units)
    )


def replay_train_unit(
    master_seed: int, spec: TrainDispersionSpec, train_index: int, unit_index: int
) -> RunInputs:
    """Reconstruct a single train-unit's draws standalone (§14.2), without the train.

    Re-draws the shared part from ``train_index`` and the per-unit part from
    ``(train_index, unit_index)`` — the train analog of ``dispersion.replay_inputs``.
    """
    shared = _sample_shared(master_seed, spec, train_index)
    return _compose_unit(master_seed, spec, train_index, unit_index, shared)


def _entry_shared_rng(master_seed: int, train_index: int) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence(entropy=master_seed ^ _ENTRY_SEED_MASK, spawn_key=(train_index,))
    )


def _entry_unit_rng(master_seed: int, train_index: int, unit_index: int) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence(
            entropy=master_seed ^ _ENTRY_SEED_MASK, spawn_key=(train_index, unit_index)
        )
    )


def _entry_offset(
    spec: TrainDispersionSpec, shared_lateral: Vec2, unit_rng: np.random.Generator
) -> Vec2:
    """One unit's 2-D lateral (⊥v) hand-off entry = shared centroid + per-unit scatter."""
    sigma_axis = spec.sigma_entry_lateral_perunit_m * _PER_AXIS_FROM_MAGNITUDE
    return (
        shared_lateral[0] + float(unit_rng.normal(0.0, sigma_axis)),
        shared_lateral[1] + float(unit_rng.normal(0.0, sigma_axis)),
    )


def _shared_entry_lateral(master_seed: int, spec: TrainDispersionSpec, train_index: int) -> Vec2:
    rng = _entry_shared_rng(master_seed, train_index)
    sigma_axis = spec.sigma_entry_lateral_shared_m * _PER_AXIS_FROM_MAGNITUDE
    return (float(rng.normal(0.0, sigma_axis)), float(rng.normal(0.0, sigma_axis)))


def sample_train_entry_offsets(
    master_seed: int, spec: TrainDispersionSpec, train_index: int
) -> tuple[Vec2, ...]:
    """The D1.1 hand-off lateral entry offsets for a train: shared centroid + per-unit scatter.

    The Φ-composed midcourse residual the terminal funnel flies from (ADR 0018): the *shared*
    lateral (centroid drift the retarget absorbs) is drawn once; each unit adds an independent
    *per-unit* lateral (the scatter the funnel must null).  A 2-D vector in the ⊥v plane, mapped
    to the orbit-normal / in-plane-⊥v axes by the JVM run.
    """
    shared_lateral = _shared_entry_lateral(master_seed, spec, train_index)
    return tuple(
        _entry_offset(spec, shared_lateral, _entry_unit_rng(master_seed, train_index, j))
        for j in range(spec.n_units)
    )


def replay_train_entry_offset(
    master_seed: int, spec: TrainDispersionSpec, train_index: int, unit_index: int
) -> Vec2:
    """Reconstruct a single unit's hand-off entry offset standalone (§14.2)."""
    shared_lateral = _shared_entry_lateral(master_seed, spec, train_index)
    return _entry_offset(
        spec, shared_lateral, _entry_unit_rng(master_seed, train_index, unit_index)
    )


@dataclass(frozen=True)
class TrainCaptureStats:
    """A train's arrival reduction split into centroid drift (retarget) vs scatter (plate)."""

    n_units: int
    centroid_m: tuple[float, float]
    centroid_retarget_m: float
    scatter_cov_m2: tuple[tuple[float, float], tuple[float, float]]
    scatter_sigma_m: float
    capture_about_centroid: float
    capture_absolute: float
    toa_centroid_drift_s: float
    toa_scatter_rms_s: float
    plate_radius_m: float = PLATE_RADIUS_M
    scatter_sigma_max_m: float = CAPTURE_SIGMA_MAX_M

    @property
    def centroid_drift_m(self) -> float:
        """The common-mode lateral shift the plane's centroid retarget must absorb."""
        return float(math.hypot(self.centroid_m[0], self.centroid_m[1]))

    @property
    def retarget_ok(self) -> bool:
        return self.centroid_drift_m <= self.centroid_retarget_m

    @property
    def scatter_sigma_ok(self) -> bool:
        """Per-axis scatter σ within the ADR 0015 plate-capture criterion (1.65 m)."""
        return self.scatter_sigma_m <= self.scatter_sigma_max_m


def summarize_train_capture(
    misses: Sequence[PlateMiss], spec: TrainDispersionSpec
) -> TrainCaptureStats:
    """Reduce a train's plate-frame arrivals into the centroid-vs-scatter verdict (ADR 0016).

    The centroid is the common-mode shift (absorbed by the plane's ±2 km retarget); each unit is
    then *re-centered* on the centroid (where the plane aims) before judging plate capture, so the
    scatter alone faces the plate.  ``capture_absolute`` (no retarget) is reported for contrast.
    """
    lateral: NDArray[np.float64] = np.array(
        [m.lateral_m for m in misses], dtype=np.float64
    ).reshape(-1, 2)
    toa: NDArray[np.float64] = np.array([m.toa_error_s for m in misses], dtype=np.float64)
    n = lateral.shape[0]

    centroid = lateral.mean(axis=0)
    toa_centroid = float(toa.mean())
    recentered = lateral - centroid
    recentered_misses = [
        PlateMiss(lateral_m=(float(d[0]), float(d[1])), toa_error_s=float(t - toa_centroid))
        for d, t in zip(recentered, toa, strict=True)
    ]

    cov = np.cov(recentered, rowvar=False) if n > 1 else np.zeros((2, 2))
    # Per-axis σ for an isotropic 2D model: E[|d|²] = 2σ² (matches the ADR 0015 Rayleigh σ).
    scatter_sigma = math.sqrt(float(np.mean(np.sum(recentered**2, axis=1))) / 2.0)
    toa_scatter_rms = math.sqrt(float(np.mean((toa - toa_centroid) ** 2)))

    return TrainCaptureStats(
        n_units=n,
        centroid_m=(float(centroid[0]), float(centroid[1])),
        centroid_retarget_m=spec.centroid_retarget_m,
        scatter_cov_m2=(
            (float(cov[0, 0]), float(cov[0, 1])),
            (float(cov[1, 0]), float(cov[1, 1])),
        ),
        scatter_sigma_m=scatter_sigma,
        capture_about_centroid=capture_fraction(recentered_misses),
        capture_absolute=capture_fraction(list(misses)),
        toa_centroid_drift_s=toa_centroid,
        toa_scatter_rms_s=toa_scatter_rms,
    )


def format_train_capture(stats: TrainCaptureStats) -> str:
    """One-screen train-mode capture report — centroid drift (retarget) vs scatter (plate)."""
    lines = [
        f"Train-mode capture — N={stats.n_units} units, centroid vs scatter (ADR 0016)",
        f"  Centroid drift: {stats.centroid_drift_m:.1f} m"
        f" vs ±{stats.centroid_retarget_m:.0f} m retarget"
        f" — {'absorbed' if stats.retarget_ok else 'EXCEEDS plane capability'}"
        f" (ToA {stats.toa_centroid_drift_s * 1e3:+.2f} ms).",
        f"  Scatter about centroid: σ {stats.scatter_sigma_m:.2f} m"
        f" vs ≤{stats.scatter_sigma_max_m} m criterion"
        f" — {'OK' if stats.scatter_sigma_ok else 'TOO WIDE'}"
        f" (ToA RMS {stats.toa_scatter_rms_s * 1e3:.2f} ms vs ≤{TOA_LIMIT_S * 1e3:g} ms).",
        f"  P(capture) about centroid: {stats.capture_about_centroid * 100:.1f}%"
        f" on the {stats.plate_radius_m:g} m plate"
        f" (absolute, no retarget: {stats.capture_absolute * 100:.1f}%).",
    ]
    return "\n".join(lines)


@dataclass(frozen=True)
class TrainEnsembleFinding:
    """The D1.1 verdict surface: capture + propellant + perigee over a train ensemble."""

    capture: TrainCaptureStats
    terminal_dv_mean_m_s: float
    terminal_dv_max_m_s: float
    midcourse_dv_m_s: float
    propellant: tuple[PropellantPoint, ...]
    perigee_mean_m: float
    perigee_min_m: float
    perigee_max_m: float

    @property
    def total_dv_worst_m_s(self) -> float:
        """Worst-unit mission Δv = the characterized midcourse + the flown terminal aim."""
        return self.midcourse_dv_m_s + self.terminal_dv_max_m_s

    @property
    def within_budget(self) -> bool:
        """The worst unit clears the <2 % line at the most conservative Isp anchor (50 s)."""
        return self.propellant[0].within_budget


def summarize_train_ensemble(
    misses: Sequence[PlateMiss],
    plans: Sequence[FeedforwardPlan],
    perigees_m: Sequence[float],
    spec: TrainDispersionSpec,
) -> TrainEnsembleFinding:
    """Assemble the D1.1 finding: the train-relative capture + propellant ledger + perigee."""
    capture = summarize_train_capture(misses, spec)
    terminal_dv = np.array([p.dv_m_s for p in plans], dtype=np.float64)
    perigee = np.array(perigees_m, dtype=np.float64)
    dv_max = float(terminal_dv.max())
    return TrainEnsembleFinding(
        capture=capture,
        terminal_dv_mean_m_s=float(terminal_dv.mean()),
        terminal_dv_max_m_s=dv_max,
        midcourse_dv_m_s=spec.midcourse_dv_m_s,
        propellant=propellant_curve(spec.midcourse_dv_m_s + dv_max),
        perigee_mean_m=float(perigee.mean()),
        perigee_min_m=float(perigee.min()),
        perigee_max_m=float(perigee.max()),
    )


def format_train_ensemble(finding: TrainEnsembleFinding) -> str:
    """One-screen D1.1 report: capture + propellant (<2 %) + perigee (deorbit diagnostic)."""
    worst = finding.propellant[0]
    lines = [
        "Rung D / D1.1 — closed-loop train ensemble (Φ-composed entry + flown terminal; ADR 0018)",
        format_train_capture(finding.capture),
        f"  Propellant: terminal aim Δv mean {finding.terminal_dv_mean_m_s:.2f}"
        f" / max {finding.terminal_dv_max_m_s:.2f} m/s; worst mission Δv"
        f" {finding.total_dv_worst_m_s:.2f} m/s (+ {finding.midcourse_dv_m_s:g} midcourse)"
        f" → {worst.fraction * 100:.2f}% @ Isp {worst.isp_s:g} s"
        f" — {'under' if finding.within_budget else 'OVER'} the 2% budget.",
        f"  Perigee (diagnostic, low=good): mean {finding.perigee_mean_m / 1e3:.1f} km"
        f" [min {finding.perigee_min_m / 1e3:.1f}, max {finding.perigee_max_m / 1e3:.1f}].",
    ]
    return "\n".join(lines)


@dataclass(frozen=True)
class PooledCaptureFinding:
    """Per-unit capture pooled over several *trains* — the statistically-honest D1 headline.

    A single train shares one common-mode draw, so pooling units *within* one train samples the
    per-unit legs only; the escape tail is driven by the shared entry leg, which needs several
    trains to sample. ``lower_bound_95`` is the one-sided Clopper–Pearson bound: the point
    estimate is what the counting yields, this is what it *establishes* (ADR 0018, 2026-07-22).
    """

    n_units: int
    n_trains: int
    escapes: int
    entry_escape_threshold_m: float
    max_entry_m: float
    core_sigma_m: float
    plate_radius_m: float = PLATE_RADIUS_M

    @property
    def capture(self) -> float:
        return 1.0 - self.escapes / self.n_units

    @property
    def lower_bound_95(self) -> float:
        return binomial_lower_bound_95(self.n_units - self.escapes, self.n_units)

    @property
    def entry_limited(self) -> bool:
        """Every escape is a hand-off entry past the funnel, not a terminal-noise event."""
        return self.escapes == 0 or self.max_entry_m >= self.entry_escape_threshold_m

    @property
    def meets_bound(self) -> bool:
        """The counting *establishes* ≥99 %, rather than merely estimating it."""
        return self.lower_bound_95 >= 0.99


def binomial_lower_bound_95(successes: int, trials: int, alpha: float = 0.05) -> float:
    """One-sided Clopper–Pearson lower bound on the success probability."""
    if trials <= 0:
        return 0.0
    if successes >= trials:
        return float(alpha ** (1.0 / trials))
    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        tail = sum(
            math.comb(trials, i) * mid**i * (1.0 - mid) ** (trials - i)
            for i in range(successes, trials + 1)
        )
        lo, hi = (mid, hi) if tail < alpha else (lo, mid)
    return lo


def summarize_pooled_capture(
    misses: Sequence[PlateMiss],
    entries: Sequence[Vec2],
    n_trains: int,
    plate_radius_m: float = PLATE_RADIUS_M,
) -> PooledCaptureFinding:
    """Reduce a multi-train batch to the counted capture + the escape/entry attribution.

    ``core_sigma_m`` is the per-axis σ of the *captured* arrivals only: the pooled set is bimodal
    (captured core plus funnel escapes), so a σ over everything is a mixture statistic that must
    not be read against the ADR 0015 criterion.
    """
    if len(misses) != len(entries):
        raise ValueError("misses and entries must be paired")
    radii = [math.hypot(*m.lateral_m) for m in misses]
    entry_mags = [math.hypot(*e) for e in entries]
    escaped = [mag for r, mag in zip(radii, entry_mags, strict=True) if r > plate_radius_m]
    captured = [m for m, r in zip(misses, radii, strict=True) if r <= plate_radius_m]
    sq = sum(m.lateral_m[0] ** 2 + m.lateral_m[1] ** 2 for m in captured)
    return PooledCaptureFinding(
        n_units=len(misses),
        n_trains=n_trains,
        escapes=len(escaped),
        entry_escape_threshold_m=min(escaped) if escaped else math.inf,
        max_entry_m=max(entry_mags) if entry_mags else 0.0,
        core_sigma_m=math.sqrt(sq / (2.0 * len(captured))) if captured else 0.0,
        plate_radius_m=plate_radius_m,
    )


def format_pooled_capture(finding: PooledCaptureFinding) -> str:
    """One-screen pooled-capture report: the point estimate vs the bound it establishes."""
    threshold = (
        f"{finding.entry_escape_threshold_m:.0f} m"
        if math.isfinite(finding.entry_escape_threshold_m)
        else "n/a"
    )
    return "\n".join(
        [
            f"  Pooled per-unit capture: {finding.capture * 100:.2f}%"
            f" ({finding.escapes} escapes / {finding.n_units} units,"
            f" {finding.n_trains} trains, {finding.plate_radius_m:g} m plate).",
            f"  One-sided 95% lower bound: {finding.lower_bound_95 * 100:.2f}%"
            f" — {'establishes' if finding.meets_bound else 'does NOT establish'} ≥99%.",
            f"  Captured-core σ (per axis): {finding.core_sigma_m:.2f} m"
            f" — the pooled σ over escapes too would be a mixture statistic, not a Gaussian σ.",
            f"  Escapes are entry-limited: {'yes' if finding.entry_limited else 'NO'}"
            f" (lowest escaping entry {threshold}; max entry"
            f" {finding.max_entry_m:.0f} m).",
        ]
    )
