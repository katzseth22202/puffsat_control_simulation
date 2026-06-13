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

# The plane's launch-window retarget capability (CONTEXT: Centroid retarget; ADR 0006/0016).
CENTROID_RETARGET_M: float = 2000.0


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
