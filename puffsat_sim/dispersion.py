"""Open-loop dispersion harness — pure sampling, RTN geometry, and statistics.

JVM-free core of the Monte Carlo / capstone harness (design doc §13, ADR 0002):
the per-run input draws, the RTN local-frame math, and the ensemble statistics.
The propagate-and-record loop that consumes these lives in
:mod:`puffsat_sim.montecarlo`; keeping this module pure lets it be unit-tested
without booting Orekit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Vec3 = tuple[float, float, float]
Basis = tuple[Vec3, Vec3, Vec3]


@dataclass(frozen=True)
class DispersionSpec:
    """Swept knobs for one ensemble: nominal values + per-input 1σ.

    Defaults match the reference mission and preset coefficients.  Injection-error
    σ are per-RTN-axis [m/s] (transverse is the dr_p/dv_a lever, §8); coefficient
    and space-weather σ are fractional (coefficient of variation) for the
    median-nominal multiplicative log-normals.
    """

    sigma_dv_radial_m_s: float = 0.1
    sigma_dv_transverse_m_s: float = 0.1
    sigma_dv_normal_m_s: float = 0.1
    cd_area_over_mass: float = 0.04
    cr_area_over_mass: float = 0.02
    sigma_cd_frac: float = 0.20
    sigma_cr_frac: float = 0.20
    f10p7: float = 150.0
    ap: float = 15.0
    sigma_f10p7_frac: float = 0.15
    sigma_ap_frac: float = 0.50


@dataclass(frozen=True)
class RunInputs:
    """One run's sampled draws plus its run index (replay needs index + master seed)."""

    run_index: int
    dv_rtn_m_s: Vec3
    cd_area_over_mass: float
    cr_area_over_mass: float
    f10p7: float
    ap: float


def _lognormal_factor(rng: np.random.Generator, cv: float) -> float:
    """Multiplicative factor with median 1 and coefficient of variation ``cv``.

    ln(factor) ~ N(0, s) with s = √(ln(1+cv²)), so the median is exp(0)=1 and the
    distribution is right-skewed — a small positive mean bias (exp(s²/2)), by
    design (ADR 0002), which the capstone reports rather than hides.
    """
    if cv <= 0.0:
        return 1.0
    s = math.sqrt(math.log(1.0 + cv * cv))
    return float(math.exp(rng.normal(0.0, s)))


def sample_run_inputs(rng: np.random.Generator, spec: DispersionSpec, run_index: int) -> RunInputs:
    """Draw one run's inputs in a fixed order (so a seed reproduces the run exactly)."""
    dv_r = float(rng.normal(0.0, spec.sigma_dv_radial_m_s))
    dv_t = float(rng.normal(0.0, spec.sigma_dv_transverse_m_s))
    dv_n = float(rng.normal(0.0, spec.sigma_dv_normal_m_s))
    cd = spec.cd_area_over_mass * _lognormal_factor(rng, spec.sigma_cd_frac)
    cr = spec.cr_area_over_mass * _lognormal_factor(rng, spec.sigma_cr_frac)
    f10p7 = spec.f10p7 * _lognormal_factor(rng, spec.sigma_f10p7_frac)
    ap = spec.ap * _lognormal_factor(rng, spec.sigma_ap_frac)
    return RunInputs(run_index, (dv_r, dv_t, dv_n), cd, cr, f10p7, ap)


def _to_vec3(a: NDArray[np.float64]) -> Vec3:
    return (float(a[0]), float(a[1]), float(a[2]))


def rtn_basis(position_m: Vec3, velocity_m_s: Vec3) -> Basis:
    """Orthonormal RTN basis (radial, transverse, normal) from a Cartesian state.

    R = r̂; N = (r×v)̂ (orbit normal); T = N×R (in-plane, toward motion).  At
    apogee the radial velocity is zero, so T coincides with the velocity direction
    — the tangential dr_p/dv_a lever.
    """
    r = np.asarray(position_m, dtype=np.float64)
    v = np.asarray(velocity_m_s, dtype=np.float64)
    r_hat = r / np.linalg.norm(r)
    n_hat = np.cross(r, v)
    n_hat = n_hat / np.linalg.norm(n_hat)
    t_hat = np.cross(n_hat, r_hat)
    return _to_vec3(r_hat), _to_vec3(t_hat), _to_vec3(n_hat)


def rtn_components(vector_m: Vec3, basis: Basis) -> Vec3:
    """Project a Cartesian vector onto an RTN basis → (radial, transverse, normal)."""
    vec = np.asarray(vector_m, dtype=np.float64)
    r_hat = np.asarray(basis[0], dtype=np.float64)
    t_hat = np.asarray(basis[1], dtype=np.float64)
    n_hat = np.asarray(basis[2], dtype=np.float64)
    return (float(vec @ r_hat), float(vec @ t_hat), float(vec @ n_hat))


def rtn_to_cartesian(components_rtn: Vec3, basis: Basis) -> Vec3:
    """Combine RTN components against an RTN basis → a Cartesian vector."""
    c_r, c_t, c_n = components_rtn
    r_hat = np.asarray(basis[0], dtype=np.float64)
    t_hat = np.asarray(basis[1], dtype=np.float64)
    n_hat = np.asarray(basis[2], dtype=np.float64)
    return _to_vec3(c_r * r_hat + c_t * t_hat + c_n * n_hat)


@dataclass(frozen=True)
class EnsembleStats:
    """Aggregate dispersion statistics for one ensemble (lengths in m, times in s)."""

    n: int
    miss_rtn_mean_m: Vec3  # bias (R,T,N) — fold into the open-loop aim
    miss_rtn_std_m: Vec3  # per-axis spread the controller must shrink
    miss_rtn_cov_m2: Basis  # 3×3 covariance (the dispersion ellipsoid)
    toa_miss_mean_s: float
    toa_miss_std_s: float
    perigee_alt_mean_m: float
    perigee_alt_std_m: float
    perigee_alt_min_m: float
    perigee_alt_max_m: float
    total_dv_mean_m_s: float  # the Rung A Δv-floor ledger (0 when control is off)
    total_dv_std_m_s: float
    total_dv_max_m_s: float
    converged_fraction: float  # fraction of runs the corrector solved (1.0 open-loop)


def summarize(
    miss_rtn_m: NDArray[np.float64],
    toa_miss_s: NDArray[np.float64],
    perigee_alt_m: NDArray[np.float64],
    total_dv_m_s: NDArray[np.float64],
    converged: NDArray[np.bool_],
) -> EnsembleStats:
    """Aggregate per-run arrays into ensemble statistics.

    ``miss_rtn_m`` is (N, 3) in the nominal-crossing RTN frame; its mean is the aim
    bias and its covariance the dispersion ellipsoid.  perigee_alt min/max bound the
    debris-disposal-safety margin (a missed PuffSat must stay low enough to deorbit).
    ``total_dv_m_s`` and ``converged`` are per-run control outcomes (all-zero / all-True
    for the open-loop capstone).
    """
    n = int(miss_rtn_m.shape[0])
    mean = miss_rtn_m.mean(axis=0)
    std = miss_rtn_m.std(axis=0, ddof=1) if n > 1 else np.zeros(3)
    cov = np.cov(miss_rtn_m, rowvar=False) if n > 1 else np.zeros((3, 3))
    return EnsembleStats(
        n=n,
        miss_rtn_mean_m=_to_vec3(mean),
        miss_rtn_std_m=_to_vec3(std),
        miss_rtn_cov_m2=(_to_vec3(cov[0]), _to_vec3(cov[1]), _to_vec3(cov[2])),
        toa_miss_mean_s=float(toa_miss_s.mean()),
        toa_miss_std_s=float(toa_miss_s.std(ddof=1)) if n > 1 else 0.0,
        perigee_alt_mean_m=float(perigee_alt_m.mean()),
        perigee_alt_std_m=float(perigee_alt_m.std(ddof=1)) if n > 1 else 0.0,
        perigee_alt_min_m=float(perigee_alt_m.min()),
        perigee_alt_max_m=float(perigee_alt_m.max()),
        total_dv_mean_m_s=float(total_dv_m_s.mean()),
        total_dv_std_m_s=float(total_dv_m_s.std(ddof=1)) if n > 1 else 0.0,
        total_dv_max_m_s=float(total_dv_m_s.max()),
        converged_fraction=float(converged.mean()),
    )
