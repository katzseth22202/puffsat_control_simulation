"""Pure importance-sampling tail P(capture) core — a Rung-D / D1.x sub-slice (ADR 0018, no JVM).

D1.1's headline P(capture) was a **16-unit empirical** ("94 % about-centroid / 100 % absolute"),
which cannot resolve the figure of merit: at the achievable ~3.2 µrad grade the per-unit arrival
scatter is σ ≈ 1.35 m against a 5 m plate, so a capture failure is a **~3.7 σ tail event**
(p ≈ 10⁻³) whose heavy-tailed physics (the C3b catch-radius cliff, the significance-gate noise
rectification) is exactly what a small empirical batch misses.  ADR 0018 decision 6 resolves it by
**importance sampling, validated by a brute-force batch** that confirms the reweighting is unbiased
— LinCov is the IS-proposal *designer* (the entry covariance + cliff location set the inflation),
never a replacement for the flown tail.

This module is the pure side: the IS proposal over the 2-D per-unit hand-off entry offset (the
per-unit *scatter* the funnel must null — :data:`puffsat_sim.train.ENTRY_LATERAL_PERUNIT_M`), the
exact likelihood-ratio weights, the weighted tail estimator with its IS variance / CI / effective
sample size, and the brute-force-vs-IS agreement reduction.  The tracker noise is sampled *fresh*
per trajectory by the JVM glue (:mod:`puffsat_sim.runs.tail_capture`), so it never enters the
weight — the estimator is unbiased for P(escape) marginalized over both entry and noise.  Truth
physics are held at nominal (drag/SRP were ≤ 0.19 m in D1.1 — 2nd-order), isolating the binding
entry×noise driver; Cr/storm importance sampling is a noted deferred extension.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.guidance import CAPTURE_SIGMA_MAX_M, PLATE_RADIUS_M, TOA_LIMIT_S
from puffsat_sim.train import ENTRY_LATERAL_PERUNIT_M

# A 2-D isotropic entry magnitude r has E[r²] = 2σ_axis², so the per-axis σ is the RMS magnitude
# over √2 (the same mapping train.py uses for the sampled entry offsets).
_PER_AXIS_FROM_MAGNITUDE: float = 1.0 / math.sqrt(2.0)

# Two-sided 95 % normal quantile for the Wald intervals (the project avoids a scipy dependency; a
# fixed 95 % CI is adequate for a tail-probability report, cf. estimation._normal_quantile).
_Z_95: float = 1.959963984540054


@dataclass(frozen=True)
class TailCaptureSpec:
    """Swept knobs for one IS tail estimate.

    ``entry_lateral_m`` is the per-unit hand-off entry RMS magnitude (D1.1's 141 m C1 nav leg);
    ``kappa`` inflates the proposal σ (the LinCov-informed designer knob — a *gentle* κ ≈ 1.35
    oversamples the escape ring without pushing draws past the ~450 m funnel edge, where the
    saturation catastrophe blows up the IS variance; aggressive κ is counterproductive);
    ``r_val_m`` is the shallow validation radius (where brute force still sees events) the IS
    reweighting is checked against; ``n_units`` drives the train read; ``target_capture`` is the
    P(capture) the lower CI bound is graded against; ``rel_error_target`` sets when the IS tail
    estimate counts as resolved; ``control_period_s`` is the C3b control clock the loop is flown on.
    """

    entry_lateral_m: float = ENTRY_LATERAL_PERUNIT_M
    kappa: float = 1.35
    plate_radius_m: float = PLATE_RADIUS_M
    toa_limit_s: float = TOA_LIMIT_S
    r_val_m: float = 3.0
    n_units: int = 16
    target_capture: float = 0.99
    rel_error_target: float = 0.5
    control_period_s: float = 1.0

    @property
    def sigma_axis_m(self) -> float:
        return self.entry_lateral_m * _PER_AXIS_FROM_MAGNITUDE


def sample_entry_nominal(
    rng: np.random.Generator, spec: TailCaptureSpec, n: int
) -> NDArray[np.float64]:
    """Draw ``n`` per-unit 2-D entry offsets from the *nominal* distribution (brute-force batch).

    The funnel-relative per-unit scatter the terminal loop must null — an isotropic 2-D Gaussian at
    ``spec.sigma_axis_m`` per axis, matching ``train.sample_train_entry_offsets``' per-unit leg.
    """
    return rng.normal(0.0, spec.sigma_axis_m, size=(n, 2))


def sample_entry_is(
    rng: np.random.Generator, spec: TailCaptureSpec, n: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Draw ``n`` entry offsets from the variance-inflated IS proposal with their weights.

    The proposal is the nominal isotropic 2-D Gaussian widened by ``spec.kappa``; the exact
    likelihood ratio ``w = p/q = κ²·exp(−(r²/2σ²)(1 − 1/κ²))`` reweights each draw back to the
    nominal measure (downweighting the inflated tail draws), so the weighted indicator is an
    unbiased estimator of the nominal tail probability.
    """
    sigma = spec.sigma_axis_m
    kappa = spec.kappa
    entries = rng.normal(0.0, kappa * sigma, size=(n, 2))
    r_sq = np.sum(entries**2, axis=1)
    weights = kappa**2 * np.exp(-(r_sq / (2.0 * sigma**2)) * (1.0 - 1.0 / kappa**2))
    return entries, weights


def escaped_flags(
    lateral_norm_m: NDArray[np.float64],
    toa_error_s: NDArray[np.float64],
    radius_m: float,
    toa_limit_s: float,
) -> NDArray[np.bool_]:
    """Per-arrival plate-escape flag (the complement of ``guidance.capture_fraction``)."""
    lateral = np.asarray(lateral_norm_m, dtype=np.float64)
    toa = np.asarray(toa_error_s, dtype=np.float64)
    return (lateral > radius_m) | (np.abs(toa) > toa_limit_s)


@dataclass(frozen=True)
class TailEstimate:
    """A weighted tail-probability estimate with its IS variance, CI, and effective sample size."""

    probability: float
    std_error: float
    ci_low: float
    ci_high: float
    ess: float
    weight_sum_ratio: float
    n: int

    @property
    def relative_error(self) -> float:
        return self.std_error / self.probability if self.probability > 0.0 else math.inf


def weighted_tail_probability(
    weights: NDArray[np.float64], escaped: NDArray[np.bool_]
) -> TailEstimate:
    """Importance-sampling estimator of P(escape) from per-draw weights and escape flags.

    ``p̂ = mean(wᵢ·yᵢ)`` is unbiased because ``E_q[w·y] = E_p[y] = p``; the variance is the sample
    variance of the per-draw contributions ``wᵢ·yᵢ`` over ``n``.  Brute force is the ``w ≡ 1`` case
    (``p̂`` = the plain fraction, the Wald variance ``p(1−p)/n``), so this serves both batches.
    """
    w = np.asarray(weights, dtype=np.float64)
    y = np.asarray(escaped, dtype=np.bool_).astype(np.float64)
    n = int(w.shape[0])
    contributions = w * y
    p_hat = float(np.clip(contributions.mean(), 0.0, 1.0))
    variance = float(contributions.var(ddof=1)) / n if n > 1 else 0.0
    std_error = math.sqrt(variance)
    half = _Z_95 * std_error
    ess = float(w.sum() ** 2 / np.sum(w**2)) if np.sum(w**2) > 0.0 else 0.0
    return TailEstimate(
        probability=p_hat,
        std_error=std_error,
        ci_low=max(0.0, p_hat - half),
        ci_high=min(1.0, p_hat + half),
        ess=ess,
        weight_sum_ratio=float(w.mean()),
        n=n,
    )


def _rayleigh_escape(radius_m: float, sigma_axis_m: float) -> float:
    """Gaussian/Rayleigh tail ``exp(−r²/2σ²)`` — what the σ ≤ 1.65 m criterion *implies*.

    The ADR 0015 capture criterion is a per-axis-σ bound on the arrival scatter, read through an
    isotropic 2-D Gaussian (Rayleigh magnitude); this is the escape probability that σ predicts,
    against which the flown tail's heaviness is measured.
    """
    if sigma_axis_m <= 0.0:
        return 0.0
    return math.exp(-(radius_m**2) / (2.0 * sigma_axis_m**2))


@dataclass(frozen=True)
class TailCaptureFinding:
    """The D1.x tail verdict: a tail-resolved per-unit P(capture) + CI, IS-validated unbiased.

    At the achievable grade the capture-failure event is *shallow* (a few %, not a deep rare
    event), so brute force is the robust plate-tail estimator (``bf_plate``) and IS is the
    validated cross-check (agreement at ``r_val``) plus the heavier-than-Gaussian quantifier and
    the tool reserved for the deeper tails a tighter grade/plate would create.
    """

    spec: TailCaptureSpec
    is_plate: TailEstimate
    is_rval: TailEstimate
    bf_plate: TailEstimate
    bf_rval: TailEstimate
    scatter_sigma_m: float

    @property
    def p_escape(self) -> float:
        """The robust per-unit plate-escape probability (brute force — see the class note)."""
        return self.bf_plate.probability

    @property
    def p_capture(self) -> float:
        """Per-unit P(capture) — also the expected fraction of an iid train that connects."""
        return 1.0 - self.bf_plate.probability

    @property
    def p_capture_ci(self) -> tuple[float, float]:
        """The per-unit P(capture) interval (the escape CI flipped: lower = 1 − escape upper)."""
        return (1.0 - self.bf_plate.ci_high, 1.0 - self.bf_plate.ci_low)

    @property
    def gaussian_escape(self) -> float:
        """The escape probability the measured arrival σ predicts under a Gaussian tail."""
        return _rayleigh_escape(self.spec.plate_radius_m, self.scatter_sigma_m)

    @property
    def tail_excess_factor(self) -> float:
        """How much heavier the flown tail is than its Gaussian extrapolation (> 1 = heavier)."""
        g = self.gaussian_escape
        return self.p_escape / g if g > 0.0 else math.inf

    @property
    def is_agrees_at_plate(self) -> bool:
        """The (fragile) IS plate estimate brackets the brute-force one — a consistency read."""
        return (
            self.is_plate.ci_low <= self.bf_plate.ci_high
            and self.bf_plate.ci_low <= self.is_plate.ci_high
        )

    @property
    def validated(self) -> bool:
        """IS and brute force agree at ``r_val`` (CIs overlap) — the reweighting is unbiased."""
        return (
            self.is_rval.ci_low <= self.bf_rval.ci_high
            and self.bf_rval.ci_low <= self.is_rval.ci_high
        )

    @property
    def tail_resolved(self) -> bool:
        """The headline (brute-force) plate estimate is resolved: events seen and tight enough."""
        return (
            self.bf_plate.probability > 0.0
            and self.bf_plate.relative_error <= self.spec.rel_error_target
        )

    @property
    def meets_criterion(self) -> bool:
        """The per-unit P(capture) *lower* CI bound clears the design target (the honest read)."""
        return self.p_capture_ci[0] >= self.spec.target_capture


def _scatter_sigma_axis_m(lateral_norm_m: NDArray[np.float64]) -> float:
    """Per-axis arrival σ for an isotropic 2-D model: ``√(⟨|d|²⟩/2)`` (train.py's convention)."""
    norms = np.asarray(lateral_norm_m, dtype=np.float64)
    return math.sqrt(float(np.mean(norms**2)) / 2.0) if norms.size else 0.0


def summarize_tail_capture(
    is_weights: NDArray[np.float64],
    is_lateral_norm_m: NDArray[np.float64],
    is_toa_error_s: NDArray[np.float64],
    bf_lateral_norm_m: NDArray[np.float64],
    bf_toa_error_s: NDArray[np.float64],
    spec: TailCaptureSpec,
) -> TailCaptureFinding:
    """Reduce an IS batch + a brute-force batch into the tail-resolved P(capture) finding.

    ``is_*`` are the flown IS-batch arrivals (with their entry weights); ``bf_*`` the nominal
    brute-force arrivals (implicit unit weight).  Both are read at the plate (the deliverable —
    brute force is the robust estimator at this shallow depth, IS the cross-check) and at ``r_val``
    (where the IS reweighting is validated against the brute-force fraction).
    """
    is_w = np.asarray(is_weights, dtype=np.float64)

    def _is(radius_m: float) -> TailEstimate:
        return weighted_tail_probability(
            is_w, escaped_flags(is_lateral_norm_m, is_toa_error_s, radius_m, spec.toa_limit_s)
        )

    def _bf(radius_m: float) -> TailEstimate:
        flags = escaped_flags(bf_lateral_norm_m, bf_toa_error_s, radius_m, spec.toa_limit_s)
        return weighted_tail_probability(np.ones(flags.shape[0], dtype=np.float64), flags)

    return TailCaptureFinding(
        spec=spec,
        is_plate=_is(spec.plate_radius_m),
        is_rval=_is(spec.r_val_m),
        bf_plate=_bf(spec.plate_radius_m),
        bf_rval=_bf(spec.r_val_m),
        scatter_sigma_m=_scatter_sigma_axis_m(bf_lateral_norm_m),
    )


def format_tail_capture(finding: TailCaptureFinding) -> str:
    """One-screen D1.x tail report — validation, the resolved P(capture), heavier-than-Gaussian."""
    spec = finding.spec
    cap_lo, cap_hi = finding.p_capture_ci
    lines = [
        "Rung D / D1.x — importance-sampling tail P(capture) (ADR 0018 decision 6)",
        f"  Setup: per-unit entry {spec.entry_lateral_m:.0f} m, achievable grade; arrival scatter"
        f" σ {finding.scatter_sigma_m:.2f} m (criterion ≤ {CAPTURE_SIGMA_MAX_M} m). IS proposal"
        f" ×{spec.kappa:g} (σ_axis → {spec.kappa * spec.sigma_axis_m:.0f} m),"
        f" ESS {finding.is_plate.ess:.0f}, Σw/N {finding.is_plate.weight_sum_ratio:.2f}.",
        f"  Validation @ r_val {spec.r_val_m:g} m: IS P(escape)"
        f" {finding.is_rval.probability * 100:.2f}% [{finding.is_rval.ci_low * 100:.2f},"
        f" {finding.is_rval.ci_high * 100:.2f}] vs brute force"
        f" {finding.bf_rval.probability * 100:.2f}% [{finding.bf_rval.ci_low * 100:.2f},"
        f" {finding.bf_rval.ci_high * 100:.2f}] (N={finding.bf_rval.n})"
        f" — {'unbiased' if finding.validated else 'INCONSISTENT'}.",
        f"  Tail @ plate {spec.plate_radius_m:g} m: P(escape) {finding.p_escape * 100:.2f}%"
        f" (±{finding.bf_plate.std_error * 100:.2f}%, N={finding.bf_plate.n})"
        f" — {'resolved' if finding.tail_resolved else 'UNDER-resolved (raise N)'};"
        f" {finding.tail_excess_factor:.1f}× the Gaussian"
        f" {finding.gaussian_escape * 100:.2f}% the σ-criterion implies.",
        f"  Per-unit P(capture): {finding.p_capture * 100:.2f}%"
        f" [{cap_lo * 100:.2f}, {cap_hi * 100:.2f}] vs target {spec.target_capture * 100:g}%"
        f" — {'MEETS' if finding.meets_criterion else 'BELOW target'}"
        " (lower bound; = expected fraction of the train that connects, the rest deorbit by §9).",
    ]
    return "\n".join(lines)
