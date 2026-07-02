"""Pure surveyor-anchored centering budget — the Tier-2 ~10 cm sizing module (no JVM).

The closed-loop Monte Carlo (Rung D / D1, ADR 0018) closed out **Tier 1**: a PuffSat captures a
**5 m** pusher plate at ≥99 % per-unit confidence.  The white paper carries a deferred **Tier 2**
tightening (CONCLUSION.md; CONTEXT.md *Surveyor-anchored centering*): shrink the plate toward
**~10 cm** by driving the *per-unit* arrival miss to centimetres.  That regime is **not** a control
problem — cm trims sit deep inside the ~475 m catch-radius funnel, so authority is not the binder
and D2/MPC does not return.  It is **knowledge/metrology-limited**, and its binding numbers are two
*hardware characterizations* a Monte Carlo cannot produce.  This module is the right rigour for
that: a **pure sizing module** (the ``tracker_budget`` / ``apogee_nav`` / ``distortion_field``
precedent) that turns those characterizations into a plate-size sensitivity curve — the ``"if we
get hardware that does X, we hit N cm"`` statement, with reference-hardware anchors for realism.

Two levers on two axes compose the per-unit arrival σ, then the committed 99 % capture criterion
(2-D Rayleigh, ``plate = (R/σ)·σ`` with ``R/σ = 3.03`` from :mod:`puffsat_sim.guidance`) turns it
into a plate radius:

* **Block shift (common-mode bias).**  A sacrificial "surveyor" PuffSat whose true crossing is read
  by an **independent one-shot hoop** (lidar/microwave trilateration — *not* the optical tracker,
  which would be circular against its own bias) pins the swarm's quasi-static optical-distortion
  *bias* to the plate.  The residual is the hoop precision **σ_hoop** (1σ), inherited by every unit.

* **Per-unit scatter (independent).**  Strobed known-pattern beacons make the swarm
  **camera-rigid**: the surveyor centroids each follower's bearing, and that knowledge is only as
  good as the camera angular grade times the intra-train link range, **σ_δ = σ_θ·v/f** (the
  committed CONTEXT.md formula — the plate-crossing cadence ``f`` sets the nearest-follower spacing
  ``v/f``, so a faster train is a closer, sharper link).  σ_θ is the RSS of photon-limited
  centroiding (reused from :mod:`puffsat_sim.tracker_budget`) and the calibrated distortion floor,
  reduced across ``N`` cameras by the equicorrelated law reused from ``distortion_field``.

The photon term is the one place a naïve read goes wrong: at the km-class intra-train range a
*gram-scale* 5 mm aperture on a naive **wide-cone 1 W CW** beacon is **photon-limited** (~4 µrad,
over the 3 µrad distortion floor), not distortion-limited, so the thesis would fail.  Two
independent levers each clear it, and the design uses both for margin:

* a **Q-switched beacon** — bright (~100 kW peak) ns pulses at known timings, a few hundred mW
  *average* — read in a gate matched to the pulse.  The measurement collects the whole pulse
  *energy* (avg ÷ rep-rate), so a 30 mJ pulse dumps ~10⁵–10⁸ photons through the tiny aperture; the
  ns gate also freezes the motion-smear term to nothing, and a narrowband filter matched to the
  line keeps the link signal-dominated against stray light.
* a **somewhat directional** beam — coarse-pointed a few degrees along the train axis toward the
  surveyor (no fine tracking), which buys ~an order of magnitude of intensity over the wide cone.

Together they put the photon term ~an order of magnitude under the distortion floor (≈0.2 µrad), so
σ_θ is set by the *calibrated distortion floor* with comfortable margin.  The surveyor's backward
look at dark sky (no forward flash/glare) keeps the common smear floor at zero — the swarm is
honestly distortion-limited.  (This reuses the Q-switched 1064 nm lineage the near-Sun extension
carries; CONTEXT.md *Near-Sun optical nav*.)

That ρ knob on the array answers the "does 3-camera voting kill the bias" question directly: at ρ=1
(identical cameras) voting does *not* touch the shared distortion; only physically diverse cameras
(ρ→0) or differential astrometry cut it.

Paper-side, reversible (ADR 0022).  The module does not change any committed Tier-1 number; it sizes
the deferred Tier-2 claim, and by construction stays *sized*, not simulated: its binders are bench
characterizations a Monte Carlo cannot produce, and its dynamics sit inside the ~475 m funnel.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field, replace

from puffsat_sim.distortion_field import DetectorBudget, array_sigma_theta_rad
from puffsat_sim.guidance import CAPTURE_SIGMA_MAX_M, PLATE_RADIUS_M
from puffsat_sim.tracker_budget import (
    INTERCEPT_SPEED_M_S,
    TrackerHardware,
    photon_sigma_theta_rad,
    photon_snr,
)

# The beacon is "somewhat directionally pointed": coarse-pointed along the train axis toward the
# surveyor (a few degrees), no fine closed-loop tracking.  It sits between the wide no-pointing LED
# cone (20°) and a tightly-steered laser — required intensity ∝ beam solid angle, so a ~5° cone buys
# ~an order of magnitude of SNR over the 20° cone for only coarse (train-axis) attitude knowledge,
# which the swarm has for free (the surveyor is the lead unit along the velocity).
SOMEWHAT_DIRECTIONAL_HALFANGLE_RAD: float = math.radians(5.0)

# The committed 99 % per-unit capture criterion as a plate-to-σ ratio: the Tier-1 5 m plate ↔
# σ ≤ 1.65 m is a 2-D Rayleigh R/σ = 3.03, and the same ratio sizes any smaller plate (ADR 0015).
R_OVER_SIGMA: float = PLATE_RADIUS_M / CAPTURE_SIGMA_MAX_M

# The Tier-2 goal and stretch (CONTEXT.md *Surveyor-anchored centering*) and the Tier-1 baseline.
TARGET_PLATE_RADIUS_M: float = 0.10
STRETCH_PLATE_RADIUS_M: float = 0.05
COMMITTED_PLATE_RADIUS_M: float = PLATE_RADIUS_M

# Intra-train geometry: the committed CONTEXT.md formula is σ_δ = σ_θ·v/f, where the plate-crossing
# cadence f sets the nearest-follower spacing (and thus the link range) v/f.  The committed band is
# "cm at 2–4 Hz"; the conservative default is the *slowest* cadence (longest link, worst scatter).
STROBE_RATE_HZ: float = 2.0
DEFAULT_LINK_RANGE_M: float = INTERCEPT_SPEED_M_S / STROBE_RATE_HZ

# Gram-scale PuffSat camera: a ~5 mm aperture (the honest gram-scale optic — a small aperture makes
# a big diffraction spot that centroiding beats down by the SNR) and the same bench-calibratable
# focal-plane distortion floor as the terminal tracker (the term that actually sets the budget).
GRAM_SCALE_APERTURE_M: float = 5.0e-3
DEFAULT_DISTORTION_FLOOR_RAD: float = TrackerHardware().distortion_floor_rad

# Realistic block-shift metrology: a rendezvous/docking-lidar-class hoop reaches ~1 cm ranging.
DEFAULT_SIGMA_HOOP_M: float = 0.01

# Q-switched beacon (the design that makes the km-class link distortion-limited, not photon-ltd):
# bright ns pulses at known timings.  ~100 kW peak, 10 Hz rep, a few hundred mW *average* → a 30 mJ,
# 300 ns pulse.  The gated measurement collects the whole pulse energy, so the aperture sees
# ~10⁵–10⁸ photons even at km range through the coarse-pointed cone; the ns gate + narrowband line
# filter keep the read signal-dominated (the surveyor looks backward at dark sky — low background).
DEFAULT_PEAK_POWER_W: float = 100.0e3
DEFAULT_REP_RATE_HZ: float = 10.0
DEFAULT_AVG_POWER_W: float = 0.30
# Conservative shot-equivalent noise floor per gated read; a narrowband filter + ns gate hold it
# here (or below) even against a sunlit background — and the link is signal-dominated regardless.
NARROWBAND_NOISE_PHOTONS: float = 100.0

# Majority/averaging cameras per unit, and the pessimistic default that identical copies share their
# distortion bias fully (ρ=1) — so voting cannot cut it; the report shows diversity (ρ=0) as upside.
DEFAULT_N_CAMERAS: int = 3
IDENTICAL_CAMERA_CORRELATION: float = 1.0
DIVERSE_CAMERA_CORRELATION: float = 0.0


@dataclass(frozen=True)
class QSwitchedBeacon:
    """A Q-switched inter-PuffSat beacon: bright, short, known-timing pulses read in a matched gate.

    The centroiding measurement integrates one pulse, so the collected signal is the pulse *energy*
    (``avg_power / rep_rate``), not the average power — that is what makes a few hundred mW of
    average power act like a ~100 kW source for the duration of the gate.  ``noise_photons`` is the
    shot-equivalent noise floor left after the narrowband line filter and the ns time-gate.
    """

    peak_power_w: float = DEFAULT_PEAK_POWER_W
    rep_rate_hz: float = DEFAULT_REP_RATE_HZ
    avg_power_w: float = DEFAULT_AVG_POWER_W
    noise_photons: float = NARROWBAND_NOISE_PHOTONS

    @property
    def pulse_energy_j(self) -> float:
        """Energy in one pulse — what the gated centroiding measurement collects: ``avg / rep``."""
        return self.avg_power_w / self.rep_rate_hz

    @property
    def pulse_duration_s(self) -> float:
        """The gated exposure (pulse width): ``pulse_energy / peak_power``."""
        return self.pulse_energy_j / self.peak_power_w

    @property
    def duty_cycle(self) -> float:
        """Fraction of time the beacon is on: ``avg / peak`` (tiny — a strobed source)."""
        return self.avg_power_w / self.peak_power_w


@dataclass(frozen=True)
class CameraModel:
    """The PuffSat's camera-rigid bearing sensor: optics, Q-switched beacon, and N-camera voting.

    ``distortion_correlation`` is the detector-to-detector correlation ρ of the focal-plane
    distortion across the ``n_cameras``: ρ=1 (identical copies) means voting/averaging cannot reduce
    the shared bias (the array collapses to a single camera on that term), ρ=0 (physically diverse,
    separately calibrated) means it averages with √N.
    """

    aperture_m: float = GRAM_SCALE_APERTURE_M
    distortion_floor_rad: float = DEFAULT_DISTORTION_FLOOR_RAD
    beacon: QSwitchedBeacon = field(default_factory=QSwitchedBeacon)
    beam_half_angle_rad: float = SOMEWHAT_DIRECTIONAL_HALFANGLE_RAD
    n_cameras: int = DEFAULT_N_CAMERAS
    distortion_correlation: float = IDENTICAL_CAMERA_CORRELATION

    @property
    def hardware(self) -> TrackerHardware:
        """The :class:`TrackerHardware` this camera + Q-switched beacon presents to the photon fn.

        The Q-switched pulse maps onto the tracker-budget beacon as *peak power held for the gated
        pulse duration* — ``signal_photons`` is power × exposure = the pulse energy — with the
        narrowband/gated noise floor.
        """
        return replace(
            TrackerHardware(),
            aperture_m=self.aperture_m,
            beam_half_angle_rad=self.beam_half_angle_rad,
            distortion_floor_rad=self.distortion_floor_rad,
            beacon_power_w=self.beacon.peak_power_w,
            exposure_s=self.beacon.pulse_duration_s,
            noise_photons=self.beacon.noise_photons,
        )


def camera_sigma_theta_rad(camera: CameraModel, range_m: float) -> float:
    """The camera-rigid angular grade after N-camera voting: RSS(photon, distortion) reduced by ρ.

    The photon-limited centroiding term is independent per camera (averages with √N); the distortion
    floor is the ρ-correlated bias.  Reuses the equicorrelated array law from ``distortion_field``
    (there is no post-impact smear here — the surveyor looks backward at dark sky, and the ns gate
    freezes any residual body rate — so the common floor is zero).
    """
    photon = photon_sigma_theta_rad(camera.hardware, range_m)
    budget = DetectorBudget(
        sigma_distortion_rad=camera.distortion_floor_rad,
        sigma_indep_other_rad=photon,
        sigma_common_rad=0.0,
    )
    return array_sigma_theta_rad(budget, camera.n_cameras, camera.distortion_correlation)


def camera_photon_sigma_theta_rad(camera: CameraModel, range_m: float) -> float:
    """The photon-limited centroiding term alone (before RSS with the distortion floor)."""
    return photon_sigma_theta_rad(camera.hardware, range_m)


def camera_snr(camera: CameraModel, range_m: float) -> float:
    """The beacon centroiding SNR at ``range_m`` (per gated pulse)."""
    return photon_snr(camera.hardware, range_m)


def per_unit_scatter_m(camera: CameraModel, link_range_m: float) -> float:
    """The independent per-unit cross-track scatter: ``σ_θ · R_link`` (bearing knowledge limit)."""
    return camera_sigma_theta_rad(camera, link_range_m) * link_range_m


def arrival_sigma_m(sigma_hoop_m: float, scatter_m: float) -> float:
    """The per-unit arrival 1σ: RSS of the hoop-pinned block-shift bias and the per-unit scatter."""
    return math.hypot(sigma_hoop_m, scatter_m)


def plate_radius_for_sigma(arrival_sigma_m: float, r_over_sigma: float = R_OVER_SIGMA) -> float:
    """Plate radius for 99 % per-unit capture at a given arrival σ: ``(R/σ)·σ``."""
    return r_over_sigma * arrival_sigma_m


def sigma_for_plate_radius(plate_radius_m: float, r_over_sigma: float = R_OVER_SIGMA) -> float:
    """The arrival σ a plate radius allows at 99 % capture — the inverse map (see above)."""
    return plate_radius_m / r_over_sigma


@dataclass(frozen=True)
class CenteringSpec:
    """One declared surveyor-anchored centering point — the metrology the plate size sizes to."""

    sigma_hoop_m: float = DEFAULT_SIGMA_HOOP_M
    camera: CameraModel = field(default_factory=CameraModel)
    link_range_m: float = DEFAULT_LINK_RANGE_M
    target_plate_radius_m: float = TARGET_PLATE_RADIUS_M
    stretch_plate_radius_m: float = STRETCH_PLATE_RADIUS_M
    r_over_sigma: float = R_OVER_SIGMA

    @property
    def strobe_rate_hz(self) -> float:
        """The plate-crossing cadence implied by the link range: ``f = v / R_link``."""
        return INTERCEPT_SPEED_M_S / self.link_range_m

    @property
    def scatter_m(self) -> float:
        return per_unit_scatter_m(self.camera, self.link_range_m)

    @property
    def arrival_sigma_m(self) -> float:
        return arrival_sigma_m(self.sigma_hoop_m, self.scatter_m)

    @property
    def plate_radius_m(self) -> float:
        return plate_radius_for_sigma(self.arrival_sigma_m, self.r_over_sigma)

    @property
    def meets_target(self) -> bool:
        return self.plate_radius_m <= self.target_plate_radius_m

    @property
    def meets_stretch(self) -> bool:
        return self.plate_radius_m <= self.stretch_plate_radius_m

    @property
    def scatter_dominates(self) -> bool:
        """Whether per-unit camera scatter, not the hoop bias, is the larger of the two legs."""
        return self.scatter_m > self.sigma_hoop_m


@dataclass(frozen=True)
class PlatePoint:
    """A (label, hoop, scatter) → plate-radius reading, used by the hoop sweep and references."""

    label: str
    sigma_hoop_m: float
    scatter_m: float
    spec: CenteringSpec

    @property
    def arrival_sigma_m(self) -> float:
        return arrival_sigma_m(self.sigma_hoop_m, self.scatter_m)

    @property
    def plate_radius_m(self) -> float:
        return plate_radius_for_sigma(self.arrival_sigma_m, self.spec.r_over_sigma)

    @property
    def meets_target(self) -> bool:
        return self.plate_radius_m <= self.spec.target_plate_radius_m

    @property
    def meets_stretch(self) -> bool:
        return self.plate_radius_m <= self.spec.stretch_plate_radius_m


@dataclass(frozen=True)
class CameraRead:
    """The camera angular grade and scatter at one (N, ρ) — the voting question as a number."""

    label: str
    n_cameras: int
    distortion_correlation: float
    sigma_theta_rad: float
    scatter_m: float


@dataclass(frozen=True)
class ReferencePoint:
    """A named, realistic hardware class — the anchors that make ``"if we get X"`` concrete.

    ``source`` names the instrument *class* (not a verified part number): the realism claim is that
    the requirement lands inside an existing class, which a datasheet citation would then confirm.
    """

    label: str
    sigma_hoop_m: float
    distortion_floor_rad: float
    source: str


# Reference hardware classes.  σ_hoop values are rendezvous/ranging-instrument class; distortion
# floors are calibrated-star-tracker (3 µrad) vs an uncalibrated micro-optic (coarse).  These are
# class-level anchors to be confirmed against datasheets, not verified part specs.
DEFAULT_REFERENCES: tuple[ReferencePoint, ...] = (
    ReferencePoint(
        "docking lidar + calibrated micro-tracker", 0.01, 3.0e-6, "rendezvous-lidar class"
    ),
    ReferencePoint(
        "mm laser rangefinder + calibrated tracker", 0.003, 3.0e-6, "laser-tracker class"
    ),
    ReferencePoint(
        "microwave trilateration + calibrated tracker", 0.03, 3.0e-6, "RF-metrology class"
    ),
    ReferencePoint("docking lidar + uncalibrated micro-optic", 0.01, 15.0e-6, "coarse-optic class"),
)

# Hoop-precision sweep: sub-cm laser-tracker to few-cm RF trilateration.
DEFAULT_HOOP_SWEEP_M: tuple[float, ...] = (0.002, 0.005, 0.01, 0.02, 0.03, 0.05)


def hoop_sweep(spec: CenteringSpec, sigma_hoops_m: Sequence[float]) -> tuple[PlatePoint, ...]:
    """Plate radius vs hoop precision at the spec's camera scatter (the block-shift requirement)."""
    scatter = spec.scatter_m
    return tuple(PlatePoint(f"σ_hoop {s * 1e2:.1f} cm", s, scatter, spec) for s in sigma_hoops_m)


def reference_points(
    spec: CenteringSpec, references: Sequence[ReferencePoint]
) -> tuple[PlatePoint, ...]:
    """Evaluate each reference hardware class → its plate radius (its own hoop and distortion)."""
    points: list[PlatePoint] = []
    for ref in references:
        camera = replace(spec.camera, distortion_floor_rad=ref.distortion_floor_rad)
        scatter = per_unit_scatter_m(camera, spec.link_range_m)
        points.append(PlatePoint(ref.label, ref.sigma_hoop_m, scatter, spec))
    return tuple(points)


def camera_reads(spec: CenteringSpec) -> tuple[CameraRead, ...]:
    """The (N, ρ) camera grades that answer the 3-camera majority-voting question numerically."""
    configs = (
        ("1 camera", 1, IDENTICAL_CAMERA_CORRELATION),
        (
            f"{spec.camera.n_cameras} identical (ρ=1)",
            spec.camera.n_cameras,
            IDENTICAL_CAMERA_CORRELATION,
        ),
        (
            f"{spec.camera.n_cameras} diverse (ρ=0)",
            spec.camera.n_cameras,
            DIVERSE_CAMERA_CORRELATION,
        ),
    )
    reads: list[CameraRead] = []
    for label, n, rho in configs:
        camera = replace(spec.camera, n_cameras=n, distortion_correlation=rho)
        sigma_theta = camera_sigma_theta_rad(camera, spec.link_range_m)
        reads.append(CameraRead(label, n, rho, sigma_theta, sigma_theta * spec.link_range_m))
    return tuple(reads)


@dataclass(frozen=True)
class CenteringFinding:
    """The Tier-2 centering sizing: nominal plate, hoop requirement, references, camera voting."""

    spec: CenteringSpec
    hoop_points: tuple[PlatePoint, ...]
    reference_points: tuple[PlatePoint, ...]
    camera_reads: tuple[CameraRead, ...]

    @property
    def arrival_sigma_m(self) -> float:
        return self.spec.arrival_sigma_m

    @property
    def plate_radius_m(self) -> float:
        return self.spec.plate_radius_m

    @property
    def meets_target(self) -> bool:
        return self.spec.meets_target

    @property
    def meets_stretch(self) -> bool:
        return self.spec.meets_stretch

    @property
    def sigma_theta_rad(self) -> float:
        """The N-camera angular grade driving the per-unit scatter at the link range."""
        return camera_sigma_theta_rad(self.spec.camera, self.spec.link_range_m)

    @property
    def photon_sigma_theta_rad(self) -> float:
        """The single-camera photon-limited centroiding term at the link range (before RSS)."""
        return camera_photon_sigma_theta_rad(self.spec.camera, self.spec.link_range_m)

    @property
    def beacon_snr(self) -> float:
        """The Q-switched beacon centroiding SNR at the link range (per gated pulse)."""
        return camera_snr(self.spec.camera, self.spec.link_range_m)

    @property
    def photon_negligible(self) -> bool:
        """Whether the photon term sits under the distortion floor — the distortion-ltd claim."""
        return self.photon_sigma_theta_rad <= self.spec.camera.distortion_floor_rad

    @property
    def improvement_over_committed(self) -> float:
        """How many times smaller than the Tier-1 5 m plate the sized plate is."""
        return COMMITTED_PLATE_RADIUS_M / self.plate_radius_m

    @property
    def max_hoop_for_target_m(self) -> float:
        """The loosest hoop precision that still reaches the 10 cm target at the spec's scatter.

        ``plate = (R/σ)·√(σ_hoop² + scatter²) ≤ target`` inverts to
        ``σ_hoop ≤ √((target/(R/σ))² − scatter²)`` — negative under the root means the scatter alone
        already overruns the target (returns 0).
        """
        allowed = sigma_for_plate_radius(self.spec.target_plate_radius_m, self.spec.r_over_sigma)
        slack = allowed**2 - self.spec.scatter_m**2
        return math.sqrt(slack) if slack > 0.0 else 0.0

    @property
    def scatter_limits_target(self) -> bool:
        """The per-unit scatter alone exceeds the 10 cm target budget — no hoop can rescue it."""
        return self.max_hoop_for_target_m <= 0.0


def centering_finding(
    spec: CenteringSpec | None = None,
    *,
    sigma_hoops_m: Sequence[float] = DEFAULT_HOOP_SWEEP_M,
    references: Sequence[ReferencePoint] = DEFAULT_REFERENCES,
) -> CenteringFinding:
    """Assemble the Tier-2 centering sizing finding (the pure runner)."""
    s = spec if spec is not None else CenteringSpec()
    return CenteringFinding(
        spec=s,
        hoop_points=hoop_sweep(s, sigma_hoops_m),
        reference_points=reference_points(s, references),
        camera_reads=camera_reads(s),
    )


def _verdict(point: PlatePoint | CenteringFinding) -> str:
    if point.meets_stretch:
        return "≤5 cm"
    if point.meets_target:
        return "≤10 cm"
    return "over"


def format_centering_budget(finding: CenteringFinding) -> str:
    """One-screen Tier-2 centering sizing — nominal plate, hoop requirement, references, voting."""
    s = finding.spec
    cam = s.camera
    beacon = cam.beacon
    plate_cm = finding.plate_radius_m * 1e2
    scatter_cm = s.scatter_m * 1e2
    hoop_cm = s.sigma_hoop_m * 1e2

    if finding.meets_stretch:
        verdict = "reaches the 5 cm stretch"
    elif finding.meets_target:
        verdict = "reaches the 10 cm target"
    else:
        verdict = "over the 10 cm target"

    photon_note = "distortion-limited" if finding.photon_negligible else "PHOTON-LIMITED"

    lines = [
        "Surveyor-anchored centering budget — Tier-2 ~10 cm sizing (paper-side, no ADR)",
        "  Metrology-limited, not control: cm trims sit deep inside the ~475 m funnel, so this"
        " sizes hardware,",
        "  it does not resurrect MPC. Plate = 3.03·σ_arrival (99 % 2-D Rayleigh, ADR 0015).",
        f"  Link: strobe f = {s.strobe_rate_hz:.1f} Hz → σ_δ = σ_θ·v/f at R_link ="
        f" {s.link_range_m / 1e3:.1f} km ({cam.n_cameras} gram-scale cameras,"
        f" {cam.aperture_m * 1e3:.0f} mm, distortion floor"
        f" {cam.distortion_floor_rad * 1e6:.0f} µrad).",
        f"  Q-switched beacon: {beacon.peak_power_w / 1e3:.0f} kW peak @"
        f" {beacon.rep_rate_hz:.0f} Hz,"
        f" {beacon.avg_power_w * 1e3:.0f} mW avg → {beacon.pulse_energy_j * 1e3:.0f} mJ /"
        f" {beacon.pulse_duration_s * 1e9:.0f} ns pulse (duty {beacon.duty_cycle:.1e}),"
        f" ±{math.degrees(cam.beam_half_angle_rad):.0f}° coarse-pointed, narrowband + ns-gated."
        f" SNR {finding.beacon_snr:.0f} → photon σ_θ {finding.photon_sigma_theta_rad * 1e6:.2f}"
        f" µrad ({photon_note}).",
        f"  Legs: block-shift bias σ_hoop {hoop_cm:.2f} cm ⊕ per-unit scatter σ_θ·R"
        f" {scatter_cm:.2f} cm → arrival σ {finding.arrival_sigma_m * 1e2:.2f} cm →"
        f" plate {plate_cm:.1f} cm — {verdict}",
        f"  ({finding.improvement_over_committed:.0f}× smaller than the 5 m Tier-1 plate).",
        "",
        "  Hoop-precision requirement (plate vs σ_hoop at the nominal camera scatter):",
    ]
    for p in finding.hoop_points:
        lines.append(
            f"    σ_hoop {p.sigma_hoop_m * 1e2:5.1f} cm → plate {p.plate_radius_m * 1e2:5.1f} cm"
            f" [{_verdict(p)}]"
        )
    if finding.scatter_limits_target:
        lines.append(
            "    (per-unit scatter alone overruns the 10 cm budget — tighten the camera/link)."
        )
    else:
        lines.append(
            f"    → 10 cm needs σ_hoop ≤ {finding.max_hoop_for_target_m * 1e2:.1f} cm here."
        )

    lines.append("")
    lines.append("  Camera voting (does 3-camera majority voting kill the distortion bias?):")
    for r in finding.camera_reads:
        lines.append(
            f"    {r.label:20} → σ_θ {r.sigma_theta_rad * 1e6:5.2f} µrad,"
            f" scatter {r.scatter_m * 1e2:5.2f} cm"
        )
    lines.append(
        "    Identical copies (ρ=1) share the bias — voting only averages the random term;"
    )
    lines.append(
        "    physically diverse cameras (ρ=0) or differential astrometry cut the distortion floor."
    )

    lines.append("")
    lines.append("  Reference hardware classes (is 'if we get X, we hit N cm' realistic?):")
    for p in finding.reference_points:
        lines.append(f"    {p.label:42} → plate {p.plate_radius_m * 1e2:5.1f} cm [{_verdict(p)}]")
    lines.append(
        "    (class-level anchors — confirm against datasheets before citing a specific part)."
    )
    return "\n".join(lines)
