"""Pure σ_θ tracker-budget gate — the D1 blocking pre-gate (ADR 0018, no JVM).

ADR 0015 made the terminal navigation a *target-side astrometric tracker* and derived the
load-bearing hardware requirement **σ_θ ≤ 10 µrad** (5 µrad design target): the homing floor
is ``σ_miss ≈ 2σ_θ²v²/a_max`` (1.46 m at 10 µrad), so the whole catch-radius / plate-capture
story rides on that angular grade.  ADR 0018 makes *achievability* of 10 µrad a **blocking**
gate before the Rung-D Monte Carlo — if a µrad-class tracker on a 1 Hz-hammered bus cannot
reach it with bench-realistic hardware, the verdict falls and there is nothing to simulate.

This module is that gate: a pure optical error budget (no Orekit — angular precision is a
focal-plane question, not an orbit-propagation one).  It does two things.

* **Error budget → achieved σ_θ.**  An RSS of named, individually-traceable sub-terms for a
  declared hardware point: photon-limited centroiding (the active laser beacon makes the
  "dim, fast target" *bright* — the SNR is enormous, so this term is negligible and the
  budget is **not** photon-limited), the post-impact smear residual left after the
  differential measurement cancels the common bus motion, the gyro inter-frame bridge, and
  the bench-calibratable focal-plane-distortion floor (the actual limit).  The achieved σ_θ
  feeds the C3b homing floor (:func:`puffsat_sim.guidance.homing_floor_m`) and the ADR 0015
  capture criterion.

* **Acquisition.**  The tracker FOV must contain the PuffSat at first lock given the C1
  hand-off delivery dispersion (FOV ≥ n·σ_lateral / R), the beacon beam must cover that same
  cone, and the FOV must hold enough reference stars for the differential solution.  The
  binding FOV is the larger of the acquisition and reference-star requirements, and a
  commodity detector must resolve the diffraction spot to Nyquist across it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from puffsat_sim.anti_drag import PEAK_THRUST_LIMIT_N
from puffsat_sim.constants import PLANCK_J_S, SPEED_OF_LIGHT_M_S
from puffsat_sim.guidance import CAPTURE_SIGMA_MAX_M, homing_floor_m

# The ADR 0015 requirement and design target for the terminal tracker angular grade.
REQUIRED_SIGMA_THETA_RAD: float = 10e-6
TARGET_SIGMA_THETA_RAD: float = 5e-6

# Terminal closing geometry (reference orbit / ADR 0014/0015): interception speed and the
# 400 mN / 25 kg actuator that set the homing floor 2σ_θ²v²/a_max.
INTERCEPT_SPEED_M_S: float = 10_780.0
TERMINAL_MASS_KG: float = 25.0
A_MAX_M_S2: float = PEAK_THRUST_LIMIT_N / TERMINAL_MASS_KG

# C1 measured midcourse delivery dispersion at the 800 km hand-off (ADR 0012/0015): the
# lateral 1σ the terminal tracker must acquire within, and the range by which angle-based
# knowledge must be actionable at the metre scale (ADR 0015).
HANDOFF_SIGMA_LATERAL_M: float = 141.0
ACQUISITION_RANGE_M: float = 300_000.0
ACQUISITION_SIGMAS: float = 3.0

# Reference-star density brighter than ~10th magnitude (≈3.5e5 stars over 4π sr): the
# differential astrometry needs a handful of these in the FOV alongside the beacon.
STAR_DENSITY_PER_SR: float = 2.8e4
MIN_REFERENCE_STARS: int = 3

# Sample the diffraction spot at Nyquist (~2 px across) so the centroid is well-posed.
NYQUIST_SAMPLES: float = 2.0


def photon_energy_j(wavelength_m: float) -> float:
    """Energy of one photon at ``wavelength_m``: ``h·c/λ``."""
    return PLANCK_J_S * SPEED_OF_LIGHT_M_S / wavelength_m


def beacon_intensity_w_sr(power_w: float, beam_half_angle_rad: float) -> float:
    """Radiant intensity toward the PuffSat: power over the beam-cone ``π·θ²`` solid angle."""
    return power_w / (math.pi * beam_half_angle_rad**2)


def diffraction_spot_rad(aperture_m: float, wavelength_m: float) -> float:
    """The diffraction spot scale ``λ/D`` — the angular size an isolated point source images to."""
    return wavelength_m / aperture_m


@dataclass(frozen=True)
class TrackerHardware:
    """One declared terminal-tracker hardware point — the inputs the budget is computed for.

    Conservative defaults: a 5 cm aperture, 1 ms exposures on a 1 W laser beacon at 1064 nm
    (ADR 0015), a beam wide enough to cover the acquisition cone, a residual post-impact
    body rate *after* the differential measurement cancels the common bus motion, a
    nav-grade gyro for the inter-frame bridge, and a bench-calibratable focal-plane
    distortion floor (the term that actually sets the budget).
    """

    aperture_m: float = 0.05
    exposure_s: float = 1.0e-3
    efficiency: float = 0.3
    wavelength_m: float = 1.064e-6
    beacon_power_w: float = 1.0
    beam_half_angle_rad: float = 2.0e-3
    noise_photons: float = 100.0
    residual_body_rate_rad_s: float = 3.0e-3
    gyro_arw_rad_sqrt_s: float = 5.8e-7
    frame_interval_s: float = 1.0
    distortion_floor_rad: float = 3.0e-6


def signal_photons(hw: TrackerHardware, range_m: float) -> float:
    """Detected beacon photons in one exposure: collected power × QE × t_exp / photon energy."""
    intensity = beacon_intensity_w_sr(hw.beacon_power_w, hw.beam_half_angle_rad)
    aperture_area = math.pi * (hw.aperture_m / 2.0) ** 2
    received_w = intensity * aperture_area / range_m**2
    return received_w * hw.efficiency * hw.exposure_s / photon_energy_j(hw.wavelength_m)


def photon_snr(hw: TrackerHardware, range_m: float) -> float:
    """Centroiding SNR of the beacon: ``N_sig / √(N_sig + noise)`` (shot-noise-dominated)."""
    n_sig = signal_photons(hw, range_m)
    return n_sig / math.sqrt(n_sig + hw.noise_photons)


def photons_for_snr(snr: float, noise_photons: float) -> float:
    """Invert ``SNR = N/√(N+noise)`` for the signal photon count ``N`` (the positive root)."""
    s2 = snr * snr
    return (s2 + math.sqrt(s2 * s2 + 4.0 * s2 * noise_photons)) / 2.0


def photon_sigma_theta_rad(hw: TrackerHardware, range_m: float) -> float:
    """Photon-limited centroiding precision of an isolated point source: ``(λ/D) / SNR``.

    An isolated bright source can be centroided well below the diffraction limit — the
    Airy scale ``λ/D`` divided by the SNR — so a bright laser beacon drives this term far
    under the budget (the "dim, fast target" worry is dissolved by making it active).
    """
    spot = diffraction_spot_rad(hw.aperture_m, hw.wavelength_m)
    return spot / photon_snr(hw, range_m)


def power_for_centroid_floor(
    hw: TrackerHardware, range_m: float, target_sigma_theta_rad: float
) -> float:
    """Minimum beacon (peak) power so the photon-limited centroid term ≤ ``target_sigma_theta_rad``.

    The inversion of :func:`photon_sigma_theta_rad`: the target σ_θ fixes the SNR
    (``(λ/D)/target``), :func:`photons_for_snr` turns that into a photon count, and
    :func:`signal_photons` is linear in power, so the power is the count divided by the
    photons-per-watt at this range and aperture.  This is the *peak* power — the power held
    during one exposure; cadence (duty cycle) only scales the thermal average, not this.
    """
    spot = diffraction_spot_rad(hw.aperture_m, hw.wavelength_m)
    required_snr = spot / target_sigma_theta_rad
    n_req = photons_for_snr(required_snr, hw.noise_photons)
    photons_per_watt = signal_photons(replace(hw, beacon_power_w=1.0), range_m)
    return n_req / photons_per_watt


def smear_sigma_theta_rad(residual_body_rate_rad_s: float, exposure_s: float) -> float:
    """Centroid blur left after differential cancellation: the residual-rate streak / √12.

    Differential astrometry (beacon vs the same-frame star field) cancels the *common* bus
    motion; what survives is the residual body rate × exposure, treated as a uniform streak
    whose centroid spread is ``length/√12`` (a uniform-blur 1σ, the conservative reading).
    """
    return residual_body_rate_rad_s * exposure_s / math.sqrt(12.0)


def gyro_bridge_sigma_theta_rad(gyro_arw_rad_sqrt_s: float, frame_interval_s: float) -> float:
    """Attitude knowledge accrued bridging between the ~1 Hz absolute frames: ``ARW·√Δt``."""
    return gyro_arw_rad_sqrt_s * math.sqrt(frame_interval_s)


@dataclass(frozen=True)
class BudgetTerm:
    """One angular-error contribution in the σ_θ budget (RSS'd as independent sources)."""

    label: str
    sigma_theta_rad: float


def error_budget(hw: TrackerHardware, range_m: float) -> tuple[BudgetTerm, ...]:
    """The four σ_θ contributions for a hardware point at the (worst, longest) design range."""
    return (
        BudgetTerm(
            "photon-limited centroiding (active beacon)",
            photon_sigma_theta_rad(hw, range_m),
        ),
        BudgetTerm(
            "post-impact smear residual (differential)",
            smear_sigma_theta_rad(hw.residual_body_rate_rad_s, hw.exposure_s),
        ),
        BudgetTerm(
            "gyro inter-frame bridge",
            gyro_bridge_sigma_theta_rad(hw.gyro_arw_rad_sqrt_s, hw.frame_interval_s),
        ),
        BudgetTerm("focal-plane distortion floor (calibrated)", hw.distortion_floor_rad),
    )


def rss_sigma_theta_rad(terms: tuple[BudgetTerm, ...]) -> float:
    """Root-sum-square of the budget terms — the achieved angular grade."""
    return math.sqrt(sum(t.sigma_theta_rad**2 for t in terms))


def required_fov_halfangle_rad(
    handoff_sigma_lateral_m: float, acquisition_range_m: float, n_sigma: float
) -> float:
    """The FOV half-angle that contains the PuffSat at first lock: ``n·σ_lateral / R``."""
    return n_sigma * handoff_sigma_lateral_m / acquisition_range_m


def reference_star_fov_halfangle_rad(n_stars: int, star_density_per_sr: float) -> float:
    """The FOV half-angle whose cone holds ``n_stars`` reference stars at the given density.

    ``n_stars / density = Ω ≈ π·θ²`` for a small cone, so ``θ = √(n / (π·density))``.
    """
    return math.sqrt(n_stars / (math.pi * star_density_per_sr))


def detector_pixels_across(fov_halfangle_rad: float, pixel_ifov_rad: float) -> int:
    """Pixels across the full FOV (``2·θ_half / IFOV``), rounded up."""
    return math.ceil(2.0 * fov_halfangle_rad / pixel_ifov_rad)


@dataclass(frozen=True)
class TrackerBudgetFinding:
    """The σ_θ gate verdict: the achieved grade, its capture consequence, and acquisition."""

    hardware: TrackerHardware
    design_range_m: float
    terms: tuple[BudgetTerm, ...]
    handoff_sigma_lateral_m: float
    acquisition_range_m: float
    acquisition_sigmas: float
    star_density_per_sr: float
    min_reference_stars: int
    required_sigma_theta_rad: float = REQUIRED_SIGMA_THETA_RAD
    target_sigma_theta_rad: float = TARGET_SIGMA_THETA_RAD
    speed_m_s: float = INTERCEPT_SPEED_M_S
    a_max_m_s2: float = A_MAX_M_S2
    capture_sigma_max_m: float = CAPTURE_SIGMA_MAX_M

    @property
    def achieved_sigma_theta_rad(self) -> float:
        return rss_sigma_theta_rad(self.terms)

    @property
    def driving_term(self) -> BudgetTerm:
        return max(self.terms, key=lambda t: t.sigma_theta_rad)

    @property
    def meets_requirement(self) -> bool:
        return self.achieved_sigma_theta_rad <= self.required_sigma_theta_rad

    @property
    def meets_target(self) -> bool:
        return self.achieved_sigma_theta_rad <= self.target_sigma_theta_rad

    @property
    def requirement_margin(self) -> float:
        """How many times under the 10 µrad requirement the achieved grade sits."""
        return self.required_sigma_theta_rad / self.achieved_sigma_theta_rad

    @property
    def homing_floor_m(self) -> float:
        """The C3b homing floor at the *achieved* grade (the realized terminal nav limit)."""
        return homing_floor_m(self.achieved_sigma_theta_rad, self.speed_m_s, self.a_max_m_s2)

    @property
    def required_grade_floor_m(self) -> float:
        """The homing floor at the bare 10 µrad requirement (ADR 0015's thin-margin reference)."""
        return homing_floor_m(self.required_sigma_theta_rad, self.speed_m_s, self.a_max_m_s2)

    @property
    def capture_floor_met(self) -> bool:
        return self.homing_floor_m <= self.capture_sigma_max_m

    @property
    def acquisition_fov_halfangle_rad(self) -> float:
        return required_fov_halfangle_rad(
            self.handoff_sigma_lateral_m, self.acquisition_range_m, self.acquisition_sigmas
        )

    @property
    def reference_star_fov_halfangle_rad(self) -> float:
        return reference_star_fov_halfangle_rad(self.min_reference_stars, self.star_density_per_sr)

    @property
    def fov_halfangle_rad(self) -> float:
        """The binding FOV — the larger of the acquisition and reference-star requirements."""
        return max(self.acquisition_fov_halfangle_rad, self.reference_star_fov_halfangle_rad)

    @property
    def beam_covers_acquisition(self) -> bool:
        """The beacon beam must illuminate the whole acquisition cone."""
        return self.hardware.beam_half_angle_rad >= self.acquisition_fov_halfangle_rad

    @property
    def pixel_ifov_rad(self) -> float:
        """Pixel scale that samples the diffraction spot at Nyquist."""
        spot = diffraction_spot_rad(self.hardware.aperture_m, self.hardware.wavelength_m)
        return spot / NYQUIST_SAMPLES

    @property
    def detector_pixels_across(self) -> int:
        return detector_pixels_across(self.fov_halfangle_rad, self.pixel_ifov_rad)


def tracker_budget_finding(
    hardware: TrackerHardware | None = None,
    *,
    design_range_m: float = ACQUISITION_RANGE_M,
    handoff_sigma_lateral_m: float = HANDOFF_SIGMA_LATERAL_M,
    acquisition_range_m: float = ACQUISITION_RANGE_M,
    acquisition_sigmas: float = ACQUISITION_SIGMAS,
    star_density_per_sr: float = STAR_DENSITY_PER_SR,
    min_reference_stars: int = MIN_REFERENCE_STARS,
) -> TrackerBudgetFinding:
    """Assemble the σ_θ gate finding for a hardware point (the pure runner).

    ``hardware=None`` uses the conservative default point.  The photon term is evaluated at
    ``design_range_m`` — the longest (worst) range the grade must hold, the acquisition range
    — since the SNR only improves as the PuffSat closes.
    """
    hw = hardware if hardware is not None else TrackerHardware()
    return TrackerBudgetFinding(
        hardware=hw,
        design_range_m=design_range_m,
        terms=error_budget(hw, design_range_m),
        handoff_sigma_lateral_m=handoff_sigma_lateral_m,
        acquisition_range_m=acquisition_range_m,
        acquisition_sigmas=acquisition_sigmas,
        star_density_per_sr=star_density_per_sr,
        min_reference_stars=min_reference_stars,
    )


def format_tracker_budget(finding: TrackerBudgetFinding) -> str:
    """One-screen σ_θ gate report — hardware point, the budget, the verdict, acquisition."""
    hw = finding.hardware
    achieved_urad = finding.achieved_sigma_theta_rad * 1e6
    req_urad = finding.required_sigma_theta_rad * 1e6
    target_urad = finding.target_sigma_theta_rad * 1e6

    if finding.meets_target:
        verdict = f"PASS — meets the {target_urad:.0f} µrad design target"
    elif finding.meets_requirement:
        verdict = f"PASS — meets the {req_urad:.0f} µrad requirement (5 µrad target not reached)"
    else:
        verdict = f"FAIL — over the {req_urad:.0f} µrad requirement (blocks Rung D)"

    snr = photon_snr(hw, finding.design_range_m)
    lines = [
        "σ_θ tracker budget — D1 blocking pre-gate (ADR 0018; requirement ADR 0015)",
        f"  Hardware point: aperture {hw.aperture_m * 1e2:.0f} cm,"
        f" exposure {hw.exposure_s * 1e3:.0f} ms, beacon {hw.beacon_power_w:g} W"
        f" @ {hw.wavelength_m * 1e9:.0f} nm, beam ±{hw.beam_half_angle_rad * 1e3:g} mrad,"
        f" η {hw.efficiency:g}",
        f"  Beacon SNR at {finding.design_range_m / 1e3:.0f} km: {snr:.0f}"
        " — active beacon, so the budget is calibration/jitter-limited, not photon-limited.",
        "  σ_θ error budget (1σ, RSS of independent sources):",
    ]
    lines += [f"    {t.label}: {t.sigma_theta_rad * 1e6:.3f} µrad" for t in finding.terms]
    lines.append(
        f"    RSS {achieved_urad:.2f} µrad (driven by {finding.driving_term.label}) → {verdict}"
    )
    lines.append(
        f"  Margin: {finding.requirement_margin:.1f}× under the {req_urad:.0f} µrad requirement;"
        f" {target_urad:.0f} µrad target {'met' if finding.meets_target else 'not met'}."
    )

    capture_verdict = "MEETS" if finding.capture_floor_met else "FAILS"
    lines.append(
        f"  Homing floor 2σ_θ²v²/a_max at the achieved grade: {finding.homing_floor_m:.3f} m"
        f" vs σ ≤ {finding.capture_sigma_max_m:g} m — {capture_verdict}"
        f" (bare {req_urad:.0f} µrad requirement → {finding.required_grade_floor_m:.2f} m,"
        " ADR 0015's thin-margin reference)."
    )

    acq_urad = finding.acquisition_fov_halfangle_rad * 1e6
    star_urad = finding.reference_star_fov_halfangle_rad * 1e6
    beam_verdict = "covered" if finding.beam_covers_acquisition else "NOT covered by the beam"
    lines.append("  Acquisition:")
    lines.append(
        f"    FOV ≥ {finding.acquisition_sigmas:g}σ·{finding.handoff_sigma_lateral_m:.0f} m"
        f" / {finding.acquisition_range_m / 1e3:.0f} km = ±{acq_urad / 1e3:.2f} mrad"
        f" (acquisition cone — {beam_verdict});"
    )
    lines.append(
        f"    FOV ≥ ±{star_urad / 1e3:.2f} mrad for {finding.min_reference_stars} reference stars"
        f" → binding FOV ±{finding.fov_halfangle_rad * 1e3:.2f} mrad, resolved to Nyquist by a"
        f" {finding.detector_pixels_across}-px detector"
        f" ({finding.pixel_ifov_rad * 1e6:.1f} µrad/px)."
    )
    return "\n".join(lines)


# --- Peak-power sizing sweep -------------------------------------------------------------
#
# The σ_θ budget above is *not* photon-limited, so beacon power has slack: the only job of the
# photons is to push the centroiding term well under the calibration floor.  This sweep sizes
# the *minimum* peak power that does that for realistic ranges and apertures, so the LED beacon
# is not over-specified.  "Peak" because the beacon strobes — cadence sets the thermal average,
# not this number; and exposure is fixed (the smear term already shows 1 ms is motion-safe), so
# the binding axes are range (R², the photon falls off) and aperture (1/D⁴: a bigger dish both
# collects more and shrinks the spot).  The default beam is a wide LED cone (no fine pointing) —
# required power ∝ the beam solid angle, so a steered narrow beam buys the (θ/θ_laser)² factor
# back at the cost of pointing (the LED-vs-laser trade).

# Photon term held this far under the distortion floor counts as "good centroiding" — past it
# the photon contribution is negligible in the RSS and σ_θ is purely calibration-limited.
PEAK_POWER_MARGIN_BELOW_FLOOR: float = 0.1

# A lensed LED beacon wide enough to need no fine pointing (~20° half-cone), and the steered
# laser beam the original budget assumed — the reference for the narrow-beam power saving.
LED_BEAM_HALFANGLE_RAD: float = 0.35
LASER_REFERENCE_BEAM_HALFANGLE_RAD: float = 2.0e-3

# Representative link ranges, close swarm to far cooperative target: surveyor↔PuffSat (100 m,
# 1 km), the acquisition range (300 km), the co-flyer (500 km) and the target (2603 km).
DEFAULT_SWEEP_RANGES_M: tuple[float, ...] = (1.0e2, 1.0e3, 3.0e5, 5.0e5, 2.603e6)
DEFAULT_SWEEP_APERTURES_M: tuple[float, ...] = (0.02, 0.05, 0.10)

# The peak power the sizing is read against (the figure the design was questioning).
REFERENCE_PEAK_POWER_W: float = 1.0e3


@dataclass(frozen=True)
class PeakPowerCell:
    """The minimum beacon peak power for good centroiding at one (range, aperture) point."""

    range_m: float
    aperture_m: float
    peak_power_w: float
    pulse_energy_j: float
    photon_sigma_theta_rad: float
    snr: float


@dataclass(frozen=True)
class PeakPowerSpec:
    """The peak-power sweep grid: ranges × apertures at a fixed beam, exposure and target."""

    ranges_m: tuple[float, ...] = DEFAULT_SWEEP_RANGES_M
    apertures_m: tuple[float, ...] = DEFAULT_SWEEP_APERTURES_M
    beam_half_angle_rad: float = LED_BEAM_HALFANGLE_RAD
    margin_below_floor: float = PEAK_POWER_MARGIN_BELOW_FLOOR
    base_hardware: TrackerHardware = field(default_factory=TrackerHardware)


@dataclass(frozen=True)
class PeakPowerSweepFinding:
    """The sized peak-power grid and the reads off it (narrow-beam saving, 1 kW reference)."""

    spec: PeakPowerSpec
    cells: tuple[PeakPowerCell, ...]
    target_sigma_theta_rad: float
    reference_peak_power_w: float = REFERENCE_PEAK_POWER_W

    @property
    def narrow_beam_factor(self) -> float:
        """Power a steered ±2 mrad laser saves over the wide LED cone: ``(θ_beam/θ_laser)²``."""
        return (self.spec.beam_half_angle_rad / LASER_REFERENCE_BEAM_HALFANGLE_RAD) ** 2

    @property
    def max_peak_power_w(self) -> float:
        return max(c.peak_power_w for c in self.cells)

    @property
    def min_peak_power_w(self) -> float:
        return min(c.peak_power_w for c in self.cells)


def peak_power_sweep(spec: PeakPowerSpec | None = None) -> PeakPowerSweepFinding:
    """Size the minimum beacon peak power for good centroiding over the range × aperture grid."""
    s = spec if spec is not None else PeakPowerSpec()
    target = s.margin_below_floor * s.base_hardware.distortion_floor_rad
    cells: list[PeakPowerCell] = []
    for range_m in s.ranges_m:
        for aperture_m in s.apertures_m:
            hw = replace(
                s.base_hardware,
                aperture_m=aperture_m,
                beam_half_angle_rad=s.beam_half_angle_rad,
            )
            power_w = power_for_centroid_floor(hw, range_m, target)
            powered = replace(hw, beacon_power_w=power_w)
            cells.append(
                PeakPowerCell(
                    range_m=range_m,
                    aperture_m=aperture_m,
                    peak_power_w=power_w,
                    pulse_energy_j=power_w * hw.exposure_s,
                    photon_sigma_theta_rad=photon_sigma_theta_rad(powered, range_m),
                    snr=photon_snr(powered, range_m),
                )
            )
    return PeakPowerSweepFinding(spec=s, cells=tuple(cells), target_sigma_theta_rad=target)


def _format_power(power_w: float) -> str:
    if power_w >= 1.0e3:
        return f"{power_w / 1e3:.2f} kW"
    if power_w >= 1.0:
        return f"{power_w:.2f} W"
    if power_w >= 1.0e-3:
        return f"{power_w * 1e3:.2f} mW"
    return f"{power_w * 1e6:.2f} µW"


def _format_range(range_m: float) -> str:
    return f"{range_m / 1e3:g} km" if range_m >= 1.0e3 else f"{range_m:g} m"


def format_peak_power_sweep(finding: PeakPowerSweepFinding) -> str:
    """One-screen peak-power sizing — the minimum power for good centroiding by range × aperture."""
    s = finding.spec
    hw = s.base_hardware
    target_urad = finding.target_sigma_theta_rad * 1e6
    floor_urad = hw.distortion_floor_rad * 1e6

    lines = [
        "LED beacon peak-power sizing — minimum peak power for good centroiding (ADR 0018)",
        f"  Target: photon-limited centroid term ≤ {target_urad:.2f} µrad"
        f" ({s.margin_below_floor:g}× the {floor_urad:.0f} µrad distortion floor) — below this the"
        " photons are negligible and σ_θ is calibration-limited, not photon-limited.",
        f"  Beacon: LED cone ±{math.degrees(s.beam_half_angle_rad):.0f}° (wide, no fine pointing),"
        f" {hw.exposure_s * 1e3:.0f} ms exposure @ {hw.wavelength_m * 1e9:.0f} nm,"
        f" η {hw.efficiency:g}.",
        "  Minimum peak power, range (rows) × aperture (cols):",
        "    range \\ aperture  " + "  ".join(f"{a * 1e2:>7.0f} cm" for a in s.apertures_m),
    ]
    for range_m in s.ranges_m:
        row = [c for c in finding.cells if c.range_m == range_m]
        powers = "  ".join(f"{_format_power(c.peak_power_w):>10}" for c in row)
        lines.append(f"    {_format_range(range_m):>16}  {powers}")

    lines.append(
        "  Required peak power ∝ beam solid angle × range² ÷ aperture⁴."
        f"  A steered ±{LASER_REFERENCE_BEAM_HALFANGLE_RAD * 1e3:g} mrad laser would need"
        f" ~{finding.narrow_beam_factor:,.0f}× less — but at the cost of pointing (LED-vs-laser)."
    )
    lines.append(
        f"  vs the {finding.reference_peak_power_w / 1e3:g} kW reference: the close-swarm beacon"
        f" links ({_format_power(finding.min_peak_power_w)}) sit orders of magnitude under it;"
        " only a wide-LED long-range link approaches it (there, narrow the beam or grow the dish)."
    )
    return "\n".join(lines)
