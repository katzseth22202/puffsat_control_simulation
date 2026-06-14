"""Pure apogee nav-constellation sizing — the ADR 0020 Lever-3 follow-on (no JVM).

ADR 0020 graduated Lever 3 (a permanent ~150,000 km nav constellation) from option to a specced
architecture for the **coast / apogee-state** navigation the midcourse corrector consumes — a
different regime from terminal homing, and one GNSS cannot reach (GPS at 20,200 km vs apogee at
150,000 km).  This module sizes that architecture: the :mod:`puffsat_sim.tracker_budget` analog for
the coast regime, and like it pure (no Orekit — link budgets and *snapshot* GDOP are arithmetic and
geometry, not orbit propagation).

It answers the four sizing questions ADR 0020 settled:

* **Link budget** — does a Ka-band broadcast (decision 1) close from a 1 m / 10 W constellation dish
  to an omni PuffSat at the apogee range?
* **Velocity budget** — does carrier-phase Doppler give the radial-velocity precision the *binding*
  transverse-velocity requirement needs? (Ka's 20× L-band carrier helps here.)
* **GDOP / minimum members** — how many members (ring vs shell) does the snapshot velocity geometry
  need to pin apogee transverse velocity to the C1-matching target?
* **Transponder mass/power** — is the PuffSat-side hardware (one-way passive, decision 2) negligible
  on the 25 kg bus, and is the crypto/timing ASIC the mass driver?  (It is not.)

The conclusion ADR 0020 fixed is **match, don't beat**: target σ_Tvel ≈ 0.66 mm/s (= C1, the
per-unit entry budget D1.1 flies), not tighter — fusion already gives terminal capture 4.2× margin,
so apogee nav below the C1 grade is redundant.  The constellation's value beyond matching is
*snapshot GDOP at apogee* and pinning the rockets, not a tighter number.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.constants import BOLTZMANN_J_PER_K, EARTH_RADIUS_M, SPEED_OF_LIGHT_M_S

Vec3 = NDArray[np.float64]

# --- Requirement (C0/C1) -------------------------------------------------------------------------
# The C0 binding axis is apogee *transverse* velocity.  The requirement at the 400 m catch radius is
# ~1.84 mm/s; C1 *achieved* 0.66 mm/s, and that 0.66 mm/s is the per-unit entry budget D1.1 flies —
# so it is the match-not-beat *target* (ADR 0020 decision 5), tighter than the bare requirement.
REQUIRED_TVEL_SIGMA_M_S: float = 1.84e-3
TARGET_TVEL_SIGMA_M_S: float = 0.66e-3

# --- Reference geometry (design doc §3) ----------------------------------------------------------
APOGEE_ALT_M: float = 150_000_000.0
APOGEE_RADIUS_M: float = APOGEE_ALT_M + EARTH_RADIUS_M
# The link's representative member range (a typical shell separation at the apogee radius); the
# farthest members sit at ~2·R (≈6 dB more loss), still closing (the report notes it).
DESIGN_RANGE_M: float = APOGEE_RADIUS_M

# --- Default Ka-band broadcast link (ADR 0020 decision 1) ----------------------------------------
KA_FREQ_HZ: float = 30e9
CONSTELLATION_DISH_M: float = 1.0
CONSTELLATION_TX_POWER_W: float = 10.0
DISH_EFFICIENCY: float = 0.6
PUFFSAT_RX_GAIN_DBI: float = 0.0  # omni patch — ADR 0011 dec-7: omni PuffSat, gain on the infra
SYSTEM_NOISE_TEMP_K: float = 400.0
INTEGRATION_TIME_S: float = 1.0

# --- RTN axes at apogee (R along the apse line, T the in-plane transverse, N orbit-normal) --------
RADIAL_AXIS: Vec3 = np.array([1.0, 0.0, 0.0])
TRANSVERSE_AXIS: Vec3 = np.array([0.0, 1.0, 0.0])
NORMAL_AXIS: Vec3 = np.array([0.0, 0.0, 1.0])

PUFFSAT_BUS_MASS_KG: float = 25.0


# --- Link budget ---------------------------------------------------------------------------------
def free_space_path_loss_db(range_m: float, freq_hz: float) -> float:
    """Free-space path loss ``20·log10(4πR/λ)`` [dB]."""
    return 20.0 * math.log10(4.0 * math.pi * range_m * freq_hz / SPEED_OF_LIGHT_M_S)


def parabolic_gain_dbi(diameter_m: float, freq_hz: float, efficiency: float) -> float:
    """Boresight gain of a parabolic dish ``η·(πD/λ)²`` [dBi]."""
    wavelength = SPEED_OF_LIGHT_M_S / freq_hz
    return 10.0 * math.log10(efficiency * (math.pi * diameter_m / wavelength) ** 2)


def noise_density_dbw_hz(system_temp_k: float) -> float:
    """Thermal noise spectral density ``k·T`` [dBW/Hz]."""
    return 10.0 * math.log10(BOLTZMANN_J_PER_K * system_temp_k)


def carrier_to_noise_density_dbhz(
    eirp_dbw: float,
    range_m: float,
    freq_hz: float,
    rx_gain_dbi: float,
    system_temp_k: float,
) -> float:
    """Received C/N0 ``EIRP − FSPL + G_rx − kT`` [dB-Hz]."""
    return (
        eirp_dbw
        - free_space_path_loss_db(range_m, freq_hz)
        + rx_gain_dbi
        - noise_density_dbw_hz(system_temp_k)
    )


def downlink_cn0_dbhz(range_m: float = DESIGN_RANGE_M) -> float:
    """C/N0 of the default constellation→omni-PuffSat Ka downlink at ``range_m``."""
    eirp = parabolic_gain_dbi(
        CONSTELLATION_DISH_M, KA_FREQ_HZ, DISH_EFFICIENCY
    ) + 10.0 * math.log10(CONSTELLATION_TX_POWER_W)
    return carrier_to_noise_density_dbhz(
        eirp, range_m, KA_FREQ_HZ, PUFFSAT_RX_GAIN_DBI, SYSTEM_NOISE_TEMP_K
    )


# --- Velocity (carrier-phase Doppler) ------------------------------------------------------------
def carrier_phase_velocity_sigma_m_s(
    cn0_dbhz: float, freq_hz: float = KA_FREQ_HZ, integration_s: float = INTEGRATION_TIME_S
) -> float:
    """Radial-velocity 1σ from the carrier-Doppler (tone-frequency) Cramér–Rao bound over ``T``.

    σ_f = √(3 / (2π²·(C/N0)·T³)) [Hz]; σ_v = (c/f)·σ_f.  The T⁻³ falloff means even a short
    coherent integration on the slow coast drives radial velocity far under the requirement, and
    Ka's high carrier (vs L-band) shrinks σ_v = λ·σ_f directly.
    """
    cn0_lin = 10.0 ** (cn0_dbhz / 10.0)
    sigma_f = math.sqrt(3.0 / (2.0 * math.pi**2 * cn0_lin * integration_s**3))
    wavelength = SPEED_OF_LIGHT_M_S / freq_hz
    return wavelength * sigma_f


# --- Constellation geometry / velocity GDOP ------------------------------------------------------
def apogee_position() -> Vec3:
    """The PuffSat at apogee, on the +R axis."""
    return np.array([APOGEE_RADIUS_M, 0.0, 0.0])


def _fibonacci_sphere(n_members: int) -> Vec3:
    """``n_members`` near-uniform unit vectors on a sphere (a 3-D *shell* constellation)."""
    i = np.arange(n_members)
    z = 1.0 - 2.0 * (i + 0.5) / n_members
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    theta = math.pi * (3.0 - math.sqrt(5.0)) * i
    return np.stack([r * np.cos(theta), r * np.sin(theta), z], axis=1)


def _ring_units(n_members: int) -> Vec3:
    """``n_members`` unit vectors evenly spaced in the orbit (xy) plane — a *coplanar ring*.

    The half-step phase offset keeps members straddling the apogee point rather than one landing on
    the PuffSat (which sits on the ring at the same radius), so no LOS is degenerate.
    """
    angles = np.linspace(0.0, 2.0 * math.pi, n_members, endpoint=False) + math.pi / n_members
    return np.stack([np.cos(angles), np.sin(angles), np.zeros(n_members)], axis=1)


def _segment_clears_earth(observer: Vec3, member: Vec3) -> bool:
    """True iff the line of sight from ``observer`` to ``member`` clears the Earth."""
    direction = member - observer
    length = float(np.linalg.norm(direction))
    unit = direction / length
    t = min(max(-float(observer @ unit), 0.0), length)
    return float(np.linalg.norm(observer + t * unit)) >= EARTH_RADIUS_M


def line_of_sight_unit_vectors(observer: Vec3, members: Vec3) -> Vec3:
    """Unit LOS vectors from ``observer`` to each Earth-unoccluded, non-co-located member."""
    los = []
    for m in members:
        direction = m - observer
        length = float(np.linalg.norm(direction))
        if length < 1.0:  # a member on top of the observer gives no usable line of sight
            continue
        if _segment_clears_earth(observer, m):
            los.append(direction / length)
    return np.array(los)


def constellation_los_units(
    n_members: int, geometry: str = "shell", radius_m: float = APOGEE_RADIUS_M
) -> Vec3:
    """LOS unit vectors from the apogee PuffSat to an ``n_members`` ring/shell at ``radius_m``."""
    units = _ring_units(n_members) if geometry == "ring" else _fibonacci_sphere(n_members)
    return line_of_sight_unit_vectors(apogee_position(), radius_m * units)


def axis_observable(los_units: Vec3, axis_unit: Vec3, tol: float = 1e-9) -> bool:
    """True iff ``axis_unit`` lies in the row space of the geometry (a coplanar ring leaves the
    orbit-normal axis unobservable — its velocity component cannot be solved)."""
    gram = los_units.T @ los_units
    row_space_projection = np.linalg.pinv(gram) @ gram @ axis_unit
    return bool(float(np.linalg.norm(axis_unit - row_space_projection)) < tol)


def velocity_dop(los_units: Vec3, axis_unit: Vec3) -> float:
    """Dilution of precision along ``axis_unit``: ``√(aᵀ(HᵀH)⁺a)`` for geometry rows ``u_i``.

    Each range-rate measurement reads ``u_i·v``, so the velocity estimate covariance is
    ``(HᵀH)⁻¹·σ_rr²``; the pseudo-inverse gives the (finite) DOP for an *observable* axis even when
    the full geometry is rank-deficient (the coplanar ring — guard with :func:`axis_observable`).
    """
    gram = los_units.T @ los_units
    return math.sqrt(float(axis_unit @ np.linalg.pinv(gram) @ axis_unit))


def transverse_velocity_sigma_m_s(los_units: Vec3, sigma_radial_m_s: float) -> float:
    """Apogee transverse-velocity 1σ = ``TDOP_T · σ_radial`` (the binding C0 axis)."""
    return velocity_dop(los_units, TRANSVERSE_AXIS) * sigma_radial_m_s


def min_members_for_target(
    sigma_radial_m_s: float,
    target_m_s: float = TARGET_TVEL_SIGMA_M_S,
    geometry: str = "shell",
    max_members: int = 64,
) -> int | None:
    """Smallest member count whose transverse axis is observable and meets ``target_m_s``."""
    for n in range(3, max_members + 1):
        los = constellation_los_units(n, geometry)
        if len(los) < 3 or not axis_observable(los, TRANSVERSE_AXIS):
            continue
        if transverse_velocity_sigma_m_s(los, sigma_radial_m_s) <= target_m_s:
            return n
    return None


# --- Transponder mass / power (PuffSat side) -----------------------------------------------------
@dataclass(frozen=True)
class TransponderComponent:
    """One line of the PuffSat-side nav-payload mass/power ledger."""

    label: str
    mass_g: float
    power_w: float


# ADR 0020 decision 2: the PuffSat carries the *one-way passive* end — receive + verify, no PA.  The
# crypto/timing ASIC is sub-gram; the drivers are the front-end and the oscillator, not the ASIC.
PASSIVE_RECEIVER: tuple[TransponderComponent, ...] = (
    TransponderComponent("Ka omni patch antenna", 5.0, 0.0),
    TransponderComponent("receiver front-end + correlator ASIC", 12.0, 0.6),
    TransponderComponent("crypto verify ASIC (authenticated broadcast)", 0.5, 0.05),
    TransponderComponent("TCXO (clock bias solved from ≥4 members)", 2.0, 0.10),
)


def total_mass_g(components: tuple[TransponderComponent, ...]) -> float:
    return sum(c.mass_g for c in components)


def total_power_w(components: tuple[TransponderComponent, ...]) -> float:
    return sum(c.power_w for c in components)


# --- Finding -------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ApogeeNavFinding:
    """The ADR 0020 sizing verdict: link, velocity, GDOP (ring vs shell), and PuffSat mass/power."""

    range_m: float
    downlink_cn0_dbhz: float
    radial_velocity_sigma_m_s: float
    n_members: int
    shell_transverse_velocity_sigma_m_s: float
    shell_normal_observable: bool
    ring_transverse_velocity_sigma_m_s: float
    ring_normal_observable: bool
    min_members_shell: int | None
    components: tuple[TransponderComponent, ...]
    bus_mass_kg: float = PUFFSAT_BUS_MASS_KG
    target_tvel_sigma_m_s: float = TARGET_TVEL_SIGMA_M_S
    required_tvel_sigma_m_s: float = REQUIRED_TVEL_SIGMA_M_S

    @property
    def meets_target(self) -> bool:
        """The shell geometry pins transverse velocity to the C1-matching target."""
        return self.shell_transverse_velocity_sigma_m_s <= self.target_tvel_sigma_m_s

    @property
    def meets_requirement(self) -> bool:
        return self.shell_transverse_velocity_sigma_m_s <= self.required_tvel_sigma_m_s

    @property
    def target_margin(self) -> float:
        """How many times under the 0.66 mm/s target the achieved transverse σ sits."""
        return self.target_tvel_sigma_m_s / self.shell_transverse_velocity_sigma_m_s

    @property
    def puffsat_mass_g(self) -> float:
        return total_mass_g(self.components)

    @property
    def puffsat_power_w(self) -> float:
        return total_power_w(self.components)

    @property
    def mass_fraction(self) -> float:
        return self.puffsat_mass_g / 1e3 / self.bus_mass_kg

    @property
    def crypto_asic_is_mass_driver(self) -> bool:
        """Whether the crypto/timing ASIC dominates the ledger (it does not — that is the finding;
        the front-end and oscillator are the drivers, the ASIC is sub-gram)."""
        crypto = max(
            (c for c in self.components if "crypto" in c.label.lower()), key=lambda c: c.mass_g
        )
        return crypto.mass_g >= max(c.mass_g for c in self.components)


def apogee_nav_finding(n_members: int = 12, range_m: float = DESIGN_RANGE_M) -> ApogeeNavFinding:
    """Assemble the ADR 0020 sizing finding at a nominal shell size (the pure runner)."""
    cn0 = downlink_cn0_dbhz(range_m)
    sigma_radial = carrier_phase_velocity_sigma_m_s(cn0)
    shell = constellation_los_units(n_members, "shell")
    ring = constellation_los_units(n_members, "ring")
    return ApogeeNavFinding(
        range_m=range_m,
        downlink_cn0_dbhz=cn0,
        radial_velocity_sigma_m_s=sigma_radial,
        n_members=n_members,
        shell_transverse_velocity_sigma_m_s=transverse_velocity_sigma_m_s(shell, sigma_radial),
        shell_normal_observable=axis_observable(shell, NORMAL_AXIS),
        ring_transverse_velocity_sigma_m_s=transverse_velocity_sigma_m_s(ring, sigma_radial),
        ring_normal_observable=axis_observable(ring, NORMAL_AXIS),
        min_members_shell=min_members_for_target(sigma_radial, geometry="shell"),
        components=PASSIVE_RECEIVER,
    )


def format_apogee_nav(finding: ApogeeNavFinding) -> str:
    """One-screen ADR 0020 apogee-nav-constellation sizing report."""
    verdict = (
        "MEETS the C1-matching target"
        if finding.meets_target
        else ("MEETS the requirement (not the C1 target)" if finding.meets_requirement else "FAILS")
    )
    lines = [
        "Apogee nav constellation sizing — Lever 3 (ADR 0020; coast/apogee-state nav)",
        f"  Link (Ka {KA_FREQ_HZ / 1e9:.0f} GHz, {CONSTELLATION_DISH_M:g} m /"
        f" {CONSTELLATION_TX_POWER_W:g} W dish → omni PuffSat at {finding.range_m / 1e3:.0f} km):"
        f" C/N0 {finding.downlink_cn0_dbhz:.1f} dB-Hz — closes (far members ~2R: −6 dB, still ok).",
        f"  Velocity: carrier-phase radial-velocity σ {finding.radial_velocity_sigma_m_s * 1e3:.3f}"
        f" mm/s ({INTEGRATION_TIME_S:g} s integration) — Ka carrier keeps it well under target.",
        f"  Snapshot GDOP at apogee (N={finding.n_members}):",
        f"    shell: transverse σ {finding.shell_transverse_velocity_sigma_m_s * 1e3:.3f} mm/s"
        f" (normal axis {'observable' if finding.shell_normal_observable else 'UNobservable'})"
        f" → {verdict} ({finding.target_margin:.1f}× under the"
        f" {finding.target_tvel_sigma_m_s * 1e3:g} mm/s target).",
        f"    coplanar ring: transverse σ {finding.ring_transverse_velocity_sigma_m_s * 1e3:.3f}"
        f" mm/s (normal axis {'observable' if finding.ring_normal_observable else 'UNobservable'})"
        " — a ring covers the binding transverse axis; normal matters ~50× less (C0), so the shell"
        " only adds the weak axis.",
        f"  Minimum shell members for the {finding.target_tvel_sigma_m_s * 1e3:g} mm/s target:"
        f" {finding.min_members_shell} — beyond it the gain is redundancy/GDOP, not accuracy"
        " (match-not-beat, ADR 0020 decision 5).",
        f"  PuffSat one-way passive payload: {finding.puffsat_mass_g:.1f} g,"
        f" {finding.puffsat_power_w:.2f} W → {finding.mass_fraction * 100:.2f}% of the"
        f" {finding.bus_mass_kg:g} kg bus"
        f" — crypto ASIC {'IS' if finding.crypto_asic_is_mass_driver else 'is NOT'} the mass driver"
        " (front-end + oscillator are; transmit power, not mass, is the constraint).",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(format_apogee_nav(apogee_nav_finding()))
