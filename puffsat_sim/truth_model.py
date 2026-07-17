"""Truth-model runner: reference-orbit checks and per-force signature reports.

Exercises the Orekit / JVM bridge and the reference orbit, then reports each
perturbation's signature against its analytic prediction as forces are added
one at a time (design-doc Stage 1 physics ladder: two-body → J2 → third-body →
SRP → drag).

These functions PRINT — they are a human-readable demo, not the verification
surface.  The assertions live in ``tests/integration/test_propagator_physics.py``.

Run with:
    make run
or:
    python -m puffsat_sim.truth_model

The orbit used matches the near-term architecture from the paper:
  - periapsis 50 km (orbit periapsis; PuffSat burns up here after impact)
  - interception at 200 km during descent, before periapsis
  - apogee  ~150 000 km altitude (recommended deployment apogee from design doc)
  - eccentricity ~0.921, period ~2.68 days
  - perigee speed ~10.91 km/s
"""

from __future__ import annotations

import math
from typing import Any

from puffsat_sim import mission, presets
from puffsat_sim.config import OrbitalConfig, PhysicsConfig
from puffsat_sim.constants import EARTH_RADIUS_M, SUN_RADIUS_M
from puffsat_sim.forces import AtmosphericDrag, Geopotential, Relativity, SolarRadiation
from puffsat_sim.forces.drag import drag_deceleration
from puffsat_sim.forces.geopotential import (
    j2_apsidal_precession_rate,
    j2_nodal_regression_rate,
)
from puffsat_sim.forces.relativity import (
    schwarzschild_apsidal_advance_per_orbit,
    schwarzschild_perigee_advance_per_orbit_m,
)
from puffsat_sim.forces.srp import srp_acceleration
from puffsat_sim.forces.third_body import lunar_tidal_ratio, solar_tidal_ratio
from puffsat_sim.orbital_math import (
    keplerian_elements,
    keplerian_period,
    perigee_speed,
    wrap_to_pi,
)

# Importing descent starts the JVM and loads Orekit data.
# All org.orekit.* imports must follow this line.
from puffsat_sim.descent import earth_model, to_absolute_date  # noqa: E402
from puffsat_sim.propagator import build_propagator

from org.orekit.bodies import CelestialBodyFactory
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit
from org.orekit.propagation.events import AltitudeDetector, EclipseDetector, EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent, StopOnDecreasing
from org.orekit.utils import Constants


def _final_pos(
    orbital_config: OrbitalConfig, physics_config: PhysicsConfig, epoch: Any, period: float
) -> tuple[float, float, float]:
    """The one-period endpoint position, for differencing two physics configs (ADR 0017)."""
    prop = build_propagator(orbital_config, physics_config)
    pv = prop.propagate(epoch.shiftedBy(period)).getPVCoordinates().getPosition()
    return float(pv.getX()), float(pv.getY()), float(pv.getZ())


def propagate_one_period(orbital_config: OrbitalConfig, physics_config: PhysicsConfig) -> None:
    """Propagate the PuffSat reference orbit for one Keplerian period."""
    a, e = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    v_perigee = perigee_speed(a, orbital_config.perigee_alt_m)

    epoch = to_absolute_date(orbital_config.epoch)
    propagator = build_propagator(orbital_config, physics_config)
    initial_pv = propagator.getInitialState().getPVCoordinates()
    r_0 = initial_pv.getPosition()
    v_0 = initial_pv.getVelocity()

    final_state = propagator.propagate(epoch.shiftedBy(period))
    pos = final_state.getPVCoordinates().getPosition()
    vel = final_state.getPVCoordinates().getVelocity()

    # Residual after one Keplerian period (should be ~floating-point zero)
    dr = math.sqrt(
        (pos.getX() - r_0.getX()) ** 2
        + (pos.getY() - r_0.getY()) ** 2
        + (pos.getZ() - r_0.getZ()) ** 2
    )
    dv = math.sqrt(
        (vel.getX() - v_0.getX()) ** 2
        + (vel.getY() - v_0.getY()) ** 2
        + (vel.getZ() - v_0.getZ()) ** 2
    )

    print("PuffSat Control Simulation — truth model: Keplerian reference orbit")
    print("  Orekit / JVM : OK")
    print()
    print("  Reference orbit (near-term architecture):")
    print(
        f"    Orbit periapsis  : {orbital_config.perigee_alt_m / 1e3:.0f} km"
        "  (burns up here; interception at 200 km during descent)"
    )
    print(f"    Apogee altitude  : {orbital_config.apogee_alt_m / 1e6:.0f} × 10³ km  (deployment)")
    print(f"    Semi-major axis  : {a / 1e3:.1f} km")
    print(f"    Eccentricity     : {e:.6f}")
    print(f"    Inclination      : {math.degrees(orbital_config.inclination_rad):.1f}°")
    print(f"    Orbital period   : {period:.1f} s  ({period / 86400:.2f} days)")
    print(f"    Perigee speed    : {v_perigee / 1e3:.3f} km/s")
    print()
    print("  One-period propagation residual (Keplerian → should be ~0):")
    print(f"    |Δr| = {dr:.3e} m")
    print(f"    |Δv| = {dv:.3e} m/s")


def propagate_to_interception(orbital_config: OrbitalConfig, physics_config: PhysicsConfig) -> None:
    """Propagate from apogee, stopping at the 200 km descent crossing (interception).

    Uses AltitudeDetector with StopOnDecreasing: the g-function is
    (altitude − 200 km), which decreases through zero as the PuffSat descends.
    Starting at apogee, the first zero-crossing is the descending one, so no
    additional filtering is needed.
    """
    frame = FramesFactory.getEME2000()
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)

    epoch = to_absolute_date(orbital_config.epoch)
    propagator = build_propagator(orbital_config, physics_config)

    # WGS84 ellipsoid in the Earth-fixed frame — used by AltitudeDetector to
    # compute geodetic altitude above the surface.
    earth = earth_model()
    detector = AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(
        StopOnDecreasing()  # type: ignore[no-untyped-call]
    )
    propagator.addEventDetector(detector)

    # Upper bound is one full period; event fires ~halfway through (apo → peri descent).
    final_state = propagator.propagate(epoch.shiftedBy(period))

    elapsed: float = final_state.getDate().durationFrom(epoch)
    pos = final_state.getPVCoordinates().getPosition()
    vel = final_state.getPVCoordinates().getVelocity()
    v: float = math.sqrt(vel.getX() ** 2 + vel.getY() ** 2 + vel.getZ() ** 2)
    geodetic = earth.transform(pos, frame, final_state.getDate())
    alt_km: float = geodetic.getAltitude() / 1e3

    print("  Propagation to interception (200 km descent crossing):")
    print(f"    Coast time from apogee : {elapsed / 3600:.3f} h  ({elapsed / 86400:.3f} days)")
    print(f"    Altitude at stop       : {alt_km:.3f} km  (event target: 200 km)")
    print(f"    Speed at interception  : {v / 1e3:.3f} km/s")


def report_j2_signatures(orbital_config: OrbitalConfig) -> None:
    """Propagate one period with J2 and compare ΔRAAN / Δω against analytic predictions.

    Checks that the NumericalPropagator + HolmesFeatherstoneAttractionModel (degree 2)
    produces nodal regression and apsidal precession rates consistent with the
    first-order J2 secular formulas.  Agreement within ~1% confirms the force model,
    integrator tolerances, and rate formulas are mutually consistent.
    """
    a, e = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    i = orbital_config.inclination_rad

    epoch = to_absolute_date(orbital_config.epoch)
    propagator = build_propagator(orbital_config, presets.j2())

    initial_state = propagator.getInitialState()
    initial_orbit = KeplerianOrbit(initial_state.getOrbit())
    raan_0: float = initial_orbit.getRightAscensionOfAscendingNode()
    omega_0: float = initial_orbit.getPerigeeArgument()

    final_state = propagator.propagate(epoch.shiftedBy(period))
    final_orbit = KeplerianOrbit(final_state.getOrbit())
    raan_f: float = final_orbit.getRightAscensionOfAscendingNode()
    omega_f: float = final_orbit.getPerigeeArgument()

    # Wrap to (-π, π]: the osculating angles straddle the 0/2π branch cut, so a
    # raw subtraction reports a small retrograde drift as ~+2π instead of ~−ε.
    d_raan_num = math.degrees(wrap_to_pi(raan_f - raan_0))
    d_omega_num = math.degrees(wrap_to_pi(omega_f - omega_0))

    rate_raan = j2_nodal_regression_rate(a, e, i)
    rate_omega = j2_apsidal_precession_rate(a, e, i)
    d_raan_analytic = math.degrees(rate_raan * period)
    d_omega_analytic = math.degrees(rate_omega * period)

    pct_raan = abs(d_raan_num - d_raan_analytic) / abs(d_raan_analytic) * 100.0
    pct_omega = abs(d_omega_num - d_omega_analytic) / abs(d_omega_analytic) * 100.0

    print("  J2 signature report (geopotential degree 2):")
    print(
        f"    ΔRAAN  numeric: {d_raan_num:+.4f}°  "
        f"analytic: {d_raan_analytic:+.4f}°  agree: {100.0 - pct_raan:.2f}%"
    )
    print(
        f"    Δω     numeric: {d_omega_num:+.4f}°  "
        f"analytic: {d_omega_analytic:+.4f}°  agree: {100.0 - pct_omega:.2f}%"
    )


def report_higher_order_gravity_signatures(orbital_config: OrbitalConfig) -> None:
    """Report the non-J2 geopotential (the truth model uses an 8×8 field).

    Beyond J2 the harmonics fall off as (Rₑ/r)^ℓ — dead at the apogee, biting only
    near perigee — yet over a pass they reach ~km scale at orbit level.  Confirms
    the higher-degree field is wired by differencing one-period endpoints: J2 vs
    8×8 (full) and J2 vs degree-8 zonal-only (tesserals dominate the difference).
    """
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    epoch = to_absolute_date(orbital_config.epoch)

    pos_j2 = _final_pos(orbital_config, presets.j2(), epoch, period)
    pos_8x8 = _final_pos(
        orbital_config, PhysicsConfig((Geopotential(degree=8, order=8),)), epoch, period
    )
    pos_8z = _final_pos(
        orbital_config, PhysicsConfig((Geopotential(degree=8, order=0),)), epoch, period
    )
    dr_full = math.sqrt(sum((p - q) ** 2 for p, q in zip(pos_j2, pos_8x8, strict=True)))
    dr_zonal = math.sqrt(sum((p - q) ** 2 for p, q in zip(pos_j2, pos_8z, strict=True)))

    print("  Higher-degree gravity signature report (truth field 8×8):")
    print(f"    One-period endpoint drift J2 → 8×8 (full)    : {dr_full / 1e3:.3f} km")
    print(f"    One-period endpoint drift J2 → degree-8 zonal: {dr_zonal / 1e3:.3f} km")


def report_third_body_signatures(orbital_config: OrbitalConfig) -> None:
    """Report third-body Sun + Moon perturbations.

    Confirms that ThirdBodyAttraction is active by comparing one-period
    propagation with J2-only vs J2+Sun+Moon, then reports the analytic tidal
    acceleration ratios at apogee (Hill approximation, design doc §2 benchmark).
    """
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    epoch = to_absolute_date(orbital_config.epoch)

    pos_j2 = _final_pos(orbital_config, presets.j2(), epoch, period)
    pos_j2_tb = _final_pos(orbital_config, presets.j2_third_body(), epoch, period)

    dr = math.sqrt(sum((a - b) ** 2 for a, b in zip(pos_j2, pos_j2_tb, strict=True)))

    moon_pct = lunar_tidal_ratio(orbital_config.apogee_alt_m) * 100.0
    sun_pct = solar_tidal_ratio(orbital_config.apogee_alt_m) * 100.0

    print("  Third-body signature report (Sun + Moon):")
    print(f"    Tidal ratio at apogee  — Moon: {moon_pct:.3f}%  Sun: {sun_pct:.3f}%")
    print(f"    One-period position drift (J2 → J2+Sun+Moon): {dr / 1e3:.3f} km")


def report_srp_signatures(orbital_config: OrbitalConfig) -> None:
    """Report SRP cannonball force model and shadow/eclipse detection.

    Confirms SRP is active by comparing one-period propagation with J2+third_body
    vs J2+third_body+SRP, then counts eclipse entry/exit events using EclipseDetector
    to verify the cylindrical shadow model is wired.  Prints the analytic SRP
    acceleration at 1 AU as the order-of-magnitude benchmark.
    """
    physics_srp = presets.j2_third_body_srp()
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    epoch = to_absolute_date(orbital_config.epoch)

    pos_tb = _final_pos(orbital_config, presets.j2_third_body(), epoch, period)
    pos_srp = _final_pos(orbital_config, physics_srp, epoch, period)
    dr = math.sqrt(sum((a - b) ** 2 for a, b in zip(pos_tb, pos_srp, strict=True)))

    # Eclipse detection: count shadow entry/exit crossings during one period.
    sun = CelestialBodyFactory.getSun()
    earth = earth_model()
    eclipse_prop = build_propagator(orbital_config, physics_srp)
    logger: EventsLogger = EventsLogger()  # type: ignore[no-untyped-call]
    # EclipseDetector defaults to StopOnIncreasing, which would halt the arc at
    # the first umbra exit; ContinueOnEvent lets it log every crossing across the
    # whole period so the count is meaningful.
    eclipse_prop.addEventDetector(
        logger.monitorDetector(
            EclipseDetector(sun, SUN_RADIUS_M, earth).withHandler(
                ContinueOnEvent()  # type: ignore[no-untyped-call]
            )
        )
    )
    eclipse_prop.propagate(epoch.shiftedBy(period))
    n_eclipse = len(logger.getLoggedEvents())

    srp_spec = next(p for p in physics_srp.perturbations if isinstance(p, SolarRadiation))
    cr_am = srp_spec.cr_area_over_mass
    accel = srp_acceleration(cr_am)

    print("  SRP signature report (cannonball):")
    print(f"    Analytic SRP accel (1 AU, Cr·A/m={cr_am}): {accel:.3e} m/s²")
    print(f"    One-period divergence (J2+3body → +SRP): {dr / 1e3:.3f} km")
    print(f"    Eclipse crossings detected in one period: {n_eclipse}")


def report_drag_signatures(orbital_config: OrbitalConfig) -> None:
    """Report NRLMSISE-00 atmospheric drag.

    Propagates from apogee to the 200 km descent crossing with the full force model
    and with the same model minus drag.  Isolating drag this way (rather than
    against the lower-fidelity j2_third_body_srp) keeps the geopotential and
    relativity identical, so the energy/speed difference is the work drag removes,
    not a conservative-potential mismatch.  The analytic deceleration at key
    altitudes gives the expected order of magnitude.

    Design doc §4: drag "bites below ~300-400 km" — shown by comparing energy at
    200 km (interception altitude) with and without drag.
    """
    physics_full = presets.full_force()
    physics_no_drag = PhysicsConfig(
        tuple(p for p in physics_full.perturbations if not isinstance(p, AtmosphericDrag))
    )
    a, _ = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    mu: float = float(Constants.WGS84_EARTH_MU)

    def _energy_and_speed_at_200km(physics_config: PhysicsConfig) -> tuple[float, float]:
        earth = earth_model()
        prop = build_propagator(orbital_config, physics_config)
        prop.addEventDetector(
            AltitudeDetector(mission.INTERCEPTION_ALT_M, earth).withHandler(
                StopOnDecreasing()  # type: ignore[no-untyped-call]
            )
        )
        epoch = to_absolute_date(orbital_config.epoch)
        state = prop.propagate(epoch.shiftedBy(period))
        pv = state.getPVCoordinates()
        pos = pv.getPosition()
        vel = pv.getVelocity()
        r = math.sqrt(pos.getX() ** 2 + pos.getY() ** 2 + pos.getZ() ** 2)
        v_sq = vel.getX() ** 2 + vel.getY() ** 2 + vel.getZ() ** 2
        energy = 0.5 * v_sq - mu / r
        return energy, math.sqrt(v_sq)

    # Both arcs stop at the same 200 km altitude (the event fixes position and time),
    # so the meaningful drag signatures are the loss of specific orbital energy and the
    # resulting speed reduction at the crossing — not the position, which barely moves.
    energy_no_drag, speed_no_drag = _energy_and_speed_at_200km(physics_no_drag)
    energy_full, speed_full = _energy_and_speed_at_200km(physics_full)
    d_energy = energy_full - energy_no_drag
    d_speed = speed_full - speed_no_drag

    drag_spec = next(p for p in physics_full.perturbations if isinstance(p, AtmosphericDrag))
    cd_am = drag_spec.cd_area_over_mass
    v_200 = math.sqrt(mu * (2.0 / (EARTH_RADIUS_M + 200_000.0) - 1.0 / a))
    v_300 = math.sqrt(mu * (2.0 / (EARTH_RADIUS_M + 300_000.0) - 1.0 / a))
    a_drag_200 = drag_deceleration(cd_am, v_200, 200_000.0)
    a_drag_300 = drag_deceleration(cd_am, v_300, 300_000.0)

    print("  Drag signature report (NRLMSISE-00):")
    print(f"    Analytic drag decel @ 200 km: {a_drag_200:.3e} m/s²")
    print(f"    Analytic drag decel @ 300 km: {a_drag_300:.3e} m/s²")
    print(f"    Orbital energy at 200 km — no drag: {energy_no_drag:.1f} J/kg")
    print(f"    Orbital energy at 200 km — with drag: {energy_full:.1f} J/kg")
    print(f"    ΔE from drag: {d_energy:.1f} J/kg  (negative = drag removed energy)")
    print(f"    Speed reduction at 200 km from drag: {d_speed * 100.0:+.3f} cm/s")


def report_relativity_signatures(orbital_config: OrbitalConfig) -> None:
    """Report the Schwarzschild relativistic perturbation.

    Isolates relativity by differencing one-period endpoints with and without it
    (J2 baseline vs J2+relativity); the shared J2 and integrator settings make the
    integration errors common-mode, so the difference is the relativity signal.
    Prints the closed-form apsidal advance and its ~cm/orbit perigee displacement.

    Relativity is negligible at the orbit-level (km-to-m) scale but ~cm per pass on
    this high-e orbit — the deferred 5 cm terminal-centering budget (design doc §7).
    """
    a, e = keplerian_elements(orbital_config.perigee_alt_m, orbital_config.apogee_alt_m)
    period = keplerian_period(a)
    epoch = to_absolute_date(orbital_config.epoch)

    pos_j2 = _final_pos(orbital_config, presets.j2(), epoch, period)
    pos_j2_rel = _final_pos(
        orbital_config, PhysicsConfig((Geopotential(degree=2), Relativity())), epoch, period
    )
    dr = math.sqrt(sum((p - q) ** 2 for p, q in zip(pos_j2, pos_j2_rel, strict=True)))

    d_pomega = schwarzschild_apsidal_advance_per_orbit(a, e)
    d_perigee_cm = schwarzschild_perigee_advance_per_orbit_m(a, e) * 100.0
    arcsec_per_orbit = math.degrees(d_pomega) * 3600.0

    print("  Relativity signature report (Schwarzschild):")
    print(
        f"    Analytic apsidal advance : {d_pomega:.3e} rad/orbit  ({arcsec_per_orbit:.3e} arcsec)"
    )
    print(f"    Perigee displacement     : {d_perigee_cm:.2f} cm/orbit")
    print(f"    One-period endpoint drift (J2 → J2+relativity): {dr:.3f} m")


def main() -> None:
    propagate_one_period(mission.NOMINAL_CONFIG, presets.two_body())
    print()
    propagate_to_interception(mission.NOMINAL_CONFIG, presets.two_body())
    print()
    report_j2_signatures(mission.NOMINAL_CONFIG)
    print()
    report_higher_order_gravity_signatures(mission.NOMINAL_CONFIG)
    print()
    report_third_body_signatures(mission.NOMINAL_CONFIG)
    print()
    report_srp_signatures(mission.NOMINAL_CONFIG)
    print()
    report_drag_signatures(mission.NOMINAL_CONFIG)
    print()
    report_relativity_signatures(mission.NOMINAL_CONFIG)


if __name__ == "__main__":
    main()
