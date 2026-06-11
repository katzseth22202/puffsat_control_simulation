# Plate-capture success criterion, GNSS-free beacon/astrometric navigation, aim-bias torque management

**Status:** accepted

## Context

A grilling session (2026-06-11) stress-tested a proposed target respec: back away from
the parked ~5 cm plate-centering toward "2 m cross/radial × 25 m along-track," delete
GNSS hardware from the PuffSat, track a PuffSat homing beacon from ground stations and
the target's pusher plate ("especially accurate in the last 50 km"), and counter
off-center plate impacts with target-side RCS.

The grill found that half the proposal had already happened — ADR 0014 re-baselined the
terminal story to a ~400 m thrust-limited catch radius + ~1 m nav-limited endpoint, with
5 cm parked behind integrator work that exists only to serve it — and that each proposed
number transformed under arithmetic:

- **"2 m radius" was accidentally a tightening.** As a 3σ bound it demands σ ≤ 0.67 m,
  below ADR 0014's 1 m floor, while deleting the sensor that achieves 1 m. Resolved by
  widening the plate: it must shadow the cabin anyway (no PuffSat may strike the vehicle
  body), and plate mass is set by pulse smoothing, not area — **M_p ≈ J·τ/2s** (arrest a
  J ≈ 2.7×10⁵ N·s hit over stroke s in time τ): ~54 t at τ = 1 s / s = 2.5 m, ~9–13 t
  with a long 10–15 m absorber stroke. A 5 m-radius disc is mass-free given that budget.
  Off-center hits tip the plate (~1.3 rad/s at 54 t; harder when lighter), so push-pull
  (tensile) struts are mandatory, not a nicety.
- **"25 m along-track" is not a position miss.** At closest approach the miss vector is
  ⊥ v_rel by definition; the along-track axis is **time-of-arrival** (25 m ↔ 2.3 ms at
  10.8 km/s) — 100× tighter than the ~1 s pulse cadence and damper reset appear to need.
- **"Accurate in the last 50 km" cannot be consumed by the PuffSat.** Thrust-limited
  lateral authority ½·a_max·(R/v)² at a_max = 0.016 m/s² is **0.17 m at 50 km range**
  (0.69 m at 100 km, 2.0 m at 170 km, 6.2 m at 300 km). Knowledge must arrive by
  ~200–300 km range to be actionable at the meter scale. Angle-based knowledge improves
  as σ_θ·R while authority dies as R²; the crossover sets a homing floor
  **σ_miss ≈ 2σ_θ²·v²/a_max** — 1.46 m at σ_θ = 10 µrad, 0.36 m at 5 µrad, 36 m
  (fails) at 50 µrad. Ground stations (instrumentation-radar ~50–100 µrad at hundreds of
  km slant range) are funnel-entry class, not terminal class.
- **Per-hit torque is too big for RCS, and a torque *engine* taxes the mission.** A
  typical 2 m-offset hit deposits J·d ≈ 5.4×10⁵ N·m·s (≈2.8°/s on a 150 t / 30 m
  vehicle). Canceling within the 1 s cadence needs ~27 kN at a 20 m front arm — easy
  *thrust* (AJ10-class chemistry) but the impulse fraction has a closed form:
  **torque propellant / delivered momentum = d̄/arm ≈ 2 m/20 m = 10%** (~9–16 t per
  boost train, a continuous sideways 27 kN burn) — in direct tension with the
  propellant-free premise. Meanwhile each *scheduled* hit can carry up to ~5×10⁵ N·m·s
  of commanded counter-torque for free by aiming it off-center on purpose.
- **C1 already measured the midcourse GNSS-delete.** Coordinator-beacon nav (no GNSS)
  reaches the 800 km hand-off at 141 m lateral — 2.8× inside the 400 m radius. But at
  400 m the **range-only fallback is dead** (386 m, margin 1.04×): in a GNSS-free
  flight, **two-way Doppler is load-bearing**, not optional.

## Decision

1. **The success criterion is plate capture, not a centering σ.** Plate radius **5 m**
   (spec: ≥ cabin radius + L·tan α_max approach-angle overhang + margin — shadowing is a
   safety requirement: any miss beyond the plate edge is a clean flyby and a burn-up,
   by design). Requirement: **≥99% per-PuffSat capture ↔ σ_lateral ≤ 1.65 m**
   (2D-Rayleigh, R/σ = 3.03), plus **ToA ≤ ~10 ms** at closest approach (derived from
   pulse cadence/damper reset, not from a position tolerance). Capture probability is an
   *economics* output (missed PuffSats deorbit themselves, paper §9), not a pass/fail
   gate; the goal is capped at 99% because the ~0.3% funnel-saturation tail (ADR 0014's
   MCC-2 case) already consumes part of any tighter budget.

2. **GNSS is deleted from the PuffSat.** Midcourse stands on the measured C1
   coordinator-beacon architecture **with Doppler required** (range-only is dead at the
   400 m radius). The honest rationale recorded for the paper: COCOM (unlocked receivers
   × thousands of expendable units), retirement of ADR 0014's flagged 10.8 km/s
   carrier-tracking stretch item, and free time transfer via two-way ranging. The
   proposal's "lower latency / higher power than GNSS" justification is **rejected**
   (GNSS is passive; latency was never its failure mode) — crosslink latency lands in
   C4's τ budget where it belongs.

3. **Terminal navigation: the PuffSat homes on angles from a target-side tracker.**
   Derived hardware requirement: **σ_θ ≤ 10 µrad** (meets 1.65 m with thin margin);
   **5 µrad is the design target** (0.36 m floor, comfortable RSS headroom). Tracker
   architecture: **beacon-vs-star-background differential astrometry in the same focal
   plane** — self-referencing (platform attitude error cancels between exposures; the
   measurement is absolute, which also kills the common-mode aim-bias failure mode
   below), ~1 ms exposures on a bright laser beacon keep worst post-impact smear
   ~3 µrad, nav-grade gyros bridge between ~1 Hz frames. A GOCE-class accelerometer on
   the tracker platform is **rejected** (wrong axis — linear, not angular — and
   ~10⁵–10⁶× over its full scale in a 1–2 m/s² pulsed environment; the ADR 0013 fiction,
   re-buried). Plume flash at each impact is handled by **scheduled gating** (the ToA
   spec makes impacts predictable to a frame time) + a narrowband filter at the beacon
   wavelength; the residual risk is detector saturation recovery — bench-testable
   hardware qualification, not a stretch item.

4. **Sim encoding (C3b/C4/Rung D deltas; everything else stands).** The ADR 0014 σ_rel
   injection becomes **range-dependent: σ_rel(R) = σ_θ·R**, swept
   **σ_θ ∈ {2, 10, 50 µrad}** (+ the constant code-differential 1 m point retained for
   ADR 0014 continuity) — "which tracker grade buys the claim" is the measured output,
   the C1 range-only pattern. New pure summarizer: **plate-frame miss decomposition**
   (2D lateral ⊥ v_rel + ToA) and the **capture-vs-plate-radius curve**; Rung D's
   headline becomes **P(capture)**. Optional measurement-dropout knob (impact outages)
   feeds C4. ZEM law, 400 m catch radius, C3a/C3b/C3c slices, MCC-2 — all unchanged.

5. **Retired outright:** the parked 5 cm centering; its numerical-fidelity prerequisite
   (`rel_tol` 1e-13 + Encke/element-formulation study — 1.65 m is measurable at the
   current 1e-10 ~cm crossing floor); the relativity-in-the-onboard-filter question.

6. **Torque management (paper-side architecture note): hybrid aim-bias + non-toxic trim
   engine.** Primary: **closed-loop aim-point biasing** — the next scheduled hit cancels
   accumulated angular momentum (zero propellant; the 5 m plate provides the aim room).
   Steady-state residual is one hit's noise: ±~1.4° vehicle jitter at ~1 Hz at
   σ_lateral = 1.65 m, <1° at the 5 µrad tracker target — the cheapest jitter lever is
   tracker grade, not propellant. Secondary: a **few-kN-class non-toxic trim engine**
   (loop start-up, dropout, end-of-train cleanup, optional jitter damping; sub-tonne
   propellant in this role, so the Isp penalty vs hypergolics is immaterial).
   **Hypergolics are excluded**: the vehicle is a passenger craft with aircraft-style
   turnaround, and hypergolic ground servicing (SCAPE-suit operations) is the
   non-starter — crewed *flight* with hypergolics is routine (Dragon/Orion), gate
   operations with them are not. Named candidates: **H₂O₂** (X-15/Mercury/Soyuz crewed
   heritage, storable, pulse-friendly) or **LOX/methane** if the engine ever grows into
   the continuous jitter-damping role. Engine-primary is rejected via the d̄/arm ≈ 10%
   momentum tax. The sim's deliverables to this ledger: the offset distribution
   (σ_lateral, capture %) and its **correlation structure** (independent offsets →
   √N random walk; common-mode tracker bias → linear accumulation — mitigated by the
   absolute star-referenced measurement and owned by the aim-bias loop).

## Considered options

- **2 m plate radius / 3σ semantics** — rejected. As 3σ it silently *tightens* the
  requirement below the achievable nav floor; as 1σ a 2 m plate captures only 39%. The
  5 m plate dissolves the dilemma and is mass-free given the smoothing mass.
- **Keep relative GNSS as the terminal baseline (ADR 0014)** — rejected in favor of the
  beacon/astrometric architecture; kept as the named comparison point in the σ sweep.
  The swap retires the carrier-tracking stretch item and COCOM burden at the cost of two
  new, *bench-testable* items (µrad tracker on a hammered vehicle; saturation recovery).
- **Ground stations as the terminal nav source** — rejected for terminal (angle error
  scales with slant range; ~tens of m), retained as funnel-entry/MCC-2 support.
- **Target-side terminal trim as baseline** (rocket translates the plate the last 1–2 m)
  — kept as fallback/margin only; it puts active GNC on a passenger vehicle during
  ascent and would add a coupled two-vehicle terminal loop to the sim scope.
- **RCS-only torque cancellation** (the original proposal, re-tested at the 30 m /
  front-RCS geometry) — rejected: ~7× per-hit shortfall at 4 kN/20 m, several °/s
  stationary wander attacking the shadowing angle budget, ~0.3 m/s lateral kick per
  cancellation into the homing corridor.
- **Torque-cancellation engine as primary (27 kN continuous)** — rejected: thrust is
  easy, impulse is not; d̄/arm ≈ 10% of all delivered momentum (~9–16 t hypergolic-class
  propellant per boost train) contradicts the propellant-free premise.
- **GOCE-quality accelerometers for tracker platform knowledge** — rejected (wrong
  measurement axis; ~10⁵–10⁶× over full scale in the pulse environment).
- **Keep "25 m along-track"** — rejected; respecced as ToA with a traceable parent
  (pulse cadence), which is *looser* (≈10 ms ↔ ~100 m) and honestly derived.

## Consequences

- Design doc edits: §2 out-of-scope note (target-RCS loop → ADR 0015 torque note),
  §3 "Two precision scales" (terminal level re-specced; 5 cm retired), §13 C3
  (plate-capture criterion + σ_θ sweep + retirements), §16.3 (crosslink note superseded:
  no PuffSat GNSS; target-side astrometric tracker), §16.4 (terminal nav grade re-keyed
  from σ_rel GNSS grades to σ_θ), §13 numerical-fidelity paragraph (retired, kept for
  the record).
- New pure code planned (with C3b): plate-frame miss summarizer + capture-vs-radius
  curve; the σ_rel(R) = σ_θ·R profile type; the rest of C3 as per ADR 0014.
- C1/C2/C0 measured results unchanged; the C1 re-read at 400 m gains one sentence:
  range-only is no longer a viable fallback (Doppler load-bearing in GNSS-free flight).
- Rung D inherits: P(capture) as the headline metric, σ_θ as a sampled axis, offset
  correlation structure as a reported output for the paper's torque ledger.
- Paper-side pins (outside this repo, flagged not assumed): β (momentum-transfer
  efficiency 1–2×, everything above scales linearly), vehicle moment of inertia for the
  torque ledger, plume-capture-driven plate sizing (may exceed 5 m), v_rel refinement
  from the actual ascent profile, the non-toxic trim-engine selection.
