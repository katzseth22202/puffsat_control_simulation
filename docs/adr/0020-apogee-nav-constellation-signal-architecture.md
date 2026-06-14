# Apogee nav constellation (Lever 3): a Ka-band authenticated-broadcast signal, match-not-beat accuracy

**Status:** accepted

## Context

ADR 0019 re-architected **terminal** relative navigation (the 800→200 km homing phase) into a fused
multi-tracker grade and recorded **Lever 3 — a permanent apogee-regime nav constellation (~150,000 km,
above GPS)** as an *option, not built*. Lever 3 addresses a **different regime**: the **coast /
apogee-state navigation** the midcourse corrector consumes. That regime is not optional — the
corrector must pin the apogee state regardless, and ADR 0019 left it untouched.

Today that coast nav is the C1 coordinator nodes (ADR 0012). The binding requirement is **apogee
transverse velocity** (C0/ADR 0011): C1 achieved **σ_Tvel ≈ 0.66 mm/s, lateral ≈ 140 m**, and that
140 m *is* the per-unit hand-off entry offset D1.1 flies (`train.ENTRY_LATERAL_PERUNIT_M`, = 0.66 mm/s
× Φ 2.15e5). GNSS cannot serve this regime: GPS sits at ~20,200 km, apogee at ~150,000 km — GNSS owns
*low* altitude only, and a single co-flyer is not a multilateration set (and can't fix its own absolute
state at apogee). So coast nav needs coordinator nodes **or** this constellation.

A grilling session (2026-06-14) took the constellation's *existence* as given (per the user) and asked
the **signal/transponder sizing** questions: same technology as GNSS? a higher-frequency radio (no
atmosphere to penetrate)? laser? is a crypto/ASIC transponder with a nanosecond-exact turnaround delay
a realistic anti-jam scheme, and does it work for the co-flyer/target rockets without being too heavy on
the PuffSat? And the accuracy question: now that we *could* pin tighter than 140 m, is it worth it?

## Decision

1. **Ka-band (~30 GHz), spread-spectrum RF broadcast — not L-band GNSS, not laser.** The method
   (TOF pseudoranging + multilateration) is altitude-agnostic, but the implementation is re-keyed for
   deep space. With no atmosphere/ionosphere between shell and PuffSat, the **dominant GNSS error term
   (ionospheric delay) vanishes** and high frequency is free; Ka buys what the *binding* axis needs —
   **velocity sensitivity** (Doppler f_d = (v/c)·f_carrier, ~20× L-band → carrier-phase radial velocity
   to sub-mm/s; transverse velocity still comes from the LOS sweep over the inbound-to-apogee arc),
   **ranging bandwidth** (range to ~30–50 m is ample for 140 m position at decent GDOP — easy headroom),
   and **antenna gain** (a 1 m Ka dish ≈ 48 dBi vs a ~6 m L-band dish) to help close the ~34 dB-harder
   link. Laser is rejected for the multilateration role: it is **point-to-point** (wrong shape for one
   signal seen simultaneously by a whole train + the rockets), needs pointing/acquisition of a fast small
   target, and serializes the fix; it stays the **terminal-phase** technology (the 1064 nm beacon + σ_θ
   tracker), at most a high-precision adjunct to the steerable rockets. Carries forward ADR 0011 dec-7:
   **omni transponder on the PuffSat, gain on the infrastructure side** (multilateration needs
   simultaneous multi-node visibility).

2. **PuffSat: one-way passive receive + authenticated broadcast.** The constellation cryptographically
   signs the nav message (Galileo-OSNMA style); the PuffSat verifies it with a **sub-gram, verify-only
   ASIC** and solves clock bias from ≥4 members — **no transmitter**. This is the right shape: the
   corrector runs onboard and needs only the PuffSat's *own* state, which passive reception delivers, and
   it defeats spoofing without two-way. The **downlink closes** (~35 dB-Hz with a modest 1 m / 10 W Ka
   dish at 150,000 km, ample for carrier-phase velocity over the slow coast); the uplink at 1 W omni is
   marginal (~25 dB-Hz), an additional reason to stay passive.

3. **Co-flyer and target rockets: the full two-way crypto-ns transponder.** The rockets are not
   mass/power constrained, so they carry the heavy end: a key-determined, calibrated **nanosecond-exact
   turnaround delay** gives cryptographic **distance-bounding** (anti-relay / anti-spoof — a party
   without the key cannot forge a valid reply) plus round-trip range, and lets the constellation/ground
   track them. **Synergy with ADR 0019 Lever 2:** the constellation pins the rockets' *absolute* states
   in the apogee regime, so the rocket→target vector the co-flyer tracker needs is known at *high*
   altitude, not only via GNSS at low altitude.

4. **Anti-jam is geometry + processing gain; crypto is anti-spoof.** Crypto authentication kills
   spoofing (one-way) and relay (two-way distance-bounding). **Jamming** (brute-force noise) is defeated
   by spread-spectrum processing gain (30–50 dB) + directional infrastructure antennas + the deep-space
   150,000 km regime itself — a jammer faces the same path loss into a spread, directional, two-way link.
   So jam risk is **not excessive here — because** it is the apogee/deep-space regime, not despite it.
   (1 ns ↔ 30 cm of range, so the turnaround delay must be *stable* to ~ns for ranging and is
   ground-calibratable; a constant delay cancels in Doppler.)

5. **Accuracy: match, don't beat — target σ_Tvel ≈ 0.66 mm/s / ~140 m at apogee (= C1).** After
   ADR 0019 the binding lever moved to the *terminal* grade, and fusion already gives terminal capture
   **100 % at 4.2× margin** (σ 0.21 m vs ≤1.65 m), with propellant under the 2 % line. Tightening the
   *entry* therefore pushes on a non-binding constraint (§16.6/§16.10: don't add capability that flips no
   verdict). The only payoff of a tighter entry is **substituting for the fusion hedge** — letting a bare
   10 µrad detector pass — and that needs the entry ~**3× tighter (~45 m / ~0.2 mm/s)**, a costlier
   constellation than the fusion it would replace. So **match 140 m** (the amount needed to keep the
   baseline funnel working, the real reason Lever 3 exists), and value the constellation for **snapshot
   GDOP at apogee** (robustness vs accumulating LOS diversity over an arc) and for **pinning the rockets**
   — not for a tighter number. Better-than-140 m from good GDOP is welcome **free, independent margin**
   (a second hedge on the entry×noise mode, in the spirit of ADR 0019's "two independent ways"); it is
   not a requirement to engineer toward.

**Mass/power pin (the PuffSat side).** The crypto+timing ASIC is sub-gram silicon at milliwatts — never
the mass driver. The drivers are the **oscillator** (TCXO ~1–3 g if solving clock-bias from ≥4 members,
else a CSAC ~35 g) and, *only if two-way*, the **PA**. A passive Ka receiver (patch antenna + front-end
+ correlator/verify ASIC + TCXO) is **~15–50 g, <1 W → ~0.2 % of the 25 kg bus**. **Mass is a non-issue;
transmit power is the real constraint**, which is exactly why the PuffSat stays one-way and the rockets
carry the two-way transponder.

## Considered options

- **Reuse L-band GNSS technology unchanged.** Rejected: lower Doppler/velocity sensitivity on the
  binding axis, no atmosphere to justify staying low, and larger antennas for the same gain.
- **Laser multilateration for the constellation.** Rejected: point-to-point and pointing-limited;
  it serializes a fix that must be seen simultaneously by the whole train + the rockets. Kept as the
  terminal-phase technology and a possible rocket-side precision adjunct.
- **Two-way crypto transponder on the PuffSat too.** Rejected for the baseline: the uplink is marginal
  (~25 dB-Hz at 1 W omni) and costs transmit power for a capability the corrector does not need (it needs
  its own state → passive receive suffices). Kept for the unconstrained rockets.
- **Spec the constellation tighter than 140 m.** Rejected: redundant margin (terminal capture already
  4.2× with fusion, propellant under budget); the substitution case needs ~3× tighter — a worse trade
  than the fusion architecture it would delete.
- **Lean on crypto for anti-jam.** Rejected as the mechanism: crypto is anti-spoof; jam immunity is
  processing gain + directional antennas + the deep-space regime.

## Consequences

- **Lever 3 graduates from "option" (ADR 0019 dec 3) to a specced architecture:** Ka-band authenticated
  broadcast, PuffSat passive one-way, rockets two-way crypto-ns transponder, **match-140 accuracy**.
- **ADR 0019 Lever 2 is strengthened:** the constellation pins the co-flyer's and target's absolute
  states in the apogee regime, extending the independently-known rocket→target vector above the GPS
  ceiling.
- **ADR 0011/0012 coordinator nodes:** the permanent constellation is their generalization; built, it
  replaces them for coast nav with *snapshot* GDOP at apogee (the C1 "min coordinator nodes / LOS
  diversity over arc" question becomes a constellation member-count + ring-vs-shell GDOP sweep).
- **Build steps (not in this ADR):** a pure sizing module (`apogee_nav.py`: link budget + Doppler/velocity
  budget + transponder mass/power, in the `tracker_budget.py` style, no JVM) and a GDOP / minimum-member
  sweep would quantify the architecture. The constellation's *existence* is assumed here (per scope);
  its deployment/cost is out of scope, deferred with the GMAT-class cross-checks toward a later rung.
- **Docs:** CONTEXT gains *Apogee nav constellation*, *Authenticated broadcast*, *Secure transponder /
  distance-bounding*; design §13 gains the Lever-3 signal spec.
- **Untouched:** A/B, C0–C2a, C3, D1, and the ADR 0019 results. The match-not-beat spec keeps D1.1's
  entry budget unchanged — this is additive infrastructure that hardens the coast-nav assumption, not a
  re-key of any measured result.
