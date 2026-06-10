# C3 terminal guidance: thrust-limited ~400 m catch radius, ZEM law, relative-GNSS 1 m endpoint, high-node trim for tails

**Status:** accepted

## Context

C3 (design doc §13, ADR 0010) is the one closed-loop dynamics slice: the deferred §6.2
fixed-step Cowell terminal phase plus the terminal burn's aim role. A grilling session
(2026-06-10, after C2) found that the quantity every C-rung slice has been judged
against — the **§9 "few km" terminal catch radius — does not survive the ADR 0004
actuator**:

- §9 sized the terminal authority from the *propellant* side ("32 m/s eats ~5 km in
  ~300 s"), silently assuming the Δv is deliverable in the window. Delivering 32 m/s in
  300 s needs ~2.7 N; the actuator is **400 mN** → a_max = 0.016 m/s² at 25 kg.
- The 600→200 km descent is **~170–180 s** (vis-viva: radial speed ~3.0 km/s at 600 km
  falling to ~1.7 km/s at 200 km), not §9's "~5 minutes".
- Thrust-limited lateral authority ½·a·t²: **~240 m from 600 km, ~450 m from the 800 km
  hand-off, ~700 m from 1000 km**. Buying 5 km needs t ≈ 790 s — a burn starting
  ~3000+ km up, where there is no drag to reject (a midcourse trim wearing a
  terminal-burn costume). The binding constraint flipped from propellant to thrust.
- Meanwhile the upstream 1σ error budget (ADR 0013: RSS 224 m) rides exactly at the
  600-km authority, and its ~3σ tail (~670 m) exceeds **any** thrust-limited radius.

Two more grill findings reshaped the sensing and the goal:

- **Terminal drag rejection is a meters-scale problem.** B3a's measured anti-drag Δv is
  0.015 m/s over the descent; even *uncompensated*, drag displaces the crossing by only
  ~1–2 m (the Δv accrues near the bottom, with no lever arm left). A tactical-grade MEMS
  accelerometer (bias ~3–10 µg) sees the ~1700 µg peak signal at SNR of several hundred;
  used as ADR 0013 reassigned it (direct force feedforward, subtracting commanded
  thrust), the feedforward residual is sub-meter. The `Cd`-parameter prediction question
  C3 inherited **dies quantitatively**: a 100% `Cd` error moves the crossing ~1–2 m,
  below the nav floor.
- **The endpoint floor is set by navigation, and the miss is *relative*** (§16.3: the
  pusher plate is on the target). Absolute single-frequency GNSS floors at ~3–10 m;
  **differential GNSS against the target** (both vehicles fly receivers; the target
  broadcasts over the same link class the coordinator nodes already use) cancels the
  common-mode errors over the ≤100 km closing baseline: code-differential ~0.5–1 m,
  carrier-differential ~0.1–0.3 m (rendezvous practice, PRISMA-class). Control is not
  the binder: nulling 1 m needs only ~11 s of remaining burn, thruster granularity
  contributes ~mm, and predictor model bias (third-body ~10⁻⁶ m/s² near perigee) is
  suppressed by re-prediction at every control step.

## Decision

1. **The catch radius is re-baselined to a working ~400 m, and C3 measures the real
   curve.** The terminal aim burn starts at the **existing 800 km regime hand-off** (B0's
   event; 800→600 km is drag-negligible pure aim time), giving ~450 m of thrust-limited
   authority; 400 m leaves margin for slew, drag feedforward, and the fact that ½at² is
   an isotropic upper bound. C3b's deliverable includes the **residual-vs-entry-offset
   curve** whose knee is the *measured* catch radius, replacing the assumption. The
   C0/C1/C2 tables re-read at 400 m with honest margins (T-vel requirement 1.84 mm/s vs
   0.66 achieved → 2.8×; nav lateral 2.8×; Cr prior 2.7×; budget RSS 1.8×) — they were
   parameterized by radius precisely for this day; nothing reruns.

2. **Terminal sensing: MEMS force feedforward + GNSS feedback; coordinator nodes are the
   GNSS-free fallback variant only.** The accelerometer feeds the *measured*
   non-gravitational force forward (no `Cd` estimate in the loop, ADR 0013), subtracting
   commanded thrust (few-% thrust knowledge ≈ tens of µg ≈ meters at worst). GNSS sets
   the aim floor. No UKF runs in the loop — terminal nav quality enters as **injected
   state noise** (the requirements-by-covariance discipline, ADR 0010/0011).

3. **The endpoint goal is 1 m, conditioned on relative GNSS.** The terminal nav-noise
   grade σ_rel is a **swept parameter** — {absolute GNSS 10 m, code-differential 1 m,
   carrier-differential 0.2 m} — so "what 1 m costs in sensor architecture" is a
   *measured output* (the C1 range-only pattern). The 1 m claim rests on
   code-differential; carrier-differential is the reported upside (its high-dynamics
   carrier tracking at ~10.8 km/s is the stretch item, flagged not assumed). Timing is
   free (GNSS time ~ns; 1 m ↔ 0.1 ms along-track). §16.3/§16.4 gain the
   target-crosslink note. **The claim pair is: catch radius ~400 m (entry, thrust-limited)
   + endpoint floor ~1 m (nav-limited)** — different numbers doing different jobs.

4. **The dumb fixed law is ZEM/ZEV guidance, zero-order-hold at the control clock.** At
   each control step: predict the crossing miss under no-further-thrust (onboard
   two-body+J2 + drag feedforward), command `a = k·ZEM/t_go²` capped at a_max, hold over
   the step. Classical, one gain, optimal for the double integrator, TDD-able against
   closed forms — read ADR 0010's "PID/LQR" as this. Control cadence swept coarsely
   {0.1, 1, 10 Hz} (expectation per the C4 bandwidth logic: 1 Hz ample, 0.1 Hz nearly
   ties — which also settles the JPype-boundary cost question). The fixed-step Cowell
   terminal integrator steps at ≤ the control step and is **equivalence-pinned against
   the proven adaptive-30 s config on an unburned descent before any burn is added**.

5. **The km-scale tails are carried by a high-node impulsive trim (MCC-2), and C3c
   measures its cost curve now.** 3σ of the upstream budget (~670 m) exceeds any
   thrust-limited terminal radius, so ~0.3% of Rung-D runs saturate without a trim.
   A2's table locates the authority: ~1.7 m/s per km of along-track at 30,000 km, dead
   at ≤5,000 km — so §16.6's "correction 2 at 800–1000 km" is the wrong node for this
   job, while the *role* is exactly ADR 0006's "observable-drift correction, once
   estimation enters at Rung C/D" (tracking has revealed the −Φδ residual by
   mid-descent; impulsive beats continuous 2× for the same Δx). C3c measures Δv-per-km
   vs node altitude at the ~km scale, reusing the **kept** `solve_two_burn_correction`
   (ADR 0006). Trim *scheduling* (always-fire vs threshold-triggered) stays a Rung-D/MPC
   question.

6. **Slices and floors.** **C3a** (= B3b): fixed-step Cowell terminal + execute B3a's
   feedforward as a real ZOH burn — known drag, open-loop; executed residual + propellant
   + ADR 0004 gates. **C3b**: close the ZEM loop — entry offsets up to saturation,
   dispersed drag (`Cd` factor + f10.7/Ap), σ_rel injection; deliverables = the
   residual-vs-entry curve (measured catch radius), endpoint floor vs σ_rel, Δv, peak
   thrust/slew vs gates. **C3c**: the authority curve vs burn-start altitude + the MCC-2
   cost curve. The 1 m endpoint is measurable at the current `rel_tol=1e-10` (~cm
   crossing floor); the **5 cm plate-centering stays parked** (needs the flagged 1e-13 +
   Encke work and plate-relative sensing — not this rung).

## Considered options

- **Keep the 5 km catch radius and let MCC-2 deliver into it** — rejected. It leaves the
  terminal claim dishonest (the burn cannot absorb what the radius promises) and hides
  the thrust limit that actually binds; the radius parameter exists to be re-baselined.
- **Raise thrust above 400 mN** — rejected. The actuator is the paper's (ADR 0004);
  changing hardware to rescue an assumption inverts the project's direction of inference.
- **PID/LQR path-tracking as written in ADR 0010** — rejected in favor of ZEM. Tracking a
  reference path is the wrong shape for "null a predicted miss"; ZEM is the classical
  terminal law, simpler to verify, and stays within "dumb fixed law."
- **Run the C1 UKF inside the terminal loop** — rejected (the UKF-in-every-trajectory
  anti-pattern, ADR 0010); injected σ_rel noise carries the knowledge quality.
- **Absolute-GNSS-only terminal nav** — kept as the swept 10 m point, not the baseline;
  it cannot support the 1 m endpoint and the miss is relative anyway.
- **Defer MCC-2 wholly to Rung D** — rejected. It is load-bearing for the tails of the
  C-rung story and its cost curve is cheap to measure with kept code; only its
  *scheduling* defers.
- **Unpark the 5 cm centering into C3** — rejected. Different floor regime (integrator +
  plate-relative sensing); C3's honest floor is ~1 m.

## Consequences

- Design doc §9 gets the thrust-limited correction (the "few km" framing was
  propellant-side); §13 C3 reframed to the C3a/C3b/C3c slices with ZEM and the claim
  pair; §16.3/§16.4 gain the target-crosslink / relative-GNSS notes.
- The C0/C1/C2 requirement tables are re-read at the 400 m working radius (2–3× honest
  margins); the C1 finding "100 m is marginal" is unchanged and now visibly close to the
  working point — the third-body-in-filter upgrade path stays the named lever.
- New pure code planned: the ZEM law + t_go/predicted-miss machinery (TDD vs double
  integrator closed forms), the σ_rel injection, trim-cost post-processing over the kept
  A2 solver. JVM side: the fixed-step terminal configuration + ZOH burn segments in
  `montecarlo.py` (B0/B1 machinery generalized).
- C4's τ-sweep rides on C3b's loop as specced; the swept control cadence here feeds it.
- Rung D inherits: trim scheduling (threshold vs always), the measured catch radius as
  the saturation boundary in the MC, and σ_rel as a sampled axis.
