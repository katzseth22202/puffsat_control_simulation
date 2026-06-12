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

## Implementation findings — C3a (2026-06-11)

C3a is built and measured (`puffsat_sim/terminal.py` pure planner + report,
`build_fixed_step_propagator_from_orbit` in `propagator.py`,
`run_terminal_feedforward` / `terminal_feedforward_report` in `montecarlo.py`;
integration test `tests/integration/test_terminal_feedforward.py`).  Nominal run
(1 Hz control clock, 1 s fixed step, ~25 s wall):

- **Equivalence pin (decision 4): 5.5 mm, ToA < 1 µs.**  The fixed-step Cowell
  (classical RK4, Cartesian, 1 s) reproduces the proven adaptive-30 s unburned descent
  from the 800 km hand-off to the crossing at the integrator-floor scale — five orders
  under the 400 m working radius.  The terminal integrator swap is settled.
- **Drag displacement at the crossing is 8.5 cm — ~20× under this ADR's ~1–2 m
  estimate.**  The naive bound (0.015 m/s × ~100 s) overcounts: drag is concentrated in
  the final seconds before the crossing, where the accumulated Δv has almost no time to
  integrate into position.  Terminal drag rejection is a *centimeters*-at-the-crossing
  problem, which further de-stresses the C3b feedback's share of the 1 m endpoint.
- **Executed residual 2 mm → rejection ~45×.**  The open-loop ZOH burn (180 one-second
  commands, drag known) cancels the displacement to below the equivalence pin itself.
  B3a's measured-only profile survives execution end-to-end: ZOH quantization,
  command-boundary alignment, and the 800→600 km uncompensated band are all sub-pin.
- **Plan ≈ B3a's profile: Δv 0.0145 m/s (B3a trapezoid: 0.015), peak thrust 15.96 mN
  (B3a: 16.7), peak slew 0.048 °/s.**  Both ADR 0004 gates PASS on the *executed* burn
  (B3a only measured the demand); the small deficits are the hold-at-tick ZOH reading a
  steeply rising profile, exactly the expected direction.  Propellant 0.0030 % of wet
  mass at the conservative Isp 50 anchor — invisible against the 2 % line.
- **Mass-depletion stays the B1 sentinel-Isp convention.**  ADR 0008 deferred real
  depletion to "B3's large anti-drag burn"; B3a/C3a falsified the premise (the burn is
  ~1.5 g, 0.006 % of wet mass), so constant-mass execution + pure Tsiolkovsky transform
  remains exact to ~6e-5.
- Planner detail that mattered: the descent's final sub-period step carries a
  disproportionate share of the impulse (exponential drag rise), so the plan holds the
  last command to the end of the span instead of truncating at the last full tick.

C3b inherits: the fixed-step terminal + maneuver-segment machinery as-is; the measured
8.5 cm/2 mm scale as the drag-feedforward floor under the ZEM loop; the 1 Hz cadence as
the validated center of the {0.1, 1, 10} Hz sweep.

## Implementation findings — C3b (2026-06-12)

C3b is built and measured (`puffsat_sim/guidance.py` pure ZEM/noise/plate-frame core,
`puffsat_sim/runs/guidance.py` tick-by-tick closed loop + one-axis sweep; integration
test `tests/integration/test_run_guidance.py`).  Sweep: ~60 closed-loop descents from
the 800 km hand-off over one shared context, seeds fixed, 1 Hz nominal clock; the aim
point is the drag-free nominal crossing.

- **Measured catch radius: 500 m — the 400 m working re-baseline (decision 1) holds
  with margin.**  Noiseless lateral entry offsets land capture-grade through 500 m
  (residual 3 mm at 91 % tick saturation) and collapse at 600 m (95 m miss, 100 %
  saturated): a clean thrust-authority cliff at the Δv ceiling a_max·t_terminal ≈
  0.016 × 246 s ≈ 3.94 m/s, bracketing the ½·a_max·t² ≈ 480 m isotropic estimate.
  Inside the radius the residual *falls* with offset (0.136 m at zero entry → 3 mm at
  500 m): large entries keep the demand above the 5 mN proportional floor and the
  significance gate so the loop tracks continuously, while the zero-entry leftover is
  onboard-model bias (two-body+J2 predictor, re-predicted each tick) plus floor
  deadband — not authority.
- **The nominal 10 µrad tracker grade is capture-grade: RMS 1.07 m, max 2.70 m over 8
  seeds (requirement σ ≤ 1.65 m; closed-form floor 2σ_θ²v²/a_max = 1.45 m).**  50 µrad
  reads RMS 3.68 m / capture 88 % — fails the 99 % line, so the grade requirement
  genuinely binds between 10 and 50 µrad; ADR 0015's σ_θ ≤ 10 µrad spec is the
  validated knee.  The measured curve is much flatter than the σ² floor: fine grades
  bottom out on the model-bias/deadband floors (0.58 m at 2 µrad vs 0.06 m floor),
  coarse grades under-run the floor by refusing to act (3.7 m vs 36 m at 50 µrad),
  paying in capture tail instead.  The constant 1 m grade reads 0.53 m, consistent.
  ToA errors ≤ 0.06 ms RMS everywhere — two orders inside ADR 0015's 10 ms window;
  arrival timing is a non-issue at this rung.
- **The raw ZEM law cannot be closed on σ_θ·R knowledge — noise discipline is
  load-bearing, worth two orders of magnitude.**  Closing the textbook law directly on
  the noise envelope rectifies zero-mean knowledge error through the slew-limited
  single engine into ~150 m RMS at 10 µrad (double-integrator harness).  Three declared
  constants in `guidance.py` fix it: a 3σ significance gate that decides *whether* to
  act, never *how much* (soft-thresholding strands n²× the floor); a 35 s track window
  holding the gate open through the endgame; a 45° firing-lag hold (burning mid-slew is
  the rectifier itself).  The injected error is Gauss–Markov (τ = 10 s, the
  gyro-bridged track), not white per-tick draws no real filter would pass.
- **Dispersed drag is absorbed silently: Cd ×0.5 / ×2 / storm F10.7 = 250, Ap = 100
  land 0.09–0.19 m with Δv ≤ 0.04 m/s.**  C3a's cm-scale drag story survives feedback
  with a *wrong* feedforward; the 5 mN floor means most low-drag feedforward ticks
  never fire and the ZEM endgame mops up.  ADR 0013's ±6.7 coefficient tolerance is
  untouched at factor-2 truth error.
- **Cadence (decision 4): 1 Hz validated; the zero-offset axis is soft.**  0.1 and
  10 Hz both stay capture-grade at the nominal grade; 10 Hz is mildly worse (0.38 m —
  more ticks, more rectification chances), 0.1 Hz mildly better by barely firing.  At
  zero entry offset "do almost nothing" is already an 8.5 cm trajectory (C3a), so this
  axis mostly measures firing eagerness; C4's τ-sweep should ride entry offsets.
- **Gates (ADR 0004) PASS — by construction, on the rails.**  The closed loop *uses*
  both rails (peak thrust = the 400 mN cap on saturated entries; the gimbal rides the
  1 °/s slew rail), and the execution machinery — cap, slew limiter, firing-lag hold —
  is what keeps the command history physical; C3a's "demand under the gate" framing
  inverts here.  Pin detail: recomputing a rail-limited step angle through acos(dot)
  carries ~1e-12 of round-off, so the gate verdict takes a matching allowance.
- **Propellant: terminal aim is now the dominant ledger line at the radius edge — and
  the 2 % claim still holds.**  Δv vs entry: 0.026 m/s at zero (feedforward + bias
  trim), 2.40 at 400 m, 3.75 at the 500 m edge (0.77 % @Isp 50), ceiling 3.94.  A
  typical entry at the C2a budget RSS (224 m) costs ~1.3 m/s.  The worst-case stack
  (B2 mission ledger 2.19 + radius-edge 3.75 ≈ 5.9 m/s) reads ~1.2 % @Isp 50 — under
  the paper's 2 % at the conservative anchor, but B2's "aim is cheap" headline gains a
  thrust-limited asterisk: cheap *inside* the funnel, ceiling-priced at its edge.

C3c inherits: the 500 m measured radius as the saturation boundary, the ~670 m 3σ
upstream tail it does not cover (MCC-2 stays load-bearing, decision 5), and the kept
A2 solver for the cost curve.  Rung D inherits the measured radius and the σ_θ grade
as sampled axes.
