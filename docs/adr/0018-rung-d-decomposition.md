# Rung D decomposition: a C-baseline feasibility gate (D1) split from the MPC-value question (D2), behind hardware-requirement and truth-validation gates

**Status:** accepted

## Context

Rung D is the Monte Carlo that produces the design's feasibility verdict (§10; headline
**P(capture)**). A grill (2026-06-13, after C4 closed the C-rung) examined how to break it
into steps, make it performant for N=10³–10⁴, reach a clear yes/no, name what is still
missing for a *true* answer, and place LinCov/NEES. Three reframings drove the decomposition:

- **The feasibility yes/no is a property of the C baseline, not of MPC.** The whole C-rung
  was built on a fixed, transparent control law precisely so a miss is attributable to
  *knowledge quality*, not controller cleverness (§16.6). §16.10 already says MPC must *beat*
  that baseline on the *same* Monte Carlo — which presupposes a C-baseline MC exists. So the
  verdict can — and should — be reached on the C baseline first; MPC is a separate value
  question.
- **The MC gives a *conditional* verdict** — feasible *given* the nav/actuator specs (the
  10 µrad terminal tracker, the C1 nav Σ, the 1°/s slew). Those are sim *inputs*, not sim
  *outputs*; the strongest honest statement converts assumed specs into **derived
  requirements** ("feasible, and here is what each subsystem must achieve").
- **The headline is a tail probability**, and the binding physics lives in the tail —
  the catch-radius cliff (C3b: capture-grade to 500 m, 95 m miss at 600 m), actuator
  saturation, the significance-gate noise rectification, the lognormal coefficient skew, the
  MCC-2 firing threshold. All break the linear/Gaussian superposition LinCov assumes.

## Decision

1. **Split Rung D into D1 (feasibility gate) and D2 (MPC value).** D1 is the full
   closed-loop MC on the **C baseline** — the A1/A3 corrector + C3b ZEM terminal + C3c
   **MCC-2** trim + finite burn — and *is* the yes/no. D2 prototypes MPC and measures it
   against the C baseline on the **same** MC; MPC earns its place only if D1 shows a
   threshold/constraint violation or only-marginal capture (§16.10 a/b). **The feasibility
   verdict does not wait on MPC.**

2. **Three gates precede D1** — two tighten the conditional, one protects confidence:
   - **σ_θ tracker budget (blocking).** A pure `tracker_budget.py` (no JVM): derive what the
     10 µrad terminal grade demands — aperture / exposure / residual jitter / SNR for a dim,
     fast target on a shaking bus — **and acquisition** (tracker FOV vs the hand-off delivery
     Σ). Converts the load-bearing terminal-nav assumption from a guess into a derived
     requirement; if it is unmeetable, the catch radius (and the verdict) falls, so it blocks.
   - **Torque-margin back-of-envelope.** Confirm the ≥1°/s slew rail the C3b loop's noise
     discipline (the 45° firing-lag hold) rides.
   - **Truth-validation gate.** Tier 1 (energy/angular-momentum conservation +
     tolerance-halving on the Orekit nominal coast) + Tier 2 (an *independent* Python
     conservative-force Cowell cross-check of the coast — the coast-dominated 99 % where a
     truth-model bug would show). The full-force **GMAT** cross-check is **Rung F**
     (deferred), run as a **headless batch script → report → compare**, not via the
     CPython-version-fragile Python API and not through conda.

3. **Train mode + swept correlation pins.** D1.0 extends `DispersionSpec` /
   `sample_run_inputs` with the ADR 0016 shared-vs-per-unit split. The correlation inputs
   ADR 0016 named as paper-side pins — coefficient bias/spread ratio, deployer systematic,
   plane launch-window flexibility (the ±2 km **centroid retarget**) — are **swept axes**, not
   point values, so the verdict carries its own sensitivity. The §16.7 "multiplicative density
   factor" gap largely collapses here: per-unit density error ≈ the Cd·(A/m) draw to first
   order (drag ∝ ρ·Cd·A/m), and the *common* density component is one shared-axis pin.

4. **Nav Σ is a swept D1 axis parameterized by node count; report the minimum coordinator
   nodes.** GDOP is demoted from a gate to a *confirmation* that a realizable geometry lands
   inside D1's feasible Σ-region (ADR 0012 kept node geometry a derived requirement, never an
   assumed constellation, so there is nothing concrete to gate on). "Minimum nodes" is set by
   the **LOS diversity accumulated over the coast arc** (range + Doppler integrated by the
   filter), not snapshot multilateration count — so it can be smaller than the ≥4 a
   single-epoch range-only fix would need.

5. **Nav error is injected from the sampled C1 Σ, not a live UKF** (ADR 0012,
   requirements-by-covariance). **NEES is the upstream C1 gate** that earned the right to
   sample from that Σ (it caught the third-body-tide q error) — it is *not* a Rung-D sizing
   tool. Live-UKF spot-checks re-enter only if a nav-marginal tail forces them.

6. **Performance: parallelism + a cheaper corrector + tail variance reduction.**
   - **Process-level parallelism** over run indices, reusing the resume sink (each worker its
     own Python+JVM; a crashed worker costs nothing).
   - Replace the per-run **FD-Jacobian Newton** (≈40–60 descents/run) with a **Φ-Jacobian
     (the C0 STM) warm-started quasi-Newton**, warm-started from the A3 nominal correction
     (≈2–3 descents/run) — **with FD-Newton fallback** on the nonlinear tail runs (near the
     A1 authority boundary), where it matters most.
   - Resolve the P(capture) tail by **importance sampling / subset simulation (B)** on the
     Cr / nav / storm drivers, **validated by a brute-force batch (A)** that confirms the
     reweighting is unbiased. **LinCov never replaces the tail MC**; it serves as the
     IS-proposal designer, the control variate (tightening the Gaussian core), and a
     pre-screen (is the core comfortably inside the catch radius before spending core-hours?).

7. **D1 deliverables (the verdict surface):** headline **P(capture)** about the train
   centroid; **centroid-drift** distribution vs the ±2 km retarget; **scatter** about the
   centroid vs the plate; **propellant** vs <2 %; **perigee** diagnostic (low = good);
   **minimum node count**; and per-axis sensitivities (nav Σ, σ_θ, train-correlation
   fraction). A pass reads: *"feasible with a dumb, transparent law, given the [derived]
   nav/actuator specs."*

## Considered options

- **Bundle MPC into the feasibility demonstration** — rejected: it conflates knowledge
  quality with controller cleverness (the §16.6 logic the C-rung exists to preserve), delays
  the yes/no, and yields a weaker claim (feasible-*if*-MPC-is-clever vs feasible-with-a-dumb-law).
- **GDOP as a D1 gate** — rejected: no concrete constellation exists to test; sweeping Σ in
  D1 carries the sensitivity and demotes GDOP to a confirmation.
- **LinCov replacing the MC** — rejected (ADR 0012): a tail probability plus
  saturation/gate/lognormal nonlinearities break superposition. LinCov screens and
  accelerates; it does not replace the tail.
- **Brute-force-only tail** — retained as the *validation batch*, rejected as the *primary*:
  ~10⁴ for ~10 % tail error is affordable but wasteful when IS reaches the same precision
  10–100× cheaper.
- **GMAT via the Python API / conda** — rejected: GMAT is not a conda/pip package, the
  bundled API is CPython-version-fragile, and a one-shot cross-check wants loose coupling
  (headless script → report → compare). Hence Rung F, not a D1 dependency.

## Consequences

- **Pure-side:** `DispersionSpec` / `sample_run_inputs` gain the shared-vs-per-unit
  structure; a new pure `tracker_budget.py`; the Tier-1/2 truth-validation checks; the IS
  estimator + control variate (all unit-testable without a JVM).
- **JVM-side:** a Rung-D `runs/` slice strings the C-rung pieces into one closed-loop run;
  the Φ-Jacobian quasi-Newton corrector and the parallel worker harness.
- **Docs:** design-doc §13 queue gains the D1/D2 + gates breakdown; §10 train-mode note is
  already present; CONTEXT gains **Rung D (D1 / D2)** and **Tracker budget**.
- **Deferred rungs after D:** **E** cylinder shape (ADR 0009); **F** GMAT full-force
  cross-check.
- **The verdict is explicitly conditional.** The σ_θ budget (plus the GDOP/torque
  confirmations) is what tightens it toward "feasible, *and* here is what each subsystem must
  achieve" — the strongest statement a sim makes without a bench.

## Implementation findings — σ_θ tracker budget (2026-06-13)

The first pre-D gate, built pure (`puffsat_sim/tracker_budget.py`, no JVM — angular precision
is a focal-plane question, not an orbit one; like C4 it has no `runs/` glue). It is a four-term
RSS error budget for a declared `TrackerHardware` point plus the acquisition geometry.

- **GATE PASS, and it even meets the 5 µrad target.** The conservative default point (5 cm
  aperture, 1 ms exposure, 1 W laser beacon @ 1064 nm, beam ±2 mrad, η 0.3, nav-grade gyro,
  bench-calibratable 3 µrad focal-plane distortion) achieves **σ_θ = 3.2 µrad RSS — 3.1×
  under the 10 µrad requirement**. So the load-bearing terminal-nav grade is a *derived*
  hardware requirement now, not an assumed input, and D1 is unblocked on this gate.
- **The budget is calibration/jitter-limited, not photon-limited.** The active beacon gives
  SNR ≈ 1670 at the 300 km design (worst, longest) range, so photon-limited centroiding
  contributes only ~0.01 µrad; the RSS is dominated by the **focal-plane distortion floor**
  (3 µrad), with the post-impact smear residual (0.87 µrad, after differential cancellation)
  and the gyro bridge (0.58 µrad) minor. This dissolves the "dim, fast target on a shaking
  bus" worry: making the target an *active* beacon converts it to a bright source, and the
  residual limits are bench-calibratable, not fundamental.
- **It closes the loop back to capture.** Homing floor `2σ_θ²v²/a_max` at the achieved grade
  is **0.15 m ≪ 1.65 m**; the bare 10 µrad requirement reproduces ADR 0015's thin-margin
  1.45 m reference — so the achievable hardware sits well clear of the criterion that drives
  the catch radius.
- **Acquisition is governed by reference-star availability, not the delivery dispersion.**
  The ±2 mrad beam covers the **±1.4 mrad** (3σ · 141 m C1 lateral / 300 km) acquisition cone;
  the **binding FOV is the ±5.8 mrad** needed for 3 reference stars (~10th-mag density), which
  a **~1100-px detector at 10.6 µrad/px** resolves to Nyquist. The narrow-FOV-vs-star-count
  tension a naive acquisition-only sizing would miss is surfaced and comfortably met.
- A coarse distortion floor (≥ ~10 µrad) flips the gate to FAIL — the blocking semantics work,
  and the gate's `meets_requirement` / `meets_target` reads are the D1 entry condition.

## Implementation findings — torque margin (2026-06-13)

The second pre-gate, built pure (`puffsat_sim/torque_margin.py`, no JVM — attitude agility is
an inertia/actuator question, not an orbit one; like the σ_θ budget it has no `runs/` glue).
It is **non-blocking** (a confirmation, not a D1 entry condition): §13 flagged the ~10× slack
over the perigee sweep as "a result to confirm, not assume," and this is that confirmation.

- **CONFIRMED.** The thrust-direction demand is the perigee LOS rate `v_p/r_p` = **0.097 °/s**,
  which sits **10.3× inside the 1 °/s C3b direction-loop rail** (`anti_drag.PEAK_SLEW_LIMIT_DEG_S`);
  the B3a *measured* descent demand (0.048 °/s) carries 21×. The rail the terminal aim rides is
  not the binding constraint.
- **The rail is deliverable on conservative pins.** A whole-body case (I ≈ 5.06 kg·m² at a 0.45 m
  gyradius, a 50 mN·m cold-gas couple = two 50 mN thrusters at 0.5 m) gives α ≈ 0.57 °/s² — it
  reaches the demand rate in **0.17 s** (inside the 1 s control period) and the full rail in 1.8 s,
  and out-torques the aero disturbance (peak drag 16.7 mN × 0.15 m CP–CM offset = 2.5 mN·m) **20×**.
  Gimballing the small nozzle beats whole-body slew a fortiori.
- **Break-even on the unpinned margins** (the inertia and actuator are paper-side pins, like the
  target inertia in ADR 0015): the actuator still reaches demand in a period up to **I ≈ 29 kg·m²**,
  and still holds against drag down to a **2.5 mN·m** control torque — wide headroom on both.

## Implementation findings — truth validation (2026-06-13)

The third pre-gate, **non-blocking**. The pure core is `puffsat_sim/truth_validation.py` (the
invariants, the independent propagation, the reductions); unlike the other two pre-gates it does
have a `runs/` glue (`runs/truth_validation.py`) because the coast it validates is flown by
Orekit. Rung D's verdict is only as trustworthy as the truth it rides on, and the coast is ~99 %
of the trajectory, so a frame / μ / J2-sign / leaking-integrator bug would hide there.

- **VALIDATED on the reference apogee→800 km coast (~32 h).** **Tier 1** — a *numerical* two-body
  coast (the analytic Kepler route bypassed) conserves specific energy to **5.3e-15** and |h| to
  **6.5e-16** (machine precision — no integrator leak); tolerance-halving (×0.1 `rel_tol`) moves
  the trajectory **0 m** (the step is max-step-bound in the benign apogee region, so the
  conservation drift is the stronger read). **Tier 2** — an independent pure-Python RK4 Cowell
  (`estimation.two_body_j2_flow`, sharing only the pinned constants) matches the Orekit J2 coast
  to **15.7 m** (1.0e-7 of the orbit scale), confirming the frame / μ / J2 / force-assembly setup
  in the dominant dynamics.
- **Scope.** Tier 1 validates integrator health on conservative dynamics; Tier 2 validates the
  perturbed-dynamics *setup*. The non-conservative forces (drag / SRP / third-body) stay validated
  by the Rung-A force-signature tests. The full-force **GMAT** cross-check remains **Rung F**
  (headless batch, not a conda/API dependency) — this gate is the in-repo confirmation, not that.
- The `rel_tol` override added to `propagator._build_numerical_propagator` is the only production
  change: it both enables tolerance-halving and forces a numerical two-body propagation (bypassing
  the `is_keplerian` analytic route) so the conservation check exercises the real integrator.

## Implementation findings — train mode (D1.0, 2026-06-13)

The first D1 sub-slice (decision 3), built **pure** (`puffsat_sim/train.py`, no JVM) and
**TDD**, faithful to the project's "pure core first, JVM glue next" rhythm: the JVM closed-loop
train ensemble (stringing the C3b ZEM + C3c MCC-2 pieces into `run_record`) is **D1.1**.

- **The shared-vs-per-unit split is a generalization of the existing sampler, not a rewrite.**
  `TrainDispersionSpec` splits each `DispersionSpec` σ into a per-train **shared** draw
  (coefficient *bias*, F10.7/Ap *drivers*, deployer *systematic*) and a per-PuffSat **per-unit**
  draw (coefficient *spread*, injection *scatter*). `sample_train` composes the **same**
  `RunInputs` the JVM `run_record` already consumes, so D1.1 wires in with no record change.
  Two limiting cases pin the decomposition: per-unit σ → 0 makes a train internally identical
  (pure common mode), shared σ → 0 recovers independent per-unit draws (the single-PuffSat
  behaviour). The bias and spread compose **multiplicatively**, so the marginal per-unit
  log-variance is `s_bias² + s_spread²` (a tested invariant) — and the §16.7 density gap closes
  as predicted: the common density component *is* the shared F10.7/Ap driver, the per-unit
  density error folds into the per-unit `Cd·(A/m)` spread, no separate axis.
- **Standalone replay survives the train structure.** Arity-distinguished spawn keys — shared on
  `(train_index,)`, unit on `(train_index, unit_index)` — let `replay_train_unit` reconstruct any
  single PuffSat bit-for-bit without the train (§14.2), the train analog of `replay_inputs`. The
  flat `run_index = train_index·n_units + unit_index` keeps the resume sink unchanged.
- **The reduction makes ADR 0016's split operational.** `summarize_train_capture` reuses the
  `guidance` plate-capture machinery: the train **centroid** is the common-mode shift (charged to
  the ±2 km centroid retarget, `retarget_ok`), and each unit is **re-centered** on the centroid
  (where the plane aims) before judging capture, so the **scatter** alone faces the plate
  (per-axis σ vs 1.65 m, `scatter_sigma_ok`). A spot-check confirms the headline: a 1.5 km
  common-mode shift + ~1 m per-unit scatter reads `retarget_ok` with **100 % capture about the
  centroid** vs **0 % absolute** — the retarget earns its keep, exactly the ADR 0016 story. The
  correlation pins (bias/spread split, systematic, retarget capability) are `TrainDispersionSpec`
  fields, so D1.1's sweep is "run ensembles across specs," carrying its own sensitivity.

## Implementation findings — D1.1 closed-loop train ensemble (2026-06-14)

The first *runnable* Rung-D result: the C-baseline closed loop flown per PuffSat
(`runs/train.py`), reduced over a train into P(capture). **Architecture (user-confirmed):
Φ-composed entry + flown terminal.** The midcourse→hand-off residual is linear in nav/coefficient
error (C0's Φ), so each unit's hand-off lateral entry offset is *sampled* from the characterized
C0/C1/C2a budget (per-unit 141 m C1 nav, shared 149 m C2a Cr-prior, as 2-D lateral magnitudes);
the C3b ZEM terminal loop is *flown* through the cliff / gate / σ_θ·R noise. The corrector-in-loop
brute-force validation, MCC-2 scheduling, and the node-count Σ sweep stay later D1.x (decisions 4/6).

- **THE FINDING — the combined entry×noise stress is the binding effect, and it tightens the
  terminal-nav requirement from the 10 µrad *ceiling* to the ~3 µrad *target*.** C3b measured the
  two stresses *separately* (entry offsets noiseless → null to mm; tracker noise at zero entry →
  1.07 m at 10 µrad) and deferred their combination to Rung D. D1.1 measures it. The capture cliff
  vs grade (16-unit train, ~155 m mean entry): **2 µrad → σ 0.65 m / 100 %; 3.2 µrad → σ 1.35 m /
  ~94 %; 5 µrad → σ 2.54 m / 81 %; 10 µrad → σ 5.77 m / 31 %.** The nominal **10 µrad requirement
  ceiling FAILS** the σ ≤ 1.65 m criterion under combined stress (σ 5.77 m); the **achievable
  3.2 µrad grade PASSES** (σ 1.35 m < 1.65 m → capture-grade by the Rayleigh criterion; the 16-unit
  empirical 94 %/100 % about-centroid/absolute is small-N). So the σ_θ budget gate's "3.1× margin
  under 10 µrad" is **load-bearing, not slack** — D1.1 shows the achievable grade is *required*.
- **Mechanism (why the combination is worse than either alone).** The hand-off→target range is
  **2603 km**, so σ_θ·R is **26 m at 10 µrad** early in the descent. C3b's significance gate
  protects the *zero-entry* case by staying silent through that noisy large-R regime (|ZEM| ≈
  noise < 3σ gate), acting only late when R — and the noise — is small. A real ~150 m entry
  **defeats that protection**: it forces the loop to fire hard early, through the 26 m noise, which
  rectifies through the slew-limited gimbal into the trajectory. A tighter grade shrinks the
  large-R noise and recovers — hence the clean cliff.
- **The terminal funnel absorbs the common-mode entry, so the centroid retarget is unstressed.**
  Because every unit homes to the same aim point, the shared (149 m) entry offset is *nulled by the
  funnel*, not left as an arrival bias: the measured arrival **centroid drift is ~0–4 m** (≪ the
  ±2 km retarget) and the per-unit misses are unbiased zero-mean noise residuals. So the ±2 km
  retarget (a *pre-launch* plane-positioning mechanism) is the backstop for common-mode the funnel
  *cannot* remove (a common ToA bias, or common-mode beyond the ~450 m funnel authority) — neither
  of which arises here. The D1.0 shared→centroid / per-unit→scatter mapping holds for the *delivery*
  budget; the in-flight funnel then removes the common-mode delivery on top.
- **The other deliverables pass comfortably.** Propellant: worst-unit mission Δv ≈ 4.4 m/s
  (midcourse 2.19 + terminal aim ≤ 2.3) → **~0.9 % @ Isp 50 s**, well under 2 %. Perigee **~65 km**
  (deorbit-good). ToA scatter **≤ 0.7 ms** — two orders inside the 10 ms window.
- **Verdict (conditional, on the C baseline — no MPC).** D1 is **feasible on the dumb, transparent
  C law, conditional on the ~3 µrad terminal-nav grade** (achievable per the σ_θ budget gate), not
  the 10 µrad ceiling. D2 (MPC) is **not triggered** — the baseline meets the criterion at the
  achievable grade. The honest D1.1 caveat: the entry σ is the C1/C2a *crossing* budget applied at
  the hand-off (a conservative proxy; the true hand-off-lateral Φ is a D1.x refinement), and the
  tail P(capture) wants the importance-sampling batch (decision 6), not a 16-unit empirical.
