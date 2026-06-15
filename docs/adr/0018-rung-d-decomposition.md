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

## Implementation findings — D1.x corrector-in-loop validation (2026-06-14)

The first D1.1 caveat, closed. D1.1 *sampled* each unit's hand-off entry from the linear C0/C1/C2a
budget instead of running the midcourse corrector; this is decision 6's **brute-force validation
batch (A)** that the sampled-entry shortcut is unbiased. Scope (user-confirmed): the **nav leg**
only — the Cr-prior predict/execute mismatch leg and the Φ-Jacobian quasi-Newton speedup stay
separate follow-ons. The pure reduction is `puffsat_sim/corrector_validation.py`; the JVM glue
`runs/corrector_validation.py` flies the **real corrector** (the C0 path — `run_record` with
`report_controller` + a predict-side nav offset, nominal coefficients/zero injection) over a batch
of **combined** nav draws sampled from the C1 nominal-cell Σ, records the interception miss, and
coasts the corrected execute state to the 800 km hand-off to measure the true lateral entry. Φ is a
minimal C0 sweep and Σ the pure C1 LinCov nominal cell — both reused. The σ-consistency band is
sample-size-aware (`max(0.15, 3/√(2N))`), so a smoke batch is judged honestly (a flat tolerance
false-fails small N on RMS sampling noise — caught when the n=8 smoke read "inconsistent" at a
19.5 % deviation that was < 1σ of sampling error).

- **VALIDATED (N=64).** Three results triangulate D1.1's keystone:
  1. **Linearity / superposition.** Per-draw `|miss − Φ·δ| = 0.01 %` of `|miss|`. C0 proved this
     *one axis at a time*; here all six axes are perturbed together at realistic C1 magnitudes and
     the residual is still ~0 — **no cross-terms**, so superposing the C0/C1/C2a budget legs (what
     D1.1's sampled entry does) is sound. Bias **14 m**, within the sample-mean noise of zero.
  2. **Magnitude (end-to-end).** Measured crossing-miss σ **147 m** ≈ `ΦΣΦᵀ` **141.3 m** ≈ D1.1's
     **141 m** (`ENTRY_LATERAL_PERUNIT_M`) proxy. The real corrector reproduces the per-unit entry
     magnitude D1.1 fed the terminal loop — the 141 m number is confirmed with the corrector in the
     loop, not just by the C0 sensitivity.
  3. **Hand-off conservatism.** The *actual* 800 km hand-off displacement is **68.9 m** —
     **2.13× smaller** than the 141 m crossing miss D1.1 fed the loop as the entry. The terminal
     loop must still null the fully-developed crossing miss (~141 m) regardless of where it starts,
     so front-loading that as a static hand-off offset makes the loop face it earlier, through the
     larger early-R (2603 km → 26 m at 10 µrad) noise — the binding D1.1 mechanism. So D1.1's entry
     proxy **over-stresses** the loop: the ~3 µrad grade it derived is **conservative** (its verdict
     is pessimistic, the safe direction).
- **Net.** D1.1's "Φ-composed sampled entry + flown terminal" architecture is validated for the nav
  leg: the corrector residual is linear and unbiased at the combined realistic draw, the fed entry
  magnitude is right, and the crossing-budget proxy is conservative at the hand-off. D1's conditional
  feasibility verdict stands and is, if anything, pessimistic on the entry stress.
- **Remaining D1.x:** the Cr-prior predict/execute mismatch leg (the shared 149 m entry — needs a
  `RunVariant` coefficient-mismatch knob); the **Φ-Jacobian warm-started quasi-Newton corrector**
  (this brute-force batch is its validation reference); nav-Σ-by-node-count → min nodes; MCC-2
  scheduling; the importance-sampling tail.

## Implementation findings — D1.x importance-sampling tail P(capture) (2026-06-15)

The last standing D1.1 caveat, closed (decision 6). D1.1's headline P(capture) was a 16-unit
empirical ("94 % about-centroid / 100 % absolute") that cannot resolve the figure of merit. This
slice resolves the per-unit capture-failure tail by importance sampling, validated against a
brute-force batch. The pure reduction is `puffsat_sim/tail_capture.py` (the IS proposal over the
2-D per-unit hand-off entry offset, the exact likelihood-ratio weights, the weighted estimator with
CI / effective sample size, and the brute-force-agreement check); the JVM glue
`runs/tail_capture.py` flies the **same C3b ZEM loop D1.1 uses** over an IS batch and a brute-force
batch at the achievable 3.2 µrad grade — tracker noise sampled fresh per trajectory (so it stays
out of the importance weight), truth physics held at nominal (the feedforward then fully rejects
drag, isolating the binding entry×noise tail driver).

- **TAIL-RESOLVED per-unit P(capture) ≈ 99.2 %** [95 % CI 98.4–99.98], N = 500 brute force at the
  3.2 µrad grade. This replaces D1.1's unresolved 16-unit "100 %" with a real number: the
  capture-failure is a **shallow ~0.8 % event**, not a deep 10⁻³ rare event. The point estimate
  clears the paper's 99 % target; the *lower* CI bound (98.4 %) sits just under it — N = 500-limited
  (4 escape events), so a larger batch tightens it (the slice's one remaining refinement).
- **The flown tail is 2.0× heavier than Gaussian — the key finding.** The arrival scatter σ is
  **1.51 m, which MEETS the ≤ 1.65 m criterion** (ADR 0015), yet the measured plate-escape (0.80 %)
  is **2.0× the 0.41 % that σ predicts under the Gaussian/Rayleigh tail the criterion assumes.** The
  catch-radius cliff + significance-gate noise rectification fatten the tail beyond Gaussian —
  exactly the nonlinearity decision 6 said the LinCov/Gaussian screen cannot see. So the σ-criterion
  is *mildly optimistic* about the tail (the absolute miss stays tiny, but "σ ≤ 1.65 m" over-states
  the capture it implies). Tail-resolution, not a covariance screen, is what surfaces this.
- **IS validated, and aggressive inflation is counterproductive.** Because the tail is shallow,
  brute force resolves it directly and IS is the validated cross-check: at a *gentle* κ = 1.35 the
  proposal is well-behaved (ESS 240 / 300, Σw/N = 1.00) and IS agrees with brute force at the
  shallow r_val = 3 m (14.8 % [9.6, 19.9] vs 14.6 % [11.5, 17.7] — CIs overlap, the reweighting is
  unbiased). An aggressive κ = 2.5 *inflated* the variance instead of reducing it: it pushed draws
  past the ~450 m funnel-authority edge into the saturation-catastrophe regime (a single 499 m miss)
  that nominal entries never reach, blowing up the weighted variance (rel error 1.0 at the plate).
  The funnel cliff makes naive variance-scaling fragile; the right designer choice oversamples the
  escape ring without crossing the funnel edge. IS is the tool reserved for the **deeper** tails a
  tighter plate / stricter target / nav-marginal grade would create, where brute force gets
  expensive; at 3.2 µrad / 5 m plate the tail is shallow enough that brute force suffices.
- **Knowledge-limited, confirmed in the tail.** The capture-failure tail is driven by the
  entry×noise rectification (nav), not control authority: the funnel catch radius (~500 m) binds
  only in the catastrophe regime nominal entries never reach. No fat dispersion tail eats the funnel
  margin — the knowledge-vs-authority conclusion holds where it is most likely to break.
- **Verdict.** D1 stays feasible on the dumb C baseline at the achievable 3.2 µrad grade (~99.2 %
  per-unit capture, the rest deorbiting by §9 as intended); the ADR 0019 fused grades (5-array σ
  0.58 m, +co-flyer σ 0.21 m) drive the tail far below the 99 % target with margin. **D2 (MPC) is
  not triggered.** The σ-criterion's 2.0× tail optimism is the slice's caution for the marginal
  single-detector grade; the fused architecture absorbs it.
- **Remaining D1.x:** the Cr-prior predict/execute mismatch leg; the Φ-Jacobian quasi-Newton
  corrector; nav-Σ-by-node-count → min nodes; MCC-2 scheduling. (Cr/storm IS is a 2nd-order
  extension per D1.1; a larger BF batch tightens the P(capture) lower bound.)

## Implementation findings — D1.x fused-grade tail: entry-limited, recorded not coded (2026-06-15)

The IS-tail slice above resolved the **single-detector** (3.2 µrad) tail. The natural follow-on —
what is the resolved tail at the **ADR 0019 fused grades** (array-only 1.62 µrad, +co-flyer
0.76 µrad), whose "far below target" was only inferred from a 16-unit run — was *measured* (a
throwaway driving the same `runs/tail_capture` flight path: a deterministic entry-magnitude sweep to
locate the catch radius plus an N=250 IS batch at κ=2.5). Recorded as a finding, **not coded**: the
result is that the fused tail leaves the regime `tail_capture.py` was built for, so a faithful
runner would need an inverted estimator (below) — deferred until a fused-grade headline is actually
needed.

- **The fused tail is entry/authority-limited, not noise-limited.** At the fused grades the arrival
  scatter is σ ≈ 0.30 m, so a capture *failure by noise* is a ~16 σ event that never happens (brute
  force saw 0/250 plate escapes). Every escape comes from the **entry exceeding the catch radius**: a
  deterministic entry sweep puts a **hard cliff at ~475 m** (entry 450 m → 0.38 m miss; 475 m →
  14.6 m; 500 m → 40.7 m; 700 m → 249.6 m — the saturated funnel passes the overage through almost
  un-attenuated), and the cliff is **identical at 1.62 and 0.76 µrad**. The catch radius is set by
  actuator authority (the ½·a·t² funnel, ~487 m in C3c), not by the tracker grade.
- **Resolved P(capture) ≈ 99.999 %** (escape 1.41e-5 [6.6e-6, 2.15e-5], relerr 0.27, ESS 74 at
  N=250), cross-validated by the semi-analytic catch-radius estimate `P(entry > 475 m)` = 1.18e-5
  (2-D Rayleigh, per-axis σ 100 m) — they agree. This is the regime IS was *reserved* for (the
  IS-tail finding above): κ = 2.5 was counterproductive at the shallow 3.2 µrad tail but is **correct
  here**, because the > 475 m region *is* the tail, and brute force is blind to it (0/250). The roles
  invert — IS becomes the headline estimator and the catch-radius analytic the validator, since BF
  cannot resolve *or* validate at this depth with feasible N.
- **The co-flyer is redundant for the capture *number*.** Array-only (1.62 µrad) and array+co-flyer
  (0.76 µrad) give the **same** entry-limited tail (same 475 m cliff, same ~1.2e-5 escape), because
  both already drowned the noise. The co-flyer's value is **robustness margin** (it holds capture-grade
  if the 3 µrad distortion floor proves optimistic), not capture probability.
- **Terminal-noise reduction has a knee at ~the array grade.** Past it, capture is governed by the
  midcourse **entry budget** (entries exceeding the funnel) and the **catch radius** (authority), not
  the tracker. Further σ_θ reduction buys nothing for capture — to go lower you raise authority
  (earlier hand-off / more thrust) or shrink the entry (midcourse nav / Cr-prior), per ADR 0021.
- **The 99.999 % is itself conservative** — it rides D1.1's 141 m entry proxy; at the
  corrector-validated *actual* hand-off displacement (68.9 m, 2.13× smaller; the corrector-validation
  finding above), `P(entry > 475 m)` ≈ 1e-21 and the fused tail vanishes. So 99.999 % is a floor, and
  the fused architecture clears the 99 % target by orders.
- **If a coded fused-grade headline is ever wanted:** extend `TailCaptureFinding` to pick the
  *resolving* estimator (BF when escapes are seen — the shallow regime; IS when BF is blind but IS
  resolves — this deep regime) and add a catch-radius validation path (entry sweep → `R_catch` →
  2-D-Rayleigh) since the r_val BF-vs-IS check degenerates (both 0) once the whole tail sits past the
  funnel edge. A `runs/tail_capture` entry point over the fused architectures would then mirror
  `runs/train.fused_train_rerun_report`.

## Implementation findings — D1.x Cr-prior mismatch leg, recorded not coded (2026-06-15)

The other deferred leg of decision 6's brute-force batch A (the corrector-validation above did the
nav leg). It flies the **real corrector planning with a prior Cr** and executing against a different
**truth Cr** (a throwaway over the existing flight primitives — `report_controller` solving against a
`physics_prior` crossing closure, then `descend` under `physics_truth`; the harness shares one physics
between predict and execute, so the predict-side coefficient mismatch is the one knob it lacks, done
out-of-band). The Cr is the SRP `cr_area_over_mass`; the mismatch sweep is the factor δ.

- **C2a validated end-to-end.** The crossing miss is **exactly linear** in δ and **injection-decoupled**
  (zero and non-zero injection give the identical slope — clean superposition), and in C2a's (T, N)
  lateral metric it reproduces the analytic **745 m/factor exactly** (149 m at the 0.2 prior). So
  C2a's `coefficient_sensitivity` (Φ · ∂Δv/∂Cr) holds against the real nonlinear corrector — the Cr
  twin of the nav-leg superposition result.
- **But the 149 m is 98.6 % along-track.** At the near-perigee crossing (e ≈ 0.92) the velocity is
  transverse, so C2a's T axis is ~along v: the δ = 0.2 miss decomposes R/T/N = 0/149/7 with along-v
  147 m. The **true plate-frame ⊥v Cr entry is only 22 m** (110 m/factor) — **6.8× smaller** than the
  (T, N) budget figure, exactly the along-track dominance of SRP that ADR 0021 records. C2a's (T, N)
  lateral is a conservative proxy that folds the along-track ToA component into "lateral."
- **The Cr-prior leg is comprehensively benign and does not feed the 475 m cliff.** It is
  **shared/common-mode** (a prior bias shifts every unit the same way): the 22 m cross-track part is
  absorbed by the ±2 km centroid retarget (90×), and the 147 m ≈ 13.6 ms along-track part by the ±5 s
  launch-window slip the same retarget commands (CONTEXT). The train split already books Cr as the
  149 m **shared** leg (`ENTRY_LATERAL_SHARED_M`), absorbed by the homing per D1.1's centroid-drift ~0;
  this leg confirms it end-to-end and adds that even the *raw* cross-track Cr is 6.8× below the (T, N)
  figure. The per-unit cliff entry (141 m) stays nav-driven; the per-unit Cr *spread* (a few % factor)
  contributes single-digit metres of ⊥v, negligible against the nav leg.

## D1 closeout — feasibility verdict (2026-06-15)

With the three pre-gates passed, the architecture decided + sized, the terminal tail resolved, and
both entry legs validated against the real corrector, the D1 feasibility gate is **closed**.

- **Verdict: D1 is FEASIBLE on the dumb C baseline (A1/A3 corrector + C3b ZEM terminal + C3c MCC-2 +
  finite burn), conditional on the fused terminal-nav grade — which the σ_θ tracker-budget gate
  (3.2 µrad per detector) and ADR 0019 fusion deliver. D2 (MPC) is NOT triggered** (no measured D1
  violation; §16.10's "MPC must beat the baseline on the same MC" has no baseline failure to beat).
- **Architecture (Seth, 2026-06-15):** zero separate coordinator nodes — coast/apogee-state nav is the
  150k Ka-band **apogee nav constellation** (ADR 0020, which generalizes the coordinator node; min 3
  shell / 4 ring members, sized to *match* the C1 grade = the 140 m entry nav leg); terminal relative
  homing is the fused **co-flyer + target array** (ADR 0019).
- **Capture story:** single-detector 3.2 µrad → **99.2 %** (noise-limited, marginal, 2× heavy tail);
  fused array (≥ ~1.6 µrad) → **99.999 %** (entry/authority-limited; the co-flyer is robustness margin,
  redundant for the number). The binding regime flips at the **knee** (arrival σ ≈ 1.05 m): below it,
  capture is governed by the midcourse entry budget vs the 475 m catch radius, not σ_θ.
- **Entry budget, both legs validated against the real corrector:** nav per-unit 141 m (the
  cliff-relevant scatter — corrector-validation: Φ-map reproduced, hand-off 2.13× conservative vs the
  crossing proxy); Cr-prior shared 149 m (common-mode, absorbed by the ±2 km / ±5 s retarget; 98 %
  along-track ToA, true ⊥v 22 m). Drag is feedforward-solved and non-binding (ADR 0021).
- **Margins, all comfortable:** catch radius 475 m vs entry RSS ~224 m ≈ 2.1×; propellant < 2 % at the
  C3b worst stack; ToA ≤ 0.7 ms vs 10 ms; perigee ~65 km (intended debris disposal, §9). Conservatisms
  stack in the safe direction (entry crossing proxy 2.1×, Cr cross-track 6.8×, fused tail rides the
  141 m proxy so the real tail ≈ 1e-21).
- **Remaining D1.x are non-blocking refinements, not gates:** the Φ-Jacobian warm-started quasi-Newton
  corrector (performance; this batch + the nav batch are its reference), MCC-2 scheduling, and a larger
  brute-force batch to tighten the single-detector lower bound. None changes the verdict; each only
  sharpens a number or speeds the MC. **Rung D's open question now is D2 (MPC value), which is not
  triggered** — so the substantive build ladder is complete at the C baseline.

## Post-closeout refinement — the entry budget is along-track-dominated (2026-06-15)

> **Superseded by the re-run below (same day).** The decomposition in this section is correct (the
> entry *is* along-track-dominated), but its **inferences** — that the 141 m proxy over-stresses the
> ⊥v cliff ~6×, that the requirement is over-conservative, and that a looser grade / single 10 µrad
> detector might pass — were **falsified** by the flown re-run. The along-track is *not* irrelevant to
> capture: it is range-observable and trips the significance gate, which is what lets the loop null
> the small ⊥v. See "Re-run executed" below. Read this section as the hypothesis, not the verdict.

Prompted by Seth's observation that the Cr leg's miss was nearly all along-track: if the *plate*-frame
capture miss is ⊥v_rel (ADR 0015) and the along-track axis is ToA, then the (T, N) "lateral" the whole
budget chain (C0/C1/C2a) uses — which excludes only radial — may be **counting along-track as
cross-track**. Measured (a throwaway flying the real corrector over N=32 C1-Σ nav draws, decomposing
the crossing miss into true ⊥v_rel vs along-v): **it is, and badly.**

- **The nav leg is along-track-dominated.** (T, N) RMS 173 m (≈ the 141 m budget metric, N=32
  sampling), but the **true ⊥v_rel cross-track entry is only 24.6 m** (14 %); along-v (ToA) is 172 m,
  and the per-axis split is R/T/N = 0 / 173 / 3 m. R = 0 *exactly* — the crossing is an altitude-
  triggered event (fixed 200 km), so the miss lives in the horizontal (T, N) plane with T = downrange
  ≈ along-velocity. The 24.6 m ⊥v is the genuine plate-frame lateral (the ~8° flight-path angle tilts a
  sliver of the downrange offset into ⊥v).
- **Why** (the deeper point): the dominant nav uncertainty is *apogee transverse velocity* (C0's
  binding lever), and a transverse-velocity error changes the orbital *period* → it surfaces as
  **along-track phase / ToA**, not cross-track. So C0's "binding requirement" is effectively a *ToA*
  sensitivity that the (T, N) convention has been booking as cross-track. Both entry legs are
  along-track-dominated (nav ⊥v 24.6 m; Cr ⊥v 22 m).
- **Implication — the closeout's margins are conservative on the binding axis.** True cross-track entry
  ≈ √(24.6² + 22²) ≈ **33 m RSS** vs the 475 m catch radius — **~14× margin, not the 2.1×** the (T, N)
  RSS (224 m) implied. And D1.1 fed a **141 m ⊥v entry when the real cross-track ZEM is ~25 m**, so it
  **over-stressed the terminal loop ~6×**. Since D1.1's binding mechanism is the entry tripping the
  significance gate and forcing early firing through σ_θ·R noise, a 6× smaller entry means the gate
  holds far better — so the **~3 µrad requirement is likely over-conservative, and the single 10 µrad
  detector may pass**, which would make ADR 0019 fusion a *hedge*, not a necessity.
- **Caveats (don't oversell):** the 172 m along-track is not free — it is ToA (~16 ms uncontrolled) —
  but the terminal loop already nulls it (D1.1: ToA ≤ 0.7 ms) and the ±5 s launch-window slip absorbs
  the common-mode part; ToA is handled, just not by tight apogee nav. And the requirement relaxation is
  a **nonlinear entry×noise question** — the linear argument above is suggestive, not proof.
- **Suggested re-run (the confirmation):** re-run D1.1 (`runs/train.py`) with the corrected per-unit
  ⊥v entry magnitude (~25 m nav, ~33 m total RSS) in place of `ENTRY_LATERAL_PERUNIT_M = 141 m`, and
  read whether the **10 µrad single detector** now clears σ ≤ 1.65 m / capture. If it does, fusion
  (ADR 0019) drops from *required* to *margin* — a real architecture simplification — and the apogee
  constellation's coast-nav grade can likewise relax. This **refines** the closeout (D1 is *more*
  feasible, possibly at a looser grade), it does not overturn it.

## Re-run executed — the "over-stress 6×" hypothesis is FALSIFIED; the 141 m proxy is ≈accurate (2026-06-15)

The suggested re-run above was flown (`runs/train.py`, controlled same-seed sweeps). Its premise — that
the (T, N) 141 m entry over-stresses the ⊥v cliff ~6×, so D1.1 is conservative and a looser grade /
single detector might pass — **does not survive contact with the flown loop.** The model was **not**
changed; the committed 141 m-isotropic entry stands, now *validated* as a ≈accurate (mildly conservative)
proxy. The reasoning chain that produced the hypothesis missed one fact: **the C3b significance gate
keys on the full 3-D ZEM, and the along-track entry is range-observable.**

- **A ⊥v-only entry is the wrong experiment — it is pessimistic, not accurate.** Replacing the 141 m
  isotropic entry with a 24.6 m ⊥v-only entry makes capture *worse*, not better: at 3.2 µrad, pooled
  σ rises 1.47 → 2.14 m and P(capture) falls 99 % → 91 % (same-seed; consistent across trains). The
  10 µrad single detector still fails either way. So the naïve "shrink the entry to the true ⊥v" move
  is an artifact, not a correction.
- **Why: the significance gate needs a signal to trip.** `guidance.significant_zem` holds fire while
  |ZEM| ≤ 3σ of the σ_θ·R noise. A *lone* small ⊥v entry (24.6 m) sits at the gate threshold (SNR ≈ 3
  at the single-detector grade), so the loop fires late with little authority → larger residual. The
  **real** entry is anisotropic — ~172 m along-v (which is *along the closing LOS, so range-observable
  to ~1 m*) + ~24.6 m ⊥v — and the large along-v makes |ZEM| large, **tripping the gate early**; the
  loop then homes on the *full* ZEM and nulls the small ⊥v cleanly.
- **Confirmed by the anisotropic re-run.** Restoring the along-v term (172 m along-v + 24.6 m ⊥v)
  recovers capture: pooled σ 2.14 → **1.35 m**, P(capture) 91 % → **100 %** at 3.2 µrad — matching
  D1.1's original 141 m-isotropic result (σ 1.47 m, ~99 %). The committed 141 m-isotropic entry trips
  the gate via a large (mis-attributed) ⊥v exactly as the real anisotropic entry trips it via the
  along-v, so it lands within ~9 % of the accurate answer and on the **conservative** side.
- **Net: D1.1's conclusions stand unchanged.** Requirement ~3.2 µrad (single detector marginal at the
  gate-tripped operating point; **fusion the robust baseline** — 5-array 1.62 µrad → σ 0.88 m / 100 %,
  +co-flyer 0.76 µrad → σ 0.32 m / 100 %, pooled). The "over-stressed 6× / over-conservative / single
  10 µrad may pass / fusion becomes a hedge" conjecture in the section above is **withdrawn**: the
  along-track is *not* irrelevant to capture (it carries the gate), the requirement does not relax, and
  no hardware loosens.
- **A reusable lesson on the entry model.** The plate miss is ⊥v, but the ⊥v *capture* is coupled to
  the along-v through the significance gate — so the entry cannot be reduced to its ⊥v magnitude. The
  (T, N) figures across the budget chain (`tracker_budget` acquisition FOV, `coeff_requirement` /
  `authority` RSS ledger, `corrector_validation` crossing proxy) therefore remain the right ones to
  carry. The genuinely accurate model is the *anisotropic* entry (along-v + ⊥v); implementing it
  (vs. the validated 141 m-isotropic proxy) is an optional D1.x refinement, not a correction —
  the outcome is ~identical.
