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
