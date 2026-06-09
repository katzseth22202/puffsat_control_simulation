# Rung C decomposition: estimation in the loop — requirements-by-covariance, not a giant ensemble; one knowledge source at a time

**Status:** accepted

## Context

Rung C drops the perfect-state assumption of Rungs A/B and asks whether the §8
controllability survives *realistic knowledge* — navigation error, model mismatch, and
latency. A grilling session (2026-06-09) decomposed it and settled several decisions a
future reader would otherwise find surprising — most importantly that Rung C is **not** a
giant Monte Carlo running a real UKF in every trajectory. Rung C populates the
predict/execute seam that already exists (ADR 0003): at A/B `predict == execute`
(full-force truth); at C the corrector's onboard `predict` finally diverges from truth.

## Decision

1. **Decompose by source of imperfect knowledge, one at a time** — C0 nav-error placeholder
   → C1 UKF state → C2 UKF coefficients → C3 closed-loop terminal → C4 latency — each
   measured against the Rung-B perfect-state baseline. Mirrors the §13 "add one new thing at
   a time" discipline so a failure is attributable to a single cause.

2. **Requirements-by-covariance, not a faithful-filter ensemble (the load-bearing compute
   decision).** The paper asks a *requirements* question — *how accurately must the PuffSat
   know its state / the coefficients?* — answered by parameterizing the estimation error as a
   (regime-dependent) covariance and sweeping it: cheap, deterministic, no filter. The real
   UKF is built only to confirm **feasibility** (its achieved covariance sits inside the
   required threshold) on a *few representative runs*. The giant N=10³–10⁴ ensemble — and even
   there, *sampling* nav error from the characterized covariance rather than re-running the
   filter — is **Rung D**. Running a UKF at 100 Hz × thousands of trajectories is the wrong
   tool for a requirements paper.

3. **C2's coefficient requirement is ~free from the A3 map.** Residual interception miss ≈
   (coefficient *estimation* error) × (A3 sensitivity); A3 already mapped required Δv « budget
   for ±2σ coefficient error, so the coefficient-knowledge *tolerance* is a post-processing of
   A3 with no new propagation. C2's only new work is the *feasibility* check (does the
   accelerometer-fed UKF hit that tolerance).

4. **C3 is the lone closed-loop *dynamics* slice, and it owns the fixed-step terminal.** The
   terminal drag disturbance is time-varying and the feedback *chases* it down the descent, so
   it cannot be collapsed into a static covariance — C3 genuinely spins up the §6.2 fixed-step
   Cowell terminal (B3b's deferred consumer) and the executed feedforward burn. **C3a (= B3b)**
   executes the feedforward under *known* drag; **C3b** closes the PID/LQR loop on *uncertain*
   drag (the §10.1 residual that C2 sizes), where the terminal burn earns its aim role. Still a
   few representative runs (nominal + worst-case drag), not an ensemble.

5. **Controller: deterministic corrector + dumb terminal feedback at C; MPC at D.** The
   midcourse keeps the A1/A3 corrector fed an *estimated* state (no replanning); the terminal
   adds a fixed PID/LQR magnitude loop. MPC — adaptive constrained replanning — is held for
   Rung D and measured against the C baseline. Same §16.6 logic that kept the corrector
   deterministic at Rung A: hold the control law fixed and transparent so a miss is
   attributable to the *knowledge quality* under study, not to controller cleverness.

6. **Latency is a lumped dead-time *per loop*, judged by bandwidth not sample rate.** Comms
   round-trip + coordinator compute on the slow outer loop; onboard sensor + valve on the fast
   inner loop. Dead-time erodes phase margin as `ω_c·τ`; the drag-rejection bandwidth is ~1 Hz
   even though the inner loop *samples* at 100 Hz, so a tens-of-ms τ costs single-digit degrees
   of phase. C4 depends on C3b (no loop to destabilize before the terminal feedback exists).

7. **Onboard fast loop + a compute-budget feasibility note (paper, not sim).** The fast
   drag-rejection loop is onboard (comms can't sit in a fast loop) and near-trivial (~10⁴
   flops/s, sub-gram MCU); the UKF (~10⁵–10⁶ flops/s on 8 states) and the Rung-D MPC (small QP
   at ~1 Hz) fit a ~Cortex-M7 or off-load to the coordinator. **50 g is ample for the classical
   first pass** (single-use → COTS, non-rad-hard parts). Neural augmentation is the deferred
   stressor. This is a paper sizing note alongside the propellant/actuator estimates, not a
   simulation step.

## Considered options

- **A giant UKF-in-the-loop Monte Carlo at C** — rejected: it answers the *distribution*
  question (Rung D), not the *requirements* question the paper asks, and is the compute blow-up
  the covariance sweep avoids.
- **One combined UKF slice (state + coefficients)** instead of split C1/C2 — rejected: position
  observability (ranging/GNSS) and coefficient observability (accelerometer/drag) are different
  measurement problems; splitting lets the loop close on nav error before the harder
  drag-estimation coupling, and lets C2 reuse A3 for its requirement.
- **MPC at Rung C** — rejected: confounds the knowledge-quality measurement with controller
  cleverness (the §16.6 argument); MPC enters at D against the C baseline.
- **Off-board fast loop** — rejected: a 100 Hz loop can't tolerate a ~13 ms comms round-trip;
  the fast loop is onboard, which is also what keeps comms off the fast loop.
- **Model latency by physical cause (comms / compute / actuation) as separate sim knobs** —
  rejected: the sim only needs the lumped τ per loop; the causes are a budget *table*, not
  separate dynamics parameters.

## Consequences

- Rung C is cheap: C0/C1/C2/C4 are requirement sweeps + A3 reuse + back-of-envelope, with
  pure-TDD building blocks (UKF math, measurement models, dead-time buffer); only C3 spins up a
  closed-loop sim, on a few runs.
- B3b and the §6.2 fixed-step Cowell terminal land in **C3a** (their real consumer); the
  terminal burn's *aim* role lands in C3b.
- The characterized nav-error covariance feeds **Rung D** (sample from it rather than re-running
  the filter every trajectory).
- TDD guidance is recorded per slice in §13: UKF math / measurement models / dead-time buffer /
  PID-LQR law / nav-injection / A3-reuse = **`/tdd`**; threshold sweeps and loop-closes findings
  = **measured** (not asserted), like the A2/B1/B3a findings.
- Deferred: the full distribution (D), MPC (D), the optimistic real-shape comparison (E), neural
  augmentation (second pass), the cm-aim terminal sub-problem.
