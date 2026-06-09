# Rung B decomposition: impulsive corrector + finite-execution layer, feedforward terminal burn as a cost baseline, regime-switch first

**Status:** accepted

## Context

Rung B replaces Rung A's impulsive Δv with a finite-burn actuator (ADR 0004 settled
the actuator *model* — proportional cold-gas, 400 mN, Isp sweep, ~1 °/s slew, mass
depletion). A grilling session (2026-06-09) decomposed Rung B into testable
end-to-end slices and settled three decisions a future reader would otherwise find
surprising. A profile during that session found the *terminal step-cap*, not lack of
parallelism, was the A3 integration-test slowness.

## Decision

1. **Impulsive corrector + finite execution layer — `predict ≠ execute` at B, by
   design.** The A1 differential corrector stays impulsive (it solves a commanded Δv);
   a new **Actuator** layer maps each Δv to a finite, mass-depleting burn. At apogee a
   ~seconds burn ≈ impulsive (the orbit is days), so this is sound — and it deliberately
   breaks the A-rung `predict == execute` identity (ADR 0003): the corrector *predicts*
   impulsive, truth *executes* finite, and the residual interception miss **is** the
   actuator-realism erosion B1 measures. A finite-burn-aware targeter is deferred (needed
   only if the erosion proves larger than tolerance, or for the terminal phase at C/D).

2. **Feedforward terminal anti-drag burn as a deliberate *cost baseline*, not a closed
   loop.** Under Rung B's perfect state the apogee corrector already nulls the aim
   including the *known* drag (A3: controllable everywhere, miss < 2 m), so a terminal
   anti-drag burn is propellant-*suboptimal for aim*. It is kept at B anyway as a
   feedforward **measurement**: instrument the known-drag descent 600 → 200 km to get the
   anti-drag Δv (the bulk of the propellant budget — §13 coarse/fine) and confirm the
   400 mN / ~5 min / ~1 °/s actuator suffices. This is the perfect-knowledge cost +
   feasibility floor that Rung C/D's *feedback* rejection — where drag is uncertain and
   the terminal burn earns its aim role — is later measured against. A closed terminal
   loop at B (cancel drag, re-target the apogee burn to the drag-free arc) is deferred to
   C/D; at perfect state it is suboptimal and overlaps MPC's job.

3. **Regime-switch the terminal propagation first (B0); defer multiprocessing to Rung
   D.** Profiling (single full-force descent): the 30 s terminal step-cap — a §6.2 interim
   guard against overshoot-below-surface on *sparse* configs — is a **~5–6× tax** per
   descent on `full_force` with no accuracy gain (identical perigee/ToA at 30 vs 120–600 s
   caps), and the **~40 s JVM cold-start**, not parallelism, bounds small-run wall-clock.
   So B0 does the §6.2 regime-switch (coast on the big adaptive step → altitude-event
   hand-off → fixed-step Cowell terminal phase), which also gives B3's continuous burn its
   stable terminal phase. Process-level multiprocessing stays a **Rung-D** tool
   (§11.2/§12): it only parallelizes the multi-scenario axis, and its small-N payoff is
   capped by the per-worker JVM cold-start.

Slices (detail in §13): **B0** regime-switch → **B1** finite actuator (reproduces the A1
null + closes the propellant ledger; measures the erosion) → **B3** feedforward terminal
anti-drag (Δv cost + 400 mN / slew feasibility) → **B2** propellant fraction-vs-Isp curve
against the <2%-of-25 kg line.

## Considered options

- **Finite-burn-aware targeter at B** — deferred: heavier, and finite ≈ impulsive at
  apogee, so measure the erosion with the impulsive corrector first and add a finite-aware
  targeter only if it bites.
- **Defer the terminal burn entirely to C/D** — rejected: B should report the
  perfect-knowledge anti-drag propellant cost (the bulk of the budget) and actuator
  feasibility *now*, as the baseline C/D improves on.
- **Closed terminal feedforward loop at B** — deferred to C/D: drag-uncertainty's job; at
  perfect state it is propellant-suboptimal and overlaps MPC.
- **Pull multiprocessing up to Rung B** — rejected: not what was slow. The cap fix is
  ~5–6× on *everything* (single tests included); the JVM cold-start caps small-N parallel
  gains; §11.2/§12 already place process-level parallelism at Rung D.
- **Interim cap-loosen instead of the full regime-switch** — rejected: B3's continuous
  burn wants the fixed-step Cowell terminal phase anyway, so do the regime-switch once.

## Consequences

- `predict ≠ execute` at B is intentional; the B1 residual is a *measured erosion* (a
  finding, like A2/ADR 0006), not a bug. If it exceeds tolerance, a finite-aware corrector
  becomes a B follow-up.
- The B2 propellant deliverable is **dominated by the terminal anti-drag Δv**, not the
  corrections — the <2% test is largely an anti-drag-cost question (expect the 50 s Isp
  anchor to fail it ~3×, per ADR 0004). B2 reuses the ADR 0004 linear model
  (`fraction = Δv/(Isp·g₀)`) / the A3 `budget_dv_m_s` helper.
- B0 changes the propagation architecture (regime-switch); nominal and dispersed runs
  share it so the interception miss stays common-mode (as the 30 s cap did).
- Dead-time, MIB / cluster geometry, cosine losses, the coarse/fine chemical+cold-gas
  split, and the cm-aim 100 Hz trim remain deferred (ADR 0004).

## Implementation notes — B0 (2026-06-09): regime-switch hand-off, confirmed by a closed-loop falsification

Decision 3 said "do the §6.2 regime-switch once" and rejected "interim cap-loosen." We
measured B0 before building it (the A2/A3 discipline). The experiment came in two acts; both
are recorded here because the first nearly led us to the wrong design.

**Act 1 — the open-loop probe argued for a simple cap-loosen.** A probe swept the single
global adaptive step-cap and an altitude hand-off across the dispersion, **including the
joint low-drag tail** (cd, F10.7, Ap each at −4σ of the §13 log-normals — where the adaptive
step near 200 km stays largest and the sub-surface-overshoot risk peaks; the earlier profile
only saw the *high*-drag corner, which self-limits):

- The tax is in the **coast**, not the terminal phase: the interim 30 s cap throttles the
  long apogee→800 km coast (which wants ≥300 s steps); the terminal leg is already step-
  limited by drag.
- **No-cap is unsafe.** A 600 s (un-capped) descent is bit-identical at nominal but throws
  `OrekitException: point is inside ellipsoid` at −2σ and below — the integrator oversteps
  the 200 km event below the surface where the atmosphere model is undefined.
- A **global 300 s cap** was bit-identical (dP 0.000 m), ~2.5× faster, safe with ~2× margin
  to the measured sub-surface cliff (~575–600 s on the −4σ corner), and *faster than* the
  hand-off (whose small-step terminal + second propagator setup costs more **open-loop**).

On open-loop evidence alone the cap-loosen looked best and the hand-off looked dominated.
**That conclusion was wrong, and the open-loop probe could not see why.**

**Act 2 — the full integration suite falsified the cap-loosen.** With the 300 s cap the two
**closed-loop** corrector tests failed with the same `point is inside ellipsoid` blow-up
(garbage state `a = −1.16e8 m, e = 0.73`). Cause: the corrector probes *large re-phasing Δv*
(e.g. the rejected ~88 m/s root on an uncorrectable tail run), and those probe orbits have
wildly varying perigee; their terminal descent oversteps at 300 s where the 30 s cap kept
every integration stage above ground. The open-loop capstone never drives the corrector, so
it never saw this. **A single global cap cannot be both fast in the long coast and safe in
the stiff terminal phase for the orbits the corrector explores** — which is precisely the
§6.2 prediction that "a single pure-adaptive scheme chokes on terminal discontinuities."

**Resolution — the regime-switch hand-off (decision 3 vindicated, with one refinement).**
B0 ships the §6.2 hand-off: coast on the big adaptive step (600 s) → **800 km** altitude
event → terminal phase on a tight cap (30 s). The terminal cap equals the old global 30 s
cap, so it is *provably* as safe as the prior code (identical terminal behaviour) for every
orbit, corrector probes included; the coast runs at 600 s, recovering the tax. Verified:
the closed-loop tests pass, the open-loop capstone is byte-identical (common-mode preserved)
and faster (N=50: 4m18s → 3m20s wall, CPU 2m48s → 0m47s). **Refinement vs decision 3:** the
terminal phase here stays *adaptive* (a capped DP853), not fixed-step Cowell. A fixed-step
terminal is deferred to **B3**, where the continuous anti-drag burn — needing a deterministic
cadence aligned to the control clock — is its first consumer. So B0 builds the coast/terminal
*seam* (the §6.2 architecture) and B3 swaps the terminal integrator at it.

**Lesson recorded:** measure the *closed-loop* path, not just the open-loop capstone, before
settling a propagation-architecture decision — the corrector's Δv probing is a stressor the
dispersion alone does not reproduce.

### Forward-compatibility with the 5 cm terminal aim (cm-trim, deferred)

The eventual goal is to center on the pusher plate to ~5 cm, not the current ~2 m / km-scale
miss. Does that reopen B0? **Mostly no, and the coast never.** Split the two senses of
"precision":

- **Numerical accuracy** is *already* ~mm end-to-end (relTol 1e-10), coast and terminal alike
  — the step-cap only ever forces steps *smaller* than the controller wants, never larger
  (bit-identical results confirm it). 5 cm is far inside current numerical accuracy
  everywhere, so accuracy alone reopens nothing.
- **Control cadence** is the real driver of the 5 cm aim, and it lives **below 800 km**: the
  100 Hz fine-trim loop (cold-gas/electric, §13 / ADR 0004) and the fixed-step deterministic
  stepping (§6.2) that lands integration nodes on the control clock. That is exactly the
  fixed-step Cowell terminal phase deferred to B3/C/D, and the terminal step there shrinks
  from 30 s toward the ~10 ms control period (or dense-output sampling off the integrator).

So the user's reasoning is sound: the extra precision is a terminal-phase concern; the coast
above 800 km only needs to deliver the right *ballpark* for the terminal aim and needs no
revisit. The **800 km hand-off boundary we built is the seam where the cm-trim machinery
plugs in** — the entire fine-aim region (drag onset ~800 km → interception 200 km) sits below
it. The regime-switch *architecture* is therefore forward-compatible with 5 cm; only the
terminal *integrator* changes (adaptive 30 s → fixed-step ~100 Hz), and that change is already
on the roadmap (B3 builds it; the cm-trim fine stage is the deferred §13 / ADR 0004 item).
