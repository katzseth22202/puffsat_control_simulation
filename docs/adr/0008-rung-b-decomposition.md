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
