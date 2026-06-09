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

## Implementation notes — B1 (2026-06-09): finite-burn execution + the erosion finding

B1 lands decision 1 (impulsive corrector + finite execution layer). The A1 corrector is
**unchanged** — it still solves an impulsive Δv; a new layer fires that Δv as a finite Orekit
`ConstantThrustManeuver` at the apogee node, and the residual interception miss is the
measured `predict ≠ execute` erosion. The pure burn kinematics (Tsiolkovsky) live in
`puffsat_sim/actuator.py` (`Actuator` / `FiniteBurn` / `plan_burn`, unit-tested without a JVM);
the maneuver construction and the finite-execute path live in `montecarlo.py`.

**The mass convention (the one design fork, settled with the user).** The propagator runs at a
fictitious **1 kg** so the lumped `Cd·(A/m)` / `Cr·(A/m)` scale drag/SRP directly — but a burn's
`a = F/m` and Tsiolkovsky depletion want the real **~25 kg**. Two paths were on the table:
(A) keep the 1 kg propagator and scale the thrust to it (`F_eff = F·m_p/m_wet ≈ 16 mN`), so the
*acceleration* and therefore the burn duration are physically exact, with real propellant
computed by Tsiolkovsky in the pure module; (B) carry the real 25 kg and rescale the drag/SRP
cross-sections by mass. We took **A**: it leaves `build.py` and the byte-identical open-loop
capstone untouched, and — because the executed burn is constant-mass (a sentinel Isp makes the
in-propagation depletion ~0) — the trajectory is **Isp-free**, exactly the ADR 0004 decision-2
shape (Δv computed once, propellant swept across Isp at B2). Option B's only added fidelity is
in-propagation mass-depletion *coupled to drag*, which is negligible at B1's ~0.2 % propellant
fraction and is precisely what **B3**'s large drag-coupled descent burn needs — so it is deferred
there, mirroring B0's deferral of the fixed-step terminal to its real consumer.

**The finding — the erosion is real, along-track, and compensable.** Seed 20260608 run 0
(commanded Δv ≈ 2.17 m/s → a ~136 s burn) lands the finite-executed crossing **~89 m** off the
impulsive null, **88.8 m of it along-track** (the dr_p/dv_a lever), with the ToA shifting
~0.18 s. The cause is geometric, not numerical: the burn cannot be *centered* on the apogee node
(deployment is at apogee = the propagation start, so the burn fires forward), so its impulse
centroid lands ~68 s late, and the dr_p/dv_a lever amplifies that small phasing offset. It is far
above the 0.7 m impulsive residual (a genuine effect, not solver noise) yet small against the
km-scale open-loop dispersion — so at the current km aim the actuator-realism erosion is
negligible and **B1's question ("can the actuator deliver the correction?") is answered yes**.
But it scales with Δv (longer burn → larger offset) and is ~1780× the eventual 5 cm aim, so it is
the concrete trigger for the deferred **finite-burn-aware targeter** (which would simply shift the
burn so its impulse-equivalent centers on the node, nulling most of the along-track erosion). The
propellant ledger for that correction is 0.44 % / 0.32 % / 0.11 % of 25 kg at Isp ∈ {50, 70, 200} s.

**Why the existing record sufficed.** With the corrector nulling its *impulsive* prediction onto
the nominal aim, the `miss_rtn_m` of a *finite*-executed run **is** the erosion by construction —
so B1 added no `RunRecord` / sink field, and the propellant ledger is the pure `plan_burn`
transform over the already-recorded `total_dv_m_s` (B2 sweeps it). Slew (~1 °/s) stays unmodelled:
a single fixed-direction apogee burn needs no slew; it binds only at B3's turning anti-drag burn.

## Implementation notes — B3a (2026-06-09): split B3, defer the fixed-step terminal, measure the cost

Decision 3 (and the B0 note) assigned the §6.2 fixed-step Cowell terminal phase to B3, "the
continuous burn being its first consumer." Building B3 surfaced a sharper reading: B3's *headline*
— anti-drag Δv, peak thrust ≤ 400 mN, peak slew ≤ 1°/s, <2% propellant — is answerable by
**instrumenting the existing adaptive descent**, with *no executed burn and no fixed-step terminal*.
So B3 splits, by the same "defer infrastructure until its consumer exists" logic that shaped B0:

- **B3a (measure, done):** a step-handler over the ephemeris of the existing descent samples the
  truth drag acceleration through 600 → 200 km and reduces it (pure `anti_drag.summarize_anti_drag`)
  to the requirement. Drag is evaluated at the propagator's 1 kg, which *is* the real a_drag (the
  lumped `Cd·(A/m)` is the real coefficient, ADR 0009); peak thrust scales by the real 25 kg.
- **B3b (execute, deferred → C/D):** build the fixed-step Cowell terminal and fly the open-loop
  feedforward burn end-to-end. Its true first consumer is the *executed / closed-loop* rejection of
  *uncertain* drag at C/D, so the fixed-step terminal moves there. B3a answers feasibility without it.

**Finding (nominal descent, conservative cannonball `Cd·(A/m)=0.04`):** peak thrust **16.7 mN**
(~24× under 400 mN), peak slew **0.048 °/s** (~21× under 1 °/s, ≈ the paper's ~0.1 °/s sweep
estimate), anti-drag Δv **0.015 m/s** → **~1.5 g @ Isp 25 s (0.006 % of 25 kg)**. The paper's
`sec:estimate_cold_gas` 374 g / 400 mN is a *deliberately stacked-pessimistic* upper bound (×10
area, ×20 surge, ×2 v², ×2 pulsing on GOCE's solar-min baseline); the physical NRLMSISE requirement
on the fast eccentric pass is ~24× / ~250× below it. **Both gates pass with enormous margin and the
<2% claim holds a fortiori — even on the conservative coefficient, so no ADR 0009 grounding is
triggered.** A no-`RunRecord`-change measurement, like the B1 erosion: the profile is reported, not
stored. The cross-check confirms the paper's feasibility floor with ~20× headroom.

## Implementation notes — B2 (2026-06-09): the "anti-drag dominates / 50 s fails ~3×" expectation, falsified

The Consequences section predicted B2 would be "dominated by the terminal anti-drag Δv, not the
corrections ... expect the 50 s Isp anchor to fail it ~3×." **Both halves are wrong, measured.**
B3a found the anti-drag Δv is **0.015 m/s** (0.006 % of 25 kg at 50 s — negligible, not dominant),
and A3 found the corrections « budget. So the mission Δv is correction-dominated and small: a
representative correctable run (correction 2.17 m/s + anti-drag → **2.19 m/s**) costs **0.45 % /
0.32 % / 0.11 %** at Isp 50 / 70 / 200 — <2 % at every anchor, ~4.5× margin even at 50 s.

The pure ledger (`propellant.propellant_curve`, ADR 0004 linear `Δv/(Isp·g₀)`) is the exact inverse
of `sweep.budget_dv_m_s`, so the 2 % line at 50 s **is** the 9.8 m/s authority budget; A3's
"controllable everywhere « budget" therefore means every correctable run is <2 % at 50 s by
construction. Why the original expectation was wrong: it took the paper's stacked-pessimistic
anti-drag estimate (374 g, ~3.6 m/s-equivalent) as the physical cost, but B3a's NRLMSISE
instrumentation showed the real anti-drag is ~250× smaller. The full propellant *distribution* over
the dispersion is deferred to **Rung D**; B2 is the deterministic ledger it feeds through. (Aside:
the N=8 sizing ensemble converged 1/8 under the *untuned* A1 plain-Newton corrector — the A3
LM-tuned corrector is controllable everywhere; that convergence rate is the A1 authority/tuning
story, not the propellant story.) **Rung B is complete: B0 / B1 / B3a / B2; B3b → C/D.**
