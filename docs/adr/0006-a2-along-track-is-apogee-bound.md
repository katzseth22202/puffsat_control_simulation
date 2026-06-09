# Rung A2 finding: along-track/phase error is apogee-bound; no useful mid-descent authority

**Status:** accepted — **supersedes [ADR 0005](0005-a2-two-burn-midcourse.md)**

## Context

ADR 0005 settled A2 as an apogee + **mid-descent (900 km)** two-burn minimum-Δv
least-norm corrector. Its premise (the A1 carry-forward): a single apogee impulse is
Δv-*inefficient* for along-track/phase error — a 2.4σ radial injection cost ~88 m/s to
null from apogee alone — and *two well-placed impulses, the second at a mid-descent node,
would null the same miss cheaply*. The phrase that drove the design was "the bulk of the
dispersion nulled cheaply" by a classical statistical-midcourse second burn.

Implementing A2 (the pure solver `solve_two_burn_correction` was committed at 1f4e7e3; a
JVM two-leg harness was prototyped) we tested that premise directly with a **node-altitude
sweep** on the exact hard run: `master_seed=20260608`, run 1 — the 2.4σ radial injection
whose ~28 km along-track miss A1 could not null sub-budget.

## Finding

**The presumption that a 900 km mid-descent burn is a game-changer is wrong.** A second
burn at mid-descent adds essentially *zero* along-track authority. The sweep (run 1, seed
20260608, `max_step=5 m/s`, `max_iter=12`):

| 2nd-burn node | converged | interception miss | total Δv | apogee / mid Δv |
|---|---|---|---|---|
| 900 km   | no   | 27.3 km | 83 m/s (capped) | 33 / 50 |
| 5 000 km  | no   | 24.5 km | 69 m/s (capped) | 9 / 59  |
| 30 000 km | yes  | 0.9 m   | 59.5 m/s        | 11 / 48 |
| 80 000 km | yes  | 0.0 m   | 55.1 m/s        | 19 / 36 |
| 140 000 km | ~ (iter-limited) | 27.7 m (from 28 km) | 30.9 m/s | 20 / 11 |

(A1's single-apogee re-phase root for this run is ~88 m/s — capped out, so A1 records it as
non-converged.)

Along-track Δv falls **monotonically as the second node rises toward apogee**; a mid-descent
node has no affordable authority at all. The mechanism: along-track miss is a phase/timing
error, correctable only by a tangential re-phase whose phase shift accrues over the
*remaining* arc. At 900 km there are ~minutes to the crossing → ~zero accumulation →
hundreds of m/s. Near apogee there is a full descent → the miss nulls for tens of m/s.

Two consequences make this an **apogee-bound** result, not a "move the node up" fix:

1. **A near-apogee second burn is not a distinct midcourse capability.** Where the second
   burn *does* help (80 000 km: 88 → 55 m/s) it sits right next to the first burn — it is
   operationally part of the apogee maneuver, not a separate mid-course phase.
2. **At Rung A there is no new information to act on mid-descent.** Under perfect state and
   perfect model, nothing is learned between the apogee burn and a 900 km node, so a
   mid-descent burn *cannot* beat the apogee burn — it can only re-phase, inefficiently. The
   classical statistical-midcourse second burn earns its place only at **Rung C/D**, once
   coast **drift becomes observable** (real new information), and/or for non-along-track
   error modes via the `dr_p/dv_a` perigee lever (§8) — not for along-track timing at Rung A.

**Conclusion for this phase: A1 accuracy dominates.** Getting the apogee injection and the
apogee correction right *is* the along-track lever; there is no cheap mid-descent rescue for
a timing error. (Budget context: the near-apogee ~31 m/s tail correction is ~6.2 % propellant
mass at 50 s Isp but ~1.56 % at 200 s — under the paper's 2 % claim; the Isp sweep, ADR 0004,
decides whether the tail is affordable.)

## Decision

- **A2 is reframed as this negative controllability result.** Its deliverable is the node
  sweep above and the apogee-bound conclusion — *not* a same-seed "A2 beats A1 on Δv" win.
- **Keep the pure two-burn solver.** `solve_two_burn_correction` + its unit tests (control.py,
  1f4e7e3) stay as the auditable tool behind this finding and a seed for the A3 least-norm
  controllability map. It is a tested library primitive, **not wired into `run_ensemble`**.
- **Do not retain the mid-descent harness glue.** The prototyped two-leg segmented-propagation
  + fixed-node knob encoded the disproven premise; it was reverted rather than committed. The
  table + recipe below make the finding reproducible without it.

## Considered options

- **Revise ADR 0005's node upward (near apogee) and keep A2 as a corrector that beats A1 on
  the tail** — rejected. A near-apogee second burn is barely separable from the apogee burn
  and carries no new information at Rung A, so it is apogee-maneuver efficiency, not a distinct
  control gain. Folding it into "A1 accuracy" is the honest reading.
- **Commit the mid-descent harness glue for reproducibility** — rejected. Production code for
  a disproven node is clutter; the sweep table + recipe reproduce the finding without it.
- **Revert the pure solver too (documentation-only finding)** — rejected. It is small, pure,
  always-green, and is both the auditable evidence and an A3 seed.

## Consequences

- ADR 0005 is superseded. Design doc §13 (A2 bullet) and §15 (acceptance) are rewritten to
  this finding; §9 and §16.6 keep the classical two-burn schedule but are annotated that the
  mid-descent correction buys ~0 along-track authority at Rung A and is a Rung C/D drift tool.
- A1, the open-loop capstone, and A3 are unaffected. `montecarlo.run_ensemble` keeps only the
  A1 single-burn control path; there is no A2 harness path.

## Reproduction recipe

Run `montecarlo.run_ensemble(DispersionSpec(), n=2, master_seed=20260608, control=ctrl)` with
`ctrl = lambda predict, target: solve_two_burn_correction(predict, target, max_step_m_s=5.0,
max_iter=12)`, where the harness applies burn 1 at apogee and burn 2 via a segmented
descending-altitude event-restart at node ∈ {900 km, 5, 30, 80, 140 Mm} (a `predict` that
maps the stacked RTN-6 Δv → 200 km-crossing position). Read `result.records[1]` (run 1). The
solver is in `control.py`; only the two-leg harness glue must be re-prototyped.
