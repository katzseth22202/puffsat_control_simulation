# Rung A2 two-burn statistical midcourse: minimum-Δv least-norm targeter

**Status:** superseded by [ADR 0006](0006-a2-along-track-is-apogee-bound.md) — the
mid-descent (900 km) node it specifies provides ~0 along-track authority; along-track/phase
error is apogee-bound. The minimum-Δv least-norm *solver* below stands and is retained; the
*node placement* and the "cheap two-burn split beats A1" premise do not.

## Context

A2 (design doc §13) adds a second coast correction to A1's single apogee impulse
(ADR 0003): correction 1 near apogee, correction 2 at a fixed mid-descent altitude
(§16.6). Still impulsive, perfect-state, perfect-model. A grilling session
(2026-06-09) settled the structure and four choices a future reader would find
surprising. The driving carry-forward from A1: a single apogee impulse is
Δv-*inefficient* for along-track/phase error — a 2.4σ radial injection cost ~88 m/s
to null from apogee alone — so A2 exists to show two well-placed impulses null the
same miss cheaply.

## Decision

A2 is **exactly two impulsive burns** (apogee + mid-descent at 900 km), nulling the
3-component interception position miss. Four decisions worth recording:

1. **Surplus DOF → minimum-Δv, via least-norm Gauss-Newton.** Two burns are 6 DOF
   against 3 position constraints — underdetermined. The spare 3 DOF are spent on Δv
   efficiency (not on a richer position+velocity target, and not on a sequential
   cheap-then-cleanup split): A2 solves for the minimum-Σ‖Δv‖² correction. This is
   the most faithful test of §13's "bulk of dispersion nulled cheaply," and the
   measurable result is A2 total Δv ≪ A1's on the same seeds.

2. **Least-norm stays pure-numpy; scipy stays deferred to A3.** The underdetermined
   minimum-norm step is `pinv(J) @ residual` (`np.linalg.lstsq`) on the same
   finite-difference Jacobian (now 3×6) — one line off A1's square `np.linalg.solve`.
   ADR 0003's scipy deferral was about trust-region robustness near the A3 *authority
   boundary*, a different concern that still belongs in A3. A1's per-step `max_step`
   cap carries over unchanged (same free-ToA spurious-far-root risk).

3. **No terminal aim at Rung A.** §9's terminal aim is a *continuous, drag-rejecting*
   600→200 km burn; under A2's perfect state + perfect model an *impulsive* terminal
   aim is vestigial — converged coast burns null everything downstream to machine
   precision, so it solves for identically zero. It becomes real as a finite burn at
   Rung B. The 800/600 km hand-offs remain as event-firing propagation checks only.

4. **Segmented event-restart for burn 2.** Burn 2 fires at a descending 900 km
   altitude event, applied by propagate-stop-apply-restart (design doc §4), not an
   Orekit `ImpulseManeuver`. Its RTN basis is re-derived from the stopped node state
   each predict call (state-dependent because perturbing burn 1 moves the node), which
   keeps the corrector a pure black-box-propagator method (ADR 0003 decision 5). RTN
   is orthonormal per burn, so minimizing Σ‖Δv‖² in the mixed (apogee-RTN, node-RTN)
   coordinates equals minimizing true EME2000 Σ‖Δv‖².

The deliverable is a **same-seed A2-vs-A1 comparison**: Δv ledger, per-burn split,
and convergence-fraction gain. No total-Δv budget gate in the solver — mapping where
Δv exceeds budget is the A3 controllability map.

## Considered options

- **Spend surplus DOF on a richer 6-state target (null position+velocity, 6×6 square
  Newton)** — rejected: at Rung A the target is a fixed point/epoch (§16.3; moving
  target is Rung D), so pinning crossing velocity solves a problem we don't have and
  *increases* Δv rather than measuring its floor. Pre-stages Rung D but premature.
- **Sequential cheap-then-cleanup (capped A1-solve at apogee, residual at mid-descent)**
  — rejected: the burn-1 cap policy is arbitrary and it yields *a* feasible split, not
  the minimum two-burn Δv that A2 is meant to measure.
- **Retain the terminal aim as a wired-but-zero placeholder node** — rejected: adds
  3 DOF and a node modelling nothing at Rung A; the B/C seam is already provided by
  the multi-action `ControlPlan` and the predict/execute split.
- **Orekit `ImpulseManeuver` for burn 2** — rejected: buries the Δv inside event
  handling and muddies the predict/execute seam; segmented restart keeps it in our
  hands, consistent with A1 and §4.
- **True minimum-propellant (L1, Σ‖Δvᵢ‖) and a total-Δv budget gate** — deferred:
  least-norm (Σ‖Δv‖²) is a pure-numpy proxy; the L1 objective and the budget-boundary
  map are A3 concerns.

## Consequences

- `control.py` grows a two-burn least-norm solver (`PredictFn` generalizes from
  RTN-3→position to RTN-6→position; `ControlPlan` carries two `ControlAction`s — the
  type already supports it). The harness closure (`montecarlo.py`) gains the
  segmented-restart mid-descent burn (the `ImpulseManeuver` TODO at the A1 call site).
- CONTEXT.md's "Differential corrector" entry will move from "square Newton for the
  apogee Δv" to "two-burn least-norm" — updated when A2 is *implemented*, not now.
- A1, the capstone, and A3 are unaffected; A2 reuses A1's FD-Jacobian/Newton machinery
  and `max_step` cap. FD cost doubles (3→6 extra arcs/iteration) — negligible at N=50,
  a factor in the N=10⁴ scipy/parallelism trigger (ADR 0003).
