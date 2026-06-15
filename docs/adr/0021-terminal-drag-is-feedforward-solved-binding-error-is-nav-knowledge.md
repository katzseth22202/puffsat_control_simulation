# Terminal drag is feedforward-solved; the binding terminal error is cross-track nav knowledge, not drag

**Status:** accepted

## Context

A natural, recurring intuition is that **terminal atmospheric drag** is the dominant
terminal-accuracy error — a 25 kg body falling through the upper atmosphere at ~10.8 km/s
toward a 200 km intercept *looks* like a drag-dominated problem. The B-rung and C-rung
findings falsify it by ~3 orders, but the result is scattered across B3a (anti-drag
feedforward), C3a (terminal feedforward), and C3b (dispersed drag) and is easy to
re-litigate. A 2026-06-15 session crystallized *why* drag cannot bind the lateral miss and,
more usefully, recorded **where the coefficient uncertainty that the intuition is reaching
for actually lives** — and the resource-allocation consequence for terminal accuracy.

This ADR is the dual of **ADR 0009**: 0009 settled the *modeling fidelity* of drag
(cannonball as conservative placeholder, cylinder deferred to Rung E). This one settles
where drag *ranks* among terminal errors and where its uncertainty enters the budget.

## Decision

1. **Terminal drag does not bind the lateral capture miss — three structural reasons, each
   measured.**
   - **It is along-track; the miss is cross-track.** Drag acts anti-parallel to velocity;
     the plate-frame capture miss is ⊥ v_rel by definition (ADR 0015). Drag projects onto
     the binding axis only through second-order trajectory-shape coupling. Its first-order
     effect is *time-of-arrival* — and ToA clears its ~10 ms budget by two orders
     (C3b: ToA ≤ 0.7 ms).
   - **It is exponentially late.** Density is negligible until the final ~100 km of descent
     (scale height 30–70 km), so the drag impulse lands in the last few seconds, when
     `½·a·t²` displacement is small. C3a measured the *uncompensated* crossing shift at
     **8.5 cm** — ~20× under the 1–2 m pre-measurement estimate.
   - **It is predictable and feedforward-cancelled.** The B3a anti-drag feedforward
     (16.7 mN peak, ~0.015 m/s Δv, both ADR 0004 gates pass ~20×) drives the executed
     residual to **2 mm** (~45× rejection, C3a) and absorbs *dispersed* drag to **≤ 0.19 m**
     (C3b). Even a **100 % coefficient error** leaves only 1–2 m — under the nav floor.

2. **The coefficient uncertainty the intuition reaches for migrates upstream — it is a
   coast-burn cost, not a terminal-drag cost.** Drag/SRP *coefficient* uncertainty does
   cost accuracy, but it enters as the **Cr-prior entry leg** of the hand-off budget
   (~149 m = the 0.2 ground-prior × the measured 745 m lateral/coefficient-factor; ADR 0013
   / C2a), **not** as terminal drag rejection. The mechanism: deploy ≈ apogee ≈ burn, the
   apogee burn Δv depends on `Cr·(A/m)` (∂Δv/∂Cr ≈ 8.8×10⁻³ m/s/factor), and no estimator
   precedes the burn — so the coefficient stays a ground prior and its error is a *coast/entry*
   dispersion. That leg is already covered ~34× by the 0.2 prior. **The "drag matters"
   instinct is half-right: it matters, but as the coast Cr-prior, already budgeted.**

3. **The binding terminal error is cross-track nav knowledge `σ_θ·R`.** The scenario is
   knowledge-limited, not disturbance-limited: catch radius is ~2.2× over the entry,
   propellant is < 2 %, and drag is cancelled — so the only quantity that moves capture is
   navigation grade (the relative-ranging grill, 2026-06-15; D1.1's entry×noise tightening
   of the requirement from the 10 µrad ceiling to the ~3 µrad target, ADR 0018). The miss is
   set by `σ_θ·R` at the *start* of terminal (26 m at the 2603 km hand-off, 10 µrad), a
   **noise** problem, not a force problem.

4. **Resource-allocation consequence: no further terminal-drag fidelity; terminal-accuracy
   budget goes to nav.** Concretely:
   - **No terminal drag/coefficient estimator in the loop.** The residual is mm-scale
     feedforward-cancelled and 100 % coefficient error → 1–2 m; estimating drag in terminal
     buys nothing the feedforward does not already deliver (ADR 0014 already settled
     no-UKF-in-the-terminal-loop — this records *why* it is also unnecessary, not just
     unwanted).
   - **The cannonball placeholder is fine for terminal.** No attitude-dependent cylinder
     drag model is needed for *accuracy*; the cylinder refinement only tightens margin in
     the optimistic direction and stays **Rung E** (ADR 0009, unchanged).
   - **Terminal-accuracy effort is allocated to nav knowledge `σ_θ·R`** — the ADR 0019
     multi-tracker / co-flyer levers (reduce σ_θ by √N averaging and R by a closer tracker)
     and the ADR 0020 apogee nav constellation (reduce the entry offset's nav leg), **not**
     drag modeling.

## Considered options

- **Build an attitude-dependent / cylinder terminal drag model for accuracy** — rejected:
  drag does not bind the lateral miss (decision 1); the cylinder only tightens in the
  optimistic direction and belongs in Rung E's clean pessimistic→optimistic A/B (ADR 0009).
- **Add a terminal drag/coefficient estimator (UKF-in-loop) to reject drag harder** —
  rejected: feedforward already drives the residual to mm / ≤ 0.19 m dispersed, and the
  worst-case 100 % coefficient error is under the nav floor; estimation buys nothing
  (decision 4).
- **Treat terminal accuracy as a disturbance-rejection problem** — rejected: it is a
  *knowledge* problem; nav grade `σ_θ·R` is the only lever that moves capture (decision 3).
- **Re-tighten the Cr ground prior to cut "drag" error** — rejected as a terminal lever:
  the Cr-prior leg is a *coast/entry* budget item already covered ~34× (decision 2); it
  belongs to the midcourse budget, not the terminal accuracy push.

## Consequences

- The recurring "drag will dominate terminal" intuition now has a durable, quantified
  answer; cite this when it resurfaces.
- Terminal-accuracy effort is formally pointed at **nav knowledge** (ADR 0019 multi-tracker /
  co-flyer for `σ_θ·R`; ADR 0020 apogee constellation for the entry-offset nav leg), and
  *away from* terminal drag modeling.
- **ADR 0009 is unaffected**: the cannonball stays the terminal placeholder; the cylinder
  stays Rung E. This ADR adds the orthogonal fact that the cannonball is fine for terminal
  *accuracy* (not just propellant), because terminal drag does not bind the miss.
- The drag/SRP **coefficient** uncertainty is budgeted in the **coast/entry** (Cr-prior leg,
  ADR 0013), already ~34× covered — so "drag" never re-enters as a terminal line item.
- Still true and unchanged: drag is intentionally over-modeled (cannonball-pessimistic), and
  a low perigee from inflated drag is the debris-disposal diagnostic ("low = good"), not a
  mission metric (ADR 0009).
