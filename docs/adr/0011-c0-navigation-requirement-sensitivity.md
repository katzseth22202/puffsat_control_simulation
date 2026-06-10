# Rung C0: the navigation requirement is a deterministic sensitivity slice — the apogee corrector cancels nav error to first order, so C0 maps Φ and a catch-radius threshold, not a sampled ensemble

**Status:** accepted

## Context

C0 is the first cash-in of the predict/execute seam (ADR 0003 decision 4): the corrector's
onboard `predict` finally diverges from truth, here by a navigation-error offset δ at the apogee
planning node. The design doc (§13) asks for a *navigation-accuracy threshold* via a "small
deterministic A3-style sweep — no filter, no ensemble," with Σ "motivated by what coordinator-node
ranging can plausibly deliver, not an arbitrary blob." A grilling session (2026-06-10) settled how
to realize that and surfaced a reframe a future reader would otherwise find surprising: **the
apogee corrector does not help against nav error**, so C0 is fundamentally a sensitivity (STM) map,
not a control experiment.

## Decision

1. **C0 is a deterministic sensitivity slice; the apogee corrector cancels nav error to first
   order.** With `P(x, u)` the truth crossing-position map, the corrector solves
   `P(x_est, u*) = target` from the estimate `x_est = x_true + δ`, then truth executes
   `P(x_true, u*)`. Because `u*` is applied *identically* in predict and execute,
   `residual = P(x_true, u*) − P(x_est, u*) ≈ −Φ·δ`, where `Φ = ∂(crossing)/∂(apogee state)` is the
   3×6 apogee→crossing STM — **independent of `u*`**. Two consequences: (i) the nav-error-induced
   miss is **uncontrollable at the apogee node by construction** — the corrector cannot reject an
   error it does not observe (it believes `x_est` is truth); the only levers are better apogee
   navigation (smaller Σ) or a shorter-lever-arm correction node (C3's terminal). (ii) The
   deliverable is **Φ**; the corrector stays in the loop only to faithfully realize the seam and
   catch second-order amplification (a large δ → large phantom `u*` → nonlinear executed arc).

2. **Per-component Φ basis, not a sampled-Σ ensemble (the load-bearing compute decision).** Because
   the residual is linear in δ, Φ *is* the complete requirement: for **any** Σ the induced
   interception-miss covariance is `Φ Σ Φᵀ`. C0 sweeps each of the 6 apogee-RTN error components
   (R/T/N position and velocity) one at a time over a log-spaced magnitude range, both signs — the
   A3 grid pattern in additive physical units. Each axis yields the sensitivity (small-offset slope
   = that column of Φ), the linearity range, and a per-axis tolerance. The structured, off-diagonal
   Σ a real filter produces is then pure post-processing (`Φ Σ Φᵀ`) and is **C1's** to supply. A
   single seeded-Gaussian draw at one representative Σ cross-checks that `Φ Σ Φᵀ` predicts the
   sample covariance — the only place randomness enters. Sampling N Gaussians per Σ-magnitude point
   (a mini-ensemble) is rejected: it reintroduces the ensemble ADR 0010 banished to Rung D, wastes
   compute given linearity, and ties the result to one assumed Σ structure.

3. **All 6 components; velocity-dominance is measured, not assumed.** `dr_p/dv_a` is ~30 km/m/s at
   the 150,000 km reference orbit — **~8× weaker** than the 0.9-Hill case the "1–2 cm/s" figure
   comes from — so apogee *position*-error sensitivity may matter more than the design-doc
   intuition imports, and even within velocity the *radial* axis is strongly along-track-amplifying
   (ADR 0003 finding 2). C0 reports the dominance ranking (which column of Φ has the largest
   induced-miss norm) as a **measured** finding.

4. **The threshold is the terminal catch radius (~few km, §9), parameterized — not the pusher-plate
   size.** A single apogee correction under imperfect apogee nav cannot hold a meter-scale plate:
   the coast lever amplifies ~1 cm/s apogee velocity error into ~300 m of crossing miss, so a
   plate-scale hit would demand ~0.03 mm/s apogee-velocity knowledge — absurd, and not the apogee
   phase's job. Plate-scale aim is **C3's** terminal sub-problem (re-navigate near perigee, short
   lever arm). C0 sizes only the *hand-off*: apogee nav good enough that the interception miss lands
   inside what the terminal burn can null (§9: "a few km"). The threshold is applied as
   post-processing, like A3's Isp/mass budget, with a 5/1/0.1 km table; this **closes the loop with
   Appendix A** (few km ÷ `dr_p/dv_a` is exactly its ~1–2 cm/s figure). The catch radius applies to
   the **lateral (T–N) miss**, since the radial component is pinned ~0 by the 200 km altitude-event
   crossing definition.

5. **Isolation: nav error is the *only* divergence.** Predict uses the **true** coefficients
   (perfect model — coefficient-knowledge error is C2); execute is **impulsive with the actuator
   off** (B1's finite burn is a separate divergence); **injection is zero** (so `x_true` is the
   nominal apogee state, the corrector chases a pure phantom, and `residual = −Φδ`;
   injection-independence is *asserted* from the first-order result, not re-measured); the **target
   stays the nominal crossing** (nav error is in own-state knowledge, never the target). Byproduct
   recorded: the phantom correction `u* = −Φ_u⁻¹ Φ_x δ` is nonzero, so nav error costs both a
   residual miss (primary metric) and wasted Δv (secondary metric).

6. **Seam: a defaulted predict-side onboard-state offset; dedicated `run_nav_sweep`, `RunInputs`
   kept clean.** The only production change is a defaulted 6-vector apogee-RTN offset threaded into
   `_run_record`'s predict closure (execute untouched); the zero default leaves A1/A2/A3/B
   byte-for-byte unchanged. The harness is a dedicated `run_nav_sweep` consuming a pure
   `NavSweepSpec` (the A3 `SweepSpec` analogue), recording `RunRecord`s (residual miss +
   phantom-Δv — no schema change) into a `NavSweepResult`. `RunInputs` stays the truth-side
   sampled-draws type — the nav offset is a deterministic predict-side knob, not a draw. Φ is
   cross-checked **two ways** (corrector-in-loop slope vs an open-loop `∂P(x,0)/∂x` finite
   difference) as the empirical proof of the decision-1 cancellation.

7. **Σ-motivation records the sensing physics; the architecture decision is C1's.** So the sweep
   brackets are physically motivated rather than arbitrary (§13), C0 anchors them to coordinator-node
   ranging at ~3000 km separation: **cross-range error = R×angle dominates** (1″ ≈ 14.5 m at
   3000 km → tens of m), range is ~cm–m (two-way time-of-flight), LOS velocity ~mm/s (carrier-phase
   Doppler), and **GNSS is sidelobe-only / non-standalone** at the 150,000 km apogee (above the
   GEO-bounded Space Service Volume). This motivates velocity sweep `1e-4–1e0 m/s` and position
   sweep `1e-1–1e4 m`. The physics favors **range-only multilateration from ≥4 coordinator nodes
   with a few-gram omni *coherent transponder* on the PuffSat** (1/R² per leg vs a retroreflector's
   1/R⁴; gain and any steering live on the capable coordinator receivers, omni on the PuffSat
   because multilateration needs *simultaneous* multi-node visibility and PuffSat-side steering would
   both fight that geometry and add mass against the few-gram ethos). **This sensing-architecture
   decision is C1's to settle and ADR** — recorded here only as the bracket-motivation C0 borrows;
   the retroreflective patch is demoted to a passive optical fallback.

## Considered options

- **A literal Σ-magnitude sweep with Gaussian sampling at each point** — rejected (decision 2):
  reintroduces the Rung-D ensemble, wasteful under linearity, Σ-structure-specific.
- **Report closed-loop miss without recognizing the cancellation** — rejected: it would read as a
  *control* result when it is a pure navigation-*sensitivity* result; the corrector cannot reduce
  the residual.
- **Threshold = pusher-plate size** — rejected (decision 4): demands absurd apogee-nav accuracy and
  conflates C0 (hand-off) with C3 (terminal aim).
- **Nav offset inside `RunInputs`** — rejected (decision 6): it is a deterministic predict-side
  knob, not a sampled draw; keeping `RunInputs` clean preserves the truth-side meaning.
- **PuffSat-side beam steering / high-gain antenna** — rejected (decision 7): does not change 1/R²,
  fights simultaneous multilateration, and adds mass against the few-gram ethos; put the gain on the
  coordinator.

## Consequences

- C0 is cheap: ~6 axes × a handful of magnitudes of short corrector solves + a 6-arc open-loop Φ
  cross-check + one Gaussian cross-check. Pure blocks (`NavSweepSpec`/grid, Φ assembly,
  `Φ Σ Φᵀ` + tolerance, seeded Gaussian sampler) are **`/tdd`**; the sensitivity values, dominance
  ranking, catch-radius table, linearity range, and phantom-Δv are **measured** findings (recorded
  in design-doc §13, like A2/B1/B3a).
- The predict-side state-offset seam is the **C1 consumer**: C1's UKF estimate-error covariance
  feeds `Φ Σ Φᵀ` through the same Φ, and its achieved covariance is checked against C0's threshold.
- The **sensing architecture** (multilateration / few-gram transponder / GNSS-at-apogee) is carried
  to **C1** and the paper's sizing notes; it is *not* settled by this ADR.
- Φ and the characterized requirement feed **Rung D** (sample nav error from the C1 covariance
  rather than re-running the filter every trajectory).

## Implementation findings (2026-06-10)

Built (`puffsat_sim/navigation.py` pure core + `montecarlo.run_nav_sweep`) and ran the default
61-cell sweep (`points_per_sign=5`, `pos 1e-1–1e4 m`, `vel 1e-4–1e0 m/s`), corrector at
`tol_m=0.01`, LM:

1. **The cancellation holds, measured.** 100% of cells converged and the residual stayed **linear
   across the entire swept range** (10 km position, 1 m/s velocity), so `−Φδ` is an excellent model
   far past the nav errors of interest — vindicating the Φ-basis / `Φ Σ Φᵀ` decision (no sampling
   needed).

2. **Velocity-dominance is overwhelming, and it's transverse.** Lateral-miss sensitivity
   ‖Φ_TN‖ = **2.15×10⁵ m per m/s** (T-vel), **4.3×10³** (N-vel), **1.3×10³** (R-vel), **0.60** (R-pos),
   **≈0** (T-pos, N-pos). Velocity beats position ~5 orders; transverse velocity beats the other
   velocity axes 50–160×. The T-pos/N-pos ≈ 0 result is the geometric fact that an along-track apogee
   displacement is a pure phase shift (same orbit, only ToA moves) — so decision 3's "measure, don't
   assume" paid off by *confirming* the design-doc intuition rather than importing it.

3. **The binding requirement closes the loop with Appendix A.** Apogee transverse velocity must be
   known to **~2.3 cm/s** (5 km catch radius) → 4.7 mm/s (1 km) → 0.47 mm/s (100 m). This both
   confirms and *refines* the paper's ~1–2 cm/s: the binding axis is the **along-track/timing
   amplification** (~215 km per m/s at the 200 km crossing), not the perigee-radial `dr_p/dv_a`
   (~30 km per m/s) — a sharper statement than Appendix A's perigee framing.

4. **Phantom Δv is real but small.** The corrector burns up to ~1 m/s (at the 1 m/s velocity cell)
   chasing the unobserved error (decision 5) — a cost, not a catastrophe.

5. **Numerical fidelity is ample here; flagged for the cm terminal.** The result is km/m-scale,
   ~5 orders above the truth model's `rel_tol=1e-10` floor (~cm at the 1.56e8 m apogee scale). float64
   representation (ULP ~35 nm at apogee) and the small-force disparity (~1e-6 of gravity, ~10 orders
   above the 1e-16 roundoff floor) are non-issues; the integrator *tolerance*, not the float, is the
   limiter — recorded as a C3 prerequisite (§13).
