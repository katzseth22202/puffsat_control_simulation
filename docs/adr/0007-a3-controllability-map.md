# Rung A3 deterministic controllability map: perfect-model coefficient sweep, LM corrector, ToA-gated authority boundary

**Status:** accepted

## Context

A3 (design doc §13) holds the targeter fixed and sweeps `Cd·(A/m)` and
`Cr·(A/m)` **deterministically** to map how required Δv grows with coefficient
error and to localize where the single-apogee targeter runs out of authority —
the controllability map, built before any noise or MPC. ADR 0003 finding 3
explicitly parked one decision here: *for a trustworthy authority map,
"non-converged" must mean "no sub-budget solution exists," not "Newton-from-zero
missed it."* A2's finding (ADR 0006) is orthogonal — it showed along-track/phase
authority is apogee-bound (A1 accuracy dominates) on the **injection** axis; A3 is
the **coefficient** axis against a fixed single-apogee targeter. A grilling
session (2026-06-09) settled the structure and five decisions a future reader
would find surprising.

## Decision

1. **Perfect-model deterministic sweep (not model mismatch).** `predict` ==
   `execute` == `full_force(swept coeffs)`, injection **zeroed**, and the corrector
   must hit the **fixed nominal aim point** (the nominal crossing at nominal
   coefficients). A larger `Cd·(A/m)` bends the trajectory off that fixed target;
   the Δv to drag it back with one apogee impulse — under *perfect knowledge* of the
   swept coefficient — is the mapped cost. This preserves Rung A's perfect-model
   line and isolates the coefficient axis. The map is therefore an **optimistic
   floor**: realistic *unknown* drag (a divergent onboard model) is Rung C's job,
   not A3's.

2. **The boundary is the interception-miss Δv; perigee is a diagnostic.** The
   recorded primitive is **raw required Δv to null the 200 km crossing miss** (ADR
   0003 decision 6 — Δv is the model-independent Rung-A primitive). The Δv **budget
   is not baked into the solver**: it is drawn as **post-processing contours** at the
   Isp anchors (2% of 25 kg at Isp ∈ {50, 70, 200} s; ADR 0004 made Isp a *reported
   sweep*). This separates two regions that carry different mission meaning:
   **over-budget** (a solution exists but costs more propellant than allowed — "buy
   Isp/mass") versus **uncontrollable** (no valid solution at all — "physically
   unreachable"). Perigee per grid point is an **overlay only** (the `dr_p/dv_a`
   cross-check, and the debris-disposal margin confirming an uncontrollable point
   still deorbits) — never a success criterion.

3. **Robust solver: in-house LM damping + budget-scale step cap + harness ToA
   gate.** A1's 2 m/s step cap is load-bearing (it suppresses the free-ToA spurious
   far root — ADR 0003 found uncapped Newton converging to a bogus ~96.9 m/s
   re-crossing a revolution later) but is **20× too tight** for A3, which must
   explore Δv up to and past the ~32–40 m/s budget. The fix has three parts:
   **(i)** add Levenberg-Marquardt `λI` damping to the finite-difference Newton step
   — it regularizes the near-singular Jacobian at the along-track wall so required-Δv
   grows *smoothly* toward the boundary instead of diverging; **(ii)** raise the step
   cap to budget-scale so genuine sub-budget solutions are reachable; **(iii)** add a
   **ToA-window gate** (sized off the open-loop ToA dispersion the capstone already
   measures) as the *physical* spurious-root discriminator — the far root crosses
   200 km a revolution off-nominal, so a solution whose ToA falls outside the window
   is rejected as "no valid local solution." Together these make `converged=False`
   mean region (b) genuinely, which is what ADR 0003 finding 3 required. `scipy`
   stays the **named escape hatch** (adopt `optimize.least_squares` only if the FD
   Jacobian proves too noisy near the wall, or bounded TRF is wanted) — not added
   speculatively, consistent with ADR 0003 decision 5 (own the convergence test so
   the boundary stays clean data, not a library status code).

4. **Factor-space grid straddling nominal; 2D deliverable, 1D tracer.** The sweep
   axis is the **multiplicative factor** on nominal (how the coefficients are
   sampled), running **below and above** nominal (drag up → lower perigee + one
   along-track sign; drag down → the opposite). The **σ-equivalent** (`factor =
   exp(k·s)`, `s = √(ln(1+cv²))`) is a reporting overlay that ties the deterministic
   axis back to the Rung-D sampling distribution, so the boundary's location reads as
   a probability ("the wall is at +4σ of drag error"). The deliverable is the **2D
   `Cd×Cr` grid** (the only way to show whether the two axes couple); **1D cuts**
   (sweep one at nominal of the other) are the build tracer bullet and first
   verification. Extent and resolution are **parameters**, widened after the first
   cut so the grid brackets the boundary.

5. **New pure `sweep.py`; reuse the harness; default-off solver and gate.** A3 is
   the same physics path as the capstone/A1 with only the inputs source and the
   control hook differing, so it reuses `_RunContext` / `_run_record` / the
   nominal-crossing setup verbatim. New pure module `sweep.py` holds `SweepSpec`,
   `grid_inputs(spec) -> tuple[RunInputs, ...]` (zero injection, swept `Cd/Cr`,
   nominal `f10p7`/`ap`), a lightweight `SweepResult` (`spec` + `records` +
   `nominal`), and `to_grid` + the σ/budget overlays. `run_sweep` lives in
   `montecarlo.py`. LM damping is a **default-off** option on
   `solve_apogee_correction` (A1's committed Newton path is the default, tests stay
   green); the ToA gate is an **optional `toa_window_s`** threaded through
   `_run_record`, default off (so `run_ensemble`/A1/the capstone are untouched).

## Considered options

- **Interpretation X — sweep truth coefficients but let `predict` use nominal
  (deliberate model mismatch)** — rejected. A cleaner "controllability against
  *unknown* drag" story, but it *is* Rung C's job and it breaks Rung A's
  perfect-model line. A3 measures the perfect-knowledge floor.
- **Bake a single Δv budget into the solver as a pass/fail gate** — rejected. Isp is
  a *reported sweep* (ADR 0004), so the budget is a family of contours, not one line;
  the solver returns raw min-Δv and the budget is applied in post-processing.
- **Keep A1's 2 m/s cap** — rejected. 20× too tight to explore the budget region;
  it would collapse the entire interesting map into a false "uncontrollable" blob.
- **Just raise the cap, no ToA gate** — rejected. Raising the cap re-exposes the
  free-ToA spurious far root; ToA (not Δv magnitude) is the physical discriminator
  between the real and the bogus root.
- **Adopt scipy now** — deferred. In-house LM is ~20 lines on the existing FD-Newton
  with no new dependency, and owning the convergence test keeps the authority
  boundary clean (ADR 0003 decision 5). scipy is the named escape hatch.
- **Overload `DispersionSpec`/`EnsembleResult`, or a full shared-core refactor of
  `run_ensemble`** — rejected/deferred. `DispersionSpec` is *the distribution*;
  mean/cov in `EnsembleStats` are meaningless on a deterministic grid. A separate
  `sweep.py` is cleaner, and the planned post-Rung-D consolidation refactor is the
  home for any unification.
- **1D cuts only** — rejected as the *deliverable* (can't reveal `Cd/Cr` coupling)
  but kept as the tracer bullet and first verification.
- **A JSONL sink/resume and a `make` target for A3** — deferred (YAGNI). ~120 fast
  deterministic points; the seed-based resume guard would need a grid-aware variant.
  Add only if the grid grows enough to need it — the same "sequential first" call A1
  made.

## Consequences

- New pure `sweep.py` (`SweepSpec`, `grid_inputs`, `SweepResult`, `to_grid` + σ/budget
  overlays). `montecarlo.py` gains `run_sweep` and an optional `toa_window_s` on
  `_run_record` (default off). `solve_apogee_correction` gains default-off LM damping
  (A1's committed tests stay green). CONTEXT.md gains the A3 terms (controllability
  map, authority boundary, σ-equivalent, `SweepSpec`/`SweepResult`) **when the code
  lands**, not now.
- No new third-party dependency — numpy only; scipy still deferred.
- A1, A2's kept solver, the open-loop capstone, and `run_ensemble` are unaffected:
  the new behavior is additive and default-off.
- The map is **perfect-model** — an optimistic controllability floor. The realistic
  unknown-drag question (does the loop close with a divergent onboard model / UKF
  estimate?) is Rung C, checked against this floor.
