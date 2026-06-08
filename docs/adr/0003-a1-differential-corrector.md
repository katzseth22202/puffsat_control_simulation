# Rung A1 differential corrector design

Rung A1 (design doc §13, "single impulsive midcourse") closes the loop for the
first time: a deterministic differential corrector fills `run_ensemble`'s reserved
`control` hook (ADR 0002 §6, §14.1), establishing the controllability ground truth
— the Δv floor and authority boundary — that MPC is later checked against. It is
perfect-state, perfect-model, impulsive. Decisions settled (grilling session
2026-06-08):

1. **The corrector nulls the EME2000 interception position miss, not perigee.**
   It drives the 3-component RTN *position* miss at the 200 km descent crossing
   (the §16.3 target the capstone already measures) to zero. A single impulsive Δv
   has exactly 3 DOF, matching the 3 position components; time-of-arrival floats
   (the crossing is an altitude event, not a clock — the moving target is Rung D).
   Perigee / `dr_p/dv_a` is demoted to an **acceptance cross-check** (§15), never
   the objective. This **supersedes the design doc's A1 wording** ("vary the
   correction Δv to null predicted *perigee* error"), which predates the 2026-06-08
   correction that made perigee a diagnostic and the interception miss the mission
   metric (ADR 0002 §5). Keeping the objective in the same frame the open-loop
   dispersion was measured in makes "did the loop close" apples-to-apples with the
   baseline it must beat.

2. **Single burn at the apogee deployment node, zero latency.** Apogee is the
   cheapest-leverage point and the cleanest analytic anchor: radial velocity is
   zero there, so RTN-Transverse coincides exactly with the velocity direction and
   the `dr_p/dv_a` cross-check is exact (it degrades as the node moves down the
   arc). Because the node *is* the initial epoch, the correction is just an
   adjustment to the run's initial velocity — it reuses
   `build_propagator_from_orbit` → `_propagate_to_interception` with no
   `ImpulseManeuver`. It is **not** degenerate with the injection: the corrector
   solves a boundary-value problem to the nominal crossing under full nonlinear
   `full_force` dynamics, nulling 3 position components — establishing the
   *best-case* Δv floor, which A2/A3/D then erode. The downstream mid-descent node
   and the burn split are deferred to A2; zero latency matches Rung A's
   idealization (latency enters at Rung B).

3. **Predictor = the truth model; finite-difference Jacobian; Newton with recorded
   non-convergence.** At Rung A the premise is perfect state *and* perfect model,
   so the corrector's internal predictor carries the run's true drawn coefficients
   (`physics_from_inputs`) — mismatch is exactly zero and it converges to machine
   precision (an onboard model that diverges from truth is the new thing Rung C
   introduces). The 3×3 Jacobian is by finite difference (3 extra arcs per
   iteration), **not** Orekit's State Transition Matrix: it keeps the corrector a
   black-box-propagator method that survives unchanged into B (finite burns) and C
   (model mismatch). Newton with a step cap, convergence at ‖position miss‖ < ~1 m
   (well under the km-scale dispersion), ~8-iteration ceiling. Non-convergence is
   **recorded as a run outcome, never thrown** — at A3/D a swept coefficient can
   push the orbit past the authority boundary, and that failure *is* the
   controllability-map signal.

4. **Two distinct propagation roles: `predict` vs `execute`.** The harness hands
   the controller a `predict` callback (the "onboard model" — swap for a divergent
   model at Rung C) and keeps its own `execute` propagation of the applied plan
   against truth (the recorded reality — an actuator model maps commanded→applied
   here at Rung B). At Rung A both run identical `full_force(true coeffs)`, so the
   recorded miss ≈ 0 to machine precision. They are made separate seams *now*,
   while identical, precisely so B and C drop in without reshaping the loop; the
   converged plan is re-propagated once through `execute` for the authoritative
   record (~1 arc on top of the 3–5 solve arcs).

5. **Own the solver (pure Newton + FD); no control library.** The corrector is a
   pure-numpy `solve_apogee_correction(predict, target_position, …) -> ControlPlan`
   parameterized by a `predict: Callable[[Vec3], Vec3]` callback (RTN Δv → EME2000
   crossing position) — unit-testable with a synthetic predict, no JVM. We rolled
   our own rather than adding `scipy.optimize.root` because the problem is trivial
   for the method (near-linear, well-conditioned, 2–4 iterations) and, decisively,
   **the corrector's non-convergence is the A3 deliverable** — owning the
   convergence test keeps the authority boundary as clean data rather than a
   library's hidden status code (a robust solver would actively mask the boundary).
   `scipy` is deferred (triggers: trust-region robustness near the A3 edge, or FD
   cost dominating at N=10⁴). `python-control` is rejected as wrong-category
   (feedback-control design, not trajectory targeting); the MPC library is a
   Rung-D decision. Swappability is already provided by the `Controller` hook and
   the `solve_*` signature — no strategy abstraction is added (YAGNI). New pure
   module `control.py` holds the value types and solver; JVM glue lives in
   `montecarlo.py`. `Controller = Callable[[PredictFn, Target, Basis], ControlPlan]`
   replaces the `Callable[..., Any]` placeholder.

6. **`RunRecord` grows by a control superset; propellant-mass and applied-Δv
   deferred.** Added: `control_log: tuple[ControlAction, ...]` (the commanded plan;
   one action at A1 — the §14.2 dense-replay "commanded Δv"), `total_dv_m_s` (the
   Rung A Δv-floor ledger), and `converged: bool` / `iterations: int` (per-run
   authority-boundary queryability). `EnsembleStats` gains Δv mean/std/max and a
   converged-fraction. Explicitly **not** added: a propellant *mass* / Isp fraction
   (Isp is a Rung-B actuator detail; Δv is the model-independent Rung-A primitive,
   and mass conversion is a reporting helper once Isp is pinned — §14.3); an
   `applied` Δv field (commanded = applied at Rung A; a Rung-B bump); and dense
   per-run trajectories (regenerated on demand via `replay_inputs`, §14.3 — never
   bulk-persisted). The open-loop capstone (`control=None`) leaves `control_log=()`,
   `total_dv_m_s=0.0`, `converged=True` — a superset, so it validates unchanged.
   `RunRecord` stays *deeply* immutable (tuples, frozen dataclasses), which is its
   entire thread-safety story.

7. **Resume sink: JSONL per-`run_index`, resume-by-complement; sequential now,
   process-sharding deferred.** Because runs are independent and seed-reproducible
   (`replay_inputs`), recovery is run-granular — no within-run integrator
   snapshots. Completed `RunRecord`s stream to a newline-delimited JSON sink keyed
   by `run_index`; on restart the harness runs only the missing indices and
   `summarize`s over the full reloaded set. Shipped against the sequential loop
   first (A1/A2 at N=50 don't need throughput). When N=10³–10⁴ demands speed,
   parallelism will be **process-based** — each worker boots its own JVM, because
   Orekit propagators and `FramesFactory`/`DataContext` caches are *not*
   thread-safe — with a **shard-per-worker / file-per-index** sink (a file existing
   = that index is done = lock-free resume). The JSONL/per-index keying is chosen
   now so that follow-on is purely additive (shards are concatenated JSONL).
   `RunRecord` needs no thread-safety: it is immutable, and process isolation means
   no cross-worker sharing regardless.

Consequences: introduces `control.py`; grows `RunRecord` / `EnsembleStats` and
replaces the `Controller = Callable[..., Any]` placeholder with a typed signature;
adds a JSONL sink and a resume path. The **predict/execute seam** (decision 4) and
the **pure-solver signature** (decision 5) are load-bearing — Rung B (actuator on
`execute`), Rung C (onboard model on `predict`), and Rung D (MPC behind the
`Controller` hook, process-parallel sink) all build directly on them. No new
third-party dependency (numpy only; scipy explicitly deferred).
