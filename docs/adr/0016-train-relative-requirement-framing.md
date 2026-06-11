# Train-relative requirement framing: the plane's centroid retarget absorbs common-mode drift; per-PuffSat scatter owns the catch radius

**Status:** accepted

## Context

A grilling session (2026-06-11, after ADR 0015) examined Seth's question: the target
rocket-plane can "fly to the formation" — tolerate shifting its interception aim by up
to ~25 km — so does that relax the per-PuffSat requirements, or is each PuffSat's
arrival error independent?

The answer is determined entirely by the error **correlation structure**. Decompose
each PuffSat's arrival error as **e_i = c + δ_i** (a common-mode shift `c` of the whole
arrival corridor + per-PuffSat independent scatter `δ_i`):

- **The plane absorbs `c` essentially for free — and not as a blind tolerance.** By
  ADR 0006's observable-drift logic, coordinator tracking has *measured* each PuffSat's
  actual trajectory days before the plane launches: the Cr-prior bias, B1 erosion, and
  density history are revealed in the tracked orbits, not residual uncertainties. The
  plane aims at the **measured** corridor via launch-time/azimuth retargeting (Earth
  rotation alone gives ~460 m/s of along-track aim authority per second of launch-window
  slip; 25 km ≈ 54 s of slip). A **±~2 km retarget spec is ample** — the correlated 3σ
  drift is km-scale at most; 25 km is ~10× overkill.
- **The plane cannot chase `δ_i`.** Arrivals are ~1 s apart; at the few-kN RCS /
  boost-steering scale (~0.2–1 m/s² lateral) the plane repositions **~0.1–0.5 m per
  arrival gap** — a static target against hundreds-of-meters scatter. (Consistent with
  ADR 0014 scoping target-side trim at the ~1 m level.) Independent scatter must be
  corrected by each PuffSat's own terminal burn; **the 400 m thrust-limited catch
  radius requirement is unchanged**.
- **The ADR 0013 budget splits across that line.** Cr prior error (149 m): mostly
  common (shared material/model bias) and revealed pre-launch regardless. C1 nav error
  (141 m): mostly independent (per-PuffSat measurement noise) + a common node-ephemeris
  part. B1 erosion (89 m): common, deterministic, software-compensable. Descent density
  (F10.7/Ap): common over the train's minutes-long span. Injection: independent per
  deployment, plus any deployer systematic (common). So the train-relative per-PuffSat
  budget collapses toward the nav-dominated ~141 m — the core 1σ barely moves, but the
  **systematics and tails** move a lot.
- **The sim as built cannot answer the question.** `DispersionSpec` samples every axis
  independently per run — correct for the single-PuffSat scope, but it makes the
  common-vs-independent split an *unmodeled assumption*, and this question is the one
  that makes it load-bearing.

## Decision

1. **Requirements are judged train-relative.** ADR 0015's criterion (σ_lateral ≤ 1.65 m
   vs the 5 m plate, ≥99% capture; ToA ≤ ~10 ms) applies to each PuffSat **about the
   train centroid**. The **train** is the time-ordered arrival sequence (~1 s spacing,
   minutes-long span) — PuffSats never station-keep; the train exists only as the
   arrival schedule the plane flies through.
2. **Centroid drift belongs to the plane: the centroid retarget.** Declared capability
   **±~2 km** (pre-launch window/azimuth adjustment aiming at the tracked corridor, plus
   slow in-train following at the ~m/s level). Justified by observable drift (ADR 0006),
   not by tolerance: the drift is measured before commitment.
3. **MCC-2 re-scopes to the independent tail only.** The ~0.3% saturation story
   (ADR 0014: 3σ ≈ 670 m > any thrust-limited radius) assumed the full budget was
   per-PuffSat. Common tail events — a 2σ Cr model bias, a density storm — shift the
   train together and are absorbed by one retarget instead of N MCC-2 burns; MCC-2
   fires correspondingly rarer. **Headline claim: the cannonball-prior bias risk that
   ADR 0009's placeholder strategy accepts converts into a free launch-window
   adjustment — systematic-error insurance.**
4. **Rung D gains a train mode.** `DispersionSpec` axes get a shared-vs-per-unit
   designation — shared draws: Cr/Cd *model bias*, F10.7/Ap drivers, deployer
   systematic; per-unit draws: nav noise, injection scatter, coefficient unit spread.
   Deliverables split accordingly: the **centroid-drift distribution** (the plane's
   load, checked against the ±2 km retarget spec) and the **scatter about the
   centroid** (the catch radius / plate's load, checked against ADR 0015). The
   bias/spread ratio for the coefficients and the deployer systematic are **declared
   modeling inputs** (paper-side pins), not silent assumptions.

## Considered options

- **Keep the inertial-point framing (status quo)** — rejected as the *claim* basis,
  retained as the conservative bound. It charges the PuffSat for systematics the plane
  demonstrably eats, overstates MCC-2 demand, and forgoes the systematic-insurance
  claim. (The single-PuffSat ensembles remain valid measurements of per-unit
  dispersion.)
- **Per-PuffSat plane chase** — rejected on arithmetic: ~0.1–0.5 m of repositioning per
  1 s arrival gap vs hundreds-of-meters scatter.
- **A 25 km retarget spec** — rejected as overkill; ±~2 km bounds the correlated 3σ
  drift and keeps the claim modest.
- **Defer the framing wholly to Rung D** — rejected: the correlation split is a
  modeling input that needs paper-side data (coefficient bias vs unit spread, deployer
  systematic); deciding it at MC-build time risks improvising it.

## Consequences

- Rung D's MC design inherits the shared-draw axis designation (pure change:
  `DispersionSpec` gains the sharing structure when Rung D is built; `sample_run_inputs`
  splits into per-train and per-unit draws) and the two-tier deliverable.
- Design doc §10.2 gains the train-mode note; §16.3 gains the train-relative framing
  bullet. CONTEXT.md gains **Train** and **Centroid retarget**; "formation" is avoided
  (implies station-keeping that does not exist).
- The C-rung measured results are untouched: catch radius, nav, coefficient tolerances
  all stand — they are per-unit quantities and the train framing only re-assigns which
  ledger their common components land in.
- Paper-side pins added: coefficient bias-vs-spread ratio, deployer systematic, the
  plane's actual launch-window flexibility (the ±2 km spec consumes ~5 s of slip).
