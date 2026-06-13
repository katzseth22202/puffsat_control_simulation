# PuffSat Simulation — Domain Context

Shared vocabulary for the PuffSat control simulation. This file names the seams
the code is organized around so reviews and AI-assisted work stay consistent.
See `puffsat_control_sim_design.md` for the physics/control specification.

## Language

**Seam**:
A boundary where one side can be exercised or substituted without the other (the
testing sense). The master seam is **pure vs JVM** — every pure module unit-tests
without booting Orekit; others include **Predict vs Execute**, the `control=` hook,
the **Actuator**'s commanded→applied map, and the **Resume sink**. Decomposition
rule (2026-06-11): *ownership* says **when** to split a file — more than one owning
ADR in its CLAUDE.md entry means more than one reason to change (line counts are the
secondary tripwire: suspicious ~400 logic lines, act by ~600 even single-owner) —
and seams say **where** to cut: a module interface must sit on a substitution
boundary; a cut that crosses a seam is wrong even if it shrinks files.
_Avoid_: "seam" for mere file cohesion (a single-owner slice is a slice, not
necessarily a seam).

**Perturbation**:
A force acting on the PuffSat beyond Earth two-body gravity, represented as a
small frozen pure-Python spec (`Geopotential`, `ThirdBody`, `SolarRadiation`,
`AtmosphericDrag`, `Relativity`). Carries only the parameters that force needs
(`Relativity` is parameter-free); constructs no Orekit objects.
_Avoid_: force flag, perturbation model (the spec is not the Orekit model).

**Force Model**:
The Orekit object that actually applies a **Perturbation** during propagation
(`HolmesFeatherstoneAttractionModel`, `DragForce`, …). Built from a Perturbation
on the JVM side, never on the pure-config side.
_Avoid_: using "force model" for the pure spec.

**PhysicsConfig**:
A pure value object holding `tuple[Perturbation, ...]` — the set of forces active
for one run. JVM-free so it stays unit-testable without booting Orekit.
`is_keplerian` ⇔ the tuple is empty.

**Preset**:
A named, content-described `PhysicsConfig` bundle, exposed as a factory function
in `presets.py`: `two_body()`, `j2()`, `j2_third_body()`, `j2_third_body_srp()`,
`full_force()`. Names describe which **Perturbations** are active — not the
design-doc build-ladder step.
_Avoid_: "rung" in preset names; it collides with the control-stage Rung A–D
ladder (see Flagged ambiguities).

**Environment**:
The JVM-side bundle of frames and bodies the **Force Models** act in (ITRF,
EME2000, WGS84 ellipsoid, Sun, Moon), built once and handed to every builder in
`forces/build.py`. Replaces today's repeated re-derivation of the ellipsoid/Sun.
_Avoid_: context (too generic), world.

**Analytic signature**:
A pure closed-form prediction of a **Perturbation**'s effect (J2 secular rates,
tidal ratios, SRP/drag acceleration), living in that force's pure module. Used
by integration tests (to assert) and by `truth_model`'s reports (to print) —
it is not a method on the **Perturbation** interface.

**Lumped coefficient** (`Cd·(A/m)`, `Cr·(A/m)`):
The single scalar that sets each non-conservative force's magnitude — drag and SRP
acceleration are `½ρv²·Cd·(A/m)` and `P₀(d₀/r)²·Cr·(A/m)`. It folds the drag/reflectivity
coefficient, cross-sectional area, and 1/mass into one number, so it is the **only**
drag/SRP lever: the propagator's spacecraft mass (a fictitious 1 kg, ADR 0008) is pure
normalization, and the UKF estimates these two scalars directly. The defaults (`0.04`,
`0.02` m²/kg) are deliberate **conservative cannonball placeholders** — a sphere-equivalent,
attitude-independent; the real ANFO-cylinder, attitude-dependent area is the optimistic
**Rung E** refinement (ADR 0009). _Avoid_: treating it as area or as ballistic coefficient
alone (it is `Cd·A/m`, mass folded in); "use the real 25 kg mass" (a no-op — the coefficient,
not the propagator mass, is the lever).

**DispersionSpec**:
The swept knobs for one Monte Carlo ensemble — nominal values plus per-input 1σ
(injection-Δv RTN axes, log-normal `Cd·(A/m)`/`Cr·(A/m)`, F10.7/Ap). Pure value
object in `dispersion.py`. It is the *distribution*, not a sample.
_Avoid_: conflating it with the per-run draws (**RunInputs**).

**RunInputs**:
One run's sampled draws (RTN injection Δv plus the four log-normal coefficients /
drivers) and its `run_index`; a run replays standalone from `master_seed` +
`run_index` (§14.2). Pure; produced by `sample_run_inputs(rng, spec, i)`.

**EnsembleResult**:
The result of `montecarlo.run_ensemble`: per-run `RunRecord`s plus the aggregate
`EnsembleStats` (mean = aim bias, covariance = dispersion ellipsoid) and the nominal
reference crossing. A pure value type (in `records.py`, with `RunRecord`) even though
the JVM loop produces it, so the **Resume sink** can serialize it without Orekit. The
Stage-1 capstone is `run_ensemble(…, control=None)`.

**RTN frame**:
The satellite-local orbital frame — Radial (outward), Transverse (in-plane, toward
motion), Normal (orbit-normal). The **Interception miss** is reported here; at apogee
the Transverse axis is the tangential `dr_p/dv_a` lever (§8).

**Interception miss**:
Crossing position − target position at the 200 km EME2000 descent crossing — three
components in the nominal-crossing **RTN frame**. The capstone's primary metric and
the **Differential corrector**'s objective (ADR 0003). _Avoid_: perigee error
(perigee is a diagnostic, not the target — see Flagged ambiguities).

**Catch radius**:
The largest lateral error at the 800 km terminal hand-off that the thrust-limited
terminal burn can null to a capture-grade hit. It is a **capability**, not a precision
target, and **larger is better**: set by control authority (½·a_max·t² with the 400 mN
**Actuator**), measured **~500 m** in C3b (ADR 0014's working re-baseline was 400 m). It
propagates *backwards* as the requirement on everything upstream — it thresholds the C0
navigation sensitivity (ADR 0011), and the C1 nav / C2 coefficient-prior margins are read
against it — so a *smaller* catch radius makes midcourse delivery, navigation, and the
coefficient prior all *harder*. Distinct from the **endpoint floor**: the achieved miss
once inside the funnel (nav-limited, ~1 m at the 10 µrad terminal tracker grade). ADR
0014's claim pair: *catch radius ~400 m (entry, thrust-limited) + endpoint floor ~1 m
(nav-limited)* — different numbers doing different jobs. _Avoid_: using it for the
achieved miss or as a precision goal (it is the funnel size, not the hit accuracy);
"reduce the catch radius" to mean "hit more precisely" (the two move oppositely — a
tighter hit is a tracker-grade question, a bigger funnel is a control-authority one).

**MCC-2** (mid-course correction, second):
The second discrete midcourse burn of the transfer — a small **impulsive, out-of-plane
trim fired at a high node during descent** (design point ~30,000 km altitude) that
tilts the orbit plane to null the *lateral* entry error which the terminal **Catch
radius** cannot. The first correction is the gross apogee-injection correction (the
**Differential corrector** at the deployment node); MCC-2 is the fine plane-tilt that
chips the 3σ lateral tail (~671 m) down inside the ~500 m funnel. Its lever (lateral m
per m/s of trim, vs node altitude) was measured in C3c (ADR 0014); per ADR 0016 it is
scoped to the **independent** tail only (common-mode tail → the plane's **Centroid
retarget**). _Avoid_: "MCC-2" for the apogee correction (that is the first/gross
correction); "trim" for the gross burn.

**Controller**:
The `control=` hook on `run_ensemble`: a callable `(predict, target, basis) ->
ControlPlan`. Rung A1 supplies the **Differential corrector**; Rung D supplies MPC;
`control=None` is the open-loop capstone. _Avoid_: "control law" (this is a targeter,
not a feedback compensator).

**Differential corrector**:
The Rung A1 targeter — Newton iteration with a finite-difference Jacobian that solves
for the apogee Δv nulling the **Interception miss**. Pure `solve_apogee_correction`
in `control.py`, parameterized by a **Predict** callback; non-convergence is a
recorded outcome (the **Authority boundary**), not an error. A3 adds two default-off
options: Levenberg-Marquardt `λI` damping (`lm=`, regularizing the near-singular
altitude-event direction) and a ToA-window gate (`passes_toa_gate`, rejecting the
spurious far-revolution root); A1/A2 keep the plain Newton path. _Avoid_: optimizer,
MPC (MPC is the Rung-D replacement).

**Predict vs Execute**:
The two propagation roles (ADR 0003). **Predict** is the onboard model handed to the
**Controller** for its internal shooting (swapped for a divergent model at Rung C);
**Execute** is the harness propagating the applied plan against truth — the recorded
reality (the **Actuator** maps commanded→applied here at Rung B). Identical
`full_force` at Rung A; at Rung B they **diverge** — the corrector predicts impulsive
while the Actuator executes a finite burn, so the residual interception miss is the
measured actuator-realism erosion (ADR 0008). At Rung **C0** the divergence moves to
predict's *starting state*: a navigation-error offset on the apogee planning state
(ADR 0011). Because the *same* correction is applied in predict and execute, the
residual is the apogee→crossing sensitivity (the 3×6 STM Φ) times the nav error — so
nav error is **uncontrollable at the apogee node**, and C0 is a sensitivity sweep, not
a control experiment.

**Actuator (finite burn)**:
The Rung B execution layer that turns a commanded impulsive Δv (from the **Differential
corrector**, which stays impulsive) into a finite, mass-depleting burn — a single
omnidirectional proportional cold-gas thruster (400 mN max, ~5 mN floor, ~1°/s slew; ADR
0004). It is the commanded→applied map *inside* **Execute**; Isp is a post-processing sweep
on the resulting Δv, not an actuator state. _Avoid_: "thruster" for the model (the model is
the proportional abstraction of a bang-bang cluster, not the hardware).

**Coordinator node**:
The paper's off-board capable asset (gain, steering, compute; `sec:coordinator_node_dry_mass_disposal`)
that ranges the PuffSats. Co-flies with the swarm on a **matched-period** neighboring orbit
(same semi-major axis — no secular along-track drift — with perigee raised to **≥200 km**,
apogee trimmed to compensate: nodes are reusable assets, not debris, and must never enter
the burn-up zone; end-of-life disposal is a ~5–15 m/s perigee-lowering burn at apogee).
Nearby through the whole high-altitude coast. In the sim (C1 on) a node is a
**known-ephemeris beacon**: its
position is a measurement-model *input*, never a filter state; its own nav error folds
into the measurement-noise `R` as an inflation term. The node *geometry* (count, LOS
angular diversity) is a **derived requirement** — C1 outputs what geometry/Doppler quality
hits the C0 threshold; the paper consumes it as a sizing result, the sim never assumes a
fixed constellation. _Avoid_: ground station (nodes are not Earth-fixed), beacon alone
(loses the compute/latency role C4 prices).

**Train**:
The time-ordered sequence of PuffSat arrivals at the interception point (~1 s spacing,
minutes-long span). PuffSats never station-keep — each flies an independent orbit
(design doc §2 keeps fleet interactions out of scope); the train exists only as the
arrival schedule the target plane flies through. ADR 0015's accuracy criterion
(σ ≤ 1.65 m vs the 5 m plate, ToA ≤ 10 ms) is judged per PuffSat **about the train
centroid** (ADR 0016). _Avoid_: formation (implies station-keeping / relative control
that does not exist).

**Centroid retarget**:
The target plane's absorption of common-mode **Train** drift: coordinator tracking
reveals the actual arrival corridor days before launch (ADR 0006's observable drift),
and the plane aims there via launch-time/azimuth adjustment — ±~2 km declared
capability (≈5 s of window slip; Earth rotation gives ~460 m/s of aim per second).
It cannot chase per-PuffSat scatter (~0.1–0.5 m of repositioning per 1 s arrival gap),
so independent scatter stays with each PuffSat's terminal burn (the catch radius), and
MCC-2 re-scopes to independent tails (ADR 0016). _Avoid_: "the plane flies to meet the
PuffSat" (true at centroid scale only, never per-PuffSat).

**Transponder (PuffSat RF link)**:
The PuffSat's only RF nav hardware — a few-gram omni **coherent turnaround**
(receive-mix-amplify MMIC + patch antenna, mW): it phase-locks to the interrogating
carrier and re-transmits at a fixed ratio, so two-way range/Doppler precision is bound
by the *interrogator's* OCXO and **the PuffSat carries no precision clock** (ADR 0011
§7, ADR 0012 §3). A ~1 g TCXO (no oven, no temperature control) covers local mixing
and ToA holdover between contacts; two-way time transfer resyncs every interrogation.
Chosen over a passive retroreflector for midcourse: 1/R² per leg vs 1/R⁴, and
multilateration needs simultaneous omni visibility. _Avoid_: "homing beacon"
(overloaded — see Flagged ambiguities); one-way beacon designs (those are what would
need an onboard precision clock).

**ControlAction / ControlPlan**:
Pure value types in `control.py`. A **ControlAction** is one commanded maneuver (node
+ RTN Δv + magnitude); a **ControlPlan** is the ordered actions a **Controller**
returns plus its `converged` / `iterations` metadata. The plan is logged into
`RunRecord.control_log`. _Avoid_: "burn" for the value type (a burn is the physical
event; these are the record).

**Resume sink**:
The newline-delimited JSON store of completed `RunRecord`s keyed by `run_index`.
Recovery is run-granular: on restart, run only the missing indices and re-summarize —
enabled by per-run seed reproducibility (`replay_inputs`, §14.2), never within-run
integrator snapshots. _Avoid_: "checkpoint" implying mid-propagation state.

**Controllability map (A3)**:
The deterministic `Cd·(A/m)` × `Cr·(A/m)` sweep (ADR 0007): hold the targeter fixed, sweep
the two lumped coefficients across a factor grid straddling nominal (injection zeroed), and
record required Δv per cell. Built in `sweep.py` (`SweepSpec`, `grid_inputs`, `to_grid` →
`SweepGrid`) and run by `montecarlo.run_sweep` → **SweepResult**. Perfect-model, so it is an
*optimistic floor* — the unknown-drag question is Rung C. The built map came out controllable
everywhere at Δv « budget and ~1D in `Cr`. _Avoid_: confusing it with the stochastic
**DispersionSpec** ensemble — A3 is a deterministic grid, not a distribution.

**SweepSpec / SweepResult**:
Pure value types in `sweep.py`. **SweepSpec** is the deterministic grid (nominal `Cd`/`Cr` +
per-axis factor range + resolution); **SweepResult** bundles the spec, the per-cell
`RunRecord`s, and a dedicated factor-(1,1) `nominal` reference run. `SweepSpec` is to A3 what
**DispersionSpec** is to the capstone — but a grid, not a distribution.

**Authority boundary**:
Where the fixed targeter runs out of control authority for a sub-budget solution — recorded
as `converged=False`, never an exception (ADR 0003 finding 3). Resolved post-hoc into
*over-budget* (a solution exists but exceeds the Isp/mass budget) vs *uncontrollable* (no
valid solution) by `classify_controllability`. For A3's coefficient axis the boundary is not
reached (controllable everywhere); it bites on the injection axis (A1).

**σ-equivalent**:
The reporting overlay (`sweep.sigma_equivalent`) mapping a deterministic sweep factor to its
σ in the Rung-D log-normal (`k = ln(factor)/s`, `s = √(ln(1+cv²))`), so the **Controllability
map**'s axis reads in σ of coefficient error — tying the A3 grid back to the capstone's
sampling distribution.

**Rung D (D1 / D2)**:
The feasibility Monte Carlo, split (ADR 0018). **D1** is the full closed-loop ensemble on
the *C baseline* (corrector + C3b ZEM + **MCC-2** trim + finite burn) — it produces the
**conditional** feasibility yes/no (headline **P(capture)** about the train centroid). **D2**
prototypes MPC and measures it against the C baseline on the *same* ensemble; MPC earns its
place only on a measured D1 threshold/constraint violation (§16.10). "Conditional" = feasible
*given* the nav/actuator specs, which the MC takes as inputs and the **Tracker budget** +
GDOP/torque analyses convert toward derived requirements. _Avoid_: treating "true MPC" as a
prerequisite for the feasibility verdict (it is a D2 value question, not the D1 gate).

**Tracker budget**:
The pure derivation (a D1 blocking gate, ADR 0018) of what the 10 µrad terminal
**Tracker grade** demands — aperture / exposure / residual jitter / SNR for a dim, fast target
on a shaking bus, plus *acquisition* (tracker FOV vs the hand-off delivery Σ). It promotes the
load-bearing terminal-nav assumption from an assumed spec to a derived hardware requirement;
if 10 µrad is unmeetable the **Catch radius** story falls, so it blocks the MC. _Avoid_:
confusing it with the **Tracker grade** itself (the budget is what *achieving* that grade
costs, not the grade).

## Relationships

- A **PhysicsConfig** contains zero or more **Perturbations**.
- Each **Perturbation** maps to one or more **Force Models** via a single
  JVM-side dispatch (`forces/build.py`).
- A **Preset** is a **PhysicsConfig** with a content-describing name.
- Each **Perturbation** owns its **Analytic signature** in the same pure module.
- The pure-spec → Force-Model split is an **internal seam**: Orekit runs live in
  the integration test suite, so no mock adapter is needed.
- **Verification vs report.** The integration tests are the verification surface
  (they assert against the **Analytic signature**). `truth_model`'s `report_*`
  functions only print — they are a demo, not a check.
- Physical constants shared by ≥2 force modules (Earth radius, μ, J2) live in one
  `constants.py`; force-specific constants (atmosphere layers, SRP P₀) live in
  that force's module.
- A **Controller** consumes a **Predict** callback and returns a **ControlPlan** of
  **ControlAction**s; the harness applies it via **Execute** and logs it into the
  **RunRecord** (`control_log` + `total_dv_m_s`).
- The **Differential corrector** is the Rung A1 **Controller**: it nulls the
  **Interception miss**, with perigee / `dr_p/dv_a` kept only as an acceptance
  cross-check.
- An **EnsembleResult**'s **RunRecord**s stream to the **Resume sink** keyed by
  `run_index`; **EnsembleStats** is recomputed from the reloaded set.
- The **Controllability map** reuses the corrector and harness: `run_sweep` shares
  `_build_context` / `_run_record` with `run_ensemble`, and its **SweepResult** is the
  deterministic counterpart of **EnsembleResult**.

## Example dialogue

> **Dev:** "Where does `f10p7` live now?"
> **Architect:** "Inside the `AtmosphericDrag` **Perturbation**, not on
> **PhysicsConfig** — it only means anything when drag is active."

## Flagged ambiguities

- **"Rung" overloaded.** The design doc has two ladders: a *physics* ladder
  (Rung 0/1/2a–2d) and a *control* ladder (Rung A–D). `truth_model.py`'s
  docstring mislabels itself "Rung A truth model"; design-doc Rung A is the
  impulsive-Δv controllability core, not truth propagation. Resolution:
  **Presets are named by content, never by rung.**

- **"Homing beacon" overloaded (resolved 2026-06-11, no architecture change).** Two
  distinct hardware roles hide under the phrase: the **Transponder** (RF, clock-free,
  midcourse range/Doppler + ToA/time transfer — ADR 0011/0012) and the terminal
  **optical beacon** (the point source ADR 0015's astrometric tracker centroids; a
  light source, no clock). A "go back to retroreflective tape" proposal was examined
  and dropped: the motivating clock worry was unfounded (the two-way coherent link
  keeps all precision timing on the interrogator; one-way designs are what would need
  a CSAC-class onboard clock — that is why TDOA was never adopted), and ADR 0011 had
  already rejected the retroreflector for midcourse on 1/R² vs 1/R⁴. Residual idea,
  noted not adopted: a corner-cube LRA (~10 g, passive) as a dead-PuffSat
  safety/debris-verification tracker — a paper-side option, no sim impact.

- **Perigee as target vs diagnostic (resolved).** The design doc's A1 wording says
  the corrector nulls "predicted perigee error," but perigee is a *diagnostic*
  (debris-disposal margin + the §8 lever), not the objective. Resolution (ADR 0003):
  the **Differential corrector** nulls the **Interception miss**; perigee /
  `dr_p/dv_a` is an acceptance cross-check only. The design doc A1 bullet is
  superseded on this point.
