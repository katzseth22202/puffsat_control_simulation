# PuffSat Simulation — Domain Context

Shared vocabulary for the PuffSat control simulation. This file names the seams
the code is organized around so reviews and AI-assisted work stay consistent.
See `puffsat_control_sim_design.md` for the physics/control specification.

## Language

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
measured actuator-realism erosion (ADR 0008).

**Actuator (finite burn)**:
The Rung B execution layer that turns a commanded impulsive Δv (from the **Differential
corrector**, which stays impulsive) into a finite, mass-depleting burn — a single
omnidirectional proportional cold-gas thruster (400 mN max, ~5 mN floor, ~1°/s slew; ADR
0004). It is the commanded→applied map *inside* **Execute**; Isp is a post-processing sweep
on the resulting Δv, not an actuator state. _Avoid_: "thruster" for the model (the model is
the proportional abstraction of a bang-bang cluster, not the hardware).

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

- **Perigee as target vs diagnostic (resolved).** The design doc's A1 wording says
  the corrector nulls "predicted perigee error," but perigee is a *diagnostic*
  (debris-disposal margin + the §8 lever), not the objective. Resolution (ADR 0003):
  the **Differential corrector** nulls the **Interception miss**; perigee /
  `dr_p/dv_a` is an acceptance cross-check only. The design doc A1 bullet is
  superseded on this point.
