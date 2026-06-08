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
The JVM-side output of `montecarlo.run_ensemble`: per-run `RunRecord`s plus the
aggregate `EnsembleStats` (mean = aim bias, covariance = dispersion ellipsoid) and
the nominal reference crossing. The Stage-1 capstone is `run_ensemble(…, control=None)`.

**RTN frame**:
The satellite-local orbital frame — Radial (outward), Transverse (in-plane, toward
motion), Normal (orbit-normal). The interception miss is reported here; at apogee
the Transverse axis is the tangential `dr_p/dv_a` lever (§8).

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
