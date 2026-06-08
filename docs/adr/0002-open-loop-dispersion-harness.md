# Open-loop dispersion capstone harness design

The Stage-1 capstone (design doc §13) characterizes where an *uncontrolled*,
perturbed PuffSat lands at the 200 km interception crossing, and how much it
disperses — sizing the midcourse Δv authority and surfacing any aim bias before a
controller exists. Per §14.1 it is also the reusable seeded Monte Carlo harness
that Rung D fills with control. Several design choices were settled (grilling
session 2026-06-07/08):

1. **Injection error is an RTN Cartesian Δv applied to the apogee deployment
   state**, not a Keplerian-element perturbation. A deployment imparts a velocity
   error, so a Δv on the apogee state is the physical model. In the
   Radial/Transverse/Normal local frame the **transverse** axis coincides with the
   velocity at apogee (radial velocity is zero at apsides), so T is exactly the
   `dr_p/dv_a` perigee lever (§8), R is timing, N is cross-track. Element
   perturbation was rejected: it would force inverting `dr_p/dv_a` to set σ — the
   wrong direction.

2. **Coefficients and space-weather drivers are median-nominal multiplicative
   log-normals** (`Cd·(A/m)`, `Cr·(A/m)`, F10.7, Ap; design doc §10.2). Median =
   nominal keeps the multiplier unbiased in log space; the resulting small *mean*
   bias is physical and is reported, not hidden by forcing mean = nominal.

3. **Per-run seeding via NumPy `SeedSequence.spawn`** — independent streams per run
   and standalone replay of any single run from its recorded seed (§14.2). All
   draws happen in a fixed order so a seed reproduces a run exactly.

4. **Pure core, JVM-thin loop.** `DispersionSpec`, `RunInputs`,
   `sample_run_inputs`, the RTN math, and the ensemble statistics are pure Python
   in `dispersion.py` (JVM-free, unit-tested like the rest of the config layer);
   only the propagate-and-record step in `montecarlo.py` touches Orekit. A small
   `build_propagator_from_orbit` seam lets the loop propagate from the perturbed
   apogee state.

5. **Primary metric: interception miss at 200 km**, decomposed in the
   nominal-crossing RTN frame (along-track ≈ timing, cross-track, radial) plus
   time-of-arrival. **Perigee is a dual diagnostic only** — the §8 lever, and a
   debris-disposal-safety margin (low is good; a *missed* PuffSat must deorbit,
   paper §9) — never an interception pass/fail. This corrects the inverted
   "perigee < 120 km = mission-killer" framing.

6. **One parameterized path (§14.1)** with a `control=None` hook, so Rung D reuses
   the harness with control switched on rather than forking a second loop.

Cost: introduces a NumPy dependency; adds `RunInputs` / `RunRecord` /
`EnsembleStats` value types; and the seam choices are load-bearing because Rung D
builds directly on them.
