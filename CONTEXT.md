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

**TrainDispersionSpec**:
The **Train**-mode dispersion (D1.0; `train.py`, ADR 0016/0018): each **DispersionSpec**
σ split into a **shared** (per-train, common-mode: coefficient *bias*, F10.7/Ap
*drivers*, deployer *systematic*) and a **per-unit** part (coefficient *spread*,
injection *scatter*). `sample_train` draws the shared part once per train and composes
`n_units` **RunInputs**; `replay_train_unit` is the standalone replay. The verdict
(`summarize_train_capture` → `TrainCaptureStats`) splits a train's arrivals into
**Centroid retarget** drift (common-mode, vs ±2 km) and scatter about the centroid
(per-unit, vs the plate). D1.1 adds the **hand-off entry offset** (`sample_train_entry_offsets`):
the Φ-composed midcourse residual (per-unit C1 nav 141 m + shared C2a Cr-prior 149 m,
2-D lateral) the C3b terminal funnel flies from — `summarize_train_ensemble` →
`TrainEnsembleFinding` (capture + propellant + perigee) is the D1.1 verdict surface.
_Avoid_: reading it as a station-keeping formation (a **Train** never station-keeps);
conflating shared bias with per-unit spread; assuming the centroid retarget absorbs the
entry offset (the funnel removes the common-mode entry in flight — the retarget is the
pre-launch backstop for common-mode beyond the funnel authority).

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
if 10 µrad is unmeetable the **Catch radius** story falls, so it blocks the MC. The **smear**
term's defensibility on a 1 Hz-hammered bus rests on a stack (grill 2026-06-15): differential
astrometry **cancels** *rigid* bus motion (common-mode), so the 3 mrad/s residual is the
*post-cancellation* rate; *structural flexure* (not common-mode) is handled by **time
separation** — the propulsion concept's Orion-style shock absorber + a small-stiff optic
(active-beacon → 5 cm aperture → high first mode) ring down ≪ the 1 s inter-impact window, so
the ToA-scheduled gate covers it; the continuous N-star plate solve is a **live distortion
monitor** for any quasi-static residual. Net: the grade is **calibration-limited, not
vibration-limited** — the **Tracker array** (ADR 0019) is the hedge if the mount's
first-mode/damping bench number comes back bad. _Avoid_: confusing it with the **Tracker
grade** itself (the budget is what *achieving* that grade costs, not the grade); reading the
smear cancellation as covering *flexure* (it covers rigid motion only — flexure is gated out
in time, not cancelled).

**Tracker array**:
N independent target-side detectors (ADR 0019) fused to beat a single **Tracker grade** by
**√N down to a common-mode floor** (correlated distortion + beacon-shape asymmetry). Independence
is the whole game: separate optics, *separately bench-calibrated* distortion maps, each doing its
own beacon-vs-its-own-star-field astrometry — then 5× 10 µrad → ~1.4 µrad. The X-pattern *spatial*
spread is for **coverage + redundancy + common-boresight rejection, not precision** (√N is
statistical, separation-agnostic). Hedge audit (2026-06-15): the fusion model banks the *entire*
3 µrad distortion in the independent bucket (**zero** assumed common-mode distortion), so the
1.62 µrad 5-array grade is the *optimistic* end — the open question is the **common-mode-distortion
tolerance** (how much correlated distortion before 1.62 µrad exceeds the 3.2 µrad capture grade),
folded into the **Differential astrometry** distortion-field study. _Avoid_: crediting the baseline
with precision, or expecting
**ranging** to sharpen the lateral (a short baseline gives ~33 mrad — angles do the lateral work).
Inter-PuffSat *cameras*/beacon cross-links and a *moving target plate* are the same trap (grill
2026-06-15): common-mode tools that sharpen the swarm's internal scatter or track the centroid
drift — both already absorbed by the **Centroid retarget** + plate width — and are blind to the
binding *per-unit* residual, which is **engine reach** past the **Catch radius** (a thrust limit,
not relative-nav knowledge). **Refined (grill 2026-06-16):** that trap holds for *passive*
cameras / move-the-plate (sub-arch A — reads only the common-mode, leaves the per-unit residual);
it is *escaped* by an **independent anchor + camera-rigid swarm** (see **Surveyor-anchored
centering**) for the *deferred cm-centering follow-on* — there the per-unit residual is optical-nav,
not engine-reach (cm trims sit deep inside the 475 m **Catch radius** funnel). Inverts near the Sun
(relative placement *is* the mission, no anchor floor), where this family is worth re-exploring as a
close-in collision sensor.

**Surveyor-anchored centering** (the deferred cm-centering follow-on, grill 2026-06-16):
The scheme for the paper's deferred **terminal cm-centering** sub-problem (§2) — *not* a change to
the committed plate-capture baseline (which stands on the 5 m plate / closed-out D1); the optional
tightening that shrinks the pusher plate from 5 m toward ~10 cm by driving the *per-unit* arrival
miss to cm. Two levers on two axes. **(1) Block shift** — a **sacrificial first PuffSat**
("surveyor") whose true crossing is read by an **independent one-shot instrumented gate** (a
**hoop**: lidar/microwave trilateration — *not* the optical tracker, which would be circular against
its own bias; *not* a hit, which dumps off-center torque + flash + shock) pins the swarm's
quasi-static optical-distortion *bias* to the plate. **(2) Per-unit scatter** — **strobed
known-pattern LED beacons** make the swarm **camera-rigid**: inter-PuffSat *bearing* gives the
binding cross-track at the short intra-train range (σ_δ = σ_θ·v/f), PnP *range* lands on the
non-binding along-track, **blink-code ID** (no cross-unit clock — cross-track is range-insensitive,
so a shared clock would buy only ID, which blink-codes give for free), **anchor-as-surveyor**
topology (the first unit looks *backward* at dark sky, maps the trailing swarm, gets hoop-pinned,
broadcasts each follower's plate-relative offset → followers never stare into the forward
flash/glare). The per-unit floor moves from geometry (easy — cm at 2–4 Hz even with a crude camera)
to **metrology**: σ_hoop (≤1 cm → 5 cm plate; ~3 cm → 10 cm plate) and the camera's *own* calibrated
distortion floor (differential-star astrometry → ~3 µrad). f is **passenger-g-bounded** (2–3 g:
ascent gravity-loss floor below, comfort ceiling above; at fixed g, f ∝ rocket mass → heavier
rockets run tighter trains → smaller plates); flash recovery (5 ms) and the ~95 ms mount ring-down
are both non-binding ≤4 Hz, and the projectile-side cameras are immune regardless. Long train →
periodic re-anchoring (~1 unit/min, <1%) if the bias drifts. **Committed claim: 10 cm plate (robust,
~50× off the 5 m); 5 cm stretch contingent on σ_hoop ≤1 cm + a calibrated camera; 2 cm dropped.**
Authority is *not* the binder (cm trims, 475 m funnel) — this stays **knowledge/metrology-limited**.
_Avoid_: reading it as a change to the committed criterion; "cameras cut per-unit scatter" *without*
the independent anchor (= the **Tracker array** trap); demanding a cross-unit clock.

**Differential astrometry** (the σ_θ distortion hedge):
Measuring the target beacon's bearing *relative to reference stars in the same FOV/exposure*
(inertial directions known to µas) rather than absolutely on the focal plane — the Gaia trick: the
focal-plane **distortion** common to beacon + nearby stars cancels, leaving only the distortion
*gradient* over the small target-to-star separation. The acquisition design already needs ≥3
reference stars in the FOV (`tracker_budget.py`); this uses them *metrically*. Centerpiece of the
terminal cross-track **hedge** program (idea menu 2026-06-15, `todos/improve_terminal_crosstrack.md`):
the binding terminal axis is the cross-track-to-target position σ_θ·R, and D1 is already feasible
*past the noise knee*, so this is **defense-in-depth on the load-bearing 3 µrad distortion-calibration
number plus extra capture margin**, not closing a capture gap. Genuine **uncredited** headroom — the
σ_θ **Tracker budget**'s 3 µrad is the *absolute* bench residual; today the differential is invoked
only to cancel bus *smear* and to justify cross-detector independence, never to cut the distortion
floor itself. **Payoff is spectrum-contingent** (no bench data — pure sim): a *low-spatial-frequency*
(smooth) residual leaves a gradient ≪3 µrad over the nearest-star separation → a several-× win; a
*high-spatial-frequency* (pixel-scale) residual makes differencing two uncorrelated errors √2 *worse*.
So the deliverable is a **sensitivity curve** (differential residual vs. assumed distortion correlation
length, break-even marked) that *outputs a bench-characterization requirement* (a paper §-worthy
result), not a point claim. **Merged with the Tracker array audit** (#1+#4 = one distortion-field
study: one field model, the spatial correlation length → per-detector differential gain, the
detector-to-detector correlation → cross-detector common-mode fraction). The VLBI-swarm
synthesized-baseline "beat" (idea #3) is **relocated to the paper's solar-collision chapter**, not
modeled for LEO — a compound-failure-only hedge whose cm cross-track baseline-knowledge requirement
re-encounters the lateral-blind wall (relative-ranging grill 2026-06-15). _Avoid_: treating distortion
as a scalar *floor* when crediting this (it must be a spatial *field*); claiming the win without naming
the assumed spectrum; the co-flyer-beacon-as-reference variant (#5) as a primary path (dominated by
stars — dense, µas-known, anchor-floor-free; keep only as a star-starved fallback).

**Co-flying tracker**:
The reused launch rocket (ADR 0019) as a *close* terminal tracker (~500 km vs the target's
2603 km → 5× less σ_θ·R), attacking the early large-R noise at its source. It tracks the **Train**
*centroid* (shrinks the common early error; per-unit scatter stays with the target **Tracker
array**). Load-bearing: the rocket→target vector must be pinned *independently* (unlocked
spaceborne GNSS — the terminal phase is low-altitude, inside the GPS volume), not by inter-rocket
ranging (which fixes range, not the long-baseline lateral). Gated on a **phasing** sim (can it stay
close + low + non-decaying at the terminal window). Hedge reframe (2026-06-15): the co-flyer is
**anchor-floored** — its `COFLYER_RELGEOM_SIGMA_M` 2 m GNSS-pinned floor dominates its σ_θ·R
(≈0.99 m at 500 km), so a 20× range cut buys only ~10 %; its real hedge value is being a
**σ_θ-independent ~2 m / 0.77 µrad backstop** (the anchor is GNSS-set, *not* distortion-set, so it
survives a bad distortion floor), not "push it closer." The only true linear-in-R lever — a *later
target hand-off* (no anchor on the target leg) — is dropped: it spends the thin ~2.2× **Catch
radius** margin *quadratically* for *linear* precision. _Avoid_: "they range to each other so they
know the relative vector" (ranging ≠ lateral); "push the co-flyer closer for ~20× gain"
(anchor-floored — ~10 %).

**Effective σ_θ**:
The *fused* relative-nav grade the PuffSat's ZEM loop actually sees (ADR 0019) — the
inverse-variance combination of every **Tracker array** / **Co-flying tracker** measurement at its
own range — as opposed to a single tracker's σ_θ. The D1.1 capture criterion reads against the
effective σ_θ; the σ_θ **Tracker budget** is *per detector*, fusion on top. _Avoid_: conflating the
per-detector σ_θ with the effective (system) σ_θ.

**Apogee nav constellation** (Lever 3):
A permanent ~150,000 km nav infrastructure (ADR 0020) that pins the **coast/apogee-state** nav the
**Differential corrector** consumes — a *different regime* from terminal homing, and one GNSS cannot
serve (GPS at 20,200 km vs apogee at 150,000 km). It generalizes the **Coordinator node** into
permanent, well-characterized *snapshot* GDOP at apogee. Signal architecture (sized pure in
`apogee_nav.py`): **Ka-band authenticated broadcast** (no atmosphere → no ionospheric term + free
high frequency for velocity sensitivity; not L-band GNSS, not laser — laser is point-to-point, wrong
for simultaneous multilateration and kept terminal-phase), **omni PuffSat / gain on the infra**
(ADR 0011 dec-7). The binding axis is apogee *transverse velocity*; a coplanar ring covers it (the
normal axis matters ~50× less, C0), a shell adds the weak axis. _Avoid_: "it's just GNSS at altitude"
(GNSS can't reach apogee, and the regime/geometry differ); conflating it with terminal **Tracker**
nav (it is the coast regime).

**Match-not-beat** (apogee nav accuracy):
ADR 0020's accuracy decision: the **Apogee nav constellation** should *match* the C1 grade
(σ_Tvel ≈ 0.66 mm/s / ~140 m — the per-unit entry budget D1.1 flies), not push tighter. After
ADR 0019 the binding lever is the *terminal* fused **Effective σ_θ** (already 4.2× capture margin),
so tighter apogee nav is redundant; the only payoff — substituting for the fusion hedge so a bare
10 µrad detector passes — needs ~3× tighter (~45 m), a worse trade than the fusion it would replace.
Better-than-140 m from good GDOP is *free* independent margin, not a target. _Avoid_: reading the
constellation's value as a tighter number (it is snapshot GDOP + pinning the rockets).

**Authenticated broadcast / secure transponder**:
The ADR 0020 anti-tamper split. **Authenticated broadcast** — the constellation cryptographically
signs the one-way nav message (OSNMA-style), verified by a sub-gram ASIC on the **passive** PuffSat
(no transmitter) — defeats *spoofing*. The **secure transponder** — a two-way, key-determined,
nanosecond-exact turnaround delay (cryptographic *distance-bounding* + round-trip range) — rides the
mass/power-unconstrained **rockets**, not the PuffSat. *Anti-jam* is a separate mechanism: spread-
spectrum processing gain + directional antennas + the deep-space regime, **not** crypto. _Avoid_:
"crypto stops jamming" (crypto is anti-spoof; geometry + processing gain are anti-jam); putting the
two-way transponder on the PuffSat (transmit power, not mass, is why it stays one-way).

**Near-Sun optical nav** (the Q-switched-laser carrier choice; paper-side, grill 2026-06-29):
The near-Sun / Parker-periapsis chapter's nav carrier — **Q-switched 1064 nm laser** behind a
**narrow monochromatic filter** with **ns time-gating**, *walking back* an earlier UV-C "solar-blind"
idea. Solar-blind is a *terrestrial* trick (Earth's ozone darkens the sky <290 nm); in space there is
no ozone, the Sun emits its UV-C straight at you, so the advantage evaporates while UV-C keeps all its
hard-tech cost — meanwhile the real background (broadband K-corona + thermal-IR off the ~1500–2000 K
heat shield) is beaten by **narrow-band × ns-gate × high peak power** at *any* line, and Q-switching is
what supplies the peak power (the committed LEO `tracker_budget.peak_power_sweep` already treats beacon
power as *peak*/strobed). 1064 nm **reuses the LEO tracker wavelength/detector lineage**. Two layers,
two jobs. **(1) Local mm lever = two-way *pulsed* laser ranging** between the converging projectiles +
the **Coordinator node** — self-clocking (no onboard precision clock, same rationale as the
**Transponder** term, which *avoids* one-way for exactly that reason), scintillation-robust (a pulse
punches through where a continuous lock drops), one pulse = range + optical two-way time-transfer
(pins the turnaround delay as a bench constant) + blink-code ID. **(2) Inertial frame = one-way
artificial-star beacons at ~1 AU** that broadcast their (Earth-pinned, km-known) positions; the
corona-blinded star tracker is replaced by matching their angle-of-arrival pattern to the broadcast
ephemerides, needing only **bearing-grade (~arcsec)** — mm rides on layer (1), not absolute attitude.
RF-vs-optical, stated precisely: coronal plasma hurts RF as ~1/f² but **path-integral**, so the **long
radial** 1 AU links are devastated (optical *mandatory*) while the **local transverse** link's bulk
slab is **common-mode, differential-cancelled** (optical a *robustness upgrade*) — sharper than the
chapter's current "tenuous, dual-frequency-removable." Two open *sizing* numbers (not forks): in-band
thermal-IR rate at the heat-shield temperature, and coronal optical angle-of-arrival jitter on the
radial path (both expected to close). Paper-side, reversible — no ADR. _Avoid_: one-way for the *mm*
lever (clock-limited — that is what the **Transponder** term avoids); letting the 1 AU beacons carry
mm *relative-range* load (GDOP/clock — they are the bearing frame only); "replace *all* RF ranging"
(overreach — optical is mandatory on the long links, an upgrade on the local one); keeping UV-C for a
"solar-blind" benefit that does not exist off-Earth.

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
