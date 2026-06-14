# CLAUDE.md — PuffSat Control Simulation

Context for AI-assisted development in this repository.

## What this project is

Closed-loop simulation to determine whether a single PuffSat can navigate from
deployment near apogee (~150 000 km altitude) to interception at perigee (200 km
altitude) with enough precision to hit a pusher plate on a target rocket.

The science it validates is in **`puffsat_control_sim_design.md`** (read that
before touching any physics or control code) and in the companion paper
*Aim Is All You Need: A Speculative White Paper on PuffSat Pulsed Propulsion*
by Seth Katz (https://doi.org/10.5281/zenodo.16741183).

The LaTeX source and design document for the paper live at
https://github.com/katzseth22202/Balloon-Pulse-Propulsion — the
`puffsat_control_sim_design.md` there is the authoritative near-term feasibility
and control algorithm specification that drives this simulation.

## What the simulation does

Three decoupled pieces in a closed loop:

1. **Truth model** — Orekit high-fidelity propagator: full geopotential, SPICE
   Sun/Moon, NRLMSISE atmospheric drag with stochastic F10.7/Ap, SRP with
   eclipse and attitude (cannonball first pass). Runs in short arcs, stops at
   each maneuver/regime boundary.

2. **UKF estimator** — estimates position, velocity, and the lumped drag/SRP
   coefficients `Cd·(A/m)` and `Cr·(A/m)`. Uses altitude-scheduled sensor
   rates: ~0.03 Hz in coast, 100 Hz in terminal.

3. **MPC controller** — discrete midcourse corrections during coast; continuous
   terminal burn 600 → 200 km for drag rejection and final aim.

Output: Monte Carlo distribution of interception miss at 200 km, perigee altitude,
and propellant consumed. The mission-killing event is **failing the 200 km
interception** (missing the pusher plate) — that is where the mission succeeds, by
transferring momentum to the target rocket (paper §2). A low perigee (~50 km) is
*intended*: it deorbits PuffSat dry mass and burns up any PuffSat that misses, for
debris disposal (paper §9, `sec:handling_space_debris`) — so burn-up is the desired
outcome on a miss, **not** a failure. Perigee is therefore a diagnostic (the §8
`dr_p/dv_a` lever, and a debris-disposal-safety margin where *low is good*), not a
success criterion. The paper claims <2% of 25 kg PuffSat mass for propellant; the
sim tests that.

## Architecture decisions (do not relitigate without reading the design doc)

- **Orekit, in-process JVM via JPype.** Not a gRPC server. Not tudatpy (good
  alternative but Orekit chosen for pedigree and event-detection framework).
  Not REBOUND (wrong strength: N-body, not GNC).
- **Segmented propagation, not one continuous arc.** Adaptive integrator (DOPRI8
  or IAS15) for coast; fixed-step Cowell for terminal phase. Event-restart at
  each maneuver/regime boundary.
- **No symplectic integrators.** Impulsive thrust + dissipative drag destroy the
  Hamiltonian structure.
- **MPC stays in Python.** Neural warm-start (second pass) and RL anomaly
  recovery (second pass) live here naturally.
- **Cannonball SRP/drag model for first pass.** Attitude-dependent cylinder model
  deferred to a later rung.
- **Classical UKF + MPC first, neural augmentation second.** The neural stuff
  (KalmanNet-style estimator, transformer warm-start, RL recovery) is a second
  pass; the first pass is clean classical baselines.

## Orekit Python specifics

The conda-forge package is `orekit_jpype` (JPype-based), not `orekit` (JCC-based).
Only `orekit_jpype` is available for aarch64/ARM.

### Initialization order (must not be changed)
The JVM boot lives in `puffsat_sim/jvm.py`; importing it starts the JVM and loads
`orekit-data.zip` (once — Python caches the module).  The one ordering rule in the
codebase: **import `puffsat_sim.jvm` before any `org.orekit` import.**
```python
import puffsat_sim.jvm   # boots the JVM; must precede any org.orekit import
from org.orekit.frames import FramesFactory
```
`jvm.py` itself does:
```python
import orekit_jpype
_VM = orekit_jpype.initVM()   # start JVM — must precede any org.orekit import
from orekit_jpype.pyhelpers import setup_orekit_curdir
setup_orekit_curdir()         # load Earth/time data from orekit-data.zip in cwd
```

`orekit-data.zip` must exist in the working directory (project root when using make).
Download once with `make data`; it fetches from `gitlab.orekit.org`.

### Type annotations
Orekit ships no `.pyi` stubs. All JVM-side objects surface as `Any`.
Suppress at the import boundary with `# type: ignore[import-untyped]`.
Use `ignore_errors = true` in the mypy overrides for `org.*` / `orekit.*`
(already configured in `pyproject.toml`). Our own code must be fully typed.

### Common JPype friction
- Java method overloads: JPype resolves by argument types; when it ambiguates,
  cast explicitly using Java primitive wrappers (`JDouble`, `JInt`).
- Java exceptions surface with cross-boundary stack traces; look for the
  `java.lang.*Exception` line in the traceback.
- `AbsoluteDate` / `TimeScalesFactory` require the data file loaded by
  `setup_orekit_curdir()` before first use.

### Orekit data file
`orekit_jpype` does NOT bundle `orekit-data.zip`; it must be downloaded separately.
If `setup_orekit_curdir()` raises `FileNotFoundError`, run `make data` from the
project root, or manually:
```python
import orekit_jpype; orekit_jpype.initVM()
from orekit_jpype.pyhelpers import download_orekit_data_curdir
download_orekit_data_curdir()   # downloads to current directory
```

## Reference orbit parameters

| Parameter | Value | Source |
|---|---|---|
| Orbit periapsis | 50 km | Below Kármán line — debris disposal by reentry |
| Interception altitude | 200 km | During descent, before periapsis; Paper §3, design doc §3 |
| Apogee altitude | ~150 000 km (from surface) | Design doc §3 recommendation |
| Semi-major axis | ~81 403 km | Derived |
| Eccentricity | ~0.921 | Derived |
| Inclination | 28.5° (nominal) | Mid-latitude launch |
| Speed at interception | ~10.78 km/s | Near escape speed at 200 km |
| Orbital period | ~2.68 days | Derived |
| Post-impact burnup | ~120 km onset, ~50 km complete | Design doc §6.3 |

## Build rungs (from design doc)

- **Rung A (now):** single-PuffSat truth propagation with force models.
  Verify reference orbit. No control yet.
- **Rung B:** add UKF estimator. Verify state estimation on truth.
- **Rung C:** add MPC controller. Verify closed-loop, single trajectory.
- **Rung D:** Monte Carlo (N=10³–10⁴). Measure perigee/miss/propellant
  distributions. This is the result.

> **Naming note:** "Rung" is overloaded in the design doc (a *physics* ladder
> 0/1/2a–2d and a *control* ladder A–D). Force-model presets are therefore named
> by **content** (`presets.two_body`, `j2`, `j2_third_body`, `j2_third_body_srp`,
> `full_force`), never by rung number. See `CONTEXT.md` and
> `docs/adr/0001-pure-perturbation-specs.md`.

## Coding conventions

- Python 3.11+; strict mypy; ruff for lint and formatting. `ruff format` is
  enforced — `make all` fails on unformatted code, so run `make format` before
  committing. (Do not hand-align inline comments; ruff collapses them.)
- Every public function must have full type annotations, including return type.
- No comments explaining *what* the code does. Only comment *why* when the
  reason is non-obvious (hidden constraint, workaround, subtle invariant).
- No orekit-specific logic in `tests/`; unit tests cover pure Python helpers.
  Integration tests (requiring a live JVM) go in a separate `tests/integration/`
  directory when they exist.
- `# type: ignore[import-untyped]` on orekit imports — do not silence other mypy
  errors with bare `# type: ignore`.
- God-file/function tripwires (ADR 0017, design doc §14.5): a module with more than
  one owning ADR in its CLAUDE.md entry gets split (lines secondary: suspicious ~400
  logic, act ~600); cut on Seams (CONTEXT.md) — ownership says when, seams say where.
  JVM glue mirrors pure slices in `puffsat_sim/runs/`. Functions: suspicious ~50
  branching logic lines (~80 linear recipes); ≥3 default-off behavior knobs → bundle
  into a value object.

## Environment

Full setup — see README.md for the detailed version.

```bash
# 1. Install Miniconda if needed (pick the right arch):
#    https://docs.conda.io/en/latest/miniconda.html

# 2. Create the conda environment (Python 3.11 + orekit_jpype + mypy + ruff + pytest)
conda env create -f environment.yml
conda activate puffsat-sim

# 3. Download orekit-data.zip once (fetches from gitlab.orekit.org, ~37 MB)
make data
```

No pip install step — flat layout, Python finds `puffsat_sim/` directly from the repo root.
`orekit-data.zip` is gitignored; `make run` guards against it being missing.

## Common tasks (Makefile)

```bash
make all       # mypy + lint + format-check + test
make run       # truth-model report: reference orbit + per-force signatures
make capstone  # open-loop dispersion capstone (smoke N=50; puffsat_sim.montecarlo)
make test      # pytest
make mypy      # strict type check
make lint      # ruff check
make format    # ruff format (auto-fix)
make format-check  # ruff format --check (gate; part of make all)
make clean     # remove caches
```

## Files

- `puffsat_control_sim_design.md` — authoritative design doc. Read before any
  physics/control work.
- `environment.yml` — conda environment (no pip).
- `pyproject.toml` — tool config only: mypy, ruff, pytest.
- `Makefile` — task runner.
- `CONTEXT.md` — domain glossary (Perturbation, Force Model, Preset, Environment).
- `docs/adr/` — architecture decision records.
- `puffsat_sim/jvm.py` — the JVM boot seam: import it before any `org.orekit` import.
- `puffsat_sim/constants.py` — single source of truth for scalar physical constants (pure).
- `puffsat_sim/config.py` — `OrbitalConfig` and `PhysicsConfig` (a `tuple[Perturbation, ...]`); pure Python.
- `puffsat_sim/presets.py` — named, content-described `PhysicsConfig` bundles (pure).
- `puffsat_sim/forces/` — one pure module per perturbation (spec + analytic signature);
  `forces/build.py` is the JVM side (`Environment` + `to_force_models()` dispatch).
- `puffsat_sim/orbital_math.py` — foundational two-body helpers only (pure).
- `puffsat_sim/orbital_plane.py` — `orbital_config_from_cities()` great-circle plane builder (pure).
- `puffsat_sim/mission.py` — reference scenario: altitudes, epoch, `NOMINAL_CONFIG` (pure, single source).
- `puffsat_sim/propagator.py` — `build_propagator()` (element-based), `build_propagator_from_orbit()` (state-based seam for the MC harness), and `build_fixed_step_propagator_from_orbit()` (fixed-step Cowell terminal for the executed burn, ADR 0014); attaches force models.
- `puffsat_sim/truth_model.py` — `make run` report runner: reference orbit + per-force signatures.
- `puffsat_sim/dispersion.py` — pure MC core: `DispersionSpec`, `RunInputs`, `sample_run_inputs`, `replay_inputs` (§14.2 seed-replay), `lognormal_factor`, RTN math, `summarize`, `EnsembleStats` (no JVM).
- `puffsat_sim/train.py` — pure train-mode dispersion core, the first Rung-D / D1 sub-slice (ADR 0016/0018, no JVM): the shared-vs-per-unit draw split — `TrainDispersionSpec` separates each `DispersionSpec` σ into a **shared** (per-train: coefficient bias / F10.7-Ap drivers / deployer systematic; the §16.7 common-density component) and a **per-unit** part (coefficient spread / injection scatter); `sample_train` composes the same `RunInputs` the JVM `run_record` consumes (so D1.1 wires in unchanged), `replay_train_unit` is the standalone §14.2 replay (shared from `train_index`, per-unit from `(train_index, unit_index)`); plus the train-relative reduction — `summarize_train_capture` → `TrainCaptureStats` splits a train's `guidance.PlateMiss` arrivals into **centroid drift** (vs the ±2 km `CENTROID_RETARGET_M`, `retarget_ok`) and **scatter about the centroid** (per-axis σ vs `guidance.CAPTURE_SIGMA_MAX_M`, `scatter_sigma_ok`), reporting `capture_about_centroid` (reusing `guidance.capture_fraction`) vs `capture_absolute`; `format_train_capture`. **D1.1 extension:** the Φ-composed hand-off entry offset — `sample_train_entry_offsets` / `replay_train_entry_offset` draw a per-unit 2-D lateral (⊥v) entry = shared centroid + per-unit scatter, magnitudes the characterized C0/C1/C2a budget legs (`ENTRY_LATERAL_PERUNIT_M` 141 m / `ENTRY_LATERAL_SHARED_M` 149 m, sampled isotropically at σ/√2; seed tree masked off the coefficient tree by `_ENTRY_SEED_MASK`); and the verdict surface `summarize_train_ensemble` → `TrainEnsembleFinding` (the D1.0 `TrainCaptureStats` + terminal-aim Δv `propellant.propellant_curve` vs <2 % + perigee diagnostic) / `format_train_ensemble`.
- `puffsat_sim/control.py` — pure Rung A1 differential corrector: `Target`/`ControlAction`/`ControlPlan` + `solve_apogee_correction` (Newton + finite-difference Jacobian; no JVM) (ADR 0003); `report_controller` is the LM-tuned report-grade config the C-rung reports share.
- `puffsat_sim/records.py` — pure result value types `RunRecord` / `EnsembleResult` (no JVM, so the resume sink stays serializable/testable) (ADR 0003).
- `puffsat_sim/sink.py` — pure JSONL resume sink: `record_to_dict`/`record_from_dict`, `append_record`/`read_records`, `plan_resume` (resume-by-complement; no JVM) (ADR 0003).
- `puffsat_sim/estimation.py` — pure C1 estimation core (ADR 0012, no JVM): owned typed UKF (sigma points, unscented transform, predict/update), two-body+J2 onboard filter dynamics, measurement models (range/LOS-Doppler/GNSS) with `NodeState` known-ephemeris beacons, white-acceleration `Q`, LinCov recursion (`run_lincov`), NEES consistency bounds.
- `puffsat_sim/nav_feasibility.py` — pure C1 sweep harness (ADR 0012, no JVM): one-axis-at-a-time `NavFeasibilitySpec` grid (range σ / Doppler σ incl. range-only / cadence / cone geometry / n_nodes / Q), per-cell LinCov along the coast arc to apogee, Σ→RTN, `ΦΣΦᵀ` catch-radius verdict via C0's Φ; `sweep_nav_feasibility` + `format_nav_feasibility`; layer-2 `validate_cell` (seeded UKF truth run judged by NEES) + `format_nav_validation`.
- `puffsat_sim/coeff_requirement.py` — pure C2a coefficient-knowledge requirement core (ADR 0013, no JVM): 1D A3 cut → per-cell Δv vectors (`cut_dv_vectors`) → `dv_gradient` at nominal → lateral sensitivity through C0's Φ → `coefficient_tolerance` vs the ground prior (`summarize_coeff_requirement`/`format_coeff_requirement`); analytic SRP-impulse cross-check (`analytic_srp_dv`) and the RSS lateral error-budget ledger (`BudgetEntry`/`rss_lateral`/`MEASURED_BUDGET`).
- `puffsat_sim/terminal.py` — pure C3a terminal feedforward core (ADR 0014, no JVM): `plan_feedforward` realizes B3a's anti-drag profile as zero-order-hold `ThrustCommand`s on the control clock (thrust `m·|a_drag|` anti-parallel to drag, saturated at the 400 mN actuator cap); `FeedforwardPlan` reports delivered Δv, saturation, peak thrust, and peak slew rate (negligible-thrust commands excluded as direction noise); `TerminalFeedforwardFinding` + `format_terminal_feedforward` render the C3a report (equivalence pin, drag displacement vs executed residual + rejection ratio, ADR 0004 gate verdicts, propellant curve).
- `puffsat_sim/guidance.py` — pure C3b terminal-guidance core (ADR 0014/0015, no JVM): ZEM law (`zem_acceleration`, `predicted_zem` on the two-body+J2 onboard model), the noise-discipline tick (`terminal_tick`: significance gate / track window / firing-lag hold — see the measured constants block), `TrackerGrade` σ_θ·R nav noise as a Gauss–Markov `NavNoiseProcess`, `slew_limited_direction`, plate-frame miss (`plate_frame_miss`/`PlateMiss`/`capture_fraction`), sweep value types (`GuidanceSweepSpec`/`GuidanceCell`/`TerminalGuidanceFinding`), `measured_catch_radius_m`, `format_terminal_guidance`.
- `puffsat_sim/authority.py` — pure C3c terminal-authority / tail-correction core (ADR 0014 decision 5/6, no JVM): the funnel model (`thrust_limited_radius_m` = ½·a_max·t², `saturation_dv_m_s` = a·t, `funnel_growth_dv_m_s` = √(2·a·entry) saturation-edge cost; all C3b-validated), the MCC-2 trim lever (`lateral_lever_m_per_m_s` central-difference ⊥v projection → `dv_per_km_m_s`), curve value types (`AuthorityPoint`/`TrimPoint`/`TailAuthoritySpec`/`TailAuthorityFinding` carrying the full C2a `budget_entries`), the crossover reads (`handoff_alt_to_cover_tail_m`, `cheapest_trim`, `trim_dv_to_cover_tail_m_s`, `tail_m`/`uncovered_tail_m`), `format_tail_authority`.
- `puffsat_sim/latency.py` — pure C4 control-loop latency core (ADR 0014, §16.8, **no JVM** — dead-time is a loop-transfer effect, so Orekit adds nothing; the only C-rung with no `runs/` glue): the per-loop dead-time budget (`LatencySource`/`ControlLoop` with `tau_s` + `phase_margin_loss_deg` = ω_c·τ, `None` for the discrete midcourse loop), `TERMINAL_LOOP` (7.3 ms) / `MIDCOURSE_LOOP` (70 ms, discrete), the `DeadTimeBuffer` (equivalent e^{-sτ}), `fly_terminal_loop` (the C3b ZEM law on a double integrator with the stale-fix buffer; noiseless by default, `rng=` for the deferred Rung-D combined-stress case), `tau_sweep`/`TauSweepPoint`, `LatencyFinding` (relative-degradation reads: `tolerated_latency_s`/`breakdown_latency_s`/`budget_margin`), `latency_finding` (the pure runner) + `format_latency`.
- `puffsat_sim/tracker_budget.py` — pure σ_θ tracker-budget gate (ADR 0018, the first Rung-D blocking pre-gate; **no JVM**, like `latency.py` — angular precision is a focal-plane question, not an orbit one, so no `runs/` glue): a four-term RSS σ_θ error budget for a declared `TrackerHardware` point — `photon_sigma_theta_rad` (active beacon → huge SNR, negligible), `smear_sigma_theta_rad` (post-impact streak after differential cancellation), `gyro_bridge_sigma_theta_rad`, and the bench-calibratable focal-plane `distortion_floor_rad` (the dominant term) — plus acquisition (`required_fov_halfangle_rad` vs the C1 hand-off Σ, `reference_star_fov_halfangle_rad` the binding FOV, `detector_pixels_across`); `TrackerBudgetFinding` ties the achieved grade to C3b's `homing_floor_m` and the ADR 0015 capture criterion (`meets_requirement`/`meets_target`/`capture_floor_met` are the D1 entry condition); `tracker_budget_finding` (pure runner) + `format_tracker_budget`.
- `puffsat_sim/tracker_fusion.py` — pure multi-tracker fusion gate (ADR 0019, the early-error reduction; **no JVM**, builds on `tracker_budget.py`): D1.1 found the combined entry×noise stress needs an **effective** ~3 µrad terminal-nav grade; this quantifies how a multi-tracker architecture recovers it from cruder, redundant 10 µrad detectors by two levers — **averaging** (`angular_sigma_theta_rad` = √(σ_common² + σ_indep²/N): the independent part — distortion⊕gyro⊕photon, separately calibrated — averages with √N, the smear `SMEAR_COMMON_SIGMA_RAD` common-mode floor does not) and **range** (`lateral_sigma_m` = σ_θ·R ⊕ rel-geom floor; a closer `coflyer` at `COFLYER_RANGE_M` 500 km vs the target's `TARGET_RANGE_M` 2603 km has proportionally less lateral error, plus the GNSS-pinned `COFLYER_RELGEOM_SIGMA_M` rocket→target vector); `Tracker` value type, `target_array(n)` / `coflyer(n)` constructors, `fuse_lateral_sigma_m` (inverse-variance) → `effective_sigma_theta_rad` at the design range read against the D1.1 `D1_CAPTURE_GRADE_SIGMA_THETA_RAD` (3.2 µrad); `TrackerFusionFinding` (`meets_capture_grade`/`margin`/`homing_floor_m`) + `tracker_fusion_finding` (pure runner) + `format_tracker_fusion`. Measured: target 5-array alone → 1.62 µrad (2.0× inside, no phasing dependency); +co-flyer → 0.76 µrad (4.2×). **Stage 2 — the Lever-2 co-flyer phasing gate** (pure verdict; the JVM run is `runs/coflyer.py`): `CoflyerPhasing` / `phasing_verdict` reduce the rocket↔centroid range and rocket altitude sampled over the 800→200 km terminal window to a two-part feasibility — `range_ok` (peak separation ≤ the `COFLYER_RANGE_M` angle-useful design range, so the σ_θ·R lever the fusion credits actually holds) and `gps_ok` (peak rocket altitude ≤ `GPS_CEILING_M` 20 200 km, inside the unlocked-spaceborne-GNSS volume that pins the rocket→target vector independently of the long baseline) → `feasible`; `format_coflyer_phasing`. **The re-key (ADR 0019 decision 4)**: `fused_tracker_grade` collapses an architecture into the C3b/D1.1 loop's `guidance.TrackerGrade` (effective σ_θ + unchanged ranging σ — the noise model already consumed a scalar σ_θ, so fusion adds no noise code); the canonical architectures `single_target_detector` / `target_array_only` / `array_with_coflyer` (`ARRAY_N_DETECTORS` 5, `COFLYER_N_DETECTORS` 3) feed the D1.1 re-run.
- `puffsat_sim/torque_margin.py` — pure torque-margin confirmation gate (ADR 0018, a Rung-D pre-gate; **non-blocking**, **no JVM** like `tracker_budget.py`/`latency.py` — attitude agility is an inertia/actuator question, not an orbit one): confirms the C3b 1 °/s direction-loop rail (`anti_drag.PEAK_SLEW_LIMIT_DEG_S`) carries margin over the thrust-direction demand — the perigee LOS rate `v_p/r_p` (`perigee_los_rate_rad_s`, ~0.1 °/s) and the B3a-measured 0.048 °/s — and that a conservative whole-body actuator (`inertia_from_gyradius`, `couple_torque_n_m`, `angular_accel_deg_s2`) reaches the demand inside one control period and out-torques the aero disturbance (`aero_disturbance_torque_n_m`); `TorqueMarginFinding` reports the rate / agility / disturbance margins and the **break-even** inertia and control torque (the paper-side pins); `torque_margin_finding` (pure runner) + `format_torque_margin`.
- `puffsat_sim/truth_validation.py` — pure truth-model validation core (ADR 0018, a Rung-D pre-gate; **non-blocking**; no JVM *in this module*, but unlike `tracker_budget.py` it does have a `runs/` glue — the coast it validates is flown by Orekit): two tiers over the apogee→hand-off coast (~99 % of the trajectory) — **Tier 1** integrator health on a *numerical* two-body coast (`conservation_drift`: `specific_energy_j_per_kg` and `angular_momentum_magnitude` are constants of motion, so `max_fractional_drift` is pure numerical leak) plus tolerance-halving; **Tier 2** an `independent_coast` (chained `estimation.two_body_j2_flow` RK4 Cowell, sharing only pinned constants) cross-checked against the Orekit J2 coast via `max_position_divergence_m`; `TruthValidationFinding` (`conservation_ok`/`convergence_ok`/`crosscheck_ok`/`validated`) + `format_truth_validation`.
- `puffsat_sim/descent.py` — truth-path kernel (ADR 0017): the §6.2 regime-switched descent every JVM harness consumes — `descend` (coast 600 s → 800 km hand-off → terminal 30 s), `propagate_to_interception` → `Crossing`, `coast_to_altitude` (the §6.2 coast generalized to any descending altitude event) / `coast_to_handoff`, `apogee_state`, `earth_model`/`to_absolute_date`.
- `puffsat_sim/montecarlo.py` — the predict/execute dispersion harness (ADR 0002; control hook ADR 0003): public surface `build_context` → `RunContext`, `run_record` + `RunVariant` (the bundled rung knobs: control / toa-gate / actuator / nav-offset), `run_ensemble` with `sink_path` checkpoint/resume, `physics_from_inputs`/`nominal_inputs`; `make capstone`.
- `puffsat_sim/runs/` — JVM glue per slice, names mirroring the pure cores across the master seam (ADR 0017): `sweep.py` (`run_sweep`, A3 / ADR 0007); `navigation.py` (`run_nav_sweep` + `nav_requirement_report`, C0 / ADR 0011); `nav_feasibility.py` (`truth_arc_to_apogee` + `run_nav_feasibility` + `nav_feasibility_report`, C1 / ADR 0012); `coeff_requirement.py` (`coeff_requirement_report`, C2a / ADR 0013); `anti_drag.py` (`instrument_anti_drag` + `sample_drag_window`, B3a / ADR 0008); `terminal.py` (`run_terminal_feedforward` + `terminal_feedforward_report`, C3a / ADR 0014); `guidance.py` (`build_guidance_context` + `run_guidance` tick loop + `run_terminal_guidance` sweep + `terminal_guidance_report`, C3b / ADR 0014/0015); `authority.py` (`build_authority_context` + `measure_authority_point` / `measure_trim_point` + `run_tail_authority` + `tail_authority_report`, C3c / ADR 0014); `truth_validation.py` (`run_truth_validation` + `truth_validation_report`, Rung-D pre-gate / ADR 0018 — flies the reference coast three ways and feeds the pure Tier 1 / Tier 2 checks); `train.py` (`run_train_dispersion` + `train_dispersion_report`, Rung D / D1.1 / ADR 0018 — flies a train's units through the C3b terminal loop with Φ-composed entry offsets + tracker noise + dispersed drag, reduces via `summarize_train_ensemble`; `trackers=` re-keys the noise grade to a fused architecture and `fused_train_rerun_report` re-runs D1.1 across them, ADR 0019 decision 4 — the 10 µrad ceiling fails, fused 5-array/+co-flyer recover capture-grade); `coflyer.py` (`run_coflyer_phasing` + `coflyer_phasing_report`, ADR 0019 Stage 2 — flies the constant-a / phase-locked launch rocket alongside the descent and samples its range-to-centroid and altitude over the 800→200 km window to feed `tracker_fusion.phasing_verdict`).
- `tests/` — pure-Python unit suite (no JVM); `tests/integration/` requires a live JVM.

## License

Copyright (c) 2026 Seth Katz. All Rights Reserved.
