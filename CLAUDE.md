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
- `puffsat_sim/propagator.py` — `build_propagator()` (element-based) and `build_propagator_from_orbit()` (state-based seam for the MC harness); attaches force models.
- `puffsat_sim/truth_model.py` — `make run` report runner: reference orbit + per-force signatures.
- `puffsat_sim/dispersion.py` — pure MC core: `DispersionSpec`, `RunInputs`, `sample_run_inputs`, RTN math, `summarize` (no JVM).
- `puffsat_sim/montecarlo.py` — JVM-side open-loop dispersion harness: `run_ensemble` (`make capstone`); §14.1 `control=None` hook for Rung D (ADR 0002).
- `tests/` — pure-Python unit suite (no JVM); `tests/integration/` requires a live JVM.

## License

Copyright (c) 2026 Seth Katz. All Rights Reserved.
