# PuffSat Control Simulation

A closed-loop orbital mechanics simulation for **PuffSat pulsed propulsion** — the concept proposed in the paper [*Aim Is All You Need: A Speculative White Paper on PuffSat Pulsed Propulsion*](https://doi.org/10.5281/zenodo.16741183) by Seth Katz.

## What this is

A PuffSat is a single-use microsatellite (~25–100 kg) that converts nearly all its mass into a high-velocity gas puff. A sequence of PuffSats deployed into a highly eccentric Earth orbit sequentially impacts a pusher plate on a target rocket, transferring momentum without requiring the target rocket to carry most of its own propellant.

This repository simulates whether that actually works. Specifically, it answers: **can a PuffSat navigate from deployment near apogee (~150 000 km altitude) to interception near perigee (200 km altitude) with enough precision to hit a pusher plate, under realistic orbital perturbations?**

The simulation is a closed-loop Monte Carlo dispersion study:

```
Truth model (Orekit, high-fidelity)
    → UKF state estimator
    → MPC trajectory controller
    → impulsive maneuver commands
    → back into truth model
```

Output metrics per Monte Carlo run: perigee altitude achieved (mission-killer if below ~120 km), interception miss distance, and cold-gas propellant consumed (the paper claims <2% of PuffSat mass).

For the full design rationale — force models, propagation strategy, controllability analysis, tooling decisions — see [`puffsat_control_sim_design.md`](puffsat_control_sim_design.md).

## Companion paper & near-term feasibility

The LaTeX source for the companion paper, along with the authoritative control-algorithm design document, lives at **[katzseth22202/Balloon-Pulse-Propulsion](https://github.com/katzseth22202/Balloon-Pulse-Propulsion)**. The `puffsat_control_sim_design.md` in that repo is the specification this simulation implements.

### Why the orbit is controllable (near-term feasibility case)

The design document establishes controllability around one pivotal choice: **apogee at ~150 000 km instead of the 0.9 Hill radius (~1.35 M km)**. That single decision:

- Reduces perigee-altitude sensitivity from ~250 km/m/s to ~30 km/m/s, making injection errors tolerable
- Cuts solar-tidal perturbations from ~49% to ~0.1% of local gravity, removing chaos from the problem
- Shrinks SRP-driven dispersions by roughly an order of magnitude through a shorter coast arc

With those margins, a 25 kg PuffSat carrying ~400 g of cold gas (Isp ~200 s, ~32 m/s total Δv) can execute two midcourse corrections (~1 m/s each, buying hundreds of kilometers of perigee adjustment) plus a terminal continuous burn for drag rejection and final aim — and still stay within the paper's <2% propellant claim.

### Control algorithm structure

The algorithm decomposes forces by feedback bandwidth across three tiers:

| Tier | Forces | Mechanism |
|---|---|---|
| Feedforward (open-loop) | Geopotential, J2+, third-body Sun/Moon | Recomputed from current state estimate |
| Slow feedback | Solar radiation pressure, coast divergence | `Cr·(A/m)` as UKF state; discrete midcourse burns |
| Fast feedback | Atmospheric drag near perigee | `Cd·(A/m)` as UKF state; continuous terminal burn |

Four clocks are strictly separated (conflating them caused errors in earlier framing):

```
Integrator step     — adaptive, non-uniform (minutes → sub-second)
UKF update          — altitude-scheduled (0.03 Hz coast → 100 Hz terminal)
MPC replan          — per-maneuver in coast; 1–10 Hz terminal
Inner tracking loop — up to 100 Hz (thrust modulation, not gimbal)
```

MPC replanning is **not compute-bound** — a warm-started QP solves in µs–ms, well inside any replan budget. The binding limits are information arrival rate and how fast the reference changes (~minute timescale in coast, ~second in terminal).

## Requirements

- **conda** (Miniconda or Anaconda) — [install Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you don't have it
- No separate Java runtime — `orekit_jpype` bundles its own JVM

## Setup

### 1. Install Miniconda (if needed)

```bash
# Linux x86_64
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh | bash -b -p ~/miniconda3
# Linux aarch64 (Apple Silicon, ARM servers)
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh | bash -b -p ~/miniconda3
# macOS (Apple Silicon)
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh | bash -b -p ~/miniconda3
# macOS (Intel)
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh | bash -b -p ~/miniconda3

# Add to shell (restart your shell or source the rc file after)
echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc   # or ~/.zshrc
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate puffsat-sim
```

This installs Python 3.11, `orekit_jpype` (the JPype-based Orekit bridge with a
bundled JVM), `mypy`, `ruff`, and `pytest`. No pip install step needed — flat
layout means Python finds `puffsat_sim/` directly from the repo root.

### 3. Download Orekit data (once)

Orekit requires Earth orientation parameters, leap seconds, and planetary
ephemerides from a ~37 MB data file that is **not** bundled in the conda
package. Download it once into the project root:

```bash
make data
```

This creates `orekit-data.zip` in the project root. It is gitignored. You only
need to do this once per clone; `make run` will remind you if it is missing.

`make data` fetches from `gitlab.orekit.org`. If that domain is blocked (e.g. in
a restricted network sandbox), allow it first:

```bash
# In a Sandbox environment only — run on your host machine:
sbx policy allow network -g gitlab.orekit.org
```

The `orekit` package from conda-forge bundles `orekit-data.zip` (Earth orientation parameters, leap seconds, etc.). `setup_orekit_curdir()` in the code locates it automatically. If you see a `FileNotFoundError` about `orekit-data.zip`, run:

```python
from orekit.pyhelpers import download_orekit_data_curdir
download_orekit_data_curdir()   # downloads to current directory; requires internet
```

## Run the Rung A truth model

Verify the Orekit / JVM bridge and the reference orbit parameters:

```bash
make run
```

or directly:

```bash
python -m puffsat_sim.truth_model
```

Expected output (numbers are exact for Keplerian propagation):

```
PuffSat Control Simulation — Rung A: Keplerian reference orbit
  Orekit / JVM : OK

  Reference orbit (near-term architecture):
    Orbit periapsis  : 50 km  (burns up here; interception at 200 km during descent)
    Apogee altitude  : 150 × 10³ km  (deployment)
    Semi-major axis  : 81403.1 km
    Eccentricity     : 0.921033
    Inclination      : 70.0°  (great circle through Tokyo and New York)
    Orbital period   : 231138.7 s  (2.68 days)
    Perigee speed    : 10.914 km/s

  One-period propagation residual (Keplerian → should be ~0):
    |Δr| = 7.856e-09 m
    |Δv| = 2.823e-13 m/s
```

## Common tasks

```bash
make all        # mypy + lint + test (CI equivalent)
make run        # Rung A truth model (reference orbit verification)
make test       # pytest
make mypy       # strict type check
make lint       # ruff check (subsumes flake8/isort/pyupgrade)
make format     # ruff format
make clean      # remove __pycache__, .mypy_cache, etc.
```

## Project layout

```
puffsat_control_simulation/
├── environment.yml                  # conda environment (orekit, mypy, ruff, pytest)
├── pyproject.toml                   # tool config: mypy, ruff, pytest
├── Makefile                         # run, test, mypy, lint, format, all, clean
├── LICENSE                          # All Rights Reserved
├── README.md
├── CLAUDE.md                        # context for AI-assisted development
├── puffsat_control_sim_design.md    # detailed design document (read this first)
├── puffsat_sim/
│   ├── __init__.py
│   ├── config.py                    # OrbitalConfig + PhysicsConfig dataclasses (pure Python)
│   ├── orbital_math.py              # Keplerian helpers + orbital_config_from_cities() (pure Python)
│   ├── propagator.py                # JVM boundary: initVM + build_propagator()
│   └── truth_model.py               # Rung A runner: reference orbit verification
└── tests/
    ├── __init__.py
    ├── test_config.py               # OrbitalConfig / PhysicsConfig unit tests
    └── test_orbital_math.py         # orbital mechanics + city-helper unit tests
```

## Tooling decisions

| Concern | Choice | Why |
|---|---|---|
| Environment | conda | `orekit_jpype` is only on conda-forge; pip-only installs are painful |
| Orbital mechanics | [Orekit](https://www.orekit.org/) (Java, driven from Python) | Most mature validated astrodynamics library; event-detection framework fits the segmented closed-loop architecture; in-process JVM via JPype, not a server |
| Python/JVM bridge | `orekit_jpype` (JPype, in-process) | Single-digit-microsecond per-call latency; no serialization; Orekit event callbacks work natively; available for all platforms including aarch64 |
| Type checking | mypy strict | Full annotations on our code; orekit JVM boundary suppressed with `ignore_errors = true` in overrides |
| Linting/formatting | ruff | Fast; handles isort, pyupgrade, bugbear, and annotations in one tool |
| Independent validation | GMAT (NASA GSFC) | Design and verify reference trajectories; cross-check Orekit truth model on the nominal trajectory |

## License

Copyright (c) 2026 Seth Katz. All Rights Reserved. See [LICENSE](LICENSE).
