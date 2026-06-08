# PuffSat Closed-Loop Control Simulation — Design Document

**Status:** Draft for review
**Date:** 2026-06-02 (revised 2026-06-03)
**Author context:** Seth Katz (paper author). Design distilled from grilling sessions on simulating the control algorithms in *Aim Is All You Need: A Speculative White Paper on PuffSat Pulsed Propulsion*.
**Intended home:** a separate implementation repository (not this LaTeX paper repo).
**2026-06-03 revision:** decomposed the build ladder into a physics stage and a control stage (Rung A detailed; B–D deferred until we reach them); recorded the in-process-JVM vs gRPC rationale and the boundary-crossing strategy; added the per-run log-normal coefficient uncertainty model and the cannonball idealization; noted pure-Python astrodynamics (boinor) as out of scope for the truth engine.

---

## 1. Purpose

Build a closed-loop simulation to answer one question: **how controllable is a PuffSat on its way to interception?** Concretely, quantify the distribution of:

1. **Interception success at 200 km** — whether the PuffSat reaches the 200 km interception altitude at the right position and time for pusher-plate impact. The orbit periapsis is 50 km by design (debris disposal; see Section 3), so the mission-killer is the PuffSat failing to reach the 200 km interception point intact and on-target, not the periapsis altitude per se.
2. **Interception miss** (along-track timing and lateral position relative to the target rocket's pusher plate).
3. **Propellant expended** for drag compensation and trajectory correction (the paper claims < 2% of a 25 kg PuffSat; this sim tests that).

This is fundamentally a **Monte Carlo dispersion study**, not a single-trajectory propagation. The result that supports or refutes the paper's thesis is a *distribution with a resolved tail*, produced by propagating many closed-loop runs under sampled uncertainty.

The integrator question that started the discussion (REBOUND IAS15 vs Cowell, symplectic applicability) turned out to be the least important decision in the stack. It is resolved in Section 6 and is not the design driver.

---

## 2. Scope

### In scope
- Single-PuffSat closed-loop guidance from deployment near apogee to interception near perigee.
- Truth dynamics + estimator (UKF) + controller (MPC, classical first) + impulsive/finite maneuvers.
- Monte Carlo over the dominant uncertainties.
- Orbit-level controllability (kilometre-to-metre scale: hitting the right perigee at the right time).

### Out of scope (first pass)
- **Terminal centimetre-scale centering** in front of the pusher plate (a nested, tighter sub-problem; follow-on once orbit-level control is shown feasible).
- **Target-rocket RCS fast loop** (the "react within fractions of a second after each impact" problem). Separate sub-system.
- **Fleet/formation interactions** (hundreds of PuffSats). PuffSats do not gravitate each other meaningfully; each orbit is independent, so the single-PuffSat study generalizes by running independent ensembles. Collision-avoidance gaps are a later concern.
- **Neural augmentation** (KalmanNet-style estimator, transformer warm-start, RL anomaly recovery). Deferred to a second pass; pass one is classical UKF + MPC to establish a clean baseline.

---

## 3. Reference orbit and the apogee decision

Constants: Earth gravitational parameter `mu = 398,600 km^3/s^2`; Earth radius `~6,378 km`; Earth Hill radius `r_Hill ≈ 1.496e6 km`.

The orbit periapsis is **50 km** — intentionally below the Kármán line (~100 km) so the PuffSat burns up on the same pass after impact, handling its own debris disposal. The **interception occurs at 200 km altitude during descent**, before the PuffSat reaches its 50 km periapsis; this altitude is the control target. For orbital mechanics, `r_p ≈ 6,428 km` (the ~150 km difference from the 200 km interception altitude changes the sensitivity formula in Section 8 by <3%, negligible). A 300 km interception fallback is noted in the paper if 200 km proves unmanageable.

Apogee is a **free design knob**, and choosing it well is one of the largest controllability levers available.

### Candidate apogees

| Quantity | 0.9 × Hill (1.346e6 km) | Lunar distance (384,400 km) | Recommendation: ~150,000 km |
|---|---|---|---|
| Semi-major axis a | 676,000 km | 195,500 km | ~78,300 km |
| Eccentricity e | 0.990 | 0.966 | ~0.916 |
| Orbital period | ~64 days | ~10 days | ~3.9 days |
| Apo→peri coast | ~32 days | ~5 days | ~2 days |
| Perigee speed v_p | 11.0 km/s | 10.9 km/s | ~10.8 km/s |
| Apogee speed v_a | ~54 m/s | ~187 m/s | ~330 m/s |
| Solar tide / Earth-g at apogee | **49%** | **1.1%** | **~0.1%** |
| Sensitivity dr_p/dv_a | **~250 km per m/s** | **~72 km per m/s** | **~30 km per m/s** |
| SRP-driven perigee dispersion | baseline | ~10–20× smaller | smaller still |

### Key finding: pull apogee in

The hypervelocity impact rationale ("max speed at the bottom") is **nearly independent of apogee**: for `e → 1`, perigee speed approaches escape speed at perigee (~11 km/s) regardless of how far out apogee sits. So pulling apogee in from 0.9 Hill costs almost nothing in impact speed (11.0 → ~10.8 km/s) while buying large reductions in:

- **Solar-tidal complexity.** At 0.9 Hill the Sun's tidal acceleration is ~49% of Earth's gravity at apogee; the orbit sits near the stability boundary and is dynamically sensitive (near-chaotic). At ~150,000 km it falls to ~0.1%, a benign highly-elliptical orbit.
- **SRP-driven dispersion.** Scales roughly as `a^2.5` (coast time ∝ `a^1.5`, leverage ∝ `a`). Pulling apogee in by 3–4× cuts SRP dispersion by ~an order of magnitude.
- **Error amplification.** `dr_p/dv_a` scales ~linearly with apogee radius, so a closer apogee makes the trajectory less sensitive to every velocity error at apogee.

Cost of pulling apogee in: separation/adjustment maneuvers at apogee cost ~3.5× more delta-v for the same perigee effect (because apogee speed is higher, lowering the aim leverage), but this is 3.5× of a small (m/s-scale) number. Also a shorter coast tightens the deployment timeline. An apogee inside lunar distance (e.g. ~100,000–200,000 km) additionally **avoids intermittent lunar close encounters**, which an apogee right at the Moon's orbit cannot reliably phase around (the Moon laps a ~10-day orbit every 27 days).

**Recommendation:** treat apogee as a tunable parameter, default it to ~150,000 km for the controllability study, and let the Monte Carlo show how perigee dispersion trades against apogee for the actual mission. Do **not** anchor on 0.9 Hill; its only advantage (lower apogee speed → more aim leverage per m/s) is outweighed by its disturbance amplification.

### Two precision scales (note for later)
- **Orbit-level** (this study): perigee altitude and interception timing, kilometre-to-metre.
- **Terminal-level** (follow-on): ~5 cm std-dev lateral centering at the pusher plate just before gasification.

---

## 4. System architecture: closed-loop, three decoupled pieces

```
        ┌──────────────────────────────────────────────────────────┐
        │  TRUTH MODEL (high fidelity)                              │
        │  full geopotential, SPICE Sun/Moon, NRLMSISE drag with    │
        │  stochastic F10.7/Ap, SRP with eclipse + attitude         │
        └──────────────┬───────────────────────────────────────────┘
                       │ true state (sampled at sensor cadence)
                       ▼
        ┌──────────────────────────┐     residual / unmodeled
        │  SENSOR + UKF estimator   │     dynamics = truth − filter
        │  estimates state +        │
        │  Cd·(A/m), Cr·(A/m)       │
        └──────────────┬───────────┘
                       │ estimated state
                       ▼
        ┌──────────────────────────┐
        │  CONTROLLER (MPC)         │  midcourse: discrete burns
        │  classical first;         │  terminal: ~0.2–1 Hz replan
        │  neural warm-start later  │
        └──────────────┬───────────┘
                       │ impulsive / finite maneuver command
                       ▼
        back into TRUTH MODEL as an event: stop, apply Δv, restart
```

The load-bearing consequence: **the truth propagator does not run as one continuous month-long shot.** It runs in short smooth arcs, stops at each maneuver/pulse and at each regime boundary, applies the velocity change instantaneously (or as a finite burn), and restarts. This event-restart structure is what makes the simulation cheap and accurate, and it dissolves the original "100 Hz integration for a month" framing.

---

## 5. The four clocks (do not conflate them)

| Clock | Rate | Notes |
|---|---|---|
| Integrator internal step | adaptive, set by accuracy | Non-uniform: minutes near apogee, sub-second through perigee. **Not** a user-chosen rate. |
| Sensor / UKF update | **altitude-scheduled: ~0.03 Hz coast → 100 Hz terminal** | The "100 Hz" from the original framing is the *terminal* rate. Best read as an accelerometer-class rate (GOCE-style drag measurement). Position fixes (GNSS/ranging) are more like 1–10 Hz. Schedule below. |
| MPC replan (outer) | per-maneuver in coast; **~1–10 Hz** terminal | Sets the reference trajectory + feedforward, not the per-sample command. A constrained solve at 100 Hz is unnecessary (the reference evolves slowly) even if a warm-started QP could do it. |
| Inner tracking loop | up to **100 Hz** terminal | Cheap linear law (PID/LQR) rejecting drag by modulating thrust *magnitude*, tracking the MPC reference between replans. Drag alone likely closes at ~10–20 Hz; 100 Hz is headroom for the cm-aim. |
| Output / logging | mode-dependent | Dense (100 Hz full state) for single runs; per-run summary for the big Monte Carlo. |

Sample/observe at any cadence cheaply by **dense-output interpolation** off the adaptive integrator (both high-order adaptive RK and IAS15 ship interpolants). A 100 Hz UKF feed costs a polynomial evaluation, not a 100 Hz force-evaluation cadence.

**The sensor clock is altitude-scheduled, not flat (settled 2026-06-03).** The required rate is set by how fast disturbances build state error and how fast the loop must respond, both negligible in coast and large in the terminal phase. Atmospheric density drives both and is roughly exponential in altitude (scale height ~30–70 km), so the principled rule is **rate ∝ local drag acceleration (or density), clamped between a coast floor and a 100 Hz terminal ceiling** — which auto-adjusts when a solar storm lifts the drag-onset altitude. A single sharp step at 600 km is the wrong shape: drag turns on around 800 km (§6.3) and the 800→600 km band (crossed in ~2 minutes) is exactly where it ramps from negligible to dominant, so a step there under-samples the onset and the filter enters the terminal burn still catching up. For the first build, approximate the exponential with a few hard-coded tiers, starting the ramp *above* 600 km:

| Regime | Altitude | Sensor rate |
|---|---|---|
| Deep coast | apogee → ~2,000 km | ~0.03–0.1 Hz (every 10–30 s; coarser near apogee fine) |
| Approach | ~2,000 → 800 km | ~1 Hz |
| Drag onset | ~800 → 600 km | ~10 Hz |
| Terminal | below 600 km | 100 Hz |
| Final cm-aim | last approach | possibly >100 Hz (deferred terminal sub-problem) |

Density-proportional scheduling is the refinement once the tiers prove the point. Note this is the *sensor* clock; the control clock above 600 km is even sparser (two discrete midcourse burns, §9, §16.6), and the MPC replan clock is separate again.

**Replan rate ≠ actuator command rate ≠ gimbal slew rate (settled 2026-06-03, do not weld these either).** The required *thrust direction* sweeps only at about the orbital rate near perigee (`v_p/r_p ≈ 0.1 deg/s`), so the direction loop is slow (~1 Hz, ~1 deg/s) and 100 Hz gimballing is both unnecessary and mechanically out of reach (gimbal/TVC bandwidth is single-digit Hz). What benefits from ~100 Hz is the thrust *magnitude* for drag rejection, modulated by valve PWM / pulse-frequency modulation of a Newton-class cold-gas or small monoprop thruster (bang-bang averaging to "continuous"), not by gimbal motion. MEMS microthrusters (µN–mN) are ~10^3 too weak for the ~0.4–few N terminal drag force; they are attitude / fine-pointing devices, not the translational drag-rejection actuator. So the terminal actuator is a *fast magnitude loop* (PWM, minimum-impulse-bit limited) plus a *slow direction loop* (gimbal or whole-body attitude); this is a Rung B input.

---

## 6. Propagation strategy

### 6.1 Formulation vs integrator (clearing the original confusion)
- **Cowell's method** is a *formulation*: integrate the second-order Cartesian ODE `r'' = -mu·r/|r|^3 + Σ a_perturb` directly. It is not an alternative to an integrator.
- **IAS15** is an *integrator* (15th-order adaptive Gauss–Radau). It can integrate Cowell's formulation. "IAS15 vs Cowell" is a category error.
- **Symplectic integrators (e.g. WHFast) do not apply.** Not mainly because drag is non-conservative, but because (a) the long-term energy-conservation benefit is irrelevant on a short dissipative + controlled arc, (b) impulsive thrust destroys the Hamiltonian structure anyway, and (c) fixed-step symplectic schemes handle a high-eccentricity perigee badly.

### 6.2 Regime-switched propagation
The orbit's perigee/apogee speed ratio is ~200:1 (≈ 42,000:1 in acceleration). One method for the whole arc is wrong in both directions: a single fixed step is absurd, and a single pure-adaptive scheme chokes on terminal discontinuities. Switch methods by regime, handing off osculating state at altitude boundaries.

| Regime | Span | Dynamics | Method |
|---|---|---|---|
| **Coast** | apogee → ~1,000–2,000 km | smooth: third-body (SPICE), SRP, J2; no drag | adaptive high-order Cowell (DOPRI8 / Dormand-Prince, or IAS15 as a cross-check); tight tolerance because injection/SRP sensitivity makes integration error masquerade as physical dispersion |
| **Terminal** | below hand-off → perigee | stiff + discontinuous: drag on, eclipse transitions, J2+ strong, control pulses | fixed-step Cowell (DOPRI8 or RK4; it is only ~minutes) **segmented at control instants with event-restart at each pulse**; machine precision pointless (drag uncertainty dominates), deterministic cadence aligns with the discrete UKF/MPC |

Regularized formulations (Sundman time-stretch, KS, EDromo, Dromo) are the elegant classical answer to `e ≈ 0.99` and would cut coast step count, but they complicate real-time events ("fire at wall-clock T", "sample at 100 Hz") because you step in a fictitious time. **Not recommended for the first build**; revisit only if coast propagation becomes a measured bottleneck (it will not — see Section 11).

### 6.3 Regime hand-off altitudes
- **~800 km:** turn the drag model on. Conservative guard band; drag is actually negligible here (~1e-8 m/s^2 at ~10 km/s through ~1e-14 kg/m^3). Cheap insurance against missing the onset.
- **~600 km:** start the terminal continuous burn (drag rejection + final aim).
- **200 km:** interception altitude — PuffSat impacts the pusher plate here, during descent.
- **~120 km:** onset of significant aerobraking on the 50 km periapsis descent. The PuffSat passes through this on every pass by design; it is not a failure criterion. It marks where thermal and structural loads accumulate after the 200 km interception.
- **50 km:** orbit periapsis. Below the Kármán line; PuffSat burns up here after impact — intentional debris disposal.

Sensor sampling ramps with altitude across these boundaries (≈30 s in coast → 100 Hz in terminal), not as a single step; the ramp starts above 600 km so the filter enters the terminal burn already at full rate. See the schedule in §5.

---

## 7. Force models: feedforward vs feedback

Do not think "drag vs everything else." Think three tiers of guidance, mapped to how each force is handled.

| Tier | Forces | Handling |
|---|---|---|
| **Feedforward (known model)** | monopole, J2 + higher geopotential, third-body Sun/Moon (SPICE), solid/ocean tides, relativity | open-loop nominal, recomputed from current estimate |
| **Slow feedback (estimated)** | SRP + slowly-accumulating coast divergence | estimate `Cr·(A/m)` as a UKF state; correct with discrete low-bandwidth midcourse burns |
| **Fast feedback (estimated)** | atmospheric drag near perigee | estimate `Cd·(A/m)` as a UKF state; reject in real time with the terminal burn |

Notes on individual forces:
- **SRP is not "easy."** It depends on A/m, attitude, and optical properties of a 25 kg sat that may spin, deploy bladders, or have its albedo deliberately altered by lasers. Over a long coast it deposits ~0.5 m/s of delta-v at the highest-leverage point (apogee). It is the **second hardest force after drag**, and at apogee possibly the first. Eclipse entry/exit is a discontinuity (event). First-build idealization: a fixed effective `Cr·(A/m)` on a sphere-equivalent (cannonball), not a cylinder; saying "cylinder" implicitly signs up for an attitude model (orientation relative to Sun and velocity), a separate subsystem deferred to a later rung.
- **Relativity** (Schwarzschild / post-Newtonian) has a tiny *acceleration* (~1e-9 of the monopole at perigee), but it accumulates coherently as a prograde apsidal precession: for this high-e orbit Δϖ = 6π·GM/(c²·a(1−e²)) ≈ 6.8e-9 rad/orbit, which displaces perigee by **~4 cm/orbit** (≈1 m at the apogee endpoint, measured). So it is **not** sub-millimetre, and "dropping it costs nothing at cm accuracy" is wrong for this orbit — at cm accuracy dropping it costs ~the whole budget. It is **negligible at the orbit-level (km-to-m) scale but right at the deferred ~5 cm terminal-centering budget**, so it is carried in the truth model (`full_force`, Orekit's `Relativity`), **truth-only** (below the estimator noise floor, not an onboard/feedback force). Its closed-form apsidal advance is a clean analytic signature. There is no separate "special relativity" force; the velocity-dependent terms live in the same PN expansion.
- **Drag** carries the real epistemic uncertainty (10–30% density error, worse in storms), driven by F10.7/Ap (forecast-uncertain), local time, season, latitude, and Cd·A/m. It bites below ~300–400 km, not 800 km.
- **"Pre-computable" ≠ "predictable trajectory."** The third-body and SRP *forces* are smooth and need no high-bandwidth feedback, but in a sensitive orbit the *trajectory* diverges from nominal, so you still need closed-loop correction of accumulated drift, just at low bandwidth.
- **Deferred or dropped (first build).** *Yarkovsky thermal recoil*: dropped — wrong regime (an asteroid / decade-scale effect; ~cm of displacement over an hours-long descent, orders below drag uncertainty). *Earth albedo + Earth IR radiation pressure*: deferred — real (~10–40% of direct SRP for a low skimmer) but below the drag/SRP uncertainty floor for a first controllability cut. *Solar-constant (TSI) variation*: dropped — the intrinsic ~0.1% bolometric fluctuation (11-yr cycle / 27-day rotation; flares are ~0.01% bolometric) is ~200× below the 10–30% `Cr·(A/m)` uncertainty and is *degenerate* with that estimated coefficient in the SRP force law (`F ∝ P₀·Cr·(A/m)·geometry`), so the UKF absorbs any slow drift. The large deterministic intensity swing (±3.4% from Earth–Sun distance) is already carried by the `(d₀/r)²` ephemeris term; the solar variability that *does* bite enters through `F10.7`/`Ap` on **drag** (§10.2, §16.7), not SRP. *Cylinder attitude-dependent area*: deferred — folded into the lumped cannonball `Cr·(A/m)` / `Cd·(A/m)` for now.

---

## 8. Controllability analysis (the governing numbers)

For a tangential velocity error at apogee, perigee radius moves by:

```
dr_p / dv_a  =  v_a · (r_a + r_p)^2 / mu
```

| Apogee | dr_p/dv_a | Implication |
|---|---|---|
| 0.9 Hill | ~250 km per m/s | 0.5 m/s SRP drift → perigee 200 → ~75 km → burn-up |
| Lunar | ~72 km per m/s | 3.5× less sensitive |
| ~150,000 km | ~30 km per m/s | ~8× less sensitive than 0.9 Hill |

**Required apogee velocity control accuracy** to arrive within a few km of nominal perigee:
- 0.9 Hill: ~1–2 cm/s (few km ÷ 250 km per m/s).
- Closer apogee relaxes this proportionally.

This accuracy requirement is the actual hard problem of the project. The simulation exists to show whether it is achievable under realistic navigation and disturbance.

---

## 9. Guidance strategy: discrete midcourse + terminal burn

The paper floats a "let drift accumulate, fix it all in a gimballed 600 → 200 km burn" strategy. **It does not work for accumulated drift, only for drag.**

- **Terminal delta-v budget is small.** 400 g cold gas, Isp ≈ 200 s, 25 kg → total delta-v ≈ 32 m/s. The 600 → 200 km descent takes ~5 minutes near perigee. Correcting Δx of lateral miss in ~300 s costs ≈ `2·Δx/t` of delta-v: ~5 km eats the entire budget. So the terminal burn can null only a **few km** of residual.
- **Letting the coast drift uncorrected** puts you tens-to-hundreds of km off (250 km per m/s). Fixing 100 km in 5 minutes would need ~670 m/s and ~2 g — infeasible.

Therefore:
1. **Discrete midcourse corrections during the coast** (event-triggered, e.g. one near apogee, one mid-descent) null the bulk of the dispersion *cheaply*, where 1 m/s buys hundreds of km of perigee. Standard statistical-midcourse + terminal guidance, as in interplanetary practice.
2. **Continuous gimballed terminal burn 600 → 200 km** handles drag rejection and the final few km of aim. It must be *handed* a near-nominal state by the midcourse corrections.

No continuous burn is needed above ~600 km (drag negligible); only a small number of discrete corrections.

---

## 10. Monte Carlo design

### 10.1 Two deliberately mismatched models
- **Truth model:** full fidelity (full geopotential, SPICE Sun/Moon, NRLMSISE drag with a stochastic F10.7/Ap draw, SRP with attitude + eclipse). Generates "reality."
- **Onboard model (UKF/MPC):** reduced, with `Cd·(A/m)` and `Cr·(A/m)` as estimated states. Flies the vehicle.
- The residual between them is the disturbance the controller must reject. One run with matched models proves almost nothing.

### 10.2 Sampled uncertainties (per run, seeded)
- Injection error at apogee (velocity, the dominant lever via `dr_p/dv_a`).
- **Lumped coefficients `Cd·(A/m)` and `Cr·(A/m)`**, each drawn once per trajectory from a multiplicative log-normal (~10–30% 1-sigma) and held fixed for the run. Rationale: a *constant* coefficient error integrates coherently over the multi-day coast and is what drives perigee dispersion; high-frequency in-flight flutter (a slow Gauss-Markov wiggle on top) largely averages out and is deferred to a later robustness pass. The drag-coefficient sigma moves the burn-up tail more than any other single input.
- Atmospheric density driver (F10.7, Ap) and any storm realization.
- Lunar/solar ephemeris error (small with SPICE, include for completeness).
- Sensor noise (the 100 Hz channel) and any dropouts/outliers.

### 10.3 Output metrics (per run)
- Interception success: PuffSat reaches 200 km altitude at the correct position and time (pass/fail vs miss-distance threshold). Periapsis altitude is not a failure criterion — 50 km is the intended orbit periapsis.
- Interception miss (along-track timing + lateral position).
- Propellant expended (vs the < 2% claim).
- Seed + draws (to enable replay).

### 10.4 Ensemble sizes
- **N = 1** (single trajectory): debugging workhorse.
- **N = 50:** smoke test for the harness and edge cases. **Not** a result — cannot estimate a ~1% burn-up tail from 50 samples.
- **N = 10^3 – 10^4:** the real controllability result, sized to resolve the low-perigee tail.

---

## 11. Tooling

### 11.1 Decision
- **Closed-loop stochastic study (this sim): Orekit** (Java core, driven from Python). It is a *library you drive*, which fits the segmented closed-loop where your own UKF/MPC injects maneuvers at each control instant. Its event-detection framework fits the regime hand-offs (altitude triggers), eclipse boundaries, and per-pulse stops directly. Most mature, best-documented, most supportable.
- **Trajectory design + maneuver optimization + independent validation oracle: GMAT** (NASA GSFC). It is an *application that owns the mission sequence*, which is the wrong shape for closed-loop-with-your-own-controller but the right shape for designing the deployment orbit, the apogee separation strategy, and the deorbit targeting (built-in differential corrector / optimizer), and for cross-checking the Orekit truth model on the nominal trajectory.
- **Two tools, two jobs:** design and verify the reference in GMAT; run the closed-loop Monte Carlo in Orekit.

### 11.2 How Orekit + Python interoperate
**In-process JVM bridge, not a server.** A JVM is started inside the Python process (`orekit_jpype`, JPype-based, recommended; the older JCC-based `orekit` is the native-shim path and is not recommended). Orekit's Java classes appear as Python objects; calls cross the Python↔JVM boundary via JNI in the same process and memory (single-digit microseconds, not a network hop). Your PyTorch MPC lives in the same Python process and is called between arcs. Monte Carlo parallelism is process-level (N workers, each its own Python + JVM, JVM cold-start amortized over many trajectories per worker).

**Why in-process rather than a gRPC/REST server.** A closed loop is chatty: every control cycle pulls state and pushes a maneuver. In-process JNI is single-digit microseconds per call; a gRPC loopback round-trip is tens to low-hundreds of microseconds (serialize + socket + deserialize), 10–100× more, compounding across cycles × trajectories. The rich Orekit objects (`SpacecraftState`, covariance, frames) would each need a protobuf schema or be flattened to bare numbers and lose their methods; and Orekit event detectors could not call back into a Python controller mid-propagation without bidirectional streaming. A server also reintroduces lifecycle, port, and client/server-version management, and forces either one JVM per worker (no savings) or a shared serialization bottleneck. The "two runtimes in one address space" fear is calibrated for embedding a native C library that can segfault and corrupt memory; Orekit is managed-code numerics that does not. That scary reputation belongs to the old JCC shim, not modern JPype.

**Boundary-crossing strategy (two modes), so 100 Hz never becomes 100 Hz of crossings:**
- *Coast:* Orekit densely propagates a segment and hands Python one dense-output interpolator per crossing; the UKF samples its 100 Hz measurements locally off that batch, dropping crossings to roughly the maneuver rate.
- *Terminal:* register the controller as an in-process Orekit step/event handler (a Python callback via JPype), so the integrator drives it at the right cadence with internal state intact across each call.

Friction to expect: Java idioms in Python (`AbsoluteDate`, `TimeScalesFactory`, `Vector3D`, JPype overload/coercion quirks); Java exceptions surface with murky cross-boundary stack traces.

### 11.3 Alternatives considered
- **tudatpy** (TU Delft, C++ core + Python): solid, native SPICE, strong estimation, pure-Python single process. Best pick if frictionless Python and a deep neural integration outweigh operational pedigree and debugging support.
- **Basilisk** (CU Boulder, C++ + SWIG): purpose-built GNC co-sim with actuator/sensor/flight-software modules and a Monte Carlo harness. Best if the project's centre of gravity becomes the GNC co-sim itself; steeper learning curve, younger.
- **REBOUND / IAS15** (Hanno Rein): excellent N-body integrator, great for an *independent coast cross-check*, but lacks native SPICE, atmosphere models, maneuver/estimation/control plumbing. Its N-body strength is wasted here (PuffSats do not gravitate each other). **Security note:** the project's marketing site `rebound.hanno-rein.de` showed signs of the polyfill.io supply-chain hack as of 2026-06-02. Pull REBOUND only from the `rebound` PyPI package or the `hannorein/rebound` GitHub repo; do not interact with popups on that site.
- **Nyx** (Rust, `nyx-space`): modern, capable, Python bindings, but largely one maintainer and thinner docs. Reach for it only if Rust-native is a hard requirement.
- **Pure-Python astrodynamics (poliastro lineage: boinor, hapsira)** (NumPy/Numba, SciPy `solve_ivp` Cowell): clean API, strong for analytic two-body, Lambert/porkchop, and plotting. **Not the truth engine.** Ships only toy perturbations (exponential-atmosphere drag, J2/J3 examples), no high-fidelity density (NRLMSISE/JB2008/DTM), no conical-shadow SRP, thin frames/time/EOP, no event-restart, no estimation framework, no operational validation. Driving Orekit from Python (Section 11.2) already gives Python ergonomics with a validated engine, so "Python-only" buys nothing here. Fine as a lightweight design-phase scratchpad or an independent analytic cross-check (same niche as GMAT in this project), not as the closed-loop truth model. The real near-peer to Orekit is tudatpy, not boinor.

### 11.4 On dropping into C++/Rust for performance
Not needed. The force evaluation and integration are already compiled (C++/Java). The real hotspot is the MPC solve, which stays in Python for the neural second pass and is accelerated with a **compiled QP backend (OSQP or qpOASES, both with Python bindings)**, not a rewrite. Touch C++ only to add a force model the library does not ship; for this stack (SPICE third-body, NRLMSISE drag, SRP + shadow, geopotential) they ship everything.

### 11.5 MPC language
**Keep the MPC in Python.** Right for the neural warm-start (KalmanNet-style estimator, transformer warm-start) and for the RL anomaly-recovery path (RL treats the sim as a black-box environment, which a Python-driven loop gives for free). The truth-tool choice does not constrain the MPC language: the segmented architecture decouples them, requiring only stop-at-event and inject-maneuver-from-outside, which all candidate tools provide.

---

## 12. Compute budget

Compute is **not** the constraint, and not in the way originally framed. Integrator order (IAS15 vs RK8) is rounding error. The bill is:

```
N_ensemble × N_perigee_passes × MPC_solve_cost
```

A single warm-started QP solve dwarfs thousands of force evaluations. At N = 10^4 runs × a few hundred MPC solves each × ~ms per QP ≈ single-digit core-hours, embarrassingly parallel across runs. A laptop-overnight or small-cluster-minutes job. The compute knobs that matter are **ensemble size** and **MPC replan rate**, both independent of the integrator.

### 12.1 The onboard replan rate is bound by usefulness, not compute (settled 2026-06-03)
Distinct from the offline Monte Carlo cost above: the *onboard* MPC replan rate is **not compute-limited** at the rates that matter. A small convex MPC (≈6 states + a few augmented, modest horizon) solves warm-started in tens of µs to low ms (OSQP/qpOASES), supporting 100 Hz–kHz; even a real-time-iteration NMPC stays sub-10 ms. Against a 1–10 Hz terminal target that is 2–3 orders of headroom. The binding limits are instead (i) *information* (no point replanning faster than new state arrives, §5) and (ii) *usefulness* (the optimal reference evolves on minute-to-hour coast / second terminal timescales, so a faster resolve returns a near-identical plan). Fast rejection is delegated to the 100 Hz inner loop so the MPC need not be fast.

**What the neural accelerator actually buys (paper's neural section).** Since compute does not bind the rate, an NPU does *not* usefully raise the replan rate. Its leverage: (a) **lower energy per solve** on a power-limited 25 kg single-use sat; (b) a **learned fast inner policy** (RL/imitation or explicit-MPC approximation) running at 100 Hz–kHz as approximate-MPC, with the true constrained solve as a slower supervisor — the accelerator *is* the fast loop; (c) a **learned warm-start** making worst-case solve time predictable enough to *guarantee* the real-time deadline; (d) the **RL anomaly-recovery** path. Compute binds only if the MPC is pushed to a full-fidelity stochastic/scenario NMPC with collision constraints — and that is exactly where the learned warm-start earns its keep. Surplus compute is better spent on horizon, model fidelity, or robustness (scenario/tube MPC) than on raw rate.

---

## 13. Build ladder

Two stages, built in order. **Stage 1** stands up and validates the truth model: physics only, perfect knowledge, no control. **Stage 2** adds guidance, then estimation, then the full stochastic study. Within each stage the discipline is the same: add exactly one force or one erosion per step, each checked against a known analytic signature, so any failure has a single suspect.

### Stage 1 — truth model (physics, perfect state, no control)

- **Rung 0 — build and "hello world."** Prove the toolchain. `orekit_jpype` installs and imports; the `orekit-data` bundle loads (leap seconds, EOP, gravity field); a date and a frame round-trip; print a state vector. No dynamics. Projects often stall here on classpath / data-file issues, so isolating it de-risks everything after.
- **Rung 1 — unperturbed two-body orbit, perfect known state.** Monopole gravity only, propagate one period of the reference orbit. Sampling/observation is just a step-handler on this propagator (no separate "sampling" rung). Checks: recover the design numbers (perigee 6,578 km, apogee per Section 3, period, v_p, v_a); conserve energy and angular momentum to integrator tolerance (validates frames, time scales, units, integrator at once); nudge apogee velocity by exactly 1 m/s and confirm perigee moves by the predicted `dr_p/dv_a` (validates the controllability lever before anything is built on it).
- **Rung 2 — perturbations, one force at a time, cumulative.** Add forces in order of decreasing magnitude and increasing messiness, validating each signature before adding the next:
  - **2a — J2, extending to an 8×8 field for the truth model.** Validate pure J2 (degree 2) first: nodal regression and apsidal precession at analytically known rates; if these match, frames/time/integrator hold up under a perturbation. The `full_force` truth model then extends the geopotential to **degree/order 8** (default, tunable). Beyond J2 the harmonics fall off as (Rₑ/r)^ℓ and only bite near perigee, but over a pass they reach **~km scale at orbit level** (measured ~1.8 km one-period endpoint drift vs J2 for the 8×8 field; tesserals dominate, ~1.8 km vs ~0.18 km zonal-only), so they belong in the orbit-level truth model — not only the deferred terminal problem. The coast does not need them; raise the degree further for the regime-switched terminal phase.
  - **2b — third-body Sun + Moon (SPICE).** The "tidal forces" item: point-mass third-body gravity (not solid-Earth tides). Validate the known long-period signature.
  - **2c — SRP, cannonball area.** Fixed effective `Cr·(A/m)` on a sphere-equivalent (Section 7). Validate eclipse entry/exit timing (conical-shadow event) and the along-track signature.
  - **2d — atmospheric drag.** The messy non-conservative one, last. Validate energy and perigee-altitude decay rate. Bites below ~300–400 km.
  - **2e — relativity (Schwarzschild), truth-only.** Conservative and tiny in acceleration but accumulates as apsidal precession: ~4 cm/orbit on this high-e orbit — negligible at the orbit-level (km-to-m) scale, but the deferred ~5 cm terminal budget (Section 7). Validate against the closed-form apsidal advance Δϖ = 6π·GM/(c²·a(1−e²)). Folded into `full_force` (which also carries the 8×8 geopotential, 2a).
  - Introduce the **per-run sampled coefficients** here (Section 10.2): `Cd·(A/m)` and `Cr·(A/m)` each drawn once per trajectory from a multiplicative log-normal and held fixed; a constant coefficient error integrates coherently over the coast and drives the dispersion, while high-frequency flutter averages out and is deferred. Dropped/deferred forces per Section 7 (Yarkovsky, Earth albedo/IR, cylinder attitude).

- **Stage 1 capstone — open-loop dispersion (settled 2026-06-07).** Before any control, run the `full_force` truth model with the Section 10.2 per-run sampling (injection error, log-normal `Cd·(A/m)`/`Cr·(A/m)`, F10.7/Ap) but **control off**, N ≈ 100–1000, summary stats only. Output: the **interception-state dispersion** at the 200 km descent crossing (covariance of position, velocity, and time-of-arrival) about the nominal aim point (Section 16.3). Sizes the midcourse Δv authority and surfaces any bias *before* the corrector is designed, validates Section 8's analytic dispersion against full nonlinear propagation, and is the baseline the closed loop must beat. Truth-only — no estimator/onboard model yet.
  - **Harness design (settled 2026-06-08):** injection error is sampled as an **RTN Cartesian Δv on the apogee deployment state** (Radial / Transverse / Normal; at apogee the radial velocity is zero so **T = the tangential `dr_p/dv_a` lever**, R = timing, N = cross-track). Coefficients/drivers are **median-nominal multiplicative log-normals** (`Cd·(A/m)`, `Cr·(A/m)`, F10.7, Ap; median = nominal, so the slight mean-bias is real and reported, not hidden). Per-run seeding uses NumPy `SeedSequence.spawn`, each run replayable standalone (§14.2). A pure `DispersionSpec` / `RunInputs` / `sample_run_inputs` core stays JVM-free and unit-testable; only the propagate-and-record step touches Orekit. **Primary metric:** interception miss decomposed in the **nominal-crossing RTN frame** (along-track ≈ timing, cross-track, radial) + time-of-arrival. **Perigee** is a *dual diagnostic only* — the §8 lever, and a **debris-disposal-safety margin** (low = good; a *missed* PuffSat must deorbit, paper §9) — never an interception pass/fail. Built as the one parameterized path (§14.1) with an optional control hook (absent here) so Rung D reuses it.

### Stage 2 — control (guidance, then estimation, then the stochastic study)

Rung A is decomposed in full below, because it is the controllability experiment and the first place real design choices bite. **Rungs B–D are headlined only and will be decomposed in detail when we reach them**, once Stage 1 and Rung A have taught us the actuator authority required, the dispersion magnitude, and the filter observability. Decomposing B–D now would be guessing ahead of those results.

Ordering rationale (control before estimation): a noisy state readout is meaningless until a controller consumes it. Establish "with perfect knowledge the controller hits the target" (Rung A) before asking "does it still hit with noisy knowledge" (Rung C). This isolates failure causes, adding under-actuation, blindness, and MPC tuning one at a time against a known-good baseline rather than all at once.

- **Rung A — controllability core: perfect state, impulsive Δv, deterministic targeter.** Decomposed below.
- **Rung B — actuator realism.** Replace impulsive Δv with finite gimballed burns: thrust magnitude, Isp, mass depletion, dead-time / latency. Still perfect state. Question: can the real ~400 g actuator deliver the needed correction in the minutes available, especially the terminal descent? *Detailed decomposition deferred.* Early inputs (§5): model a fast magnitude loop (PWM / pulse-frequency, minimum-impulse-bit limited) plus a slow direction loop (gimbal/attitude rate-limited ~1 deg/s); cold-gas micro-pulses above 800 km, sustained pulsed thrust below; MEMS reserved for attitude, not translation. **100 Hz modulation feasibility (settled 2026-06-03):** clean for cold gas (no combustion transient; fast solenoid/piezo valves ms-class, or a *proportional* flow valve for smooth thrust with no MIB quantization), but marginal for chemical mono/biprop (ignition + chamber-pressure + catalyst-bed dynamics cap clean pulse-mode at ~tens of Hz). This is not a problem: drag rejection closes at ~10–20 Hz (chemical's comfortable range), and only the deferred cm-aim wants true 100 Hz, which is a cold-gas/electric fine stage anyway. Favor a **coarse/fine split**: chemical (or higher-thrust) carries the steady gross anti-drag burn (Isp matters, most Δv), cold-gas micro-pulses do the fast 100 Hz trim (small Δv, low-Isp penalty negligible). Watch-item: Appendix A assumes Isp ≈ 200 s, optimistic for true cold gas (N₂ ~70 s, He ~165 s); pin per-device Isp when sizing this rung, since the cold-gas portion may run below 200 s.

  **Settled 2026-06-08 — first-cut Rung B actuator model (ADR 0004).** A *single omnidirectional proportional cold-gas thruster*, not the coarse/fine chemical+cold-gas split above (that remains a deferred higher-fidelity option):
  - **Paradigm:** proportional throttle is the modeling *abstraction*; the physical realization is a bang-bang cluster averaging to "continuous" (§5). No MIB quantization in the model — the cluster's residual minimum-impulse-bit is a deferred refinement, and lands exactly where the cm-aim 100 Hz trim wants it.
  - **Thrust:** 400 mN max (GOCE-derived, Appendix A); ~5 mN floor (proportional ~80:1 turndown; the cluster duty-cycles below it); continuous between, ms-class response (instantaneous vs. the ~5-min terminal burn).
  - **Direction:** omnidirectional, direction-slew rate-limited to ~1 deg/s (§5). The required sweep is only ~0.1 deg/s near perigee, so the limit carries ~10× margin and is expected slack — a result to confirm, not assume. Thruster geometry / cosine losses / gimbal-vs-whole-body realization deferred until the slew limit is shown to bind.
  - **Isp:** a *reported sweep*, not a single value — Isp is a post-processing lever on the propellant ledger (`propellant ≈ Δv/(Isp·g₀)`), not on controllability, so the A-rung Δv is computed once and transformed. Anchors: 50 s (conservative — CO₂/argon/refrigerant-class storable, dense + self-pressurizing), 70 s (N₂ realistic), 200 s (Appendix A optimistic). At 50 s the §9 ~32 m/s budget costs ~6.4% of 25 kg (vs ~1.6% at 200 s), so the curve shows the paper's <2% propellant claim failing ~3× at the conservative anchor — a deliberate finding, with higher mass fraction the accepted trade.
  - **Mass depletion:** modeled (Tsiolkovsky; Orekit maneuver models deplete mass via Isp natively). ≤~6.4% acceleration rise across a burn at worst; makes total Δv → propellant mass a first-class output feeding the Isp sweep.
  - **Dead-time:** deferred to Rung C (§16.8). Valve/command latency is ms-class — negligible vs. minute-to-day maneuver timescales — and the latency that bites (off-board coordinator-node sense→compute→uplink→actuate) has no source under Rung B's perfect state; carried as a ~0-default parameter so the seam exists.
- **Rung C — estimation in the loop.** Feed the controller a noisy state instead of truth: a Gaussian placeholder first (covariance motivated by what coordinator-node ranging can plausibly deliver, not an arbitrary blob), then the UKF estimating `Cd·(A/m)` and `Cr·(A/m)`. Question: does the loop still close with realistic knowledge? *Detailed decomposition deferred.*
- **Rung D — MPC + Monte Carlo.** Swap the deterministic targeter for MPC (Python; neural second pass later); run N = 50 as a smoke test, then N = 10^3 – 10^4 for the result with a resolved perigee tail. *Detailed decomposition deferred.*

#### Rung A decomposition (controllability core)

Goal: with perfect knowledge and idealized thrust, is there enough control authority to hit the perigee target, and at what Δv? This answers the headline question in its purest form; B–D each erode the margin Rung A establishes.

- **A1 — single impulsive midcourse.** On the nominal forces from Stage 1, inject a known apogee error (e.g. 0.5 m/s). One impulsive midcourse correction at the apogee deployment node, solved by a differential corrector (vary the correction Δv to null the predicted **interception miss** — the 3-component RTN position miss at the 200 km EME2000 crossing, §16.3), restores the nominal crossing. Confirm the resulting perigee shift matches the `dr_p/dv_a` prediction (an **acceptance cross-check**, not the objective — perigee is a diagnostic, not the target) and the propellant ledger closes (Δv applied = Δv commanded). This is the old "Rung 1" check, now with perfect state and no filter. (Corrector objective settled by ADR 0003; supersedes the earlier "null predicted perigee error" wording.)
- **A2 — two-burn statistical midcourse + terminal aim.** Add the Section 9 structure: one correction near apogee, one mid-descent, plus a small terminal aim, all still impulsive and perfect-state. Schedule is **fixed and always-executed** (§16.6, settled): correction 1 near apogee, correction 2 at a fixed mid-descent altitude (~800–1,000 km); threshold/adaptive triggering is deferred to Rung D. Confirm the bulk of dispersion is nulled cheaply during coast (1 m/s buys hundreds of km) and only a few km of residual reaches the terminal phase. Regime hand-offs (800 / 600 km) and eclipse events fire correctly.
- **A3 — deterministic coefficient sweep.** Hold the targeter fixed; sweep `Cd·(A/m)` and `Cr·(A/m)` across their ranges deterministically (not yet random). Map how required Δv grows with coefficient error and where the dumb targeter runs out of authority. This is the controllability map, and it localizes the burn-up boundary before any noise or MPC is added.

Targeter choice (settled 2026-06-03): a deterministic differential corrector, not MPC, for Rung A. It yields clean controllability ground truth, a Δv floor and an authority boundary. MPC enters at Rung D and is checked against this baseline; were MPC to lead at Rung A, a miss could not be attributed between an uncontrollable orbit and a mistuned controller.

---

## 14. Engineering practices

1. **One parameterized code path, not a separate "single mode."** The single-trajectory run is the ensemble harness with `N=1`, noise flags off, stochastic inputs frozen to nominal. Knobs: `N`, `noise_on`, `mismatch_on`, `seed`, over a single scenario definition. Prevents the quick-mode and the real-MC from drifting apart.
2. **Seed per run, replayable.** Any ensemble run (including a failing run #7,432 at N=10^4) must replay in single-trajectory mode from its seed, with dense logging. Single-trajectory and Monte Carlo share one draw mechanism.
3. **Mode-dependent logging.** Single runs: full 100 Hz state + covariance + residuals + event log + commanded/applied Δv. Big MC: per-run summary only (perigee achieved, miss, propellant, pass/fail, seed) or storage drowns.
4. **Apogee as a tunable parameter** (Section 3), so the dispersion-vs-apogee trade is a sweep, not a recompile.

---

## 15. Acceptance / sanity checks (Stage 1 and Rung A)

- Conservative-only invariants conserved to integrator tolerance.
- Design orbit parameters recovered from initial conditions.
- `dr_p/dv_a` measured numerically matches the closed-form value for the chosen apogee.
- Known injected apogee error nulled by midcourse correction to the target perigee.
- Propellant ledger closes; maneuvers fire at the correct events; eclipse and altitude events fire at the correct boundaries.

---

## 16. Open questions and decisions to make

Each item below lists options with trade-offs and a lean. None is settled; these are the decisions to make before or during the first build.

### 16.1 Apogee altitude (most consequential)
Sets disturbance amplification, coast duration, and dynamical regime. See Section 3.
- **A. 0.9 × Hill (~1.35M km).** Maximum aim leverage per m/s (54 m/s apogee, so a tiny launcher delta-v redirects far). Costs: solar tide ~49% of local g (near-chaotic), ~32-day coast, largest SRP dispersion, `dr_p/dv_a ~250 km` per m/s (errors amplified), requires ~1-2 cm/s apogee velocity accuracy.
- **B. Lunar distance (~384,400 km).** Solar tide drops to ~1.1%. Costs: intermittent lunar close encounters that a ~10-day orbit cannot reliably phase around (Moon laps every 27 days).
- **C. ~150,000 km (recommended default).** Solar tide ~0.1%, no lunar encounters, ~2-day coast, `dr_p/dv_a ~30 km` per m/s, perigee impact still ~10.8 km/s. Costs ~more launcher delta-v for separation (still m/s-scale, small).
- **D. Sweep it.** Make apogee a parameter; let the Monte Carlo quantify perigee dispersion vs apogee, then pick.
- **Lean:** D for the study, defaulting to C as the working point. The decision driver is how much aim leverage the launcher actually needs versus how much disturbance amplification the control can reject. Do not anchor on A.

### 16.2 Interception altitude
The orbit periapsis is fixed at **50 km** (debris disposal by reentry; see Section 3). This decision sets the altitude at which the PuffSat impacts the pusher plate during descent.
- **A. 200 km (paper primary).** Maximum PuffSat-velocity benefit (~10.8 km/s at impact) and lowest target-rocket propellant fraction. Costs: strongest and most variable drag during the terminal phase, tightest control.
- **B. 300 km (paper fallback).** Atmospheric density roughly an order of magnitude lower, much more benign terminal phase. Costs: marginally higher target propellant fraction and slightly less velocity benefit.
- **Lean:** run both (interception altitude is a cheap parameter). Lead with 200 km since that is the paper's claim; the 300 km contrast shows the control difficulty trade-off.

### 16.3 Interception / target model (sets the miss metric)
- **A. Fixed perigee point + epoch.** PuffSat must arrive at a specified position at a specified time. Simplest; decouples PuffSat guidance from target motion. Miss = position + timing error there.
- **B. Prescribed target ascent trajectory.** PuffSat rendezvous with a moving target following a known rocket ascent profile (the paper notes the optimal interception point moves with the accelerating rocket, so successive PuffSats have slightly different elements). Miss = relative position/velocity at closest approach. Requires a representative ascent profile as input.
- **C. Full co-sim with target dynamics + its RCS.** Most realistic; pulls in the out-of-scope target-RCS fast loop.
- **Lean:** A through Rung A (prove you can hit a fixed point in space-time at all), then B for the real controllability result at Rung D (the target is moving; that is the actual mission). Defer C. Define the miss in the target's pusher-plate frame eventually, but start inertial.
- **Decision (settled 2026-06-07):** A for the open-loop capstone and through Rung A — the target is a **fixed interception point + epoch in inertial EME2000, defined as the nominal trajectory's own 200 km descent crossing.** The capstone reports the **interception-state dispersion** (covariance of crossing position, velocity, and time-of-arrival) about that nominal point: its mean is the open-loop aim point (and surfaces any bias), its spread sizes the midcourse Δv authority. The corrector (Rung A1) and later MPC null *crossing-state − target-state in EME2000*, one shared miss definition. Moving target (B) deferred to Rung D; pusher-plate frame deferred to the terminal sub-problem.

### 16.4 Sensor / navigation model
Drives UKF observability. Note a strong altitude dependence.
- **A. Accelerometer only (GOCE-style, 100 Hz).** Measures non-gravitational specific force (drag + SRP), so it directly anchors the `Cd·(A/m)` and `Cr·(A/m)` estimates. Does not observe absolute position; must be fused with a positioning source.
- **B. Position fix only.** GNSS (~1-10 Hz) is fully available near perigee but essentially unavailable at a 150,000 km apogee (above the constellation, weak spillover signals, poor geometry). Coordinator-node ranging/angles are available throughout (the paper's coordinator nodes track PuffSat positions) and are the realistic apogee anchor.
- **C. Fused suite (recommended realistic case).** Accelerometer at 100 Hz (drag observability) + coordinator-node relative measurements (anchoring throughout, including apogee) + GNSS when near Earth. Matches the paper's architecture.
- **Lean:** perfect state at Rung A to establish controllability, then a Gaussian placeholder and the fused suite C at Rung C to debug the filter. Flag that apogee navigation (coordinator-node relative measurements, GNSS unavailable) versus perigee navigation (GNSS available) is itself a sub-study: navigation quality is not uniform around the orbit.

### 16.5 Maneuver model (impulsive vs finite burn)
- **A. Impulsive delta-v throughout.** Instantaneous velocity jumps, event-restart at each. Simplest; fine for brief midcourse corrections.
- **B. Finite burn for the terminal phase.** The 600 to 200 km terminal burn lasts ~5 minutes, is continuous and gimballed, and is coupled to drag during the descent. Modeling it as a single impulse loses the drag-rejection-during-descent physics and misstates propellant use.
- **Lean / settled (2026-06-03):** hybrid. Impulsive throughout Rung A (midcourse plus terminal aim, to establish controllability ground truth); finite gimballed burn introduced at Rung B, when actuator authority and drag-rejection-during-descent become the focus.

### 16.6 Midcourse correction schedule
- **A. Fixed pre-planned schedule.** Corrections at predetermined trajectory points (e.g. one near apogee, one at a fixed descent altitude), always executed. Simple, deterministic, Δv budgetable a priori; matches the classical statistical-midcourse-correction paradigm.
- **B. Threshold / MPC triggered.** A correction fires when predicted perigee error exceeds a threshold. More adaptive, fewer wasted burns, but more logic, harder to bound Δv a priori, and its trigger depends on navigation quality.
- **Decision (settled 2026-06-03): A for all of Rung A.** Two reasons beyond simplicity. (1) Rung A measures controllability against a *fixed, transparent* control law; an adaptive scheduler would confound the A3 coefficient sweep — you could not tell whether the authority boundary reflects the orbit or the scheduler's cleverness. (2) Threshold-triggering's value (fire only when needed) is entangled with estimation quality: at Rung A's perfect state the trigger is trivially perfect, which overstates how well it works once navigation error enters at Rung C, so evaluating it here gives a misleadingly optimistic read. Move to B only at Rung D, where MPC subsumes it naturally and estimation error is in the loop to test it honestly. Concrete Rung A schedule: correction 1 near apogee (cheapest leverage), correction 2 at a fixed mid-descent altitude (~800–1,000 km, above the terminal hand-off), plus the terminal aim; all always-execute. These are fixed trajectory points (timing set by geometry), which is what "fixed schedule" means here — not fire-on-realized-error.

### 16.7 Atmospheric density model and uncertainty representation
The drag model is the hard, feedback-demanding force; how its *uncertainty* is represented is what makes the Monte Carlo meaningful.
- **Model choice:** NRLMSISE-00 (standard empirical thermosphere), JB2008 (Jacchia-Bowman, often more accurate for drag work), or DTM. Orekit supports several.
- **Uncertainty representation:** either (i) perturb the space-weather drivers F10.7 and Ap within their forecast uncertainty per run, or (ii) apply a multiplicative density bias plus process noise on top of a nominal model. Optionally inject storm realizations.
- **Truth vs filter:** the truth model should use one realization; the onboard filter uses a nominal model and estimates `Cd·(A/m)` to absorb the mismatch.
- **Lean:** JB2008 (or NRLMSISE-00 if simpler to wire first) for truth; represent uncertainty by perturbing F10.7/Ap plus a multiplicative density factor; nominal model in the filter with `Cd·(A/m)` estimated. This is the core of the controllability result.
- **Implemented (2026-06-07):** NRLMSISE-00 driven by **constant per-run F10.7/Ap** via a custom `NRLMSISE00InputParameters` provider (`forces/_space_weather.py`), set from the `AtmosphericDrag` spec — replacing the calendar-tied CSSI data so the Monte Carlo samples space weather per run by constructing the spec with different `f10p7` / `ap`. Higher F10.7/Ap inflate density monotonically (verified: ρ(200 km) 1.6→2.8→4.2 ×10⁻¹⁰ kg/m³ for F10.7 70/150/250). Still to add: multiplicative density factor and storm realizations.

### 16.8 Control-loop latency
The paper has off-board coordinator nodes that "account for communications and PuffSat actuator latencies." That is a dead-time in the loop (measure → coordinator compute → uplink → actuate), which erodes stability margin, most in the fast terminal phase.
- **A. Zero latency.** First-pass simplification.
- **B. Constant delay.** A fixed dead-time on the command path.
- **C. Delay + jitter.** Stochastic latency, sampled per step.
- **Lean (sharpened 2026-06-08, ADR 0004):** A (zero latency) through **Rung B** — Rung B is still perfect-state, the actuator valve/command dead-time is ms-class (negligible vs. minute-to-day maneuver timescales), and the dead-time that actually bites (measure → coordinator compute → uplink → actuate) has no source until estimation enters. B (constant delay) therefore enters at **Rung C** with the estimator/coordinator loop, carried as a ~0-default parameter until then; C (delay plus jitter) only if the terminal phase shows sensitivity to it. Worth surfacing because a fast terminal loop with unmodeled dead-time can go unstable in a way the paper's architecture explicitly tries to avoid.

### 16.9 First deliverable
- **A. Rung 0 then Rung 1** (build / hello-world, then the unperturbed orbit with its sanity and leverage checks). Smallest artifact, de-risks the most (toolchain, frames, units, ephemeris, plus the `dr_p/dv_a` leverage check) before any force or control.
- **B. Jump straight into perturbations or control.** Riskier; debugs toolchain and physics (or control) at once.
- **Lean / settled (2026-06-03):** A. Start at Rung 0, proceed force-by-force through Stage 1 (2a→2d), then Rung A.

---

## Appendix A — Key numbers

| Quantity | Value |
|---|---|
| mu (Earth) | 398,600 km^3/s^2 |
| Earth radius | ~6,378 km |
| Earth Hill radius | ~1.496e6 km |
| Orbit periapsis (50 km alt) r_p | ~6,428 km |
| Interception altitude | 200 km (during descent, before periapsis) |
| v_at_interception (any high-e apogee) | ~10.8–11.0 km/s |
| Terminal delta-v budget (400 g, Isp 200 s, 25 kg) | ~32 m/s |
| Terminal descent 600 → 200 km | ~5 minutes |
| Post-impact aerobraking onset | ~120 km (PuffSat burns up at ~50 km; intentional) |
| Leverage dr_p/dv_a = v_a·(r_a+r_p)^2/mu | ~250 (0.9 Hill), ~72 (lunar), ~30 (150k km) km per m/s |
| Required apogee velocity accuracy (few km perigee) | ~1–2 cm/s at 0.9 Hill, looser at closer apogee |

## Appendix B — Paper cross-references
- Orbit / deployment: `sec:leo_orbit_details`, `sec:formation_challenges_current_missions`
- Guidance / estimation: `sec:neural_navigation` (UKF, MPC, KalmanNet, RL warm-start)
- Drag / propellant: `sec:formation_challenges_current_missions`, `sec:estimate_cold_gas` (GOCE-derived ~400 mN, ~400 g)
- Debris / perigee targeting: `sec:handling_space_debris`
- Coordinator nodes (measurement/computation off-board the PuffSats): `sec:coordinator_node_dry_mass_disposal`
