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
- **Target-rocket RCS fast loop** (the "react within fractions of a second after each impact" problem). Separate sub-system. (ADR 0015 records its architecture note: **aim-bias primary** — the next scheduled hit cancels accumulated angular momentum — plus a non-toxic kN-class trim engine; the sim's deliverable to that ledger is the offset distribution and its correlation structure, not a coupled loop.)
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
- **Terminal-level (re-specced 2026-06-11, ADR 0015):** plate *capture*, not cm-centering — **σ_lateral ≤ 1.65 m against a 5 m plate radius (≥99% capture) + time-of-arrival ≤ ~10 ms** at closest approach. The historical ~5 cm centering target is **retired** (and its `rel_tol` 1e-13 / Encke fidelity prerequisite with it); capture probability is an economics output (misses burn up by design), not a pass/fail gate.

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
- **SRP is not "easy."** It depends on A/m, attitude, and optical properties of a 25 kg sat that may spin, deploy bladders, or have its albedo deliberately altered by lasers. Over a long coast it deposits ~0.5 m/s of delta-v at the highest-leverage point (apogee). It is the **second hardest force after drag**, and at apogee possibly the first. Eclipse entry/exit is a discontinuity (event). First-build idealization: a fixed effective `Cr·(A/m)` on a sphere-equivalent (cannonball), not a cylinder; saying "cylinder" implicitly signs up for an attitude model (orientation relative to Sun and velocity), a separate subsystem deferred to a later rung (the optimistic ANFO-cylinder model is **Rung E**, ADR 0009; the lumped `Cr·(A/m)` / `Cd·(A/m)` are deliberate *conservative* placeholders, and the propagator's 1 kg mass is normalization — the coefficient, not the mass, is the only drag/SRP lever).
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
1. **Discrete midcourse corrections during the coast** (event-triggered, e.g. one near apogee, one mid-descent) null the bulk of the dispersion *cheaply*, where 1 m/s buys hundreds of km of perigee. Standard statistical-midcourse + terminal guidance, as in interplanetary practice. (Caveat from Rung A — ADR 0006: the **mid-descent** burn buys ~0 *along-track* authority; along-track/phase error is apogee-bound. The mid-descent correction earns its place for **perigee/radial** trim via the `dr_p/dv_a` lever and, once estimation enters at Rung C/D, for **observable coast drift** — not for along-track timing.)
2. **Continuous gimballed terminal burn 600 → 200 km** handles drag rejection and the final few km of aim. It must be *handed* a near-nominal state by the midcourse corrections.

No continuous burn is needed above ~600 km (drag negligible); only a small number of discrete corrections.

> **Corrected 2026-06-10 (ADR 0014):** the "few km" terminal authority above is
> propellant-side and does not survive the 400 mN actuator (ADR 0004) — delivering
> 32 m/s in the window needs ~2.7 N, and the 600→200 km descent is ~170–180 s
> (vis-viva), not ~5 min. The binding constraint is **thrust**: lateral authority
> ½·a_max·t² ≈ 240 m from 600 km, ~450 m starting at the 800 km hand-off. The working
> catch radius is therefore **~400 m** (aim burn starting at the hand-off), the
> km-scale dispersion tails are carried by a **high-node impulsive trim** (~30,000 km,
> where A2's measured table shows affordable authority — not the mid-descent node), and
> the terminal burn handles drag (a *meters*-scale effect, B3a) plus the last few
> hundred meters of aim down to a relative-GNSS-limited ~1 m endpoint.

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
- **Train mode (ADR 0016, Rung D):** for train-relative claims the axes split into **shared draws** (coefficient *model bias*, F10.7/Ap drivers, deployer systematic — one draw per train) vs **per-unit draws** (nav noise, injection scatter, coefficient unit spread). Deliverables split accordingly: centroid drift (absorbed by the plane's centroid retarget, vs the ±~2 km spec) and scatter about the centroid (charged to the catch radius/plate, vs ADR 0015). The single-PuffSat ensembles elsewhere in this section sample all axes independently — valid for per-unit dispersion, silent on the split.

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
- **Rung B — actuator realism.** Replace impulsive Δv with finite gimballed burns: thrust magnitude, Isp, mass depletion, dead-time / latency. Still perfect state. Question: can the real ~400 g actuator deliver the needed correction in the minutes available, especially the terminal descent? **Decomposed below (B0–B2, ADR 0008).** Early inputs (§5): model a fast magnitude loop (PWM / pulse-frequency, minimum-impulse-bit limited) plus a slow direction loop (gimbal/attitude rate-limited ~1 deg/s); cold-gas micro-pulses above 800 km, sustained pulsed thrust below; MEMS reserved for attitude, not translation. **100 Hz modulation feasibility (settled 2026-06-03):** clean for cold gas (no combustion transient; fast solenoid/piezo valves ms-class, or a *proportional* flow valve for smooth thrust with no MIB quantization), but marginal for chemical mono/biprop (ignition + chamber-pressure + catalyst-bed dynamics cap clean pulse-mode at ~tens of Hz). This is not a problem: drag rejection closes at ~10–20 Hz (chemical's comfortable range), and only the deferred cm-aim wants true 100 Hz, which is a cold-gas/electric fine stage anyway. Favor a **coarse/fine split**: chemical (or higher-thrust) carries the steady gross anti-drag burn (Isp matters, most Δv), cold-gas micro-pulses do the fast 100 Hz trim (small Δv, low-Isp penalty negligible). Watch-item: Appendix A assumes Isp ≈ 200 s, optimistic for true cold gas (N₂ ~70 s, He ~165 s); pin per-device Isp when sizing this rung, since the cold-gas portion may run below 200 s.

  **Settled 2026-06-08 — first-cut Rung B actuator model (ADR 0004).** A *single omnidirectional proportional cold-gas thruster*, not the coarse/fine chemical+cold-gas split above (that remains a deferred higher-fidelity option):
  - **Paradigm:** proportional throttle is the modeling *abstraction*; the physical realization is a bang-bang cluster averaging to "continuous" (§5). No MIB quantization in the model — the cluster's residual minimum-impulse-bit is a deferred refinement, and lands exactly where the cm-aim 100 Hz trim wants it.
  - **Thrust:** 400 mN max (GOCE-derived, Appendix A); ~5 mN floor (proportional ~80:1 turndown; the cluster duty-cycles below it); continuous between, ms-class response (instantaneous vs. the ~5-min terminal burn).
  - **Direction:** omnidirectional, direction-slew rate-limited to ~1 deg/s (§5). The required sweep is only ~0.1 deg/s near perigee, so the limit carries ~10× margin and is expected slack — a result to confirm, not assume. Thruster geometry / cosine losses / gimbal-vs-whole-body realization deferred until the slew limit is shown to bind.
  - **Isp:** a *reported sweep*, not a single value — Isp is a post-processing lever on the propellant ledger (`propellant ≈ Δv/(Isp·g₀)`), not on controllability, so the A-rung Δv is computed once and transformed. Anchors: 50 s (conservative — CO₂/argon/refrigerant-class storable, dense + self-pressurizing), 70 s (N₂ realistic), 200 s (Appendix A optimistic). At 50 s the §9 ~32 m/s budget costs ~6.4% of 25 kg (vs ~1.6% at 200 s), so the curve shows the paper's <2% propellant claim failing ~3× at the conservative anchor — a deliberate finding, with higher mass fraction the accepted trade.
  - **Mass depletion:** modeled (Tsiolkovsky; Orekit maneuver models deplete mass via Isp natively). ≤~6.4% acceleration rise across a burn at worst; makes total Δv → propellant mass a first-class output feeding the Isp sweep.
  - **Dead-time:** deferred to Rung C (§16.8). Valve/command latency is ms-class — negligible vs. minute-to-day maneuver timescales — and the latency that bites (off-board coordinator-node sense→compute→uplink→actuate) has no source under Rung B's perfect state; carried as a ~0-default parameter so the seam exists.
- **Rung C — estimation in the loop.** Feed the controller a noisy state instead of truth: a Gaussian placeholder first (covariance motivated by what coordinator-node ranging can plausibly deliver, not an arbitrary blob), then the UKF estimating `Cd·(A/m)` and `Cr·(A/m)`. Question: does the loop still close with realistic knowledge? **Decomposed below (C0–C4, ADR 0010).**
- **Rung D — MPC + Monte Carlo.** Swap the deterministic targeter for MPC (Python; neural second pass later); run N = 50 as a smoke test, then N = 10^3 – 10^4 for the result with a resolved perigee tail. *Detailed decomposition deferred.*
- **Rung E — high-fidelity shape + attitude (optimistic capstone; ADR 0009).** A–D run the *conservative cannonball* (lumped `Cd·(A/m)` / `Cr·(A/m)`); Rung E re-runs the **Rung D Monte Carlo** with the real ANFO-cylinder shape model — same trajectory and seeds, a clean A/B — to bracket the truth and quantify the margin the conservative bound leaves. Narrative: *"even cannonball-pessimistic it works — here is the optimistic, more accurate model, and the comparison."* **E1 (committed):** attitude-dependent area at the favorable orientation (face-on into the flow / projected-to-Sun), *assuming* good pointing — an Orekit panel/box drag+SRP model + attitude provider swapped into D's harness (a force-model swap, not a new harness). Deliverable: the **pessimistic-vs-optimistic comparison**, the paper's follow-up to the GOCE-ANFO appendix (`sec:estimate_cold_gas`). **E2 (optional):** closed-loop attitude *pointing feasibility* — proves the assumed orientation holds against the destabilizing aero torque on a slender body (the predictable reviewer question); a **light** torque-margin / static-stability analysis first (the attitude analog of the paper's drag back-of-napkin, appendix-grade), the **heavy** rigid-body attitude-dynamics + control sim only if the margin is thin. **Mass/coefficient note (ADR 0009):** the propagator's 1 kg is *normalization* — the lumped coefficient is the only drag/SRP lever — so Rung E changes the *coefficient/shape model*, never the propagator mass.

#### Rung A decomposition (controllability core)

Goal: with perfect knowledge and idealized thrust, is there enough control authority to hit the perigee target, and at what Δv? This answers the headline question in its purest form; B–D each erode the margin Rung A establishes.

- **A1 — single impulsive midcourse.** On the nominal forces from Stage 1, inject a known apogee error (e.g. 0.5 m/s). One impulsive midcourse correction at the apogee deployment node, solved by a differential corrector (vary the correction Δv to null the predicted **interception miss** — the 3-component RTN position miss at the 200 km EME2000 crossing, §16.3), restores the nominal crossing. Confirm the resulting perigee shift matches the `dr_p/dv_a` prediction (an **acceptance cross-check**, not the objective — perigee is a diagnostic, not the target) and the propellant ledger closes (Δv applied = Δv commanded). This is the old "Rung 1" check, now with perfect state and no filter. (Corrector objective settled by ADR 0003; supersedes the earlier "null predicted perigee error" wording.)
- **A2 — two-burn midcourse: along-track is apogee-bound (finding).** A2 added a second coast burn to A1's apogee impulse to test the A1 carry-forward belief that a **mid-descent** burn (900 km, §16.6) nulls the along-track/phase tail cheaply. **Settled 2026-06-09, then overturned by evidence — ADR 0006 supersedes ADR 0005:**
  - **The 900 km "game-changer" presumption is wrong.** A node-altitude sweep on the hard run (seed 20260608, run 1: the 2.4σ radial injection → ~28 km along-track miss A1 cannot null sub-budget) shows a mid-descent second burn adds **~0 along-track authority** — at 900 km / 5000 km the miss stays ~25 km and the Δv just hits the cap. Authority is concentrated **near apogee** and falls monotonically as the node rises (80 000 km nulls at 55 m/s vs A1's ~88; 140 000 km trends to ~31, iteration-limited). Sweep table + recipe in ADR 0006.
  - **Mechanism.** Along-track miss is a phase/timing error; the only Δv-efficient fix is a tangential re-phase whose shift accrues over the *remaining* arc — ~zero at 900 km, a full descent near apogee. And at Rung A's perfect state there is **no new information** between the apogee burn and a mid-descent node, so a mid-descent burn cannot beat the apogee burn. The near-apogee second burn that *does* help is operationally part of the apogee maneuver, not a distinct midcourse capability.
  - **Conclusion: A1 accuracy dominates this phase.** There is no cheap mid-descent rescue for a timing error; the apogee injection + apogee correction *is* the along-track lever. The classical statistical-midcourse second burn (§9, §16.6) earns its place only at **Rung C/D**, where coast drift becomes *observable*, and/or for non-along-track modes via the `dr_p/dv_a` lever (§8).
  - **Code disposition.** The pure minimum-Δv least-norm solver (`solve_two_burn_correction`: `lstsq` on a 3×6 finite-difference Jacobian, RTN-6→position; the ADR 0005 mechanics) is **kept** as the auditable probe and an A3 seed — a tested library primitive, **not** wired into `run_ensemble`. The mid-descent two-leg harness glue was **reverted** (it encoded the disproven node). Both A2 burns stay impulsive/perfect-state; the vestigial §9 terminal aim still nulls to zero at Rung A and becomes a finite drag-rejecting burn at Rung B.
- **A3 — deterministic coefficient sweep: controllable everywhere, the wall is conditioning not authority (implemented, ADR 0007).** Hold the targeter fixed; sweep `Cd·(A/m)` and `Cr·(A/m)` deterministically (factor-space, straddling nominal, injection zeroed) and map required Δv vs coefficient error — the controllability map, before any noise or MPC. **Settled and built 2026-06-09 (ADR 0007); verified live:**
  - **The perfect-model map is controllable everywhere at Δv « budget.** A transverse apogee Δv has large authority over the along-track crossing (~2×10⁵ m per m/s), so 0.5–2× coefficient dispersion is cheaply nulled by the single apogee burn — far under even the conservative 50 s-Isp budget (~9.8 m/s). On the *coefficient* axis the paper's <2% propellant claim holds comfortably; the over-budget / uncontrollable regions belong to the **injection** axis (A1/ADR 0003) and to Rung C's unknown drag, not here.
  - **The map is ~1D in `Cr`.** Drag is negligible at a 150,000 km apogee; the small required-Δv gradient is driven by SRP (`Cr·(A/m)`) over the multi-day coast. The `Cd` axis is flat.
  - **The apparent authority boundary was corrector conditioning, not physics (resolves ADR 0003 finding 3).** The first sweep read all-non-converged off-nominal with ~0 Δv — but the 200 km target is an *altitude event*, so the radial/altitude crossing component is pinned and that Jacobian direction is near-singular (cond ~10⁷). The fix is the in-house LM damping (ADR 0007 decision 3): `λI` regularizes the pinned direction so the corrector nulls the well-conditioned in-plane/normal miss cheaply (plain Newton instead wastes ~2.5 m/s fighting it); `lm_lambda` re-tuned to `1e-6`. A **2-DOF retarget** (drop the pinned radial) is the parked structural alternative.
  - **Code.** Pure `sweep.py` (`SweepSpec`/`grid_inputs`/`to_grid` + σ-equivalent & Isp-budget overlays); the LM-damped + ToA-gated `solve_apogee_correction`; the JVM `run_sweep` harness. No resume sink / `make` target yet (deterministic, fast once converged — ~4 min for 3×3).

Targeter choice (settled 2026-06-03): a deterministic differential corrector, not MPC, for Rung A. It yields clean controllability ground truth, a Δv floor and an authority boundary. MPC enters at Rung D and is checked against this baseline; were MPC to lead at Rung A, a miss could not be attributed between an uncontrollable orbit and a mistuned controller.

#### Rung B decomposition (actuator realism)

Goal: with perfect knowledge, does the real ~400 mN proportional cold-gas actuator deliver the needed correction in the time available — and what does it cost in propellant? Each slice erodes the perfect-impulse margin Rung A established. Decomposed and the surprising calls settled 2026-06-09 (**ADR 0008**); a session profile found the terminal step-cap (not parallelism) was the A3 integration-test slowness.

- **B0 — terminal-propagation regime-switch (do first; implemented 2026-06-09, ADR 0008).** The §6.2 fragility: the adaptive integrator oversteps the 200 km event below the surface (`point is inside ellipsoid`) whenever drag is too weak to force a small step there. The §6.2 fix is the hand-off: **coast on the big adaptive step (600 s) → 800 km altitude event → terminal phase on a tight cap (30 s)**. The terminal cap equals the prior global 30 s cap (provably safe for every orbit), while the coast runs at 600 s — recovering the tax it cost. Verified: byte-identical perigee/ToA (common-mode preserved), capstone N=50 4m18s→3m20s wall, and the closed-loop corrector tests pass. **Measurement arc (full detail in the ADR 0008 impl-note):** an open-loop probe first argued for a simpler one-constant cap-loosen (a global 300 s cap — bit-identical, ~2.5× faster, safe to a ~575–600 s sub-surface cliff on the −4σ low-drag tail), but the **full integration suite falsified it** — the corrector probes large re-phasing Δv whose wild-perigee orbits overstep at 300 s where 30 s held, exactly §6.2's "single pure-adaptive scheme chokes on terminal discontinuities." A single global cap cannot be both fast in the coast and safe in the terminal for the orbits the corrector explores; the hand-off decouples them. **Refinement:** the terminal phase stays *adaptive* (capped DP853) here; the *fixed-step* Cowell terminal defers to **B3**, where the continuous burn needs the deterministic control-clock cadence. **Lesson:** measure the closed-loop path, not just the open-loop capstone, before settling a propagation-architecture call. Check (met): dispersed `full_force` descents — open- and closed-loop — reach 200 km without overshoot, faster, same perigee. (Process-level multiprocessing stays a **Rung D** tool, §11.2/§12 — it only parallelizes the multi-scenario axis, and the ~40 s per-worker JVM cold-start caps its small-N payoff.)
  - **5 cm forward-compat:** the eventual cm-trim aim does *not* reopen the coast. Numerical accuracy is already ~mm everywhere (relTol 1e-10; the cap only ever forces *smaller* steps), so 5 cm is well inside current accuracy. The driver of the 5 cm aim is **control cadence below 800 km** — the 100 Hz fine-trim loop (§5 / ADR 0004) and fixed-step stepping on the control clock — which plugs into this same 800 km hand-off seam (B3 builds the fixed-step terminal; its step shrinks 30 s → ~10 ms). The coast above 800 km only needs the right ballpark for the terminal aim, and is untouched. So the regime-switch *architecture* is forward-compatible; only the terminal *integrator* swaps.
- **B1 — finite mass-depleting burn (the Actuator; implemented 2026-06-09, ADR 0008).** Map a commanded impulsive Δv (the *unchanged* A1 corrector) → a finite Orekit `ConstantThrustManeuver` at the apogee node: 400 mN, Isp, Tsiolkovsky mass depletion. **This deliberately breaks the A-rung `predict == execute` identity** (the corrector predicts impulsive, truth executes finite); the residual interception miss *is* the actuator-realism erosion B1 measures. **Mass convention (ADR 0008):** the burn thrust is scaled to the propagator's fictitious 1 kg (`F·m_p/m_wet`), so `a=F/m` and the burn duration match the real 25 kg / 400 mN actuator while the lumped Cd·(A/m) / Cr·(A/m) keep scaling drag/SRP; a sentinel Isp keeps the executed arc constant-mass, so the *trajectory* is Isp-free and propellant is the pure Tsiolkovsky transform (ADR 0004 decision 2). **Finding (seed 20260608 run 0, Δv≈2.17 m/s → ~136 s burn):** the erosion is **~89 m, almost entirely along-track** (the dr_p/dv_a lever; ToA shifts ~0.18 s) — the burn centroid landing ~68 s past the apogee node. Small vs the km-scale open-loop dispersion and far above the 0.7 m impulsive residual, so a real measured second-order effect, *compensable* (a finite-aware targeter that centers the burn would null most of it). At the current km aim it is negligible; for the deferred 5 cm aim it is ~1780× too large, so that is where the finite-aware targeter (deferred, ADR 0008) earns its place. Propellant for that correction: 0.44 % / 0.32 % / 0.11 % of 25 kg at Isp ∈ {50, 70, 200} s. Slew (~1 °/s) stays unmodelled: a single fixed-direction apogee burn needs no slew (it binds only at B3's turning anti-drag burn). Real-mass depletion *coupled to descent drag* is deferred to B3 (its large anti-drag burn is the first consumer).
- **B3 — feedforward terminal anti-drag burn (actuator authority). Split into B3a (measure, done) / B3b (execute, deferred).** Under perfect state the apogee burn already nulls the aim including the *known* drag (A3), so the terminal burn is a deliberate propellant-suboptimal **cost baseline**, not a closed loop (that is C/D's job).
  - **B3a — measure (implemented 2026-06-09, ADR 0008 B3a impl-note).** Instrument the known-drag descent 600 → 200 km — sample the truth drag acceleration and reduce it to anti-drag Δv = ∫(drag accel)·dt, peak thrust = max │drag force│, peak direction-slew. **This needs neither an executed burn nor the fixed-step terminal** — it is a step-handler over the *existing* adaptive descent (ephemeris-sampled), so the §6.2 fixed-step Cowell terminal is *not* built here (see B3b). Pure reductions in `puffsat_sim/anti_drag.py` (`summarize_anti_drag`); the JVM `montecarlo.instrument_anti_drag` samples the descent. **Finding (nominal descent, conservative cannonball `Cd·(A/m)=0.04`):** peak thrust **16.7 mN** (~24× under the 400 mN limit), peak slew **0.048 °/s** (~21× under 1 °/s, and ~the paper's ~0.1 °/s sweep estimate), anti-drag Δv **0.015 m/s** → propellant **~1.5 g @ Isp 25 s (0.006 % of 25 kg)**. The paper's `sec:estimate_cold_gas` 374 g / 400 mN is a *deliberately stacked-pessimistic* upper bound (×10 area, ×20 surge, ×2 v², ×2 pulsing on GOCE's solar-min baseline); the physical NRLMSISE requirement on the fast eccentric pass is ~24× under the thrust limit and ~250× under the propellant figure, so **both gates pass with enormous margin and the <2% claim holds a fortiori** — even on the conservative coefficient (ADR 0009: no grounding triggered).
  - **B3b — execute (deferred, likely to C/D).** Build the §6.2 **fixed-step Cowell terminal phase** (deferred from B0) and actually fly the open-loop feedforward anti-drag burn to validate the cost end-to-end. Its real first consumer is the *executed / closed-loop* drag rejection at C/D (uncertain drag), so the fixed-step terminal moves there; B3a already answers the feasibility/cost question without it.
- **B2 — propellant ledger + Isp sweep (the headline deliverable; implemented 2026-06-09).** Aggregate the mission Δv (B1 apogee corrections + B3a anti-drag) → propellant fraction = Δv/(Isp·g₀) at Isp ∈ {50, 70, 200} s → the fraction-vs-Isp curve against the <2%-of-25 kg line. Pure ledger in `puffsat_sim/propellant.py` (`propellant_curve`), the ADR 0004 linear model — the **exact inverse of the A3 `budget_dv_m_s`**, so the ledger and the controllability budget share one 2% line. **Finding (reverses the prior expectation):** the anti-drag Δv does **not** dominate — B3a measured it at 0.015 m/s (0.006 % at 50 s, negligible) — and the corrections are « budget (A3: controllable everywhere at Δv « the 9.8 m/s 50 s-budget). So a representative correctable run (correction 2.17 m/s + anti-drag → mission Δv **2.19 m/s**) costs **0.45 % / 0.32 % / 0.11 %** at Isp 50 / 70 / 200 — all « 2 %, ~4.5× margin even at the conservative 50 s. The earlier "expect 50 s to fail ~3×" (ADR 0004/0008, which assumed a large anti-drag-dominated Δv from the paper's stacked-pessimistic 374 g) is **falsified**: B3a's NRLMSISE anti-drag is ~250× smaller, so <2% holds at every Isp. Because the 2 % line at 50 s *is* the authority budget, every correctable run sits under it by construction. The full propellant *distribution* over the dispersion is the **Rung-D** deliverable; B2 is the deterministic ledger Rung D feeds through. **Rung B is now complete (B0/B1/B3a/B2; B3b deferred to C/D).**

Order: B0 → B1 → B3a → B2 (B3b's executed burn + the fixed-step Cowell terminal move to C/D, their real consumer). Dead-time, MIB / cluster geometry, cosine losses, the coarse/fine chemical+cold-gas split, and the cm-aim 100 Hz trim remain deferred (ADR 0004).

#### Rung C decomposition (estimation in the loop)

Decomposed 2026-06-09 (ADR 0010). Rung C drops the perfect-state assumption and asks whether the §8 controllability survives *realistic knowledge*. It populates the predict/execute seam that already exists (ADR 0003): at A/B `predict == execute` (full-force truth); at C the corrector's onboard `predict` finally diverges from truth. Three independent sources of imperfect knowledge are added **one at a time** against the known-good Rung-B baseline, and — the load-bearing decision — each is led by a cheap **requirement** sweep, with the real filter built only to confirm **feasibility** on a few runs. **The giant Monte Carlo (a UKF in every trajectory) is Rung D, not C** — a requirements paper sizes *how good knowledge must be*, which is a covariance sweep, not a faithful-filter ensemble.

- **C0 — Gaussian nav-error placeholder (the navigation requirement).** Inject a Gaussian navigation error of covariance Σ at the apogee planning node; the corrector plans from the perturbed estimate and executes against truth; sweep Σ → the **navigation-accuracy threshold** for the interception. No filter, no ensemble (a small deterministic A3-style sweep). The requirement is velocity-dominant (the `dr_p/dv_a` lever; Appendix A's ~1–2 cm/s apogee-velocity figure). **TDD:** the nav-error injection (seeded Gaussian perturbation of the planning state) is pure → use **`/tdd`**; the threshold is a **measured** sweep, not asserted.
  - **Implemented 2026-06-10 (ADR 0011; `puffsat_sim/navigation.py` + `montecarlo.run_nav_sweep`).** Reframed during grilling: because the apogee corrector applies the *same* correction in predict and execute, it **cancels** nav error to first order — the residual interception miss is the apogee→crossing sensitivity **Φ** (a 3×6 STM) times the nav error, **uncontrollable at the apogee node**. So C0 is a deterministic per-component Φ-basis sweep (not a sampled ensemble), thresholded against the **terminal catch radius** (§9, parameterized), with Σ motivated by coordinator-node ranging — range-only multilateration + a few-gram omni transponder (decision 7, carried to C1). **Finding (61-cell sweep, 100% converged, residual linear across the full swept 10 km / 1 m/s range — `−Φδ` is rock-solid):** the requirement is **overwhelmingly transverse-velocity-dominated**. Lateral-miss sensitivity ‖Φ_TN‖ = **2.15×10⁵ m per m/s** (T-vel) vs 4.3×10³ (N-vel), 1.3×10³ (R-vel), 0.60 (R-pos), ≈0 (T-pos/N-pos — an along-track apogee displacement is a pure phase shift: same orbit, only ToA moves). Velocity beats position by ~5 orders; transverse velocity beats the other velocity axes 50–160×. **Binding navigation requirement:** apogee transverse velocity known to **~2.3 cm/s** for a 5 km catch radius (→ 4.7 mm/s @ 1 km, 0.47 mm/s @ 100 m) — confirming *and refining* Appendix A's ~1–2 cm/s (the binding axis is the along-track/timing amplification, not the perigee-radial lever). Secondary: phantom Δv (the corrector burning fuel chasing the unobserved error) up to ~1 m/s at the 1 m/s cell. Φ feeds **C1** (UKF covariance checked via `ΦΣΦᵀ`) and **Rung D** (sample nav error from it). *Numerical note:* trustworthy to the metre here — the truth model's `rel_tol=1e-10` floor (~cm at apogee scale) sits ~5 orders under this km/m-scale finding; a cm-level terminal end-state would need the tolerance tightening flagged in C3.
- **C1 — UKF state estimation (requirement reuse + feasibility).** *Requirement:* C0's sweep already gives the position/velocity-knowledge bar. *Feasibility:* build the classical UKF estimating position/velocity from the realistic measurement suite (coordinator-node ranging at apogee, GNSS + accelerometer near perigee — §16.4 fused suite) and confirm on a **few representative runs** that its steady-state error covariance sits inside the C0 threshold. The two navigation regimes — apogee (coordinator ranging, GNSS unavailable) vs perigee (GNSS) — are the §16.4 sub-study. **TDD:** the UKF math (sigma points, unscented transform, predict/update) and the measurement models `h(x)` (ranging, GNSS) are pure and check against analytic/linear systems → **strongly `/tdd`**; "the filter achieves the required covariance on truth" is a **measured** few-run finding.
  - **Approach settled 2026-06-10 (ADR 0012; grilling session).** C0's reframe (the binding axis is *transverse velocity*, 2.3 cm/s) turns the architecture question from "how well does multilateration fix position" into "how is transverse velocity observed" — directly (two-way carrier Doppler, mm/s LOS) and/or indirectly (ranges filtered through the dynamics over the days-long apogee arc); **both are modeled, and "is Doppler load-bearing?" is a measured output** (range-only sweep point), with cadence swept coarsely (0.003/0.03/0.3 Hz). **Coordinator nodes**: co-flying matched-`a` orbits (no secular drift) with **perigee ≥200 km** (reusable assets, never in the burn-up zone; apogee trimmed to keep the period), modeled as **known-ephemeris beacons** (inputs, not states — node nav error inflates `R`); node **geometry is a derived GDOP-style requirement** ("N nodes / this LOS spread / this Doppler quality"), never an assumed constellation. **Filter**: owned, fully typed UKF (no FilterPy — untyped, dormant, conditioning risk; verified against the analytic linear KF), 6-state Cartesian EME2000 (dimension parameterized for C2's +2), **pure-Python two-body+J2 onboard dynamics** (realistic flight architecture; the truth−filter gap is `Q`'s to absorb, and `Q` is a swept knob, never quietly tuned). **Method, mirroring C0's compute decision: a LinCov sweep engine** (the covariance recursion is measurement-value-independent — pinned to the reference it sweeps the 5-axis grid with zero truth propagation) **validated by a few seeded UKF runs on Orekit truth with NEES consistency**; Σ at the apogee node → `Φ Σ Φᵀ` vs the catch radius. GNSS: position-fix `h(x)`, *reported not swept* (descent hand-off with/without GNSS — "GNSS optional" is the stronger paper claim if it holds; ~30–50 g unlocked spaceborne receiver, COCOM-limited COTS can't track 10.8 km/s); accelerometer **deferred to C2** (it observes coefficients, not position). Sizing notes banked in the ADR: clock-free coherent-transponder link budget (coordinator-OCXO-bound, ~10× Doppler margin at 3000 km), GNSS mass, node constellation principle. **LinCov is also a Rung-D asset** — screen + control variate that shrinks the MC (plausibly 10–100× fewer runs), never replaces it (saturation, non-Gaussian drivers, tail-probability headlines). Layout: pure `puffsat_sim/estimation.py` + a `montecarlo.run_nav_feasibility` seam.
  - **Implemented 2026-06-10 (ADR 0012 findings; `puffsat_sim/estimation.py` + `nav_feasibility.py` + `montecarlo.run_nav_feasibility`).** Φ re-derived live (137 s, matches ADR 0011); 15-cell LinCov envelope ~50 s; seeded NEES validation against full-force truth. **Findings:** (1) **the NEES layer earned its place on the first run — the grilled SRP-scale `q = 5e-8 m/s²` was fiction**: claimed σ_Tvel ≈ 6 µm/s vs *actual* error 0.09 m/s (above the 2.3 cm/s requirement itself; NEES ~8×10⁸) — the truth−filter gap at apogee is dominated by **third-body tidal acceleration** (~3×10⁻⁵ m/s², three orders above SRP), and a LinCov-only C1 would have shipped the fiction. (2) The q-ladder localizes the consistency crossover at the tidal scale; the validated nominal is **q = 1e-4** (NEES 5.18, mildly *pessimistic* — the safe direction). (3) **C1 verdict: requirement met with ~35× margin at the honest q** — σ_Tvel = 0.66 mm/s vs 2.3 cm/s, lateral 1σ miss 141 m vs the 5 km catch radius, **all 15 cells MEET**, every swept degradation (100 m ranges, 5° cone, 3 nodes, 0.003 Hz) clears by ≥ an order; the 1 km radius clears ~7×. (4) **Doppler is load-bearing at the honest q** (the apparent range-only tie was a q=5e-8 artifact): range-only degrades ~2.7× to 1.8 mm/s and goes mildly NEES-optimistic (9.4 vs 6.19) — still meets ~13×, but the with-Doppler nominal is the validated configuration. (5) The **100 m catch radius is the marginal case**; the named upgrade path is third-body in the onboard filter dynamics (recovering ~2 orders of q), not better sensors — banked for C2/C3.
- **C2 — UKF coefficient estimation (requirement ≈ free + feasibility).** *Requirement:* largely a **post-processing of the A3 controllability map** — the residual interception miss ≈ (coefficient *estimation* error) × (A3 sensitivity), and A3 already mapped required Δv « budget for ±2σ coefficient error, so the coefficient-knowledge *tolerance* falls out with **no new propagation**. *Feasibility:* confirm the UKF (accelerometer giving drag observability) estimates `Cd·(A/m)` / `Cr·(A/m)` to within that tolerance, a few runs. **TDD:** the A3-reuse residual-miss computation is pure → **`/tdd`**; the coefficient-estimation feasibility reuses the C1 UKF (already TDD'd) plus a **measured** run.
  - **Reframed 2026-06-10 (ADR 0013; grilling session).** Two premises in the original feasibility leg are broken. (1) **Timing:** deploy ≈ apogee ≈ burn (B1) and along-track authority is apogee-bound (ADR 0006), so the coefficient knowledge the corrector acts on is the **ground prior, full stop** — no estimator can improve it before the only high-authority burn. The requirement verdict is therefore **tolerance vs prior** (area few %, material `Cr` ~10–20%), not "can a filter converge." (2) **Instrument:** the §16.4 "GOCE-style" accelerometer is not a PuffSat payload, and a flyable MEMS unit (bias instability ~1–10 µg) can **never observe SRP at 0.01 µg** — `Cr` is unobservable by accelerometer; terminal drag (~1700 µg peak) is trivially observable but is better used as a **direct force measurement** (the product `Cd·(A/m)·ρv²` that terminal control actually needs) → the accelerometer is **reassigned to C3** as a drag-disturbance sensor, and the `Cd`-parameter prediction question moves to C3 with it. **C2a** = the pure tolerance computation: residual miss ≈ `Φ_lat,vel · (∂Δv/∂c) · δc` → tolerance = catch radius / sensitivity, with `∂Δv/∂c` **measured** from 1D A3 cuts (the per-cell Δv *vectors* in `control_log`; never recorded by A3's scalar map) and cross-checked against the analytic SRP impulse (`Cr·(A/m)·P₀·t_coast` ≈ 0.011 m/s per 100% → ≤ ~2.3 km lateral per doubling → tolerance ~factor-2 at 5 km, prior covers ~10×); plus an explicit RSS **error-budget ledger** (C1 nav 141 m, B1 erosion 89 m, coefficient-prior residual, corrector floor) vs the catch radius. **C2b** = LinCov coefficient augmentation (+1/+2 states, SRP/drag in filter dynamics, `σ_Cr(t)` along the coast) — quantifies ADR 0006's MCC-2 "observable-drift correction" contingency for the day Rung D tightens the radius. Expected verdict: C2 closes as an A2-pattern finding — *in-flight coefficient estimation is unnecessary for the midcourse because the burn precedes any possible estimate.*
  - **C2a implemented 2026-06-10 (ADR 0013 findings; `puffsat_sim/coeff_requirement.py` + `montecarlo.coeff_requirement_report`).** Measured (172 s wall; Φ from the minimal C0 sweep + two 3-point 1D cuts under the unchanged A3/C0 LM corrector, all converged): `‖∂Δv/∂Cr‖` = **8.8×10⁻³ m/s per 1.0 factor** (analytic SRP impulse 0.0105; **measured/analytic 0.83** — the envelope bounds from above as derived), lateral sensitivity **745 m per factor** → **tolerance ±6.7 factor units @ 5 km; the 0.2 prior is COVERED ~34×**. `Cd`: 3.4 µm/s per factor (the corrector's tolerance floor), tolerance ~6×10⁴ — **quantitatively unconstrained**; A3's "~1D in `Cr`" is now a number. **Error budget: RSS 224 m vs 5 km (~22× headroom)** — `Cr` prior 149 m ≈ C1 nav 141 m > B1 erosion 89 m, no dominant residual. Radius scaling: the prior covers down to ~150 m catch radius (6.7× @ 1 km; not covered @ 100 m — but there the C1 nav lateral alone already exceeds the radius, so **the coefficient prior is never the binding constraint at any feasible radius**). **Verdict (the A2-pattern finding, confirmed): in-flight coefficient estimation is unnecessary for the midcourse.** C2b (LinCov `σ_Cr(t)`, the MCC-2 contingency quantifier) stays specced, not scheduled — built only when a tightened Rung-D radius or C3's terminal predictor consumes it.
- **C3 — closed-loop terminal drag rejection (the one closed-loop *dynamics* slice).** Where **B3b** and the §6.2 **fixed-step Cowell terminal phase** (both deferred from B0/B3) finally have their consumer.
  - **Approach settled 2026-06-10 (ADR 0014; grilling session).** The grill's headline: **§9's "few km" catch radius does not survive the 400 mN actuator** — it was sized propellant-side (32 m/s "deliverable" in the window needs ~2.7 N), and the 600→200 km descent is ~170–180 s (vis-viva), not ~5 min. Thrust-limited lateral authority ½at² ≈ **240 m from 600 km / ~450 m from the 800 km hand-off**; the binding constraint flipped from propellant to thrust. **Re-baseline: the aim burn starts at the existing 800 km hand-off and the working catch radius is ~400 m** (C3b's measured residual-vs-entry curve replaces the assumption); the C0/C1/C2 tables re-read at 400 m with honest 2–3× margins (T-vel 1.84 mm/s vs 0.66 achieved; nav 2.8×; Cr prior 2.7×; RSS 1.8×). **Sensing:** terminal drag rejection is a *meters*-scale problem (B3a's 0.015 m/s anti-drag Δv → ~1–2 m crossing displacement even uncompensated), so MEMS feedforward (bias 3–10 µg vs 1700 µg peak, commanded thrust subtracted) + GNSS feedback suffice; the inherited `Cd`-parameter prediction question **dies quantitatively** (100% error → 1–2 m, under the nav floor); no UKF in the loop — σ_rel state-noise injection. **Endpoint goal 1 m, conditioned on relative GNSS** (the miss is *relative*, §16.3): σ_rel swept {absolute 10 m, code-differential 1 m, carrier-differential 0.2 m} so "what 1 m costs in sensors" is a measured output; the claim pair is **catch radius ~400 m (entry, thrust-limited) + endpoint floor ~1 m (nav-limited)**. **Law: ZEM/ZEV** (`a = k·ZEM/t_go²` capped, ZOH; read this rung's "PID/LQR" as ZEM), cadence swept {0.1, 1, 10 Hz}; fixed-step Cowell ≤ control step, equivalence-pinned vs the proven adaptive-30 s config on an unburned descent first. **Tails: a high-node impulsive MCC-2 trim** (3σ of the 224 m budget ≈ 670 m exceeds any thrust-limited radius; A2's table puts trim authority at ~1.7 m/s per km @30,000 km, dead ≤5,000 km — the §16.6 "correction 2 at 800–1000 km" node is wrong for this job, the *role* is ADR 0006's observable-drift correction); **C3c** measures the authority + MCC-2 cost curves with the kept A2 solver; trim *scheduling* defers to Rung D. Slices: **C3a** fixed-step terminal + executed B3a feedforward → **C3b** ZEM closed loop (dispersed drag + σ_rel; the measured catch radius) → **C3c** curves. 1 m is measurable at `rel_tol=1e-10` (~cm floor); the **5 cm centering stays parked**.
  - **Re-specced 2026-06-11 (ADR 0015; grilling session).** The success criterion becomes **plate capture**: σ_lateral ≤ 1.65 m vs a 5 m plate radius (≥99%, 2D-Rayleigh R/σ = 3.03) + **ToA ≤ ~10 ms** at closest approach (the "along-track" axis is time-of-arrival; the miss vector is ⊥ v_rel by definition). **GNSS is deleted from the PuffSat**: midcourse stands on the measured C1 coordinator-beacon suite (Doppler required — range-only is dead at the 400 m radius, 386 m ≈ 1.04×), and terminal nav becomes a **target-side astrometric tracker** (beacon-vs-star-field differential astrometry; **σ_θ ≤ 10 µrad required / 5 µrad target**) feeding the PuffSat's ZEM loop over the crosslink. The homing floor is **σ_miss ≈ 2σ_θ²v²/a_max** (1.46 m @ 10 µrad; knowledge ∝ R loses to authority ∝ R², so plate knowledge inside ~50 km is unusable by the PuffSat — 0.17 m of authority left). C3b's σ_rel sweep is re-keyed: **σ_rel(R) = σ_θ·R, σ_θ ∈ {2, 10, 50 µrad}** (+ the constant 1 m point for ADR 0014 continuity); new pure plate-frame miss summarizer (2D ⊥ v_rel + ToA, capture-vs-radius curve; Rung D headline = P(capture)). The **5 cm park and the 1e-13/Encke prerequisite are retired**.
  - **C3a (= B3b):** build the fixed-step Cowell terminal and *execute* the B3a feedforward anti-drag profile as a real continuous burn on the control-clock cadence — **known drag, open-loop feedforward**; measure the executed residual (validates end-to-end what B3a only measured).
  - **C3a DONE 2026-06-11 (measured; ADR 0014 implementation findings).** Equivalence pin **5.5 mm** / <1 µs ToA (fixed-step Cowell, RK4-Cartesian @1 s, vs the proven adaptive-30 s on the unburned descent). Uncompensated drag displaces the crossing **8.5 cm — ~20× under the ~1–2 m estimate above** (drag concentrates in the final seconds before the crossing, where it has no time to integrate into position); the executed 180-command ZOH burn leaves **2 mm (rejection ~45×)**, below the pin itself. Plan ≈ B3a's profile (Δv 0.0145 m/s vs 0.015 trapezoid; peak 15.96 mN; slew 0.048 °/s) — **both ADR 0004 gates PASS on the *executed* burn**; propellant 0.0030 % of wet mass @Isp 50. Terminal drag rejection is a *centimeters*-at-the-crossing problem; C3b inherits the fixed-step + maneuver-segment machinery and the validated 1 Hz cadence. Mass depletion stays the B1 sentinel-Isp convention (ADR 0008's "large anti-drag burn" premise falsified at ~1.5 g).
  - **C3b:** add the PID/LQR magnitude **feedback** on the accelerometer residual, now rejecting **uncertain** drag (truth `Cd`/density ≠ the onboard model's — the §10.1 disturbance that C2 sizes); measure whether the feedback recovers the miss. **This is where the terminal burn earns its aim role** (§7 / ADR 0008).
  Unlike C0/C1/C2 this **cannot** be a static-covariance sweep: the drag disturbance is time-varying and the feedback *chases* it down the descent, so C3 genuinely spins up the closed-loop terminal integrator — but still on a **handful of representative runs** (nominal + a few worst-case drag realizations), not an ensemble. **TDD:** the fixed-step integrator configuration and the PID/LQR law (vs a synthetic disturbance) are pure/small → **`/tdd`**; "does the loop close / the executed residual" is **measured**.
  - **Numerical-fidelity prerequisite (flagged from C0, 2026-06-10; binds the parked cm-aim, not the few-km drag rejection).** Pushing the terminal aim to **cm-level fidelity** requires tightening the truth-model integrator tolerance first. The current `rel_tol=1e-10` (`puffsat_sim/propagator.py`, DOP853) gives a per-step local tolerance of ~1.5 cm at the 1.56×10⁸ m apogee scale and a ~cm–dm global crossing-position floor — ~5 orders under the km/m aim (so C0/A/B and the few-km C3 drag rejection are unaffected), but it would **swamp a cm-level end-state signal**. The fix is **not** more than float64: representation is hugely adequate (ULP ~35 nm at apogee, ~1.5 nm at the 200 km crossing), the small-force disparity is safe (drag/SRP ~10⁻⁵–10⁻⁶ of gravity, ~10 orders above the 10⁻¹⁶ roundoff floor), and Orekit's two-part `AbsoluteDate` already handles the long-span/fine-step *time* disparity to sub-ns. The fix is **tightening `rel_tol` to ~1e-13** (float64's usable floor is ~1e-14 before roundoff dominates → sub-mm global), and likely switching the eccentric coast to an **Encke / element-based formulation** (equinoctial / DROMO) so the slowly-varying integrand keeps truncation+roundoff well under cm. **Verify before trusting any cm result:** the Rung-1 energy/angular-momentum conservation drift bounds the floor (`ε_E·a` ≈ 8 mm at `rel_tol=1e-10`), and a `rel_tol` 1e-10-vs-1e-12 crossing-difference plus a reverse apogee→perigee→apogee round-trip closure size it empirically in minutes. **TDD:** the conservation/round-trip checks are pure assertions → **`/tdd`**; the tolerance/formulation choice is a **measured** convergence study. **Retired 2026-06-11 (ADR 0015):** the cm-aim died with the plate-capture respec — 1.65 m is measurable at the current 1e-10 floor. Kept for the record; do not schedule.
- **C4 — control-loop latency (§16.8).** Model the dead-time as a lumped delay **per loop**, because the comms lag and the actuation lag live on different loops: comms round-trip + coordinator compute on the **slow outer loop** (midcourse/replan), onboard sensor + valve on the **fast inner loop** (C3b's terminal tracking). What decides whether it bites is the loop **bandwidth**, not the sample rate: dead-time erodes phase margin as `ω_c·τ`, and the drag-rejection bandwidth is ~1 Hz (drag varies over the ~3-min descent) even though the inner loop *samples* at 100 Hz — so a tens-of-ms τ costs single-digit degrees of phase. Deliverable: a budget table (distance/c ≈ 6.7 ms one-way at 2000 km + ~ms compute + ms-class valve), the `ω_c·τ` phase-margin check, and a τ-sweep on the C3b loop confirming insensitivity up to τ ≫ the budget. Jitter (§16.8 C) only if the sweep shows sensitivity. C4 **depends on C3b** (no loop to destabilize before the terminal feedback exists). **TDD:** the dead-time mechanism (a buffer holding the command stale for τ) is pure → **`/tdd`**; the budget table + τ-sweep tolerance is a **measured** back-of-envelope.

**Onboard-compute feasibility note (a paper appendix, not a sim rung).** The fast drag-rejection loop runs *onboard* (it must — comms can't sit in a fast loop) and is near-trivial: read the accelerometer, subtract the B3a feedforward (a table lookup), one PID/LQR update → ~10⁴ flops/s, a sub-gram MCU. The UKF (~10⁵–10⁶ flops/s on 8 states at the sensor cadence) and the Rung-D MPC (a small QP at ~1 Hz) fit a ~Cortex-M7 (~grams, hundreds of mW) **or** off-load to the coordinator. So **50 g of compute is ample for the classical first pass**; PuffSats being single-use and burn-up lets them fly COTS (non-rad-hard) parts. The deferred **neural** augmentation (KalmanNet / RL, CLAUDE.md second pass) is the item that would first stress a tiny budget. This belongs in the paper as a sizing note alongside the propellant/actuator estimates, not as a simulation step.

**Controller at C (vs D):** the **deterministic A1/A3 corrector** stays (fed an *estimated* state, no replanning) for the midcourse, plus the **dumb fixed terminal feedback** (C3b's PID/LQR); **MPC is held for Rung D**. Same rationale as §16.6 for Rung A: hold the control law fixed and transparent so a miss is attributable to the *knowledge quality* under study, not to controller cleverness — then Rung D introduces MPC and measures its value against this C baseline. **Success criterion:** the interception miss stays under threshold despite nav error (within C1's covariance) + model mismatch (within C2's tolerance) + latency (within C4's budget) — i.e., the §8 controllability survives realistic knowledge.

Order: C0 → C1 → C2 → C3 (C3a → C3b) → C4. The big Monte Carlo *distribution* (and MPC) are Rung D; the optimistic real-shape comparison is Rung E.

### Implementation queue (consolidated 2026-06-11)

One checklist for what remains; the named ADRs hold the authoritative specs — this list only
points. Everything above C3a is **done and unaffected** by the ADR 0015/0016 respec (both ADRs
state the measured A/B/C0–C2a results stand; the gate is green as of this entry).

- [x] **C3a (= B3b)** — fixed-step Cowell terminal arc, equivalence-pinned vs the proven
  adaptive-30 s config on an unburned descent; execute the B3a feedforward as a real burn
  (ADR 0014). **DONE 2026-06-11:** pin 5.5 mm; drag displacement 8.5 cm (~20× under the
  estimate); executed residual 2 mm (~45× rejection); both gates PASS (ADR 0014
  implementation findings).
- [x] **C3b** — ZEM closed loop under dispersed drag; σ_rel re-keyed to **σ_rel(R) = σ_θ·R**,
  σ_θ ∈ {2, 10, 50 µrad} + the constant 1 m ADR 0014 continuity point; new pure
  **plate-frame miss summarizer** (2D lateral ⊥ v_rel + ToA) and capture-vs-plate-radius
  curve (ADR 0014, 0015). **DONE 2026-06-12:** measured catch radius 500 m; 10 µrad grade
  capture-grade (RMS 1.07 m vs σ ≤ 1.65 m), binds by 50 µrad; **noise discipline
  load-bearing (2 orders)** — 3σ gate + 35 s track window + 45° firing-lag hold tame the
  σ_θ·R rectification (ADR 0014 C3b findings).
- [x] **C3c** — thrust-authority + MCC-2 cost curves with the kept A2 solver; MCC-2 sized
  against the *independent* tail only (ADR 0014, 0016). **DONE 2026-06-13:** ½·a_max·t²
  funnel validated (800 km → 487 m model vs 500 m measured); the 671 m tail (full C2a
  budget 224 m, 3σ) is reachable by funnel growth at only ~942 km burn-start but at the
  ceiling Δv ~4.63 m/s, vs a **0.12 m/s** high-node impulsive trim (**~38× cheaper** —
  MCC-2 vindicated); the out-of-plane trim lever stays finite at low nodes (unlike A2's
  along-track lever), scheduling defers to Rung D (ADR 0014 C3c findings).
- [x] **C4** — dead-time budget table + τ-sweep on the C3b loop (`puffsat_sim/latency.py`, pure,
  no JVM; §16.8). Terminal budget **7.3 ms** (sensor+crosslink+compute+valve), midcourse **70 ms**
  but discrete (no phase loop); ω_c·τ = **2.64°** at 1 Hz — negligible. Noiseless homing sweep flat
  through ~1 s (1 tick ≈ 136× the budget), degrades at 2 ticks; budget flat at 100 Hz too (the ~1 s
  tolerance is dynamics-set, not cadence-set). Read as relative degradation (double-integrator
  stand-in, ~2× hot vs C3b's Orekit loop); jitter/measurement-dropout knob unbuilt (no fragility
  shown); combined offset×noise tail stress deferred to Rung D. **Closes the C-rung.**
- [ ] **Rung D** — the feasibility Monte Carlo, decomposed (**ADR 0018**). **Pre-D gates
  (blocking):** σ_θ **tracker budget** (pure: what 10 µrad demands + acquisition FOV vs the
  hand-off Σ) — **DONE 2026-06-13 (`puffsat_sim/tracker_budget.py`, pure, no JVM): GATE
  PASS.** A conservative point (5 cm aperture, 1 ms exposure, 1 W beacon @ 1064 nm, beam
  ±2 mrad, η 0.3, nav-grade gyro, bench-calibratable 3 µrad focal-plane distortion) achieves
  **σ_θ = 3.2 µrad RSS — 3.1× under the 10 µrad requirement, inside the 5 µrad target**. The
  active beacon gives SNR ≈1670 at 300 km, so the budget is **calibration/jitter-limited
  (the distortion floor dominates), not photon-limited** — the "dim, fast target" worry is
  dissolved by making the target active. Homing floor at the achieved grade **0.15 m ≪ 1.65 m**
  (the bare 10 µrad requirement → 1.45 m, ADR 0015's thin-margin reference). Acquisition: the
  ±2 mrad beam covers the **±1.4 mrad** (3σ·141 m / 300 km) acquisition cone; the **binding FOV
  is reference-star availability** (±5.8 mrad for 3 stars at ~10th-mag density), resolved to
  Nyquist by a **~1100-px detector at 10.6 µrad/px**. The load-bearing terminal-nav grade is
  now a *derived* hardware requirement, not a guess. *Remaining pre-D gates:* torque-margin
  (the 1°/s slew rail); truth-validation **Tier 1** (conservation / tolerance-halving on the
  nominal coast) + **Tier 2** (independent Python conservative-force coast cross-check).
  **D1 — feasibility gate on the C baseline**
  (corrector + C3b ZEM + C3c MCC-2 + finite burn): train-mode `DispersionSpec`
  (shared-vs-per-unit, correlation pins **swept**); nav Σ a swept axis parameterized by node
  count → **minimum node count**; nav error sampled from the C1 Σ (not a live UKF); Φ-Jacobian
  warm-started quasi-Newton (FD-Newton fallback) + process parallelism; **importance-sampling
  tail** + brute-force validation batch, LinCov as pre-screen / control-variate / IS-designer
  (never replaces the tail). Headline **P(capture)** about the train centroid + centroid-drift
  (vs ±2 km) + scatter (vs plate) + propellant (<2 %) + perigee diagnostic + per-axis
  sensitivities → the **conditional** yes/no. **D2 — MPC value:** MPC vs the C baseline on the
  same MC, robustness only on a measured violation (§16.10) (ADR 0012, 0015, 0016, **0018**).
- **Deferred rungs:** **E** cylinder shape comparison (ADR 0009); **F** GMAT full-force
  cross-check of the nominal trajectory — the Tier-3 truth validation, run headless
  batch-script → report → compare, *not* the conda/Python-API path (ADR 0018).
- **Parked, build only on consumption:** C2b LinCov `σ_Cr(t)` (ADR 0013).
- **Retired, do not schedule:** 5 cm centering; `rel_tol` 1e-13 / Encke study;
  relativity-in-filter (ADR 0015).

---

## 14. Engineering practices

1. **One parameterized code path, not a separate "single mode."** The single-trajectory run is the ensemble harness with `N=1`, noise flags off, stochastic inputs frozen to nominal. Knobs: `N`, `noise_on`, `mismatch_on`, `seed`, over a single scenario definition. Prevents the quick-mode and the real-MC from drifting apart.
2. **Seed per run, replayable.** Any ensemble run (including a failing run #7,432 at N=10^4) must replay in single-trajectory mode from its seed, with dense logging. Single-trajectory and Monte Carlo share one draw mechanism.
3. **Mode-dependent logging.** Single runs: full 100 Hz state + covariance + residuals + event log + commanded/applied Δv. Big MC: per-run summary only (perigee achieved, miss, propellant, pass/fail, seed) or storage drowns.
4. **Apogee as a tunable parameter** (Section 3), so the dispersion-vs-apogee trade is a sweep, not a recompile.
5. **One slice, one module — on both sides of the master seam (ADR 0017).** A module is split when more than one ADR owns it (its CLAUDE.md entry is the test; logic lines are the secondary tripwire — suspicious ~400, act by ~600), and the cut must land on a Seam (CONTEXT.md): ownership says *when*, seams say *where*. JVM-side glue mirrors the pure slice modules in `puffsat_sim/runs/`. Functions: suspicious at ~50 branching logic lines (~80 for linear recipes); ≥3 default-off behavior knobs is the hard tripwire — bundle into a value object (`RunVariant` precedent).

---

## 15. Acceptance / sanity checks (Stage 1 and Rung A)

- Conservative-only invariants conserved to integrator tolerance.
- Design orbit parameters recovered from initial conditions.
- `dr_p/dv_a` measured numerically matches the closed-form value for the chosen apogee.
- Known injected apogee error nulled by midcourse correction to the nominal **interception crossing** (the position miss, not perigee — ADR 0003; perigee is the `dr_p/dv_a` acceptance cross-check, §8).
- A2 (finding, ADR 0006): the node-altitude sweep shows a **mid-descent** second burn adds ~0 along-track authority — authority is apogee-bound, so for this phase A1 accuracy dominates. The sweep table is the recorded result (not an "A2 beats A1" pass/fail); the disproven "cheap two-burn split" claim is retired with ADR 0005.
- A3 (finding, ADR 0007): the deterministic `Cd×Cr` perfect-model sweep is **controllable everywhere** at Δv « the propellant budget — coefficient dispersion is cheaply nulled by the apogee burn (transverse authority over the along-track crossing is large). The map is ~1D in `Cr` (drag negligible at apogee); the apparent authority boundary was corrector conditioning (the 200 km altitude event pins the radial direction → near-singular Jacobian), resolved by LM damping — closing ADR 0003 finding 3.
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
- **Crosslink note (ADR 0014):** because the miss is ultimately *relative* (the plate is on the target), C3's 1 m terminal endpoint rests on **differential GNSS against the target** — the target rocket flies a receiver and broadcasts over the same link class the coordinator nodes already use; common-mode GNSS errors cancel over the ≤100 km closing baseline (code-differential ~0.5–1 m, carrier-differential ~0.1–0.3 m). The target model stays fixed-point (A) in the sim; only the *nav-noise grade* σ_rel changes.
- **Crosslink note superseded (ADR 0015):** the PuffSat flies **no GNSS**. The relative-nav source is a **target-side tracker**: beacon-vs-star-background differential astrometry (self-referencing — platform attitude error cancels; absolute, killing common-mode aim bias) at **σ_θ ≤ 10 µrad required / 5 µrad target**, plus two-way ranging (ToA + time transfer for free). σ_rel becomes **range-dependent, σ_θ·R** — lateral knowledge improves as the PuffSat closes, but authority dies as R², setting the homing floor σ_miss ≈ 2σ_θ²v²/a_max. The target model still stays fixed-point (A) in the sim; the differential-GNSS grades survive only as the retained comparison point in C3b's sweep.
- **Train-relative framing (ADR 0016):** the fixed point (A) is re-read as the **train centroid**. The plane retargets pre-launch to the *measured* arrival corridor (ADR 0006's observable drift; ±~2 km ≈ 5 s of launch-window slip), absorbing the common-mode budget components (Cr prior bias, B1 erosion, density) and the *systematic* tails — the cannonball-prior-bias risk (ADR 0009) becomes a free retarget. Per-PuffSat independent scatter (~nav-dominated) still owns the 400 m catch radius; the plane repositions only ~0.1–0.5 m per 1 s arrival gap, so it can never chase individual PuffSats; MCC-2 re-scopes to independent tails only.

### 16.4 Sensor / navigation model
Drives UKF observability. Note a strong altitude dependence.
- **A. Accelerometer only (GOCE-style, 100 Hz).** Measures non-gravitational specific force (drag + SRP), so it directly anchors the `Cd·(A/m)` and `Cr·(A/m)` estimates. Does not observe absolute position; must be fused with a positioning source.
- **B. Position fix only.** GNSS (~1-10 Hz) is fully available near perigee but essentially unavailable at a 150,000 km apogee (above the constellation, weak spillover signals, poor geometry). Coordinator-node ranging/angles are available throughout (the paper's coordinator nodes track PuffSat positions) and are the realistic apogee anchor.
- **C. Fused suite (recommended realistic case).** Accelerometer at 100 Hz (drag observability) + coordinator-node relative measurements (anchoring throughout, including apogee) + GNSS when near Earth. Matches the paper's architecture.
- **Lean:** perfect state at Rung A to establish controllability, then a Gaussian placeholder and the fused suite C at Rung C to debug the filter. Flag that apogee navigation (coordinator-node relative measurements, GNSS unavailable) versus perigee navigation (GNSS available) is itself a sub-study: navigation quality is not uniform around the orbit.
- **Sharpened 2026-06-10 (ADR 0012):** the C1 suite is two-way range + two-way carrier Doppler from co-flying matched-`a` coordinator nodes (known-ephemeris beacons, perigee ≥200 km), with node geometry a *derived* requirement and "is Doppler load-bearing?" a *measured* output; GNSS is the perigee-regime position fix, reported with/without (~30–50 g unlocked spaceborne receiver — COCOM-limited COTS cannot track 10.8 km/s); the accelerometer enters at **C2** (it observes the coefficients, not position, and is blind during the coast). The apogee/perigee sub-study lands as C1's LinCov-swept apogee regime + reported descent hand-off. **Measured (ADR 0012 findings):** Doppler *is* load-bearing at the NEES-honest q — range-only degrades σ_Tvel ~2.7× (0.66 → 1.8 mm/s) and goes mildly NEES-optimistic, though it still meets the C0 requirement with ~13× margin.
- **Accelerometer reality (ADR 0013):** option A's "GOCE-style" instrument is kg-class with power/thermal demands a PuffSat cannot host; a flyable MEMS/tactical unit (bias instability ~1–10 µg) sits 2–3 orders above the SRP signal (9.1×10⁻⁸ m/s² ≈ 0.01 µg at nominal `Cr·(A/m)`), so **`Cr` is never accelerometer-observable**. Terminal drag (peak ~1700 µg, B3a) is trivially observable — and the accelerometer's flight role there is **direct drag-force sensing for C3** (it measures the product `Cd·(A/m)·ρv²` that terminal control needs, sidestepping the density/coefficient conflation), not coefficient estimation.
- **Terminal nav grade (ADR 0014):** the terminal aim floor is GNSS-bound, and the miss is *relative* (§16.3 crosslink note) — so the C3 sweep axis is the relative-nav grade **σ_rel ∈ {absolute GNSS 10 m, code-differential 1 m, carrier-differential 0.2 m}**, with the 1 m endpoint claim resting on code-differential (carrier tracking at ~10.8 km/s is the flagged stretch for the 0.2 m upside). MEMS drag feedforward (sub-meter residual after subtracting commanded thrust) and GNSS feedback never strain each other; coordinator nodes are the GNSS-free reported fallback only.
- **GNSS deleted (ADR 0015):** the PuffSat flies no GNSS receiver. Midcourse is the measured C1 beacon suite with **Doppler load-bearing** (range-only is dead at the 400 m radius); terminal is the §16.3 target-side astrometric tracker + two-way ranging, and the C3b sweep axis is re-keyed from σ_rel grades to **tracker grade σ_θ ∈ {2, 10, 50 µrad}** (σ_rel(R) = σ_θ·R; the constant 1 m point retained for comparison). Rationale recorded: COCOM (unlocked receivers × thousands of expendable units), retirement of the 10.8 km/s carrier-tracking stretch item, free two-way time transfer — *not* "latency/power" (GNSS is passive; crosslink latency lands in C4's τ budget). New bench-testable hardware items in exchange: a µrad-class tracker on a 1 Hz-hammered vehicle (mitigated by star-referenced self-calibration, ~1 ms exposures, gyro bridging) and detector saturation recovery through scheduled, narrowband-filtered impact flashes.

### 16.5 Maneuver model (impulsive vs finite burn)
- **A. Impulsive delta-v throughout.** Instantaneous velocity jumps, event-restart at each. Simplest; fine for brief midcourse corrections.
- **B. Finite burn for the terminal phase.** The 600 to 200 km terminal burn lasts ~5 minutes, is continuous and gimballed, and is coupled to drag during the descent. Modeling it as a single impulse loses the drag-rejection-during-descent physics and misstates propellant use.
- **Lean / settled (2026-06-03):** hybrid. Impulsive throughout Rung A (midcourse plus terminal aim, to establish controllability ground truth); finite gimballed burn introduced at Rung B, when actuator authority and drag-rejection-during-descent become the focus.

### 16.6 Midcourse correction schedule
- **A. Fixed pre-planned schedule.** Corrections at predetermined trajectory points (e.g. one near apogee, one at a fixed descent altitude), always executed. Simple, deterministic, Δv budgetable a priori; matches the classical statistical-midcourse-correction paradigm.
- **B. Threshold / MPC triggered.** A correction fires when predicted perigee error exceeds a threshold. More adaptive, fewer wasted burns, but more logic, harder to bound Δv a priori, and its trigger depends on navigation quality.
- **Decision (settled 2026-06-03): A for all of Rung A.** Two reasons beyond simplicity. (1) Rung A measures controllability against a *fixed, transparent* control law; an adaptive scheduler would confound the A3 coefficient sweep — you could not tell whether the authority boundary reflects the orbit or the scheduler's cleverness. (2) Threshold-triggering's value (fire only when needed) is entangled with estimation quality: at Rung A's perfect state the trigger is trivially perfect, which overstates how well it works once navigation error enters at Rung C, so evaluating it here gives a misleadingly optimistic read. Move to B only at Rung D, where MPC subsumes it naturally and estimation error is in the loop to test it honestly. Concrete Rung A schedule: correction 1 near apogee (cheapest leverage), correction 2 at a fixed mid-descent altitude (~800–1,000 km, above the terminal hand-off), plus the terminal aim; all always-execute. These are fixed trajectory points (timing set by geometry), which is what "fixed schedule" means here — not fire-on-realized-error. (A2 implements correction 1 + correction 2 only; the impulsive terminal aim is vestigial under Rung A's perfect state — it nulls to zero — and becomes a finite drag-rejecting burn at Rung B. A2's finding, ADR 0006 superseding ADR 0005: correction 2 at a **mid-descent** node adds ~0 *along-track* authority — that authority is apogee-bound — so at Rung A correction 2 does not earn its place for timing error; its role is perigee/radial trim and, at Rung C/D, observable-drift correction.)

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
- **Measured (C4, 2026-06-13; `puffsat_sim/latency.py`, pure, no JVM — dead-time is a loop-transfer effect, so the Orekit physics add nothing the stability question needs):** **Per-loop budget, split onto the loops the latency actually lives on.** The **terminal inner loop** carries only the fast local chain — tracker exposure 1 ms + ≤100 km crosslink 0.33 ms + onboard compute 1 ms + ms-class valve 5 ms = **7.3 ms**. The big comms round-trip + coordinator UKF/corrector + uplink (**70 ms**) rides the **midcourse outer loop, a discrete replan with no continuous phase loop to erode** (its latency only shifts *when* an impulsive correction is computed against an hours-scale timeline). **Phase check (the deliverable):** at the ~1 Hz drag-rejection bandwidth the 7.3 ms budget erodes ω_c·τ = **2.64°** of phase — single-digit, negligible against a 30–60° margin. **Noiseless τ-sweep on the C3b ZEM loop (confirmation):** the dead-time is a buffer holding the nav fix stale for τ; at the validated 1 Hz cadence the budget is 0.73 % of the control period (sub-tick → reproduces the zero-delay miss byte-for-byte), the homing miss stays flat through a full **1 s (1 tick ≈ 136× the budget)** and only degrades at 2 ticks. Cross-checked at the 100 Hz inner-sample rate, the 7.3 ms budget (0.7 ticks) is still flat (≈1× baseline) — **the absolute dead-time the loop tolerates (~1 s) is set by the closing dynamics, not the cadence**, and the budget sits ≫100× inside it, so the "1 control period" tolerance at 1 Hz is a cadence coincidence, not a fundamental limit. The sweep is read as **relative degradation, not absolute capture**: it is a pure double-integrator stand-in (running ~2× hot vs C3b's Orekit loop — no plate-frame coast, no feedforward, gravity-free geometry), so the absolute miss is not the figure of merit, but dead-time is a loop-timing effect independent of the floor. **C (delay + jitter) and the measurement-dropout knob are not built** — the sweep showed no fragility within the budget. The **combined entry-offset × tracker-noise stress at the dispersion tail** (C3b measured catch radius and nav floor on *separate* axes) is the **Rung-D full-MC question by construction.**

### 16.9 First deliverable
- **A. Rung 0 then Rung 1** (build / hello-world, then the unperturbed orbit with its sanity and leverage checks). Smallest artifact, de-risks the most (toolchain, frames, units, ephemeris, plus the `dr_p/dv_a` leverage check) before any force or control.
- **B. Jump straight into perturbations or control.** Riskier; debugs toolchain and physics (or control) at once.
- **Lean / settled (2026-06-03):** A. Start at Rung 0, proceed force-by-force through Stage 1 (2a→2d), then Rung A.

### 16.10 Controller at Rung D: MPC formulation & library

Rungs A–C deliberately hold a **fixed, transparent control law** — the deterministic A1/A3 corrector for the midcourse plus a dumb PID/LQR terminal feedback (ADR 0010 decision 5) — so a miss is attributable to the *knowledge quality* under study, not to controller cleverness (the §16.6 logic). **Rung D introduces MPC and measures its value against this C baseline.** This subsection records the *criteria and options*; the choice is **deferred to Rung D** — picking a formulation or library before the terminal dynamics (C3) and the uncertainty bounds (C1/C2/C4) exist would be premature commitment. It graduates to an ADR when chosen.

- **(a) When does MPC earn its place over the C baseline?** MPC (adaptive constrained replanning) must beat the corrector + terminal-feedback baseline on the *same* Rung-D Monte Carlo by a margin that justifies its complexity and compute: (i) lower interception miss and/or propellant at equal authority; (ii) it handles constraints the fixed law handles only implicitly — thrust ceiling/floor, the <2% propellant budget, keep-out/plume; (iii) it subsumes threshold-triggered scheduling and anomaly recovery (§16.6 B — the "fire only when needed" value that is only honestly testable once estimation error is in the loop). If the C baseline already meets the miss threshold under the Rung-D dispersion, MPC's value is marginal and it stays a second-pass refinement.
- **(b) When does *robustness* earn its place over nominal MPC? (the C→D bridge).** The C-rung characterizes the uncertainty MPC must survive — C1 a navigation covariance Σ, C2 a coefficient-knowledge tolerance, C3 a terminal drag-disturbance bound, C4 a latency budget. Robust MPC (tube / scenario-tree / min-max) explicitly accounts for these; nominal (certainty-equivalent) MPC does not. **Criterion: robustness earns its place only if nominal MPC violates the interception-miss threshold or a constraint under the characterized C-rung uncertainty** (measured on the Rung-D MC). Otherwise nominal MPC + the terminal feedback suffices and the robust machinery is unjustified complexity.
- **(c) Formulation options (surveyed, not chosen).** Nominal NMPC; **tube MPC** (a robust invariant tube around a nominal plan — clean if the dynamics linearize acceptably over the horizon); **scenario / multi-stage MPC** (a tree of sampled uncertainty realizations — a natural fit for the C-rung's sampled Σ / coefficient draws); **min-max MPC** (worst-case — usually too conservative for a propellant-tight mission). The terminal phase is nonconvex (drag + finite burn), pointing at sequential-convex (SCP/GuSTO) or direct NMPC; the midcourse may linearize to a convex QP.
- **(d) Library options (surveyed, not chosen).** Mirrors ADR 0003's "own the simple solver, defer the library": (i) **`do-mpc`** (CasADi-based Python; supports multi-stage *scenario* robust MPC out of the box — the fastest path to the Rung-D study); (ii) **CasADi** directly (DIY NLP, maximum control); (iii) **`acados`** (C-generated SQP/RTI, embedded-grade — a *flight-code* concern, not the sim's, named so it isn't conflated); (iv) **OSQP / cvxpy** (if a leg linearizes to a convex QP); (v) **SCP / GuSTO** for the nonconvex terminal. `python-control` is rejected for the same reason as ADR 0003 (feedback-compensator design, not trajectory optimization).
- **Lean / deferred to Rung D:** prototype in Python (`do-mpc` or CasADi) for the study; pick the robust formulation *empirically* by criterion (b) — start nominal, add tube/scenario only where the Rung-D dispersion shows a threshold/constraint violation. Neural warm-start / RL anomaly recovery (CLAUDE.md second pass) layer on top of whichever classical MPC is chosen. Cross-refs: ADR 0003 decision 5, ADR 0010 decision 5, §16.6.

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
