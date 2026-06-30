# Conclusion — Is PuffSat Pulsed Propulsion Navigationally Feasible?

This document summarizes what the closed-loop simulation in this repository set out
to answer, the architecture it converged on, and the verdict it reached. It is the
top-level read for the companion speculative white paper
[*Aim Is All You Need*](https://doi.org/10.5281/zenodo.16741183); the supporting
detail lives in [`puffsat_control_sim_design.md`](puffsat_control_sim_design.md),
the domain glossary [`CONTEXT.md`](CONTEXT.md), and the architecture decision records
in [`docs/adr/`](docs/adr/).

## The question

Can a single-use PuffSat — deployed near apogee (~150 000 km altitude) on a highly
eccentric Earth orbit — navigate down to an interception near perigee (200 km
altitude) with enough precision to **hit a pusher plate on a target rocket**, under
realistic orbital perturbations, while spending less than ~2% of its mass on
propellant? A sequence of PuffSats then transfers momentum to the rocket by
sequential plate impacts; any PuffSat that misses cleanly flies by and burns up on
reentry (debris disposal by design, paper §9).

## Verdict (three confidence tiers)

The repository establishes feasibility at three deliberately distinct levels of
proof. Keeping them separate is the whole point — only the first is a flown result.

| Tier | Claim | Status |
|---|---|---|
| **1 — Simulated & closed out** | A PuffSat captures a **5 m pusher plate at ≥99% per-unit confidence** under full dispersions (LEO scenario). | Closed-loop Monte Carlo; Rung D / D1 closed ([ADR 0018](docs/adr/0018-rung-d-decomposition.md)). |
| **2 — Argued / sized, not closed-loop** | A surveyor-anchored, camera-rigid extension drives the per-unit miss to a **metrology-limited ~10 cm plate**. | Analytic architecture argument; deferred follow-on, no ADR ([CONTEXT.md](CONTEXT.md): *Surveyor-anchored centering*). |
| **3 — Architectural extension sketch** | The same design pattern extends to a **near-Sun / Parker-class** trajectory, swapping star references for artificial beacons. | Paper-side, reversible, **two open sizing numbers, not simulated in this repo** ([CONTEXT.md](CONTEXT.md): *Near-Sun optical nav*). |

The headline result of the simulation is **Tier 1**. Tiers 2 and 3 are honest
extrapolations carried for the speculative paper, labeled as such.

## The architecture (Tier 1 — the simulated LEO design)

### Reference orbit and the apogee decision

The pivotal choice (design doc §3) is **apogee at ~150 000 km** rather than near the
Hill radius. Pulling apogee in cuts perigee-altitude sensitivity from ~250 to
~30 km per m/s, drops solar-tidal perturbation from ~49% to ~0.1% of local gravity
(removing chaos), and shrinks SRP dispersion by an order of magnitude through a
shorter coast. The reference orbit: e ≈ 0.921, a ≈ 81 403 km, period ≈ 2.68 days,
interception at 200 km during descent, perigee 50 km (below the Kármán line, for
intentional burn-up).

### Closed-loop structure

Three decoupled pieces in a loop: an Orekit high-fidelity **truth model**, a
**UKF state estimator** (position, velocity, and the lumped drag/SRP coefficients
`Cd·(A/m)`, `Cr·(A/m)`), and a **controller** — a differential corrector for the
discrete midcourse burns plus a continuous terminal burn. Feasibility is a property
of this "C baseline" controller; full MPC (D2) was scoped but **never triggered**,
because the baseline meets the requirement with margin ([ADR 0018](docs/adr/0018-rung-d-decomposition.md)).

### Navigation, by regime

- **Coast / apogee state** is pinned by the **Apogee nav constellation** — a permanent
  ~150 000 km **Ka-band authenticated-broadcast** infrastructure that the passive PuffSat
  receives one-way ([ADR 0020](docs/adr/0020-apogee-nav-constellation-signal-architecture.md)).
  This is **not GNSS** — GPS reaches only ~20 200 km, and the apogee regime and geometry
  differ. It is sized "match-not-beat" the C1 nav grade (transverse-velocity σ well inside
  requirement at a minimal 3-member geometry). The corrector consumes this state.
- **Terminal homing** is angular: the PuffSat homes on an **active optical beacon**
  (LED, or a **Q-switched 1064 nm laser** for strobed peak power) **centroided against a
  reference-star field via differential astrometry** — the "Gaia trick," which cancels
  common-mode focal-plane distortion. The figure of merit is the tracker angular grade
  **σ_θ** ([ADR 0015](docs/adr/0015-plate-capture-criterion-gnss-free-navigation.md)). The
  σ_θ budget is **calibration-limited, not photon- or vibration-limited** (~3.2 µrad
  achievable; [ADR 0018](docs/adr/0018-rung-d-decomposition.md)).
- **Real GNSS appears only at low altitude**, in the terminal phase, as the co-flyer's
  ~2 m relative-geometry anchor — never at apogee.

### Terminal sensor: array (load-bearing) + co-flyer (hedge)

The **target-side tracker array** (5 independent, separately bench-calibrated detectors)
is the load-bearing terminal sensor; fusion buys √N down to a common-mode floor, reaching
~1.62 µrad ([ADR 0019](docs/adr/0019-multi-tracker-terminal-navigation.md)). The
**co-flying tracker** (the reused launch rocket, a closer vantage) is a **σ_θ-independent
redundant hedge** that adds margin (to ~0.76 µrad) — it is *not* required for the 5 m
capture verdict, and it carries the relative-sensing role in the Tier-2/3 extensions.

### Terminal guidance and drag

Guidance is a zero-effort-miss (ZEM) law against a thrust-limited **catch-radius funnel**
of ~475 m (½·a_max·t² with the 400 mN actuator;
[ADR 0014](docs/adr/0014-c3-terminal-guidance-architecture.md)). A small out-of-plane
**MCC-2** trim nulls residual lateral entry error the funnel cannot. Terminal **drag is
feedforward-solved**, not feedback-rejected — it is along-track, exponentially late, and
cancels to ~mm; the binding terminal error is **cross-track navigation knowledge, not drag**
([ADR 0021](docs/adr/0021-terminal-drag-is-feedforward-solved-binding-error-is-nav-knowledge.md)).

### Train-relative framing

PuffSats arrive as a **train** (~1 s spacing), never station-keeping. Common-mode train
drift is absorbed by the target plane's **centroid retarget** (±~2 km launch-time aim);
only **independent per-unit scatter** must be caught by the terminal burn
([ADR 0016](docs/adr/0016-train-relative-requirement-framing.md)). Capture is judged
per-unit *about the centroid*.

## Tier 1 result — the numbers

The success criterion ([ADR 0015](docs/adr/0015-plate-capture-criterion-gnss-free-navigation.md)):
**5 m plate radius**, **≥99% per-unit capture ↔ σ_lateral ≤ 1.65 m** (2-D Rayleigh,
R/σ = 3.03), plus **ToA ≤ ~10 ms** at closest approach.

The closed-loop Monte Carlo ([ADR 0018](docs/adr/0018-rung-d-decomposition.md)) found,
flying the real terminal loop through the catch-radius cliff, the significance gate, and
σ_θ·R measurement noise:

- **Combined entry × noise stress sets the binding terminal-nav requirement at ~3.2 µrad**
  — tighter than the 10 µrad ceiling each effect implied alone. A **single detector** at the
  achievable 3.2 µrad grade is **marginal** (P(capture) ≈ 99.2%, σ ≈ 1.5 m), so the
  **committed architecture is the fused tracker array**: 5-array → σ ≈ 0.9 m / 100% capture;
  + co-flyer → σ ≈ 0.3 m / 100%. The importance-sampled tail at the fused grade reaches
  **P(capture) ≈ 99.999%**.
- **Propellant** stays under the ~2% claim across the stack (~0.9% typical, ~1.2% worst case
  at Isp 50 s).
- **ToA** is ≤ ~0.7 ms — two orders inside the 10 ms requirement.
- **Perigee** lands ~65 km, a diagnostic where *low is good* (debris disposal), not a target.
- **Entry margin:** the per-unit hand-off entry budget (~224 m on the conservative (T, N)
  convention) sits ~2.1× inside the 475 m catch radius; the truly binding cross-track entry is
  far smaller, so the margin is conservative.

Both upstream entry legs — the C1 nav residual ([ADR 0012](docs/adr/0012-c1-sensing-architecture-lincov-feasibility.md))
and the C2 coefficient-prior residual ([ADR 0013](docs/adr/0013-c2-coefficient-knowledge-prior-bound.md)) —
were validated end-to-end against the real corrector. **Verdict: D1 feasible on the dumb C
baseline, conditional on the fused terminal-nav grade. MPC is not needed.**

## Tier 2 — the deferred ~10 cm tightening (argued, not simulated)

The 5 m plate can in principle be shrunk toward **~10 cm** by driving the *per-unit* arrival
miss to centimeters with a **Surveyor-anchored centering** scheme ([CONTEXT.md](CONTEXT.md)):
a sacrificial "surveyor" PuffSat read by an **independent** lidar/microwave hoop pins the
swarm's quasi-static optical-distortion bias to the plate, while strobed known-pattern LED
beacons make the swarm camera-rigid (bearing → binding cross-track).

This is **knowledge/metrology-limited, not control-limited** — cm trims sit deep inside the
475 m funnel, so it does not resurrect MPC. Its binding numbers (the hoop precision σ_hoop and
the camera's calibrated distortion floor) are **bench/hardware characterizations a simulation
cannot produce**; a closed-loop Monte Carlo here would only sharpen a distribution, not falsify
the architecture. It is therefore carried as a **deferred extension**. The committed claim:
**10 cm robust** (~50× off the 5 m baseline); **5 cm contingent** on σ_hoop ≤ 1 cm and a
calibrated camera. (If a notch more rigor is ever wanted, the right move is a *pure sizing
module* outputting plate-size-vs-metrology sensitivity — matching the
`tracker_budget.py` / `apogee_nav.py` / `distortion_field.py` precedent — not an MC.)

## Tier 3 — near-Sun / Parker extension (sketch, open numbers)

The same architectural pattern — an off-board nav infrastructure anchoring the hard end of a
highly eccentric orbit — extends to a **near-Sun / Parker-class** trajectory. There the corona
blinds the star tracker, so the natural star field is replaced by **one-way artificial-star
beacons at ~1 AU** (broadcasting Earth-pinned positions, bearing-grade), with the millimeter
relative lever carried by **two-way pulsed (Q-switched 1064 nm) laser ranging** between the
converging projectiles and the coordinator node — self-clocking, scintillation-robust, reusing
the LEO tracker wavelength lineage ([CONTEXT.md](CONTEXT.md): *Near-Sun optical nav*).

This layer is **the most speculative**: paper-side, reversible, **no ADR**, with **two still-open
sizing numbers** (in-band thermal-IR rate off the heat shield; coronal angle-of-arrival jitter on
the radial path — both expected to close) and **no simulation in this repo**. It is offered as a
plausibility argument, not a verified result.

## Scope and honest caveats

- The feasibility verdict is **conditional** — feasible *given* the nav and actuator specs the
  Monte Carlo takes as inputs. The σ_θ tracker budget, apogee-constellation GDOP, and torque
  margin convert those inputs toward derived hardware requirements
  ([ADR 0018](docs/adr/0018-rung-d-decomposition.md)).
- The drag/SRP shape model is a **conservative cannonball** placeholder; the optimistic
  attitude-dependent cylinder (Rung E) can only improve a passing number
  ([ADR 0009](docs/adr/0009-shape-fidelity-cannonball-placeholder.md)).
- An **independent external cross-check** (GMAT, Rung F) is deferred — the in-house two-tier
  truth validation (machine-precision conservation + an independent RK4 Cowell cross-check)
  suffices for this speculative verdict; GMAT is reviewer-defensibility insurance.

## Architecture decision records cited

[0003](docs/adr/0003-a1-differential-corrector.md) differential corrector ·
[0004](docs/adr/0004-rung-b-actuator-model.md) finite-burn actuator ·
[0009](docs/adr/0009-shape-fidelity-cannonball-placeholder.md) shape-fidelity ·
[0011](docs/adr/0011-c0-navigation-requirement-sensitivity.md) nav sensitivity ·
[0012](docs/adr/0012-c1-sensing-architecture-lincov-feasibility.md) C1 sensing ·
[0013](docs/adr/0013-c2-coefficient-knowledge-prior-bound.md) coefficient prior ·
[0014](docs/adr/0014-c3-terminal-guidance-architecture.md) terminal guidance ·
[0015](docs/adr/0015-plate-capture-criterion-gnss-free-navigation.md) plate-capture criterion ·
[0016](docs/adr/0016-train-relative-requirement-framing.md) train-relative framing ·
[0018](docs/adr/0018-rung-d-decomposition.md) Rung D / D1 closeout ·
[0019](docs/adr/0019-multi-tracker-terminal-navigation.md) multi-tracker fusion ·
[0020](docs/adr/0020-apogee-nav-constellation-signal-architecture.md) apogee nav constellation ·
[0021](docs/adr/0021-terminal-drag-is-feedforward-solved-binding-error-is-nav-knowledge.md) terminal drag is feedforward-solved.
</content>
</invoke>
