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
| **1 — Simulated & closed out** | Closed-loop terminal control delivers a PuffSat onto a **5 m pusher plate at ≥99% per-unit capture, directly counted** (0 escapes / 320 units over 8 independent trains → one-sided 95% LB **99.07%**) at the committed fused grade under the **implemented** dispersions (LEO scenario). That figure credits the ±2 km centroid retarget the design specifies; the flown loop does not implement it, and *as flown* the same batch gives **99.06% (3 escapes / 320, LB 97.6%)** — see the numbers below, which name both conventions. The target rocket is scored in a relative plate frame against a placeholder trajectory but is still *flown to* as a fixed aim point (see caveats). | Closed-loop Monte Carlo; Rung D / D1 closed ([ADR 0018](docs/adr/0018-rung-d-decomposition.md)). |
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
  requirement at a minimal 3-member geometry — 4 members when the receiver clock bias is
  solved from the signal rather than held). The corrector consumes this state.
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
  achievable 3.2 µrad grade is **marginal**, and this is the **directly-flown result**:
  brute-force Monte Carlo (N = 500) gives **P(capture) ≈ 99.2%, arrival σ ≈ 1.51 m**. So the
  **committed architecture is the fused tracker array**: fusion drives the arrival scatter down
  (5-array → σ ≈ 0.9 m; + co-flyer → σ ≈ 0.3 m).
- **What the counting establishes at the fused grade — and the entry convention it depends on**
  (measured 2026-07-22, 8 independent trains × 40 units = **320** at the 5-array 1.62 µrad grade;
  reproducible via `runs/train.pooled_train_capture_report`). This **supersedes an earlier "zero
  escapes in N = 192"**, which was a **train-count artefact**: a train shares one common-mode draw,
  so pooling units *within* a train never samples the shared entry leg. The result splits cleanly:

  | Entry convention | Result | One-sided 95% LB |
  |---|---|---|
  | **As flown** — funnel nulls shared ⊕ per-unit entry | 3 escapes / 320 → **99.06%**, core σ 0.53 m | 97.59% — does *not* establish ≥99% |
  | **Retarget credited** — funnel nulls per-unit entry only | 0 escapes / 320 → **100%**, core σ 0.60 m | **99.07% — establishes ≥99%** |

  The per-train **shared** entry leg (149 m) is *specified* to be absorbed by the ±2 km centroid
  retarget ([ADR 0016](docs/adr/0016-train-relative-requirement-framing.md)), leaving the funnel only
  the **per-unit** leg (141 m). The flown loop **does not implement that retarget** — it makes the
  funnel null both legs — so the as-flown row is the *conservative* one, and the gap is
  implementation-vs-spec, not a hidden assumption.
- **All escapes are entry-cliff events, not noise events.** The lowest escaping hand-off entry was
  **491 m** against the ~475 m analytic catch radius; the retarget-credited batch never exceeded
  344 m and saw none. Arrivals are therefore **bimodal** — a captured core at σ ≈ 0.5 m plus
  tens-of-metres outliers — so a pooled σ over the whole set is a mixture statistic and **must not**
  be read against the 1.65 m criterion.
- The often-quoted **≈99.999%** is the deep-tail figure on the *retarget-credited* convention, and
  is additionally an **importance-sampling + catch-radius-analytic extrapolation** — the committed
  `tail_capture` runner defaults to the single-detector grade, so it is recorded from analysis, not
  reproduced end-to-end. The directly-counted claims are the single-detector **99.2%**, the as-flown
  fused **99.06%**, and the retarget-credited fused **≥99% at 95% confidence**.
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
the architecture. It is therefore carried as a **deferred extension**, now *sized* rather than
hand-waved by a pure sizing module ([`centering_budget.py`](puffsat_sim/centering_budget.py),
[ADR 0022](docs/adr/0022-surveyor-anchored-centering-tier2-sizing.md)), matching the
`tracker_budget.py` / `apogee_nav.py` / `distortion_field.py` precedent.

The plate is the **RSS of two legs**, `plate = 3.03·√(σ_hoop² + (σ_θ·v/f)²)`. The km-class
intra-train link is kept **distortion-limited, not photon-limited**, by a **Q-switched,
coarse-pointed beacon** (bright ns pulses — ~100 kW peak, few-hundred-mW average — read in a
matched gate; photon term ~17× under the 3 µrad floor). The committed claim: **10 cm robust**
(σ_hoop ~1 cm rendezvous-lidar class ⊕ ~3 µrad calibrated camera → **5.8 cm nominal, 87× off the
5 m**; 10 cm tolerates σ_hoop ≤ 2.9 cm); **5 cm contingent on tightening *both* legs** — a
mm-class hoop **and** a smaller scatter (4 Hz train, diverse cameras, or a lower distortion
floor), since the scatter leg alone sets a ~4.9 cm floor at the nominal point. It remains
**argued/sized, not simulated**: nothing new is simulatable for it with the current architecture,
so the simulated frontier stays the 5 m plate and the right next rigor for 10 cm is a **bench
test**.

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
- **The moving target is modeled by a relative plate frame with one placeholder trajectory**
  ([ADR 0023](docs/adr/0023-moving-target-terminal-frame.md)). The plate-frame miss is now
  measured ⊥ the **closing** velocity `v_rel = v − v_target` against a constant-velocity target;
  the historical fixed-point frame is the `v_target = 0` case. For the nominal **mirror-ascending**
  target the closing speed is only ~3 km/s — interception at 200 km sits ~8° from horizontal, so a
  co-moving ascending rocket closes gently — and capture is *unaffected* (the fixed-point result is
  mildly conservative, not optimistic, for this model). Scaling the mirror velocity sweeps a family
  of co-directional targets whose closing speed rises monotonically as the target *slows*, and the
  fixed frame is exactly its `v_target = 0` endpoint — so the committed numbers bound that whole
  family from the conservative side. **They do not bound geometries outside it:** a head-on or
  high-cross-track launch would raise `|v_rel|` past the PuffSat speed and could erode margin, so
  the real launch trajectory is a mission-design input still to be supplied. **Still deferred, and
  this is the substantive gap:** the loop is *scored* relative but still *flown* to a fixed aim
  point — no target lead is modeled — so ADR 0023 is a relative-frame sensitivity study, not an
  end-to-end powered-target simulation. Also deferred: a dispersed/powered target (its *position*
  dispersion is largely nulled by beacon homing; the residual is target *velocity* uncertainty),
  and wiring the moving frame into the D1.1 train and tail-capture runs (which today use the
  fixed frame).
- **Terminal velocity knowledge is idealized.** The onboard ZEM state receives **exact truth
  velocity**; only position carries σ_θ·R measurement noise, and the hand-off residual is
  displaced in position with nominal velocity. Terminal miss at ~10.8 km/s depends on both
  position *and* velocity error, so the absence of any range-rate / clock / velocity-estimate
  error is a simplification that favors the result.
- **The dominant σ_θ floor is flown as zero-mean noise but argued to be a persistent bias.**
  The guidance loop injects the calibration floor as a 10 s zero-mean Gauss–Markov process, which
  the √N array averaging then reduces; but the same floor is argued elsewhere to be a quasi-static
  focal-plane **distortion bias** (the reason there is "no hardware walk-back"). The fusion benefit
  (5-array 1.62 µrad; +co-flyer 0.76 µrad) is therefore an **optimistic bound** whose realization
  depends on the bench-measured distortion correlation length
  ([`distortion_field.py`](puffsat_sim/distortion_field.py)); the co-flyer, as a physically
  independent platform, is the hedge that survives fully-correlated distortion.
- **Lower-order navigation-model simplifications** sit under the C1 apogee-nav numbers. They do
  not threaten the *knowledge-limited* conclusion — that is a regime finding, robust to the exact
  σ's — but they should be tightened before those σ's are quoted as hardware requirements:
  (i) the C1 LinCov places beacon nodes at rigid RTN offsets carrying the PuffSat's exact inertial
  velocity, so the observability geometry is reconstructed each epoch rather than propagated as
  realizable orbits; (ii) the filter process noise is a white-acceleration PSD but is parameterized
  as a physical acceleration, so `q`-vs-tidal-acceleration comparisons are only order-of-magnitude;
  (iii) NEES consistency is averaged over temporally correlated epochs of a single trajectory, so
  it is an indicative check, not a full Monte-Carlo consistency test; (iv) the passive one-way
  apogee-Doppler result solves three velocity components with the receiver clock held, so its
  thermal precision omits TCXO drift and the **match-not-beat** margin there is the softest of the
  sizing numbers. These are lower-order relative to the terminal-nav story and are carried as
  refinements for the expert review.
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
