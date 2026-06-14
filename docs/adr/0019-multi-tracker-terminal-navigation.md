# Multi-tracker terminal navigation: a target-side detector array + a co-flying close tracker recover capture-grade nav from 10 µrad detectors

**Status:** accepted

## Context

D1.1 (the first runnable Rung-D result) surfaced the binding feasibility effect: the **combined
entry × tracker-noise stress** — which C3b had measured *separately* (entries noiseless → null to
mm; tracker noise at zero entry → 1.07 m) and deferred to Rung D — fails the σ ≤ 1.65 m capture
criterion at the nominal **10 µrad** single-target-tracker grade (scatter σ 5.77 m, **31 %
capture**), and needs **~3 µrad**. The mechanism is the **early large-R noise**: the hand-off→target
range is 2603 km, so σ_θ·R = **26 m** at 10 µrad early in the descent, and a real ~150 m entry
offset forces the loop to fire hard *early*, through that noise, defeating the significance gate
that protects the zero-entry case. A tighter grade shrinks the early noise and recovers (3.2 µrad →
capture-grade; the grade cliff is sharp).

A single conservative detector *already* achieves 3.2 µrad (the σ_θ budget gate, ADR 0018), so D1
closes on the C baseline — but the dominant term is a **3 µrad bench-calibratable focal-plane
distortion floor**, the one new bench-testable risk the gate flagged. A grilling session
(2026-06-14) asked: if that 3 µrad floor proves optimistic on a hammered, vibrating vehicle, can we
recover capture-grade relative nav from cruder, redundant **10 µrad** detectors instead of betting
on one precise calibration — and can we attack the early error at its source? The answer is a
**multi-tracker terminal-nav architecture**, resolved lever by lever against the existing docs.

## Decision

1. **Target-side detector array.** Put N independent 10 µrad detectors on the target (an X-pattern
   across a ~30 m airframe). Fusing them buys **σ_θ/√N down to a common-mode floor** — *iff* each is
   a genuinely independent measurement: separate optics, **separately bench-calibrated distortion
   maps**, each doing its own beacon-vs-its-own-star-field astrometry (ADR 0015's "absolute
   star-referenced measurement" is what keeps the attitude term independent). 5 detectors → 3.2/√5
   ≈ **1.4 µrad**. The common-mode floor √N cannot cross is **correlated distortion** (shared cal /
   manufacturing batch) + **beacon-shape asymmetry** (all detectors image the *same* beacon), so the
   design rules are *separate calibrations* and *a clean symmetric point beacon*. The **spatial
   spread is for coverage + redundancy + common-boresight rejection, not precision** — the √N is
   statistical and separation-agnostic (5 co-located independent detectors give the same √5), so the
   G650-scale baseline does **not** tighten σ_θ. A **shared clock is fine** (distortion is a spatial
   error, independent of timing; the clock only couples into the already-met ToA and the small
   smear term), even helpful for ranging. **Ranging is a red herring** for the binding lateral axis:
   a 30 m baseline at 100 km+ gives ~33 mrad from range-differencing — three orders worse than a
   single 10 µrad angle — so the angles do all the lateral work.

2. **Co-flying launch-rocket close tracker.** Reuse the already-on-orbit launch rocket (it deployed
   the PuffSats) as a **close** tracker — a small apogee maneuver (raise perigee, lower apogee) keeps
   it near the descending train. At ~500 km vs the target's 2603 km, σ_θ·R is **5× smaller**, which
   attacks the early error at its *source* (a bigger lever than the array's √N). It tracks the
   **train centroid** — it shrinks the *common* early error; per-unit scatter stays owned by the
   target array. **Load-bearing condition:** the measurement is PuffSat-relative-to-*rocket*, but the
   miss is relative to the *target*, so the rocket→target vector must be known **independently of
   long-baseline angle tracking** — otherwise that ~2100 km baseline's lateral error (≈ 21 m at
   10 µrad) eats the whole close-range advantage. This is satisfied because the terminal phase is
   **low altitude** (200–800 km, deep inside the GPS constellation): **unlocked spaceborne GNSS**
   pins the rocket→target vector to ~metres (the 150,000 km apogee is irrelevant *here*). Note
   inter-rocket *ranging alone* fixes the range but **not** the lateral — independent absolute
   position is required, not just two-way ranging. The **gating unknown is phasing**: can the rocket
   stay ~≤ 500 km from the descending centroid *and* low enough for GPS *and* high enough not to
   decay through the 800→200 km window, given the maneuver + ~32 h drift? That is a simulation
   question (Stage 2 below), and it gates this lever — not the architecture.

3. **High-altitude nav infrastructure (optional, midcourse).** A permanent nav constellation in the
   apogee regime (~150,000 km, above GPS) would pin the **apogee state** better than an ad-hoc
   co-flying node, shrinking the **entry offset** (the forcing) rather than the terminal noise. It
   generalizes the C1 coordinator nodes into permanent, well-characterized GDOP. It is the secondary
   lever (terminal noise dominates the D1.1 miss) and is **recorded as an option, not built now** —
   the entry offset is already a swept D1.1 spec field, so its effect is explored by a sweep.

4. **What this re-keys downstream (bounded).** The terminal-nav grade of ADR 0015 is re-read as a
   **per-detector** σ_θ; the **system** grade the ZEM loop sees is the **fused effective σ_θ**. The
   C3b noise model (`runs/guidance`) is re-keyed from a single-tracker σ_θ·R to the fused grade, and
   D1.1 is re-run at it. The σ_θ budget gate finding (3.2 µrad *per detector*) **stands** — fusion
   sits on top of it. **A/B and the C0–C2a requirement results are untouched.**

## Considered options

- **Just build the single 3.2 µrad tracker (the σ_θ-gate path).** Retained as the *baseline* — it
  already closes D1. The array is the **hedge** against the 3 µrad distortion floor proving
  optimistic; it is not a correction of a wrong baseline.
- **Credit the X-pattern baseline with precision (stereo / triangulation).** Rejected: the √N is
  statistical (separation-agnostic), and stereo range from a 30 m baseline at 100 km+ is negligible.
- **Use ranging / multilateration for the lateral.** Rejected: short baseline → ~33 mrad, useless
  for the binding axis; ranging only sharpens range/ToA, already met by two orders.
- **Establish the rocket→target relative vector by inter-rocket ranging.** Rejected for the lateral:
  ranging fixes range, not the lateral of the long baseline; independent GNSS/ground position is
  required.
- **Bundle MPC (D2) to recover capture.** Rejected (ADR 0018 §16.6/§16.10): the verdict must be
  reached on the dumb C law; MPC earns its place only on a measured D1 violation, which the fused
  grade removes.

## Consequences

- **Pure-side:** a new `tracker_fusion.py` gate (fused effective σ_θ from a set of `Tracker`s,
  building on `tracker_budget.py`), unit-testable without a JVM.
- **JVM-side:** a `runs/coflyer.py` phasing-feasibility run (Stage 2, the Lever-2 gate); the C3b
  noise re-key in `runs/guidance` and the D1.1 re-run (`runs/train`).
- **Docs:** ADR 0015 terminal grade re-read as per-detector (system = fused); ADR 0011/0012
  coordinator nodes generalize to a reusable co-flying rocket + an optional permanent infra; design
  §13/§16.4 gain the multi-tracker revision; CONTEXT gains *Tracker array*, *Co-flying tracker*,
  *Effective σ_θ*.
- **The conditional verdict is strengthened.** D1's "feasible conditional on ~3 µrad" now has the
  condition reachable **two independent ways** — one tighter detector (σ_θ gate) *or* fused cruder
  10 µrad detectors (this ADR) — so the load-bearing terminal-nav assumption is robust to the
  distortion floor, not a single point of failure.

## Implementation findings — multi-tracker fusion gate (Stage 1, 2026-06-13)

The pure gate (`puffsat_sim/tracker_fusion.py`, no JVM, building on `tracker_budget.py`), TDD.
Both levers were quantified against the D1.1 capture-grade (3.2 µrad effective).

- **The target-side array buys √N down to the smear common-mode floor.** A `Tracker`'s σ_θ splits
  into the **independent** part (distortion ⊕ gyro ⊕ photon — separately bench-calibrated, so it
  averages as σ_indep/√N) and the **common** part (`SMEAR_COMMON_SIGMA_RAD`, the same beacon imaged
  by every detector, which √N cannot cross). Five 10 µrad detectors fuse to **1.62 µrad** at the
  2603 km target range — **2.0× inside** the requirement — with **no phasing/baseline dependency**
  (the √N is statistical, separation-agnostic; the X-pattern spread is for coverage and
  common-boresight rejection, not precision). This confirms decision 1's "ranging is a red herring"
  framing: the angles do all the lateral work.
- **The close co-flyer attacks the early error at its source.** Adding the rocket at the
  `COFLYER_RANGE_M` (500 km) design range — 5× closer than the target — drops the inverse-variance
  fused grade to **0.76 µrad** (**4.2× inside**), a bigger lever than the array's √N because
  σ_θ·R scales with range. The credit is real only if the rocket→target vector is pinned
  independently of the long baseline (`COFLYER_RELGEOM_SIGMA_M`, the GNSS-pinned floor); whether the
  rocket can *hold* that 500 km range and stay in the GNSS volume is the Stage 2 gate.

## Implementation findings — co-flyer phasing gate (Stage 2, 2026-06-14)

The JVM run (`runs/coflyer.py`) for Lever 2's gating unknown — phasing — feeding the pure
`phasing_verdict`. The rocket orbit is **constructed directly** at the constant-semi-major-axis
elements (perigee +100 km / apogee −100 km, co-phased in mean anomaly), since the maneuver that
realizes them is a separate propulsion detail and the phasing question is purely geometric. Physics
is J2 (a ~32 h coast geometry, consistent with the truth-validation coast), flown alongside the
nominal descent and sampled across the 800→200 km terminal window.

- **The constant-a maneuver holds the period exactly, so there is no secular drift.** +perigee /
  −apogee by the same amount preserves `perigee_alt + apogee_alt` → `a` → period (an integration-test
  invariant), and the rocket keeps the train's mean anomaly. Over the window the rocket↔centroid
  separation peaks at only **125 km** — **well inside the 500 km angle-useful design range** the
  fusion credits — so the σ_θ·R advantage Stage 1 banked is real, not assumed.
- **The rocket stays in the GNSS volume and aloft through interception.** Its altitude over the
  window is **295–879 km**, far below the **20 200 km GPS ceiling** (so an unlocked spaceborne
  receiver pins the rocket→target vector) and above the 200 km crossing (raised perigee → it does
  not decay through the window). Verdict: **PHASING-FEASIBLE** — the Lever-2 close-tracker credit
  holds, and the conditional-verdict hedge of decision 2 stands.

## Implementation findings — C3b noise re-key + D1.1 re-run (2026-06-14)

Decision 4 made runnable. The re-key is **one pure function** — `tracker_fusion.fused_tracker_grade`
collapses a multi-tracker architecture into the C3b loop's `TrackerGrade` (effective σ_θ at the
target design range + the unchanged along-LOS ranging σ). The C3b noise model (`guidance.NavNoiseProcess`)
already consumed a scalar σ_θ·R, so fusion needed **no new noise code**; `runs/train.run_train_dispersion`
gains a `trackers=` path that flies the train at the fused grade, and `fused_train_rerun_report`
re-runs D1.1 across the architectures over one shared hand-off context.

- **The re-run closes the conditional verdict.** Flown at the **legacy single-tracker 10 µrad
  ceiling**, D1.1 fails exactly as before (scatter σ **5.21 m** vs ≤1.65 m, **50 % capture**). The
  fused architectures recover capture-grade with margin that tracks the effective grade:
  single target detector (σ_θ gate, 3.18 µrad) σ **1.17 m** / 100 %; **target 5-array** (Lever 1,
  1.62 µrad) σ **0.58 m** / 100 %; **target 5-array + co-flyer** (Levers 1+2, 0.76 µrad) σ **0.21 m**
  / 100 %. Scatter σ falls ~linearly with the effective grade — the σ_θ·R early-noise mechanism D1.1
  identified, now shrunk at its source. So D1's "feasible conditional on ~3 µrad" is **delivered**:
  cruder detectors fused reach the grade, two independent ways.
- **The cost stays inside budget.** Worst-unit mission Δv lands 0.77–1.03 % @ Isp 50 across all
  architectures (well under 2 %), and perigee holds ~65 km (deorbit-good). The co-flyer's tighter
  scatter does not buy lower Δv (the funnel still pays for each unit's entry); its value is capture
  margin against the distortion-floor risk, exactly the decision-2 hedge.
- **What did *not* change (decision 4's "bounded").** The C3b sweep's σ_θ axis is now re-read as the
  *system* (fused) grade — its capture characterization (3.2 µrad capture-grade, 10 µrad fails)
  **stands** as the reference curve fusion sits on. The σ_θ budget gate (3.2 µrad per detector) and
  the A/B + C0–C2a results are untouched. **ADR 0019 is complete** (all three levers: array + co-flyer
  built, high-altitude infra recorded as an option).
