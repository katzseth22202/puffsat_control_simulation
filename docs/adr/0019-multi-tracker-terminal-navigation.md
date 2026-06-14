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
