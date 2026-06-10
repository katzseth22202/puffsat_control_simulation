# C2 coefficient knowledge is prior-bound at the burn; accelerometer demoted to C3 force sensing

**Status:** accepted

## Context

C2 (design doc §13, ADR 0010) was specced as "requirement ≈ free + feasibility": the
coefficient-knowledge *tolerance* falls out of the A3 controllability map with no new
propagation, and the *feasibility* leg confirms "the UKF (accelerometer giving drag
observability) estimates `Cd·(A/m)` / `Cr·(A/m)` to within that tolerance." A grilling
session (2026-06-10, after C1) pushed on both legs and found one gap and two broken
premises:

- **The load-bearing number was never recorded.** A3's finding ("controllable everywhere,
  ~1D in `Cr`, the small gradient is SRP over the coast") is qualitative; the per-cell Δv
  values were not persisted (the grid-resume sink was deferred, ADR 0007). The C2 chain
  rule needs `∂Δv/∂c` numerically.
- **There is no tracking arc before the only high-authority burn.** Deploy ≈ apogee ≈
  burn (B1: "can't center the burn"), and along-track authority is apogee-bound (ADR
  0006). So the coefficient knowledge the corrector acts on is the **ground prior, full
  stop** — no estimator, however good, can improve it in time. "Confirm the UKF estimates
  to within tolerance" silently assumed the estimate exists when it is needed; it cannot.
- **The accelerometer premise is fiction.** §16.4's "GOCE-style" electrostatic
  accelerometer is a kg-class instrument with serious power/thermal demands — not a
  PuffSat payload. A realistic MEMS/tactical-grade accelerometer has bias instability
  ~1–10 µg; SRP at the nominal `Cr·(A/m) = 0.02` is 9.1×10⁻⁸ m/s² ≈ **0.01 µg** — 2–3
  orders below the floor, so **no flyable accelerometer ever observes `Cr`**. Terminal
  drag (peak ~0.017 m/s² ≈ 1700 µg, B3a) is trivially observable — but only during the
  one-shot descent (first perigee *is* the interception; there are no calibration passes),
  and what terminal control needs is the *force* (the product `Cd·(A/m)·ρv²`, which the
  accelerometer measures directly, density conflation and all), not the parameter.

The analytic envelope says the requirement is loose: the SRP Δv-gradient over the
~1.34-day coast is `Cr·(A/m)·P₀·t ≈ 0.011 m/s per 100% Cr error`; through Φ's transverse
amplification (2.15×10⁵ m per m/s, ADR 0011) that bounds the lateral miss at ≤ ~2.3 km
per *doubling* of `Cr` — so the 5 km catch radius tolerates a factor-~2 coefficient
error, and a ground prior (manufactured balloon area to a few %, material `Cr` to
~10–20%) covers it by ~an order. But the envelope is an upper bound on sensitivity built
from an unrecorded gradient; C2a measures it.

## Decision

1. **C2a is the pure tolerance computation, and it gates everything downstream.** The
   chain: the corrector burns `Δv(ĉ)` believing ĉ while truth flies c, so the residual
   interception miss ≈ `Φ_vel · [Δv(ĉ) − Δv(c)] ≈ Φ_lat,vel · G · δc` with
   `G = ∂Δv/∂c` — giving **tolerance = catch radius / ‖Φ_lat,vel · G‖** (factor units,
   matching A3's axes). `G` is *measured* from 1D A3 cuts (3-point `Cr` cut and `Cd` cut
   via the existing `run_sweep`; the per-cell Δv **vectors** come from
   `RunRecord.control_log`, which A3's scalar map never used), and cross-checked against
   the analytic SRP impulse `Cr·(A/m)·P₀·t_coast`. Both axes are measured; `Cd` is
   expected flat at apogee (A3) and its terminal role belongs to C3.

2. **The requirement verdict is "does the ground prior cover the tolerance"** — not "can
   a filter converge." Because deploy ≈ apogee ≈ burn, prior knowledge is what the
   corrector acts on; the comparison is tolerance vs prior σ (default 0.2 of nominal,
   conservative for area-few-% + material-`Cr`-~10–20%). If the prior covers it, C2
   closes as a documented finding in the A2 pattern: *in-flight coefficient estimation is
   unnecessary for the midcourse, because the only burn that could use a better estimate
   happens before any estimate can exist.*

3. **The accelerometer is demoted from C2 and reassigned to C3 as a direct
   drag-disturbance sensor.** It cannot observe `Cr` (ever, on a PuffSat); it observes
   terminal drag superbly but as a *force measurement* feeding feedforward/feedback —
   no coefficient estimate in the loop. The `Cd`-*parameter* requirement (prediction over
   the remaining descent) is assessed against C3's predictor needs, in C3. §16.4 gets the
   honest accelerometer-class note.

4. **C2b (second slice): a LinCov coefficient-augmentation layer quantifies the
   contingency.** Augment the C1 filter with coefficient state(s) (the dimension is
   already parameterized, ADR 0012), put SRP/drag in the filter dynamics, and read
   `σ_Cr(t)` along the coast arc — what tracking alone pins the coefficient to *during*
   the descent. That number converts ADR 0006's named contingency ("observable-drift
   correction" at a near-apogee MCC-2, where the orbit lingers and authority is still
   high) from words to a quantity, against the day Rung D tightens the catch radius.
   Built after C2a's verdict, not before.

5. **C2a ships an explicit lateral-miss error budget.** An RSS ledger of the known
   contributions against the catch radius — C1 nav-induced lateral (141 m at the honest
   q), B1 finite-burn erosion (~89 m), the coefficient-prior residual (this ADR), the
   corrector convergence floor — establishing the budget discipline Rung D will need.
   Allocation is reported, not enforced: the ledger shows headroom, it does not gate.

## Considered options

- **Augmented-UKF feasibility run as the primary deliverable (the original spec)** —
  rejected. The comparison is against the prior (decision 2); a filter demonstration
  cannot change what the corrector knows at the burn. The filter-side question that *is*
  real (what does tracking pin `Cr` to during the coast) is C2b's LinCov layer, which
  needs no per-run UKF.
- **Accelerometer-based coefficient observability** — rejected on instrument physics:
  0.01 µg signal vs ~1–10 µg MEMS bias instability. Carrying it forward would repeat the
  C1 lesson (the SRP-scale q fiction) at the hardware level.
- **Trust the analytic envelope without a live cut** — rejected. Six corrector runs
  (~3 min) buy a measured `∂Δv/∂c` and the Δv *vector* the chain rule actually wants;
  ADR 0007's unrecorded per-cell values are exactly the gap being closed.
- **Re-run the full A3 grid with a sink to bank all Δv vectors** — deferred. The 1D cuts
  answer C2a; the grid-resume sink stays deferred until a wider grid needs it (ADR 0007).
- **Skip the error budget** — rejected. Cheap, and the catch radius is now being consumed
  by three measured residuals from different rungs; Rung D needs the ledger anyway.

## Consequences

- New pure `puffsat_sim/coeff_requirement.py` (Δv-vector cut → gradient, sensitivity,
  tolerance, analytic SRP cross-check, error-budget ledger, report formatting), TDD'd.
- `montecarlo.py` gains a thin `coeff_requirement_report` seam reusing `run_sweep` (1D
  cuts), `run_nav_sweep` (Φ), and `_c0_controller` — no new physics path. Following the
  B2 precedent, no slow integration test: the verdict is a measured finding and the seam
  is glue over already-integration-covered harnesses.
- Design doc §13 C2 is reframed (requirement vs prior; feasibility leg replaced by the
  C2b LinCov layer); §16.4 gets the accelerometer-class reality note.
- The accelerometer moves to C3's slice spec as a force sensor; C3 inherits the
  `Cd`-parameter prediction question.
- If C2a's measured tolerance lands tight (it should not, per the envelope), the recorded
  fallback is the C2b-quantified MCC-2 route — not better instruments.
