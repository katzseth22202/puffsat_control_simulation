# Rung C1: navigation feasibility is a LinCov sweep against C0's Φ — derived node geometry, range + Doppler from a clock-free coherent transponder, an owned typed UKF validated by NEES on a few truth runs

**Status:** accepted

## Context

C0 (ADR 0011) measured the navigation requirement: the binding axis is **apogee transverse
velocity at ~2.3 cm/s** (5 km catch radius), with position tolerances km-scale or unbounded.
C1 must settle the sensing architecture ADR 0011 decision 7 carried forward, and answer
*feasibility*: can a realistic measurement suite, filtered by a realistic onboard estimator,
deliver a covariance Σ at the apogee correction node with `Φ Σ Φᵀ` inside the catch radius?
A grilling session (2026-06-10) walked the decision tree: node geometry, measurement suite,
filter architecture, and the sweep engine. The reframe a future reader needs: C0 moved the
requirement from *position* to *velocity*, so the architecture question is no longer "how well
can multilateration fix position" but "how is the **transverse velocity** observed" — LOS
Doppler directly, and/or ranges filtered through the dynamics over the slow apogee arc.

## Decision

1. **Coordinator nodes co-fly on matched-period orbits with perigee ≥ 200 km; in the sim they
   are known-ephemeris beacons.** Nodes deploy with the swarm on neighboring orbits with the
   **same semi-major axis** (no secular along-track drift; the formation re-converges every
   orbit) and small e/i/RAAN differences that *provide* the ~10³ km apogee baselines. Perigee
   is raised to **≥200 km with apogee trimmed to keep `a`** — nodes are reusable capable assets
   (dish, OCXO, compute) and must never enter the burn-up zone (a node on the PuffSat's 50 km
   perigee would be destroyed on its first pass); raising perigee *alone* would cost ~320 s of
   period mismatch ≈ 144 km/orbit of apogee drift. The perigee lever is cheap at apogee
   (~30 km per m/s → 150 km ≈ 5 m/s), so longevity-minded nodes sit higher (~400–600 km), and
   end-of-life disposal is the same burn reversed. In the filter a node is an **input, never a
   state**: its ephemeris error folds into the measurement noise `R` as an inflation term —
   estimating node states inside every PuffSat's filter contradicts the paper's asymmetry
   (gain, steering, compute, and self-navigation live on the capable node).

2. **Node geometry is a derived requirement, not an assumed constellation.** Formations that
   are "nearby with good geometry" *everywhere* on an e ≈ 0.92 orbit do not exist (speed runs
   0.45→10.9 km/s, separations breathe ~24×; a 1° plane difference is ~2700 km cross-track at
   apogee, ~110 km at perigee). They don't need to: coordinator ranging must win only in the
   **high-altitude coast where GNSS is unavailable** — exactly where slow dynamics keep a
   co-deployed cluster coherent for days — and the binding C0 requirement lives at the apogee
   node, where the cluster is tightest. C1 therefore **sweeps geometry** (node count × LOS
   angular diversity) and outputs the GDOP-style requirement — "N nodes / this spread / this
   Doppler quality hits 2.3 cm/s transverse" — as a paper sizing result; geometry quality vs
   orbit position (including the mid-descent hand-off, where ADR 0006 says little authority
   lives anyway) is *reported*, not required.

3. **Measurement suite: two-way range + two-way carrier Doppler per visible node; whether
   Doppler is load-bearing is a measured output.** Two routes to transverse velocity: direct
   (carrier Doppler, ~mm/s LOS, needs LOS diversity) and indirect (ranges filtered through the
   dynamics — ~10 m fixes differenced over ~10³ s already give ~1.4 cm/s, and the apogee arc is
   *days*). Both `h(x)` are trivially pure, so C1 models both and includes a **range-only sweep
   point** (Doppler σ → ∞): if range-only clears the bar, the paper claims less hardware.
   Cadence is swept coarsely (**0.003 / 0.03 / 0.3 Hz**) since cadence × arc is what the
   indirect route trades against the direct one. *Link-budget sizing note (recorded for the
   paper):* at S-band over 3000 km (FSPL ≈ 169 dB), 10 mW into the PuffSat omni and a ~25 dBi
   coordinator dish give ~30 dB carrier-loop SNR (10 Hz loop) → thermal two-way Doppler
   ~0.01 mm/s at 60 s; the **binding term is the coordinator oscillator** (OCXO 1×10⁻¹² →
   ~0.15 mm/s two-way), and the **coherent turnaround means the PuffSat carries no clock** —
   one receive-mix-amplify MMIC chain + patch antenna at mW power is how "few grams" stays
   honest. Multi-node service is time-multiplexed turnaround (trivial at these cadences).
   mm/s-class LOS Doppler therefore has ~10× margin, but Doppler σ stays a swept axis
   (~0.1–10 mm/s) so the conclusion rests on the measurement, not the link budget.

4. **GNSS: modeled as a position fix, reported not swept; its necessity is a measured claim.**
   Near perigee the PuffSat is under the constellation looking up — easy availability — and
   `h(x)` is the cheapest of all (direct position + noise). Hardware is honestly **tens of
   grams, not few**: COTS receivers are COCOM-limited (~515 m/s) and the PuffSat crosses 200 km
   at ~10.8 km/s — 20× over, beyond even standard LEO space receivers — so it must be an
   unlocked spaceborne receiver with high-dynamics tracking loops, ~30–50 g installed with
   antenna (~0.15% of the 25 kg PuffSat; same scale as the §13 50 g compute note). C1 reports
   the descent covariance hand-off **with and without GNSS**: if the coordinator-only suite
   meets the threshold end-to-end, the paper gets "GNSS optional" — stronger than assuming it.
   The terminal *relative*-nav question (plate-relative aim; absolute GNSS can't answer it
   alone) is **C3's**, flagged forward.

5. **An owned, fully typed UKF — no FilterPy, in the loop or as oracle.** The UKF core (sigma
   points, unscented transform, predict/update) is ~150 lines of NumPy; FilterPy is untyped and
   dormant (no release since ~2018) — an `ignore_errors` island in exactly the pure core we
   keep strictly typed — and its plain-Cholesky UKF is known to fail on ill-conditioned
   covariances like ours (state magnitudes span ~10⁸ m to ~10² m/s). Owning it keeps the
   square-root/Joseph-form option open and makes the TDD'd filter part of the paper's
   verification story (ADR 0003 "own the solver"). Verification needs no library: on a linear
   system the UKF must reproduce the **analytic Kalman filter exactly**, plus one nonlinear
   sanity case. Start standard-form with symmetrization; go square-root only if tests force it.

6. **Filter state: 6-state Cartesian position/velocity in EME2000, dimension parameterized;
   dynamics: pure-Python two-body + J2; `Q` is a swept knob, never quietly tuned.** Cartesian
   because the Gaussian assumption behaves better there than in near-parabolic elements at
   e ≈ 0.92 and every `h(x)` is Cartesian-natural; parameterized dimension so C2 appends the
   two coefficient states without surgery. The predict step propagates sigma points with the
   cheap **onboard** model (RK4 substeps of two-body+J2) — the realistic flight architecture
   (no flight filter runs Orekit), fully TDD-able, and it makes the truth−filter model gap a
   *measured* quantity absorbed by `Q` at C1 (omitted SRP ~5×10⁻⁸ m/s², third-body tidal
   residual) and by the estimated coefficients at C2. `Q` (white unmodeled-acceleration) is an
   explicit sweep axis reported with the envelope.

7. **Two-layer method: a LinCov sweep engine, validated by a few seeded UKF runs on truth with
   NEES (the load-bearing compute decision, mirroring ADR 0011 decision 2).** The Kalman
   covariance recursion does not depend on measurement *values* — pinned to the reference
   trajectory it is deterministic (classical linear covariance analysis), reusing the same
   TDD'd machinery with no truth propagation. **Sweep engine:** for each (range σ, Doppler σ,
   cadence, geometry, `Q`) cell, run the covariance recursion along the coast arc → Σ at the
   apogee node → `Φ Σ Φᵀ` vs the catch radius → the feasibility envelope. **Validation:** at a
   handful of cells (nominal + the envelope edge), generate seeded synthetic measurements from
   an Orekit truth arc, run the real UKF, and check **NEES consistency** — the filter's actual
   error statistically inside its claimed Σ. This is what LinCov cannot certify (sigma-point
   spread through real nonlinearity over 33–333 s gaps at e ≈ 0.92, the truth−filter gap `Q`
   must absorb); if NEES fails, the envelope is fiction and `Q` is retuned. Same epistemic
   shape as C0: a measured requirement envelope + a validated-filter existence proof, no
   ensemble.

8. **LinCov is a Rung-D accelerant, not a replacement for the Monte Carlo.** LinCov breaks
   exactly where D's question lives: actuator saturation (400 mN / 5 mN floor), propellant
   depletion, corrector non-convergence, and regime switches destroy superposition; log-normal
   coefficients and F10.7/Ap drive non-Gaussian density response; and the headline result is a
   **tail probability** (P(miss > catch radius), P(propellant > 2%)) — precisely where Gaussian
   projections lie. But the C1 engine makes D *small*: ADR 0010 decision 2 (sample nav error
   from Σ instead of running the filter per trajectory) is already LinCov-inside-MC; add LinCov
   as the **screen** (which axes matter, where the cliffs are) and as a **control variate**
   (regress MC results on the LinCov prediction; only the nonlinear residual needs samples,
   plausibly 10–100× fewer for the same confidence). MC stays the certifying instrument.

9. **Regime slicing: apogee swept, GNSS reported, accelerometer deferred to C2.** The binding
   C0 requirement lives at the apogee node, so the apogee/coast regime is C1's subject. The
   accelerometer observes non-gravitational specific force — the *coefficients*, C2's
   measurement problem by ADR 0010's own split — and senses ~nothing during the coast (drag ≈ 0,
   SRP at the noise floor of any few-gram MEMS unit); pulling it into C1 would re-couple what
   ADR 0010 deliberately decoupled.

10. **Layout: one pure module `puffsat_sim/estimation.py`** (sigma points, UKF predict/update,
    measurement models `range`/`doppler`/`gnss_fix`, node-ephemeris helpers, LinCov runner),
    with seeded validation runs entering through a `montecarlo.run_nav_feasibility` JVM seam
    mirroring `run_nav_sweep`. Split the LinCov runner out only if the module outgrows itself.

## Considered options

- **FilterPy (or pykalman / Stone Soup) in the loop** — rejected (decision 5): untyped + dormant
  in the strictly-typed pure core, plain-Cholesky conditioning risk, and it would move "the
  filter is correct" off our test suite. Even as a test oracle it adds a dependency weaker than
  the analytic linear-KF equivalence check.
- **Orekit inside the filter predict** — rejected (decision 6): 13 JPype propagations per cycle,
  unrealistic flight architecture, and filter dynamics = truth dynamics tests less (the filter
  can never be wrong about the physics).
- **Node states estimated in the filter** — rejected (decision 1): blows up the state vector
  6-per-node and contradicts the capable-node asymmetry.
- **Assume a constellation geometry** — rejected (decision 2): deriving the geometry requirement
  is the stronger paper claim and the sim never owns an arbitrary constellation.
- **Full UKF truth runs per sweep cell** — rejected (decision 7): the covariance recursion is
  measurement-value-independent; running truth arcs across a 5-axis grid re-imports the compute
  blow-up ADR 0010 banished.
- **LinCov replacing the Rung-D Monte Carlo** — rejected (decision 8): saturation/discrete
  events, non-Gaussian drivers, and tail-probability headlines; LinCov screens and
  variance-reduces instead.
- **Accelerometer in the C1 suite** — rejected (decision 9): it observes coefficients (C2), not
  position, and is blind during the coast.
- **Keplerian/equinoctial filter state** — rejected (decision 6): near-parabolic e ≈ 0.92 makes
  element-space Gaussians badly behaved at exactly the regime C1 studies; Cartesian + RK4
  substeps is boring and testable. (Element formulations stay on the table for the *truth*
  integrator at C3 — different problem.)
- **Node perigee at the PuffSat's 50 km / perigee raised without trimming apogee** — rejected
  (decision 1): the first destroys the node on its first pass; the second breaks the matched
  period that keeps the formation drift-free.

## Consequences

- C1 stays cheap: the sweep is pure Python (no truth propagation); only the handful of NEES
  validation runs touch the JVM. Pure blocks — sigma points, UT, predict/update, `h(x)` models,
  LinCov recursion, NEES bounds, node-ephemeris helpers — are **`/tdd`**; the feasibility
  envelope, geometry requirement, range-only verdict, and with/without-GNSS hand-off are
  **measured** findings (recorded in design-doc §13, like C0).
- The **C0 seam closes**: C1's Σ feeds `Φ Σ Φᵀ` through C0's measured Φ against the same
  catch-radius table.
- **C2 consumes** the parameterized state dimension and the whole UKF; the accelerometer model
  enters there. **C3 consumes** the GNSS-regime hand-off covariance and the flagged terminal
  relative-nav question. **Rung D consumes** Σ (sampling nav error, ADR 0010) and the LinCov
  engine (screen + control variate).
- **Paper sizing notes banked here**: the clock-free coherent-transponder link budget (OCXO-bound,
  ~10× Doppler margin at 3000 km), the GNSS mass/COCOM reality (~30–50 g, unlocked receiver),
  and the node constellation principle (matched `a`, perigee ≥200 km, ~5–15 m/s
  raise/disposal burns).

## Implementation findings (2026-06-10)

Built across three commits — `73f88bc` (pure `estimation.py`: owned UKF, two-body+J2 flow,
measurement models, LinCov, NEES), `b8d9a9d` (pure `nav_feasibility.py`: one-axis-at-a-time
sweep, cone-geometry node family, Σ→RTN, `Φ Σ Φᵀ` verdict), `ad3186b` (`validate_cell` +
the `montecarlo` truth-arc seam) — then run for real: Φ re-derived live (a 13-cell
`points_per_sign=1` C0 sweep, 137 s; lateral column norms match ADR 0011's full sweep), the
15-cell LinCov envelope (~50 s), and the seeded NEES validation against full-force truth.

1. **The validation layer earned its place immediately: the SRP-scale q was fiction.** At the
   grilled nominal `q = 5e-8 m/s²` the filter *claimed* σ_Tvel ≈ 6 µm/s while its *actual*
   error against full-force truth was **0.09 m/s — above the 2.3 cm/s C0 requirement itself**
   (average NEES ~8×10⁸ vs bounds [5.81, 6.19]). The truth−filter gap at apogee is dominated
   not by SRP but by the **third-body tidal acceleration** (~3×10⁻⁵ m/s², Moon+Sun — three
   orders above SRP). A LinCov-only C1 would have shipped this fiction; decision 7's layer 2
   caught it on the first run.

2. **The q-ladder localizes the consistency crossover at the tidal scale, and the verdict
   survives there.** Measured (12 h arc, nominal suite): q = 5e-8 → NEES 8e8; 5e-7 → 3×10⁶;
   5e-6 → 7×10³; **3e-5 → 39** (claimed 0.40 mm/s ≈ actual 0.42 mm/s — magnitudes agree, the
   residual NEES is the colored-vs-white signature of a deterministic tide under a white Q);
   **1e-4 → 5.2, honest** (just under the lower bound — mildly pessimistic; actual error
   0.08 mm/s, ~8× under claimed). `NavFeasibilitySpec.nominal_q_accel_m_s2` now carries the
   validated **1e-4**, with 5e-8 (the trap) and 3e-5 (the crossover) as the swept points.

3. **C1 verdict: the navigation requirement is met with ~35× margin at the honest q.** All 15
   envelope cells MEET the 5 km catch radius; at the validated nominal, σ_Tvel = 0.66 mm/s vs
   the 2.3 cm/s requirement and lateral 1σ miss = 141 m vs 5 km. Every swept degradation
   (100 m ranges, 5° LOS cone, 3 nodes, 0.003 Hz cadence) stays under by ≥ an order. The
   1 km catch radius (4.7 mm/s) also clears ~7×.

4. **Doppler is load-bearing at the honest q — the apparent range-only tie was a q=5e-8
   artifact.** At the fictional q the process noise was so tight that velocity knowledge
   accrued through the dynamics (the grilling's "route 2") regardless of measurement type,
   and the two cells tied. The canonical run at the validated q = 1e-4 separates them:
   with-Doppler claims σ_Tvel = 0.66 mm/s and is NEES-honest (5.18, mildly pessimistic);
   range-only degrades ~2.7× to 1.8 mm/s and goes mildly OPTIMISTIC (NEES 9.4 vs upper
   bound 6.19 — mild, not the 10⁸ pathology; its actual T-vel error 1.33 mm/s stays under
   its claim). Range-only still MEETS with ~13× margin (lateral 386 m vs 5 km), so the
   paper claim *can* fall back on ranging alone, but a range-only flight would want its own
   slightly higher q — the with-Doppler nominal is the validated configuration.

5. **The 100 m catch radius is the marginal case and names the upgrade path.** Its tolerance
   (0.47 mm/s, ADR 0011) sits *at* the honest-q claimed σ (0.66 mm/s fails it; the 3e-5
   crossover's 0.40 mm/s grazes it). If plate-scale hand-off ever tightens to ~100 m, the
   lever is **third-body in the onboard filter dynamics** (recovering ~2 orders of q), not
   better sensors — recorded for C2/C3, consistent with decision 6's "Q absorbs the gap"
   stance holding only down to the km-scale radii C0 actually requires.

6. **Mechanics that held up**: the equatorial pure reference vs the 28.5°-inclined truth
   orientation is immaterial (everything is judged in RTN — pinned by integration test);
   truth arcs cache per cadence; the full default report (Φ + sweep + 2 validations) runs in
   ~4 min wall.
