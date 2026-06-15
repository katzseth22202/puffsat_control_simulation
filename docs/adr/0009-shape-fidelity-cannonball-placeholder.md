# Shape fidelity: cannonball as conservative placeholder, mass/coefficient decoupling, optimistic cylinder deferred to Rung E

**Status:** accepted

## Context

B1's finite-burn work runs the propagator at a fictitious 1 kg (so the lumped
`Cd·(A/m)` / `Cr·(A/m)` scale drag/SRP directly; ADR 0008). That kept prompting a
recurring question — *is the propagation pessimistic without the real 25 kg mass, and
should drag/SRP use the real ANFO-cylinder shape instead of a cannonball?* A grilling
session (2026-06-09) settled the sequencing and recorded four things a future reader —
or the recurring "just use the real mass" instinct — would otherwise re-derive. The
paper's `sec:estimate_cold_gas` (GOCE→ANFO cold-gas estimate) is the external anchor.

## Decision

1. **The 1 kg propagator mass is normalization, not pessimism — the lumped coefficient
   is the only drag/SRP lever.** Orekit drag/SRP acceleration is `½ρv²·(crossSection/mass)`;
   the sim passes `crossSection = Cd·(A/m)` with `mass = 1`, which is *algebraically
   identical* to `crossSection = Cd·(A/m)·25` with `mass = 25` — the mass cancels. So
   switching the propagator to the real 25 kg changes drag/SRP by **nothing** unless the
   coefficient is also rescaled; it is a no-op for the perturbations. Only the *burn*
   (`a = F/m`) ever needed the real mass, and B1 handled that by scaling the thrust to the
   fictitious mass (ADR 0008). **Do not "use the real mass for realism."**

2. **The cannonball coefficients are deliberate conservative placeholders.**
   `Cd·(A/m)=0.04` and `Cr·(A/m)=0.02` m²/kg are ungrounded but conservative. Against the
   paper's GOCE→ANFO cylinder (`sec:estimate_cold_gas`: D=1 m, L=4 m, ρ=840 kg/m³ → a
   geometrically-similar 25 kg PuffSat of ~0.63 m² total surface, r≈0.106 m, L≈0.85 m):
   0.04 is ~13× the *streamlined face-on* value (~0.003), ~2.5× broadside (~0.016), and
   ~0.7× the paper's *total-surface-area* drag proxy (~0.056); SRP's 0.02 *exceeds* the
   worst-case broadside-aluminum (Cr≈1.9 → ~0.014). The pessimism is in the directions
   that matter — more drag = more anti-drag propellant *and* more terminal perturbation to
   reject — so a pass under it is *a fortiori* a pass for the real, lower-drag, oriented
   cylinder. The one place the pessimism flips: more drag → lower perigee → burn-up looks
   *easier*; that is the debris-disposal diagnostic ("low = good"), not a mission-kill
   metric, so it is harmless (just don't over-claim deorbit margin off an inflated-drag
   run). The cannonball is also architecturally coherent: it is exactly what the UKF
   estimates as `Cd·(A/m)` / `Cr·(A/m)`, whereas "cylinder" implies an *attitude subsystem*
   (orientation vs. flow and Sun), not a parameter.

3. **B3 grounds the drag coefficient only on a gate failure.** B3 keeps 0.04 and
   cross-checks its measured peak drag force / anti-drag Δv against the paper's 374 g /
   400 mN; it re-grounds the coefficient in the GOCE-ANFO geometry **only if B3 fails
   either gate** — propellant > 2% of 25 kg *or* peak thrust > 400 mN. A back-of-envelope
   (NRLMSISE at 200 km, ρ≈2.5×10⁻¹⁰, v≈10.8 km/s, `Cd·A = 0.04×25 = 1.0 m²`) gives a peak
   drag force of only ~15 mN — ~30× under the 400 mN surge — so a comfortable pass under
   the pessimistic placeholder is the expected, and strongest, result.

4. **The optimistic cylinder is the final rung (Rung E), after D.** A–D run
   cannonball-pessimistic; **Rung E** re-runs the Rung D Monte Carlo with the shape model
   (same trajectory and seeds — a clean A/B), producing the pessimistic-vs-optimistic
   comparison that becomes the paper's follow-up to the GOCE-ANFO appendix.

   > **Amendment 2026-06-15 (after D1 closeout): E1 downgraded from committed to optional.**
   > D1 closed feasible with drag confirmed **non-gating** (ADR 0021: terminal drag is
   > feedforward-solved and the coefficient's only entry is the coast Cr-prior leg, which the
   > Cr-mismatch leg showed is shared/absorbed and 98 % along-track). The cannonball is
   > pessimistic in *every* GNC-relevant direction (terminal rejection, Cr-prior entry,
   > anti-drag propellant), so the lower-drag cylinder can only *improve* margins that already
   > pass — it cannot change the verdict. **E1 therefore survives only as optional paper polish**
   > (the pessimistic→optimistic comparison numbers), not a gate. **The one place the cannonball
   > is *optimistic*, not pessimistic — burn-up / debris disposal** (more drag → easier burn-up;
   > the "low = good" flip already noted in decision 2) — is the only residual that wants the
   > low-drag case, and it is a **reentry-heating diagnostic, not the Rung-E GNC re-run** and not
   > a mission gate (burn-up is the desired outcome on a miss). If the paper claims deorbit
   > margin, check it at the realistic (high-β, low-drag) shape — separately — rather than off the
   > inflated-drag run.

   - **E1 (optional — downgraded 2026-06-15): attitude-dependent area at favorable pointing.** Model the cylinder
     face-on into the flow / its projected area to the Sun, *assuming* good pointing (an
     Orekit panel/box `DragSensitive`+`RadiationSensitive` + an attitude provider swapped
     into D's harness — a force-model swap, not a new harness). Deliverable: the comparison
     numbers. Assuming optimal pointing is honest *because* this is labelled the optimistic
     bound.
   - **E2 (optional): closed-loop attitude pointing feasibility.** Proves the assumed
     orientation is *achievable* against the destabilizing aerodynamic torque on a slender
     body — the predictable reviewer question. Light form first: a torque-margin /
     static-stability analysis (disturbance torque at peak dynamic pressure vs. available
     restoring/control torque — the attitude analog of the paper's drag back-of-napkin),
     appendix-grade and likely sufficient. Heavy form (full rigid-body attitude-dynamics +
     control sim) only if the margin is thin.

## Considered options

- **Switch the propagator to the real 25 kg "for realism"** — rejected: a no-op for
  drag/SRP (decision 1); only the burn needed it, and B1 already scaled the thrust.
- **Ground the drag coefficient now, before B3** — rejected: a pass under the pessimistic
  placeholder is the stronger result; grounding is a contingent follow-up (decision 3).
- **Model the cylinder/attitude now, or fold it into Rung C/D** — rejected: it is a
  subsystem the cannonball folds into the UKF-estimated lumped coefficient, and the
  refinement only *tightens margin in the optimistic direction*, so it belongs last as a
  separate rung for the clean pessimistic→optimistic comparison.
- **Commit E2 (closed-loop attitude control)** — deferred to optional: not needed for the
  propellant claim or the comparison numbers; kept on the roadmap because a reviewer will
  ask whether the favorable orientation holds.

## Consequences

- The recurring "use the real mass" instinct now has a durable answer; cite this when it
  resurfaces. CONTEXT.md gains a **Lumped coefficient** term making the decoupling explicit.
- B3 scope stays narrow: instrument drag at the placeholder, the coefficient is one
  documented lever, grounding is contingent on a gate failure (not a B3 prerequisite).
- The roadmap gains **Rung E** (E1 committed, E2 optional); the paper's GOCE-ANFO appendix
  (`sec:estimate_cold_gas`) is the pessimistic endpoint of E1's comparison.
  **Superseded by the 2026-06-15 amendment (decision 4): E1 downgraded to optional after D1
  closed feasible with drag non-gating (ADR 0021).** Rung E is now optional paper polish; the
  only residual that wants the low-drag shape is the burn-up / debris-disposal claim, which is a
  separate reentry-heating diagnostic, not the GNC re-run.
- The cannonball deferral — previously only a §7 note + a parked-decision memory — is now a
  formal record.
- Still deferred within shape fidelity: optical properties beyond a single `Cr`,
  deployed-structure (fins / bladder / nitrogen-tank) area growth, Earth albedo/IR — all
  below the first-cut floor (§7).
