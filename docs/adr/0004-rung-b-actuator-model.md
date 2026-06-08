# Rung B actuator model: proportional abstraction + reported Isp sweep

**Status:** accepted

## Context

Rung B replaces Rung A's impulsive Δv with a finite-burn actuator (§13, §16.5);
A2 itself stays impulsive (the ADR 0003 corrector is unaffected). A grilling
session (2026-06-08) settled the first-cut actuator model and two choices a future
reader would otherwise find surprising.

## Decision

First-cut Rung B actuator: a **single omnidirectional proportional cold-gas
thruster**. The full parameter set lives in §13 (Rung B bullet). Two decisions
are worth recording here because each was a real fork:

1. **Proportional throttle as the modeling abstraction, over bang-bang hardware.**
   We model continuously-variable thrust (400 mN max, ~5 mN floor, no
   quantization) even though a real cold-gas system is on/off valves. Justified by
   §5: a bang-bang cluster (or PWM'd valve) averages to "continuous." Modeling the
   effective proportional behavior keeps the actuator parameter set minimal — no
   minimum-impulse-bit (MIB) bookkeeping — and the MIB residual becomes a deferred
   refinement that lands exactly where the cm-aim 100 Hz trim needs it.

2. **Isp as a reported sweep, not a single fixed value.** Isp is a post-processing
   lever on the propellant ledger (`propellant ≈ Δv/(Isp·g₀)`), not on
   controllability: the A-rung Δv result is computed once and transformed across
   Isp ∈ {50, 70, 200} s. 50 s is the conservative anchor (CO₂/argon/refrigerant-
   class storable propellant — realistic, dense, self-pressurizing); 70 s is N₂;
   200 s is the Appendix A optimistic figure. The sim exists in part to test the
   paper's "<2% of 25 kg" propellant claim, and Isp is the single lever that
   decides it — so the trade must be *shown on a curve*, not pre-decided.

Mass depletion is modeled (Tsiolkovsky, Orekit-native). Dead-time is deferred to
Rung C (§16.8 — perfect state through Rung B gives loop latency no source).

## Considered options

- **Pulsed / PWM with explicit MIB quantization** — rejected for the first cut:
  premature quantization bookkeeping. Reintroduced as the bang-bang-cluster
  residual when the cm-aim demands it.
- **Hardcode Isp = 50 s** — rejected: pre-decides the propellant answer the sim is
  built to measure (it would just report the <2% claim failing).
- **Hardcode Isp = 200 s** — rejected: Appendix A optimistic for true cold gas;
  hides the Isp-vs-mass-fraction trade.
- **Commit now to gimbal vs. whole-body vs. cluster for direction** — deferred: the
  ~1 deg/s slew limit has ~10× margin over the ~0.1 deg/s required sweep, so the
  realization only matters once the limit is shown to bind (it shouldn't).

## Consequences

- The propellant deliverable becomes a **fraction-vs-Isp curve** against the <2%
  line, not a single number; expect the conservative (50 s) anchor to fail it ~3×.
- A2 and the ADR 0003 corrector are unchanged (still impulsive, perfect-state).
- Rung B adds a finite-burn maneuver with mass depletion and a ~1 deg/s direction
  slew; everything else (MIB / cluster geometry, cosine losses, the coarse/fine
  chemical+cold-gas split, and control-loop latency) is explicitly deferred.
