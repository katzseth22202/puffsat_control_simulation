# Surveyor-anchored centering: a pure sizing module for the deferred ~10 cm Tier-2 claim

**Status:** accepted

## Context

The closed-loop Monte Carlo (Rung D / D1, ADR 0018) closed out **Tier 1**: a PuffSat captures a
**5 m** pusher plate at ≥ 99 % per-unit confidence on the dumb C baseline, conditional on the fused
terminal-nav grade. The white paper carries a deferred **Tier 2** tightening (CONCLUSION.md;
CONTEXT.md *Surveyor-anchored centering*, grill 2026-06-16): shrink the plate toward **~10 cm** by
driving the *per-unit* arrival miss to centimetres with a surveyor-anchored, camera-rigid swarm.

That regime is **not** a control problem — cm trims sit deep inside the ~475 m catch-radius funnel,
so authority is not the binder and D2/MPC does not return (ADR 0021). It is
**knowledge/metrology-limited**, and its binding numbers are two *hardware characterizations* a Monte
Carlo cannot produce (the hoop precision σ_hoop and the camera's calibrated distortion floor).
CONCLUSION.md Tier 2 already named the right next rigor for it: *"a pure sizing module outputting
plate-size-vs-metrology sensitivity — matching the `tracker_budget.py` / `apogee_nav.py` /
`distortion_field.py` precedent — not an MC."* This ADR records building exactly that, plus the
beacon-design decision the sizing forced.

## Decision

1. **Build `centering_budget.py` as a pure sizing module (no JVM), following the
   `tracker_budget` / `apogee_nav` / `distortion_field` precedent.** It composes the per-unit
   arrival σ from the two levers and applies the committed 99 % 2-D Rayleigh capture criterion
   (`plate = 3.03·σ_arrival`, R/σ from ADR 0015) to output a plate-size-vs-metrology sensitivity
   curve with reference-hardware anchors. It changes **no** committed Tier-1 number.

2. **The plate is the RSS of two legs**, `plate = 3.03·√(σ_hoop² + (σ_θ·v/f)²)`:
   - **Block shift (common-mode bias)** — a sacrificial "surveyor" PuffSat read by an *independent*
     one-shot hoop (lidar/microwave trilateration, **not** the optical tracker, which would be
     circular against its own bias) pins the swarm's quasi-static optical-distortion bias to the
     plate. Residual = the hoop precision **σ_hoop** (1σ), inherited by every unit.
   - **Per-unit scatter (independent)** — the camera-rigid bearing solution, **σ_θ·v/f** (the
     committed CONTEXT.md formula: the plate-crossing cadence `f` sets the nearest-follower spacing
     `v/f`, so a faster train is a closer, sharper link). σ_θ is the RSS of photon-limited
     centroiding and the calibrated distortion floor, reduced across N cameras by the equicorrelated
     ρ law reused from `distortion_field`.

3. **The intra-train link is made distortion-limited by a Q-switched, somewhat-directional beacon —
   not photon-limited.** A naive *wide-cone 1 W CW* beacon on a 5 mm gram-scale aperture at the
   km-class link is **photon-limited** (~4 µrad, over the 3 µrad floor), which would break the
   thesis. Two independent levers each clear it, and the design uses both for margin:
   - a **Q-switched beacon** — bright ns pulses at known timings (~100 kW peak, 10 Hz rep, a few
     hundred mW *average* → a 30 mJ / 300 ns pulse), read in a gate matched to the pulse. The
     measurement collects the whole pulse *energy* (avg ÷ rep), so ~10⁵ photons reach the tiny
     aperture even at km range; the ns gate also freezes the motion-smear term, and a narrowband
     line filter keeps the read signal-dominated against stray light (the surveyor looks backward at
     dark sky). This reuses the Q-switched 1064 nm lineage of the near-Sun extension.
   - a **somewhat directional** beam — coarse-pointed a few degrees along the train axis toward the
     surveyor (no fine tracking), buying ~an order of magnitude of intensity over the wide cone.

   Together they put the photon term at **≈ 0.18 µrad, ~17× under the distortion floor**, so σ_θ is
   set by the *calibrated distortion floor* with comfortable margin. (This also corrected a draft
   bug: the link cadence is the committed 2–4 Hz band, `v/f ≈ 5.4 km` at 2 Hz — not a 1 Hz /
   10.8 km link.)

4. **Finding: 10 cm is robust with existing-instrument-class hardware; 5 cm needs both legs
   tightened.** At the nominal point — σ_hoop 1 cm (rendezvous-lidar class) ⊕ scatter 1.62 cm
   (3 µrad calibrated camera × 5.4 km) → arrival 1.9 cm → **5.8 cm plate (87× off the 5 m
   baseline)**; the 10 cm target tolerates σ_hoop ≤ **2.9 cm**. **5 cm is not reached at nominal**:
   the scatter leg alone (1.62 cm) sets a ~4.9 cm floor — a mm-class hoop by itself lands at
   **4.99 cm** (the module's rangefinder reference point), formally under 5 cm but with zero
   margin — so *robust* 5 cm requires a mm-class hoop **and** a smaller scatter — a 4 Hz train
   (→ 0.8 cm), diverse cameras (ρ=0 → 0.94 cm), or a lower distortion floor. This is a **refinement** of the earlier "σ_hoop ≤ 1 cm → 5 cm" shorthand, which credited
   only the hoop and ignored the co-binding scatter leg. **2 cm dropped.**

5. **The distortion floor is thermally robust because thermal distortion is smooth (low
   spatial-frequency), which is exactly where differential-star astrometry cancels it.** From
   `distortion_field.py`: the nearest-star separation is Δθ ≈ 2.99 mrad and the break-even
   correlation length L\* ≈ 2.54 mrad (L > L\* → differencing helps; L < L\* → √2 worse). Thermal
   distortion is a low-order optical figure change (focus / astigmatism / coma; thermal diffusion
   smooths out high frequencies) → L ≫ L\* → the differential residual is ≪ the absolute floor
   (10× at L = 30 mrad). So the division of labor is: **bench-calibrate the static high-frequency
   fabrication pattern** (thermally stable → calibrate-once works), and **differential-astrometry
   out the smooth, drifting thermal part** (frozen within a ms/ns exposure). Active thermal
   management is a **backup / margin** item, not the primary lever. Near the Sun this becomes a
   **reflective** heat-reject / narrowband front element (dumps the broadband load back out before
   it reaches the sensitive optics; only ~1 % absorbed, and its common-path smooth distortion is
   itself cancelled by the star cross-check) plus artificial-star beacons (the corona blinds natural
   stars) — the in-band coronal/heat-shield IR the filter *cannot* reject is a background/noise axis
   handled by ns time-gating, and remains one of the open near-Sun numbers.

6. **What "physically diverse" (ρ→0) means, and what does not qualify.** The ρ knob is
   **detector-to-detector** correlation across *physically separate detection chains*. Two beacons
   at different wavelengths imaged on **one** camera do **not** give N=2 diverse detectors: the
   dominant distortion terms — detector geometric (pixel-grid) distortion and the thermal figure —
   are **achromatic**, so both colors see the same field (ρ ≈ 1). What two colors *do* decorrelate
   is the **chromatic** distortion (lateral color / chromatic aberration), which they also let you
   *measure and subtract* as a targeted self-calibration. So color diversity is a self-cal for the
   chromatic component, **not** a substitute for spatial diversity (the array) or the co-flyer
   against the achromatic detector/thermal floor. (Corollary for differential astrometry: measure
   the reference in the beacon's band to avoid a lateral-color residual.)

## Consequences

- **The Tier-2 10 cm claim is now *sized*, not hand-waved**: reachable with rendezvous-lidar-class
  metrology and a calibrated star-tracker-class camera, with margin, on legible reference-hardware
  anchors. It remains **argued/sized (Tier 2)**, not simulated — see the boundary below.
- **Nothing new is simulatable for the 10 cm claim with the current architecture.** Its binders are
  bench characterizations a Monte Carlo cannot produce, and its dynamics sit trivially inside the
  already-simulated 475 m funnel — a closed-loop surveyor MC would only re-derive the Rayleigh
  distribution the sizing already applies. The **simulated frontier stays the 5 m plate** (Tier 1,
  D1 closed out); the right next rigor for 10 cm is a **bench test**, not a sim. This is consistent
  with the ADR 0018 project-closeout (falsification test: sharpen-only items are retired).
- **Reversibility.** Paper-side; the module imports the canonical constants from
  `tracker_budget` / `guidance` / `distortion_field` and changes none. A live re-key would be
  warranted only if a bench measurement contradicts an input (e.g. the distortion floor's spatial
  spectrum falls below the differential break-even, or the hoop class cannot reach ~cm).
- **Docs:** CONTEXT.md *Surveyor-anchored centering* updated (two-leg RSS, Q-switched beacon, refined
  5 cm contingency, sized-not-simulated boundary); CONCLUSION.md Tier 2 updated to point at the built
  module. No change to any ADR 0015–0021 committed number.
