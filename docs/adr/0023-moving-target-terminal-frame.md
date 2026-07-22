# The terminal plate frame is relative to a moving target, not a fixed interception point

**Status:** accepted

## Context

Through Rungs C and D the terminal loop scored capture against a **fixed point** in inertial
space: the PuffSat's own nominal drag-free 200 km crossing stood in for the target rocket's
pusher plate (`build_guidance_context`). A design review (2026-07) correctly flagged that the
mission target is a *moving* rocket, and the design doc (§ Rung D) always called for the
prescribed moving target. With a fixed point:

- the plate-frame miss is measured ⊥ the **PuffSat's** inertial velocity, not ⊥ the **closing**
  velocity `v_rel = v − v_target`;
- the time-of-arrival lever uses the PuffSat speed, not the closing speed;
- the target's own motion between the nominal meet time and the PuffSat's actual crossing time
  is ignored.

`PlateMiss` already *documented* itself as "⊥ v_rel … relative velocity" — the intent was
always relative kinematics; the implementation had simply hard-coded `v_target = 0`. So this is
a correction of a modeling shortcut, not a new subsystem.

This ADR does **not** attempt a full launch-to-200 km ascent model. It installs the relative
kinematics and one physically-grounded target model, and records what the fixed-point
approximation was actually costing.

## Decision

1. **The plate frame is relative.** `plate_frame_miss` takes `target_velocity_m_s` (default
   `(0,0,0)`), models the target constant-velocity through its nominal crossing
   (`p_target(t) = target_position_m + v_target·(t − target_toa_s)`), and decomposes the
   crossing about the closing velocity `v_rel = v − v_target`. `v_target = 0` reproduces the
   historical fixed-point frame **byte-identically**, so every committed C3b/D1 number is
   unchanged and the moving target is strictly additive.

2. **The target is constant-velocity near the meet.** Over the ~30 s terminal window the
   rocket's ascent acceleration is second-order; a straight-line target isolates the
   first-order effect (closing geometry) that the fixed point dropped. A dispersed / powered
   target trajectory is deferred (see below).

3. **The nominal target model is a mirror-ascending rocket.** `reflect_radial_velocity` mirrors
   the PuffSat's crossing velocity across the local horizontal (radial sign flipped, horizontal
   and speed preserved) — a clean single-parameter stand-in for an ascending launch rocket at
   200 km, deliberately fast so the closing speed is on the pessimistic side of a real, slower
   rocket. The actual target trajectory is a mission-design input to be supplied later; this is
   the placeholder that makes the machinery exercisable.

4. **The homing flight is unchanged; only the scoring is relative.** `v_target` affects the
   plate frame, not the flown trajectory (the loop homes on the target *position*), so one flown
   set is re-scored in both frames — the fixed-vs-moving delta is purely geometry, not a
   different flight. This is why the crossing state is now retained on `GuidanceRun`.

## Result (measured, mirror-ascending target, nominal grade, N = 48)

| Quantity | Value |
|---|---|
| PuffSat speed at 200 km | 10.79 km/s |
| **Closing speed** `\|v_rel\|` | **3.03 km/s** (0.28× PuffSat) |
| σ_lateral, fixed frame | 0.79 m |
| σ_lateral, moving frame | 0.58 m (0.73×) |
| Capture (5 m plate) | 100 % (fixed 98 %) |

**The closing speed is *low*, and the fixed-point approximation was mildly conservative, not
optimistic — for this target model.** The reason is geometric: interception at 200 km sits well
above the 50 km perigee, where the descending trajectory is only ~8° below horizontal (flight
path angle γ ≈ −8.1° from vis-viva on the reference orbit). A mirror-ascending target therefore
climbs at ~8° with nearly the same horizontal velocity, so the two bodies are close to
co-moving horizontally and `|v_rel| ≈ 2·v_radial ≈ 3 km/s`, not the ~18 km/s a naive head-on
picture suggests. The gentler closing velocity slightly *tightens* the plate-frame miss.

**This result is entirely target-model-dependent** — it is the mirror corner, not a universal
answer. A head-on or high-cross-track launch geometry would raise `|v_rel|` toward and past the
PuffSat speed and could erode the margin. The closing geometry is now the explicit, single knob;
pinning it down is a mission-design / expert question, and the code is ready to consume the real
trajectory the moment it exists.

## Consequences

- "Interception" in the simulation now means terminal control to a **moving** target under the
  mirror model, not a fixed point. The Tier-1 conclusion is unchanged (and mildly strengthened)
  for this model; the honest caveat becomes *"sensitive to the launch geometry"* rather than
  *"the target is not moving at all."*
- **Deferred (non-blocking):** (a) a dispersed / powered target trajectory and the target's own
  nav uncertainty — note the terminal tracker homes on the target beacon, so target *position*
  dispersion is largely nulled by homing; the residual is target *velocity* uncertainty feeding
  the plate frame; (b) a closing-geometry sweep (head-on / cross-track) once the real launch
  azimuth-vs-orbit relationship is known; (c) wiring the moving frame into the D1.1 train
  ensemble and the tail-capture estimator (today they run the fixed frame, which this ADR shows
  is conservative for the mirror model).
- Reversible: the entire change is gated behind `target_velocity_m_s`, default zero.
