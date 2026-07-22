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

2. **The target is constant-velocity near the meet.** The straight line is *not* asked to
   cover the terminal flight: `v_target` enters only the frame direction `v_rel` (one velocity,
   evaluated at the meet) and the target's displacement over `toa_s − target_toa_s` plus the
   closest-approach coast `dt` — together of order 10⁻² s, over which even a 15 m/s² ascent
   displaces the target by millimetres. The approximation is therefore excellent *for the
   scoring*; what is deferred is a dispersed / powered target trajectory and, more importantly,
   feeding target motion to the guidance (decision 4).

3. **The nominal target model is a mirror-ascending rocket.** `reflect_radial_velocity` mirrors
   the PuffSat's crossing velocity across the local horizontal (radial sign flipped, horizontal
   and speed preserved) — a clean single-parameter stand-in for an ascending launch rocket at
   200 km. It is **not** a pessimistic corner: because interception is nearly horizontal, a
   *faster* co-directional target closes *more gently*, so the mirror is the **minimum**-closing
   member of the scaled family `v_target = s · v_mirror` (see Result). The actual target
   trajectory is a mission-design input to be supplied later; this is the placeholder that makes
   the machinery exercisable.

4. **The homing flight is unchanged; only the scoring is relative.** `v_target` affects the
   plate frame, not the flown trajectory — the loop still homes on the *fixed* aim point
   `target_position_m` — so one flown set is re-scored in both frames and the fixed-vs-moving
   delta is purely geometry, not a different flight. This is why the crossing state is now
   retained on `GuidanceRun`. It is also the **main limitation** of this ADR: a real moving
   target requires the aim point to *lead* the target, and no lead is modeled here. This is a
   relative-frame sensitivity study, not an end-to-end powered-target simulation.

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
path angle γ ≈ −8.5° from vis-viva on the reference orbit). A mirror-ascending target therefore
climbs at ~8° with nearly the same horizontal velocity, so the two bodies are close to
co-moving horizontally and `|v_rel| ≈ 2·v_radial ≈ 3 km/s`, not the ~18 km/s a naive head-on
picture suggests. The gentler closing velocity slightly *tightens* the plate-frame miss.

### The two rows bracket a family

Scaling the mirror velocity, `v_target = s · v_mirror`, sweeps a one-parameter family of
co-directional targets. Closing speed rises **monotonically as the target slows** (Kepler
check on the reference crossing):

| `s` | 1.00 (mirror) | 0.75 | 0.50 | 0.25 | 0.00 |
|---|---|---|---|---|---|
| `\|v_rel\|` km/s | 3.19 | 3.86 | 5.85 | 8.24 | 10.78 |

Both ends of this family are **already measured**: `s = 1` is the moving row (σ 0.58 m, 100 %)
and `s = 0` is exactly the fixed-point frame (σ 0.79 m, 98 %). So for any co-directional
ascending target — including a realistically *slower* one than the mirror — the committed
fixed-frame numbers are the conservative end. That, not "the mirror is deliberately fast", is
why the fixed frame is safe to keep.

**This result is target-model-dependent and the family is not exhaustive.** A head-on or
high-cross-track launch geometry leaves the scaled family entirely and raises `|v_rel|` *past*
the PuffSat speed (a fully head-on mirror reaches ~21.6 km/s), which the fixed frame does **not**
bound and which could erode the margin. The closing geometry is now the explicit, single knob;
pinning it down is a mission-design / expert question, and the code is ready to consume the real
trajectory the moment it exists.

## Consequences

- Terminal capture in the simulation is now *scored* against a **moving** target under the
  mirror model, not a fixed point — but it is still *flown* to a fixed aim point (decision 4).
  The Tier-1 conclusion is unchanged (and mildly strengthened) for this model; the honest caveat
  becomes *"scored relative to a placeholder target, sensitive to the launch geometry, no target
  lead in the guidance"* rather than *"the target is not moving at all."*
- **Deferred (non-blocking):** (a) a dispersed / powered target trajectory and the target's own
  nav uncertainty — note the terminal tracker homes on the target beacon, so target *position*
  dispersion is largely nulled by homing; the residual is target *velocity* uncertainty feeding
  the plate frame; (b) a closing-geometry sweep (head-on / cross-track) once the real launch
  azimuth-vs-orbit relationship is known; (c) wiring the moving frame into the D1.1 train
  ensemble and the tail-capture estimator (today they run the fixed frame, which this ADR shows
  is conservative for the mirror model).
- Reversible: the entire change is gated behind `target_velocity_m_s`, default zero.
