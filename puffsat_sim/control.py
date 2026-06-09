"""Rung A1 differential corrector — pure targeter (no JVM).

The corrector solves for the apogee correction Δv that nulls the **interception
miss**: the 3-component position miss at the 200 km EME2000 descent crossing
(ADR 0003).  It is a 3×3 root-find — Newton with a finite-difference Jacobian —
parameterized by a ``predict`` callback (correction Δv in RTN → crossing position
in EME2000).  Keeping it a black-box-propagator method (rather than harvesting an
Orekit STM) lets the same solver survive into Rung B (finite burns) and Rung C
(onboard-model mismatch); keeping it pure lets it be unit-tested with a synthetic
``predict``, no Orekit boot.  Non-convergence is a recorded outcome (the Rung A3
authority boundary), never an exception.

The JVM glue that builds the ``predict`` closure and executes the returned plan
lives in :mod:`puffsat_sim.montecarlo`.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat_sim.dispersion import Vec3

# Correction Δv in the apogee RTN frame [m/s] → crossing position in EME2000 [m].
# The RTN→EME2000 conversion lives inside the closure, so the solver works in RTN
# input space and the solved Δv is already RTN (transverse = the dr_p/dv_a lever).
PredictFn = Callable[[Vec3], Vec3]

# Two-burn correction: the stacked (apogee RTN, mid-descent RTN) Δv → crossing
# position.  Each 3-vector is in its own node's RTN frame; the closure re-derives the
# mid-descent basis from the stopped node state (ADR 0005).
Vec6 = tuple[float, float, float, float, float, float]
TwoBurnPredictFn = Callable[[Vec6], Vec3]


@dataclass(frozen=True)
class Target:
    """What the corrector aims at: the nominal crossing position in EME2000 [m].

    Velocity / time-of-arrival are not carried — Rung A1 nulls position only (a
    single impulsive Δv has 3 DOF); later rungs extend this type.
    """

    position_m: Vec3


@dataclass(frozen=True)
class ControlAction:
    """One commanded maneuver: an impulsive Δv at a node (the record, not the burn)."""

    node_label: str
    elapsed_s: float  # time from epoch at which it is applied (0.0 = apogee node)
    dv_rtn_m_s: Vec3
    dv_mag_m_s: float


@dataclass(frozen=True)
class ControlPlan:
    """A Controller's output: the ordered actions plus convergence metadata.

    Deeply immutable (``actions`` is a tuple of frozen ``ControlAction``) so it can
    be logged into a frozen ``RunRecord`` and shared across processes safely.
    """

    actions: tuple[ControlAction, ...]
    converged: bool
    iterations: int

    @property
    def total_dv_m_s(self) -> float:
        return float(sum(a.dv_mag_m_s for a in self.actions))


# A Controller maps the run's prediction model and target to a plan.  Rung A1
# supplies ``solve_apogee_correction``; Rung D supplies MPC; ``control=None`` is the
# open-loop capstone.
Controller = Callable[[PredictFn, Target], ControlPlan]


def passes_toa_gate(converged: bool, toa_miss_s: float, toa_window_s: float | None) -> bool:
    """Whether a converged plan is the intended local root, not the spurious far one (ADR 0007).

    With the step cap raised to budget scale (A3 decision 3ii), Newton can re-converge on a
    far, high-Δv orbit that re-crosses 200 km a revolution off-nominal.  ToA — not Δv
    magnitude — is the physical discriminator: a converged solution whose crossing falls
    outside ±``toa_window_s`` of the nominal time of arrival is rejected as "no valid local
    solution".  ``toa_window_s=None`` disables the gate (A1/A2/the capstone keep their
    verdict); a non-converged plan is never rescued.
    """
    if toa_window_s is None or not converged:
        return converged
    return abs(toa_miss_s) <= toa_window_s


def _vec3_of(a: NDArray[np.float64]) -> Vec3:
    return (float(a[0]), float(a[1]), float(a[2]))


def _vec6_of(a: NDArray[np.float64]) -> Vec6:
    return (float(a[0]), float(a[1]), float(a[2]), float(a[3]), float(a[4]), float(a[5]))


def _finite_difference_jacobian(
    predict: PredictFn, x: NDArray[np.float64], base_position: NDArray[np.float64], step: float
) -> NDArray[np.float64]:
    """Forward-difference ∂(crossing position)/∂(correction Δv), reusing ``base_position``."""
    jac = np.empty((3, 3), dtype=np.float64)
    for j in range(3):
        xp = x.copy()
        xp[j] += step
        perturbed = np.asarray(predict(_vec3_of(xp)), dtype=np.float64)
        jac[:, j] = (perturbed - base_position) / step
    return jac


def _newton_step(jac: NDArray[np.float64], residual: NDArray[np.float64]) -> NDArray[np.float64]:
    """A1's committed square Newton step; raises ``LinAlgError`` on a singular Jacobian."""
    return np.asarray(np.linalg.solve(jac, residual), dtype=np.float64)


def _lm_step(
    jac: NDArray[np.float64], residual: NDArray[np.float64], lm_lambda: float
) -> NDArray[np.float64]:
    """Levenberg-Marquardt step: solve ``(JᵀJ + λI) δ = Jᵀr`` (ADR 0007 decision 3i).

    ``λ`` is scaled to the Jacobian (``lm_lambda × max diag(JᵀJ)``) so the damping is
    scale-invariant.  Near the along-track wall the Jacobian goes near-singular and the
    plain Newton step diverges; the ``λI`` term regularizes the near-null direction so
    required-Δv grows *smoothly* toward the authority boundary instead of blowing up.
    ``JᵀJ + λI`` is SPD for any ``λ>0``, so this never raises unless the Jacobian is
    exactly zero (no authority at all), which the caller records as non-convergence.
    """
    jtj = jac.T @ jac
    damping = lm_lambda * float(np.max(np.diag(jtj)))
    step = np.linalg.solve(jtj + damping * np.eye(3), jac.T @ residual)
    return np.asarray(step, dtype=np.float64)


def solve_apogee_correction(
    predict: PredictFn,
    target: Target,
    *,
    tol_m: float = 1.0,
    max_iter: int = 8,
    fd_step_m_s: float = 1.0e-3,
    max_step_m_s: float = 2.0,
    lm: bool = False,
    lm_lambda: float = 1.0e-3,
) -> ControlPlan:
    """Solve for the apogee correction Δv (RTN) nulling the interception miss.

    Newton on a finite-difference Jacobian, starting from zero correction.  Converges
    when the EME2000 position residual falls below ``tol_m``; a singular Jacobian or
    exhausting ``max_iter`` returns ``converged=False`` (the authority boundary), not
    an exception.  ``iterations`` is the number of Newton steps taken.

    ``max_step_m_s`` caps each Newton step at the *physical* correction scale (apogee
    corrections are O(0.1–1 m/s)).  This is load-bearing, not just guarding an
    ill-conditioned Jacobian: the 200 km target is an altitude event with free
    time-of-arrival, so a far, high-Δv orbit can re-cross at the same inertial point
    much later — a spurious root the cap keeps Newton from wandering into.  A run that
    genuinely needs more than the cap is past the authority boundary and reads as
    non-converged, which is the honest A3 outcome.

    ``lm`` switches the square Newton step for a Levenberg-Marquardt damped step
    (default off, so A1/A2/the capstone keep their committed Newton path).  A3's
    controllability sweep turns it on to walk the near-singular along-track wall
    smoothly toward the boundary rather than diverging on it (ADR 0007 decision 3).
    """
    target_pos = np.asarray(target.position_m, dtype=np.float64)
    x = np.zeros(3, dtype=np.float64)
    iterations = 0
    converged = False

    while True:
        base_pos = np.asarray(predict(_vec3_of(x)), dtype=np.float64)
        residual = base_pos - target_pos
        if float(np.linalg.norm(residual)) < tol_m:
            converged = True
            break
        if iterations >= max_iter:
            break

        jac = _finite_difference_jacobian(predict, x, base_pos, fd_step_m_s)
        try:
            step = _lm_step(jac, residual, lm_lambda) if lm else _newton_step(jac, residual)
        except np.linalg.LinAlgError:
            break  # singular geometry → out of authority, recorded as non-converged

        step_norm = float(np.linalg.norm(step))
        if step_norm > max_step_m_s:
            step = step * (max_step_m_s / step_norm)
        x = x - step
        iterations += 1

    dv_rtn = _vec3_of(x)
    action = ControlAction(
        node_label="apogee",
        elapsed_s=0.0,
        dv_rtn_m_s=dv_rtn,
        dv_mag_m_s=float(np.linalg.norm(x)),
    )
    return ControlPlan(actions=(action,), converged=converged, iterations=iterations)


def _finite_difference_jacobian_6(
    predict: TwoBurnPredictFn,
    x: NDArray[np.float64],
    base_position: NDArray[np.float64],
    step: float,
) -> NDArray[np.float64]:
    """Forward-difference 3×6 ∂(crossing position)/∂(stacked Δv), reusing ``base_position``."""
    jac = np.empty((3, 6), dtype=np.float64)
    for j in range(6):
        xp = x.copy()
        xp[j] += step
        perturbed = np.asarray(predict(_vec6_of(xp)), dtype=np.float64)
        jac[:, j] = (perturbed - base_position) / step
    return jac


def solve_two_burn_correction(
    predict: TwoBurnPredictFn,
    target: Target,
    *,
    tol_m: float = 1.0,
    max_iter: int = 8,
    fd_step_m_s: float = 1.0e-3,
    max_step_m_s: float = 2.0,
) -> ControlPlan:
    """Solve for the minimum-Δv two-impulse correction nulling the interception miss (ADR 0005).

    Two burns give 6 DOF against the 3 position constraints, so the system is
    underdetermined; each Gauss-Newton step is the **minimum-norm** least-squares step
    (``np.linalg.lstsq``), so starting from zero the solver walks the least-Σ‖Δv‖² path.
    Unlike A1's square ``solve``, ``lstsq`` never raises on a rank-deficient lever — an
    unreachable target simply leaves the residual above ``tol_m`` after ``max_iter``,
    which is recorded as ``converged=False`` (the authority boundary), never thrown.

    ``max_step_m_s`` caps the full 6-vector step at the physical correction scale, the
    same free-ToA spurious-far-root guard as A1.  Returns two ``ControlAction``s
    ("apogee", "midcourse"); the midcourse ``elapsed_s`` is a placeholder the harness
    re-stamps with the real 900 km event time (execution is altitude-event-driven).
    """
    target_pos = np.asarray(target.position_m, dtype=np.float64)
    x = np.zeros(6, dtype=np.float64)
    iterations = 0
    converged = False

    while True:
        base_pos = np.asarray(predict(_vec6_of(x)), dtype=np.float64)
        residual = base_pos - target_pos
        if float(np.linalg.norm(residual)) < tol_m:
            converged = True
            break
        if iterations >= max_iter:
            break

        jac = _finite_difference_jacobian_6(predict, x, base_pos, fd_step_m_s)
        step, *_ = np.linalg.lstsq(jac, residual, rcond=None)

        step_norm = float(np.linalg.norm(step))
        if step_norm > max_step_m_s:
            step = step * (max_step_m_s / step_norm)
        x = x - step
        iterations += 1

    actions = (
        ControlAction("apogee", 0.0, _vec3_of(x[:3]), float(np.linalg.norm(x[:3]))),
        ControlAction("midcourse", 0.0, _vec3_of(x[3:]), float(np.linalg.norm(x[3:]))),
    )
    return ControlPlan(actions=actions, converged=converged, iterations=iterations)
