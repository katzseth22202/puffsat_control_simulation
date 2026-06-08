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


def _vec3_of(a: NDArray[np.float64]) -> Vec3:
    return (float(a[0]), float(a[1]), float(a[2]))


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


def solve_apogee_correction(
    predict: PredictFn,
    target: Target,
    *,
    tol_m: float = 1.0,
    max_iter: int = 8,
    fd_step_m_s: float = 1.0e-3,
    max_step_m_s: float = 2.0,
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
            step = np.linalg.solve(jac, residual)
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
