"""Perturbations — pure-Python force specs (no JVM dependency).

Each perturbation is a small frozen spec carrying only the parameters its force
needs, plus that force's pure analytic signature.  The Orekit force models are
built from these specs on the JVM side in :mod:`puffsat_sim.forces.build`.
"""

from __future__ import annotations

from puffsat_sim.forces.drag import AtmosphericDrag
from puffsat_sim.forces.geopotential import Geopotential
from puffsat_sim.forces.srp import SolarRadiation
from puffsat_sim.forces.third_body import ThirdBody

Perturbation = Geopotential | ThirdBody | SolarRadiation | AtmosphericDrag

__all__ = [
    "AtmosphericDrag",
    "Geopotential",
    "Perturbation",
    "SolarRadiation",
    "ThirdBody",
]
