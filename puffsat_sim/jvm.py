"""Orekit JVM lifecycle — the single boot seam.

Importing this module starts the JVM and loads ``orekit-data.zip`` from the
current working directory.  Every module that imports ``org.orekit.*`` must
import this module first, so the JVM is guaranteed booted before any Java class
is referenced.  Python's module cache makes the boot happen exactly once, no
matter how many modules import it.

The one ordering rule in the codebase is therefore::

    import puffsat_sim.jvm   # before any org.orekit import
"""
from __future__ import annotations

from typing import Any, Final

import orekit_jpype

_VM: Final[Any] = orekit_jpype.initVM(vmargs="--enable-native-access=ALL-UNNAMED")

from orekit_jpype.pyhelpers import setup_orekit_curdir

setup_orekit_curdir()  # loads orekit-data.zip from the current working directory
