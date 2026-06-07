# Perturbations are pure specs, separate from their Orekit force-model builders

We considered making each perturbation a single self-contained class that both
holds its parameters *and* builds its Orekit force model (a `.force_models(env)`
method per force). We rejected that: building a force model imports `org.orekit`,
so constructing a config would boot the JVM and require `orekit-data.zip`, which
would force `tests/test_config.py` and the preset/scenario unit tests to run
inside the JVM.

Instead, a **Perturbation** is a small frozen pure-Python spec
(`forces/geopotential.py`, `forces/third_body.py`, `forces/srp.py`,
`forces/drag.py`); a single JVM-side dispatch (`forces/build.py`,
`to_force_models`) turns each spec into its Orekit **Force Model**. This keeps
`PhysicsConfig` JVM-free and fast to unit-test. The spec→model split is an
internal seam — Orekit runs live in the integration suite, so no mock adapter is
needed. The cost is that one force's wiring is split across two files (its pure
module and the dispatch arm in `build.py`) rather than living in one class.
