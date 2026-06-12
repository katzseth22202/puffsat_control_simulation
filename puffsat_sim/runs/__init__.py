"""JVM-side run modules — one per slice, names mirroring the pure cores (ADR 0017).

Each module here is the JVM glue for the same-named pure module across the master
seam: it builds the Orekit runs, feeds the pure core, and formats the slice's
report.  The shared machinery lives in :mod:`puffsat_sim.descent` (truth-path
kernel) and :mod:`puffsat_sim.montecarlo` (the predict/execute harness surface).
"""
