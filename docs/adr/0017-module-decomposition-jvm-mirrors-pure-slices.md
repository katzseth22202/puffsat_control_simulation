# 0017 — Module decomposition: ownership says when, seams say where; JVM glue mirrors pure slices

Date: 2026-06-11
Status: accepted

## Context

`montecarlo.py` reached 931 lines with **five owning ADRs** (0002, 0003, 0012, 0013,
0014) cited in its CLAUDE.md entry — five reasons to change, the definition of a god
file. The cause is an asymmetry: the pure side has an unwritten rule (every rung slice
gets its own content-named module — `navigation.py`, `nav_feasibility.py`,
`coeff_requirement.py`, `terminal.py`) and is healthy because of it, while the JVM side
had no rule, so every slice's JVM glue accreted into the one harness file. The same
disease shows one level down: `_run_record` carries four default-off behavior knobs
from four different rungs (`control` A1, `toa_window_s` A3, `actuator` B1,
`nav_offset_rtn6` C0). A grilling session (2026-06-11) settled the rules and the split
before C3b adds the next ~150 lines of JVM glue.

## Decision

1. **When to split — ownership, not lines.** A module is a god file when its CLAUDE.md
   entry cites more than one owning ADR. Logic lines are the secondary tripwire:
   suspicious at ~400, act by ~600 even with a single owner. (This is why
   `nav_feasibility.py` at 479/one-owner stands and `montecarlo.py` at 931/five fell.)

2. **Where to cut — on a Seam.** "Seam" is now a CONTEXT.md term: a substitution
   boundary where one side can be exercised without the other (master seam: pure vs
   JVM). A module interface must sit on a seam; a cut that crosses one is wrong even
   if it shrinks files. A single-owner slice that is not a seam is acceptable.

3. **JVM glue mirrors pure slices in `puffsat_sim/runs/`.** One JVM module per slice,
   named identically to its pure core across the master seam
   (`nav_feasibility.py` ↔ `runs/nav_feasibility.py`). The truth-path kernel — the
   regime-switched descent every harness consumes (`descend`, `coast_to_handoff`,
   `propagate_to_interception`, `Crossing`, `apogee_state`) — is a real substitution
   seam and becomes `descent.py`. `montecarlo.py` keeps its name, its
   `make capstone` entry point, and exactly the ADR 0002/0003 predict/execute
   ensemble harness.

4. **The harness surface goes public.** `run_record`, `build_context`, `RunContext`
   lose their underscores: the corrector harness is a seam with first-party consumers
   (`runs/*`), not an implementation detail. The narrow-facade alternative was
   rejected as an abstraction with one implementation.

5. **Function rule.** Suspicious at ~50 logic lines (docstrings excluded) for
   branching logic, ~80 for linear single-path recipes (e.g.
   `run_terminal_feedforward` — splitting those manufactures once-called helpers with
   the whole physics context as parameters). The hard tripwire is **≥3 default-off
   behavior knobs** on one function: bundle them into a frozen value object.
   `run_record`'s four knobs become **`RunVariant`** (control, actuator,
   nav_offset_rtn6, toa_window_s); `run_ensemble`'s public signature is unchanged.
   Rung D knobs (train mode, MPC) get a home instead of a fifth parameter.

6. **`report_controller` moves to pure `control.py`.** The LM-tuned report-grade
   corrector config (`_c0_controller`) is pure and consumed by the C0/C1/C2a reports;
   it never belonged in the JVM file.

## Considered options

- **Hard line cap (300/500) as the primary rule** — rejected: it would force
  make-work splits of cohesive single-owner files (`nav_feasibility.py`,
  `truth_model.py`) while the real signal is multiple owners.
- **Flat three-way split** (`instruments.py` + `reports.py`) — rejected: each new
  file carries 2–3 owning ADRs, violating the when-rule on day one.
- **Kernel-only extraction** — rejected: leaves the god file diagnosed and unfixed,
  and C3b lands on it next.
- **Narrow facade over private harness internals** — rejected: one implementation,
  six first-party callers; altitude the project avoids elsewhere.

## Consequences

- C3b's ZEM-loop glue lands as growth of `runs/terminal.py` (or a sibling), not on
  the harness.
- `truth_model.py` (445 lines, single owner) stands, with a note: its three
  near-clone `report_*_signatures` functions get table-driven the next time any of
  them changes — not before.
- The thresholds live as design-doc §14 practice 5 and a CLAUDE.md convention line;
  CONTEXT.md owns the Seam vocabulary; this ADR owns the shape and the why.
- Mirror-naming means a pure module and its JVM glue share a basename; imports
  disambiguate by package path (`from puffsat_sim.runs.nav_feasibility import …`),
  and the convention is to import functions, not the module, when shadowing could
  confuse.
