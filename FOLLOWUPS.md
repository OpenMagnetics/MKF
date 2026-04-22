# MKF / PyMKF follow-up plan

Short working doc covering the follow-ups we accumulated over the last session.
Written when we stopped instead of kept pushing, so the next session can pick
these up without reconstructing context. **Do NOT treat this as permanent
documentation** — delete once the work is done.

## 1. MKF: delete `add_transformer_turn_variants`

- **Location**: `src/advisers/CoreAdviser.cpp:1214` (function body), `:2550` and
  `:2818` (call sites).
- **What it does**: after the power-application filter chain culls to
  `maximumNumberResults`, it loops over the top-5 magnetics and appends a copy
  with `N ← ceil(1.3·N)` at the end of the result vector. Purely cosmetic —
  users see their "top N" plus 5 alternative choices with more turns.
- **Why delete**: after the sweep-minimum refactor in
  `MagneticFilterCoreDcAndSkinLosses` (commit d65ad108), the loss filter already
  explores N upward from the saturation floor and tracks the loss-optimal point
  per core. The +30% variant is strictly worse at the loss-optimal operating
  point, so appending it is misleading UX.
- **Why NOT delete yet**: the DAB regression showed top-10 turns values are
  identical pre- and post-Change 2/3. That could mean either "sweep correctly
  converges on the initial N" or "sweep doesn't run meaningfully". Without a
  case where the sweep visibly moves the result, I don't want to drop the
  fallback. **Decision depends on finding a test case where the sweep moves
  the answer.**
- **Concrete plan**:
  1. Find a test where the saturation-floor N is clearly suboptimal — try a
     larger core (E 55/28/21 at 1 kW 100 kHz DAB, or any EE core where window
     fill allows 2×N_start).
  2. Log inside `MagneticFilterCoreDcAndSkinLosses::evaluate_magnetic` to verify
     `bestNumberTurnsPrimary != currentNumberTurns` on that case — i.e. the
     sweep actually picks a better N.
  3. Once verified, delete the function and both call sites.
  4. Remove the `[adviser][magnetic-adviser]` tag add-transformer-turn-variants
     references in any test, if any.
- **Risk**: low. Removing an appended vector tail only shrinks the result set.
  Max-results changes from `N + min(N,5)` to `N`.

## 2. MKF: port loss-optimal target to `refine_gaps_for_saturation`

- **Location**: `src/advisers/CoreAdviser.cpp:1642`.
- **Current behavior**: iterates gap length to drive B_peak toward
  `0.7·Bsat` — the same "saturation-edge floor as target" pattern the transformer
  branch had before Change 1. For energy-storing inductors the loss-optimal B
  is below the saturation cap (typically 0.1–0.2 T for ferrites), but this
  function pushes gap design to the cap.
- **Why hold back**: changing the inductor target requires care — the flux
  density in an energy-storing inductor has both DC and AC components, and the
  optimization is over `ΔB` while saturation constrains `B_peak`. Simpler to
  leave it targeting 0.7·Bsat for now and revisit after the transformer path
  has a track record.
- **Concrete plan**:
  1. Add an `iGSE`-aware `B_target` calculation for the inductor branch using
     material Steinmetz data. Target the `ΔB` that minimizes
     `P_fe(ΔB) + P_cu(N(gap(B_peak)))` subject to `B_peak ≤ 0.7·Bsat`.
  2. Replace the hardcoded `0.7·Bsat` target in the refine loop with
     `min(B_target_opt, 0.7·Bsat)`.
  3. Regression-check on flyback / buck / boost / PFC adviser tests.
- **Risk**: medium. Inductor design is where the hardcoded 0.7 has historically
  "just worked"; changing the target point will shift gap sizes and may expose
  latent issues elsewhere.

## 3. PyMKF: replace string-sentinel exceptions with real Python exceptions

- **Scope**: 56 occurrences across 9 files in `/home/alf/OpenMagnetics/PyMKF/src/`.
  ```
  advisers.cpp   (3)
  bobbin.cpp     (5)
  core.cpp       (13)
  losses.cpp     (6)
  settings.cpp   (2)
  simulation.cpp (14)
  utils.cpp      (3)
  winding.cpp    (3)
  wire.cpp       (7)
  ```
- **Pattern to remove**: every function wraps its body in
  ```cpp
  try { ... }
  catch (const std::exception &exc) {
      json exception;
      exception["data"] = "Exception: " + std::string{exc.what()};
      return exception;
  }
  ```
  which hides errors from Python callers (they get a dict with a string inside
  instead of a raised exception). Bugs 1 and 2 in `PYOM_BUG_REPORT.md` both
  surface this way.
- **Fix**: remove the `try { ... } catch` scaffolding, let `std::exception`
  propagate. pybind11 auto-converts to `RuntimeError` (or specific Python
  exception types for `std::invalid_argument` → `ValueError`,
  `std::out_of_range` → `IndexError`, etc.).
- **Why not do it in this session**: mechanical, file-scale edit across 56
  sites. Better to do as a focused pass with fresh context and a grep-verify
  step, so the next session can execute the plan cleanly.
- **Concrete plan**:
  1. For each file, `Grep` the pattern `exception\["data"\] = "Exception:`
     and Read enough context to see the wrapping try/catch.
  2. Edit to remove the try block opener + catch block closer, de-indent the
     body by 4 spaces.
  3. Final verify: `grep -c 'exception\["data"\] = "Exception:' src/*.cpp`
     should be `0`.
  4. Rebuild PyMKF from source and run the bug report reproducers. Bugs 1 and
     2 should now raise Python `RuntimeError` with a meaningful message
     instead of returning a dict sentinel.
- **Risk**: low-medium. The catch blocks never did anything useful; removing
  them preserves the happy-path behavior exactly and only changes error-path
  behavior from "return sentinel dict" to "raise Python exception".

## 4. PyMKF: add `py::arg()` names to binding definitions

- **Scope**: every `m.def` call in `PyMKF/src/*.cpp` that currently passes the
  function pointer without `py::arg("name")` markers.
- **Why**: python stubs and `help()` show `calculate_leakage_inductance(arg0,
  arg1, arg2)` and `calculate_maxwell_capacitance_matrix(arg0, arg1)` with no
  parameter names — callers can't discover whether the first arg is a magnetic
  or a coil (Bug 2 from the report).
- **Concrete plan**:
  1. Walk each `m.def(...)` in simulation.cpp, losses.cpp, etc.
  2. Add `py::arg("param_name")` for each positional arg so stubs generate
     `calculate_leakage_inductance(magnetic: json, frequency: float,
     source_winding_index: int)` instead of `(arg0, arg1, arg2)`.
  3. Do this in the same commit as item 3 (exceptions) — both touch every
     pybind11 binding file.
- **Risk**: zero. Parameter names don't change ABI.

## 5. PyMKF `calculate_maxwell_capacitance_matrix`: accept magnetic or coil

- **Related**: Bug 1 in the report.
- **Actual signature** (confirmed in `PyMKF/src/simulation.cpp:291`):
  `calculate_maxwell_capacitance_matrix(coil_json, capacitance_among_windings_json)`
  — takes a **coil**, not a magnetic; takes a **capacitance matrix dict**, not
  an operating point.
- **Bug report misreads the signature** as `(magnetic, operating_point)`. Not a
  code bug in PyMKF — it's a documentation/discoverability bug that item 4
  (py::arg names) fixes.
- **Fix**: item 4 (add `py::arg("coil")`, `py::arg("capacitance_among_windings")`)
  makes the expected arguments self-explanatory. Then either update the bug
  report, or file a new issue against the FEM test harness that calls this
  function wrong.
- **Additional nice-to-have**: raise a specific `InvalidArgumentException`
  when `coil["bobbin"]` is missing, rather than letting nlohmann::json throw a
  bare `key 'bobbin' not found`. That's a one-line guard in the C++
  `StrayCapacitance::calculate_maxwell_capacitance_matrix` entry point.

## 6. MKF: Bug 3 — `calculate_core_losses` silent zero on multi-winding toroid

- **Reproducer**: `toroidal_inductor_round_wire_multilayer_with_insulation.json`
  → `coreLosses = 0.0, magneticFluxDensityAcPeak = 0.0`.
- **Hypothesis from the report**: `calculate_core_losses` may internally skip
  when `magneticFluxDensity` is not pre-populated in the operating point; the
  multi-winding case may produce that pre-population via a different code
  path than the single-winding case.
- **What to check**:
  1. Does `MagneticSimulator::calculate_core_losses(operatingPoint, magnetic)`
     early-return when `magneticFluxDensity` is unset?
  2. Does the multi-winding operating point have `magneticFluxDensity` set at
     the point of the call?
  3. Is the primary winding's excitation the one driving the flux, or does
     the code try the secondary which has no current in this setup?
- **Concrete plan**:
  1. Add a focused Catch2 test that replicates the reproducer and binary-searches
     the code path until the zero falls out.
  2. Fix at the point where the zero originates — likely a missing `process_`
     call on the multi-winding input, or a wrong winding index when
     `get_primary_excitation` returns the first winding that has a voltage
     instead of the first with a current.
- **Risk**: unknown until the root cause is identified.

## 7. MKF: Bug 5 — duplicated harmonic frequencies in `calculate_winding_losses`

- **Reproducer** (from the report):
  ```
  proximityEffectLosses.harmonicFrequencies: [0.0, 100000.0, 300000.0, 300000.0]
  proximityEffectLosses.lossesPerHarmonic: [0.0, 0.28890, 0.06162, 0.00524]
  ```
- **Root cause hypothesis**: the winding-loss code iterates harmonics separately
  per winding and appends into a shared output vector without deduplicating
  on frequency. The two 300 kHz entries are likely from two windings with the
  same harmonic, presented as independent contributions but schema-wise
  collapsed into one list.
- **Concrete plan**:
  1. Find where `harmonicFrequencies` and `lossesPerHarmonic` are populated in
     `WindingProximityEffectLosses` / `WindingSkinEffectLosses`
     (`src/physical_models/Winding*Losses.cpp`).
  2. Decide the correct schema: either
     (a) deduplicate by frequency and sum the loss contributions (simplest,
         but loses per-winding attribution), or
     (b) keep per-winding arrays and surface them as a 2D structure
         `[[winding, freq, loss], ...]`.
  3. Document whichever choice is made in the `WindingLossesOutput` schema.
- **Risk**: low. Output-shape change is additive or dedup-only; callers that
  already sum `lossesPerHarmonic` keep working.

## 8. MKF: Inputs — the proper Bug 4 fix

- The current Bug 4 commit (d68efdd1) keeps mutating `excitation.frequency`
  but does it by a smarter (f·amplitude) rule. That works, but the underlying
  problem — downstream consumers reading `excitation.frequency` when they
  really want `processed.effectiveFrequency` — is still there.
- **Proper fix** (deferred):
  1. Stop mutating `excitation.frequency` entirely in `process_operating_point`.
  2. Ensure `processed.effectiveFrequency` and `processed.acEffectiveFrequency`
     are populated for every excitation.
  3. Audit every consumer of `excitation.get_frequency()` in MKF — decide per
     call-site whether it wants the user-supplied fundamental (keep as-is) or
     the loss-driving frequency (switch to `processed.effectiveFrequency`).
     Expected sites: MagneticFilter\*.cpp, CoreLosses.cpp, WindingLosses.cpp,
     MagnetizingInductance.cpp, SaturationFilter\*.
  4. Keep the new point fix (d68efdd1) as a safety net until the audit is
     done, then remove the whole heuristic.
- **Risk**: medium-high. Touches many consumers.

## 9. FEM: already-failing PFC adviser test

- `Test_Magnetic_Adviser_PFC_Only_Current` (TestMagneticAdviser.cpp:1908) fails
  on HEAD before any of my changes — 94 cores survive until the Saturation
  filter then all get rejected. Not a regression, but flagged here so it
  doesn't look like one.
- **Concrete plan**: separate investigation. Likely the saturation filter for
  inductors-from-CSV operating points has a bug in how it reads
  `processed.peak` vs recomputes from waveform.

---

## Suggested order of execution

1. **Item 3 (PyMKF exceptions)** + **Item 4 (py::arg names)** — one commit each
   or one bundled commit. Mechanical, high-impact, unblocks the FEM team's
   debugging.
2. **Item 5 (capacitance docstring + bobbin guard)** — small, closes Bug 1.
3. **Item 6 (Bug 3 core losses)** — needs investigation but highest-value bug
   fix.
4. **Item 7 (Bug 5 duplicate harmonics)** — small scope, well-defined.
5. **Item 1 (delete add_transformer_turn_variants)** — only after a test case
   proves the sweep is doing useful work.
6. **Item 8 (proper Bug 4)** — only after items 3–5 settle; requires broader
   audit.
7. **Item 2 (inductor loss-optimal target)** — largest scope, lowest urgency.
8. **Item 9 (PFC adviser sat filter)** — can happen anytime, orthogonal.
