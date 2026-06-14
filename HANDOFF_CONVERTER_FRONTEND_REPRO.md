# Handoff: DMC / LLC / CLLC failures against MKF main (b784917a)

**Context.** The OpenMagnetics WebFrontend Playwright suite regressed from 51 →
93 failures after rebuilding the WASM against MKF main `b784917a` +
WebLibMKF `78cdec1`. The frequency-0 cluster was fixed, but three converter-model
clusters now fail. This doc + `tests/TestConverterFrontendRepro.cpp` capture the
exact frontend inputs so the fixes can be developed and verified natively.

Run the repro tests (hidden by default):

```bash
cd build && ninja MKF_tests && ./MKF_tests "[converter-frontend-repro]"
```

All three currently FAIL on main; they assert the desired post-fix behaviour.

---

## 1. DMC — `[MISSING_DATA] 'peakCurrent' is required` (≈22 frontend tests)

**Symptom.** Every DMC wizard flow (down to the trivial layout test) dies at
analytical/simulated with:

> `Exception: [MISSING_DATA] DifferentialModeChoke: 'peakCurrent' is required but
> was not provided. Set it on the converter inputs (no default is substituted).`

**Path.** `propose_dmc_design` → `DifferentialModeChoke::propose_design()` →
`process_operating_points()` → `DifferentialModeChoke.cpp:139`:

```cpp
double peakCurrent = require_input(get_peak_current(), "DifferentialModeChoke", "peakCurrent");  // 20% margin if not specified
```

**Root cause.** The trailing comment `// 20% margin if not specified` is a relic:
DMC used to derive `peakCurrent` from the operating current when it was absent.
Commit `8bb5443c (fix(converters,physics): throw on missing data)` replaced that
derivation with a hard `require_input`. The frontend DMC wizard (a CM/DM choke
sized from impedance/attenuation, not a converter with a defined peak current)
never supplies `peakCurrent` — it provides `operatingCurrent` (3 A here).

**Captured input.** `singlePhaseBalanced`, `operatingCurrent: 3`, no `peakCurrent`.

**Decision needed (MKF owner):**
- **Option A (recommended) — derive in DMC.** Restore the documented default:
  `peakCurrent = operatingCurrent * (1 + rippleFraction)` with an explicit,
  named ripple fraction (the old code implied 20%). This keeps the "no silent
  fallback" spirit if the derivation is documented and the ripple is a real,
  named model parameter rather than a magic constant. The repro test
  `DMC_FrontendRepro_HelpMode_DerivesPeakCurrent` asserts this (no throw).
- **Option B — require it at the boundary.** Keep the throw, and have the
  frontend DMC wizard compute and pass `peakCurrent`. This pushes a sizing
  decision (what ripple?) into the UI layer, which is arguably the wrong place
  for a magnetics model parameter. If chosen, the repro test should instead
  assert the throw, and a matching frontend change is required.

Note line 142 already computes `currentRipple = peakCurrent - operatingCurrent`
and falls back to `operatingCurrent * 0.2` when negative — so the 20% ripple
assumption already exists downstream; Option A just moves it ahead of the throw.

---

## 2. LLC — Nielsen TDA solver non-convergence (≈3–9 frontend tests)

**Symptom.** `expectNoConsoleErrors` catches the stderr line:

> `[LLC] Nielsen TDA solver did not fully converge: residual=4.32014 A
> (Vi=400V, Vo=399.84V, fsw=100000Hz) — analytical waveform is bounded but may be
> imperfect for this op point.`

**Path.** `Llc::process()` → `process_operating_points()` → `solve_steady_state()`;
warning at `Llc.cpp:1017`, gate `i_thresh = max(0.5, 0.02*|Iload_reflected|)`.

**Root cause (already documented in-code).** The comment at `Llc.cpp:1004-1014`
states this is a *known analytical-model accuracy limitation*, not a fabrication
bug: the Newton iterate is bounded and the PtP reference-design tests validate
the waveforms against hardware. "Tightening the TDA solver to drive this residual
down is tracked separately." This handoff supplies a concrete failing design for
that tracked work.

**Captured input.** 400 V → 48 V / 10.42 A, fr=100 kHz, Q=0.4, Ln=5, full bridge,
integrated resonant inductor. Residual ~4.32 A vs. gate floor 0.5 A.

**Decision needed:** This is genuinely MKF solver work (tighten the TDA /
steady-state Newton so the residual falls below the gate for mainstream designs).
Until then the frontend test correctly surfaces a real stderr signal — per the
project's no-allowlist policy, the fix belongs in the solver, not in a
`KNOWN_NOISE` entry. Repro: `LLC_FrontendRepro_NielsenTDA_Converges`.

---

## 3. CLLC — steady-state solver non-convergence (≈8 frontend tests)

**Symptom.**

> `[CLLC] steady-state solver did not fully converge: residual=2.43625 A
> (Vi=400V, Vo=400V)` (and a second corner at residual=1.38611 A)

**Path.** `CllcConverter::process_operating_points({1.0}, 200e-6)` →
steady-state Newton; warning at `Cllc.cpp:1359`.

**Captured input.** 400 V → 400 V / 8.25 A (n=1, unity step), fr=120 kHz, Q=0.4,
Ln=5, Lm=200 µH, forward power flow. Residual ~2.44 A / 1.39 A across corners.

**Decision needed:** Same class as LLC — the CLLC steady-state solver doesn't
converge below the gate for this standard 1:1 / 400→400 design. MKF solver work.
Repro: `CLLC_FrontendRepro_SteadyState_Converges`.

---

## Summary table

| Converter | Frontend tests | Entry point | Failure | Owner / fix |
|-----------|---------------:|-------------|---------|-------------|
| DMC  | ~22 | `DifferentialModeChoke::propose_design` | throws `[MISSING_DATA] peakCurrent` | MKF: derive peakCurrent (Option A) **or** frontend supplies it (Option B) |
| LLC  | ~3–9 | `Llc::process` | Nielsen TDA residual 4.32 A > 0.5 A gate | MKF: tighten steady-state solver |
| CLLC | ~8 | `CllcConverter::process_operating_points` | steady-state residual 2.44 A > 0.5 A gate | MKF: tighten steady-state solver |

DMC is the highest-leverage and the only one with a frontend-vs-MKF policy
question; LLC/CLLC are MKF analytical-solver accuracy work already tracked
in-code. Inputs are embedded verbatim in `tests/TestConverterFrontendRepro.cpp`.
