# Next Session Handoff — Fix the 21 Failing Tests, One at a Time

**Repo state at handoff:** `main` at `4fbd6371` (SPICE-config sweep
complete, `extras` map → `std::optional<double>` migration landed).

**Fast-suite ground truth:**
```
test cases:   2209 |  2156 passed |  21 failed |  23 skipped |  9 failed as expected
assertions: 347783 | 347753 passed |  21 failed
runtime: 3027 s (~50 min) — fast subset only, excludes [adviser][slow][ptp][refdesign]
```

99.0% case pass / 99.99% assertion pass. The remaining **21 failures
are pre-existing** (verified by stashing the last commit and re-running
at HEAD~1 for the affected suites).

The actual list (enumerated from a 50-min `-r compact` run captured to
`/tmp/mkf_full.log` at handoff time) is **further down in §"The actual
21 failing test cases (ground-truth list)"** — the speculative
categories A-E I wrote first turned out to be misleading (those were
[adviser]/[ptp]-tagged failures that DON'T appear in the fast suite).
**The fast-suite failures are concentrated in core-losses, winding-
losses, temperature, and fracpole — not the adviser pattern I
expected.** Read the ground-truth list before anything else.

---

## ⚠️ Working protocol for this session

**Run ONLY the failing tests until each is individually green.** Do NOT
run the smoke suite, full fast suite, or any broad filter while fixing
these — every wide run wastes ~50 min that should be spent on the
narrow fix. Use the `Test_*` name as the Catch2 filter:

```bash
cd /home/alf/OpenMagnetics/MKF/build
./MKF_tests "Test_Specific_Failing_Test_Name"
```

For each failing test:
1. Reproduce it in isolation with the exact `Test_*` filter.
2. Diagnose root cause (read the test, read the implementation it
   exercises, look at the actual vs expected value).
3. Fix the underlying issue (per CLAUDE.md §"Never route around broken
   things" — do NOT tag with `[!shouldfail]`, do NOT relax the
   assertion to make it pass, do NOT skip).
4. Re-run JUST that test. Confirm green.
5. Commit + push the fix.
6. Move to the next test.

Only after **all 21 are green individually** should you run the broader
suite to confirm no new regressions, in this order:
1. Each affected topology's `[xxx-topology]` tag.
2. `~[adviser]~[slow]~[ptp]~[refdesign]` (the fast suite — but only at
   the very end).

---

## The actual 21 failing test cases (ground-truth list)

Already enumerated during this session — 30 unique `file:line` assertion
sites in 10 source locations, distributed across 21 `TEST_CASE` runs.
Listed below in priority order (cheap-to-debug first, hardest last).

### 1. `Test_Vienna_missing_operating_points_fails_run_checks`
- **File:** `tests/TestVienna.cpp:619`
- **Symptom:** `unexpected exception: 'Vienna: 'operatingPoints' must be a non-empty array'`
- **Cause:** `Vienna::Vienna(const json&)` calls `validate_vienna_spec_shape(j)` which throws on empty `operatingPoints` BEFORE the test calls `run_checks(true)`. Test expects the ctor to succeed and `run_checks` to throw — but the ctor pre-empts it.
- **Fix options:** (a) move the empty-operatingPoints check from `validate_vienna_spec_shape` into `run_checks`; (b) rewrite the test to expect the ctor throw. Option (a) is cleaner — keep `validate_vienna_spec_shape` for malformed JSON (wrong types) and let `run_checks` handle empty business-logic states.

### 2. `Test_Fracpole_Evaluate_Debug`
- **File:** `tests/TestCircuitSimulatorInterface.cpp:2857`
- **Symptom:** `R_high > R_low * 10 for: 0.00000000016167893 > 0.00147055906061926` — i.e. R_high (~1.6e-10) is FAR smaller than R_low (~1.5e-3), inverting the expected stage ordering.
- **Cause:** The 12-stage fracpole network's R values are being computed in the wrong order or with a wrong scaling. The test name has `[debug]` — possibly a half-finished investigation.
- **Investigation seed:** read `emit_fracpole_winding_spice` in `src/processors/CircuitSimulatorInterface.cpp` around the R-stage generation. Compare the ordering against `Test_Fracpole_Evaluate_Smoke` if that one passes.

### 3. Three `TestTemperature` failures
- **`Temperature: concentric_flyback_rectangular_column`** — `tests/TestTemperature.cpp:3664` — `tempsByType.at("core") <= 605.73 for: 676.79378...` — core temp 671 K vs gate 605.73 K (≈12 % over).
- **`Temperature: concentric_transformer_contiguous_rectangular_wire`** — `tests/TestTemperature.cpp:3777` — `tempsByType.at("core") >= 283.82 for: 270.933...` — core temp 271 vs gate 283.82 (≈5 % under).
- **`Temperature: BuckInductor T134_77_27 from MAS file`** — `tests/TestTemperature.cpp:3889` — `'[json.exception.parse_error.101] parse error at line 1, column 1: attempting to parse an empty input'` — the MAS file load is returning empty content.
- **Cause:** Two are thermal-model drift (the bounds are tight ±5-12 %) — likely a thermal-network parameter shifted in some commit. The third is a file-loading bug — `T134_77_27` MAS file path or load function returns empty input.
- **Investigation seed:** Bisect when `Test_Buck_Inductor*` started failing for the empty-input case (`git log -p tests/TestTemperature.cpp | head -200`). For the two ± bound failures, check if the thermal-network resistance values in `src/physical_models/Temperature.cpp` have moved recently.

### 4. `Test_Kool_Mu_Ultra_60`
- **File:** `tests/TestCoreLosses.cpp:2712`
- **Symptom:** `Catch::Matchers::WithinAbs(...) for: 0.23200000000000001 is within 0.0232 of 0.26913806865459983` — actual `B_offset` is 0.269 vs expected 0.232 ± 0.023.
- **Cause:** Kool Mu Ultra 60 material parameters drifted, or the new test data was authored against a different model coefficient. ~16 % offset error.

### 5. `Test_Core_Losses_Hoganas`
- **File:** `tests/TestCoreLosses.cpp:2855`
- **Symptom:** `'[json.exception.parse_error.101] parse error at line 1, column 1: attempting to parse an empty input'` — same empty-MAS-file pattern as the Buck temperature test above. Likely the same root cause.

### 6. Mass core-losses formula blowup — `TestCoreLosses.cpp:700`
- **File:** `tests/TestCoreLosses.cpp:700` (assertion site inside the helper `test_core_losses_*`; called from MANY TEST_CASEs starting at line 1742).
- **Symptom:** **Multiple failures with catastrophic numerical blowup:**
  - `138214 within 79389 of 49618` — 2.8× off (recoverable)
  - `343943 within 99200 of 62000` — 5.5× off
  - `1.4e+27` value
  - `6.6e+42` value
  - `1.8e+33` value
  - `2.8e+233` value
  - `4.7e+97` value
  - `6.7e+115` value
  - **Three `inf` results**
- **Cause:** One or more core-loss formulas overflow at certain operating-point ranges. The pattern (clusters around specific value scales) suggests a `pow(B, n)` or `exp(...)` blowing up when an intermediate goes negative or exceeds a representable range. Inf appearances point at division-by-zero somewhere.
- **Investigation seed:** Grep the helper to identify which model is computing each value: `awk 'NR>=600 && NR<=750' tests/TestCoreLosses.cpp`. Then check the model classes in `src/physical_models/CoreLosses*.cpp` — almost certainly one of Steinmetz / iGSE / Magnetics / Micrometals has a bad input-range guard.
- **Triage:** Run the suite filtered to ONE model at a time (each `TEST_CASE` uses a specific `CoreLossesModels` enum) to isolate which model is bad: `MKF_tests "[steinmetz-core-losses-model]"`, `"[igse-core-losses-model]"`, etc.

### 7. Mass winding-losses tolerance drift — `TestWindingLosses.cpp:61` and `:98`
- **File:** `tests/TestWindingLosses.cpp:61` (2 fails) and `:98` (6+ fails). Both are inside helper functions called from MANY TEST_CASEs starting at line 108.
- **Symptom:** Mix of mild and catastrophic drift:
  - Line 61: `0.000442 within 0.000110 of 0.000573` — ~30 % off (mild)
  - Line 61: `0.01035 within 0.0026 of 0.01295` — ~25 % off
  - Line 98: `0.34495 within 0.086 of 0.25370` — ~36 %
  - Line 98: `1421 within 426 of 9.74` — ~146× off (catastrophic)
  - Line 98: `1438 within 432 of 9.30` — ~155×
  - Line 98: `219 within 55 of 164` — ~33 %
  - Line 98: `289 within 87 of 7.98` — ~36×
  - Line 98: `513 within 154 of 93.5` — ~5.5×
- **Cause:** Two distinct issues — the 25-36 % drifts are tolerance/coefficient shifts; the 36-155× catastrophic results are a wrong-formula regression for a subset of inputs (probably a unit confusion or missing reflection factor for some winding geometries).
- **Investigation seed:** `tests/TestWindingLosses.cpp:60-100` to see which helper and which inputs trigger each. Then look at `src/physical_models/WindingOhmicLosses*.cpp` for recent changes.
- **Note:** Lines 61 and 98 are in TEST_CASEs `Test_Winding_Losses_One_Turn_Round_Sinusoidal_Stacked` (line 108) and tests around line 195+ — the smoke-tagged ones explicitly tagged `[!mayfail]` (which is a milder XFAIL, will fail-as-expected, but NOT the `[!shouldfail]` we should leave alone). Confirm the failing ones don't carry `[!mayfail]` — if they do, they're already accepted-as-broken and shouldn't be in the 21 count.

---

## Quick-repro per failure (no 50-min suite needed)

```bash
cd /home/alf/OpenMagnetics/MKF/build
./MKF_tests "Vienna: missing operating points fails run_checks"
./MKF_tests "Test_Fracpole_Evaluate_Debug"
./MKF_tests "Temperature: concentric_flyback_rectangular_column"
./MKF_tests "Temperature: concentric_transformer_contiguous_rectangular_wire"
./MKF_tests "Temperature: BuckInductor T134_77_27 from MAS file"
./MKF_tests "Test_Kool_Mu_Ultra_60"
./MKF_tests "Test_Core_Losses_Hoganas"
# For the mass TestCoreLosses.cpp:700 failures, run each model tag:
./MKF_tests "[steinmetz-core-losses-model]"
./MKF_tests "[igse-core-losses-model]"
./MKF_tests "[magnetics-core-losses-model]"
./MKF_tests "[micrometals-core-losses-model]"
# For winding-losses, both lines come from many TEST_CASEs — run the tag:
./MKF_tests "[winding-losses]"
```

---

## (Archived) Speculative categories from earlier in this session

The five categories below were what I predicted the failures would look
like based on targeted suite runs (`[llc-topology]`, `[psfb-topology]`,
…) without proper exclusion of `[adviser]` / `[ptp]` / `[refdesign]`.
**They are NOT the fast-suite failures** — those are the ones in the
ground-truth list above. The categories here are still informative for
the `[adviser]` / `[ptp]` / `[refdesign]` sweeps once the fast-suite
ones are clean (a separate workstream). Skim, don't act on them yet.

## 5 failure categories observed during the SPICE sweep (slow suites)

### Category A — `MagneticAdviser::get_advised_magnetic_from_converter` returns 0 results

**Symptom:**
```
REQUIRE( results.size() >= 1 )
with expansion:
  0 >= 1
```

**Test files seen failing this way:**
- `tests/TestMagneticAdviserConverterVariants.cpp` — many topologies
  (PSFB, SSF/TSF/ACF forward family, SRC, IsolatedBuck, SEPIC, …).
- `tests/TestMagneticAdviser.cpp:3649` (Push-Pull), `:3681`
  (IsolatedBuck), `:4179`, `:5052` — adviser integration tests.

**Hypothesis (root cause):**
The pattern is consistent across topologies, so it's a single bug in
`MagneticAdviser::get_advised_magnetic_from_converter` or its filter
chain — not per-topology. Investigation seeds:
- One log line we kept seeing: `"[MagneticAdviser] No magnetics found
  with toroids enabled. Retrying without toroids..."` — the adviser is
  going through both passes and ending up with `results.empty()`.
- The `Test_CoreAdviser_LLC_From_Frontend_Inputs` test in
  `TestCoreAdviser.cpp:2404` documents the same shape: "After Saturation
  filter, all 163 cores are eliminated (163 → 0)" — strongly suggests
  the saturation filter is rejecting every core for the integrated-from-
  converter operating point. Likely the converter→adviser bridge
  computes saturation differently than the adviser standalone path
  expects (units? RMS vs peak? B vs flux?).

**Where to start:** `src/advisers/MagneticAdviser.cpp ::
get_advised_magnetic_from_converter` and the `MagneticFilterSaturation`
chain. Run one specific test in this category first, instrument the
filter to print how many cores survive each stage, find the stage that
eats everything.

### Category B — Vienna ctor throws BEFORE run_checks

**Test:** `tests/TestVienna.cpp:619` —
`Test_Vienna_missing_operating_points_fails_run_checks`.

**Symptom:**
```
unexpected exception with message:
  Vienna: 'operatingPoints' must be a non-empty array
```

**Cause:** `Vienna::Vienna(const json&)` calls
`validate_vienna_spec_shape(j)` which throws on empty `operatingPoints`
BEFORE the test can call `run_checks()`. The test asserts that
`run_checks(true)` throws, but the ctor already threw.

**Fix options:**
1. Move the empty-operatingPoints check out of `validate_vienna_spec_shape`
   into `run_checks` (where the test expects it).
2. Update the test to expect `Vienna(j)` to throw, not `run_checks`.

Option 1 is more consistent with the file's documented run_checks
contract; option 2 is the smaller diff. Discuss with the user.

### Category C — LLC `bad_optional_access` (excluded from fast suite by `[adviser]`)

**Test:** `tests/TestCoreAdviser.cpp:2404` —
`Test_CoreAdviser_LLC_From_Frontend_Inputs`.

**Symptom:**
```
unexpected exception with message: 'bad optional access'
```

**Tag:** `[adviser][core-adviser][standard-cores][llc-topology][debug]`
— excluded from the fast suite I ran. Shows up when you run any LLC
suite without the `~[adviser]` exclusion.

**Cause:** Standard `std::bad_optional_access` — somewhere in the
adviser code path triggered by this scenario, `std::optional::value()`
is called on a `nullopt`. Run with gdb or `--abort` to get a stack.
Likely related to Category A (the saturation filter eliminating all
cores leaves the adviser holding an empty container it then tries to
`.value()` on).

### Category D — LLC PtP SIGSEGV (excluded from fast suite by `[ptp]`)

**Test:** `tests/TestLlcReferenceDesignsPtp.cpp:144`.

**Symptom:** `SIGSEGV - Segmentation violation signal` during one of
the LLC reference-design analytical-vs-SPICE comparisons.

**Tag:** `[ptp][refdesign]` — excluded from fast suite.

**Hypothesis:** Likely a null pointer or out-of-bounds in the Nielsen
TDA convergence path for one specific design point (the LLC `[llc-topology]`
suite also shows a non-SIGSEGV LLC fail —
`"[LLC] Nielsen TDA solver did not fully converge: residual=1.6854 A"`
at `Vi=200 V Vo=200.16 V fsw=100 kHz` — the same operating point may
crash in the PtP path because it relies on a converged solution).

Run under `valgrind` or `-fsanitize=address` to localise.

### Category E — ACF PtP NRMSE drift (excluded from fast suite by `[ptp]`)

**Test:** `tests/TestActiveClampForwardReferenceDesignsPtp.cpp:292`.

**Symptom:**
```
CHECK( nrmse < s.tol_nrmse )
with expansion:
  0.21315870588131403 < 0.14999999999999999
```

i.e. one ACF reference design shows ~21% NRMSE between analytical and
SPICE, vs a 15% gate.

**Pre-existing** (confirmed during the SPICE sweep — stashed my changes,
re-ran at HEAD~1, same 21% NRMSE). Predates this session.

**Hypothesis:** Either the analytical model is wrong for the
specific design point, or the SPICE netlist doesn't faithfully model
something the analytical model expects (clamp-cap settling time,
dead-time positioning, reset asymmetry). Compare the analytical
waveform shape to the SPICE waveform with `dump_waveforms_csv` (the
test already calls it on line 294-295) and inspect side-by-side.

---

## What you should NOT touch (or touch with extreme care)

- The `SpiceSimulationConfig` struct / `extras` migration just landed
  (`4fbd6371`). Don't revert it. If you find a test that breaks because
  of it, the migration is byte-equivalent at the netlist level (verified
  per-topology at commit time) — the symptom you're seeing is almost
  certainly pre-existing and was previously hidden by other failures.
- Vienna Phase 3 (`3ad41c39` … `f2fb21cc`) — full feature, tests under
  `[phase3]`, all green.
- The wizard fixes (PSFB/PSHB hang, PSHB secondary `Vin/(4·Vout)`, CSV
  importer hardening) — all landed, all green for their unit tests.

---

## Commit + push protocol per fix

Each test fix gets its own commit. Commit-message shape:

```
fix(<area>): <what> for <test>

<root cause analysis: where the bug lives>
<fix: what changed and why>

After this commit:
  Test_X_Y: 1/1 pass (was bad_optional_access)
  [xxx-topology]: <pass/total> (same as HEAD on everything else)
```

Push directly to `main` (this series' style — small, focused, verified
commits). Use the SSH key per `~/.claude/CLAUDE.md`:
```bash
git -c core.sshCommand="/usr/bin/ssh -i ~/.ssh/hephaestus_om" push origin main
```

---

## When all 21 are green

Run the fast suite one final time to confirm:
```bash
./MKF_tests "~[adviser]~[slow]~[ptp]~[refdesign]" -r compact 2>&1 | tail -5
```

Expected: `test cases: 2209 | 2209 passed | 23 skipped | 9 failed as
expected`. The "failed as expected" 9 are `[!shouldfail]`-tagged tests
that pass-by-being-expected-to-fail; leave those alone (per
`feedback_shouldfail_tag.md`).

Then run the slower suites at your discretion:
```bash
./MKF_tests "[adviser]" -r compact 2>&1 | tail -5
./MKF_tests "[ptp]"     -r compact 2>&1 | tail -5
./MKF_tests "[refdesign]" -r compact 2>&1 | tail -5
```

These will surface additional failures (mostly more Category A adviser
issues). Fix them in the same one-at-a-time discipline.

---

## Reference: this session's SPICE-config sweep commits (FYI, don't revert)

```
4fbd6371  refactor(spice-config): collapse extras map → std::optional<double>
bc20eed1  mop up ACF + Vienna + PSFB + AHB
db742ee4  topology-specific extras map (no-fallback throws)
cda7754b  ssf+tsf+ibb full
7714e448  pfc
1b76d55b  cmc+dmc
6f795fef  isolated-buck  ← contains the std::defaultfloat lesson (read this)
dc0a797c  src
5a58b6da  ssf+tsf+ibb partial (PULSE)
dd5dcb4d  ahb
a37e262e  psfb+pshb registry split + full
69586e2e  psfb+pshb partial (PULSE)
ef6ecf44  cllc+clllc
ab61883f  llc
5e5e4910  cllc+boost test gates (physics fix sync)
```

`6f795fef`'s commit message has the **most-useful debugging lesson** of
the series — read it before doing any more SPICE-stream-format work.

End of handoff.
