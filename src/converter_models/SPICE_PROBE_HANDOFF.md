# SPICE Probe Correctness — Agent Handoff

**Repo:** `OpenMagnetics/MKF`  
**Directory:** `src/converter_models/`  
**Written:** 2026-05-24  
**Status:** Two real bugs remain. Everything else has already been fixed.

---

## Background: §8a.5 probe rules

Three correctness rules for every SPICE extractor
(`simulate_and_extract_topology_waveforms` in each topology `.cpp`):

| Rule | Must be | Must NOT be |
|------|---------|-------------|
| (A) Input voltage | `v(vin_dc)` — the DC supply rail | Tank node (`vab`, `v(sw)`, `pri_trafo_in`) |
| (B) Input current | `i(Vq1_sense)` where `Vq1_sense vin_dc q*_drain 0` | `i(Vdc)` (RC snubber spikes), `vpri_sense#branch` (bipolar tank, averages to 0) |
| (C) Winding voltage | Per-winding differential probe (`Evpri_w`, `Evsec_w_o*`, or bare node where winding bottom is netlist reference) | Lumped tank node or post-rectifier DC bus |

**Golden reference:** `AsymmetricHalfBridge.cpp` (post-fix AHB) — all three rules clean. Use its SW1-mode circuit and extractor as a template.

---

## Current status per topology

### Already clean — do not touch

| Topology | File | Notes |
|----------|------|-------|
| AHB | `AsymmetricHalfBridge.cpp` | Golden reference |
| Boost | `Boost.cpp` | Clean |
| Flyback | `Flyback.cpp` | Clean |
| Buck | `Buck.cpp` | `Vq1_sense vin_dc q1_drain 0` ✓ |
| Isolated Buck | `IsolatedBuck.cpp` | `Vq1_sense vin_dc q1_drain 0` ✓ |
| Isolated Buck-Boost | `IsolatedBuckBoost.cpp` | `Vq1_sense vin_dc q1_drain 0` ✓ |
| LLC | `Llc.cpp` | FB: `Vq1_sense`+`Vq3_sense` on `vin_dc` ✓; HB: `Vq1_sense vin_dc qhi_drain 0` ✓ |
| CLLC | `Cllc.cpp` | `Vq1_sense vin_p qa_drain 0` (`vin_p` = positive DC rail) ✓ |
| DAB | `Dab.cpp` | `Vq1_sense vin_dc1 q1_drain 0` + `Vq3_sense` ✓ |
| Push-Pull | `PushPull.cpp` | `Vct_sense vin_dc center_tap 0` — correct for center-tap architecture ✓ |
| Single-Switch Forward | `SingleSwitchForward.cpp` | `Vq1_sense vin_dc q1_drain 0` ✓; bare `pri_in` winding probe intentional (winding bottom = GND) |
| Two-Switch Forward | `TwoSwitchForward.cpp` | `Vq1_sense vin_dc q1_drain 0` ✓ |
| Active Clamp Forward | `ActiveClampForward.cpp` | `Vq1_sense vin_dc q1_drain 0` ✓ |
| CLLLC | `Clllc.cpp` | `Vq1_sense` on `vin_*` ✓ |
| SRC | `Src.cpp` | Check on next rebuild |
| Weinberg | `Weinberg.cpp` | Check on next rebuild |
| Zeta | `Zeta.cpp` | Check on next rebuild |

> **Note:** CLAUDE.md §8a.5 listed LLC, DAB, CLLC, PushPull, Buck-family, and all forwards as "buggy". That section is **outdated** — the code already has proper `Vq1_sense vin_dc q*_drain 0` everywhere. Only PSFB and PSHB have the remaining issue described below.

---

## Bug 1 — PSFB BEHAVIORAL_PULSE mode: wrong `Vq1_sense` placement

**File:** `PhaseShiftedFullBridge.cpp`  
**Lines:** ~788–791 (BEHAVIORAL_PULSE branch inside `generate_ngspice_circuit`)

### What the code does

```cpp
// BEHAVIORAL_PULSE branch (fast path)
circuit << "Vleg_A leg_a_src 0 PULSE(0 " << Vin << " ...)\n";
circuit << "Vq1_sense leg_a_src mid_A 0\n";   // ← BUG: not on vin_dc
```

### What it should do

```cpp
// Should follow the SW1 branch pattern (lines ~820):
circuit << "Vq1_sense vin_dc qa_drain 0\n";
```

### Why it matters

In BEHAVIORAL_PULSE mode, `Vdc` is decoupled from the bridge (the leg sources are independent PULSE sources). So `i(Vdc) = 0` by construction, and the current `i(Vq1_sense)` measures the leading-leg's current drawn from the behavioral source — which approximates, but does not equal, the true DC bus draw from `vin_dc`. The RC snubbers and body diodes connected to `vin_dc` are not included. This causes the displayed input-current waveform to be an approximation, with possible DC offset errors in operating-point calculations.

The SW1 (switch-model) branch at line ~820 is already correct: `Vq1_sense vin_dc qa_drain 0`.

### How to reproduce

1. Open the PSFB Wizard in the WebFrontend (`http://localhost:5173/wizards`).
2. Choose "I know the design I want", enter any reasonable PSFB params.
3. Click **Get SPICE Code** → download the `.sp` netlist.
4. Inspect the netlist — look for `Vq1_sense`. If it reads `Vq1_sense leg_a_src mid_A 0`, the BEHAVIORAL mode fired.
5. The `.save` line will include `i(Vq1_sense)` but that current flows through the behavioral source terminal, not `vin_dc`.

Alternatively via C++ test:
```bash
cd /home/alf/OpenMagnetics/MKF
./build/tests/OpenMagneticsTest --gtest_filter="*Psfb*spice*"
```
Compare `input_current` waveform peak against `output_power / input_voltage` — a significant discrepancy reveals the bug.

### How to fix

In the BEHAVIORAL_PULSE branch of `generate_ngspice_circuit`:

1. Change the netlist topology so `vin_dc` connects through a real source and `Vq1_sense` sits between `vin_dc` and the leading-leg node — analogous to the SW1 branch. The behavioral voltage sources become dependent on `vin_dc` rather than being independent.
2. OR: accept that BEHAVIORAL mode is an approximation and document it clearly, but ensure the SW1 path is always used when probe accuracy matters (e.g., make SW1 the default for `Get SPICE Code`).

### How to validate when fixed

1. **Unit test:** PSFB simulated-mode waveforms. Assert:
   - `input_voltage` waveform is flat at `Vin` (DC rail constant).
   - `avg(input_current) ≈ output_power / Vin` within 5 % (power balance).
   - Peak `input_current` is non-negative (draws from bus, never pushes).
2. **Waveform shape test:** Input current must be pulsed at twice the switching frequency (both leading-leg high-side switches). If it shows `2·Fs` pulses, the probe is on the right node.
3. **BEHAVIORAL vs SW1 agreement:** run both modes for the same operating point; `avg(input_current)` must agree within 3 %.

---

## Bug 2 — PSHB BEHAVIORAL_PULSE mode: wrong `Vq1_sense` placement

**File:** `PhaseShiftedHalfBridge.cpp`  
**Lines:** ~729–741 (BEHAVIORAL_PULSE branch)

### What the code does

```cpp
// BEHAVIORAL_PULSE branch
circuit << "Vbridge_a leg_a_src 0 PWL(...)\n";
circuit << "Vq1_sense leg_a_src bridge_a 0\n";   // ← BUG: not on vin_dc
```

### What it should do

Match the SW1 branch (lines ~772):
```cpp
circuit << "Vq1_sense vin_dc qa_drain 0\n";
```

### Why it matters

Identical reasoning to PSFB Bug 1. PSHB's half-bridge behavioral model has `vin_dc` decoupled from the leg sources, so `i(Vq1_sense)` from `leg_a_src` → `bridge_a` measures the behavioral leg current, not the actual DC bus draw. The half-bridge uses the DC link capacitor mid-point, so the current imbalance between upper and lower halves makes this approximation worse than in PSFB.

### How to reproduce

Same approach as PSFB:
1. PSHB Wizard → "I know the design" → Get SPICE Code.
2. Inspect netlist for `Vq1_sense leg_a_src bridge_a 0`.
3. Or run: `./build/tests/OpenMagneticsTest --gtest_filter="*Pshb*spice*"` and check power balance on `input_current`.

### How to fix

Same pattern as PSFB fix. The SW1 branch at line ~772 is already correct.

### How to validate when fixed

Same checklist as PSFB, plus:
- Half-bridge specific: `avg(input_current)` must equal approximately `output_power / Vin` (the split-cap mid-point should not introduce a DC offset in the averaged current).
- Check that `input_current` is non-negative (half-bridge draws from bus in each half-cycle).

---

## No missing WASM functions (as of 2026-05-24)

All `simulate_*_ideal_waveforms` functions are present in the built WASM binary AND wired in `src/stores/taskQueue.js`. The following were previously thought missing but are now confirmed present:

| Function in WASM | Topology |
|-----------------|---------|
| `simulate_clllc_ideal_waveforms` | CLLLC |
| `simulate_psfb_ideal_waveforms` | PSFB |
| `simulate_pshb_ideal_waveforms` | PSHB |
| `simulate_src_ideal_waveforms` | SRC |
| `simulate_weinberg_ideal_waveforms` | Weinberg |
| `simulate_zeta_ideal_waveforms` | Zeta |
| `simulate_ahb_ideal_waveforms` | AHB |
| `simulate_cuk_ideal_waveforms` | Cuk |
| `simulate_vienna_ideal_waveforms` | Vienna |

To verify for a future rebuild:
```bash
strings src/assets/js/libMKF.wasm.wasm | grep "^simulate_" | sort
```
All entries from `taskQueue.js`'s `mkf.simulate_*` calls must appear in that list.

---

## WASM rebuild workflow

After any C++ change:

```bash
# 1. Edit source in MKF (this repo)
# 2. Rebuild — MKF source is symlinked into WebLibMKF build tree
ninja -j5 -C /home/alf/OpenMagnetics/WebLibMKF/build/ libMKF.wasm.js

# 3. Deploy to BOTH locations (they must stay in sync)
cp /home/alf/OpenMagnetics/WebLibMKF/build/libMKF.wasm.js  \
   /home/alf/OpenMagnetics/WebFrontend/src/assets/js/
cp /home/alf/OpenMagnetics/WebLibMKF/build/libMKF.wasm.wasm \
   /home/alf/OpenMagnetics/WebFrontend/src/assets/js/
cp /home/alf/OpenMagnetics/WebLibMKF/build/libMKF.wasm.js  \
   /home/alf/OpenMagnetics/WebFrontend/MagneticBuilder/src/assets/js/
cp /home/alf/OpenMagnetics/WebLibMKF/build/libMKF.wasm.wasm \
   /home/alf/OpenMagnetics/WebFrontend/MagneticBuilder/src/assets/js/

# 4. Verify functions present
strings /home/alf/OpenMagnetics/WebFrontend/src/assets/js/libMKF.wasm.wasm \
  | grep "^simulate_" | sort

# 5. Run Playwright (headless ONLY)
cd /home/alf/OpenMagnetics/WebFrontend
npx playwright test --project=site
npx playwright test --project=wizards-light
```

**-fexceptions required:** WASM must be built with `-fexceptions`. Without it, all `catch` blocks are no-ops and C++ exceptions surface as integer pointers in JS (native tests pass, WASM crashes silently). Verify in `WebLibMKF/CMakeLists.txt`.

---

## Test coverage for SPICE paths

SPICE is tested via the `spice: true` capability flag in `WebFrontend/tests/utils/catalog.js`. Every wizard with `spice: true` gets a **Group G** battery test that clicks "Get SPICE Code" and verifies a non-empty `.sp` file is returned. PFC (`spice: false`) and Vienna (`spice: false`) are intentionally excluded (PFC SPICE not yet wired; Vienna SPICE is single-phase emulation only).

To add a SPICE-specific assertion test (e.g., power balance check), add it to the `makeBatterySpec` G-group in `tests/utils/battery.js` or as a standalone spec in `tests/wizards/<key>.spec.js`.

---

## Key contacts / files

| What | Where |
|------|-------|
| Probe rules reference | `CLAUDE.md` §8a.5 (WebFrontend) |
| Golden AHB implementation | `AsymmetricHalfBridge.cpp` SW1 branch |
| PSFB behavioral bug | `PhaseShiftedFullBridge.cpp` ~line 791 |
| PSHB behavioral bug | `PhaseShiftedHalfBridge.cpp` ~line 741 |
| libMKF embind exports | `/home/alf/OpenMagnetics/WebLibMKF/src/libMKF.cpp` |
| taskQueue SPICE calls | `WebFrontend/src/stores/taskQueue.js` (search `simulate_`) |
| Wizard test catalog | `WebFrontend/tests/utils/catalog.js` |
