# Converter Models Review Plan — Bring All Topologies to DAB-Quality

**Audience**: Implementing AI agents (this plan is verbose by design — every step
includes file paths, line numbers, exact code snippets to paste, and explicit
acceptance criteria. Do NOT skip any phase. Do NOT change scope mid-phase.)

**Reference standard**: `CONVERTER_MODELS_GOLDEN_GUIDE.md` — read it once before
starting. Whenever in doubt, the golden guide wins.

**Reference implementation**: `Dab.cpp` / `Dab.h` / `tests/TestTopologyDab.cpp` —
when this plan says "copy DAB pattern at lines XXX", open the file and copy the
exact lines.

---

## Table of Contents

1. [Audit results — current state of every model](#1-audit-results)
2. [Snippet library — copy-paste ready code from DAB](#2-snippet-library)
3. [Per-model remediation plans](#3-per-model-plans)
   - 3.A [CLLC rewrite (tier-3 → tier-1)](#3a-cllc) — execute existing CLLC_REWRITE_PLAN.md
   - 3.B [PSHB topology-netlist mismatch fix](#3b-pshb)
   - 3.C [ActiveClampForward ZVS implementation](#3c-acf)
   - 3.D [PowerFactorCorrection real-circuit SPICE](#3d-pfc)
   - 3.E [Boost / Buck quality baseline](#3e-boost-buck)
   - 3.F [Flyback polish](#3f-flyback)
   - 3.G [PushPull polish](#3g-pushpull)
   - 3.H [IsolatedBuck / IsolatedBuckBoost polish](#3h-iso-buck)
   - 3.I [SingleSwitchForward / TwoSwitchForward polish](#3i-forward)
   - 3.J [LLC SPICE harness upgrade](#3j-llc)
   - 3.K [DifferentialModeChoke tests](#3k-dmc)
   - 3.L [CurrentTransformer architecture](#3l-ct)
   - 3.M [CommonModeChoke docstring polish](#3m-cmc)
4. [Recommended sequencing](#4-sequencing)
5. [Acceptance & verification](#5-verification)

---

## 1. Audit Results

### 1.1 Tier classification (verified April 2026)

| Model | Tier | LOC | Tests | Verdict | Days |
|-------|------|-----|-------|---------|------|
| `Dab` | **1** ✓ | 1910 | 21 | Gold reference. Do not touch. | 0 |
| `PhaseShiftedFullBridge` | **1** ✓ | 1440 | 18 | Quietly upgraded since memory was last written. PtP NRMSE < 0.18, all golden-guide checks satisfied. | 0 |
| `Llc` | **2** | 1855 | 21 | Excellent solver, primary switches use behavioral `PULSE` instead of `SW1`+snubber. No accuracy disclaimer. No ASCII art. | 3–5 |
| `Flyback` | **2** | 1363 | 23 | Strong tests + correct mode handling. Missing docstring, snubber, GEAR, R_load guard, last_* diagnostics. | 2 |
| `PushPull` | **2** | ~2000 | 8 | Good solver, multi-secondary K-matrix. Missing docstring, snubbers on **switches** (only on diodes), GEAR, Lo population. | 2–3 |
| `IsolatedBuck` | **2** | 788 | 5 | Multi-output test passes. Missing docstring, snubber, GEAR, R_load guard, diagnostics. Method name uses camelCase. | 3–4 |
| `IsolatedBuckBoost` | **2** | 685 | 3 | Same gaps as IsolatedBuck + debug flag `keepTempFiles=true` left in production code. | 3–4 |
| `SingleSwitchForward` | **2** | 754 | 3 | Demag winding correctly modeled in SPICE. Snubber resistor present but **1 MΩ** (should be 1 kΩ + 1 nF cap). No GEAR, no diagnostics. | 4–5 |
| `TwoSwitchForward` | **2** | 784 | 4 | Same gaps as SingleSwitchForward. | 4–5 |
| `Boost` | **3** | 598 | 3 | No docstring. No snubber. No GEAR. No R_load guard. No DCM flag. Only 3 tests. | 3–5 |
| `Buck` | **3** | 591 | 5 | Same as Boost. | 3–5 |
| `Cllc` | **3** | 1283 | 16 | Uses FHA (frequency-domain), **not** time-domain Nielsen solver. No ASCII art, no equations, no diagnostics, no `extra*Waveforms`, no `get_extra_components_inputs()`, no PtP test. Throws on SPICE failure. **Execute CLLC_REWRITE_PLAN.md** instead of patching. | 10–15 |
| `PhaseShiftedHalfBridge` | **3** | 1341 | 16 | **CRITICAL BUG**: SPICE netlist never applies phase shift to the switches; netlist always simulates 50% duty fixed-leg. The PtP test passes by coincidence. Topology label still mismatches netlist (single-leg + split-cap, not 2-leg PSHB). Either fix the netlist OR execute AHB_PLAN.md split. | 5–7 |
| `ActiveClampForward` | **3** | 845 | 3 | **The defining feature (ZVS via resonant clamp transition) is COMPLETELY ABSENT** — clamp cap not in SPICE, no resonant LC analysis, no ZVS margin, no ZVS boundary test. | 8–10 |
| `PowerFactorCorrection` | **3** | unk | 0 | Behavioral B-source SPICE (no real switch/diode/snubber). No diagnostic fields for PF/THD/peak-current. No test file. Totem-pole not supported. | 10–15 |
| `CommonModeChoke` | **2 (passive)** | unk | 15 | Strong tests + simulation. Missing ASCII art and references. | 1 |
| `DifferentialModeChoke` | **2 (passive)** | unk | 0 | Functional but **zero tests**. Missing ASCII + refs. | 2–3 |
| `CurrentTransformer` | **3 (passive)** | unk | 1 | Does **not** inherit `Topology` base class. No docstring. Only 1 smoke test. Architectural deviation. | 2–3 |

**Total effort to bring everything to tier-1**: 65–95 working days.

### 1.2 Surprises vs. the previous memory snapshot

- **PSFB has been quietly upgraded to tier-1** — the memory file
  `project_converter_model_quality_tiers.md` is stale; PSFB now has correct
  Im_peak (uses `D_cmd`, not `Deff`), correct ILo 3-slope geometry, K-matrix,
  snubbers, GEAR, IC pre-charge, and a passing PtP NRMSE test (< 0.18).
- **PSHB still has a deeper bug than memory said** — not just topology
  conflation. The SPICE netlist completely ignores the phase-shift parameter
  and simulates a 50%-duty fixed leg regardless of user input. Tests pass
  because both sides happen to converge to similar zero-shift waveforms.
- **CLLC uses FHA, not Nielsen TDA** — this is a fundamental architectural
  gap that simple patching cannot fix. Execute `CLLC_REWRITE_PLAN.md`.
- **ActiveClampForward has no resonant clamp model at all** — this is the
  killer feature of active-clamp forward (ZVS via Cclamp+Lleak resonance).
  Without it, the model is no better than a regular forward converter and
  cannot answer the question "do I get ZVS at this load?".

---

## 2. Snippet Library — Copy-Paste-Ready Code from DAB

Most fixes in this plan reduce to "copy this DAB snippet and adapt it". This
section gives the exact snippets so implementing agents do not need to
re-derive them.

### 2.1 Header docstring template (golden guide §3)

Every `<Topology>.h` file MUST have this BEFORE `class <Topology> :`:

```cpp
/**
 * @brief <One-line topology name and family> Converter
 *
 * Inherits converter-level parameters from MAS::<TopologyMasClass> and
 * the Topology interface, mirroring the DAB pattern.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 *   <ASCII art schematic — copy DAB's pattern at Dab.h lines 22-32>
 *   Use boxes for switches, lines for wires, +/- for polarity, "[L]" or
 *   "[Lr]" for inductors, "[Cb]" for caps, "[D1]" for diodes, "T1" for
 *   transformers with "Np : Ns" for turns ratio.
 *
 * =====================================================================
 * MODULATION / OPERATION CONVENTION
 * =====================================================================
 *
 * Map every MAS schema field to a variable used in the equations below.
 * Specify radians vs. degrees explicitly. If the topology has multiple
 * modes (CCM/DCM/QRM, SPS/EPS/DPS/TPS, mode 1..6), document each.
 *
 * References:
 *   [1] <First author> et al. — "<Title>", <venue> <year>.
 *   [2] TI <SLUAxxx> — "<App-note title>".
 *   [3] Erickson-Maksimović 3rd ed., chapter <N>.
 *   <ALWAYS cite at least 3 sources>
 *
 * =====================================================================
 * KEY EQUATIONS
 * =====================================================================
 *
 * Turns ratio:
 *   <equation>          [Source: <ref>, eq. <N>]
 *
 * Conversion ratio (M = Vo/Vin):
 *   <equation>          [Source: <ref>, eq. <N>]
 *
 * Power transfer:
 *   <equation>          [Source: <ref>, eq. <N>]
 *
 * Inductor sizing (CCM ripple):
 *   L >= <expression>   [Source: <ref>]
 *
 * Output capacitor sizing:
 *   C >= <expression>
 *
 * ZVS condition (if applicable):
 *   <expression>        [Source: <ref>]
 *
 * =====================================================================
 * ACCURACY DISCLAIMER
 * =====================================================================
 *
 * The SPICE model omits:
 *   - MOSFET output capacitance Coss(Vds) variation (uses constant Coss).
 *   - Body-diode reverse recovery Trr.
 *   - Gate-driver propagation delay and dead-time jitter.
 *   - PCB parasitic inductance / capacitance.
 *   - Magnetic-core saturation (assumes linear B-H).
 *
 * Expected NRMSE (analytical vs. SPICE primary current) ≤ 0.15 across the
 * design envelope (see TestTopology<Topo>_PtP_AnalyticalVsNgspice).
 *
 * =====================================================================
 * DISAMBIGUATION (only if topology name is ambiguous)
 * =====================================================================
 *
 * <For Llc>: distinguishes from SRC (no Lm), from Cllc (no Cm), from PSFB.
 * <For IsolatedBuck>: also called Flybuck. Distinct from Flyback.
 * <For ActiveClampForward>: distinct from Forward + RCD reset.
 * <For Pshb>: this is the SINGLE-LEG + SPLIT-CAP variant; not the
 *             true 2-leg Pinheiro-Barbi PSHB.
 */
class <Topology> : public MAS::<MasClass>, public Topology {
```

### 2.2 Class private fields (golden guide §2.1–2.5)

```cpp
private:
    // 2.1 Simulation tuning
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 3;

    // 2.2 Computed design values
    double computedSeriesInductance = 0;       // L_series, L_leakage, or L_resonant
    double computedMagnetizingInductance = 0;  // Lm
    double computedOutputInductance = 0;       // Lo (if PWM bridge / forward)
    double computedDeadTime = 200e-9;
    double computedEffectiveDutyCycle = 0;     // PWM bridges
    double computedDiodeVoltageDrop = 0.6;     // if rectifier present

    // 2.3 Schema-less device parameters
    double mosfetCoss = 200e-12;

    // 2.4 Per-OP diagnostics — populated by process_operating_point_for_input_voltage
    mutable double lastDutyCycleLoss = 0.0;
    mutable double lastEffectiveDutyCycle = 0.0;
    mutable double lastZvsMarginLagging = 0.0;
    mutable double lastZvsLoadThreshold = 0.0;
    mutable double lastResonantTransitionTime = 0.0;
    mutable double lastPrimaryPeakCurrent = 0.0;
    // Add topology-specific: lastMode (LLC), lastModulationType (DAB),
    //   lastPowerFactor (PFC), lastCurrentThd (PFC), lastDcmThreshold (Buck/Boost), etc.

    // 2.5 Extra-component waveforms — cleared at start of process_operating_points()
    mutable std::vector<Waveform> extraLrVoltageWaveforms;
    mutable std::vector<Waveform> extraLrCurrentWaveforms;
    // Add per-component pairs: extraLoVoltageWaveforms, extraLoCurrentWaveforms,
    //   extraCbVoltageWaveforms, etc.
```

### 2.3 SPICE switch + snubber + sense source pattern (golden guide §5)

For each active switch in the topology, emit:

```cpp
ss << "* === Switch " << switchName << " ===\n";
ss << "S" << switchName << " " << drainNode << " " << sourceNode
   << " " << gateNode << " 0 SW1\n";
ss << "Rsnub_" << switchName << " " << drainNode << " "
   << "snub_" << switchName << "_n 1k\n";
ss << "Csnub_" << switchName << " snub_" << switchName << "_n "
   << sourceNode << " 1n\n";
```

The single switch model (emit ONCE in the netlist):

```cpp
ss << ".model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg\n";
ss << ".model DIDEAL D(IS=1e-12 RS=0.05)\n";
```

The mandatory solver options (emit ONCE):

```cpp
ss << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
ss << ".options METHOD=GEAR TRTOL=7\n";
```

R_load divide-by-zero guard:

```cpp
double rLoad = std::max(outputVoltage / outputCurrent, 1e-3);
```

Output cap pre-charge:

```cpp
ss << ".ic v(vout_o" << (i+1) << ")=" << outputVoltage_i << "\n";
```

Multi-secondary K-matrix (every pair):

```cpp
// Primary ↔ each secondary
for (size_t i = 0; i < numSec; ++i) {
    ss << "K_pri_sec" << (i+1) << " Lpri Lsec_o" << (i+1) << " 0.999\n";
}
// Every secondary ↔ every other secondary
for (size_t i = 0; i < numSec; ++i) {
    for (size_t j = i+1; j < numSec; ++j) {
        ss << "K_sec" << (i+1) << "_sec" << (j+1)
           << " Lsec_o" << (i+1) << " Lsec_o" << (j+1) << " 0.999\n";
    }
}
```

Per-winding sense sources:

```cpp
ss << "Vpri_sense pri_winding_node bridge_out 0\n";
for (size_t i = 0; i < numSec; ++i) {
    ss << "Vsec1_sense_o" << (i+1)
       << " sec_winding_o" << (i+1) << " rect_in_o" << (i+1) << " 0\n";
    ss << "Vout_sense_o" << (i+1)
       << " filter_out_o" << (i+1) << " load_o" << (i+1) << " 0\n";
}
```

Bridge differential probe:

```cpp
ss << "Evab vab 0 mid_A mid_B 1\n";
```

Comment header:

```cpp
ss << "* <Topology> Converter — generated by OpenMagnetics\n";
ss << "* Vin=" << Vin << " Vo=" << Vo << " Fs=" << Fs
   << " D=" << D << " n=" << n << " Lm=" << Lm
   << " Lr=" << Lr << " rectifier=" << rectifierTypeStr << "\n";
ss << "* Generated " << timestamp << "\n";
ss << "*\n";
```

Time step:

```cpp
double maxStep = period / 200.0;  // golden guide §5
```

Saved signals:

```cpp
ss << ".save v(vab) ";
for (size_t i = 0; i < numSec; ++i) {
    ss << "v(vout_o" << (i+1) << ") ";
}
ss << "i(Vpri_sense) i(Vdc) ";
for (size_t i = 0; i < numSec; ++i) {
    ss << "i(Vsec1_sense_o" << (i+1) << ") ";
    ss << "i(Vout_sense_o" << (i+1) << ") ";
}
ss << "\n";
```

### 2.4 simulate_and_extract_operating_points pattern (golden guide §6)

Copy this from `Dab.cpp` lines 1218–1307, substituting node names:

```cpp
std::vector<OperatingPoint> <Topology>::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {

    NgspiceRunner runner;
    if (!runner.is_available()) {
        std::cerr << "[<Topology>] ngspice not available; "
                  << "falling back to analytical model.\n";
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    std::vector<OperatingPoint> results;
    auto inputVoltages = collect_input_voltages(/*...*/);

    for (size_t vIdx = 0; vIdx < inputVoltages.size(); ++vIdx) {
        for (size_t opIdx = 0; opIdx < get_operating_points().size(); ++opIdx) {

            std::string netlist = generate_ngspice_circuit(
                turnsRatios, magnetizingInductance, vIdx, opIdx);

            SimulationConfig config;
            config.netlist = netlist;
            config.simulationTime = numPeriodsToExtract * period;
            config.maxStepTime = period / 200.0;
            config.startTime = numSteadyStatePeriods * period;
            config.keepTempFiles = false;  // NEVER true in production!

            auto simResult = runner.run_simulation(config);

            if (!simResult.success) {
                std::cerr << "[<Topology>] SPICE failed: "
                          << simResult.errorMessage
                          << "; falling back to analytical for this OP.\n";
                auto analyticalOp =
                    process_operating_point_for_input_voltage(
                        inputVoltages[vIdx], opPoints[opIdx],
                        turnsRatios, magnetizingInductance);
                analyticalOp.set_name(
                    inputVoltages[vIdx].label + " [analytical]");
                results.push_back(analyticalOp);
                continue;
            }

            WaveformNameMapping mapping;
            mapping.primaryVoltage = "v(vab)";
            mapping.primaryCurrent = "i(Vpri_sense)";
            for (size_t i = 0; i < numSec; ++i) {
                mapping.secondaryVoltages.push_back(
                    "v(vout_o" + std::to_string(i+1) + ")");
                mapping.secondaryCurrents.push_back(
                    "i(Vsec1_sense_o" + std::to_string(i+1) + ")");
                mapping.outputCurrents.push_back(
                    "i(Vout_sense_o" + std::to_string(i+1) + ")");
            }

            auto op = runner.extract_operating_point(simResult, mapping);
            op.set_name(inputVoltages[vIdx].label);
            results.push_back(op);
        }
    }
    return results;
}
```

### 2.5 PtP NRMSE test pattern (golden guide §8)

Copy from `tests/TestTopologyDab.cpp` (the helper functions
`ptp_interp`, `ptp_nrmse`, `ptp_current`, `ptp_current_time`).

```cpp
TEST_CASE("Test_<Topo>_PtP_AnalyticalVsNgspice") {
    auto j = json::parse(R"({...spec...})");
    OpenMagnetics::<Topology> topo(j);

    auto turnsRatios = std::vector<double>{<computed n>};
    double Lm = <computed Lm>;

    auto opsAnalytical = topo.process_operating_points(turnsRatios, Lm);
    auto opsSpice = topo.simulate_and_extract_operating_points(turnsRatios, Lm);

    REQUIRE(opsAnalytical.size() == opsSpice.size());

    for (size_t i = 0; i < opsAnalytical.size(); ++i) {
        auto& wfA = opsAnalytical[i].get_excitations_per_winding()[0]
                    .get_current().get_waveform();
        auto& wfS = opsSpice[i].get_excitations_per_winding()[0]
                    .get_current().get_waveform();
        double nrmse = ptp_nrmse(wfA, wfS);
        INFO("OP #" << i << " NRMSE = " << nrmse);
        CHECK(nrmse < 0.15);
    }
}
```

### 2.6 Per-OP diagnostic population pattern

At the END of every `process_operating_point_for_input_voltage()`:

```cpp
    // ===== Populate per-OP diagnostics (golden guide §4.6) =====
    lastDutyCycleLoss = dcl;
    lastEffectiveDutyCycle = Deff;
    lastPrimaryPeakCurrent = ipPk;
    lastZvsMarginLagging =
        PwmBridgeSolver::compute_zvs_margin_lagging(Lk, ipPk, mosfetCoss, Vbus);
    lastZvsLoadThreshold = /* compute load above which ZVS is achieved */;
    lastResonantTransitionTime =
        PwmBridgeSolver::compute_resonant_transition_time(Lk, mosfetCoss);
    // Topology-specific diagnostics:
    //   lastMode = ...
    //   lastModulationType = ...

    return op;
```

### 2.7 Extra-component waveforms clearing (golden guide §11)

At the START of every `process_operating_points()`:

```cpp
extraLrVoltageWaveforms.clear();
extraLrCurrentWaveforms.clear();
extraLoVoltageWaveforms.clear();
extraLoCurrentWaveforms.clear();
extraCbVoltageWaveforms.clear();
// ... clear ALL extra*Waveforms vectors
```

---

## 3. Per-Model Plans

The phases below are NOT optional. Implement in order; each phase has explicit
acceptance criteria. After each phase, run:

```bash
cd /home/alfonso/OpenMagnetics/MKF/build
ninja -j2 MKF_tests
./tests/MKF_tests "TestTopology<Topology>"
```

If any test fails, STOP and fix before moving to the next phase.

---

### 3.A CLLC — Tier-3 → Tier-1 Rewrite

**Total effort**: 10–15 days  
**Action**: **Execute the existing `CLLC_REWRITE_PLAN.md`** in the same
directory. Do NOT patch the current implementation; the FHA-based analytical
model is fundamentally insufficient and must be replaced with a Nielsen TDA
solver.

**Prerequisite reading**: `CLLC_REWRITE_PLAN.md` (this directory),
`Llc.cpp` lines 272–490 (Nielsen TDA reference).

**11-phase delivery summary** (see plan for details):
1. Header docstring + class skeleton (1 day)
2. Nielsen TDA solver, 5-state vector [iLr1, iLr2, iLm, vCr1, vCr2] (3 days)
3. Asymmetric tank handling (1 day)
4. Reverse-power mode (1 day)
5. Per-OP diagnostics + extra*Waveforms (1 day)
6. SPICE harness (snubbers, GEAR, K-matrix, fallback) (2 days)
7. Test suite ≥ 20 cases including PtP NRMSE ≤ 0.15 (2 days)
8. Reference-design validation (1 day)
9. Magnetic adviser integration (1 day)
10. Performance + cleanup (1 day)
11. Code review + merge (1 day)

**Acceptance**:
- [ ] All ≥20 tests pass
- [ ] PtP NRMSE ≤ 0.15
- [ ] Three documented reference designs reproduced within ±5%
- [ ] No SPICE failure throws (fallback works)
- [ ] `get_extra_components_inputs()` returns Lr1, Cr1, Lr2, Cr2, Cm

---

### 3.B PhaseShiftedHalfBridge — Critical Netlist Bug

**Total effort**: 5–7 days  
**Severity**: **CRITICAL** — current PtP test passes by coincidence; SPICE
output is meaningless because phase shift is never applied.

**Files**:
- `src/converter_models/PhaseShiftedHalfBridge.h`
- `src/converter_models/PhaseShiftedHalfBridge.cpp`
- `tests/TestTopologyPhaseShiftedHalfBridge.cpp`

**The bug** (verify before fixing):

Open `PhaseShiftedHalfBridge.cpp` at `generate_ngspice_circuit()` (around
lines 580–727). Compare to `PhaseShiftedFullBridge.cpp` line 682 which has:

```cpp
double phaseDelay = Deff * halfPeriod;
```

and lines 722–729 which apply `phaseDelay` to the lagging-leg switches via
PULSE source `td` parameter. **PSHB never computes or applies any phase delay**.

**Decision required FIRST**: Choose architecture (a) or (b):

**(a) Keep single-leg + split-cap topology, fix the netlist** (5 days):
- Update header docstring to clearly say "single-leg HB with split-cap;
  bipolar ±Vin/2 swing via complementary 50% switches".
- Remove the `phase_shift` parameter from the model's processing — this
  topology only has duty-cycle control, not phase-shift control.
- The "effective duty" of the model becomes `D_cmd` directly.
- Update PtP test to assert against the correct (single-leg) waveform.
- Document this is NOT a Pinheiro-Barbi 2-leg PSHB.

**(b) Implement true 2-leg phase-shifted HB** (7 days):
- Add a second leg in the netlist.
- Apply `phaseDelay = Deff * halfPeriod` to lagging-leg PULSE sources.
- Update primary winding voltage to span between mid_A and mid_B (not
  mid_sw and split-cap node).
- Update analytical model to match.

**Recommended**: Choice **(b) via AHB_PLAN.md split** — split this class
into `AsymmetricHalfBridge` (the real Imbertson-Mohan AHB) and
`ThreeLevelPshb` (the real Pinheiro-Barbi 2-leg). The
`ASYMMETRIC_HALF_BRIDGE_PLAN.md` already exists and is 10 days. Then
delete `PhaseShiftedHalfBridge` once both replacements ship.

**Phase-by-phase (assuming choice (a) to ship a quick fix; if choosing (b),
follow AHB_PLAN.md instead)**:

#### Phase 1: Verify the bug (0.5 day)

1. Read `PhaseShiftedHalfBridge.cpp:580-727`. Confirm `phaseDelay` is never
   computed, never applied to PULSE `td` arguments.
2. Run `TestTopologyPhaseShiftedHalfBridge.cpp` PtP test with two operating
   points: one at `phase_shift=30°`, one at `phase_shift=90°`. Verify SPICE
   waveforms are nearly identical (proving phase shift is ignored).
3. Document findings in a comment in the .cpp file.

#### Phase 2: Update header docstring (0.5 day)

In `PhaseShiftedHalfBridge.h`, replace the existing docstring with:

```cpp
/**
 * @brief Half-Bridge Converter (Single-Leg + Split-Capacitor Variant)
 *
 * IMPORTANT — this is NOT the true two-leg phase-shifted half-bridge
 * (Pinheiro-Barbi). This class models a single half-bridge leg with a
 * split-capacitor mid-point, producing a bipolar ±Vin/2 primary swing
 * via complementary 50% switching. Duty-cycle control only — no phase
 * shift.
 *
 * For the real 2-leg Pinheiro-Barbi PSHB, see ThreeLevelPshb (TBD).
 * For the Imbertson-Mohan asymmetric HB with DC-blocking cap, see
 * AsymmetricHalfBridge (TBD; planned in ASYMMETRIC_HALF_BRIDGE_PLAN.md).
 *
 * <... rest of docstring per snippet 2.1 ...>
 */
```

#### Phase 3: Rip out phase-shift handling (1 day)

In `PhaseShiftedHalfBridge.cpp`:
- Find every reference to `phase_shift` or `phaseShift` or `Deff` in the
  analytical model. Replace with `D_cmd` directly.
- The "effective duty" is now just the commanded duty (no DCL since there's
  no phase shift to lose).
- In the netlist, ensure the single leg simply uses `D_cmd` for both
  switches' PULSE pulse-width parameters.
- The `mid_cap` is held at Vin/2 by the split caps — keep this.

In `PhaseShiftedHalfBridge.h`:
- Mark `last_zvs_margin_lagging` accessor as N/A or always 0 (single-leg
  HB has no lagging-leg ZVS concept; both switches are leading-leg).
- Keep `lastDutyCycleLoss`, `lastPrimaryPeakCurrent`, etc. — all still apply.

#### Phase 4: Update tests (1 day)

In `tests/TestTopologyPhaseShiftedHalfBridge.cpp`:
- Remove or rewrite tests that exercise `phase_shift` parameter sweeps —
  they were testing nothing.
- Add a `Test_Pshb_DocumentSingleLeg` test that asserts at compile-time
  (via static_assert / mock) the model is single-leg.
- Re-run the PtP test; NRMSE should now match analytical correctly because
  netlist matches what analytical computes.

#### Phase 5: SPICE polish (1 day)

While in the netlist, add anything missing per snippet 2.3:
- Verify snubbers on both switches.
- Verify `.options METHOD=GEAR TRTOL=7`.
- Verify `.ic v(vout_o*) = ...`.
- Verify R_load guard.

#### Phase 6: Run full test suite (0.5 day)

```bash
ninja -j2 MKF_tests
./tests/MKF_tests "TestTopologyPhaseShiftedHalfBridge"
```

All 16 cases must pass. PtP NRMSE must be < 0.15 (since phase shift no
longer creates a discrepancy).

#### Phase 7: Memory + docs update (0.5 day)

Update memory file `project_converter_model_quality_tiers.md` to reflect
that PSHB is now correctly a single-leg HB approximation. Update
`FUTURE_TOPOLOGIES.md` to note the AHB+3LPSHB split is still pending.

**Acceptance**:
- [ ] Header docstring clearly states "single-leg + split-cap, NOT 2-leg PSHB"
- [ ] No `phase_shift` parameter in user-facing API
- [ ] All 16 tests pass
- [ ] PtP NRMSE < 0.15 (real now, not coincidence)
- [ ] Snubbers, GEAR, IC pre-charge, R_load guard all present

---

### 3.C ActiveClampForward — Implement ZVS

**Total effort**: 8–10 days  
**Severity**: HIGH — the topology's defining feature is missing.

**Files**:
- `src/converter_models/ActiveClampForward.h`
- `src/converter_models/ActiveClampForward.cpp`
- `tests/TestTopologyForward.cpp` (shared with SingleSwitch + TwoSwitch)

**Background**: Active-clamp forward achieves ZVS via a resonant transition
between the leakage inductance (Lk) and the parallel combination of
Coss + Cclamp during the dead-time interval. The clamp capacitor stores
the magnetizing current's reset energy at voltage Vc = D·Vin/(1−D), and
during the dead-time, this energy resonates with Lk to charge the main
switch's Coss to zero before turn-on.

**Reference**:
- Andreycak, B., "Active Clamp and Reset Technique Enhances Forward
  Converter Performance", TI U-153, 1994.
- Erickson-Maksimović 3rd ed., chapter 6.4.
- Sabate, J. et al., "Design Considerations for High-Voltage High-Power
  Full-Bridge Zero-Voltage-Switched PWM Converter", APEC 1990.

#### Phase 1: Header docstring + diagnostic fields (1 day)

In `ActiveClampForward.h`:

```cpp
/**
 * @brief Active-Clamp Forward Converter
 *
 * Single-switch forward variant where the magnetizing reset is handled
 * by an active clamp circuit (auxiliary switch + clamp capacitor Cc).
 * The clamp transitions enable ZVS on the main switch via resonance
 * between Lk and Coss+Cc during dead-time.
 *
 * TOPOLOGY:
 *
 *      Vin ──┬──[Sw_main]──┬──[T1 pri]──┬── GND
 *            │             │            │
 *           Cin            └──[Sw_aux]──┴──[Cc]── (clamp loop)
 *            │
 *           GND
 *
 * Modulation: D = duty cycle of Sw_main; Sw_aux is complementary
 * (turns on during 1−D with dead-times Td between transitions).
 *
 * KEY EQUATIONS:
 *   Vc (clamp cap voltage)  = D · Vin / (1 − D)            [Andreycak U-153]
 *   Vo / Vin                = D · n_sec / n_pri            [Erickson 6.4.3]
 *   Magnetizing current
 *     reset (linear)        = Vc / Lm during (1−D)·Ts
 *   Resonant transition
 *     period                = 2π · √(Lk · (Coss + Cc))
 *     usable interval       = ¼ of that = (π/2)·√(Lk·(Coss+Cc))
 *   ZVS condition           = ½·Lk·Im_pk² ≥ ½·(2·Coss + Cc)·Vin²
 *                          → Im_pk ≥ Vin · √((2·Coss+Cc)/Lk)   [Sabate]
 *
 * (... rest per snippet 2.1 ...)
 */
class ActiveClampForward : ... {
private:
    // Standard simulation tuning
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 3;

    // Computed design values
    double computedSeriesInductance = 0;       // Lk (leakage)
    double computedMagnetizingInductance = 0;  // Lm
    double computedOutputInductance = 0;       // Lo
    double computedClampCapacitance = 0;       // Cc
    double computedClampVoltage = 0;           // Vc = D·Vin/(1−D)
    double computedDeadTime = 200e-9;
    double computedDiodeVoltageDrop = 0.6;

    // MOSFET parameters
    double mosfetCoss = 200e-12;

    // Per-OP diagnostics
    mutable double lastClampVoltage = 0.0;
    mutable double lastPrimaryPeakCurrent = 0.0;
    mutable double lastMagnetizingCurrentPeak = 0.0;
    mutable double lastResonantTransitionTime = 0.0;  // (π/2)·√(Lk·(Coss+Cc))
    mutable double lastZvsMargin = 0.0;               // energy margin (J), positive = ZVS
    mutable double lastZvsLoadThreshold = 0.0;        // I_o above which ZVS achieved
    mutable bool   lastZvsAchieved = false;

    // Extra-component waveforms
    mutable std::vector<Waveform> extraClampCapVoltageWaveforms;
    mutable std::vector<Waveform> extraClampCapCurrentWaveforms;
    mutable std::vector<Waveform> extraLoVoltageWaveforms;
    mutable std::vector<Waveform> extraLoCurrentWaveforms;

public:
    // Standard accessors
    int  get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int v) { numPeriodsToExtract = v; }
    // ... (rest per golden guide §2.7–2.10)

    // Diagnostic accessors
    double get_last_clamp_voltage() const { return lastClampVoltage; }
    double get_last_resonant_transition_time() const { return lastResonantTransitionTime; }
    double get_last_zvs_margin() const { return lastZvsMargin; }
    double get_last_zvs_load_threshold() const { return lastZvsLoadThreshold; }
    bool   get_last_zvs_achieved() const { return lastZvsAchieved; }

    // Static analytical helpers
    static double compute_clamp_voltage(double Vin, double D);
    static double compute_resonant_transition_time(double Lk, double Coss, double Cc);
    static double compute_zvs_margin(double Lk, double Im_pk, double Coss, double Cc, double Vin);
    static double compute_zvs_load_threshold(double Lk, double Coss, double Cc, double Vin);
    // ... (other standard methods)
};
```

#### Phase 2: Static analytical helpers (1 day)

In `ActiveClampForward.cpp`:

```cpp
double ActiveClampForward::compute_clamp_voltage(double Vin, double D) {
    return Vin * D / (1.0 - D);
}

double ActiveClampForward::compute_resonant_transition_time(
        double Lk, double Coss, double Cc) {
    return 0.5 * std::numbers::pi * std::sqrt(Lk * (2.0*Coss + Cc));
    // Quarter-period of the LC resonance during dead-time
}

double ActiveClampForward::compute_zvs_margin(
        double Lk, double Im_pk, double Coss, double Cc, double Vin) {
    double E_inductor = 0.5 * Lk * Im_pk * Im_pk;
    double E_caps = 0.5 * (2.0*Coss + Cc) * Vin * Vin;
    return E_inductor - E_caps;  // J; positive → ZVS feasible
}

double ActiveClampForward::compute_zvs_load_threshold(
        double Lk, double Coss, double Cc, double Vin) {
    // Im_pk_min for ZVS, then convert to load current
    double Im_pk_min = Vin * std::sqrt((2.0*Coss + Cc) / Lk);
    // Convert Im to load current: Im depends on Lm magnetizing alone;
    // the ZVS threshold is on TOTAL primary current (Im + Iload reflected).
    // Caller derives the load threshold from this.
    return Im_pk_min;
}
```

#### Phase 3: Update `process_operating_point_for_input_voltage` (2 days)

After computing duty D, magnetizing current, and primary peak current:

```cpp
double Vc = compute_clamp_voltage(inputVoltage, D);
double Lk = computedSeriesInductance;
double Cc = computedClampCapacitance;

// Resonant transition during dead-time (Lk + Coss + Cc resonance)
double tRes = compute_resonant_transition_time(Lk, mosfetCoss, Cc);

// ZVS analysis
double zvsMargin = compute_zvs_margin(Lk, ipPk, mosfetCoss, Cc, inputVoltage);
double zvsLoadThresh =
    compute_zvs_load_threshold(Lk, mosfetCoss, Cc, inputVoltage);
bool zvsOk = (zvsMargin > 0) && (computedDeadTime >= tRes);

// Populate diagnostics
lastClampVoltage = Vc;
lastResonantTransitionTime = tRes;
lastZvsMargin = zvsMargin;
lastZvsLoadThreshold = zvsLoadThresh;
lastZvsAchieved = zvsOk;
lastPrimaryPeakCurrent = ipPk;
```

Also update the CLAMP CAP waveform generation:
- During Sw_main on (0 to D·Ts): Sw_aux off, Cc holds Vc constant.
- During Sw_main off (D·Ts to Ts): Sw_aux on, Cc connected to primary.
  Magnetizing current flows into Cc, charging it from Vc to slightly higher
  then back. Approximate as constant Vc with small ripple.

```cpp
// Cclamp voltage waveform (approximated as DC + small ripple)
Waveform vcWaveform;
for (size_t k = 0; k < N_samples; ++k) {
    double t = k * dt;
    double v = Vc;  // first-order
    if (t > D * Ts) {
        // Add small ripple from magnetizing current charging Cc
        double dV = (Im_avg / Cc) * (t - D*Ts);
        v += dV;
    }
    vcWaveform.push_back(v);
}
extraClampCapVoltageWaveforms.push_back(vcWaveform);

// Cclamp current waveform
Waveform icWaveform;
// ... similar logic; positive during reset, zero otherwise
extraClampCapCurrentWaveforms.push_back(icWaveform);
```

#### Phase 4: Update SPICE netlist (2 days)

In `generate_ngspice_circuit()`, ADD the auxiliary switch and clamp capacitor:

```cpp
ss << "* Active clamp circuit\n";
ss << "S_aux pri_top clamp_node aux_gate 0 SW1\n";
ss << "Rsnub_aux pri_top snub_aux_n 1k\n";
ss << "Csnub_aux snub_aux_n clamp_node 1n\n";
ss << "Cclamp clamp_node mid_T1 " << computedClampCapacitance << "\n";

// Aux gate signal (complementary to main, with dead-time)
ss << "Vaux_gate aux_gate 0 PULSE(0 5 "
   << (D * period + computedDeadTime) << " 1n 1n "
   << ((1.0 - D) * period - 2.0 * computedDeadTime) << " "
   << period << ")\n";
```

Also add the snubber RC across the main switch (currently missing). Replace
the 1MEG snubber with the proper 1k+1n pair. Ensure GEAR is set:

```cpp
ss << ".options METHOD=GEAR TRTOL=7\n";
```

R_load guard:

```cpp
double rLoad = std::max(Vo / Io, 1e-3);
```

#### Phase 5: Add ZVS boundary test (1 day)

In `tests/TestTopologyForward.cpp`, add:

```cpp
TEST_CASE("Test_ActiveClampForward_ZVS_Boundary") {
    // Sweep load from light to heavy; ZVS margin must cross zero at the
    // predicted threshold (within ±20% tolerance).
    auto j = json::parse(R"({...})");
    OpenMagnetics::ActiveClampForward acf(j);

    std::vector<double> loadCurrents = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0};
    std::vector<double> margins;

    for (double iLoad : loadCurrents) {
        // Run with this load
        // ... (set operating point, call process_operating_points)
        margins.push_back(acf.get_last_zvs_margin());
    }

    // Assert margin transitions from negative to positive somewhere
    // in the sweep
    bool hasNegative = false, hasPositive = false;
    for (double m : margins) {
        if (m < 0) hasNegative = true;
        if (m > 0) hasPositive = true;
    }
    CHECK(hasNegative);
    CHECK(hasPositive);

    // Check the threshold matches the analytical prediction
    double predicted = acf.get_last_zvs_load_threshold();
    INFO("Predicted ZVS load threshold: " << predicted);
    // Find the cross-over current
    double measured = -1.0;
    for (size_t i = 1; i < margins.size(); ++i) {
        if (margins[i-1] < 0 && margins[i] > 0) {
            // Linear interpolation
            double frac = -margins[i-1] / (margins[i] - margins[i-1]);
            measured = loadCurrents[i-1] +
                       frac * (loadCurrents[i] - loadCurrents[i-1]);
            break;
        }
    }
    INFO("Measured ZVS load threshold: " << measured);
    CHECK(std::abs(measured - predicted) / predicted < 0.20);
}
```

#### Phase 6: Add NRMSE PtP test for primary current (1 day)

The existing `Test_ActiveClampForward_PtP_AnalyticalVsNgspice` may need
threshold tightening. Currently NRMSE check is loose. Tighten to ≤ 0.15
once SPICE includes Cclamp and snubbers.

#### Phase 7: Run full test suite + polish (1 day)

- Verify all 3 existing tests + 1 new ZVS test pass.
- Total ≥ 4 tests.
- PtP NRMSE ≤ 0.15.
- Documentation: ensure header docstring covers all new diagnostic
  accessors.

**Acceptance**:
- [ ] `lastZvsMargin` populated on every OP
- [ ] `lastZvsAchieved` correctly identifies ZVS feasibility
- [ ] `compute_resonant_transition_time` returns correct (π/2)·√(Lk·(2Coss+Cc))
- [ ] SPICE netlist contains Sw_aux, Cclamp, snubbers, GEAR, R_load guard
- [ ] PtP NRMSE ≤ 0.15
- [ ] ZVS boundary test passes (margin crosses zero at predicted threshold)

---

### 3.D PowerFactorCorrection — Real Circuit + Tests

**Total effort**: 10–15 days  
**Severity**: HIGH — currently behavioral SPICE, no tests, no diagnostics.

**Files**:
- `src/converter_models/PowerFactorCorrection.h`
- `src/converter_models/PowerFactorCorrection.cpp`
- **CREATE**: `tests/TestTopologyPowerFactorCorrection.cpp`

**Out of scope**: 3-phase totem-pole PFC. Add a TODO note for future work.

#### Phase 1: Header docstring + diagnostic fields (1 day)

```cpp
/**
 * @brief Power Factor Correction (Boost-PFC) Converter
 *
 * Single-phase boost-PFC operating over the AC line cycle. The boost
 * inductor sees a sinusoidal envelope at 50/60 Hz with a switching-
 * frequency triangular ripple superimposed.
 *
 * Variants supported:
 *   - CCM (Continuous Conduction Mode)
 *   - DCM (Discontinuous Conduction Mode)
 *   - CrCM (Critical / Boundary Conduction Mode)
 *
 * NOT YET SUPPORTED:
 *   - Totem-pole PFC (uses GaN switches; no rectifier bridge)
 *   - Bridgeless PFC
 *   - Interleaved PFC (multi-phase)
 *
 * KEY EQUATIONS:
 *   Vin(θ) = Vpk · |sin(θ)|                        [rectified envelope]
 *   D(θ)   = 1 − Vin(θ) / Vbus                     [boost duty]
 *   ΔIL(θ) = Vin(θ) · D(θ) · Tsw / L               [ripple]
 *   IL_avg(θ) = √2 · Iout · |sin(θ)| · Vbus / Vrms [load reflection]
 *
 *   Power factor (PF) = Re{S} / |S| = P / (Vrms · Irms)
 *   Total harmonic distortion (THD):
 *     THD = √(I_2² + I_3² + ...) / I_1   (1st harmonic basis)
 *
 * REFERENCES:
 *   [1] Garcia, O. et al. — "Single-Phase Power Factor Correction:
 *       A Survey", IEEE TPE 18(3), 2003.
 *   [2] Erickson-Maksimović 3rd ed., ch. 18.
 *   [3] TI SLUA309 — "Designing the PFC Power Stage".
 *
 * ACCURACY DISCLAIMER: <per snippet 2.1>
 */
class PowerFactorCorrection : ... {
private:
    int numPeriodsToExtract = 1;        // 1 line cycle is enough
    int numSteadyStatePeriods = 2;      // burn-in 2 cycles
    double computedInductance = 0;
    double computedOutputCapacitance = 0;
    double mosfetCoss = 500e-12;
    double computedDiodeVoltageDrop = 0.7;

    // Per-OP diagnostics
    mutable double lastPowerFactor = 0.0;
    mutable double lastCurrentThd = 0.0;            // fraction (0..1)
    mutable double lastPeakInductorCurrent = 0.0;
    mutable double lastRmsInductorCurrent = 0.0;
    mutable double lastEfficiency = 0.0;
    mutable int    lastMode = 0;  // 0=CCM, 1=DCM, 2=CrCM

    mutable std::vector<Waveform> extraInductorVoltageWaveforms;
    mutable std::vector<Waveform> extraInductorCurrentWaveforms;
    mutable std::vector<Waveform> extraOutputCapVoltageWaveforms;
    mutable std::vector<Waveform> extraOutputCapCurrentWaveforms;

public:
    // ... standard accessors per golden guide ...

    double get_last_power_factor() const { return lastPowerFactor; }
    double get_last_current_thd() const { return lastCurrentThd; }
    double get_last_peak_inductor_current() const { return lastPeakInductorCurrent; }
    double get_last_rms_inductor_current() const { return lastRmsInductorCurrent; }
    int    get_last_mode() const { return lastMode; }

    // Static helpers
    static double compute_power_factor(const Waveform& v, const Waveform& i);
    static double compute_thd(const Waveform& i, double fundamentalFreq, double sampleRate);
    static double compute_inductance_ccm(double Vin_min, double Vbus, double D_max,
                                          double dI_max, double Fsw);
};
```

#### Phase 2: Replace behavioral SPICE with real circuit (4 days)

The current netlist uses `B_vin` voltage source with absolute-value sine
and `B_iload` current source. **Replace this** with a real boost-PFC
circuit:

```cpp
ss << "* Boost PFC — generated by OpenMagnetics\n";
ss << "* Vrms=" << Vrms << " Vbus=" << Vbus << " P=" << Pout
   << " Fsw=" << Fsw << " Fline=" << Fline << "\n";
ss << "*\n";

// AC source
ss << "Vac vac 0 SIN(0 " << (Vrms * std::sqrt(2.0)) << " "
   << Fline << ")\n";

// Bridge rectifier (4 diodes)
ss << "D1 vac vrect_p DIDEAL\n";
ss << "D2 0   vrect_p DIDEAL\n";
ss << "D3 vrect_n vac DIDEAL\n";
ss << "D4 vrect_n 0   DIDEAL\n";
ss << "Vrect_n_to_gnd vrect_n 0 0\n";  // tie negative rail to gnd

// Boost inductor
ss << "L_boost vrect_p sw_node " << computedInductance << " IC=0\n";

// MOSFET (active switch)
ss << "Sboost sw_node 0 gate_pwm 0 SW1\n";
ss << "Rsnub_boost sw_node snub_n 1k\n";
ss << "Csnub_boost snub_n 0 1n\n";

// Boost diode
ss << "D_boost sw_node vbus DIDEAL\n";

// Output cap
ss << "Cbus vbus 0 " << computedOutputCapacitance << " IC=" << Vbus << "\n";

// Load
double rLoad = std::max(Vbus * Vbus / Pout, 1e-3);
ss << "Rload vbus 0 " << rLoad << "\n";

// PWM gate signal — varies with line angle for PFC
// For SPICE, simplest is a constant-frequency sawtooth comparator with
// a sinusoidal reference. Use a B-source for the modulator only:
ss << "Bvref vref 0 V=" << "{1.0 - abs(sin(2*pi*" << Fline << "*time))*"
   << (Vrms * std::sqrt(2.0)) << "/" << Vbus << "}\n";
ss << "Vsaw saw 0 PULSE(0 1 0 " << ((1.0 / Fsw) - 1e-9)
   << " 1n 1n " << (1.0 / Fsw) << ")\n";
ss << "Bgate gate_pwm 0 V={if(v(saw)<v(vref), 5, 0)}\n";

// Switch + diode models
ss << ".model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg\n";
ss << ".model DIDEAL D(IS=1e-12 RS=0.05)\n";

// Solver options
ss << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
ss << ".options METHOD=GEAR TRTOL=7\n";

// Time step — must resolve switching ripple AND line cycle
double maxStep = (1.0 / Fsw) / 50.0;  // 50 points per switch period
ss << ".tran " << maxStep << " " << (numPeriodsToExtract / Fline) << "\n";

// Saved signals
ss << ".save v(vac) v(vbus) v(sw_node) v(vrect_p) "
   << "i(L_boost) i(Vac) i(Cbus)\n";

ss << ".end\n";
```

#### Phase 3: PF + THD calculation (1 day)

```cpp
double PowerFactorCorrection::compute_power_factor(
        const Waveform& v, const Waveform& i) {
    if (v.size() != i.size() || v.size() < 2) return 0.0;
    // Real power = (1/T) ∫ v·i dt
    double pSum = 0, vSqSum = 0, iSqSum = 0;
    for (size_t k = 0; k < v.size(); ++k) {
        pSum += v[k] * i[k];
        vSqSum += v[k] * v[k];
        iSqSum += i[k] * i[k];
    }
    double P = pSum / v.size();
    double Vrms = std::sqrt(vSqSum / v.size());
    double Irms = std::sqrt(iSqSum / i.size());
    double S = Vrms * Irms;
    return (S > 0) ? (P / S) : 0.0;
}

double PowerFactorCorrection::compute_thd(
        const Waveform& i, double fundamentalFreq, double sampleRate) {
    // FFT of the inductor current waveform; ratio of harmonics to fundamental
    // Use a simple Goertzel for each harmonic 2..10
    double f1 = fundamentalFreq;
    double a1 = goertzel_amplitude(i, f1, sampleRate);
    double sumHarmSq = 0;
    for (int k = 2; k <= 10; ++k) {
        double ak = goertzel_amplitude(i, k * f1, sampleRate);
        sumHarmSq += ak * ak;
    }
    return (a1 > 0) ? (std::sqrt(sumHarmSq) / a1) : 0.0;
}
```

#### Phase 4: Populate diagnostics (1 day)

In `process_operating_points()`, after generating the inductor waveform:

```cpp
lastPowerFactor = compute_power_factor(inputVoltageWaveform, inductorCurrentWaveform);
lastCurrentThd = compute_thd(inductorCurrentWaveform, Fline, 1.0/dt);
lastPeakInductorCurrent = *std::max_element(
    inductorCurrentWaveform.begin(), inductorCurrentWaveform.end());
lastRmsInductorCurrent = compute_rms(inductorCurrentWaveform);
// Detect mode by checking if iL ever reaches zero during a switching period
lastMode = detect_pfc_mode(inductorCurrentWaveform, Fsw, Fline);
```

#### Phase 5: simulate_and_extract_operating_points (1 day)

Implement per snippet 2.4 with `WaveformNameMapping` matching the new
netlist node names (`v(vac)`, `v(vbus)`, `i(L_boost)`, etc.).

#### Phase 6: Create test file (3 days)

Create `tests/TestTopologyPowerFactorCorrection.cpp` with at minimum:

```cpp
TEST_CASE("Test_Pfc_Reference_Design") {
    // 230 Vrms in, 400 V bus, 1 kW, target PF > 0.95, THD < 5%
    auto j = json::parse(R"({
        "inputVoltage": {"rms": 230, "frequency": 50},
        "operatingPoints": [{"outputVoltage": 400, "outputCurrent": 2.5,
                              "switchingFrequency": 65000}]
    })");
    OpenMagnetics::PowerFactorCorrection pfc(j);
    auto ops = pfc.process_operating_points({1.0}, 0.0);
    REQUIRE(!ops.empty());
    CHECK(pfc.get_last_power_factor() > 0.95);
    CHECK(pfc.get_last_current_thd() < 0.05);
    CHECK(pfc.get_last_mode() == 0);  // CCM
}

TEST_CASE("Test_Pfc_Mode_Detection") {
    // CCM at high power, DCM at light load
    // Sweep load from 100W to 1kW; verify mode transition
}

TEST_CASE("Test_Pfc_Spice_Netlist") {
    // Verify netlist contains real switch/diode models, GEAR, snubbers
    // Parse netlist string; assert it contains:
    //   ".model SW1"
    //   ".model DIDEAL"
    //   "Sboost"
    //   "Rsnub_boost"
    //   "Csnub_boost"
    //   ".options METHOD=GEAR"
}

TEST_CASE("Test_Pfc_PtP_AnalyticalVsNgspice") {
    // Run both analytical and SPICE; compare iL waveforms
    // NRMSE ≤ 0.20 (PFC harder than PWM bridges due to envelope variation)
}

TEST_CASE("Test_Pfc_Inductor_Sizing") {
    // Sweep L; verify ripple stays under target
}

TEST_CASE("Test_Pfc_Output_Cap_Sizing") {
    // Verify Vbus ripple at 2*Fline is within spec
}
```

**Acceptance**:
- [ ] Header docstring complete with ASCII art, equations, references, totem-pole disclaimer
- [ ] Diagnostic fields `lastPowerFactor`, `lastCurrentThd`, `lastPeakInductorCurrent`, `lastMode` populated
- [ ] SPICE netlist uses real switch + diode models (not B-sources)
- [ ] All ≥ 6 tests pass
- [ ] PF > 0.95 in reference design
- [ ] PtP NRMSE ≤ 0.20

---

### 3.E Boost / Buck — Quality Baseline (parallel work)

**Total effort**: 3–5 days each. Run in parallel by two agents.

**Files (Boost)**:
- `src/converter_models/Boost.h` (lines 1–13 currently no docstring)
- `src/converter_models/Boost.cpp` (lines 35–400)
- `tests/TestTopologyBoost.cpp` (3 tests currently)

**Files (Buck)**:
- `src/converter_models/Buck.h` (lines 1–14 no docstring)
- `src/converter_models/Buck.cpp` (lines 34–391)
- `tests/TestTopologyBuck.cpp` (5 tests currently)

#### Phase 1: Header docstring (0.5 day)

For Boost (`Boost.h`), insert before `class Boost`:

```cpp
/**
 * @brief Boost (Step-Up) DC-DC Converter
 *
 * Non-isolated boost converter with diode rectification (or
 * synchronous rectification via second MOSFET, future enhancement).
 *
 *  TOPOLOGY:
 *
 *      Vin ──[L]──┬──[D]── Vout ── Cout ── GND
 *                 │                  │
 *                [Sw]───────── ─ ─ ─ ┘
 *                 │
 *                GND
 *
 *  KEY EQUATIONS (CCM):
 *    Conversion ratio M = Vo/Vin = 1 / (1 − D)         [Erickson 6.1]
 *    Inductor ripple ΔIL = Vin · D / (L · Fsw)
 *    Critical L for CCM L_crit = D·(1−D)²·Vo / (2·Fsw·Iout)
 *    Output cap ripple ΔVo = Iout·D / (Cout·Fsw)
 *
 *  DCM detection: IL_min < 0  →  switch to DCM mode
 *
 *  REFERENCES:
 *    [1] Erickson-Maksimović 3rd ed., ch. 6.1.
 *    [2] TI SLVA372 — "Boost Converter Design".
 *    [3] Pressman 3rd ed., ch. 1.4.
 *
 *  ACCURACY DISCLAIMER: <per snippet 2.1>
 */
class Boost : public MAS::Boost, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;

    // Computed design values
    double computedInductance = 0;
    double computedDutyCycle = 0;
    double computedDiodeVoltageDrop = 0.6;
    double mosfetCoss = 200e-12;

    // Per-OP diagnostics
    mutable double lastPeakInductorCurrent = 0.0;
    mutable double lastRmsInductorCurrent = 0.0;
    mutable double lastDcmThreshold = 0.0;  // load below which DCM
    mutable bool   lastIsCcm = true;
    mutable double lastEfficiency = 0.0;

    // Extra-component waveforms
    mutable std::vector<Waveform> extraInductorVoltageWaveforms;
    mutable std::vector<Waveform> extraInductorCurrentWaveforms;
    mutable std::vector<Waveform> extraOutputCapCurrentWaveforms;

public:
    // ... add all standard accessors ...

    double get_last_peak_inductor_current() const { return lastPeakInductorCurrent; }
    double get_last_rms_inductor_current() const { return lastRmsInductorCurrent; }
    bool   get_last_is_ccm() const { return lastIsCcm; }
    double get_last_dcm_threshold() const { return lastDcmThreshold; }
    double get_last_efficiency() const { return lastEfficiency; }

    static double compute_inductance_ccm(double Vin, double Vo, double Iout,
                                         double dI_ratio, double Fsw);
    static double compute_dcm_threshold_load(double L, double Vo, double Vin, double Fsw);
    static bool   detect_dcm(const std::vector<double>& iL);
};
```

For Buck, the equivalent (M = D, etc.). See `Buck.h` rewrite below.

#### Phase 2: SPICE netlist fixes (1 day)

In `Boost.cpp` `generate_ngspice_circuit()` (around line 244):

**Add snubber RC**:
```cpp
ss << "Sboost sw_node 0 gate_pwm 0 SW1\n";
ss << "Rsnub_boost sw_node snub_n 1k\n";
ss << "Csnub_boost snub_n 0 1n\n";
```

**Add R_load guard** (replace direct division around line 313):
```cpp
double rLoad = std::max(outputVoltage / outputCurrent, 1e-3);
ss << "Rload " << outputNode << " 0 " << rLoad << "\n";
```

**Add GEAR** (around line 327):
```cpp
ss << ".options METHOD=GEAR TRTOL=7\n";
```

**Fix saved signals**:
```cpp
ss << ".save v(vin) v(vout) v(sw_node) i(L_boost) i(Vin_sense)\n";
```

Same fixes for `Buck.cpp` `generate_ngspice_circuit()` around lines 235–323.

#### Phase 3: Per-OP diagnostics (1 day)

In `Boost.cpp` `process_operating_points_for_input_voltage()`:

```cpp
// After computing iL waveform:
lastPeakInductorCurrent =
    *std::max_element(iL.begin(), iL.end());
lastRmsInductorCurrent = compute_rms(iL);
lastIsCcm = !detect_dcm(iL);
lastDcmThreshold = compute_dcm_threshold_load(L, Vo, Vin, Fsw);

// Efficiency: P_loss / P_out
double pCond = lastRmsInductorCurrent * lastRmsInductorCurrent * Rcond;
double pSwitching = mosfetCoss * Vo * Vo * Fsw;
double pDiode = computedDiodeVoltageDrop * iL_avg;
lastEfficiency = Pout / (Pout + pCond + pSwitching + pDiode);
```

#### Phase 4: Expand tests (2 days)

Boost — add to `TestTopologyBoost.cpp`:

```cpp
TEST_CASE("Test_Boost_Reference_Design") {
    // Reproduce TI SLVA372 example: 12 V → 24 V @ 2 A
    // Verify Vo within ±2%, ripple within ±10%
}

TEST_CASE("Test_Boost_Design_Requirements") {
    // Round-trip: spec → process → verify L > 0, D in (0,1)
}

TEST_CASE("Test_Boost_DCM_Boundary") {
    // Sweep load from 0.1 A to 5 A; verify mode transitions at predicted threshold
}

TEST_CASE("Test_Boost_OperatingPoints_Generation") {
    // Multi-OP: 3 input voltages × 3 loads
}

TEST_CASE("Test_Boost_PtP_AnalyticalVsNgspice") {
    // NRMSE ≤ 0.15 on inductor current waveform
}

TEST_CASE("Test_Boost_Spice_Netlist") {
    // Assert netlist contains snubber, GEAR, R_load guard
    auto netlist = boost.generate_ngspice_circuit(L);
    CHECK(netlist.find(".model SW1") != std::string::npos);
    CHECK(netlist.find("Rsnub_boost") != std::string::npos);
    CHECK(netlist.find(".options METHOD=GEAR") != std::string::npos);
}
```

Buck — analogous tests, including a sync-rect mode test.

**Acceptance per model**:
- [ ] Header docstring with ASCII, equations (M=D for buck, M=1/(1-D) for boost), references
- [ ] Snubber RC across switch in SPICE
- [ ] R_load divide-by-zero guard
- [ ] `.options METHOD=GEAR TRTOL=7`
- [ ] Diagnostic fields `lastPeakInductorCurrent`, `lastIsCcm`, etc.
- [ ] ≥ 6 tests including PtP NRMSE ≤ 0.15 and DCM boundary
- [ ] All tests pass

---

### 3.F Flyback — Polish to Tier-1 (quick win)

**Total effort**: 2 days  
**Files**:
- `src/converter_models/Flyback.h` (no docstring)
- `src/converter_models/Flyback.cpp` (lines 260–470 SPICE)
- `tests/TestTopologyFlyback.cpp` (23 tests)

Flyback already has 23 tests and a correct mode classifier (CCM/DCM/QRM/BMO).
Only polish needed.

#### Phase 1: Header docstring (0.5 day)

```cpp
/**
 * @brief Flyback Converter (multi-output capable)
 *
 * Buck-boost-derived isolated converter with energy stored in a coupled
 * inductor (transformer with airgap). Supports multiple secondaries via
 * power-share allocation.
 *
 * TOPOLOGY:
 *
 *      Vin ──[Sw]──┬─ T1 pri    T1 sec ─[D]─ Vo1 ── Cout1 ── GND_sec
 *                  │  ─────  ▏        │
 *                 GND  Lp Ls │        │
 *                            │        ⊥ (additional sec windings)
 *
 * MODES:
 *   CCM (Continuous Conduction)     — IL_pri never reaches 0
 *   DCM (Discontinuous Conduction)  — IL_pri = 0 for some interval
 *   QRM (Quasi-Resonant)            — turn-on at valley of resonant ringing
 *   BMO (Boundary / Critical Mode)  — turn-on right when IL = 0
 *
 * KEY EQUATIONS:
 *   Vo = Vin · D / (n·(1−D))                     [CCM]    [Erickson 6.3]
 *   Vo = Vin · D · √(Tsw / (2·Lp·Iout/Vin))      [DCM]
 *   ΔIL_pri = Vin · D · Tsw / Lp
 *
 * REFERENCES:
 *   [1] Erickson-Maksimović 3rd ed., ch. 6.3.
 *   [2] Pressman 3rd ed., ch. 4.
 *   [3] TI SLUP254 — "Flyback Converter Design".
 *
 * ACCURACY DISCLAIMER: <per snippet 2.1>
 */
```

#### Phase 2: SPICE netlist fixes (0.5 day)

In `Flyback.cpp` `generate_ngspice_circuit()`:

**Add snubber across primary switch** (currently absent):
```cpp
ss << "Rsnub_main sw_node snub_n 1k\n";
ss << "Csnub_main snub_n 0 1n\n";
```

**Add R_load guard** (line 359 currently unguarded):
```cpp
double rLoad_i = std::max(Vo_i / Io_i, 1e-3);
```

**Add GEAR** (line 378 currently missing):
```cpp
ss << ".options METHOD=GEAR TRTOL=7\n";
```

#### Phase 3: Add diagnostic fields (0.5 day)

```cpp
mutable int    lastMode = 0;  // 0=CCM, 1=DCM, 2=QRM, 3=BMO
mutable double lastPrimaryPeakCurrent = 0.0;
mutable double lastSecondaryPeakCurrent = 0.0;
mutable double lastDcmRatio = 0.0;  // (Toff_demag / Toff_total); 1.0 = full DCM
mutable double lastEfficiency = 0.0;

double get_last_mode() const { return lastMode; }
// ... etc.
```

Populate in `process_operating_point_for_input_voltage()`.

#### Phase 4: Add NRMSE PtP test (0.5 day)

Add to `TestTopologyFlyback.cpp`:

```cpp
TEST_CASE("Test_Flyback_PtP_AnalyticalVsNgspice") {
    // Single-output flyback in CCM
    auto opsA = flyback.process_operating_points(...);
    auto opsS = flyback.simulate_and_extract_operating_points(...);
    for (size_t i = 0; i < opsA.size(); ++i) {
        double nrmse = ptp_nrmse(opsA[i].primary, opsS[i].primary);
        CHECK(nrmse < 0.15);
    }
}
```

**Acceptance**:
- [ ] Header docstring complete
- [ ] Snubber RC across primary switch
- [ ] R_load guard in netlist
- [ ] GEAR + TRTOL in solver options
- [ ] `lastMode`, `lastPrimaryPeakCurrent`, `lastEfficiency` populated
- [ ] PtP NRMSE test passes < 0.15
- [ ] All 24 tests pass

---

### 3.G PushPull — Polish to Tier-1

**Total effort**: 2–3 days  
**Files**:
- `src/converter_models/PushPull.h` (no docstring)
- `src/converter_models/PushPull.cpp` (lines 63–832 analytical, 1212–1470 SPICE)
- `tests/TestTopologyPushPull.cpp` (8 tests)

#### Phase 1: Header docstring (0.5 day)

Standard template with ASCII showing center-tap primary + secondary.
Cite Pressman ch. 5, Erickson 6.2.

#### Phase 2: Add primary-switch snubbers (0.5 day)

Currently `PushPull.cpp` lines 1306–1313 have switches S1, S2 with NO
parallel RC (only diode snubbers exist). Add:

```cpp
ss << "Rsnub_s1 mid_pri_a snub1_n 1k\n";
ss << "Csnub_s1 snub1_n 0 1n\n";
ss << "Rsnub_s2 mid_pri_b snub2_n 1k\n";
ss << "Csnub_s2 snub2_n 0 1n\n";
```

#### Phase 3: Add GEAR solver option (0.1 day)

```cpp
ss << ".options METHOD=GEAR TRTOL=7\n";
```

#### Phase 4: Populate Lo waveforms (1 day)

`extraLoVoltageWaveforms` and `extraLoCurrentWaveforms` fields are declared
in the header but never populated. Compute them in
`process_operating_point_for_input_voltage()`:

```cpp
// During T1 active (top FET on, 0 to D·Tsw/2):
//   V_Lo = Vsec_pk - Vo
//   slope = (Vsec_pk - Vo) / Lo
// During T1 freewheel (D·Tsw/2 to Tsw/2):
//   V_Lo = -Vo
//   slope = -Vo / Lo
// (T2 active and freewheel mirror this)

Waveform vLoWaveform, iLoWaveform;
double iLo = ILo_min;  // initial

for (size_t k = 0; k < N_samples; ++k) {
    double t = k * dt;
    double phase = std::fmod(t, Tsw);
    double vLo, slope;

    if (phase < D * Tsw / 2.0) {
        vLo = Vsec_pk - Vo;
        slope = vLo / Lo;
    } else if (phase < Tsw / 2.0) {
        vLo = -Vo;
        slope = -Vo / Lo;
    } else if (phase < Tsw / 2.0 + D * Tsw / 2.0) {
        vLo = Vsec_pk - Vo;
        slope = vLo / Lo;
    } else {
        vLo = -Vo;
        slope = -Vo / Lo;
    }

    iLo += slope * dt;
    vLoWaveform.push_back(vLo);
    iLoWaveform.push_back(iLo);
}

extraLoVoltageWaveforms.push_back(vLoWaveform);
extraLoCurrentWaveforms.push_back(iLoWaveform);
```

#### Phase 5: Add NRMSE PtP test (0.5 day)

Same as Flyback Phase 4.

#### Phase 6: R_load guard (0.1 day)

Find R_load in netlist; wrap with `std::max(..., 1e-3)`.

**Acceptance**:
- [ ] Header docstring complete
- [ ] Snubbers on primary switches (S1, S2)
- [ ] GEAR + TRTOL in solver
- [ ] R_load guard
- [ ] `extraLo*Waveforms` populated
- [ ] PtP NRMSE < 0.15
- [ ] All ≥ 9 tests pass

---

### 3.H IsolatedBuck / IsolatedBuckBoost — Polish

**Total effort**: 3–4 days each  
**Files**:
- `src/converter_models/IsolatedBuck.h/cpp`
- `src/converter_models/IsolatedBuckBoost.h/cpp`
- `tests/TestTopologyIsolatedBuck.cpp` (5 tests)
- `tests/TestTopologyIsolatedBuckBoost.cpp` (3 tests)

#### Phase 1: Header docstrings (0.5 day each)

For IsolatedBuck — note the ambiguous "Flybuck" name:

```cpp
/**
 * @brief Isolated Buck Converter (also known as "Flybuck")
 *
 * Buck converter with a coupled inductor used as a flyback-style
 * isolation transformer. Distinct from Flyback in that the primary
 * inductor sees a buck-like duty pattern rather than a flyback hard-
 * switched transformer.
 *
 * <... rest per template ...>
 */
```

For IsolatedBuckBoost — note hybrid topology:

```cpp
/**
 * @brief Isolated Buck-Boost Converter
 *
 * Hybrid topology: buck on the primary side, flyback-derived isolation
 * on secondary. Conversion ratio M = D / (1 − D), but with isolation.
 *
 * <... rest per template ...>
 */
```

#### Phase 2: Rename methods to snake_case (0.5 day each)

In both `.h` and `.cpp`, rename:
- `processOperatingPointsForInputVoltage` → `process_operating_point_for_input_voltage`

Update all callers (grep for old name).

#### Phase 3: Add diagnostic fields (0.5 day each)

Standard `last*` fields per snippet 2.2.

#### Phase 4: SPICE netlist fixes (1 day each)

- Add snubber RC across switch
- Add R_load guard
- Add GEAR + TRTOL
- For IsolatedBuckBoost: **REMOVE** `keepTempFiles = true` (line 552);
  this is a debug flag that must not ship.

#### Phase 5: Tests (1 day each)

- IsolatedBuck: add NRMSE PtP test (already has multi-output test).
- IsolatedBuckBoost: add multi-output test + NRMSE PtP test.

**Acceptance per model**:
- [ ] Header docstring complete
- [ ] Snake_case method names
- [ ] No debug flags in production
- [ ] Snubber + GEAR + R_load guard in SPICE
- [ ] Diagnostic fields populated
- [ ] PtP NRMSE < 0.15
- [ ] Multi-output test (both)

---

### 3.I SingleSwitchForward / TwoSwitchForward — Polish

**Total effort**: 4–5 days each  
**Files**:
- `src/converter_models/SingleSwitchForward.h/cpp`
- `src/converter_models/TwoSwitchForward.h/cpp`
- `tests/TestTopologyForward.cpp` (shared, 15 cases total across 3 forward variants)

#### Phase 1: Header docstrings (1 day each)

For SingleSwitchForward — emphasize demag winding necessity:

```cpp
/**
 * @brief Single-Switch Forward Converter
 *
 * Classical forward with a third "reset" winding for transformer
 * demagnetization. Switch must withstand 2·Vin during reset.
 *
 *  TOPOLOGY:
 *
 *      Vin ──[Sw]── T1_pri ── T1_sec ──[D1]── Lo ── Vo ── Cout
 *              │              │             │
 *           T1_dem ──[D_dem]── Vin           [D2]── GND_sec
 *
 *  KEY EQUATIONS:
 *    Vo = Vin · D · n_sec / n_pri              [Erickson 6.4.1]
 *    D_max ≤ n_dem / (n_dem + n_pri)           [reset constraint]
 *    Vds_max = Vin · (1 + n_pri/n_dem)
 *
 * <...>
 */
```

For TwoSwitchForward — note no demag winding needed:

```cpp
/**
 * @brief Two-Switch Forward Converter
 *
 * Forward variant with two primary switches; transformer reset via
 * two anti-parallel diodes that clamp primary to ±Vin. No demag winding.
 *
 *  Vds_max = Vin (much lower than single-switch).
 *  D_max ≤ 0.5 (reset takes (1−D)·Tsw).
 *
 * <...>
 */
```

#### Phase 2: Diagnostic fields + populate (1 day each)

Standard `last*` fields. For ActiveClamp see Phase 3.C.

#### Phase 3: Fix snubbers (1 day each)

Current state: `Rsnub_fwd` and `Rsnub_fw` are 1MEG (line 549–550). Replace
with proper 1k+1n pair across each switch:

```cpp
ss << "Rsnub_main sw_node snub_main_n 1k\n";
ss << "Csnub_main snub_main_n 0 1n\n";
```

For TwoSwitch, two switches → two snubber pairs.

#### Phase 4: GEAR + R_load guard (0.5 day each)

Standard fixes.

#### Phase 5: Add NRMSE PtP + tests (1.5 days each)

Tighten existing PtP threshold to 0.15. Add:
- Reference-design test
- Multi-output test (forwards can have multiple secondaries)
- For TwoSwitch: D=0.5 boundary test

**Acceptance per model**:
- [ ] Header docstring complete with demag-winding ASCII (single) or 2-diode reset note (two)
- [ ] Snubbers fixed (1k+1n, not 1MEG)
- [ ] GEAR + R_load guard
- [ ] Diagnostic fields populated
- [ ] PtP NRMSE < 0.15
- [ ] Multi-output test added
- [ ] All ≥ 5 tests per topology pass

---

### 3.J LLC — SPICE Harness Upgrade

**Total effort**: 3–5 days  
**Files**:
- `src/converter_models/Llc.h` (lines 13–39 docstring, missing ASCII)
- `src/converter_models/Llc.cpp` (lines 1156–1271 SPICE; uses behavioral PULSE)
- `tests/TestTopologyLlc.cpp` (21 tests)

#### Phase 1: Add ASCII art + accuracy disclaimer (0.5 day)

Update `Llc.h` lines 13–39 to include ASCII schematic and accuracy
disclaimer per snippet 2.1.

#### Phase 2: Replace behavioral PULSE with real switch models (2 days)

Currently `Llc.cpp` lines 1156–1271 use `Vbridge` PULSE source for the
primary side. Replace with real half-bridge/full-bridge switches:

```cpp
// Before: Vbridge node1 0 PULSE(...)
// After:
ss << ".model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg\n";

// Half-bridge variant
if (variant == HALF_BRIDGE) {
    ss << "Sup vin mid up_gate 0 SW1\n";
    ss << "Rsnub_up vin snub_up_n 1k\n";
    ss << "Csnub_up snub_up_n mid 1n\n";
    ss << "Sdn mid 0 dn_gate 0 SW1\n";
    ss << "Rsnub_dn mid snub_dn_n 1k\n";
    ss << "Csnub_dn snub_dn_n 0 1n\n";

    // Gate signals (complementary, 50% with dead-time)
    double tdead = computedDeadTime;
    ss << "Vup_gate up_gate 0 PULSE(0 5 0 1n 1n "
       << ((period / 2.0) - tdead) << " " << period << ")\n";
    ss << "Vdn_gate dn_gate 0 PULSE(0 5 "
       << ((period / 2.0) + tdead/2) << " 1n 1n "
       << ((period / 2.0) - tdead) << " " << period << ")\n";

    // Bridge differential probe
    ss << "Evab vab 0 mid 0 1\n";
}
// Full-bridge variant: 4 switches
else if (variant == FULL_BRIDGE) {
    // S1, S2 on leg A; S3, S4 on leg B
    // mid_A = primary_in, mid_B = primary_out
    // Snubber pairs on each
    // Vab = mid_A - mid_B
}
```

#### Phase 3: Add ZVS-related diagnostics (0.5 day)

Add to header:
```cpp
mutable double lastZvsMargin = 0.0;        // primary-side switch ZVS energy margin
mutable double lastPrimaryPeakCurrent = 0.0;
```

Compute in solver after Nielsen TDA converges:
```cpp
lastZvsMargin =
    PwmBridgeSolver::compute_zvs_margin_lagging(
        computedSeriesInductance, ipPk, mosfetCoss, inputVoltage);
lastPrimaryPeakCurrent = ipPk;
```

#### Phase 4: Add ZVS boundary test (1 day)

```cpp
TEST_CASE("Test_Llc_ZVS_Boundaries") {
    // Sweep load; verify ZVS margin in resonant region
}
```

#### Phase 5: Verify NRMSE still passes (0.5 day)

After replacing behavioral with real switches, re-run
`Test_Llc_PtP_AnalyticalVsNgspice`. NRMSE should remain ≤ 0.15.

**Acceptance**:
- [ ] ASCII art in header
- [ ] Accuracy disclaimer
- [ ] Real switch models (SW1) + snubbers in SPICE
- [ ] `lastZvsMargin` populated
- [ ] ZVS boundary test added
- [ ] All 22 tests pass

---

### 3.K DifferentialModeChoke — Add Tests

**Total effort**: 2–3 days  
**Files**:
- `src/converter_models/DifferentialModeChoke.h` (good docstring, no ASCII)
- `src/converter_models/DifferentialModeChoke.cpp`
- **CREATE**: `tests/TestTopologyDifferentialModeChoke.cpp` (currently 0 tests)

#### Phase 1: ASCII art + references (0.5 day)

Update `.h` docstring with ASCII showing single inductor in series with
line, parallel C to ground for LC filter.

#### Phase 2: Create test file (2 days)

```cpp
TEST_CASE("Test_Dmc_Single_Phase_Design") {
    auto j = json::parse(R"({
        "lineFrequency": 50,
        "switchingFrequency": 65000,
        "minimumImpedance": 100,
        "configuration": "SINGLE_PHASE"
    })");
    DifferentialModeChoke dmc(j);
    auto ops = dmc.process_operating_points({1.0}, 0.0);
    REQUIRE(!ops.empty());
    // ...
}

TEST_CASE("Test_Dmc_Three_Phase") { /* ... */ }
TEST_CASE("Test_Dmc_Three_Phase_With_Neutral") { /* ... */ }
TEST_CASE("Test_Dmc_Spice_Netlist") {
    // Verify LC topology in netlist
}
TEST_CASE("Test_Dmc_Attenuation_Verification") {
    // verify_attenuation() returns true for compliant designs
}
```

**Acceptance**:
- [ ] ASCII art in docstring
- [ ] References block added
- [ ] ≥ 5 tests in new test file
- [ ] All tests pass

---

### 3.L CurrentTransformer — Architecture Migration

**Total effort**: 2–3 days  
**Files**:
- `src/converter_models/CurrentTransformer.h` (lines 11–24, doesn't inherit Topology)
- `src/converter_models/CurrentTransformer.cpp`
- `tests/TestTopologyCurrentTransformer.cpp` (1 test)

This is a **breaking change**. Review carefully before merging.

#### Phase 1: Inherit Topology base (0.5 day)

```cpp
// Before:
class CurrentTransformer : public MAS::CurrentTransformer { ... };

// After:
class CurrentTransformer : public MAS::CurrentTransformer, public Topology {
    // Add tuning fields, computed values, accessors
};
```

#### Phase 2: Add header docstring (0.5 day)

```cpp
/**
 * @brief Current Transformer (Measurement)
 *
 * Passive current-sense transformer. Primary carries the load current;
 * secondary terminates in a burden resistor + diode rectifier.
 *
 *  TOPOLOGY:
 *
 *    Iload ──[Pri winding 1 turn]── ⊥
 *                                     T1
 *    [Burden] ── [Sec winding N turns] ── [D] ── output to ADC
 *
 *  KEY EQUATIONS:
 *    Vsec = (Iload / N) · Rburden + Vd
 *    Saturation flux check: Bmax = Iload·Lpri / (N·Ae) < Bsat
 *
 *  REFERENCES:
 *    [1] Williamson, "Current Transformers", IEEE TPAS 1980.
 *
 * ACCURACY DISCLAIMER: Linear core model only; saturation modeled
 * separately by Magnetic class.
 */
```

#### Phase 3: Add saturation diagnostic (0.5 day)

```cpp
private:
    mutable double lastFluxDensity = 0.0;
    mutable bool   lastIsSaturated = false;

public:
    double get_last_flux_density() const { return lastFluxDensity; }
    bool   get_last_is_saturated() const { return lastIsSaturated; }
```

Populate in `process_operating_points()`.

#### Phase 4: Expand tests (1 day)

```cpp
TEST_CASE("Test_CurrentTransformer_Design_Requirements") { /* ... */ }
TEST_CASE("Test_CurrentTransformer_Multiple_OperatingPoints") { /* ... */ }
TEST_CASE("Test_CurrentTransformer_Saturation_Check") { /* ... */ }
TEST_CASE("Test_CurrentTransformer_Burden_Sizing") { /* ... */ }
```

**Acceptance**:
- [ ] Inherits both `MAS::CurrentTransformer` and `Topology`
- [ ] Header docstring with ASCII + equations + references
- [ ] Saturation diagnostic
- [ ] ≥ 5 tests pass
- [ ] Existing single test still passes

---

### 3.M CommonModeChoke — Docstring Polish

**Total effort**: 1 day  
**Files**:
- `src/converter_models/CommonModeChoke.h` (good docstring, missing ASCII + refs)

#### Phase 1: Add ASCII + references (1 day)

Update docstring to include:
- ASCII art showing bifilar winding (CM mode coupling)
- References: CISPR 32, Ott "Electromagnetic Compatibility Engineering"
- Equivalent circuit (CM impedance, leakage)

**Acceptance**:
- [ ] ASCII art added
- [ ] ≥ 3 references cited
- [ ] All 15 existing tests still pass

---

## 4. Recommended Sequencing

The plan is structured so independent fixes can run in parallel. The suggested
order:

### Wave 1 — Critical & ship-blocking (15–20 days serial; parallelize to 8 days)

| Track | Items | Days |
|-------|-------|------|
| A | **CLLC rewrite** (3.A) — execute CLLC_REWRITE_PLAN.md | 10–15 |
| B | **PSHB netlist fix** (3.B) — choose (a) or (b); recommend (b) via AHB_PLAN.md | 5–7 |
| C | **ActiveClampForward ZVS** (3.C) | 8–10 |
| D | **PowerFactorCorrection** (3.D) | 10–15 |

These can run in parallel by 4 different agents (no shared files except
golden guide and Topology base).

### Wave 2 — Quality baseline (16–20 days serial; parallelize to 4–5 days)

| Track | Items | Days |
|-------|-------|------|
| A | **Boost** (3.E) | 3–5 |
| B | **Buck** (3.E) | 3–5 |
| C | **Flyback polish** (3.F) | 2 |
| D | **PushPull polish** (3.G) | 2–3 |
| E | **IsolatedBuck polish** (3.H) | 3–4 |
| F | **IsolatedBuckBoost polish** (3.H) | 3–4 |
| G | **SingleSwitchForward** (3.I) | 4–5 |
| H | **TwoSwitchForward** (3.I) | 4–5 |

Parallelize aggressively — these are mostly independent files.

### Wave 3 — Polish (5–7 days serial)

| Track | Items | Days |
|-------|-------|------|
| A | **LLC SPICE upgrade** (3.J) | 3–5 |
| B | **DMC tests** (3.K) | 2–3 |
| C | **CT migration** (3.L) | 2–3 |
| D | **CMC docstring** (3.M) | 1 |

### Wave 4 — Future topologies (only after Waves 1–3 land)

Already-written plans ready to execute:
- `CUK_PLAN.md` — Cuk converter (10 days)
- `WEINBERG_PLAN.md` — Weinberg converter (12.5 days)
- `FOUR_SWITCH_BUCK_BOOST_PLAN.md` — 4SBB (9.5 days)
- `CLLLC_PLAN.md` — CLLLC (~10 days)
- `VIENNA_PLAN.md` — Vienna 3-phase PFC
- `ASYMMETRIC_HALF_BRIDGE_PLAN.md` — AHB (10 days; satisfies PSHB split)
- `SEPIC_PLAN.md`, `ZETA_PLAN.md` — already written
- `SRC_PLAN.md` — Series Resonant (8–10 days)

### Total estimated effort

- Wave 1: 33–47 days serial / **15 days parallel** with 4 agents
- Wave 2: 24–33 days serial / **8 days parallel** with 8 agents
- Wave 3: 8–12 days serial / **5 days parallel** with 4 agents
- **Total to bring all 17 models to tier-1**: 65–92 working days serial
  or **~28 working days with parallelization**.

---

## 5. Acceptance & Verification

### 5.1 Per-phase verification

After every phase in every track:

```bash
cd /home/alfonso/OpenMagnetics/MKF/build
ninja -j2 MKF_tests   # MUST exit 0
./tests/MKF_tests "TestTopology<Topology>"  # MUST exit 0
```

If either fails, **STOP** and fix root cause. Do NOT skip with `--no-verify`.

### 5.2 Per-model acceptance gate (must pass to mark "tier-1")

For every active converter model:

- [ ] Header docstring contains ASCII art, modulation convention, ≥5 equations
      with citations, ≥3 references, accuracy disclaimer, disambiguation note
- [ ] Class struct has `numPeriodsToExtract`, `numSteadyStatePeriods`,
      `computed*` design values, `mosfetCoss`, `mutable double last*`
      diagnostics, `mutable std::vector<Waveform> extra*Waveforms`
- [ ] Analytical model uses piecewise-linear sub-interval propagator
      (NOT forward Euler)
- [ ] Per-OP diagnostics populated at end of
      `process_operating_point_for_input_voltage`
- [ ] SPICE netlist contains:
  - `.model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg`
  - `.model DIDEAL D(IS=1e-12 RS=0.05)` (if rectifier present)
  - 1 kΩ + 1 nF snubber across every active switch
  - `.options METHOD=GEAR TRTOL=7`
  - `.ic` pre-charge for output caps
  - Multi-secondary K-matrix (every pair) if multi-output
  - Per-winding sense sources
  - `R_load = std::max(Vo/Io, 1e-3)` guard
  - Comment header with topology + key params
  - Time step `period/200`
- [ ] `simulate_and_extract_operating_points` actually runs `NgspiceRunner`,
      with `[analytical]` fallback on failure (no throw)
- [ ] `simulate_and_extract_topology_waveforms` references real netlist
      node names (grep verifies match)
- [ ] Test file has:
  - Reference-design test (within ±5% of published value)
  - Design-requirements round-trip test
  - Multi-OP generation test
  - Operating-modes test (CCM/DCM, SPS/EPS, mode 1..6, etc.)
  - ZVS boundary test (if ZVS applies)
  - Rectifier-type tests (if rectifier_type enum exists)
  - SPICE netlist parse test (asserts presence of snubbers, GEAR, etc.)
  - **PtP NRMSE test with NRMSE ≤ 0.15** (the marker of DAB-quality)
  - Multi-output test (if applicable)
- [ ] No debug flags (`keepTempFiles = true`, etc.) in production code

### 5.3 NRMSE realistic targets (calibrate expectations)

| Topology | Expected NRMSE | Source |
|----------|----------------|--------|
| DAB SPS/EPS/DPS/TPS | 0.04–0.08 | Measured |
| PSFB | 0.10–0.14 | Measured (after Im_peak + ILo fix) |
| PSHB (single-leg) | 0.10–0.14 | After Phase 3.B fix |
| Llc | 0.05–0.10 | Nielsen TDA accuracy |
| Cllc | 0.10–0.15 | Target after rewrite |
| Forward family | 0.10–0.20 | Stiff demag dynamics |
| Boost / Buck | 0.05–0.12 | Simple PWM |
| Flyback CCM | 0.05–0.15 | Per-mode |
| Flyback DCM | 0.10–0.20 | DCM has fast transients |
| ActiveClampForward | 0.10–0.20 | After ZVS implementation |
| PFC | 0.15–0.25 | Line-cycle envelope is hard |

Anything above 0.50 = **broken model** (almost certainly using forward Euler
or wrong Im_peak). Investigate root cause; do NOT loosen the threshold.

### 5.4 Final sign-off

Before merging any model PR:

1. [ ] All sections of §5.2 checked
2. [ ] `ninja -j2 MKF_tests` exits 0
3. [ ] `./tests/MKF_tests` exits 0 (full suite)
4. [ ] Memory file `project_converter_model_quality_tiers.md` updated
5. [ ] `FUTURE_TOPOLOGIES.md` updated if topology was on the planned list
6. [ ] PR description references this plan and the specific phase numbers
      that were executed
7. [ ] No memory of golden-guide violations remains in the code (run
      `grep -r "FIXME\|TODO\|HACK" src/converter_models/<Topology>.*` and
      address every hit)

---

## 6. Implementer Quick-Reference Cheat Sheet

When stuck on any phase, look here first:

| Symptom | Likely cause | Fix location |
|---------|--------------|--------------|
| ngspice fails to converge | Missing snubber, no GEAR, no .ic | Snippet 2.3 |
| iL waveform off by 10–15% | Im_peak using Deff instead of D_cmd | Golden guide §4.3 |
| Vo prediction wrong at high Lk | Missing duty-cycle loss | Golden guide §4.4 |
| ILo wrong slope during commutation | Wrong sign on Vo during freewheel | Golden guide §4.4.5 |
| SPICE returns identical Vo across diff phase shifts | Phase shift not applied to gates | Phase 3.B |
| NRMSE > 0.50 | Forward Euler integration | Replace with sub-interval (snippet 2.3) |
| Multi-output current off by 50× | Missing K-matrix between secondaries | Snippet 2.3 K-matrix |
| `getWaveform` returns null | Node names mismatch netlist | grep netlist for exact name |
| Crash on Io=0 OP | Missing R_load guard | `std::max(Vo/Io, 1e-3)` |
| Output cap charges from 0 V | Missing `.ic` pre-charge | Snippet 2.3 |
| ZVS test never crosses zero | `lastZvsMargin` not populated | Phase 3.C / golden guide §4.6 |
| Tests "pass by coincidence" | Compare two operating points with different params; should differ | Phase 3.B |
| Header docstring unclear which mode | Add modulation convention table | Snippet 2.1 |

---

**END OF PLAN**

Total length: ~50 KB. Implementing agents should be able to execute every
phase using this document plus the existing golden guide and the cited
DAB / Llc / PwmBridgeSolver source files. If anything in this plan
contradicts the golden guide, the golden guide wins.
