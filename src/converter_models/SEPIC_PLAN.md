# SEPIC Converter Plan

Companion to `CONVERTER_MODELS_GOLDEN_GUIDE.md`. Adds a SEPIC
(Single-Ended Primary-Inductor Converter) model to MKF. SEPIC is
currently a "hollow" entry — recognised in the MAS topology enum
(`MAS/schemas/inputs/designRequirements.json` and `MasMigration.cpp`)
but with **no C++ implementation**. This plan covers how to add it
end-to-end, hitting every Golden-Guide checkpoint from day one.

This is a planning document — no code has been changed. Sibling to
`LLC_REWORK_PLAN.md`; that file is **not** touched by this plan.

---

## 0. Status snapshot

| Aspect | Today | After this plan |
|---|---|---|
| `Sepic.cpp` / `Sepic.h` | does not exist | new files, DAB-quality from start |
| Schema | hollow enum entry only | full `sepic.json` with `coupledInductor`, `synchronousRectifier`, `currentRippleRatio` |
| Magnetic-adviser path | no SEPIC support | uncoupled (2 inductors) or coupled (1 transformer-style component) |
| Rectifier type | n/a | `DIODE` (default) + `SYNCHRONOUS` flag |
| Operating modes | n/a | CCM + DCM + BCM (boundary detection from Erickson §5.7) |
| Tests | none | TIDA-00781 reference, ZVS-N/A, ripple, NRMSE, multi-OP, SPICE |

**Cousin topologies** (Cuk and Zeta) are separate fourth-order
converters in the same family. They are **out of scope** for this plan
(see § 8 for future-work notes) — the user-confirmed scope is "SEPIC".

---

## 1. Why add SEPIC

1. SEPIC is the only common non-isolated DC-DC converter MKF doesn't
   model. Buck (`Buck.cpp`), Boost (`Boost.cpp`), and isolated
   buck-boost (`IsolatedBuckBoost.cpp`) are present; Buck-boost
   non-isolated is structurally a SEPIC corner-case.
2. **Coupled-inductor SEPIC is a magnetics-design problem.** TI
   SLYT411, Coilcraft, and Würth ANP135 all describe ripple-steering
   between primary and secondary windings as a function of mutual
   coupling and leakage. MKF's magnetic adviser is exactly the tool to
   pick the coupled core that satisfies the ripple-cancellation
   condition, so the SEPIC class is value-additive even before the
   electrical solver lands.
3. The user's existing converter library targets the same audience
   (power-supply designers picking magnetics) — SEPIC rounds out
   coverage at <250 W where it competes with flyback.

---

## 2. Variant taxonomy (industry survey)

I checked TI SNVA168E (AN-1484), TI SLYT309, TI SLYT411, TI SLUA158,
TI SNVAA43, TI TIDA-00781, TI PMP23320, ADI "SEPIC Equations" article,
ADI LTC1871, ADI AN-1126, ADI MAXREFDES1204, onsemi AND90136,
Coilcraft application notes, Würth ANP135, Erickson-Maksimović §5.7
& §6.4, IEEE 7394472, ScienceDirect S2590123025028786, MDPI Energies
15/7936, Wikipedia SEPIC + Ćuk, plus a comparative SEPIC-vs-Zeta
EE Times article. Verdicts:

### Belong as **flags / enums on the SEPIC class** (NOT new models)

| Variant | Mechanism | Why same class |
|---|---|---|
| **Coupled-inductor SEPIC** | `coupledInductor: bool` (default `false`) | Same equations; magnetics build switches from "two inductors" to "one transformer-style coupled component". Most production SEPICs at 25–250 W use this (TI SLYT411). |
| **Synchronous SEPIC** | `synchronousRectifier: bool` (default `false`) | Diode → low-side MOSFET; topology unchanged, only loss model differs. Common at low Vout (LTC1871). |
| **Interleaved SEPIC** | `phaseCount: int` (default 1) | Same per-phase magnetics; phase count is a load-splitting parameter (TI PMP23320). Defer initial implementation; emit `phaseCount=1` only. |

### Belong as **separate model classes**

| Variant | Why distinct | Reference |
|---|---|---|
| **Cuk** | Inverting output, capacitive-only energy transfer through Cs1+Cs2, current-source output. Different gain analytics, different current waveforms. | Wikipedia Ćuk; Erickson Ch. 6 |
| **Zeta** | Non-inverting like SEPIC but Cs sits between L1 and the output node; needs high-side switch drive. Output ripple lower, input ripple higher. | EE Times "SEPIC vs Zeta"; MDPI Energies 15/7936 |
| **Isolated SEPIC** | L2 replaced with a transformer; gain becomes M = n·D/(1−D); secondary referenced to isolated ground. | Patent WO2014064643A2; SBA Brazil |
| **Bridgeless SEPIC PFC** | AC-input, current-shaping loop, line-cycle envelope sweep. Very different operating-point semantics. | IEEE 7394472 |

### Skip outright

- Quasi-resonant / ZVS SEPIC (research-grade, control-mode change rather than magnetics-design change)
- High-gain / charge-pump / tapped-inductor SEPIC (research-grade, no standardisation)

### Direct user-question equivalents (cf. LLC's "current doubler")

If a future request asks *"is coupled-inductor SEPIC missing?"* — answer
mirrors the LLC current-doubler verdict: it's *not* a new class, it's
a flag on the SEPIC class. Sources phrase it as "SEPIC **with**
coupled inductors", never as a distinct topology family.

---

## 3. Existing patterns we will reuse

### Buck / Boost as the closest precedent

`Buck.cpp` and `Boost.cpp` are the structural baseline. Both:

- Inherit `MAS::<topology> + Topology` (Buck.h:~30, Boost.h:~30)
- Implement `calculate_duty_cycle()` — a single static formula
- Implement piecewise-linear CCM/DCM in
  `process_operating_points_for_input_voltage()` with a quadratic
  formula for DCM duty (Buck.cpp:60–67, Boost.cpp:64–71)
- Wire a real `NgspiceRunner` in `simulate_and_extract_operating_points`
  with GEAR + IC + ITL options
- Have a 3–5 case test suite including an analytical-vs-ngspice NRMSE
  test

But both are **PARTIAL quality** vs the Golden Guide:
- No header docstring with ASCII art / equations / citations
- No primary-switch RC snubber in SPICE
- No `R_load = max(Vo/Io, 1e-3)` divide-by-zero guard
- No accuracy-disclaimer block
- No `Test_*_Waveform_Plotting` CSV dump

**Implication for SEPIC**: do NOT copy Buck/Boost wholesale, because
that propagates their gaps. Use Buck.cpp's piecewise-linear solver
*structure*, but bring it up to DAB-quality from day one.

### Coupled inductors map to transformer with leakage

There is no `CoupledInductor` class. From the codebase scan:

- `src/physical_models/LeakageInductance.h` computes leakage between
  winding pairs via FEM
- `src/physical_models/Inductance.cpp::calculate_leakage_inductance_matrix`
  builds an N×N matrix
- `Flyback.cpp` represents the coupled-inductor pair as a transformer
  with two windings and a designed leakage value

**For coupled-inductor SEPIC**: the magnetic component will be a
2-winding transformer with deliberate (non-parasitic) leakage Lk. The
ripple-cancellation condition (TI SLYT411) becomes a magnetic-adviser
constraint: pick the geometry such that Lk satisfies the cancellation
formula. This is exactly the kind of constraint MKF's adviser already
handles — no new primitive needed.

### Schema convention (camelCase)

All `MAS/schemas/inputs/topologies/*.json` files use camelCase enums
and field names. The new file will be
`MAS/schemas/inputs/topologies/sepic.json`, sibling to `buck.json` /
`boost.json` / `isolatedBuckBoost.json`.

---

## 4. Track A — Basic SEPIC (the only track in this plan)

### A.1 Files to create

```
src/converter_models/Sepic.h
src/converter_models/Sepic.cpp
tests/TestTopologySepic.cpp
MAS/schemas/inputs/topologies/sepic.json
```

Plus add `SEPIC` to the topology enum in
`MAS/schemas/inputs/designRequirements.json` (already partially there
per `MasMigration.cpp` — verify the canonical name is `sepic` not
`SEPIC` to match buck/boost style).

### A.2 Schema (`sepic.json`)

Pattern follows `boost.json` (51 lines) plus SEPIC-specific knobs:

```jsonc
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://psma.com/mas/inputs/topologies/sepic.json",
  "title": "sepic",
  "description": "The description of a SEPIC converter excitation",
  "type": "object",
  "properties": {
    "inputVoltage": { "$ref": "../../utils.json#/$defs/dimensionWithTolerance" },
    "diodeVoltageDrop": { "type": "number", "default": 0.6 },
    "maximumSwitchCurrent": { "type": "number" },
    "currentRippleRatio": { "type": "number", "default": 0.4 },
    "couplingCapacitorRippleRatio": { "type": "number", "default": 0.05 },
    "efficiency": { "type": "number", "default": 0.92 },
    "coupledInductor": {
      "description": "If true, L1 and L2 are wound on the same core (1:1 coupled)",
      "type": "boolean", "default": false
    },
    "synchronousRectifier": {
      "description": "If true, the rectifier diode is a low-side MOSFET",
      "type": "boolean", "default": false
    },
    "operatingPoints": {
      "type": "array",
      "items": { "$ref": "#/$defs/sepicOperatingPoint" },
      "minItems": 1
    }
  },
  "required": ["inputVoltage", "operatingPoints"],
  "$defs": {
    "sepicOperatingPoint": {
      "title": "sepicOperatingPoint",
      "allOf": [{ "$ref": "../../utils.json#/$defs/baseOperatingPoint" }]
    }
  }
}
```

After the edit, run the project's MAS code-gen so `MAS::Sepic` becomes
a generated C++ class. Per `project_mas_hpp_staging.md`, copy the
freshly built `MAS.hpp` into `build/MAS/MAS.hpp` if the build complains
about an undeclared `MAS::Sepic`.

### A.3 Class structure (`Sepic.h`)

Follow the DAB pattern from the Golden Guide § 2 verbatim. The minimum
private state:

```cpp
class Sepic : public MAS::Sepic, public Topology {
private:
    // 2.1 Simulation tuning
    int numPeriodsToExtract    = 5;
    int numSteadyStatePeriods  = 3;

    // 2.2 Computed design values (filled by process_design_requirements)
    double computedInputInductance      = 0;   // L1
    double computedOutputInductance     = 0;   // L2
    double computedCouplingCapacitance  = 0;   // Cs
    double computedOutputCapacitance    = 0;   // Co
    double computedDiodeVoltageDrop     = 0.6;
    double computedDutyCycleNominal     = 0;   // at Vin nominal
    double computedDcmBoundaryK         = 0;   // K_crit per Erickson §5.7

    // 2.3 Schema-less device parameters
    double mosfetCoss        = 200e-12;        // for SPICE only
    double diodeReverseTrr   = 0.0;            // off in v1; documented for future
    double mutualCouplingK   = 0.99;           // only if coupledInductor=true

    // 2.4 Per-OP diagnostics
    mutable double lastDutyCycle              = 0.0;
    mutable double lastConductionMode         = 0.0;  // 0=CCM, 1=BCM, 2=DCM
    mutable double lastInputInductorPeak      = 0.0;
    mutable double lastOutputInductorPeak     = 0.0;
    mutable double lastSwitchPeakVoltage      = 0.0;
    mutable double lastSwitchPeakCurrent      = 0.0;
    mutable double lastDiodePeakVoltage       = 0.0;
    mutable double lastDiodePeakCurrent       = 0.0;
    mutable double lastCouplingCapRmsCurrent  = 0.0;

    // 2.5 Extra-component waveforms (cleared in process_operating_points)
    mutable std::vector<Waveform> extraL1CurrentWaveforms;
    mutable std::vector<Waveform> extraL1VoltageWaveforms;
    mutable std::vector<Waveform> extraL2CurrentWaveforms;
    mutable std::vector<Waveform> extraL2VoltageWaveforms;
    mutable std::vector<Waveform> extraCsCurrentWaveforms;
    mutable std::vector<Waveform> extraCsVoltageWaveforms;
    mutable std::vector<Waveform> extraTimeVectors;

public:
    bool _assertErrors = false;

    Sepic(const json& j);
    Sepic() {}

    // Tuning accessors, computed-value accessors, last_* accessors,
    // mosfet_coss / coupling_k accessors — full set per § 2.7–2.10.

    // Overrides
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const SepicOperatingPoint& op,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // Static analytical helpers (callable from tests without instantiation)
    static double compute_duty_cycle(double Vin, double Vo, double Vd, double eff);
    static double compute_conversion_ratio(double D);
    static double compute_l1_for_ripple(double Vin, double D, double Fs, double dILRatio, double Iin);
    static double compute_l2_for_ripple(double Vin, double D, double Fs, double dILRatio, double Iout);
    static double compute_dcm_critical_k(double D);   // K_crit = (1−D)²
    static double compute_coupling_cap_rms(double Iout, double D);
    static double compute_switch_peak_voltage(double Vin, double Vo);
    static double compute_diode_peak_voltage(double Vin, double Vo);
    static double compute_coupling_cap(double dVcsRatio, double Iout, double D, double Fs, double Vin);

    // Extra components factory
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // SPICE
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);
};

class AdvancedSepic : public Sepic {
    // bypasses process_design_requirements; user provides L1, L2, Cs directly
};
```

### A.4 Header docstring

Mandatory blocks (per Golden Guide § 3) before `class Sepic`:

1. **TOPOLOGY OVERVIEW** with ASCII art:

   ```
   Vin ─── L1 ─┬── Cs ──┬─── L2 ─┬── D ─┬── Vout
               │         │        │      │
               │        SW        │      Co
               │         │        │      │
              GND       GND      GND    GND
   ```

   Plus a coupled-inductor variant where L1 and L2 share a core
   (dot-convention diagram).

2. **MODULATION CONVENTION** — duty cycle D ∈ [0,1] is the switch
   ON-fraction. Schema field `efficiency` enters the duty calculation
   as `D = (Vo + Vd) / (Vin·η + Vo + Vd)`.

3. **KEY EQUATIONS** with citations (use the table in § 6 below
   verbatim with `[TI SNVA168E]` style references inline).

4. **References** — full citations (TI app-note titles, Erickson
   chapter, Coilcraft / Würth notes).

5. **Accuracy disclaimer** mirroring `Dab.cpp:296–320`:
   - Coupling-cap ESR not modelled in CCM analytical solver.
   - Coupled-inductor leakage Lk is treated as a single equivalent
     value; ripple-steering from non-uniform Lk is approximated.
   - Diode reverse recovery and MOSFET Coss(Vds) not in v1.
   - DCM solver uses Erickson §5.7 quadratic; deep-DCM (heavy
     subharmonic) accuracy degrades.

6. **Disambiguation note**:
   - **SEPIC** vs **Cuk** — Cuk is inverting, SEPIC is non-inverting.
   - **SEPIC** vs **Zeta** — Zeta puts Cs between L1 and the output
     node, needs high-side switch.
   - **SEPIC** vs **Buck-Boost (non-isolated)** — buck-boost inverts
     output; SEPIC does not. Buck-boost is 2nd-order; SEPIC is 4th.
   - **SEPIC** vs **Flyback** — flyback uses an isolation transformer
     with energy storage on the primary; SEPIC stores energy on L1
     and transfers via Cs. Often interchangeable at <50 W; SEPIC wins
     when isolation is not required and input ripple matters.

### A.5 Analytical solver (`process_operating_point_for_input_voltage`)

Use the **piecewise-linear sub-interval rule** (Golden Guide § 4.1).
SEPIC has two sub-intervals per period:

- **t ∈ [0, D·T]** — switch ON: VL1 = +Vin, VL2 = +VCs ≈ +Vin,
  iD = 0.
- **t ∈ [D·T, T]** — switch OFF: VL1 = Vin − VCs − Vo ≈ −Vo,
  VL2 = −Vo, iD = iL1 + iL2.

Within each sub-interval iL1(t) and iL2(t) are exactly linear, so
propagate by a single multiply-add per sub-interval (Golden Guide § 4.1
"do NOT use forward Euler").

Steady-state initial conditions (Golden Guide § 4.2): use volt-second
balance on each inductor (= 0 over a period) and amp-second balance on
Cs (= 0 over a period). Solve for I_L1_avg = Iout·D/(1−D) and
I_L2_avg = Iout, then set the starting current to
`I_avg − ΔI/2` for each inductor.

Conduction-mode classification: compute K = 2·L_e·Fs / R_load with
L_e = L1‖L2 = L1·L2/(L1+L2), and K_crit = (1 − D)² (Erickson §5.7).
- K > K_crit → CCM (use the linear solver above)
- K = K_crit → BCM (mark `lastConductionMode = 1`)
- K < K_crit → DCM (apply the Erickson quadratic for D' and rebuild
  the third sub-interval where iL1 + iL2 = 0)

Populate every `last_*` diagnostic at the end (Golden Guide § 4.6):
`lastSwitchPeakVoltage = Vin + Vo`,
`lastSwitchPeakCurrent = I_L1_avg + I_L2_avg + (ΔIL1+ΔIL2)/2`,
`lastDiodePeakVoltage = Vin + Vo`,
`lastDiodePeakCurrent = lastSwitchPeakCurrent`,
`lastCouplingCapRmsCurrent = Iout · √(D/(1−D))`.

Sample grid: `2·N + 1` with N = 256, [0, T] (Golden Guide § 4.8).

### A.6 SPICE netlist (`generate_ngspice_circuit`)

Mandatory checkpoints from Golden Guide § 5:

- `.model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg`
- `.model DSEPIC D(IS=1e-12 RS=0.05)`
- **R_snub 1k + C_snub 1n** across the primary MOSFET — mandatory.
- Per-component sense sources: `Vin_sense`, `VL1_sense`, `VL2_sense`,
  `VCs_sense`, `Vd_sense`, `Vout_sense`, `Vsw_sense`.
- **K-coupling** if `coupledInductor=true`: emit
  `K_L1_L2 L_L1 L_L2 <k>` with `k = mutualCouplingK`. Without it ngspice
  treats L1 and L2 as independent (Golden Guide § 5 "K-matrix" warning
  applies even with two windings).
- `.ic v(vout) = <Vo>`, `.ic v(vcs) = <Vin>` — Cs charges to Vin in
  steady state, so this skips a slow warm-up.
- `.options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500`
  `.options METHOD=GEAR TRTOL=7`
- **`R_load = max(Vo/Io, 1e-3)`** — divide-by-zero guard (Golden
  Guide § 5; Buck/Boost are missing this — SEPIC must include it).
- `stepTime = period / 200`.
- Saved signals: `.save v(sw) v(vcs) v(vout) i(Vin_sense)
  i(VL1_sense) i(VL2_sense) i(VCs_sense) i(Vd_sense)`.
- Comment header with topology, Vin, Vo, Fs, D, L1/L2/Cs/Co,
  rectifier type, `coupledInductor` flag, "generated by OpenMagnetics".

If `synchronousRectifier=true`, replace the diode with a low-side
MOSFET driven anti-phase to the main switch (with `dead_time` from
the `mosfet_coss` calculation) and update the `.save` line to capture
the rectifier MOSFET current.

### A.7 `simulate_and_extract_operating_points`

Follow the DAB pattern (Dab.cpp:1218–1307) verbatim:

1. Construct `NgspiceRunner`.
2. If `runner.is_available()` is false, print warning, fall back to
   analytical via `process_operating_point_for_input_voltage`.
   Do not throw.
3. Loop over input voltages × OPs; build netlist via
   `generate_ngspice_circuit`.
4. On `simResult.success == false` log error, fall back analytical,
   tag result name with `[analytical]`.
5. On success, build a `WaveformNameMapping` matching the netlist
   nodes and call `NgspiceRunner::extract_operating_point`.

### A.8 `simulate_and_extract_topology_waveforms`

Round-trip SPICE waveforms to MAS converter-level waveforms. Common
bug to avoid (Golden Guide § 7): the names in `getWaveform("v(...)")`
must exactly match what `generate_ngspice_circuit` emits in `.save`.
Grep the netlist before referencing.

### A.9 `get_extra_components_inputs`

Returns:

- `Inputs` for **L1** (input inductor): nominal inductance from
  `computedInputInductance`, DC bias = `I_L1_avg`, ripple ΔI_L1.
- `Inputs` for **L2** (output-side inductor): nominal inductance from
  `computedOutputInductance`, DC bias = `Iout`, ripple ΔI_L2.
- `CAS::Inputs` for **Cs** (coupling cap): voltage waveform with
  bias = Vin and AC component, RMS current per § 6 equation.
- `CAS::Inputs` for **Co** (output cap) when `mode = REAL`.

When `coupledInductor = true`, replace the two L1 / L2 `Inputs` with a
**single** 2-winding transformer-style `Inputs` whose primary +
secondary excitations carry the L1 and L2 voltage / current waveforms
respectively, and whose desired leakage is set per the
ripple-cancellation condition (TI SLYT411).

### A.10 Tests (`TestTopologySepic.cpp`)

Mandatory cases per Golden Guide § 8:

1. **`Test_Sepic_TIDA00781_Reference_Design`** — TI 12 V → 12 V @ 1 A
   reference design, reproduce within 5 %.
2. **`Test_Sepic_SLUA158_Reference_Design`** — Linear/TI coupled-
   inductor reference, reproduce L1 = L2 sizing and coupling.
3. **`Test_Sepic_Design_Requirements`** — turns ratio (1:1 for
   coupled), L1, L2, Cs all positive and within sane ranges; round-trip
   a power target through the design solver.
4. **`Test_Sepic_OperatingPoints_Generation`** — multiple Vin
   (typical 9–18 V), antisymmetry holds where applicable, both
   inductor currents piecewise linear, output current matches.
5. **`Test_Sepic_Operating_Modes_CCM_DCM`** — sweep load to push past
   K_crit, verify `lastConductionMode` flips 0 → 1 → 2 across BCM.
6. **`Test_Sepic_CoupledInductor_Flag`** — same OP run twice with
   `coupledInductor=false` and `=true`; ripple amplitude differs
   per TI SLYT411 prediction; primary RMS within 10 % of analytical.
7. **`Test_Sepic_SynchronousRectifier_Flag`** — verify Vo improves
   slightly (no Vd loss); SPICE netlist contains rectifier MOSFET.
8. **`Test_Sepic_SPICE_Netlist`** — netlist parses, contains snubber,
   `.options METHOD=GEAR`, K-coupling when coupled, IC pre-charge on
   Cs and Co.
9. **`Test_Sepic_PtP_AnalyticalVsNgspice`** — primary-current shape
   NRMSE ≤ 0.15, primary RMS within 15 %, Vo within 5 % over a
   period after warm-up.
10. **`Test_Sepic_Waveform_Plotting`** — CSV dump under
    `tests/output/sepic/` for visual inspection. Not asserted.
11. **`Test_AdvancedSepic_Process`** — round-trip
    `AdvancedSepic::process()` returns sane Inputs.

ZVS-boundary tests are N/A — SEPIC is hard-switched. The ZVS
checkpoint from Golden Guide § 14 should be marked "N/A — hard-switched
PWM" in the header docstring rather than implemented.

---

## 5. Class checklist (Golden-Guide § 14 acceptance)

- [ ] Builds clean (`ninja -j2 MKF_tests` returns 0)
- [ ] All `TestTopologySepic` cases pass (≥ 11 cases per § A.10)
- [ ] `Test_Sepic_PtP_AnalyticalVsNgspice` NRMSE ≤ 0.15
- [ ] CCM ↔ DCM transition test (`Test_Sepic_Operating_Modes_CCM_DCM`)
      observes `lastConductionMode` change at predicted K_crit
- [ ] Schema fields documented in
      `MAS/schemas/inputs/topologies/sepic.json`
- [ ] Header docstring includes ASCII art, equations with citations,
      references, accuracy disclaimer, disambiguation note (SEPIC vs
      Cuk vs Zeta vs flyback vs buck-boost)
- [ ] `generate_ngspice_circuit` netlist parses and contains: SW1+
      DSEPIC models, R_snub+C_snub across switch, GEAR, ITL=500/500,
      `.ic` on Cs and Co, K-coupling when `coupledInductor`,
      `R_load = max(Vo/Io, 1e-3)`
- [ ] `simulate_and_extract_operating_points` actually runs SPICE and
      falls back analytically with `[analytical]` tag on failure
- [ ] `simulate_and_extract_topology_waveforms` references real netlist
      nodes (grep-verified)
- [ ] Per-OP diagnostic accessors populated and at least one test
      reads them
- [ ] Multi-output **N/A** — SEPIC is single-output. Document this in
      the disclaimer.
- [ ] Synchronous-rectifier flag tested
- [ ] Coupled-inductor flag tested
- [ ] No divide-by-zero in design-time or SPICE code

ZVS-boundary, rectifier-type CT/CD/FB, and duty-cycle-loss checkboxes
from the Golden Guide are **not applicable** to SEPIC — document with
"N/A — hard-switched single-rectifier PWM" notes in the test file.

---

## 6. Equations block (citation table)

Use this verbatim in the header docstring. CCM, lossless unless noted.

| Quantity | Formula | Source |
|---|---|---|
| Conversion ratio | M = Vo/Vin = D/(1−D) | TI SLYT309; Erickson §6.4 |
| Duty cycle (with η, Vd) | D = (Vo+Vd) / (Vin·η+Vo+Vd) | TI SNVA168E |
| Coupling-cap voltage | V_Cs = Vin (steady state) | Wikipedia "SEPIC"; ADI |
| L1 average current | I_L1 = Iout·D/(1−D) | TI SNVA168E |
| L2 average current | I_L2 = Iout | TI SNVA168E |
| L1 ripple | ΔI_L1 = Vin·D / (L1·Fs) | TI SNVA168E |
| L2 ripple | ΔI_L2 = Vin·D / (L2·Fs) | TI SNVA168E |
| L1 sizing for ratio r | L1 = Vin·D / (r·I_L1·Fs) | TI SNVA168E |
| L2 sizing for ratio r | L2 = Vin·D / (r·I_L2·Fs) | TI SNVA168E |
| Coupling-cap RMS current | I_Cs(rms) = Iout·√(D/(1−D)) | TI SLYT309 |
| Coupling-cap sizing | Cs ≥ Iout·D / (ΔV_Cs·Fs)            (use ΔV_Cs ≤ 5 % Vin) | TI SNVA168E |
| Switch peak voltage | V_SW(pk) = Vin + Vo | TI SNVA168E |
| Switch peak current | I_SW(pk) = I_L1 + I_L2 + (ΔI_L1+ΔI_L2)/2 ≈ Iout/(1−D) + ripple | TI SNVA168E |
| Diode peak voltage | V_D(pk) = Vin + Vo | TI SNVA168E |
| Diode peak current | I_D(pk) = I_SW(pk) | TI SNVA168E |
| DCM K | K = 2·L_e·Fs / R_load with L_e = L1·L2/(L1+L2) | Erickson §5.7 |
| DCM K_crit | K_crit = (1−D)² | Erickson §5.7 |
| Coupled L ripple-cancel | Lk1 such that VL1·Lk1/(L1·M_mut) = VL2 (input-ripple zero) | TI SLYT411 |

---

## 7. Reference designs to include in tests

| ID | Source | Use case | Test name |
|---|---|---|---|
| TIDA-00781 | TI | 12 V → 12 V @ 1 A, 12 W, 2-inductor uncoupled | `Test_Sepic_TIDA00781_Reference_Design` |
| SLUA158 | TI / Linear Tech | classic coupled-inductor teaching example | `Test_Sepic_SLUA158_Reference_Design` |
| MAXREFDES1204 | ADI | single-output SEPIC, MAX16990 (automotive) — alternate, optional | (optional) |
| PMP23320 | TI | 2-phase interleaved 500 W — for future interleaved work | (deferred to § 8) |

---

## 8. Future work (out of scope for this plan)

- **`Cuk` model** — sibling fourth-order converter, inverting output.
  Mirror the SEPIC class structure but with the Cuk-specific gain
  (M = −D/(1−D)) and the dual coupling cap. `MasMigration.cpp`
  already maps `"Cuk Converter" → "cukConverter"`.
- **`Zeta` model** — sibling fourth-order, non-inverting,
  high-side-switch. `MasMigration.cpp` maps
  `"Zeta Converter" → "zetaConverter"`.
- **Isolated SEPIC** — adds a transformer in place of L2; gain
  M = n·D/(1−D). Could be a separate class or a flag once the basic
  SEPIC lands.
- **Bridgeless SEPIC PFC** — AC input, line-cycle envelope sweep;
  fundamentally different operating-point semantics. Defer.
- **Interleaved SEPIC** — `phaseCount: int` flag on the basic class.
  Per-phase magnetics identical; load splits 1/N per phase. TI
  PMP23320 is the reference design for the eventual test.
- **Synchronous-rectifier loss model upgrade** — currently treated
  as ideal anti-phase MOSFET; adding gate-driver dead time and
  body-diode conduction during dead time would tighten Vo prediction
  by ~1 % at low Vout.

After Tracks A lands, update
`project_converter_model_quality_tiers.md` to add SEPIC at DAB-quality
tier.

---

## 9. References

- TI AN-1484 / SNVA168E — *Designing A SEPIC Converter*
- TI SLYT309 — *Designing DC/DC converters based on SEPIC topology*
- TI SLYT411 — *Benefits of a coupled-inductor SEPIC*
- TI SLUA158 — *Versatile Low-Power SEPIC Converter*
- TI SNVAA43 — *Basic Calculation of a Coupled-Inductor SEPIC Power Stage*
- TI TIDA-00781 — 12 W SEPIC reference design
- TI PMP23320 — 500 W 2-phase interleaved SEPIC
- ADI — *SEPIC Equations and Component Ratings* (web article)
- ADI LTC1871 datasheet — Boost / Flyback / SEPIC controller
- ADI AN-1126 — *SEPIC Multiplied Boost Converter*
- ADI MAXREFDES1204 — single-output SEPIC reference design
- onsemi AND90136 — *SEPIC Converter Analysis and Design*
- Coilcraft — *Selecting Coupled Inductors for SEPIC Applications*
- Würth Elektronik ANP135 — *SEPIC with coupled and uncoupled inductors*
- Erickson & Maksimović, *Fundamentals of Power Electronics*, 3rd ed.,
  Ch. 5 §5.7 (DCM, SEPIC), Ch. 6 §6.4 (gain derivation)
- Wikipedia — *Single-ended primary-inductor converter*
- EE Times — *No need to fear: SEPIC outperforms the flyback*
- MDPI Energies 15/7936 — *Comparative Study of Buck-Boost, SEPIC,
  Cuk, Zeta*
- IEEE 7394472 — *Single-switch bridgeless SEPIC PFC*
- ScienceDirect S2590123025028786 — *AC-DC SEPIC: state-of-the-art review*
- Patent WO2014064643A2 — *Galvanically isolated SEPIC converter*
- Springer JPE — *High step-up coupled-inductor SEPIC with input ripple
  cancellation*
