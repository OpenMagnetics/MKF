# 4-Switch Buck-Boost (4SBB) Converter Plan

Companion to `CONVERTER_MODELS_GOLDEN_GUIDE.md`. Adds a 4-switch
buck-boost (also called **non-inverting buck-boost**, **H-bridge
buck-boost**, **cascaded buck-boost**) model to MKF. 4SBB is
ubiquitous in laptops, USB-PD chargers, drones, single-cell-Li
chargers, and automotive 12 V↔48 V mild-hybrid rails — every
application where Vin can be **above OR below** Vout.

Distinctive headline: **the single inductor sees three different
excitation patterns** (buck region, boost region, transition region)
selected by the controller across the Vin sweep. The magnetic adviser
must size for the worst-case region, which is non-trivial — the
adviser already does Vin sweeps, but it has never had to swap
solver branches inside a single excitation. This is the
"three-region inductor problem".

This is a planning document — no code changed. The existing plans
(`LLC_REWORK_PLAN.md`, `CLLC_REWRITE_PLAN.md`, `SEPIC_PLAN.md`,
`ZETA_PLAN.md`, `CUK_PLAN.md`, `WEINBERG_PLAN.md`,
`FUTURE_TOPOLOGIES.md`) are not touched by this plan.

---

## 0. Status snapshot

| Aspect | Today | After this plan |
|---|---|---|
| `FourSwitchBuckBoost.cpp` / `.h` | does not exist | new files, DAB-quality from start |
| Schema | nothing | full `fourSwitchBuckBoost.json` with `controlMode`, `transitionMode`, `bidirectional`, `phaseCount` |
| Magnetic-adviser path | no 4SBB support | sweeps Vin and picks worst-case L from {buck @ Vin_max, boost @ Vin_min, buck-boost @ Vin ≈ Vo} |
| Three-region solver | n/a | `process_operating_point_for_input_voltage` branches on Vo/Vin into BUCK / BOOST / BUCK_BOOST |
| Control modes | n/a | `peakCurrent` (TI LM5176-style), `peakBuckPeakBoost` (ADI LT8390-style), `averageCurrent` (USB-PD), `voltageMode` |
| Operating modes (CCM/DCM) | n/a | CCM default; DCM detection per region |
| Tests | none | TIDA-01411 + LM5176EVM-HP + LT8390 references, three-region transition test, worst-case L test, NRMSE, SPICE |

**Single inductor** is canonical (~99 % of commercial parts).
Coupled-inductor variants exist but are not standard for 4SBB and are
out of scope. Bidirectional 4SBB and multi-phase 4SBB are flags on
the same class.

---

## 1. Why add 4SBB

1. **Most-requested missing topology in the buck-boost family.** Every
   modern laptop / USB-PD adapter / drone / single-cell-Li charger
   uses one. TI sells LM5175 / LM5176 / LM5177 / TPS55288 / BQ25710 /
   BQ25713; ADI sells LT8390 / LT8391 / LT8392; Renesas ships
   ISL81801 / ISL78224 / ISL78226; MPS ships MPQ4214 / MP4247.
2. **The three-region inductor problem is unique among MKF
   topologies.** Buck.cpp and Boost.cpp each see one excitation
   pattern; SEPIC / Zeta / Cuk see one (fourth-order PWM). 4SBB is
   the first model where the SAME magnetic component sees buck-shape,
   boost-shape, and buck-boost-shape excitation across the Vin sweep
   — the magnetic adviser must select cores that survive all three.
3. **Saturation is dominated by boost @ Vin_min** (DC bias =
   `Iout·Vo/Vin_min`); **core loss is dominated by buck-boost mode**
   (highest harmonic content). These are different magnetic-design
   constraints — the model must report both so the adviser doesn't
   pick a core that survives one and fails the other.
4. **Production volume** is in the kW class for automotive 12 V↔48 V
   mild-hybrid (Renesas ISL78224 / ISL78226), 100–500 W for laptop
   USB-PD, 50–500 W for drones. All overlap with the magnetics
   range MKF already targets.

---

## 2. Variant taxonomy (industry survey)

Sources: TI SLVA535B (canonical 4SBB equations), TI LM5175 / LM5176 /
LM5177 datasheets, TI SLVAFJ3 (layout), TI SNVA826 (parallel 4SBB),
TI TIDA-01411 reference, TI PMP20042, TI TPS55288, TI BQ25710 / BQ25713,
ADI LT8390 / LT8391 / LT8392 datasheets, ADI MAXREFDES1255, ADI
"H-Bridge Buck-Boost Layout" Analog Dialogue article, Renesas
ISL81801 + ISL78224 + ISL78226, MPS MPQ4214 / MP4247, Restrepo et al.
IEEE 8986423 (operation modes), Yan et al. IEEE 9219154 (MPC smooth
transition), "Three-Mode Peak Current" Springer chapter, Wiley ADC2.223
(RHPZ-elimination modulation), EDN article, Erickson-Maksimović
3rd ed. Ch. 5 + 6, Wikipedia.

### Belong as **flags / enums on the 4SBB class** (NOT new models)

| Variant | Mechanism | Why same class |
|---|---|---|
| **Bidirectional 4SBB** | `bidirectional: bool` (default `false`) | Same hardware, sign of I_L flips; same L sizing per region. Used in V2G / battery / 12 V↔48 V mild-hybrid. References: ISL78224, ISL81801, LM5177. |
| **Multi-phase / interleaved** | `phaseCount: int` (default 1) | Per-phase magnetics identical; load splits 1/N per phase. Used > 200 W (automotive 12 V↔48 V). References: ISL78224 (4-phase), ISL78226. |
| **Control mode** | `controlMode: enum` (`peakCurrent`, `peakBuckPeakBoost`, `averageCurrent`, `voltageMode`) | Determines duty assignments in the buck-boost region (Mode 1 vs Mode 2). The model needs this to compute waveforms, not the magnetics adviser. |
| **Buck-boost transition mode** | `transitionMode: enum` (`simultaneous`, `splitPwm`) | Mode 1 (simultaneous Q1+Q3 / Q2+Q4) vs Mode 2 (sequential split-PWM, LM5176/LT8390 default). Affects ripple amplitude in the transition region. |
| **Synchronous vs diode** | always `synchronous` | 4SBB is intrinsically synchronous (4 MOSFETs); a diode-on-the-low-side variant is rare and not worth a flag. |

### Belong as **separate model classes**

| Variant | Why distinct | Reference |
|---|---|---|
| **Three-level / flying-capacitor 4SBB** | Six switches + flying cap; different commutation, different magnetics | research-grade; defer |

### Skip outright

- **Coupled-inductor 4SBB** — not a standard variant (SEPIC / Cuk
  territory). Skip.
- **Quasi-resonant ZVS 4SBB** — research-only.
- **GaN / SiC variants** — same topology, only switch model differs.

### Direct user-question equivalents

If a future request asks *"is bidirectional 4SBB missing?"* — answer
mirrors the prior verdicts: it's a flag on the class, not a separate
class. Sources phrase it as "4SBB **with** bidirectional power flow",
never as a distinct topology.

---

## 3. Existing patterns we will reuse

### Buck.cpp + Boost.cpp as the structural baseline

The buck-region and boost-region analytical solvers in 4SBB are
**literally Buck.cpp's and Boost.cpp's solvers**. The 4SBB class
should call them directly (or via static helpers extracted from
Buck.cpp / Boost.cpp). The novel piece is the buck-boost (transition)
region.

But Buck.cpp / Boost.cpp are themselves **partial-quality** vs the
Golden Guide (per the SEPIC plan audit): no header docstring, no
primary-switch RC snubber, no `R_load = max(Vo/Io, 1e-3)` divide-by-
zero guard. If `FourSwitchBuckBoost.cpp` calls the Buck/Boost solvers,
that gap becomes 4SBB's gap. **Recommended approach**: extract the
Buck and Boost piecewise-linear solvers into static helpers in
`PwmBridgeSolver.h` (which already exists per the SEPIC plan), with
DAB-quality polish (R_load guard, Coss-aware switch model). Both
Buck/Boost AND 4SBB call those helpers. This is a slight refactor of
Buck/Boost but keeps the codebase DRY.

### DAB.cpp as the ngspice-wrapper baseline

For the SPICE side, 4SBB has **four MOSFETs + bootstrap drivers**,
which is the most complex SPICE netlist of any non-isolated converter
in MKF. Use `Dab.cpp:1218–1307` as the simulate-and-extract pattern,
with the Dab-style switch model and snubber discipline.

### No coupled-inductor primitive needed

4SBB uses a single inductor — no transformer, no coupled-inductor
mapping. The magnetic component is a straight power inductor (same as
Buck or Boost). No new primitive in MKF's magnetic library.

### Schema convention (camelCase)

The new file is `MAS/schemas/inputs/topologies/fourSwitchBuckBoost.json`,
sibling to `buck.json` / `boost.json` / `sepic.json` / `zeta.json`.

---

## 4. Track A — Basic 4SBB (the only track in this plan)

### A.1 Files to create

```
src/converter_models/FourSwitchBuckBoost.h
src/converter_models/FourSwitchBuckBoost.cpp
tests/TestTopologyFourSwitchBuckBoost.cpp
MAS/schemas/inputs/topologies/fourSwitchBuckBoost.json
```

Plus add `fourSwitchBuckBoost` to the topology enum in
`MAS/schemas/inputs/designRequirements.json` (verify the canonical
camelCase name matches the buck/boost/sepic/zeta style).

### A.2 Schema (`fourSwitchBuckBoost.json`)

```jsonc
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://psma.com/mas/inputs/topologies/fourSwitchBuckBoost.json",
  "title": "fourSwitchBuckBoost",
  "description": "The description of a 4-switch buck-boost converter excitation",
  "type": "object",
  "properties": {
    "inputVoltage": { "$ref": "../../utils.json#/$defs/dimensionWithTolerance" },
    "maximumSwitchCurrent": { "type": "number" },
    "currentRippleRatio": { "type": "number", "default": 0.4 },
    "outputVoltageRippleRatio": { "type": "number", "default": 0.01 },
    "efficiency": { "type": "number", "default": 0.94 },
    "controlMode": {
      "description": "How the controller assigns duty across the three regions",
      "$ref": "#/$defs/fsbbControlMode",
      "default": "peakCurrent"
    },
    "transitionMode": {
      "description": "Buck-boost-region modulation strategy",
      "$ref": "#/$defs/fsbbTransitionMode",
      "default": "splitPwm"
    },
    "bidirectional": {
      "description": "If true, controller allows reverse power flow",
      "type": "boolean", "default": false
    },
    "phaseCount": {
      "description": "Number of interleaved phases (default 1)",
      "type": "integer", "default": 1, "minimum": 1
    },
    "transitionHysteresisRatio": {
      "description": "Width of the buck-boost band as |1 − Vo/Vin| (default 0.15)",
      "type": "number", "default": 0.15
    },
    "operatingPoints": {
      "type": "array",
      "items": { "$ref": "#/$defs/fsbbOperatingPoint" },
      "minItems": 1
    }
  },
  "required": ["inputVoltage", "operatingPoints"],
  "$defs": {
    "fsbbControlMode": {
      "title": "fsbbControlMode",
      "type": "string",
      "enum": ["peakCurrent", "peakBuckPeakBoost", "averageCurrent", "voltageMode"]
    },
    "fsbbTransitionMode": {
      "title": "fsbbTransitionMode",
      "type": "string",
      "enum": ["simultaneous", "splitPwm"]
    },
    "fsbbOperatingPoint": {
      "title": "fsbbOperatingPoint",
      "allOf": [{ "$ref": "../../utils.json#/$defs/baseOperatingPoint" }]
    }
  }
}
```

After the edit, run the project's MAS code-gen so
`MAS::FourSwitchBuckBoost` becomes a generated C++ class. Per
`project_mas_hpp_staging.md`, copy the freshly built `MAS.hpp` into
`build/MAS/MAS.hpp` if the build complains.

### A.3 Class structure (`FourSwitchBuckBoost.h`)

```cpp
class FourSwitchBuckBoost : public MAS::FourSwitchBuckBoost,
                            public Topology {
public:
    enum class Region { BUCK, BOOST, BUCK_BOOST };

private:
    // 2.1 Simulation tuning
    int numPeriodsToExtract    = 5;
    int numSteadyStatePeriods  = 3;

    // 2.2 Computed design values
    double computedInductance               = 0;   // single inductor L
    double computedOutputCapacitance        = 0;   // Co
    double computedInputCapacitance         = 0;   // Cin
    double computedDutyCycleBuckNominal     = 0;   // D_buck @ Vin_nom
    double computedDutyCycleBoostNominal    = 0;   // D_boost @ Vin_nom
    double computedTransitionVinLow         = 0;   // boundary into BUCK_BOOST
    double computedTransitionVinHigh        = 0;   // boundary out of BUCK_BOOST
    double computedWorstCaseInductancePeak  = 0;   // I_L_peak in worst region
    double computedWorstCaseInductanceRms   = 0;   // I_L_rms in worst region

    // 2.3 Schema-less device parameters
    double mosfetCoss          = 200e-12;
    double bootstrapDeadTime   = 100e-9;
    double minOnTime           = 80e-9;            // limits dead-zone width
    double minOffTime          = 80e-9;

    // 2.4 Per-OP diagnostics
    mutable Region lastRegion                = Region::BUCK;
    mutable double lastDutyCycleBuck         = 0.0;
    mutable double lastDutyCycleBoost        = 0.0;
    mutable double lastConductionMode        = 0.0;  // 0=CCM, 1=BCM, 2=DCM
    mutable double lastInductorPeakCurrent   = 0.0;
    mutable double lastInductorAvgCurrent    = 0.0;
    mutable double lastInductorRmsCurrent    = 0.0;
    mutable double lastInductorRippleAmpl    = 0.0;
    mutable double lastInductorBSatHeadroom  = 0.0;  // I_sat_core / I_pk
    mutable double lastInputCurrentRipple    = 0.0;
    mutable double lastOutputVoltageRipple   = 0.0;
    mutable double lastQ1RmsCurrent          = 0.0;
    mutable double lastQ2RmsCurrent          = 0.0;
    mutable double lastQ3RmsCurrent          = 0.0;
    mutable double lastQ4RmsCurrent          = 0.0;

    // 2.5 Extra-component waveforms (cleared in process_operating_points)
    mutable std::vector<Waveform> extraInductorCurrentWaveforms;
    mutable std::vector<Waveform> extraInductorVoltageWaveforms;
    mutable std::vector<Waveform> extraSw1NodeVoltageWaveforms;  // SW1 (Q1/Q2 mid)
    mutable std::vector<Waveform> extraSw2NodeVoltageWaveforms;  // SW2 (Q3/Q4 mid)
    mutable std::vector<Waveform> extraTimeVectors;

public:
    bool _assertErrors = false;

    FourSwitchBuckBoost(const json& j);
    FourSwitchBuckBoost() {}

    // Tuning, computed, last_*, mosfet_coss, etc. accessors per Golden
    // Guide § 2.7–2.10. The lastRegion accessor returns the enum.

    Region get_last_region() const { return lastRegion; }

    // Overrides
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const FourSwitchBuckBoostOperatingPoint& op,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // Static analytical helpers
    static Region classify_region(double Vin, double Vo, double hysteresisRatio,
                                  double minOnTime, double minOffTime, double Fs);
    static double compute_buck_duty(double Vin, double Vo, double Vd, double eff);
    static double compute_boost_duty(double Vin, double Vo, double Vd, double eff);
    static double compute_buck_boost_duties(double Vin, double Vo, double eff,
                                            double& D_buck_out, double& D_boost_out,
                                            FsbbTransitionMode mode);
    static double compute_inductor_buck_region(double Vin, double Vo, double Iout,
                                               double Fs, double rippleRatio);
    static double compute_inductor_boost_region(double Vin, double Vo, double Iout,
                                                double Fs, double rippleRatio);
    static double compute_worst_case_inductor(double VinMin, double VinMax, double Vo,
                                              double Iout, double Fs, double rippleRatio);
    static double compute_inductor_dc_bias_boost(double Vin, double Vo, double Iout);  // = Iout·Vo/Vin
    static double compute_dcm_critical_l_buck(double Vo, double R_load, double Fs);
    static double compute_dcm_critical_l_boost(double D_boost, double R_load, double Fs);

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

class AdvancedFourSwitchBuckBoost : public FourSwitchBuckBoost {
    // bypasses process_design_requirements; user provides L directly
};
```

### A.4 Header docstring

Mandatory blocks (per Golden Guide § 3) before the class:

1. **TOPOLOGY OVERVIEW** with ASCII art:

   ```
   Vin ──┬── Q1 ──┬── L ──┬── Q3 ──┬── Vout
         │   (HS) │       │   (HS) │
         │        │       │        │
         Cin     Q2      [SW2]    Q4         Co
         │   (LS) │       │   (LS) │         │
         │        │       │        │         │
        GND ─── GND ───── GND ─── GND ───── GND
                          [SW1]
   ```

   With the three operating-region diagrams stacked as separate ASCII
   panels:

   ```
   BUCK region (Vin >> Vo):  Q3 always ON, Q4 always OFF, Q1/Q2 PWM
   BOOST region (Vin << Vo): Q1 always ON, Q2 always OFF, Q3/Q4 PWM
   BUCK_BOOST region (Vin ≈ Vo): all four switches modulate
   ```

   For the buck-boost region, document Mode 1 (simultaneous,
   `transitionMode = simultaneous`) vs Mode 2 (split-PWM,
   `transitionMode = splitPwm` — LM5176/LT8390 default).

2. **MODULATION CONVENTION** —

   - In BUCK: `D_buck = Vo / (Vin·η)`, applied to Q1.
   - In BOOST: `D_boost = 1 − Vin·η/Vo`, applied to Q4.
   - In BUCK_BOOST (split-PWM): both D_buck and D_boost applied
     simultaneously per Tsw, with `M = D_buck / (1 − D_boost)`.

3. **KEY EQUATIONS** — use the table in § 6 below verbatim with
   `[TI SLVA535B]` and `[TI LM5176]` style references inline.

4. **References** — TI SLVA535B / LM5176 / LM5175, ADI LT8390,
   Restrepo IEEE 8986423, Erickson-Maksimović Ch. 5 + 6.

5. **Accuracy disclaimer** mirroring `Dab.cpp:296–320`:
   - Bootstrap-cap dynamics on the high-side gates not modelled.
   - Coss(Vds) variation not represented; SPICE uses linear C.
   - Body-diode reverse recovery omitted.
   - Mode-transition glitches at the BUCK ↔ BUCK_BOOST and
     BUCK_BOOST ↔ BOOST boundaries are NOT simulated; the analytical
     solver picks one region per OP based on the controller
     hysteresis. Real controllers (LM5176, LT8390) blend smoothly,
     introducing additional ripple in a ~5 % Vin band that is
     under-predicted by the model.
   - Inductor saturation: peak current at the BUCK ↔ BUCK_BOOST mode
     boundary during a load step can exceed steady-state predictions
     by ~20 % (per ADI app note); add this margin in the
     `computedWorstCaseInductancePeak` if `bidirectional == true` or
     phaseCount > 1.

6. **Disambiguation note**:
   - **4SBB** vs **2-switch buck-boost (inverting)** — the 2-switch
     version inverts the output and is 2nd-order; 4SBB is non-
     inverting and uses two half-bridges with one inductor.
   - **4SBB** vs **SEPIC / Cuk / Zeta** — all four are
     non-inverting buck-boost. SEPIC / Cuk / Zeta are 4th-order
     (two L + two C); 4SBB is 2nd-order (one L) but uses four
     active switches. SEPIC/Cuk/Zeta have RHPZ in CCM; 4SBB has
     RHPZ only in the boost and buck-boost regions.
   - **4SBB** vs **DAB** — DAB is isolated (transformer); 4SBB is
     non-isolated. Both have 4 switches, but DAB switches modulate
     phase shift, 4SBB modulates duty.

### A.5 Analytical solver (`process_operating_point_for_input_voltage`)

This is the **most complex part** and the section that justifies a
separate model. The algorithm:

1. **Classify the region** via `classify_region(Vin, Vo, hysteresis,
   minOnTime, minOffTime, Fs)`:
   - If `Vo / Vin < 1 − hysteresis` AND duty-cycle constraint
     `D_buck > Fs · minOnTime` is satisfied → `BUCK`.
   - If `Vo / Vin > 1 + hysteresis` AND `D_boost > Fs · minOnTime` →
     `BOOST`.
   - Else → `BUCK_BOOST` (the "dead zone" where neither pure mode
     respects min on-time).
2. **Dispatch**:

   ```cpp
   switch (region) {
   case Region::BUCK:        return solve_buck(...);
   case Region::BOOST:       return solve_boost(...);
   case Region::BUCK_BOOST:  return solve_buck_boost(...);
   }
   ```

3. **`solve_buck`** — call into the static helper extracted from
   Buck.cpp's piecewise-linear CCM/DCM solver. Switch waveforms:
   Q1 PWM at D_buck, Q2 anti-phase, Q3 ON, Q4 OFF.
   Inductor sees +Vin−Vo during ON, −Vo during OFF.

4. **`solve_boost`** — call into the static helper extracted from
   Boost.cpp's solver. Switch waveforms: Q1 ON, Q2 OFF, Q4 PWM at
   D_boost, Q3 anti-phase. Inductor sees +Vin during ON,
   −(Vo−Vin) during OFF.

5. **`solve_buck_boost`** — Zeta/Sepic-style two-or-four sub-interval
   solver, branching on `transitionMode`:
   - **`simultaneous`** (Mode 1): Q1+Q3 ON during D, Q2+Q4 ON during
     1−D. Inductor sees +Vin during ON, −Vo during OFF (flyback-
     like). M = D / (1−D), but limited to Vo ≈ Vin operation.
     `D = Vo / (Vin + Vo)`. High ripple — flag in disclaimer.
   - **`splitPwm`** (Mode 2, LM5176/LT8390 default): four
     sub-intervals per Tsw. Q1 turns OFF before Q4 turns ON, etc.
     Effective duty `M = D_buck / (1 − D_boost)`. Solved with two
     duties simultaneously by enforcing `Vo = Vin · D_buck / (1 −
     D_boost)` and `D_buck + D_boost = 1 + transitionDelta`
     (controller-specific; default 0.05 per LM5176 datasheet § 8.3).

   Sample grid: `2·N + 1` with N = 256, [0, T] (Golden Guide § 4.8).

6. **Conduction-mode classification per region** — compute K = 2·L·Fs /
   R_load, K_crit per region:
   - Buck: `K_crit_buck = 1 − D_buck`.
   - Boost: `K_crit_boost = D_boost · (1 − D_boost)²`.
   - Buck-boost (Mode 1): `K_crit_bb = (1 − D)²`.
   Erickson § 5.

7. **Populate diagnostics** at the end:

   ```
   lastRegion                 = region
   lastDutyCycleBuck          = D_buck (or 1.0 in BOOST, 0 if undefined)
   lastDutyCycleBoost         = D_boost (or 0 in BUCK)
   lastInductorPeakCurrent    = max(|i_L|)
   lastInductorAvgCurrent     = mean(|i_L|)
   lastInductorRmsCurrent     = sqrt(mean(i_L²))
   lastInductorRippleAmpl     = ΔI_L per region equation
   lastQ1..Q4RmsCurrent       = computed per-switch RMS
   lastInductorBSatHeadroom   = (provided I_sat) / lastInductorPeakCurrent
                                — only populated when magnetic is supplied
   lastConductionMode         = 0/1/2 per K vs K_crit
   ```

8. **`process_design_requirements`** computes the worst-case L:

   ```cpp
   double L_buck  = compute_inductor_buck_region(VinMax, Vo, Iout, Fs, rippleRatio);
   double L_boost = compute_inductor_boost_region(VinMin, Vo, Iout, Fs, rippleRatio);
   double L = std::max(L_buck, L_boost);
   computedInductance = L;
   computedWorstCaseInductancePeak = compute_inductor_dc_bias_boost(VinMin, Vo, Iout)
                                     + compute_buck_boost_ripple(VinMin, Vo, L, Fs) / 2;
   ```

   This is the heart of the magnetic-design contribution: the adviser
   gets a worst-case sizing that respects all three regions across
   the full Vin sweep.

### A.6 SPICE netlist (`generate_ngspice_circuit`)

Mandatory checkpoints from Golden Guide § 5:

- **Four MOSFETs**: Q1 + Q2 (buck leg), Q3 + Q4 (boost leg). Use
  `.model M_NFET NMOS VTO=2.5 KP=20 LAMBDA=0.01` for low-side; for
  high-side Q1 and Q3, use a level-shift / bootstrap subcircuit:

  ```
  V_drive_Q1 vg1_drv 0 PULSE(...)
  E_bootstrap vg1 sw1 vg1_drv 0 1
  Cboot1 vboot1 sw1 100n IC=12
  Dboot1 vcc_aux vboot1 DZBOOT
  ```

  Same pattern for Q3.
- **`.model SW_RECT SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg`** — only if
  any low-side body-diode is preferred in DCM.
- **R_snub 1k + C_snub 1n across each switching node** (SW1 and SW2)
  — mandatory.
- Per-component sense sources: `Vin_sense`, `VL_sense`,
  `Vsw1_sense` (between Q2 drain and L), `Vsw2_sense` (between Q4
  drain and L), `Vout_sense`, plus `Vq1_sense`..`Vq4_sense` for
  per-switch RMS extraction.
- `.ic v(vout) = <Vo>`, `.ic v(sw1) = 0`, `.ic v(sw2) = 0`.
- `.options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500`
- `.options METHOD=GEAR TRTOL=7`
- **`R_load = max(Vo/Io, 1e-3)`** — divide-by-zero guard.
- `stepTime = period / 200`.
- Saved signals: `.save v(sw1) v(sw2) v(vout) i(Vin_sense)
  i(VL_sense) i(Vq1_sense) i(Vq2_sense) i(Vq3_sense) i(Vq4_sense)`.
- Comment header with topology, Vin, Vo, Fs, region, D_buck, D_boost,
  L, Cin, Co, transitionMode, controlMode, "generated by
  OpenMagnetics".

**Region-specific gate-drive emission**:

- In BUCK: emit Q1 PWM at D_buck, Q2 PWM at (1−D_buck) with dead
  time, Q3 = +12 V (always-on bootstrap simulated by a slow PULSE),
  Q4 = 0 V (always off).
- In BOOST: Q1 = +12 V (always-on bootstrap), Q2 = 0 V (always off),
  Q3 PWM at (1−D_boost), Q4 PWM at D_boost.
- In BUCK_BOOST (split-PWM): all four PWM with appropriate phase
  alignment per LM5176 datasheet § 8.3.4.

This region-aware netlist emission is the SPICE-side analogue of the
analytical solver branching, and the most-likely place for bugs. Add
a `Test_FSBB_SPICE_Netlist_Per_Region` test that asserts the gate-
drive PULSE statements for each region.

### A.7 `simulate_and_extract_operating_points`

Follow the DAB pattern (`Dab.cpp:1218–1307`) verbatim, with the
analytical fallback tagged `[analytical]` on failure. Per OP, the
netlist must reflect the analytical region selection — DO NOT emit a
buck-region netlist for an OP that the analytical solver classified
as BUCK_BOOST, or the SPICE-vs-analytical NRMSE will diverge.

### A.8 `simulate_and_extract_topology_waveforms`

Same hygiene rules as SEPIC / Zeta plans: grep the netlist for exact
node names before referencing them.

### A.9 `get_extra_components_inputs`

Returns:

- `Inputs` for **L** (single inductor): nominal inductance from
  `computedInductance`, **DC bias = max DC bias across regions**
  (= `Iout · Vo / Vin_min` from boost region), ripple computed from
  the worst-case region. The `Inputs.processedData` should
  carry per-region waveforms tagged via the
  `extraInductorVoltageWaveforms` vector — one entry per OP, indexed
  by the OP's classified region.

- `CAS::Inputs` for **Cin** (input cap): pulsed-current ripple in
  BUCK region, smaller in BOOST.
- `CAS::Inputs` for **Co** (output cap): pulsed-current ripple in
  BOOST region, smaller in BUCK.

When `phaseCount > 1`, divide the per-phase L value and current
proportionally; otherwise unchanged.

### A.10 Tests (`TestTopologyFourSwitchBuckBoost.cpp`)

Mandatory cases per Golden Guide § 8 + 4SBB-specific tests:

1. **`Test_FSBB_TIDA01411_Reference_Design`** — TI 9–36 V → 12 V @ 8 A
   reference design (LM5175 controller, 350 kHz). Reproduce within
   5 %. Sweep Vin from 9 to 36 V to cover all three regions.
2. **`Test_FSBB_LM5176EVM_Reference_Design`** — TI 9–55 V → 24 V @ 5 A
   reference (LM5176, 200 kHz).
3. **`Test_FSBB_LT8390_Demo_Reference_Design`** — ADI 4–60 V → 12 V @
   7 A demo board DC2825A (LT8390, 300 kHz, peakBuckPeakBoost mode).
4. **`Test_FSBB_Design_Requirements`** — turns ratio (1.0 — single
   inductor, no transformer), L positive, computed L equals
   max(L_buck, L_boost). Round-trip a power target through the
   design solver.
5. **`Test_FSBB_Region_Classification`** — sweep Vo/Vin from 0.3 to
   3.0 and verify `lastRegion` is BUCK below `1 − hysteresis`, BOOST
   above `1 + hysteresis`, and BUCK_BOOST in between. Verify the
   dead-zone width grows with Fs (per minOnTime constraint).
6. **`Test_FSBB_OperatingPoints_Generation_Buck_Region`** — multiple
   Vin in the buck region; inductor current is triangular at D_buck;
   peak = Iout + ΔIL/2.
7. **`Test_FSBB_OperatingPoints_Generation_Boost_Region`** — multiple
   Vin in the boost region; inductor current is triangular at D_boost;
   **peak = Iout/(1−D_boost) + ΔIL/2 — verify it's the worst case**.
8. **`Test_FSBB_OperatingPoints_Generation_BuckBoost_Region`** — Vin
   ≈ Vo; inductor current shows the four-sub-interval split-PWM
   pattern (Mode 2). Compare against Mode 1 (simultaneous) by toggling
   `transitionMode` and verify the ripple amplitude differs by
   the predicted factor.
9. **`Test_FSBB_WorstCase_Inductor_Sizing`** — verify that
   `process_design_requirements()` picks
   `L = max(L_buck@VinMax, L_boost@VinMin)` and that the resulting
   inductor survives the **boost @ VinMin** peak current and the
   **buck-boost mode** core-loss check. **This is the headline
   magnetic-design test.**
10. **`Test_FSBB_Operating_Modes_CCM_DCM_Per_Region`** — sweep load
    in each region; verify CCM → DCM transition at predicted K_crit
    per region.
11. **`Test_FSBB_ControlMode_Variants`** — same OP run with all four
    `controlMode` enum values; verify only the buck-boost-region
    duties differ; buck-only and boost-only OPs are control-mode-
    independent.
12. **`Test_FSBB_TransitionMode_Simultaneous_vs_SplitPwm`** — Mode 1
    (`simultaneous`) shows higher ripple and lower efficiency than
    Mode 2 (`splitPwm`) at the same OP per Restrepo IEEE 8986423.
13. **`Test_FSBB_Bidirectional_Symmetry`** — same OP with
    `bidirectional=true` and reversed sign of Iout produces inductor
    current of equal magnitude with opposite sign; magnetics sizing
    unchanged.
14. **`Test_FSBB_PhaseCount_Interleaved`** — `phaseCount=2`: per-phase
    inductor current is half, ripple cancels at 2·Fs at the input/
    output. Verify against TI SNVA826.
15. **`Test_FSBB_SPICE_Netlist_Per_Region`** — emit the netlist for a
    BUCK OP, a BOOST OP, and a BUCK_BOOST OP; assert each contains
    the expected gate-drive PULSE statements (e.g., for a BUCK OP,
    Q3's gate is held high and Q4's gate is held low; for BUCK_BOOST,
    all four switches have non-trivial PWM).
16. **`Test_FSBB_SPICE_Netlist`** — generic netlist hygiene: contains
    snubbers on both SW nodes, `.options METHOD=GEAR`, ITL=500/500,
    `R_load = max(Vo/Io, 1e-3)`, IC pre-charge on Co.
17. **`Test_FSBB_PtP_AnalyticalVsNgspice_Buck_Region`** —
    inductor-current shape NRMSE ≤ 0.15, RMS within 15 %, Vo within
    5 % over a period after warm-up.
18. **`Test_FSBB_PtP_AnalyticalVsNgspice_Boost_Region`** — same
    targets in boost region.
19. **`Test_FSBB_PtP_AnalyticalVsNgspice_BuckBoost_Region`** — relax
    NRMSE to 0.20 in this region (the analytical model is less
    accurate per the disclaimer); still assert period-averaged Vo
    within 5 %.
20. **`Test_FSBB_DeadZone_Width_vs_Frequency`** — sweep Fs and
    confirm the buck-boost band widens at higher Fs (because
    `minOnTime / Tsw` grows). Per LM5176 E2E forum.
21. **`Test_FSBB_Waveform_Plotting`** — CSV dump under
    `tests/output/four_switch_buck_boost/<region>/` per region for
    visual inspection. Not asserted.
22. **`Test_AdvancedFourSwitchBuckBoost_Process`** — round-trip
    `AdvancedFourSwitchBuckBoost::process()` returns sane Inputs.

ZVS-boundary tests are N/A — 4SBB is hard-switched (though some
controllers like LT8390 add quasi-resonant ZVS on the boost leg; not
modelled in v1).

Tests (5), (8), (9), (12), and (19) are the headline 4SBB-specific
tests that justify a separate model from Buck and Boost. They prove
the three-region solver works and that the magnetic-adviser pass
picks the worst-case core across the full Vin sweep.

---

## 5. Class checklist (Golden-Guide § 14 acceptance)

- [ ] Builds clean (`ninja -j2 MKF_tests` returns 0)
- [ ] All `TestTopologyFourSwitchBuckBoost` cases pass (≥ 22 cases per § A.10)
- [ ] `Test_FSBB_PtP_AnalyticalVsNgspice_Buck_Region` and
      `_Boost_Region` NRMSE ≤ 0.15
- [ ] `_BuckBoost_Region` NRMSE ≤ 0.20 (relaxed per disclaimer)
- [ ] Region-classification test confirms the three regions and
      hysteresis-band behaviour
- [ ] `Test_FSBB_WorstCase_Inductor_Sizing` confirms
      `L = max(L_buck@VinMax, L_boost@VinMin)` and peak current
      sized for boost @ Vin_min
- [ ] Schema fields documented in
      `MAS/schemas/inputs/topologies/fourSwitchBuckBoost.json`
- [ ] Header docstring includes ASCII art (with all four switches and
      both switching nodes labelled), three-region diagrams,
      equations with citations, references, accuracy disclaimer,
      disambiguation note (4SBB vs SEPIC/Cuk/Zeta vs DAB vs 2-switch
      buck-boost)
- [ ] `generate_ngspice_circuit` netlist parses and contains: four
      MOSFET subcircuits with bootstrap on Q1+Q3, R_snub+C_snub on
      both SW nodes, GEAR, ITL=500/500, `.ic` on Co, region-aware
      gate-drive PULSE
- [ ] `simulate_and_extract_operating_points` actually runs SPICE and
      falls back analytically with `[analytical]` tag on failure;
      netlist region matches analytical-solver region per OP
- [ ] `simulate_and_extract_topology_waveforms` references real
      netlist nodes (grep-verified)
- [ ] Per-OP diagnostic accessors populated and at least one test
      reads them (especially `lastRegion` and the per-switch RMS
      currents)
- [ ] Multi-output **N/A** — single-output. Document.
- [ ] `bidirectional`, `phaseCount`, `controlMode`, `transitionMode`
      flags all tested
- [ ] No divide-by-zero in design-time or SPICE code

ZVS-boundary, rectifier-type CT/CD/FB, and duty-cycle-loss checkboxes
from the Golden Guide are **not applicable** to 4SBB — document with
"N/A — hard-switched single-inductor synchronous PWM" notes in the
test file.

---

## 6. Equations block (citation table)

Use this verbatim in the header docstring. CCM, lossless unless
noted.

### Buck region (Vo/Vin < 1 − hysteresis)

| Quantity | Formula | Source |
|---|---|---|
| Conversion ratio | M = Vo/Vin = D_buck | TI SLVA535B Eq. 1 |
| Duty (lossy) | D_buck = (Vo+Vd) / (Vin·η) | TI SLVA535B |
| Inductor average current | I_L = Iout | TI SLVA535B |
| Inductor ripple | ΔI_L = Vo·(1−D_buck) / (L·Fs) = (Vin−Vo)·D_buck / (L·Fs) | TI SLVA535B Eq. 3 |
| Switch peak current | I_SW(pk) = Iout + ΔI_L/2 | TI SLVA535B |
| L_min (buck) | L = Vo·(Vin_max−Vo) / (K_ind·Fs·Vin_max·Iout), K_ind ∈ [0.2, 0.4] | TI SLVA535B Eq. 5 |
| Worst-case ripple | at **Vin = Vin_max** | TI SLVA535B § 3 |

### Boost region (Vo/Vin > 1 + hysteresis)

| Quantity | Formula | Source |
|---|---|---|
| Conversion ratio | M = Vo/Vin = 1/(1−D_boost) | TI SLVA535B Eq. 2 |
| Duty (lossy) | D_boost = 1 − Vin·η/(Vo+Vd) | TI SLVA535B |
| **Inductor average current** | **I_L = Iout/(1−D_boost) = Iout·Vo/Vin** ← **largest DC bias** | TI SLVA535B Eq. 4 |
| Inductor ripple | ΔI_L = Vin·D_boost / (L·Fs) | TI SLVA535B Eq. 4 |
| Switch peak current | I_SW(pk) = Iout/(1−D_boost) + ΔI_L/2 | TI SLVA535B |
| L_min (boost) | L = Vin_min²·(Vo−Vin_min) / (K_ind·Fs·Iout·Vo²) | TI SLVA535B Eq. 8 |
| Worst-case ripple | at **Vin = Vin_min** (also where saturation is worst) | TI SLVA535B § 3 |

### Buck-boost region (Vin ≈ Vo)

| Quantity | Mode 1 (`simultaneous`) | Mode 2 (`splitPwm`) | Source |
|---|---|---|---|
| Conversion ratio | M = D / (1−D) | M = D_buck / (1−D_boost) | Restrepo IEEE 8986423; LM5176 § 8.3 |
| Duty | D = Vo / (Vin+Vo) | D_buck + D_boost = 1 + δ (δ ≈ 0.05 per LM5176) | LM5176 d.s. § 8.3.4 |
| Inductor avg current | I_L ≈ Iout / (1−D) | I_L ≈ Iout · max(1, Vo/Vin) | Restrepo IEEE 8986423 |
| Inductor ripple | ΔI_L = Vin·D / (L·Fs) (large) | ΔI_L ≈ buck-mode + boost-mode summed (smaller) | LM5176 § 8.3.4 |
| Efficiency | lower (flyback-like) | higher (LM5176/LT8390 default) | EDN article |

### DCM critical inductance per region

| Region | K_crit | Source |
|---|---|---|
| Buck | (1 − D_buck) | Erickson § 5 |
| Boost | D_boost · (1 − D_boost)² | Erickson § 5 |
| Buck-boost (Mode 1) | (1 − D)² | Erickson § 5 |

### Worst-case design (the headline)

| Constraint | Region & condition | Comment |
|---|---|---|
| **Saturation peak** | boost @ Vin_min | I_pk = Iout·Vo/Vin_min + ΔI_L/2 |
| **Copper / RMS** | boost @ Vin_min | RMS ≈ I_L_avg with 5 % triangular term |
| **Core loss** | buck-boost @ Vin ≈ Vo | highest harmonic content |
| **L sizing** | L = max(L_buck@VinMax, L_boost@VinMin) | TI SLVA535B § 3 |

---

## 7. Reference designs to include in tests

| ID | Vin | Vout | Iout | Fsw | Controller | Modes covered | Source |
|---|---|---|---|---|---|---|---|
| TI TIDA-01411 | 9–36 V | 12 V | 8 A | 350 kHz | LM5175 | all 3 | ti.com/tool/TIDA-01411 |
| TI LM5176EVM-HP | 9–55 V | 24 V | 5 A | 200 kHz | LM5176 | all 3 | TI SNVU547 |
| TI LM5175 HD-EVM | 5–42 V | 12 V | 12 A | 300 kHz | LM5175 | all 3 | TI SNVU439/SNVU465 |
| ADI MAXREFDES1255 | 8–36 V | 12 V | 5 A | 440 kHz | MAX17761 | all 3 | analog.com |
| ADI LT8390 demo DC2825A | 4–60 V | 12 V | 7 A | 300 kHz | LT8390 | all 3 (peakBuckPeakBoost) | analog.com |
| Renesas ISL81801EVAL2Z | 32–80 V | 48 V | 5 A | 200 kHz | ISL81801 | mostly buck-boost | renesas.com |
| TI PMP20042 | 5–20 V | 5/9/15/20 V | 5 A | 600 kHz | LM5175 | all 3 (USB-PD) | ti.com/tool/PMP20042 |
| TI BQ25713 / PMP40441 | 5–20 V | 1S–4S Li | up to 8 A | 800 kHz–1.2 MHz | BQ25713 | all 3 (laptop charger) | ti.com |

---

## 8. Future work (out of scope for this plan)

- **Three-level / flying-capacitor 4SBB** — six switches + flying
  cap; different commutation. Research-grade. Defer.
- **Quasi-resonant ZVS extensions** — LT8390 quasi-resonant on the
  boost leg is omitted in v1 (analytical solver assumes
  hard-switching). Adding this would require Coss-aware switch model
  and ringing analysis.
- **Mode-transition glitch modelling** — real controllers blend
  smoothly between regions over a ~5 % Vin band; v1 picks one region
  per OP. A v2 polish item is to model the transient overlap.
- **Inductor saturation derate at mode-edge load steps** — peak
  current at the BUCK ↔ BUCK_BOOST mode boundary during a load step
  can exceed steady-state predictions by ~20 % (per ADI app note).
  v1 reports the steady-state worst case; a v2 item is to add the
  load-step derating factor as an optional `loadStepDerateRatio`
  schema field.
- **MPC controller mode** — Yan et al. IEEE 9219154 propose a Model
  Predictive Controller for smooth transition. Out of scope for v1.

After Track A lands, update
`project_converter_model_quality_tiers.md` to add 4SBB at
DAB-quality tier, and update `FUTURE_TOPOLOGIES.md` to remove 4SBB
from the "in progress" list.

---

## 9. References

- **TI SLVA535B** — *Basic Calculations of a 4-Switch Buck-Boost
  Power Stage* (canonical 4SBB equations)
- **TI LM5175** datasheet — synchronous 4-switch buck-boost
  controller
- **TI LM5176** datasheet — wide-Vin 4-switch buck-boost controller
- **TI LM5177** datasheet — bidirectional 4-switch buck-boost
- **TI SLVAFJ3** — *Layout Optimization of 4-Switch Buck-Boost*
- **TI SNVA826** — *Power Sharing Between Two Parallel 4-Switch
  Buck-Boost Converters* (interleaved)
- **TI SNVU547 / SNVU439 / SNVU465** — LM5176/LM5175 EVM user guides
- **TI TIDA-01411** — 9–36 V → 12 V @ 8 A reference design
- **TI PMP20042** — USB-PD 4SBB reference design
- **TI TPS55288** — USB-PD buck-boost converter
- **TI BQ25710 / BQ25713** — laptop NVDC charger
- **ADI LT8390 / LT8391 / LT8392** datasheets — peak-buck/peak-boost
  controller family
- **ADI MAXREFDES1255** — 8–36 V → 12 V @ 5 A 4SBB reference
- **ADI Analog Dialogue** — *4-Switch Buck-Boost Layout: Single vs
  Dual Hot Loop*
- **ADI Technical Article** — *Alternating Control & Optimizing
  Bandwidth in H-Bridge Buck-Boost*
- **Renesas ISL81801** datasheet (FN8937)
- **Renesas ISL78224 / ISL78226** — 4-phase 12 V/48 V bidirectional
  PWM controllers
- **MPS MPQ4214 / MP4247** datasheets — buck-boost converters
- **Restrepo et al.**, *Development of a 4-Switch Buck-Boost
  Converter: Operation Modes, Control & Experimental Results*,
  IEEE 8986423
- **Yan et al.**, *FSBB Based on MPC With Smooth Mode Transition*,
  IEEE 9219154
- **Springer chapter** — *Three-Mode Peak Current Control Strategy
  for FSBB*
- **Wiley ADC2.223** — *A novel modulation for FSBB to eliminate
  RHPZ*
- **EDN** — *FSBB delivers highest efficiency in buck or boost mode*
- **Erickson & Maksimović**, *Fundamentals of Power Electronics*,
  3rd ed., Ch. 5 (DCM), Ch. 6 (cascade connection); the textbook
  covers 4SBB only as a buck+boost cascade — the full three-region
  equation set is from app-notes (TI/ADI/Renesas), not the
  textbook
- **Wikipedia** — *Buck–boost converter*
