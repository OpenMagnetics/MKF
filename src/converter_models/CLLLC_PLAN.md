# CLLLC Bidirectional Resonant Converter Plan

Companion to `CONVERTER_MODELS_GOLDEN_GUIDE.md`. Adds a CLLLC
(Capacitor-Inductor-Inductor-Inductor-Capacitor) bidirectional
resonant converter to MKF. CLLLC is the **symmetric, fully-bidirectional**
cousin of LLC and CLLC — the resonant tank is symmetric across the
transformer, which makes the forward and reverse gain curves identical.
This is what unlocks true V2G in EV onboard chargers and bidirectional
ESS power conversion.

Distinctive headline: **the resonant tank is symmetric** (Lr2 = Lr1/n²,
Cr2 = Cr1·n²) so the same set of equations applies in both power-flow
directions. The model must support reversing power flow as a control
input, and the magnetic adviser must size two integrated leakage
inductances (Lr1 and Lr2) plus the magnetizing inductance Lm — three
magnetic constraints jointly satisfied by a single planar transformer
with controlled leakage on both sides.

This is a planning document — no code changed. The existing plans
(`LLC_REWORK_PLAN.md`, `CLLC_REWRITE_PLAN.md`, `SEPIC_PLAN.md`,
`ZETA_PLAN.md`, `CUK_PLAN.md`, `WEINBERG_PLAN.md`,
`FOUR_SWITCH_BUCK_BOOST_PLAN.md`, `FUTURE_TOPOLOGIES.md`) are not
touched by this plan.

---

## 0. Status snapshot

| Aspect | Today | After this plan |
|---|---|---|
| `Cllllc.cpp` / `.h` | does not exist | new files, DAB-quality from start |
| Schema | nothing | full `cllllcResonant.json` with `bridgeType`, `powerFlowDirection`, `controlStrategy`, `tankSymmetryRatio` |
| Magnetic-adviser path | no CLLLC support | sizes the planar transformer with **two** controlled leakages (Lr1 + Lr2) plus Lm |
| State variables in solver | 3 (LLC) | 5 (i_Lr1, v_Cr1, i_Lr2, v_Cr2, i_Lm) |
| Operating modes | n/a | LLC's Nielsen 6 modes **mirrored** for bidirectional → effectively 12 modes |
| Power-flow direction | unidirectional only (LLC, CLLC) | **forward + reverse** with a single mode-symmetric solver |
| Control strategy | n/a | `pfm` (variable frequency), `psm` (phase-shift), `hybridPfmPsm` (default per TI TIDM-02002) |
| Tests | none | TI TIDM-02002 (6.6 kW) + TIDM-02013 (7.4 kW) + KIT 20 kW + V2G reverse-direction symmetry test, ZVS boundary, NRMSE, SPICE |

**Important disambiguation up front** (and in the header docstring):
- **LLC**: one Lr, one Cr (primary side only), parallel Lm. Forward only.
- **CLLC**: one Lr (primary), TWO Cr (one each side), parallel Lm.
  Bidirectional but **asymmetric** reverse gain. Already in MKF as
  `Cllc.cpp` (with rewrite plan in `CLLC_REWRITE_PLAN.md`).
- **CLLLC**: TWO Lr (one each side), TWO Cr (one each side), parallel
  Lm. Bidirectional and **symmetric** if Lr2 = Lr1/n² and Cr2 = Cr1·n².
  This plan.

---

## 1. Why add CLLLC

1. **It is the canonical V2G-capable resonant topology.** Every modern
   EV OBC at 6.6 / 11 / 22 kW uses CLLLC because the symmetric tank
   gives identical forward (grid → battery) and reverse (battery →
   grid) gain curves. CLLC's asymmetric reverse gain forces the
   controller into a different operating region in V2G mode, which
   either de-rates the converter or requires a more complex control.
   CLLLC eliminates this asymmetry by construction.
2. **The magnetics are different from LLC and CLLC.** LLC has one
   leakage; CLLC has one leakage; CLLLC has **two** controlled
   leakages (one each side of the transformer). Modern designs
   (Navitas 3-φ CLLLC, Wolfspeed 6.6 kW bi-OBC, IET PEL 2025 Ansari)
   integrate both leakages as planar transformer magnetic shunts on
   each side. The magnetic adviser must size both Lr1 and Lr2
   simultaneously, with the symmetry constraint Lr2 = Lr1/n² as a
   joint optimization target.
3. **Production volume** is in EV OBCs (6.6 / 11 / 22 kW), bidirectional
   ESS DC/DC (5–250 kW), DC fast chargers (50–350 kW modular), and
   server PSUs (3 kW Open Compute / Open Rack bidirectional). CLLLC is
   the dominant topology in this range; LLC and CLLC cover only a
   subset.
4. **`LLC_REWORK_PLAN.md` § 7 already flagged CLLLC as deferred future
   work** for `Cllllc.cpp`, distinct from the existing `Cllc.cpp`. This
   plan is the execution of that flag.

---

## 2. Variant taxonomy (industry survey)

Sources: TI TIDM-02002 (6.6 kW CLLLC C2000 reference), TI Power Tips
135 (CLLLC ESS control), TI TIDM-02013 (7.4 kW OBC), TI PMP41042 test
report, MDPI Electronics 12(7):1605 (PFM+PSM hybrid control), MDPI
Energies 18(14):3815 (V2G fixed-frequency CLLLC), MDPI Appl Sci
10(22):8144 (symmetric isolated bidirectional), Infineon AN-2024-06
(CLLC modelling — applies in part), IEEE 6972179 (CLLLC time-domain
analysis), IEEE 9352422 (CLLLC ASSA — Advanced State-Space Analysis),
IEEE 10362332 (asymmetric CLLLC OBC), Navitas SiC 3-φ CLLLC 20 kW
(2024), KIT 20 kW CLLLC, Wiley CTA 4181 (30 kW V2G CLLLC), IET PEL
2025 Ansari (FEA inserted-shunt planar transformer), Typhoon-HIL
CLLLC reference, EDN Power Tips 135, Virginia Tech CPES V2G thesis.

### Belong as **flags / enums on the CLLLC class** (NOT new models)

| Variant | Mechanism | Why same class |
|---|---|---|
| **Power-flow direction** | `powerFlowDirection: enum {forward, reverse}` per OP, default `forward` | Same hardware; software flips which bridge is PWM-driving and which is synchronous-rectifying. Symmetric tank → identical equations with state-variable swap. |
| **Bridge type** | `bridgeTypePrimary` / `bridgeTypeSecondary: enum {halfBridge, fullBridge}` | Both sides almost always full-bridge (8 switches total) in production; half-bridge variants exist in academic <1 kW prototypes (MDPI Energies 16/16/5928). |
| **Control strategy** | `controlStrategy: enum {pfm, psm, hybridPfmPsm, fixedFrequencyPhaseShift}` | Hybrid PFM+PSM is the industry default (TI TIDM-02002, MDPI 12(7):1605). The model needs this to compute waveforms. |
| **Tank symmetry** | `tankSymmetryRatio: double = 1.0` (a · b ≡ Lr2/(n²·Lr1) · Cr2·n²/Cr1) | 1.0 = symmetric (default). Asymmetric CLLLC (a ≠ 1) optimises one direction (IEEE 10362332); model still works, just unequal fwd/rev gain. |
| **Synchronous rectification** | always on | Universal in CLLLC because bidirectional requires active switches both sides; no diode-rectifier variant. |

### Belong as **separate model classes**

| Variant | Why distinct | Reference |
|---|---|---|
| **Three-phase CLLLC** | 12 switches per side (3 legs each), 3-φ transformer with integrated leakage on each phase. Distinct magnetics. | Navitas SiC 3-φ CLLLC 20 kW (2024); IEEE 10792930 |
| **Modular / interleaved CLLLC** | N parallel cells with phase-shifted carriers; relevant > 50 kW. Per-cell magnetics same as basic CLLLC, but cell-balancing analysis is its own model. | DC fast charger architecture papers |

### Skip outright

- **Asymmetric CLLLC as separate class** — it's the same hardware with
  `tankSymmetryRatio ≠ 1.0`; covered by the flag. (IEEE 10362332.)
- **CLLC as separate class** — already covered by `Cllc.cpp`. CLLLC ≠
  CLLC.
- **Diode-rectifier CLLLC** — does not exist in production
  (bidirectional excludes it).

### Direct user-question equivalents

If a future request asks *"is V2G CLLLC missing?"* — answer: it's not
a separate class, it's `powerFlowDirection: reverse` on the existing
CLLLC class. The symmetric tank means the equations are identical,
just with primary and secondary state variables swapped.

---

## 3. Existing patterns we will reuse

### `Llc.cpp` as the closest precedent (and `Cllc.cpp` for the second-Cr template)

The existing `Llc.cpp` already implements:
- Nielsen-style time-domain piecewise solver (canonical in MKF)
- Half-wave antisymmetry steady-state convergence
- Mode 1–6 classification via `classify_mode()`
- LIP (Load Independent Point) diagnostic
- HB and FB primary bridge support
- Multi-secondary K-matrix in SPICE

`Cllc.cpp` extends this with the second resonant capacitor on the
secondary side. CLLLC extends `Cllc.cpp` further by adding the second
resonant inductor and making the entire model symmetric and
bidirectional.

**Recommended approach**: do NOT inherit from `Cllc` (it's mid-rework
per `CLLC_REWRITE_PLAN.md` and the inheritance would couple the two
plans). Instead, implement `Cllllc` as a sibling class at the same
level — copy and adapt the Nielsen solver from `Llc.cpp` (which is
the gold standard among MKF resonant models), extend the state vector
from 3 to 5, and apply the symmetric-tank invariants explicitly.

### DAB.cpp as the bidirectional-bridge precedent

`Dab.cpp` is the only existing MKF model with two active bridges and
power-flow direction as a control input. The CLLLC SPICE wrapper
should reuse:
- `Dab.cpp` per-winding sense-source pattern
- `Dab.cpp` multi-secondary K-matrix discipline (CLLLC has two
  K-coupling pairs if Lr1/Lr2 are integrated as transformer leakage)
- `Dab.cpp` bidirectional rectification: at any instant, the
  "secondary" bridge is in synchronous-rectifier mode following the
  current zero-crossing, vs the "primary" bridge driving PWM. Software
  flips this on direction reversal.

### LLC rework adds prerequisites

The LLC rework (`LLC_REWORK_PLAN.md`) introduces ZVS-margin diagnostics
(`lastZvsMarginLagging`, `lastZvsLoadThreshold`,
`lastResonantTransitionTime`) that CLLLC must populate **twice** (once
for each bridge). Land the LLC rework first if possible; otherwise
land both ZVS-diagnostic patterns at once (LLC + CLLLC) since the math
is identical.

### Schema convention (camelCase)

The new file is `MAS/schemas/inputs/topologies/cllllcResonant.json`,
sibling to `llcResonant.json` and `cllcResonant.json`.

---

## 4. Track A — Basic CLLLC (the only track in this plan)

### A.1 Files to create

```
src/converter_models/Cllllc.h
src/converter_models/Cllllc.cpp
tests/TestTopologyCllllc.cpp
MAS/schemas/inputs/topologies/cllllcResonant.json
```

Plus add `cllllcResonant` to the topology enum in
`MAS/schemas/inputs/designRequirements.json` (verify camelCase matches
the existing `llcResonant` and `cllcResonant`).

### A.2 Schema (`cllllcResonant.json`)

```jsonc
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://psma.com/mas/inputs/topologies/cllllcResonant.json",
  "title": "cllllcResonant",
  "description": "The description of a CLLLC bidirectional resonant converter excitation",
  "type": "object",
  "properties": {
    "highVoltageBusVoltage": {
      "description": "The HV-side bus voltage (typically 380–800 V for OBC)",
      "$ref": "../../utils.json#/$defs/dimensionWithTolerance"
    },
    "lowVoltageBusVoltage": {
      "description": "The LV-side bus voltage (typically 200–500 V for OBC)",
      "$ref": "../../utils.json#/$defs/dimensionWithTolerance"
    },
    "minSwitchingFrequency": { "type": "number" },
    "maxSwitchingFrequency": { "type": "number" },
    "primaryResonantFrequency": {
      "description": "f_r1 = 1/(2π√(Lr1·Cr1)). Optional; computed if absent.",
      "type": "number"
    },
    "inductanceRatioK": {
      "description": "K = Lm/Lr1, the magnetizing-to-series-inductance ratio",
      "type": "number"
    },
    "qualityFactor": { "type": "number", "default": 0.4 },
    "tankSymmetryRatio": {
      "description": "(Lr2/(n²·Lr1)) · (Cr2·n²/Cr1). 1.0 = perfectly symmetric.",
      "type": "number", "default": 1.0
    },
    "efficiency": { "type": "number", "default": 0.97 },
    "bridgeTypePrimary": {
      "description": "HV-side bridge",
      "$ref": "#/$defs/cllllcBridgeType",
      "default": "fullBridge"
    },
    "bridgeTypeSecondary": {
      "description": "LV-side bridge",
      "$ref": "#/$defs/cllllcBridgeType",
      "default": "fullBridge"
    },
    "controlStrategy": {
      "description": "How the controller modulates frequency vs phase-shift",
      "$ref": "#/$defs/cllllcControlStrategy",
      "default": "hybridPfmPsm"
    },
    "integratedResonantInductors": {
      "description": "If true, Lr1 and Lr2 are integrated as transformer leakage on each side",
      "type": "boolean", "default": true
    },
    "operatingPoints": {
      "type": "array",
      "items": { "$ref": "#/$defs/cllllcOperatingPoint" },
      "minItems": 1
    }
  },
  "required": [
    "highVoltageBusVoltage", "lowVoltageBusVoltage",
    "minSwitchingFrequency", "maxSwitchingFrequency", "operatingPoints"
  ],
  "$defs": {
    "cllllcBridgeType": {
      "title": "cllllcBridgeType",
      "type": "string",
      "enum": ["halfBridge", "fullBridge"]
    },
    "cllllcControlStrategy": {
      "title": "cllllcControlStrategy",
      "type": "string",
      "enum": ["pfm", "psm", "hybridPfmPsm", "fixedFrequencyPhaseShift"]
    },
    "cllllcPowerFlowDirection": {
      "title": "cllllcPowerFlowDirection",
      "type": "string",
      "enum": ["forward", "reverse"]
    },
    "cllllcOperatingPoint": {
      "title": "cllllcOperatingPoint",
      "allOf": [
        { "$ref": "../../utils.json#/$defs/baseOperatingPoint" },
        {
          "properties": {
            "powerFlowDirection": {
              "$ref": "#/$defs/cllllcPowerFlowDirection",
              "default": "forward"
            },
            "phaseShiftDegrees": {
              "description": "Inter-bridge phase shift; only used by psm or hybridPfmPsm",
              "type": "number"
            }
          }
        }
      ]
    }
  }
}
```

After the schema edit, regenerate MAS so `MAS::CllllcResonant` becomes
a generated C++ class. Per `project_mas_hpp_staging.md`, copy the
freshly built `MAS.hpp` into `build/MAS/MAS.hpp` if the build
complains.

### A.3 Class structure (`Cllllc.h`)

```cpp
class Cllllc : public MAS::CllllcResonant, public Topology {
public:
    enum class PowerFlowDirection { FORWARD, REVERSE };
    enum class OperatingRegion { ABOVE_RESONANCE, AT_RESONANCE, BELOW_RESONANCE };

private:
    // 2.1 Simulation tuning
    int numPeriodsToExtract    = 5;
    int numSteadyStatePeriods  = 4;   // larger than LLC (5 state vars converge slower)

    // 2.2 Computed design values (filled by process_design_requirements)
    double computedPrimarySeriesInductance     = 0;   // Lr1
    double computedSecondarySeriesInductance   = 0;   // Lr2 = Lr1 / n²  (symmetric)
    double computedPrimaryResonantCapacitance  = 0;   // Cr1
    double computedSecondaryResonantCapacitance= 0;   // Cr2 = Cr1 · n²  (symmetric)
    double computedMagnetizingInductance       = 0;   // Lm
    double computedTurnsRatio                  = 0;   // n = Np/Ns
    double computedDeadTime                    = 200e-9;
    double computedInductanceRatioK            = 0;   // K = Lm / Lr1
    double computedQualityFactor               = 0;
    double computedPrimaryResonantFrequency    = 0;   // f_r1
    double computedMagnetizingResonantFreqFwd  = 0;   // f_m1
    double computedMagnetizingResonantFreqRev  = 0;   // f_m2 (= f_m1 if symmetric)

    // 2.3 Schema-less device parameters
    double mosfetCossPrimary   = 200e-12;
    double mosfetCossSecondary = 400e-12;             // typically larger LV-side device

    // 2.4 Per-OP diagnostics
    mutable PowerFlowDirection lastPowerFlowDirection = PowerFlowDirection::FORWARD;
    mutable OperatingRegion lastOperatingRegion       = OperatingRegion::AT_RESONANCE;
    mutable int lastModeForward                       = 1;    // Nielsen 1–6
    mutable int lastModeReverse                       = 1;    // Nielsen 1–6, mirrored

    // ZVS diagnostics — populated for BOTH bridges
    mutable double lastZvsMarginPrimaryLagging   = 0.0;
    mutable double lastZvsMarginSecondaryLagging = 0.0;
    mutable double lastZvsLoadThresholdPrimary   = 0.0;
    mutable double lastZvsLoadThresholdSecondary = 0.0;
    mutable double lastResonantTransitionTime    = 0.0;

    mutable double lastPrimaryPeakCurrent       = 0.0;
    mutable double lastSecondaryPeakCurrent     = 0.0;
    mutable double lastPrimaryRmsCurrent        = 0.0;
    mutable double lastSecondaryRmsCurrent      = 0.0;
    mutable double lastMagnetizingPeakCurrent   = 0.0;
    mutable double lastCr1PeakVoltage           = 0.0;
    mutable double lastCr2PeakVoltage           = 0.0;
    mutable double lastCurrentSharingRatio      = 1.0;        // i_Lr1 / (i_Lr2/n) — should be 1.0 if symmetric
    mutable double lastSteadyStateResidual      = 0.0;

    // 2.5 Extra-component waveforms (cleared in process_operating_points)
    mutable std::vector<Waveform> extraLr1CurrentWaveforms;
    mutable std::vector<Waveform> extraLr1VoltageWaveforms;
    mutable std::vector<Waveform> extraLr2CurrentWaveforms;
    mutable std::vector<Waveform> extraLr2VoltageWaveforms;
    mutable std::vector<Waveform> extraCr1CurrentWaveforms;
    mutable std::vector<Waveform> extraCr1VoltageWaveforms;
    mutable std::vector<Waveform> extraCr2CurrentWaveforms;
    mutable std::vector<Waveform> extraCr2VoltageWaveforms;
    mutable std::vector<Waveform> extraLmCurrentWaveforms;
    mutable std::vector<Waveform> extraTimeVectors;

public:
    bool _assertErrors = false;

    Cllllc(const json& j);
    Cllllc() {}

    // Tuning, computed, last_*, mosfet_coss accessors per Golden
    // Guide § 2.7–2.10.

    PowerFlowDirection get_last_power_flow_direction() const { return lastPowerFlowDirection; }
    OperatingRegion    get_last_operating_region() const     { return lastOperatingRegion; }

    // Overrides
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const CllllcResonantOperatingPoint& op,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // Static analytical helpers
    static double compute_primary_resonant_frequency(double Lr1, double Cr1);
    static double compute_magnetizing_resonant_frequency_fwd(double Lr1, double Lm, double Cr1);
    static double compute_magnetizing_resonant_frequency_rev(double Lr2, double Lm, double n, double Cr2);
    static double compute_fha_gain(double f_n, double K, double Q, double a, double b);
    static double compute_inductance_ratio_K(double Lm, double Lr1);
    static double compute_quality_factor(double Lr1, double Cr1, double R_load_reflected);
    static double compute_zvs_lm_max(double T_dead, double C_oss_pri, double Fs);
    static double compute_turns_ratio(double V_HV, double V_LV);
    static double compute_symmetric_lr2(double Lr1, double n);   // = Lr1 / n²
    static double compute_symmetric_cr2(double Cr1, double n);   // = Cr1 · n²

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

class AdvancedCllllc : public Cllllc {
    // bypasses process_design_requirements; user provides Lr1, Lr2, Cr1, Cr2, Lm directly
};
```

### A.4 Header docstring

Mandatory blocks (per Golden Guide § 3) before the class:

1. **TOPOLOGY OVERVIEW** with ASCII art:

   ```
        HV-side full-bridge        Tank        LV-side full-bridge
        (4 SiC/GaN)                            (4 SiC/GaN)
   
       Q1 ──┬── Q3       Lr1   Cr1     Cr2  Lr2       Q5 ──┬── Q7
            │                  ╔═══╗                       │
       HV+──┤                  ║   ║                       ├──LV+
            │             ┌────╣Lm ╠────┐                  │
       HV-──┤            (n)   ║   ║   (1)                 ├──LV-
            │                  ╚═══╝                       │
       Q2 ──┴── Q4    L───C───            ───C───L         Q6 ──┴── Q8
                       (Lr1=transformer leakage)
                       (Lr2=transformer leakage / n²)
   ```

   Symmetric design constraint: **Lr2 = Lr1 / n²**, **Cr2 = Cr1 · n²**.
   Both Lr's are typically integrated as planar transformer leakage
   (IET PEL 2025 Ansari, Navitas 3-φ CLLLC).

2. **MODULATION CONVENTION** —
   - Forward direction: HV bridge drives PWM at variable Fs; LV bridge
     is in synchronous-rectifier mode (gate signals follow current
     zero-crossing).
   - Reverse direction: roles flip — LV bridge drives PWM, HV bridge
     synchronous-rectifies.
   - Hybrid PFM+PSM: PFM near f_r1 for high efficiency; PSM at the
     edges of the gain range for wide regulation.
   - Phase shift φ ∈ [0, π/2]: between leading and lagging legs of the
     active bridge (intra-bridge PSM, LM5176-style — not DAB-style
     inter-bridge phase shift, which is a separate variant).

3. **KEY EQUATIONS** — use the table in § 6 below verbatim with
   `[TI TIDM-02002]`, `[MDPI 12(7):1605]`, `[IEEE 6972179]` style
   inline citations.

4. **References** — TI TIDM-02002 / TIDM-02013, MDPI Electronics
   12(7):1605, MDPI Energies 18(14):3815, IEEE 6972179 (time-domain),
   IEEE 9352422 (ASSA), IET PEL 2025 Ansari, Navitas SiC 3-φ CLLLC.

5. **Accuracy disclaimer** mirroring `Dab.cpp:296–320`:
   - FHA gain accuracy degrades farther from f_r1 than for LLC because
     two LC poles interact (forward and reverse). Time-domain solver
     is the source of truth; FHA is for design synthesis only.
   - Coss(Vds) variation not modelled in v1; SPICE uses linear C.
   - Body-diode reverse recovery omitted.
   - Direction-reversal mode transition near gain = 1 is numerically
     tricky; v1 picks one direction per OP from the `powerFlowDirection`
     field and does NOT model the reversal transient.
   - Transformer parasitic capacitance (Cps inter-winding) not in v1;
     becomes relevant > 500 kHz.
   - Magnetic shunt-leakage tuning assumed ideal: Lr1 and Lr2 reach
     their target values exactly. Real planar transformers have ±10 %
     tolerance on each leakage independently.

6. **Disambiguation note**:
   - **CLLLC** vs **LLC** — LLC has one Lr (primary side) and one Cr;
     forward only.
   - **CLLLC** vs **CLLC** — CLLC has one Lr and TWO Cr; bidirectional
     but **asymmetric** reverse gain. Already in MKF as `Cllc.cpp`.
   - **CLLLC** vs **DAB** — DAB has no resonant tank, just transformer
     leakage; CLLLC is resonant. Both bidirectional, both 8 switches.
   - **CLLLC** vs **SRC / PRC** — SRC is series-resonant only (Lr+Cr,
     no Lm parallel branch); PRC is parallel; CLLLC has both plus a
     mirrored secondary tank.

### A.5 Analytical solver (`process_operating_point_for_input_voltage`)

This is the **most complex** section because of the 5-state symmetric
solver and the bidirectional control. The algorithm:

1. **Read direction** from the OP: `powerFlowDirection ∈ {forward, reverse}`.
2. **Reflect everything to the active side**. If direction is forward,
   work in the HV-side reference frame: Lr_active = Lr1, Cr_active =
   Cr1, Lr_passive = Lr2 · n², Cr_passive = Cr2 / n², V_in = V_HV,
   V_out = V_LV · n. If reverse, swap and work in the LV-side frame.
   This trick reduces the bidirectional model to a unidirectional
   solver with a "passive-side reflected" tank.
3. **Build sub-intervals** via a Nielsen-style time-domain propagator,
   extended to 5 state variables: i_Lr_active, v_Cr_active,
   i_Lr_passive (reflected), v_Cr_passive (reflected), i_Lm.
4. **Within each sub-interval**, all five variables evolve as a
   coupled LCLCL system. Use the closed-form modal solution:
   - Decompose the 5×5 state matrix into eigenmodes (typically 2
     resonant pairs + 1 magnetizing pole).
   - Each sub-interval is exact-sinusoid for the resonant pairs and
     exact-linear for the magnetizing pole.
   - Propagate all five variables by `expm(A·dt)` on the modal-
     decomposed form (much cheaper than re-evaluating expm every
     step).
5. **Steady-state convergence** by half-wave antisymmetry:
   - i_Lr_active(T/2) = −i_Lr_active(0)
   - v_Cr_active(T/2) = −v_Cr_active(0)
   - same for the passive-side reflected variables and i_Lm
   - Newton iteration (5×5 Jacobian) — slower than LLC's 3×3 but the
     same algorithm.
6. **Mode classification**: extend `classify_mode()` from `Llc.cpp` to
   the symmetric tank. Forward and reverse modes are mirrored: if
   forward mode is 3, reverse mode in the same OP would be 3' (the
   mirror). v1: report `lastModeForward` and `lastModeReverse`
   separately; the active direction's mode matches the analytical
   solver's chosen mode.
7. **Operating-region classification**:
   - `f_s > f_r1` → ABOVE_RESONANCE
   - `f_s = f_r1 ± 5 %` → AT_RESONANCE
   - `f_m1 < f_s < f_r1` → BELOW_RESONANCE
8. **ZVS check on BOTH bridges** (the headline new diagnostic):
   - Primary (active bridge): `lastZvsMarginPrimaryLagging` =
     i_Lm_at_switch − i_Lr_active_reflected_at_switch.
   - Secondary (synchronous-rectifying bridge): same formula but
     evaluated at the LV-side switching event.
   - Both must be > 0 for full soft-switching. Lm sizing constraint:
     `Lm ≤ T_dead / (16 · C_oss · f_s)`.
9. **Symmetric-tank current-sharing diagnostic**:
   `lastCurrentSharingRatio` = max(i_Lr1) / (max(i_Lr2) / n).
   Should be 1.0 ± 5 % for symmetric design; deviations flag a
   `tankSymmetryRatio` mismatch.
10. **Sample grid**: `2·N + 1` with N = 256 (default). Use N = 512
    for above-resonance with high f_s/f_r1 ratio (CLLLC current is
    near-sinusoidal, finer grid helps the SPICE NRMSE comparison).
11. **Populate diagnostics** at the end (Golden Guide § 4.6 + CLLLC
    extras):

    ```
    lastPowerFlowDirection      = direction
    lastOperatingRegion         = region per Fs vs (f_m1, f_r1)
    lastModeForward / Reverse   = classify_mode(...)
    lastZvsMarginPrimaryLagging = computed (per Huang SLUP263 § 6 + CLLLC mirror)
    lastZvsMarginSecondaryLagging = computed
    lastZvsLoadThresholdPrimary / Secondary = computed
    lastResonantTransitionTime  = T_dead or Coss-resonance time
    lastPrimaryPeakCurrent      = max(|i_Lr_active|)
    lastSecondaryPeakCurrent    = max(|i_Lr_passive_reflected|) / n
    lastCr1PeakVoltage          = max(|v_Cr1|)
    lastCr2PeakVoltage          = max(|v_Cr2|)
    lastCurrentSharingRatio     = i_Lr1 vs i_Lr2/n ratio
    lastSteadyStateResidual     = Newton final residual
    ```

### A.6 SPICE netlist (`generate_ngspice_circuit`)

Mandatory checkpoints from Golden Guide § 5:

- **Eight active switches**: Q1–Q4 (HV bridge) and Q5–Q8 (LV bridge).
  Use `.model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg` for both sides
  (or per-side models if Coss differs significantly).
- **Two K-coupling pairs** if Lr1 and Lr2 are integrated as transformer
  leakage:
  - `K1 L_pri L_lr1 0.998` — primary winding to its leakage
  - `K2 L_sec L_lr2 0.998` — secondary winding to its leakage
  - `K_pri_sec L_pri L_sec <coupling>` — main transformer coupling
  Without all three, ngspice treats Lr1 and Lr2 as independent and
  produces nonsensical resonance.
- **R_snub 1k + C_snub 1n across each switch** (8 snubbers total).
- Per-component sense sources: `V_pri_bridge_sense`,
  `V_sec_bridge_sense`, `V_Lr1_sense`, `V_Lr2_sense`, `V_Cr1_sense`,
  `V_Cr2_sense`, `V_Lm_sense`, `V_Vout_sense`.
- `.ic` pre-charge on all output caps and on `v_Cr1` and `v_Cr2` to
  their steady-state values (Cr1 charges to roughly half of V_HV·D
  per cycle; Cr2 mirrored).
- `.options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500`
- `.options METHOD=GEAR TRTOL=7`
- **`R_load = max(V_LV / I_LV, 1e-3)`** — divide-by-zero guard on the
  passive (output) side. In reverse direction, the guard applies on
  V_HV / I_HV instead.
- `stepTime = period / 500` (finer than LLC's `period / 200` because
  CLLLC has higher-frequency resonant content above f_r1).
- Saved signals: `.save v(hv_pos) v(lv_pos) v(vab_pri) v(vab_sec)
  v(vcr1) v(vcr2) i(V_pri_bridge_sense) i(V_sec_bridge_sense)
  i(V_Lr1_sense) i(V_Lr2_sense) i(V_Lm_sense)`.
- Comment header with topology, V_HV, V_LV, Fs, f_r1,
  `powerFlowDirection`, `controlStrategy`, `tankSymmetryRatio`,
  "generated by OpenMagnetics".

**Direction-aware gate-drive emission**:

- In `forward` direction: HV-side Q1–Q4 driven with PWM at Fs, with
  φ phase shift if `controlStrategy ∈ {psm, hybridPfmPsm}`. LV-side
  Q5–Q8 driven as synchronous rectifier, gates following the current
  zero-crossing of i_Lr2.
- In `reverse` direction: roles flip.
- The synchronous-rectifying side gate signals are derived from the
  analytical solver's predicted i_Lr2 zero-crossing times. v1: emit
  these as fixed `PULSE` statements at the predicted switch times
  (open-loop). A v2 polish item is to make the SR gates respond to
  the actual SPICE-simulated current zero-crossing (closed-loop).

### A.7 `simulate_and_extract_operating_points`

Follow the DAB pattern (`Dab.cpp:1218–1307`) verbatim, with the
analytical fallback tagged `[analytical]` on failure. Per OP, the
netlist must reflect the analytical solver's chosen direction —
forward netlist for forward OP, reverse netlist for reverse OP.

### A.8 `simulate_and_extract_topology_waveforms`

Same hygiene rules as the SEPIC / Zeta / 4SBB plans: grep the netlist
for exact node names before referencing.

### A.9 `get_extra_components_inputs`

Returns a richer set than LLC (which has Lr + Cr only):

- `Inputs` for the **transformer with two integrated leakages**
  (the canonical case when `integratedResonantInductors == true`).
  This is a 2-winding transformer where the magnetic adviser is asked
  to deliver a primary-side leakage Lr1 and a secondary-side leakage
  Lr2 = Lr1/n² simultaneously. Both leakages are first-class design
  parameters; geometry must satisfy both jointly.
- `Inputs` for **Lr1 (discrete)** and **Lr2 (discrete)** when
  `integratedResonantInductors == false` — two separate inductors
  outside the transformer.
- `CAS::Inputs` for **Cr1** (HV-side resonant cap): voltage waveform
  ±V_Cr1_pk, RMS current per § 6 equation. Voltage rating margin
  1.5× peak for safety (CLLLC Cr can ring during direction-reversal
  transients).
- `CAS::Inputs` for **Cr2** (LV-side resonant cap): mirrored. Higher
  RMS current at low V_LV (`I_Cr2_rms ≈ I_Lr2_rms`).
- `CAS::Inputs` for the **HV bus cap** when `mode = REAL`.
- `CAS::Inputs` for the **LV bus cap** when `mode = REAL`.

### A.10 Tests (`TestTopologyCllllc.cpp`)

Mandatory cases per Golden Guide § 8 + CLLLC-specific tests:

1. **`Test_Cllllc_TIDM02002_Reference_Design`** — TI 6.6 kW reference
   (380–600 V ↔ 280–450 V, 300–700 kHz). Reproduce within 5 %.
2. **`Test_Cllllc_TIDM02013_Reference_Design`** — TI 7.4 kW OBC.
3. **`Test_Cllllc_KIT_20kW_Reference_Design`** — KIT 20 kW academic
   reference (400–800 V ↔ 200–500 V, 100–300 kHz).
4. **`Test_Cllllc_Forward_Reverse_Symmetry`** — same OP run with
   `powerFlowDirection=forward` and `=reverse`; verify forward
   inductor-current waveform mirrors reverse (within 1 %) when
   `tankSymmetryRatio = 1.0`. **THIS IS THE HEADLINE TEST.**
5. **`Test_Cllllc_Asymmetric_Tank`** — set `tankSymmetryRatio = 0.8`
   and verify forward/reverse gain differs by the predicted FHA
   amount (per IEEE 10362332). Verifies the symmetric assumption is
   not silently baked into the solver.
6. **`Test_Cllllc_Design_Requirements`** — turns ratio, Lr1, Lr2,
   Cr1, Cr2, Lm all positive; symmetric design returns
   Lr2 ≈ Lr1/n², Cr2 ≈ Cr1·n². Round-trip a power target.
7. **`Test_Cllllc_OperatingPoints_Generation`** — multiple Fs values
   spanning above/at/below resonance; antisymmetry holds; current is
   piecewise sinusoidal (5-state Nielsen modes).
8. **`Test_Cllllc_Operating_Modes_Forward`** — Nielsen 6 modes in
   forward direction.
9. **`Test_Cllllc_Operating_Modes_Reverse`** — same 6 modes mirrored
   in reverse direction.
10. **`Test_Cllllc_ZVS_Boundaries_Both_Bridges`** — sweep load and
    verify BOTH `lastZvsMarginPrimaryLagging` and
    `lastZvsMarginSecondaryLagging` cross zero at the predicted
    threshold. **This is the bidirectional ZVS check that LLC and
    CLLC don't have.**
11. **`Test_Cllllc_Lm_ZVS_Sizing`** — verify that the design solver
    picks Lm satisfying `Lm ≤ T_dead / (16·C_oss·f_s)` per Huang
    SLUP263 + Infineon AN-2024-06.
12. **`Test_Cllllc_ControlStrategy_Variants`** — same OP run with
    `pfm`, `psm`, `hybridPfmPsm`, `fixedFrequencyPhaseShift` and
    verify the corresponding switch waveforms (only Fs varies in PFM,
    only φ varies in PSM, both vary in hybrid, only φ in FF-PS).
13. **`Test_Cllllc_BridgeType_Variants`** — full-bridge × full-bridge
    (default), half-bridge primary × full-bridge secondary, etc.
    Verify gain factor = 0.5 for each half-bridge side.
14. **`Test_Cllllc_LIP_Diagnostic`** — load-independent point
    detection (mirrors LLC's LIP test). For symmetric CLLLC, the LIP
    is at `f_s = f_r1` exactly.
15. **`Test_Cllllc_CurrentSharing_Symmetric`** — verify
    `lastCurrentSharingRatio` ∈ [0.95, 1.05] when `tankSymmetryRatio
    = 1.0`. Verify it deviates predictably when ratio ≠ 1.0.
16. **`Test_Cllllc_SPICE_Netlist`** — netlist parses, contains 8
    snubbers, three K-coupling statements, `.options METHOD=GEAR`,
    ITL=500/500, `.ic` on Cr1 + Cr2 + bus caps,
    `R_load = max(V/I, 1e-3)`.
17. **`Test_Cllllc_SPICE_Netlist_Forward_Reverse`** — generate netlist
    for forward and reverse OPs; assert the gate-drive PULSE
    statements differ correctly (HV bridge active in forward, LV
    bridge active in reverse).
18. **`Test_Cllllc_PtP_AnalyticalVsNgspice_Forward`** — primary-current
    NRMSE ≤ 0.20 (relaxed from LLC's 0.15 per the disclaimer about
    FHA accuracy), period-averaged V_LV within 5 %.
19. **`Test_Cllllc_PtP_AnalyticalVsNgspice_Reverse`** — same targets
    in reverse.
20. **`Test_Cllllc_Waveform_Plotting`** — CSV dump under
    `tests/output/cllllc/<direction>/` per direction.
21. **`Test_AdvancedCllllc_Process`** — round-trip
    `AdvancedCllllc::process()` returns sane Inputs.
22. **`Test_Cllllc_IntegratedLeakage_vs_Discrete`** — compare
    `integratedResonantInductors=true` (single transformer with two
    leakage winding inputs) vs `=false` (two discrete inductors); both
    produce equivalent steady-state operating points.

ZVS-boundary tests are **mandatory and headline** — CLLLC is the
first MKF model where TWO ZVS conditions must be simultaneously
satisfied. This test (10) is the audit-evidence that the bidirectional
soft-switching analysis works.

Tests (4), (5), (10), (15), and (17) are the headline CLLLC-specific
tests that justify a separate model from LLC and CLLC. They prove the
symmetric-tank bidirectional solver works correctly.

---

## 5. Class checklist (Golden-Guide § 14 acceptance)

- [ ] Builds clean (`ninja -j2 MKF_tests` returns 0)
- [ ] All `TestTopologyCllllc` cases pass (≥ 22 cases per § A.10)
- [ ] `Test_Cllllc_PtP_AnalyticalVsNgspice_Forward` and `_Reverse`
      NRMSE ≤ 0.20
- [ ] `Test_Cllllc_Forward_Reverse_Symmetry` confirms identical
      waveforms in symmetric mode
- [ ] `Test_Cllllc_ZVS_Boundaries_Both_Bridges` confirms BOTH bridges'
      ZVS-margin diagnostics change sign at predicted thresholds
- [ ] Schema fields documented in
      `MAS/schemas/inputs/topologies/cllllcResonant.json`
- [ ] Header docstring includes ASCII art (with all 8 switches and
      symmetric tank labelled), equations with citations, references,
      accuracy disclaimer, disambiguation note (CLLLC vs LLC vs CLLC
      vs DAB vs SRC)
- [ ] `generate_ngspice_circuit` netlist parses and contains: 8
      MOSFETs with snubbers, three K-coupling statements (K1, K2,
      K_pri_sec), GEAR, ITL=500/500, `.ic` on Cr1+Cr2+bus caps,
      direction-aware gate-drive PULSE
- [ ] `simulate_and_extract_operating_points` actually runs SPICE and
      falls back analytically with `[analytical]` tag on failure;
      netlist direction matches OP's `powerFlowDirection`
- [ ] `simulate_and_extract_topology_waveforms` references real
      netlist nodes (grep-verified)
- [ ] Per-OP diagnostic accessors populated and at least one test
      reads them (especially `lastZvsMarginPrimary/Secondary` and
      `lastCurrentSharingRatio`)
- [ ] Multi-output **N/A** — single primary + single secondary
      (CLLLC for multi-output is not common in industry). Document.
- [ ] `powerFlowDirection`, `controlStrategy`, `bridgeTypePrimary/
      Secondary`, `tankSymmetryRatio`, `integratedResonantInductors`
      all tested
- [ ] No divide-by-zero in design-time or SPICE code

ZVS-boundary checkbox is **the headline acceptance** for this model.
Rectifier-type CT/CD/FB and duty-cycle-loss checkboxes from the
Golden Guide are N/A — document with "N/A — synchronous-rectifier
resonant" notes in the test file.

---

## 6. Equations block (citation table)

Use this verbatim in the header docstring. Symmetric design unless
noted.

### Resonant frequencies

| Quantity | Formula | Source |
|---|---|---|
| Primary series resonance | f_r1 = 1 / (2π·√(Lr1·Cr1)) | TI TIDM-02002 |
| Secondary series resonance | f_r2 = 1 / (2π·√(Lr2·Cr2)) = f_r1 (sym) | MDPI 12(7):1605 |
| Magnetizing-loaded fwd | f_m1 = 1 / (2π·√((Lr1+Lm)·Cr1)) | MDPI 12(7):1605 |
| Magnetizing-loaded rev | f_m2 = 1 / (2π·√((Lr2+Lm/n²)·Cr2)) = f_m1 (sym) | MDPI 12(7):1605 |
| LIP (load-independent point) | f_LIP = f_r1 (symmetric) | TI Power Tips 135 |

### FHA gain (normalized)

| Symbol | Definition | Source |
|---|---|---|
| n | turns ratio N_p/N_s | — |
| K | Lm / Lr1 | Infineon AN-2024-06 |
| Q | (1/R_load_reflected)·√(Lr1/Cr1) | MDPI 12(7):1605 |
| a | Lr2 / (n²·Lr1); = 1 if symmetric | IEEE 10362332 |
| b | Cr2·n² / Cr1; = 1 if symmetric | IEEE 10362332 |
| f_n | f_s / f_r1 | — |

Gain expression: M(f_n, K, Q, a, b) — see MDPI 12(7):1605 § 3 for the
closed form. In symmetric design (a = b = 1) it reduces to:

```
M(f_n, K, Q) = 1 / sqrt( (1 + K·(1 − 1/f_n²))² + Q²·(f_n − 1/f_n)² )
```

(structurally similar to LLC FHA but with both Lr's contributing).

### Time-domain state vector (5 variables)

| Variable | Description |
|---|---|
| i_Lr1 | primary-side resonant inductor current |
| v_Cr1 | primary-side resonant cap voltage |
| i_Lr2 | secondary-side resonant inductor current |
| v_Cr2 | secondary-side resonant cap voltage |
| i_Lm | magnetizing current |

### Design constraints

| Quantity | Formula | Source |
|---|---|---|
| Symmetric Lr2 | Lr2 = Lr1 / n² | MDPI 12(7):1605 |
| Symmetric Cr2 | Cr2 = Cr1 · n² | MDPI 12(7):1605 |
| ZVS Lm max (both bridges) | Lm ≤ T_dead / (16·C_oss·f_s) | Huang SLUP263 § 6; Infineon AN-2024-06 |
| Cr1 RMS current | I_Cr1_rms = I_Lr1_rms (series) | TI TIDM-02002 |
| Cr2 RMS current | I_Cr2_rms = I_Lr2_rms = n · I_Cr1_rms (sym) | MDPI 12(7):1605 |
| Cr1 peak voltage | V_Cr1_pk = V_HV/2 + I_Lr1_pk·X_Cr1 (above res) | MDPI 12(7):1605 |

### Soft-switching boundaries

| Boundary | Formula | Source |
|---|---|---|
| Primary ZVS | i_Lm(t_switch) > i_Lr1_load_reflected(t_switch) | Huang SLUP263 (extended for CLLLC) |
| Secondary ZVS | i_Lm·n²/n(t_switch_LV) > i_Lr2_load(t_switch_LV) | mirrored — IEEE 6972179 |

---

## 7. Reference designs to include in tests

| ID | Power | V_HV | V_LV | Fs | Controller | Source |
|---|---|---|---|---|---|---|
| TI TIDM-02002 / PMP41042 | 6.6 kW | 380–600 V | 280–450 V | 300–700 kHz (500 nom) | TMS320F280049C | ti.com tidueg2c |
| TI TIDM-02013 | 7.4 kW | 400 V | 250–450 V | ~500 kHz | C2000 | ti.com tiduf18a |
| KIT 20 kW CLLLC | 20 kW | 400–800 V | 200–500 V | 100–300 kHz | — | KIT 1000123579 |
| Wiley CTA 4181 (V2G) | 30 kW | 800 V | 250–450 V | — | — | onlinelibrary cta.4181 |
| Navitas SiC 3-φ CLLLC | 20 kW | 800 V | 250–450 V | ~150 kHz | — | Navitas 2024 |

**Excluded** (despite earlier mentions in `LLC_REWORK_PLAN.md` and
`FUTURE_TOPOLOGIES.md`): Wolfspeed CRD-22DD12N is **CLLC**, not
CLLLC. Use it only as bidirectional-resonant baseline, not CLLLC
validation. Same for Infineon UG-2020-31 (REF-DAB11KIZSICSYS) —
that's CLLC.

---

## 8. Future work (out of scope for this plan)

- **Three-phase CLLLC** — 12 switches per side, 3-φ transformer with
  integrated leakage. Reference: Navitas SiC 3-φ CLLLC 20 kW.
  Distinct enough to be a separate class (`ThreePhaseCllllc`).
- **Modular / interleaved CLLLC** — N parallel cells with phase-shifted
  carriers, > 50 kW. Per-cell magnetics are basic CLLLC; cell-balancing
  is the new analysis. Defer.
- **Direction-reversal transient model** — v1 picks one direction per
  OP and does not simulate the reversal. Adding the transient would
  require switched-state-machine handling near the gain = 1 crossing.
- **Closed-loop synchronous-rectifier gating** — v1 emits SR gate
  PULSE statements at predicted current zero-crossings; v2 should
  make them respond to the actual SPICE-simulated zero-crossing.
- **Transformer parasitic capacitance** — Cps (inter-winding) becomes
  relevant > 500 kHz; v1 omits.
- **Coupled with PFC stage** — full OBC simulation (PFC + CLLLC
  cascade); deferred to a separate "OBC" composite-model.

After Track A lands, update `project_converter_model_quality_tiers.md`
to add CLLLC at DAB-quality tier, and update `FUTURE_TOPOLOGIES.md`
to remove CLLLC from the "in progress" list.

---

## 9. References

- **TI TIDM-02002** — *6.6 kW Bidirectional CLLLC Resonant Converter
  for OBC*, user guide tidueg2c
- **TI TIDM-02013** — *7.4 kW Bidirectional OBC* (CLLLC stage),
  tiduf18a
- **TI Power Tips 135** — *Control scheme of a bidirectional CLLLC
  resonant converter in an ESS* (ssztd96)
- **TI PMP41042** — CLLLC test report (tidt367)
- **MDPI Electronics 12(7):1605** — *Bidirectional CLLLC Resonant
  Converter with Hybrid PFM+PSM Control*
- **MDPI Energies 18(14):3815** — *V2G Fixed-Frequency CLLLC*
- **MDPI Appl Sci 10(22):8144** — *Symmetric Isolated Bidirectional
  Converter*
- **Infineon AN-2024-06** — *Operation and modelling analysis of a
  bidirectional CLLC converter* (applies in part to CLLLC)
- **Infineon UG-2020-31** — REF-DAB11KIZSICSYS 11 kW CLLC reference
  (note: CLLC, not CLLLC; use as baseline only)
- **IEEE 6972179** — *Steady-state analysis of bidirectional CLLLC
  in time domain*
- **IEEE 9352422** — *CLLLC ASSA — Advanced State-Space Analysis*
- **IEEE 10362332** — *Asymmetric CLLLC OBC* (Sharma, Missouri S&T)
- **Navitas SiC 3-φ CLLLC** — 20 kW SiC-based three-phase CLLLC for
  high-power OBC (2024)
- **KIT 20 kW CLLLC** — Karlsruhe Institute analysis publication
  1000123579
- **Wiley CTA 4181** — *30 kW V2G CLLLC* (2025)
- **IET PEL 2025 Ansari** — *FEA-based inserted-shunt planar
  transformer for CLLLC*
- **Typhoon-HIL** — CLLLC reference circuit documentation
- **EDN Power Tips 135** — *Control scheme of bidirectional CLLLC
  in ESS*
- **VT CPES** — Bidirectional resonant V2G thesis (Virginia Tech)
- **Huang SLUP263** — TI *Designing an LLC Resonant Half-Bridge
  Power Converter* (ZVS condition extended for CLLLC)
