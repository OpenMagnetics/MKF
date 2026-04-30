# Vienna Rectifier Plan

Companion to `CONVERTER_MODELS_GOLDEN_GUIDE.md`. Adds a Vienna
rectifier (3-phase, 3-level, unidirectional boost-type PFC) to MKF.
Vienna is the dominant front-end topology for **5–50 kW EV DC fast
chargers, telecom rectifiers, aviation power supplies, and high-end
welding / UPS rectifier stages**. Originally proposed by Kolar & Zach
in 1994; modern industrial implementations use SiC and reach
98 %+ efficiency at 10–30 kW per stage.

Distinctive headline: this is the **first MKF model with three
distinct magnetic components** (one boost inductor per phase, all
identical, all carrying line-frequency AC at 120° phase shift). It
is also the **first 3-level topology**, the **first 3-phase topology**,
and the **first line-cycle-envelope topology** in the MKF catalogue.
Three new patterns at once.

This is a planning document — no code changed. The existing plans
(`LLC_REWORK_PLAN.md`, `CLLC_REWRITE_PLAN.md`, `SEPIC_PLAN.md`,
`ZETA_PLAN.md`, `CUK_PLAN.md`, `WEINBERG_PLAN.md`,
`FOUR_SWITCH_BUCK_BOOST_PLAN.md`, `CLLLC_PLAN.md`,
`FUTURE_TOPOLOGIES.md`) are not touched by this plan.

---

## 0. Status snapshot

| Aspect | Today | After this plan |
|---|---|---|
| `Vienna.cpp` / `.h` | does not exist | new files, DAB-quality from start |
| Schema | nothing | full `vienna.json` with `viennaVariant`, `switchType`, `synchronousRectifier`, `phaseCount`, `lineFrequency`, `samplingStrategy` |
| Magnetic-adviser path | no Vienna support | sizes **three identical** boost inductors with **zero DC bias**, AC excitation at line + switching frequencies |
| Inductor return type | `Inputs` (single component) | **`std::vector<Inputs>` of length 3** — the first MKF model returning multiple magnetic components |
| Operating-point semantics | single steady-state OP | **line-cycle envelope** (1× peak-of-line + 6× 30°-sector samples, default) |
| Three-level switching | n/a | switch V-stress = Vdc/2 (the 3-level advantage); SPICE netlist emits split DC bus with neutral midpoint |
| Output structure | single Vout | split DC bus: Vdc+ / Vdc− / N (neutral midpoint) |
| Tests | none | TIDA-010257 + STDES-VIENNARECT + Microchip MSCSICPFC-REF5 references, line-cycle envelope, neutral-point balancing, three-correlated-inductor sizing, NRMSE, SPICE |

**Vienna I** is the canonical industrial variant (single bidirectional
switch per leg, 3 switches total, 6 fast rectifier diodes). Vienna II
(two switches per leg) is a `viennaVariant` flag. Bidirectional
Vienna does **not exist** — the diode bridge is intrinsically
unidirectional; bidirectional 3-phase 3-level needs ANPC, which is
a separate topology (out of scope, see § 8).

---

## 1. Why add Vienna

1. **Dominant 3-phase PFC in the 5–50 kW range.** Every modern EV DC
   fast charger front-end (per stage), telecom rectifier (5–30 kW),
   aviation 270 V_DC supply (3–25 kW), high-power welding supply, and
   UPS rectifier uses Vienna I (TI TIDA-010257, ST STDES-VIENNARECT,
   Microchip MSCSICPFC-REF5, Infineon REF_11KW_PFC_SIC_QD). Below
   ~3 kW single-phase totem-pole wins; above ~50 kW interleaved
   Vienna or ANPC win. Vienna owns the meaty middle.
2. **First 3-phase topology in MKF.** All existing MKF models are
   single-phase (DC-DC, single-phase PFC, or transformer-coupled
   single-phase). Adding Vienna establishes the patterns that future
   3-phase models (3-phase DAB, 3-phase Vienna ANPC, three-phase
   CLLLC, motor drives) will reuse: abc/αβ/dq coordinate handling,
   per-phase operating points, balanced/unbalanced excitation.
3. **First multi-component magnetic adviser request.** Vienna requires
   sizing **three identical correlated inductors** at once. The
   existing MKF API returns a single `Inputs` per converter — Vienna
   needs `std::vector<Inputs>`. This pattern unlocks future topologies
   with multi-magnetic builds (3-phase DAB transformer with three
   single-phase or one 3-leg core, interleaved buck with N inductors,
   etc.).
4. **First zero-DC-bias inductor problem.** Every existing MKF
   converter sizes inductors with a DC bias (buck Iout, boost Iout,
   SEPIC L1=Iin, etc.). Vienna's boost inductors carry **pure AC
   line current** with average ≈ 0. Core selection shifts from
   gapped ferrite (good for high DC bias) to **powder cores**
   (Kool Mu, sendust, High-Flux — soft saturation under transient,
   AC-loss-optimised). The magnetic adviser must know the difference
   to recommend appropriate core families.
5. **First line-cycle-envelope simulation.** Vienna's operating point
   sweeps over each 50/60 Hz line cycle; the inductor sees the
   sector-dependent duty cycle d(t) = |sin(ωt − φ_k)|. v1 should
   sample at peak-of-line + 30°/60° sectors (cheap, 6 samples per
   line cycle); the full 20 ms × kHz-switching simulation is too
   expensive for SPICE per OP.

---

## 2. Variant taxonomy (industry survey)

Sources: Kolar & Zach 1994 (PCIM Nürnberg, the original Vienna
paper), Kolar APEC 2018 plenary "Vienna Rectifier and Beyond" (ETH-
PES), Hartmann ETH dissertation 19755 (2011), Friedli & Kolar IEEE
TPEL "The Essence of Three-Phase PFC Rectifier Systems Part II", TI
TIDUCJ0B (TIDM-1000 user guide), TI TIDA-010257 (10 kW Vienna PFC),
TI TIDA-010210 (11 kW 3L-ANPC bidirectional comparator), ST UM2975
(STDES-VIENNARECT 15 kW), ST UM3011 (STDES-30KWVRECT 30 kW),
Infineon REF_11KW_PFC_SIC_QD app note, Avnet Vienna-for-EV-chargers
article, Microchip Vienna reference + MSCSICPFC-REF5 30 kW (DS50002952B),
Vincotech 3-level PFC benchmark (2016), MDPI Electronics 14(5):936
(2025 NP-balancing review), Wisconsin Zhu MS thesis (GaN Vienna),
IEEE 10609055 (10 kW Vienna inductor design alternatives).

### Belong as **flags / enums on the Vienna class** (NOT new models)

| Variant | Mechanism | Why same class |
|---|---|---|
| **Vienna II** (two switches per leg) | `viennaVariant: enum {viennaI, viennaII}` (default `viennaI`) | Reduced conduction loss at > 30 kW; same magnetics, only switch-loss model differs. References: Kolar APEC 2018, Vincotech 2016. |
| **Switch type** | `switchType: enum {singleMosfetIn4DiodeBridge, backToBackMosfet, tType}` (default `tType`) | All three implement the bidirectional switch per leg differently; magnetics unchanged. Modern SiC designs use T-type (TIDA-010257, ST 15 kW). |
| **Synchronous-rectifier diode replacement** | `synchronousRectifier: bool` (default `false`) | Replaces the 6 fast rectifier diodes with MOSFETs for higher efficiency at high power. Same magnetics. |
| **Phase / interleave count** | `phaseCount: int` (default 1, must be ≥ 1) | Multi-channel paralleling for > 50 kW. Per-channel magnetics identical; ripple cancels at `phaseCount × Fsw` at the DC bus. |
| **Line frequency** | `lineFrequency: number, default 50` (Hz) | 50 (EU) vs 60 (US) Hz; affects line-cycle envelope sample placement. |

### Belong as **separate model classes**

| Variant | Why distinct | Reference |
|---|---|---|
| **2-level 3-phase boost rectifier** | 6 IGBTs + 6 diodes; switch V-stress = Vdc (not Vdc/2); no neutral midpoint. Different topology, not a Vienna variant. | Friedli & Kolar Part II |
| **3-level ANPC (Active Neutral-Point Clamped)** | Bidirectional 3-phase 3-level, used in V2G fast chargers. TIDA-010210 explicitly compares Vienna vs 3L-ANPC. Different switch arrangement, different control. | TI TIDA-010210 |
| **Vienna for single-phase** | Single-phase Vienna doesn't exist as a topology — single-phase PFC uses totem-pole or basic boost. Already covered by `PowerFactorCorrection.cpp`. | n/a |

### Skip outright

- **Vienna III** — academic curiosity; modular variants exist but no production volume.
- **Bidirectional Vienna** — does NOT exist. The diode bridge is intrinsically unidirectional. If a user asks for "bidirectional Vienna", point them at 3L-ANPC.
- **Hysteresis-control variants / one-cycle control** — research only.

### Direct user-question equivalents

If a future request asks *"is bidirectional Vienna missing?"* — answer:
it's not Vienna at all. It's 3L-ANPC, a separate topology. Vienna's
unidirectionality is intrinsic.

If asked *"is GaN Vienna a separate model?"* — no. Same topology,
just different switch device characteristics. It affects the SPICE
loss model and the achievable Fsw, not the magnetics.

---

## 3. Existing patterns we will reuse

### `PowerFactorCorrection.cpp` for the line-cycle envelope idea

The existing `PowerFactorCorrection.cpp` is the only MKF model that
already deals with line-cycle excitation. Read it first to understand
how MKF currently represents a PFC operating point — likely as a
peak-of-line single-OP simplification, with the line-cycle envelope
deferred. Verify whether `PowerFactorCorrection.cpp` covers totem-pole
PFC (per `FUTURE_TOPOLOGIES.md` Tier-3 note) or only basic boost PFC.

If `PowerFactorCorrection.cpp` already has a line-cycle sampling
helper, reuse it. Otherwise, build the helper as part of this plan
and design it so totem-pole PFC and future PFC variants can call it.

### `Boost.cpp` for the per-phase boost equations

Each phase of Vienna is a boost converter from `V_ph` (instantaneous
phase voltage) to `Vdc/2` (per-half DC bus). When `V_ph > 0`, the
phase boosts to Vdc+; when `V_ph < 0`, to Vdc−. The per-phase
sub-interval solver is **literally** Boost.cpp's piecewise-linear
solver, evaluated at the line-cycle sample point.

But Boost.cpp is **partial-quality** (per the SEPIC plan audit) — no
header docstring, no primary-switch RC snubber, no `R_load = max(Vo/
Io, 1e-3)` divide-by-zero guard. Per the 4SBB plan recommendation,
extract Boost.cpp's solver into a `PwmBridgeSolver.h` static helper
(with DAB-quality polish) so Vienna and Boost both call it. This is
the same refactor 4SBB needs; do it once for both.

### `CommonModeChoke.cpp` for the AC input filter

Vienna typically has an EMI common-mode choke at the AC input, sized
for switching-frequency CM noise. This is **separate** from the boost
inductors and is already covered by `CommonModeChoke.cpp`. The Vienna
class should NOT model the CM choke; it should document that the user
adds it as a separate magnetic per the EMI-filter design (CISPR 11
Class A / IEC 61000-3-2).

### Schema convention (camelCase)

The new file is `MAS/schemas/inputs/topologies/vienna.json`, sibling
to `boost.json` / `powerFactorCorrection.json`.

### Multi-component magnetics — NEW pattern

This is the design decision that needs care. Two options:

**Option A (recommended)**: extend the `process_operating_points`
return type. Vienna returns a list of length 3 — one
`OperatingPoint` per phase inductor. The downstream magnetic adviser
iterates the list and sizes each. This is the cleanest extension of
the existing API.

**Option B**: return a single `OperatingPoint` with a list of
winding excitations, all assumed to share one core. This is the
3-phase coupled-inductor case (one 3-leg core), which is rare in
production and academic-only. Defer to a future flag
`coupledThreeLeg: bool`.

v1: implement Option A. Document Option B as future work in § 8.

---

## 4. Track A — Basic Vienna I (the only track in this plan)

### A.1 Files to create

```
src/converter_models/Vienna.h
src/converter_models/Vienna.cpp
tests/TestTopologyVienna.cpp
MAS/schemas/inputs/topologies/vienna.json
```

Plus add `vienna` to the topology enum in
`MAS/schemas/inputs/designRequirements.json` (verify camelCase).

### A.2 Schema (`vienna.json`)

```jsonc
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://psma.com/mas/inputs/topologies/vienna.json",
  "title": "vienna",
  "description": "The description of a Vienna rectifier (3-phase 3-level PFC) excitation",
  "type": "object",
  "properties": {
    "lineToLineVoltage": {
      "description": "AC line-to-line RMS voltage (e.g. 400 V for EU, 480 V for US)",
      "$ref": "../../utils.json#/$defs/dimensionWithTolerance"
    },
    "lineFrequency": {
      "description": "Line frequency (50 or 60 Hz)",
      "type": "number", "default": 50
    },
    "outputDcVoltage": {
      "description": "Total DC bus voltage Vdc = Vdc+ + |Vdc-| (typical 700–800 V)",
      "type": "number"
    },
    "switchingFrequency": { "type": "number" },
    "currentRippleRatio": {
      "description": "ΔI_L,pp / I_pk at peak-of-line (typical 0.2–0.3)",
      "type": "number", "default": 0.25
    },
    "neutralPointVoltageRippleRatio": {
      "description": "(Vdc+ − |Vdc-|) / Vdc tolerance (typical 0.02)",
      "type": "number", "default": 0.02
    },
    "powerFactor": { "type": "number", "default": 0.99 },
    "thdLimit": { "type": "number", "default": 0.05 },
    "efficiency": { "type": "number", "default": 0.97 },
    "viennaVariant": {
      "description": "Switch arrangement",
      "$ref": "#/$defs/viennaVariant",
      "default": "viennaI"
    },
    "switchType": {
      "description": "How each per-leg bidirectional switch is realised",
      "$ref": "#/$defs/viennaSwitchType",
      "default": "tType"
    },
    "synchronousRectifier": {
      "description": "If true, fast rectifier diodes replaced with MOSFETs",
      "type": "boolean", "default": false
    },
    "phaseCount": {
      "description": "Number of interleaved Vienna channels in parallel (1 = single Vienna)",
      "type": "integer", "default": 1, "minimum": 1
    },
    "samplingStrategy": {
      "description": "How the line-cycle envelope is sampled",
      "$ref": "#/$defs/viennaSamplingStrategy",
      "default": "peakOfLinePlusSectors"
    },
    "operatingPoints": {
      "type": "array",
      "items": { "$ref": "#/$defs/viennaOperatingPoint" },
      "minItems": 1
    }
  },
  "required": [
    "lineToLineVoltage", "outputDcVoltage", "switchingFrequency", "operatingPoints"
  ],
  "$defs": {
    "viennaVariant": {
      "title": "viennaVariant",
      "type": "string",
      "enum": ["viennaI", "viennaII"]
    },
    "viennaSwitchType": {
      "title": "viennaSwitchType",
      "type": "string",
      "enum": ["singleMosfetIn4DiodeBridge", "backToBackMosfet", "tType"]
    },
    "viennaSamplingStrategy": {
      "title": "viennaSamplingStrategy",
      "type": "string",
      "enum": [
        "peakOfLineOnly",
        "peakOfLinePlusSectors",
        "fullLineCycle"
      ]
    },
    "viennaOperatingPoint": {
      "title": "viennaOperatingPoint",
      "allOf": [
        { "$ref": "../../utils.json#/$defs/baseOperatingPoint" },
        {
          "properties": {
            "outputPower": { "type": "number" },
            "lineCurrentImbalanceRatio": {
              "description": "max(I_a, I_b, I_c) / min(I_a, I_b, I_c) for unbalanced grid",
              "type": "number", "default": 1.0
            }
          }
        }
      ]
    }
  }
}
```

### A.3 Class structure (`Vienna.h`)

```cpp
class Vienna : public MAS::Vienna, public Topology {
public:
    enum class Phase { A, B, C };
    enum class LineSector { S1, S2, S3, S4, S5, S6 };  // 6 × 60° sectors per cycle

private:
    // 2.1 Simulation tuning
    int numPeriodsToExtract     = 3;     // switching periods at each sample point
    int numSteadyStatePeriods   = 2;
    int lineCycleSampleCount    = 6;     // sectors sampled when peakOfLinePlusSectors

    // 2.2 Computed design values
    double computedBoostInductance       = 0;   // single value (all 3 phases identical)
    double computedDcBusCapacitancePlus  = 0;   // C+
    double computedDcBusCapacitanceMinus = 0;   // C− (= C+ in symmetric design)
    double computedModulationIndex       = 0;   // M
    double computedSwitchVoltageStress   = 0;   // = Vdc/2
    double computedDiodeVoltageStress    = 0;   // = Vdc/2
    double computedLinePeakCurrent       = 0;   // I_pk per phase
    double computedSwitchingDeadTime     = 200e-9;

    // 2.3 Schema-less device parameters
    double mosfetCoss        = 200e-12;
    double diodeReverseTrr   = 0.0;
    double bootstrapDeadTime = 100e-9;

    // 2.4 Per-OP diagnostics — populated for each line-cycle sample
    mutable std::vector<double> lastInductorPeakCurrentPerSector;     // size 6
    mutable std::vector<double> lastInductorRmsCurrentPerSector;      // size 6
    mutable std::vector<double> lastInductorRippleAmplPerSector;      // size 6
    mutable std::vector<double> lastDutyCyclePerSector;               // size 6
    mutable double lastInductorPeakCurrentWorstCase = 0.0;            // max over sectors
    mutable double lastInductorRmsCurrentLineCycle  = 0.0;            // line-cycle RMS
    mutable double lastSwitchPeakCurrent            = 0.0;
    mutable double lastSwitchRmsCurrent             = 0.0;
    mutable double lastDiodePeakCurrent             = 0.0;
    mutable double lastDiodeAvgCurrent              = 0.0;
    mutable double lastNeutralPointVoltageRipple    = 0.0;
    mutable double lastNeutralPointDcOffset         = 0.0;            // |Vdc+| - |Vdc-|
    mutable double lastTotalHarmonicDistortion      = 0.0;            // input-current THD
    mutable double lastPowerFactor                  = 0.0;
    mutable double lastConductionMode               = 0.0;            // 0=CCM, 2=DCM (rare)

    // 2.5 Extra-component waveforms — ONE PER PHASE × ONE PER SECTOR.
    // Cleared in process_operating_points.
    mutable std::vector<std::vector<Waveform>> extraInductorCurrentWaveforms;  // [phase][sector]
    mutable std::vector<std::vector<Waveform>> extraInductorVoltageWaveforms;
    mutable std::vector<std::vector<Waveform>> extraSwitchingNodeVoltageWaveforms;
    mutable std::vector<Waveform> extraDcBusPlusVoltageWaveforms;     // [sector]
    mutable std::vector<Waveform> extraDcBusMinusVoltageWaveforms;
    mutable std::vector<Waveform> extraNeutralCurrentWaveforms;       // i_N at midpoint
    mutable std::vector<Waveform> extraTimeVectors;

public:
    bool _assertErrors = false;

    Vienna(const json& j);
    Vienna() {}

    // Accessors per Golden Guide § 2.7–2.10.

    // Static analytical helpers
    static double compute_modulation_index(double V_phase_peak, double Vdc);
    static double compute_line_peak_current(double P, double V_phase_rms, double eff, double pf);
    static double compute_inductor_for_ripple(double V_phase_peak, double duty_at_peak,
                                              double Fsw, double rippleRatio, double I_pk);
    static double compute_dc_bus_cap(double P, double Vdc, double npRippleRatio,
                                     double lineFreq);
    static double compute_switch_rms(double I_pk, double M);   // Kolar 1994 § III
    static double compute_diode_avg(double I_pk);              // = I_pk / π
    static double compute_duty_at_phase_angle(double M, double angle, Phase phase);
    static double compute_phase_voltage_at_angle(double V_phase_peak, double angle, Phase phase);
    static LineSector classify_sector(double angle);

    // Returns a list of OperatingPoints, one per phase inductor.
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    // Per-sector worker — solves one switching cycle at a fixed line-cycle phase angle.
    OperatingPoint process_operating_point_for_sector(
        double lineAngleRad,
        Phase phase,
        const ViennaOperatingPoint& op);

    // Extra components factory — returns 3 inductor Inputs + 2 cap Inputs.
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // SPICE
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0,
        double lineAngleRad = M_PI / 2.0);   // default = peak-of-line
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);
};

class AdvancedVienna : public Vienna {
    // bypasses process_design_requirements; user provides L directly
};
```

### A.4 Header docstring

Mandatory blocks (per Golden Guide § 3) before the class:

1. **TOPOLOGY OVERVIEW** with ASCII art (Vienna I, T-type switch):

   ```
                         L_a                D1_a
       Va ──/\/\/\──/\/\/\── ●────|>|──┬── Vdc+ ──┬──
            CM                │         │          │
                            S_a (T)     │          C+
                              │         │          │
                              │         │          ●── Neutral N
                              │         │          │
                              │         │          C−
                              │         │          │
                            S_a (T)     │          │
                              │         │          │
       Vb ── (idem L_b, D1_b, S_b) ──── ┤          │
                                        │          │
       Vc ── (idem L_c, D1_c, S_c) ──── ┴── Vdc− ──┴──
   ```

   Three-level structure: each phase node switches between {Vdc+, N,
   Vdc−} relative to the neutral. Switch blocking voltage = Vdc/2
   (the 3-level advantage).

   Document the three switch realisations (`viennaSwitchType`):
   - **`tType`** (modern SiC default): two anti-series MOSFETs
     between phase node and neutral N.
   - **`backToBackMosfet`**: same idea, slightly different physical
     packaging.
   - **`singleMosfetIn4DiodeBridge`** (Kolar 1994 classic): single
     MOSFET inside a 4-diode bridge per phase.

2. **MODULATION CONVENTION** — average-current control with
   sinusoidal current shaping. Per phase:
   `I_ref(t) = I_pk · sin(ωt − φ_phase)`, with the inner current loop
   forcing `I_phase(t) → I_ref(t)`.
   Per-phase duty `d_k(t) = 1 − |sin(ωt − φ_k)| · M` (where M is the
   modulation index `√2·V_phase / (Vdc/2)`).

3. **KEY EQUATIONS** — use the table in § 6 below verbatim with
   `[Kolar 1994]`, `[TI TIDUCJ0B]`, `[ST UM2975]` citations inline.

4. **References** — Kolar & Zach 1994, Kolar APEC 2018, Friedli &
   Kolar TPEL Part II, TI TIDUCJ0B, TI TIDA-010257, ST UM2975, MDPI
   Electronics 14(5):936 (NP balancing review).

5. **Accuracy disclaimer** mirroring `Dab.cpp:296–320`:
   - **Line-cycle envelope is sampled, not simulated continuously.**
     v1 takes 6 samples (every 30°) plus peak-of-line; THD is
     interpolated from sector duties, not extracted from a full
     line-cycle FFT. Accuracy: ±5 % on per-inductor RMS, ±10 % on
     THD, vs a full 20 ms line-cycle simulation.
   - **Three-phase grid is assumed balanced** unless
     `lineCurrentImbalanceRatio ≠ 1.0` is provided. Sequence-component
     analysis (positive/negative/zero) under fault is NOT implemented.
   - **Neutral-point voltage is assumed balanced** by an outer control
     loop; v1 reports `lastNeutralPointDcOffset` as a diagnostic but
     does not model the NP-balancing controller. References: MDPI
     14(5):936.
   - Coss(Vds) variation, body-diode reverse recovery, and gate-driver
     delay omitted.
   - DCM is rare in Vienna (only at < 5 % load); v1 implements CCM
     only and flags `lastConductionMode = 2` if K < K_crit.
   - 3-phase coupled inductor (one 3-leg core for all phases) NOT
     supported in v1; defaults to three independent powder cores.

6. **Disambiguation note**:
   - **Vienna** vs **2-level 3-phase boost rectifier** — 2-level uses
     6 IGBTs + 6 diodes, switch V-stress = Vdc (full bus). Different
     topology, separate model.
   - **Vienna** vs **3L-ANPC** (Active Neutral-Point Clamped) — ANPC
     is bidirectional, used for V2G fast charging. Vienna is
     unidirectional (intrinsic — diode bridge). TI TIDA-010210
     compares them directly.
   - **Vienna I** vs **Vienna II** — single switch per leg vs two
     switches per leg. Same topology, different conduction-loss
     trade-off; covered by `viennaVariant` flag.
   - **Vienna** vs **Single-phase PFC / totem-pole** — three-phase vs
     single-phase. `PowerFactorCorrection.cpp` handles single-phase.

### A.5 Analytical solver (`process_operating_points`)

This is the **headline new pattern** — line-cycle envelope sampling
+ multi-component magnetic output. The algorithm:

1. **Validate the OP**: line-to-line voltage in plausible range,
   `Vdc > √2 · V_LL` (boost requirement satisfied), `M = √2 · V_phase
   / (Vdc/2) ≤ 1` (no over-modulation). Throw with a clear error
   if any violated.

2. **Compute line-cycle invariants**:
   ```
   V_phase_peak = √2 · V_LL / √3
   I_pk         = √2 · P / (3 · V_phase_rms · η · PF)
   M            = V_phase_peak / (Vdc / 2)
   ```

3. **Choose sample angles based on `samplingStrategy`**:
   - `peakOfLineOnly`: 1 sample at ωt = π/2 (peak of phase A).
   - `peakOfLinePlusSectors` (default): 7 samples — peak of A, peak
     of B (ωt = π/2 + 2π/3), peak of C (ωt = π/2 − 2π/3), and
     one each at the 30° points within sector 1 (ωt = 0, π/6,
     π/3, π/2). Captures both the peak-current case (peak-of-phase)
     and the peak-ripple case (which occurs near 30° / 60° within
     each sector).
   - `fullLineCycle`: 36 samples (every 10°). Expensive — only used
     for THD validation tests.

4. **For each sample angle ωt and each phase k ∈ {A, B, C}**:
   - Compute the instantaneous phase voltage `V_phase_k(ωt)`.
   - Compute the instantaneous duty `d_k(ωt) = 1 − |sin(ωt − φ_k)| · M`.
     (When `V_phase_k > 0`: phase boosts to Vdc+; when negative: to
     Vdc−. The Vienna switch routes the current via the diode bridge.)
   - Solve **one switching cycle** of the per-phase boost via the
     Boost.cpp piecewise-linear helper (CCM):
     - Sub-interval 1: switch ON, V_L = V_phase_k.
     - Sub-interval 2: switch OFF, V_L = V_phase_k − Vdc/2 (sign
       depends on V_phase_k sign).
   - Extract per-phase: peak current, ripple amplitude, RMS within
     this switching cycle.

5. **Aggregate over the line cycle**:
   - `lastInductorPeakCurrentPerSector[s] = max over phases k`.
   - `lastInductorPeakCurrentWorstCase = max over s`.
   - `lastInductorRmsCurrentLineCycle = sqrt(integral of i_L² over
     line cycle / line period)`. Use sector samples weighted by
     sector arc length.
   - Switch RMS via Kolar 1994 closed form
     (`I_sw_rms = I_pk · sqrt((1 − 8·√2·M / (3π)) / 2)`).
   - Diode average via `I_D_avg = I_pk / π`.

6. **Build three `OperatingPoint` objects** (one per phase). Each
   carries:
   - The phase-A / phase-B / phase-C inductor's voltage and current
     waveforms across all sample sectors.
   - The same DC bus voltage and current.
   - The phase-specific line-frequency excitation pattern (used by
     the magnetic adviser to compute AC core loss at the line
     fundamental).

7. **`process_design_requirements`** computes:
   - `computedBoostInductance = max(L_for_ripple_constraint,
                                    L_for_thd_constraint)`.
   - `computedDcBusCapacitancePlus = computedDcBusCapacitanceMinus
     = compute_dc_bus_cap(P, Vdc, npRippleRatio, lineFreq)`.

8. **Populate diagnostics** at the end:

   ```
   lastInductorPeakCurrentPerSector  = per-sector max over phases
   lastInductorRmsCurrentLineCycle   = sector-weighted RMS
   lastInductorRippleAmplPerSector   = per-sector ΔI_L
   lastDutyCyclePerSector            = per-sector d_k
   lastSwitchPeakCurrent             = same as inductor peak (series)
   lastSwitchRmsCurrent              = Kolar closed form
   lastDiodePeakCurrent              = same as inductor peak
   lastDiodeAvgCurrent               = I_pk / π
   lastNeutralPointVoltageRipple     = computed from diode current at 6× line freq
   lastNeutralPointDcOffset          = 0 in balanced case
   lastTotalHarmonicDistortion       = THD via FFT of sampled current
   lastPowerFactor                   = computed from cos(φ)
   lastConductionMode                = 0 (CCM) or 2 (DCM, rare)
   ```

### A.6 SPICE netlist (`generate_ngspice_circuit`)

The SPICE netlist captures **one switching cycle at a fixed line-cycle
angle** (default = peak-of-phase-A, ωt = π/2). To generate a full
line-cycle envelope, the test harness calls `generate_ngspice_circuit`
6× with different `lineAngleRad` values and post-processes.

Mandatory checkpoints from Golden Guide § 5:

- **Three "frozen" DC sources** representing the instantaneous
  three-phase voltages at the chosen line angle:
  ```
  V_a a_in 0 DC <V_phase_peak · sin(lineAngleRad)>
  V_b b_in 0 DC <V_phase_peak · sin(lineAngleRad − 2π/3)>
  V_c c_in 0 DC <V_phase_peak · sin(lineAngleRad − 4π/3)>
  ```
  (No 50/60 Hz sinusoidal source — the line frequency is much slower
  than the switching cycle.)
- **Three boost inductors** L_a, L_b, L_c with sense sources
  `V_La_sense`, `V_Lb_sense`, `V_Lc_sense`.
- **Three bidirectional active switches** per
  `viennaSwitchType` choice. For `tType` default:
  ```
  S_a1 a_node N_node SW1
  S_a2 N_node a_node SW2  (anti-series MOSFETs for bidirection)
  ```
  Same for phases B and C.
- **R_snub 1k + C_snub 1n across each switch (6 snubbers total)** —
  mandatory.
- **Six fast rectifier diodes** D1..D6 forming the bridge to Vdc±.
  `.model DVIENNA D(IS=1e-12 RS=0.05)`.
- **Two output capacitors** C_plus (Vdc+ to N) and C_minus
  (N to Vdc−), with `.ic v(vdc_plus) = <Vdc/2>` and
  `.ic v(vdc_minus) = <-Vdc/2>` — pre-charging avoids slow start-up.
- **Neutral midpoint node** explicitly named `N`. The split-DC-bus
  node is the 3-level signature; do not omit it.
- **Two load resistors** `R_load_plus = max(Vdc/2 / I_load, 1e-3)`
  from Vdc+ to N, and `R_load_minus` from N to Vdc−. The split keeps
  the neutral defined.
- `.options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500`
- `.options METHOD=GEAR TRTOL=7`
- `stepTime = 1 / (Fsw · 200)`.
- Saved signals: `.save v(a_node) v(b_node) v(c_node) v(vdc_plus)
  v(vdc_minus) v(N) i(V_La_sense) i(V_Lb_sense) i(V_Lc_sense)
  i(V_load_plus_sense) i(V_load_minus_sense)`.
- Comment header with topology, V_LL, Vdc, Fsw, lineAngle, M, P,
  `viennaVariant`, `switchType`, `synchronousRectifier`,
  "generated by OpenMagnetics".

**Per-phase gate-drive emission** at the chosen line angle:

For each phase k, compute the steady-state duty `d_k = 1 − |sin
(lineAngleRad − φ_k)| · M`. Emit the gate signal as a PULSE of
period `1 / Fsw` with on-fraction `d_k`. The phase relationships
between the three gate signals are preserved (all three switch at
the same Fsw, just with different duties).

### A.7 `simulate_and_extract_operating_points`

Follow the DAB pattern (`Dab.cpp:1218–1307`) verbatim. Per OP, run
SPICE at each line-cycle sample angle (6× for default
`peakOfLinePlusSectors`), aggregate the per-sector results, and
return one `OperatingPoint` per phase. Analytical fallback tagged
`[analytical]` on failure.

### A.8 `simulate_and_extract_topology_waveforms`

Same hygiene rules as prior plans: grep the netlist for exact node
names. The waveform set per phase is `(time, V_phase_node, I_L_phase,
V_switching_node)` — the same triple Boost.cpp emits, replicated 3×
with phase suffixes.

### A.9 `get_extra_components_inputs`

Returns:

- **Three `Inputs` for L_a, L_b, L_c** (the three boost inductors).
  Each:
  - Same nominal inductance `computedBoostInductance`.
  - DC bias = 0 (line-frequency AC).
  - AC excitation: line-frequency fundamental at `I_pk` amplitude
    plus switching-frequency ripple at `ΔI_L` amplitude.
  - Recommended core class: powder (Kool Mu, sendust, or
    high-flux). Document this in the `Inputs` metadata so the
    magnetic adviser can pre-filter the core catalogue.
- **`CAS::Inputs` for C+** (DC bus positive cap):
  - DC bias = Vdc/2.
  - Ripple = 6 × line-frequency component (NP imbalance) +
    switching-frequency ripple.
- **`CAS::Inputs` for C−** (mirrored, same as C+).
- **`CAS::Inputs` for C_neutral** if the user has an explicit
  neutral cap (typically not separate; same as C+ and C−).

When `phaseCount > 1` (interleaved Vienna), divide the per-channel
inductance by the number of channels in parallel only if the
inductors are physically merged; otherwise return `phaseCount × 3`
inductor inputs (one per channel × phase).

### A.10 Tests (`TestTopologyVienna.cpp`)

Mandatory cases per Golden Guide § 8 + Vienna-specific tests:

1. **`Test_Vienna_TIDA010257_Reference_Design`** — TI 380–400 V_LL ↔
   650 V × 10 kW × 40 kHz reference. Reproduce within 5 %.
2. **`Test_Vienna_TIDM1000_Reference_Design`** — TI 400 V_LL ↔ 800 V
   × 10 kW × 50 kHz.
3. **`Test_Vienna_STDES_VIENNARECT_Reference_Design`** — ST 400 V_LL
   ↔ 800 V × 15 kW × 70 kHz (UM2975).
4. **`Test_Vienna_MSCSICPFC_REF5_Reference_Design`** — Microchip
   30 kW SiC reference (DS50002952B).
5. **`Test_Vienna_Design_Requirements`** — V_LL, Vdc, P round-trip
   through design solver. Verify `Vdc > √2 · V_LL`, `M ≤ 1`, three
   identical inductances returned.
6. **`Test_Vienna_OperatingPoints_Generation_Returns_Three_Inductors`**
   — verify the return value is **a list of length 3**, one per
   phase, with phase-shifted excitations (120° between A/B/C).
   **THIS IS THE HEADLINE TEST — it proves the multi-magnetic API.**
7. **`Test_Vienna_LineCycle_Envelope_Sampling`** — verify
   `peakOfLinePlusSectors` returns exactly 6 sector samples;
   `peakOfLineOnly` returns 1; `fullLineCycle` returns 36. Compare
   the line-cycle RMS computed from the 6-sample weighted sum to the
   36-sample weighted sum; expect agreement within 5 %.
8. **`Test_Vienna_Inductor_Zero_DC_Bias`** — assert the DC component
   of each phase inductor's current waveform (averaged over the
   line cycle) is < 1 % of `I_pk`. **This is the audit-evidence that
   the model correctly handles the AC-only excitation pattern that
   drives core selection.**
9. **`Test_Vienna_Inductor_Worst_Case_Ripple_Sector`** — verify
   `lastInductorPeakCurrentPerSector` peaks near the 30°/60° points
   within each sector (not at peak-of-phase, where `I_pk` is largest
   but ripple is smallest).
10. **`Test_Vienna_Switch_Voltage_Stress_Equals_Vdc_Half`** — assert
    `computedSwitchVoltageStress ≈ Vdc/2` (the 3-level advantage).
    Negate the test for a hypothetical 2-level rectifier, where the
    stress would be Vdc.
11. **`Test_Vienna_Switch_RMS_Kolar_Formula`** — verify
    `lastSwitchRmsCurrent` matches Kolar 1994 § III closed-form
    expression within 2 %.
12. **`Test_Vienna_Diode_Avg_Current`** — verify
    `lastDiodeAvgCurrent ≈ I_pk / π`.
13. **`Test_Vienna_NeutralPoint_Balancing`** — assert
    `lastNeutralPointDcOffset` is < 1 % of Vdc/2 in the balanced
    case. With `lineCurrentImbalanceRatio = 1.1`, verify the offset
    grows predictably (per MDPI 14(5):936).
14. **`Test_Vienna_THD_Sanity`** — assert `lastTotalHarmonicDistortion
    < thdLimit` (default 5 %) per IEC 61000-3-2.
15. **`Test_Vienna_PowerFactor_Sanity`** — assert
    `lastPowerFactor ≥ powerFactor` field (default 0.99).
16. **`Test_Vienna_Variant_Vienna_I_vs_II`** — same OP run with both
    `viennaVariant` settings; verify magnetics unchanged but
    switch-RMS differs (Vienna II has lower per-switch conduction
    loss).
17. **`Test_Vienna_SwitchType_Variants`** — three `switchType` choices
    produce equivalent steady-state OPs but different SPICE
    netlists.
18. **`Test_Vienna_SynchronousRectifier_Flag`** — assert efficiency
    improves by the predicted Vd loss; SPICE netlist contains
    rectifier MOSFETs in place of diodes.
19. **`Test_Vienna_PhaseCount_Interleaved`** — `phaseCount = 2`:
    per-channel current is half, ripple cancels at `2 × Fsw` at the
    DC bus. Returns 6 inductor `Inputs` (3 per channel).
20. **`Test_Vienna_SPICE_Netlist_Three_Phase_Sources`** — netlist
    contains three frozen-DC sources at 120° phase shift, six fast
    diodes, three bidir switches per `switchType`, two output caps
    with `.ic`, neutral midpoint node `N`, snubbers on all switches.
21. **`Test_Vienna_PtP_AnalyticalVsNgspice_PeakOfLine`** — single-
    cycle simulation at peak-of-phase-A; per-inductor current NRMSE
    ≤ 0.15 vs analytical, RMS within 15 %, period-averaged Vdc
    within 5 %.
22. **`Test_Vienna_PtP_AnalyticalVsNgspice_30Degrees`** — same at
    30° within sector 1 (where ripple is largest). NRMSE ≤ 0.15.
23. **`Test_Vienna_LineCycle_FullSweep_THD_Validation`** — long
    test: simulate the 6 sector samples in SPICE, post-process the
    aggregated current to extract THD via FFT, compare to
    `lastTotalHarmonicDistortion`. Expect agreement within 1 %.
24. **`Test_Vienna_Waveform_Plotting`** — CSV dump under
    `tests/output/vienna/<sector>/` per sector.
25. **`Test_AdvancedVienna_Process`** — round-trip
    `AdvancedVienna::process()` returns sane Inputs.

ZVS-boundary tests are **N/A** — Vienna is hard-switched. Document.

Tests (6), (7), (8), (10), and (13) are the headline Vienna-specific
tests that justify a separate model:
- (6) proves the multi-magnetic API works.
- (7) proves the line-cycle envelope sampling.
- (8) proves the zero-DC-bias inductor handling.
- (10) proves the 3-level voltage-stress advantage.
- (13) proves the neutral-point balancing diagnostic.

---

## 5. Class checklist (Golden-Guide § 14 acceptance)

- [ ] Builds clean (`ninja -j2 MKF_tests` returns 0)
- [ ] All `TestTopologyVienna` cases pass (≥ 25 cases per § A.10)
- [ ] `Test_Vienna_PtP_AnalyticalVsNgspice_PeakOfLine` and `_30Degrees`
      NRMSE ≤ 0.15
- [ ] `Test_Vienna_OperatingPoints_Generation_Returns_Three_Inductors`
      confirms multi-magnetic API
- [ ] `Test_Vienna_Inductor_Zero_DC_Bias` confirms AC-only excitation
- [ ] `Test_Vienna_LineCycle_Envelope_Sampling` confirms the
      6-sample sampling matches the 36-sample reference within 5 %
- [ ] Schema fields documented in
      `MAS/schemas/inputs/topologies/vienna.json`
- [ ] Header docstring includes ASCII art (with all three phases,
      all three switches per chosen `switchType`, six diodes, split
      DC bus with neutral midpoint), equations with citations,
      references, accuracy disclaimer, disambiguation note (Vienna
      vs 2L 3φ boost vs 3L-ANPC vs single-phase PFC vs Vienna II)
- [ ] `generate_ngspice_circuit` netlist parses and contains: three
      frozen-DC sources, six fast diodes, three bidir switches per
      `switchType`, six snubbers, two output caps with `.ic`, neutral
      midpoint node `N`, GEAR, ITL=500/500, divide-by-zero-guarded
      load resistors
- [ ] `simulate_and_extract_operating_points` runs SPICE at each
      sector sample and aggregates; falls back analytically with
      `[analytical]` tag on failure
- [ ] `simulate_and_extract_topology_waveforms` references real
      netlist nodes (grep-verified)
- [ ] Per-OP diagnostic accessors populated and at least one test
      reads them (especially the per-sector arrays)
- [ ] **Multi-output: `std::vector<Inputs>` of length 3 returned by
      `get_extra_components_inputs`** — the new pattern
- [ ] `viennaVariant`, `switchType`, `synchronousRectifier`,
      `phaseCount`, `samplingStrategy` all tested
- [ ] Neutral-point balancing diagnostic populated and tested
- [ ] Powder-core recommendation propagated through the magnetic-
      adviser hint
- [ ] No divide-by-zero in design-time or SPICE code

ZVS-boundary, rectifier-type CT/CD/FB, and duty-cycle-loss checkboxes
from the Golden Guide are **not applicable** to Vienna — document
with "N/A — hard-switched 3-phase 3-level PFC" notes.

---

## 6. Equations block (citation table)

Use this verbatim in the header docstring. CCM, balanced 3-phase,
unity power factor unless noted.

### Line-cycle invariants

| Quantity | Formula | Source |
|---|---|---|
| Phase peak voltage | V_phase_peak = √2 · V_LL / √3 | basic 3-phase |
| Modulation index | M = V_phase_peak / (Vdc/2) (≤ 1) | Kolar 1994 |
| Boost requirement | Vdc ≥ √2 · V_LL · 1.05 | Kolar 1994 |
| Per-phase line peak current | I_pk = √2 · P / (3 · V_phase_rms · η · PF) | basic |
| Per-phase line RMS current | I_rms = I_pk / √2 | basic |

### Per-switching-cycle (sampled at line angle ωt)

| Quantity | Formula | Source |
|---|---|---|
| Per-phase voltage at angle | V_k(ωt) = V_phase_peak · sin(ωt − φ_k) | basic |
| Per-phase duty | d_k(ωt) = 1 − \|sin(ωt − φ_k)\| · M | TI TIDUCJ0B |
| Inductor ripple | ΔI_L,pp = V_phase_peak · sin(ωt − φ_k) · (1 − d_k) · d_k · T_sw / L | basic boost |
| Inductor peak current | I_L,pk(ωt) = I_pk · sin(ωt − φ_k) + ΔI_L,pp(ωt) / 2 | derived |
| **Inductor DC bias** | **I_L,dc = 0** | line-frequency AC only |
| Switch peak voltage | V_sw = Vdc / 2 | 3-level advantage [Kolar 1994] |
| Diode peak voltage | V_D = Vdc / 2 | 3-level advantage |

### Component-level RMS / averages

| Quantity | Formula | Source |
|---|---|---|
| Switch RMS | I_sw,rms = I_pk · √((1 − 8√2·M / (3π)) / 2) | Kolar 1994 § III |
| Switch peak current | I_sw,pk = I_pk + ΔI_L_max / 2 | derived |
| Diode average current | I_D,avg = I_pk / π | classical bridge |
| Diode peak current | I_D,pk = I_sw,pk | diode + switch series |
| Inductor RMS over line | I_L,rms = √( I_pk² / 2 + ΔI_rms² ) | line + switching |

### Design

| Quantity | Formula | Source |
|---|---|---|
| L for ripple ratio r at peak-of-line | L ≥ V_phase_peak · (1 − M)·M · T_sw / (r · I_pk) | basic boost |
| DC bus cap C+ = C− | C ≥ P / (2π · f_line · 6 · ΔV_NP · Vdc) | 6× line-freq NP ripple |
| Neutral-point ripple frequency | 6 × f_line | 3-phase rectifier 3rd-harmonic |

---

## 7. Reference designs to include in tests

| ID | V_LL | Vdc | P | Fsw | Switch | Controller | Source |
|---|---|---|---|---|---|---|---|
| TI TIDA-010257 | 380–400 V | 650 V | 10 kW | 40 kHz | T-type SiC | C2000 F2800137 | ti.com |
| TI TIDM-1000 | 400 V | 800 V | 10 kW | 50 kHz | classic Si | C2000 F28377D | TIDUCJ0B |
| ST STDES-VIENNARECT | 400 V | 800 V | 15 kW | 70 kHz | T-type SiC | STM32 | UM2975 |
| ST STDES-30KWVRECT | 400 V | 800 V | 30 kW | 70 kHz | T-type SiC | STM32 | UM3011 |
| Microchip MSCSICPFC-REF5 | 400 V | 800 V | 30 kW | 70 kHz | SiC | dsPIC | DS50002952B |
| Infineon REF_11KW_PFC_SIC_QD | 400 V | 800 V | 11 kW | 65 kHz | SiC Q-DPAK | XMC | Infineon app note |
| Kolar / Miniböck PCIM '01 | 400 V | 800 V | 10 kW | 400 kHz | classic | analog hysteresis | PCIM 2001 |

---

## 8. Future work (out of scope for this plan)

- **3L-ANPC (Active Neutral-Point Clamped) bidirectional 3-level
  rectifier** — separate model class. TIDA-010210 reference. Used
  for V2G fast charging where Vienna's unidirectionality is a
  blocker. Different switch arrangement, different control. Belongs
  in a future `Anpc.cpp` or `ThreePhaseAnpc.cpp` plan.
- **2-level 3-phase boost rectifier** — separate model. 6 IGBTs +
  6 diodes; switch V-stress = Vdc (full bus). Different topology.
  Could share line-cycle envelope helpers with Vienna.
- **3-phase coupled inductor** (one 3-leg core for all three
  phases) — `coupledThreeLeg: bool` flag on the existing class. v1
  uses three independent powder cores per industry default.
- **Full line-cycle SPICE simulation** — v1 samples at 6 sector
  points; full 20 ms × kHz-switching simulation is too expensive.
  Future v2 could add a `samplingStrategy = fullLineCycleSpice`
  path for premium THD / EMI validation.
- **Sequence-component (positive / negative / zero) analysis under
  fault / unbalance** — v1 handles balanced 3-phase only. Adding
  unbalance handling would require dq-frame transformation.
- **Neutral-point balancing controller modelling** — v1 reports the
  diagnostic but assumes the controller does its job. Adding an
  actual NP-balancing simulation (zero-sequence injection) is a v2
  polish item per MDPI 14(5):936.
- **EMI filter design** — common-mode and differential-mode filter
  sizing for IEC 61000-6-3 / CISPR 11 Class A. Belongs in a future
  composite-EMI-filter helper, not in the Vienna class itself.
- **3-phase DAB / 3-phase CLLLC** — once Vienna establishes the
  3-phase patterns (multi-magnetic API, abc/αβ/dq helpers,
  per-phase OP semantics), these higher-power resonant 3-phase
  topologies become tractable.

After Track A lands, update `project_converter_model_quality_tiers.md`
to add Vienna at DAB-quality tier, and update `FUTURE_TOPOLOGIES.md`
to remove Vienna from the "in progress" list.

---

## 9. References

- **Kolar & Zach 1994** — *A novel three-phase three-switch three-
  level unity power factor PWM rectifier*, PCIM '94 Nürnberg,
  pp. 125–138 (the original Vienna paper)
- **Kolar APEC 2018 plenary** — *Vienna Rectifier and Beyond*,
  ETH-PES
- **Hartmann ETH Diss. 19755 (2011)** — *Ultra-Compact and
  Ultra-Efficient Three-Phase PWM Rectifier Systems*
- **Friedli & Kolar TPEL** — *The Essence of Three-Phase PFC
  Rectifier Systems Part II*
- **TI TIDUCJ0B** — *Vienna Rectifier-Based Three-Phase PFC
  Reference Design* (TIDM-1000 user guide)
- **TI TIDA-010257** — 10 kW Vienna PFC reference design
- **TI TIDA-010210** — 11 kW 3L-ANPC bidirectional comparator
- **ST UM2975** — *15 kW Three-Level Three-Phase Vienna Rectifier
  with Digital PFC Control* (STDES-VIENNARECT)
- **ST UM3011** — STDES-30KWVRECT 30 kW user manual
- **Infineon REF_11KW_PFC_SIC_QD** — *11 kW PFC SiC bidirectional
  three-/single-phase inverter app note*
- **Microchip MSCSICPFC-REF5** — *3-Phase 30 kW Vienna PFC
  Reference Design* (DS50002952B)
- **Microchip Vienna 3-Phase PFC Reference Design**
- **Vincotech 2016** — *Benchmark of High-Efficient 3-Level PFC
  Topologies*
- **Avnet** — *Designing Vienna Rectifiers for EV Chargers*
- **MDPI Electronics 14(5):936 (2025)** — *Review of Neutral-Point
  Voltage Balancing in Three-Level Converters*
- **Wisconsin Zhu MS thesis** — *Vienna Rectifier with GaN Devices*
- **IEEE 10609055** — *Part 2: Design Alternatives for Power
  Inductor of the 10 kW 3-Phase Vienna Rectifier*
- **Wikipedia** — *Vienna rectifier*
