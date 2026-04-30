# Zeta Converter Plan

Companion to `CONVERTER_MODELS_GOLDEN_GUIDE.md` and sibling to
`SEPIC_PLAN.md`. Adds a Zeta converter model to MKF. Zeta is the
**dual of SEPIC** — same four reactive elements, same gain
M = D/(1−D), but the LC filter is moved from the input side to the
output side. The consequences for magnetics design and control
behaviour are non-trivial enough to deserve its own model.

This is a planning document — no code changed. The existing plans
(`LLC_REWORK_PLAN.md`, `CLLC_REWRITE_PLAN.md`, `SEPIC_PLAN.md`,
`FUTURE_TOPOLOGIES.md`) are not touched by this plan.

---

## 0. Status snapshot

| Aspect | Today | After this plan |
|---|---|---|
| `Zeta.cpp` / `Zeta.h` | does not exist | new files, DAB-quality from start |
| Schema | hollow enum entry only (`MasMigration.cpp`: `"Zeta Converter" → "zetaConverter"`) | full `zeta.json` with `coupledInductor`, `synchronousRectifier`, `currentRippleRatio` |
| Magnetic-adviser path | no Zeta support | uncoupled (2 inductors) or coupled (1 transformer-style component) |
| Switch drive | n/a | **HIGH-SIDE** (PFET default; bootstrap NFET option) — opposite of SEPIC |
| Rectifier type | n/a | `DIODE` (default) + `SYNCHRONOUS` flag |
| Operating modes | n/a | CCM + BCM + DCM (multi-sub-mode DCM1/DCM2 per MDPI 21/22/7434) |
| Tests | none | TI PMP9581 reference, RHPZ-free Bode sanity, ripple, NRMSE, SPICE |

The cousin topology **Cuk** is intentionally **out of scope** for this
plan. Once Zeta lands, Cuk is the third sibling in the fourth-order
family — see `FUTURE_TOPOLOGIES.md` for the recommended order.

---

## 1. Why add Zeta

1. **It is structurally distinct from SEPIC, not a flag on it.** The
   reactive-element rearrangement changes which inductor sees the
   DC-link ripple and which sees the load ripple. The gate drive
   moves from low-side to high-side. The control-loop dynamics are
   fundamentally different: **Zeta has no right-half-plane zero in
   CCM**, so a Zeta loop can be compensated like a buck (high
   bandwidth) while SEPIC suffers RHPZ-limited bandwidth (TI SLYT372,
   ADI "Low Output-Ripple Zeta"). A magnetic-adviser pass that
   optimises for output-ripple targets will pick different magnetics
   for Zeta than for SEPIC.
2. **The DC-bias allocation between L1 and L2 is the opposite of
   SEPIC.** In Zeta, L1 (input-side / magnetizing) carries the large
   `Iout·D/(1−D)` DC bias and L2 (output-side / filter) carries the
   smaller `Iout` DC bias. The adviser must swap which inductor gets
   the high-saturation core class — this is not derivable by a
   downstream user from the SEPIC model.
3. **Real production volume.** PV / MPPT charge controllers, LED
   drivers, automotive 12 V → 13.2 V rails (TI PMP9581 family),
   medical / sensitive-load rails where RHPZ-free compensation is a
   hard requirement.
4. **`MasMigration.cpp` already maps the schema name.** The slot is
   reserved, just hollow.

---

## 2. Variant taxonomy (industry survey)

Sources used: TI SLYT372 (the canonical Zeta app note), TI SLVAFJ6
(Zeta equations Pt 4), TI SNVA608 (LM5085 P-FET Zeta), TI PMP9581_REVB
+ REVC reference designs, ADI "Low Output Voltage Ripple Zeta DC/DC
Converter Topology" article, ADI LT8711 datasheet, MPS Zeta tutorial,
Electronics Weekly "Zeta: that other DC/DC topology",
Erickson-Maksimović 3rd ed. Ch. 6, MDPI Sensors 21/22/7434
(SEPIC/Cuk/Zeta DCM modelling), IEEE 8343387 (bridgeless Zeta PFC),
WSEAS 2010 isolated-Zeta paper, IEEE 6506052 (asymmetric interleaved),
Coilcraft "Guide to Coupled Inductors", ICnavigator SEPIC/Zeta note.

### Belong as **flags / enums on the Zeta class** (NOT new models)

| Variant | Mechanism | Why same class |
|---|---|---|
| **Coupled-inductor Zeta** | `coupledInductor: bool` (default `false`) | Same equations; magnetics build switches from "two inductors" to "one transformer-style coupled component". Ripple-cancellation condition is identical to coupled-SEPIC (1:1, k ≈ 0.95–0.99). Reference: TI PMP9581_REVB; Coilcraft. |
| **Synchronous Zeta** | `synchronousRectifier: bool` (default `false`) | Diode → low-side MOSFET; topology unchanged, only loss model differs. Reference: ADI LT8711, TI PMP9581 family. |
| **High-side switch type** | `highSideSwitch: enum {PFET, NFET_BOOTSTRAP}` (default `PFET`) | Affects SPICE netlist (driver, dead-time) and Vin range only. PFET simpler, ≤40 V. Bootstrap NFET extends to ~60 V (LT8711). |
| **Interleaved phases** | `phaseCount: int` (default 1) | Same per-phase magnetics, multi-phase ripple cancel (asymmetric, IEEE 6506052). Defer initial implementation; emit `phaseCount=1` only. |

### Belong as **separate model classes**

| Variant | Why distinct | Reference |
|---|---|---|
| **Isolated Zeta** | L1 replaced with a transformer; gives flyback-like isolation but with non-inverting CCM-output filter. Used in PV / LED-driver PFC where isolation + low output ripple are both required. | WSEAS 2010 isolated-Zeta paper |
| **Bridgeless Zeta PFC** | AC-input, current-shaping loop, line-cycle envelope sweep. Eliminates the input bridge. Used in LED drivers > 25 W (IEC 61000-3-2). | IEEE 8343387 |

### Skip outright

- High-gain Zeta with charge-pump / voltage-multiplier cells (research-only — MDPI Electronics 11/3/483)
- Three-level / quadratic Zeta (research-only — KIPE 2016)

### Direct user-question equivalents

If a future request asks *"is coupled-inductor Zeta missing?"* — answer
mirrors the SEPIC and LLC verdicts: it's a flag on the class, not a
separate class. Sources phrase it as "Zeta **with** coupled inductors",
never as a distinct topology.

---

## 3. Existing patterns we will reuse

### SEPIC as the closest precedent

`SEPIC_PLAN.md` defines an analogous fourth-order, single-output,
single-rectifier, hard-switched PWM converter. **Zeta inherits 80 % of
that structure unchanged** — class skeleton, sub-interval analytical
solver, conduction-mode classification, SPICE wrapper, multi-output
warn-once (N/A — single output), test catalogue. The actual Zeta
differences are:

- High-side switch (SPICE driver model + dead-time semantics)
- Roles of L1 and L2 are swapped in the magnetic adviser
- Output capacitor sees small triangular ripple instead of large
  pulsed ripple — different sizing constraint
- DCM has multiple sub-modes (DCM1/DCM2) per MDPI Sensors 21/22/7434

Build Zeta **after** SEPIC lands so the SEPIC infrastructure is
proven; do not block on it.

### Buck.cpp / Boost.cpp as the structural baseline

Same caveat as in `SEPIC_PLAN.md`: copy the *structure* of
`Buck.cpp`'s piecewise-linear CCM/DCM solver, but bring it up to
DAB-quality from day one (header docstring, primary-switch RC snubber,
`R_load = max(Vo/Io, 1e-3)` divide-by-zero guard, accuracy-disclaimer
block, `Test_Zeta_Waveform_Plotting` CSV dump).

### Coupled inductors map to transformer with leakage

Same as SEPIC: there is no `CoupledInductor` class. A coupled-inductor
Zeta uses a 2-winding transformer-style magnetic with deliberate
leakage Lk; ripple cancellation is identical to coupled-SEPIC (TI
SLYT411 condition applies because L1 and L2 see the same voltage
during both subintervals — this is true for both SEPIC and Zeta).

### Schema convention (camelCase)

The new file is `MAS/schemas/inputs/topologies/zeta.json`, sibling to
`buck.json` / `boost.json` / `sepic.json` (when SEPIC lands).

---

## 4. Track A — Basic Zeta (the only track in this plan)

### A.1 Files to create

```
src/converter_models/Zeta.h
src/converter_models/Zeta.cpp
tests/TestTopologyZeta.cpp
MAS/schemas/inputs/topologies/zeta.json
```

Plus add `zeta` to the topology enum in
`MAS/schemas/inputs/designRequirements.json` (already partially there
per `MasMigration.cpp:"Zeta Converter" → "zetaConverter"` — verify
the canonical lowercase name `zeta` matches the buck/boost/sepic
style).

### A.2 Schema (`zeta.json`)

Mirrors `sepic.json` plus a `highSideSwitch` enum:

```jsonc
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://psma.com/mas/inputs/topologies/zeta.json",
  "title": "zeta",
  "description": "The description of a Zeta converter excitation",
  "type": "object",
  "properties": {
    "inputVoltage": { "$ref": "../../utils.json#/$defs/dimensionWithTolerance" },
    "diodeVoltageDrop": { "type": "number", "default": 0.6 },
    "maximumSwitchCurrent": { "type": "number" },
    "currentRippleRatio": { "type": "number", "default": 0.4 },
    "couplingCapacitorRippleRatio": { "type": "number", "default": 0.05 },
    "outputVoltageRippleRatio": { "type": "number", "default": 0.01 },
    "efficiency": { "type": "number", "default": 0.92 },
    "coupledInductor": {
      "description": "If true, L1 and L2 are wound on the same core (1:1 coupled)",
      "type": "boolean", "default": false
    },
    "synchronousRectifier": {
      "description": "If true, the rectifier diode is a low-side MOSFET",
      "type": "boolean", "default": false
    },
    "highSideSwitch": {
      "description": "The type of high-side switch device",
      "$ref": "#/$defs/zetaHighSideSwitch"
    },
    "operatingPoints": {
      "type": "array",
      "items": { "$ref": "#/$defs/zetaOperatingPoint" },
      "minItems": 1
    }
  },
  "required": ["inputVoltage", "operatingPoints"],
  "$defs": {
    "zetaHighSideSwitch": {
      "title": "zetaHighSideSwitch",
      "type": "string",
      "enum": ["pfet", "nfetBootstrap"],
      "default": "pfet"
    },
    "zetaOperatingPoint": {
      "title": "zetaOperatingPoint",
      "allOf": [{ "$ref": "../../utils.json#/$defs/baseOperatingPoint" }]
    }
  }
}
```

After the edit, run the project's MAS code-gen so `MAS::Zeta` becomes
a generated C++ class. Per `project_mas_hpp_staging.md`, copy the
freshly built `MAS.hpp` into `build/MAS/MAS.hpp` if the build complains
about an undeclared `MAS::Zeta`.

### A.3 Class structure (`Zeta.h`)

Mirror `Sepic.h` (when it exists) with three Zeta-specific deltas
flagged with `// Zeta:` comments. The minimum private state:

```cpp
class Zeta : public MAS::Zeta, public Topology {
private:
    // 2.1 Simulation tuning
    int numPeriodsToExtract    = 5;
    int numSteadyStatePeriods  = 3;

    // 2.2 Computed design values (filled by process_design_requirements)
    double computedInputInductance      = 0;   // Zeta: L1 = magnetizing
    double computedOutputInductance     = 0;   // Zeta: L2 = output filter
    double computedCouplingCapacitance  = 0;   // Cc
    double computedOutputCapacitance    = 0;   // Co
    double computedDiodeVoltageDrop     = 0.6;
    double computedDutyCycleNominal     = 0;
    double computedDcmBoundaryK         = 0;   // K_crit = (1−D)² per Erickson §5/6

    // 2.3 Schema-less device parameters
    double mosfetCoss        = 200e-12;
    double diodeReverseTrr   = 0.0;
    double mutualCouplingK   = 0.99;           // only if coupledInductor=true
    double bootstrapDeadTime = 100e-9;         // Zeta: only if highSideSwitch=NFET_BOOTSTRAP

    // 2.4 Per-OP diagnostics
    mutable double lastDutyCycle              = 0.0;
    mutable double lastConductionMode         = 0.0;  // 0=CCM, 1=BCM, 2=DCM1, 3=DCM2
    mutable double lastInputInductorPeak      = 0.0;  // L1 peak (large DC bias)
    mutable double lastOutputInductorPeak     = 0.0;  // L2 peak (small DC bias)
    mutable double lastSwitchPeakVoltage      = 0.0;  // Vin + Vout
    mutable double lastSwitchPeakCurrent      = 0.0;  // = lastDiodePeakCurrent
    mutable double lastDiodePeakVoltage       = 0.0;
    mutable double lastDiodePeakCurrent       = 0.0;
    mutable double lastCouplingCapRmsCurrent  = 0.0;
    mutable double lastOutputVoltageRipple    = 0.0;  // Zeta: small (L2 in series)
    mutable double lastInputCurrentRipple     = 0.0;  // Zeta: large (Q1 pulsed input)

    // 2.5 Extra-component waveforms (cleared in process_operating_points)
    mutable std::vector<Waveform> extraL1CurrentWaveforms;
    mutable std::vector<Waveform> extraL1VoltageWaveforms;
    mutable std::vector<Waveform> extraL2CurrentWaveforms;
    mutable std::vector<Waveform> extraL2VoltageWaveforms;
    mutable std::vector<Waveform> extraCcCurrentWaveforms;
    mutable std::vector<Waveform> extraCcVoltageWaveforms;
    mutable std::vector<Waveform> extraTimeVectors;

public:
    bool _assertErrors = false;

    Zeta(const json& j);
    Zeta() {}

    // Tuning, computed, last_*, mosfet_coss, coupling_k, bootstrap_dead_time
    // accessors — full set per Golden Guide § 2.7–2.10.

    // Overrides
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const ZetaOperatingPoint& op,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // Static analytical helpers
    static double compute_duty_cycle(double Vin, double Vo, double Vd, double eff);
    static double compute_conversion_ratio(double D);
    static double compute_l1_for_ripple(double Vin, double D, double Fs, double dILRatio, double IL1);
    static double compute_l2_for_ripple(double Vin, double D, double Fs, double dILRatio, double Iout);
    static double compute_dcm_critical_k(double D);   // K_crit = (1−D)²
    static double compute_coupling_cap_rms(double Iout, double D);
    static double compute_switch_peak_voltage(double Vin, double Vo, double Vd, double dVcc);
    static double compute_diode_peak_voltage(double Vin, double Vo, double dVcc);
    static double compute_coupling_cap(double dVccRatio, double Iout, double D, double Fs, double Vin);
    static double compute_output_cap_for_ripple(double dVoRatio, double Vo, double dIL2, double Fs);

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

class AdvancedZeta : public Zeta {
    // bypasses process_design_requirements; user provides L1, L2, Cc directly
};
```

### A.4 Header docstring

Mandatory blocks (per Golden Guide § 3) before `class Zeta`:

1. **TOPOLOGY OVERVIEW** with ASCII art:

   ```
   Vin ──┤ Q1 ├──┬── L1 ──┬── Cc ──┬── L2 ──┬── Vout
          (HS)   │        │        │        │
                 │        │       D1        Co
                 │        │   (cathode      │
                 │        │   to L1/Cc)     │
                GND      GND      GND      GND
   ```

   Plus a coupled-inductor variant where L1 and L2 share a core
   (dot-convention diagram). Mark Q1 as **HIGH-SIDE** explicitly —
   this is the single most-confused detail in Zeta tutorials.

2. **MODULATION CONVENTION** — duty cycle D ∈ [0,1] is the high-side
   switch ON-fraction. With `efficiency`:
   D = (Vo + Vd) / (Vin·η + Vo + Vd) (TI SLVAFJ6 Eq. 6).

3. **KEY EQUATIONS** with citations (use the table in § 6 below
   verbatim with `[TI SLYT372]` / `[TI SLVAFJ6]` style references
   inline).

4. **References** — TI SLYT372 / SLVAFJ6 / SNVA608, ADI Zeta article,
   Erickson-Maksimović Ch. 6, MDPI Sensors 21/22/7434.

5. **Accuracy disclaimer** mirroring `Dab.cpp:296–320`:
   - High-side gate-driver bootstrap behaviour modelled as ideal
     square-wave drive; bootstrap-cap ripple not represented.
   - Cc ESR not modelled in CCM analytical solver.
   - DCM solver uses Erickson § 5/6 quadratic for first-mode DCM1;
     deep DCM2 (where iL1 + iL2 cross zero simultaneously inside the
     OFF interval) is approximated, not exactly solved.
   - Coupled-inductor leakage Lk treated as single equivalent value.
   - Diode reverse recovery and MOSFET Coss(Vds) not in v1.

6. **Disambiguation note**:
   - **Zeta** vs **SEPIC** — SEPIC has the LC filter on the **input**
     side (low-side switch, low input ripple, RHPZ in CCM); Zeta has
     it on the **output** side (high-side switch, low output ripple,
     no RHPZ in CCM). Same gain M = D/(1−D).
   - **Zeta** vs **Cuk** — Cuk is inverting; Zeta is non-inverting.
   - **Zeta** vs **Buck-Boost (non-isolated)** — buck-boost inverts
     output and is 2nd-order; Zeta does not invert and is 4th-order.
   - **Zeta** vs **Flyback** — flyback uses an isolated transformer
     with energy storage; Zeta is non-isolated. The "isolated Zeta"
     variant (out of scope) uses a transformer in place of L1.

### A.5 Analytical solver (`process_operating_point_for_input_voltage`)

Use the **piecewise-linear sub-interval rule** (Golden Guide § 4.1).
Zeta has two sub-intervals per period in CCM:

- **t ∈ [0, D·T]** — Q1 ON, D1 OFF: VL1 = +Vin, VL2 = +Vin (both
  inductors see Vin during ON; Cc voltage = Vin in steady state).
  Diode is reverse-biased.
- **t ∈ [D·T, T]** — Q1 OFF, D1 ON: VL1 = −Vo, VL2 = −Vo. Diode
  conducts iL1 + iL2.

Within each sub-interval iL1(t) and iL2(t) are exactly linear, so
propagate by a single multiply-add per sub-interval (Golden Guide
§ 4.1).

Steady-state initial conditions (Golden Guide § 4.2):

- Volt-second balance on each inductor over a period → 0.
- Amp-second balance on Cc over a period → 0.
- Solve for `I_L1_avg = Iout·D/(1−D)` and `I_L2_avg = Iout`.
  (Zeta L1 is the magnetizing inductor — *opposite role* from SEPIC
   where L1 is the input filter.)
- Set the starting current to `I_avg − ΔI/2` for each inductor.
- Set initial `Vcc(0) = Vin` (matches steady state).

Conduction-mode classification: compute K = 2·L_e·Fs / R_load with
L_e = L1‖L2 = L1·L2/(L1+L2), and K_crit = (1 − D)² (Erickson §5/6).

- K > K_crit → **CCM** (use the linear solver above)
- K = K_crit → **BCM** (mark `lastConductionMode = 1`)
- K < K_crit → **DCM** — Zeta has two sub-modes (MDPI Sensors
  21/22/7434):
  - **DCM1**: only one inductor reaches zero current first.
  - **DCM2**: both inductors reach zero simultaneously within the
    OFF interval.
  v1: detect DCM1 via Erickson quadratic; flag DCM2 as
  `lastConductionMode = 3` and use a perturbed CCM solution as
  fallback (note this in the accuracy disclaimer).

Populate every `last_*` diagnostic at the end (Golden Guide § 4.6):

```
lastSwitchPeakVoltage    = Vin + Vo + ΔVcc/2          # TI SLVAFJ6 Eq. 7
lastSwitchPeakCurrent    = I_L1_avg + I_L2_avg + (ΔIL1+ΔIL2)/2
lastDiodePeakVoltage     = Vin + Vo + ΔVcc/2          # TI SLVAFJ6 Eq. 8
lastDiodePeakCurrent     = lastSwitchPeakCurrent
lastCouplingCapRmsCurrent= Iout · √(D/(1−D))
lastOutputVoltageRipple  = ΔIL2 / (8·Fs·Co)           # tiny — L2 filters
lastInputCurrentRipple   = I_L1_avg + ΔIL1/2          # large — pulsed input
```

Sample grid: `2·N + 1` with N = 256, [0, T] (Golden Guide § 4.8).

### A.6 SPICE netlist (`generate_ngspice_circuit`)

Mandatory checkpoints from Golden Guide § 5:

- **`.model M_PFET PMOS VTO=-2.5 KP=20 LAMBDA=0.01`** — for the
  default high-side PFET. If `highSideSwitch == NFET_BOOTSTRAP`,
  emit an NMOS model with a charge-pump driver subcircuit.
- **`.model SW_RECT SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg`** — only
  if `synchronousRectifier == true`.
- **`.model DZETA D(IS=1e-12 RS=0.05)`** — diode (default rectifier).
- **R_snub 1k + C_snub 1n across the high-side switch** — mandatory.
- Per-component sense sources: `Vin_sense`, `VL1_sense`, `VL2_sense`,
  `VCc_sense`, `VD_sense`, `Vout_sense`, `Vsw_sense`.
- **K-coupling** if `coupledInductor=true`: emit
  `K_L1_L2 L_L1 L_L2 <k>` with `k = mutualCouplingK`. Without it
  ngspice treats L1 and L2 as independent, breaking the
  ripple-cancellation analysis.
- `.ic v(vout) = <Vo>`, `.ic v(vcc) = <Vin>` — Cc charges to Vin in
  steady state.
- `.options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500`
- `.options METHOD=GEAR TRTOL=7`
- **`R_load = max(Vo/Io, 1e-3)`** — divide-by-zero guard.
- `stepTime = period / 200`.
- Saved signals: `.save v(vsw) v(vcc) v(vout) i(Vin_sense)
  i(VL1_sense) i(VL2_sense) i(VCc_sense) i(VD_sense)`.
- Comment header with topology, Vin, Vo, Fs, D, L1/L2/Cc/Co,
  rectifier type, `coupledInductor` flag, `highSideSwitch` choice,
  "generated by OpenMagnetics".

**High-side-drive specifics** (Zeta-only, not in SEPIC):

- For PFET: gate driver references VSW node to source (Vin); use a
  voltage-controlled voltage source `Egate vg vsw VALUE = ...` to
  hold the gate at `vsw − Vgs` during ON.
- For NFET_BOOTSTRAP: emit a bootstrap capacitor `Cboot vboot vsw 100n`
  charged from a 12 V auxiliary rail through `Dboot`. Drive the gate
  via `vboot` reference. Add `bootstrapDeadTime` to the Q1 ON edge to
  let Cboot recharge.

### A.7 `simulate_and_extract_operating_points`

Follow the DAB pattern (`Dab.cpp:1218–1307`) verbatim, with the
`simulate_and_extract_*` invariant: real `NgspiceRunner` + analytical
fallback tagged with `[analytical]` on failure.

### A.8 `simulate_and_extract_topology_waveforms`

Round-trip SPICE waveforms to MAS converter-level waveforms. Same
golden-guide hygiene as SEPIC: grep the netlist for exact node names
before referencing them in `getWaveform("v(...)")` / `getWaveform("i(...)")`.

### A.9 `get_extra_components_inputs`

Returns:

- `Inputs` for **L1** (input-side / magnetizing inductor): nominal
  inductance from `computedInputInductance`, **DC bias =
  `I_L1_avg = Iout·D/(1−D)`** (large), ripple ΔI_L1.
- `Inputs` for **L2** (output-side / filter inductor): nominal
  inductance from `computedOutputInductance`, **DC bias = Iout**
  (small), ripple ΔI_L2.
- `CAS::Inputs` for **Cc** (coupling cap): voltage waveform with
  bias = Vin and AC component, RMS current per § 6 equation.
- `CAS::Inputs` for **Co** (output cap) when `mode = REAL`. Sized for
  `outputVoltageRippleRatio` constraint per § 6.

When `coupledInductor = true`, replace the two L1 / L2 `Inputs` with a
**single** 2-winding transformer-style `Inputs`. **Note**: the Zeta
coupled-inductor configuration has the SAME ripple-cancellation
condition as the SEPIC one (TI SLYT411) — the inductors see identical
voltage waveforms during both subintervals — so the magnetic-adviser
calculation can reuse the SEPIC helper.

### A.10 Tests (`TestTopologyZeta.cpp`)

Mandatory cases per Golden Guide § 8:

1. **`Test_Zeta_PMP9581_Reference_Design`** — TI 12 V → 13.2 V @ 4 A
   reference design (uncoupled, REVC), reproduce within 5 %.
2. **`Test_Zeta_PMP9581_CoupledInductor_Reference_Design`** — TI
   PMP9581_REVB (coupled-inductor variant) — verifies
   ripple-cancellation prediction.
3. **`Test_Zeta_LM5085_PFET_Reference_Design`** — TI SNVA608, 12 V
   automotive → 12 V @ 1–2 A, P-FET high-side switch.
4. **`Test_Zeta_Design_Requirements`** — turns ratio (1:1 for
   coupled), L1, L2, Cc, Co all positive; round-trip a power target
   through the design solver. Verify the L1 / L2 DC-bias allocation
   is **opposite of SEPIC**.
5. **`Test_Zeta_OperatingPoints_Generation`** — multiple Vin
   (typical 9–18 V automotive), antisymmetry of inductor ripple,
   both inductor currents piecewise linear, output current matches.
6. **`Test_Zeta_Operating_Modes_CCM_BCM_DCM`** — sweep load to push
   past K_crit, verify `lastConductionMode` flips 0 → 1 → 2 across
   BCM. Add a deeper-load case where DCM2 triggers (mode 3).
7. **`Test_Zeta_CoupledInductor_Flag`** — same OP run twice with
   `coupledInductor=false` and `=true`; ripple amplitude differs
   per TI SLYT411 prediction.
8. **`Test_Zeta_SynchronousRectifier_Flag`** — verify Vo improves
   slightly (no Vd loss); SPICE netlist contains rectifier MOSFET.
9. **`Test_Zeta_HighSideSwitch_PFET_vs_NFET_Bootstrap`** — both
   switch types produce equivalent steady-state Vo; netlist contains
   PMOS model in PFET case and bootstrap subcircuit in NFET case.
10. **`Test_Zeta_OutputRipple_LowerThan_SEPIC`** — when SEPIC plan
    has landed: same OP, same magnetics, compare
    `lastOutputVoltageRipple` of Zeta vs SEPIC. Expect Zeta < SEPIC.
    (This is the headline magnetic-adviser distinction.)
11. **`Test_Zeta_InputRipple_HigherThan_SEPIC`** — dual to (10);
    expect Zeta > SEPIC.
12. **`Test_Zeta_RHPZ_Free_Sanity`** — confirm the analytical solver
    produces a Bode plot or equivalent transfer-function snapshot
    with no RHPZ. Reference: TI SLYT372 § 4. (Optional if a Bode
    helper isn't already in the codebase; otherwise mark a TODO.)
13. **`Test_Zeta_SPICE_Netlist`** — netlist parses, contains snubber,
    `.options METHOD=GEAR`, K-coupling when coupled, IC pre-charge on
    Cc and Co, PMOS or NFET-bootstrap subcircuit per choice.
14. **`Test_Zeta_PtP_AnalyticalVsNgspice`** — primary-current shape
    NRMSE ≤ 0.15, primary RMS within 15 %, Vo within 5 % over a
    period after warm-up.
15. **`Test_Zeta_Waveform_Plotting`** — CSV dump under
    `tests/output/zeta/` for visual inspection. Not asserted.
16. **`Test_AdvancedZeta_Process`** — round-trip
    `AdvancedZeta::process()` returns sane Inputs.

ZVS-boundary tests are N/A — Zeta is hard-switched. Mark explicitly.

Tests (10) and (11) are the two distinctive Zeta tests that justify
having a separate model from SEPIC. They are the audit-evidence that
the magnetic-adviser pass picks different magnetics for Zeta vs SEPIC.

---

## 5. Class checklist (Golden-Guide § 14 acceptance)

- [ ] Builds clean (`ninja -j2 MKF_tests` returns 0)
- [ ] All `TestTopologyZeta` cases pass (≥ 16 cases per § A.10)
- [ ] `Test_Zeta_PtP_AnalyticalVsNgspice` NRMSE ≤ 0.15
- [ ] CCM ↔ BCM ↔ DCM transition test observes `lastConductionMode`
      change at predicted K_crit
- [ ] `Test_Zeta_OutputRipple_LowerThan_SEPIC` and
      `Test_Zeta_InputRipple_HigherThan_SEPIC` confirm the dual
      relationship with SEPIC
- [ ] Schema fields documented in
      `MAS/schemas/inputs/topologies/zeta.json`
- [ ] Header docstring includes ASCII art (with HIGH-SIDE switch
      explicitly labelled), equations with citations, references,
      accuracy disclaimer, disambiguation note (Zeta vs SEPIC vs Cuk
      vs flyback vs buck-boost)
- [ ] `generate_ngspice_circuit` netlist parses and contains: PMOS or
      NFET-bootstrap subcircuit, R_snub+C_snub across switch, GEAR,
      ITL=500/500, `.ic` on Cc and Co, K-coupling when
      `coupledInductor`, `R_load = max(Vo/Io, 1e-3)`
- [ ] `simulate_and_extract_operating_points` actually runs SPICE and
      falls back analytically with `[analytical]` tag on failure
- [ ] `simulate_and_extract_topology_waveforms` references real
      netlist nodes (grep-verified)
- [ ] Per-OP diagnostic accessors populated and at least one test
      reads them (especially the new
      `lastOutputVoltageRipple` / `lastInputCurrentRipple` pair)
- [ ] Multi-output **N/A** — single-output. Document.
- [ ] Synchronous-rectifier flag tested
- [ ] Coupled-inductor flag tested
- [ ] High-side-switch type both branches tested
- [ ] No divide-by-zero in design-time or SPICE code

ZVS-boundary, rectifier-type CT/CD/FB, and duty-cycle-loss checkboxes
from the Golden Guide are **not applicable** to Zeta — document with
"N/A — hard-switched single-rectifier PWM" notes in the test file.

---

## 6. Equations block (citation table)

Use this verbatim in the header docstring. CCM, lossless unless noted.

| Quantity | Formula | Source |
|---|---|---|
| Conversion ratio | M = Vo/Vin = D/(1−D) | Erickson Ch. 6; TI SLYT372 |
| Duty (with η, Vd) | D = (Vo+Vd) / (Vin·η+Vo+Vd) | TI SLVAFJ6 Eq. 6 |
| Coupling-cap voltage | V_Cc = Vin (steady state) | TI SLYT372; ADI Zeta |
| Coupling-cap voltage ripple | ΔV_Cc = Iout·D / (Cc·Fs) | TI SLYT372 |
| L1 average current | I_L1 = Iout·D/(1−D) | TI SLYT372 |
| L2 average current | I_L2 = Iout | TI SLYT372 |
| L1 ripple | ΔI_L1 = Vin·D / (L1·Fs) | TI SLYT372 |
| L2 ripple | ΔI_L2 = Vin·D / (L2·Fs) | TI SLYT372 |
| L1 sizing for ratio r | L1 = Vin·D / (r·I_L1·Fs) | TI SNVA168E (analogous) |
| L2 sizing for ratio r | L2 = Vin·D / (r·I_L2·Fs) | TI SNVA168E (analogous) |
| Coupling-cap RMS current | I_Cc(rms) = Iout·√(D/(1−D)) | TI SLYT372 |
| Coupling-cap sizing | Cc ≥ Iout·D / (ΔV_Cc·Fs) | TI SLYT372 |
| **Output cap sizing (Zeta)** | Co ≥ ΔI_L2 / (8·Fs·ΔVo) | Erickson Ch. 6 (buck output filter) |
| Switch peak voltage | V_SW(pk) = Vin + Vo + ΔV_Cc/2 | TI SLVAFJ6 Eq. 7 |
| Switch peak current | I_SW(pk) = I_L1 + I_L2 + (ΔI_L1+ΔI_L2)/2 = Iout/(1−D) + ripple | TI SLYT372 |
| Diode peak voltage | V_D(pk) = Vin + Vo + ΔV_Cc/2 | TI SLVAFJ6 Eq. 8 |
| Diode peak current | I_D(pk) = I_SW(pk) | TI SLYT372 |
| DCM K | K = 2·L_e·Fs / R_load with L_e = L1·L2/(L1+L2) | Erickson §5/6 |
| DCM K_crit | K_crit = (1−D)² | Erickson §5/6; MDPI 21/22/7434 |
| DCM2 sub-mode boundary | (multi-condition; see MDPI 21/22/7434 § III.B) | MDPI Sensors 21/22/7434 |
| Coupled L ripple-cancel | Lk1 such that VL1·Lk1/(L1·M_mut) = VL2 (input-ripple zero) — same as SEPIC | TI SLYT411 |

**Distinctive Zeta property (NOT in the table; document in disclaimer)**:
the small-signal control-to-output transfer function has **no
right-half-plane zero in CCM** (TI SLYT372 § 4; Erickson Ch. 12). This
gives Zeta significantly higher achievable closed-loop bandwidth than
SEPIC for the same magnetics.

---

## 7. Reference designs to include in tests

| ID | Vin | Vout | Iout | Fsw | Magnetics | Test name |
|---|---|---|---|---|---|---|
| TI PMP9581_REVC | 5–14 V | 13.2 V | 4 A (5 A > 10 V) | ~600 kHz | uncoupled | `Test_Zeta_PMP9581_Reference_Design` |
| TI PMP9581_REVB | 5–14 V | 13.2 V | 4 A | ~600 kHz | coupled-inductor | `Test_Zeta_PMP9581_CoupledInductor_Reference_Design` |
| TI LM5085 (SNVA608) | 9–18 V automotive | 12 V | 1–2 A | constant-on (≈300 kHz) | uncoupled, P-FET | `Test_Zeta_LM5085_PFET_Reference_Design` |
| ADI LT8711 typical app | 4.5–40 V | 12 V | up to 4 A | 200 kHz–750 kHz | coupled, NFET-bootstrap | (optional, alternate) |

---

## 8. Future work (out of scope for this plan)

- **Cuk model** — sibling fourth-order, inverting output, capacitive
  energy transfer through Cs1+Cs2. `MasMigration.cpp` already maps
  `"Cuk Converter" → "cukConverter"`. After Zeta lands, the third
  sibling is a low-marginal-cost addition — copy-paste of Sepic /
  Zeta structure with sign flip.
- **Isolated Zeta** — replace L1 with a transformer; gain becomes
  M = n·D/(1−D). PV / LED-driver use case. Could be a separate class
  or a flag once the basic Zeta lands.
- **Bridgeless Zeta PFC** — AC input, line-cycle envelope sweep. LED
  drivers > 25 W (IEEE 8343387). Defer.
- **Interleaved Zeta** — `phaseCount: int` flag on the basic class
  (asymmetric ripple cancel, IEEE 6506052). Defer.
- **DCM2 exact solver** — Zeta DCM has multiple sub-modes (MDPI
  Sensors 21/22/7434). v1 uses approximation; an exact solver is a
  v2 polish item.
- **High-side gate-driver loss model upgrade** — currently treated
  as ideal; adding bootstrap-cap dynamics and dead-time analysis
  would tighten Vo prediction.

After Track A lands, update `project_converter_model_quality_tiers.md`
to add Zeta at DAB-quality tier, and update `FUTURE_TOPOLOGIES.md`
to remove Zeta from the "planned" list.

---

## 9. References

- **TI SLYT372** — *Designing DC/DC Converters Based on ZETA Topology*
  (the canonical Zeta app note)
- **TI SLVAFJ6** — *How to Approach a Power-Supply Design Pt 4: Zeta
  Equations*
- **TI SNVA608** — *LM5085: Designing Non-Inverting Buck-Boost (Zeta)
  with a Buck P-FET*
- **TI PMP9581_REVC** — Zeta reference design, uncoupled
- **TI PMP9581_REVB** — Zeta reference design, coupled-inductor
- **TI SLYT411** — *Benefits of a coupled-inductor SEPIC* (the
  ripple-cancellation condition, applies identically to coupled Zeta)
- **TI SNVA168E** — *Designing A SEPIC Converter* (analogous sizing
  formulae)
- **ADI** — *The Low Output Voltage Ripple Zeta DC/DC Converter
  Topology* (web article)
- **ADI LT8711** datasheet — multi-topology controller (Zeta, SEPIC,
  inverting buck-boost)
- **MPS** — *Zeta Converters* tutorial
- **Electronics Weekly** — *Zeta: that other DC/DC topology*
- **Erickson & Maksimović**, *Fundamentals of Power Electronics*,
  3rd ed., Ch. 5 (DCM), Ch. 6 (gain derivation), Ch. 12 (small-signal,
  RHPZ analysis)
- **MDPI Sensors 21/22/7434** — *SEPIC/Cuk/Zeta DCM modelling*
  (PMC8618931); covers DCM1/DCM2 sub-modes
- **IEEE 8343387** — *Bridgeless Zeta PFC for LED driver*
- **WSEAS 2010** — *Isolated Zeta Converter principles & design*
- **IEEE 6506052** — *Asymmetric multiphase interleaved DC-DC ripple
  cancel*
- **Coilcraft** — *Guide to Coupled Inductors*
- **ICnavigator** — *Single-Inductor Buck-Boost (SEPIC/Zeta) CV/CC
  limits*
