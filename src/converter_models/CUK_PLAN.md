# Cuk Converter — Topology Add Plan

> Purpose: add a new `Cuk` topology to `src/converter_models/` (none exists
> today), wired into the MAS schema, with full ngspice harness and tests,
> meeting **DAB-quality** as defined by `CONVERTER_MODELS_GOLDEN_GUIDE.md`.
>
> Read `CONVERTER_MODELS_GOLDEN_GUIDE.md` first; this plan layers on top.
> Sibling references:
> - `Boost.{h,cpp}` — single-inductor PWM, simplest analytical pattern
> - `Flyback.{h,cpp}` — transformer-isolated PWM, useful for the isolated-Cuk variant
> - `IsolatedBuckBoost.{h,cpp}` — transformer + Lo, the closest 2-magnetic shape
> - `DifferentialModeChoke.{h,cpp}` — coupled-inductor magnetic pattern (zero-ripple Cuk maps here)
> - `CONVERTER_MODELS_GOLDEN_GUIDE.md`, `CLLC_REWRITE_PLAN.md`
>
> Cuk is **not** a resonant converter; it's a 4th-order PWM converter.
> Solver is analytical piecewise-linear (no Newton, no matrix exponential).
> The hard part is **SPICE convergence** (4th-order series LC with hard
> switching) and the **magnetic structure** (single inductor vs. coupled
> 2-winding vs. isolation transformer).

---

## 1. Why Cuk in MKF

- Continuous input AND output current → the lowest-EMI non-resonant
  topology. Customers asking for PV / LED / audio / EMI-sensitive
  designs are the natural users.
- Coupled-inductor zero-ripple variant exercises a magnetic-design
  problem MKF currently has no canonical example for: a 2-winding
  magnetic with a *constrained* leakage inductance (neither pure
  transformer nor pure inductor).
- Isolated Cuk gives a flyback-class topology with **no DC bias** in
  the transformer (Ca, Cb block DC) — design-target distinct from
  flyback.
- Buck, Boost, Flyback, IsolatedBuckBoost, SEPIC-like coverage is
  already in MKF; Cuk fills the inverting / 4th-order gap.

---

## 2. Variants in scope

From the research dossier (Cuk PESC 1977; Erickson-Maksimović 3rd ed.
§2.4, §5, §11, §15.5; Cuk patents US4,257,087 and US4,355,352; TI
LM2611; MDPI Sustainability/WEVJ/Energies 2022–2023; Nature Sci. Rep.
2022):

1. **Non-isolated Cuk (V1)** — original Slobodan Cuk topology. Two
   *separate* inductors L1 and L2, coupling cap C1, switch S, diode D.
   Inverting (Vo < 0). CCM operation.
2. **Coupled-inductor Cuk (V2 — "zero-ripple")** — L1 and L2 wound on
   a common core; mutual M chosen so input ripple ΔIL1 → 0 (or output
   ripple, or both with extra trim winding). 2-winding magnetic with
   constrained leakage.
3. **Isolated Cuk (V3)** — splits C1 into Ca (primary) + Cb (secondary)
   around an isolation transformer. Vo polarity becomes free; turns
   ratio adds an extra gain knob. 1 transformer + 2 separate inductors
   (or a single integrated-magnetics core).
4. **Synchronous Cuk (V4)** — diode → MOSFET; same equations, just
   netlist swap. Implemented as a flag on V1.
5. **Bidirectional Cuk (V5)** — both switches active, regen capable.
   Implemented as a flag on V1 with an extra OP `powerFlow` field.

**Out of scope for the initial add** (defer to a follow-up):

- DCM (DICM) gain expressions wired into the design solver — V1 will
  enforce strict CCM in `process_design_requirements`; DCM is
  detected and flagged but not designed for.
- Tapped-inductor / high-gain Cuk (extra turns ratio in L2).
- Three-phase / interleaved Cuk.
- DVM (discontinuous-cap-voltage mode) — exotic, rarely used.

---

## 3. Reference designs as test fixtures

| # | Source | Vin | Vo | Pout | fsw | D | L1 | L2 | C1 | Cin | Co | η | Use |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Erickson-Maksimović 3rd ed. §2.4 worked example | 25 V | -25 V | 25 W | 100 kHz | 50 % | 100 µH | 100 µH | 1 µF | 4.7 µF | 10 µF | ideal | unit-test for CCM equations |
| 2 | Simon Bramble LT3757-based | 10 V | -5 V | 5 W | 300 kHz | 33 % | 47 µH | 47 µH | film, 15 V | small | 3.3 µF MLCC | 85 % | small-power smoke test |
| 3 | TI LM2611 typical app (1.4 MHz) | 3.3 V | -5 V | 0.5 W | 1.4 MHz | ~60 % | 10 µH | 10 µH | 0.22-1 µF | 4.7 µF | 10-22 µF | ~80 % | very-high-fsw SPICE convergence test |
| 4 | MDPI Sustainability 15 (2023) "New Cuk DC-DC" | 60 V | -120 V | 200 W | 50 kHz | 67 % | 1 mH | 0.5 mH | 4.7 µF film | 10 µF | 47 µF | 96 % | mid-power CCM reference |
| 5 | MDPI WEVJ 14 (2023) Flicker-Free LED isolated Cuk | 90-270 Vac | 50 Vdc | 50 W | 65 kHz | varies | integrated | integrated | 470 nF + 470 nF film (Ca, Cb) | EMI filter | 220 µF | ~92 % | isolated variant reference |
| 6 | Nature Sci. Rep. 2022 (s41598-022-17801-z) — bidirectional tapped Cuk | 48 V ↔ 380 V | bidir | 1 kW | 20 kHz | 0.7 fwd | 500 µH (tapped) | 200 µH | 6.8 µF film | 47 µF | 470 µF | 95 % buck / 91 % boost | bidirectional reference (out of initial scope, but spec'd for V5) |

URLs (for citation in header docstring):

- Cuk patent US4,257,087 (zero-ripple integrated magnetics): https://patents.google.com/patent/US4257087
- Cuk patent US4,355,352: https://patents.google.com/patent/US4355352
- Cuk "A Zero Ripple Technique Applicable To Any DC Converter" (~1983): https://static.elitesecurity.org/uploads/3/6/3617667/azrtatad.pdf
- Erickson-Maksimović "Fundamentals of Power Electronics" 3rd ed. (2020), §2.4, §5, §8, §11, §15.5 (textbook, no URL)
- Simon Bramble Cuk walkthrough: http://www.simonbramble.co.uk/dc_dc_converter_design/inverting_converter/inverting_dc_to_dc_converter_design.html
- TI LM2611: https://www.ti.com/product/LM2611
- TI SLUP084 (RHP-zero primer): https://e2e.ti.com/cfs-file/__key/communityserver-discussions-components-files/196/slup084.pdf
- MDPI Sustainability 2023: https://www.mdpi.com/2071-1050/15/11/8515
- MDPI WEVJ 2023 (LED isolated Cuk): https://www.mdpi.com/2032-6653/14/3/75
- MDPI Energies 2022 (SEPIC/Cuk/Zeta MPPT): https://www.mdpi.com/1996-1073/15/21/7936
- Nature Sci. Rep. 2022 (bidirectional tapped Cuk): https://www.nature.com/articles/s41598-022-17801-z

---

## 4. Equations (canonical block for the new header docstring)

Notation: `D` = duty cycle, `Tsw = 1/fsw`, `Iin = IL1,avg`,
`Iout = IL2,avg` (positive when flowing out of the negative output
terminal toward ground), `R = |Vo|/Iout`. All voltages magnitudes
unless explicitly signed; in the topology output `Vo < 0`.

### 4.1 CCM steady state (ideal, lossless)

```
(E1)  Conversion ratio:        M(D) = Vo/Vin = -D/(1-D)
(E2)  Coupling cap voltage:    VC1 = Vin + |Vo| = Vin/(1-D)
(E3)  Switch peak voltage:     VS,peak = VC1 = Vin + |Vo|
(E4)  Diode peak rev voltage:  VD,peak = VC1 = Vin + |Vo|
(E5)  Average inductor curr:   IL1,avg = Iin = Iout·D/(1-D)
                               IL2,avg = Iout
(E6)  Switch peak current:     IS,peak ≈ IL1,avg + IL2,avg + (ΔIL1+ΔIL2)/2
(E7)  Diode peak current:      ID,peak ≈ IL1,avg + IL2,avg + (ΔIL1+ΔIL2)/2
(E8)  Switch RMS current:      IS,rms ≈ (Iin+Iout)·√D
(E9)  Diode RMS current:       ID,rms ≈ (Iin+Iout)·√(1-D)
(E10) C1 RMS current:          IC1,rms ≈ √(D·(1-D))·(Iin+Iout)   ← film/MLCC required
```

### 4.2 Sizing (CCM)

```
(S1)  L1   ≥ Vin·D / (ΔIL1·fsw)            ; ΔIL1 = 20-40% of IL1,avg
(S2)  L2   ≥ |Vo|·(1-D) / (ΔIL2·fsw)        ; ΔIL2 = 20-40% of IL2,avg
(S3)  C1   ≥ Iout·D / (ΔVC1·fsw)            ; ΔVC1 ≤ 5% of VC1, FILM cap only
(S4)  Cin  ≥ ΔIL1 / (8·fsw·ΔVin)
(S5)  Co   ≥ ΔIL2 / (8·fsw·ΔVo)
```

### 4.3 CCM/DCM boundary

```
(B1)  Le        = L1·L2 / (L1 + L2)            (parallel combination)
(B2)  K(D)      = 2·Le·fsw / R
(B3)  Kcrit(D)  = (1-D)²
(B4)  CCM if K(D) > Kcrit(D);  DCM (DICM) otherwise
(B5)  DCM gain magnitude:  |M(D,K)| = D / √K
```

### 4.4 Small-signal poles / zeros (control-loop limits)

```
(P1)  ω0,a² ≈ 1/(L1·C1)                       (input-side LC)
(P2)  ω0,b² ≈ 1/(L2·Co)                       (output-side LC)
(P3)  ωRHP  ≈ R·(1-D)² / (D·L2)               (dominant RHP zero)
      Loop bandwidth limit:  fc ≤ ωRHP/(2π·5)
```

Cuk has **two** RHP zeros in the vc/d transfer function (Erickson §8).
The dominant one gates loop bandwidth; the second is typically a
decade higher.

### 4.5 Coupled-inductor zero-ripple (V2)

Define `n_eff = √(L2/L1)` and `k = M/√(L1·L2)`.

```
(Z1)  Zero input ripple:  n_eff · k = 1
                          ⇔ M = L1
                          ⇔ leakage referred to L1 winding = 0
```

Implementation: single E-core or pot core, two windings, asymmetric
gap (single outer leg) — symmetric gapping yields equal leakage on
both sides and prevents zero ripple. References: Cuk patent
US4,257,087; Erickson §15.5.

### 4.6 Isolated Cuk (V3)

```
(I1)  Vo/Vin               = -(N2/N1) · D/(1-D) = -n · D/(1-D)
(I2)  VCa,avg              = Vin/(1-D)
(I3)  VCb,avg              = |Vo|/D
(I4)  VS,peak              = Vin + VCa = Vin/(1-D)
(I5)  VD,peak              = |Vo| + VCb = |Vo|/D
(I6)  Transformer DC bias  = 0   (Ca, Cb block DC — automatic flux reset)
```

The transformer is sized for low magnetizing current (ungapped or
minimally gapped), **not** for energy storage. This is the
load-bearing distinction vs. flyback.

---

## 5. Solver: piecewise-linear analytical (no TDA)

Cuk is a PWM converter; its sub-intervals are LTI and the inductor
currents are exact-linear within each. No matrix exponential or Newton
needed (unlike LLC/CLLC). Pattern to follow: `Boost.cpp:35–88` and
`IsolatedBuckBoost::processOperatingPointsForInputVoltage`.

### 5.1 Sub-intervals (CCM, V1)

Per switching period:

| Sub-interval | Length | S | D | vL1 | vL2 | iC1 |
|---|---|---|---|---|---|---|
| ON  | D·Tsw    | closed | open  | +Vin     | -VC1+Vo  | -IL2 |
| OFF | (1-D)·Tsw | open  | closed | Vin-VC1 | +Vo      | +IL1 |

(With `Vo < 0`, "+Vo" in vL2 row reads as a negative number; the
slopes are correct.)

Steady-state algebra → `M = -D/(1-D)` and `VC1 = Vin/(1-D)`. Inductor
currents propagated by single multiply-add per sub-interval. Time
grid: `2·N_samples + 1` with `N_samples = 256`, span `[0, Tsw]`
exactly.

### 5.2 DCM detection

After CCM solve, compute `K(D) = 2·Le·fsw / R` and compare against
`(1-D)²`. If DCM:

- **For V1 process_operating_point_for_input_voltage**: still emit
  CCM-style waveforms but tag `lastOperatingMode = DCM` and populate
  `lastDcmK` and `lastDcmRatio = K/Kcrit`. Issue a one-time warning.
  Full DCM waveform construction is deferred (out of initial scope).
- **For run_checks**: warn (not error) if any OP is DCM.

### 5.3 Per-OP diagnostics (mandatory per Guide §2.4 / §4.6)

```cpp
mutable double lastDutyCycle = 0;
mutable double lastConversionRatio = 0;            // M(D) signed
mutable double lastCouplingCapVoltage = 0;         // VC1
mutable double lastSwitchPeakVoltage = 0;
mutable double lastSwitchPeakCurrent = 0;
mutable double lastDiodePeakReverseVoltage = 0;
mutable double lastDiodePeakCurrent = 0;
mutable double lastSwitchRmsCurrent = 0;
mutable double lastDiodeRmsCurrent = 0;
mutable double lastCouplingCapRmsCurrent = 0;      // drives film-cap selection
mutable double lastInputInductorRipple = 0;        // ΔIL1
mutable double lastOutputInductorRipple = 0;       // ΔIL2
mutable double lastRhpZeroFrequency = 0;           // fRHP
mutable double lastRecommendedLoopBandwidth = 0;   // fRHP/5
mutable int    lastOperatingMode = 0;              // 0 = CCM, 1 = DCM
mutable double lastDcmK = 0;
```

### 5.4 Coupled-inductor (V2) extension

Add `coupledInductor: bool` to schema. When true, the topology returns
a **single** 2-winding magnetic (DesignRequirements with isolation
sides covering both inductor windings on the same core). The leakage
constraint `M = L1` (E1 in §4.5) is published via
`get_extra_components_inputs(IDEAL/REAL)` rather than embedded in the
DesignRequirements (MKF doesn't have a "constrained leakage" field).

### 5.5 Isolated (V3) extension

Add `isolated: bool` and `turnsRatio: double` to schema. When isolated:

- DesignRequirements returns a 2-winding transformer (mirror
  `Flyback.cpp` pattern).
- Lm is sized for low magnetizing current: target
  `ΔILm,pp ≤ 5% of (Iin + n·Iout)`.
- `get_extra_components_inputs` emits **L1, L2, Ca, Cb** as
  ancillary components (4 entries, not 2 like CLLC).

---

## 6. SPICE netlist (generate_ngspice_circuit)

Cuk is famously a convergence troublemaker — 4th-order LC with hard
switching → ringing during commutation. The netlist below is built to
the Guide §5 contract, plus **mandatory IC pre-charge** (skipping
this guarantees divergence).

### 6.1 Topology — non-isolated V1

```
* Cuk Converter (non-isolated, CCM) — Generated by OpenMagnetics
*
* Vin─Vin_sense─L1─node_A──Vc1─C1──node_B──L2─out_sense─Vout(neg)
*                  │                  │
*                  S to GND          D to GND (cathode at node_B)
*                  │                  │
*                snubber           snubber
*                  │                  │
*                 GND                GND
*
* Iout flows out of Vout(neg) to GND through Rload
* Output sign: V(vout) is NEGATIVE
```

### 6.2 Mandatory netlist contents (Guide §5 + Cuk specifics)

- `Vin vin_p 0 {Vin}` and zero-volt `Vin_sense` for input current.
- `L1 vin_sense node_A {L1} ic={IL1_avg}`  ← **IC mandatory**
- `Vc1_sense node_A node_C 0`  (sense source for C1 current)
- `C1 node_C node_B {C1} ic={Vin/(1-D)}`  ← **IC mandatory**, voltage from `node_A` to `node_B` is `+VC1` so polarity is `node_C` (= `node_A` after sense) referenced to `node_B`.
- `L2 node_B out_sense {L2} ic={IL2_avg}`  ← **IC mandatory**
- Switch S: `S1 node_A 0 pwm 0 SW1` (drain at `node_A`, source at GND).
- Diode D: `D1 0 node_B DIDEAL` (anode at GND, cathode at `node_B`).
  In synchronous mode (V4): replace D1 with `S2 0 node_B pwm_inv 0 SW1`.
- Snubber across S: `Rsnub_S node_A 0 1k` and `Csnub_S node_A snub_S_int 1n`.
- Snubber across D: `Rsnub_D 0 node_B 1k` and `Csnub_D 0 snub_D_int 1n`.
- Output cap: `Co out_sense 0 {Co} ic={Vo}` (Vo is negative).
- Load: `Rload out_sense 0 {max(|Vo|/Iout, 1e-3)}`.
- Sense sources (zero-volt, exposing branch currents to .save):
  `Vin_sense`, `Vl1_sense` (or use the existing inline V), `Vc1_sense`,
  `Vl2_sense`, `Vd_sense`, `Vs_sense`, `Vout_sense`.
- Bridge-voltage probe (Guide §5 mandates one for every topology):
  `Esw sw_probe 0 VALUE={V(node_A)}` — exposes the switch-node
  voltage as `v(sw_probe)`.
- Models (verbatim, copy from `Boost.cpp:288–328`):
  ```
  .model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg
  .model DIDEAL D(IS=1e-12 RS=0.05)
  ```
- PWM source: `Vpwm pwm 0 PULSE(0 5 0 10n 10n {D*Tsw} {Tsw})`.
- Solver options (mandatory; copy from `Boost.cpp:327`, **add GEAR**):
  ```
  .options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500
  .options METHOD=GEAR TRTOL=7
  ```
- Time step: `period/200` baseline. Force `period/500` if `fsw > 1 MHz`
  (LM2611 fixture).
- `.save v(sw_probe) v(node_A) v(node_B) v(out_sense)
   i(Vin_sense) i(Vl1_sense) i(Vc1_sense) i(Vl2_sense)
   i(Vd_sense) i(Vs_sense) i(Vout_sense)`.
- `.tran {Tsw/200} {totalTime} {steadyStateStart} UIC` — UIC enables
  `.ic` directives.
- Comment header: topology, Vin/Vo/Pout/fsw/D/L1/L2/C1/Co, "generated
  by OpenMagnetics".
- **Add small parasitic series R** to L1, L2 (DCR ≈ 50 mΩ) and ESR
  to C1, Co (≈ 5 mΩ). Purely lossless reactive networks make ngspice
  unhappy. Either bake into L/C lines as separate R, or use ngspice
  `.model L1mod IND` with `RP=`.

### 6.3 Coupled-inductor V2 netlist diff

Replace the two separate L1, L2 lines with:

```
Lp vin_sense node_A {L1} ic={IL1_avg}
Ls node_B   out_sense {L2} ic={IL2_avg}
K1 Lp Ls {k}    ; k computed from zero-ripple condition or user-set
```

`k` defaults to `√(L1/L2)` (unity n_eff·k); user can override.

### 6.4 Isolated V3 netlist diff

```
* Primary side
L1 vin_sense node_A {L1} ic={IL1_avg}
Vca_sense node_A node_Ca 0
C_a node_Ca pri_p {Ca} ic={Vin/(1-D)}
S1 pri_p 0 pwm 0 SW1
* Transformer (coupled inductors)
Lpri pri_p pri_n {Lpri} ic=0
Lsec sec_p sec_n {Lsec} ic=0
K1 Lpri Lsec 0.9999
Vpri_n pri_n 0 0          ; primary return tied to GND
* Secondary side
Vcb_sense sec_n node_Cb 0
C_b node_Cb node_B {Cb} ic={|Vo|/D}
D1 sec_p node_B DIDEAL
L2 node_B out_sense {L2} ic={IL2_avg}
Co out_sense 0 {Co} ic={Vo}
Rload out_sense 0 {max(|Vo|/Iout, 1e-3)}
```

(Polarity of secondary winding tags determines whether Vo is positive
or negative; pick to match user convention.)

### 6.5 simulate_and_extract_operating_points

Follow Guide §6 verbatim:

1. Construct `NgspiceRunner`.
2. If `runner.is_available()` is false, warn and fall back to analytical.
3. For each (Vin, OP), build netlist, run, on `success==false` log
   error, fall back to `process_operating_point_for_input_voltage`,
   tag name with `[analytical]`.
4. On success, build `WaveformNameMapping` and call
   `NgspiceRunner::extract_operating_point`.

Mapping (V1):
- Primary winding (= L1): `voltage = node_A`, `current = vl1_sense#branch`
- For coupled-inductor V2 the same magnetic carries both windings → emit
  both Lp and Ls excitations on a single OperatingPoint.
- For isolated V3, primary = transformer Lpri, secondary = Lsec; L1
  and L2 emerge via `get_extra_components_inputs`.

### 6.6 simulate_and_extract_topology_waveforms

```cpp
wf.set_input_voltage(getWaveform("vin_sense"));   // = Vin
wf.set_input_current(getWaveform("vin_sense#branch"));
wf.get_mutable_output_voltages().push_back(getWaveform("out_sense"));
wf.get_mutable_output_currents().push_back(getWaveform("vout_sense#branch"));
```

---

## 7. Schema additions (`MAS/schemas/inputs/topologies/cuk.json`)

New file. Pattern off `boost.json` (single inductor) + `flyback.json`
(turns ratio + isolation):

```jsonc
{
    "title": "cuk",
    "description": "Inverting Cuk converter (non-isolated, coupled-inductor, or isolated)",
    "type": "object",
    "properties": {
        "inputVoltage":         { "$ref": "../../utils.json#/$defs/dimensionWithTolerance" },
        "diodeVoltageDrop":     { "type": "number" },
        "maximumSwitchCurrent": { "type": "number" },
        "currentRippleRatio":   { "type": "number" },
        "efficiency":           { "type": "number", "default": 0.95 },

        "variant": {
            "enum": ["nonIsolated", "coupledInductor", "isolated"],
            "default": "nonIsolated"
        },
        "synchronous":          { "type": "boolean", "default": false },
        "bidirectional":        { "type": "boolean", "default": false },

        "couplingCapacitance":  { "type": "number" },     // C1 (or Ca for isolated)
        "couplingCapacitanceSecondary": { "type": "number" },  // Cb, isolated only
        "outputCapacitance":    { "type": "number" },
        "inputCapacitance":     { "type": "number" },
        "couplingCoefficient":  { "type": "number", "default": 0.9999 },
                                 // V2: zero-ripple condition; V3: leakage spec

        "operatingPoints": {
            "type": "array",
            "items": { "$ref": "#/$defs/cukOperatingPoint" },
            "minItems": 1
        }
    },
    "required": ["inputVoltage", "diodeVoltageDrop", "operatingPoints"],
    "$defs": {
        "cukOperatingPoint": {
            "title": "cukOperatingPoint",
            "allOf": [
                { "$ref": "../../utils.json#/$defs/baseOperatingPoint" },
                {
                    "type": "object",
                    "properties": {
                        "powerFlow": {
                            "enum": ["forward", "reverse"],
                            "default": "forward"
                        }
                    }
                }
            ]
        }
    }
}
```

After schema edit, regenerate `MAS.hpp` and copy to `build/MAS/MAS.hpp`
(per user-memory staging gotcha).

Add `"CUK_CONVERTER"` to the `Topologies` enum in `MAS/schemas/inputs/`
(see how BUCK_CONVERTER and FLYBACK_CONVERTER are listed).

---

## 8. Class structure (`Cuk.h` / `Cuk.cpp`)

Mirror the IsolatedBuckBoost pattern (it's the closest 2-magnetic
sibling):

```cpp
namespace OpenMagnetics {
using namespace MAS;

class Cuk : public MAS::Cuk, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 20;

    // 2.2 Computed design values
    double computedInputInductance       = 0;   // L1
    double computedOutputInductance      = 0;   // L2
    double computedCouplingCapacitance   = 0;   // C1
    double computedOutputCapacitance     = 0;   // Co
    double computedDutyCycle             = 0;
    double computedDiodeVoltageDrop      = 0.6;

    // 2.3 Schema-less device parameters
    double mosfetCoss = 200e-12;

    // 2.4 Per-OP diagnostics (see §5.3)
    mutable double lastDutyCycle = 0;
    /* ... full list from §5.3 ... */
    mutable int    lastOperatingMode = 0;       // 0=CCM, 1=DCM

    // 2.5 Extra-component waveforms (cleared in process_operating_points)
    mutable std::vector<Waveform> extraL2VoltageWaveforms;
    mutable std::vector<Waveform> extraL2CurrentWaveforms;
    mutable std::vector<Waveform> extraC1VoltageWaveforms;
    mutable std::vector<Waveform> extraC1CurrentWaveforms;
    // For isolated V3: also extraCaWaveforms, extraCbWaveforms

public:
    bool _assertErrors = false;

    Cuk(const json& j);
    Cuk() {};

    // 2.7 Tuning + 2.8 computed accessors + 2.9 device params + 2.10 diagnostics
    /* boilerplate per Guide §2 — copy IsolatedBuckBoost shape */

    // 2.11 Topology interface
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // 2.12 Worker
    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const CukOperatingPoint& opPoint,
        double l1,
        double l2);

    // 2.13 Static analytical helpers
    static double compute_duty_cycle(double Vin, double Vo,
                                     double diodeDrop, double efficiency);
    static double compute_conversion_ratio(double D);          // = -D/(1-D)
    static double compute_coupling_cap_voltage(double Vin, double D);
    static double compute_l1_min(double Vin, double D, double dIL1, double fsw);
    static double compute_l2_min(double Vo,  double D, double dIL2, double fsw);
    static double compute_c1_min(double Iout, double D, double dVC1, double fsw);
    static double compute_dcm_K(double L1, double L2, double fsw, double R);
    static double compute_rhp_zero_freq(double R, double D, double L2);

    // 2.14 Extra-components factory
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // 2.15 SPICE
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

class AdvancedCuk : public Cuk {
    /* desiredL1, desiredL2, desiredC1, desiredCo, desiredTurnsRatios */
    /* Inputs process(); */
};

void from_json(const json& j, AdvancedCuk& x);
void to_json(json& j, const AdvancedCuk& x);
}
```

The "primary" magnetic returned by `process_design_requirements`
defaults to **L1** (input inductor) — single isolation side, no turns
ratio. L2, C1, Co are emitted via `get_extra_components_inputs`. For
V2 (coupled-inductor) the primary is the 2-winding magnetic; for V3
(isolated) the primary is the transformer, and L1/L2/Ca/Cb all flow
through `get_extra_components_inputs`.

---

## 9. Tests (`tests/TestCuk.cpp`) — DAB-quality bar

Mirror DAB and Boost test sets. Required `TEST_CASE`s per Guide §8:

| Test | Purpose | Tolerance |
|---|---|---|
| `Test_Cuk_Erickson_Reference_Design` | reproduce Erickson §2.4 worked example (25 V → -25 V, 100 kHz, D=0.5) | ±5% on D, M, VC1, IL1, IL2, ΔIL1, ΔIL2 |
| `Test_Cuk_Bramble_Reference_Design` | LT3757-based Cuk (10 V → -5 V, 300 kHz) | ±5% on key scalars |
| `Test_Cuk_Design_Requirements` | turns-ratio (1.0 for non-isolated), L1/L2/C1/Co positive and within range | round-trip a power target |
| `Test_Cuk_OperatingPoints_Generation` | multiple Vin (min/nom/max), antisymmetry doesn't apply (asymmetric duty), but: vL1·D + vL1·(1-D) = 0 (volt-second balance), iL1 piecewise linear, IL1,avg = D·Iout/(1-D) | algebraic |
| `Test_Cuk_CCM_DCM_Boundary` | sweep load; `lastOperatingMode` flips at K(D) = (1-D)² | sign change |
| `Test_Cuk_Conversion_Ratio` | `lastConversionRatio` matches `-D/(1-D)` for D ∈ {0.2, 0.5, 0.7} | exact |
| `Test_Cuk_RHP_Zero_Frequency` | `lastRhpZeroFrequency` matches Erickson §8 expression | exact |
| `Test_Cuk_Diagnostics_Populated` | every diagnostic in §5.3 is non-zero after a successful OP | non-null |
| `Test_Cuk_SPICE_Netlist` | netlist parses, contains snubbers, **METHOD=GEAR**, IC pre-charge for L1/L2/C1/Co, Esw probe, sense sources, no diodes if synchronous=true | string `find` |
| `Test_Cuk_DivideByZero_Guard` | Iout=0 → no crash; Rload clamps to 1e-3 | no throw |
| `Test_Cuk_PtP_AnalyticalVsNgspice_LowFsw` | Erickson 100 kHz fixture; primary-current (= IL1) NRMSE | **≤ 0.15** |
| `Test_Cuk_PtP_AnalyticalVsNgspice_HighFsw` | LM2611 1.4 MHz fixture | **≤ 0.15** |
| `Test_Cuk_Coupled_Inductor_ZeroRipple` | V2 with k = √(L1/L2); ΔIL1 in SPICE drops by ≥ 10× vs uncoupled | ratio test |
| `Test_Cuk_Isolated_TransformerNoDcBias` | V3; transformer DC-bias current ≤ 1% of IL1,avg | ratio |
| `Test_Cuk_Isolated_ConversionRatio` | V3; `M = -n·D/(1-D)` for n ∈ {0.5, 1.0, 2.0} | ±5% |
| `Test_Cuk_Synchronous_Variant` | V4; netlist contains S2 (no diode), simulation runs | smoke |
| `Test_Cuk_ExtraComponents_L2_C1_Co` | `get_extra_components_inputs(IDEAL/REAL)` returns 3 entries (or 4 for V3) | count + sanity |
| `Test_Cuk_Waveform_Plotting` | CSV dump under `tests/output/`; visual inspection | not asserted |
| `Test_AdvancedCuk_Process` | round-trip user-specified L1/L2/C1/Co | round-trip equal |

NRMSE helpers (`ptp_nrmse`, `ptp_interp`, `ptp_current`,
`ptp_current_time`) extracted to a shared header
`tests/TestTopologyHelpers.h` (this is also a CLLC plan deliverable —
do it once, both plans benefit).

---

## 10. Phased delivery

| # | Scope | Files | Tests added | Effort |
|---|---|---|---|---|
| **P1** | MAS schema (`cuk.json`, `Topologies` enum entry), regenerate MAS.hpp, stub `Cuk.h`/`Cuk.cpp` (constructors + run_checks only, throws on the rest), `from_json`/`to_json` for `AdvancedCuk` | schema, `Cuk.h`, `Cuk.cpp` | none | 0.5 d |
| **P2** | Static analytical helpers (`compute_duty_cycle`, `compute_conversion_ratio`, `compute_l1_min`, `compute_l2_min`, `compute_c1_min`, `compute_dcm_K`, `compute_rhp_zero_freq`), unit tests for each | `Cuk.cpp`, tests | `Conversion_Ratio`, `RHP_Zero_Frequency` | 0.5 d |
| **P3** | `process_design_requirements` (sizes L1, L2, C1, Co from current-ripple-ratio or maximum-switch-current; mirror `Boost::process_design_requirements`); `process_operating_point_for_input_voltage` (CCM piecewise-linear waveforms for L1, L2, plus diagnostics); `process_operating_points` (loop over input voltages × OPs) | `Cuk.cpp` | `Erickson_Reference_Design`, `Bramble_Reference_Design`, `Design_Requirements`, `OperatingPoints_Generation`, `Diagnostics_Populated` | 1.5 d |
| **P4** | DCM detection & flag (no DCM waveform construction yet); CCM/DCM boundary test | `Cuk.cpp`, tests | `CCM_DCM_Boundary` | 0.5 d |
| **P5** | SPICE netlist (V1 non-isolated): full Guide §5 contract + IC pre-charge for L1/L2/C1/Co + GEAR + snubbers + Esw probe + sense sources. `simulate_and_extract_operating_points` with analytical fallback. `simulate_and_extract_topology_waveforms`. | `Cuk.cpp` | `SPICE_Netlist`, `DivideByZero_Guard` | 1.5 d |
| **P6** | NRMSE acceptance: tighten until ≤ 0.15 for low-fsw and high-fsw fixtures. Failure-discovery phase — expect to fix 1-2 SPICE convergence bugs (ic mismatch, missing parasitic R) | tests + minor SPICE fixes | `PtP_AnalyticalVsNgspice_LowFsw`, `PtP_AnalyticalVsNgspice_HighFsw` | 1 d |
| **P7** | `get_extra_components_inputs` returning L2, C1, Co ancillary `Inputs` (IDEAL + REAL modes) | `Cuk.cpp` | `ExtraComponents_L2_C1_Co` | 0.5 d |
| **P8** | V2 coupled-inductor variant: schema flag, 2-winding DesignRequirements, K1 Lp Ls statement in netlist, zero-ripple condition publishing | `Cuk.cpp` | `Coupled_Inductor_ZeroRipple` | 1 d |
| **P9** | V3 isolated variant: transformer DesignRequirements (mirror Flyback), Ca/Cb sizing, isolated netlist (Lpri/Lsec/K1, Ca/Cb, no DC bias check) | `Cuk.cpp` | `Isolated_TransformerNoDcBias`, `Isolated_ConversionRatio` | 1.5 d |
| **P10** | V4 synchronous + V5 bidirectional variant: schema flags, netlist diff (replace D1 with S2 driven by `~pwm`), bidirectional `powerFlow` field | `Cuk.cpp` | `Synchronous_Variant` | 1 d |
| **P11** | AdvancedCuk class, header docstring rewrite (ASCII art, equations §4, references, accuracy disclaimer per Guide §3), waveform plotting test | `Cuk.h`, `Cuk.cpp`, tests | `AdvancedCuk_Process`, `Waveform_Plotting` | 1 d |

**Total ~10.5 working days** for V1+V2+V3+V4 coverage.

Dependencies: P1 → P2 → P3; P4 ⊥ P3; P5 needs P3; P6 needs P5; P7 ⊥
P5; P8/P9/P10 need P5; P11 last.

---

## 11. Acceptance criteria (per Guide §15)

- [ ] `ninja -j2 MKF_tests` builds clean.
- [ ] All `TestCuk` cases pass.
- [ ] Both `*_PtP_AnalyticalVsNgspice_*` tests run and NRMSE ≤ 0.15.
- [ ] `lastOperatingMode` flips CCM ↔ DCM at K(D) = (1-D)² boundary
      within ±5%.
- [ ] Schema documented; `Topologies::CUK_CONVERTER` added.
- [ ] Header has ASCII art, equations §4 with citations,
      accuracy disclaimer (Coss neglected, body-diode reverse-recovery
      neglected), **polarity convention note** (Vo < 0 in non-isolated),
      RHP-zero-aware loop-bandwidth note.
- [ ] Netlist parses, contains snubbers, GEAR, IC pre-charge for **all
      four** reactive elements (L1, L2, C1, Co), Esw probe, sense
      sources, no rectifier diodes when `synchronous=true`.
- [ ] `simulate_and_extract_operating_points` runs SPICE for real,
      falls back to analytical with `[analytical]` tag on failure.
- [ ] `simulate_and_extract_topology_waveforms` references only nodes
      the netlist actually emits.
- [ ] `get_extra_components_inputs` supports `IDEAL` and `REAL` for
      L2, C1, Co (and Ca, Cb in V3).
- [ ] V2 coupled-inductor variant achieves ≥ 10× ripple reduction in
      SPICE vs. uncoupled.
- [ ] V3 isolated variant achieves ≤ 1% transformer DC bias.

---

## 12. Common Cuk implementation mistakes (avoid these)

From the research dossier — explicit pitfalls flagged in the literature:

1. **Polarity sign error.** Vo is *negative*. Wiring the load to GND
   with the wrong end of Co inverts feedback. Internally model
   `Vo = -|Vo|` and document the convention in the header.
2. **Drifting into DCM at light load.** M departs from -D/(1-D),
   loop opens up. Either size for full-load CCM at min duty, or
   detect+flag DCM and warn (P4 does the latter; full DCM design is
   deferred).
3. **Wrong C1 technology.** IC1,rms ≈ √(D(1-D))·(Iin+Iout) — easily
   0.5–3 A RMS. Electrolytic dies in seconds. Always document
   "polypropylene film or X7R MLCC stack, 50–100% V derating" in the
   recommended-component output.
4. **Cold-start inrush.** Without `.ic` SPICE dumps Vin into a series
   LC tank → IL peak >> rated. Pre-charge L1, L2, C1, Co — non-negotiable.
5. **Loop bandwidth too high.** Two RHP zeros + complex pole pair →
   bandwidth > ωRHP/5 gives non-minimum-phase ringing on every load
   step. Publish `lastRecommendedLoopBandwidth = fRHP/5`.
6. **K = 1.0 for transformer in V3 SPICE** — singular factorization;
   always 0.999..0.9999.
7. **Treating isolated-Cuk transformer like a flyback transformer.**
   No DC bias (Ca, Cb block it) → ungapped or minimally gapped, low
   magnetizing current. Test this.
8. **Forgetting the MOSFET stress is Vin+|Vo|** (not max). Surface in
   `lastSwitchPeakVoltage` so it shows up in design reports.
9. **Coupled-inductor V2: gap on the wrong leg.** Symmetric gapping
   yields equal leakage on both sides → ripple appears on both. The
   zero-ripple condition is asymmetric — needs single-outer-leg gap.
   Document in the V2 fixture.
10. **No clamp on D**. As Vo → ∞, D → 1; the loop runs away. Add
    `D ≤ 0.9` in `process_design_requirements`; throw if any OP needs more.
11. **No parasitic R in SPICE**. Lossless reactive networks are
    badly conditioned in ngspice. Bake DCR ≈ 50 mΩ into L1/L2 and
    ESR ≈ 5 mΩ into C1/Co.

---

## 13. Magnetic structure summary (the bit MKF cares about most)

| Variant | "Primary" magnetic returned by topology | Extra components via `get_extra_components_inputs` |
|---|---|---|
| V1 non-isolated | L1 (single inductor, 1 isolation side) | L2, C1, Co |
| V2 coupled-inductor | 2-winding magnetic with constrained `M = L1` (zero-ripple); 1 isolation side covering both windings | C1, Co |
| V3 isolated | Transformer (2 windings, **0 DC bias**, 2 isolation sides) | L1, L2, Ca, Cb, Co |

V2 is the most novel from MKF's perspective — it requires publishing a
*leakage constraint* (`L_leak1 = 0` or `M = L1`) alongside the
DesignRequirements. Today `Topology::create_isolation_sides` doesn't
have a "constrained leakage" knob; investigate whether to:

- (a) add a new optional field `desired_leakage_inductance` to
  DesignRequirements (clean but cross-cutting),
- (b) only publish via `get_extra_components_inputs` (less canonical
  but fully self-contained in `Cuk.cpp`).

Recommend (b) for the initial add; option (a) becomes a follow-up if
SEPIC / Zeta / Forward-with-coupled-inductor land later.

---

## 14. Checklist for the implementing agent

Before each phase:

- [ ] Re-read `CONVERTER_MODELS_GOLDEN_GUIDE.md`.
- [ ] Verify the latest `Boost.cpp` SPICE pattern (closest single-inductor sibling).
- [ ] Verify the latest `Flyback.cpp` SPICE pattern (closest transformer sibling for V3).
- [ ] Verify the latest `IsolatedBuckBoost.cpp` for the "main magnetic + extra components" pattern.
- [ ] Verify the latest `DifferentialModeChoke.cpp` for coupled-inductor magnetic returns (V2).
- [ ] Build with `ninja -j2` in the build dir, not `cmake --build` (per user memory).
- [ ] After MAS schema edits, copy `MAS.hpp` from source to `build/MAS/MAS.hpp` (per user memory).
- [ ] Mark each phase as complete only when its NRMSE / sign-change / build-clean criteria are met.

After P11 (final):

- [ ] Update CLAUDE.md / user-memory entry on converter quality tiers
      to include Cuk as DAB-quality.
- [ ] Open follow-up issues for: full DCM waveform construction;
      tapped/high-gain Cuk; three-phase / interleaved Cuk; SEPIC and
      Zeta as Cuk siblings.
- [ ] Consider extracting shared 4th-order PWM helpers (Cuk, SEPIC,
      Zeta) into `PwmFourthOrderSolver.h` analogous to
      `PwmBridgeSolver.h`, when a second 4th-order topology lands.

---

## 15. Reference: equations cheat-sheet (paste into header docstring)

```
TOPOLOGY OVERVIEW (non-isolated Cuk, V1)
────────────────────────────────────────
            L1               C1                L2
   Vin ──┳──/\/\─────┳──┤├──┳──/\/\──┳── Vo (negative)
         │           │      │        │
         │          S│     D│       Co
         │           │      │        │
        Cin         GND    GND      GND ── Rload
         │
        GND

KEY EQUATIONS (CCM, ideal)
──────────────────────────
M(D)        = -D / (1 - D)                                     — Erickson §2.4
VC1         = Vin / (1 - D) = Vin + |Vo|                       — Erickson §2.4
IL1,avg     = D · Iout / (1 - D)            (= Iin)
IL2,avg     = Iout
VS,peak     = VD,peak = Vin + |Vo|
IS,peak     = ID,peak ≈ IL1,avg + IL2,avg + (ΔIL1+ΔIL2)/2
IC1,rms     ≈ √(D·(1-D)) · (Iin + Iout)                        — film-cap RMS

L1,min      ≥ Vin · D / (ΔIL1 · fsw)
L2,min      ≥ |Vo| · (1-D) / (ΔIL2 · fsw)
C1,min      ≥ Iout · D / (ΔVC1 · fsw)
Co,min      ≥ ΔIL2 / (8 · fsw · ΔVo)

K(D)        = 2 · Le · fsw / R,    Le = L1·L2/(L1+L2)
Kcrit(D)    = (1 - D)²
CCM if K > Kcrit, else DCM

ωRHP        ≈ R · (1 - D)² / (D · L2)                          — Erickson §8
fc,max      ≤ ωRHP / (2π · 5)

ZERO-RIPPLE (V2, coupled-inductor)
n_eff·k = 1, where n_eff = √(L2/L1), k = M/√(L1·L2)            — Cuk patent US4,257,087

ISOLATED (V3)
M           = -n · D / (1 - D)        (n = N2/N1)
VCa         = Vin / (1 - D)
VCb         = |Vo| / D
Transformer DC bias = 0 (Ca, Cb block DC)

REFERENCES
[1] Slobodan Cuk, "A New Optimum Topology Switching DC-to-DC Converter", IEEE PESC 1977.
[2] Cuk & Middlebrook, "A general unified approach to modelling switching-converter power stages", IEEE PESC 1976.
[3] Cuk, "A 'Zero' Ripple Technique Applicable To Any DC Converter", ~1983.
[4] Erickson & Maksimović, "Fundamentals of Power Electronics" 3rd ed. (2020), §2.4, §5, §8, §11, §15.5.
[5] Cuk patent US4,257,087 (1979) — integrated-magnetics zero-ripple Cuk.
[6] Cuk patent US4,355,352 (1982).
[7] TI LM2611 datasheet (SNOS965F).
[8] TI SLUP084 — Right-Half-Plane Zero primer.

ACCURACY DISCLAIMER
SPICE model omits MOSFET Coss variation with Vds, body-diode reverse
recovery, gate-driver propagation delay, and core-loss resistance.
Predicted η is an upper bound; expect 2-5% loss vs. measured at full load.

POLARITY CONVENTION
Vo < 0 in V1 / V2 (non-isolated). Iout reported as a positive scalar
flowing OUT of the negative output terminal toward GND. In V3
(isolated), Vo polarity is set by the secondary winding orientation
and may be positive.
```
