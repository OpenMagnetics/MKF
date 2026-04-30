# CLLC Bidirectional Resonant Converter — Rewrite Plan

> Purpose: bring `src/converter_models/Cllc.{h,cpp}` and `tests/TestCllc.cpp`
> up to **DAB-quality** as defined by `CONVERTER_MODELS_GOLDEN_GUIDE.md`.
> The current implementation is FHA-only with a diode rectifier and cannot
> physically run reverse power flow. This document is self-contained: a
> future agent should be able to execute it end-to-end without re-deriving.
>
> Read `CONVERTER_MODELS_GOLDEN_GUIDE.md` first; this plan layers on top.
> Reference siblings: `Llc.cpp` (closest topology — TDA solver to mirror),
> `Dab.cpp` (gold standard for diagnostics, NRMSE tests, SPICE quality).

---

## 1. Gap analysis (current Cllc.{h,cpp})

| # | Issue | Severity | Location |
|---|---|---|---|
| 1 | Waveform model is **pure FHA**: `Ip_pk·sin(ωs·t)` plus a free triangular Lm ramp. No real sub-interval propagation. Cannot reproduce sub-resonant freewheel, super-resonant truncation, asymmetric-tank gain, or any mode-dependent shape. | **Critical** | `Cllc.cpp:362–425` |
| 2 | **Diode-only secondary bridge** in SPICE → physically impossible to run reverse mode. Reverse "support" today is only a string suffix on the OP name. | **Critical** | `Cllc.cpp:661–665` |
| 3 | No piecewise propagator → `Im_peak = Vin·(T/2−td)/(2·Lm)` ignores that during freewheel `iLm = iLr` (Lm participates in resonance). Predicts wrong magnetizing waveform sub-resonant. | **Critical** | `Cllc.cpp:304` |
| 4 | SPICE missing snubbers (1 kΩ‖1 nF), `METHOD=GEAR`, IC pre-charge, sense sources `Vpri_sense`/`Vsec1_sense_o1`, Evab VCVS probe (Guide §5). | **High** | `Cllc.cpp:680–696` |
| 5 | Throws on simulation failure instead of analytical fallback (Guide §6). | **High** | `Cllc.cpp:752, 845` |
| 6 | No `get_extra_components_inputs` — Cr1, Cr2, Lr1, Lr2 cannot round-trip through MAS as ancillary magnetics. | **High** | missing |
| 7 | No diagnostics (`lastZvsMargin*`, `lastMode`, `lipFrequency`, `lastSteadyStateResidual`) — Guide §2.4, §4.6. | **High** | missing |
| 8 | No analytical-vs-ngspice NRMSE test → cannot certify ≤ 0.15 (Guide §8 marker). | **High** | missing |
| 9 | `R_load = Vout/Iout` with no divide-by-zero guard. | Medium | `Cllc.cpp:591` |
| 10 | Q convention undocumented (3 conventions in literature). | Medium | `Cllc.h:20–33` |
| 11 | Hardcoded magic constants (`Cout=100µF`, `Rdc_sec=1G`, `Rsnub=1MEG`, `Vgnd_sec` shorting `vout_n` to ground). | Medium | `Cllc.cpp:657–680` |
| 12 | Method misnamed (Guide says singular: `process_operating_point_for_input_voltage`). | Low | `Cllc.h:151` |
| 13 | Asymmetric a, b enter sizing but not the FHA gain function — gain returns the symmetric expression even when `symmetricDesign=false`. | Medium | `Cllc.cpp:189–217` |
| 14 | Header lacks ASCII art, accuracy disclaimer, mode table, asymmetric-Q disambiguation (Guide §3). | Medium | `Cllc.h:35–53` |
| 15 | No half-bridge variant (Guide §9 conventions allow it via a bridgeType field). | Medium | missing |

---

## 2. Variants in scope

From the research dossier (Infineon AN-2024-06, KIT 20 kW thesis, Wei TDA,
Min IEEE TPE 2021, Bartecka MDPI Energies 2024, Rezayati thesis 2023, TI
TIDM-02002, Lin 3.3 kW, Jung-Kim TPE 2013):

1. **Full-bridge symmetric CLLC** (`a=b=1`) — primary baseline; gain
   reduces to LLC under FHA. Most common (V2G, OBC).
2. **Full-bridge asymmetric CLLC** (`a≠1`, `b≠1`, `a·b≈1` so `fr1=fr2`) —
   different forward/reverse gain envelopes for wide-range battery V.
3. **Half-bridge CLLC** — sub-1 kW (server PSU, V2L). Bridge-voltage
   factor 0.5 like LLC.
4. **CLLLC** (extra discrete Lr2 inductor on secondary, vs. CLLC where
   Lr2 = secondary leakage) — TI TIDM-02002, KIT 20 kW. Same equations,
   different physical layout for SPICE / extra-components.
5. **Reverse power flow** for any of the above — secondary is inverter,
   primary is active rectifier. Identical TDA matrices, just swap source
   / load and turns ratio.

**Out of scope for this rewrite:** three-phase / interleaved CLLC.

---

## 3. Reference designs as test fixtures

| Fixture | Vin | Vout | Pout | fr | n | Lr | Lm | Cr | Use |
|---|---|---|---|---|---|---|---|---|---|
| Infineon REF-DAB11KIZSICSYS (UG-2020-31) | 750 V | 600 V | 11 kW | 73 kHz | 1.25 | ~36 µH | ~160 µH | ~132 nF | already in `TestCllc.cpp:807–849` — keep, tighten to ±5% |
| KIT 20 kW CLLLC (asymmetric) | mid-bus | mid-bus | 20 kW | fr1=146 kHz, fr2=117.6 kHz | 7/9 ≈ 0.78 | LR1=LR2=15 µH | LP=15 µH | CR1=79.2 nF, CR2=114.4 nF | **NEW** — asymmetric reference |
| Lin 3.3 kW symmetric | 400 V | 240–420 V | 3.3 kW | 200 kHz | 1.0 | symmetric | — | symmetric | **NEW** — symmetric V2G |
| TI TIDM-02002 CLLLC | 380–600 V | 280–450 V | 6.6 kW | 500 kHz nom | — | symmetric CLLLC | — | — | **NEW** — high-fr SiC |
| Small-power 400 V→48 V 500 W 200 kHz | 400 V | 48 V | 480 W | 200 kHz | 8.33 | — | — | — | already in `TestCllc.cpp:93–120` — keep for SPICE smoke tests |

URLs (for citation in header docstring):

- Infineon AN-2024-06: https://www.infineon.com/dgdl/Infineon-Operation-and-modeling-analysis-of-a-bidirectiona-CLLC-converter-ApplicationNotes-v01_00-EN.pdf
- Infineon UG-2020-31 / REF-DAB11KIZSICSYS: https://www.infineon.com/dgdl/Infineon-UG-2020-31_REF_DAB11KIZSICSYS-UserManual-v01_01-EN.pdf
- TI TIDM-02002 / TIDUEG2C: https://www.ti.com/tool/TIDM-02002
- Bartecka MDPI Energies 2024: https://www.mdpi.com/1996-1073/17/1/55
- KIT 20 kW CLLLC: https://publikationen.bibliothek.kit.edu/1000123579/151512204
- Jung & Kim IEEE TPE 28(4) 2013: https://ui.adsabs.harvard.edu/abs/2013ITPE...28.1741J
- Min et al. IEEE TPE 36(6) 2021 (asymmetric): https://ieeexplore.ieee.org/document/9241436
- Rezayati PhD thesis 2023 (state-plane): https://theses.hal.science/tel-04094415
- Wei et al. TDA: https://www.researchgate.net/publication/335499530
- Liu et al. SR via Lr voltage IEEE TPE: https://ieeexplore.ieee.org/iel7/63/4359240/09497655.pdf

---

## 4. Equations (canonical block for the new header docstring)

Notation: `ωs = 2π·fs`, `ωr = 1/√(Lr·Cr)`, `ωn = ωs/ωr`, `n = Np/Ns`,
`k = Lm/Lr`, **`Q = √(Lr/Cr) / Ro`** (Infineon convention — make explicit).

### 4.1 Sizing (Infineon AN-2024-06 §2.3, 11-step procedure)

```
Step 1   n        = Vin_nom / Vout_nom                       (full-bridge)
                  = Vin_nom / (2·Vout_nom)                   (half-bridge)
Step 2   Mg_min   = n·Vout_min / Vin_max
         Mg_max   = n·Vout_max / Vin_min
Step 3   pick k ∈ [3, 6], Q ∈ [0.2, 0.5] so M(ωn,Q,k) covers
         [Mg_min, Mg_max] within [fs_min, fs_max]
Step 4   Ro       = (8·n²/π²) · Vout² / Pout                 (FHA-equivalent AC load)
Step 5   Cr       = 1 / (2π·Q·fr·Ro)
Step 6   Lr       = 1 / ((2π·fr)²·Cr)
Step 7   Lm       = k · Lr
Step 8   sym :    L2 = Lr/n²,           C2 = n²·Cr           (a=b=1)
         asym:    L2 = a·Lr/n²,         C2 = n²·b·Cr,  a·b≈1
Step 9   ZVS:     Lm ≤ td / (16 · Coss_eq · fs_max)         (Infineon AN Eq 28)
```

### 4.2 FHA voltage gain

Symmetric (`a = b = 1`):

```
M(ωn, Q, k) = 1 / √( (1 + 1/k - 1/(k·ωn²))²
                   + Q²·(ωn - 1/ωn)² )                       (Infineon AN Eq 2)
```

Asymmetric: closed-form per Min IEEE TPE 2021 — `M_fwd(ωn) ≠ M_rev(ωn)`.
The full expression has 6 free parameters (Q, k, ωn, a, b, plus
direction). Implement as a complex-impedance numerical evaluation —
already done correctly in `Cllc.cpp:189–217` for the symmetric case.

### 4.3 Per-sub-interval natural frequencies

```
Power-delivery (rectifier on, Lm clamped to ±n·Vout):
   ω_P = 1 / √( (Lr1 + Lr2/n²) · (Cr1 · n²·Cr2) / (Cr1 + n²·Cr2) )
   At symmetric a=b=1 this equals  ωr  exactly.

Freewheel (rectifier off, Lm in series):
   ω_F = ω_P / √(1 + k)   ≈   1 / (2π·√((Lr+Lm)·Cr))   (= fm)
```

### 4.4 LIP / mode boundaries

- LIP frequency (load-independent point): `f_LIP = ω_P / (2π)`. At
  `fs = f_LIP` the FHA gain is ≈ 1 regardless of load. Use as the
  natural design centre.
- Sub-resonant (boost):     `fs ≤ f_LIP`,  modes 1–3 (Nielsen-style)
- Super-resonant (buck):    `fs > f_LIP`,  modes 4–6
- Capacitive boundary:      `fs > f_capacitive(Q,k)` — forbidden, ZVS lost

---

## 5. Time-domain solver (the heart of the rewrite)

Mirror the LLC pattern (`Llc.cpp:272–490`) with a wider state vector:

```cpp
struct CllcStateVector {
    double iLr1;   // primary resonant inductor current
    double vCr1;   // primary resonant capacitor voltage
    double iLm;    // magnetizing current
    double iLr2;   // secondary resonant inductor current (referred to sec)
    double vCr2;   // secondary resonant capacitor voltage
};

enum class CllcSubState {
    P_POS,   // power delivery, Vp = +Vin, Vs = +n·Vout (rectifier ON)
    P_NEG,   // power delivery, Vp = -Vin, Vs = -n·Vout
    F_POS,   // freewheel,      Vp = +Vin, rectifier OFF (iLr ≈ iLm)
    F_NEG    // freewheel,      Vp = -Vin, rectifier OFF
};
```

### 5.1 State-space matrices

For each sub-state `i`, the system is `dx/dt = Aᵢ·x + Bᵢ·u`, with
`u = [Vin, n·Vout]`. The Aᵢ are LTI; propagation per sub-interval is
**closed form** via matrix exponential `x(t+Δt) = expm(Aᵢ·Δt)·x(t) +
(integral term)` — never Euler.

Recommended derivation reference: Wei et al., IEEE PEAC TDA paper. The
2-stage (P, F) reduction collapses the 5-state problem to a chain of
2×2 / 3×3 oscillator solves with linear coupling to `iLm` — closed-form
per stage, no numerical integration required.

For symmetric `a=b=1` the cap chain `Cr1` and `n²·Cr2` collapses into a
single equivalent `Cr_eq = Cr1·(n²·Cr2)/(Cr1+n²·Cr2) = Cr1/2` and the
inductor chain `Lr1 + Lr2/n²` into `Lr_eq = 2·Lr1`. Verify symmetric
TDA against this collapsed form as a sanity test before tackling
asymmetric.

### 5.2 Sub-interval enumeration

Per half-cycle the ordered sequence is determined by `fs` vs `f_LIP`
and the load:

| Region | fs | Sequence (one half-cycle) |
|---|---|---|
| Sub-resonant boost | `fs ≤ fm` | `[P_POS, F_POS]` (rectifier completes its half-sine before the next switching edge) |
| Mid-resonant | `fm < fs ≤ f_LIP` | `[P_POS, F_POS]` with shorter `F_POS` |
| At resonance | `fs = f_LIP` | `[P_POS]` (no freewheel, exact half-sine) |
| Super-resonant buck | `fs > f_LIP` | `[P_POS]` truncated; rectifier turns off hard at next switching edge |
| Capacitive (forbidden) | `fs > f_cap` | rectifier polarity reverses — flag and refuse |

### 5.3 Steady state via 2D Newton on antisymmetry

`(vCr1(0), iLr1(0))` are the two unknowns. Other initial conditions
follow:

- `iLm(0) = iLr1(0)` at the start of `F` if the previous sub-interval
  ended in freewheel (else derived from charge balance).
- `vCr2(0) = -n²·vCr1(0)` for symmetric tanks (mirror symmetry).
- `iLr2(0) = -n·iLr1(0)` at antisymmetric start (full-bridge symmetric).

Residual `F(x0) = x(Ts/2) + x(0)` (half-wave antisymmetry: state at
`Ts/2` should equal `-x(0)`). Solve `F(x0) = 0` with damped Newton;
populate `lastSteadyStateResidual = ‖F(x0)‖₂` for diagnostics.

### 5.4 Diagnostics to populate per OP

```cpp
mutable int    lastMode;                  // 1 = sub-res, 2 = at res, 3 = super-res
mutable double computedLipFrequency;      // f_LIP = ω_P / (2π)
mutable double lastZvsMarginPrimary;      // iLr(switching edge) - iLm_threshold  (>0 = ZVS)
mutable double lastZvsMarginSecondary;    // for reverse mode
mutable double lastResonantTransitionTime;// dead-time fraction the bridge node takes
mutable double lastSteadyStateResidual;
mutable double lastPrimaryPeakCurrent;
mutable double lastResonantCapPeakVoltage;
mutable std::vector<int> lastSubStateSequence;
```

### 5.5 Reverse mode

Use the **same** solver, swap `Vin ↔ n·Vout`, and pass `direction =
REVERSE` to the SPICE driver so it gates the secondary bridge as the
inverter. The state-space matrices are reciprocal so no new code is
needed beyond a sign flip on `u`.

---

## 6. SPICE netlist (generate_ngspice_circuit)

### 6.1 Topology

```
* Forward-mode CLLC full-bridge bidirectional
*
* Vin─┬─S1───node_A─Cr1─Lr1─Lpri─node_B─S3──┐
*     │                  ║K              │
*     S2  (snubber)    Lm Lsec   (snub) S4  GND_pri
*     │                  ║              │
*     └──────────────────┴──────────────┘
*
*     Lsec─Lr2─Cr2─node_C─Sa─┬─out_node_o1─Cout─Rload─┐
*                            │                         │
*                            Sb         (snub)        GND_sec
*                            │                         │
*                  node_D─Sc─┴─Sd──────────────────────┘
*
* PWM: pwm_p1 → S1, S4   (pwm_p2 = ¬pwm_p1 with td)
* PWM: pwm_p2 → S2, S3
* SR (forward): pwm_s1 = pwm_p1, pwm_s2 = pwm_p2  (synchronous)
* SR (reverse): pwm_s* drives the secondary as inverter, primary becomes rectifier
```

### 6.2 Mandatory netlist contents (Guide §5 + CLLC specifics)

- 4× active full-bridge primary switches (`SW` model
  `VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg`).
- **4× active synchronous-rectifier switches on secondary** (NOT
  diodes — required for reverse). Gating: external clock 180° offset
  from primary (forward) or B-source ideal-diode emulator
  `B_Sa G S V=if(V(D,S)>0.1, 5, 0)` (works across all modes).
- Snubber `1k + 1nF` across **every** switch (8 total).
- Discrete `Lr1`, `Cr1` on primary; discrete `Lr2`, `Cr2` on secondary.
  Treat as separate components even when integrated; let the user
  later set `Lr1 = transformer leakage` via `extra_components`.
- Coupled inductors `K Lpri Lsec 0.9999` (never 1.0).
  `Lm = Lpri` by construction. If user sets `discreteLm=true`, emit
  `Lm` as a separate parallel inductor with its own K to a dummy
  winding.
- Sense sources: `Vpri_sense`, `Vsec1_sense_o1`, `Vout_sense_o1` (all
  zero-volt).
- Bridge-voltage probe: `Evab vab 0 VALUE={V(node_A) - V(node_B)}`.
- IC pre-charge: `.ic v(out_node_o1)=Vo v(cr1)=0 v(cr2)=0`.
- `R_load = max(Vo/Io, 1e-3)`.
- Solver options (verbatim):
  ```
  .options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500
  .options METHOD=GEAR TRTOL=7
  ```
- Time step: `period/200` baseline; allow `period/500` override when
  `fs > 1.5·fr` (super-resonant ringing).
- `.save v(vab) v(node_a) v(node_b) v(node_c) v(node_d)
  i(Vpri_sense) i(Vsec1_sense_o1) i(Vout_sense_o1)
  v(out_node_o1) v(cr1) v(cr2)`.
- Comment header: topology, Vin/Vout/Pout/fs/n/Lr/Lm/Cr, "generated by
  OpenMagnetics".

### 6.3 simulate_and_extract_operating_points

Follow Guide §6:

1. Construct `NgspiceRunner`.
2. If `runner.is_available()` is false, warn and fall back to analytical.
3. For each (Vin, OP), build netlist, run, on `success==false` log
   error, fall back to `process_operating_point_for_input_voltage`,
   tag name with `[analytical]`.
4. On success, build `WaveformNameMapping` and call
   `NgspiceRunner::extract_operating_point`.

DAB lines 1218–1307 are the canonical pattern.

### 6.4 simulate_and_extract_topology_waveforms

Output one entry per simulation with `v(vab)` as input voltage,
`i(Vpri_sense)` as input current, `v(out_node_o1)` as output voltage,
`i(Vout_sense_o1)` as output current. Clean up the existing
`v(pri_trafo_in)` references — they will not exist in the new netlist.

---

## 7. Schema additions (`MAS/schemas/inputs/topologies/cllcResonant.json`)

Add (all optional, defaults match existing behavior):

```jsonc
"bridgeType": {
    "enum": ["fullBridge", "halfBridge"],
    "default": "fullBridge"
},
"integratedResonantInductor1": {
    "type": "boolean",
    "default": true,
    "description": "If true, Lr1 is the transformer primary leakage"
},
"integratedResonantInductor2": {
    "type": "boolean",
    "default": true,
    "description": "If true, Lr2 is the transformer secondary leakage"
},
"resonantInductorRatio": {
    "type": "number",
    "default": 1.0,
    "description": "a = n²·L2/L1 (1.0 for symmetric, ~0.95 for asymmetric)"
},
"resonantCapacitorRatio": {
    "type": "number",
    "default": 1.0,
    "description": "b = C2/(n²·C1) (1.0 for symmetric, ~1.052 for asymmetric)"
}
```

Keep `symmetricDesign` as a derived alias (`a==1 && b==1`). Migrate
existing fixtures so `symmetricDesign=false` sets `a=0.95, b=1.052`.

After schema edit, regenerate `MAS.hpp` and copy to
`build/MAS/MAS.hpp` (per the staging gotcha in user memory).

---

## 8. Tests (`tests/TestCllc.cpp`) — DAB-quality bar

Keep the 16 existing tests (smoke-level). **Add per Guide §8:**

| New test | Purpose | NRMSE / tolerance |
|---|---|---|
| `Test_Cllc_Infineon_Reference_Design` | Tighten existing `InfineonExample_Parameters` | ±5 % on n, Ro, C1, L1, Lm |
| `Test_Cllc_KIT_20kW_Asymmetric_Reference_Design` | Asymmetric tank validation | ±10 % |
| `Test_Cllc_Operating_Modes` | sub-resonant / at-fr / super-resonant | `lastMode` & substate sequence match expected |
| `Test_Cllc_ZVS_Boundaries_Forward` | sweep Lm; `lastZvsMarginPrimary` crosses zero at Eq 28 threshold | sign change |
| `Test_Cllc_ZVS_Boundaries_Reverse` | sweep load in reverse mode | sign change |
| `Test_Cllc_BridgeType_HB_vs_FB` | same fr/Q/k → HB has half the Vp swing | exact |
| `Test_Cllc_Reverse_Mode_Power_Balance` | reverse simulation w/ active SR; Vin terminal absorbs power | sign of P_in < 0 |
| `Test_Cllc_SPICE_Netlist` | snubbers, METHOD=GEAR, IC pre-charge, **8** active switches, Evab probe, sense sources | string `find` |
| `Test_Cllc_PtP_AnalyticalVsNgspice_Symmetric_AtResonance` | primary-current NRMSE | **≤ 0.15** |
| `Test_Cllc_PtP_AnalyticalVsNgspice_Symmetric_BelowResonance` | same | **≤ 0.15** |
| `Test_Cllc_PtP_AnalyticalVsNgspice_Symmetric_AboveResonance` | same | **≤ 0.15** |
| `Test_Cllc_PtP_AnalyticalVsNgspice_Asymmetric_AtResonance` | KIT 20 kW fixture | **≤ 0.15** |
| `Test_Cllc_PtP_AnalyticalVsNgspice_Reverse` | reverse direction NRMSE | **≤ 0.15** |
| `Test_Cllc_ExtraComponents_Cr1_Cr2_Lr1_Lr2` | `get_extra_components_inputs(IDEAL/REAL)` returns 4 ancillary `Inputs` (or fewer if integrated) | count + sanity |
| `Test_Cllc_GainCurve_FHA_vs_TDA` | at low Q, FHA and TDA average gains agree | within 5 % |
| `Test_Cllc_LIP_Frequency` | `lipFrequency = ω_P/(2π)` matches design | exact |
| `Test_Cllc_Waveform_Plotting` | CSV dump under `tests/output/` | not asserted |
| `Test_AdvancedCllc_Process` | round-trip user-specified Lr1/Cr1/Lr2/Cr2 | round-trip equal |
| `Test_Cllc_DivideByZero_Guard` | Iout=0 doesn't crash netlist generation | no throw |

NRMSE helpers (`ptp_nrmse`, `ptp_interp`, `ptp_current`,
`ptp_current_time`) extracted from `TestTopologyDab.cpp` to a shared
header `tests/TestTopologyHelpers.h`. (Currently each test file copies
them — DAB / PSFB / PSHB / LLC; this is the right time to share.)

---

## 9. Phased delivery

| # | Scope | Files | Tests added | Effort |
|---|---|---|---|---|
| **P1** | Header rewrite (ASCII art, modes, equations, references, Q convention disambiguation), schema v2, MAS regen, `run_checks` for new fields | `Cllc.h`, schema, MAS regen, fixtures | none (existing keep) | 0.5 d |
| **P2** | TDA solver: `CllcStateVector`, `CllcSubState` enum, `propagate_substate`, `find_next_event`, steady-state Newton, mode detection. Replace `process_operating_point_for_input_voltage` body. | `Cllc.cpp` | `Operating_Modes`, `LIP_Frequency` | 2 d |
| **P3** | Diagnostics population (`lastMode`, `lastZvsMargin*`, `lipFrequency`, `lastSteadyStateResidual`). ZVS condition check Eq 28. | `Cllc.h`, `Cllc.cpp` | `ZVS_Boundaries_Forward/Reverse` | 0.5 d |
| **P4** | SPICE rewrite: active SR bridge, snubbers, GEAR, IC pre-charge, sense sources, Evab probe, divide-by-zero guard, analytical fallback on failure | `Cllc.cpp` | `SPICE_Netlist`, `DivideByZero_Guard` | 1.5 d |
| **P5** | `get_extra_components_inputs` returning Cr1/Cr2/Lr1/Lr2 ancillary `Inputs` (IDEAL + REAL modes), with leakage absorption when `integratedResonantInductor*=true` | `Cllc.cpp` | `ExtraComponents_*` | 0.5 d |
| **P6** | NRMSE acceptance: tighten until ≤ 0.15 for all 5 PtP variants. Failure-discovery phase — expect to fix 2–3 bugs in the TDA solver and SPICE. | tests + minor TDA fixes | `PtP_AnalyticalVsNgspice_*` (5) | 2 d |
| **P7** | Asymmetric tank gain (a, b through FHA + TDA), KIT 20 kW fixture, half-bridge variant, AdvancedCllc rewrite | `Cllc.cpp`, fixtures | `KIT_20kW_Asymmetric`, `BridgeType_HB_vs_FB`, `AdvancedCllc_Process` | 1.5 d |
| **P8** | Reverse mode physical: SR drive, reverse netlist, OP power balance test | `Cllc.cpp`, tests | `Reverse_Mode_Power_Balance`, `PtP_AnalyticalVsNgspice_Reverse` | 1 d |

**Total ~10 working days.**

Dependencies: P2 → P3 → P6; P4 ⊥ P2; P5 ⊥ P4; P7 ⊥ P6; P8 needs P4 and P7.

---

## 10. Acceptance criteria (per Guide §15)

- [ ] `ninja -j2 MKF_tests` builds clean.
- [ ] All `TestCllc` cases pass (16 existing + ~14 new).
- [ ] At least 3 `*_PtP_AnalyticalVsNgspice` tests run and NRMSE ≤ 0.15
      (symmetric at-res, below-res, above-res).
- [ ] `lastZvsMarginPrimary` changes sign at Lm threshold predicted by
      Infineon AN Eq 28.
- [ ] Schema v2 documented; `symmetricDesign` migrated to `a/b` ratios.
- [ ] Header has ASCII art, equations with citations, accuracy
      disclaimer, **Q-convention disambiguation note**, mode table.
- [ ] Netlist parses, contains 8 SW switches (no rectifier diodes),
      snubbers, GEAR, IC pre-charge, Evab probe.
- [ ] `simulate_and_extract_operating_points` runs SPICE for real,
      falls back to analytical with `[analytical]` tag on failure.
- [ ] `simulate_and_extract_topology_waveforms` references only nodes
      the netlist actually emits.
- [ ] `get_extra_components_inputs` supports `IDEAL` and `REAL`.
- [ ] `Test_Cllc_Reverse_Mode_Power_Balance` passes — physical reverse,
      not just label change.

---

## 11. Common CLLC implementation mistakes to avoid

From the research dossier — explicit pitfalls flagged in the literature:

1. **FHA-only ZVS verification** — FHA is necessary, not sufficient.
   Verify Lm-max with the actual `iLm(t_switch)` and Coss charge balance.
2. **Reverse-mode label-only swap** — secondary MOSFETs must be
   re-gated; body-diode-only conduction loses ZVS and overheats.
3. **Mis-defining Q** — three conventions in circulation. Pin one
   (Infineon: `Q = √(Lr/Cr)/Ro`) and document it in the header.
4. **Ignoring magnetizing current in primary current waveform** —
   the triangular Lm component is essential for ZVS prediction.
5. **Symmetric assumption with wide V_battery** — gain envelope
   becomes asymmetric at high SOC; fixed `a=b=1` over-widens fs range.
6. **K = 1.0 for transformer in SPICE** — singular factorization;
   always 0.999..0.9999.
7. **Diode rectifier in "bidirectional" simulation** — physically
   impossible to push power back. Use active SR.
8. **Dead-time longer than half-cycle slack** at high fs.
9. **k too high (Lm too low) "to guarantee ZVS"** — increases
   circulating current and conduction loss; sweet spot 3.5–6.
10. **No clamp on capacitive-region operation** — controller must
    hard-limit fs above the resonance point at light load.

---

## 12. Checklist for the implementing agent

Before each phase:

- [ ] Re-read `CONVERTER_MODELS_GOLDEN_GUIDE.md`.
- [ ] Verify the latest `Llc.cpp` TDA solver pattern (the closest
      sibling to copy from).
- [ ] Verify the latest `Dab.cpp` SPICE pattern (gold-standard
      netlist).
- [ ] Build with `ninja -j2` in the build dir, not `cmake --build`
      (per user memory).
- [ ] After MAS schema edits, copy `MAS.hpp` from source to
      `build/MAS/MAS.hpp` (per user memory).
- [ ] Mark each phase as complete only when its NRMSE / sign-change /
      build-clean criteria are met.

After P8 (final):

- [ ] Delete `process_operating_points_for_input_voltage` (misnamed;
      replaced by `process_operating_point_for_input_voltage` in P2).
- [ ] Update CLAUDE.md / user-memory entry on converter quality tiers
      to mark CLLC as DAB-quality.
- [ ] Open a follow-up issue for three-phase / interleaved CLLC if
      ever needed (out of scope here).
