# Asymmetric Half-Bridge (AHB) — Topology Add Plan

> Purpose: add a new `AsymmetricHalfBridge` topology to
> `src/converter_models/` (none exists today), wired into the MAS
> schema, with full ngspice harness and tests, meeting **DAB-quality**
> as defined by `CONVERTER_MODELS_GOLDEN_GUIDE.md`.
>
> Read `CONVERTER_MODELS_GOLDEN_GUIDE.md` first; this plan layers on top.
> Sibling references:
> - `PhaseShiftedHalfBridge.{h,cpp}` — closest two-switch primary; **header docstring already documents the AHB-vs-Pinheiro-Barbi disambiguation** (lines 21-42 of the .h file). The MKF `Pshb` class is the Pinheiro-Barbi 3-level half-bridge — AHB is a *separate* topology, not a rewrite.
> - `PhaseShiftedFullBridge.{h,cpp}` — closest **rectifier-type-aware** pattern (CT / CD / FB)
> - `ActiveClampForward.{h,cpp}` — closest **DC-blocking-cap + asymmetric-magnetizing-ramp** pattern
> - `SingleSwitchForward.{h,cpp}` / `TwoSwitchForward.{h,cpp}` — closest **forward-class secondary** pattern (Lo + rectifier + Co)
> - `Dab.{h,cpp}` — gold standard for diagnostics, NRMSE tests, SPICE quality
> - `PwmBridgeSolver.h` — primary-leg PWM scaffold to reuse
> - `CONVERTER_MODELS_GOLDEN_GUIDE.md`, `CLLC_REWRITE_PLAN.md`, `CUK_PLAN.md`, `WEINBERG_PLAN.md`, `FOUR_SWITCH_BUCK_BOOST_PLAN.md`
>
> AHB is a non-resonant, isolated, single-leg PWM converter with a
> DC-blocking cap and **asymmetric** D/(1-D) complementary drive.
> Solver is analytical piecewise-linear (no Newton, no matrix
> exponential). The distinctive challenges for MKF are: (a) the
> **asymmetric primary voltage** (+(1-D)·Vin / -D·Vin); (b) the
> **non-monotonic gain curve** Vo = 2·D·(1-D)·Vin/n peaking at D=0.5;
> (c) **transient DC-bias risk** in the transformer when V_Cb lags
> Vin during a transient; (d) the existing PSHB header's
> disambiguation, which the AHB header must mirror so future readers
> never conflate them.

---

## 1. Why AHB in MKF

- **Existing PSHB header already promises this disambiguation** — the
  `PhaseShiftedHalfBridge.h` lines 21-42 explicitly call out that the
  MKF `Pshb` class is Pinheiro-Barbi 3-level, *not* AHB. Adding the
  proper AHB class closes that loop.
- **Sweet-spot topology for 30-500 W AC/DC** — server 5/12/24 V
  bricks, telecom 48 V, 30-90 W AC adapters, single-stage LED
  drivers. ON NCP1605, ST L6591, TI UCC28251, Power Integrations
  LinkSwitch use this topology in commodity volume.
- **Fills the gap between Forward and PSHB.** Single-switch forward
  (uses ½ B-H curve) + two-switch forward (1 ½ B-H + RCD or tertiary
  reset) → AHB (full B-H, automatic flux balance, ZVS for free, only
  2 switches). The next step up is PSFB (4 switches, 1-5 kW). AHB
  is the right answer at 100-500 W.
- **Exercises a novel magnetic case for MKF**: transient DC bias
  when V_Cb lags Vin during input transients. Adviser must check
  *transient* flux excursion, not just steady state — a check no
  other topology in MKF currently demands.

---

## 2. Variants in scope

From the research dossier (Imbertson-Mohan IEEE TIA 1993; Pressman
Ch. 16; Erickson-Maksimović 2e §6.3; TI SLUP223; ON AN-4153; ST
AN2852; TI SLUA121; TI UCC28251; Renesas ISL6745):

1. **Classic AHB (V1)** — single leg (Q1 HS, Q2 LS), DC-blocking
   cap Cb in series with primary, **center-tapped (CT)** secondary
   + Lo + Co. Sweet spot 100-300 W telecom 48 V output.
   **Baseline.**
2. **AHB with full-bridge (FB) secondary (V2)** — same primary,
   secondary uses 4-diode FB rectifier with single secondary winding.
   Sweet spot 30-100 W AC adapters where copper utilization matters.
3. **AHB with current-doubler (CD) secondary (V3)** — single
   secondary winding feeding two output inductors L1 + L2 sharing
   Iout/2. Halves secondary copper loss. Preferred for 5-12 V /
   ≥20 A server outputs. TI SLUA121.
4. **Synchronous-rectifier AHB (V4)** — diodes → SR FETs (any of
   V1/V2/V3). Self-driven SR variant uses the asymmetric primary
   voltage to drive SRs directly.
5. **AHB-Flyback (V5)** — Lo collapsed to zero; Lm itself is the
   energy-storage element. Single-stage 30-90 W AC/DC adapters.
   ST L6591 reference.
6. **Multi-output AHB (V6)** — multiple secondaries on the same
   transformer (cross-regulation worse than forward; warn-once like
   DAB §4.7).

**Out of scope for the initial add** (defer to follow-ups):

- Two-transformer AHB (NTU 2013) — splits flux excursion across two
  cores; high-power variant.
- AHB-LLC hybrid (datacenter PSU, 2022+ IEEE TPEL) — adds a small
  Cr to enable DCX-mode; really a different topology.
- GaN-AHB at 1-2 MHz — no topology change, just scaling.
- DCM design (V1 enforces strict CCM in `process_design_requirements`;
  DCM detected and flagged but not designed for — same approach as
  Cuk, Weinberg, 4SBB plans).

---

## 3. Reference designs as test fixtures

| # | Source | Vin | Vout | Pout | fsw | n | Lo | Lm | Llk | Cb | Use |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **TI SLUP223 worked example** (Steigerwald) | 100 V | 5 V | 100 W (20 A) | 200 kHz | 10:1 (CT) | 1 µH | 50 µH | 1 µH | 0.47 µF | **primary fixture** — full numerics, CT secondary |
| 2 | **ON AN-4153** (FPS-based) | 90-264 Vac (rectified ~125-375 V) | 12 V | 60 W | 100 kHz | depends on D_max | – | – | – | – | low-power AC adapter, FB rectifier |
| 3 | **ST AN2852 / EVL6591-90WADP** | 90-264 Vac | 19 V | 90 W | 100 kHz | – | – | – | – | – | mid-power AC adapter, full BoM published |
| 4 | **Renesas ISL6745EVAL1Z** | 36-72 V | 12 V | 120 W (10 A) | 250 kHz | – | – | – | – | – | telecom 48 V brick |
| 5 | **Imbertson-Mohan IEEE TIA 1993** | 200 V | 50 V | 250 W | 100 kHz | 4:1 | – | – | – | – | foundational reference, original ZVS demo |
| 6 | **TI SLUA121 CD example** | 48 V | 3.3 V | 60 W (18 A) | 200 kHz | 8:1 | 2 × 1.5 µH | – | – | – | V3 current-doubler reference |
| 7 | **MDPI Energies 11(7):1721 (2018)** AHB-Flyback | 90-264 Vac | 24 V | 65 W | 65 kHz | – | (no Lo) | – | – | – | V5 flyback variant reference |
| 8 | **Power Integrations LinkSwitch app circuits** | 90-264 Vac | 5/12/19 V | 30-150 W | 66-132 kHz | – | – | – | – | – | volume AC-adapter fixtures |

URLs (for citation in header docstring):

- Imbertson & Mohan IEEE TIA 29(1) 1993: https://doi.org/10.1109/28.195897
- ON Semi AN-4153: https://www.onsemi.com/pub/collateral/an-4153.pdf
- ST AN2852 / EVL6591-90WADP: https://www.st.com/resource/en/application_note/an2852-evl659190wadp-90-w-acdc-asymmetrical-halfbridge-adapter-using-l6591-and-l6563-stmicroelectronics.pdf
- TI SLUA121 Current-Doubler Rectifier: https://www.ti.com/lit/an/slua121/slua121.pdf
- TI UCC28251 in AHB mode (E2E): https://e2e.ti.com/support/power-management-group/power-management/f/power-management-forum/1104129/ucc28251-asymmetric-half-bridge-operation-for-off-line-converters-400vin-to-isolated-24vout
- Renesas ISL6745EVAL1Z: https://www.renesas.com/en/document/gde/isl6745eval1z-user-guide
- AHB ZVS-PWM springer chapter: https://link.springer.com/chapter/10.1007/978-3-319-96178-1_10
- High-Efficiency AHB IEEE TPEL 2022: https://ieeexplore.ieee.org/document/9796986/
- DTU IECON 2022 AHB datacenter: https://orbit.dtu.dk/files/262634420/IECON_AHB_datacenter.pdf
- Two-Transformer AHB (NTU 2013): https://scholars.lib.ntu.edu.tw/bitstreams/99078c69-7ddf-40a5-971b-c877a9a8351c/download
- AHB-Flyback MDPI Energies: https://www.mdpi.com/1996-1073/11/7/1721
- Resonant AHB Flyback MDPI Apps: https://www.mdpi.com/2076-3417/12/13/6685
- AHB modeling and simulation: https://www.researchgate.net/publication/334329298_Modeling_and_Simulation_of_the_Asymmetrical_Half_Bridge_Converter

---

## 4. Equations (canonical block for the new header docstring)

Notation: `Vin` = input bus voltage; `Vo` = output voltage; `n = Np/Ns`
(half winding for CT, full winding for FB / CD); `D` = duty of Q1
(complementary 1-D for Q2); `Tsw = 1/fsw`; `td` = dead-time;
`R = Vo/Iout`; `Vd` = rectifier voltage drop; `Llk` = primary leakage;
`Lm` = primary magnetizing inductance; `Coss` = MOSFET output
capacitance (per FET).

### 4.1 Primary voltage and DC-blocking cap

```
(E1)  V_Cb  = (1−D)·Vin                  (steady-state DC voltage on Cb — Imbertson-Mohan eq. 3)
(E2)  v_pri = +(1−D)·Vin   during Q1 ON  (D·Tsw)
(E3)  v_pri = −D·Vin       during Q2 ON  ((1−D)·Tsw)
(E4)  Volt-second balance: (1−D)·Vin·D·Tsw + (−D·Vin)·(1−D)·Tsw = 0   ✓ AUTOMATIC
       (the Cb voltage tracks D continuously; controller does not
        need active flux-balance compensation)
```

### 4.2 Conversion ratio

```
(E5)  Vo (CT or FB rectifier, CCM)  = 2·D·(1−D)·Vin / n     ← non-monotonic in D
(E6)  Vo (CD rectifier, CCM)        = D·(1−D)·Vin / n        (factor of 2 absorbed by doubler)
(E7)  Maximum gain at D = 0.5       : Vo,max = Vin / (2·n)
```

The gain curve is **non-monotonic** — same Vo achievable at two
different D values (one above, one below 0.5). Controllers operate
on **one branch only**, conventionally `D ∈ [0, 0.5]` (so that
increasing D increases Vo).

### 4.3 Sizing (CCM)

```
(S1)  n        = 2·D_max·(1−D_max)·Vin_min / (Vo + Vd)        (target D_max ≤ 0.5 at Vin_min)
(S2)  Lo       ≥ Vo · (1 − 2·D·(1−D)) / (ΔILo · fsw)
                  (at D = 0.5: Lo ≥ 0.5·Vo / (ΔILo·fsw))
(S3)  Cb       ≥ I_pri,pk · D / (fsw · ΔV_Cb)                  (ΔV_Cb ≤ 5 % V_Cb — ON AN-4153 eq. 16)
(S4)  Lm       ≤ ((1−D)·Vin·D·Tsw) / (2·Im,target)             (sized to provide ZVS assist current)
(S5)  Co       ≥ ΔILo / (8·fsw·ΔVo)
(S6)  ZVS:    0.5·(Llk + Lm,refl)·I_pri,min²  ≥  2·Coss·Vin²   (Imbertson-Mohan; same form as PSFB)
(S7)  Dead-time: td  ≥  π·√(Llk · 2·Coss)                      (resonant ZVS transition)
```

### 4.4 Switch and rectifier stresses

```
(E8)  V_Q1,pk = V_Q2,pk = Vin            (each blocks the full bus when the other conducts)
(E9)  I_Q1,rms = I_pri,Q1 · √D
(E10) I_Q2,rms = I_pri,Q2 · √(1−D)        (UNEQUAL when D ≠ 0.5 — different switch RMS)
(E11) I_pri,rms ≈ (Iout/n) · √(D² + (1−D)²)    (worst at D = 0.5: 1/√2)
(E12) V_D,pk (CT/FB) = Vin/n
(E13) V_D,pk (CD)    = Vin/n
```

### 4.5 Asymmetric flux excursion (the magnetic-design crux)

```
(F1)  λ⁺ = (1−D)·Vin · D·Tsw         (positive volt-seconds, Q1 ON)
(F2)  λ⁻ = D·Vin · (1−D)·Tsw         = λ⁺   (steady-state symmetry)
(F3)  ΔΦ = λ⁺ / Np = D·(1−D)·Vin·Tsw / Np
(F4)  Bpk = ΔΦ / (2·Ae) = D·(1−D)·Vin·Tsw / (2·Np·Ae)
```

Steady-state ΔΦ scales as `D·(1−D)` → peaks at D = 0.5. Core
utilization at D = 0.5 ≈ 70-75 % of B-H curve, degrading to ≈ 30 %
at D = 0.1 (vs PSFB ~85 % at any D). This is the price of single-leg
simplicity.

### 4.6 Transient DC-bias risk

When `Vin` steps (load transient, brown-out recovery), `V_Cb` cannot
change instantaneously — its time-constant is `R_dcr · L_m`. During
the transient, the volt-second balance (E4) is temporarily violated
and a transient DC magnetizing current appears. **Peak transient B
can exceed steady-state Bpk (E4)** and saturate the core even if
steady-state design is fine. The MKF magnetic adviser must check
**transient** flux excursion, not just steady state — this is the
load-bearing novel constraint AHB introduces.

### 4.7 ZVS condition (per Imbertson-Mohan 1993)

```
(Z1)  E_stored,Llk + E_stored,Lm,refl  ≥  E_charge,Coss
       0.5·(Llk + Lm,refl)·I_pri(t_sw)²  ≥  2·Coss·Vin²
(Z2)  Dead-time matched to resonant period:  td ≈ π·√(Llk · 2·Coss)
```

ZVS extension via Lm sizing (S4): make Im,target large enough that
even at minimum load the assist current keeps I_pri above the ZVS
threshold.

---

## 5. Solver: piecewise-linear analytical (no TDA)

AHB is a PWM converter; sub-intervals are LTI, primary current is
exact-linear within each. Pattern to follow: combine
`SingleSwitchForward.cpp` (forward-class secondary with rectifier +
Lo) and `ActiveClampForward.cpp` (DC-blocking cap + asymmetric
magnetizing ramp). Per Guide §4.

### 5.1 Sub-intervals (CCM, V1 classic, CT secondary)

Per switching period, the canonical sequence is:

| Sub-int | Length | Q1 | Q2 | v_pri | i_Lm slope | rectifier conducting | V_Lo |
|---|---|---|---|---|---|---|---|
| **A** Q1 ON | D·Tsw      | ON  | OFF | +(1-D)·Vin | +(1-D)·Vin/Lm | top diode (Ds_a) | (1-D)·Vin/n − Vo |
| **B** dead-time | td   | OFF | OFF | resonant Llk-Coss ringing | non-linear briefly | one diode reaches OFF | freewheel via Lo via both diodes |
| **C** Q2 ON | (1-D)·Tsw  | OFF | ON  | −D·Vin     | −D·Vin/Lm     | bottom diode (Ds_b) | D·Vin/n − Vo |
| **D** dead-time | td   | OFF | OFF | mirror of B | mirror | mirror | freewheel |

For analytical waveforms ignore B/D's resonant transition (it's a
small fraction of Tsw and contributes negligibly to RMS / peak); the
SPICE netlist captures it correctly.

`Im` propagates linearly through A and C with **asymmetric slopes**
(reuse `ActiveClampForward.cpp` magnetizing-ramp computation).

`I_Lo` (output inductor current) sees three slope regions per cycle
(rectifier-on, freewheel, rectifier-on, freewheel) — same shape as
forward.

Time grid: `2·N_samples + 1` with `N_samples = 256`, span `[0, Tsw]`
exactly. Per Guide §4.8.

### 5.2 Rectifier-type-aware secondaries

Mirror PSFB / PSHB §4.5:

- **CT (V1)**: two secondary windings (`Secondary 0a`, `Secondary 0b`).
  Each carries unipolar pulses on alternate half-cycles. Reverse
  voltage = 2·Vsec_pk.
- **CD (V3)**: one secondary winding, bipolar 3-level voltage; two
  output inductors L1, L2 emitted via `extraLoCurrentWaveforms` and
  `extraLo2CurrentWaveforms`.
- **FB (V2)**: one secondary winding, bipolar 3-level voltage; one
  output inductor.

### 5.3 Per-OP diagnostics (mandatory per Guide §2.4 / §4.6)

```cpp
mutable double lastDutyCycle = 0;
mutable double lastConversionRatio = 0;
mutable double lastDcBlockingCapVoltage = 0;     // V_Cb = (1-D)·Vin
mutable double lastDcBlockingCapRipple = 0;      // ΔV_Cb
mutable double lastPrimaryPeakVoltagePositive = 0;   // (1-D)·Vin
mutable double lastPrimaryPeakVoltageNegative = 0;   // −D·Vin
mutable double lastSwitchPeakVoltageQ1 = 0;       // = Vin
mutable double lastSwitchPeakVoltageQ2 = 0;       // = Vin
mutable double lastSwitchRmsCurrentQ1 = 0;        // asymmetric: ≠ Q2 if D ≠ 0.5
mutable double lastSwitchRmsCurrentQ2 = 0;
mutable double lastZvsMargin = 0;                 // (Llk + Lm,refl)·I_pri²/2 − 2·Coss·Vin²
mutable double lastResonantTransitionTime = 0;    // π·√(Llk·2·Coss)
mutable double lastSteadyStateFluxExcursion = 0;  // ΔΦ steady state
mutable double lastTransientFluxExcursionEstimate = 0;  // worst-case during Vin step
mutable double lastMagnetizingCurrentRipple = 0;
mutable double lastOutputInductorRipple = 0;
mutable int    lastOperatingMode = 0;             // 0 = CCM, 1 = DCM
mutable int    lastRectifierType = 0;             // 0=CT, 1=CD, 2=FB
```

### 5.4 Multi-output (V6) note

Cross-regulation in AHB is **worse than forward** because the
secondary winding voltage is asymmetric — different positive vs
negative half-cycle amplitude. Print the warn-once message at the
start of `process_operating_points` when secondaries > 1, similar
to `Dab.cpp:855-906` pattern. Document explicitly in the header
docstring that multi-output AHB is approximate unless per-secondary
leakage is provided.

### 5.5 DC-bias transient check

When user supplies `inputVoltageStepRange` schema field (optional),
compute the worst-case transient `B` excursion under a Vin step:

```
(T1)  ΔV_Cb,transient ≈ ΔVin · (1 − exp(−τ/(R_dcr·Lm)))     (settles after τ)
(T2)  Bpk,transient   ≈ Bpk,ss + (D·ΔVin·Tsw) / (2·Np·Ae)    (worst-case adder)
```

Surface as `lastTransientFluxExcursionEstimate` so the magnetic
adviser sizes the core for transient peak, not steady-state peak.

---

## 6. SPICE netlist (generate_ngspice_circuit)

AHB has 2 active switches, a DC-blocking cap, a transformer, and a
rectifier-type-dependent secondary. The netlist follows Guide §5; the
DC-blocking cap **MUST be IC pre-charged to V_Cb = (1-D)·Vin** —
without this, ngspice will saturate the transformer in the first
cycle.

### 6.1 Topology — V1 classic AHB with CT secondary

```
* Asymmetric Half-Bridge (V1 classic, CT secondary) — Generated by OpenMagnetics
*
*    Vin ──┬── Q1 ──┬── A ── Cb ── Np ── B ─────  GND_pri
*          │        │                │
*          │       SW1               │
*          │        │                │
*          ├── Q2 ──┘                Llk (lumped on primary side)
*          │        │
*         Cin       │
*          │       GND_pri
*         GND_pri
*
*    Secondary (CT FW rectifier):
*         Ns_a ─── Ds_a ───┬── Lo ── out_node ── Co ── 0
*                          │
*         Ns_b ─── Ds_b ───┤
*                          │
*                         Rload ── 0
```

### 6.2 Mandatory netlist contents (Guide §5 + AHB specifics)

- `Vin vin_p 0 {Vin}`; zero-volt `Vin_sense` for input current.
- 2 active switches:
  ```
  Q1 vin_p mid_A pwm1 0 SW1
  Q2 mid_A 0     pwm2 0 SW1
  ```
- **DC-blocking cap with mandatory IC**:
  ```
  Vcb_sense mid_A cb_in 0
  C_b cb_in pri_a {Cb} ic={(1-D)*Vin}
  ```
- Primary winding with leakage (use lumped Llk on primary side):
  ```
  Llk pri_a pri_a_int {Llk}
  Lpri pri_a_int 0 {Lm}
  ```
  (For coupling, use `K_xfmr Lpri Lsec_a Lsec_b 0.99` — not 1.0.)
- CT secondary (V1):
  ```
  Lsec_a sec_ct sec_a {Lsec_half}
  Lsec_b sec_ct sec_b {Lsec_half}
  K_xfmr Lpri Lsec_a Lsec_b 0.99
  Vsec1_sense_o1 sec_ct sec_ct_sense 0
  Ds_a sec_a node_lo DIDEAL
  Ds_b sec_b node_lo DIDEAL
  Lo node_lo lo_out {Lo} ic={Iout}
  Vout_sense_o1 lo_out out_node_o1 0
  Co out_node_o1 0 {Co} ic={Vo}
  Rload out_node_o1 0 {max(Vo/Iout, 1e-3)}
  ```
- FB secondary (V2): single secondary winding, 4 diodes.
- CD secondary (V3): single secondary winding, 2 diodes, two output
  inductors L1, L2. Mirror PSFB / PSHB CD wiring exactly.
- **Snubbers** across each switch (mandatory for ZVS dead-time
  ringing):
  ```
  Rsnub_Q1 vin_p snub_q1 1k
  Csnub_Q1 snub_q1 mid_A 1n
  ; same for Q2
  ```
- Bridge-voltage probe (Guide §5):
  `Evab vab 0 VALUE={V(mid_A)}`.
- Models (verbatim from `Boost.cpp`):
  ```
  .model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg
  .model DIDEAL D(IS=1e-12 RS=0.05)
  ```
- PWM sources (complementary with dead-time `td`):
  ```
  Vpwm1 pwm1 0 PULSE(0 5 0          {td} {td} {D*Tsw - td} {Tsw})
  Vpwm2 pwm2 0 PULSE(0 5 {D*Tsw}    {td} {td} {(1-D)*Tsw - td} {Tsw})
  ```
- Solver options (mandatory; **GEAR is non-negotiable for ZVS
  resonant transitions**):
  ```
  .options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500
  .options METHOD=GEAR TRTOL=7
  ```
- Time step: `period/200` baseline; `period/500` if `td < period/40`
  (i.e. dead-time is short relative to switching period; need finer
  resolution to capture ZVS transition).
- IC pre-charge (**all four mandatory**):
  `.ic v(out_node_o1)={Vo} v(cb_in)={(1-D)*Vin} v(pri_a)=0 i(Lo)={Iout} i(Lpri)=0`
- `.tran {stepTime} {totalTime} {steadyStateStart} UIC`.
- Saved signals:
  ```
  .save v(vab) v(mid_A) v(pri_a) v(cb_in) v(out_node_o1)
        i(Vin_sense) i(Vcb_sense) i(Vsec1_sense_o1) i(Vout_sense_o1)
        i(Lpri) i(Lo)
  ```
- Comment header: topology, rectifierType, Vin/Vo/Pout/fsw/D/n/Lo/Lm/Cb,
  "generated by OpenMagnetics".
- Parasitic R: bake DCR ≈ 50 mΩ into Lo and Lpri, ESR ≈ 5 mΩ into Co
  and Cb — same rationale as Cuk/Weinberg/4SBB plans.

### 6.3 simulate_and_extract_operating_points

Follow Guide §6 verbatim:

1. Construct `NgspiceRunner`.
2. If `runner.is_available()` is false, warn and fall back to analytical.
3. For each (Vin, OP), build netlist for the configured `rectifierType`,
   run, on `success==false` log error, fall back to
   `process_operating_point_for_input_voltage`, tag name with `[analytical]`.
4. On success, build `WaveformNameMapping` and call
   `NgspiceRunner::extract_operating_point`.

Mapping: primary winding via `Vin_sense` and `mid_A`. Secondary
windings depend on rectifier type — emit per-winding excitations
matching the netlist nodes.

### 6.4 simulate_and_extract_topology_waveforms

```cpp
wf.set_input_voltage(getWaveform("vin_sense"));
wf.set_input_current(getWaveform("vin_sense#branch"));
wf.get_mutable_output_voltages().push_back(getWaveform("out_node_o1"));
wf.get_mutable_output_currents().push_back(getWaveform("vout_sense_o1#branch"));
```

---

## 7. Schema additions (`MAS/schemas/inputs/topologies/asymmetricHalfBridge.json`)

New file. Pattern off `phaseShiftFullBridge.json` (rectifier-type-aware)
+ `flyback.json` (turns ratio + isolation):

```jsonc
{
    "title": "asymmetricHalfBridge",
    "description": "Asymmetric Half-Bridge converter (Imbertson-Mohan 1993, single-leg + DC-blocking cap)",
    "type": "object",
    "properties": {
        "inputVoltage":         { "$ref": "../../utils.json#/$defs/dimensionWithTolerance" },
        "diodeVoltageDrop":     { "type": "number", "default": 0.6 },
        "maximumSwitchCurrent": { "type": "number" },
        "currentRippleRatio":   { "type": "number" },
        "efficiency":           { "type": "number", "default": 0.95 },
        "useLeakageInductance": { "type": "boolean", "default": true },

        "rectifierType": {
            "enum": ["centerTapped", "currentDoubler", "fullBridge"],
            "default": "centerTapped"
        },
        "synchronousRectifier": { "type": "boolean", "default": false },
        "variant": {
            "enum": ["classic", "ahbFlyback"],
            "default": "classic"
        },

        "outputInductance":          { "type": "number" },
        "outputCapacitance":         { "type": "number" },
        "dcBlockingCapacitance":     { "type": "number" },
        "couplingCoefficient":       { "type": "number", "default": 0.99 },
        "deadTime":                  { "type": "number", "default": 100e-9 },

        "inputVoltageStepRange": {
            "$ref": "../../utils.json#/$defs/dimensionWithTolerance",
            "description": "Optional Vin step magnitude for transient flux-excursion check"
        },

        "operatingPoints": {
            "type": "array",
            "items": { "$ref": "#/$defs/asymmetricHalfBridgeOperatingPoint" },
            "minItems": 1
        }
    },
    "required": ["inputVoltage", "diodeVoltageDrop", "operatingPoints"],
    "$defs": {
        "asymmetricHalfBridgeOperatingPoint": {
            "title": "asymmetricHalfBridgeOperatingPoint",
            "$ref": "../../utils.json#/$defs/baseOperatingPoint"
        }
    }
}
```

After schema edit, regenerate `MAS.hpp` and copy to `build/MAS/MAS.hpp`
(per user-memory staging gotcha).

Add `"ASYMMETRIC_HALF_BRIDGE_CONVERTER"` to the `Topologies` enum
(alongside PHASE_SHIFTED_HALF_BRIDGE_CONVERTER, etc.).

---

## 8. Class structure (`AsymmetricHalfBridge.h` / `AsymmetricHalfBridge.cpp`)

Mirror the `PhaseShiftedFullBridge` pattern (closest rectifier-type-aware
sibling) with magnetizing-ramp logic from `ActiveClampForward`:

```cpp
namespace OpenMagnetics {
using namespace MAS;

class AsymmetricHalfBridge : public MAS::AsymmetricHalfBridge, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 20;

    // 2.2 Computed design values
    double computedTurnsRatio                   = 0;
    double computedOutputInductance             = 0;   // Lo
    double computedMagnetizingInductance        = 0;   // Lm
    double computedLeakageInductance            = 1e-6;// Llk
    double computedDcBlockingCapacitance        = 0;   // Cb
    double computedOutputCapacitance            = 0;   // Co
    double computedDutyCycle                    = 0;
    double computedDeadTime                     = 100e-9;
    double computedDiodeVoltageDrop             = 0.6;

    // 2.3 Schema-less device parameters
    double mosfetCoss = 200e-12;

    // 2.4 Per-OP diagnostics — see §5.3

    // 2.5 Extra-component waveforms (cleared in process_operating_points)
    mutable std::vector<Waveform> extraLoVoltageWaveforms;
    mutable std::vector<Waveform> extraLoCurrentWaveforms;
    mutable std::vector<Waveform> extraLo2VoltageWaveforms;    // CD only
    mutable std::vector<Waveform> extraLo2CurrentWaveforms;    // CD only
    mutable std::vector<Waveform> extraCbVoltageWaveforms;
    mutable std::vector<Waveform> extraCbCurrentWaveforms;

public:
    bool _assertErrors = false;

    AsymmetricHalfBridge(const json& j);
    AsymmetricHalfBridge() {};

    // 2.7-2.10 boilerplate accessors per Guide §2 (copy PhaseShiftFullBridge.h shape)

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
        const AsymmetricHalfBridgeOperatingPoint& opPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // 2.13 Static analytical helpers
    static double compute_duty_cycle(double Vin, double Vo, double n,
                                     double diodeDrop, double efficiency);
    static double compute_dc_blocking_cap_voltage(double Vin, double D);
    static double compute_conversion_ratio(double D, double n);          // 2D(1-D)/n
    static double compute_turns_ratio(double Vin_min, double Vo, double D_max,
                                       double diodeDrop);
    static double compute_lo_min(double Vo, double D, double dILo, double fsw);
    static double compute_lm_min_for_zvs(double Vin, double D, double Tsw,
                                          double Im_target);
    static double compute_cb_min(double Iprimary_pk, double D, double dVCb,
                                  double fsw);
    static double compute_zvs_energy_balance(double Llk, double Lm_refl,
                                              double Ipri, double Coss, double Vin);
    static double compute_dead_time(double Llk, double Coss);
    static double compute_steady_state_flux_excursion(double Vin, double D,
                                                       double Tsw, double Np, double Ae);
    static double compute_transient_flux_excursion(double Vin, double dVin, double D,
                                                    double Tsw, double Np, double Ae,
                                                    double Lm, double Rdcr);

    // 2.14 Extra-components factory (returns one Inputs per ancillary magnetic)
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

class AdvancedAsymmetricHalfBridge : public AsymmetricHalfBridge {
    /* desiredTurnsRatios, desiredMagnetizingInductance, desiredLeakageInductance,
       desiredOutputInductance, desiredDcBlockingCapacitance, desiredOutputCapacitance */
    /* Inputs process(); */
};

void from_json(const json& j, AdvancedAsymmetricHalfBridge& x);
void to_json(json& j, const AdvancedAsymmetricHalfBridge& x);
}
```

The "primary" magnetic returned by `process_design_requirements` is
the **2-winding transformer** (or 3-winding for CT). Lo (and Lo2 for
CD), Cb, Co flow through `get_extra_components_inputs`.

---

## 9. Tests (`tests/TestAsymmetricHalfBridge.cpp`) — DAB-quality bar

Mirror DAB and PSFB / PSHB test sets. Required `TEST_CASE`s per
Guide §8:

| Test | Purpose | Tolerance |
|---|---|---|
| `Test_AHB_SLUP223_Reference_Design` | reproduce TI SLUP223 fixture (100 V → 5 V/20 A, 200 kHz, 10:1 CT) | ±5 % on D, M, V_Cb, ZVS margin |
| `Test_AHB_AN4153_Reference_Design` | ON Semi AN-4153 60 W AC adapter | ±10 % |
| `Test_AHB_ST_AN2852_Reference_Design` | ST EVL6591 90 W AHB | ±10 % |
| `Test_AHB_Renesas_ISL6745_Reference_Design` | telecom 36-72 V → 12 V/120 W | ±5 % |
| `Test_AHB_Imbertson_Mohan_Original` | reproduce IEEE TIA 1993 250 W example | ±10 % |
| `Test_AHB_Design_Requirements` | turns ratio, Lo, Lm, Cb, Co positive and within sane ranges; round-trip a power target | algebraic |
| `Test_AHB_OperatingPoints_Generation` | multiple Vin (min/nom/max), volt-sec balance: λ⁺ + λ⁻ = 0; primary current piecewise linear | algebraic |
| `Test_AHB_NonMonotonic_Gain_Curve` | sweep D from 0.1 to 0.9; verify M = 2D(1-D)/n peaks at D=0.5 | exact |
| `Test_AHB_DC_Blocking_Cap_Voltage` | `lastDcBlockingCapVoltage = (1-D)·Vin` for D ∈ {0.2, 0.4, 0.6} | ±2 % |
| `Test_AHB_Asymmetric_Switch_RMS` | I_Q1,rms ≠ I_Q2,rms for D ≠ 0.5; equal at D=0.5 | ratio |
| `Test_AHB_ZVS_Boundaries` | sweep Llk; `lastZvsMargin` crosses zero at energy-balance threshold | sign change |
| `Test_AHB_Steady_State_Flux_Excursion` | ΔΦ = D·(1-D)·Vin·Tsw / Np matches (E2)/(F3) | exact |
| `Test_AHB_Transient_Flux_Excursion` | with `inputVoltageStepRange` set, transient B exceeds steady-state B | inequality |
| `Test_AHB_RectifierType_CT` | V1; correct CT secondary waveforms (2 windings, alternating) | shape |
| `Test_AHB_RectifierType_CD` | V3; 1 secondary, 2 output inductors L1+L2 each carrying Iout/2 | algebraic |
| `Test_AHB_RectifierType_FB` | V2; 1 secondary, 4 diodes, single Lo | shape |
| `Test_AHB_DCM_Detection` | sweep load; `lastOperatingMode` flips at boundary | sign change |
| `Test_AHB_Synchronous_Variant` | V4; netlist contains 2 SR FETs (no diodes); simulation runs | smoke |
| `Test_AHB_Flyback_Variant` | V5 ahbFlyback (Lo = 0); Lm energy-storage check | algebraic |
| `Test_AHB_Multiple_Outputs_Warn_Once` | 2 outputs → warn-once message printed; cross-regulation flagged | string find |
| `Test_AHB_SPICE_Netlist` | netlist parses, contains 2 SW switches, snubbers, METHOD=GEAR, **IC pre-charge for Cb at (1-D)·Vin**, IC for Lo, Co; rectifier-type-aware structure | string find |
| `Test_AHB_DivideByZero_Guard` | Iout=0 → no crash; Rload clamps to 1e-3 | no throw |
| `Test_AHB_PtP_AnalyticalVsNgspice_AtD05` | SLUP223 fixture at D=0.5; primary-current NRMSE | **≤ 0.15** |
| `Test_AHB_PtP_AnalyticalVsNgspice_LightLoad` | D≈0.3; NRMSE | **≤ 0.15** |
| `Test_AHB_PtP_AnalyticalVsNgspice_HighDuty` | D≈0.45; NRMSE | **≤ 0.15** |
| `Test_AHB_PsHB_Disambiguation` | a fixture *labelled* PSHB (Pinheiro-Barbi) should NOT parse as AHB; verify schema rejects PSHB-only fields | throws |
| `Test_AHB_ExtraComponents_Lo_Lo2_Cb_Co` | `get_extra_components_inputs(IDEAL/REAL)` returns Lo, Cb, Co (V1) or Lo, Lo2, Cb, Co (V3) | count + sanity |
| `Test_AHB_Waveform_Plotting` | CSV dump under `tests/output/`; visual inspection | not asserted |
| `Test_AdvancedAHB_Process` | round-trip user-specified turnsRatios/Lm/Llk/Lo/Cb/Co | round-trip equal |

NRMSE helpers (`ptp_nrmse`, `ptp_interp`, `ptp_current`,
`ptp_current_time`) extracted to a shared header
`tests/TestTopologyHelpers.h` (also a CLLC, Cuk, Weinberg, 4SBB plan
deliverable — do it once, all five plans benefit).

---

## 10. Phased delivery

| # | Scope | Files | Tests added | Effort |
|---|---|---|---|---|
| **P1** | MAS schema (`asymmetricHalfBridge.json`, `Topologies` enum entry), regenerate MAS.hpp, stub `AsymmetricHalfBridge.{h,cpp}` (constructors + run_checks only, throws on the rest), `from_json`/`to_json` for AdvancedAHB | schema, `.h`, `.cpp` | none | 0.5 d |
| **P2** | Static analytical helpers (`compute_duty_cycle`, `compute_dc_blocking_cap_voltage`, `compute_conversion_ratio`, `compute_turns_ratio`, `compute_lo_min`, `compute_lm_min_for_zvs`, `compute_cb_min`, `compute_zvs_energy_balance`, `compute_dead_time`, `compute_steady_state_flux_excursion`, `compute_transient_flux_excursion`); unit tests | `.cpp`, tests | `NonMonotonic_Gain_Curve`, `DC_Blocking_Cap_Voltage`, `Conversion_Ratio` | 0.5 d |
| **P3** | `process_operating_point_for_input_voltage` (CCM piecewise-linear waveforms for primary, Lm asymmetric ramp, Lo three-region, CT secondary; uses ActiveClampForward magnetizing-ramp pattern); `process_operating_points` (loop over Vin × OPs); diagnostics population | `.cpp` | `OperatingPoints_Generation`, `Asymmetric_Switch_RMS`, `Steady_State_Flux_Excursion`, `Diagnostics_Populated`, `SLUP223_Reference_Design`, `Imbertson_Mohan_Original` | 1.5 d |
| **P4** | Rectifier-type-aware secondaries (CT / CD / FB) with extra-component waveform emission for Lo / Lo2 / Cb | `.cpp` | `RectifierType_CT`, `RectifierType_CD`, `RectifierType_FB` | 1 d |
| **P5** | Transient flux-excursion check (when `inputVoltageStepRange` set); ZVS boundary diagnostic; DCM detection | `.cpp`, tests | `ZVS_Boundaries`, `Transient_Flux_Excursion`, `DCM_Detection` | 1 d |
| **P6** | `process_design_requirements` (size n, Lm, Lo, Cb, Co; verify ZVS at min load) | `.cpp` | `Design_Requirements`, `AN4153_Reference_Design`, `ST_AN2852_Reference_Design`, `Renesas_ISL6745_Reference_Design` | 1 d |
| **P7** | SPICE netlist (V1 classic CT secondary): full Guide §5 contract + IC pre-charge for **all four reactive elements** (Lo, Cb, Co, Lpri) + GEAR + snubbers + Vab probe + sense sources. `simulate_and_extract_operating_points` with analytical fallback. `simulate_and_extract_topology_waveforms`. | `.cpp` | `SPICE_Netlist`, `DivideByZero_Guard` | 1.5 d |
| **P8** | NRMSE acceptance: tighten until ≤ 0.15 for D≈0.3, 0.45, 0.5 fixtures. Failure-discovery phase — expect to fix 1-2 SPICE convergence bugs (IC mismatch on Cb, Llk too small, dead-time too short) | tests + minor SPICE fixes | `PtP_AnalyticalVsNgspice_AtD05`, `PtP_AnalyticalVsNgspice_LightLoad`, `PtP_AnalyticalVsNgspice_HighDuty` | 1.5 d |
| **P9** | `get_extra_components_inputs` returning Lo, Cb, Co (V1) or Lo, Lo2, Cb, Co (V3) (IDEAL + REAL modes) | `.cpp` | `ExtraComponents_Lo_Lo2_Cb_Co` | 0.5 d |
| **P10** | V4 synchronous-rectifier variant (netlist diff) and V5 AHB-Flyback variant (Lo=0, Lm storage) | `.cpp`, tests | `Synchronous_Variant`, `Flyback_Variant` | 0.5 d |
| **P11** | Multi-output (V6) with warn-once cross-regulation message | `.cpp`, tests | `Multiple_Outputs_Warn_Once` | 0.5 d |
| **P12** | AdvancedAHB class, header docstring rewrite (ASCII art, equations §4, references, accuracy disclaimer per Guide §3, **PSHB / Pinheiro-Barbi disambiguation block mirroring PhaseShiftedHalfBridge.h lines 21-42**), waveform plotting, PSHB disambiguation test | `.h`, `.cpp`, tests | `PsHB_Disambiguation`, `AdvancedAHB_Process`, `Waveform_Plotting` | 1 d |

**Total ~10 working days** for V1+V2+V3+V4+V5+V6 coverage with all
three rectifier types and NRMSE-validated.

Dependencies: P1 → P2 → P3; P4 needs P3; P5 ⊥ P4; P6 needs P3+P4; P7
needs P3+P4; P8 needs P7; P9 ⊥ P7; P10/P11 ⊥ P8; P12 last.

---

## 11. Acceptance criteria (per Guide §15)

- [ ] `ninja -j2 MKF_tests` builds clean.
- [ ] All `TestAsymmetricHalfBridge` cases pass.
- [ ] Three `*_PtP_AnalyticalVsNgspice_*` tests run and NRMSE ≤ 0.15
      (D=0.3, D=0.45, D=0.5 — covers the gain-curve branch).
- [ ] Non-monotonic gain curve: M peaks at D=0.5 with M = 1/(2n);
      M(0.3) < M(0.5) and M(0.7) < M(0.5).
- [ ] `lastZvsMargin` changes sign at the predicted Llk threshold.
- [ ] Asymmetric switch RMS test: I_Q1,rms / I_Q2,rms = √(D/(1-D))
      within 5 % at D=0.3.
- [ ] Transient flux excursion exceeds steady-state when
      `inputVoltageStepRange` is set.
- [ ] Schema documented; `Topologies::ASYMMETRIC_HALF_BRIDGE_CONVERTER`
      added.
- [ ] Header has ASCII art, equations §4 with citations, accuracy
      disclaimer, **PSHB / Pinheiro-Barbi disambiguation block**
      mirroring PSHB header lines 21-42.
- [ ] Netlist parses, contains 2 SW switches, snubbers, METHOD=GEAR,
      **IC pre-charge for Cb at (1-D)·Vin** (and Lo, Co, Lpri),
      Vab probe, sense sources, rectifier-type branch (CT/CD/FB).
- [ ] `simulate_and_extract_operating_points` runs SPICE for real,
      falls back to analytical with `[analytical]` tag on failure.
- [ ] `simulate_and_extract_topology_waveforms` references only
      nodes the netlist actually emits.
- [ ] `get_extra_components_inputs` supports IDEAL and REAL for
      Lo, Cb, Co (and Lo2 in V3).
- [ ] Multi-output emits the warn-once message and tags
      cross-regulation as approximate.

---

## 12. Common AHB implementation mistakes (avoid these)

From the research dossier — explicit pitfalls flagged in the
literature:

1. **Cb sized too small** → ΔV_Cb > 5 % of V_Cb → flux excursion
   asymmetry → core saturation at D extremes. Surface
   `lastDcBlockingCapRipple` so it's visible to the user; enforce
   ΔV_Cb ≤ 5 % default.
2. **Ignoring transient DC offset** when Vin steps. V_Cb's RC
   time-constant is large (R_dcr · Lm); transformer can saturate
   during brown-out recovery even if steady-state design is fine.
   Check `lastTransientFluxExcursionEstimate` against core Bsat.
3. **Wrong dead-time** → ZVS lost at light load. Size Lm to provide
   ZVS assist current at min load (S4).
4. **Confusing AHB with PSHB / Pinheiro-Barbi.** Different switch
   count (2 vs 4), different primary voltage (asymmetric 2-level
   vs 3-level), different conversion ratio. The MKF PSHB header
   lines 21-42 already calls this out — the AHB header MUST mirror
   the same disambiguation block.
5. **Treating V_Cb as constant** in small-signal — V_Cb is a slow
   state variable that introduces an extra pole and zero. Document
   in header docstring.
6. **Light-load DCM not modeled correctly** — gain becomes
   load-dependent. Either size for full-load CCM at min duty, or
   detect+flag DCM (P5 does the latter).
7. **Cross-regulation in multi-output** — warn-once.
8. **Asymmetric switch ratings** — at D ≠ 0.5, Q1 and Q2 have
   different RMS currents → using identical FETs leaves headroom on
   one and stresses the other. Surface `lastSwitchRmsCurrentQ1` and
   `lastSwitchRmsCurrentQ2` separately.
9. **Wrong gain-curve branch** — controller designers sometimes
   land on D > 0.5 expecting Vo to keep increasing; the curve is
   non-monotonic. Document explicitly: nominal operating range is
   `D ∈ [0, 0.5]`.
10. **Cold-start inrush** — without `.ic` for Cb at (1-D)·Vin, first
    cycle saturates transformer. Pre-charge **all four** reactive
    elements (Cb, Co, Lo, Lpri) — non-negotiable for SPICE
    convergence.
11. **K = 1.0 in SPICE transformer** — singular factorization;
    always 0.99 (intentional leakage) or 0.999 (minimal leakage).
12. **No parasitic R in SPICE** — lossless reactive networks badly
    conditioned in ngspice. Bake DCR ≈ 50 mΩ into Lo/Lpri, ESR ≈
    5 mΩ into Cb/Co.
13. **Self-driven SR phasing wrong** in V4 — SR FETs derived from
    primary asymmetric voltage need explicit phase relationship;
    sloppy implementation causes secondary shoot-through.

---

## 13. Magnetic structure summary (the bit MKF cares about most)

| Object | Windings | Coupling | DC bias | Gap | MKF mapping |
|---|---|---|---|---|---|
| **Main transformer** | 2 (V2 FB) or 3 (V1 CT) or 2 (V3 CD) | k = 0.99 (Llk usable for ZVS) | 0 (steady state); transient bias possible | typically minimal gap | 2/3-winding transformer; mirror `Flyback` for CT pattern + `PhaseShiftFullBridge` for FB/CD |
| **Output inductor Lo** | 1 | n/a | Iout (continuous) | gapped (energy storage) | single-inductor magnetic; mirror `Buck` |
| **Output inductor Lo2** (V3 CD only) | 1 | n/a | Iout/2 | gapped | mirror `Buck` |
| **DC-blocking cap Cb** | n/a (capacitor) | n/a | (1-D)·Vin | n/a | non-magnetic, lives in netlist + ancillary `Inputs` |

**`process_design_requirements` returns**: the **main transformer**
by default (2 isolation sides). Lo, Lo2 (V3), Cb, Co flow through
`get_extra_components_inputs`.

**`get_extra_components_inputs`** emits:
- V1 CT or V2 FB: `[Lo, Cb, Co]` (3 entries)
- V3 CD: `[Lo, Lo2, Cb, Co]` (4 entries)
- V5 AHB-Flyback: `[Cb, Co]` (no Lo — Lm carries the storage)

**Multi-output**: V6 adds more secondaries to the same transformer
(plus more rectifier diodes / SR FETs and per-secondary Lo + Co).
Use the same per-output indexing as DAB (`output_voltages[]`,
`output_currents[]`, `o<N>` netlist suffix). Print warn-once for
cross-regulation per Guide §4.7.

**Novel constraint MKF doesn't have today**: the **transient flux
excursion check** under Vin step. Today's adviser sizes the core for
steady-state Bpk; AHB can saturate during transients. Investigate
either:

- (a) extend DesignRequirements with a `transientPeakFluxDensity`
  field (cleaner; benefits all topologies that have this risk —
  Forward, ACF, AHB),
- (b) bake the transient check into AHB-specific
  `process_design_requirements` and surface as a diagnostic
  (less canonical but self-contained).

Recommend (b) for the initial add. (a) becomes a follow-up if
Forward / ACF want to surface the same check.

---

## 14. Checklist for the implementing agent

Before each phase:

- [ ] Re-read `CONVERTER_MODELS_GOLDEN_GUIDE.md`.
- [ ] **Re-read PhaseShiftedHalfBridge.h lines 21-42** to absorb the
      disambiguation block — the AHB header must mirror it
      symmetrically (PSHB ≠ AHB; AHB ≠ split-cap HB; AHB ≠ HB-LLC).
- [ ] Verify the latest `PhaseShiftedFullBridge.cpp` for
      rectifier-type-aware (CT / CD / FB) pattern.
- [ ] Verify the latest `ActiveClampForward.cpp` for
      DC-blocking-cap + asymmetric-magnetizing-ramp pattern.
- [ ] Verify the latest `SingleSwitchForward.cpp` /
      `TwoSwitchForward.cpp` for forward-class secondary pattern.
- [ ] Verify `PwmBridgeSolver.h` — its scaffold can be extended for
      the AHB asymmetric duty (non-equal half-cycles).
- [ ] Build with `ninja -j2` in the build dir, not `cmake --build`
      (per user memory).
- [ ] After MAS schema edits, copy `MAS.hpp` from source to
      `build/MAS/MAS.hpp` (per user memory).
- [ ] Mark each phase as complete only when its NRMSE / sign-change
      / build-clean criteria are met.

After P12 (final):

- [ ] Update CLAUDE.md / user-memory entry on converter quality
      tiers to include AHB as DAB-quality.
- [ ] Open follow-up issues for: two-transformer AHB; AHB-LLC
      hybrid (datacenter PSU); GaN-AHB at 1-2 MHz;
      `transientPeakFluxDensity` field generalization to
      Forward / ACF.
- [ ] Update `FUTURE_TOPOLOGIES.md` to mark AHB as completed.
- [ ] Consider adding a cross-reference from
      `PhaseShiftedHalfBridge.h` disambiguation block to the new
      `AsymmetricHalfBridge.h` so readers can find the proper
      Imbertson-Mohan AHB.

---

## 15. Reference: equations cheat-sheet (paste into header docstring)

```
TOPOLOGY OVERVIEW (V1 classic AHB, CT secondary)
────────────────────────────────────────────────
    Vin ──┬── Q1 ──┬── A ── Cb ── Np ── B ─── GND_pri
          │        │
          │       SW1
          │        │
          ├── Q2 ──┘
          │
         Cin
          │
         GND_pri

    CT secondary:
          Ns_a ─── Ds_a ──┬── Lo ── Vout+
                          │              ┌─ Co ─┐
          Ns_b ─── Ds_b ──┘              │      │
                                       Rload    │
                                        │       │
                                       0       0

MODULATION CONVENTION
─────────────────────
Q1 ON for D·Tsw, Q2 ON for (1−D)·Tsw, complementary with dead-time td
Operating range: D ∈ [0, 0.5] (chosen branch of non-monotonic gain curve)

| Sub-int     | Length   | Q1 | Q2 | v_pri      | i_Lm slope    | rectifier   |
|-------------|----------|----|----|------------|---------------|-------------|
| A           | D·Tsw    | ON | OFF| +(1-D)·Vin | +(1-D)·Vin/Lm | top diode   |
| B (dead)    | td       | OFF| OFF| ringing    | non-linear    | freewheel   |
| C           | (1-D)·Tsw| OFF| ON | -D·Vin     | -D·Vin/Lm     | bottom diode|
| D (dead)    | td       | OFF| OFF| ringing    | non-linear    | freewheel   |

KEY EQUATIONS (CCM, ideal)
──────────────────────────
DC-blocking cap voltage : V_Cb = (1−D)·Vin                     — Imbertson-Mohan eq. 3
Primary voltage         : v_pri = +(1-D)·Vin during D·Tsw,
                          v_pri = -D·Vin    during (1-D)·Tsw    ← ASYMMETRIC, automatic flux balance
Conversion ratio (CT/FB): Vo = 2·D·(1-D)·Vin / n               — non-monotonic, peak at D=0.5
Conversion ratio (CD)   : Vo = D·(1-D)·Vin / n
Maximum gain at D=0.5   : Vo,max = Vin / (2·n)
Switch peak voltage     : V_Q,pk = Vin (each switch)
Switch RMS              : I_Q1,rms = I_pri,Q1 · √D ;  I_Q2,rms = I_pri,Q2 · √(1-D)
                          (UNEQUAL when D ≠ 0.5)
ZVS condition           : 0.5·(Llk + Lm,refl)·I_pri²  ≥  2·Coss·Vin²
Dead-time               : td ≥ π·√(Llk · 2·Coss)

Steady-state flux excursion: ΔΦ = D·(1-D)·Vin·Tsw / Np
Bpk                       : D·(1-D)·Vin·Tsw / (2·Np·Ae)
Core utilization at D=0.5  ≈ 70-75% of B-H curve, degrades to 30% at D=0.1
Transient B (Vin step)    : may exceed steady-state Bpk during V_Cb settle
                            time-constant τ ≈ R_dcr·Lm

Sizing — Cb (ripple ≤ 5 % V_Cb):
  Cb ≥ I_pri,pk · D / (fsw · ΔV_Cb)                            — ON AN-4153 eq. 16

Sizing — Lo (CCM):
  Lo ≥ Vo · (1 − 2·D·(1-D)) / (ΔILo · fsw)
       (at D=0.5: Lo ≥ 0.5·Vo / (ΔILo · fsw))

Sizing — Lm for ZVS (assist current target):
  Lm ≤ ((1-D)·Vin·D·Tsw) / (2·Im,target)

DISAMBIGUATION (mirror PhaseShiftedHalfBridge.h lines 21-42)
────────────────────────────────────────────────────────────
This is NOT the Pinheiro-Barbi 3-level half-bridge (4 switches, two
legs, split-cap input, 3-level ±Vin/2 primary). The MKF
`PhaseShiftedHalfBridge` (`Pshb`) class implements the Pinheiro-Barbi
version. AHB is a separate topology with 2 switches, single-leg, and
2-level asymmetric primary voltage.

This is NOT a symmetric split-cap half-bridge (50% duty each, primary
swing ±Vin/2, requires current-mode control to prevent flux walking).
AHB's flux balance is automatic via the DC-blocking cap.

This is NOT a Half-Bridge LLC (resonant tank, frequency-modulated, no
PWM duty). AHB is PWM-controlled with naturally-resonant ZVS dead-time
transitions.

REFERENCES
[1] P. Imbertson & N. Mohan, "Asymmetrical Duty Cycle Permits Zero
    Switching Loss in PWM Circuits with No Conduction Loss Penalty,"
    IEEE Trans. Industry Applications, 29(1):121-125, 1993.
[2] R. L. Steigerwald, "Designing High-Efficiency Asymmetric Half-Bridge
    PWM Converters," TI SLUP223 (Power Supply Design Seminar).
[3] ON Semi AN-4153, "Designing Asymmetric PWM Half-Bridge DC/DC
    Converter with FPS."
[4] STMicro AN2852, "EVL6591-90WADP - 90 W AC-DC AHB adapter."
[5] R. P. Severns & G. E. Bloom, "Modern DC-to-DC Switchmode Power
    Converter Circuits," (sometimes attributed); also Pressman,
    Switching Power Supply Design 3e Ch. 16.
[6] Erickson & Maksimović, "Fundamentals of Power Electronics" 2e §6.3
    (DC analysis), §11.1.4 (small-signal).
[7] TI SLUA121, "The Current-Doubler Rectifier: An Alternative
    Rectification Technique."

ACCURACY DISCLAIMER
SPICE model omits MOSFET Coss variation with V_DS, body-diode
reverse recovery, gate-driver propagation delay, transformer
inter-winding capacitance, and core-loss resistance. Predicted η is
an upper bound; expect 2-4 % loss vs. measured at full load.
Multi-output AHB cross-regulation is approximate unless per-secondary
leakage inductance is provided (warn-once message printed in
process_operating_points).

POLARITY / FLOW CONVENTION
Power flows Vin → Vo (unidirectional). Vo > 0. Iout positive when
flowing from Vo through Rload to GND. Operating range D ∈ [0, 0.5]
(chosen branch of non-monotonic gain curve).
```
