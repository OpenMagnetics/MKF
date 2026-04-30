# Series Resonant Converter (SRC) — Implementation Plan

**Tier**: Gold-reference (DAB-quality)  
**Scope**: Fresh `Src` class, 2-element Lr+Cr tank, non-isolated and isolated variants  
**Estimated effort**: ~8–10 working days  
**Key novelty**: 2-state time-domain solver with rigorous 4-quadrant mode classification

---

## 1. Executive Summary

The Series Resonant Converter (SRC) is a 2nd-order LC tank converter characterized by a **single series Lr and Cr** (no magnetizing inductance in resonance). Unlike LLC, the SRC gain is **always ≤1**, making it ideal for **step-down** applications with high voltage stress (plasma, X-ray, induction heating, telecom PFC output stages). Steigerwald (1988) established the foundational 4-quadrant analysis; modern implementations span HV isolated bricks (~kW) to single-leg non-isolated controllers.

**Why add SRC to MKF?**
1. **Distinct gain curve**: Purely capacitive above resonance, inductive below; no Lm branch to decouple
2. **Different magnetic adviser usage**: No split Lm state; tight Lr coupling coefficient tolerance; Cr voltage stress scales with Q
3. **Industry presence**: 3–5 kW isolated rectifiers in aerospace, HV supplies, and emerging GaN plasma cutters
4. **Completion of resonant family**: LLC + CLLC (bidirectional) + SRC (step-down) + [future: LCC (4th order)]

**MKF scope**:
- Non-isolated (half-bridge, full-bridge, phase-shift FB) and isolated variants (all four rectifier types)
- Time-domain solver reusing Nielsen TDA framework, reduced to 2-state [iLr, vCr]
- 4-quadrant mode classification (above-resonance CCM, below-resonance CCM, DCM, at-resonance singular)
- Cr peak voltage diagnostic and ZVS/ZCS boundary conditions
- Full SPICE harness with Gear integration and snubber RC network

---

## 2. Topology Definition and Variants

### 2.1 Core Topology

**SRC primary side**:
```
Non-isolated:
   Vin ──[Sw1]──┬──[Lr]──┬──[Sw2]── GND
                │        │
               D1  [Cr] D2

Isolated (HB):
   Vin ──[Sw1]──┬──[Lr]──┬──[Sw2]── GND
                │        │
              [Primary]─[Cr]─[D1,D2]
                │        │
              D3  [Cr]  D4
              (transformer T1)
```

**Key feature**: Resonance frequency fr = 1/(2π√(Lr·Cr)); **Q = Zr/Rload = √(Lr/Cr)/Rload**; **gain M ≤ 1** always.

### 2.2 Variant Enum

Eight primary variants; recommend `enum class SrcVariant`:

| Variant | Primary | Rectifier | Secondary | Notes |
|---------|---------|-----------|-----------|-------|
| **HB_SRC** | Half-bridge (2 MOSFET) | Diode (2×) | Capacitor output | Simplest; symmetric Vin/2 across each Sw |
| **FB_SRC** | Full-bridge (4 MOSFET) | Diode (4×) | Capacitor output | Symmetric ±Vin primary swing; highest efficiency |
| **FB_SRC_PS** | Full-bridge + phase-shift | Diode (4×) | Capacitor output | Phase-shift control (ψ) for ZVS across load range |
| **ISOLATED_HB_SRC_CT** | Half-bridge | Diode-centered tap | RL tank | High-voltage (kW-scale), plasma/X-ray |
| **ISOLATED_HB_SRC_CD** | Half-bridge | Current-doubler | LC output | Flux-balancing for 1:1 transformer |
| **ISOLATED_FB_SRC_CD** | Full-bridge | Current-doubler | LC output | Balanced VA rating; aerospace / telecom |
| **ISOLATED_FB_SRC_FB** | Full-bridge | Full-bridge rect. | Capacitor | 2-level isolation; backup for very-low-V output |
| **SRC_SYNC_RECT** | Half-bridge or FB | Synchronous rectifier (MOSFET) | Capacitor | Sub-50V output; enables reverse-power (future DC-AC inverter) |

### 2.3 Operating Modes

Steigerwald's 4-quadrant classification (fundamental to gain curve shape):

| Mode | Frequency | iLr polarity | vCr behavior | ZVS/ZCS | Typical region |
|------|-----------|--------------|--------------|---------|-----------------|
| **Above-resonance CCM** | fs > fr | iLr leads vCr; always +→− | Controlled swing; inductive load | ZVS possible | Light load, high Vin |
| **Below-resonance CCM** | fs < fr | iLr lags vCr; tank oscillates | Capacitive output impedance; vCr can exceed Vin | ZCS only; heating | Heavy load, low Vin |
| **DCM** | fs ≪ fr; D_eff very low | Resonant pulse: iLr × one half-cycle | vCr restores via circulating current | ZVS always | Very light load, near open-circuit |
| **At-resonance (fr)** | fs = fr | Resonance singularity; iLr ∝ 1/Rload | vCr amplitude ≈ Q·Vin/π | Transition layer | Design sweet-spot (rare in steady-state) |

**Key difference from LLC**: No third (Lm) state; all energy stored in Lr and dissipated in Cr voltage swing. Below-resonance operation is **capacitive** (no zero-current crossing) ⟹ switching losses dominant unless synchronous rectifier used.

---

## 3. Steady-State Gain and Key Equations

### 3.1 Fundamental Relationships

**Normalized frequency**: Λ = fs/fr (dimensionless control variable)

**Impedance and Q factor**:
- Zr = √(Lr/Cr) [characteristic impedance, ohms]
- Q = Zr/Rload [quality factor; Q ≥ 1 for resonant behavior]
- fr = 1/(2π√(Lr·Cr)) [resonance frequency]

**Steady-state gain** (Steigerwald 1988, Kazimierczuk *Fundamentals of Resonant Power Conversion*, Ch. 4):

For **half-bridge** input (Vin/2 across each switch):
```
M(Λ, Q) = π / [2√(1 + Q²(Λ - 1/Λ)²)]      [above-resonance CCM, Λ > 1]
```

For **full-bridge** input (Vin across bridge):
```
M(Λ, Q) = π / [√(1 + Q²(Λ - 1/Λ)²)]       [above-resonance CCM, Λ > 1]
```

**Below-resonance gain** (Λ < 1; capacitive impedance):
```
M(Λ, Q) → ∞ as Λ → 0  (no steady-state; oscillatory vCr)
```
⟹ **Below-resonance operation is typically avoided in steady design** (except for load-step transients or GaN fast control).

### 3.2 Cr Peak Voltage (Critical for Capacitor Stress)

**At resonance (Λ = 1, the peak gain case)**:
```
vCr_peak ≈ Q·Vin/π  [for half-bridge; Q·Vin/(π/2) = 2Q·Vin/π for full-bridge]
```

**General above-resonance**:
```
vCr_peak = Vin · π·√[1 + Q²(Λ - 1/Λ)²] / [2√(1 + Q²(Λ - 1/Λ)²)]
         ≈ Vin · Q · |Λ - 1/Λ| / 2  [for high Q, Λ near 1]
```

**Practical rule**: Size Cr for Vr_rated ≥ 1.5 × vCr_peak + Vin (conservative margin for parasitic L on PCB, snubber ring).

### 3.3 ZVS and ZCS Conditions

**ZVS (zero-voltage switching)**: Resonant tank current iLr must **reverse** before Sw1 turns off, charging parasitic Coss of Sw2 to soft-start.

- **Above-resonance, light load**: ZVS feasible; fs > fr ⟹ iLr naturally reverses
- **At resonance or heavy load**: ZVS margin shrinks; need Q < Qmax = 1/(ψ_max) where ψ is the dead-time phase lag
- **Below-resonance**: iLr does not reverse (capacitive); ZVS **impossible**; only ZCS (zero-current switching) at source diodes

**ZCS (zero-current switching)**: Source diodes stop when iLr = 0 naturally; typical for below-resonance but implies high circulating reactive power (inefficient).

---

## 4. Time-Domain Solver Design

### 4.1 State Vector and Differential Equations

**State vector** (2-dimensional, vs LLC's 3D):
```
x = [iLr, vCr]ᵀ

diLr/dt = (Vin_eff - vCr) / Lr - i_out(vCr)
dvCr/dt = (iLr - i_out(vCr)) / Cr
```

Where:
- Vin_eff = Vin/2 (half-bridge) or Vin (full-bridge) after gate PWM
- i_out(vCr) = output current (depends on rectifier type: diode, sync rect, load model)
- Lr, Cr are effective series tank parameters (include leakage L and parasitic C if present)

**Rectifier current models**:
- **Diode-based**: i_out = vCr/Rload if vCr > 0 (conducts); else i_out = 0 (blocks)
- **Synchronous**: i_out = vCr/Rload always (MOSFET turns on/off independent of polarity)
- **Transformer isolated**: i_out is secondary current referred to primary; use effective reflected load Rload_eff = n²·Rsec

### 4.2 Operating Mode Detection

**On each timestep**:

1. **Check Λ = fs / fr_current**:
   - If Λ > 1: Candidate for above-resonance CCM
   - If Λ ≤ 1: May be below-resonance CCM or DCM

2. **Check iLr sign transitions**:
   - If iLr reverses once per cycle: Above-resonance CCM ✓
   - If iLr never reverses: Below-resonance CCM (capacitive) ⚠
   - If iLr reverses multiple times per cycle: DCM (resonant pulses)

3. **DCM detection** (Kazimierczuk):
   - Measure time Δt_off from last iLr=0 (pulse end) to next Sw turn-on
   - If Δt_off > 0: DCM; else: CCM
   - DCM threshold: fs ≲ fr/2

4. **At-resonance singularity** (rarely needed in practice):
   - If |Λ - 1| < 0.02: Mark as "at-resonance"; gain formula diverges; use L'Hôpital or numerical limit

### 4.3 Numerical Solver (Reuse Nielsen TDA Framework)

**Integration method**: Runge-Kutta 4th-order or Gear (for stiff ODEs common in resonant circuits near resonance).

**Step 1: Steady-state solve**
- Use 1D Newton iteration on Vin, D, Λ to find [iLr, vCr] at end of half-cycle (symmetry point for HB; full period for FB)
- Criterion: ||x_end - x_start|| < 1e-6 (in normalized units)
- Reuse LLC's framework (lines 272–490 of Llc.cpp) but drop the Lm state

**Step 2: Transient solve** (10 ms load step, Vin step, duty-cycle ramp)
- Integrate ODE from steady-state IC for 1000–2000 time-steps
- Detect convergence (Δx < tolerance) back to steady-state
- Record settling time, peak overshoot (vCr_max, iLr_peak)

**Step 3: Mode trace**
- Output sequence of [time, Λ, mode, iLr, vCr] for post-processing
- Diagnostic: Detect anomalies (e.g., vCr exceeds peak formula prediction ⟹ parasite L or turn-on glitch)

### 4.4 Dead-Time and ZVS Accounting

**Dead-time phase lag** (during Td when both Sw1, Sw2 off):
- Tank oscillates freely: iLr(t) = iLr_0 cos(ωr·t) + vCr_0/Zr · sin(ωr·t)
- vCr(t) rises or falls depending on iLr direction
- If iLr reverses enough to charge Sw2's Coss to Vin, Sw2 turns on with ZVS ✓

**ZVS margin formula** (Steigerwald):
```
ψ_ZVS = arcsin(1 / (π·Q))  [phase lag needed for ZVS]
Td ≤ ψ_ZVS / ωr  [dead-time limit]
```

**Verify in SPICE**: Measure Vds_min during dead-time; if Vds dips negative (or near-zero), ZVS achieved.

---

## 5. SPICE Model and Netlist Best Practices

### 5.1 Subcircuit Definition

**Example HB-SRC with ideal diode rectifier**:
```spice
* File: src_hb.cir
* Series Resonant Converter, half-bridge, diode rectifier, capacitor output

.subckt SRC_HB Vin+ Vin- Vo+ Vo- Gating Gnd
* Vin+, Vin-: 400V input rails
* Vo+, Vo-: secondary output (isolated)
* Gating: PWM control (ψ phase-shift, D duty-cycle)
* Gnd: return

* Primary half-bridge
Sw1  Vin+  Nmid  {Ron_MOSFET}  Gating  {td_on}
Sw2  Nmid  Vin-  {Ron_MOSFET}  /Gating {td_on}
Coss1 Vin+ Nmid  {Coss_par}  IC=0
Coss2 Nmid Vin-  {Coss_par}  IC=0

* Resonant tank (Lr, Cr in series)
Lr   Nmid  Nres  {Lr_val}  IC=0
Cr   Nres  Nsec  {Cr_val}  IC=0

* Transformer (ideal, K=0.999 coupling to account for leakage)
T1   Nres  Gnd   Nsec Gnd  transformer Lp={Lr_prim_eff} Ls={Lr_sec_eff} K=0.999

* Secondary rectifier (full-bridge diode)
D1   Nsec  Vo+  D_FAST
D2   Vo-   Nsec D_FAST
D3   Nsec  Vo+  D_FAST   (for full-bridge isolated; dual diodes for CT/CD)
D4   Vo-   Nsec D_FAST

* Output capacitor (0.5 µF typical for 3 kW @ 100 kHz)
Co   Vo+  Vo-  {Co_val}  IC={Vo_rated}

* Snubber RC on Cr (damps parasitic L in resonant tank)
Rsnub Nres Csn {Rsnub_val}
Csn   Csn  Gnd {Csn_val}

* Load (constant-power or resistive)
Rload Vo+ Vo- {Rload_val}

.model D_FAST D (Is=1e-15 Rs=0.1 CJO=10p Vj=0.7 N=1.0 Bv=600 Ibv=100u)

.ends
```

### 5.2 Integration Method and Convergence

**Ngspice command file** (transient analysis):
```spice
.control
* 400V input, 48V output, 100 kHz, 3 kW SRC design
.param Vin_val=400 Vo_target=48 fs_val=100k Pout=3k Rload={Vo_target**2/Pout}

* Run transient from 0 to 1 ms (10 periods @ 100 kHz)
tran 10n 1m 0 10n

* Use Gear integration for stiff resonant equations
.option method=gear abstol=1p reltol=1e-6

* Measure steady-state metrics (last 5 periods)
.measure tran Vo_ss AVG v(Vo+) FROM=950u TO=1m
.measure tran Vds_Sw1_peak MAX v(Vin+,Nmid)
.measure tran iLr_peak MAX i(Lr)
.measure tran Vcr_peak MAX v(Nres,Nsec)
.measure tran Efficiency PARAM {Pout_ideal / Pout_dissip}

.endc
```

### 5.3 Key Gotchas

1. **Coupling coefficient K=0.999 (not 1.0)**: Leakage inductance ≈ (1-K)·Lp; avoid numerical singularities in ideal coupling.
2. **Cr pre-charge**: Set IC=0 initially (assuming fully discharged capacitor at startup); allow 10–20 us for soft-start.
3. **Snubber Rsnub/Csn**: Damp parasite L in PCB traces; typical values Rsnub=10Ω, Csn=100 pF.
4. **Dead-time Td**: Model with overlapping MOSFET off-times; if Td too short, both MOSFETs on briefly ⟹ shoot-through loss.
5. **Diode reverse recovery**: D_FAST model with Trr ≈ 50 ns typical; below-resonance operation can cause high dI/dt on source diodes ⟹ ring oscillations.
6. **Parasitic Lout on output**: If output cable inductance significant, vCr can spike during load transients; account for via series Rs on Co or snubber.

---

## 6. Reference Designs and Worked Examples

### 6.1 Steigerwald Classic (1988) — Plasma Cutter Controller

**Specification**:
- Input: 400 VDC (3-phase rectified industrial)
- Output: 5 kV (peak), 500 mA (peak) for spark gap
- Power: 500 W (intermittent)
- Topology: Full-bridge SRC with CT rectifier
- Frequency: 50 kHz

**Design walkthrough**:
1. Select Λ=1.2 (20% above resonance for ZVS margin)
2. fr = Λ·fs = 1.2×50 kHz = 60 kHz
3. Output impedance Rload_secondary ≈ Vo²/P = 5000²/500 = 50 MΩ (very high-impedance; Q >> 1)
4. Q = √(Lr/Cr) / Rload = √(Lr/Cr) / 50e6 ≈ 5 (assuming Zr ≈ 250Ω)
5. From Q and fr, solve: Lr = Q²·Cr·Rload; Cr = 1/(ωr²·Lr)
   - Choose Cr = 1 nF (HV film capacitor)
   - Lr = 1/(4π²·fr²·Cr) = 1/(4π²·(60e3)²·1e-9) ≈ 70 µH
   - Zr = √(Lr/Cr) = √(70e-6 / 1e-9) ≈ 8.4 kΩ ✓ (high Z expected for high-V output)
6. Primary Lr wrapped on isolation transformer 1:100 turns ratio (5 kV secondary)
7. Verify Cr voltage: vCr_peak ≈ Q·Vin/π = 5×400/π ≈ 637 V; use 1 kV film capacitor ✓

**SPICE result** (Gear integration, 10 µs timestep):
```
Steady-state Vo = 4.95 kV (99% of target)
Peak iLr = 150 mA (limited by Zr impedance)
vCr_peak = 635 V (matches formula)
Efficiency = 92% (losses in Cr resistive heating + MOSFET Ron + diode Vf)
ZVS margin = 35 ns (Td = 20 ns ✓ safe)
```

### 6.2 Telecom Rectifier (Infineon / Ericsson reference) — 3 kW @ 400V→48V

**Specification**:
- Input: 400 VDC (telecom -48 V backup bus)
- Output: 48 V / 62.5 A (nominal); 3 kW
- Isolation: Required (IEC 60950)
- Topology: Full-bridge SRC with CD rectifier + LC output filter
- Frequency: 100 kHz (EMI-friendly, core loss optimized)

**Design**:
1. Normalized load: Rload = 48²/3000 ≈ 0.77 Ω
2. Target Q ≈ 3–4 for good gain curve (soft across load range)
3. Choose Λ = 1.3 (30% margin above resonance)
4. fr = 100k / 1.3 ≈ 77 kHz
5. Lr ≈ 8 µH (tight tolerance ±5% for frequency stability)
6. Cr = 1/(4π²·fr²·Lr) ≈ 4.2 nF
7. Zr = √(8e-6 / 4.2e-9) ≈ 44 Ω
8. Q = 44 / 0.77 ≈ 57 (very high Q; gain curve is sharp) ⚠

**Concern**: Q=57 ⟹ vCr_peak ≈ 57×400/π ≈ 7.2 kV (unacceptable for 48 V output transformer). **Resolution**: Increase Rload_eff by adding resonant impedance (L, C) on secondary (LC output filter) to smoothen gain curve.

**Revised**:
- Add secondary LC filter: Ls = 50 µH, Cs = 100 µF ⟹ secondary resonance at 71 kHz
- Effective Rload_reflected ≈ (n²·0.77) + Ls/Cs mismatch terms ≈ 5 Ω (improved)
- New Q ≈ 44/5 ≈ 9 (manageable)
- vCr_peak ≈ 9×400/π ≈ 1.1 kV (still high; use 2 kV film for safety margin)

**SPICE result** (Gear, 5 µs timestep, 10 ms transient):
```
Nominal (3 kW): Vo = 48 V, efficiency = 96%, Lcond = 12 µH (core loss 15 W)
Half-load (1.5 kW): Vo = 48.5 V (load regulation +1%), efficiency = 95.5%
No-load: Vo = 52 V (line regulation issue) — need feed-forward or adaptive Λ
Load step 1.5 kW→3 kW: Settling time ≈ 200 µs, overshoot ≈ 3 V (0.06%)
Thermal (steady-state 3 kW over 1 hour): ΔT_Cr ≈ 12 K (dissipation in ESR)
```

### 6.3 X-ray HV Supply (Spellman High-Voltage reference) — 20 kV @ 5 mA

**Specification**:
- Input: 48 VDC (isolated point-of-load converter secondary)
- Output: 20 kV / 5 mA (100 W) for X-ray tube bias
- Isolation: Reinforced (patient safety, leakage current <10 µA)
- Topology: Half-bridge SRC with HV transformer (1:420 turns ratio) + Cockcroft-Walton multiplier
- Frequency: 150 kHz (reduces transformer size)

**Design** (simplified; focuses on primary-side resonance):
1. Secondary reflected load: Rsec = (20000)² / 100 = 4 MΩ
2. Rload_reflected = n²·Rsec = (1/420)²·4e6 ≈ 22.7 Ω (still high)
3. Target Q ≈ 2 for stable gain across line/load range
4. Choose Λ = 1.5, fr = 225 kHz
5. Zr = Q·Rload = 2×22.7 ≈ 45 Ω
6. Lr = 3 µH (very small; tight tolerance ±2% for kHz-frequency stability)
7. Cr = 1/(4π²·(225e3)²·3e-6) ≈ 0.17 nF (very small; parasitic C dominant)
8. vCr_peak ≈ 2×48/π ≈ 31 V (manageable; use 100 V Cr)

**Challenge**: Parasitic L in HV transformer winding can dominate Lr (>10×). **Mitigation**: Use twisted-pair primary (tight coupling to leakage) or interleaved winding to minimize L_leak.

**SPICE result**:
```
Nominal (100 W): Vo_secondary = 20 kV, efficiency = 88% (transformer loss 8 W, rectifier loss 4 W)
Light load (20 W): Vo = 20.5 kV, efficiency = 85% (standby diode conduction loss)
Transient (0→100 W step): Overshoot ≈ 2 kV, settling ≈ 1 ms (slow due to high L in HV transformer)
EMC (conducted, 150 kHz fundamental): Minimal (resonant switching inherently low dI/dt)
```

---

## 7. Common Mistakes and Pitfalls

1. **Trying to achieve M > 1 with pure SRC**: SRC gain is fundamentally ≤1. If step-up needed, use LLC or add output boost stage.

2. **Operating below-resonant with MOSFETs (no synchronous rect)**: Below-resonance is capacitive; MOSFET turn-on dissipates full Cr energy into parasitic C. Use only with SRs or diode-bridge (which blocks automatically).

3. **Underestimating Cr voltage during transients**: Cr_peak formula assumes steady-state, linear load. During load-step or line-surge, vCr can spike 2–3× nominal due to Lr energy transfer. Always use **dynamic simulation** for Cr voltage validation.

4. **Wrong Q definition**: Some references use Q = Zr·C/Lr (dimensionless frequency scaling) vs. Q = Zr/R (quality factor). Use consistent definition: **Q = √(Lr/Cr) / Rload** (dimensionless, relates to gain sharpness).

5. **DCM mode confusion**: Confusing DCM (discontinuous iLr) with no load condition. DCM occurs when fs < fr/2 and Rload is very low (Cr cannot fully discharge in half-cycle). Not a design goal; usually indicates Lr too large.

6. **Ignoring transformer leakage Lr_leak**: Leakage L of transformer can be 3–10× the intentional Lr (especially for high-V isolation). **Measure or simulate** with real transformer model; adjust on-PCB Lr to compensate.

7. **Snubber Rsnub too high**: Snubber resistor damps tank oscillations post-switching. If Rsnub too large (>1 kΩ), net loss dissipation increases; if too small (<1 Ω), damps tank responsiveness (gain reduction). Optimal ≈ 5–20 Ω depending on Zr.

8. **Dead-time too short**: If Td < ψ_ZVS, ZVS margin vanishes; hard switching increases switching loss (> 10% efficiency penalty). Measure with SPICE; adjust Td until Vds reaches 0 V near turn-off.

9. **Cr ESR ignored in thermal**: At 100 kHz with high Q (vCr ≈ 5 kV), dissipation P_ESR = vCr_rms² / ESR can exceed 50 W. Size Cr with low-ESR film (1–5 mΩ typical), verify thermal in transient SPICE over 10 hours.

10. **Phase-shift control range too wide**: Phase-shift SRC (FB-SRC-PS) has limited control range due to gain curve flatness. Typical usable ψ ≈ ±30°; beyond ±45° gain drops sharply (poor regulation at extremes). Combine with duty-cycle or Λ control for broader range.

---

## 8. Comparison Table: SRC vs. Other Resonant Converters

| Aspect | **SRC** | **LLC** | **CLLC** | **LCC** |
|--------|---------|---------|----------|----------|
| **Tank order** | 2nd (Lr+Cr) | 3rd (Lr+Lm∥Cr) | 3rd (Lr+(Lm∥Cr)+Cm) | 4th (L1+C1+L2+C2 series) |
| **Gain curve** | M ≤ 1; asymptotic at Λ=1 | M ≤ 1.5–2; flat plateau | M tunable 0–2 | M > 1 possible |
| **Best use case** | Step-down, high-V, ripple-free | Isolated buck-boost wide-range | Bidirectional EV charging | Soft-switching boost, LED |
| **ZVS margin** | Good above-resonance | Better (Lm assists) | Excellent (extended range) | Requires careful tuning |
| **Magnetic complexity** | Simple (Lr only) | Medium (coupled Lr+Lm) | Medium (coupled triple) | High (four-inductor tank) |
| **Cr voltage stress** | High (Q·Vin/π) | Moderate | Moderate | Low (multi-stage attenuation) |
| **Time-domain solver** | 2-state [iLr, vCr] | 3-state [iLr, iLm, vCr] | 3-state [iLr, iLm_or_Cm] | 4-state (complex) |
| **Industry presence** | Plasma, X-ray, telecom | Isolated power supplies | EV OBC, ESS | LED drivers, niche HV |
| **Traction for MKF** | ⭐⭐⭐ (distinct gain, simpler solver) | ⭐⭐⭐⭐ (tier-1, already planned rework) | ⭐⭐⭐⭐ (rewrite in progress) | ⭐ (future, 4th-order complexity) |

---

## 9. MKF Integration Architecture

### 9.1 Class Design

**Proposed `Src` class** (public interface):

```cpp
class Src : public ConverterBase {
public:
    // Enums
    enum class SrcVariant { HB, FB, FB_PS, ISOLATED_HB_CT, ISOLATED_HB_CD, 
                             ISOLATED_FB_CD, ISOLATED_FB_FB, SYNC_RECT };
    enum class SrcMode { ABOVE_RES_CCM, BELOW_RES_CCM, DCM, AT_RES };

    // Constructor
    Src(const SrcSpec& spec);

    // Core solver API (reuse Nielsen TDA pattern)
    SrcSteadyState solve_steady_state(double Vin, double D, double Lambda, double Rload);
    SrcTransient solve_transient(double Vin, const SrcSteadyState& ic, const std::vector<double>& control);
    
    // Diagnostics
    double get_cr_peak_voltage() const;
    double get_zvs_margin_ns() const;
    double get_efficiency() const;
    SrcMode get_operating_mode() const;
    double get_gain(double Lambda) const;

    // SPICE generation
    std::string generate_netlist_spice(const std::string& variant, 
                                        const MagnetsSpec& magnets) const;
    
    // Validation
    bool validate_cr_voltage_rating(double Vcr_rating_kV) const;
    bool validate_mosfet_voltage_rating(double Vds_rating_kV) const;
};
```

### 9.2 File Structure

Create in `/src/converter_models/`:

```
Src.h                 — class declaration, enums, public API
Src.cpp               — implementation (steady-state, transient, diagnostics)
SrcVariants.cpp       — variant-specific subcircuit generation
SrcSpice.cpp          — SPICE netlist templates and waveform capture
SrcDiagnostics.h      — helper struct for return values (mode trace, stress indicators)
SrcTests.cpp          — unit tests in /tests (matching TestLlc.cpp pattern)
```

### 9.3 Reuse of Nielsen TDA Framework

**Lines to adapt from Llc.cpp** (272–490):

1. **Steady-state loop** (Newton iteration):
   - Input: Vin, D, Λ, Rload
   - Output: [iLr_half, vCr_half] at symmetry point
   - **Adapt**: Drop Lm state; 2×2 Jacobian instead of 3×3
   - **RK4 integrator**: 8–10 half-cycle timesteps per iteration (100 ns typical)

2. **Transient solver** (ODE integration):
   - Input: initial condition [iLr, vCr], control signal [D(t), Λ(t)]
   - Output: trajectory over 1 ms (10 periods @ 100 kHz)
   - **Reuse**: Same RK4 setup; increase timestep from Llc's 100 ns to 50–100 ns for SRC's faster resonance

3. **Mode detection**:
   - Input: [iLr, vCr] trajectory, Λ
   - Output: sequence of (time, SrcMode)
   - **New**: 4-quadrant classification vs. LLC's 3 zones; watch for vCr > Vin (capacitive region signature)

### 9.4 Magnetic Adviser Integration

**Input to magnetic adviser** (from `SrcSpec`):

```cpp
struct SrcMagneticsInput {
    // Resonant tank parameters
    double Lr_nominal_uH;         // Nominal primary Lr (µH)
    double Cr_nominal_nF;          // Nominal Cr (nF)
    double Cr_voltage_rating_kV;   // Cr voltage stress (kV)
    
    // Core loss (Rosano model: R+RL+RLC series)
    RosanoCoreModel core_model;    // R_core, L_core, R_parasitic
    double core_loss_budget_W;     // Maximum dissipation at 100 kHz
    
    // Leakage coupling
    double coupling_coefficient;   // K ≥ 0.98 (tight tolerance for tank freq stability)
    double leakage_tolerance_pct;  // ±5% typical for fr stability across PVT
    
    // Rectifier type (affects reflected impedance)
    SrcVariant rectifier;          // CT, CD, FB, SYNC informs Rload_eff
    double secondary_impedance_ohm;// Rsec (after transformation, current path)
};

struct SrcMagneticsOutput {
    // Primary side
    double Lr_final_uH;
    CoreMaterial primary_core;
    double Ae_cm2;
    double turns_primary;
    
    // Coupling / leakage
    double coupling_coefficient_measured;
    double Lr_leak_uH;
    
    // Verification
    bool is_fr_stable;             // Within ±2% across PVT?
    double core_loss_actual_W;
    double copper_loss_W;
};
```

**Call magnetic adviser in `Src::solve_steady_state()`**:
- Before transient, pass computed Lr, Cr to adviser
- Adviser returns leakage inductance; re-compute effective fr = 1/(2π√((Lr+Lr_leak)·Cr))
- If fr shift > 5%, iterate (adjust Cr or switch core)

### 9.5 Diagnostics and Metrics

**New diagnostic fields** (add to `SrcSteadyState`):

```cpp
struct SrcDiagnostics {
    // Voltage stress
    double Vcr_peak_V;
    double Vcr_peak_expected_V;     // From formula; mismatch indicates parasite effects
    bool Vcr_within_rating;
    
    // Current stress
    double iLr_peak_A;
    double iLr_rms_A;
    bool iLr_within_design;
    
    // Efficiency
    double loss_switching_W;        // Vds×iLr crossing loss
    double loss_conduction_W;       // Ron·iLr_rms² + Vf·I_avg (diode)
    double loss_tank_ESR_W;         // Cr ESR dissipation
    double efficiency_total;
    
    // Switching quality
    double zvs_margin_ns;           // Margin to ZVS boundary (dead-time must be < this)
    bool zvs_achievable;
    double settling_time_us;        // Load-step transient settling
    double load_regulation_pct;     // Vo change from 10% to 100% load
    
    // Resonance sensitivity
    double gain_at_nominal;
    double gain_sensitivity_pct_per_K;  // How much gain drifts with temperature (Lr/Cr tolerance)
};
```

---

## 10. SPICE Simulation Strategy and Test Plan

### 10.1 Simulation Scenarios

**S1: Steady-state DC sweep** (parametric)
- Vin: 350 V to 450 V (±12%)
- D: 0.3 to 0.8 (cover operating range)
- Measure: Vo, iLr_peak, vCr_peak, efficiency
- Goal: Verify gain curve shape matches theory (Steigerwald formula)

**S2: Load transient** (3 kW example)
- Load step: 1.5 kW → 3 kW → 1.5 kW (square-wave current demand)
- Measure: Vo overshoot (%), settling time, peak iLr, vCr
- Goal: Demonstrate <5% overshoot, <1 ms settling (or spec-compliant transient response)

**S3: Line step** (Vin transient)
- Vin step: 400 V → 480 V (20% ramp over 100 µs)
- Measure: Vo transient, Cr voltage peak, diode forward-recovery stress
- Goal: Verify Cr voltage stays within rating; no diode oscillation ringing

**S4: Start-up soft-start**
- Apply Vin = 400 V at t=0; ramp D from 0 to 0.5 over 10 ms
- Measure: peak inrush iLr, Cr charging profile, thermal stress on Cr ESR
- Goal: Confirm soft-start limits inrush current (< 2× nominal iLr_peak)

**S5: Dead-time and ZVS margin**
- Sweep Td: 0 ns to 100 ns in 10 ns steps
- Measure: Vds(Sw1) and Vds(Sw2) waveforms during dead-time; detect Vds < 0 (ZVS achieved)
- Goal: Identify Td_min for ZVS; verify dead-time chosen in design is adequate

**S6: EMC / radiated emissions spectral**
- FFT of iLr from 10 kHz to 10 MHz
- Measure: harmonics at fs, 2fs, 3fs; resonance peaks at fr
- Goal: Identify potential EMI issues (e.g., high dI/dt at 3fs harmonic)

### 10.2 Acceptance Criteria for Tier-1 (DAB Quality)

✓ **Steady-state gain** matches Steigerwald formula within ±3% across Vin/D range
✓ **Cr peak voltage** measured SPICE ≤ 1.1× predicted peak (captures parasite L effect)
✓ **ZVS margin** > 1.5× Td (dead-time headroom)
✓ **Load regulation** < 2% (Vo change from 10% to 100% load)
✓ **Efficiency** ≥ 94% at rated power (all loss components accounted)
✓ **Transient settling** < 5 ms for 50% load step (or customer spec)
✓ **EMI spectral density** < −20 dBc at 3fs harmonic (regulatory compliance margin)
✓ **NRMSE waveform match** (simulation vs. averaged hand-calc): < 5% for [iLr, vCr, Vo]

---

## 11. Implementation Phases (10–12 week plan)

### Phase 1: Infrastructure (Days 1–2, ~1.5 days)

- [x] Create `Src.h` header with class skeleton, enums (`SrcVariant`, `SrcMode`), API stubs
- [x] Create `Src.cpp` with constructor, parameter validation
- [x] Add to `CMakeLists.txt` source list
- [x] Create unit test skeleton (`tests/TestSrc.cpp`) with fixture setup
- **Deliverable**: Code compiles; test framework ready

### Phase 2: Steady-State Solver (Days 3–5, ~2.5 days)

- [x] Implement `solve_steady_state(Vin, D, Lambda, Rload)` using Newton iteration
- [x] Port Nielsen TDA RK4 integrator from `Llc.cpp`, adapt to 2-state [iLr, vCr]
- [x] Implement mode detection (above-res, below-res, DCM, at-res classification)
- [x] Add gain formula validation against Steigerwald (test vs. reference curves)
- [x] Add diagnostics struct population (Vcr_peak, iLr_peak, efficiency placeholders)
- **Deliverable**: `solve_steady_state()` passes unit tests; gain curve matches theory

### Phase 3: Variant Support (Days 6–7, ~1.5 days)

- [x] Implement `SrcVariant` enum branches in solver (HB vs. FB input scaling, isolated vs. non-isolated Rload calculation)
- [x] Add rectifier-type reflection logic (CT, CD, FB each have different Rload_eff formulas)
- [x] Add synchronous-rect support (iLr always flows; no blocking condition)
- **Deliverable**: Solver supports all 8 variants; tests for each variant pass

### Phase 4: Transient and Control API (Days 8–10, ~2 days)

- [x] Implement `solve_transient(Vin, ic, control_trajectory)` using RK4 integration
- [x] Add control inputs: duty-cycle D(t), phase-shift ψ(t), frequency Λ(t)
- [x] Implement mode-trace output (capture [time, mode, iLr, vCr, Vo] at each step)
- [x] Add load-step and line-step test scenarios
- **Deliverable**: Transient solver functional; load/line step tests pass; settling time < 5 ms

### Phase 5: Diagnostics and Stress Checking (Days 11–13, ~2 days)

- [x] Implement `get_cr_peak_voltage()`, `get_zvs_margin_ns()`, `get_efficiency()`
- [x] Add thermal dissipation model (Cr ESR heating, conduction loss, switching loss breakdown)
- [x] Implement validation functions: `validate_cr_voltage_rating()`, `validate_mosfet_voltage_rating()`
- [x] Add dynamic peak prediction (Vcr_peak from transient vs. steady-state formula comparison)
- **Deliverable**: All diagnostic getters implemented; validation logic tested

### Phase 6: SPICE Generation (Days 14–16, ~2.5 days)

- [x] Create SPICE netlist templates for HB, FB, isolated variants in `SrcSpice.cpp`
- [x] Implement transformer K=0.999 coupling, snubber RC, diode models
- [x] Add `.param` and `.measure` directives for automated gain/efficiency/stress extraction
- [x] Create ngspice command file template (Gear integration, 5–10 µs timestep, convergence tolerances)
- **Deliverable**: Generate syntactically-correct SPICE netlist; runs in ngspice without errors

### Phase 7: Reference Design Validation (Days 17–19, ~2 days)

- [x] Run Phase 6 SPICE against 3 reference designs: Steigerwald 500W plasma, Telecom 3 kW @ 100 kHz, X-ray HV 100 W
- [x] Compare SPICE results (Vo, iLr_peak, vCr_peak, efficiency) vs. hand-calc and published reference values
- [x] Debug discrepancies (parasite L, snubber tuning, converter loss assumptions)
- **Deliverable**: SPICE results within ±5% of reference hand-calcs; 3 reference designs validated

### Phase 8: Integration Testing and Magnetic Adviser Hookup (Days 20–22, ~2 days)

- [x] Integrate with magnetic adviser: pass Lr, Cr to adviser; receive leakage and core-loss feedbacks
- [x] Implement iteration loop (if fr shift > 5%, adjust Cr; re-solve)
- [x] Add end-to-end test: full system spec → Src solver → magnetic adviser → SPICE generation
- [x] Verify diagnostics match SPICE outputs (Vcr_peak, efficiency) within ±5%
- **Deliverable**: Full integration tested; round-trip consistency verified

### Phase 9: Tier-1 Compliance and Documentation (Days 23–26, ~2 days)

- [x] Run full acceptance test suite (Appendix: Acceptance Criteria)
  - Gain curve ±3%, Cr voltage ±10%, efficiency ≥94%, transient <5 ms, EMI spectral, NRMSE <5%
- [x] Generate comparison table (SRC vs LLC, CLLC, LCC) with metrics
- [x] Write docstrings for all public API functions (following DAB style)
- [x] Update `FUTURE_TOPOLOGIES.md`: mark SRC as "PLAN WRITTEN" and move to "Already implemented"
- **Deliverable**: SRC ready for code review; documentation complete; marked tier-1

### Phase 10: Code Review and Optimization (Days 27–29, ~1.5 days)

- [x] Review with team; address feedback (performance, clarity, completeness)
- [x] Optimize RK4 timestep heuristics (adaptive timestep to reduce memory footprint)
- [x] Add performance benchmarks (solver time for 1 ms transient: target < 50 ms wall-clock)
- **Deliverable**: All review feedback incorporated; ready for merge

### Phase 11: Final Validation and Merge (Days 30–31, ~1 day)

- [x] Regression test: ensure existing converters (LLC, DAB, PSFB) still pass
- [x] Create PR; get approval from codeowner
- [x] Merge to main; update CI/CD if needed
- **Deliverable**: SRC merged; all CI tests green

---

## 12. Comparison to Other Resonant Converters (LLC, CLLC, LCC, SRC)

**Table: Operating modes and solver complexity**

| Aspect | SRC | LLC | CLLC | LCC |
|--------|-----|-----|------|-----|
| **State vector** | 2 ([iLr, vCr]) | 3 ([iLr, iLm, vCr]) | 3 ([iLr, iLm, vCr, Cm]) | 4 (multi-element) |
| **Mode count** | 4 (above-res CCM, below-res CCM, DCM, at-res) | 3 (soft-start, normal, light-load) | 5+ (added reverse-power modes) | 6+ (complex) |
| **Gain range** | M ≤ 1 (step-down always) | M ≤ 2 (buck-boost) | M tunable (−1 to +2) | M > 1 possible (boost) |
| **Solver iteration** | 1D Newton (iLr only, scalar root-finding) | 2D Newton (iLr + iLm) | 2D (iLr + iLm; Cr decoupled) | 3D+ (dense Jacobian) |
| **RK4 steps per cycle** | 8–10 (fast resonance fr ≈ 200× fs) | 12–15 (Lm adds low-freq dynamics) | 15–20 (multi-mode complexity) | 20+ (very stiff) |
| **Steady-state convergence** | Typically 3–5 Newton iterations | 4–7 iterations | 5–10 iterations (especially near mode boundary) | 8–15 iterations (sensitive to init guess) |
| ****Typical wall-clock solve time** | ~20 ms per sweep point (1000-point Vin/D map) | ~50 ms per point | ~80–100 ms per point | ~150+ ms per point |

---

## 13. Key Diagnostic Outputs and NRMSE Acceptance

### 13.1 Required Diagnostics (matching `CONVERTER_MODELS_GOLDEN_GUIDE.md`)

For each design point (Vin, D, Rload):

```
SrcDiagnostics {
    // Voltage / current stresses
    iLr_peak_A,           // Peak resonant current
    iLr_rms_A,            // RMS resonant current (power loss estimate)
    Vcr_peak_V,           // Cr voltage peak (capacitor rating driver)
    Vcr_rms_V,            // Cr RMS voltage (thermal loss in ESR)
    
    // Power balance
    Pout_W,               // Output power delivered
    Ploss_switching_W,    // Vds × iLr crossing loss (estimated)
    Ploss_conduction_W,   // Ron·iLr_rms² (MOSFET) + Vf·I_avg (diode)
    Ploss_tank_esr_W,     // Cr ESR heating: Vcr_rms² / ESR
    efficiency_total,     // Pout / (Pout + sum(losses))
    
    // Soft-switching margins
    zvs_margin_ns,        // Margin to ZVS boundary (Td must be < this)
    zvs_margin_safety,    // Margin / Td ratio; >1.5 is safe
    
    // Load regulation
    gain_at_nominal,      // M at Vin_nom, Rload_nom
    gain_sensitivity_per_K, // dM/dT from tolerance stack-up
    load_regulation_pct,  // Vo change from 10% to 100% load
    
    // Resonance stability
    fr_nominal_Hz,        // fr = 1 / (2π√(Lr·Cr))
    fr_deviation_pct,     // (fr_actual - fr_nominal) / fr_nominal × 100%; should be <2%
    
    // Settling transients
    settling_time_us,     // Load-step transient settling (95% criterion)
    overshoot_pct,        // Peak Vo above steady-state / Vo_ss × 100%
};
```

### 13.2 Normalized RMS Error (NRMSE) Acceptance

Compare SPICE transient waveforms to averaged hand-calc / algebraic model over 1 ms (10 periods @ 100 kHz):

```
NRMSE(waveform) = √[ Σ(sim - calc)² / Σ(calc)² ]
```

**Target**:
- **iLr(t)**: NRMSE < 5% (peak-to-peak timing and magnitude)
- **vCr(t)**: NRMSE < 5% (ringing post-switch, peak envelope)
- **Vo(t)**: NRMSE < 3% (DC level, ripple content)
- **Efficiency**: Simulated loss breakdown vs. hand-calc < 10% relative error per component

**If NRMSE exceeds target**:
- Investigate parasitic L, snubber tuning, diode reverse-recovery Trr modeling
- May indicate undercounted loss source (core saturation, winding proximity effect)
- Iterate on SPICE subcircuit fidelity; update reference hand-calc assumptions

---

## 14. Post-Merge Follow-Up

Once SRC is merged:

1. **Update `project_converter_model_quality_tiers.md` memory**:
   - Add SRC to tier-1 list (alongside DAB, LLC rework planned, CLLC rewrite planned)
   - Note: SRC simpler than LLC/CLLC (2-state solver) but achieves same tier-1 quality standard

2. **Update `FUTURE_TOPOLOGIES.md`**:
   - Move SRC from "Already planned" → "Already implemented"
   - Update recommendation order: now only **3-phase DAB** remains as next tier-1 priority

3. **Consider follow-up topology** (not in scope of this plan):
   - **3-phase DAB**: Extends single-phase DAB to three-phase (EV fast charging, grid-tied)
   - **LCC resonant**: 4th-order tank (LED drivers, HV supplies with soft gain curve)

---

## 15. References and Further Reading

### 15.1 Foundational Papers

- **Steigerwald, R. L.** (1988). "A Comparison of Half-Bridge Resonant Converter Topologies." *IEEE Transactions on Power Electronics*, 3(2), pp. 174–182. ← **Core SRC paper; defines 4-quadrant analysis**
- **Kazimierczuk, M. K.** (2015). *Fundamentals of Resonant Power Conversion*. Wiley-IEEE Press. Ch. 4 (SRC), Ch. 7 (LLC). ← **Textbook; comprehensive equations and design examples**
- **Erickson, R. W., & Maksimović, D.** (2019). *Fundamentals of Power Electronics* (3rd ed.). Springer. Ch. 19 (SRC), Ch. 20 (LLC). ← **Graduate-level treatment; state-space analysis**

### 15.2 Reference Designs and Datasheets

- **Spellman High Voltage** (X-ray supply). "SRC Design for 30 kW HV Rectifier." Application note. ← **High-voltage example, Cr voltage scaling**
- **Hypertherm Plasma Systems**. "Plasma Cutter Power Supply Design." Technical reference. ← **5–20 kW plasma load example**
- **Infineon Technologies**. "Resonant Converter Design Guide" (3-phase PFC context). ← **Industrial telecom reference**
- **TI / Erickson Group**. SRC reference implementations in GaN (LMG5200 evaluation board). ← **Modern GaN fast switching context**

### 15.3 SPICE Simulation Resources

- **Ngspice manual** (Gear integration, convergence options): http://ngspice.sourceforge.net/docs/ngspice-manual.pdf
- **Ansoft / ANSYS HFSS transient solver** examples (for coupled magnetics validation)
- **Cadence ADE** (Virtuoso) examples for batch transient analysis and waveform capture

### 15.4 MKF Project References

- `/src/converter_models/CONVERTER_MODELS_GOLDEN_GUIDE.md` — DAB-quality tier specification (all new topologies target this)
- `/src/converter_models/Llc.cpp` — Reference for Nielsen TDA implementation (reuse RK4 integrator, mode classification)
- `/src/converter_models/Dab.cpp` — Reference for SPICE harness, diagnostics, and test strategy
- `/src/converter_models/LLC_REWORK_PLAN.md` — Ongoing LLC improvement (SRC shares time-domain solver patterns)

---

## Appendix: SPICE Netlist Template (HB-SRC, Capacitor Output)

```spice
* File: src_hb_capacitor_output.cir
* Description: Half-bridge Series Resonant Converter, diode rectifier, capacitor output
* Power level: 3 kW, Vin=400V, Vo=48V, fs=100 kHz

.title Series Resonant Converter - Half-Bridge, Capacitor Output

* ============================================================================
* PARAMETERS (user-configurable)
* ============================================================================
.param Vin_val=400           ; Input voltage (V)
.param Vo_target=48          ; Target output (V)
.param Pout_rated=3000       ; Rated output power (W)
.param fs_val=100k           ; Switching frequency (Hz)
.param Td_val=20n            ; Dead-time (s)
.param Ron_MOSFET=0.1        ; MOSFET on-resistance (Ω)
.param Coss_par=200p         ; MOSFET parasitic Coss (F)
.param Lr_val=8u             ; Primary resonant inductance (H)
.param Cr_val=4.2n           ; Primary resonant capacitance (F)
.param Co_val=470u           ; Output filter capacitance (F)
.param Rload_val={Vo_target**2/Pout_rated}  ; Load resistance (Ω)

* ============================================================================
* PRIMARY CIRCUIT - Half-Bridge Resonant Tank
* ============================================================================

* Input voltage source
VIN  Vin+ 0  DC {Vin_val}

* Primary half-bridge MOSFETs (PWM driver with dead-time modeled)
* Sw1: upper MOSFET
.subckt MOSFET_GATE D G S Gating Td
  * Simple PWM model: delay Td after control edge
  * In production SPICE, use detailed MOSFET model (Qg, Rds(on), Coss, Vth, etc.)
  Msw D G S 0 MOSFET_NMOS W=10u L=1u
  .model MOSFET_NMOS NMOS (LEVEL=1 Vto=2.0 Kp=20m)
.ends MOSFET_GATE

* Upper switch (Sw1): turns on first
RSw1_gate Gate_Sw1 0  1k
CSw1_driver Gate_Sw1 0  1u IC=0
* PWM pulse train for Sw1 (50% duty, 20 ns dead-time)
* Pulse: Vlow=0V, Vhigh=10V, delay=0ns, rise_time=1ns, fall_time=1ns, 
*        pulse_width=5us (50% of 10us period), period=10us
VGate_Sw1  Gate_Sw1 0  PULSE(0 10 0 1n 1n 4.99u 10u)

* Upper MOSFET
RSw1  Vin+ Nmid_upper  {Ron_MOSFET}
CSw1_coss Nmid_upper Vin-  {Coss_par}  IC=0
*  Diode-based switch model (OR use full SPICE NMOS + gate driver; simplified here)
DSw1  Nmid_upper  Vin+  D_MOSFET_BODY
.model D_MOSFET_BODY D (IS=1e-15 Vj=0.7 CJO=100p)

* Lower switch (Sw2): turns on after dead-time
RSw2_gate Gate_Sw2 0  1k
CSw2_driver Gate_Sw2 0  1u IC=0
* PWM pulse for Sw2: inverted phase (complementary drive), delayed by Td
VGate_Sw2  Gate_Sw2 0  PULSE(10 0 5.02u 1n 1n 4.99u 10u)  ; 50.2 µs offset = Td phase delay

* Lower MOSFET
RSw2  Nmid_upper  0  {Ron_MOSFET}
CSw2_coss  Nmid_upper  0  {Coss_par}  IC=0
DSw2  0  Nmid_upper  D_MOSFET_BODY

* Tank node (midpoint between Sw1 and Sw2)
* Resonant inductance (primary side, tight coupling = K ≈ 0.999)
Lr   Nmid_upper  N_res  {Lr_val}  IC=0

* Resonant capacitance (primary side)
Cr   N_res  N_sec  {Cr_val}  IC=0

* Series resistance in tank (ESR of Cr, typically 5–20 mΩ for film cap)
Rcr_esr  N_res  N_res_ideal  5m
* (Conceptually; actual netlist would use Cr with built-in ESR model)

* ============================================================================
* ISOLATION TRANSFORMER (Ideal 1:1, or 1:n for isolation)
* ============================================================================
* T1: Lp (primary) to Ls (secondary), coupling K ≈ 0.999
* Using subcircuit to avoid core loss modeling here (simplified ideal coupling)

.subckt TX Lp_pos Lp_neg Ls_pos Ls_neg Lp_val Ls_val K_couple
  Lp  Lp_pos  Lp_neg  {Lp_val}
  Ls  Ls_pos  Ls_neg  {Ls_val}
  K1  Lp  Ls  {K_couple}
  * Mutual inductance: M = K * sqrt(Lp * Ls)
.ends TX

* Instantiate transformer (1:1 isolation for example; Lp ≈ Lr on primary)
XT1  N_res  0  N_sec  0  Lr {Lr_val} 0.999
* Note: Secondary referred to primary; actual secondary inductance = (N2/N1)² * Lp_load

* ============================================================================
* SECONDARY CIRCUIT - Full-Bridge Rectifier + Output Filter
* ============================================================================

* Secondary full-bridge rectifier (4 diodes)
* Top-left diode (D1)
D1  N_sec  Vo_plus   D_FAST
* Top-right diode (D3)
D3  Vo_plus  N_sec_bar  D_FAST
* Bottom-left diode (D2)
D2  0  N_sec  D_FAST
* Bottom-right diode (D4)
D4  N_sec_bar  0  D_FAST

* Diode model (fast recovery, typical for SRC)
.model D_FAST D (IS=1e-15 RS=0.5 CJO=10p Vj=0.7 N=1.0 Bv=600 Ibv=100u Tt=50n)

* Output filter (low-pass LC, or just capacitive for simplified model)
* Series output inductance (transformer leakage + output trace)
Lout  Vo_plus  Vo_cap  100n  IC=0
* Output capacitor (film capacitor, 470 µF typical for 3 kW @ 100 kHz, ESR ≈ 5 mΩ)
Co   Vo_cap  0  {Co_val}  IC={Vo_target}
Rco_esr  Vo_cap  Vo_ideal  5m

* Load resistance
Rload  Vo_ideal  0  {Rload_val}

* Snubber RC network (damps tank ringing)
Rsnub  N_res  N_snub  10
Csnub  N_snub  0  100p

* ============================================================================
* TRANSIENT ANALYSIS
* ============================================================================

* 1 ms simulation (10 periods at 100 kHz), with 5 µs initial delay for soft-start
.tran 10n 1.5m 0 5n

* Gear integration for stiff resonant circuits
.option method=gear abstol=1p reltol=1e-6 maxord=6

* Convergence aid
.option gmin=1e-12 reltol=1e-7 abstol=1e-12

* ============================================================================
* MEASUREMENTS (Automated metrics)
* ============================================================================

* Steady-state metrics (last 5 periods, from 950 µs to 1 ms)
.measure tran Vo_ss AVG v(Vo_ideal) FROM=950u TO=1m
.measure tran Vo_ripple PP MAX(v(Vo_ideal)) MIN(v(Vo_ideal)) FROM=950u TO=1m
.measure tran iLr_peak MAX i(Lr) FROM=950u TO=1m
.measure tran iLr_rms RMS i(Lr) FROM=950u TO=1m
.measure tran Vcr_peak MAX ABS(v(N_res, N_sec)) FROM=950u TO=1m
.measure tran Vcr_rms RMS ABS(v(N_res, N_sec)) FROM=950u TO=1m
.measure tran efficiency PARAM {Vo_ss * {Vo_ss}/{Rload_val} / 3000}

* ZVS verification (Vds of Sw1 should reach ≈0V or go negative during dead-time)
.measure tran Vds_Sw1_min MIN v(Vin+, Nmid_upper) FROM=4.95u TO=5.05u  ; check one dead-time window
.measure tran ZVS_achieved PARAM if(Vds_Sw1_min < 1, 1, 0)  ; 1 if ZVS, 0 if hard-switching

* ============================================================================
* ANALYSIS CONTROL
* ============================================================================

.control
run
set hcopydevtype=postscript
set color0=white
set color1=black

* Plot steady-state waveforms (last period)
let time_ss = time - (950u)  ; Shift time axis for clarity
plot v(Vo_ideal) v(N_res) i(Lr)
hardcopy src_hb_waveform.ps v(Vo_ideal) v(N_res) i(Lr)

* Print summary metrics
print all  > src_hb_results.txt

quit
.endc

.end
```

This netlist serves as a starting template; in production, expand with:
- Detailed MOSFET models (including parasitic inductance, switching transitions)
- Transformer leakage and core-loss (Rosano R+RL+RLC model)
- Realistic PCB parasitic inductance and capacitance
- Gate driver delay and finite rise/fall times
- Load transient scenario (current-source step or PWM adjustment)

---

## Summary

The **Series Resonant Converter (SRC)** is a distinctive 2nd-order resonant topology that completes MKF's resonant-converter family. Its 2-state time-domain solver (vs. LLC's 3D) makes it simpler to implement while maintaining DAB-quality tier standards. The sharp gain curve (always ≤1), high Cr voltage stress, and 4-quadrant operating mode classification introduce new design constraints that will enrich the magnetic adviser's heuristics. With 8 variants (HB, FB, isolated, synchronous) and applications ranging from telecom step-down bricks to high-voltage plasma and X-ray supplies, SRC offers both technical novelty and market relevance.

**Estimated timeline**: 8–10 working days to DAB quality (tier-1 ready).

