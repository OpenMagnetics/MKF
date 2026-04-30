# Weinberg Converter — Topology Add Plan

> Purpose: add a new `Weinberg` topology to `src/converter_models/`
> (none exists today), wired into the MAS schema, with full ngspice
> harness and tests, meeting **DAB-quality** as defined by
> `CONVERTER_MODELS_GOLDEN_GUIDE.md`.
>
> Read `CONVERTER_MODELS_GOLDEN_GUIDE.md` first; this plan layers on top.
> Sibling references:
> - `Llc.{h,cpp}` — closest 2-magnetic structure (Lr + xfmr); reuse the
>   "main magnetic + ancillary" pattern
> - `Flyback.{h,cpp}` — closest center-tap-secondary rectifier pattern
> - `IsolatedBuckBoost.{h,cpp}` — multi-magnetic DesignRequirements
> - `PushPull.{h,cpp}` — closest **drive scheme** (push-pull, two
>   switches, 180° phase, flux-walk concerns)
> - `ActiveClampForward.{h,cpp}` — closest **multi-element extras** pattern
> - `CONVERTER_MODELS_GOLDEN_GUIDE.md`, `CLLC_REWRITE_PLAN.md`, `CUK_PLAN.md`
>
> Weinberg is a **non-resonant, isolated, current-fed push-pull-derivative
> PWM converter** with mandatory switch overlap above D = 0.5. Solver is
> analytical piecewise-linear (no Newton, no matrix exponential). The
> hard parts are: (a) the dual-magnetic structure (input 1:1 coupled L
> + main 4-winding center-tapped transformer); (b) the three operating
> regimes (D < 0.5 / D = 0.5 / D > 0.5); (c) energy-recovery diode
> modeling; (d) flux-walk safety in the main transformer.

---

## 1. Why Weinberg in MKF

- **Spacecraft heritage.** ESA Cluster, ENVISAT, Galileo GSTB-V2/A,
  reportedly NASA ISS DC-DC unit. The space-power user community has
  been a recurring ask in MKF issues (hi-rel, radiation-tolerant
  magnetic design).
- **Continuous input AND output current** with isolation — unique
  combination in MKF's catalogue. Push-pull has pulsating input;
  full-bridge has pulsating output until you add a big Lo; isolated
  Cuk has continuous I/O but a stressed coupling cap.
- **Boost-capable isolated topology.** MKF currently has Flyback and
  IsolatedBuckBoost as boost-capable isolated converters; Weinberg
  adds the high-power (0.5–5 kW) regime where Flyback fails on
  energy-storage limits.
- **Exercises a magnetic-design problem MKF lacks**: the input 1:1
  coupled inductor (closely coupled, **gapped**, equal DC bias on
  both windings) — neither a pure transformer (it stores energy) nor
  a pure inductor (its 1:1 coupling is load-bearing).
- **Disambiguation value.** "Weinberg" is often confused with
  Watkins-Johnson and plain current-fed push-pull. Adding it cleanly
  documents the differences (Guide §3 disambiguation requirement).

---

## 2. Variants in scope

From the research dossier (Weinberg ESA/ESTEC 1974; Weinberg &
Schreuders IEEE PESC ~1985; Larico-Barbi IEEE TIE 2012; Yadav-Gowrishankara
IEEE 2016; Bhandari et al. J. Power Electron. 2025; Patent CN112821760B;
Mazumder UIC SiC paper):

1. **Classic Weinberg (V1)** — 1:1 input coupled inductor + push-pull
   primary (center-tapped, 2 switches) + main isolation transformer +
   center-tapped (CT) full-wave-rectifier secondary + energy-recovery
   diode D3 to Vin. Unidirectional. **Baseline.**
2. **Bridge Weinberg (V2)** — replace push-pull with H-bridge primary
   (4 switches). Halves switch V_DS stress (Vin instead of 2·Vin/(1−D));
   used at high-bus-voltage spacecraft (300–600 V).
3. **Bidirectional Weinberg (V3)** — active rectifier on the secondary
   (4 MOSFETs replacing CT diodes); runs as charger one way,
   discharger the other. Modern Li-ion satellite buses.
4. **Multi-output Weinberg (V4)** — multiple secondaries on T (each
   with its own rectifier). Cluster/ENVISAT-style multi-bus.
5. **ZVS/ZCS soft-switched Weinberg (V5)** — active-clamp snubber adds
   resonant turn-on/off. Patent CN112821760B; IEEE 8397710. Implemented
   as a flag.
6. **FB-rectifier secondary (V6)** — secondary full-bridge instead of
   CT full-wave. Rectifier-type switch like PSFB / Flyback.

**Out of scope for the initial add** (defer to follow-ups):

- Three-phase / interleaved Weinberg (Larico-Barbi 2012) — separate
  topology, materially different magnetics.
- All-SiC parallel Weinberg with current sharing — control-system
  problem, not a magnetic-model problem.
- Non-isolated Weinberg variant (1:1 main "transformer" used purely
  as coupled inductor).
- DCM design (V1 enforces strict CCM in `process_design_requirements`;
  DCM detected and flagged but not designed for — same approach as
  the Cuk plan).

---

## 3. Reference designs as test fixtures

| # | Source | Vin | Vout | Pout | fsw | D | n (Np/Ns) | L1 | η | Use |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Weinberg & Schreuders, IEEE PESC ~1985 | ~50 V | ~150–300 V | 1.5 kW | 30–50 kHz | 0.6–0.7 | 1:3 | – | >93 % | foundational reference, V1 baseline |
| 2 | Galileo GSTB-V2/A BDR (ESA SP-589, 2005) | 30–50 V (battery) | ~50 V (bus) | several kW | mid-kHz | – | – | – | >96 % | space-heritage reference (numerical specs scarce) |
| 3 | Larico & Barbi, IEEE TIE 59(2):888–896, 2012 (3-phase Weinberg) | 120 V | 75 V | 735 W | 42 kHz | 0.33 / 0.66 | 1:1 | – | high | proves equations on a published 3-phase build (out of initial scope) |
| 4 | Yadav & Gowrishankara, IEEE 2016 (fault-tolerant) | 42 V | 100 V | 500 W | 100 kHz | 0.7 | 1:1.4 | – | – | mid-power CCM smoke test |
| 5 | Bhandari et al., J. Power Electron. 25:591, 2025 (low-ripple BDR) | 50 V | 100 V | 5 kW BDR module | – | – | – | – | 95–97 % | high-power reference, recent |
| 6 | "Practical BDR with Weinberg," IJRTE 9(1), 2020 | 24–32 V | 50 V | – | 100 kHz | 0.6–0.7 | 1:2 | – | – | low-power smoke test |
| 7 | Mazumder et al., All-SiC parallel Weinberg unit | 100 V | 270 V | 5 kW (parallel) | 100 kHz (SiC) | – | – | – | high | high-fsw SiC reference |
| 8 | NASA NTRS 20020083039 (ISS Weinberg) | 120 V bus | 28 V | several hundred W | – | – | – | – | – | space-heritage reference |

URLs (for citation in header docstring):

- Weinberg ESA SP-589 (2005): https://ui.adsabs.harvard.edu/abs/2005ESASP.589E..11W
- Weinberg & Schreuders IEEE PESC: https://ieeexplore.ieee.org/document/254756/
- Weinberg & Schreuders, semantic scholar mirror: https://www.semanticscholar.org/paper/A-High-Power-High-Voltage-DC-DC-Converter-for-Space-Weinberg-Schreuders/a802d90fc951753245dc892c118b7bad97c2aa40
- Larico & Barbi IEEE TIE 2012: https://ieeexplore.ieee.org/document/5764531/
- Yadav & Gowrishankara IEEE 2016: https://ieeexplore.ieee.org/document/7516476/
- Bhandari et al. J. Power Electron. 2025: https://link.springer.com/article/10.1007/s43236-025-00989-4
- IJRTE 9(1) 2020 Practical BDR Weinberg: https://www.ijrte.org/wp-content/uploads/papers/v9i1/F1124038620.pdf
- Patent CN112821760B (ZVS-ZCS Weinberg): https://patents.google.com/patent/CN112821760B/en
- Mazumder UIC SiC parallel Weinberg: https://mazumder.lab.uic.edu/wp-content/uploads/sites/504/2021/03/27.pdf
- Marsala & Olivieri (Weinberg vs current-fed push-pull): https://ieeexplore.ieee.org/document/573330/
- ZVS-ZCS DC-DC Weinberg IEEE 8397710: https://ieeexplore.ieee.org/document/8397710
- Modernized Weinberg charge-discharge IEEE 8434944: https://ieeexplore.ieee.org/document/8434944
- Improved low-ripple Weinberg IEEE 8927072: https://ieeexplore.ieee.org/document/8927072
- NASA NTRS ISS Weinberg: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20020083039.pdf
- Watkins-Johnson disambiguation (ESPC 2017): https://www.e3s-conferences.org/articles/e3sconf/pdf/2017/04/e3sconf_espc2017_14004.pdf

---

## 4. Equations (canonical block for the new header docstring)

Notation: `D` = duty of each switch (each switch is ON for `D·Tsw`
within its half-period, with 180° phase between Q1 and Q2);
`(2D−1)` = overlap fraction (both switches simultaneously ON when
`D > 0.5`); `n = Np/Ns` (each half of CT primary vs each half of CT
secondary); `Tsw = 1/fsw`; `R = Vout/Iout`; subscripts `L1` for input
coupled inductor, `Lm` for main transformer magnetizing inductance,
`Llk` for main transformer leakage.

### 4.1 Three operating regimes

| Regime | D range | Switch state | M = Vout/Vin (ideal, lossless) |
|---|---|---|---|
| **Buck-like** | 0 < D < 0.5 | non-overlap | `M = D / n` |
| **Boundary** | D = 0.5 | continuous alternation, no overlap | `M = 1 / n` |
| **Boost-like** | 0.5 < D < 1 | overlap fraction `(2D−1)` | `M = 1 / (2·n·(1−D))` |

The two regimes meet at `D = 0.5, M = 1/n`. Most BDR designs operate
in the **boost regime** (D nominal 0.6–0.7).

### 4.2 Steady-state equations (CCM, 1:1 input coupled-L, n = Np/Ns)

```
(E1)  M (boost, D > 0.5)        = Vout/Vin = 1 / (2·n·(1 − D))
(E2)  M (buck,  D < 0.5)        = Vout/Vin = D / n
(E3)  Volt-sec balance (boost)  : Vin·D·Tsw = n·Vout·(1−D)·Tsw  (per primary half)
(E4)  Iin_avg                    = Iout / (η · M)
(E5)  Switch peak voltage       : V_Q,pk = 2·Vin / (1 − D)  ≈ 2·n·Vout
                                  (clamped by D3 in practice)
(E6)  Switch peak current       : I_Q,pk = Iin + ΔI_L1/2     (non-overlap)
                                          = Iin/2 + ΔI_L1/2  (during overlap, shared)
(E7)  Switch RMS current        : I_Q,rms ≈ Iin · √D  (approx; depends on overlap)
(E8)  Diode peak rev voltage    : V_D,pk = 2·Vout (CT FW rectifier)
(E9)  Diode avg current         : I_D,avg = Iout/2  (per diode, FW rectifier)
(E10) Output cap RMS current    : I_Co,rms ≈ Iout · √((1−D)/D)
                                  (much smaller than push-pull because Iout is continuous)
(E11) Energy-recovery D3 avg I  : I_D3,avg ≈ Iin · (2D−1) · k_leak
                                  (depends on leakage L_lk and clamp action)
(E12) RHP zero (boost)          : ω_RHPZ ≈ R · (1 − D)² / (n² · L1)
                                  Loop bandwidth limit: fc ≤ ω_RHPZ / (2π · 5)
```

### 4.3 Sizing (CCM, boost regime)

```
(S1)  L1   ≥ Vin · (2D − 1) · Tsw / (4 · ΔI_L1)
        keeps ΔI_L1 ≤ 20 % of Iin/2 per winding
(S2)  Lm   ≥ Vin · D · Tsw / (4 · ΔI_mag)        (per primary half)
        sized for ΔI_mag ≤ 10 % of reflected Iout
(S3)  Co   ≥ Iout · D / (ΔVo · fsw)               (continuous output current)
(S4)  n    chosen at D_nom = 0.6–0.7 to satisfy (E1) at nominal Vin/Vout
(S5)  D_max ≤ 0.85 (avoid M → ∞ singularity); D_min ≥ 0.5 + margin if pure-boost regime is required
```

### 4.4 CCM/DCM boundary

The Weinberg enters DCM when the L1 magnetizing current goes to zero
between switching events (light load). In CCM `M` is load-independent;
in DCM it becomes load-dependent (boost-like
`1 + something/(R·Tsw/L1)`).

```
(B1)  K_crit (boost) = 2 · L1 · fsw · (1 − D)² / R
        CCM if K_crit < 1, else DCM
(B2)  DCM gain (boost): M_DCM = 1/(2·n·(1−D)) · f(K, D)   [load-dependent]
```

For BDR application, designers force CCM at all expected loads by
sizing `L1` for `Iout_min`. V1 enforces CCM in design; DCM detected
and flagged in operating-point analysis only.

### 4.5 Flux balance (main transformer)

```
(F1)  Asymmetric Q1/Q2 timing → flux walk in T → saturation
(F2)  Mitigation: peak-current-mode control OR DC-blocking series cap on primary
       OR strict on-time matching (typ. <50 ns)
(F3)  Volt-sec per half-cycle (each primary half):  Vin · D · Tsw - n·Vout·(1−D)·Tsw = 0
       (this is (E3); flux balance is automatic if D matches at both switches)
```

### 4.6 Energy-recovery diode D3

D3 clamps the leakage spike from `Llk` of the main transformer back
to `Vin`. The "new energy-transfer principle" of Weinberg's 1974
paper is that this energy is *recovered* (not dissipated in a
snubber). In SPICE, D3 is modeled as a real diode from the primary
center-tap (or from each switch drain) to `Vin`.

---

## 5. Solver: piecewise-linear analytical (no TDA)

Weinberg is a PWM converter; its sub-intervals are LTI and the
inductor currents are exact-linear within each. Pattern to follow:
`Boost.cpp`, `IsolatedBuckBoost::processOperatingPointsForInputVoltage`,
and the secondary-rectifier shape from `Flyback.cpp`.

### 5.1 Sub-intervals (CCM, V1, boost regime D > 0.5)

Per switching period, the canonical sequence is:

| Sub-interval | Length | Q1 | Q2 | D3 | Notes |
|---|---|---|---|---|---|
| **A** Q1 only | (1−D)·Tsw / 2 | ON | OFF | OFF | Iin flows through L1a → Q1 → Np_a → diode_pos; output gets `Vin/n` |
| **B** Overlap | (2D−1)·Tsw / 2 | ON | ON | OFF | primary shorted; L1 stores energy; output sees only L1's reflected current |
| **C** Q2 only | (1−D)·Tsw / 2 | OFF | ON | OFF | mirror of A; diode_neg conducts |
| **D** Overlap | (2D−1)·Tsw / 2 | ON | ON | OFF | mirror of B |

(For D < 0.5 buck regime, replace overlap with **dead-time** (both
OFF); the canonical sequence becomes A, dead, C, dead.)

Within each sub-interval:

- **iL1a, iL1b** propagate linearly from a known initial condition:
  `iL1(t+Δt) = iL1(t) + (vL1 / L1) · Δt` per winding.
- **iLm** propagates linearly under primary applied voltage (which
  alternates +Vin / 0 / −Vin depending on sub-interval).
- **iout** is the rectified positive of `n · (iL1a + iL1b − iLm)` —
  the secondary sees the difference of input-coupled-L current and
  the magnetizing current.

Steady-state: by symmetry (Q1/Q2 mirror), only half a switching
period needs to be solved. Volt-sec balance gives `M(D)` exactly per
(E1)/(E2). Time grid: `2·N_samples + 1` with `N_samples = 256`,
span `[0, Tsw]` exactly. Per Guide §4.8.

### 5.2 Per-OP diagnostics (mandatory per Guide §2.4 / §4.6)

```cpp
mutable double lastDutyCycle = 0;
mutable int    lastOperatingRegime = 0;     // 0 = buck (D<0.5), 1 = boundary, 2 = boost (D>0.5)
mutable int    lastOperatingMode = 0;       // 0 = CCM, 1 = DCM
mutable double lastConversionRatio = 0;
mutable double lastOverlapFraction = 0;     // (2D-1) clamped to [0,1]
mutable double lastSwitchPeakVoltage = 0;   // 2·Vin/(1-D)
mutable double lastSwitchPeakCurrent = 0;
mutable double lastSwitchRmsCurrent = 0;
mutable double lastDiodePeakReverseVoltage = 0;
mutable double lastDiodeAverageCurrent = 0;
mutable double lastEnergyRecoveryAvgCurrent = 0;  // I_D3,avg
mutable double lastInputInductorRipple = 0;       // ΔI_L1
mutable double lastMagnetizingCurrentRipple = 0;  // ΔI_mag
mutable double lastFluxImbalanceMargin = 0;       // (D_Q1 - D_Q2) / D — should be < 1%
mutable double lastRhpZeroFrequency = 0;
mutable double lastRecommendedLoopBandwidth = 0;  // fRHP/5
```

### 5.3 Multi-magnetic structure

Like LLC (Lr + xfmr) and IsolatedBuckBoost (Lo + xfmr). Default
"primary" magnetic returned by `process_design_requirements` is the
**main center-tapped transformer** (4 windings: Np_a, Np_b, Ns_a,
Ns_b, where Np_a + Np_b = full primary, Ns_a + Ns_b = full secondary).

The **input 1:1 coupled inductor L1** (2 windings, gapped, equal DC
bias) is emitted via `get_extra_components_inputs`. So is the output
filter cap if present. So is D3's path inductance if user models it.

Schema flag `mainMagnetic` = `inputCoupledInductor` | `transformer`
lets users flip which one is the "primary" magnetic to design for.

---

## 6. SPICE netlist (generate_ngspice_circuit)

Push-pull current-fed topologies with overlap are **convergence
gremlins** — both switches ON simultaneously short the primary
through the input coupled-L; without leakage modelling the matrix is
singular. The netlist below is built to the Guide §5 contract, plus
**mandatory leakage inductance** in both magnetics.

### 6.1 Topology — V1 classic Weinberg

```
* Weinberg Converter (classic, CCM, boost regime) — Generated by OpenMagnetics
*
*  Vin─┬───L1a─┬─Np_a─Lleak_p─┬─D3 cathode (clamp to Vin)
*      │  k=1  │              │
*      └───L1b─┴─Np_b─Lleak_p─┘
*               │              │
*              Q1             Q2
*               │              │
*              GND            GND
*
*  Secondary (CT full-wave):
*       Ns_a ──── D_pos ──┬── L_out (optional) ── out_node ── Co ── R_load ── 0
*                          │
*       Ns_b ──── D_neg ──┘
*
* Drive:
*   Q1 gated by pwm1 = PULSE(0,5,0,..., D·Tsw, Tsw)
*   Q2 gated by pwm2 = PULSE(0,5, Tsw/2, ..., D·Tsw, Tsw)
*   When D > 0.5 the two pulses overlap by (2D−1)·Tsw/2 each per period.
```

### 6.2 Mandatory netlist contents (Guide §5 + Weinberg specifics)

- `Vin vin_p 0 {Vin}`; zero-volt `Vin_sense` for Iin extraction.
- **Input coupled inductor** (1:1, two windings on same magnetic):
  ```
  L1a vin_sense node_PriCT_a {L1} ic={Iin/2}     ; IC mandatory
  L1b vin_sense node_PriCT_b {L1} ic={Iin/2}
  K_in L1a L1b 0.999                               ; never 1.0
  ```
- **Main transformer** (4 windings, center-tapped pri AND CT sec):
  ```
  Lpri_a node_PriCT_a node_DrainQ1 {Lpri_half}
  Lpri_b node_PriCT_b node_DrainQ2 {Lpri_half}
  Lsec_a node_SecCT   node_DiodePos {Lsec_half}
  Lsec_b node_SecCT   node_DiodeNeg {Lsec_half}
  K_xfmr Lpri_a Lpri_b Lsec_a Lsec_b 0.99          ; 4-way coupling, leave headroom for leakage
  ```
  Add explicit `Lleak_pri_a`, `Lleak_pri_b` series with each primary
  half (typ 2 % of Lpri_half) — leakage is *intentionally usable*
  for ZCS in modern variants and is essential for SPICE convergence
  during overlap.
- **Switches** (2 active, push-pull):
  ```
  Q1 node_DrainQ1 0 pwm1 0 SW1
  Q2 node_DrainQ2 0 pwm2 0 SW1
  ```
  In V2 (bridge) replace with 4 active switches; in V3 (bidirectional)
  also replace secondary diodes.
- **Energy-recovery diode D3** (the key element):
  ```
  D3a node_DrainQ1 vin_p DIDEAL
  D3b node_DrainQ2 vin_p DIDEAL
  ```
  (Two diodes — one per switch drain — return the leakage spike to
  Vin. Some Weinberg variants use a single D3 from primary center-tap
  to Vin instead; pick by schema flag `energyRecoveryStyle`.)
- **Secondary rectifier** (CT full-wave, V1):
  ```
  D_pos node_DiodePos out_node DIDEAL
  D_neg node_DiodeNeg out_node DIDEAL
  ```
  In V6 (FB rectifier) emit 4 diodes. In V3 (bidirectional) emit 2 or 4
  active SR MOSFETs.
- **Optional Lo + Co** (most Weinberg designs omit Lo because the input
  L1 already smooths the secondary current via reflection; emit only
  if `outputInductance` is set):
  ```
  Lo out_node out_sense {Lo} ic={Iout}
  Co out_sense 0 {Co} ic={Vout}
  ```
- **Load**: `Rload out_sense 0 {max(Vout/Iout, 1e-3)}`.
- **Sense sources** (zero-volt, all branches Guide §5 requires):
  `Vin_sense`, `Vpri_sense_a`, `Vpri_sense_b`, `Vsec1_sense_o1` (per
  output), `Vout_sense_o1`.
- **Bridge-voltage probe** (Guide §5 mandates one):
  `Evab vab 0 VALUE={V(node_DrainQ1) - V(node_DrainQ2)}`.
- **Snubbers** (push-pull leakage spikes are notorious — mandatory):
  `Rsnub_Q1 node_DrainQ1 snubQ1 100`, `Csnub_Q1 snubQ1 0 1n` (per switch).
  Same for D_pos, D_neg.
- **Models** (verbatim, copy from `Boost.cpp`):
  ```
  .model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg
  .model DIDEAL D(IS=1e-12 RS=0.05)
  ```
- **PWM sources** with explicit overlap:
  ```
  Vpwm1 pwm1 0 PULSE(0 5 0      10n 10n {D*Tsw} {Tsw})
  Vpwm2 pwm2 0 PULSE(0 5 {Tsw/2} 10n 10n {D*Tsw} {Tsw})
  ```
  When D > 0.5 these pulses naturally overlap by (2D−1)·Tsw/2 each.
- **Solver options** (mandatory, GEAR is non-negotiable for overlap
  modes):
  ```
  .options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500
  .options METHOD=GEAR TRTOL=7
  ```
- **Time step**: `period/200` baseline; `period/500` if `fsw > 200 kHz`
  or active-clamp ZVS variant.
- **IC pre-charge** (skipping this guarantees 50–100 cycles of inrush):
  `.ic v(out_sense)={Vout} i(L1a)={Iin/2} i(L1b)={Iin/2}`.
- `.tran {stepTime} {totalTime} {steadyStateStart} UIC`.
- **Saved signals**:
  ```
  .save v(vab) v(node_DrainQ1) v(node_DrainQ2) v(out_node) v(out_sense)
        i(Vin_sense) i(Vpri_sense_a) i(Vpri_sense_b)
        i(Vsec1_sense_o1) i(Vout_sense_o1)
        i(L1a) i(L1b) i(Lpri_a) i(Lpri_b)
        i(D3a) i(D3b)
  ```
- **Comment header**: topology, Vin/Vout/Pout/fsw/D/n/L1/Lm, "generated
  by OpenMagnetics".
- **Parasitic R**: bake DCR ≈ 50 mΩ into L1a/L1b and ESR ≈ 5 mΩ into
  Co — purely lossless reactive networks make ngspice unhappy. Same
  reason as Cuk plan §6.

### 6.3 simulate_and_extract_operating_points

Follow Guide §6 verbatim:

1. Construct `NgspiceRunner`.
2. If `runner.is_available()` is false, warn and fall back to analytical.
3. For each (Vin, OP), build netlist, run, on `success==false` log
   error, fall back to `process_operating_point_for_input_voltage`,
   tag name with `[analytical]`.
4. On success, build `WaveformNameMapping` and call
   `NgspiceRunner::extract_operating_point`.

Mapping (V1) — primary "magnetic" is the main transformer; emit 4
winding excitations (Np_a, Np_b, Ns_a, Ns_b). The input coupled L1
emerges through `get_extra_components_inputs`.

### 6.4 simulate_and_extract_topology_waveforms

```cpp
wf.set_input_voltage(getWaveform("vin_sense"));
wf.set_input_current(getWaveform("vin_sense#branch"));
wf.get_mutable_output_voltages().push_back(getWaveform("out_sense"));
wf.get_mutable_output_currents().push_back(getWaveform("vout_sense_o1#branch"));
```

---

## 7. Schema additions (`MAS/schemas/inputs/topologies/weinberg.json`)

New file. Pattern off `flyback.json` (turns ratio + isolation) +
`pushPull.json` (push-pull drive) + `cllcResonant.json` (variant
flag):

```jsonc
{
    "title": "weinberg",
    "description": "Weinberg current-fed push-pull-derivative DC-DC converter (ESA Weinberg 1974)",
    "type": "object",
    "properties": {
        "inputVoltage":         { "$ref": "../../utils.json#/$defs/dimensionWithTolerance" },
        "diodeVoltageDrop":     { "type": "number" },
        "maximumSwitchCurrent": { "type": "number" },
        "currentRippleRatio":   { "type": "number" },
        "efficiency":           { "type": "number", "default": 0.95 },

        "variant": {
            "enum": ["classic", "bridge", "bidirectional", "softSwitched"],
            "default": "classic"
        },
        "rectifierType": {
            "enum": ["centerTapped", "fullBridge"],
            "default": "centerTapped"
        },
        "energyRecoveryStyle": {
            "enum": ["perSwitch", "centerTap", "none"],
            "default": "perSwitch"
        },
        "synchronousRectifier": { "type": "boolean", "default": false },
        "mainMagnetic": {
            "enum": ["transformer", "inputCoupledInductor"],
            "default": "transformer",
            "description": "Which of the two magnetics is the 'main' design target"
        },

        "inputCoupledInductance":   { "type": "number" },        // L1
        "outputInductance":         { "type": "number" },        // optional Lo
        "outputCapacitance":        { "type": "number" },
        "couplingCoefficientInput": { "type": "number", "default": 0.999 },
        "couplingCoefficientMain":  { "type": "number", "default": 0.99 },

        "operatingPoints": {
            "type": "array",
            "items": { "$ref": "#/$defs/weinbergOperatingPoint" },
            "minItems": 1
        }
    },
    "required": ["inputVoltage", "diodeVoltageDrop", "operatingPoints"],
    "$defs": {
        "weinbergOperatingPoint": {
            "title": "weinbergOperatingPoint",
            "$ref": "../../utils.json#/$defs/baseOperatingPoint"
        }
    }
}
```

After schema edit, regenerate `MAS.hpp` and copy to `build/MAS/MAS.hpp`
(per user-memory staging gotcha).

Add `"WEINBERG_CONVERTER"` to the `Topologies` enum (alongside
BUCK_CONVERTER, BOOST_CONVERTER, etc.).

---

## 8. Class structure (`Weinberg.h` / `Weinberg.cpp`)

Mirror the `IsolatedBuckBoost` / `Llc` pattern (closest 2-magnetic
siblings):

```cpp
namespace OpenMagnetics {
using namespace MAS;

class Weinberg : public MAS::Weinberg, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;   // Weinberg is slow to settle (large L1)

    // 2.2 Computed design values (populated by process_design_requirements)
    double computedTurnsRatio              = 0;   // n = Np/Ns (each half)
    double computedInputCoupledInductance  = 0;   // L1
    double computedMagnetizingInductance   = 0;   // Lm of main xfmr
    double computedOutputInductance        = 0;   // Lo (optional)
    double computedOutputCapacitance       = 0;
    double computedDutyCycle               = 0;
    double computedDiodeVoltageDrop        = 0.6;
    double computedDeadTime                = 200e-9;  // for buck regime D < 0.5

    // 2.3 Schema-less device parameters
    double mosfetCoss = 200e-12;
    double leakageInductanceFraction = 0.02;     // Llk = 2% of Lpri_half

    // 2.4 Per-OP diagnostics — see §5.2

    // 2.5 Extra-component waveforms (cleared in process_operating_points)
    mutable std::vector<Waveform> extraL1aVoltageWaveforms;
    mutable std::vector<Waveform> extraL1aCurrentWaveforms;
    mutable std::vector<Waveform> extraL1bVoltageWaveforms;
    mutable std::vector<Waveform> extraL1bCurrentWaveforms;
    mutable std::vector<Waveform> extraLoVoltageWaveforms;     // if Lo present
    mutable std::vector<Waveform> extraLoCurrentWaveforms;

public:
    bool _assertErrors = false;

    Weinberg(const json& j);
    Weinberg() {};

    // 2.7-2.10 boilerplate accessors per Guide §2 (copy IsolatedBuckBoost shape)

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
        const WeinbergOperatingPoint& opPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        double inputCoupledInductance);

    // 2.13 Static analytical helpers
    static double compute_duty_cycle(double Vin, double Vout, double n,
                                     double diodeDrop, double efficiency);
    static double compute_conversion_ratio_boost(double D, double n);
    static double compute_conversion_ratio_buck(double D, double n);
    static double compute_switch_peak_voltage(double Vin, double D);
    static double compute_overlap_fraction(double D);   // max(0, 2D−1)
    static double compute_l1_min(double Vin, double D, double dIL1, double fsw);
    static double compute_lm_min(double Vin, double D, double dImag, double fsw);
    static double compute_dcm_K_boost(double L1, double D, double R, double fsw);
    static double compute_rhp_zero_freq(double R, double D, double n, double L1);
    static int    detect_operating_regime(double D);   // 0 buck, 1 boundary, 2 boost

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

class AdvancedWeinberg : public Weinberg {
    /* desiredTurnsRatios, desiredMagnetizingInductance, desiredInputCoupledInductance,
       desiredOutputInductance, desiredOutputCapacitance */
    /* Inputs process(); */
};

void from_json(const json& j, AdvancedWeinberg& x);
void to_json(json& j, const AdvancedWeinberg& x);
}
```

The "primary" magnetic returned by `process_design_requirements`
defaults to the **main transformer** (4-winding center-tapped pri+sec,
isolation = 2 sides). The input coupled inductor L1 (2-winding,
1 isolation side) flows through `get_extra_components_inputs`. User
can flip via `mainMagnetic` schema flag (§7).

---

## 9. Tests (`tests/TestWeinberg.cpp`) — DAB-quality bar

Mirror DAB and Llc test sets. Required `TEST_CASE`s per Guide §8:

| Test | Purpose | Tolerance |
|---|---|---|
| `Test_Weinberg_Weinberg_Reference_Design` | reproduce Weinberg & Schreuders 1985 fixture (50 V → 150 V, 1.5 kW, 50 kHz, D≈0.6, 1:3) | ±5 % on D, M, V_Q,pk, I_Q,pk |
| `Test_Weinberg_BDR_Reference_Design` | Yadav-Gowrishankara fixture (42 V → 100 V, 500 W, 100 kHz) | ±5 % |
| `Test_Weinberg_HighPower_Reference_Design` | Bhandari et al. 2025 fixture (50 V → 100 V, 5 kW BDR module) | ±10 % |
| `Test_Weinberg_Design_Requirements` | turns ratio, Lm, L1, Co positive and within sane ranges; round-trip a power target | algebraic |
| `Test_Weinberg_OperatingPoints_Generation` | multiple Vin (min/nom/max), volt-sec balance holds per primary half, iL1a == iL1b by symmetry | algebraic |
| `Test_Weinberg_Operating_Regimes` | sweep D from 0.3 to 0.85; `lastOperatingRegime` flips buck → boundary → boost; M matches (E1)/(E2) | exact within 5 % |
| `Test_Weinberg_Overlap_Fraction` | `lastOverlapFraction = max(0, 2D-1)` exactly | exact |
| `Test_Weinberg_CCM_DCM_Boundary` | sweep load; `lastOperatingMode` flips at K_crit | sign change |
| `Test_Weinberg_Switch_Peak_Voltage` | `lastSwitchPeakVoltage = 2·Vin/(1−D)` for D ∈ {0.5, 0.6, 0.7, 0.8} | ±5 % |
| `Test_Weinberg_Flux_Balance_Symmetric_Drive` | with identical D for Q1/Q2, magnetizing flux returns to zero each cycle | ‖∫Vp dt‖ < 1 % |
| `Test_Weinberg_Flux_Walk_Asymmetric_Drive` | with 1 % D mismatch, flux walks linearly (sanity test the mismatch is detectable in the model) | non-zero drift |
| `Test_Weinberg_RHP_Zero_Frequency` | `lastRhpZeroFrequency` matches (E12) | exact |
| `Test_Weinberg_Diagnostics_Populated` | every diagnostic in §5.2 is non-zero after a successful OP | non-null |
| `Test_Weinberg_RectifierType_CT_vs_FB` | CT-FW vs FB rectifier; output ripple frequency 2·fsw both ways; secondary RMS scales correctly | ratio |
| `Test_Weinberg_SPICE_Netlist` | snubbers, METHOD=GEAR, IC pre-charge, **2 sw + energy-recovery diodes** (V1) or **4 sw** (V2), 4-way K coupling on main xfmr, 2-way K on input L | string `find` |
| `Test_Weinberg_DivideByZero_Guard` | Iout=0 → no crash; Rload clamps to 1e-3 | no throw |
| `Test_Weinberg_PtP_AnalyticalVsNgspice_BoostRegime` | Weinberg-Schreuders fixture (D≈0.6); primary-side current NRMSE | **≤ 0.15** |
| `Test_Weinberg_PtP_AnalyticalVsNgspice_BuckRegime` | D < 0.5 fixture; primary current NRMSE | **≤ 0.15** |
| `Test_Weinberg_PtP_InputCoupledInductor_Symmetry` | NRMSE between iL1a and iL1b waveforms — should be near zero in steady state | ≤ 0.05 |
| `Test_Weinberg_ExtraComponents_L1_Lo_Co` | `get_extra_components_inputs(IDEAL/REAL)` returns 2 (V1 no Lo) or 3 (V1 with Lo) entries; reverses for `mainMagnetic = inputCoupledInductor` | count + sanity |
| `Test_Weinberg_Bridge_Variant_SwitchVoltage` | V2; `lastSwitchPeakVoltage = Vin` (not 2·Vin) | exact |
| `Test_Weinberg_Bidirectional_PowerBalance` | V3; reverse direction Vin terminal absorbs power | sign check |
| `Test_Weinberg_Synchronous_Variant` | V3 / synchronousRectifier=true; netlist contains 4 SR MOSFETs (no diodes); simulation runs | smoke |
| `Test_Weinberg_Watkins_Johnson_Disambiguation` | a fixture *labelled* W-J should NOT parse as Weinberg (verify the schema rejects W-J specs that don't have `inputCoupledInductance`) | throws |
| `Test_Weinberg_Waveform_Plotting` | CSV dump under `tests/output/`; visual inspection | not asserted |
| `Test_AdvancedWeinberg_Process` | round-trip user-specified turnsRatios/Lm/L1/Lo/Co | round-trip equal |

NRMSE helpers (`ptp_nrmse`, `ptp_interp`, `ptp_current`,
`ptp_current_time`) extracted to a shared header
`tests/TestTopologyHelpers.h` (also a CLLC and Cuk plan deliverable
— do it once, all three plans benefit).

---

## 10. Phased delivery

| # | Scope | Files | Tests added | Effort |
|---|---|---|---|---|
| **P1** | MAS schema (`weinberg.json`, `Topologies` enum entry), regenerate MAS.hpp, stub `Weinberg.h`/`Weinberg.cpp` (constructors + run_checks only, throws on the rest), `from_json`/`to_json` for `AdvancedWeinberg` | schema, `Weinberg.h`, `Weinberg.cpp` | none | 0.5 d |
| **P2** | Static analytical helpers (`compute_duty_cycle`, `compute_conversion_ratio_*`, `compute_switch_peak_voltage`, `compute_overlap_fraction`, `compute_l1_min`, `compute_lm_min`, `compute_dcm_K_boost`, `compute_rhp_zero_freq`, `detect_operating_regime`); unit tests | `Weinberg.cpp`, tests | `Operating_Regimes`, `Overlap_Fraction`, `Switch_Peak_Voltage`, `RHP_Zero_Frequency` | 0.5 d |
| **P3** | `process_design_requirements` (size n, Lm, L1, Co); `process_operating_point_for_input_voltage` (CCM piecewise-linear waveforms for main xfmr + input coupled L, with overlap handling); `process_operating_points` (loop) | `Weinberg.cpp` | `Weinberg_Reference_Design`, `BDR_Reference_Design`, `Design_Requirements`, `OperatingPoints_Generation`, `Diagnostics_Populated` | 2 d |
| **P4** | DCM detection & flag; flux-balance / flux-walk math + tests | `Weinberg.cpp`, tests | `CCM_DCM_Boundary`, `Flux_Balance_Symmetric_Drive`, `Flux_Walk_Asymmetric_Drive` | 1 d |
| **P5** | SPICE netlist (V1 classic): full Guide §5 contract + IC pre-charge + GEAR + snubbers + 4-way K coupling on main xfmr + 2-way K on input L + energy-recovery D3 + Vab probe + sense sources. `simulate_and_extract_operating_points` with analytical fallback. `simulate_and_extract_topology_waveforms`. | `Weinberg.cpp` | `SPICE_Netlist`, `DivideByZero_Guard` | 2 d |
| **P6** | NRMSE acceptance: tighten until ≤ 0.15 for boost-regime and buck-regime fixtures + iL1a/iL1b symmetry. Failure-discovery phase — expect to fix 2–3 SPICE convergence bugs (leakage too low, IC mismatch, wrong K-matrix structure for 4-winding xfmr) | tests + minor SPICE fixes | `PtP_AnalyticalVsNgspice_BoostRegime`, `PtP_AnalyticalVsNgspice_BuckRegime`, `PtP_InputCoupledInductor_Symmetry` | 1.5 d |
| **P7** | `get_extra_components_inputs` returning L1, Co (and Lo if present) ancillary `Inputs` (IDEAL + REAL modes); honour `mainMagnetic` flag flip | `Weinberg.cpp` | `ExtraComponents_L1_Lo_Co` | 0.5 d |
| **P8** | Rectifier-type-aware secondary (CT FW vs FB); test both | `Weinberg.cpp`, tests | `RectifierType_CT_vs_FB` | 1 d |
| **P9** | V2 bridge variant: schema flag, 4-switch primary netlist, switch-V_DS halving | `Weinberg.cpp` | `Bridge_Variant_SwitchVoltage` | 1 d |
| **P10** | V3 bidirectional + V5 synchronous-rectifier variants: schema flags, active-SR netlist, reverse-flow OP | `Weinberg.cpp` | `Bidirectional_PowerBalance`, `Synchronous_Variant` | 1.5 d |
| **P11** | High-power reference design (Bhandari et al. 2025), AdvancedWeinberg class, header docstring rewrite (ASCII art, equations §4, references, accuracy disclaimer per Guide §3, **Watkins-Johnson disambiguation**), waveform plotting | `Weinberg.h`, `Weinberg.cpp`, tests | `HighPower_Reference_Design`, `Watkins_Johnson_Disambiguation`, `AdvancedWeinberg_Process`, `Waveform_Plotting` | 1 d |

**Total ~12.5 working days** for V1+V2+V3+V5 coverage with both
rectifier types.

Dependencies: P1 → P2 → P3 → P4; P5 needs P3; P6 needs P5; P7 ⊥ P5;
P8 needs P5; P9/P10 need P5+P8; P11 last.

---

## 11. Acceptance criteria (per Guide §15)

- [ ] `ninja -j2 MKF_tests` builds clean.
- [ ] All `TestWeinberg` cases pass.
- [ ] At least 2 `*_PtP_AnalyticalVsNgspice_*` tests run and NRMSE ≤ 0.15
      (boost and buck regimes).
- [ ] `lastOperatingRegime` correctly identifies buck / boundary /
      boost across D ∈ {0.3, 0.5, 0.7}.
- [ ] `lastOperatingMode` flips CCM ↔ DCM at K_crit boundary within
      ±5 %.
- [ ] Symmetric-drive flux balance test passes (drift < 1 %).
- [ ] Asymmetric-drive flux walk test passes (drift detected).
- [ ] Schema documented; `Topologies::WEINBERG_CONVERTER` added.
- [ ] Header has ASCII art, equations §4 with citations,
      accuracy disclaimer (Coss neglected, body-diode reverse-recovery
      neglected, leakage spike clamped only by D3), **Watkins-Johnson
      disambiguation note**, **regime table** (buck/boundary/boost).
- [ ] Netlist parses, contains snubbers, GEAR, IC pre-charge for L1
      and Co, **4-way K coupling on main transformer**, **2-way K on
      input coupled inductor**, energy-recovery diodes per
      `energyRecoveryStyle`, Evab probe, sense sources, rectifier-type
      branch (CT FW vs FB).
- [ ] `simulate_and_extract_operating_points` runs SPICE for real,
      falls back to analytical with `[analytical]` tag on failure.
- [ ] `simulate_and_extract_topology_waveforms` references only nodes
      the netlist actually emits.
- [ ] `get_extra_components_inputs` supports `IDEAL` and `REAL` for
      L1 and Co.
- [ ] V2 bridge variant achieves switch-V_DS halving (Vin instead of
      2·Vin/(1−D)).
- [ ] V3 bidirectional variant Vin terminal absorbs power in reverse
      mode.

---

## 12. Common Weinberg implementation mistakes (avoid these)

From the research dossier — explicit pitfalls flagged in the
literature:

1. **Asymmetric Q1/Q2 timing → flux walk in main xfmr.** Mandatory:
   peak-current-mode control or DC-blocking series cap on primary
   or strict on-time matching (< 50 ns). Surface
   `lastFluxImbalanceMargin` so it shows up in design reports.
2. **Wrong input coupled-L turns ratio.** Must be **exactly 1:1**
   with `k > 0.99`. Any asymmetry → one half handles more current →
   thermal runaway. Hard-code n = 1:1 in V1; expose as schema field
   only if a future variant needs it.
3. **Inadequate energy-recovery clamp.** Without D3, leakage Llk
   causes Q1/Q2 V_DS to ring to >2·Vin/(1−D) → avalanche. The D3
   path returns leakage energy to Vin — the original "new
   energy-transfer principle" of Weinberg's 1974 paper. Schema
   `energyRecoveryStyle` must default to `perSwitch` (always
   present); the `none` value is exposed only for academic-comparison
   experiments.
4. **Treating it as plain current-fed push-pull (Watkins-Johnson)
   and missing the overlap requirement.** W-J uses a *single* input
   inductor; Weinberg uses a *coupled 1:1* inductor that
   intrinsically supports overlap without needing a freewheel
   switch. Drop overlap and the topology stops working. Document in
   header docstring per Guide §3.
5. **Ignoring leakage ringing → no snubber.** Even with D3, the
   1–3 % leakage produces 3–5× resonant V_DS overshoot at hard
   turn-off. Always add an RC or RCD snubber; or move to V5
   (active-clamp ZVS).
6. **L1 saturation at minimum-Vin / max load.** L1 sees full Iin
   (not Iin/2) plus ΔI_L1; designers often size for nominal and
   forget worst-case battery-low. Size for `Iin_max =
   Iout / (η · M_min)` at `Vin_min` and `D_max`.
7. **Forgetting the RHP zero** when designing the voltage loop.
   Boost mode (D > 0.5) Weinberg has the same `ω_RHPZ ≈ R(1−D)²/(n²·L1)`
   limitation as any boost-derived topology. Surface
   `lastRecommendedLoopBandwidth = fRHP/5`.
8. **Cold-start inrush.** Without `.ic`, SPICE rings from zero state
   for 50–100 cycles. Pre-charge L1a, L1b, Co — non-negotiable.
9. **K = 1.0 in SPICE** for either magnetic — singular factorization;
   always 0.999 for input L1 (very tight coupling) and 0.99 for main
   xfmr (some leakage required).
10. **No parasitic R in SPICE**. Push-pull current-fed with overlap
    is *the* convergence-killer category in ngspice. Mandatory:
    Llk on every transformer winding (2 % of self), DCR on every
    inductor (50 mΩ), ESR on every cap (5 mΩ).
11. **Confusing CT-FW rectifier with FB rectifier**. CT-FW uses 2
    diodes (or SR MOSFETs) and a center-tapped secondary winding;
    FB uses 4 diodes and a single secondary winding. Output ripple
    frequency is 2·fsw in *both* cases (because of push-pull drive),
    but the secondary winding's RMS current and number of windings
    differ. Branch the netlist on `rectifierType` (Guide §4.5
    rectifier-type-aware secondaries).

---

## 13. Magnetic structure summary (the bit MKF cares about most)

| Object | Windings | Coupling | DC bias | Gap | MKF mapping |
|---|---|---|---|---|---|
| **Input coupled L (L1)** | 2 (1:1) | k = 0.999 (very tight) | Iin/2 each, OPPOSING (cancels in 1:1 ideal) | gapped (stores energy) | 2-winding magnetic with constrained turns ratio = 1:1; mirror `DifferentialModeChoke` shape |
| **Main transformer** | 4 (Np_a, Np_b, Ns_a, Ns_b) — both pri and sec center-tapped | k = 0.99 (some leakage usable for ZCS) | 0 (symmetric push-pull drive) | ungapped | 4-winding transformer; mirror `PushPull` for drive scheme + `Flyback` for CT-secondary rectifier handling |

**`process_design_requirements` returns**: the **main transformer**
by default (4-winding object, 2 isolation sides). `mainMagnetic`
schema flag flips it to the input coupled L if the user wants to
design that one.

**`get_extra_components_inputs` returns**: whichever magnetic isn't
"main", plus Co, plus optional Lo. So: usually `[L1, Co]` (V1, no
Lo), or `[L1, Lo, Co]` (V1 with Lo), or `[Transformer, Co]` if
`mainMagnetic = inputCoupledInductor`.

**Multi-output**: V4 multi-output Weinberg adds more secondaries to
the same transformer. Use the same per-output indexing convention as
DAB / IsolatedBuckBoost (`output_voltages[]`, `output_currents[]` in
the OP, `o<N>` suffix in netlist nodes). Initial scope is single-output.

**Novel constraint MKF doesn't have today**: the 1:1 input coupled
inductor needs a **constrained turns ratio = 1:1 with high coupling**
guaranteed. The current `Topology::create_isolation_sides` doesn't
have a "turns-ratio constraint" knob. Investigate either:

- (a) add an optional `turnsRatioConstraint` field to
  DesignRequirements (clean but cross-cutting),
- (b) bake the 1:1 into `get_extra_components_inputs` and rely on
  the user not to override it (V1 default; less canonical).

Recommend (b) for the initial add, matching the Cuk plan §13
recommendation. (a) becomes a follow-up if more topologies need it.

---

## 14. Checklist for the implementing agent

Before each phase:

- [ ] Re-read `CONVERTER_MODELS_GOLDEN_GUIDE.md`.
- [ ] Verify the latest `PushPull.cpp` SPICE pattern (closest **drive
      scheme** sibling — center-tapped primary, 2 switches, 180°
      phase, flux-walk concerns).
- [ ] Verify the latest `Llc.cpp` for the **multi-magnetic + ancillary**
      structure pattern (Lr + xfmr).
- [ ] Verify the latest `Flyback.cpp` for the CT-secondary rectifier
      pattern.
- [ ] Verify the latest `IsolatedBuckBoost.cpp` for the
      "main magnetic + extra components" emission pattern.
- [ ] Verify the latest `DifferentialModeChoke.cpp` for the
      **constrained-turns-ratio 2-winding magnetic** pattern (input L1).
- [ ] Verify the latest `ActiveClampForward.cpp` for energy-recovery
      diode + clamp-cap modelling (closest analogue to D3).
- [ ] Build with `ninja -j2` in the build dir, not `cmake --build`
      (per user memory).
- [ ] After MAS schema edits, copy `MAS.hpp` from source to
      `build/MAS/MAS.hpp` (per user memory).
- [ ] Mark each phase as complete only when its NRMSE / sign-change /
      build-clean criteria are met.

After P11 (final):

- [ ] Update CLAUDE.md / user-memory entry on converter quality tiers
      to include Weinberg as DAB-quality.
- [ ] Open follow-up issues for: full DCM design; three-phase
      Weinberg (Larico-Barbi 2012); paralleled Weinberg with current
      sharing (Mazumder); non-isolated Weinberg variant; multi-output
      Weinberg.
- [ ] Add a separate `CurrentFedPushPull.{h,cpp}` (Watkins-Johnson)
      topology if requested — explicitly distinct from Weinberg.

---

## 15. Reference: equations cheat-sheet (paste into header docstring)

```
TOPOLOGY OVERVIEW (V1 classic Weinberg)
───────────────────────────────────────
            L1a               Np_a         D_pos       Lo (opt)
   Vin ──┳──[~~~]───┳────────[==o==]────┬───┤▶├────┬───[~~~]──┳── Vout+
         │  k=0.999 │                   │          │           │
         │  L1b     │                   │ T 4-w    │          Co
         │  [~~~]   │                   │ k=0.99   │           │
         └──────────┴──Np_b──[==o==]────┘          │           │
                       │      │                    │          GND
                      Q1     Q2  ─── center-tap   │
                       │      │                    │
                     D3a/D3b clamp to Vin          │
                       │      │                   D_neg
                      GND    GND  ◀── Ns_b ────────┘

MODULATION CONVENTION
─────────────────────
D = duty of EACH switch (each ON for D·Tsw within its half-period)
Q1 phase = 0;  Q2 phase = Tsw/2 (180° offset)
Overlap fraction = max(0, 2D − 1)

| Regime          | D range        | Switch state         | M = Vout/Vin (ideal) |
|-----------------|----------------|----------------------|----------------------|
| Buck-like       | 0 < D < 0.5    | non-overlap (dead)   | D / n                |
| Boundary        | D = 0.5        | continuous, no overlap| 1 / n                |
| Boost-like      | 0.5 < D < 1    | overlap (2D−1)·Tsw/2 | 1 / (2·n·(1 − D))    |

KEY EQUATIONS (CCM, 1:1 input coupled-L, n = Np/Ns)
───────────────────────────────────────────────────
M (boost)        = Vout/Vin = 1 / (2·n·(1 − D))                         — Weinberg PESC 1985
M (buck)         = Vout/Vin = D / n
V_Q,pk           = 2·Vin / (1 − D)  ≈ 2·n·Vout                          — clamped by D3
I_Q,pk           = Iin + ΔI_L1/2  (non-overlap)
                 = Iin/2 + ΔI_L1/2  (during overlap)
V_D,pk           = 2·Vout                          (CT FW rectifier)
I_Co,rms         ≈ Iout · √((1−D)/D)
ω_RHPZ (boost)   ≈ R · (1 − D)² / (n² · L1)                             — Erickson §8

L1,min           ≥ Vin · (2D − 1) · Tsw / (4 · ΔI_L1)
Lm,min           ≥ Vin · D · Tsw / (4 · ΔI_mag)
Co,min           ≥ Iout · D / (ΔVo · fsw)

K_crit (boost)   = 2 · L1 · fsw · (1 − D)² / R
CCM if K_crit < 1, else DCM

DISAMBIGUATION
──────────────
This is NOT the Watkins-Johnson topology (single input inductor, 1
switch, M = (2D−1)/(D·n)). The Weinberg's defining feature is the
1:1 INPUT COUPLED INDUCTOR with energy-recovery diode D3, which
intrinsically handles the overlap interval without needing a
freewheel switch. See Marsala & Olivieri, IEEE 1996.

This is NOT a plain current-fed push-pull. The plain version uses a
single input inductor and needs an explicit clamp/snubber.

REFERENCES
[1] A. H. Weinberg, "A Boost Regulator with a New Energy-Transfer
    Principle", ESA/ESTEC Spacecraft Power Conditioning Electronics
    Seminar, 1974.
[2] A. H. Weinberg & L. Schreuders, "A High-Power High-Voltage DC-DC
    Converter for Space Applications", IEEE PESC ~1985.
[3] A. Weinberg, "Modelling Technique for DC-DC Converters Using
    Weinberg Topology", IEEE 1996.
[4] R. C. Larico & I. Barbi, "Three-Phase Weinberg Isolated DC-DC
    Converter", IEEE Trans. Ind. Electron., 59(2):888–896, 2012.
[5] Yadav & Gowrishankara, "Analysis and Development of Weinberg
    Converter with Fault-Tolerant Feature", IEEE 2016.
[6] Bhandari et al., "Low-output-current-ripple Weinberg converter
    for BDR in GEO satellites", J. Power Electron. 25:591, 2025.
[7] Erickson & Maksimović, "Fundamentals of Power Electronics" 3rd
    ed. (2020), §8 (RHP zero in boost-derived topologies).
[8] Marsala & Olivieri, "Comparison Between Two Current-Fed
    Push-Pull DC-DC Converters", IEEE 1996.
[9] ESA SP-589 (2005) — Galileo GSTB-V2/A BDR.

ACCURACY DISCLAIMER
SPICE model omits MOSFET Coss variation with V_DS, body-diode
reverse recovery, gate-driver propagation delay, transformer
inter-winding capacitance, and core-loss resistance. Predicted η is
an upper bound; expect 2–4 % loss vs. measured at full load.
Energy-recovery diode D3 is modeled as ideal — actual leakage-spike
overshoot may exceed 2·Vin/(1−D) by 10–20 % unless an explicit
RCD snubber is added.

POLARITY / FLOW CONVENTION
Power flows Vin → Vout in V1 (unidirectional). Vout > 0. Iout
positive when flowing from Vout through Rload to GND. In V3
(bidirectional) `powerFlow` field selects direction; secondary SR
MOSFETs gated accordingly.
```
