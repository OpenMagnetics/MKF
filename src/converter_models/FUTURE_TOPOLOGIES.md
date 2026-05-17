# Future Topologies for MKF

A prioritised list of converter topologies that could be added to MKF
beyond the ones currently implemented or already planned. Sibling to
the existing plan documents in this directory:

- `CONVERTER_MODELS_GOLDEN_GUIDE.md` — the DAB-quality spec
- `LLC_REWORK_PLAN.md` — **DONE** (Tracks A + B): brought Llc to DAB-quality with all four rectifier variants (CT / FB / CD / VD)
- `CLLC_REWRITE_PLAN.md` — rewrite of the Cllc model
- `SEPIC_PLAN.md` — add SEPIC from scratch
- `ZETA_PLAN.md` — add Zeta from scratch (sibling to SEPIC)
- `CUK_PLAN.md` — add Cuk from scratch (third fourth-order sibling)
- `WEINBERG_PLAN.md` — add Weinberg converter
- `FOUR_SWITCH_BUCK_BOOST_PLAN.md` — add 4SBB / FSBB from scratch
  (USB-PD, Li-ion chargers; first envelope-sweep topology in MKF)
- `CLLLC_PLAN.md` — bidirectional symmetric-tank resonant converter
  (EV OBC, ESS, server PSU)
- `VIENNA_PLAN.md` — 3-phase 3-level PFC rectifier (EV DC fast
  charging, telecom)
- `ASYMMETRIC_HALF_BRIDGE_PLAN.md` — Imbertson-Mohan AHB
  (single-leg + DC-blocking cap; 30-500 W AC adapters / telecom bricks)
- `SRC_PLAN.md` — Series Resonant Converter (2-element Lr+Cr tank; step-down, high-V, plasma/X-ray)

Ranking criterion: **magnetic-design distinctiveness × industry pull**.
A topology that exercises the magnetic adviser in a new way *and* has
real production volume scores higher than a topology that's just
"another buck variant".

---

## Already implemented in MKF

| Class | Notes |
|---|---|
| `Boost.cpp` / `Buck.cpp` | DAB-quality (header, 3 ref designs, NRMSE, snubber, R_load guard, GEAR) |
| `IsolatedBuck.cpp` / `IsolatedBuckBoost.cpp` | |
| `Flyback.cpp` | **DAB-quality** — 3 ref designs (PMP30817 1.2 W, LM5180EVM-DUAL 3 W, TIDA-00709 33 W) with iPri NRMSE ≤ 0.15 vs SPICE; snubber split into ideal-K=1 (100 kΩ, no DC bias) and real-magnetic (1 kΩ, Lk-spike damped) paths; shape-match OP uses SPICE-settled Iout; 256-shift NRMSE alignment |
| `SingleSwitchForward.cpp` / `TwoSwitchForward.cpp` / `ActiveClampForward.cpp` | |
| `PushPull.cpp` | **DAB-quality** — 3 ref designs (SN6501 1.2 W, SN6505B 3.3 W, SN6507 5 W) with iPri NRMSE ≤ 0.15 vs SPICE; shape-match OP uses SPICE-settled Iout; 256-shift NRMSE alignment |
| `PhaseShiftedFullBridge.cpp` | rectifier-type-aware (CT/CD/FB); 3 ref designs + ZVS + multi-output + NgSpice PtP gates (current iPri NRMSE ≤ 0.20; commutation-ramp refinement pending to reach DAB-quality ≤ 0.15). FW decay model τ = Lr/(n²·R_d_diode) per Vlatkovic 1992 §III asymmetric-commutation. |
| `PhaseShiftedHalfBridge.cpp` | **DAB-quality** — rectifier-type-aware (CT/CD/FB); 3 ref designs (NRMSE ≤ 0.10), ZVS test, physics-derived freewheel τ = Lr/(2·RON), current-dependent diode drop matched to SPICE |
| `Dab.cpp` | gold reference, all four DAB modulations |
| `Llc.cpp` | DAB-quality (Tracks A + B): ZVS diagnostics, ITL=500/500, multi-output warn, 41 test cases / 1335 assertions; rectifier-type-aware (CT / FB / CD / VD) with per-variant netlist + extras (CD: 2× MAS::Inputs for Lo; VD: 2× CAS::Inputs for stacked Co) |
| `Cllc.cpp` | resonant; rewrite plan exists (`CLLC_REWRITE_PLAN.md`) |
| `Sepic.cpp` / `Zeta.cpp` / `Cuk.cpp` | DAB-quality (header, 3 ref designs, NRMSE, snubber, GEAR) |
| `Weinberg.cpp` | DAB-quality (header, 3 ref designs, NRMSE, snubber, GEAR) |
| `FourSwitchBuckBoost.cpp` | DAB-quality (3 ref designs spanning BUCK/BOOST/BUCK_BOOST regions, NRMSE, VoltSec, ConvPort) |
| `AsymmetricHalfBridge.cpp` | DAB-quality (3 ref designs, ZVS test, VoltSec_AllWindings for CT/FB/CD) |
| `Clllc.cpp` | Phase A (analytical 4-state ODE solver, fwd + rev, above-resonance) **plus Phase B-1** (8-MOSFET SPICE netlist with 3 K-couplings, GEAR/ITL=500/500, direction-aware bridges; `simulate_and_extract_*` DAB-pattern wrappers with analytical fallback; NRMSE-vs-analytical PtP gates fwd + rev) **plus Phase B-2** (`get_extra_components_inputs`: Cr1+Cr2 CAS::Inputs always, HV+LV bus-cap CAS::Inputs in REAL mode, discrete Lr1+Lr2 MAS::Inputs when `integratedResonantInductors == false`; Lr voltage waveforms produced by central-difference `L·di/dt`). `AdvancedClllc::process` + reference designs still pending (see `CLLLC_PLAN.md`) |
| `Vienna.cpp` | pragmatic Phase 1+2 (peak-of-line analytical) — no SPICE, no DAB-quality tests |
| `Src.cpp` | pragmatic Phase 1+2 (FHA solver) — no SPICE, no DAB-quality tests |
| `PowerFactorCorrection.cpp` | basic PFC (verify totem-pole coverage) |
| `CommonModeChoke.cpp` / `DifferentialModeChoke.cpp` / `CurrentTransformer.cpp` | passive magnetic components |

## Already planned (have plan docs)

- **SEPIC** — `SEPIC_PLAN.md`
- **Zeta** — `ZETA_PLAN.md`
- **Cuk** — `CUK_PLAN.md`
- **Weinberg** — `WEINBERG_PLAN.md`
- **4-switch buck-boost** — `FOUR_SWITCH_BUCK_BOOST_PLAN.md`
- **CLLLC** — `CLLLC_PLAN.md`
- **Vienna rectifier** — `VIENNA_PLAN.md`
- **Asymmetric Half-Bridge** — `ASYMMETRIC_HALF_BRIDGE_PLAN.md`
- **Series Resonant Converter** — `SRC_PLAN.md`
- **LLC rework** — `LLC_REWORK_PLAN.md` **(DONE)**
- **CLLC rewrite** — `CLLC_REWRITE_PLAN.md`

---

## Tier 1 — high value, clear gap, broad industry use

### Cuk converter — **PLAN WRITTEN** (`CUK_PLAN.md`)

Sibling of SEPIC; inverting output, capacitive energy transfer through
Cs1+Cs2. Already a `MasMigration.cpp` enum entry
(`"Cuk Converter" → "cukConverter"`). Magnetic adviser handles it the
same way as SEPIC (uncoupled or coupled). Conversion ratio
M = −D/(1−D). Plan exists at `CUK_PLAN.md`.

### 4-switch buck-boost (4SBB) — **PLAN WRITTEN** (`FOUR_SWITCH_BUCK_BOOST_PLAN.md`)

Ubiquitous in laptops, USB-PD, drones, single-cell-Li chargers. Single
inductor + four MOSFETs operating in three regions (buck / buck-boost /
boost) selected by the controller depending on Vin/Vo. The inductor
sees three different excitation patterns over a Vin sweep — non-trivial
for the magnetic adviser to size for the worst case. No existing MKF
model covers this.

References: TI LM5176 / LM5175 datasheets; TI TIDA-01411 reference
design; ADI LT8390.

### CLLLC — **PLAN WRITTEN** (`CLLLC_PLAN.md`)

Bidirectional EV OBC / server PSU. Symmetric Lr1+Cr1+Lm+Cr2+Lr2 tank,
distinct from existing `Cllc.cpp` (which is asymmetric). Wolfspeed
22 kW CRD-22DD12N is the canonical reference design.

References: Wolfspeed CRD-22DD12N; MDPI Electronics 12(7):1605;
Infineon UG-2020-31 11 kW CLLC reference.

---

## Tier 2 — medium value, real production use, narrower scope

### Asymmetric Half-Bridge converter (AHB) — **PLAN WRITTEN** (`ASYMMETRIC_HALF_BRIDGE_PLAN.md`)

Distinct from PSHB. Same primary topology but with a DC-blocking
coupling cap; transformer reset is asymmetric. `PhaseShiftedHalfBridge.cpp`
header docstring already calls out the disambiguation; the AHB version
itself is not implemented. Common in 100–500 W server / telecom.

References: Imbertson-Mohan IEEE TIA 1993; TI SLUP223; ON AN-4153;
ST AN2852; Erickson Ch. 6.

### Vienna rectifier (3-phase PFC) — **PLAN WRITTEN** (`VIENNA_PLAN.md`)

Three boost inductors with a central neutral leg; common in EV
chargers and 3-phase telecom rectifiers. Magnetics adviser would
handle three correlated boost inductors. PFC mode (line-cycle envelope
sweep) inherits from existing `PowerFactorCorrection.cpp`.

References: Kolar & Zach 1994 (the original Vienna paper); Infineon
"3-phase PFC" app notes; TI TIDA-010210.

### Series Resonant Converter (SRC) — **PLAN WRITTEN** (`SRC_PLAN.md`)

LLC's no-Lm cousin (just Lr+Cr in series). Historical but still used
in plasma cutters, induction heating, X-ray HV supplies. Distinct
gain analytics from LLC (no parallel branch → no Lm state, 2-state
time-domain solver, 4-quadrant mode classification). Fresh class,
DAB-quality tier.

References: Steigerwald 1988 (the SRC analysis paper);
Erickson-Maksimović Ch. 19; Kazimierczuk *Fundamentals of Resonant
Power Conversion*.

### 3-phase DAB

High-power EV / grid-tied. Two 3-leg bridges, three-phase transformer.
Substantially different magnetics — 3-phase transformer with ∠120°
phase displacement. Distinct from existing `Dab.cpp` (single-phase).

References: De Doncker / Divan 1991 (original DAB paper); IEEE TPEL
papers on 3-phase DAB; Infineon EVAL_5K5W_3PH_DAB.

---

## Tier 3 — specialty, deferrable

| Topology | Why deferred | Future trigger |
|---|---|---|
| **LCC resonant** | 4th-order tank, mostly HV / LED-driver (Infineon ICL5102) | only if MKF wants HV-LED market |
| **Totem-pole PFC** | check whether existing `PowerFactorCorrection.cpp` covers it; if not, modern GaN-PFC topology worth its own class | confirm `PowerFactorCorrection.cpp` coverage |
| **Class-E / Class-D resonant inverter** | induction heating, wireless power Tx; niche but distinct magnetics | wireless-power customer ask |
| **Multi-phase / interleaved Buck** | `phaseCount` flag on existing `Buck.cpp` rather than a new class | demand for high-current point-of-load |
| **Bridgeless SEPIC PFC** | already noted in `SEPIC_PLAN.md` § 8; AC-input semantics differ | LED-driver / EV-charger ask |
| **Modular Multilevel Converter (MMC)** | grid-tied, very high power; out of MKF's typical SMPS scope | grid / HVDC ask |
| **Solid-state transformer** | research-grade; no standardisation | research ask |

---

## Recommendation order (after current plans land)

Plans written so far: SEPIC, Zeta, Cuk, Weinberg, 4SBB, CLLLC, Vienna,
AHB, **SRC (latest)**.  Remaining recommended order, in priority:

1. **3-phase DAB** (high-power EV / grid-tied, distinct 3-phase magnetics)

After each addition, update `project_converter_model_quality_tiers.md`
(memory) so future sessions know which models are DAB-quality vs
partial.
