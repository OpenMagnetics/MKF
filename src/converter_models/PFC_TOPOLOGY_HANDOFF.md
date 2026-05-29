# PFC Topology Expansion — Handoff

Branch: `pfc-topology-expansion` (off `main` @ 6b14ef1f). **Not pushed.**
Two commits:

- `678778ed` fix(inputs): bound volt-second integration to one switching period
- `944152fc` feat(pfc): review fixes + bridgeless/semi/interleaved/SEPIC/Ćuk variants

> ⚠️ The working tree also contains **pre-existing, unrelated, uncommitted**
> changes that are **NOT mine** and were deliberately left unstaged:
> `src/advisers/{CoilAdviser,CoreAdviserDataset,CoreAdviserGapping,CoreAdviserPipeline,MagneticAdviser}.cpp`,
> `src/advisers/MagneticFilterInternal.h`, `tests/TestCoilAdviser.cpp`, and the
> untracked `CAPABILITIES.md`, `HANDOFF_FAST_PATH_GAPPING.md`,
> `src/converter_models/SPICE_PROBE_HANDOFF.md`, `tests/TestWireAdviserHangRepro.cpp`.
> Do not fold these into PFC commits — they belong to someone else's in-flight work.

---

## 1. Current PFC variant support matrix

`PfcTopologyVariants` (MAS schema, 10 values). Status after this work:

| Variant | Sizing/waveforms | Switching netlist + PtP | Notes |
|---|---|---|---|
| BOOST | ✅ | ✅ | baseline |
| BRIDGELESS | ✅ | ✅ | inductor ≡ boost (series, rectified side) |
| SEMI_BRIDGELESS | ✅ | ✅ | inductor ≡ boost |
| INTERLEAVED_BOOST | ✅ (per-phase) | ✅ (per-phase cell) | `numberOfPhases` ∈ {2,3} |
| TOTEM_POLE | ✅ (bipolar waveform) | ✅ (bipolar 4-switch stage) | CCM requires `wideBandgapSwitch=true` |
| SEPIC | ✅ (analytical, CCM/CrCM) | ❌ throws | buck-boost class; DCM rejected |
| CUK | ✅ (analytical, CCM/CrCM) | ❌ throws | buck-boost class; DCM rejected |
| BUCK | ❌ throws | ❌ | discontinuous input current — no series-L path |
| BUCK_BOOST | ❌ throws | ❌ | discontinuous input current |
| VIENNA | ❌ throws | ❌ | 3-phase/3-level; own model — see VIENNA_PLAN.md |

"Throws" everywhere is a *loud, specific* exception (CLAUDE.md: no silent
fallbacks). Every variant either works or refuses with a precise reason.

### Key design split (read this before extending)
- **Boost family** (boost/bridgeless/semi-bridgeless/interleaved): the series
  inductor sits on the rectified side and sees `+Vin` during ON, `Vin−Vout`
  during OFF → identical boost sizing `L = Vin·D/(ΔI·fsw)` and the SAME
  switching netlist (`generate_ngspice_switching_circuit`). Interleaved emits
  ONE per-phase cell (`per_phase_power` = Pout/N, per-phase L).
- **Totem-pole** (bridgeless, bipolar): same boost sizing, but the inductor is
  on the AC side and carries SIGNED current. `generate_ngspice_switching_circuit`
  has a `bipolar` branch (selected for TOTEM_POLE) emitting a genuine 4-switch
  stage: a floating line source `B_vac` between `vac` and the LF-leg midpoint
  `mid_lf`; an HF leg (`S_hf_hi`/`S_hf_lo` + rectifying diodes `D_hf_hi`/
  `D_hf_lo`) where ONE switch is the active boost switch and the OPPOSITE leg's
  diode rectifies — roles swap by line polarity; and a line-frequency LF leg
  (`S_lf_hi`/`S_lf_lo`) routing the AC return to gnd (positive half) or vbus
  (negative half). Only one switch + one diode conduct at any instant, so it
  converges like the boost. The controller is REUSED verbatim (same
  `derive_pfc_controller_tuning`, type-II EAs, warm-start) but runs in the
  rectified domain: it senses `|i_L|` (`abs(I(Vl_sense))`), the multiplier
  reference is on `|vline|`, and a polarity selector (`pol_pos`) ROUTES the
  single PWM command to the active HF switch + drives the LF leg. This is how
  industry totem-pole average-current controllers actually work — no separate
  bipolar controller was needed.
- **Buck-boost class** (SEPIC/Ćuk): `is_buck_boost_class()` in the header.
  Series input inductor (continuous input current) but buck-boost conversion:
  `D = (Vout+Vd)/(Vin+Vout+Vd)` (see `calculate_duty_cycle`), L1 sees `+Vin`
  ON (so the boost L formula is shared) but `−(Vout+Vd)` OFF. Handled in all
  THREE synthesis loops (`process_operating_points`,
  `simulate_and_extract_waveforms`, `simulate_and_extract_topology_waveforms`).
- Anything with **discontinuous input current** (buck, buck-boost) has no
  shared series-inductor path → not supported.

---

## 2. Deferred work — what's left, why, and how to do it

### 2a. Totem-pole switching PtP  ✅ DONE (2026-05-29)
Implemented as the `bipolar` branch of `generate_ngspice_switching_circuit`
(see the "Totem-pole" bullet in §1). Key outcome vs the original plan: a
**separate bipolar controller was NOT needed** — the average-current loop runs
in the rectified domain (senses `|i_L|`, multiplier on `|vline|`) and a
line-polarity selector routes the single PWM command to the active HF switch,
exactly as industry totem-pole controllers do. The genuine 4-switch stage
(floating source + HF leg with role-swapped active-switch/diode + LF polarity
leg) converged on the first attempt with the existing PFC solver settings (only
one switch + one diode conduct at a time, so it's boost-like numerically).
PtP: `TestPfcReferenceDesignsPtp.cpp` "Totem-pole 100 W" — envelope gate
compares against the SIGNED sine; input-port check uses a zero-mean bipolar
line (`check_pfc_switching_ports(..., bipolarInput=true)`, fed `v(vline)`).
Sim wall-time ≈ 10 s, bus lands ~382 V (−4.5 %, inside ±6 %), all gates pass.

### 2b. SEPIC/Ćuk switching PtP  (highest value remaining)
4th-order plant (L1, L2, coupling cap, Cout). Needs its own converging
controller. Analytical sizing is done and correct for CCM/CrCM; DCM needs the
`Le = L1‖L2` model (currently rejected in `calculate_inductance_dcm`).

### 2c. Buck / buck-boost PFC
Different topology (discontinuous input current, inductor not in series with
the line). Low industry use. Only do if specifically needed.

### 2d. Vienna
3-phase, 3-level. Keep it a separate model (VIENNA_PLAN.md) — do NOT shoehorn
into the PFC variant branch; the code intentionally redirects.

---

## 3. ⚠️ Caveat the next agent must know: PtP bus-regulation gate is fragile

The PFC PtP bus-regulation gate (`TestPfcReferenceDesignsPtp.cpp`,
`run_ptp_gates`, "Gate 1") measures a **still-settling 3-cycle transient** —
the voltage loop crosses over at ~fline/5 (≈10 Hz) and barely corrects in 3
line cycles, so the bus mean sits ~5% below nominal (NCP1654 lands 379 V vs
400 V target, just inside ±6%). It is **hypersensitive to sub-percent L
shifts**: a 0.17% inductance change moved the bus mean 379→370 V (tripping
±6%). That's why the interleaved PtP pins its per-phase L to the reference
value (3.30 mH) rather than the model's 0.17%-different computed value.

If you tighten or extend PtP, the right levers are **more line cycles** (better
settling, slower) or **feed-forward in the controller** — not nudging the gate
tolerance. Don't "fix" a failing bus gate by loosening it.

---

## 4. Pre-existing issues observed (NOT introduced here — do not attribute)

- `Test_MagneticAdviser_Inductor_Only_Toroidal_Cores` (TestMagneticAdviser.cpp
  ~line 1325) fails on a clean baseline too (verified by reverting all my
  files). It's an interaction with the committed "reject ferrite toroids for
  energy-storing topologies" change (7b1517b5) and/or the uncommitted adviser
  work — an inductor forced onto toroidal cores that all get rejected.
- Running the **whole `[adviser]` suite together** shows ~tens of failures +
  an occasional SIGSEGV. This is documented test-isolation flakiness (`[!mayfail]`
  tags reference "global settings state when run after Buck/Boost tests") and
  reproduces on baseline. Use focused tags (`[pfc-topology]`,
  `[from-converter]`, `[forward-topology]`, `[converter-model]`) for signal.

---

## 5. Where things live

- Model: `src/converter_models/PowerFactorCorrection.{h,cpp}`
  - `validate_topology_variant()` — the per-variant gate (accept/throw).
  - `is_buck_boost_class()` (header) — SEPIC/Ćuk classifier.
  - `calculate_duty_cycle()` — topology-aware duty.
  - `calculate_inductance_{ccm,dcm,crcm}()` — sizing (DCM throws for SEPIC/Ćuk).
  - `generate_ngspice_switching_circuit()` — `switch(variant)` allows unipolar
    boost family; uses `per_phase_power` so interleaved emits one cell.
  - three synthesis loops (all topology-aware now).
- Controller: `PfcControllerDesign.{h,cpp}`, `PfcControllerSubcircuits.h`.
- Tests: `tests/TestTopologyPowerFactorCorrection.cpp` (analytical/variant),
  `tests/TestPfcSwitchingNetlist.cpp` (netlist well-formed/rejection),
  `tests/TestPfcReferenceDesignsPtp.cpp` (`[slow]` ngspice PtP, harness with
  `variant`/`phases` fields + per-phase gating).

---

## 6. Verification commands

```bash
ninja -C build -j5 MKF_tests        # always -j5
build/MKF_tests "[pfc-topology]"                 # 30 cases incl 6 slow PtP
build/MKF_tests "[pfc-topology]~[slow]"          # fast subset (24)
build/MKF_tests "[totem-pole]"                   # totem-pole netlist + PtP only
build/MKF_tests "[converter-model]"              # 535 cases — shared-loop regression
build/MKF_tests "[pfc]"                          # PFC adviser path
```

All green as of this handoff: pfc-topology 30/30 (467 assertions, all 6 PtP
incl totem-pole), converter-model 535/535 (8407), pfc adviser 4/4.

---

## 7. Suggested next step
2a (totem-pole switching PtP) is **done** — every boost-family AND bridgeless
variant now has a converging switching PtP. The remaining switching-PtP gap is
**2b (SEPIC/Ćuk)**: a 4th-order plant (L1, L2, coupling cap, Cout) that needs
its own converging controller — a genuinely harder problem than totem-pole was
(totem-pole reused the boost loop verbatim; SEPIC/Ćuk cannot). Analytical
sizing for SEPIC/Ćuk is already done and correct for CCM/CrCM. If you take it,
expect the real work to be controller convergence, not the power stage. As with
totem-pole, keep the variant analytical-only rather than shipping a
non-converging circuit.
