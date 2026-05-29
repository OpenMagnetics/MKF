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
| TOTEM_POLE | ✅ (bipolar waveform) | ❌ throws | **needs bipolar power stage + controller** |
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

### 2a. Totem-pole switching PtP  (highest value remaining)
Totem-pole is analytically complete (bipolar inductor current/voltage already
synthesized in `process_operating_points`, `bipolar` flag). The gap is the
**switching netlist + PtP**. Why it's a real project, not an increment:
- No input bridge → the inductor sees a **signed sine**, current is bipolar.
- Needs a **bipolar power stage**: bipolar source `B_vin = vpk*sin(ωt)`,
  a totem-pole HF leg (high/low switch, complementary PWM, roles swap by
  half-cycle), and a line-frequency polarity leg.
- Needs a **bipolar rewrite of the average-current controller**: the
  multiplier reference `i_ref = G_mul·vea·vin/vrms_ff` must follow signed
  `sin` (not `|sin|`), and the current EA must regulate signed current. The
  existing controller (`PfcControllerDesign.{h,cpp}`) is unipolar and its
  convergence took heavy soft-start/anti-windup/IC-warm-start tuning — budget
  similar effort.
- Where to hook: `generate_ngspice_switching_circuit` currently has a
  `switch(variant)` that allows only the unipolar boost family; TOTEM_POLE
  falls into the `default: throw`. Add a `generate_ngspice_switching_circuit_totempole`
  (or branch) + a bipolar tuning path.
- PtP: extend `TestPfcReferenceDesignsPtp.cpp`. The gates compare the inductor
  envelope; for bipolar you'll compare against a signed-sine reference, not
  `|sin|`.

### 2b. SEPIC/Ćuk switching PtP
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
build/MKF_tests "[pfc-topology]"                 # 28 cases incl 5 slow PtP
build/MKF_tests "[pfc-topology]~[slow]"          # fast subset (23)
build/MKF_tests "[converter-model]"              # 533 cases — shared-loop regression
build/MKF_tests "[pfc]"                          # PFC adviser path
```

All green as of this handoff: pfc-topology 28/28 (415 assertions, all 5 PtP),
converter-model 533/533 (8355), pfc adviser 4/4.

---

## 7. Suggested next step
Take **2a (totem-pole switching PtP)** as a dedicated session — it's the last
boost-family gap and the highest-value remaining item. Start from the boost
netlist in `generate_ngspice_switching_circuit`, fork a bipolar power stage,
and adapt the controller multiplier/current-EA for signed current. Expect
real ngspice convergence iteration; keep totem-pole analytical-only (it
already works) if the bipolar netlist won't converge, rather than ship a
half-working circuit.
