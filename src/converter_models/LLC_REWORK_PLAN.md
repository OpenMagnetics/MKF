# LLC Rework Plan

Companion to `CONVERTER_MODELS_GOLDEN_GUIDE.md`. Brings the existing
`Llc.{h,cpp}` model up to DAB-quality and adds the missing rectifier
variants (most notably the **current doubler** the user asked about).

This is a planning document — no code has been changed. A future agent
or human can execute Tracks A and B end-to-end from this file alone.

---

## 0. Status snapshot

| Aspect | Current state | After this plan |
|---|---|---|
| Quality tier | Tier-2 (partial) | Tier-1 (DAB-quality) |
| Bridge variants | HB + FB ✅ | unchanged |
| Rectifier variants | CT only | **CT + FB + CD + VD** |
| ZVS diagnostics | none | `lastZvsMargin*` populated + tested |
| Header docstring | equations only | + ASCII art, disclaimer, disambiguation |
| SPICE primary snubbers | missing | RC across each MOSFET |
| SPICE solver options | `ITL1=300 ITL4=100` | `ITL1=500 ITL4=500` per guide |
| Multi-output warn-once | missing | mirrors `Dab.cpp:855–906` |
| Required tests | 21 cases, 3 missing | + `ZVS_Boundaries`, `RectifierType_*`, `Waveform_Plotting` |

The new `Cllllc` class identified during the survey (symmetric
Lr1+Cr1+Lm+Cr2+Lr2 tank, distinct from existing `Cllc`) is **out of
scope** for this rework. See § 7 for future-work notes.

---

## 1. Why this rework

The PSFB / PSHB rework already followed the Golden Guide. LLC is the
next-largest gap in the converter library and is the model the
magnetic-adviser path exercises most heavily for resonant designs. Two
specific drivers:

1. **Rectifier variants** — every authoritative LLC reference (TI
   SLUP263, IEEE 7930935, KIPE LLC-CD paper, MDPI Energies 17(24):6262)
   describes CT, FB, CD, and VD as four interchangeable secondary
   options on the *same* tank. PSFB / PSHB already implement
   `PsfbRectifierType` with three of these. LLC must catch up.
2. **DAB-quality polish** — ZVS-margin diagnostics, primary snubbers,
   ASCII art, accuracy disclaimer, `Waveform_Plotting` test. These are
   inexpensive and unblock the magnetic adviser from making safe LLC
   design decisions across the load range.

---

## 2. Variant taxonomy (industry survey)

Conclusions from a literature pass (TI SLUP263 / SLUA733, onsemi
AN-4151, Infineon LLC + 5.5 kW 3-phase + 11 kW CLLC app notes,
MDPI Electronics 8(1):3 / 12(7):1605, MDPI Energies 17(24):6262, IEEE
7930935 / 6980705 / 9124392, Wolfspeed CRD-22DD12N, Plexim PLECS demo).

### Belong as **rectifier_type enum values** on the existing class

| Variant | Diodes | Extra magnetics | Typical use |
|---|---|---|---|
| `CENTER_TAPPED` (default) | 2 | none | sub-1 kW, mid-V output |
| `FULL_BRIDGE`             | 4 | none | wide V_in, server PSU |
| `CURRENT_DOUBLER`         | 2 | 2 output inductors L1, L2 | high-I / low-V (≤12 V, ≥50 A) |
| `VOLTAGE_DOUBLER`         | 2 | 2 output capacitors | high-V / low-I (≥200 V) |

Tank, gain, and ZVS conditions are identical across the four — only
the secondary network changes. This mirrors `PsfbRectifierType` in
`PhaseShiftedFullBridge.cpp:34–61`.

### Belongs as **a new class** (deferred)

| Variant | Why distinct | Reference |
|---|---|---|
| `Cllllc` | Symmetric Lr1+Cr1+Lm+Cr2+Lr2 tank, fully bidirectional, different gain analytics | Wolfspeed CRD-22DD12N (22 kW OBC); MDPI Electronics 12(7):1605 |

### Skip

- 3-phase / interleaved LLC — model later as a `phaseCount` field on the
  same class (Infineon EVAL_5K5W_3PH_LLC_SIC2 = "3× HB-LLC" with
  current-sharing layered on top, not a new tank).
- LCC — 4th-order tank, only useful for HV/LED-driver applications.
- LCLC, ISOP-stacked — niche, no commercial pull.

---

## 3. Track A — DAB-quality polish

Mandatory checkpoints from `CONVERTER_MODELS_GOLDEN_GUIDE.md` § 3 / 4.6 /
5 / 8 that LLC currently fails.

### A.1 Header docstring (`Llc.h:13–39`)

Add three blocks before the existing equations section:

1. **ASCII schematic** — at minimum two diagrams: HB-LLC + CT, FB-LLC +
   FB-rectifier. Use `Dab.cpp:22–32` as the template. The four rectifier
   variants don't need separate ASCII; a textual table is sufficient.
2. **Accuracy disclaimer** — copy the structure from `Dab.cpp:296–320`
   verbatim and adapt:
   - Coss is fixed at design value; SPICE uses linear C, not Coss(Vds).
   - Body-diode reverse recovery is omitted (D model has RS but no Trr).
   - Gate-driver propagation delay ignored; dead time enforced ideally.
   - Magnetic core nonlinearity (μ vs B) not represented in the
     analytical solver; SPICE inductors are linear.
3. **Disambiguation** — one paragraph distinguishing:
   - **LLC** — series Lr + Cr, parallel Lm; resonant; passive rectifier
   - **SRC** — series Lr + Cr only, no parallel Lm
   - **DAB** — active rectifier on both sides, no resonant tank
   - **PSHB** — PWM half-bridge, not resonant; the Pinheiro-Barbi PSHB
     vs the AHB labelling caveat applies here too
   - **CLLC / CLLLC** — bidirectional resonant, distinct models

### A.2 Per-OP diagnostics (`Llc.h:59–62` + `Llc.cpp` solver)

Add to the private section of `Llc`:

```cpp
mutable double lastZvsMarginLagging = 0.0;     // A; >0 ⇒ ZVS achieved
mutable double lastZvsLoadThreshold = 0.0;     // A; load above which ZVS lost
mutable double lastResonantTransitionTime = 0.0; // s; commutation interval
mutable double lastPrimaryPeakCurrent = 0.0;   // A; peak iL_pri over period
```

Public read accessors `get_last_zvs_margin_lagging()` etc. — same prefix
convention as DAB.

Populate these at the end of
`process_operating_point_for_input_voltage` (around `Llc.cpp:1000` after
the steady-state solver converges):

- `lastPrimaryPeakCurrent` = max(|iL_pri|) across the sampled period.
  Already trivially available from the existing `Ipri_full` array.
- `lastZvsMarginLagging` = `Im_at_switching_instant − I_load_reflected`.
  Use the magnetizing current at the lagging-leg switching event minus
  the reflected load current at that instant. Reference: Huang
  SLUP263 § 6 ("ZVS condition"); same formula DAB / PSFB use.
- `lastZvsLoadThreshold` = `Vbus·dead_time / (4·Lm)` — the maximum
  load current the magnetizing energy can sustain ZVS for. Per
  Huang SLUP263 eq. (12).
- `lastResonantTransitionTime` = the dead-time interval, or the
  Coss-resonance time if Coss is provided.

### A.3 SPICE fixes (`Llc.cpp:1238–1271`)

1. **Primary-bridge snubbers**. After each `SW<i>` declaration, add:

   ```spice
   R_snub<i> <drain_node> <snub<i>_mid> 1k
   C_snub<i> <snub<i>_mid> <source_node> 1n
   ```

   Mandatory per Golden Guide § 5; both PSFB and DAB emit these.

2. **Solver options**. Change line 1270 from
   `.options ... ITL1=300 ITL4=100`
   to
   `.options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500`.

   The 500/500 limit is the guide-specified value and matches DAB.

### A.4 Multi-output warn-once (`Llc.cpp` `process_operating_points`)

Mirror `Dab.cpp:855–906`. At the top of `process_operating_points`,
after determining `numSecondaries`:

```cpp
if (numSecondaries > 1) {
    static thread_local bool warned = false;
    if (!warned) {
        std::cerr << "[Llc] Multi-output: per-secondary currents are "
                  << "approximated via load-share projection. Provide "
                  << "per-secondary leakage inductance for accurate "
                  << "individual current waveforms.\n";
        warned = true;
    }
}
```

### A.5 New tests (`tests/TestTopologyLlc.cpp`)

1. **`Test_Llc_ZVS_Boundaries`** — sweep load current at fixed Lm and
   verify `lastZvsMarginLagging` crosses zero at the predicted threshold
   `Vbus·dt / (4·Lm)`. Tolerance: ±10 %.
2. **`Test_Llc_Waveform_Plotting`** — write a CSV of (t, vab, ipri,
   isec_o1, vout_o1) to `tests/output/llc/waveforms_<rectType>.csv` for
   each rectifier variant. Not asserted, just generated for visual
   inspection. Use the existing `dump_waveform_csv` helper if present;
   otherwise inline a simple writer (≤ 30 lines).

---

## 4. Track B — `rectifier_type` enum (the "current doubler" answer)

### B.1 Schema change (`MAS/schemas/inputs/topologies/llcResonant.json`)

Note the file is camelCase (`llcResonant.json`), not snake_case.

Add to the top-level `properties`:

```jsonc
"rectifierType": {
    "description": "The type of secondary rectifier",
    "$ref": "#/$defs/llcRectifierType"
}
```

Add to `$defs`:

```jsonc
"llcRectifierType": {
    "description": "The type of secondary rectifier for LLC",
    "title": "llcRectifierType",
    "type": "string",
    "enum": [
        "centerTapped",
        "fullBridge",
        "currentDoubler",
        "voltageDoubler"
    ]
}
```

Do **not** make the field `required` — keep it optional and have the
C++ side default to `CENTER_TAPPED` via `get_rectifier_type().value_or(...)`.
This preserves backward compatibility with every existing test fixture
and persisted JSON file. (User confirmed this preference.)

After the schema edit, regenerate the C++ types so `MAS::LlcResonant`
exposes `get_rectifier_type()` returning
`std::optional<MAS::LlcRectifierType>`. Follow whatever code-gen step
the project already uses for the existing `psfbRectifierType` enum
(check `MAS/CMakeLists.txt` or the build script).

### B.2 Header (`Llc.h`)

Add four new mutable waveform vectors near the existing extras
(`Llc.h:68–72`):

```cpp
// CURRENT_DOUBLER only — two output inductors L1 / L2.
mutable std::vector<Waveform> extraLoCurrentWaveforms;
mutable std::vector<Waveform> extraLoVoltageWaveforms;
mutable std::vector<Waveform> extraLo2CurrentWaveforms;
mutable std::vector<Waveform> extraLo2VoltageWaveforms;

// VOLTAGE_DOUBLER only — Vout-cap voltage.
mutable std::vector<Waveform> extraVoutCapVoltageWaveforms;
```

Clear them at the top of `process_operating_points` alongside the
existing `extra*Waveforms.clear()` calls.

### B.3 Analytical secondary branching (`Llc.cpp:1013–1065`)

The existing function builds CT-shape secondary waveforms unconditionally.
Wrap the existing block in a `switch (rectType)` and add the three new
cases. Reference: `PhaseShiftedFullBridge.cpp:439–513` is the closest
pattern in the codebase.

```cpp
LlcRectifierType rectType =
    get_rectifier_type().value_or(LlcRectifierType::CENTER_TAPPED);

// Sinusoidal-ish secondary current shape comes from the Nielsen solver.
// Extract the rectified envelope first; then shape Vsec / iLo per type.
std::vector<double> Isec_envelope = ...;  // |reflected_pri_current|

switch (rectType) {
case LlcRectifierType::CENTER_TAPPED:
    // existing code path — two half-windings, unipolar pulses
    break;
case LlcRectifierType::FULL_BRIDGE:
    // one secondary winding, bipolar 3-level Vsec, full-wave rectified
    // current at output. Reverse diode voltage = Vsec_pk (not 2·Vsec_pk).
    break;
case LlcRectifierType::CURRENT_DOUBLER:
    // one secondary winding. Each output inductor carries Iout/2 DC
    // with ripple at 2·Fs (ripple cancels at the load summing node).
    // Emit extraLo*, extraLo2*. Reverse diode voltage = Vsec_pk.
    break;
case LlcRectifierType::VOLTAGE_DOUBLER:
    // one secondary winding, capacitor stack provides 2·Vsec_pk.
    // Emit extraVoutCap*. Reverse diode voltage = 2·Vsec_pk.
    break;
}
```

Each branch must populate the same set of MAS waveform fields the
existing CT code does (output_voltage, output_current, secondary
voltage / current per winding) so downstream MAS consumers see a
uniform shape.

### B.4 SPICE branching (`Llc.cpp:1245–1266`)

Branch the diode + output filter section. Reuse the netlist patterns
PSFB already emits (`PhaseShiftedFullBridge.cpp:660–730` is the
template).

Per-variant netlist deltas:

- **CT** — keep existing `D1_o<i> / D2_o<i>` + center-tap node.
- **FB** — emit four diodes `Dh1, Dh2, Dl1, Dl2` between
  `vsec_pos_o<i>` / `vsec_neg_o<i>` and the output rails. No center-tap
  node.
- **CD** — emit two diodes + two output inductors `L_o1_o<i>`,
  `L_o2_o<i>`. Each inductor needs its own initial-condition
  `.ic v(...) = ...` if pre-charge is desired. Sense source on each
  inductor for current readout: `V_Lo1_sense_o<i>`, `V_Lo2_sense_o<i>`.
- **VD** — emit two diodes + two output capacitors `C_o1_o<i>`,
  `C_o2_o<i>` with `.ic` setting the midpoint to Vo/2.

Sense sources (`Vsec_sense_o<i>`, `Vout_sense_o<i>`) keep the same
naming so `simulate_and_extract_*` doesn't need a parallel branch.

Update the `.save` directive to include the new sense nodes
(`i(V_Lo1_sense_o<i>)` etc.) when CD is selected.

### B.5 `get_extra_components_inputs` (`Llc.cpp:1482–1605`)

Currently returns the resonant cap (`CAS::Inputs`) and optionally the
external Lr (`Inputs`). Branch on `rectType`:

- **CT** / **FB** — unchanged
- **CD** — append two `Inputs` for L1 and L2. Average current = Iout/2,
  ripple amplitude derived from `Vo·(1−duty)/(Fs·Lo)` per
  `PhaseShiftedFullBridge.cpp:670–680`.
- **VD** — append two `CAS::Inputs` for the output capacitor stack.

### B.6 New tests

- **`Test_Llc_RectifierType_FullBridge`** — Vo within 5 % of analytical;
  secondary RMS = 0.707·Ipri_pk/n; reverse diode voltage = Vsec_pk
  (verify against the netlist max).
- **`Test_Llc_RectifierType_CurrentDoubler`** — Lo1 + Lo2 each carry
  ≈ Iout/2 DC; output ripple cancels at 2·Fs (FFT-based check or
  peak-to-peak of `Vout − Vo_avg` < threshold).
- **`Test_Llc_RectifierType_VoltageDoubler`** — Vo ≈ 2·Vsec_pk;
  output current ≤ 0.5·Iout_CT for the same secondary RMS.

---

## 5. File touchpoints summary

| File | Track A | Track B |
|---|---|---|
| `src/converter_models/Llc.h` | docstring, `lastZvs*` members | `extraLo*`, `extraLo2*`, `extraVoutCap*` |
| `src/converter_models/Llc.cpp` | populate `lastZvs*`, primary snubbers, `ITL=500/500`, multi-output warn-once | rectifier branching in waveform gen, SPICE, and `get_extra_components_inputs` |
| `MAS/schemas/inputs/topologies/llcResonant.json` | — | add `rectifierType` field + `llcRectifierType` enum |
| `tests/TestTopologyLlc.cpp` | `ZVS_Boundaries`, `Waveform_Plotting` | `RectifierType_FullBridge`, `RectifierType_CurrentDoubler`, `RectifierType_VoltageDoubler` |

Out of scope: any new file under `src/converter_models/`. No `Cllllc.*`
in this rework.

---

## 6. Verification

Per `feedback_build.md` (use `ninja -j2`, not `cmake --build`):

```bash
cd /home/alfonso/OpenMagnetics/MKF/build
ninja -j2 MKF_tests
./MKF_tests "[Llc]"
./MKF_tests "Test_Llc_ZVS_Boundaries"
./MKF_tests "Test_Llc_RectifierType_*"
./MKF_tests "Test_Llc_Waveform_Plotting"
./MKF_tests "Test_Llc_PtP_AnalyticalVsNgspice"   # NRMSE ≤ 0.15
```

If `MAS.hpp` complains about an undeclared `LlcRectifierType` after the
schema edit, copy the freshly-built `MAS.hpp` from the source tree into
`build/MAS/MAS.hpp` (per the `project_mas_hpp_staging` memory note).

Acceptance per Golden Guide § 14:

- [ ] Build clean
- [ ] All `TestTopologyLlc` cases pass (existing 21 + 5 new = 26)
- [ ] `Test_Llc_PtP_AnalyticalVsNgspice` NRMSE ≤ 0.15 still holds for
      every rectifier variant
- [ ] `lastZvsMarginLagging` changes sign at the predicted Lm in
      `Test_Llc_ZVS_Boundaries`
- [ ] `rectifierType` field present in `llcResonant.json`; default
      behaviour (no field) still produces CT result
- [ ] Header docstring includes ASCII art, accuracy disclaimer,
      LLC/SRC/DAB/PSHB disambiguation
- [ ] SPICE netlist contains primary-bridge `R_snub` + `C_snub` and
      `ITL1=500 ITL4=500`
- [ ] Schema regeneration emits a `MAS::LlcRectifierType` enum
      consumed by `Llc::get_rectifier_type()`

Update `project_converter_model_quality_tiers.md` once Tracks A + B
land: LLC moves from Tier-2 to "DAB-quality with all four rectifier
variants supported."

---

## 7. Future work (not part of this rework)

- **`Cllllc` class** — symmetric tank Lr1+Cr1+Lm+Cr2+Lr2; bidirectional
  EV/server applications. Mirror `Cllc.{cpp,h}` structure. References:
  Wolfspeed CRD-22DD12N (22 kW OBC), MDPI Electronics 12(7):1605.
  Estimated effort similar to the original `Cllc` introduction.
- **3-phase / interleaved LLC** — add a `phaseCount` field on the
  existing `Llc` class (defaults to 1). Each phase shares the tank
  formulae; current-sharing analysis goes in a helper. Reference:
  Infineon EVAL_5K5W_3PH_LLC_SIC2.
- **Synchronous-rectifier variants** — could become a `synchronous: bool`
  flag orthogonal to `rectifierType`. SR replaces diodes with MOSFET
  drops; analytical Vo gains slightly, SPICE becomes more complex.

---

## 8. References

- TI SLUP263 — Hong Huang, *Designing an LLC Resonant Half-Bridge Power
  Converter*
- TI SLUA733 — *LLC Design for UCC29950*
- onsemi AN-4151 — *Half-Bridge LLC Resonant Converter Design*
- Infineon — *Resonant LLC Converter: Operation and Design*
- Infineon EVAL_5K5W_3PH_LLC_SIC2 — 5.5 kW 3-phase interleaved LLC
- Infineon — *Operation and modelling analysis of a bidirectional CLLC*
- Wolfspeed CRD-22DD12N — 22 kW Bidirectional CLLLC reference design
- IEEE 7930935 — Nigsch & Schlenk, *Detailed analysis of a current-doubler
  rectifier for an LLC resonant converter*
- IEEE 6980705 — *LLC resonant converter with current-doubler rectification*
- IEEE 9124392 — *Three-Phase Interleaved LLC for Telecom and Server*
- KIPE — *Mathematical Analysis of LLC Series Resonant Converter with
  Current Doubler Rectifier*
- MDPI Electronics 8(1):3 — *Resonant Converter with Voltage-Doubler or
  Full-Bridge Rectifier for Wide-Output Voltage*
- MDPI Energies 17(24):6262 — *A New Voltage-Doubler Rectifier for
  High-Efficiency LLC Resonant Converters*
- MDPI Electronics 12(7):1605 — *Bidirectional CLLLC Resonant Converter*
- ITACOIL — *LLC and LCC Resonant Topologies*
