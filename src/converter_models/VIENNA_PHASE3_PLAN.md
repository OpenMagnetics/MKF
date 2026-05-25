# Vienna Rectifier — Phase 3+ Implementation Plan

**Repo:** `OpenMagnetics/MKF`
**File:** `src/converter_models/Vienna.{cpp,h}`
**Written:** 2026-05-25
**Status:** Item 1 (fullLineCycle) landed in commit (TBD).
Items 2-6 below are next.

This document is the spec for finishing Vienna Phase 3+. Each item is a
self-contained workstream — pick one per session/PR. Order is suggested
(low risk → high risk), but not load-bearing.

---

## Reference papers (read once, then keep handy)

- **[K1994]** J. W. Kolar, F. C. Zach, "A Novel Three-Phase Three-Switch
  Three-Level Unity Power Factor PWM Rectifier", PCIM 1994 —
  original Vienna paper, all base formulas.
- **[H2011]** R. Hartmann, "Ultra-Compact and Ultra-Efficient Three-Phase
  PWM Rectifier Systems for More Electric Aircraft", ETH dissertation
  19755, 2011 — closed-form switch RMS / diode RMS for Vienna I, §3.2.
- **[FK2014]** T. Friedli, J. W. Kolar, "The Essence of Three-Phase PFC
  Rectifier Systems Part II", IEEE TPEL 29 (2), 543-560 — Vienna II
  variant, alternative switch arrangements, comparative loss model.
- **[KMG2019]** J. Kolar et al., "Vienna Rectifier and Beyond", APEC
  2018-2019 plenary slides — modern SiC implementations, sync rec.
- **[TIDA-010257]** TI 10 kW Vienna PFC ref design — concrete
  interleaving (phaseCount > 1) design example.

All exist as PDFs; URLs in `VIENNA_PLAN.md`.

---

## Item 2 — `samplingStrategy = peakOfLinePlusSectors`

### Goal

Six DPWM sector samples (every 60° of the line cycle) plus the peak-of-
line sample, in a single OperatingPoint with 3 windings × 7 samples each.
Adviser sees worst-case across sectors. Lighter than `fullLineCycle`
(7 snapshots vs. a 4096-point envelope) and more accurate than
`peakOfLineOnly` (catches sector edges where the duty / ripple shape
deviates from the sinusoidal model).

### Physics

For each line angle θ ∈ {30°, 90°, 150°, 210°, 270°, 330°} (sector
mid-points) plus the peak (90°, already covered by `peakOfLineOnly`),
recompute the local boost duty d(θ) = 1 − M·|sin θ| and the local
switching-period triangular ripple ΔI_pp(θ) = V_phase·|sin θ|·d(θ)/(L·Fsw).
Emit one switching-period waveform per sector per phase.

### Code sketch

```cpp
case PEAK_OF_LINE_PLUS_SECTORS: {
    const double sectorAngles[6] = {
        M_PI/6, M_PI/2, 5*M_PI/6, 7*M_PI/6, 3*M_PI/2, 11*M_PI/6,
    };
    for (int ph = 0; ph < 3; ++ph) {
        for (double θ : sectorAngles) {
            double Vphase_θ = V_peak * std::sin(θ + phaseOffsets[ph]);
            double duty_θ   = 1.0 - std::abs(Vphase_θ) / Vhalf;
            double ΔI_θ     = std::abs(Vphase_θ) * duty_θ * T_sw / L;
            // Emit a switching-period waveform with these params, tagged
            // "Phase A sector 1" / "Phase A sector 2" / ...
        }
    }
}
```

### Decisions / open questions

- **Output shape**: 3 windings × 7 sectors = 21 excitations per OP, OR
  21 separate operating points each with 3 windings? The latter aligns
  with how `process_operating_points` produces multiple OPs for input-
  voltage sweeps (PSFB/PSHB pattern) and the adviser handles `vector<OP>`
  natively. Recommended.
- **Naming**: `"Phase A @ 30°"`, etc. Wizard label needs to handle this.

### Test

`Test_Vienna_PeakOfLinePlusSectors_EmitsSevenSectorsPerPhase`
- Replay default JSON with `"samplingStrategy": "peakOfLinePlusSectors"`.
- Assert 7 sector-OPs returned, each carrying 3 windings.
- Assert ΔI_pp at θ=30° is roughly 0.5·ΔI_pp(peak) (since sin(30°)=0.5).

### Estimate

~150 lines + 50 lines test. Half a session.

---

## Item 3 — `viennaVariant = viennaII`

### Goal

Vienna II uses **two switches per phase leg** (bidirectional clamp made
of back-to-back MOSFETs) instead of Vienna I's single switch + bidir
diode bridge. Differences: switch V-stress is still Vdc/2, but switch
RMS / average current and switching loss split differently between the
two devices per leg. Diode count drops from 6 fast + 6 line to 6 line.

### Physics — switch currents per leg

Vienna II [FK2014 §IV.B]:

- Each of the two MOSFETs per leg conducts only during ONE half-cycle of
  the line current, gated by the sign of phase current.
- Per-switch RMS = `I_pk · sqrt(1/8 − M/(3·π))` (half the conduction
  window of Vienna I → factor sqrt(1/2) below the Vienna I form).
- Per-switch average = `I_pk · (1/(2π) · ∫₀^π sin(t)·(1 − M·sin(t)) dt)`.
- No fast-rectifier diode current (rectification is via the MOSFET body
  diode of the OPPOSITE-polarity device when both are off — sync rec
  variant covered in Item 5).

### Code changes

- New helper: `compute_switch_rms_viennaII(I_pk, M)`,
  `compute_switch_avg_viennaII(I_pk, M)`.
- Branch in `process_operating_point_for_input_voltage`:
  ```cpp
  bool isViennaII = (get_vienna_variant() == ViennaVariant::VIENNA_II);
  lastSwitchRmsCurrent = isViennaII
      ? compute_switch_rms_viennaII(I_pk, M)
      : compute_switch_rms(I_pk, M);
  ```
- `run_checks`: drop `if (viennaVariant != VIENNA_I) throw`.

### Test

`Test_Vienna_ViennaII_SwitchRmsClosedForm`
- Compare against [FK2014] Table III data points (5 kW, 10 kW design
  examples, V_LL = 400, V_dc = 800).

### Estimate

~250 lines + 100 lines test. One session.

---

## Item 4 — Alternative switch types

### `backToBackMosfet`

Same physical arrangement as `viennaII` but with explicit
back-to-back MOSFETs (drain-drain or source-source) per leg. Reuses the
Vienna II formulas from Item 3. The MAS schema distinguishes them for
SPICE / BOM purposes but the analytical model is identical.

### `singleMosfetIn4DiodeBridge`

A single MOSFET inside a four-diode bridge per leg. The MOSFET sees the
full inductor current (both polarities) → conduction loss model uses
the FULL line-current RMS, not the half-cycle form:

- MOSFET RMS = `I_pk · sqrt(1/4 − 2M/(3π))` (same as Vienna I single-
  switch leg, but with 4 forward diode drops in series in the
  conduction path).
- 4 diodes carry MOSFET-RMS each.

Higher conduction loss (4× Vf) but lowest switch count. Used when MOSFET
cost dominates BOM.

### Test

`Test_Vienna_SingleMosfetIn4DiodeBridge_DiodeCount`
`Test_Vienna_BackToBackMosfet_MatchesViennaII`

### Estimate

~150 lines + 50 lines test (mostly delegating to Item 3 formulas).
Half a session.

---

## Item 5 — `synchronousRectifier`

### Goal

Replace the boost diodes with active MOSFETs operated in synchronous
rectification mode. Lower conduction loss (Rds·I² instead of Vf·I_avg)
at the cost of additional gate-drive complexity.

### Physics

- Diode conduction loss model:  `P_diode = Vf · I_diode_avg`
- Sync-rec MOSFET loss model:   `P_swSR = Rds_on · I_diode_rms²`
- Crossover: sync rec wins when `Rds·I_rms < Vf`, i.e. at high current.
- I_diode_avg already computed (`I_pk / π`).
- I_diode_rms needs a new helper:
  `I_diode_rms = I_pk · sqrt(1/4 + M/(3π))` [H2011 §3.2.3].

### Code changes

- New diagnostic field: `lastDiodeRmsCurrent`.
- New helper: `compute_diode_rms(I_pk, M)`.
- `run_checks`: drop the `synchronousRectifier` gate.
- The user-supplied Rds_on isn't part of the schema today — either add
  it (preferred; mirror PSFB's MOSFET Rds_on slot) or compute the
  diagnostic with `Rds_on = 0` and let the downstream loss model
  multiply by a user-supplied value.

### Test

`Test_Vienna_SyncRec_DiodeRmsClosedForm`
- Spot-check `I_diode_rms` at the design point against [H2011] Table 3.2.

### Estimate

~150 lines + 50 lines test. Half a session.

---

## Item 6 — `phaseCount > 1` (interleaved channels)

### Goal

N parallel boost stages per phase, gated 360°/N out of phase on the
switching scale (NOT line scale). Each channel carries I_phase/N average
and has its own inductor with peak-to-peak ripple ΔI_pp/N. By ripple
cancellation at the SUM node, the input-side observable ripple drops by
roughly N² (in the ideal case).

### Physics — per-channel inductor sizing

Inductor for N-channel interleaving:

- Each channel inductor: L_ch = L_single / N (same di/dt requirement,
  N× the switching events per line cycle on the sum node).
- Per-channel I_avg = I_pk_phase / N
- Per-channel ΔI_pp = V_phase·d·T_sw / L_ch = N·(V_phase·d·T_sw / L_single)
  — so each channel's ripple is N× larger than the single-channel case,
  but the SUM ripple at the converter input cancels to ~1/N² of single.

The adviser sizes ONE channel inductor — it's smaller and carries less
DC than the single-phase equivalent, but with higher ΔB.

### Output shape

3 phases × N channels = 3N windings per operating point. Naming:
`"Phase A ch 1"`, `"Phase A ch 2"`, ..., `"Phase C ch N"`.

### Code changes

- `process_operating_point_for_input_voltage`: loop `for (ch = 0; ch < N;
  ++ch)` inside the per-phase loop, each channel offset by `2π·ch/N` of
  the switching period (NOT line period).
- The line-cycle envelope is the same across channels of a given phase;
  only the switching ripple phase shifts. This is invisible in the
  current sub-sampled visualisation but visible in FFT.
- `process_design_requirements`: divide computed L by N.
- `run_checks`: drop the `phaseCount > 1` gate.

### Test

`Test_Vienna_PhaseCount_3_EmitsNineWindings`
`Test_Vienna_PhaseCount_2_InductorIsHalfSingleChannel`

### Estimate

~200 lines + 100 lines test. One session.

---

## Estimated total for items 2-6

- Item 2: 0.5 session
- Item 3: 1 session
- Item 4: 0.5 session
- Item 5: 0.5 session
- Item 6: 1 session

**Total: ~3.5 sessions, sequential.**

Items can run in parallel if multiple agents pick them up — they touch
overlapping but non-conflicting code paths (each adds its own
`run_checks` exemption, new helpers, and one branch in
`process_operating_point_for_input_voltage`).

---

## Cross-cutting work (do once, after items 2-6)

### Full 3-phase SPICE netlist

Replace the current single-phase-boost emulation in
`Vienna::generate_ngspice_circuit` with a true 3-phase circuit:
3 phase voltage sources at 120° offset, 3 boost inductors, T-type
switch tree per phase (or viennaII back-to-back per phase), 6-diode
front-end rectifier, split DC bus with neutral-point ground reference.
~400 lines. The single-phase emulator stays for fast adviser pre-
screening; the full netlist is gated on `sampling = fullLineCycle` or
`phaseCount > 1`.

### Neutral-point voltage balance

Vienna's split bus needs a NP voltage balancer (in real designs an
outer control loop; in the model a constraint that the upper and lower
half-bus capacitor average currents match). Currently assumed balanced.
Required for accurate cap RMS / lifetime sizing.

### DC-bus capacitor sizing

Add to `process_design_requirements`: target ripple voltage on Cdc/2
based on third-harmonic-of-line ripple injected by the 6-pulse
rectification. Standard formula:
`Cdc ≥ I_dc·D / (Δv_dc · 6·F_line)`.

---

## Wizard wiring (already complete)

`WebFrontend/src/components/Wizards/ViennaWizard.vue` already exposes
the controls for all of the above (`viennaVariant`, `switchType`,
`samplingStrategy`, `synchronousRectifier`, `phaseCount`). Today they
pass through to MKF, which throws on anything beyond the Phase 1+2
defaults. As each item lands above, the wizard option becomes live with
no further frontend change.
