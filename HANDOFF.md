# MKF Converter-Models DAB-Quality Push — Handoff

## Goal

Bring all non-DAB-quality converter models in MKF up to the DAB-quality
standard defined in `src/converter_models/CONVERTER_MODELS_GOLDEN_GUIDE.md`
§15 (primary-current NRMSE ≤ 0.15 vs SPICE on 3 reference designs, ZVS test,
rectifier-aware, multi-output, etc.). User directive: **correct physics >
passing tests; SPICE is the source of truth, analytical must match SPICE;
look online if needed**.

## Hard rules (CLAUDE.md)

- No silent shortcuts / fallbacks / magic-number defaults — **throw**.
- NEVER force-push to main. NEVER `git add -A/.` — explicit paths only.
- NEVER `git stash`.
- **Concurrent compile rule**: `pgrep -fa "ninja -j"` before any build;
  don't build if another ninja is anywhere on the box.
- **Don't touch** other-session dirty files: `Cllc.cpp`, `TestCllc.cpp`.
- One commit per logical phase; explicit `git add <paths>`;
  push via `git -c core.sshCommand="/usr/bin/ssh -i ~/.ssh/hephaestus_om" push origin main`.
- Build: `cd build && ninja -j6 MKF_tests` foreground.

## This session — accomplishments

### SRC SPICE buildout (May 2026) — Phase-1 DAB-quality

Brought SRC from "Phase 1+2 analytical only, 0 SPICE references" to
"DAB-quality tier-1 (FB-diode rectifier scope)" in one session.

**Final state — 3 PtP designs, all green (69 assertions):**

| Design | Pin (W) | Loss % | NRMSE | Peak ratio (spice/anal) |
|---|---|---|---|---|
| Telecom-500W-FB-AboveFr (400→48V/10A, fr=100kHz, fs=110kHz, Q=1) | 505 | 18.7 | 8.3 % | 1.045 |
| Telecom-250W-HB-AboveFr (400→48V/5A, fr=100kHz, fs=110kHz, Q=1)  | 253 | 17.9 | 8.4 % | 1.096 |
| Telecom-1k4W-FB-AboveFr (400→48V/30A, fr=100kHz, fs=120kHz, Q=0.5)| 1516 | 25.4 | 12.6 % | 0.975 |

All NRMSE ≤ 0.15 gate. All peak ratios within [0.5, 1.7] band (much
tighter than CLLC's ~1.35; SRC has no Lm-branch artifact in FHA).
Combined wall time < 1 s.

**Changes:**
- `src/converter_models/Src.h` — added `numPeriodsToExtract`,
  `numSteadyStatePeriods`, `extra{Cap,Ind}{Voltage,Current}Waveforms`
  state; declared 3 SPICE methods.
- `src/converter_models/Src.cpp` — implemented
  `generate_ngspice_circuit` (behavioural-PULSE bridge, HB+FB,
  FB-diode rectifier, separate Lr, pairwise K matrix at K=0.999),
  `simulate_and_extract_operating_points`,
  `simulate_and_extract_topology_waveforms`. Mirrors LLC's
  BEHAVIORAL_PULSE branch verbatim; ~280 LOC.
- `tests/TestSrcReferenceDesignsPtp.cpp` — NEW. 3 ref designs, same
  gate pattern as `TestCllcReferenceDesignsPtp.cpp` (Vin sanity,
  power balance, NRMSE, peak amplitude, wall time).

**Documented deferrals:**
- **FIXME-src-1** (rectifier coverage): `SRC_PLAN.md` §6 designs
  (Steigerwald plasma, Telecom 3 kW CD, X-ray HV CT) use CT/CD
  rectifiers. Phase-1 SPICE is FB-diode only; analytical solver
  already throws on non-FB rectifiers. Adding CT/CD requires
  branches in both `process_operating_point_for_input_voltage`
  and `generate_ngspice_circuit`. Until then the PtP suite uses
  FB-diode variants.
- **`get_extra_components_inputs`** for Cr/Lr: state members are
  in place but `process_operating_point_for_input_voltage` does
  not yet populate `extra*Waveforms`, and `get_extra_components_inputs`
  itself is not yet implemented. Not blocking PtP suite. Mirror
  `Llc.cpp:1934-2143` when needed.

### CLLC polish (May 2026) — committed & pushed

After the prior CLLC DAB-quality session (commit `5f075d73`), a 6-item polish
list was tackled. Snapshot of final state for all 6 PtP designs:

| Design | Pin (W) | Loss % | NRMSE | Peak ratio (spice/anal) |
|---|---|---|---|---|
| Telecom-500W (at fr) | 578 | 17.6 | 5.6 % | 1.345 |
| Telecom-500W-BelowFr | 560 | 15.4 | 4.4 % | 1.354 |
| Telecom-250W (above fr) | 371 | **34.5** | 4.6 % | 1.357 |
| Infineon-11kW | 11491 | 0.9 | 5.4 % | 1.352 |
| KIT-20kW-Asymmetric | 21464 | 11.3 | 4.5 % | 1.367 |
| Telecom-500W-REVERSE | n/a | n/a | **8.8 %** | 0.895 |

All 132 assertions in 6 cases pass. Highlights:

- **`98f7d470`** `test(cllc): add SPICE loss-breakdown diagnostic for KIT-20kW`
  — `tests/TestCllcKitLossBreakdown.cpp`, tagged `[.diag][cllc-loss-breakdown]`.
- **`556301d1`** `fix(cllc): add primary body diodes + shrink Cout — phantom snubber loss 23 kW → 0.5 kW`
  — `.model DBODY`, 4 antiparallel diodes on primary only (NOT SR), Cout
  100 µF → 10 µF. KIT phantom loss 23 060 W → 589 W. NRMSE 15.4 % → 4.5 %.
- **`c6b0fccb`** `test(cllc): add peak-amplitude gate to PtP suite`
  — `ptp_nrmse` is mean-subtracted AND RMS-normalized, so amplitude scale
  errors register as only ~4-5 % NRMSE: shape-only. Added explicit
  `tol_peak_ratio_min/max` per `RefDesignSpec`, default [0.5, 1.7].
  Catches the FIXME-item1 regression band.
- **`41db0ff0`** `test(cllc): T_common resampling + diagnostics for reverse-mode PtP`
  — applied the same time-horizon resampling fix that forward had; only
  moves reverse NRMSE 8.88 → 8.85 % though (residual is real shape gap).
- **`2260bca1`** `docs(cllc): document FIXME-item3 (Telecom-250W above-fr SR snubber loss)`
  — diagnosed cause (SR sync-gated with primary PWM; above-fr the SR
  opens with non-zero secondary current that gets dumped into the 1k
  snubbers, 6× under-delivery on secondary). Attempted SR-body-diode
  fix turned reverse into 171× explosion; reverted.

### Deferred (documented in-file as FIXME-itemN)

- **#1** KIT primary peak ±53.8 A vs analytical ±38.7 A (systematic 1.35-1.39×
  across ALL 5 forward designs). Stiff-Vout experiment proved diagnosis
  (drops peak to ±33.5 A, collapses Pin 16× → 1.34 kW). Means at the
  design fsw the tank only naturally delivers ~7 % of Pnom into Vnom; the
  Cout+Rload sim hides this by letting Vout drift. **Fix belongs in
  `Cllc.cpp` solver** (`process_operating_point_for_input_voltage` or
  `CllcResonantParameters` derivation), not the testbench.
- **#2** Telecom-500W-REVERSE 8.85 % NRMSE. Likely reverse-mode TDA seed-bias
  (`Cllc.cpp:1158-1168`); under 15 % gate.
- **#3** Telecom-250W 34 % loss. Requires current-direction-aware SR gating
  (B-source conditional on i_sec polarity). Under 60 % gate.

### Other items (already complete or non-issues)

- **#4** DAB body-diode fix: **already in place**. `Dab.cpp:1065-1074` and
  `:1162-1171` have antiparallel `D1..D8 DIDEAL` diodes on every switch.
  All 2444 DAB assertions pass.
- **#5** This HANDOFF.md update.
- **#6** `-Wmaybe-uninitialized` at `Cllc.cpp:1567` — gcc 13/14 false
  positive from `std::optional<Cooling>` storage byte-wise traversal
  during `set_conditions(conditions)` copy. The optional is explicitly
  set to `nullopt`. Suppressed with a localized `#pragma GCC diagnostic
  push/pop`.

---

## Previous session — accomplishments

### Committed & pushed
- **`477773a4`** `feat(psfb): replace MOSFET-RON freewheel tau with diode-dynamic-R model`
  - PSFB freewheel τ changed from `Lr/(2·RON_mosfet)` ≈ 250 µs to
    `Lr/(n²·R_d_diode)` ≈ 1.7 µs (Vlatkovic 1992 §III asymmetric-rectifier-
    commutation; reflected diode dynamic resistance, not MOSFET RON).
  - Improvement scales as 1/n²: Telecom (n=20) gains most, EV-Aux (n=5) barely
    affected.
  - NRMSE (3 ref designs):
    - Telecom 12V/50A : 22.9 % → **15.8 %** (−7.1 pp)
    - Server  24V/50A : 15.4 % → **13.5 %** (−1.9 pp)
    - EV-Aux  48V/21A : 17.0 % → **16.6 %** (−0.4 pp)
    - Single-OP 12V/50A: 11.3 % → 16.2 % (+4.9 pp — exposes residual
      active-interval commutation undershoot that the old flat-FW model was
      coincidentally compensating for)
    - Multi-output uniform: ~0.22 → 12.6 %
    - Multi-output 2:1 spread: → 12.5 %
  - Gates tightened: 0.30/0.25 → 0.20 (0.40 → 0.25 for 2:1 spread).
  - 168/168 PSFB assertions / 24 tests pass.
  - Files: `src/converter_models/PhaseShiftedFullBridge.cpp` (freewheel
    block at ~lines 478-510), `tests/TestPhaseShiftedFullBridgeReference
    DesignsPtp.cpp` (3 gates), `tests/TestTopologyPhaseShiftedFullBridge.cpp`
    (3 gates), `src/converter_models/FUTURE_TOPOLOGIES.md` (status row).

### Diagnosed but not fixed (deferred)

**PSFB last 1-2 pp to DAB-quality (≤ 0.15)** — three ref designs cluster at
13.5–16.6 %, mean 15.3 %. Residual gap is in the **active-interval iSec
winding amplitude**:
- Analytical iSec winding range: ±54 A
- SPICE iSec winding range: ±67 A (24 % bigger)
- Hypothesis: SPICE's `iSec_winding = iLo + Im·n` (correct transformer KCL)
  while analytical's `iSec_winding = sign(Vpri) · iLo` (omits Im·n).
- For n=20, Im·n ≈ 8 A — explains part of 13 A gap. Remaining ~5 A is real
  PSFB rectifier asymmetry during active interval (current redistributes
  between diodes due to Lr, making transient peaks higher than steady-state
  ILo). Requires Vlatkovic-style multi-state ZVS-transition model
  (Sadhana 2025 Eq. 11a/12 — at least one full session).

**Two candidate fixes for next session:**
1. **Quick (~20 min)**: add `iSec += Im·n` to secondary winding model
   in PSFB.cpp lines ~595-608 (the `else` branch for FB/CD rectifier).
   Affects winding 1 RMS/copper-loss; may or may not improve gated winding-0
   NRMSE.
2. **Deep (multi-session)**: full Vlatkovic 6-state ZVS model with sinusoidal
   LC during transitions (Sadhana 2025 Eq. 11a/12).

## Topology queue — priority order

1. **PSFB** — at 13.5-16.6 % NRMSE, gate 0.20. Within 1-2 pp of DAB-quality.
   *Defer further work* unless user wants the iSec KCL fix or full Vlatkovic.
2. **Flyback** — DAB-quality (commit `733ca7f3`).
3. **PushPull** — DAB-quality (commit `6ecab010`).
4. **IsolatedBuck** — DAB-quality (this session). NRMSE Design1 12.0% / Design2 3.3% / Design3 2.3% (gate 0.15). Root cause was double: (a) old analytical primary current was a symmetric triangle that ignored the off-time secondary-reflected step; (b) SPICE used K=0.99 coupling which broke secondary-diode forward-bias for low-V secondaries (Vout_sec ≤ Vd_diode). Fix: SPICE K=1 (Kpri_sec and Ksec_sec, all coupled inductors), analytical primary current rewritten as PIECEWISE 4-sample waveform `(0, D·T, D·T, T) → (i_mag_min, i_mag_max, i_mag_max − I_off, i_mag_min − I_off)` where `I_off = Σ Iout_sec/((1−D)·n)` from Cout_sec charge balance. Magnetizing current `i_mag = i_pri + Σ i_sec/n` is the continuous symmetric-triangle component; winding currents step at switching events (intrinsic to K=1 coupled-L). All 24 IsolatedBuck test cases / 237 assertions pass.
5. **IsolatedBuckBoost** — DAB-quality (this session). NRMSE Design1 1.57% / Design2 0.82% / Design3 1.03% (gate 0.15). Root cause was topology wiring: prior `Dpri pri_in→vpri_rect` with `Cpri IC=+Vout` left the primary diode permanently reverse-biased because V(pri_in) goes NEGATIVE during OFF in flyback operation. Fix per TI SNVAA84 "Fly-Buck-Boost" (negative-Vout variant of Fly-Buck): reversed Dpri to `vpri_rect → pri_in`, `Cpri IC = -|Vout|`, `.nodeset v(vpri_out)=-|Vout|`; analytical primary-current average now `(Iout_pri + Σ Iout_sec/n) / (1-D)` per Erickson §5.2 buck-boost charge balance (was missing /(1-D)); primary voltage min corrected to `-(Vout_pri + Vd)`. Continuous Bimag probe gives clean triangle (no piecewise needed unlike IsolatedBuck). All 20 IsolatedBuckBoost test cases / 156 assertions pass.
6. **Forward × 3 (SSF / TSF / ACF)** — DAB-quality with caveat (this session). NRMSE:
   - SSF: UC3845-50W 11.6 % / NCP1252-100W 13.0 % / Erickson-25W 11.6 % (gate 0.15).
   - TSF: SLUP089-100W 14.3 % / UCC2897-50W 13.7 % / Erickson-25W 13.6 % (gate 0.15).
   - ACF: UCC2897A-100W **26.0 %** / Erickson-50W 11.8 % / AN1023-200W **26.8 %** (gate 0.15 Erickson, 0.30 others).
   Two fixes:
   - **Volt-second `/2` removal**: SSF/TSF/ACF `process_operating_points_for_input_voltage` had `t1 = period/2 · (Vout+Vd)/(Vin/n)` (left over from push-pull/half-bridge thinking).  Forward has NO inherent /2 — volt-second balance on Lout: `Vout = Vin/n·D − Vd ⇒ D = (Vout+Vd)/(Vin/n)`.  Matches SPICE netlist `tOn = period * dutyCycle`.  Drops SSF/TSF NRMSE 74-76% → 11-14% and ACF 89-98% → 46-64%.
   - **Shape-match OP for ACF**: open-loop SPICE drives D=Dmax fixed; Vout settles below nominal.  `build_for_shape_match` in TestActiveClampForwardReferenceDesignsPtp.cpp uses `Vout_for_t1 = Dmax·Vin/n − Vd` (forces analytical t1 = SPICE tOn), `Iout_spice = Vout_spice² / Rload_nom`, and recomputed `current_ripple_ratio` from physical Lout at the SPICE-settled OP.  Vin range collapsed to nominal to avoid Vin_min triggering the `t1 > T/2` guard at the elevated Vout_for_t1.  Drops Erickson to 11.8% (DAB-quality).
   - **Residual on UCC2897A and AN1023** (high N or high Iout): SPICE primary current collapses to ~0 during OFF instead of ramping +Im/2 → −Im/2 through the active-clamp leg.  The magnetizing current commutates via the secondary-side rectifier (small Lsec=Lm/n² and K=0.9999 leakage favor that path over the active-clamp Cclamp/S_clamp loop).  Analytical is correct; SPICE netlist needs fix (tighter K, body-diode on S_clamp, or Cclamp-referenced-to-Vin standard ACF topology).  Gate held at 0.30 to keep the regression honest — documented in test header.
7. **PFC** — SKIPPED this session per user decision. Status: NCP1654-100W envelope NRMSE **7.5 %** (DAB ✓); UCC28180-360W **22.2 %** (gate 30 %); L4981-1000W **51.4 %** (gate 60 %). All 3 PtP tests PASS at current gates. Root cause investigated and is **not** the di/dt ZC cusp (cusp width = `2L·Iin_pk/Vpk ≈ 25 µs` on L4981, <5 % contribution). Dominant error is SPICE controller behavior: iL stays near zero for ~500 µs after each line ZC (vs designed `τ = 1/(2π·fc_i) ≈ 32 µs` for `fc_i = fsw/10 = 5 kHz`), then catches up and overshoots peak (iL_peak 7.5 A vs ideal 6.15 A, +22 %). Physical bound check at t=7.4 ms post-ZC: `i_max ≈ 4.18 A`, SPICE `iL_avg = 0.32 A` — well below di/dt limit, so it's the controller not the inductor. Next agent who picks this up should investigate `PfcControllerDesign.cpp` (likely `ic_vc_i`, FF clamp, or current-EA topology) and/or `PfcControllerSubcircuits.h` to find why the 5 kHz current loop responds at 500 µs not 32 µs. Use the env-gated CSV dump added this session (see Useful invocations).
8. CLLC rewrite (NB: `Cllc.cpp` + `TestCllc.cpp` are other-session dirty —
   do not touch until coordinated)
9. CLLLC — **BLOCKED**: `AdvancedClllc::process` at
   `src/converter_models/Clllc.cpp:1573` is a stub that throws "not yet
   implemented. Depends on Clllc::process_operating_points." — underlying
   MKF logic doesn't exist yet, only the wrapper. Cannot work on CLLLC
   until CLLC (§8) is implemented.
10. Vienna (NB: SPICE netlist missing — `Src.cpp` 0 SPICE refs, `Vienna.cpp` only docstring; will require ground-up SPICE rewrite before DAB-quality applies).
11. SRC (same — analytical only, no SPICE).

**Campaign status**: DAB-quality push is now effectively stalled. CLLC/CLLLC blocked, Vienna/SRC need SPICE rewrites, PFC needs controller deep-dive. Resume requires either (a) other-session coordination on CLLC, (b) commitment to multi-session Vienna/SRC SPICE buildout, or (c) PFC controller investigation in `PfcControllerDesign.cpp`.

## Already DAB-quality — do NOT touch
`Llc.cpp`, `Dab.cpp`, `AsymmetricHalfBridge.cpp`, `Sepic.cpp`, `Zeta.cpp`,
`Cuk.cpp`, `Weinberg.cpp`, `FourSwitchBuckBoost.cpp`, `Boost.cpp`, `Buck.cpp`,
`PhaseShiftedHalfBridge.cpp`, `Flyback.cpp`, `PushPull.cpp`, `IsolatedBuck.cpp`,
`IsolatedBuckBoost.cpp`, `SingleSwitchForward.cpp`, `TwoSwitchForward.cpp`,
`ActiveClampForward.cpp` (caveat: ACF still has 2/3 designs at 0.26-0.27 NRMSE pending SPICE netlist fix; see Topology queue §6).

## Test infrastructure (committed prior sessions, still in use)

- `tests/WaveformDumpHelpers.h` — `dump_waveforms_csv(tag, analytical_op,
  spice_op, N=1024)`. Env-gated via `MKF_DUMP_WAVEFORMS=1`. Writes CSV to
  `/tmp/mkf_wf_<tag>_<ms>.csv` with columns `time_s, analytical_w<N>_i_A,
  spice_w<N>_i_A, analytical_w<N>_v_V, spice_w<N>_v_V` per winding, plus
  `analytical_im_A, spice_im_A` if magnetizing current is set on excitation
  (currently NOT set in PSFB — would need `set_magnetizing_current(...)` on
  primary excitation before dump captures it).
- Wired into `tests/TestTopologyPhaseShiftedFullBridge.cpp:668` as
  `dump_waveforms_csv("psfb_singleop", analytical[0], sim[0])`.

## Important: NRMSE is shape-only

`tests/PtpHelpers.h:198 ptp_nrmse()` does:
1. Subtract mean from both signals
2. Normalize both to unit AC-RMS
3. Try 64 phase shifts, return minimum

So **DC offset and amplitude scale are factored out** — only the shape
matters. This is why iSec mean-magnitude gaps don't directly affect NRMSE,
but ripple-ratio differences (range/scale relative shape) do.

## Reference / spec files
- `src/converter_models/CONVERTER_MODELS_GOLDEN_GUIDE.md` §15 — DAB-quality
  spec
- `src/converter_models/FUTURE_TOPOLOGIES.md` — tier table (PSFB now annotated
  with diode-dynamic-R FW model status)

## Useful invocations

```bash
# Build (always check no ninja is running first)
pgrep -fa "ninja -j" || echo free
cd build && ninja -j6 MKF_tests 2>&1 | tail -5

# Dump waveforms for a single PSFB test
cd build && MKF_DUMP_WAVEFORMS=1 timeout 120 ./MKF_tests \
    "Test_Psfb_PtP_AnalyticalVsNgspice" 2>&1 | tail -8

# Run full PSFB suite
cd build && timeout 600 ./MKF_tests "[psfb-topology]" 2>&1 | tail -8

# Run Flyback PtP only
cd build && timeout 600 ./MKF_tests "[flyback-topology][ptp]" 2>&1 \
    | grep -iE "nrmse|fail|PtP"

# Push to OpenMagnetics with correct SSH key
git -c core.sshCommand="/usr/bin/ssh -i ~/.ssh/hephaestus_om" push origin main

# Inspect latest waveform dump CSV
ls -t /tmp/mkf_wf_*.csv | head -1
```

## Next session recommended start

Campaign effectively stalled (see queue §7 PFC notes and §8-11). Options for next agent:

1. **CLLC — DAB-quality REACHED on all 6 PtP designs** (commit `5f075d73`).
   Per-design NRMSE: Telecom-500W 3.88 %, Telecom-500W-BelowFr 5.44 %,
   Telecom-250W 4.77 %, **Infineon-11kW 10.43 %** (was 15.75 %),
   KIT-20kW-Asymmetric 12.99 %, Telecom-500W-REVERSE 11.69 %. All under
   the 0.15 DAB-quality gate. Full CLLC suite green (33 cases / 353 assert).

   **What's already done this campaign**:
   - `0f5bde38` 4-state asymmetric TDA per Sun 2020 TPEL — KIT amplitude
     fix (±301 A → physically plausible). State (iLr1, iLm, vCr1, vCr2_pri);
     shifted system block-antidiagonal, closed-form 2×2 eigendecomp.
   - `254770da` damped Picard 4-state solver (α=0.4, ≤200 iters) with
     multi-start over 6 physically-motivated seeds + sanity + degeneracy
     filters. Replaces 4-D Newton that trapped on low-amplitude basins.
   - `5f075d73` **route ALL CLLC tanks through 4-state TDA** (symmetric
     too). 3-state collapsed path retired. F sub-state explicitly models
     dead-time / commutation interval — that's the missing physics that
     was costing Infineon ~5 pp. Infineon gate tightened 0.16 → 0.15.

   **Residual amplitude gap on KIT** (analytical ±38.7 A vs SPICE ±54 A,
   ~28 % low): Picard map is biased toward low-amplitude basin (lossless
   tank → continuum of fixed-points from zero up to load-driven). NRMSE
   is shape-only so still passes gate, but a future amplitude-aware
   improvement would replace the Picard solve with multi-cycle transient
   simulation alternating bridge polarity (no antisymmetry trick), letting
   small numerical dissipation settle the steady state. Not required for
   DAB-quality; document only.

   **Test infrastructure ready** (unchanged): env-gated CSV dump
   (`MKF_DUMP_WAVEFORMS=1`) writes 9-column waveform dumps to
   `/tmp/mkf_wf_cllc_<name>_<ms>.csv`; T_common shared-horizon NRMSE
   removes the period-mismatch artefact; per-test diagnostics print
   t_span, i_range, mean for both signals.

2. **PFC controller deep-dive** — investigate why `fc_i = 5 kHz` current loop
   responds at ~500 µs instead of ~32 µs. Start in
   `src/converter_models/PfcControllerDesign.cpp` (warm-start IC for `vc_i`,
   FF clamp `vrms_ff`, current-EA topology) and
   `src/converter_models/PfcControllerSubcircuits.h`. Use env-gated CSV dump:
   ```bash
   cd build && MKF_DUMP_WAVEFORMS=1 timeout 300 ./MKF_tests \
       "PFC reference design PtP — L4981 1000 W" 2>&1 | tail -8
   # Inspect: awk -F, 'NR==1; NR%50==0' /tmp/mkf_wf_pfc_L4981-1000W.csv | head -40
   ```
3. **Vienna/SRC SPICE buildout** — multi-session, ground-up. Vienna 3-phase
   PFC is non-trivial (6 active switches, 3-phase symmetric model).
