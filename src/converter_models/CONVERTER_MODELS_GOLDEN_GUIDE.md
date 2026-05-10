# Converter Models — Golden Guide

This file is the **mandatory** specification every converter model in
`src/converter_models/` must follow. The DAB model
(`Dab.cpp` / `Dab.h` / `TestTopologyDab.cpp`) is the reference
implementation. Any new topology and any rework of an existing one
must clear every checkbox below before it is considered "DAB-quality".

If you are an AI agent improving an existing model or adding a new one,
read this guide end-to-end first, then implement against it. Push back on
the user only if a checkbox is genuinely impossible for the topology
(very rare); otherwise just satisfy it.

---

## 1. File layout

Every converter contributes exactly the following files:

```
src/converter_models/<Topology>.h
src/converter_models/<Topology>.cpp
tests/TestTopology<Topology>.cpp
```

Optionally, an `Advanced<Topology>` sub-class exposed through the same
header (for users who want to bypass the design solver and provide
turns-ratio / Lm / Lr / Lo directly).

Shared helpers live in `src/converter_models/PwmBridgeSolver.h`,
`src/converter_models/Topology.h`, `src/converter_models/Topology.cpp`.

---

## 2. Class structure (mirror the DAB pattern)

```cpp
class <Topology> : public MAS::<TopologyMasClass>, public Topology {
private:
    // 2.1 Simulation tuning
    int numPeriodsToExtract = 5;        // post-steady-state periods to sample
    int numSteadyStatePeriods = 3;      // burn-in periods before sampling

    // 2.2 Computed design values, populated by process_design_requirements()
    double computedSeriesInductance = 0;     // resonant or leakage
    double computedMagnetizingInductance = 0;
    double computedOutputInductance = 0;     // if applicable
    double computedDeadTime = 200e-9;
    double computedEffectiveDutyCycle = 0;   // if applicable
    double computedDiodeVoltageDrop = 0.6;   // if rectifier present

    // 2.3 Schema-less device parameters (Coss, etc) with sensible defaults
    double mosfetCoss = 200e-12;        // Drain-source output capacitance

    // 2.4 Per-OP diagnostic outputs (mutable — populated by process)
    mutable double lastDutyCycleLoss = 0.0;
    mutable double lastEffectiveDutyCycle = 0.0;
    mutable double lastZvsMarginLagging = 0.0;
    mutable double lastZvsLoadThreshold = 0.0;
    mutable double lastResonantTransitionTime = 0.0;
    mutable double lastPrimaryPeakCurrent = 0.0;
    // (Add more if the topology needs them — modulation type, mode index,
    //  sub-interval angles, etc. DAB has ~6, LLC has ~5.)

    // 2.5 Extra-component waveforms — one entry per OP, cleared in
    //     process_operating_points().
    mutable std::vector<Waveform> extraLrVoltageWaveforms;
    mutable std::vector<Waveform> extraLrCurrentWaveforms;
    // Add per-component pairs (Lo, Lo2, Cb, etc.) as the topology requires.

public:
    bool _assertErrors = false;

    // 2.6 Constructors
    <Topology>(const json& j);
    <Topology>() {};

    // 2.7 Tuning accessors
    int  get_num_periods_to_extract() const;
    void set_num_periods_to_extract(int v);
    int  get_num_steady_state_periods() const;
    void set_num_steady_state_periods(int v);

    // 2.8 Computed-value accessors (one get/set per private double)
    double get_computed_series_inductance() const;
    void   set_computed_series_inductance(double v);
    /* …repeat for every computed* field… */

    // 2.9 Schema-less device parameters
    double get_mosfet_coss() const;
    void   set_mosfet_coss(double v);

    // 2.10 Per-OP diagnostic accessors (read-only)
    double get_last_duty_cycle_loss() const;
    double get_last_effective_duty_cycle() const;
    double get_last_zvs_margin_lagging() const;
    double get_last_zvs_load_threshold() const;
    double get_last_resonant_transition_time() const;
    double get_last_primary_peak_current() const;

    // 2.11 Topology interface (override these)
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // 2.12 Per-input-voltage worker (called from process_operating_points)
    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const <TopologyMasClass>OperatingPoint& opPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // 2.13 Static analytical helpers (topology-specific equations).
    //      These exist as STATIC functions so tests can call them without
    //      instantiating the class. Examples: compute_power, compute_phase_shift,
    //      compute_zvs_margin, compute_turns_ratio.
    static double compute_<thing>(double a, double b, ...);

    // 2.14 Extra-components factory (returns one Inputs per ancillary
    //      magnetic — Lr, Lo, Cb, etc.).
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // 2.15 SPICE simulation
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
```

---

## 3. Header docstring requirements

The class docstring (just before `class <Topology> :`) must include:

- [ ] **TOPOLOGY OVERVIEW** — ASCII art of the schematic (DAB lines 22–32).
- [ ] **MODULATION CONVENTION** — explicit mapping between MAS schema
      fields and the variables used in equations. Specify radians vs
      degrees. If the topology has multiple modes (DAB SPS/EPS/DPS/TPS,
      LLC mode-1..6) document each one.
- [ ] **KEY EQUATIONS** — at minimum: turns-ratio, conversion ratio,
      power transfer, inductor / capacitor sizing, ZVS condition. Each
      equation tagged with the reference paper / app-note that introduces
      it (Sabate APEC 1990, Andreycak U-136A, etc.).
- [ ] **References list** — full citations (author, title, venue, year).
- [ ] **Accuracy disclaimer** — what the SPICE model omits (Coss
      variation, body-diode reverse recovery, gate-driver delay).
      DAB lines 296–320.
- [ ] **Disambiguation note** if the topology name is ambiguous (PSHB
      must call out the AHB vs Pinheiro-Barbi distinction; LLC vs SRC
      vs DAB; flyback vs flybuck etc).

---

## 4. Analytical model (process_operating_point_for_input_voltage)

### 4.1 The sub-interval rule

Every PWM converter must use a **piecewise-linear sub-interval solver**.
Do NOT use forward Euler integration. Within each sub-interval the
primary voltage is constant, so iL(t) and Im(t) are EXACTLY linear and
should be propagated by a single multiply-add per sub-interval.

For DAB this is `dab_build_subintervals` / `dab_propagate_period`. For
PSFB / PSHB / AHB use `PwmBridgeSolver::build_first_half_cycle`.
Resonant converters (LLC, CLLC) use a Nielsen-style time-domain
propagator. The principle is the same: enumerate the breakpoints, treat
each segment as exact-linear (or exact-sinusoid for resonant), build
initial conditions from half-wave antisymmetry.

### 4.2 Initial conditions

Derive from **half-wave antisymmetry**: `iL(Ts/2) = −iL(0)`, applied to
the integrated voltage over the first half-cycle. DAB
`dab_initial_conditions` is the canonical implementation. Never start
the solver at "t = 0, iL = 0" and rely on warm-up periods.

### 4.3 Magnetizing current

`Im(t)` integrates the primary winding voltage divided by Lm. It is
ALWAYS linear within a sub-interval and is added on top of the
reflected load current to produce the primary winding current.
Warning: `Im` should clamp at ±Im_peak to prevent drift due to
floating-point integration error.

**Critical bug to avoid (PSFB / PSHB):** `Im_peak` is computed from
the *commanded* duty cycle `D_cmd`, NOT the post-DCL effective duty
cycle `Deff`. The primary winding sees +Vbus during BOTH the
commutation interval (length `dcl·Thalf`) AND the active power
transfer interval (length `(D_cmd − dcl)·Thalf`). Total +Vbus time per
half-cycle = `D_cmd · Thalf`. Therefore:

```
Im_peak = Vbus · D_cmd / (4 · Fs · Lm)
```

Using `Deff` instead under-estimates Im_peak by a factor of
`Deff/D_cmd` — typically ~10 % off and a major contributor to
analytical-vs-SPICE NRMSE.

### 4.4 Duty-cycle loss (PWM bridges only)

Compute via `PwmBridgeSolver::compute_duty_cycle_loss(Vbus, Lk, Io, n, Fs)`.
Use the EFFECTIVE duty cycle (D_cmd − ΔD) for all downstream waveform
generation. Iterate the design step (turns ratio) so n converges to a
self-consistent value.

**Per-section duty usage** (PSFB / PSHB):

| Quantity                | Use D_cmd or Deff?                       |
|-------------------------|------------------------------------------|
| `Im_peak`               | **D_cmd** (Vbus appears for full D_cmd)  |
| `Vo` (output voltage)   | **Deff**  (rectifier delivers Deff)      |
| Output inductor ripple  | **Deff** (active interval = Deff·Thalf)  |
| `t_active` boundary     | **D_cmd** (Vp goes to 0 at D_cmd·Thalf)  |
| `t_dcl` boundary        | **dcl_duty = ΔD**                        |
| Sub-interval Vp values  | Vbus during [0, D_cmd·Thalf], 0 after    |
| ZVS Ip_pk (lagging leg) | end of active = ILo_max/n + Im_peak      |

### 4.4.5 Output-inductor (Lo) waveform — PSFB / PSHB

The output inductor sees three slope regions per half-cycle:

| Sub-interval                        | Vo applied to Lo | Slope            |
|-------------------------------------|------------------|------------------|
| Commutation `[0, t_dcl]`            | −Vo (freewheel)  | `-Vo / Lo`       |
| Active `[t_dcl, t_active]`          | Vsec_pk − Vo     | `(Vsec_pk-Vo)/Lo`|
| Freewheel `[t_active, Thalf]`       | −Vo              | `-Vo / Lo`       |

**Steady-state geometry:**
- **Trough of ILo** at `t = t_dcl` (end of commutation, just before
  rectifier picks up).
- **Peak of ILo** at `t = t_active` (end of active interval).
- ΔILo (peak-to-peak) = `Vo · (1 − Deff) / (Fs · Lo)` — uses `Deff`
  because that's the true active fraction.
- ILo at `t = 0` (start of half-cycle) = `ILo_min + Vo · t_dcl / Lo`,
  i.e. partway between trough and peak.

The output inductor's ripple repeats every half-cycle (rectified
shape), so the load sees ripple at `2 · Fs`.

A wrong-direction ramp during commutation (treating it like the active
interval) is one of the largest sources of analytical-vs-SPICE NRMSE.

### 4.5 Rectifier-type-aware secondaries

If the schema has a `rectifierType` enum (CT / CD / FB), the secondary
waveforms must differ between the three modes:

- **CENTER_TAPPED**: emit two separate winding excitations per output
  (`Secondary 0a`, `Secondary 0b`). Each carries unipolar pulses on
  alternate half-cycles. Reverse voltage = 2·Vsec_pk.
- **CURRENT_DOUBLER**: one secondary winding, bipolar 3-level
  voltage; two output inductors L1 and L2 emitted via
  `extraLoCurrentWaveforms` and `extraLo2CurrentWaveforms`.
- **FULL_BRIDGE**: one secondary winding, bipolar 3-level voltage;
  one output inductor.

### 4.6 Per-OP diagnostics (populate every time)

At the end of `process_operating_point_for_input_voltage`, populate
`lastZvsMarginLagging`, `lastZvsLoadThreshold`,
`lastResonantTransitionTime`, `lastDutyCycleLoss`,
`lastEffectiveDutyCycle`, `lastPrimaryPeakCurrent` from the converged
solution. Tests rely on these for ZVS-boundary checks.

### 4.7 Multi-output

If the topology supports multiple outputs, follow the DAB load-share
projection (Dab.cpp lines 855–906). Print an explicit warning
exactly once per process() call when the user has more than one
output, stating that per-output secondary currents are approximate
unless per-secondary leakage inductance is provided.

### 4.8 Sample grid

Use `2*N_samples + 1` points where `N_samples = 256`. Time vector
spans `[0, T]` exactly. Primary uses the same grid as secondaries.
Extra-component waveforms (Lo, Lr, Cb) reuse the same grid.

---

## 5. SPICE netlist (generate_ngspice_circuit)

### 5.0 The winding-port vs converter-port rule (CRITICAL)

**Every signal that ends up in `excitations_per_winding[i].voltage` /
`excitations_per_winding[i].current` MUST be measured at the magnetic
component's *winding ports* — i.e. the two terminals of the inductor
or transformer winding inside the magnetic — and NEVER at a
downstream converter node such as a rectified output capacitor, a
filtered output rail, an averaged DC bus, or anything past the
secondary-side rectifier.**

The output rails (rectified Vout, post-LC Iout, DC link) are also
worth capturing, but they belong to the **separate** converter-level
waveform stream returned by
`simulate_and_extract_topology_waveforms` (consumed by the frontend's
*Converter Waveforms* section), and they MUST be tested independently
of winding signals.

Why this matters:

1. The MKF magnetic-design pipeline (`MagneticFilterAreaProduct`,
   `CoreAdviser`, `WindingLossesAdviser`, …) consumes
   `excitations_per_winding`. It assumes each entry is a true winding
   v(t) / i(t) so that `|V·I|`, volt-second product, RMS current and
   peak flux are physically meaningful. Feeding a rectified-output
   node into this stream collapses `|V·I|` by 100×–1000× (the
   secondary "voltage" looks like a flat DC), defeats the
   area-product filter, and lets `CoreAdviser` advance sub-mm³ cores
   that never could have built the winding. This bug class struck
   PushPull (commit `4ffd3c28`) and IsolatedBuck (commit `22bfee9e`)
   independently.
2. The frontend's *Magnetic* view shows the winding waveforms; the
   *Converter* view shows the rectified rails. Conflating the two
   strands silently breaks both.
3. Volt-second balance must hold on every true winding signal in
   steady state (`|avg(v_winding)| < ε · |v_winding|_peak`). It
   intentionally does **not** hold on rectified outputs (Vout averages
   to Vo, by design). The cross-topology
   `Test_<Topo>_VoltSecondBalance_AllWindings` (§8) is the gate.

Mechanical rules in the netlist:

- **Primary winding voltage** = differential across the primary
  inductor's terminals. If `Lpri A B`, emit
  `Bvpri_diff vpri_diff 0 V=V(A)-V(B)` and `.save v(vpri_diff)`.
  NEVER `.save v(A)` — that is node-to-GND and is wrong whenever B
  is not GND (PushPull center tap, IsolatedBuck output node, any
  bridge midpoint).
- **Primary winding current** = `i(Vpri_sense)` of a 0-V series
  source on the primary winding branch, with sign matched to the
  passive convention against the differential voltage probe (so
  that `avg(V·I) > 0` on average over a switching period — the
  winding absorbs energy from the source).
- **Secondary winding voltage** = differential across the secondary
  inductor's terminals (`Bvsec<N>_diff vsec<N>_diff 0
  V=V(<dot>)-V(<other>)`). NEVER `v(sec<N>_rect)`,
  `v(out_node_o<N>)`, or any node past the rectifier.
- **Secondary winding current** = `i(Vsec<N>_sense)` of a 0-V
  series source in the secondary winding branch (between the
  winding and the rectifier), passive-convention to the
  differential voltage probe.
- **Output-rail voltage / current** (Vout, ILout) — emit
  `.save v(out_node_o<N>) i(Vout_sense_o<N>)` so they are
  available, but they go to
  `simulate_and_extract_topology_waveforms`'s
  `output_voltages` / `output_currents` slots, NOT into
  `excitations_per_winding`.

`waveformMapping` in `simulate_and_extract_operating_points` MUST
list only winding-port signals for each winding. Any output-rail
signal listed there is a bug, not a feature.

Verify per topology with the **passive-sign sanity test**:
`avg(V_winding · I_winding) > 0` and within ±20 % of the expected
per-winding processed power. A negative `avg(V·I)` means the probe
polarity is reversed (flip the sign in the `B` source or swap the
current sense's terminal order).

### 5.1 Checklist



- [ ] **Switch model**: `.model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg`
      (DAB pattern). VH=0.8 reduces switching-event chatter.
- [ ] **Diode model**: `.model DIDEAL D(IS=1e-12 RS=0.05)`. Slightly
      softer than `RS=0` for solver stability.
- [ ] **Snubber RC across every active switch**: 1 kΩ in series with
      1 nF. Mandatory for converging multi-position switching events.
- [ ] **Per-winding sense sources**: `Vpri_sense`, `Vsec1_sense_oN`,
      `Vsec2_sense_oN`, `Vout_sense_oN`. These zero-volt sources expose
      branch currents to the saved-signal list.
- [ ] **Bridge-voltage probe**: a VCVS named `Evab` that exposes the
      differential bridge output `v(vab) = v(mid_A) − v(mid_C)` (or the
      equivalent for the topology's primary winding).
- [ ] **Multi-secondary K-matrix**: when there are N secondaries, emit
      EVERY pair of K statements:
      `K_<i>: L_pri ↔ L_sec_o<i>` for i = 1..N, and
      `K_<i,j>: L_sec_o<i> ↔ L_sec_o<j>` for every i < j.
      Without this, SPICE treats secondaries as independent transformers
      sharing only the primary, producing currents off by 50× or more.
- [ ] **Output-cap pre-charge**: every `vout_node_oN` has an `.ic`
      directive setting it to its target Vo. Skipping this forces the
      simulator to charge the cap from 0 V via the rectifier, adding
      seconds of warm-up.
- [ ] **Solver options** (mandatory, exact strings):
```
.options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500
.options METHOD=GEAR TRTOL=7
```
- [ ] **Rectifier-type-aware emission**: if the topology has a
      `rectifierType` field, branch the netlist to emit the correct
      diode count and output filter. Hard-coding a 4-diode FB rectifier
      regardless of `rectifierType` is a bug.
- [ ] **R_load divide-by-zero guard**: `R_load = max(Vo / Io, 1e-3)`.
      A zero `Io` in the operating point would otherwise crash ngspice.
- [ ] **Time step**: `stepTime = period / 200` (DAB pattern). Use
      `period / 500` only if the topology has very fast resonant ringing.
- [ ] **Saved signals — winding-port stream (consumed by
      `excitations_per_winding`):**
      `v(vpri_diff) i(Vpri_sense)
       v(vsec1_diff) i(Vsec1_sense_o1)
       v(vsec2_diff) i(Vsec2_sense_o2) …`
      One differential voltage probe + one current sense per winding.
      See **§5.0**: NEVER substitute `v(<rect_node>)` or
      `v(<out_node>)` here. Sign must be passive
      (`avg(V·I) > 0`).
- [ ] **Saved signals — converter-port stream (consumed by
      `simulate_and_extract_topology_waveforms`):**
      `v(vab) v(<midpoints>) i(Vdc)
       v(out_node_o1) i(Vout_sense_o1)
       v(out_node_o2) i(Vout_sense_o2) …`
      Match what `simulate_and_extract_topology_waveforms` expects.
      These are for the *Converter Waveforms* frontend section and
      are tested independently of the winding signals.
- [ ] **Comment header**: a `*` comment block at the top with topology
      name, key parameters (Vin, Vo, Fs, Deff/φ/D3, n, L*, rectifier
      type), and the literal "generated by OpenMagnetics" tag.

---

## 6. simulate_and_extract_operating_points (mandatory pattern)

This function fills `excitations_per_winding`, the **winding-port
stream** (see §5.0). The `WaveformNameMapping` it builds MUST
reference only the differential voltage probes (`v(v<winding>_diff)`)
and the per-winding current senses — never an output-rail node.

DO NOT just call `process_operating_points()` and label the result
"SPICE". The function must:

1. Construct an `NgspiceRunner`.
2. If `runner.is_available()` is false, print a warning and fall back to
   analytical generation. Do not throw.
3. For each (input voltage, operating point) pair, build the netlist via
   `generate_ngspice_circuit`, populate `SimulationConfig`, run the
   simulation.
4. On `simResult.success == false`, log the error message, fall back to
   `process_operating_point_for_input_voltage` for that OP, and tag the
   result name with `[analytical]`.
5. On success, build a `WaveformNameMapping` matching the netlist nodes
   and call `NgspiceRunner::extract_operating_point`.

DAB lines 1218–1307 are the canonical pattern; use them verbatim with
node-name substitutions for the topology.

---

## 7. simulate_and_extract_topology_waveforms

This function fills the **converter-port stream** (see §5.0):
input voltage / current at the DC source, output voltage / current
at every rectified rail, and bridge-midpoint waveforms. It feeds
the frontend's *Converter Waveforms* view and is tested
independently of the winding signals (see
`Test_<Topo>_ConverterPortWaveforms` in §8).

Required to round-trip the SPICE waveforms back through MAS for
converter-level (not winding-level) validation. Common bugs:

- **Wrong node names.** The names in `getWaveform("v(...)")` /
  `getWaveform("i(...)")` MUST match what `generate_ngspice_circuit`
  emits in `.save`. The PSFB and PSHB models had `v(pri_in)` and
  `i(vout_sense)` from copy-pasting an old version; the netlist
  emitted `v(vab)` and `i(Vsec1_sense_o1)`. Always grep the netlist for
  the exact names before referencing them in extraction.
- **Multi-output**: emit one entry per secondary in
  `wf.get_mutable_output_voltages()` and `…output_currents()`, indexed
  with the same `o<N>` suffix as the netlist.

---

## 8. Tests (tests/TestTopology<Topology>.cpp)

Mirror the DAB test set. The minimum required `TEST_CASE`s are:

- [ ] **Three industry reference designs** (mandatory for DAB-quality).
      Pick three independently-published designs that span the
      topology's realistic power / voltage envelope — typically:
      1. **Low corner** — small commercial design (TI/ON/Infineon
         eval board, ≤ 50 W or low Vin).
      2. **Mid corner** — typical application (TI/ADI/STMicro app
         note, 100 W – 1 kW industrial / automotive).
      3. **High corner** — high-power or academic reference (IEEE
         paper, vendor reference design, ≥ 1 kW or high Vin).
      Cite each source fully (vendor + part number + page, or
      author + venue + year + page). For **each** of the three
      designs, add **two** test cases:
      - `Test_<Topo>_RefDesign<N>_Values` — instantiate the design
         in MAS, run the design solver + analytical OP, assert the
         published numbers (Vo, Io, Ipri_rms, Lm/Lo/Lr, turns
         ratio, ZVS margin where stated) are reproduced within
         **5 %**.
      - `Test_<Topo>_RefDesign<N>_PtP` — run the same operating
         point through SPICE and assert primary-current
         `ptp_nrmse` ≤ **0.15** vs the analytical waveform.

      Total: **6 reference-design tests** per converter (3 designs ×
      {Values, PtP}). This is the single most effective gate
      against silent regressions and is the deciding criterion for
      promoting a model to "DAB-quality".
- [ ] `Test_<Topo>_Design_Requirements` — verify turns-ratio, Lm, Lr
      / Lo, are positive and within sane ranges; round-trip a power
      target through the design solver.
- [ ] `Test_<Topo>_OperatingPoints_Generation` — multiple input
      voltages, antisymmetry holds, primary current is piecewise
      linear, secondary current reflects through n.
- [ ] `Test_<Topo>_Operating_Modes` — each modulation/operation mode
      (DAB: SPS/EPS/DPS/TPS, LLC: mode-1..6, PSFB: light load /
      full load, AHB: D=0.3/0.5).
- [ ] `Test_<Topo>_ZVS_Boundaries` — sweep Lk or load and verify
      `lastZvsMarginLagging` crosses zero at the predicted threshold.
- [ ] `Test_<Topo>_RectifierType_*` — for converters with a rectifier
      enum, run CT / CD / FB and verify Vo, output-inductor ripple,
      secondary RMS.
- [ ] `Test_<Topo>_DutyCycleLoss` (PWM bridges) — verify that the
      design solver converges and the analytical Vo matches Sabate's
      formula.
- [ ] `Test_<Topo>_SPICE_Netlist` — netlist parses, contains snubbers,
      contains `.options METHOD=GEAR`, contains the right rectifier
      diode count.
- [ ] `Test_<Topo>_VoltSecondBalance_AllWindings` — for **every**
      winding (primary + every secondary), assert
      `|avg(v_winding) · T_period| / (peak(|v_winding|) · T_period)
      < 0.02` (≤ 2 %), in **both** the analytical and the
      simulated path. This is the gate for §5.0: it catches any
      probe that accidentally references a converter-port node
      (rectified output, DC bus, etc.) instead of the winding
      terminals. Skip the SPICE half only if `NgspiceRunner` is
      unavailable; never relax the bound.
- [ ] `Test_<Topo>_ConverterPortWaveforms` — independently of the
      winding-port test, assert the *converter-port* stream from
      `simulate_and_extract_topology_waveforms`:
      `avg(v(out_node_o<N>)) ≈ Vo<N>` within ±5 %,
      `avg(i(Vout_sense_o<N>)) ≈ Io<N>` within ±5 %, output ripple
      within the design spec. Output rails MUST NOT volt-second
      balance to zero — they are DC by design.
- [ ] `Test_<Topo>_SpicePowerSanity` — for the primary winding,
      assert `avg(V·I) > 0` (passive sign) and within ±20 % of
      expected processed power. Catches polarity flips in either
      the differential voltage probe or the current sense.
- [ ] `Test_<Topo>_PtP_AnalyticalVsNgspice` — primary-current shape
      NRMSE (DAB helper `ptp_nrmse`) ≤ **0.15** for all DAB-quality
      converters. The threshold is non-negotiable: it is the marker
      that distinguishes "reference-grade" from "approximate" models.
      Realistic NRMSE numbers we measured:
      - DAB SPS / EPS / DPS / TPS: **0.04–0.08** ✓
      - PSFB (with corrected Im_peak + ILo dynamics): **0.10–0.14** ✓
      - PSHB (with same corrections): **0.10–0.14** ✓
      - Naive PSFB/PSHB (Im_peak using Deff, wrong ILo direction
        during commutation): 0.40–0.50 ✗ — DO NOT ship like this.
- [ ] `Test_<Topo>_Multiple_Outputs_PtP_*` — uniform / moderate-spread
      / large-spread V₂ vectors. Skip individual current accuracy if
      leakage is unspecified, but still assert primary RMS sanity.
- [ ] `Test_<Topo>_Waveform_Plotting` — CSV dump under `tests/output/`
      for visual inspection. Not asserted, just generated.
- [ ] `Test_Advanced<Topo>_Process` — round-trip
      `Advanced<Topo>::process()` returns sane Inputs.

Use the helpers from `TestTopologyDab.cpp` (`ptp_interp`, `ptp_nrmse`,
`ptp_current`, `ptp_current_time`) — copy them or extract to a shared
header.

---

## 9. Naming conventions

- **Topology variable in the user's code**: lowercase abbreviation
  (`dab`, `psfb`, `pshb`, `llc`, `cllc`).
- **Class name**: `Pascal` (Dab, Psfb, Pshb, Llc, Cllc, Flyback).
- **MAS schema class**: `MAS::<TopologyName>` — match the JSON schema
  exactly (DualActiveBridge, PhaseShiftFullBridge, LlcResonant).
- **Operating point class**: `<TopologyName>OperatingPoint` (DabOp,
  PsfbOp, LlcOp). Multi-secondary topologies use a single OP class with
  `output_voltages: [...]` and `output_currents: [...]`.
- **Waveform names in SPICE**:
  - Bridge differential: `vab`.
  - Primary current sense: `Vpri_sense`.
  - Per-output: `out_node_o<N>`, `vsec1_sense_o<N>`,
    `vsec2_sense_o<N>`, `vout_sense_o<N>`.
- **Diagnostic accessors**: prefix `last_` (`get_last_zvs_margin`).
  Computed (design-time) values: prefix `computed_`
  (`get_computed_series_inductance`).

---

## 10. JSON schema integration

Every converter's input parameters live in `MAS/schemas/inputs/topologies/<topology>.json`. Required fields:

```jsonc
{
  "input_voltage":      "DimensionWithTolerance",     // nominal/min/max
  "operating_points":   "array of <Topo>OperatingPoint",
  "efficiency":         "optional double",
  "use_leakage_inductance": "optional bool",          // PWM bridges
  "<inductance>":       "optional double",            // series, output, etc.
  "rectifier_type":     "optional enum",              // PWM bridges with secondaries
  "<modulation>":       "optional enum"               // DAB: SPS/EPS/DPS/TPS
}
```

The `<Topo>OperatingPoint` schema must include:

```jsonc
{
  "ambient_temperature": "double (Celsius)",
  "switching_frequency": "double (Hz)",
  "output_voltages":     "array of double",            // per-secondary
  "output_currents":     "array of double",            // same length
  "phase_shift" or
   "inner_phase_shift1/2/3": "double (degrees)"        // topology-specific
}
```

Schema-less device parameters (Coss, snubber R/C, gate-driver dead time)
live as private member fields with `get_*`/`set_*` accessors. They
default to typical values and are not exposed through JSON.

---

## 11. Pitfalls and anti-patterns

- ❌ **Forward Euler integration of iL.** Use the closed-form
     piecewise-linear propagator instead.
- ❌ **Trapezoidal-ramp primary current with no duty-cycle loss.**
     The model becomes ~10 % wrong on Vo prediction at high Lk.
- ❌ **Hard-coded rectifier in SPICE regardless of rectifier_type
     field.** Always branch on the enum.
- ❌ **`simulate_and_extract_operating_points` that just calls
     `process_operating_points`.** Wire SPICE through `NgspiceRunner`
     for real.
- ❌ **`getWaveform("v(some_node)")` referencing a node that the
     netlist never declares.** Always grep the netlist for the exact
     name before referencing it.
- ❌ **Divide-by-zero in `R_load = Vo/Io`.** Always
     `max(Vo/Io, 1e-3)`.
- ❌ **Single primary-current waveform for multi-output.** Print the
     load-share approximation warning. Document the bound on per-output
     accuracy (DAB says 20–100 % depending on V₂ spread).
- ❌ **No K-coupling between secondaries in multi-output SPICE.** N
     secondaries need N×(N+1)/2 K statements.
- ❌ **Setting Im(0) = 0** instead of deriving it from half-wave
     antisymmetry. Drifts over the warm-up periods.
- ❌ **Forgetting to clear `extra*Waveforms` at the start of
     `process_operating_points`.** The vector accumulates across
     multiple calls.
- ❌ **Topology label conflation** (PSHB vs AHB, LLC HB vs PSHB,
     SRC vs LLC). Add a Disambiguation note in the header docstring.
- ❌ **Saving a converter-port node as the winding voltage.**
     `.save v(sec1_rect)` / `.save v(out_node_o1)` then feeding it
     into `excitations_per_winding[k].voltage` is wrong even when
     `|V·I|` happens to look "close to Pout" — it collapses the
     bipolar swing, breaks volt-second balance on the magnetic
     model, and silently invalidates every downstream filter
     (`AreaProduct`, `WindingLosses`, core-loss integration). See
     **§5.0**. The winding stream and the converter-port stream
     are distinct and must stay distinct.
- ❌ **Saving `v(node)` as a winding voltage when the other
     winding terminal is not GND.** Always emit a
     `Bv<name>_diff <name>_diff 0 V=V(<dot>)-V(<other>)`
     behavioral source and save *that*. Bug class confirmed in
     PushPull (`4ffd3c28`) and IsolatedBuck (`22bfee9e`); the
     symptom is `Test_MagneticAdviserFromConverter_<Topo>`
     returning 0 cores.
- ❌ **Reversed passive sign on the winding probe.**
     `avg(V_winding · I_winding) < 0` means the topology is
     reporting energy *leaving* the winding on average — a probe
     bug, not physics. Flip the sign of the `B` source or swap
     the current sense's terminal order.

---

## 12. Reference list (cite full)

When introducing equations, cite the source. Acceptable forms:

- TI app notes: `TI TIDU248`, `TI SLUA107`.
- IEEE: `Sabate et al., APEC 1990, pp. 275–284`.
- Textbooks: `Erickson-Maksimović 3rd ed., ch. 6`.

A converter must cite at least one published reference per major
equation block (turns ratio, conversion ratio, ZVS, magnetizing
inductance sizing, output inductor sizing). DAB and LLC have ~5
references each; that's the bar.

---

## 13. The Topology base class contract

Every model inherits both `MAS::<schema_class>` and `Topology`. The
`Topology` base provides:

- `Topology::collect_input_voltages(...)` — extract
  nominal/min/max from a `DimensionWithTolerance` into a labelled vector.
- `Topology::create_isolation_sides(numSecondaries, hasDemag)` —
  populate the isolation-sides array for the design requirements.
- `Topology::run_checks_common(...)` — template helper for shared
  sanity checks. Each model's `run_checks` should call it first, then
  add topology-specific checks.
- `ExtraComponentsMode::IDEAL / REAL` — the IDEAL mode returns the
  nominal Lr, REAL mode subtracts the actual leakage inductance of the
  designed transformer. Every converter's `get_extra_components_inputs`
  must support both.

Use these helpers — do not re-implement them.

---

## 14. Outstanding split: PSHB → AsymmetricHalfBridge + ThreeLevelPSHB

The current `Pshb` class implements a **single-leg + split-cap**
half-bridge with bipolar ±Vin/2 primary swing. This is a hybrid that
matches neither of the two physically-grounded HB families:

- **AsymmetricHalfBridge (AHB)** — Imbertson-Mohan IEEE TIA 1993,
  Pressman ch. 16. Single leg, complementary D / 1−D switches, DC-
  blocking cap Cb in series with primary, **2-level asymmetric**
  primary voltage `+(1−D)·Vin / −D·Vin`, conversion ratio
  `Vo = 2·D·(1−D)·Vin/n`. Suited to 100 W–1 kW.

- **Three-Level Pinheiro-Barbi (3LPSHB)** — Pinheiro & Barbi IEEE
  TPE 8(4) 1993. Two stacked half-bridge legs sharing split caps,
  **3-level** primary voltage `+Vin/2 / 0 / −Vin/2`, conversion ratio
  `Vo = (Vin/2) · D_eff / n`. Suited to 1–5 kW high-Vin telecom/EV.

The current `Pshb` model has the right amplitude (Vin/2) but lacks
the correct circuit topology for either family. Deferred follow-up:

**a) Add MAS schemas:**
```
MAS/schemas/inputs/topologies/asymmetricHalfBridge.json
MAS/schemas/inputs/topologies/threeLevelPhaseShiftHalfBridge.json
```
Each with its own OperatingPoint type carrying the relevant duty /
phase-shift parameters.

**b) Add new C++ classes** in `src/converter_models/`:
```
class AsymmetricHalfBridge : public MAS::AsymmetricHalfBridge, public Topology
class ThreeLevelPshb        : public MAS::ThreeLevelPshb, public Topology
```
Both follow this Golden Guide. The AHB primary has 2-level asymmetric
shape with a DC-blocking cap; the 3LPSHB has a true 4-switch + 2-clamp
diode netlist generating the 3-level Vp swing physically.

**c) Either keep the current `Pshb` as-is** (clearly labelled "single-
leg HB approximation, not a real Pinheiro-Barbi or AHB") **or delete
it** once both replacements ship.

**d) Update users:** check `MAS/data/conformance/` for any test
fixtures referencing `phaseShiftedHalfBridge` and migrate them.

This split is required for true DAB-quality coverage of the half-
bridge family. Until it lands, `Pshb` stays in the codebase as an
acceptable approximation but is documented in its header docstring as
NOT a real Pinheiro-Barbi or AHB topology.

## 15. Acceptance criteria for a "DAB-quality" converter

Final pre-merge checklist:

- [ ] Builds clean (`ninja -j2 MKF_tests` returns 0).
- [ ] All `TestTopology<Topo>` cases pass.
- [ ] **≥ 3 published industry reference designs** reproduced —
      each within 5 % on the published numbers (Vo, Ipri_rms, Lm,
      Lo, Lr, turns ratio, ZVS margin) **and** each with a paired
      `*_PtP` test whose primary-current NRMSE ≤ 0.15. The three
      designs must span the topology's power / voltage envelope
      (low / mid / high corner) and come from independent sources
      (vendor app note + textbook + IEEE paper, or three independent
      vendors). This is the deciding criterion for "DAB-quality".
- [ ] At least one ZVS-boundary test where `lastZvsMarginLagging`
      changes sign at the predicted Lk.
- [ ] Schema fields documented in `MAS/schemas/inputs/topologies/`.
- [ ] Header docstring includes ASCII art, equations, references,
      accuracy disclaimer, optional disambiguation note.
- [ ] `generate_ngspice_circuit` netlist parses and contains snubbers,
      GEAR, IC pre-charge, multi-secondary K-matrix.
- [ ] `simulate_and_extract_operating_points` actually runs SPICE.
- [ ] `simulate_and_extract_topology_waveforms` references real netlist
      nodes.
- [ ] Per-OP diagnostic accessors populated and at least one test
      reads them.
- [ ] Multi-output supported (warn-once on approximation; emit K-matrix
      in SPICE; per-secondary sense sources).
- [ ] Rectifier-type-aware secondary waveforms (if applicable).

---

## 16. Cross-stack synchronization (MAS ↔ MKF ↔ WebLibMKF ↔ WebFrontend)

A converter is **not done** when its `.cpp` builds. A converter ships
across four codebases and breaks silently if any of them fall out of
sync. The following rules are non-negotiable and apply to every new
topology.

### 16.1 The four moving parts

```
            ┌─────────────────── MAS ──────────────────┐
            │  schemas/inputs/topologies/<X>.json      │   schema
            │  schemas/inputs/designRequirements.json  │   topology enum
            │  schemas/inputs.json::supportedTopologies│   map key
            │  scripts/check-mas-hpp.sh                │   regen list
            │  MAS.hpp                                 │   generated C++
            └──────────────────────────────────────────┘
                              │  inherits
                              ▼
            ┌─────────────── MKF ──────────────────────┐
            │  src/converter_models/<X>.{h,cpp}        │   class : public MAS::<X>, Topology
            │  src/advisers/CoreAdviser.cpp            │   is_energy_storing_topology, routing
            │  src/advisers/CoilAdviser.cpp            │   sub_application branches
            │  src/advisers/MagneticFilter.cpp         │   topology switches
            │  src/processors/Inputs.cpp               │   per-topology magnetizing-current synth
            │  tests/TestTopology<X>.cpp               │   DAB-quality bar
            └──────────────────────────────────────────┘
                              │  exposed via embind
                              ▼
            ┌────────── WebLibMKF ─────────────────────┐
            │  src/libMKF.cpp                          │   simulate_<x>_ideal_waveforms
            └──────────────────────────────────────────┘
                              │  loaded by Vite
                              ▼
            ┌────────── WebFrontend ───────────────────┐
            │  src/components/Wizards/<X>Wizard.vue    │   getTopology() returns SCHEMA STRING
            │  src/stores/state.js                     │   Wizards routing enum
            │  src/stores/mas.js                       │   topology readers
            │  src/stores/taskQueue.js                 │   topology routing
            │  WebSharedComponents/assets/js/defaults  │   per-topology MAS defaults
            └──────────────────────────────────────────┘
```

### 16.2 The MAS class is the single source of truth

The MKF C++ class **must inherit from the MAS-generated class**:

```cpp
class MyTopology : public MAS::MyTopology, public Topology { ... };
```

This is the rule that flybuck/CMC/DMC/PFC violated for years. The
penalty is duplication: when the C++ class declares its own
`_inputVoltage` field instead of inheriting from
`MAS::MyTopology::input_voltage`, every read/write must be plumbed
manually through the wizard, the libMKF binding, and back. Any drift
between the schema and the duplicated C++ field is a silent
correctness bug.

**Concretely:**
- Do **not** declare private fields that shadow MAS-generated members.
- Use the inherited `get_*`/`set_*` accessors; do not invent your own.
- Schema-less device parameters (Coss, snubber values, dead-time
  defaults) are the **only** legitimate private fields.
- If you need a field the MAS class does not have, **add it to the
  schema first**, regen `MAS.hpp`, then use the inherited accessor.

### 16.3 Naming consistency rule

Three names must agree across the stack — or the schema validation,
JS routing, and serialization all desync:

| Layer | Token | Example |
|---|---|---|
| `designRequirements.json::topology` enum value | `<x>Converter` | `flybackConverter`, `boostConverter` |
| `inputs.json::supportedTopologies` map key | `<x>` (no `Converter` suffix) | `flyback`, `boost` |
| `MAS::Topologies` C++ enum | `<X>_CONVERTER` | `FLYBACK_CONVERTER`, `BOOST_CONVERTER` |
| `MAS::SupportedTopologies` C++ field | `<x>` | `flyback`, `boost` |
| `<X>Wizard.vue::getTopology()` return value | `<x>Converter` (matches enum) | `flybackConverter` |
| `defaults.js` topology default value | `<x>Converter` | `flybackConverter` |

A historical pitfall: PSFB used `phaseShiftFullBridge` as the map key
but `phaseShiftedFullBridgeConverter` as the enum value, with no
matching C++ class. The 2026 cleanup standardized on the
"shifted" spelling everywhere.

A second historical pitfall: every wizard's `getTopology()` returned
human-readable strings (`'Flyback Converter'`, `'CommonModeChoke'`)
instead of the schema enum string. Round-tripped MAS files were
schema-invalid. Wizards must return the **exact** enum string.

### 16.4 Sub-application is an enum, not a string

`utils.json::subApplication` is an enum. The MKF code uses
`SubApplication::COMMON_MODE_NOISE_FILTERING` etc. Frontend writers
must emit the enum string verbatim (`"commonModeNoiseFiltering"`,
`"differentialModeNoiseFiltering"`, `"powerFiltering"`,
`"transforming"`, `"isolation"`). String-comparing in MKF for new
sub-applications is allowed only via the `SubApplication::*` enum,
never against raw string literals.

### 16.5 The `is_energy_storing_topology` switch

`CoreAdviser::is_energy_storing_topology` (around line 102) gates gap
selection. Every new topology added to the `Topologies` enum **must**
get an explicit case. Do not rely on the `default: return false;`
branch — silent miss-classification produces ungapped designs for
buck/boost-class topologies and gapped designs for transformer-class
topologies.

Truth table guidance:
- Storing energy in a magnetic (single-winding inductor, flyback,
  forward-class with output inductor, PFC, DMC): **true**.
- Pure transformer (push-pull, full-bridge, LLC, CLLC, DAB, CMC,
  current transformer): **false**.

### 16.6 Frontend wizard ↔ MAS ↔ libMKF binding

A new topology is wired end-to-end only when **all four** of these
exist and agree:

1. `MAS/schemas/inputs/topologies/<x>.json` — design-input shape.
2. `MAS::<X>` C++ class — generated by quicktype (added to
   `check-mas-hpp.sh` `-S` list).
3. `WebLibMKF/src/libMKF.cpp::simulate_<x>_ideal_waveforms` —
   embind entry point that constructs the MKF C++ class, calls
   `process_design_requirements()` + `process_operating_points()`,
   and writes to `inputs.operatingPoints[]`.
4. `WebFrontend/src/components/Wizards/<X>Wizard.vue` — UI that
   builds the per-topology block in
   `inputs.converterInformation.supportedTopologies.<x>` **and**
   sets `designRequirements.topology = "<x>Converter"`.

A topology missing any of these is "schema only" and will fail
silently in the adviser pipeline.

### 16.7 Add-a-topology checklist

Before opening a PR adding a new converter:

**MAS:**
- [ ] Created `schemas/inputs/topologies/<x>.json`.
- [ ] Added `<x>Converter` to the `topology` enum in
      `schemas/inputs/designRequirements.json`.
- [ ] Added `<x>` map key under `supportedTopologies` in
      `schemas/inputs.json`.
- [ ] Added `-S ./schemas/inputs/topologies/<x>.json` to
      `scripts/check-mas-hpp.sh`.
- [ ] Re-ran quicktype regeneration; committed `MAS.hpp`.
- [ ] `scripts/check-mas-hpp.sh` exits clean.

**MKF:**
- [ ] `class <X> : public MAS::<X>, public Topology` — inherits, does
      not duplicate fields.
- [ ] `set_topology(Topologies::<X>_CONVERTER)` called in
      `process_design_requirements`.
- [ ] `set_application` and `set_sub_application` set explicitly.
- [ ] Added a case to `CoreAdviser::is_energy_storing_topology`.
- [ ] Added a case to any `MagneticFilter` topology switches.
- [ ] If sub-application is CMC/DMC/etc., added matching branches in
      `CoilAdviser` (insulation, repetitions, layout heuristics).

**WebLibMKF:**
- [ ] Added `simulate_<x>_ideal_waveforms` embind entry.
- [ ] Wired the MKF class constructor + `process_*` calls.

**WebFrontend:**
- [ ] Created `<X>Wizard.vue` extending `ConverterWizardBase`.
- [ ] `getTopology()` returns the **schema enum string** (not a
      human-readable label).
- [ ] Wizard writes `inputs.converterInformation.supportedTopologies.<x>`.
- [ ] Updated `stores/state.js::Wizards` if a routing enum exists.
- [ ] Updated `defaults.js` per-topology MAS defaults.
- [ ] Updated `stores/mas.js` and `stores/taskQueue.js` topology
      readers to recognize the new enum value.
- [ ] Verified MAS round-trip: load default, save, reload — no diff.

**Tests:**
- [ ] `MKF/tests/TestTopology<X>.cpp` covers the DAB-quality bar
      (§15).
- [ ] At least one MCP end-to-end test from wizard → core advisor.

### 16.8 Removing a topology

If a topology has no converter model, drop it everywhere:
- Remove enum value from `designRequirements.json::topology`.
- Remove map key from `inputs.json::supportedTopologies`.
- Remove `-S` flag from `check-mas-hpp.sh`.
- Delete schema file.
- Regen `MAS.hpp`.
- Grep MKF/WebLibMKF/WebFrontend for the enum value and the
  PascalCase class name; remove all references.

The 2026 cleanup removed `flybuck` (duplicate of `isolatedBuck`),
`invertingBuckBoostConverter`, `SEPIC`, `cukConverter`,
`zetaConverter`, `weinbergConverter`, `halfBridgeConverter`,
`fullBridgeConverter` (no converter model). Any future topology
without a converter model should be removed, not left as a stub.

### 16.9 Audit before adding

Before adding a new topology, run a cross-reference audit:
- `git grep "<x>Converter"` in MAS, MKF, WebLibMKF, WebFrontend.
- `git grep "Topologies::<X>"` in MKF.
- `git grep "MAS::<X>"` in MKF, WebLibMKF.
- `git grep "<x>Wizard"` in WebFrontend.

If grep returns hits in fewer than four codebases, the topology is
not actually wired end-to-end. Treat partial wirings as a
correctness bug, not a feature.
- [ ] No divide-by-zero in design-time or SPICE code.
