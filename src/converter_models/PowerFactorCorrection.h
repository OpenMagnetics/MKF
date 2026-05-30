#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "processors/NgspiceRunner.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief Structure holding PFC simulation waveforms for analysis
 */
struct PfcSimulationWaveforms {
    std::vector<double> time;
    double switchingFrequency;
    double lineFrequency;

    // Input signals
    std::vector<double> inputVoltage;      // Rectified AC input
    std::vector<double> inputCurrent;      // Input current from AC line

    // Inductor signals
    std::vector<double> inductorVoltage;   // Voltage across boost inductor
    std::vector<double> inductorCurrent;   // Current through boost inductor
    std::vector<double> currentEnvelope;   // Sinusoidal current envelope (ideal reference)
    std::vector<double> currentRipple;     // Current ripple amplitude vs time

    // Output signals
    std::vector<double> outputVoltage;     // DC bus voltage
    std::vector<double> outputCurrent;     // Load current

    // Metadata
    std::string operatingPointName;
    double powerFactor;       // Calculated PF
    double efficiency;        // Calculated efficiency
    double currentThd;        // Total harmonic distortion of input current
};

/**
 * @brief Power Factor Correction (PFC) boost inductor converter model.
 *
 * MAS::PowerFactorCorrection owns the spec fields (input voltage, output power,
 * switching frequency, etc.). MKF inherits and adds the engineering layer:
 * inductance sizing for CCM/DCM/CrCM, operating-point waveform synthesis with
 * envelope-plus-ripple shape, and ngspice circuit generation/extraction.
 *
 * Per the repository "no fallbacks / no defaults / no silent shortcuts — throw"
 * rule (CLAUDE.md), every optional MAS field that this engineering layer needs
 * unconditionally is exposed through a `*_required()` accessor that throws
 * `std::runtime_error` if the optional is unset. The schema-level defaults
 * (`currentRippleRatio=0.3`, `efficiency=0.95`, `lineFrequency=50`,
 * `diodeVoltageDrop=0.6`) are NOT applied silently at the read site — they
 * must be either explicitly set by the caller or written into the JSON.
 */
class PowerFactorCorrection : public Topology, public MAS::PowerFactorCorrection {
public:
    bool _assertErrors = false;

    PowerFactorCorrection() = default;
    PowerFactorCorrection(const json& j);

    MAS::Topologies topology_kind() const override { return MAS::Topologies::POWER_FACTOR_CORRECTION; }

    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;

    /**
     * @param turnsRatios Not used for single-winding inductor (pass empty vector)
     * @param magnetizingInductance The calculated inductance (pass <=0 to recompute from spec)
     */
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;

    /// CCM boost-PFC: L = Vin_peak_min · D / (ΔI · f_sw), worst-case at the
    /// minimum-Vrms peak-of-line (where ΔI is maximum). Derived from
    /// V_L = Vin during the switch ON-time, dt_on = D·Tsw, ΔI = V_L·dt_on/L.
    /// (Note: the buck-converter formula has an extra (1−D); boost does not.)
    double calculate_inductance_ccm();
    /// DCM boost-PFC: L = Vin²·D² / (2·P·f_sw).
    double calculate_inductance_dcm();
    /// CrCM/TCM: boundary between CCM and DCM (peak ripple = 2× envelope average).
    double calculate_inductance_crcm();

    double calculate_duty_cycle(double vinPeak, double vout);
    double calculate_peak_current(double vinPeak, double inductance);

    /**
     * @brief Determine actual mode given an inductance value.
     * @return "Continuous Conduction Mode", "Critical Conduction Mode", or
     *         "Discontinuous Conduction Mode"
     */
    std::string determine_actual_mode(double inductance);

    // ------------------------------------------------------------------
    // Required-field accessors — throw if the underlying MAS optional is
    // unset. Naming convention: `*_required()`. Use these instead of the
    // raw `MAS::PowerFactorCorrection::get_*()` accessors anywhere the
    // engineering math assumes the field is populated.
    // ------------------------------------------------------------------
    double get_line_frequency_required() const {
        auto v = MAS::PowerFactorCorrection::get_line_frequency();
        if (!v.has_value())
            throw std::runtime_error(
                "PowerFactorCorrection: lineFrequency is required (set via "
                "set_line_frequency() or the JSON field 'lineFrequency')");
        return v.value();
    }
    double get_current_ripple_ratio_required() const {
        auto v = MAS::PowerFactorCorrection::get_current_ripple_ratio();
        if (!v.has_value())
            throw std::runtime_error(
                "PowerFactorCorrection: currentRippleRatio is required (set via "
                "set_current_ripple_ratio() or the JSON field 'currentRippleRatio')");
        return v.value();
    }
    double get_efficiency_required() const {
        auto v = MAS::PowerFactorCorrection::get_efficiency();
        if (!v.has_value())
            throw std::runtime_error(
                "PowerFactorCorrection: efficiency is required (set via "
                "set_efficiency() or the JSON field 'efficiency')");
        return v.value();
    }
    /**
     * @brief Required diode forward drop. Only meaningful for variants that
     *        have a boost diode in the conduction path: BOOST and
     *        INTERLEAVED_BOOST. TOTEM_POLE replaces the boost diode with a
     *        synchronously-rectified low-side MOSFET, so the conduction drop
     *        is RDS(on)·I (not Vf); callers MUST NOT invoke this for
     *        TOTEM_POLE — use `effective_diode_voltage_drop()` instead, which
     *        applies the variant-specific policy in one place.
     */
    double get_diode_voltage_drop_required() const {
        auto v = MAS::PowerFactorCorrection::get_diode_voltage_drop();
        if (!v.has_value())
            throw std::runtime_error(
                "PowerFactorCorrection: diodeVoltageDrop is required (set via "
                "set_diode_voltage_drop() or the JSON field 'diodeVoltageDrop')");
        return v.value();
    }

    /**
     * @brief Conduction-path forward drop seen by the duty-cycle / off-time
     *        formulas, taking the topology variant into account:
     *          - TOTEM_POLE     → 0 (synchronous rectifier, no diode Vf)
     *          - BOOST          → get_diode_voltage_drop_required()
     *          - INTERLEAVED_BOOST → get_diode_voltage_drop_required()
     *        Centralising the policy here avoids two-site drift (we used to
     *        branch on variant in calculate_duty_cycle and again in
     *        process_operating_points, which made the TOTEM_POLE-suppresses-
     *        the-required-field invariant invisible).
     */
    double effective_diode_voltage_drop() const {
        return (get_topology_variant_or_default() == PfcTopologyVariants::TOTEM_POLE)
                   ? 0.0
                   : get_diode_voltage_drop_required();
    }

    /**
     * @brief True for the buck-boost-class PFC variants whose INPUT inductor
     *        L1 sits in series with the line (continuous input current) but
     *        whose conversion ratio is buck-boost, not boost: SEPIC and Ćuk.
     *
     * For these the duty relation is Vout/Vin = D/(1−D) ⇒
     *   D = (Vout+Vd)/(Vin+Vout+Vd)   (vs boost's D = 1 − Vin/(Vout+Vd))
     * and L1 sees +Vin during the switch ON-time (so the inductance-sizing
     * formula L = Vin·D/(ΔI·fsw) is shared with boost) but −(Vout+Vd) during
     * OFF (vs boost's Vin−Vout), which the operating-point voltage waveform
     * accounts for. BUCK and BUCK_BOOST are NOT in this class: their input
     * current is discontinuous (inductor not in series with the line), so
     * they remain unsupported.
     */
    bool is_buck_boost_class() const {
        auto v = get_topology_variant_or_default();
        return v == PfcTopologyVariants::SEPIC || v == PfcTopologyVariants::CUK;
    }

    // ------------------------------------------------------------------
    // Topology-variant helpers (PfcTopologyVariants enum, MAS schema).
    //
    // Currently implemented variants:
    //   - BOOST              (default if topologyVariant unset)
    //   - TOTEM_POLE         bridgeless, sync-rectified low side; the boost
    //                        inductor sees AC voltage and bipolar current
    //                        (true sine, not |sin|). CCM totemPole requires
    //                        wideBandgapSwitch=true (GaN/SiC) per
    //                        Erickson §17 and ON Semi AND8016.
    //   - INTERLEAVED_BOOST  N parallel boost cells phase-shifted; the spec
    //                        models per-phase magnetics (one inductor at
    //                        Pout/N), so process_design_requirements returns
    //                        the per-phase inductance and the user replicates
    //                        the resulting magnetic across N phases.
    //
    // All other PfcTopologyVariants values throw with a "not yet implemented"
    // message at the first engineering-layer call (run_checks /
    // calculate_inductance_*). VIENNA in particular is a 3-phase 3-level
    // topology with its own dedicated VIENNA_PLAN.md; do NOT enable it as
    // a PFC variant code branch.
    // ------------------------------------------------------------------
    PfcTopologyVariants get_topology_variant_or_default() const {
        auto v = MAS::PowerFactorCorrection::get_topology_variant();
        return v.has_value() ? v.value() : PfcTopologyVariants::BOOST;
    }
    int64_t get_number_of_phases_or_default() const {
        auto v = MAS::PowerFactorCorrection::get_number_of_phases();
        return v.has_value() ? v.value() : 1;
    }
    bool get_wide_bandgap_switch_or_default() const {
        auto v = MAS::PowerFactorCorrection::get_wide_bandgap_switch();
        return v.has_value() ? v.value() : false;
    }
    /**
     * @brief Validate the configured topologyVariant + (numberOfPhases,
     *        wideBandgapSwitch, mode) combination. Throws std::runtime_error
     *        on any unsupported variant or missing companion field.
     *        Idempotent and cheap; called from every public engineering
     *        entry point (run_checks, calculate_inductance_*,
     *        process_operating_points, generate_ngspice_circuit, …).
     */
    void validate_topology_variant() const;

    /**
     * @brief Vrms to use for analytical-simulation / netlist generation.
     *
     * Policy (centralised so the three sim paths agree):
     *   - If `nominal` is set → use it (the user-specified operating point).
     *   - Else                → use the MINIMUM of the range (worst-case for
     *                           inductance ripple / current stress).
     * Throws if neither is populated.
     */
    double get_vrms_for_simulation() const;

    /**
     * @brief Mode as the legacy long string used by MKF logic and JSON I/O.
     *        Throws if the `mode` field is unset (per CLAUDE.md: no silent
     *        default to CCM — the caller must pick a mode explicitly).
     */
    std::string get_mode_string() const {
        auto m = MAS::PowerFactorCorrection::get_mode();
        if (!m.has_value()) {
            throw std::runtime_error(
                "PowerFactorCorrection: mode is required (set via set_mode() "
                "or the JSON field 'mode'). Valid values: "
                "continuousConductionMode, discontinuousConductionMode, "
                "criticalConductionMode, transitionMode.");
        }
        switch (m.value()) {
            case PfcModes::CONTINUOUS_CONDUCTION_MODE:    return "Continuous Conduction Mode";
            case PfcModes::DISCONTINUOUS_CONDUCTION_MODE: return "Discontinuous Conduction Mode";
            case PfcModes::CRITICAL_CONDUCTION_MODE:      return "Critical Conduction Mode";
            case PfcModes::TRANSITION_MODE:               return "Transition Mode";
        }
        throw std::runtime_error(
            "PowerFactorCorrection::get_mode_string: unknown PfcModes enum value");
    }

    void set_num_periods_to_extract(int value) { _numberOfPeriods = value; }
    int  get_num_periods_to_extract() const    { return _numberOfPeriods; }

    /**
     * @brief Suggest a bulk (hold-up) capacitance value from the standard
     *        energy-balance formula:
     *
     *            C_bus = 2 · P_out · Δt_holdup / (V_bus² − V_bus_min²)
     *
     *        Derived from the energy stored in C_bus dropping from V_bus to
     *        V_bus_min during the hold-up interval Δt_holdup with the load
     *        drawing P_out (Erickson §17.2; TI SLUA754). Typical sizing for
     *        a 400 V bus / 20 ms hold-up: 1–4 µF/W.
     *
     * @param vBus        Nominal bus voltage [V] (> 0).
     * @param vBusMin     Minimum acceptable bus voltage during hold-up [V]
     *                    (> 0, < vBus — typically 0.7·vBus).
     * @param pOut        Output power [W] (> 0).
     * @param holdupTime  Hold-up time [s] (> 0; one mains cycle = 0.02 s
     *                    for 50 Hz, 0.0167 s for 60 Hz, common spec values
     *                    are 10, 20 or 30 ms).
     * @return            Required C_bus [F]. Throws on invalid arguments.
     */
    static double suggested_bulk_capacitance(double vBus, double vBusMin,
                                              double pOut, double holdupTime);

    /**
     * @brief Generate SPICE netlist for PFC boost converter simulation.
     */
    std::string generate_ngspice_circuit(double inductance,
                                         double dcResistance = 0.1,
                                         double simulationTime = 0.02,
                                         double timeStep = 1e-8);

    /**
     * @brief Run analytical PFC waveform synthesis and return raw arrays.
     */
    PfcSimulationWaveforms simulate_and_extract_waveforms(double inductance,
                                                          double dcResistance = 0.1,
                                                          int numberOfCycles = 1);

    /**
     * @brief Wrap simulated waveforms in OperatingPoint(s).
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance,
                                                                       double dcResistance = 0.1);

    /**
     * @brief §5.1 converter-port stream — line-frequency view of the PFC's
     *        external ports (NOT the magnetic-component winding port).
     *
     * Sweep of `collect_input_voltages(get_input_voltage(), ...)` × the analytical
     * waveform synthesis already implemented by `simulate_and_extract_waveforms`.
     * Returned per ConverterWaveforms:
     *   - `input_voltage`  = rectified line voltage  Vpk·|sin(ω_line·t)|
     *                        line-cycle MEAN  ≈ 2·√2·Vrms / π  ≈ 0.9·Vrms
     *                        line-cycle RMS   ≈ Vrms
     *   - `input_current`  = inductor / line current envelope (with switching ripple)
     *                        line-cycle MEAN ≈ 2·√2·Iin_rms / π
     *   - `output_voltages = [Vbus]`  flat DC bus  (we model the regulated bus
     *                        as a constant — the PFC analytical model does not
     *                        currently include the bulk-cap state variable, so
     *                        2·f_line ripple on Vbus is NOT represented here).
     *   - `output_currents = [Iload]`  flat DC load current = Pout / Vbus
     *
     * The PFC's "input" is fundamentally AC, so the §5.1 DC-stream gate
     * (`ConverterPortChecks::check_dc_ports`) does NOT apply to the input port.
     * Use `ConverterPortChecks::check_pfc_ports` instead, which validates
     *   - input_voltage line-cycle MEAN ≈ 0.9·Vrms and RMS ≈ Vrms
     *   - output_voltage / output_current MEAN against nameplate.
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        double inductance,
        size_t numberOfCycles = 1);

    /**
     * @brief Generate a real switching boost-PFC ngspice netlist driven by
     *        a native average-current-mode controller (built from
     *        OPAMP_IDEAL / B-source primitives — no vendor IP). Unlike
     *        `generate_ngspice_circuit`, which emits a behavioural-source
     *        synthesis of the analytical model, this method produces a
     *        full switching circuit with a real boost stage (B_vin
     *        rectified-sine source → L1 → ideal switch S1 → ideal diode D1
     *        → Cbus → Rload) and the six-block analog controller documented
     *        in `PfcControllerDesign.h` (voltage error amp, RMS² feed-
     *        forward, multiplier, current error amp, sawtooth, PWM
     *        comparator). Initial conditions warm-start every state
     *        variable to its analytical steady-state value so the
     *        simulation begins in the operating regime instead of having
     *        to settle from cold start.
     *
     * Throws for non-BOOST variants (TOTEM_POLE / INTERLEAVED_BOOST require
     * topology-specific switching netlists and are out of scope here), and
     * for any missing required spec field (no fallbacks).
     *
     * @param inductance       Boost inductance [H]
     * @param numberOfLineCycles Number of line cycles to simulate (typ. 3 —
     *                           Cbus envelope settles within ~1 cycle from
     *                           warm-start ICs; the last cycle is the one
     *                           the test consumes).
     * @return                 Complete ngspice netlist string (transient
     *                         analysis only; no .control block).
     */
    std::string generate_ngspice_switching_circuit(double inductance,
                                                    int numberOfLineCycles = 3);

    /**
     * @brief Build the DCM SEPIC/Ćuk switching netlist (voltage-mode
     *        variable-on-time control). Routed to by
     *        generate_ngspice_switching_circuit when the variant is SEPIC/CUK
     *        and the mode is DISCONTINUOUS_CONDUCTION_MODE.
     *
     * CCM SEPIC is intractable for a switching PtP (the Cc-L2 resonance limits
     * average-current control to ~22% envelope error — see PFC_TOPOLOGY_HANDOFF
     * §2b). DCM is resonance-free and self-shaping: a slow type-II voltage loop
     * sets the on-time magnitude and DCM does the current shaping. A
     * variable-on-time feed-forward Ton ∝ 1/sqrt(1+Vin/Vo) cancels the residual
     * DCM-COT distortion, giving unity power factor (Shen 2018 Electronics
     * Letters; Lin et al. 2023 Electronics 12(8):1807). No current loop / no
     * multiplier.
     *
     * @param inductance         Input inductor L1 [H] (= 2·Le for the coupled
     *                           L1=L2 SEPIC; from calculate_inductance_dcm).
     * @param numberOfLineCycles Number of line cycles to simulate.
     * @return                   Complete ngspice netlist string.
     */
    std::string generate_ngspice_switching_circuit_dcm_sepic(
        double inductance, int numberOfLineCycles = 3);

    /**
     * @brief Run the switching boost-PFC netlist via NgspiceRunner and
     *        return an OperatingPoint whose primary excitation contains the
     *        boost-inductor current waveform (along with vbus and
     *        vin_rect). Convenience wrapper over
     *        `generate_ngspice_switching_circuit` + NgspiceRunner. Throws
     *        if ngspice is unavailable or the simulation fails.
     *
     * @param inductance              Boost-inductor value [H].
     * @param numberOfLineCycles      How many full line cycles to simulate
     *                                (the warm-start transient consumes
     *                                the first ~1-2 cycles).
     * @param trimToLastLineCycle     If true (default), the returned
     *                                OperatingPoint contains only the
     *                                *final* full line period of every
     *                                waveform — the steady-state portion
     *                                that downstream consumers (Painter,
     *                                THD/PF analysers, PtP comparators)
     *                                actually want.  Pass false to keep
     *                                the full simulation history including
     *                                the soft-start transient.
     */
    OperatingPoint simulate_with_ngspice_switching(double inductance,
                                                    int numberOfLineCycles = 3,
                                                    bool trimToLastLineCycle = true);

    // ------------------------------------------------------------------
    // Computed diagnostics — populated by process_design_requirements()
    // and process_operating_points(). These mirror the values shown in
    // the wizard UI under the "Diagnostics" panel. All are populated on
    // every successful process_* call (zeros mean the C++ path skipped
    // its assignment — a bug, not a normal state).
    // ------------------------------------------------------------------
    double get_computed_inductance()          const { return computedInductance; }
    std::string get_computed_actual_mode()    const { return computedActualMode; }
    double get_last_duty_cycle_peak()         const { return lastDutyCyclePeak; }
    double get_last_peak_inductor_current()   const { return lastPeakInductorCurrent; }
    double get_last_inductor_ripple()         const { return lastInductorRipple; }
    double get_last_line_rms_current()        const { return lastLineRmsCurrent; }
    double get_last_input_power()             const { return lastInputPower; }

private:
    // MKF-only field — number of mains periods to synthesize per operating
    // point. Not part of the MAS schema (it's a simulation knob, not a
    // physical spec).
    int _numberOfPeriods = 2;

    // Computed diagnostics members (mutable so const-context process_* can
    // assign them when the topology models declare process methods as const
    // — Buck/Boost convention). PFC's process_* are non-const so plain
    // doubles would also work; using `mutable` keeps the pattern uniform
    // across topologies.
    mutable double computedInductance         = 0.0;   // L from calculate_inductance_<mode> [H]
    mutable std::string computedActualMode    = "";    // "Continuous Conduction Mode" / "Critical..." / "Discontinuous..."
    mutable double lastDutyCyclePeak          = 0.0;   // D at minimum-Vin peak-of-line — worst case
    mutable double lastPeakInductorCurrent    = 0.0;   // I_L peak at line peak, minimum Vin [A]
    mutable double lastInductorRipple         = 0.0;   // ΔI_L at the worst-case sample [A]
    mutable double lastLineRmsCurrent         = 0.0;   // I_in_rms at minimum Vin [A]
    mutable double lastInputPower             = 0.0;   // P_in_avg (per phase for interleaved) [W]

    // PFC is single-shot: process_operating_points populates the last_* fields
    // at the worst-case operating point (minimum Vin, line peak) only — there
    // is no V_in × OP loop because the line-cycle envelope is the operating
    // point. To keep the diagnostics surface uniform with the other topologies,
    // expose a single-row perOp vector named "Worst-case".
public:
    mutable std::vector<std::string> perOpName;
    mutable std::vector<double>      perOpDutyCyclePeak;
    mutable std::vector<double>      perOpPeakInductorCurrent;
    mutable std::vector<double>      perOpInductorRipple;
    mutable std::vector<double>      perOpLineRmsCurrent;
    mutable std::vector<double>      perOpInputPower;

    const std::vector<std::string>& get_per_op_name()                  const { return perOpName; }
    const std::vector<double>&      get_per_op_duty_cycle_peak()       const { return perOpDutyCyclePeak; }
    const std::vector<double>&      get_per_op_peak_inductor_current() const { return perOpPeakInductorCurrent; }
    const std::vector<double>&      get_per_op_inductor_ripple()       const { return perOpInductorRipple; }
    const std::vector<double>&      get_per_op_line_rms_current()      const { return perOpLineRmsCurrent; }
    const std::vector<double>&      get_per_op_input_power()           const { return perOpInputPower; }
};

} // namespace OpenMagnetics
