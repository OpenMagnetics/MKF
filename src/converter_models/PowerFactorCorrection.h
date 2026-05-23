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

    /// CCM: L = Vin_min·D·(1−D) / (ΔI·f_sw), worst-case at minimum input voltage.
    double calculate_inductance_ccm();
    /// DCM: L = Vin²·D² / (2·P·f_sw).
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
    double get_diode_voltage_drop_required() const {
        auto v = MAS::PowerFactorCorrection::get_diode_voltage_drop();
        if (!v.has_value())
            throw std::runtime_error(
                "PowerFactorCorrection: diodeVoltageDrop is required (set via "
                "set_diode_voltage_drop() or the JSON field 'diodeVoltageDrop')");
        return v.value();
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
     * @brief Mode as the legacy long string used by MKF logic and JSON I/O.
     */
    std::string get_mode_string() const {
        auto m = MAS::PowerFactorCorrection::get_mode();
        if (!m.has_value()) return "Continuous Conduction Mode";
        switch (m.value()) {
            case PfcModes::CONTINUOUS_CONDUCTION_MODE:    return "Continuous Conduction Mode";
            case PfcModes::DISCONTINUOUS_CONDUCTION_MODE: return "Discontinuous Conduction Mode";
            case PfcModes::CRITICAL_CONDUCTION_MODE:      return "Critical Conduction Mode";
            case PfcModes::TRANSITION_MODE:               return "Transition Mode";
        }
        return "Continuous Conduction Mode";
    }

    void set_num_periods_to_extract(int value) { _numberOfPeriods = value; }
    int  get_num_periods_to_extract() const    { return _numberOfPeriods; }

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

private:
    // MKF-only field — number of mains periods to synthesize per operating
    // point. Not part of the MAS schema (it's a simulation knob, not a
    // physical spec).
    int _numberOfPeriods = 2;
};

} // namespace OpenMagnetics
