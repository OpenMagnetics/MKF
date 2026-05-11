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
 * Default fallbacks (currentRippleRatio=0.3, efficiency=0.95, lineFrequency=50,
 * diodeVoltageDrop=0.6) are applied at the read site via value_or() so they
 * stay close to the formulas that consume them, instead of being silently
 * stamped on the MAS object during construction.
 */
class PowerFactorCorrection : public Topology, public MAS::PowerFactorCorrection {
public:
    bool _assertErrors = false;

    PowerFactorCorrection() = default;
    PowerFactorCorrection(const json& j);

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
    // Compatibility accessors with explicit defaults for fields that are
    // optional in MAS but used unconditionally by MKF formulas.
    // ------------------------------------------------------------------
    double get_line_frequency_or_default() const {
        return MAS::PowerFactorCorrection::get_line_frequency().value_or(50.0);
    }
    double get_current_ripple_ratio_or_default() const {
        return MAS::PowerFactorCorrection::get_current_ripple_ratio().value_or(0.3);
    }
    double get_efficiency_or_default() const {
        return MAS::PowerFactorCorrection::get_efficiency().value_or(0.95);
    }
    double get_diode_voltage_drop_or_default() const {
        return MAS::PowerFactorCorrection::get_diode_voltage_drop().value_or(0.6);
    }

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

private:
    // MKF-only field — number of mains periods to synthesize per operating
    // point. Not part of the MAS schema (it's a simulation knob, not a
    // physical spec).
    int _numberOfPeriods = 2;
};

} // namespace OpenMagnetics
