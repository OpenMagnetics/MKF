#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "processors/NgspiceRunner.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief DMC configuration types based on number of phases.
 * 
 * Unlike CMC which uses coupled inductors, DMC configurations use:
 * - SINGLE_PHASE: One inductor between line and load
 * - THREE_PHASE: Three separate inductors (L1, L2, L3), no neutral filtering
 * - THREE_PHASE_WITH_NEUTRAL: Four inductors including neutral line filtering
 */
enum class DmcConfiguration {
    SINGLE_PHASE,           ///< 1 winding
    THREE_PHASE,            ///< 3 windings: L1 + L2 + L3
    THREE_PHASE_WITH_NEUTRAL ///< 4 windings: L1 + L2 + L3 + N
};

/**
 * @brief Result of attenuation verification via ngspice simulation.
 */
struct DmcAttenuationResult {
    double frequency;           ///< Test frequency in Hz
    double requiredAttenuation; ///< Required attenuation in dB
    double measuredAttenuation; ///< Measured attenuation in dB (from ngspice)
    double theoreticalAttenuation; ///< Theoretical attenuation from LC filter model
    bool passed;                ///< Whether the requirement was met
    std::string message;        ///< Human-readable result message
};

/**
 * @brief Structure holding DMC simulation waveforms for analysis
 */
struct DmcSimulationWaveforms {
    std::vector<double> time;
    double frequency;
    
    // Input/Output signals
    std::vector<double> inputVoltage;    // Noise source voltage
    std::vector<double> outputVoltage;   // Filtered output voltage
    std::vector<double> inductorCurrent; // Current through DMC
    
    // Metadata
    std::string operatingPointName;
    double dmAttenuation;  // Attenuation in dB at test frequency
};

/**
 * @brief Differential Mode Choke (DMC) converter model for EMI filter applications.
 * 
 * Differential Mode Chokes are single inductors used to attenuate differential mode
 * noise (noise between line and neutral or between lines). Unlike CMCs, DMCs present
 * impedance to differential signals and are typically used in LC filter configurations.
 * 
 * Key design considerations:
 * 
 * - **Inductance**: The primary specification. Determines the filter cutoff frequency
 *   when combined with capacitors: f_c = 1/(2π√LC)
 * 
 * - **Saturation Current**: The inductor must not saturate under peak current
 *   including ripple. Core selection depends on energy storage: E = ½LI²
 * 
 * - **DC Resistance**: Lower DCR reduces I²R losses and improves filter Q-factor
 * 
 * - **Self-Resonant Frequency**: Must be above the frequencies being filtered
 * 
 * The topology generates Inputs suitable for MagneticAdviser with:
 * - SubApplication::DIFFERENTIAL_MODE_NOISE_FILTERING  
 * - INTERFERENCE_SUPPRESSION application
 * - Single winding design
 * - Inductance requirements for filter calculation
 */
class DifferentialModeChoke {
public:
    bool _assertErrors = false;

    DifferentialModeChoke() = default;
    DifferentialModeChoke(const json& j);

    /**
     * @brief Process the DMC specification into a complete Inputs structure.
     * 
     * @return Inputs structure with design requirements and operating points
     */
    Inputs process();

    /**
     * @brief Generate design requirements for DMC optimization.
     * 
     * Sets up:
     * - Single winding
     * - INTERFERENCE_SUPPRESSION application
     * - DIFFERENTIAL_MODE_NOISE_FILTERING sub-application
     * - Magnetizing inductance requirements
     */
    DesignRequirements process_design_requirements();

    /**
     * @brief Generate operating points representing DMC excitation.
     * 
     * Models the current through the choke with ripple superimposed
     * on the DC operating current.
     * 
     * @return Vector of operating points
     */
    std::vector<OperatingPoint> process_operating_points();

    // Accessors for DMC parameters
    void set_input_voltage(DimensionWithTolerance value) { _inputVoltage = value; }
    DimensionWithTolerance get_input_voltage() const { return _inputVoltage; }

    void set_operating_current(double value) { _operatingCurrent = value; }
    double get_operating_current() const { return _operatingCurrent; }

    void set_peak_current(std::optional<double> value) { _peakCurrent = value; }
    std::optional<double> get_peak_current() const { return _peakCurrent; }

    void set_minimum_inductance(std::optional<double> value) { _minimumInductance = value; }
    std::optional<double> get_minimum_inductance() const { return _minimumInductance; }

    void set_minimum_impedance(std::optional<std::vector<ImpedanceAtFrequency>> value) { _minimumImpedance = value; }
    std::optional<std::vector<ImpedanceAtFrequency>> get_minimum_impedance() const { return _minimumImpedance; }

    void set_switching_frequency(std::optional<double> value) { _switchingFrequency = value; }
    std::optional<double> get_switching_frequency() const { return _switchingFrequency; }

    void set_maximum_dc_resistance(std::optional<double> value) { _maximumDcResistance = value; }
    std::optional<double> get_maximum_dc_resistance() const { return _maximumDcResistance; }

    void set_ambient_temperature(double value) { _ambientTemperature = value; }
    double get_ambient_temperature() const { return _ambientTemperature; }

    void set_configuration(DmcConfiguration value) { _configuration = value; }
    DmcConfiguration get_configuration() const { return _configuration; }

    void set_filter_capacitance(std::optional<double> value) { _filterCapacitance = value; }
    std::optional<double> get_filter_capacitance() const { return _filterCapacitance; }

    void set_line_frequency(double value) { _lineFrequency = value; }
    double get_line_frequency() const { return _lineFrequency; }

    /**
     * @brief Get the number of windings based on configuration.
     */
    int get_number_of_windings() const {
        switch (_configuration) {
            case DmcConfiguration::SINGLE_PHASE: return 1;
            case DmcConfiguration::THREE_PHASE: return 3;
            case DmcConfiguration::THREE_PHASE_WITH_NEUTRAL: return 4;
            default: return 1;
        }
    }

    /**
     * @brief Generate an ngspice circuit for DMC LC filter testing.
     * 
     * Creates a test circuit with:
     * - DM noise source (differential between line and neutral)
     * - LC low-pass filter (DMC + capacitor)
     * - Load resistor
     * 
     * @param inductance Inductance in Henries
     * @param frequency Test frequency for AC analysis
     * @return SPICE netlist string
     */
    std::string generate_ngspice_circuit(double inductance, double frequency = 150000);

    /**
     * @brief Simulate DMC and extract waveforms.
     * 
     * @param inductance Inductance value
     * @param frequencies Test frequencies
     * @return Simulation waveforms for analysis
     */
    std::vector<DmcSimulationWaveforms> simulate_and_extract_waveforms(
        double inductance,
        const std::vector<double>& frequencies);

    /**
     * @brief Simulate and extract operating points from simulation.
     * 
     * @param inductance Inductance value
     * @return Vector of OperatingPoints from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);

    /**
     * @brief Verify that a design meets the attenuation requirements using ngspice.
     * 
     * This method can be called from the web frontend to verify a proposed design.
     * It simulates the LC filter circuit at the specified frequencies and compares
     * the measured attenuation against requirements.
     * 
     * @param inductance The inductance value of the proposed DMC design in Henries
     * @param capacitance The filter capacitance in Farads (optional, auto-calculated if not provided)
     * @return Vector of verification results for each frequency requirement
     */
    std::vector<DmcAttenuationResult> verify_attenuation(
        double inductance,
        std::optional<double> capacitance = std::nullopt);

    /**
     * @brief Calculate the required inductance to meet attenuation requirements.
     * 
     * Given a target attenuation at a specific frequency and a filter capacitance,
     * calculates the minimum inductance needed.
     * 
     * @param targetAttenuation Target attenuation in dB
     * @param frequency Frequency in Hz
     * @param capacitance Filter capacitance in Farads
     * @return Required inductance in Henries
     */
    static double calculate_required_inductance(
        double targetAttenuation,
        double frequency,
        double capacitance);

    /**
     * @brief Propose a design based on the specifications.
     * 
     * Calculates inductance and capacitance values that will meet the
     * attenuation/impedance requirements.
     * 
     * @return JSON object with proposed design parameters
     */
    json propose_design();

private:
    DimensionWithTolerance _inputVoltage;
    double _operatingCurrent = 1.0;
    std::optional<double> _peakCurrent;
    std::optional<double> _minimumInductance;
    std::optional<std::vector<ImpedanceAtFrequency>> _minimumImpedance;
    std::optional<double> _switchingFrequency;
    std::optional<double> _maximumDcResistance;
    double _ambientTemperature = 25.0;
    DmcConfiguration _configuration = DmcConfiguration::SINGLE_PHASE;
    std::optional<double> _filterCapacitance;
    double _lineFrequency = 50.0;  ///< Mains frequency (50 or 60 Hz)
};

} // namespace OpenMagnetics
