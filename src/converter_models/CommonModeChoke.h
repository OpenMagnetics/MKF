#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "processors/NgspiceRunner.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief CMC configuration types based on number of phases.
 */
enum class CmcConfiguration {
    SINGLE_PHASE,           ///< 2 windings: Line + Neutral
    THREE_PHASE,            ///< 3 windings: L1 + L2 + L3
    THREE_PHASE_WITH_NEUTRAL ///< 4 windings: L1 + L2 + L3 + N
};

/**
 * @brief Structure holding CMC simulation waveforms for analysis
 */
struct CmcSimulationWaveforms {
    std::vector<double> time;
    double frequency;
    
    // Input signals
    std::vector<double> inputVoltage;              // Common mode noise source
    
    // CMC line currents
    std::vector<std::vector<double>> windingCurrents;  // Current through each winding
    
    // LISN output (measured noise)
    std::vector<double> lisnVoltage;              // Voltage at LISN measurement point
    
    // Metadata
    std::string operatingPointName;
    double commonModeAttenuation;   // Attenuation in dB at test frequency
    double commonModeImpedance;     // Impedance in Ohms at test frequency
    double theoreticalImpedance;    // Theoretical Z = 2*pi*f*L in Ohms
};

/**
 * @brief Common Mode Choke (CMC) converter model for EMI filter applications.
 * 
 * Common Mode Chokes are coupled inductors wound on toroidal cores used to 
 * suppress common mode noise while passing differential signals. The design
 * process focuses on:
 * 
 * - **Impedance Requirements**: CMCs must meet minimum impedance specifications
 *   at specified frequencies. Higher impedance = better common mode rejection.
 * 
 * - **Coupling Coefficient**: For effective common mode rejection, windings
 *   should be tightly coupled (k ≈ 1). This is achieved through bifilar
 *   winding on toroidal cores, minimizing leakage inductance.
 * 
 * - **Self-Resonant Frequency**: The CMC must operate below its SRF to avoid
 *   impedance collapse. Typically f_op < 0.25 * SRF.
 * 
 * - **Saturation**: The magnetizing inductance must not saturate under
 *   differential mode (unbalanced) currents.
 * 
 * Supported configurations:
 * - **Single-phase (2 windings)**: Line + Neutral - typical for AC mains
 * - **Three-phase (3 windings)**: L1 + L2 + L3 - for 3-phase systems without neutral
 * - **Three-phase with neutral (4 windings)**: L1 + L2 + L3 + N - full 3-phase+N filtering
 * 
 * The topology generates Inputs suitable for MagneticAdviser with:
 * - SubApplication::COMMON_MODE_NOISE_FILTERING
 * - INTERFERENCE_SUPPRESSION application
 * - Multiple identical windings with PRIMARY isolation side
 * - Minimum impedance requirements for filter calculation
 */
class CommonModeChoke {
public:
    bool _assertErrors = false;

    CommonModeChoke() = default;
    CommonModeChoke(const json& j);

    /**
     * @brief Process the CMC specification into a complete Inputs structure.
     * 
     * @return Inputs structure with design requirements and operating points
     */
    Inputs process();

    /**
     * @brief Generate design requirements for CMC optimization.
     * 
     * Sets up:
     * - N windings with identical turns (all 1:1 turns ratios)
     * - INTERFERENCE_SUPPRESSION application
     * - COMMON_MODE_NOISE_FILTERING sub-application
     * - Minimum impedance requirements
     * - Toroidal core preference
     * 
     * N = 2 for single-phase, 3 for three-phase, 4 for three-phase+neutral
     */
    DesignRequirements process_design_requirements();

    /**
     * @brief Generate operating points representing CMC excitation.
     * 
     * For CMCs:
     * - Common mode: All windings carry in-phase current (additive MMF)
     * - Differential mode: Windings carry phase-shifted currents (canceling MMF)
     * 
     * For three-phase systems, currents are 120° phase shifted.
     * The operating point models nominal differential mode current flow.
     * 
     * @return Vector of operating points
     */
    std::vector<OperatingPoint> process_operating_points();

    /**
     * @brief Get the number of windings based on configuration.
     */
    size_t get_number_of_windings() const;

    // Accessors for CMC parameters

    void set_configuration(CmcConfiguration value) { _configuration = value; }
    CmcConfiguration get_configuration() const { return _configuration; }

    void set_operating_voltage(DimensionWithTolerance value) { _operatingVoltage = value; }
    DimensionWithTolerance get_operating_voltage() const { return _operatingVoltage; }

    void set_operating_current(double value) { _operatingCurrent = value; }
    double get_operating_current() const { return _operatingCurrent; }

    void set_line_frequency(double value) { _lineFrequency = value; }
    double get_line_frequency() const { return _lineFrequency; }

    void set_neutral_current(std::optional<double> value) { _neutralCurrent = value; }
    std::optional<double> get_neutral_current() const { return _neutralCurrent; }

    void set_minimum_impedance(std::vector<ImpedanceAtFrequency> value) { _minimumImpedance = value; }
    std::vector<ImpedanceAtFrequency> get_minimum_impedance() const { return _minimumImpedance; }

    void set_line_impedance(double value) { _lineImpedance = value; }
    double get_line_impedance() const { return _lineImpedance.value_or(50.0); }

    void set_maximum_dc_resistance(std::optional<double> value) { _maximumDcResistance = value; }
    std::optional<double> get_maximum_dc_resistance() const { return _maximumDcResistance; }

    void set_maximum_leakage_inductance(std::optional<double> value) { _maximumLeakageInductance = value; }
    std::optional<double> get_maximum_leakage_inductance() const { return _maximumLeakageInductance; }

    void set_ambient_temperature(double value) { _ambientTemperature = value; }
    double get_ambient_temperature() const { return _ambientTemperature; }

    /**
     * @brief Generate an ngspice circuit for CMC EMI testing with LISN.
     * 
     * Creates a test circuit with:
     * - Common mode noise source (represents switching noise)
     * - LISN (Line Impedance Stabilization Network) per CISPR 16
     * - CMC between noise source and LISN
     * - Measurement points for attenuation calculation
     * 
     * @param inductance Magnetizing inductance of each winding
     * @param frequency Test frequency for AC analysis
     * @return SPICE netlist string
     */
    std::string generate_ngspice_circuit(double inductance, double frequency = 150000);

    /**
     * @brief Simulate CMC and extract waveforms.
     * 
     * @param inductance Magnetizing inductance
     * @param frequencies Test frequencies
     * @return Simulation waveforms for analysis
     */
    std::vector<CmcSimulationWaveforms> simulate_and_extract_waveforms(
        double inductance,
        const std::vector<double>& frequencies);

    /**
     * @brief Simulate and extract operating points from simulation.
     * 
     * @param inductance Magnetizing inductance
     * @return Vector of OperatingPoints from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);

private:
    CmcConfiguration _configuration = CmcConfiguration::SINGLE_PHASE;  ///< Default to single-phase
    DimensionWithTolerance _operatingVoltage;
    double _operatingCurrent;
    double _lineFrequency;  ///< Mains/line frequency in Hz (50 or 60 Hz typically) - REQUIRED
    std::optional<double> _neutralCurrent;  ///< Current in neutral (for 4-winding config)
    std::vector<ImpedanceAtFrequency> _minimumImpedance;
    std::optional<double> _lineImpedance;
    std::optional<double> _maximumDcResistance;
    std::optional<double> _maximumLeakageInductance;
    double _ambientTemperature = 25.0;
};

} // namespace OpenMagnetics
