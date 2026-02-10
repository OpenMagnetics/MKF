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
 * PFC converters shape the input current to follow the input voltage waveform,
 * achieving near-unity power factor. The inductor operates with a triangular
 * current ripple superimposed on a half-sinusoidal envelope.
 * 
 * Key design considerations:
 * 
 * - **Operating Mode**:
 *   - CCM (Continuous Conduction Mode): Current never reaches zero
 *   - DCM (Discontinuous Conduction Mode): Current reaches zero each cycle  
 *   - CrCM/TCM (Critical/Transition Mode): Borderline between CCM and DCM
 * 
 * - **Inductance Calculation**:
 *   - CCM: L = Vin_min * D * (1-D) / (ΔI * f_sw)
 *   - DCM: L = Vin² * D² / (2 * P * f_sw)
 *   - CrCM: L determines the variable frequency
 * 
 * - **Core Selection**: High frequency operation (50-150 kHz typical) with
 *   significant AC flux swing requires low-loss core materials
 * 
 * - **Worst Case Operating Point**: Maximum current occurs at minimum input
 *   voltage and peak of the AC line (90° phase)
 * 
 * The topology generates Inputs suitable for MagneticAdviser with:
 * - POWER_CONVERSION application
 * - Single winding design
 * - Operating points at multiple AC line phases
 */
class PowerFactorCorrection : public Topology {
public:
    bool _assertErrors = false;

    PowerFactorCorrection() = default;
    PowerFactorCorrection(const json& j);

    /**
     * @brief Run validation checks on PFC parameters.
     */
    bool run_checks(bool assert = false) override;

    /**
     * @brief Generate design requirements for PFC inductor.
     * 
     * Sets up:
     * - Single winding
     * - POWER_CONVERSION application
     * - Calculated inductance based on ripple requirements
     */
    DesignRequirements process_design_requirements() override;

    /**
     * @brief Generate operating points representing PFC excitation.
     * 
     * Creates operating points at various phases of the AC line to capture:
     * - Peak current (90° phase at minimum Vin)
     * - Maximum flux swing
     * - RMS current for loss calculations
     * 
     * @param turnsRatios Not used for single-winding inductor (pass empty vector)
     * @param magnetizingInductance The calculated inductance
     * @return Vector of operating points
     */
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) override;

    /**
     * @brief Calculate the required inductance for CCM operation.
     * 
     * L = Vin_min * D_max * (1 - D_max) / (ΔI * f_sw)
     * where D_max occurs at minimum input voltage
     */
    double calculate_inductance_ccm();

    /**
     * @brief Calculate the required inductance for DCM operation.
     */
    double calculate_inductance_dcm();

    /**
     * @brief Calculate the required inductance for CrCM/TCM operation.
     */
    double calculate_inductance_crcm();

    /**
     * @brief Calculate duty cycle at a given input voltage and output voltage.
     */
    double calculate_duty_cycle(double vinPeak, double vout);

    /**
     * @brief Calculate peak inductor current at a given operating point.
     */
    double calculate_peak_current(double vinPeak, double inductance);

    // Accessors for PFC parameters
    void set_input_voltage(DimensionWithTolerance value) { _inputVoltage = value; }
    DimensionWithTolerance get_input_voltage() const { return _inputVoltage; }

    void set_output_voltage(double value) { _outputVoltage = value; }
    double get_output_voltage() const { return _outputVoltage; }

    void set_output_power(double value) { _outputPower = value; }
    double get_output_power() const { return _outputPower; }

    void set_line_frequency(double value) { _lineFrequency = value; }
    double get_line_frequency() const { return _lineFrequency; }

    void set_switching_frequency(double value) { _switchingFrequency = value; }
    double get_switching_frequency() const { return _switchingFrequency; }

    void set_current_ripple_ratio(double value) { _currentRippleRatio = value; }
    double get_current_ripple_ratio() const { return _currentRippleRatio; }

    void set_efficiency(double value) { _efficiency = value; }
    double get_efficiency() const { return _efficiency; }

    void set_mode(std::string value) { _mode = value; }
    std::string get_mode() const { return _mode; }

    void set_diode_voltage_drop(double value) { _diodeVoltageDrop = value; }
    double get_diode_voltage_drop() const { return _diodeVoltageDrop; }

    void set_ambient_temperature(double value) { _ambientTemperature = value; }
    double get_ambient_temperature() const { return _ambientTemperature; }

    /**
     * @brief Generate SPICE netlist for PFC boost converter simulation.
     * 
     * Creates a complete boost converter circuit including:
     * - Rectified AC input voltage source (half-sinusoid)
     * - Boost inductor model
     * - MOSFET switch with PWM control
     * - Boost diode
     * - Output capacitor and load
     * 
     * @param inductance The boost inductor value
     * @param dcResistance DC resistance of the inductor
     * @param simulationTime Total simulation time
     * @param timeStep Time step for transient analysis
     * @return Complete SPICE netlist string
     */
    std::string generate_ngspice_circuit(double inductance,
                                         double dcResistance = 0.1,
                                         double simulationTime = 0.02,
                                         double timeStep = 1e-8);

    /**
     * @brief Run ngspice simulation and extract waveforms.
     * 
     * Simulates the PFC circuit and returns detailed waveform data
     * for inductor current, input/output voltages, etc.
     * 
     * @param inductance The boost inductor value
     * @param dcResistance DC resistance of the inductor
     * @param numberOfCycles Number of line frequency cycles to simulate
     * @return PfcSimulationWaveforms structure with all waveform data
     */
    PfcSimulationWaveforms simulate_and_extract_waveforms(double inductance,
                                                          double dcResistance = 0.1,
                                                          int numberOfCycles = 1);

    /**
     * @brief Simulate and extract operating points from ngspice results.
     * 
     * Processes simulation waveforms to extract operating points suitable
     * for magnetic design analysis.
     * 
     * @param inductance The boost inductor value
     * @param dcResistance DC resistance of the inductor
     * @return Vector of operating points extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance,
                                                                       double dcResistance = 0.1);

private:
    DimensionWithTolerance _inputVoltage;
    double _outputVoltage;
    double _outputPower;
    double _lineFrequency = 50.0;
    double _switchingFrequency;
    double _currentRippleRatio = 0.3;
    double _efficiency = 0.95;
    std::string _mode = "Continuous Conduction Mode";
    double _diodeVoltageDrop = 0.6;
    double _ambientTemperature = 25.0;
};

} // namespace OpenMagnetics
