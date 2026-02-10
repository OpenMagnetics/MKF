#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

using namespace MAS;

namespace OpenMagnetics {


/**
 * @brief Structure holding topology-level waveforms for Isolated Buck-Boost converter validation
 * 
 * These waveforms are used to validate that the simulation matches expected
 * converter behavior, not for magnetic component analysis.
 */
struct IsolatedBuckBoostTopologyWaveforms {
    // Time base
    std::vector<double> time;
    double frequency;
    
    // Input side signals
    std::vector<double> inputVoltage; // v(vin_dc) - DC input voltage
    std::vector<double> primaryVoltage; // v(pri_in) - primary winding voltage
    
    // Output side signals (one per secondary winding)
    std::vector<std::vector<double>> secondaryWindingVoltages; // v(sec_N_in) - secondary winding voltages
    std::vector<std::vector<double>> outputVoltages; // v(vout_N) - DC output voltages
    
    // Currents
    std::vector<double> primaryCurrent; // i(vpri_sense) - primary winding current
    std::vector<std::vector<double>> secondaryCurrents; // i(vsec_sense_N) - secondary winding currents
    
    // Metadata
    std::string operatingPointName;
    double inputVoltageValue;
    std::vector<double> outputVoltageValues; // One per secondary
    double dutyCycle;
};


class IsolatedBuckBoost : public MAS::IsolatedBuckBoost, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;

public:
    bool _assertErrors = false;

    IsolatedBuckBoost(const json& j);
    IsolatedBuckBoost() {
    };

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, IsolatedBuckBoostOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency);

    /**
     * @brief Generate an ngspice circuit for this Isolated Buck-Boost converter
     * 
     * @param turnsRatios Turns ratios for secondary windings
     * @param magnetizingInductance Magnetizing inductance in H
     * @param inputVoltageIndex Which input voltage to use (0=nom, 1=min, 2=max)
     * @param operatingPointIndex Which operating point to simulate
     * @return SPICE netlist string
     */
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);
    
    /**
     * @brief Simulate the Isolated Buck-Boost converter and extract operating points
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);
    
    /**
     * @brief Simulate and extract topology-level waveforms for converter validation
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<IsolatedBuckBoostTopologyWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

};

class AdvancedIsolatedBuckBoost : public IsolatedBuckBoost {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;

protected:
public:


    AdvancedIsolatedBuckBoost() = default;
    ~AdvancedIsolatedBuckBoost() = default;

    AdvancedIsolatedBuckBoost(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, AdvancedIsolatedBuckBoost & x);
void to_json(json & j, const AdvancedIsolatedBuckBoost & x);


inline void from_json(const json & j, AdvancedIsolatedBuckBoost& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<IsolatedBuckBoostOperatingPoint>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedIsolatedBuckBoost & x) {
    j = json::object();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredInductance"] = x.get_desired_inductance();
}
} // namespace OpenMagnetics