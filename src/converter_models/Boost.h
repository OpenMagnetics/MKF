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
 * @brief Structure holding topology-level waveforms for Boost converter validation
 * 
 * These waveforms are used to validate that the simulation matches expected
 * converter behavior, not for magnetic component analysis.
 */
struct BoostTopologyWaveforms {
    // Time base
    std::vector<double> time;
    double frequency;
    
    // Input side signals
    std::vector<double> inputVoltage;           // v(vin_dc) - DC input voltage
    std::vector<double> switchNodeVoltage;      // v(sw) - switch node voltage
    
    // Output side signals  
    std::vector<double> inductorVoltage;        // v(vin_dc) - v(sw) - voltage across inductor
    std::vector<double> outputVoltage;          // v(vout) - DC output voltage
    
    // Currents
    std::vector<double> inductorCurrent;        // i(vl_sense) - inductor current
    
    // Metadata
    std::string operatingPointName;
    double inputVoltageValue;
    double outputVoltageValue;
    double dutyCycle;
};


class Boost : public MAS::Boost, public Topology {
private:
    int numPeriodsToExtract = 5;  // Number of periods to extract from simulation

public:
    bool _assertErrors = false;

    Boost(const json& j);
    Boost() {
    };

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) override;

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, BoostOperatingPoint outputOperatingPoint, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);
    double resolve_efficiency();

    /**
     * @brief Generate an ngspice circuit for this Boost converter
     * 
     * Uses the calculated design parameters (inductance, duty cycle)
     * to create a SPICE netlist that can be simulated.
     * 
     * @param inductance Inductance in H
     * @param inputVoltageIndex Which input voltage to use (0=nom, 1=min, 2=max)
     * @param operatingPointIndex Which operating point to simulate
     * @return SPICE netlist string
     */
    std::string generate_ngspice_circuit(
        double inductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);
    
    /**
     * @brief Simulate the Boost converter and extract operating points from waveforms
     * 
     * @param inductance Inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);
    
    /**
     * @brief Simulate and extract topology-level waveforms for converter validation
     * 
     * @param inductance Inductance in H
     * @return Vector of BoostTopologyWaveforms for each operating condition
     */
    std::vector<BoostTopologyWaveforms> simulate_and_extract_topology_waveforms(double inductance);
};

class AdvancedBoost : public Boost {
private:
    double desiredInductance;

protected:
public:
    bool _assertErrors = false;

    AdvancedBoost() = default;
    ~AdvancedBoost() = default;

    AdvancedBoost(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

};


void from_json(const json & j, AdvancedBoost & x);
void to_json(json & j, const AdvancedBoost & x);


inline void from_json(const json & j, AdvancedBoost& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<BoostOperatingPoint>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedBoost & x) {
    j = json::object();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredInductance"] = x.get_desired_inductance();
}
} // namespace OpenMagnetics