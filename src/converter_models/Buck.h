#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"


namespace OpenMagnetics {
using namespace MAS;


class Buck : public MAS::Buck, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;  // Number of periods to extract from simulation

public:
    bool _assertErrors = false;

    Buck(const json& j);
    Buck() {
    };

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, const BuckOperatingPoint& outputOperatingPoint, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    /**
     * @brief Generate an ngspice circuit for this Buck converter
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
     * @brief Simulate the Buck converter and extract operating points from waveforms
     * 
     * @param inductance Inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);
    
    /**
     * @brief Simulate and extract operating points from topology waveforms
     * 
     * Runs simulation and extracts operating points from the resulting waveforms.
     * Similar to simulate_and_extract_operating_points but returns raw OperatingPoint
     * objects for further processing.
     * 
     * @param inductance Inductance in H
     * @return Vector of OperatingPoints extracted from simulation waveforms
     */
    std::vector<OperatingPoint> simulate_and_extract_topology_waveforms(double inductance);
};

class AdvancedBuck : public Buck {
private:
    double desiredInductance;

protected:
public:
    AdvancedBuck() = default;
    ~AdvancedBuck() = default;

    AdvancedBuck(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

};


void from_json(const json & j, AdvancedBuck & x);
void to_json(json & j, const AdvancedBuck & x);


inline void from_json(const json & j, AdvancedBuck& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<BuckOperatingPoint>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedBuck & x) {
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