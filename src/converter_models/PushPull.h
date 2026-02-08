#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

using namespace MAS;

namespace OpenMagnetics {

class PushPull : public MAS::PushPull, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;

public:
    bool _assertErrors = false;

    PushPull(const json& j);
    PushPull() {
    };
    
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    bool run_checks(bool assert = false) override;

    Inputs process();
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, PushPullOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, double outputInductance);
    double get_output_inductance(double mainSecondaryTurnsRatio);
    double get_maximum_duty_cycle();

    /**
     * @brief Generate an ngspice circuit for this Push-Pull converter
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
     * @brief Simulate the Push-Pull converter and extract operating points
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);
    
    /**
     * @brief Simulate and extract operating points from topology waveforms
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

};

class AdvancedPushPull : public PushPull {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;
    std::optional<double> desiredOutputInductance;

protected:
public:
    bool _assertErrors = false;

    AdvancedPushPull() = default;
    ~AdvancedPushPull() = default;

    AdvancedPushPull(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    std::optional<double> get_desired_output_inductance() const { return desiredOutputInductance; }
    void set_desired_output_inductance(std::optional<double> value) { this->desiredOutputInductance = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, AdvancedPushPull & x);
void to_json(json & j, const AdvancedPushPull & x);

inline void from_json(const json & j, AdvancedPushPull& x) {
    x.set_current_ripple_ratio(j.at("currentRippleRatio").get<double>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_duty_cycle(get_stack_optional<double>(j, "dutyCycle"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_drain_source_voltage(get_stack_optional<double>(j, "maximumDrainSourceVoltage"));
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<PushPullOperatingPoint>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
    x.set_desired_output_inductance(get_stack_optional<double>(j, "desiredOutputInductance"));
}

inline void to_json(json & j, const AdvancedPushPull & x) {
    j = json::object();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["dutyCycle"] = x.get_duty_cycle();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumDrainSourceVoltage"] = x.get_maximum_drain_source_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredInductance"] = x.get_desired_inductance();
    j["desiredOutputInductance"] = x.get_desired_output_inductance();
}
} // namespace OpenMagnetics