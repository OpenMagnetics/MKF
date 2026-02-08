#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"
#include "converter_models/ForwardConverterUtils.h"

using namespace MAS;

namespace OpenMagnetics {

class ActiveClampForward : public MAS::Forward, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;
    
public:
    bool _assertErrors = false;

    ActiveClampForward(const json& j);
    ActiveClampForward() {
    };
    
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);
    double get_total_reflected_secondary_current(ForwardOperatingPoint forwardOperatingPoint, std::vector<double> turnsRatios, double rippleRatio=1);

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, ForwardOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, double mainOutputInductance);
    double get_output_inductance(double mainSecondaryTurnsRatio, size_t outputIndex);
    double get_maximum_duty_cycle();

    /**
     * @brief Generate an ngspice circuit for this Active Clamp Forward converter
     */
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);
    
    /**
     * @brief Simulate the Active Clamp Forward converter and extract operating points
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);
    
    /**
     * @brief Simulate and extract topology-level waveforms for converter validation
     */
    std::vector<OperatingPoint> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

};

class AdvancedActiveClampForward : public ActiveClampForward {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;
    std::optional<std::vector<double>> desiredOutputInductances;

protected:
public:
    bool _assertErrors = false;

    AdvancedActiveClampForward() = default;
    ~AdvancedActiveClampForward() = default;

    AdvancedActiveClampForward(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

    std::optional<std::vector<double>> get_desired_output_inductances() const { return desiredOutputInductances; }
    void set_desired_output_inductances(std::optional<std::vector<double>> value) { this->desiredOutputInductances = value; }
};


void from_json(const json & j, AdvancedActiveClampForward & x);
void to_json(json & j, const AdvancedActiveClampForward & x);

inline void from_json(const json & j, AdvancedActiveClampForward& x) {
    x.set_current_ripple_ratio(j.at("currentRippleRatio").get<double>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_duty_cycle(get_stack_optional<double>(j, "dutyCycle"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<ForwardOperatingPoint>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
    x.set_desired_output_inductances(get_stack_optional<std::vector<double>>(j, "desiredOutputInductance"));
}

inline void to_json(json & j, const AdvancedActiveClampForward & x) {
    j = json::object();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["dutyCycle"] = x.get_duty_cycle();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredInductance"] = x.get_desired_inductance();
    j["desiredOutputInductances"] = x.get_desired_output_inductances();
}
} // namespace OpenMagnetics