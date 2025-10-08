#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"

using namespace MAS;

namespace OpenMagnetics {


class TwoSwitchForward : public MAS::Forward, public Topology {
public:
    bool _assertErrors = false;

    TwoSwitchForward(const json& j);
    TwoSwitchForward() {
    };

    bool run_checks(bool assert = false);

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);
    double get_total_reflected_secondary_current(ForwardOperatingPoint forwardOperatingPoint, std::vector<double> turnsRatios, double rippleRatio=1);

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, ForwardOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, double mainOutputInductance);
    double get_output_inductance(double mainSecondaryTurnsRatio, size_t outputIndex);
    double get_maximum_duty_cycle();

};

class AdvancedTwoSwitchForward : public TwoSwitchForward {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;
    std::optional<std::vector<double>> desiredOutputInductances;

protected:
public:
    bool _assertErrors = false;

    AdvancedTwoSwitchForward() = default;
    ~AdvancedTwoSwitchForward() = default;

    AdvancedTwoSwitchForward(const json& j);

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


void from_json(const json & j, AdvancedTwoSwitchForward & x);
void to_json(json & j, const AdvancedTwoSwitchForward & x);

inline void from_json(const json & j, AdvancedTwoSwitchForward& x) {
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

inline void to_json(json & j, const AdvancedTwoSwitchForward & x) {
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