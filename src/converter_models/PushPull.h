#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {


class PushPull : public MAS::PushPull {
public:
    bool _assertErrors = false;

    PushPull(const json& j);
    PushPull() {
    };

    bool run_checks(bool assert = false);

    Inputs process();
    DesignRequirements process_design_requirements();
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, PushPullOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, double outputInductance);
    double get_output_inductance(double mainSecondaryTurnsRatio);
    double calculate_duty_cycle();

};

class AdvancedPushPull : public PushPull {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;
    double desiredOutputInductance;

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

    const double & get_desired_output_inductance() const { return desiredOutputInductance; }
    double & get_mutable_desired_output_inductance() { return desiredOutputInductance; }
    void set_desired_output_inductance(const double & value) { this->desiredOutputInductance = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
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
    x.set_desired_output_inductance(j.at("desiredOutputInductance").get<double>());
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