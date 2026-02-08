#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"

using namespace MAS;

namespace OpenMagnetics {


class IsolatedBuck : public MAS::IsolatedBuck, public Topology {
public:
    bool _assertErrors = false;

    IsolatedBuck(const json& j);
    IsolatedBuck() {
    };

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, IsolatedBuckOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency);
    std::vector<OperatingPoint> simulate_and_extract_topology_waveforms(const std::vector<double>& turnsRatios, double magnetizingInductance);

};

class AdvancedIsolatedBuck : public IsolatedBuck {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;

protected:
public:
    bool _assertErrors = false;

    AdvancedIsolatedBuck() = default;
    ~AdvancedIsolatedBuck() = default;

    AdvancedIsolatedBuck(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, AdvancedIsolatedBuck & x);
void to_json(json & j, const AdvancedIsolatedBuck & x);


inline void from_json(const json & j, AdvancedIsolatedBuck& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<IsolatedBuckOperatingPoint>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedIsolatedBuck & x) {
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