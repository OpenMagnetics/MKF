#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {


class IsolatedBuck : public MAS::IsolatedBuck {
public:
    bool _assertErrors = false;

    IsolatedBuck(const json& j);
    IsolatedBuck() {
    };

    bool run_checks(bool assert = false);

    Inputs process();
    DesignRequirements process_design_requirements();
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, IsolatedBuckOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency);

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
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<IsolatedBuckOperatingPoint>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedIsolatedBuck & x) {
    j = json::object();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["inputVoltage"] = x.get_input_voltage();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredInductance"] = x.get_desired_inductance();
}
} // namespace OpenMagnetics