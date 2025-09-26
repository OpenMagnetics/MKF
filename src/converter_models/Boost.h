#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {


class Boost : public MAS::Boost {
public:
    bool _assertErrors = false;

    Boost(const json& j);
    Boost() {
    };

    bool run_checks(bool assert = false);

    Inputs process();
    DesignRequirements process_design_requirements();
    std::vector<OperatingPoint> process_operating_points(double magnetizingInductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, BoostOperatingPoint outputOperatingPoint, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency);

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
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<BoostOperatingPoint>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedBoost & x) {
    j = json::object();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["inputVoltage"] = x.get_input_voltage();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredInductance"] = x.get_desired_inductance();
}
} // namespace OpenMagnetics