#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {

class FlybackOperatingPoint : public MAS::FlybackOperatingPoint{
public:
    double resolve_switching_frequency(double inputVoltage, double diodeVoltageDrop, std::optional<double> inductance = std::nullopt, std::optional<std::vector<double>> turnsRatios = std::nullopt, double efficiency = 0.85);
    FlybackModes resolve_mode(std::optional<double> currentRippleRatio = std::nullopt);
};


class Flyback : public MAS::Flyback {
private:
    std::optional<double> maximumDrainSourceVoltage = 600;
    std::vector<OpenMagnetics::FlybackOperatingPoint> operatingPoints;
    std::optional<double> maximumDutyCycle = 0.5;
    double efficiency = 1;

public:
    const std::vector<OpenMagnetics::FlybackOperatingPoint> & get_operating_points() const { return operatingPoints; }
    std::vector<OpenMagnetics::FlybackOperatingPoint> & get_mutable_operating_points() { return operatingPoints; }
    void set_operating_points(const std::vector<OpenMagnetics::FlybackOperatingPoint> & value) { this->operatingPoints = value; }


protected:
public:
    bool _assertErrors = false;

    Flyback(const json& j);
    Flyback() {
        maximumDrainSourceVoltage = 600;
        maximumDutyCycle = 0.5;
        efficiency = 1;
    };

    bool run_checks(bool assert = false);

    // According to Worked Example (7), pages 135-144 — Designing the Flyback Transformer of Switching Power Supplies A - Z (Second Edition) by Sanjaya Maniktala
    Inputs process();
    DesignRequirements process_design_requirements();
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, OpenMagnetics::FlybackOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, std::optional<FlybackModes> customMode=std::nullopt, std::optional<double> customDutyCycle=std::nullopt, std::optional<double> customDeadTime=std::nullopt);
    static double get_total_input_power(std::vector<double> outputCurrents, std::vector<double> outputVoltages, double efficiency, double diodeVoltageDrop);
    static double get_total_input_power(double outputCurrent, double outputVoltage, double efficiency, double diodeVoltageDrop);
    static double get_minimum_output_reflected_voltage(double maximumDrainSourceVoltage, double maximumInputVoltage, double safetyMargin=0.85);

};

class AdvancedFlyback : public Flyback {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;
    std::vector<std::vector<double>> desiredDutyCycle;
    std::optional<std::vector<double>> desiredDeadTime;

protected:
public:
    bool _assertErrors = false;

    AdvancedFlyback() = default;
    ~AdvancedFlyback() = default;

    AdvancedFlyback(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const std::vector<std::vector<double>> & get_desired_duty_cycle() const { return desiredDutyCycle; }
    std::vector<std::vector<double>> & get_mutable_desired_duty_cycle() { return desiredDutyCycle; }
    void set_desired_duty_cycle(const std::vector<std::vector<double>> & value) { this->desiredDutyCycle = value; }

    std::optional<std::vector<double>> get_desired_dead_time() const { return desiredDeadTime; }
    void set_desired_dead_time(std::optional<std::vector<double>> value) { this->desiredDeadTime = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, FlybackOperatingPoint & x);
void to_json(json & j, const FlybackOperatingPoint & x);

inline void from_json(const json & j, FlybackOperatingPoint& x) {
    x.set_output_voltages(j.at("outputVoltages").get<std::vector<double>>());
    x.set_output_currents(j.at("outputCurrents").get<std::vector<double>>());
    x.set_switching_frequency(get_stack_optional<double>(j, "switchingFrequency"));
    x.set_mode(get_stack_optional<FlybackModes>(j, "mode"));
    x.set_ambient_temperature(j.at("ambientTemperature").get<double>());
}

inline void to_json(json & j, const FlybackOperatingPoint & x) {
    j = json::object();
    j["outputVoltages"] = x.get_output_voltages();
    j["outputCurrents"] = x.get_output_currents();
    j["switchingFrequency"] = x.get_switching_frequency();
    j["mode"] = x.get_mode();
    j["ambientTemperature"] = x.get_ambient_temperature();
}

void from_json(const json & j, Flyback & x);
void to_json(json & j, const Flyback & x);

inline void from_json(const json & j, Flyback& x) {
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_maximum_drain_source_voltage(get_stack_optional<double>(j, "maximumDrainSourceVoltage"));
    x.set_maximum_duty_cycle(get_stack_optional<double>(j, "maximumDutyCycle"));
    x.set_current_ripple_ratio(j.at("currentRippleRatio").get<double>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<FlybackOperatingPoint>>());
    x.set_efficiency(j.at("efficiency").get<double>());
}

inline void to_json(json & j, const Flyback & x) {
    j = json::object();
    j["inputVoltage"] = x.get_input_voltage();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["maximumDrainSourceVoltage"] = x.get_maximum_drain_source_voltage();
    j["maximumDutyCycle"] = x.get_maximum_duty_cycle();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["operatingPoints"] = x.get_operating_points();
    j["efficiency"] = x.get_efficiency();
}

void from_json(const json & j, AdvancedFlyback & x);
void to_json(json & j, const AdvancedFlyback & x);

inline void from_json(const json & j, AdvancedFlyback& x) {
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
    x.set_desired_dead_time(get_stack_optional<std::vector<double>>(j, "desiredDeadTime"));
    x.set_desired_duty_cycle(j.at("desiredDutyCycle").get<std::vector<std::vector<double>>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<FlybackOperatingPoint>>());
    x.set_efficiency(j.at("efficiency").get<double>());
    x.set_current_ripple_ratio(std::numeric_limits<double>::quiet_NaN());
}

inline void to_json(json & j, const AdvancedFlyback & x) {
    j = json::object();
    j["inputVoltage"] = x.get_input_voltage();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["desiredInductance"] = x.get_desired_inductance();
    j["desiredDutyCycle"] = x.get_desired_duty_cycle();
    j["desiredDeadTime"] = x.get_desired_dead_time();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["operatingPoints"] = x.get_operating_points();
    j["efficiency"] = x.get_efficiency();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
}
} // namespace OpenMagnetics