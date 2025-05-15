#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/InputsWrapper.h"
#include "constructive_models/MagneticWrapper.h"


namespace OpenMagnetics {

class Topology {
private:
protected:
public:
    Topology() = default;
    ~Topology() = default;
};

class Flyback : public Topology {
public:

    enum class Modes : int {
        ContinuousConductionMode,
        DiscontinuousConductionMode,
        QuasiResonantMode,
        BoundaryModeOperation,
    };

    class FlybackOperatingPoint {
    private:
        std::vector<double> outputVoltages;
        std::vector<double> outputCurrents;
        std::optional<double> switchingFrequency;
        std::optional<Flyback::Modes> mode;
        double ambientTemperature;

    public:
        FlybackOperatingPoint() = default;
        ~FlybackOperatingPoint() = default;

        const std::vector<double> & get_output_voltages() const { return outputVoltages; }
        std::vector<double> & get_mutable_output_voltages() { return outputVoltages; }
        void set_output_voltages(const std::vector<double> & value) { this->outputVoltages = value; }
        
        const std::vector<double> & get_output_currents() const { return outputCurrents; }
        std::vector<double> & get_mutable_output_currents() { return outputCurrents; }
        void set_output_currents(const std::vector<double> & value) { this->outputCurrents = value; }
    
        std::optional<double> get_switching_frequency() const { return switchingFrequency; }
        void set_switching_frequency(std::optional<double> value) { this->switchingFrequency = value; }

        std::optional<Flyback::Modes> get_mode() const { return mode; }
        void set_mode(std::optional<Flyback::Modes> value) { this->mode = value; }
    
        const double & get_ambient_temperature() const { return ambientTemperature; }
        double & get_mutable_ambient_temperature() { return ambientTemperature; }
        void set_ambient_temperature(const double & value) { this->ambientTemperature = value; }

        double resolve_switching_frequency(double inputVoltage, double diodeVoltageDrop, std::optional<double> inductance = std::nullopt, std::optional<std::vector<double>> turnsRatios = std::nullopt, double efficiency = 0.85);
        Flyback::Modes resolve_mode(std::optional<double> currentRippleRatio = std::nullopt);
    };

private:

    DimensionWithTolerance inputVoltage;
    double diodeVoltageDrop;
    std::optional<double> maximumDrainSourceVoltage = 600;
    std::optional<double> maximumDutyCycle = 0.5;
    double currentRippleRatio;
    std::vector<FlybackOperatingPoint> operatingPoints;
    double efficiency = 1;

protected:
public:
    bool _assertErrors = false;

    Flyback() = default;
    ~Flyback() = default;

    Flyback(const json& j);

    const DimensionWithTolerance & get_input_voltage() const { return inputVoltage; }
    DimensionWithTolerance & get_mutable_input_voltage() { return inputVoltage; }
    void set_input_voltage(const DimensionWithTolerance & value) { this->inputVoltage = value; }
    
    const double & get_diode_voltage_drop() const { return diodeVoltageDrop; }
    double & get_mutable_diode_voltage_drop() { return diodeVoltageDrop; }
    void set_diode_voltage_drop(const double & value) { this->diodeVoltageDrop = value; }

    std::optional<double> get_maximum_duty_cycle() const { return maximumDutyCycle; }
    void set_maximum_duty_cycle(std::optional<double> value) { this->maximumDutyCycle = value; }

    std::optional<double> get_maximum_drain_source_voltage() const { return maximumDrainSourceVoltage; }
    void set_maximum_drain_source_voltage(std::optional<double> value) { this->maximumDrainSourceVoltage = value; }

    const double & get_current_ripple_ratio() const { return currentRippleRatio; }
    double & get_mutable_current_ripple_ratio() { return currentRippleRatio; }
    void set_current_ripple_ratio(const double & value) { this->currentRippleRatio = value; }

    const std::vector<FlybackOperatingPoint> & get_operating_points() const { return operatingPoints; }
    std::vector<FlybackOperatingPoint> & get_mutable_operating_points() { return operatingPoints; }
    void set_operating_points(const std::vector<FlybackOperatingPoint> & value) { this->operatingPoints = value; }
    
    const double & get_efficiency() const { return efficiency; }
    double & get_mutable_efficiency() { return efficiency; }
    void set_efficiency(const double & value) { this->efficiency = value; }

    bool run_checks(bool assert = false);

    // According to Worked Example (7), pages 135-144 — Designing the Flyback Transformer of Switching Power Supplies A - Z (Second Edition) by Sanjaya Maniktala
    InputsWrapper process();
    InputsWrapper process(MagneticWrapper magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, Flyback::FlybackOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, std::optional<Flyback::Modes> customMode=std::nullopt, std::optional<double> customDutyCycle=std::nullopt, std::optional<double> customDeadTime=std::nullopt);
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

    InputsWrapper process();
    InputsWrapper process(MagneticWrapper magnetic);

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


void from_json(const json & j, Flyback::FlybackOperatingPoint & x);
void to_json(json & j, const Flyback::FlybackOperatingPoint & x);

inline void from_json(const json & j, Flyback::FlybackOperatingPoint& x) {
    x.set_output_voltages(j.at("outputVoltages").get<std::vector<double>>());
    x.set_output_currents(j.at("outputCurrents").get<std::vector<double>>());
    x.set_switching_frequency(get_stack_optional<double>(j, "switchingFrequency"));
    x.set_mode(get_stack_optional<Flyback::Modes>(j, "mode"));
    x.set_ambient_temperature(j.at("ambientTemperature").get<double>());
}

inline void to_json(json & j, const Flyback::FlybackOperatingPoint & x) {
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
    x.set_operating_points(j.at("operatingPoints").get<std::vector<Flyback::FlybackOperatingPoint>>());
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
    x.set_operating_points(j.at("operatingPoints").get<std::vector<AdvancedFlyback::FlybackOperatingPoint>>());
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