#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {

class Topology {
private:
protected:
public:
    Topology() = default;
    ~Topology() = default;
};

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

/** ------------------
 *  Inverter topology
 *  ------------------
 * */

class MyInverter  : public MAS::Inverter {
  private:
    std::complex<double>  Zgrid;
    std::complex<double>  Zfilter;
    struct ABCVoltages {
        double Va;
        double Vb;
        double Vc;
    };
    struct PwmSignals {
        bool Sa;  // Leg A upper switch
        bool Sb;  // Leg B upper switch
        bool Sc;  // Leg C upper switch
    };
    struct Inductor1State {
        std::complex<double> vL1;  // Voltage across L1
        std::complex<double> iL1;  // Current through L1
    };

  public:
    bool _assertErrors = false;

    MyInverter() = default;
    MyInverter(const json& j);

    const std::vector<TwoLevelInverterOperatingPoint>& get_operating_points() const { return operatingPoints; }
    void set_operating_points(const std::vector<TwoLevelInverterOperatingPoint>& value) {
        this->operatingPoints = value;
    }

    bool run_checks(bool assert = false);
    compute_inverter_reference(const InverterLoad& load,
                                    const InverterOperatingPoint& op,
                                    const DimensionWithTolerance& vdc,
                                    double power);
    
    // ---- Impedances ----
    Impedance compute_load_impedance(const InverterLoad& load);
    Impedance compute_filter_impedance(const InverterDownstreamFilter& filter);

    
    SignalDescriptor compute_switching_node_voltage(const Modulation& modulation,
                                                    const Signal& inverter_voltage_reference);


    // ---- Currents ----
    SignalDescriptor compute_inductor_current();

    Inputs process();
    DesignRequirements process_design_requirements();
    std::vector<OperatingPoint> process_operating_points();
};


inline void from_json(const json& j, MAS::InverterOperatingPoint& x) {
    x.set_fundamental_frequency(j.at("fundamentalFrequency").get<double>());
    x.set_output_power(get_stack_optional<double>(j, "outputPower"));
    x.set_operating_temperature(get_stack_optional<double>(j, "operatingTemperature"));
    x.set_power_factor(get_stack_optional<double>(j, "powerFactor"));
    x.set_current_phase_angle(get_stack_optional<double>(j, "currentPhaseAngle"));
    x.set_load(j.at("load").get<MAS::InverterLoad>());
}

inline void to_json(json& j, const MAS::InverterOperatingPoint& x) {
    j = json::object();
    j["fundamentalFrequency"] = x.get_fundamental_frequency();
    j["outputPower"] = x.get_output_power();
    j["operatingTemperature"] = x.get_operating_temperature();
    j["powerFactor"] = x.get_power_factor();
    j["currentPhaseAngle"] = x.get_current_phase_angle();
    j["load"] = x.get_load();
}

inline void from_json(const json& j, MyInverter& x) {
    x.set_dc_bus_voltage(j.at("dcBusVoltage").get<DimensionWithTolerance>());
    x.set_vdc_ripple(get_stack_optional<DimensionWithTolerance>(j, "vdcRipple"));
    x.set_inverter_leg_power_rating(j.at("inverterLegPowerRating").get<double>());
    x.set_line_rms_current(j.at("lineRmsCurrent").get<DimensionWithTolerance>());
    x.set_downstream_filter(get_stack_optional<MAS::InverterDownstreamFilter>(j, "downstreamFilter"));
    x.set_modulation(get_stack_optional<MAS::Modulation>(j, "modulation"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<MAS::InverterOperatingPoint>>());
}

inline void to_json(json& j, const MyInverter& x) {
    j = json::object();
    j["dcBusVoltage"] = x.get_dc_bus_voltage();
    j["vdcRipple"] = x.get_vdc_ripple();
    j["inverterLegPowerRating"] = x.get_inverter_leg_power_rating();
    j["lineRmsCurrent"] = x.get_line_rms_current();
    j["downstreamFilter"] = x.get_downstream_filter();
    j["modulation"] = x.get_modulation();
    j["operatingPoints"] = x.get_operating_points();
}

inline void from_json(const json& j, MAS::InverterDownstreamFilter& x) {
    // filterTopology
    x.set_filter_topology(j.at("filterTopology").get<MAS::FilterTopologies>());

    // inductor (required)
    x.set_inductor(j.at("inductor").get<MAS::Inductor>());

    // capacitor (optional, required only for LC or LCL)
    if (j.contains("capacitor")) {
        x.set_capacitor(j.at("capacitor").get<MAS::Capacitor>());
    }

    // inductor2 (optional, required only for LCL)
    if (j.contains("inductor2")) {
        x.set_inductor2(j.at("inductor2").get<MAS::Inductor2>());
    }
}

inline void to_json(json& j, const MAS::InverterDownstreamFilter& x) {
    j = json::object();

    j["filterTopology"] = x.get_filter_topology();
    j["inductor"] = x.get_inductor();

    if (x.get_capacitor().has_value()) {
        j["capacitor"] = x.get_capacitor().value();
    }

    if (x.get_inductor2().has_value()) {
        j["inductor2"] = x.get_inductor2().value();
    }
}

inline void from_json(const json& j, MAS::InverterLoad& x) {
    x.set_load_type(j.at("loadType").get<MAS::LoadType>());
    x.set_resistance(get_stack_optional<double>(j, "resistance"));
    x.set_inductance(get_stack_optional<double>(j, "inductance"));
    x.set_phase_voltage(get_stack_optional<double>(j, "phaseVoltage"));
    x.set_grid_frequency(get_stack_optional<double>(j, "gridFrequency"));
    x.set_grid_resistance(get_stack_optional<double>(j, "gridResistance"));
    x.set_grid_inductance(get_stack_optional<double>(j, "gridInductance"));
}

inline void to_json(json& j, const MAS::InverterLoad& x) {
    j = json::object();
    j["loadType"] = x.get_load_type();
    j["resistance"] = x.get_resistance();
    j["inductance"] = x.get_inductance();
    j["phaseVoltage"] = x.get_phase_voltage();
    j["gridFrequency"] = x.get_grid_frequency();
    j["gridResistance"] = x.get_grid_resistance();
    j["gridInductance"] = x.get_grid_inductance();
}

inline void from_json(const json& j, MAS::Modulation& x) {
    x.set_switching_frequency(j.at("switchingFrequency").get<double>());
    x.set_pwm_type(j.at("pwmType").get<MAS::PwmType>());
    x.set_modulation_strategy(j.at("modulationStrategy").get<MAS::ModulationStrategy>());
    x.set_modulation_depth(j.at("modulationDepth").get<double>());
    x.set_rise_time(get_stack_optional<double>(j, "riseTime"));
    x.set_deadtime(get_stack_optional<double>(j, "deadtime"));
    x.set_third_harmonic_injection_coefficient(
        get_stack_optional<double>(j, "thirdHarmonicInjectionCoefficient"));
}

inline void to_json(json& j, const MAS::Modulation& x) {
    j = json::object();
    j["switchingFrequency"] = x.get_switching_frequency();
    j["pwmType"] = x.get_pwm_type();
    j["modulationStrategy"] = x.get_modulation_strategy();
    j["modulationDepth"] = x.get_modulation_depth();
    j["riseTime"] = x.get_rise_time();
    j["deadtime"] = x.get_deadtime();
    j["thirdHarmonicInjectionCoefficient"] = x.get_third_harmonic_injection_coefficient();
}

inline void from_json(const json& j, MAS::ModulationStrategy& x) {
    std::string s = j.get<std::string>();
    if (s == "SPWM")    x = MAS::ModulationStrategy::SPWM;
    else if (s == "SVPWM") x = MAS::ModulationStrategy::SVPWM;
    else if (s == "THIPWM") x = MAS::ModulationStrategy::THIPWM;
    else throw std::invalid_argument("Invalid ModulationStrategy: " + s);
}

inline void to_json(json& j, const MAS::ModulationStrategy& x) {
    switch (x) {
        case MAS::ModulationStrategy::SPWM:   j = "SPWM"; break;
        case MAS::ModulationStrategy::SVPWM:  j = "SVPWM"; break;
        case MAS::ModulationStrategy::THIPWM: j = "THIPWM"; break;
    }
}

inline void from_json(const json& j, MAS::PwmType& x) {
    std::string s = j.get<std::string>();
    if (s == "sawtooth")   x = MAS::PwmType::SAWTOOTH;
    else if (s == "triangular") x = MAS::PwmType::TRIANGULAR;
    else throw std::invalid_argument("Invalid PwmType: " + s);
}

inline void to_json(json& j, const MAS::PwmType& x) {
    switch (x) {
        case MAS::PwmType::SAWTOOTH:   j = "sawtooth"; break;
        case MAS::PwmType::TRIANGULAR: j = "triangular"; break;
    }
}

inline void from_json(const json& j, MAS::LoadType& x) {
    std::string s = j.get<std::string>();
    if (s == "grid")   x = MAS::LoadType::GRID;
    else if (s == "R-L") x = MAS::LoadType::R_L;
    else throw std::invalid_argument("Invalid LoadType: " + s);
}

inline void to_json(json& j, const MAS::LoadType& x) {
    switch (x) {
        case MAS::LoadType::GRID: j = "grid"; break;
        case MAS::LoadType::R_L:  j = "R-L"; break;
    }
}

inline void from_json(const json& j, MAS::FilterTopologies& x) {
    std::string s = j.get<std::string>();
    if (s == "L")       x = MAS::FilterTopologies::L;
    else if (s == "LC") x = MAS::FilterTopologies::LC;
    else if (s == "LCL")x = MAS::FilterTopologies::LCL;
    else throw std::invalid_argument("Invalid FilterTopologies: " + s);
}

inline void to_json(json& j, const MAS::FilterTopologies& x) {
    switch (x) {
        case MAS::FilterTopologies::L:   j = "L"; break;
        case MAS::FilterTopologies::LC:  j = "LC"; break;
        case MAS::FilterTopologies::LCL: j = "LCL"; break;
    }
}


} // namespace OpenMagnetics