#pragma once
#include "Constants.h"
#include "constructive_models/CoreWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class TemperatureNode {
public:
    TemperatureNode() = default;
    virtual ~TemperatureNode() = default;
private:
    std::string name;
    std::vector<double> coordinates;

public:

    const std::string & get_name() const { return name; }
    std::string & get_mutable_name() { return name; }
    void set_name(const std::string & value) { this->name = value; }

    const std::vector<double> & get_coordinates() const { return coordinates; }
    std::vector<double> & get_mutable_coordinates() { return coordinates; }
    void set_coordinates(const std::vector<double> & value) { this->coordinates = value; }
};

class TemperatureNodeConnection {
public:
    enum class ThermalResistanceConnectionTypes : int { SERIES, PARALLEL };
    TemperatureNodeConnection() = default;
    virtual ~TemperatureNodeConnection() = default;

private:
    TemperatureNode* firstNode;
    TemperatureNode* secondNode;
    std::vector<double> thermalResistances;
    ThermalResistanceConnectionTypes connectionType;

public:

    TemperatureNode* get_first_node() { return firstNode; } 
    void set_first_node(TemperatureNode* & value) { this->firstNode = value; }

    TemperatureNode* get_second_node() { return secondNode; }
    void set_second_node(TemperatureNode* & value) { this->secondNode = value; }

    const std::vector<double> & get_thermal_resistances() const { return thermalResistances; }
    std::vector<double> & get_mutable_thermal_resistances() { return thermalResistances; }
    void set_thermal_resistances(const std::vector<double> & value) { this->thermalResistances = value; }

    ThermalResistanceConnectionTypes get_layer_purpose() const { return connectionType; }
    void set_layer_purpose(ThermalResistanceConnectionTypes value) { this->connectionType = value; }
};

class Temperature {
private:
protected:
public:
    Temperature() = default;
    ~Temperature() = default;

    static double calculate_temperature_from_core_thermal_resistance(CoreWrapper core, double totalLosses);
    static double calculate_temperature_from_core_thermal_resistance(double thermalResistance, double totalLosses);
};


} // namespace OpenMagnetics