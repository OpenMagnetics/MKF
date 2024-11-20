#pragma once
#include "Constants.h"
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
    TemperatureNode firstNode;
    TemperatureNode secondNode;
    std::vector<double> thermalResistances;
    ThermalResistanceConnectionTypes connectionType;

public:

    const TemperatureNode & get_first_node() const { return firstNode; }
    TemperatureNode & get_mutable_first_node() { return firstNode; }
    void set_first_node(const TemperatureNode & value) { this->firstNode = value; }

    const TemperatureNode & get_second_node() const { return secondNode; }
    TemperatureNode & get_mutable_second_node() { return secondNode; }
    void set_second_node(const TemperatureNode & value) { this->secondNode = value; }

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
};


} // namespace OpenMagnetics