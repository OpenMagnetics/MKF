#pragma once

#include "InputsWrapper.h"
#include "InsulationMaterialWrapper.h"
#include "Defaults.h"
#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using json = nlohmann::json;

namespace OpenMagnetics {

class WireFieldPoint : public FieldPoint {
    public:
    WireFieldPoint() = default;
    virtual ~WireFieldPoint() = default;

    private:
    double length;
    double area;
    double lossDensity;

    public:
    const double & get_length() const { return length; }
    double & get_mutable_length() { return length; }
    void set_length(const double & length) { this->length = length; }

    const double & get_area() const { return area; }
    double & get_mutable_area() { return area; }
    void set_area(const double & area) { this->area = area; }

    const double & get_lossDensity() const { return lossDensity; }
    double & get_mutable_lossDensity() { return lossDensity; }
    void set_lossDensity(const double & lossDensity) { this->lossDensity = lossDensity; }
};

class WireMagneticStrengthFieldOutput : public WindingWindowMagneticStrengthFieldOutput {
    public:
    WireMagneticStrengthFieldOutput() = default;
    virtual ~WireMagneticStrengthFieldOutput() = default;

    private:
    std::vector<WireFieldPoint> data;

    public:
    const std::vector<WireFieldPoint> & get_data() const { return data; }
    std::vector<WireFieldPoint> & get_mutable_data() { return data; }
    void set_data(const std::vector<WireFieldPoint> & value) { this->data = value; }
};


class WireWrapper : public Wire {
    public:
        WireWrapper(Wire wire) {
            set_type(wire.get_type());
            if (wire.get_coating())
                set_coating(wire.get_coating().value());
            if (wire.get_conducting_diameter())
                set_conducting_diameter(wire.get_conducting_diameter().value());
            if (wire.get_conducting_height())
                set_conducting_height(wire.get_conducting_height().value());
            if (wire.get_conducting_width())
                set_conducting_width(wire.get_conducting_width().value());
            if (wire.get_manufacturer_info())
                set_manufacturer_info(wire.get_manufacturer_info().value());
            if (wire.get_material())
                set_material(wire.get_material().value());
            if (wire.get_name())
                set_name(wire.get_name().value());
            if (wire.get_number_conductors())
                set_number_conductors(wire.get_number_conductors().value());
            if (wire.get_outer_diameter())
                set_outer_diameter(wire.get_outer_diameter().value());
            if (wire.get_outer_height())
                set_outer_height(wire.get_outer_height().value());
            if (wire.get_outer_width())
                set_outer_width(wire.get_outer_width().value());
            if (wire.get_standard())
                set_standard(wire.get_standard().value());
            if (wire.get_standard_name())
                set_standard_name(wire.get_standard_name().value());
            if (wire.get_strand())
                set_strand(wire.get_strand().value());
            if (wire.get_conducting_area())
                set_conducting_area(wire.get_conducting_area().value());
            if (wire.get_edge_radius())
                set_edge_radius(wire.get_edge_radius().value());
        }

        WireWrapper(WireRound wire) {
            set_type(wire.get_type());
            if (wire.get_coating())
                set_coating(wire.get_coating().value());
            set_conducting_diameter(wire.get_conducting_diameter());
            if (wire.get_manufacturer_info())
                set_manufacturer_info(wire.get_manufacturer_info().value());
            if (wire.get_material())
                set_material(wire.get_material().value());
            if (wire.get_name())
                set_name(wire.get_name().value());
            if (wire.get_number_conductors())
                set_number_conductors(wire.get_number_conductors().value());
            if (wire.get_outer_diameter())
                set_outer_diameter(wire.get_outer_diameter().value());
            if (wire.get_standard())
                set_standard(wire.get_standard().value());
            if (wire.get_standard_name())
                set_standard_name(wire.get_standard_name().value());
            if (wire.get_conducting_area())
                set_conducting_area(wire.get_conducting_area().value());
        }


        WireWrapper() = default;
        virtual ~WireWrapper() = default;


        static WireRound convert_from_wire_to_strand(WireWrapper wire);
        static std::optional<InsulationWireCoating> resolve_coating(const WireWrapper& wire);
        std::optional<InsulationWireCoating> resolve_coating();
        static WireRound resolve_strand(const WireWrapper& wire);
        WireRound resolve_strand();
        static WireMaterial resolve_material(WireWrapper wire);
        static WireMaterial resolve_material(WireRound wire);
        WireMaterial resolve_material();
        static InsulationMaterialWrapper resolve_coating_insulation_material(WireWrapper wire);
        static InsulationMaterialWrapper resolve_coating_insulation_material(WireRound wire);
        InsulationMaterialWrapper resolve_coating_insulation_material();

        // Thought for enamelled round wires
        static double get_filling_factor_round(double conductingDiameter, int grade = 1, WireStandard standard = WireStandard::IEC_60317, bool includeAirInCell = false);
        static double get_outer_diameter_round(double conductingDiameter, int grade = 1, WireStandard standard = WireStandard::IEC_60317);
        // Thought for insulated round wires
        static double get_filling_factor_round(double conductingDiameter, int numberLayers, double thicknessLayers, WireStandard standard = WireStandard::IEC_60317, bool includeAirInCell = false);
        static double get_outer_diameter_round(double conductingDiameter, int numberLayers, double thicknessLayers, WireStandard standard = WireStandard::IEC_60317);
        // Thought for enamelled litz wires with or without serving
        static double get_serving_thickness_from_standard(int numberLayers, double outerDiameter);
        static double get_filling_factor_served_litz(double conductingDiameter, int numberConductors, int grade = 1, int numberLayers= 1, WireStandard standard = WireStandard::IEC_60317, bool includeAirInCell = false);
        static double get_outer_diameter_bare_litz(double conductingDiameter, int numberConductors, int grade = 1, WireStandard standard = WireStandard::IEC_60317);
        static double get_outer_diameter_served_litz(double conductingDiameter, int numberConductors, int grade = 1, int numberLayers= 1, WireStandard standard = WireStandard::IEC_60317);
        // Thought for insulated litz wires
        static double get_filling_factor_insulated_litz(double conductingDiameter, int numberConductors, int numberLayers, double thicknessLayers, int grade = 1, WireStandard standard = WireStandard::IEC_60317, bool includeAirInCell = false);
        static double get_outer_diameter_insulated_litz(double conductingDiameter, int numberConductors, int numberLayers, double thicknessLayers, int grade = 1, WireStandard standard = WireStandard::IEC_60317);
        // Thought for enamelled reactangular wires
        static double get_filling_factor_rectangular(double conductingWidth, double conductingHeight, int grade = 1, WireStandard standard = WireStandard::IEC_60317);
        static double get_conducting_area_rectangular(double conductingWidth, double conductingHeight, WireStandard standard = WireStandard::IEC_60317);
        static double get_outer_width_rectangular(double conductingWidth, int grade = 1, WireStandard standard = WireStandard::IEC_60317);
        static double get_outer_height_rectangular(double conductingHeight, int grade = 1, WireStandard standard = WireStandard::IEC_60317);

        void set_nominal_value_conducting_diameter(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_conducting_diameter(aux);
        }

        void set_nominal_value_conducting_height(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_conducting_height(aux);
        }

        void set_nominal_value_conducting_width(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_conducting_width(aux);
        }

        void set_nominal_value_outer_diameter(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_outer_diameter(aux);
        }

        void set_nominal_value_outer_height(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_outer_height(aux);
        }

        void set_nominal_value_outer_width(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_outer_width(aux);
        }

        int get_equivalent_insulation_layers(double voltageToInsulate);

        std::vector<std::vector<double>> calculate_current_density_distribution(SignalDescriptor current, double frequency, double temperature, size_t numberPoints=10);
        double calculate_dc_current_density(OperatingPointExcitation excitation);
        double calculate_dc_current_density(SignalDescriptor current);
        double calculate_effective_current_density(OperatingPointExcitation excitation, double temperature);
        double calculate_effective_current_density(SignalDescriptor current, double temperature);
        double calculate_effective_current_density(double rms, double frequency, double temperature);
        double calculate_effective_conducting_area(double frequency, double temperature);
        double calculate_conducting_area();

        static int calculate_number_parallels_needed(InputsWrapper inputs, WireWrapper& wire, double maximumEffectiveCurrentDensity, size_t windingIndex = 0);
        static int calculate_number_parallels_needed(OperatingPointExcitation excitation, double temperature, WireWrapper& wire, double maximumEffectiveCurrentDensity);
        static int calculate_number_parallels_needed(SignalDescriptor current, double temperature, WireWrapper& wire, double maximumEffectiveCurrentDensity);
        static int calculate_number_parallels_needed(double rms, double effectiveFrequency, double temperature, WireWrapper& wire, double maximumEffectiveCurrentDensity);

        double get_maximum_outer_width();
        double get_maximum_outer_height();
        double get_maximum_outer_dimension();
        double get_maximum_conducting_width();
        double get_maximum_conducting_height();
        double get_minimum_conducting_dimension();


        double get_coating_thickness();
        double get_coating_dielectric_strength();
        static double get_coating_thickness(WireWrapper wire);
        static double get_coating_dielectric_strength(WireWrapper wire);
        std::string encode_coating_label();
        static std::string encode_coating_label(WireWrapper wire);
        static std::optional<InsulationWireCoating> decode_coating_label(std::string label);
        double get_coating_relative_permittivity();
        static double get_coating_relative_permittivity(WireWrapper wire);

        void cut_foil_wire_to_section(Section section);
        void cut_planar_wire_to_section(Section section);
        double get_relative_cost();

        static WireWrapper get_equivalent_wire(WireWrapper oldWire, WireType newWireType, double effectiveFrequency=200000, double temperature=Defaults().ambientTemperature);
        static WireWrapper get_wire_for_frequency(double effectiveFrequency, double temperature=Defaults().ambientTemperature, bool exact=false);
        static WireWrapper get_wire_for_conducting_area(double conductingArea, double temperature=Defaults().ambientTemperature, bool exact=false);

};
} // namespace OpenMagnetics
