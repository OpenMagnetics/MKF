#pragma once

#include "constructive_models/InsulationMaterial.h"
#include "processors/Inputs.h"
#include "Defaults.h"
#include "json.hpp"
#include "spline.h"

#include <MAS.hpp>
#include <vector>

using namespace MAS;

using json = nlohmann::json;

namespace OpenMagnetics {
    
inline std::map<std::string, tk::spline> wireCoatingThicknessProportionInterps;
inline std::map<std::string, tk::spline> wireFillingFactorInterps;
inline std::map<std::string, tk::spline> wirePackingFactorInterps;
inline std::map<std::string, tk::spline> wireConductingAreaProportionInterps;
inline std::map<std::string, double> minWireConductingDimensions;
inline std::map<std::string, double> maxWireConductingDimensions;
inline std::map<std::string, int64_t> minLitzWireNumberConductors;
inline std::map<std::string, int64_t> maxLitzWireNumberConductors;

class WireSolidInsulationRequirements {
    public:
    WireSolidInsulationRequirements() = default;
    virtual ~WireSolidInsulationRequirements() = default;

    private:
    std::optional<int64_t> _minimumNumberLayers;
    std::optional<int64_t> _maximumNumberLayers;
    std::optional<int64_t> _minimumGrade;
    std::optional<int64_t> _maximumGrade;
    double minimumBreakdownVoltage = 0;

    public:
    std::optional<int64_t> get_minimum_number_layers() const { return _minimumNumberLayers; }
    void set_minimum_number_layers(std::optional<int64_t> value) { this->_minimumNumberLayers = value; }

    std::optional<int64_t> get_minimum_grade() const { return _minimumGrade; }
    void set_minimum_grade(std::optional<int64_t> value) { this->_minimumGrade = value; }

    std::optional<int64_t> get_maximum_number_layers() const { return _maximumNumberLayers; }
    void set_maximum_number_layers(std::optional<int64_t> value) { this->_maximumNumberLayers = value; }

    std::optional<int64_t> get_maximum_grade() const { return _maximumGrade; }
    void set_maximum_grade(std::optional<int64_t> value) { this->_maximumGrade = value; }

    double get_minimum_breakdown_voltage() const { return minimumBreakdownVoltage; }
    void set_minimum_breakdown_voltage(const double & value) { this->minimumBreakdownVoltage = value; }
};


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

    const double & get_loss_density() const { return lossDensity; }
    double & get_mutable_loss_density() { return lossDensity; }
    void set_loss_density(const double & lossDensity) { this->lossDensity = lossDensity; }
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


class Wire : public MAS::Wire {
    public:
        Wire(MAS::Wire wire) {
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

        Wire(WireRound wire) {
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

        Wire() = default;
        virtual ~Wire() = default;

        static WireRound convert_from_wire_to_strand(Wire wire);
        static std::optional<InsulationWireCoating> resolve_coating(const Wire& wire);
        std::optional<InsulationWireCoating> resolve_coating();
        static WireRound resolve_strand(const Wire& wire);
        WireRound resolve_strand();
        static WireMaterial resolve_material(Wire wire);
        static WireMaterial resolve_material(WireRound wire);
        WireMaterial resolve_material();
        static InsulationMaterial resolve_coating_insulation_material(Wire wire);
        static InsulationMaterial resolve_coating_insulation_material(WireRound wire);
        static InsulationMaterial resolve_coating_insulation_material(InsulationWireCoating coating);
        InsulationMaterial resolve_coating_insulation_material();

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

        static double calculate_outer_height(Wire wire, OpenMagnetics::DimensionalValues preferredValue = OpenMagnetics::DimensionalValues::NOMINAL);
        double calculate_outer_height(OpenMagnetics::DimensionalValues preferredValue = OpenMagnetics::DimensionalValues::NOMINAL);
        static double calculate_outer_width(Wire wire, OpenMagnetics::DimensionalValues preferredValue = OpenMagnetics::DimensionalValues::NOMINAL);
        double calculate_outer_width(OpenMagnetics::DimensionalValues preferredValue = OpenMagnetics::DimensionalValues::NOMINAL);
        static double calculate_outer_diameter(Wire wire, OpenMagnetics::DimensionalValues preferredValue = OpenMagnetics::DimensionalValues::NOMINAL);
        double calculate_outer_diameter(OpenMagnetics::DimensionalValues preferredValue = OpenMagnetics::DimensionalValues::NOMINAL);

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

        void set_nominal_value_conducting_area(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_conducting_area(aux);
        }
        void set_bare_coating();

        int get_equivalent_insulation_layers(double voltageToInsulate);

        std::vector<std::vector<double>> calculate_current_density_distribution(SignalDescriptor current, double frequency, double temperature, size_t numberPoints=10);
        double calculate_dc_current_density(OperatingPointExcitation excitation);
        double calculate_dc_current_density(SignalDescriptor current);
        double calculate_effective_current_density(OperatingPointExcitation excitation, double temperature);
        double calculate_effective_current_density(SignalDescriptor current, double temperature);
        double calculate_effective_current_density(double rms, double frequency, double temperature);
        double calculate_effective_conducting_area(double frequency, double temperature);
        double calculate_conducting_area();
        double calculate_outer_area();

        static int calculate_number_parallels_needed(Inputs inputs, Wire& wire, double maximumEffectiveCurrentDensity, size_t windingIndex = 0);
        static int calculate_number_parallels_needed(OperatingPointExcitation excitation, double temperature, Wire& wire, double maximumEffectiveCurrentDensity);
        static int calculate_number_parallels_needed(SignalDescriptor current, double temperature, Wire& wire, double maximumEffectiveCurrentDensity);
        static int calculate_number_parallels_needed(double rms, double effectiveFrequency, double temperature, Wire& wire, double maximumEffectiveCurrentDensity);

        double get_maximum_outer_width();
        double get_maximum_outer_height();
        double get_maximum_outer_dimension();
        double get_maximum_conducting_width();
        double get_maximum_conducting_height();
        double get_minimum_conducting_dimension();
        double get_minimum_outer_dimension();


        double get_coating_thickness();
        double get_coating_dielectric_strength();
        double get_coating_thermal_conductivity();  // NEW: Thermal conductivity for thermal modeling
        static double get_coating_thickness(Wire wire);
        static double get_coating_dielectric_strength(Wire wire);
        static double get_coating_thermal_conductivity(Wire wire);  // NEW: Static version
        // NEW: Overloads that work directly with InsulationWireCoating (for thermal modeling)
        static double get_coating_thickness(const InsulationWireCoating& coating);
        static double get_coating_thermal_conductivity(const InsulationWireCoating& coating);
        std::string encode_coating_label();
        static std::string encode_coating_label(Wire wire);
        static std::optional<InsulationWireCoating> decode_coating_label(std::string label);
        double get_coating_relative_permittivity();
        static double get_coating_relative_permittivity(Wire wire);

        void cut_foil_wire_to_section(Section section);
        void cut_planar_wire_to_section(Section section);
        double get_relative_cost();
        std::string get_reference();

        static Wire get_equivalent_wire(Wire oldWire, WireType newWireType, double effectiveFrequency=200000, double temperature=Defaults().ambientTemperature);
        static Wire get_wire_for_frequency(double effectiveFrequency, double temperature=Defaults().ambientTemperature, bool exact=false);
        static Wire get_wire_for_conducting_area(double conductingArea, double temperature=Defaults().ambientTemperature, bool exact=false);
        static Wire get_wire_for_dc_resistance_per_meter(double dcResistancePerMeter, double temperature=Defaults().ambientTemperature);
        
        static Wire create_quick_litz_wire(double conductingDiameter, int64_t numberStrands);
        static Wire create_quick_rectangular_wire(double conductingWidth, double conductingHeight);

};

bool operator==(Wire lhs, Wire rhs);

inline bool operator==(Wire lhs, Wire rhs) {
    if (lhs.get_type() != rhs.get_type()) {
        return false;
    }
    if (lhs.get_name() && rhs.get_name()) {
        return lhs.get_name().value() == rhs.get_name().value();
    }
    if (lhs.get_conducting_diameter() && rhs.get_conducting_diameter()) {
        if (resolve_dimensional_values(lhs.get_conducting_diameter().value()) != resolve_dimensional_values(rhs.get_conducting_diameter().value())) {
            return false;
        }
    }
    if (lhs.get_conducting_height() && rhs.get_conducting_height()) {
        if (resolve_dimensional_values(lhs.get_conducting_height().value()) != resolve_dimensional_values(rhs.get_conducting_height().value())) {
            return false;
        }
    }
    if (lhs.get_conducting_width() && rhs.get_conducting_width()) {
        if (resolve_dimensional_values(lhs.get_conducting_width().value()) != resolve_dimensional_values(rhs.get_conducting_width().value())) {
            return false;
        }
    }
    if (lhs.get_outer_diameter() && rhs.get_outer_diameter()) {
        if (resolve_dimensional_values(lhs.get_outer_diameter().value()) != resolve_dimensional_values(rhs.get_outer_diameter().value())) {
            return false;
        }
    }
    if (lhs.get_outer_height() && rhs.get_outer_height()) {
        if (resolve_dimensional_values(lhs.get_outer_height().value()) != resolve_dimensional_values(rhs.get_outer_height().value())) {
            return false;
        }
    }
    if (lhs.get_outer_width() && rhs.get_outer_width()) {
        if (resolve_dimensional_values(lhs.get_outer_width().value()) != resolve_dimensional_values(rhs.get_outer_width().value())) {
            return false;
        }
    }
    if (lhs.get_number_conductors() && rhs.get_number_conductors()) {
        if (lhs.get_number_conductors().value() != rhs.get_number_conductors().value()) {
            return false;
        }
    }
    if (lhs.get_standard() && rhs.get_standard()) {
        if (lhs.get_standard().value() != rhs.get_standard().value()) {
            return false;
        }
    }
    if (lhs.get_standard_name() && rhs.get_standard_name()) {
        if (lhs.get_standard_name().value() != rhs.get_standard_name().value()) {
            return false;
        }
    }

    return true;
}
} // namespace OpenMagnetics