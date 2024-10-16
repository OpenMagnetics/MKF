#pragma once

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


class InsulationMaterialWrapper : public InsulationMaterial {
private:
    std::vector<std::pair<double, double>> _available_thicknesses;
public:
    InsulationMaterialWrapper(const json& j) {
        from_json(j, *this);
        extract_available_thicknesses();
    }
    InsulationMaterialWrapper() = default;
    virtual ~InsulationMaterialWrapper() = default;

    InsulationMaterialWrapper(InsulationMaterial insulationMaterial) {
        if (insulationMaterial.get_aliases()) {
            set_aliases(insulationMaterial.get_aliases().value());
        }
        if (insulationMaterial.get_composition()) {
            set_composition(insulationMaterial.get_composition().value());
        }
        if (insulationMaterial.get_dielectric_constant()) {
            set_dielectric_constant(insulationMaterial.get_dielectric_constant().value());
        }
        set_dielectric_strength(insulationMaterial.get_dielectric_strength());

        if (insulationMaterial.get_manufacturer()) {
            set_manufacturer(insulationMaterial.get_manufacturer().value());
        }
        if (insulationMaterial.get_melting_point()) {
            set_melting_point(insulationMaterial.get_melting_point().value());
        }
        set_name(insulationMaterial.get_name());

        if (insulationMaterial.get_resistivity()) {
            set_resistivity(insulationMaterial.get_resistivity().value());
        }
        if (insulationMaterial.get_specific_heat()) {
            set_specific_heat(insulationMaterial.get_specific_heat().value());
        }
        if (insulationMaterial.get_temperature_class()) {
            set_temperature_class(insulationMaterial.get_temperature_class().value());
        }
        if (insulationMaterial.get_thermal_conductivity()) {
            set_thermal_conductivity(insulationMaterial.get_thermal_conductivity().value());
        }
    }

    void extract_available_thicknesses();
    std::vector<std::pair<double, double>> get_available_thicknesses();
    std::pair<double, double> get_thicker_tape();
    std::pair<double, double> get_thinner_tape();
    double get_thicker_tape_thickness();
    double get_thinner_tape_thickness();

    static double get_dielectric_strength_by_thickness(InsulationMaterial materialData, double thickness);
    double get_dielectric_strength_by_thickness(double thickness);

};
} // namespace OpenMagnetics