#include "InitialPermeability.h"

#include "Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "spline.h"
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

double InitialPermeability::get_initial_permeability(CoreMaterialDataOrNameUnion material,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency) {
    CoreMaterial material_data;
    // If the material is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(material)) {
        material_data = OpenMagnetics::find_core_material_by_name(std::get<std::string>(material));
    }
    else {
        material_data = std::get<OpenMagnetics::CoreMaterial>(material);
    }

    auto initial_permeability_data = material_data.get_permeability().get_initial();

    double initial_permeability_value = 1;

    if (std::holds_alternative<PermeabilityPoint>(initial_permeability_data)) {
        auto permeability_point = std::get<PermeabilityPoint>(initial_permeability_data);
        initial_permeability_value = permeability_point.get_value();

        if (permeability_point.get_modifiers()) {
            InitialPermeabilitModifier modifiers = (*permeability_point.get_modifiers())["default"];
            if ((*modifiers.get_method()) == "magnetics") {
                auto temperature_factor = modifiers.get_temperature_factor();
                if (temperature && temperature_factor) {
                    double permeability_variation_due_to_temperature =
                        temperature_factor->get_a() + temperature_factor->get_b().value() * (temperature.value()) +
                        temperature_factor->get_c().value() * pow(temperature.value(), 2) +
                        temperature_factor->get_d().value() * pow(temperature.value(), 3) +
                        temperature_factor->get_e().value() * pow(temperature.value(), 4);
                    initial_permeability_value *= (1 + permeability_variation_due_to_temperature);
                }

                auto frequency_factor = modifiers.get_frequency_factor();
                if (frequency && frequency_factor) {
                    double permeability_variation_due_to_frequency =
                        frequency_factor->get_a() + frequency_factor->get_b() * (frequency.value()) +
                        frequency_factor->get_c() * pow(frequency.value(), 2) + frequency_factor->get_d() * pow(frequency.value(), 3) +
                        frequency_factor->get_e().value() * pow(frequency.value(), 4);

                    initial_permeability_value *= (1 + permeability_variation_due_to_frequency);
                }

                if (magneticFieldDcBias) {
                    auto magnetic_field_dc_bias_factor = modifiers.get_magnetic_field_dc_bias_factor();
                    double permeability_variation_due_to_magnetic_field_dc_bias =
                        0.01 / (magnetic_field_dc_bias_factor.get_a() +
                                magnetic_field_dc_bias_factor.get_b() *
                                    pow(roundFloat(fabs(magneticFieldDcBias.value()), 3), magnetic_field_dc_bias_factor.get_c()));

                    initial_permeability_value *= permeability_variation_due_to_magnetic_field_dc_bias;
                }
            }
        }
    }
    else {
        auto permeability_points = std::get<std::vector<PermeabilityPoint>>(initial_permeability_data);
        bool has_temperature_requirement = temperature.has_value();
        bool has_temperature_dependency = permeability_points[0].get_temperature().has_value();

        if (has_temperature_dependency) {
            int n = permeability_points.size();
            std::vector<double> x, y;

            double temperature_point;
            if (!has_temperature_requirement)
                temperature_point = 25;
            else
                temperature_point = temperature.value();

            for (int i = 0; i < n; i++) {
                if (x.size() == 0 || (*permeability_points[i].get_temperature()) != x.back()) {
                    x.push_back(*permeability_points[i].get_temperature());
                    y.push_back(permeability_points[i].get_value());
                }
            }

            tk::spline interp(x, y, tk::spline::cspline_hermite);
            initial_permeability_value = std::max(1., interp(temperature_point));
        }
        else {
            throw std::invalid_argument("Invalid material permeability");
        }
    }

    if (material_data.get_curie_temperature() && temperature) {
        if ((temperature.value()) > (*material_data.get_curie_temperature())) {
            initial_permeability_value = 1;
        }
    }

    if (std::isnan(initial_permeability_value)) {
        throw std::runtime_error("Initial Permeability must be a number, not NaN");
    }

    return initial_permeability_value;
}

} // namespace OpenMagnetics
