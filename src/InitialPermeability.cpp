#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <streambuf>
#include <magic_enum.hpp>
#include "InitialPermeability.h"
#include "Utils.h"
#include <libInterpolate/Interpolate.hpp>


namespace OpenMagnetics {

    double InitialPermeability::get_initial_permeability(CoreMaterialDataOrNameUnion material,
                                                         double* temperature,
                                                         double* magneticFieldDcBias,
                                                         double* frequency) {
        CoreMaterial material_data;
        // If the material is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(material)) {
            material_data = OpenMagnetics::find_data_by_name<OpenMagnetics::CoreMaterial>(std::get<std::string>(material));
        }
        else {
            material_data = std::get<OpenMagnetics::CoreMaterial>(material);
        }

        auto initial_permeability_data = material_data.get_permeability().get_initial();

        double initial_permeability_value;

        if (std::holds_alternative<PermeabilityPoint>(initial_permeability_data)) {
            auto permeability_point = std::get<PermeabilityPoint>(initial_permeability_data);
            initial_permeability_value = permeability_point.get_value();

            if (permeability_point.get_modifiers()) {
                MagneticsPermeabilityMethodData modifiers = (*permeability_point.get_modifiers())["default"];
                if ((*modifiers.get_method()) == "Magnetics") {
                    if (temperature) {
                        auto temperature_factor = (*modifiers.get_temperature_factor());
                        double permeability_variation_due_to_temperature = temperature_factor.get_a() + temperature_factor.get_b() * (*temperature) + temperature_factor.get_c() * pow(*temperature, 2) + temperature_factor.get_d() * pow(*temperature, 3) + temperature_factor.get_e() * pow(*temperature, 4);
                        initial_permeability_value *= (1 + permeability_variation_due_to_temperature);
                    }

                    if (frequency) {
                        auto frequency_factor = (*modifiers.get_frequency_factor());
                        double permeability_variation_due_to_frequency = frequency_factor.get_a() + frequency_factor.get_b() * (*frequency) + frequency_factor.get_c() * pow(*frequency, 2) + frequency_factor.get_d() * pow(*frequency, 3) + frequency_factor.get_e() * pow(*frequency, 4);

                        initial_permeability_value *= (1 + permeability_variation_due_to_frequency);
                    }

                    if (magneticFieldDcBias) {
                        auto magnetic_field_dc_bias_factor = modifiers.get_magnetic_field_dc_bias_factor();
                        double permeability_variation_due_to_magnetic_field_dc_bias = 0.01 / (magnetic_field_dc_bias_factor.get_a() + magnetic_field_dc_bias_factor.get_b() * pow(roundFloat(*magneticFieldDcBias, 3), magnetic_field_dc_bias_factor.get_c()));

                        initial_permeability_value *= permeability_variation_due_to_magnetic_field_dc_bias;
                    }
                }
            }
        }
        else {
            auto permeability_points = std::get<std::vector<PermeabilityPoint>>(initial_permeability_data);
            bool has_temperature_requirement = temperature;
            bool has_temperature_dependency = permeability_points[0].get_temperature().has_value();

            if (has_temperature_dependency) {
                int n = permeability_points.size();
                _1D::CubicSplineInterpolator<double> interp;
                _1D::CubicSplineInterpolator<double>::VectorType xx(n), yy(n);

                double temperature_point;
                if (!has_temperature_requirement)
                    temperature_point = 25;
                else
                    temperature_point = *temperature;

                for(int i = 0; i < n; i++) {
                    xx(i) = *permeability_points[i].get_temperature();
                    yy(i) = permeability_points[i].get_value();
                }

                interp.setData( xx, yy );
                initial_permeability_value = std::max(1., interp(temperature_point));
            }
            else {
                throw std::invalid_argument( "Invalid material permeability" );
            }
        }

        if (material_data.get_curie_temperature() && temperature) {
            if ((*temperature) > (*material_data.get_curie_temperature())) {
                initial_permeability_value = 1;
            }
        }

        if (std::isnan(initial_permeability_value))  {
            throw std::runtime_error("Initial Permeability must be a number, not NaN");
        }

        return initial_permeability_value;
    };

}
