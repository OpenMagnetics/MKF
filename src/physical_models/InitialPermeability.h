#pragma once
#include "Constants.h"

#include "constructive_models/CoreWrapper.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

class InitialPermeability {
    private:
    protected:
    public:
        static double get_initial_permeability(std::string coreMaterialName,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt,
                                        std::optional<double> magneticFluxDensity = std::nullopt);

        static double get_initial_permeability(CoreMaterial coreMaterial,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt,
                                        std::optional<double> magneticFluxDensity = std::nullopt);
        static double has_temperature_dependency(CoreMaterial coreMaterial);
        static double has_frequency_dependency(CoreMaterial coreMaterial);
        static double has_magnetic_field_dc_bias_dependency(CoreMaterial coreMaterial);
        static double get_initial_permeability_temperature_dependent(CoreMaterial coreMaterial, double temperature);
        static double get_initial_permeability_frequency_dependent(CoreMaterial coreMaterial, double frequency);
        static double get_initial_permeability_magnetic_field_dc_bias_dependent(CoreMaterial coreMaterial, double magneticFieldDcBias);
        static double calculate_frequency_for_initial_permeability_drop(CoreMaterial coreMaterial, double percentageDrop, double maximumError = 0.01);
        static std::vector<size_t> get_only_temperature_dependent_indexes(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_temperature_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_temperature_dependent_points(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_frequency_dependent_indexes(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_frequency_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_frequency_dependent_points(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_magnetic_field_dc_bias_dependent_indexes(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_magnetic_field_dc_bias_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_magnetic_field_dc_bias_dependent_points(CoreMaterial coreMaterial);
        static double get_initial_permeability_formula(CoreMaterial coreMaterial,
                                                       std::optional<double> temperature,
                                                       std::optional<double> magneticFieldDcBias,
                                                       std::optional<double> frequency,
                                                       std::optional<double> magneticFluxDensity);
};

} // namespace OpenMagnetics