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

double InitialPermeability::get_initial_permeability(std::string coreMaterialName,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency) {
    CoreMaterial coreMaterial = CoreWrapper::resolve_material(coreMaterialName);
    return get_initial_permeability(coreMaterial,temperature, magneticFieldDcBias, frequency);
}

double InitialPermeability::get_initial_permeability(CoreMaterial coreMaterial,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();

    double initialPermeabilityValue = 1;

    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        auto permeabilityPoint = std::get<PermeabilityPoint>(initialPermeabilityData);
        initialPermeabilityValue = permeabilityPoint.get_value();

        if (permeabilityPoint.get_modifiers()) {
            InitialPermeabilitModifier modifiers = (*permeabilityPoint.get_modifiers())["default"];
            if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MAGNETICS) {
                auto temperatureFactor = modifiers.get_temperature_factor();
                if (temperature && temperatureFactor) {
                    double permeabilityVariationDueToTemperature =
                        temperatureFactor->get_a() + temperatureFactor->get_b().value() * (temperature.value()) +
                        temperatureFactor->get_c().value() * pow(temperature.value(), 2) +
                        temperatureFactor->get_d().value() * pow(temperature.value(), 3) +
                        temperatureFactor->get_e().value() * pow(temperature.value(), 4);
                    initialPermeabilityValue *= (1 + permeabilityVariationDueToTemperature);
                }

                auto frequencyFactor = modifiers.get_frequency_factor();
                if (frequency && frequencyFactor) {
                    double permeabilityVariationDueToFrequency =
                        frequencyFactor->get_a() + frequencyFactor->get_b() * (frequency.value()) +
                        frequencyFactor->get_c() * pow(frequency.value(), 2) + frequencyFactor->get_d() * pow(frequency.value(), 3) +
                        frequencyFactor->get_e().value() * pow(frequency.value(), 4);

                    initialPermeabilityValue *= (1 + permeabilityVariationDueToFrequency);
                }

                if (magneticFieldDcBias) {
                    auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
                    double permeabilityVariationDueToMagneticFieldDcBias =
                        0.01 / (magneticFieldDcBiasFactor.get_a() +
                                magneticFieldDcBiasFactor.get_b() *
                                    pow(roundFloat(fabs(magneticFieldDcBias.value()), 3), magneticFieldDcBiasFactor.get_c()));

                    initialPermeabilityValue *= permeabilityVariationDueToMagneticFieldDcBias;
                }
            }
        }
    }
    else {
        auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
        bool hasTemperatureRequirement = temperature.has_value();
        bool hasTemperatureDependency = permeabilityPoints[0].get_temperature().has_value();

        if (hasTemperatureDependency) {
            int n = permeabilityPoints.size();
            std::vector<double> x, y;

            double temperaturePoint;
            if (!hasTemperatureRequirement)
                temperaturePoint = 25;
            else
                temperaturePoint = temperature.value();

            for (int i = 0; i < n; i++) {
                if (x.size() == 0 || (*permeabilityPoints[i].get_temperature()) != x.back()) {
                    x.push_back(*permeabilityPoints[i].get_temperature());
                    y.push_back(permeabilityPoints[i].get_value());
                }
            }

            tk::spline interp(x, y, tk::spline::cspline_hermite);
            initialPermeabilityValue = std::max(1., interp(temperaturePoint));
        }
        else {
            throw std::invalid_argument("Invalid material permeability");
        }
    }

    if (coreMaterial.get_curie_temperature() && temperature) {
        if ((temperature.value()) > (*coreMaterial.get_curie_temperature())) {
            initialPermeabilityValue = 1;
        }
    }

    if (std::isnan(initialPermeabilityValue)) {
        throw std::runtime_error("Initial Permeability must be a number, not NaN");
    }

    return initialPermeabilityValue;
}

} // namespace OpenMagnetics
