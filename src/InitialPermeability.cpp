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
                                                     std::optional<double> frequency,
                                                     std::optional<double> magneticFluxDensity) {
    CoreMaterial coreMaterial = CoreWrapper::resolve_material(coreMaterialName);
    return get_initial_permeability(coreMaterial,temperature, magneticFieldDcBias, frequency, magneticFluxDensity);
}

double InitialPermeability::get_initial_permeability(CoreMaterial coreMaterial,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency,
                                                     std::optional<double> magneticFluxDensity) {
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
                    double a = magneticFieldDcBiasFactor->get_a();
                    double b = magneticFieldDcBiasFactor->get_b();
                    double c = magneticFieldDcBiasFactor->get_c();
                    double H = magneticFieldDcBias.value();
                    double permeabilityVariationDueToMagneticFieldDcBias = 0.01 / (a + b * pow(roundFloat(fabs(H), 3), c));

                    initialPermeabilityValue *= permeabilityVariationDueToMagneticFieldDcBias;
                }
            }
            else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MICROMETALS) {

                auto frequencyFactor = modifiers.get_frequency_factor();
                if (frequency && frequencyFactor) {
                    double a = frequencyFactor->get_a();
                    double b = frequencyFactor->get_b();
                    double c = frequencyFactor->get_c();
                    double d = frequencyFactor->get_d();
                    double f = frequency.value();
                    initialPermeabilityValue = 1.0 / (a + b * pow(f, c)) + d;
                }

                auto temperatureFactor = modifiers.get_temperature_factor();
                if (temperature && temperatureFactor) {
                    double permeabilityVariationDueToTemperature;
                    if (temperatureFactor->get_b()) {
                        double a = temperatureFactor->get_a();
                        double b = temperatureFactor->get_b().value();
                        double c = temperatureFactor->get_c().value();
                        double d = temperatureFactor->get_d().value();
                        double e = temperatureFactor->get_e().value();
                        double T = temperature.value();
                        permeabilityVariationDueToTemperature = (a + c * T + e * pow(T, 2)) / (1 + b * T + d * pow(T, 2));
                    }
                    else {
                        double a = temperatureFactor->get_a();
                        double T = temperature.value();
                        permeabilityVariationDueToTemperature = a * (T - 20) * 0.0001;
                    }

                    initialPermeabilityValue *= (1 + permeabilityVariationDueToTemperature * 0.01);
                }

                if (magneticFieldDcBias) {
                    auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
                    double a = magneticFieldDcBiasFactor->get_a();
                    double b = magneticFieldDcBiasFactor->get_b();
                    double c = magneticFieldDcBiasFactor->get_c();
                    double d = magneticFieldDcBiasFactor->get_d().value();
                    double H = magneticFieldDcBias.value();
                    double permeabilityVariationDueToMagneticFieldDcBias = 1.0 / (a + b * pow(roundFloat(fabs(H), 3), c)) + d;

                    initialPermeabilityValue *= permeabilityVariationDueToMagneticFieldDcBias * 0.01;
                }

                if (magneticFluxDensity) {
                    auto magneticFluxDensityFactor = modifiers.get_magnetic_flux_density_factor();
                    double a = magneticFluxDensityFactor->get_a();
                    double b = magneticFluxDensityFactor->get_b();
                    double c = magneticFluxDensityFactor->get_c();
                    double d = magneticFluxDensityFactor->get_d();
                    double e = magneticFluxDensityFactor->get_e();
                    double f = magneticFluxDensityFactor->get_f();
                    double B = magneticFluxDensity.value();
                    double permeabilityVariationDueToMagneticFluxDensity = 1.0 / (1.0 / (a + b * pow(B, c)) + 1.0 / (d * pow(B, e)) + 1.0 / f);

                    initialPermeabilityValue *= permeabilityVariationDueToMagneticFluxDensity * 0.01;
                }
            }
            else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::FAIR_RITE) {
                auto temperatureFactor = modifiers.get_temperature_factor();
                if (temperature && temperatureFactor) {
                    double permeabilityVariationDueToTemperature;
                    double a = temperatureFactor->get_a();
                    double T = temperature.value();
                    permeabilityVariationDueToTemperature = a * T;

                    initialPermeabilityValue *= (1 + permeabilityVariationDueToTemperature * 0.01);
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
            std::vector<double> x, y, distinctTemperatureValues;

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

            if (x.size() > 1) {
                tk::spline interp(x, y, tk::spline::cspline_hermite);
                initialPermeabilityValue = std::max(1., interp(temperaturePoint));
            }
            else {
                initialPermeabilityValue = permeabilityPoints[0].get_value();
            }
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
