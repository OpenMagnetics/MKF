#include "physical_models/InitialPermeability.h"
#include "physical_models/AmplitudePermeability.h"

#include "support/Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "spline.h"
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include <cfloat>


std::map<std::string, std::variant<double, tk::spline>> initialPermeabilityMagneticFieldDcBiasInterps;
std::map<std::string, std::variant<double, tk::spline>> initialPermeabilityFrequencyInterps;
std::map<std::string, std::variant<double, tk::spline>> initialPermeabilityTemperatureInterps;

namespace OpenMagnetics {

double InitialPermeability::get_initial_permeability(std::string coreMaterialName,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency,
                                                     std::optional<double> magneticFluxDensity) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_initial_permeability(coreMaterial, temperature, magneticFieldDcBias, frequency, magneticFluxDensity);
}


double InitialPermeability::get_initial_permeability(CoreMaterial coreMaterial, OperatingPoint operatingPoint) {
    if (operatingPoint.get_excitations_per_winding().size() == 0) {
        throw std::runtime_error("Operating point is missing excitations");
    }
    double temperature = operatingPoint.get_conditions().get_ambient_temperature();
    auto frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();

    std::optional<double> magneticFieldDcBias = std::nullopt;
    if (operatingPoint.get_excitations_per_winding()[0].get_magnetic_field_strength()) {
        if (operatingPoint.get_excitations_per_winding()[0].get_magnetic_field_strength()->get_processed()) {
            magneticFieldDcBias = operatingPoint.get_excitations_per_winding()[0].get_magnetic_field_strength()->get_processed()->get_offset();
        }
    }
    std::optional<double> magneticFluxDensity = std::nullopt;
    if (operatingPoint.get_excitations_per_winding()[0].get_magnetic_flux_density()) {
        if (operatingPoint.get_excitations_per_winding()[0].get_magnetic_flux_density()->get_processed()) {
            magneticFluxDensity = operatingPoint.get_excitations_per_winding()[0].get_magnetic_flux_density()->get_processed()->get_peak();
        }
    }
    return get_initial_permeability(coreMaterial, temperature, magneticFieldDcBias, frequency, magneticFluxDensity);
}

double InitialPermeability::get_initial_permeability(std::string coreMaterialName, OperatingPoint operatingPoint) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_initial_permeability(coreMaterial, operatingPoint);
}

double InitialPermeability::has_frequency_dependency(CoreMaterial coreMaterial) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        return false;
    }
    else {
        auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
        for (auto point : permeabilityPoints) {
            if (point.get_frequency()) {
                if (point.get_frequency().value() > Constants().quasiStaticFrequencyLimit) {
                    return true;
                }
            }
        }
        return false;
    }
}

double InitialPermeability::has_temperature_dependency(CoreMaterial coreMaterial) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    double firstTemperatureValue = DBL_MAX;
    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        return false;
    }
    else {
        auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
        for (auto point : permeabilityPoints) {
            if (point.get_temperature()) {
                if (firstTemperatureValue == DBL_MAX) {
                    firstTemperatureValue = point.get_temperature().value();
                }
                else if (firstTemperatureValue != point.get_temperature().value()) {
                    return true;
                }
            }
        }
        return false;
    }
}

double InitialPermeability::has_magnetic_field_dc_bias_dependency(CoreMaterial coreMaterial) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        return false;
    }
    else {
        auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
        for (auto point : permeabilityPoints) {
            if (point.get_magnetic_field_dc_bias()) {
                if (point.get_magnetic_field_dc_bias().value() > 0) {
                    return true;
                }
            }
        }
        return false;
    }
}

std::map<std::string, std::string> InitialPermeability::get_initial_permeability_equations(PermeabilityPoint permeabilityPoint) {
    std::map<std::string, std::string> equations;

    if (permeabilityPoint.get_modifiers()) {
        InitialPermeabilitModifier modifiers = (*permeabilityPoint.get_modifiers())["default"];
        if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MAGNETICS) {
            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperatureFactor) {
                equations["temperatureFactor"] = "1 + (a + b*T + c*T^2 + d*T^3 + e*T^4)";
            }

            auto frequencyFactor = modifiers.get_frequency_factor();
            if (frequencyFactor) {
                equations["frequencyFactor"] = "1 + (a + b*f + c*f^2 + d*f^3 + e*f^4)";
            }

            auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
            if (magneticFieldDcBiasFactor) {
                equations["magneticFieldDcBiasFactor"] = "0.01 / (a + b*(H^c))";
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MICROMETALS) {

            auto frequencyFactor = modifiers.get_frequency_factor();
            if (frequencyFactor) {
                equations["frequencyFactor"] = "(1.0 / (a + b * f^c) + d) / mu_ini";
            }

            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperatureFactor) {
                if (temperatureFactor->get_b()) {
                    equations["temperatureFactor"] = "1 + (a + c * T + e * T^2) / (1 + b * T + d * T^2) * 0.01";
                }
                else {
                    equations["temperatureFactor"] = "1 + (a * (T - 20) * 0.0001) * 0.01";
                }
            }

            auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
            if (magneticFieldDcBiasFactor) {
                equations["magneticFieldDcBiasFactor"] = "(1.0 / (a + b * H^c) + d) * 0.01";
            }

            auto magneticFluxDensityFactor = modifiers.get_magnetic_flux_density_factor();
            if (magneticFluxDensityFactor) {
                equations["magneticFluxDensityFactor"] = "(1.0 / (1.0 / (a + b * B^c) + 1.0 / (d * B^e) + 1.0 / f)) * 0.01";
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::FAIR_RITE) {
            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperatureFactor) {
                equations["temperatureFactor"] = "(1 + a * T * 0.01)";
            }

        }
    }
    return equations;

}

double InitialPermeability::get_initial_permeability_formula(CoreMaterial coreMaterial,
                                                             std::optional<double> temperature,
                                                             std::optional<double> magneticFieldDcBias,
                                                             std::optional<double> frequency,
                                                             std::optional<double> magneticFluxDensity) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    auto permeabilityPoint = std::get<PermeabilityPoint>(initialPermeabilityData);
    double initialPermeabilityValue = permeabilityPoint.get_value();

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
    return initialPermeabilityValue;
}

double get_closes_temperature_to_default_in_permeability_points(std::vector<PermeabilityPoint> points) {
    double closestDistance = DBL_MAX;
    double closestTemperature = DBL_MAX;
    for (auto point : points) {
        if (point.get_temperature()) {
            double distance = fabs(point.get_temperature().value() - Defaults().ambientTemperature);
            if (distance < closestDistance) {
                closestDistance = distance;
                closestTemperature = point.get_temperature().value();
            }
        }
    }

    return closestTemperature;
}

double get_minimum_frequency_in_permeability_points(std::vector<PermeabilityPoint> points) {
    auto order_frequency = []( const auto &e1, const auto &e2 ) {
        if (!e1.get_frequency()) {
            return false;
        }
        if (!e2.get_frequency()) {
            return true;
        }
        return e1.get_frequency().value() < e2.get_frequency().value();
    };

    auto minimumFrequencyPoint = std::min_element( std::begin( points ), std::end( points ), order_frequency );
    if (minimumFrequencyPoint->get_frequency()) {
        return minimumFrequencyPoint->get_frequency().value();
    }
    else {
        return Constants().quasiStaticFrequencyLimit;
    }
}

double get_maximum_frequency_in_permeability_points(std::vector<PermeabilityPoint> points) {
    auto order_frequency = []( const auto &e1, const auto &e2 ) {
        if (!e1.get_frequency()) {
            return false;
        }
        if (!e2.get_frequency()) {
            return true;
        }
        return e1.get_frequency().value() < e2.get_frequency().value();
    };

    auto maximumFrequencyPoint = std::max_element( std::begin( points ), std::end( points ), order_frequency );
    if (maximumFrequencyPoint->get_frequency()) {
        return maximumFrequencyPoint->get_frequency().value();
    }
    else {
        return DBL_MAX;
    }
}

double get_minimum_magnetic_field_dc_bias_in_permeability_points(std::vector<PermeabilityPoint> points) {
    auto minimum_magnetic_field_dc_bias = []( const auto &e1, const auto &e2 ) {
        if (!e1.get_magnetic_field_dc_bias()) {
            return false;
        }
        if (!e2.get_magnetic_field_dc_bias()) {
            return true;
        }
        return e1.get_magnetic_field_dc_bias().value() < e2.get_magnetic_field_dc_bias().value();
    };

    auto minimumMagneticFieldDcBiasPoint = std::min_element( std::begin( points ), std::end( points ), minimum_magnetic_field_dc_bias );
    if (minimumMagneticFieldDcBiasPoint->get_magnetic_field_dc_bias()) {
        return minimumMagneticFieldDcBiasPoint->get_magnetic_field_dc_bias().value();
    }
    else {
        return 0;
    }
}

std::vector<size_t> InitialPermeability::get_only_temperature_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints) {
    std::vector<size_t> temperatureIndexes;
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
    for (size_t index = 0; index < permeabilityPoints.size(); ++index) {
        auto point = permeabilityPoints[index];
        if (point.get_frequency()) {
            if (point.get_frequency().value() > minimumFrequency) {
                continue;
            }
        }
        if (point.get_magnetic_field_dc_bias()) {
            if (point.get_magnetic_field_dc_bias().value() > minimumMagneticFieldDcBias) {
                continue;
            }
        }
        temperatureIndexes.push_back(index);
    }
    return temperatureIndexes;
}

std::vector<size_t> InitialPermeability::get_only_temperature_dependent_indexes(CoreMaterial coreMaterial) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    return get_only_temperature_dependent_indexes(permeabilityPoints);
}

std::vector<PermeabilityPoint> InitialPermeability::get_only_temperature_dependent_points(CoreMaterial coreMaterial) {
    std::vector<PermeabilityPoint> temperaturePoints;
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
    for (auto point : permeabilityPoints) {
        if (point.get_frequency()) {
            if (point.get_frequency().value() > minimumFrequency) {
                continue;
            }
        }
        if (point.get_magnetic_field_dc_bias()) {
            if (point.get_magnetic_field_dc_bias().value() > minimumMagneticFieldDcBias) {
                continue;
            }
        }
        temperaturePoints.push_back(point);
    }
    return temperaturePoints;
}

std::vector<size_t> InitialPermeability::get_only_frequency_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints) {
    std::vector<size_t> frequencyIndexes;
    double defaultTemperature = get_closes_temperature_to_default_in_permeability_points(permeabilityPoints);
    double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
    for (size_t index = 0; index < permeabilityPoints.size(); ++index) {
        auto point = permeabilityPoints[index];
        if (point.get_temperature()) {
            if (point.get_temperature().value() != defaultTemperature) {
                continue;
            }
        }
        if (point.get_magnetic_field_dc_bias()) {
            if (point.get_magnetic_field_dc_bias().value() > minimumMagneticFieldDcBias) {
                continue;
            }
        }
        frequencyIndexes.push_back(index);
    }
    return frequencyIndexes;
}

std::vector<size_t> InitialPermeability::get_only_frequency_dependent_indexes(CoreMaterial coreMaterial) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    return get_only_frequency_dependent_indexes(permeabilityPoints);
}

std::vector<PermeabilityPoint> InitialPermeability::get_only_frequency_dependent_points(CoreMaterial coreMaterial) {
    std::vector<PermeabilityPoint> frequencyPoints;
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    double defaultTemperature = get_closes_temperature_to_default_in_permeability_points(permeabilityPoints);
    double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
    for (auto point : permeabilityPoints) {
        if (point.get_temperature()) {
            if (point.get_temperature().value() != defaultTemperature) {
                continue;
            }
        }
        if (point.get_magnetic_field_dc_bias()) {
            if (point.get_magnetic_field_dc_bias().value() > minimumMagneticFieldDcBias) {
                continue;
            }
        }
        frequencyPoints.push_back(point);
    }
    return frequencyPoints;
}

std::vector<size_t> InitialPermeability::get_only_magnetic_field_dc_bias_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints) {
    std::vector<size_t>  magneticFieldDcBiasIndexes;
    double defaultTemperature = get_closes_temperature_to_default_in_permeability_points(permeabilityPoints);
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    for (size_t index = 0; index < permeabilityPoints.size(); ++index) {
        auto point = permeabilityPoints[index];
        if (point.get_temperature()) {
            if (point.get_temperature().value() != defaultTemperature) {
                continue;
            }
        }
        if (point.get_frequency()) {
            if (point.get_frequency().value() > minimumFrequency) {
                continue;
            }
        }
         magneticFieldDcBiasIndexes.push_back(index);
    }
    return  magneticFieldDcBiasIndexes;
}


std::vector<size_t> InitialPermeability::get_only_magnetic_field_dc_bias_dependent_indexes(CoreMaterial coreMaterial) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    return get_only_magnetic_field_dc_bias_dependent_indexes(permeabilityPoints);
}

std::vector<PermeabilityPoint> InitialPermeability::get_only_magnetic_field_dc_bias_dependent_points(CoreMaterial coreMaterial) {
    std::vector<PermeabilityPoint>  magneticFieldDcBiasPoints;
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    double defaultTemperature = get_closes_temperature_to_default_in_permeability_points(permeabilityPoints);
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    for (auto point : permeabilityPoints) {
        if (point.get_temperature()) {
            if (point.get_temperature().value() != defaultTemperature) {
                continue;
            }
        }
        if (point.get_frequency()) {
            if (point.get_frequency().value() > minimumFrequency) {
                continue;
            }
        }
         magneticFieldDcBiasPoints.push_back(point);
    }
    return  magneticFieldDcBiasPoints;
}

double InitialPermeability::get_initial_permeability_temperature_dependent(CoreMaterial coreMaterial, double temperature) {
    double initialPermeabilityValue = 1;
    auto permeabilityPoints = get_only_temperature_dependent_points(coreMaterial);
    if (permeabilityPoints.size() == 0) {
        throw std::runtime_error("No temperature dependent points for material: " + coreMaterial.get_name());
    }

    if (!initialPermeabilityTemperatureInterps.contains(coreMaterial.get_name())) {
        int n = permeabilityPoints.size();
        std::vector<double> x, y;

        for (int i = 0; i < n; i++) {
            if (x.size() == 0 || (*permeabilityPoints[i].get_temperature()) != x.back()) {
                x.push_back(*permeabilityPoints[i].get_temperature());
                y.push_back(permeabilityPoints[i].get_value());
            }
        }

        if (x.size() > 1) {
            tk::spline interp(x, y, tk::spline::cspline_hermite);
            initialPermeabilityTemperatureInterps[coreMaterial.get_name()] = interp;
        }
        else {
            initialPermeabilityTemperatureInterps[coreMaterial.get_name()] = permeabilityPoints[0].get_value();
        }
    }

    if (std::holds_alternative<double>(initialPermeabilityTemperatureInterps[coreMaterial.get_name()])) {
        initialPermeabilityValue = std::get<double>(initialPermeabilityTemperatureInterps[coreMaterial.get_name()]);
    }
    else {
        auto value = std::get<tk::spline>(initialPermeabilityTemperatureInterps[coreMaterial.get_name()])(temperature);
        initialPermeabilityValue = std::max(1., value);
    }


    if (coreMaterial.get_curie_temperature()) {
        if (temperature > (*coreMaterial.get_curie_temperature())) {
            initialPermeabilityValue = 1;
        }
    }

    if (std::isnan(initialPermeabilityValue)) {
        throw std::runtime_error("Initial Permeability must be a number, not NaN");
    }

    return initialPermeabilityValue;
}

double InitialPermeability::get_initial_permeability_frequency_dependent(CoreMaterial coreMaterial, double frequency) {
    double initialPermeabilityValue = 1;
    auto permeabilityPoints = get_only_frequency_dependent_points(coreMaterial);
    if (permeabilityPoints.size() == 0) {
        throw std::runtime_error("No frequency dependent points for material: " + coreMaterial.get_name());
    }

    if (!initialPermeabilityFrequencyInterps.contains(coreMaterial.get_name())) {
        int n = permeabilityPoints.size();
        std::vector<double> x, y;

        for (int i = 0; i < n; i++) {
            if (x.size() == 0 || (*permeabilityPoints[i].get_frequency()) != x.back()) {
                x.push_back(*permeabilityPoints[i].get_frequency());
                y.push_back(permeabilityPoints[i].get_value());
            }
        }

        if (x.size() > 1) {
            tk::spline interp(x, y, tk::spline::cspline_hermite);
            initialPermeabilityFrequencyInterps[coreMaterial.get_name()] = interp;
        }
        else {
            initialPermeabilityFrequencyInterps[coreMaterial.get_name()] = permeabilityPoints[0].get_value();
        }
    }

    if (std::holds_alternative<double>(initialPermeabilityFrequencyInterps[coreMaterial.get_name()])) {
        initialPermeabilityValue = std::get<double>(initialPermeabilityFrequencyInterps[coreMaterial.get_name()]);
    }
    else {
        auto value = std::get<tk::spline>(initialPermeabilityFrequencyInterps[coreMaterial.get_name()])(frequency);
        initialPermeabilityValue = std::max(1., value);
    }


    if (std::isnan(initialPermeabilityValue)) {
        throw std::runtime_error("Initial Permeability must be a number, not NaN");
    }

    return initialPermeabilityValue;
}

double InitialPermeability::get_initial_permeability_magnetic_field_dc_bias_dependent(CoreMaterial coreMaterial, double magneticFieldDcBias) {
    double initialPermeabilityValue = 1;
    auto permeabilityPoints = get_only_magnetic_field_dc_bias_dependent_points(coreMaterial);
    if (permeabilityPoints.size() == 0) {
        throw std::runtime_error("No magnetic field dc bias dependent points for material: " + coreMaterial.get_name());
    }

    if (!initialPermeabilityMagneticFieldDcBiasInterps.contains(coreMaterial.get_name())) {
        int n = permeabilityPoints.size();
        std::vector<double> x, y;

        for (int i = 0; i < n; i++) {
            if (x.size() == 0 || (*permeabilityPoints[i].get_magnetic_field_dc_bias()) != x.back()) {
                x.push_back(*permeabilityPoints[i].get_magnetic_field_dc_bias());
                y.push_back(permeabilityPoints[i].get_value());
            }
        }

        if (x.size() > 1) {
            tk::spline interp(x, y, tk::spline::cspline_hermite);
            initialPermeabilityMagneticFieldDcBiasInterps[coreMaterial.get_name()] = interp;
        }
        else {
            initialPermeabilityMagneticFieldDcBiasInterps[coreMaterial.get_name()] = permeabilityPoints[0].get_value();
        }
    }

    if (std::holds_alternative<double>(initialPermeabilityMagneticFieldDcBiasInterps[coreMaterial.get_name()])) {
        initialPermeabilityValue = std::get<double>(initialPermeabilityMagneticFieldDcBiasInterps[coreMaterial.get_name()]);
    }
    else {
        auto value = std::get<tk::spline>(initialPermeabilityMagneticFieldDcBiasInterps[coreMaterial.get_name()])(magneticFieldDcBias);
        initialPermeabilityValue = std::max(1., value);
    }

    if (std::isnan(initialPermeabilityValue)) {
        throw std::runtime_error("Initial Permeability must be a number, not NaN");
    }

    return initialPermeabilityValue;
}

double InitialPermeability::calculate_frequency_for_initial_permeability_drop(CoreMaterial coreMaterial, double percentageDrop, double maximumError) {
    bool hasFrequencyDependency = has_frequency_dependency(coreMaterial);
    if (!hasFrequencyDependency) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    double maximumFrequency = get_maximum_frequency_in_permeability_points(permeabilityPoints);
    double initialPermeabilityValueReference = get_initial_permeability_frequency_dependent(coreMaterial, minimumFrequency);
    double initialPermeabilityAbsolute = initialPermeabilityValueReference * (1 - percentageDrop);
    double currentFrequency = (maximumFrequency + minimumFrequency) / 2;

    size_t timeout = 100;
    while (timeout > 0) {
        auto currentInitialPermeability = get_initial_permeability_frequency_dependent(coreMaterial, currentFrequency);
        if (fabs(currentInitialPermeability - initialPermeabilityAbsolute) / initialPermeabilityAbsolute < maximumError) {
            break;
        }
        else {
            if (currentInitialPermeability > initialPermeabilityAbsolute) {
                minimumFrequency = currentFrequency;
                currentFrequency = (maximumFrequency + currentFrequency) / 2;
            }
            else {
                maximumFrequency = currentFrequency;
                currentFrequency = (currentFrequency + minimumFrequency) / 2;
            }
        }

        timeout--;
    }

    return currentFrequency;
}


double InitialPermeability::get_initial_permeability(CoreMaterial coreMaterial,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency,
                                                     std::optional<double> magneticFluxDensity) {
    auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();

    double initialPermeabilityValue = 1;

    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        initialPermeabilityValue = get_initial_permeability_formula(coreMaterial, temperature, magneticFieldDcBias, frequency, magneticFluxDensity);
    }
    else {
        double initialPermeabilityValueReference = 1;
        bool hasTemperatureRequirement = temperature.has_value();
        bool hasFrequencyRequirement = frequency.has_value();
        bool hasMagneticFieldDcBiasRequirement = magneticFieldDcBias.has_value();
        bool hasMagneticFluxDensityRequirement = magneticFluxDensity.has_value();
        bool hasTemperatureDependency = has_temperature_dependency(coreMaterial);
        bool hasFrequencyDependency = has_frequency_dependency(coreMaterial);
        bool hasMagneticFieldDcBiasDependency = has_magnetic_field_dc_bias_dependency(coreMaterial);
        double temperatureFactor = 1;
        double frequencyFactor = 1;
        double magneticFieldDcBiasFactor = 1;
        double saturationFactor = 1;

        if (hasTemperatureDependency) {
            initialPermeabilityValueReference = get_initial_permeability_temperature_dependent(coreMaterial, Defaults().ambientTemperature);
        }
        else if (hasFrequencyDependency) {
            auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
            auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
            double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
            initialPermeabilityValueReference = get_initial_permeability_frequency_dependent(coreMaterial, minimumFrequency);
        }
        else if (hasMagneticFieldDcBiasDependency) {
            auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
            auto permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
            double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
            initialPermeabilityValueReference = get_initial_permeability_magnetic_field_dc_bias_dependent(coreMaterial, minimumMagneticFieldDcBias);
        }

        if (hasTemperatureDependency && hasTemperatureRequirement) {
            double initialPermeabilityValueTemperatureDependent = get_initial_permeability_temperature_dependent(coreMaterial, temperature.value());
            temperatureFactor = initialPermeabilityValueTemperatureDependent / initialPermeabilityValueReference;
        }

        if (hasFrequencyDependency && hasFrequencyRequirement) {
            double initialPermeabilityValueFrequencyDependent = get_initial_permeability_frequency_dependent(coreMaterial, frequency.value());
            frequencyFactor = initialPermeabilityValueFrequencyDependent / initialPermeabilityValueReference;
        }

        if (hasMagneticFieldDcBiasDependency && hasMagneticFieldDcBiasRequirement) {
            double initialPermeabilityValueMagneticFieldDcBiasDependent = get_initial_permeability_magnetic_field_dc_bias_dependent(coreMaterial, magneticFieldDcBias.value());
            magneticFieldDcBiasFactor = initialPermeabilityValueMagneticFieldDcBiasDependent / initialPermeabilityValueReference;
        }

        // Add variation due to saturation
        if (!hasMagneticFieldDcBiasDependency && (hasMagneticFieldDcBiasRequirement || hasMagneticFluxDensityRequirement)) {
            double amplitudePermeability = 1;
            double auxTemperatude = defaults.ambientTemperature;
            if (temperature) {
                auxTemperatude = temperature.value();
            }
            if(hasMagneticFieldDcBiasRequirement) {
                if (magneticFieldDcBias.value() <= 0) {
                    amplitudePermeability = initialPermeabilityValueReference;
                }
                else {
                    amplitudePermeability = AmplitudePermeability::get_amplitude_permeability(coreMaterial, std::nullopt, magneticFieldDcBias.value(), auxTemperatude);        
                }
            }
            else {
                if (magneticFluxDensity.value() <= 0) {
                    amplitudePermeability = initialPermeabilityValueReference;
                }
                else {
                    amplitudePermeability = AmplitudePermeability::get_amplitude_permeability(coreMaterial, magneticFluxDensity.value(), std::nullopt, auxTemperatude);        
                }
            }
            if (amplitudePermeability < initialPermeabilityValueReference) {
                saturationFactor = amplitudePermeability / initialPermeabilityValueReference;
            }
        }

        initialPermeabilityValue = initialPermeabilityValueReference * temperatureFactor * frequencyFactor * magneticFieldDcBiasFactor * saturationFactor;
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
