#include "CoreWrapper.h"
#include "InitialPermeability.h"
#include "CoreMaterialCrossReferencer.h"
#include "Defaults.h"
#include "Settings.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <string>
#include <vector>
#include <execution>
#include <magic_enum_utility.hpp>
#include <list>
#include <cmrc/cmrc.hpp>

CMRC_DECLARE(data);


namespace OpenMagnetics {


std::map<std::string, std::map<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters, double>> CoreMaterialCrossReferencer::get_scorings(bool weighted){
    std::map<std::string, std::map<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters, double>> swappedScorings;
    for (auto& [filter, aux] : _scorings) {
        auto filterConfiguration = _filterConfiguration[filter];

        double maximumScoring = (*std::max_element(aux.begin(), aux.end(),
                                     [](const std::pair<std::string, double> &p1,
                                        const std::pair<std::string, double> &p2)
                                     {
                                         return p1.second < p2.second;
                                     })).second; 
        double minimumScoring = (*std::min_element(aux.begin(), aux.end(),
                                     [](const std::pair<std::string, double> &p1,
                                        const std::pair<std::string, double> &p2)
                                     {
                                         return p1.second < p2.second;
                                     })).second; 

        for (auto& [name, scoring] : aux) {
            if (filterConfiguration["log"]){
                if (filterConfiguration["invert"]) {
                    if (weighted) {
                        swappedScorings[name][filter] = _weights[filter] * (1 - (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring)));
                    }
                    else {
                        swappedScorings[name][filter] = 1 - (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                    }
                }
                else {
                    if (weighted) {
                        swappedScorings[name][filter] = _weights[filter] * (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                    }
                    else {
                        swappedScorings[name][filter] = (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                    }
                }
            }
            else {
                if (filterConfiguration["invert"]) {
                    if (weighted) {
                        swappedScorings[name][filter] = _weights[filter] * (1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring));
                    }
                    else {
                        swappedScorings[name][filter] = 1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                    }
                }
                else {
                    if (weighted) {
                        swappedScorings[name][filter] = _weights[filter] * (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                    }
                    else {
                        swappedScorings[name][filter] = (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                    }
                }
            }
        }
    }
    return swappedScorings;
}

void normalize_scoring(std::vector<std::pair<CoreMaterial, double>>* rankedCoreMaterials, std::vector<double>* newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    double maximumScoring = *std::max_element(newScoring->begin(), newScoring->end());
    double minimumScoring = *std::min_element(newScoring->begin(), newScoring->end());

    for (size_t i = 0; i < (*rankedCoreMaterials).size(); ++i) {
        auto mas = (*rankedCoreMaterials)[i].first;
        auto scoring = (*newScoring)[i];
        if (std::isnan(scoring)) {
            scoring = maximumScoring;
        }
        else {
            scoring = std::max(0.0001, scoring);
        }
        minimumScoring = std::max(0.0001, minimumScoring);
        if (maximumScoring != minimumScoring) {

            if (filterConfiguration["log"]){
                if (filterConfiguration["invert"]) {
                    (*rankedCoreMaterials)[i].second = (*rankedCoreMaterials)[i].second + weight * (1 - (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring)));
                }
                else {
                    (*rankedCoreMaterials)[i].second = (*rankedCoreMaterials)[i].second + weight * (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                }
            }
            else {
                if (filterConfiguration["invert"]) {
                    (*rankedCoreMaterials)[i].second = (*rankedCoreMaterials)[i].second + weight * (1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring));
                }
                else {
                    (*rankedCoreMaterials)[i].second = (*rankedCoreMaterials)[i].second + weight * (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                }
            }
        }
        else {
            (*rankedCoreMaterials)[i].second = (*rankedCoreMaterials)[i].second + 1;
        }
    }
    sort((*rankedCoreMaterials).begin(), (*rankedCoreMaterials).end(), [](std::pair<CoreMaterial, double>& b1, std::pair<CoreMaterial, double>& b2) {
        return b1.second > b2.second;
    });
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterInitialPermeability::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;
    OpenMagnetics::InitialPermeability initialPermeability;

    double referenceInitialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(referenceCoreMaterial, temperature, std::nullopt, std::nullopt);

    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(coreMaterial, temperature, std::nullopt, std::nullopt);

        double scoring = fabs(referenceInitialPermeabilityValueWithTemperature - initialPermeabilityValueWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY, scoring);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterRemanence::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceRemanenceWithTemperature = CoreWrapper::get_remanence(referenceCoreMaterial, temperature);

    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double remanenceWithTemperature = CoreWrapper::get_remanence(coreMaterial, temperature);

        double scoring = fabs(referenceRemanenceWithTemperature - remanenceWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::REMANENCE, scoring);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::REMANENCE]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterCoerciveForce::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceCoerciveForceWithTemperature = CoreWrapper::get_coercive_force(referenceCoreMaterial, temperature);

    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double coerciveForceWithTemperature = CoreWrapper::get_coercive_force(coreMaterial, temperature);

        double scoring = fabs(referenceCoerciveForceWithTemperature - coerciveForceWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::COERCIVE_FORCE, scoring);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::COERCIVE_FORCE]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterSaturation::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceSaturationWithTemperature = CoreWrapper::get_magnetic_flux_density_saturation(referenceCoreMaterial, temperature);

    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double saturationWithTemperature = CoreWrapper::get_magnetic_flux_density_saturation(coreMaterial, temperature);

        double scoring = fabs(referenceSaturationWithTemperature - saturationWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::SATURATION, scoring);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::SATURATION]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterCurieTemperature::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceCurieTemperatureWithTemperature = CoreWrapper::get_curie_temperature(referenceCoreMaterial);

    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double curieTemperatureWithTemperature = CoreWrapper::get_curie_temperature(coreMaterial);

        double scoring = fabs(referenceCurieTemperatureWithTemperature - curieTemperatureWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE, scoring);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE]);
    }
    return *unfilteredCoreMaterials;
}

double CoreMaterialCrossReferencer::MagneticCoreFilterVolumetricLosses::calculate_average_volumetric_losses(CoreMaterial coreMaterial, double temperature, std::map<std::string, std::string> models) {
    if (models.find("coreLosses") == models.end()) {
        models["coreLosses"] = magic_enum::enum_name(Defaults().reluctanceModelDefault);
    }

    OperatingPointExcitation excitation;
    SignalDescriptor magneticFluxDensity;
    Processed magneticFluxDensityProcessed;
    magneticFluxDensityProcessed.set_label(WaveformLabel::SINUSOIDAL);
    magneticFluxDensityProcessed.set_offset(0);
    magneticFluxDensityProcessed.set_duty_cycle(0.5);
        std::shared_ptr<CoreLossesModel> coreLossesModelForMaterial = nullptr;

    try {
        auto availableMethodsForMaterial = CoreLossesModel::get_methods(coreMaterial);
        for (auto& [modelName, coreLossesModel] : _coreLossesModels) {
            if (std::find(availableMethodsForMaterial.begin(), availableMethodsForMaterial.end(), modelName) != availableMethodsForMaterial.end()) {
                coreLossesModelForMaterial = coreLossesModel;
            }
        }

        if (coreLossesModelForMaterial == nullptr) {
            throw std::runtime_error("No model found for material: " + coreMaterial.get_name());
        }

        double averageVolumetricLosses = 0;
        for (auto magneticFluxDensityPeak : _magneticFluxDensities) {
            magneticFluxDensityProcessed.set_peak(magneticFluxDensityPeak);
            magneticFluxDensityProcessed.set_peak_to_peak(magneticFluxDensityPeak * 2);
            magneticFluxDensity.set_processed(magneticFluxDensityProcessed);
            for (auto frequency : _frequencies) {
                magneticFluxDensity.set_waveform(InputsWrapper::create_waveform(magneticFluxDensityProcessed, frequency));
                excitation.set_frequency(frequency);
                excitation.set_magnetic_flux_density(magneticFluxDensity);
                double coreVolumetricLosses = coreLossesModelForMaterial->get_core_volumetric_losses(coreMaterial, excitation, temperature);
                averageVolumetricLosses += coreVolumetricLosses;
            }
        }
        averageVolumetricLosses /= _magneticFluxDensities.size() * _frequencies.size();
        return averageVolumetricLosses;
    }
    catch(const std::runtime_error& re)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterVolumetricLosses::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, std::map<std::string, std::string> models, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceVolumetricLossesWithTemperature = calculate_average_volumetric_losses(referenceCoreMaterial, temperature, models);

    OperatingPointExcitation excitation;
    SignalDescriptor magneticFluxDensity;
    Processed magneticFluxDensityProcessed;
    magneticFluxDensityProcessed.set_label(WaveformLabel::SINUSOIDAL);
    magneticFluxDensityProcessed.set_offset(0);
    magneticFluxDensityProcessed.set_duty_cycle(0.5);

    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;
        double volumetricLossesWithTemperature = calculate_average_volumetric_losses(coreMaterial, temperature, models);

        double scoring = fabs(referenceVolumetricLossesWithTemperature - volumetricLossesWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES, scoring);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterResistivity::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceResistivityWithTemperature = CoreWrapper::get_resistivity(referenceCoreMaterial, temperature);

    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double resistivityWithTemperature = CoreWrapper::get_resistivity(coreMaterial, temperature);

        double scoring = fabs(referenceResistivityWithTemperature - resistivityWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::RESISTIVITY, scoring);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::RESISTIVITY]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterDensity::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceDensityWithTemperature = CoreWrapper::get_density(referenceCoreMaterial);

    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double densityWithTemperature = CoreWrapper::get_density(coreMaterial);

        double scoring = fabs(referenceDensityWithTemperature - densityWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::DENSITY, scoring);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::DENSITY]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::get_cross_referenced_core_material(CoreMaterial referenceCoreMaterial, double temperature, size_t maximumNumberResults) {
    return get_cross_referenced_core_material(referenceCoreMaterial, temperature, _weights, maximumNumberResults);
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::get_cross_referenced_core_material(CoreMaterial referenceCoreMaterial, double temperature, std::map<CoreMaterialCrossReferencerFilters, double> weights, size_t maximumNumberResults) {
    auto defaults = Defaults();
    _weights = weights;

    std::vector<std::pair<CoreMaterial, double>> coreMaterials;

    for (auto [name, coreMaterial] : coreMaterialDatabase){
        if (name != referenceCoreMaterial.get_name()) {
            if (!_onlyManufacturer || coreMaterial.get_manufacturer_info().get_name() == _onlyManufacturer.value()) {
                coreMaterials.push_back({coreMaterial, 0.0});
            }
        }
    }

    auto filteredCoreMaterials = apply_filters(&coreMaterials, referenceCoreMaterial, temperature, weights, maximumNumberResults);
    return filteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::apply_filters(std::vector<std::pair<CoreMaterial, double>>* coreMaterials, CoreMaterial referenceCoreMaterial, double temperature, std::map<CoreMaterialCrossReferencerFilters, double> weights, size_t maximumNumberResults) {
    MagneticCoreFilterInitialPermeability filterInitialPermeability;
    MagneticCoreFilterRemanence filterRemanence;
    MagneticCoreFilterCoerciveForce filterCoerciveForce;
    MagneticCoreFilterSaturation filterSaturation;
    MagneticCoreFilterCurieTemperature filterCurieTemperature;
    MagneticCoreFilterVolumetricLosses filterVolumetricLosses;
    MagneticCoreFilterResistivity filterResistivity;
    MagneticCoreFilterDensity filterDensity;

    filterInitialPermeability.set_scorings(&_scorings);
    filterInitialPermeability.set_filter_configuration(&_filterConfiguration);
    filterRemanence.set_scorings(&_scorings);
    filterRemanence.set_filter_configuration(&_filterConfiguration);
    filterCoerciveForce.set_scorings(&_scorings);
    filterCoerciveForce.set_filter_configuration(&_filterConfiguration);
    filterSaturation.set_scorings(&_scorings);
    filterSaturation.set_filter_configuration(&_filterConfiguration);
    filterCurieTemperature.set_scorings(&_scorings);
    filterCurieTemperature.set_filter_configuration(&_filterConfiguration);
    filterVolumetricLosses.set_scorings(&_scorings);
    filterVolumetricLosses.set_filter_configuration(&_filterConfiguration);
    filterResistivity.set_scorings(&_scorings);
    filterResistivity.set_filter_configuration(&_filterConfiguration);
    filterDensity.set_scorings(&_scorings);
    filterDensity.set_filter_configuration(&_filterConfiguration);

    std::vector<std::pair<CoreMaterial, double>> rankedCoreMaterials = *coreMaterials;

    magic_enum::enum_for_each<CoreMaterialCrossReferencerFilters>([&] (auto val) {
        CoreMaterialCrossReferencerFilters filter = val;
        switch (filter) {
            case CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY: 
                rankedCoreMaterials = filterInitialPermeability.filter_core_materials(&rankedCoreMaterials, referenceCoreMaterial, temperature, weights[CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY]);
                break;
            case CoreMaterialCrossReferencerFilters::REMANENCE: 
                rankedCoreMaterials = filterRemanence.filter_core_materials(&rankedCoreMaterials, referenceCoreMaterial, temperature, weights[CoreMaterialCrossReferencerFilters::REMANENCE]);
                break;
            case CoreMaterialCrossReferencerFilters::COERCIVE_FORCE: 
                rankedCoreMaterials = filterCoerciveForce.filter_core_materials(&rankedCoreMaterials, referenceCoreMaterial, temperature, weights[CoreMaterialCrossReferencerFilters::COERCIVE_FORCE]);
                break;
            case CoreMaterialCrossReferencerFilters::SATURATION: 
                rankedCoreMaterials = filterSaturation.filter_core_materials(&rankedCoreMaterials, referenceCoreMaterial, temperature, weights[CoreMaterialCrossReferencerFilters::SATURATION]);
                break;
            case CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE: 
                rankedCoreMaterials = filterCurieTemperature.filter_core_materials(&rankedCoreMaterials, referenceCoreMaterial, temperature, weights[CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE]);
                break;
            case CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES: 
                rankedCoreMaterials = filterVolumetricLosses.filter_core_materials(&rankedCoreMaterials, referenceCoreMaterial, temperature, _models, weights[CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES]);
                break;
            case CoreMaterialCrossReferencerFilters::RESISTIVITY: 
                rankedCoreMaterials = filterResistivity.filter_core_materials(&rankedCoreMaterials, referenceCoreMaterial, temperature, weights[CoreMaterialCrossReferencerFilters::RESISTIVITY]);
                break;
            case CoreMaterialCrossReferencerFilters::DENSITY: 
                rankedCoreMaterials = filterDensity.filter_core_materials(&rankedCoreMaterials, referenceCoreMaterial, temperature, weights[CoreMaterialCrossReferencerFilters::DENSITY]);
                break;
        }    
        std::string filterString = std::string{magic_enum::enum_name(filter)};
        logEntry("There are " + std::to_string(rankedCoreMaterials.size()) + " after filtering by " + filterString + ".");
    });

    if (rankedCoreMaterials.size() > maximumNumberResults) {
        rankedCoreMaterials = std::vector<std::pair<CoreMaterial, double>>(rankedCoreMaterials.begin(), rankedCoreMaterials.end() - (rankedCoreMaterials.size() - maximumNumberResults));
    }

    return rankedCoreMaterials;
}

} // namespace OpenMagnetics
