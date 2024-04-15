#include "CoreWrapper.h"
#include "InitialPermeability.h"
#include "CoreCrossReferencer.h"
#include "MagneticSimulator.h"
#include "Reluctance.h"
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


std::map<std::string, std::map<CoreCrossReferencer::CoreCrossReferencerFilters, double>> CoreCrossReferencer::get_scorings(bool weighted){
    std::map<std::string, std::map<CoreCrossReferencer::CoreCrossReferencerFilters, double>> swappedScorings;
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

void normalize_scoring(std::vector<std::pair<CoreWrapper, double>>* rankedCores, std::vector<double>* newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    double maximumScoring = *std::max_element(newScoring->begin(), newScoring->end());
    double minimumScoring = *std::min_element(newScoring->begin(), newScoring->end());

    for (size_t i = 0; i < (*rankedCores).size(); ++i) {
        auto mas = (*rankedCores)[i].first;
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
                    (*rankedCores)[i].second = (*rankedCores)[i].second + weight * (1 - (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring)));
                }
                else {
                    (*rankedCores)[i].second = (*rankedCores)[i].second + weight * (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                }
            }
            else {
                if (filterConfiguration["invert"]) {
                    (*rankedCores)[i].second = (*rankedCores)[i].second + weight * (1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring));
                }
                else {
                    (*rankedCores)[i].second = (*rankedCores)[i].second + weight * (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                }
            }
        }
        else {
            (*rankedCores)[i].second = (*rankedCores)[i].second + 1;
        }
    }
    sort((*rankedCores).begin(), (*rankedCores).end(), [](std::pair<CoreWrapper, double>& b1, std::pair<CoreWrapper, double>& b2) {
        return b1.second > b2.second;
    });
}

std::vector<std::pair<CoreWrapper, double>> CoreCrossReferencer::MagneticCoreFilterPermeance::filter_core(std::vector<std::pair<CoreWrapper, double>>* unfilteredCores, CoreWrapper referenceCore, InputsWrapper inputs, std::map<std::string, std::string> models, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }
    if (models.find("gapReluctance") == models.end()) {
        models["gapReluctance"] = magic_enum::enum_name(Defaults().reluctanceModelDefault);
    }

    std::vector<double> newScoring;
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(models);
    std::vector<std::pair<CoreWrapper, double>> filteredCoresWithScoring;
    double referenceReluctance = 0;

    if (inputs.get_operating_points()[0].get_excitations_per_winding().size() > 0) {
        for (auto operatingPoint : inputs.get_operating_points()) {
            referenceReluctance += reluctanceModel->get_core_reluctance(referenceCore, &operatingPoint).get_core_reluctance();
        }
        referenceReluctance /= inputs.get_operating_points().size();
    }
    else {
        referenceReluctance = reluctanceModel->get_core_reluctance(referenceCore).get_core_reluctance();
    }

    std::list<size_t> listOfIndexesToErase;
    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        CoreWrapper core = (*unfilteredCores)[coreIndex].first;
        double reluctance = 0;
        if (inputs.get_operating_points()[0].get_excitations_per_winding().size() > 0) {
            for (auto operatingPoint : inputs.get_operating_points()) {
                reluctance += reluctanceModel->get_core_reluctance(core, &operatingPoint).get_core_reluctance();
            }
            reluctance /= inputs.get_operating_points().size();
        }
        else {
            reluctance = reluctanceModel->get_core_reluctance(core).get_core_reluctance();
        }

        if (fabs(referenceReluctance - reluctance) / referenceReluctance < limit) {
            double scoring = fabs(referenceReluctance - reluctance);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencer::CoreCrossReferencerFilters::PERMEANCE, scoring);
        }
        else {
            listOfIndexesToErase.push_back(coreIndex);
        }

    }

    for (size_t i = 0; i < (*unfilteredCores).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoresWithScoring.push_back((*unfilteredCores)[i]);
        }
    }

    if (filteredCoresWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCores: " + std::to_string((*unfilteredCores).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoresWithScoring.size() > 0) {
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencer::CoreCrossReferencerFilters::PERMEANCE]);
    }
    return filteredCoresWithScoring;
}

std::vector<std::pair<CoreWrapper, double>> CoreCrossReferencer::MagneticCoreFilterWindingWindowArea::filter_core(std::vector<std::pair<CoreWrapper, double>>* unfilteredCores, CoreWrapper referenceCore, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }

    std::vector<double> newScoring;
    std::vector<std::pair<CoreWrapper, double>> filteredCoresWithScoring;
    if (!referenceCore.get_winding_windows()[0].get_area()) {
        throw std::runtime_error("Winding windong is missing area");
    }
    double referenceWindingWindowArea = referenceCore.get_winding_windows()[0].get_area().value();

    std::list<size_t> listOfIndexesToErase;
    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        CoreWrapper core = (*unfilteredCores)[coreIndex].first;
        if (!core.get_winding_windows()[0].get_area()) {
            throw std::runtime_error("Winding windong is missing area");
        }
        double windingWindowArea = core.get_winding_windows()[0].get_area().value();

        if (fabs(referenceWindingWindowArea - windingWindowArea) / referenceWindingWindowArea < limit) {
            double scoring = fabs(referenceWindingWindowArea - windingWindowArea);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencer::CoreCrossReferencerFilters::WINDING_WINDOW_AREA, scoring);
        }
        else {
            listOfIndexesToErase.push_back(coreIndex);
        }

    }

    for (size_t i = 0; i < (*unfilteredCores).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoresWithScoring.push_back((*unfilteredCores)[i]);
        }
    }

    if (filteredCoresWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCores: " + std::to_string((*unfilteredCores).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoresWithScoring.size() > 0) {
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencer::CoreCrossReferencerFilters::WINDING_WINDOW_AREA]);
    }
    return filteredCoresWithScoring;
}

std::vector<std::pair<CoreWrapper, double>> CoreCrossReferencer::MagneticCoreFilterDimensions::filter_core(std::vector<std::pair<CoreWrapper, double>>* unfilteredCores, CoreWrapper referenceCore, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }

    std::vector<double> newScoring;
    std::vector<std::pair<CoreWrapper, double>> filteredCoresWithScoring;
    double referenceDepth = referenceCore.get_depth();
    double referenceHeight = referenceCore.get_height();
    double referenceWidth = referenceCore.get_width();

    std::list<size_t> listOfIndexesToErase;
    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        CoreWrapper core = (*unfilteredCores)[coreIndex].first;
        double depth = core.get_depth();
        double height = core.get_height();
        double width = core.get_width();

        if (fabs(referenceDepth - depth) / referenceDepth < limit && fabs(referenceHeight - height) / referenceHeight < limit && fabs(referenceWidth - width) / referenceWidth < limit) {
            double scoring = fabs(referenceDepth - depth) + fabs(referenceHeight - height) + fabs(referenceWidth - width);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencer::CoreCrossReferencerFilters::DIMENSIONS, scoring);
        }
        else {
            listOfIndexesToErase.push_back(coreIndex);
        }

    }

    for (size_t i = 0; i < (*unfilteredCores).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoresWithScoring.push_back((*unfilteredCores)[i]);
        }
    }

    if (filteredCoresWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCores: " + std::to_string((*unfilteredCores).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoresWithScoring.size() > 0) {
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencer::CoreCrossReferencerFilters::DIMENSIONS]);
    }
    return filteredCoresWithScoring;
}

std::pair<double, double> CoreCrossReferencer::MagneticCoreFilterCoreLosses::calculate_average_core_losses_and_magnetic_flux_density(CoreWrapper core, int64_t numberTurns, InputsWrapper inputs, std::map<std::string, std::string> models) {
    auto defaults = Defaults();
    if (models.find("coreLosses") == models.end()) {
        models["coreLosses"] = magic_enum::enum_name(Defaults().reluctanceModelDefault);
    }

    auto temperature = inputs.get_maximum_temperature();

    OperatingPointExcitation excitation;
    SignalDescriptor magneticFluxDensity;
    Processed magneticFluxDensityProcessed;
    magneticFluxDensityProcessed.set_label(WaveformLabel::SINUSOIDAL);
    magneticFluxDensityProcessed.set_offset(0);
    magneticFluxDensityProcessed.set_duty_cycle(0.5);
    std::shared_ptr<CoreLossesModel> coreLossesModelForMaterial = nullptr;

    auto reluctanceModelName = defaults.reluctanceModelDefault;
    if (models.find("reluctance") != models.end()) {
        std::string modelNameStringUpper = models["reluctance"];
        std::transform(modelNameStringUpper.begin(), modelNameStringUpper.end(), modelNameStringUpper.begin(), ::toupper);
        reluctanceModelName = magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(modelNameStringUpper).value();
    }
    auto coreLossesModelName = defaults.coreLossesModelDefault;
    if (models.find("coreLosses") != models.end()) {
        std::string modelNameStringUpper = models["coreLosses"];
        std::transform(modelNameStringUpper.begin(), modelNameStringUpper.end(), modelNameStringUpper.begin(), ::toupper);
        coreLossesModelName = magic_enum::enum_cast<OpenMagnetics::CoreLossesModels>(modelNameStringUpper).value();
    }
    auto coreTemperatureModelName = defaults.coreTemperatureModelDefault;
    if (models.find("coreTemperature") != models.end()) {
        std::string modelNameStringUpper = models["coreTemperature"];
        std::transform(modelNameStringUpper.begin(), modelNameStringUpper.end(), modelNameStringUpper.begin(), ::toupper);
        coreTemperatureModelName = magic_enum::enum_cast<OpenMagnetics::CoreTemperatureModels>(modelNameStringUpper).value();
    }

    OpenMagnetics::MagneticSimulator magneticSimulator;
    magneticSimulator.set_core_losses_model_name(coreLossesModelName);
    magneticSimulator.set_core_temperature_model_name(coreTemperatureModelName);
    magneticSimulator.set_reluctance_model_name(reluctanceModelName);

    try {

        auto availableMethodsForMaterial = CoreLossesModel::get_methods(core.resolve_material());
        for (auto& [modelName, coreLossesModel] : _coreLossesModels) {
            if (std::find(availableMethodsForMaterial.begin(), availableMethodsForMaterial.end(), modelName) != availableMethodsForMaterial.end()) {
                coreLossesModelForMaterial = coreLossesModel;
            }
        }

        if (coreLossesModelForMaterial == nullptr) {
            throw std::runtime_error("No model found for material: " + core.resolve_material().get_name());
        }
        double averageCoreLosses = 0;
        double maximumMagneticFluxDensitySaturationPeak = 0;

        if (inputs.get_operating_points()[0].get_excitations_per_winding().size() > 0) {
            Magnetic magnetic;
            magnetic.set_core(core);
            Coil coil;
            coil.set_bobbin("Dummy");
            CoilFunctionalDescription coilFunctionalDescription;
            coilFunctionalDescription.set_number_turns(numberTurns);
            coilFunctionalDescription.set_wire("Dummy");
            coil.set_functional_description({coilFunctionalDescription});
            magnetic.set_coil(coil);
            for (auto operatingPoint : inputs.get_operating_points()) {
                auto coreLossesOutput = magneticSimulator.calculate_core_losses(operatingPoint, magnetic);
                averageCoreLosses += coreLossesOutput.get_core_losses();
                double magneticFluxDensityPeak = coreLossesOutput.get_magnetic_flux_density()->get_processed()->get_peak().value();
                maximumMagneticFluxDensitySaturationPeak = std::max(maximumMagneticFluxDensitySaturationPeak, magneticFluxDensityPeak);
            }
            averageCoreLosses /= inputs.get_operating_points().size();
        }
        else {
            for (auto magneticFluxDensityPeak : _magneticFluxDensities) {
                magneticFluxDensityProcessed.set_peak(magneticFluxDensityPeak);
                magneticFluxDensityProcessed.set_peak_to_peak(magneticFluxDensityPeak * 2);
                magneticFluxDensity.set_processed(magneticFluxDensityProcessed);
                for (auto frequency : _frequencies) {
                    magneticFluxDensity.set_waveform(InputsWrapper::create_waveform(magneticFluxDensityProcessed, frequency));
                    excitation.set_frequency(frequency);
                    excitation.set_magnetic_flux_density(magneticFluxDensity);
                    double coreVolumetricLosses = coreLossesModelForMaterial->get_core_volumetric_losses(core.resolve_material(), excitation, temperature);
                    averageCoreLosses += coreVolumetricLosses * core.get_processed_description()->get_effective_parameters().get_effective_volume();
                }
            }
            averageCoreLosses /= _magneticFluxDensities.size() * _frequencies.size();
        }
        
        return {averageCoreLosses, maximumMagneticFluxDensitySaturationPeak};
    }
    catch(const std::runtime_error& re)
    {
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
}

std::vector<std::pair<CoreWrapper, double>> CoreCrossReferencer::MagneticCoreFilterCoreLosses::filter_core(std::vector<std::pair<CoreWrapper, double>>* unfilteredCores, CoreWrapper referenceCore, int64_t referenceNumberTurns, InputsWrapper inputs, std::map<std::string, std::string> models, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }
    std::vector<double> newScoring;

    auto [referenceCoreLossesWithTemperature, referenceCoreMagneticFluxDensitySaturationPeak] = calculate_average_core_losses_and_magnetic_flux_density(referenceCore, referenceNumberTurns, inputs, models);

    OperatingPointExcitation excitation;
    SignalDescriptor magneticFluxDensity;
    Processed magneticFluxDensityProcessed;
    magneticFluxDensityProcessed.set_label(WaveformLabel::SINUSOIDAL);
    magneticFluxDensityProcessed.set_offset(0);
    magneticFluxDensityProcessed.set_duty_cycle(0.5);

    std::vector<std::pair<CoreWrapper, double>> filteredCoresWithScoring;
    std::list<size_t> listOfIndexesToErase;

    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        CoreWrapper core = (*unfilteredCores)[coreIndex].first;
        auto temperature = inputs.get_maximum_temperature();

        auto magneticFluxDensitySaturationPeak = core.get_magnetic_flux_density_saturation(temperature, true);

        auto [coreLossesWithTemperature, coreMagneticFluxDensitySaturationPeak] = calculate_average_core_losses_and_magnetic_flux_density(core, referenceNumberTurns, inputs, models);

        if ((coreMagneticFluxDensitySaturationPeak < magneticFluxDensitySaturationPeak) && (coreLossesWithTemperature < referenceCoreLossesWithTemperature)) {
            double scoring = 0;
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencer::CoreCrossReferencerFilters::CORE_LOSSES, scoring);
        }
        else if ((coreMagneticFluxDensitySaturationPeak < magneticFluxDensitySaturationPeak) && (fabs(referenceCoreLossesWithTemperature - coreLossesWithTemperature) / referenceCoreLossesWithTemperature < limit) || limit >= 1) {
            double scoring = fabs(referenceCoreLossesWithTemperature - coreLossesWithTemperature);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencer::CoreCrossReferencerFilters::CORE_LOSSES, scoring);
        }
        else {
            listOfIndexesToErase.push_back(coreIndex);
        }
    }


    for (size_t i = 0; i < (*unfilteredCores).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoresWithScoring.push_back((*unfilteredCores)[i]);
        }
    }

    if (filteredCoresWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCores: " + std::to_string((*unfilteredCores).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoresWithScoring.size() > 0) {
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencer::CoreCrossReferencerFilters::PERMEANCE]);
    }
    return filteredCoresWithScoring;
}

std::vector<std::pair<CoreWrapper, double>> CoreCrossReferencer::get_cross_referenced_core(CoreWrapper referenceCore, int64_t referenceNumberTurns, InputsWrapper inputs, size_t maximumNumberResults) {
    return get_cross_referenced_core(referenceCore, referenceNumberTurns, inputs, _weights, maximumNumberResults);
}

std::vector<std::pair<CoreWrapper, double>> CoreCrossReferencer::get_cross_referenced_core(CoreWrapper referenceCore, int64_t referenceNumberTurns, InputsWrapper inputs, std::map<CoreCrossReferencerFilters, double> weights, size_t maximumNumberResults) {
    auto defaults = Defaults();
    _weights = weights;

    std::vector<std::pair<CoreWrapper, double>> cores;
    std::string referenceName = referenceCore.get_name().value();
    for (auto core : coreDatabase){
        if (referenceName != core.get_name().value()) {
            if (!_onlyManufacturer || core.get_manufacturer_info()->get_name() == _onlyManufacturer.value()) {
                if (!core.get_processed_description()) {
                    core.process_data();
                    core.process_gap();
                }
                cores.push_back({core, 0.0});
            }
        }
    }

    if (!referenceCore.get_processed_description()) {
        referenceCore.process_data();
        referenceCore.process_gap();
    }

    double limit = 0;
    std::vector<std::pair<CoreWrapper, double>> filteredCores;

    while (limit <= 1 && filteredCores.size() < maximumNumberResults) {
        limit += 0.25;
        filteredCores = apply_filters(&cores, referenceCore, referenceNumberTurns, inputs, weights, maximumNumberResults, limit);
    }

    return filteredCores;
}

std::vector<std::pair<CoreWrapper, double>> CoreCrossReferencer::apply_filters(std::vector<std::pair<CoreWrapper, double>>* cores, CoreWrapper referenceCore, int64_t referenceNumberTurns, InputsWrapper inputs, std::map<CoreCrossReferencerFilters, double> weights, size_t maximumNumberResults, double limit) {
    MagneticCoreFilterPermeance filterPermeance;
    MagneticCoreFilterCoreLosses filterVolumetricLosses;
    MagneticCoreFilterWindingWindowArea filterWindingWindowArea;
    MagneticCoreFilterDimensions filterDimensions;

    filterPermeance.set_scorings(&_scorings);
    filterPermeance.set_filter_configuration(&_filterConfiguration);
    filterVolumetricLosses.set_scorings(&_scorings);
    filterVolumetricLosses.set_filter_configuration(&_filterConfiguration);
    filterWindingWindowArea.set_scorings(&_scorings);
    filterWindingWindowArea.set_filter_configuration(&_filterConfiguration);
    filterDimensions.set_scorings(&_scorings);
    filterDimensions.set_filter_configuration(&_filterConfiguration);

    std::vector<std::pair<CoreWrapper, double>> rankedCores = *cores;

    magic_enum::enum_for_each<CoreCrossReferencerFilters>([&] (auto val) {
        CoreCrossReferencerFilters filter = val;
        switch (filter) {
            case CoreCrossReferencerFilters::DIMENSIONS: 
                rankedCores = filterDimensions.filter_core(&rankedCores, referenceCore, weights[CoreCrossReferencerFilters::DIMENSIONS], limit);
                std::cout << "DIMENSIONS rankedCores: " << rankedCores.size() << std::endl;
                break;
            case CoreCrossReferencerFilters::WINDING_WINDOW_AREA: 
                rankedCores = filterWindingWindowArea.filter_core(&rankedCores, referenceCore, weights[CoreCrossReferencerFilters::WINDING_WINDOW_AREA], limit);
                std::cout << "WINDING_WINDOW_AREA rankedCores: " << rankedCores.size() << std::endl;
                break;
            case CoreCrossReferencerFilters::PERMEANCE: 
                rankedCores = filterPermeance.filter_core(&rankedCores, referenceCore, inputs, _models, weights[CoreCrossReferencerFilters::PERMEANCE], limit);
                std::cout << "PERMEANCE rankedCores: " << rankedCores.size() << std::endl;
                break;
        }    
        std::string filterString = std::string{magic_enum::enum_name(filter)};
        logEntry("There are " + std::to_string(rankedCores.size()) + " after filtering by " + filterString + ".");
    });

    // We leave core losses for the last one, as it is the most computationally costly
    rankedCores = filterVolumetricLosses.filter_core(&rankedCores, referenceCore, referenceNumberTurns, inputs, _models, weights[CoreCrossReferencerFilters::CORE_LOSSES], limit);
    std::cout << "CORE_LOSSES rankedCores: " << rankedCores.size() << std::endl;
        logEntry("There are " + std::to_string(rankedCores.size()) + " after filtering by " + std::string{magic_enum::enum_name(CoreCrossReferencerFilters::CORE_LOSSES)} + ".");

    if (rankedCores.size() > maximumNumberResults) {
        rankedCores = std::vector<std::pair<CoreWrapper, double>>(rankedCores.begin(), rankedCores.end() - (rankedCores.size() - maximumNumberResults));
    }

    return rankedCores;
}

} // namespace OpenMagnetics
