#include "constructive_models/Core.h"
#include "physical_models/InitialPermeability.h"
#include "advisers/CoreCrossReferencer.h"
#include "processors/MagneticSimulator.h"
#include "physical_models/Reluctance.h"
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


std::map<std::string, std::map<CoreCrossReferencerFilters, double>> CoreCrossReferencer::get_scorings(bool weighted){
    std::map<std::string, std::map<CoreCrossReferencerFilters, double>> swappedScorings;
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
        minimumScoring = std::max(0.0001, minimumScoring);

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
std::map<std::string, std::map<CoreCrossReferencerFilters, double>> CoreCrossReferencer::get_scored_values(){
    std::map<std::string, std::map<CoreCrossReferencerFilters, double>> swappedScoredValues;
    for (auto& [filter, aux] : _scoredValues) {
        auto filterConfiguration = _filterConfiguration[filter];

        for (auto& [name, scoredValue] : aux) {
            swappedScoredValues[name][filter] = scoredValue;
        }
    }
    return swappedScoredValues;
}

void normalize_scoring(std::vector<std::pair<Core, double>>* rankedCores, std::vector<double>* newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
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
    sort((*rankedCores).begin(), (*rankedCores).end(), [](std::pair<Core, double>& b1, std::pair<Core, double>& b2) {
        return b1.second > b2.second;
    });
}

std::vector<std::pair<Core, double>> CoreCrossReferencer::MagneticCoreFilterPermeance::filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, Inputs inputs, std::map<std::string, std::string> models, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }
    if (models.find("gapReluctance") == models.end()) {
        models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
    }

    std::vector<double> newScoring;
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(models);
    std::vector<std::pair<Core, double>> filteredCoresWithScoring;
    double referencePermeance = 0;

    if (inputs.get_operating_points()[0].get_excitations_per_winding().size() > 0) {
        for (auto operatingPoint : inputs.get_operating_points()) {
            referencePermeance += 1.0 / reluctanceModel->get_core_reluctance(referenceCore, operatingPoint).get_core_reluctance();
        }
        referencePermeance /= inputs.get_operating_points().size();
    }
    else {
        referencePermeance = 1.0 / reluctanceModel->get_core_reluctance(referenceCore).get_core_reluctance();
    }
    add_scored_value("Reference", CoreCrossReferencerFilters::PERMEANCE, referencePermeance);


    std::list<size_t> listOfIndexesToErase;
    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        Core core = (*unfilteredCores)[coreIndex].first;

        if ((*_validScorings).contains(CoreCrossReferencerFilters::PERMEANCE)) {
            if ((*_validScorings)[CoreCrossReferencerFilters::PERMEANCE].contains(core.get_name().value())) {
                if ((*_validScorings)[CoreCrossReferencerFilters::PERMEANCE][core.get_name().value()]) {
                    newScoring.push_back((*_scorings)[CoreCrossReferencerFilters::PERMEANCE][core.get_name().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(coreIndex);
                }
                continue;
            }
        }

        double reluctance = 0;
        if (inputs.get_operating_points()[0].get_excitations_per_winding().size() > 0) {
            for (auto operatingPoint : inputs.get_operating_points()) {
                reluctance += 1.0 / reluctanceModel->get_core_reluctance(core, operatingPoint).get_core_reluctance();
            }
            reluctance /= inputs.get_operating_points().size();
        }
        else {
            reluctance = 1.0 / reluctanceModel->get_core_reluctance(core).get_core_reluctance();
        }

        if (fabs(referencePermeance - reluctance) / referencePermeance < limit) {
            double scoring = fabs(referencePermeance - reluctance);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencerFilters::PERMEANCE, scoring);
            add_scored_value(core.get_name().value(), CoreCrossReferencerFilters::PERMEANCE, reluctance);
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
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencerFilters::PERMEANCE]);
    }
    return filteredCoresWithScoring;
}

std::vector<std::pair<Core, double>> CoreCrossReferencer::MagneticCoreFilterWindingWindowArea::filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }

    std::vector<double> newScoring;
    std::vector<std::pair<Core, double>> filteredCoresWithScoring;
    if (!referenceCore.get_winding_windows()[0].get_area()) {
        throw std::runtime_error("Winding window is missing area");
    }
    double referenceWindingWindowArea = referenceCore.get_winding_windows()[0].get_area().value();
    add_scored_value("Reference", CoreCrossReferencerFilters::WINDING_WINDOW_AREA, referenceWindingWindowArea);

    std::list<size_t> listOfIndexesToErase;
    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        Core core = (*unfilteredCores)[coreIndex].first;

        if ((*_validScorings).contains(CoreCrossReferencerFilters::WINDING_WINDOW_AREA)) {
            if ((*_validScorings)[CoreCrossReferencerFilters::WINDING_WINDOW_AREA].contains(core.get_name().value())) {
                if ((*_validScorings)[CoreCrossReferencerFilters::WINDING_WINDOW_AREA][core.get_name().value()]) {
                    newScoring.push_back((*_scorings)[CoreCrossReferencerFilters::WINDING_WINDOW_AREA][core.get_name().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(coreIndex);
                }
                continue;
            }
        }
        if (!core.get_winding_windows()[0].get_area()) {
            throw std::runtime_error("Winding window is missing area");
        }
        double windingWindowArea = core.get_winding_windows()[0].get_area().value();

        if (fabs(referenceWindingWindowArea - windingWindowArea) / referenceWindingWindowArea < limit) {
            double scoring = fabs(referenceWindingWindowArea - windingWindowArea);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencerFilters::WINDING_WINDOW_AREA, scoring);
            add_scored_value(core.get_name().value(), CoreCrossReferencerFilters::WINDING_WINDOW_AREA, windingWindowArea);
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
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencerFilters::WINDING_WINDOW_AREA]);
    }
    return filteredCoresWithScoring;
}

std::vector<std::pair<Core, double>> CoreCrossReferencer::MagneticCoreFilterEffectiveArea::filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }

    std::vector<double> newScoring;
    std::vector<std::pair<Core, double>> filteredCoresWithScoring;
    if (!referenceCore.get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    double referenceEffectiveArea = referenceCore.get_processed_description()->get_effective_parameters().get_effective_area();
    add_scored_value("Reference", CoreCrossReferencerFilters::EFFECTIVE_AREA, referenceEffectiveArea);

    std::list<size_t> listOfIndexesToErase;
    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        Core core = (*unfilteredCores)[coreIndex].first;

        if ((*_validScorings).contains(CoreCrossReferencerFilters::EFFECTIVE_AREA)) {
            if ((*_validScorings)[CoreCrossReferencerFilters::EFFECTIVE_AREA].contains(core.get_name().value())) {
                if ((*_validScorings)[CoreCrossReferencerFilters::EFFECTIVE_AREA][core.get_name().value()]) {
                    newScoring.push_back((*_scorings)[CoreCrossReferencerFilters::EFFECTIVE_AREA][core.get_name().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(coreIndex);
                }
                continue;
            }
        }
        if (!core.get_processed_description()) {
            throw std::runtime_error("Core is not processed");
        }
        double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();

        if (fabs(referenceEffectiveArea - effectiveArea) / referenceEffectiveArea < limit) {
            double scoring = fabs(referenceEffectiveArea - effectiveArea);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencerFilters::EFFECTIVE_AREA, scoring);
            add_scored_value(core.get_name().value(), CoreCrossReferencerFilters::EFFECTIVE_AREA, effectiveArea);
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
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencerFilters::EFFECTIVE_AREA]);
    }
    return filteredCoresWithScoring;
}

std::vector<std::pair<Core, double>> CoreCrossReferencer::MagneticCoreFilterEnvelopingVolume::filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }

    std::vector<double> newScoring;
    std::vector<std::pair<Core, double>> filteredCoresWithScoring;
    double referenceDepth = referenceCore.get_depth();
    double referenceHeight = referenceCore.get_height();
    double referenceWidth = referenceCore.get_width();
    add_scored_value("Reference", CoreCrossReferencerFilters::ENVELOPING_VOLUME, std::max(referenceDepth, std::max(referenceHeight, referenceWidth)));

    std::list<size_t> listOfIndexesToErase;
    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        Core core = (*unfilteredCores)[coreIndex].first;
        if ((*_validScorings).contains(CoreCrossReferencerFilters::ENVELOPING_VOLUME)) {
            if ((*_validScorings)[CoreCrossReferencerFilters::ENVELOPING_VOLUME].contains(core.get_name().value())) {
                if ((*_validScorings)[CoreCrossReferencerFilters::ENVELOPING_VOLUME][core.get_name().value()]) {
                    newScoring.push_back((*_scorings)[CoreCrossReferencerFilters::ENVELOPING_VOLUME][core.get_name().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(coreIndex);
                }
                continue;
            }
        }
        double depth = core.get_depth();
        double height = core.get_height();
        double width = core.get_width();

        if (fabs(referenceDepth - depth) / referenceDepth < limit && fabs(referenceHeight - height) / referenceHeight < limit && fabs(referenceWidth - width) / referenceWidth < limit) {
            double scoring = fabs(referenceDepth - depth) + fabs(referenceHeight - height) + fabs(referenceWidth - width);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencerFilters::ENVELOPING_VOLUME, scoring);
            add_scored_value(core.get_name().value(), CoreCrossReferencerFilters::ENVELOPING_VOLUME, std::max(depth, std::max(height, width)));
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
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencerFilters::ENVELOPING_VOLUME]);
    }
    return filteredCoresWithScoring;
}

std::pair<double, double> CoreCrossReferencer::MagneticCoreFilterCoreLosses::calculate_average_core_losses_and_magnetic_flux_density(Core core, int64_t numberTurns, Inputs inputs, std::map<std::string, std::string> models) {
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
        from_json(models["reluctance"], reluctanceModelName);
    }
    auto coreLossesModelName = defaults.coreLossesModelDefault;
    if (models.find("coreLosses") != models.end()) {
        from_json(models["coreLosses"], coreLossesModelName);
    }
    auto coreTemperatureModelName = defaults.coreTemperatureModelDefault;
    if (models.find("coreTemperature") != models.end()) {
        from_json(models["coreTemperature"], coreTemperatureModelName);
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
                break;
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
            Winding winding;
            winding.set_number_turns(numberTurns);
            winding.set_wire("Dummy");
            coil.set_functional_description({winding});
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
                    magneticFluxDensity.set_waveform(Inputs::create_waveform(magneticFluxDensityProcessed, frequency));
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

std::vector<std::pair<Core, double>> CoreCrossReferencer::MagneticCoreFilterCoreLosses::filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, int64_t referenceNumberTurns, Inputs inputs, std::map<std::string, std::string> models, double weight, double limit) {
    if (weight <= 0) {
        return *unfilteredCores;
    }
    std::vector<double> newScoring;

    auto temperature = inputs.get_maximum_temperature();
    auto [referenceCoreLossesWithTemperature, referenceCoreMagneticFluxDensitySaturationPeak] = calculate_average_core_losses_and_magnetic_flux_density(referenceCore, referenceNumberTurns, inputs, models);
    auto referenceMaterialMagneticFluxDensitySaturationPeak = referenceCore.get_magnetic_flux_density_saturation(temperature, true);
    add_scored_value("Reference", CoreCrossReferencerFilters::CORE_LOSSES, referenceCoreLossesWithTemperature);
    add_scored_value("Reference", CoreCrossReferencerFilters::SATURATION, 1 + (referenceCoreMagneticFluxDensitySaturationPeak - referenceMaterialMagneticFluxDensitySaturationPeak) / referenceMaterialMagneticFluxDensitySaturationPeak);

    OperatingPointExcitation excitation;
    SignalDescriptor magneticFluxDensity;
    Processed magneticFluxDensityProcessed;
    magneticFluxDensityProcessed.set_label(WaveformLabel::SINUSOIDAL);
    magneticFluxDensityProcessed.set_offset(0);
    magneticFluxDensityProcessed.set_duty_cycle(0.5);

    std::vector<std::pair<Core, double>> filteredCoresWithScoring;
    std::list<size_t> listOfIndexesToErase;

    for (size_t coreIndex = 0; coreIndex < (*unfilteredCores).size(); ++coreIndex){
        Core core = (*unfilteredCores)[coreIndex].first;
        if ((*_validScorings).contains(CoreCrossReferencerFilters::CORE_LOSSES)) {
            if ((*_validScorings)[CoreCrossReferencerFilters::CORE_LOSSES].contains(core.get_name().value())) {
                if ((*_validScorings)[CoreCrossReferencerFilters::CORE_LOSSES][core.get_name().value()]) {
                    newScoring.push_back((*_scorings)[CoreCrossReferencerFilters::CORE_LOSSES][core.get_name().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(coreIndex);
                }
                continue;
            }
        }

        auto materialMagneticFluxDensitySaturationPeak = core.get_magnetic_flux_density_saturation(temperature, true);

        auto [coreLossesWithTemperature, coreMagneticFluxDensitySaturationPeak] = calculate_average_core_losses_and_magnetic_flux_density(core, referenceNumberTurns, inputs, models);
        add_scored_value(core.get_name().value(), CoreCrossReferencerFilters::CORE_LOSSES, coreLossesWithTemperature);
        add_scored_value(core.get_name().value(), CoreCrossReferencerFilters::SATURATION, 1 + (coreMagneticFluxDensitySaturationPeak - materialMagneticFluxDensitySaturationPeak) / materialMagneticFluxDensitySaturationPeak);

        if ((coreMagneticFluxDensitySaturationPeak < materialMagneticFluxDensitySaturationPeak) && (coreLossesWithTemperature < referenceCoreLossesWithTemperature)) {
            double scoring = 0;
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencerFilters::CORE_LOSSES, scoring);
        }
        else if ((coreMagneticFluxDensitySaturationPeak < materialMagneticFluxDensitySaturationPeak) && ((fabs(referenceCoreLossesWithTemperature - coreLossesWithTemperature) / referenceCoreLossesWithTemperature < limit) || limit >= 1)) {
            double scoring = fabs(referenceCoreLossesWithTemperature - coreLossesWithTemperature);
            newScoring.push_back(scoring);
            add_scoring(core.get_name().value(), CoreCrossReferencerFilters::CORE_LOSSES, scoring);
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
        normalize_scoring(&filteredCoresWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreCrossReferencerFilters::CORE_LOSSES]);
    }
    return filteredCoresWithScoring;
}

std::vector<std::pair<Core, double>> CoreCrossReferencer::get_cross_referenced_core(Core referenceCore, int64_t referenceNumberTurns, Inputs inputs, size_t maximumNumberResults) {
    return get_cross_referenced_core(referenceCore, referenceNumberTurns, inputs, _weights, maximumNumberResults);
}

std::vector<std::pair<Core, double>> CoreCrossReferencer::get_cross_referenced_core(Core referenceCore, int64_t referenceNumberTurns, Inputs inputs, std::map<CoreCrossReferencerFilters, double> weights, size_t maximumNumberResults) {
    _weights = weights;

    if (coreDatabase.empty()) {
        load_cores();
    }

    std::vector<std::pair<Core, double>> cores;
    std::string referenceShapeName = referenceCore.get_shape_name();
    std::string referenceMaterialName = referenceCore.get_material_name();

    if (!referenceCore.get_name()) {
        referenceCore.set_name("Custom");
    }

    MaximumDimensions maximumDimensions;
    bool useMaximumDimensions = false;
    if (inputs.get_design_requirements().get_maximum_dimensions()) {
        maximumDimensions = inputs.get_design_requirements().get_maximum_dimensions().value();
        useMaximumDimensions = true;
    }

    for (auto core : coreDatabase){
        if (referenceShapeName != core.get_shape_name() || referenceMaterialName != core.get_material_name()) {
            if (!_onlyManufacturer || core.get_manufacturer_info()->get_name() == _onlyManufacturer.value()) {
                if (!_onlyReferenceMaterial || referenceMaterialName == core.get_material_name()) {
                    if (!core.get_processed_description()) {
                        core.process_data();
                        core.process_gap();
                    }
                    if (!useMaximumDimensions || core.fits(maximumDimensions, false)) {
                        cores.push_back({core, 0.0});
                    }
                }
            }
        }
    }

    if (!referenceCore.get_processed_description()) {
        referenceCore.process_data();
        referenceCore.process_gap();
    }

    double limit = 0;
    std::vector<std::pair<Core, double>> filteredCores;

    while (limit <= _limit && filteredCores.size() < maximumNumberResults) {
        if (limit < 1) {
            limit += 0.25;
        }
        else if (limit < 10) {
            limit += 2.5;
        }
        else {
            limit += 25;
        }
        filteredCores = apply_filters(&cores, referenceCore, referenceNumberTurns, inputs, weights, maximumNumberResults, limit);
    }

    return filteredCores;
}

std::vector<std::pair<Core, double>> CoreCrossReferencer::apply_filters(std::vector<std::pair<Core, double>>* cores, Core referenceCore, int64_t referenceNumberTurns, Inputs inputs, std::map<CoreCrossReferencerFilters, double> weights, size_t maximumNumberResults, double limit) {
    MagneticCoreFilterPermeance filterPermeance;
    MagneticCoreFilterCoreLosses filterVolumetricLosses;
    MagneticCoreFilterWindingWindowArea filterEffectiveArea;
    MagneticCoreFilterEffectiveArea filterWindingWindowArea;
    MagneticCoreFilterEnvelopingVolume filterEnvelopingVolume;

    filterPermeance.set_scorings(&_scorings);
    filterPermeance.set_valid_scorings(&_validScorings);
    filterPermeance.set_scored_value(&_scoredValues);
    filterPermeance.set_filter_configuration(&_filterConfiguration);
    filterVolumetricLosses.set_scorings(&_scorings);
    filterVolumetricLosses.set_valid_scorings(&_validScorings);
    filterVolumetricLosses.set_scored_value(&_scoredValues);
    filterVolumetricLosses.set_filter_configuration(&_filterConfiguration);
    filterWindingWindowArea.set_scorings(&_scorings);
    filterWindingWindowArea.set_valid_scorings(&_validScorings);
    filterWindingWindowArea.set_scored_value(&_scoredValues);
    filterWindingWindowArea.set_filter_configuration(&_filterConfiguration);
    filterEffectiveArea.set_scorings(&_scorings);
    filterEffectiveArea.set_valid_scorings(&_validScorings);
    filterEffectiveArea.set_scored_value(&_scoredValues);
    filterEffectiveArea.set_filter_configuration(&_filterConfiguration);
    filterEnvelopingVolume.set_scorings(&_scorings);
    filterEnvelopingVolume.set_valid_scorings(&_validScorings);
    filterEnvelopingVolume.set_scored_value(&_scoredValues);
    filterEnvelopingVolume.set_filter_configuration(&_filterConfiguration);

    std::vector<std::pair<Core, double>> rankedCores = *cores;

    magic_enum::enum_for_each<CoreCrossReferencerFilters>([&] (auto val) {
        CoreCrossReferencerFilters filter = val;
        switch (filter) {
            case CoreCrossReferencerFilters::ENVELOPING_VOLUME: 
                rankedCores = filterEnvelopingVolume.filter_core(&rankedCores, referenceCore, weights[CoreCrossReferencerFilters::ENVELOPING_VOLUME], limit);
                break;
            case CoreCrossReferencerFilters::WINDING_WINDOW_AREA: 
                rankedCores = filterWindingWindowArea.filter_core(&rankedCores, referenceCore, weights[CoreCrossReferencerFilters::WINDING_WINDOW_AREA], limit);
                break;
            case CoreCrossReferencerFilters::EFFECTIVE_AREA: 
                rankedCores = filterEffectiveArea.filter_core(&rankedCores, referenceCore, weights[CoreCrossReferencerFilters::EFFECTIVE_AREA], limit);
                break;
            case CoreCrossReferencerFilters::PERMEANCE: 
                rankedCores = filterPermeance.filter_core(&rankedCores, referenceCore, inputs, _models, weights[CoreCrossReferencerFilters::PERMEANCE], limit);
                break;
        }    
        std::string filterString = to_string(filter);
        logEntry("There are " + std::to_string(rankedCores.size()) + " after filtering by " + filterString + ".", "Core Cross Referencer", 2);
    });

    if (rankedCores.size() > (1.1 * maximumNumberResults)) {
        rankedCores = std::vector<std::pair<Core, double>>(rankedCores.begin(), rankedCores.end() - (rankedCores.size() - (1.1 * maximumNumberResults)));
    }

    // We leave core losses for the last one, as it is the most computationally costly
    rankedCores = filterVolumetricLosses.filter_core(&rankedCores, referenceCore, referenceNumberTurns, inputs, _models, weights[CoreCrossReferencerFilters::CORE_LOSSES], limit);
        logEntry("There are " + std::to_string(rankedCores.size()) + " after filtering by " + to_string(CoreCrossReferencerFilters::CORE_LOSSES) + ".", "Core Cross Referencer", 2);

    if (rankedCores.size() > maximumNumberResults) {
        rankedCores = std::vector<std::pair<Core, double>>(rankedCores.begin(), rankedCores.end() - (rankedCores.size() - maximumNumberResults));
    }

    return rankedCores;
}

} // namespace OpenMagnetics
