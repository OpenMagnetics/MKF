#include "BobbinWrapper.h"
#include "CoreLosses.h"
#include "CoreTemperature.h"
#include "CoreAdviser.h"
#include "MagnetizingInductance.h"
#include "WindingSkinEffectLosses.h"
#include "WindingOhmicLosses.h"
#include "NumberTurns.h"
#include "MagneticEnergy.h"
#include "WireWrapper.h"
#include "BobbinWrapper.h"
#include "Defaults.h"
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


namespace OpenMagnetics {


std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> CoreAdviser::get_scorings(bool weighted){
    std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> swappedScorings;
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

void normalize_scoring(std::vector<std::pair<MasWrapper, double>>* masMagneticsWithScoring, std::vector<double>* newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    double maximumScoring = *std::max_element(newScoring->begin(), newScoring->end());
    double minimumScoring = *std::min_element(newScoring->begin(), newScoring->end());

    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i) {
        auto mas = (*masMagneticsWithScoring)[i].first;
        auto scoring = (*newScoring)[i];
        if (maximumScoring != minimumScoring) {

            if (filterConfiguration["log"]){
                if (filterConfiguration["invert"]) {
                    (*masMagneticsWithScoring)[i].second += weight * (1 - (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring)));
                }
                else {
                    (*masMagneticsWithScoring)[i].second += weight * (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                }
            }
            else {
                if (filterConfiguration["invert"]) {
                    (*masMagneticsWithScoring)[i].second += weight * (1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring));
                }
                else {
                    (*masMagneticsWithScoring)[i].second += weight * (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                }
            }
        }
        else {
            (*masMagneticsWithScoring)[i].second += 1;
        }
    }
    sort((*masMagneticsWithScoring).begin(), (*masMagneticsWithScoring).end(), [](std::pair<MasWrapper, double>& b1, std::pair<MasWrapper, double>& b2) {
        return b1.second > b2.second;
    });
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterAreaProduct::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    auto defaults = Defaults();
    std::vector<std::pair<MasWrapper, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::map<std::string, double> materialScaledMagneticFluxDensities;
    std::map<std::string, double> bobbinFillingFactors;

    double primaryAreaFactor = 1;
    if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
        primaryAreaFactor = 0.5;
    }

    std::vector<double> areaProductRequiredPreCalculations;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        auto voltageWaveform = InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_voltage().value().get_waveform().value();
        auto currentWaveform = InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_current().value().get_waveform().value();
        double frequency = InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency();
        if (voltageWaveform.get_data().size() != currentWaveform.get_data().size()) {
            voltageWaveform = InputsWrapper::calculate_sampled_waveform(voltageWaveform, frequency);
            currentWaveform = InputsWrapper::calculate_sampled_waveform(currentWaveform, frequency);
        }

        std::vector<double> voltageWaveformData = voltageWaveform.get_data();
        std::vector<double> currentWaveformData = currentWaveform.get_data();

        double powerMean = 0;
        for (size_t i = 0; i < voltageWaveformData.size(); ++i)
        {
            powerMean += fabs(voltageWaveformData[i] * currentWaveformData[i]);
        }
        powerMean /= voltageWaveformData.size();
        areaProductRequiredPreCalculations.push_back(powerMean / (primaryAreaFactor * 2 * frequency * defaults.maximumCurrentDensity));
    }

    std::map<std::string, double> scaledMagneticFluxDensitiesPerMaterial;
    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));
    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLosses();

    double magneticFluxDensityReference = 0.18;
    double frequencyReference = 100000;
    OperatingPointExcitation operatingPointExcitation;
    SignalDescriptor magneticFluxDensity;
    Processed processed;
    operatingPointExcitation.set_frequency(frequencyReference);
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);
    processed.set_peak(magneticFluxDensityReference);
    processed.set_peak_to_peak(2 * magneticFluxDensityReference);
    magneticFluxDensity.set_processed(processed);
    operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);

    std::list<size_t> listOfIndexesToErase;

    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = magnetic.get_core();

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::AREA_PRODUCT)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT][magnetic.get_manufacturer_info().value().get_reference().value()]) {

                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }

        double bobbinFillingFactor;
        if (core.get_winding_windows().size() == 0)
            continue;
        auto windingWindow = core.get_winding_windows()[0];
        auto windingColumn = core.get_columns()[0];
        if (!bobbinFillingFactors.contains(core.get_shape_name())) {
            if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
                bobbinFillingFactor = OpenMagnetics::BobbinWrapper::get_filling_factor(windingWindow.get_width().value(), core.get_winding_windows()[0].get_height().value());
            }
            else {
                bobbinFillingFactor = 1;
            }
            bobbinFillingFactors[core.get_shape_name()] = bobbinFillingFactor;
        }
        else {
            bobbinFillingFactor = bobbinFillingFactors[core.get_shape_name()];
        }
        double areaProductCore = windingWindow.get_area().value() * windingColumn.get_area();
        double maximumAreaProductRequired = 0;

        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            double temperature = inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature();
            double frequency = InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency();
            auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
            double wireAirFillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(2 * skinDepth);
            double windingWindowUtilizationFactor = wireAirFillingFactor * bobbinFillingFactor;
            double magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
            if (!materialScaledMagneticFluxDensities.contains(core.get_material_name())) {
                auto coreLossesMethods = core.get_available_core_losses_methods();
                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
                    double referenceCoreLosses = coreLossesModelSteinmetz->get_core_losses(core, operatingPointExcitation, temperature).get_core_losses();
                    auto aux = coreLossesModelSteinmetz->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                    magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
                }
                else {
                    double referenceCoreLosses = coreLossesModelProprietary->get_core_losses(core, operatingPointExcitation, temperature).get_core_losses();
                    auto aux = coreLossesModelProprietary->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                    magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
                }
                materialScaledMagneticFluxDensities[core.get_material_name()] = magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
            }
            else {
                magneticFluxDensityPeakAtFrequencyOfReferenceLosses = materialScaledMagneticFluxDensities[core.get_material_name()];
            }

            double areaProductRequired = areaProductRequiredPreCalculations[operatingPointIndex] / (windingWindowUtilizationFactor * magneticFluxDensityPeakAtFrequencyOfReferenceLosses);
            if (std::isnan(magneticFluxDensityPeakAtFrequencyOfReferenceLosses) || magneticFluxDensityPeakAtFrequencyOfReferenceLosses == 0) {
                throw std::runtime_error("magneticFluxDensityPeakAtFrequencyOfReferenceLosses cannot be 0 or NaN");
                break;
            }
            if (std::isnan(areaProductRequired)) {
                break;
            }
            if (std::isinf(areaProductRequired) || areaProductRequired == 0) {
                throw std::runtime_error("areaProductRequired cannot be 0 or NaN");
            }

            maximumAreaProductRequired = std::max(maximumAreaProductRequired, areaProductRequired);
        }

        if (areaProductCore >= maximumAreaProductRequired * defaults.coreAdviserThresholdValidity) {
            double scoring = fabs(areaProductCore - maximumAreaProductRequired);
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::AREA_PRODUCT, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT]);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterEnergyStored::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    auto defaults = Defaults();
    OpenMagnetics::MagneticEnergy magneticEnergy(models);    
    std::vector<std::pair<MasWrapper, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    double requiredMagneticEnergy = magneticEnergy.required_magnetic_energy(inputs).get_nominal().value();
    MagnetizingInductanceOutput magnetizingInductanceOutput;

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){  
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::ENERGY_STORED)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }

        bool validMagnetic = true;
        double totalStorableMagneticEnergy = 0;
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            auto operatingPoint = inputs.get_operating_point(operatingPointIndex);
            totalStorableMagneticEnergy = std::max(totalStorableMagneticEnergy, magneticEnergy.get_core_maximum_magnetic_energy(static_cast<CoreWrapper>(magnetic.get_core()), &operatingPoint));
            if (totalStorableMagneticEnergy >= requiredMagneticEnergy * defaults.coreAdviserThresholdValidity) {
                magnetizingInductanceOutput.set_maximum_magnetic_energy_core(totalStorableMagneticEnergy);
                magnetizingInductanceOutput.set_method_used(models["gapReluctance"]);
                magnetizingInductanceOutput.set_origin(ResultOrigin::SIMULATION);
                mas.get_mutable_outputs()[operatingPointIndex].set_magnetizing_inductance(magnetizingInductanceOutput);
            }
            else {
                validMagnetic = false;
                break;
            }
        }
        if (validMagnetic) {
            double scoring = totalStorableMagneticEnergy;
            newScoring.push_back(scoring);
            (*unfilteredMasMagnetics)[masIndex].first = mas;
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::ENERGY_STORED, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }
    }


    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED]);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterWindingWindowArea::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    auto defaults = Defaults();
    std::vector<std::pair<MasWrapper, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    double primaryCurrentRms = 0;
    double frequency = 0;
    double temperature = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        primaryCurrentRms = std::max(primaryCurrentRms, InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_current().value().get_processed().value().get_rms().value());
        frequency = std::max(frequency, InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency());
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }

    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLosses();
    auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
    double wireAirFillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(2 * skinDepth);
    double estimatedWireConductingArea = std::numbers::pi * pow(skinDepth, 2);
    double estimatedWireTotalArea = estimatedWireConductingArea / wireAirFillingFactor;
    double necessaryWireCopperArea = primaryCurrentRms / defaults.maximumCurrentDensity;
    double estimatedParallels = ceil(necessaryWireCopperArea / estimatedWireConductingArea);

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = magnetic.get_core();

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }

        double primaryNumberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();
        double estimatedNeededWindingArea = primaryNumberTurns * estimatedParallels * estimatedWireTotalArea;
        double coreWindingWindowArea = core.get_winding_windows()[0].get_area().value();

        if (coreWindingWindowArea >= estimatedNeededWindingArea * defaults.coreAdviserThresholdValidity) {
            double manufacturabilityRelativeCost;
            if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
                double estimatedNeededLayers = (primaryNumberTurns * estimatedParallels * (2 * skinDepth / wireAirFillingFactor)) / core.get_winding_windows()[0].get_height().value();
                manufacturabilityRelativeCost = estimatedNeededLayers;
            }
            else {
                double layerLength = 2 * std::numbers::pi * (core.get_winding_windows()[0].get_radial_height().value() - skinDepth);
                double estimatedNeededLayers = (primaryNumberTurns * estimatedParallels * (2 * skinDepth / wireAirFillingFactor)) / layerLength;
                if (estimatedNeededLayers > 1) {
                    manufacturabilityRelativeCost = estimatedNeededLayers * 2;
                }
                else {
                    manufacturabilityRelativeCost = estimatedNeededLayers;
                }
            }
            if (core.get_functional_description().get_number_stacks() > 1) {
                manufacturabilityRelativeCost *= 2;  // Because we need a custom bobbin
            }

            double scoring = manufacturabilityRelativeCost;
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }
    }


    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA]);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterCoreLosses::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    auto defaults = Defaults();
    std::vector<std::pair<MasWrapper, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;


    std::vector<double> powerMeans(inputs.get_operating_points().size(), 0);
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        auto voltageWaveform = InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_voltage().value().get_waveform().value();
        auto currentWaveform = InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_current().value().get_waveform().value();
        double frequency = InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency();

        if (voltageWaveform.get_data().size() != currentWaveform.get_data().size()) {
            voltageWaveform = InputsWrapper::calculate_sampled_waveform(voltageWaveform, frequency);
            currentWaveform = InputsWrapper::calculate_sampled_waveform(currentWaveform, frequency);
        }
        std::vector<double> voltageWaveformData = voltageWaveform.get_data();
        std::vector<double> currentWaveformData = currentWaveform.get_data();
        for (size_t i = 0; i < voltageWaveformData.size(); ++i)
        {
            powerMeans[operatingPointIndex] += fabs(voltageWaveformData[i] * currentWaveformData[i]);
        }
        powerMeans[operatingPointIndex] /= voltageWaveformData.size();
    }


    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(models);
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));

    OpenMagnetics::MagnetizingInductance magnetizingInductance(models["gapReluctance"]);
    auto windingOhmicLosses = OpenMagnetics::WindingOhmicLosses();
    size_t numberTimeouts = 0;

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = magnetic.get_core();

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::CORE_LOSSES)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::CORE_LOSSES].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::CORE_LOSSES][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::CORE_LOSSES][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }

        std::string shapeName = core.get_shape_name();
        if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("R ", 0) == 0) || (shapeName.rfind("T ", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
            auto bobbin = OpenMagnetics::BobbinWrapper::create_quick_bobbin(core);
            magnetic.get_mutable_coil().set_bobbin(bobbin);
            auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();

            if ((windingWindows[0].get_width().value() < 0) || (windingWindows[0].get_width().value() > 1)) {
                throw std::runtime_error("Something wrong happened in bobbins 1:   windingWindows[0].get_width(): " + std::to_string(bool(windingWindows[0].get_width())) +
                                         " windingWindows[0].get_width().value(): " + std::to_string(windingWindows[0].get_width().value()) + 
                                         " shapeName: " + shapeName
                                         );
            }
        }

        auto currentNumberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();
        OpenMagnetics::NumberTurns numberTurns(currentNumberTurns);
        std::vector<double> totalLossesPerOperatingPoint;
        std::vector<CoreLossesOutput> coreLossesPerOperatingPoint;
        std::vector<WindingLossesOutput> windingLossesPerOperatingPoint;
        double currentTotalLosses = DBL_MAX;
        double coreLosses = DBL_MAX;
        CoreLossesOutput coreLossesOutput;
        double ohmicLosses = DBL_MAX;
        WindingLossesOutput windingLossesOutput;
        double newTotalLosses = DBL_MAX;
        auto previousNumberTurnsPrimary = currentNumberTurns;

        size_t iteration = 10;

        CoilWrapper coil = magnetic.get_coil();

        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            auto operatingPoint = inputs.get_operating_point(operatingPointIndex);
            double temperature = operatingPoint.get_conditions().get_ambient_temperature();
            OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
            do {
                currentTotalLosses = newTotalLosses;
                auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
                coil.get_mutable_functional_description()[0].set_number_turns(numberTurnsCombination[0]);
                // coil = CoilWrapper(coil);
                coil.try_wind(true);

                auto aux = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint);
                auto magnetizingInductance = aux.first;
                auto magneticFluxDensity = aux.second;

                if (!check_requirement(inputs.get_design_requirements().get_magnetizing_inductance(), magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                    coil.get_mutable_functional_description()[0].set_number_turns(previousNumberTurnsPrimary);
                    // coil = CoilWrapper(coil);
                    coil.try_wind(true);
                    break;
                }
                else {
                    previousNumberTurnsPrimary = numberTurnsCombination[0];
                }

                if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("R ", 0) == 0) || (shapeName.rfind("T ", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
                    if (!coil.get_turns_description()) {
                        newTotalLosses = coreLosses;
                        break;
                    }
                }

                excitation.set_magnetic_flux_density(magneticFluxDensity);
                auto coreLossesMethods = core.get_available_core_losses_methods();
                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
                    coreLossesOutput = coreLossesModelSteinmetz->get_core_losses(core, excitation, temperature);
                    coreLosses = coreLossesOutput.get_core_losses();
                }
                else {
                    coreLossesOutput = coreLossesModelProprietary->get_core_losses(core, excitation, temperature);
                    coreLosses = coreLossesOutput.get_core_losses();
                    if (coreLosses < 0) {
                        break;
                    }
                }

                if (coreLosses < 0) {
                    throw std::runtime_error("Something wrong happend in core losses calculation for magnetic: " + magnetic.get_manufacturer_info().value().get_reference().value());
                }

                if (!coil.get_turns_description()) {
                    newTotalLosses = coreLosses;
                    break;
                }

                if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("R ", 0) == 0) || (shapeName.rfind("T ", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
                    windingLossesOutput = windingOhmicLosses.calculate_ohmic_losses(coil, operatingPoint, temperature);
                    ohmicLosses = windingLossesOutput.get_winding_losses();
                    newTotalLosses = coreLosses + ohmicLosses;
                    if (ohmicLosses < 0) {
                        throw std::runtime_error("Something wrong happend in ohmic losses calculation for magnetic: " + magnetic.get_manufacturer_info().value().get_reference().value() + " ohmicLosses: " + std::to_string(ohmicLosses));
                    }

                }
                else {
                    newTotalLosses = coreLosses;
                    break;
                }

                if (newTotalLosses == DBL_MAX) {
                    throw std::runtime_error("Too large losses");
                }

                iteration--;
                if (iteration <=0) {
                    numberTimeouts++;
                    break;
                }
            }
            while(newTotalLosses < currentTotalLosses * defaults.coreAdviserThresholdValidity);
 

            if (coreLosses < DBL_MAX && coreLosses > 0) {
                magnetic.set_coil(coil);

                currentTotalLosses = newTotalLosses;
                totalLossesPerOperatingPoint.push_back(currentTotalLosses);
                coreLossesPerOperatingPoint.push_back(coreLossesOutput);
                windingLossesPerOperatingPoint.push_back(windingLossesOutput);
            }
        }


        if (totalLossesPerOperatingPoint.size() < inputs.get_operating_points().size()) {
            listOfIndexesToErase.push_back(masIndex);
        }
        else {

            double meanTotalLosses = 0;
            for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
                meanTotalLosses += totalLossesPerOperatingPoint[operatingPointIndex];
            }
            if (meanTotalLosses > DBL_MAX / 2) {
                throw std::runtime_error("Something wrong happend in core losses calculation for magnetic: " + magnetic.get_manufacturer_info().value().get_reference().value());
            }
            meanTotalLosses /= inputs.get_operating_points().size();
            double maximumPowerMean = *max_element(powerMeans.begin(), powerMeans.end());

            if (meanTotalLosses < maximumPowerMean * defaults.coreAdviserMaximumPercentagePowerCoreLosses / defaults.coreAdviserThresholdValidity) {
                for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
                    mas.get_mutable_outputs()[operatingPointIndex].set_core_losses(coreLossesPerOperatingPoint[operatingPointIndex]);
                    mas.get_mutable_outputs()[operatingPointIndex].set_winding_losses(windingLossesPerOperatingPoint[operatingPointIndex]);
                }
                double scoring = meanTotalLosses;
                newScoring.push_back(scoring);
                (*unfilteredMasMagnetics)[masIndex].first = mas;
                if (magnetic.get_manufacturer_info().value().get_reference().value() == "P 9/5 - N30 - Ungapped") {
                }
                add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::CORE_LOSSES, scoring, firstFilter);
            }
            else {
                listOfIndexesToErase.push_back(masIndex);
            }
        }
    }



    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::CORE_LOSSES]);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterCoreTemperature::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    auto defaults = Defaults();
    std::vector<std::pair<MasWrapper, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));

    OpenMagnetics::MagnetizingInductance magnetizingInductance(models["gapReluctance"]);

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = magnetic.get_core();

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }
        auto winding = magnetic.get_coil();
        bool validMagnetic = true;
        double meanTemperature = 0;
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            double temperature = inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature();
            auto operatingPoint = inputs.get_operating_point(operatingPointIndex);
            OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
            double coreLosses;
            if (!mas.get_outputs()[operatingPointIndex].get_core_losses()) {
                auto magneticFluxDensity = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;
                excitation.set_magnetic_flux_density(magneticFluxDensity);
                auto coreLossesMethods = core.get_available_core_losses_methods();
                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
                    coreLosses = coreLossesModelSteinmetz->get_core_losses(core, excitation, temperature).get_core_losses();
                }
                else {
                    coreLosses = coreLossesModelProprietary->get_core_losses(core, excitation, temperature).get_core_losses();
                }
                CoreLossesOutput coreLossesOutput;
                coreLossesOutput.set_core_losses(coreLosses);
                coreLossesOutput.set_method_used(models["coreLosses"]);
                coreLossesOutput.set_origin(ResultOrigin::SIMULATION);
                mas.get_mutable_outputs()[operatingPointIndex].set_core_losses(coreLossesOutput);
            }
            else {
                coreLosses = mas.get_outputs()[operatingPointIndex].get_core_losses().value().get_core_losses();
            }

            auto coreTemperatureModel = OpenMagnetics::CoreTemperatureModel::factory(CoreTemperatureModels::MANIKTALA);
            auto coreTemperature = coreTemperatureModel->get_core_temperature(core, coreLosses, temperature);
            double calculatedTemperature = coreTemperature.get_maximum_temperature();
            meanTemperature += calculatedTemperature;

            if (calculatedTemperature < defaults.coreAdviserMaximumCoreTemperature / defaults.coreAdviserThresholdValidity) {
                auto coreLossesOutput = mas.get_mutable_outputs()[operatingPointIndex].get_core_losses().value();
                coreLossesOutput.set_temperature(calculatedTemperature);
                mas.get_mutable_outputs()[operatingPointIndex].set_core_losses(coreLossesOutput);
            }
            else {
                validMagnetic = false;
                break;
            }
        }
        if (validMagnetic) {
            meanTemperature /= inputs.get_operating_points().size();
            double scoring = meanTemperature;
            newScoring.push_back(scoring);
            (*unfilteredMasMagnetics)[masIndex].first = mas;
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }
    }


    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE]);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterDimensions::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    std::vector<double> newScoring;

    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = magnetic.get_core();
        double volume = core.get_width() * core.get_height() * core.get_depth();

        double scoring = volume;
        newScoring.push_back(scoring);
        add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::DIMENSIONS, scoring, firstFilter);
    }

    if ((*unfilteredMasMagnetics).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string((*unfilteredMasMagnetics).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredMasMagnetics).size() > 0) {
        normalize_scoring(unfilteredMasMagnetics, &newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::DIMENSIONS]);
    }
    return (*unfilteredMasMagnetics);
}


CoilWrapper get_dummy_coil(InputsWrapper inputs) {
    double frequency = 0; 
    double temperature = 0; 
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        frequency = std::max(frequency, InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency());
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }
    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLosses();
    auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);  // TODO material hardcoded

    // Set round wire with diameter to two times the skin depth 
    WireWrapper wire;
    wire.set_nominal_value_conducting_diameter(skinDepth * 2);
    wire.set_nominal_value_outer_diameter(skinDepth * 2.05); // Hardcoded
    wire.set_material("copper");
    wire.set_type(WireType::ROUND);
    CoilFunctionalDescription primaryCoilFunctionalDescription;
    primaryCoilFunctionalDescription.set_isolation_side(IsolationSide::PRIMARY);
    primaryCoilFunctionalDescription.set_name("primary");
    primaryCoilFunctionalDescription.set_number_parallels(1);
    primaryCoilFunctionalDescription.set_number_turns(1);
    primaryCoilFunctionalDescription.set_wire(wire);

    CoilWrapper coil;
    coil.set_bobbin("Dummy");
    coil.set_functional_description({primaryCoilFunctionalDescription});
    return coil;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::get_advised_core(InputsWrapper inputs, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::get_advised_core(InputsWrapper inputs, std::vector<CoreWrapper>* cores, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults);
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults){
    std::string file_path = __FILE__;
    auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/../../MAS/data/cores_stock.ndjson");
    std::ifstream ndjsonFile(inventory_path);
    std::string jsonLine;
    std::vector<CoreWrapper> cores;
    while (std::getline(ndjsonFile, jsonLine)) {
        json jf = json::parse(jsonLine);
        CoreWrapper core(jf, false, true, false);
        if (_includeToroids || core.get_type() != CoreType::TOROIDAL) {
            cores.push_back(core);
        }
    }
    return get_advised_core(inputs, weights, &cores, maximumNumberResults);
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::get_advised_core(InputsWrapper inputs, std::vector<CoreWrapper>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults, maximumNumberCores);
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, std::vector<CoreWrapper>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {

    std::vector<std::pair<MasWrapper, double>> results;


    for(size_t i = 0; i < (*cores).size(); i += maximumNumberCores) {
        auto last = std::min((*cores).size(), i + maximumNumberCores);
        std::vector<CoreWrapper> partialCores;
        partialCores.reserve(last - i);

        move((*cores).begin() + i, (*cores).begin() + last, std::back_inserter(partialCores));
        
        auto partialResult = get_advised_core(inputs, weights, &partialCores, maximumNumberResults);
        std::move(partialResult.begin(), partialResult.end(), std::back_inserter(results));
    }

    sort(results.begin(), results.end(), [](std::pair<MasWrapper, double>& b1, std::pair<MasWrapper, double>& b2) {
        return b1.second > b2.second;
    });

    if (results.size() > maximumNumberResults) {
        results = std::vector<std::pair<MasWrapper, double>>(results.begin(), results.end() - (results.size() - maximumNumberResults));
    }
    return results;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, std::vector<CoreWrapper>* cores, size_t maximumNumberResults){
    auto defaults = Defaults();
    _weights = weights;
 
    CoreAdviserFilters firstFilter = CoreAdviserFilters::AREA_PRODUCT;
    double maxWeight = 0;

    for (auto& pair : weights) {
        auto filter = pair.first;
        auto weight = pair.second;
        if (weight > maxWeight) {
            maxWeight = weight;
            firstFilter = filter;
        }
    }
    size_t maximumMagneticsAfterFiltering = defaults.coreAdviserMaximumMagneticsAfterFiltering;
    std::vector<std::pair<MasWrapper, double>> masMagnetics;
    bool needToAddStacks = false;

    if (firstFilter == CoreAdviserFilters::CORE_LOSSES || firstFilter == CoreAdviserFilters::CORE_TEMPERATURE) {
        masMagnetics = create_mas_dataset(inputs, cores, true);
        logEntry("We start the search with " + std::to_string(masMagnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.");
        std::string firstFilterString = std::string{magic_enum::enum_name(firstFilter)};
        logEntry("We include stacks of cores in our search because the most important selectd filter is " + firstFilterString + ".");
        auto filteredMasMagnetics = apply_filters(&masMagnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMasMagnetics.size() >= maximumNumberResults) {
            return filteredMasMagnetics;
        }

    }
    else {
        needToAddStacks = true;
        masMagnetics = create_mas_dataset(inputs, cores, false);
        logEntry("We start the search with " + std::to_string(masMagnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.");
        logEntry("We don't include stacks of cores in our search.");
        auto filteredMasMagnetics = apply_filters(&masMagnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMasMagnetics.size() >= maximumNumberResults) {
            return filteredMasMagnetics;
        }
    }

    // masMagnetics = create_mas_dataset(inputs, cores, true);
    if (needToAddStacks) {
        expand_mas_dataset_with_stacks(inputs, cores, &masMagnetics);
    }

    // auto masMagnetics = create_mas_dataset(inputs, cores, true);
    logEntry("First attempt produced not enough results, so now we are searching again with " + std::to_string(masMagnetics.size()) + " magnetics, including up to " + std::to_string(defaults.coreAdviserMaximumNumberStacks) + " cores stacked when possible.");
    weights[CoreAdviserFilters::CORE_LOSSES] = 1;
    maximumMagneticsAfterFiltering = masMagnetics.size();
    auto filteredMasMagnetics = apply_filters(&masMagnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    return filteredMasMagnetics;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::create_mas_dataset(InputsWrapper inputs, std::vector<CoreWrapper>* cores, bool includeStacks) {
    auto defaults = Defaults();
    std::vector<std::pair<MasWrapper, double>> masMagnetics;
    CoilWrapper coil = get_dummy_coil(inputs);

    MagneticWrapper magnetic;
    MasWrapper mas;
    OutputsWrapper outputs;
    for (size_t i = 0; i < inputs.get_operating_points().size(); ++i)
    {
        mas.get_mutable_outputs().push_back(outputs);
    }
    magnetic.set_coil(std::move(coil));
    for (auto& core : *cores){
        if (!_includeToroids && core.get_type() == CoreType::TOROIDAL) {
            continue;
        }

        core.process_data();
        if (!core.process_gap()) {
            continue;
        }

        if (includeStacks && (core.get_shape_family() == CoreShapeFamily::E || core.get_shape_family() == CoreShapeFamily::T || core.get_shape_family() == CoreShapeFamily::U)) {
            for (size_t i = 0; i < defaults.coreAdviserMaximumNumberStacks; ++i)
            {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                core.scale_to_stacks(1 + i);
                // core.process_gap();
                magnetic.set_core(core);
                MagneticManufacturerInfo magneticmanufacturerinfo;
                if (i!=0) {
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stack" );
                }
                magnetic.set_manufacturer_info(magneticmanufacturerinfo);
                mas.set_magnetic(magnetic);
                masMagnetics.push_back(std::pair<MasWrapper, double>{mas, 0});
            }
        }
        else {
            magnetic.set_core(core);
            MagneticManufacturerInfo magneticmanufacturerinfo;
            magneticmanufacturerinfo.set_reference(core.get_name().value());
            magnetic.set_manufacturer_info(magneticmanufacturerinfo);
            mas.set_magnetic(magnetic);
            masMagnetics.push_back(std::pair<MasWrapper, double>{mas, 0});
        }
    }

    return masMagnetics;
}

void CoreAdviser::expand_mas_dataset_with_stacks(InputsWrapper inputs, std::vector<CoreWrapper>* cores, std::vector<std::pair<MasWrapper, double>>* masMagnetics) {
    auto defaults = Defaults();
    CoilWrapper coil = get_dummy_coil(inputs);

    MagneticWrapper magnetic;
    MasWrapper mas;
    OutputsWrapper outputs;
    for (size_t i = 0; i < inputs.get_operating_points().size(); ++i)
    {
        mas.get_mutable_outputs().push_back(outputs);
    }
    magnetic.set_coil(std::move(coil));
    for (auto& core : *cores){
        if (!_includeToroids && core.get_type() == CoreType::TOROIDAL) {
            continue;
        }

        if (core.get_shape_family() == CoreShapeFamily::E || core.get_shape_family() == CoreShapeFamily::T || core.get_shape_family() == CoreShapeFamily::U) {

            core.process_data();
            if (!core.process_gap()) {
                continue;
            }

            for (size_t i = 1; i < defaults.coreAdviserMaximumNumberStacks; ++i)
            {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                core.scale_to_stacks(1 + i);
                // core.process_gap();
                MagneticManufacturerInfo magneticmanufacturerinfo;
                if (i!=0) {
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stack" );
                }
                magnetic.set_manufacturer_info(magneticmanufacturerinfo);
                magnetic.set_core(core);
                mas.set_magnetic(magnetic);
                (*masMagnetics).push_back(std::pair<MasWrapper, double>{mas, 0});
            }
        }

    }
}

void add_initial_turns(std::vector<std::pair<MasWrapper, double>> *masMagneticsWithScoring, InputsWrapper inputs) {
    MagnetizingInductance magnetizingInductance;
    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i){

        CoreWrapper core = (*masMagneticsWithScoring)[i].first.get_magnetic().get_core();
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        double numberTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs, DimensionalValues::MINIMUM);
        (*masMagneticsWithScoring)[i].first.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(numberTurns);
    }
}

void correct_windings(std::vector<std::pair<MasWrapper, double>> *masMagneticsWithScoring, InputsWrapper inputs) {
    MagnetizingInductance magnetizingInductance;
    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i){

        CoilWrapper coil = CoilWrapper((*masMagneticsWithScoring)[i].first.get_magnetic().get_coil());
        OpenMagnetics::NumberTurns numberTurns(coil.get_number_turns(0), inputs.get_design_requirements());
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();

        (*masMagneticsWithScoring)[i].first.set_inputs(inputs);
        if ((*masMagneticsWithScoring)[i].first.get_magnetic().get_core().get_functional_description().get_type() != CoreType::TOROIDAL) {
            (*masMagneticsWithScoring)[i].first.get_mutable_magnetic().get_mutable_coil().set_bobbin(BobbinWrapper::create_quick_bobbin((*masMagneticsWithScoring)[i].first.get_magnetic().get_core()));
        }
        for (size_t windingIndex = 1; windingIndex < numberTurnsCombination.size(); ++windingIndex) {
            auto winding = coil.get_functional_description()[0];
            winding.set_number_turns(numberTurnsCombination[windingIndex]);
            winding.set_isolation_side(get_isolation_side_from_index(windingIndex));
            winding.set_name(std::string{magic_enum::enum_name(get_isolation_side_from_index(windingIndex))});
            (*masMagneticsWithScoring)[i].first.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description().push_back(winding);
        }
    }
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::apply_filters(std::vector<std::pair<MasWrapper, double>>* masMagnetics, InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
    MagneticCoreFilterAreaProduct filterAreaProduct;
    MagneticCoreFilterEnergyStored filterEnergyStored;
    MagneticCoreFilterWindingWindowArea filterWindingWindowArea;
    MagneticCoreFilterCoreLosses filterCoreLosses;
    MagneticCoreFilterCoreTemperature filterCoreTemperature;
    MagneticCoreFilterDimensions filterDimensions;

    filterAreaProduct.set_scorings(&_scorings);
    filterAreaProduct.set_valid_scorings(&_validScorings);
    filterAreaProduct.set_filter_configuration(&_filterConfiguration);
    filterEnergyStored.set_scorings(&_scorings);
    filterEnergyStored.set_valid_scorings(&_validScorings);
    filterEnergyStored.set_filter_configuration(&_filterConfiguration);
    filterWindingWindowArea.set_scorings(&_scorings);
    filterWindingWindowArea.set_valid_scorings(&_validScorings);
    filterWindingWindowArea.set_filter_configuration(&_filterConfiguration);
    filterCoreLosses.set_scorings(&_scorings);
    filterCoreLosses.set_valid_scorings(&_validScorings);
    filterCoreLosses.set_filter_configuration(&_filterConfiguration);
    filterCoreTemperature.set_scorings(&_scorings);
    filterCoreTemperature.set_valid_scorings(&_validScorings);
    filterCoreTemperature.set_filter_configuration(&_filterConfiguration);
    filterDimensions.set_scorings(&_scorings);
    filterDimensions.set_valid_scorings(&_validScorings);
    filterDimensions.set_filter_configuration(&_filterConfiguration);

    CoreAdviserFilters firstFilter = CoreAdviserFilters::AREA_PRODUCT;
    double maxWeight = 0;

    for (auto& pair : weights) {
        auto filter = pair.first;
        auto weight = pair.second;
        if (weight > maxWeight) {
            maxWeight = weight;
            firstFilter = filter;
        }
    }

    std::vector<std::pair<MasWrapper, double>> masMagneticsWithScoring;

    switch (firstFilter) {
        case CoreAdviserFilters::AREA_PRODUCT: 

            masMagneticsWithScoring = filterAreaProduct.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::AREA_PRODUCT], true);
            break;
        case CoreAdviserFilters::ENERGY_STORED: 

            masMagneticsWithScoring = filterEnergyStored.filter_magnetics(masMagnetics, inputs, _models, weights[CoreAdviserFilters::ENERGY_STORED], true);
            break;
        case CoreAdviserFilters::WINDING_WINDOW_AREA: 

            masMagneticsWithScoring = filterWindingWindowArea.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::WINDING_WINDOW_AREA], true);
            break;
        case CoreAdviserFilters::CORE_LOSSES: 

            add_initial_turns(masMagnetics, inputs);
            masMagneticsWithScoring = filterCoreLosses.filter_magnetics(masMagnetics, inputs, _models, weights[CoreAdviserFilters::CORE_LOSSES], true);
            break;
        case CoreAdviserFilters::CORE_TEMPERATURE: 

            add_initial_turns(masMagnetics, inputs);
            masMagneticsWithScoring = filterCoreTemperature.filter_magnetics(masMagnetics, inputs, _models, weights[CoreAdviserFilters::CORE_TEMPERATURE], true);
            break;
        case CoreAdviserFilters::DIMENSIONS: 

            masMagneticsWithScoring = filterDimensions.filter_magnetics(masMagnetics, weights[CoreAdviserFilters::DIMENSIONS], true);
            break;
    }


    std::string firstFilterString = std::string{magic_enum::enum_name(firstFilter)};
    logEntry("There are " + std::to_string(masMagneticsWithScoring.size()) + " magnetics after the first filter, which was " + firstFilterString + ".");

    if (masMagneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
        masMagneticsWithScoring = std::vector<std::pair<MasWrapper, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumMagneticsAfterFiltering));
        logEntry("There are " + std::to_string(masMagneticsWithScoring.size()) + " after culling by the score on the first filter.");
    }

    if (firstFilter != CoreAdviserFilters::CORE_LOSSES && firstFilter != CoreAdviserFilters::CORE_TEMPERATURE) {
        add_initial_turns(&masMagneticsWithScoring, inputs);
    }
    

    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        if (filter != firstFilter) {
            switch (filter) {
                case CoreAdviserFilters::AREA_PRODUCT: 

                    masMagneticsWithScoring = filterAreaProduct.filter_magnetics(&masMagneticsWithScoring, inputs, weights[CoreAdviserFilters::AREA_PRODUCT]);
                    break;
                case CoreAdviserFilters::ENERGY_STORED: 

                    masMagneticsWithScoring = filterEnergyStored.filter_magnetics(&masMagneticsWithScoring, inputs, _models, weights[CoreAdviserFilters::ENERGY_STORED]);
                    break;
                case CoreAdviserFilters::WINDING_WINDOW_AREA: 

                    masMagneticsWithScoring = filterWindingWindowArea.filter_magnetics(&masMagneticsWithScoring, inputs, weights[CoreAdviserFilters::WINDING_WINDOW_AREA]);
                    break;
                case CoreAdviserFilters::CORE_LOSSES: 

                    masMagneticsWithScoring = filterCoreLosses.filter_magnetics(&masMagneticsWithScoring, inputs, _models, weights[CoreAdviserFilters::CORE_LOSSES]);
                    break;
                case CoreAdviserFilters::CORE_TEMPERATURE: 

                    masMagneticsWithScoring = filterCoreTemperature.filter_magnetics(&masMagneticsWithScoring, inputs, _models, weights[CoreAdviserFilters::CORE_TEMPERATURE]);
                    break;
                case CoreAdviserFilters::DIMENSIONS: 

                    masMagneticsWithScoring = filterDimensions.filter_magnetics(&masMagneticsWithScoring, weights[CoreAdviserFilters::DIMENSIONS]);
                    break;
            }    
            std::string filterString = std::string{magic_enum::enum_name(filter)};
            logEntry("There are " + std::to_string(masMagneticsWithScoring.size()) + " after filtering by " + filterString + ".");
        }
    });

    if (masMagneticsWithScoring.size() > maximumNumberResults) {
        if (_uniqueCoreShapes) {
            std::vector<std::pair<MasWrapper, double>> masMagneticsWithScoringAndUniqueShapes;
            std::vector<std::string> useShapes;
            for (size_t masIndex = 0; masIndex < masMagneticsWithScoring.size(); ++masIndex){
                MasWrapper mas = masMagneticsWithScoring[masIndex].first;
                auto core = mas.get_magnetic().get_core();
                if (std::find(useShapes.begin(), useShapes.end(), core.get_shape_name()) != useShapes.end()) {
                    continue;
                }
                else {
                    masMagneticsWithScoringAndUniqueShapes.push_back(masMagneticsWithScoring[masIndex]);
                    useShapes.push_back(core.get_shape_name());
                }

                if (masMagneticsWithScoringAndUniqueShapes.size() == maximumNumberResults) {
                    break;
                }
            }
            masMagneticsWithScoring.clear();
            masMagneticsWithScoring = masMagneticsWithScoringAndUniqueShapes;
        }
        else {
            masMagneticsWithScoring = std::vector<std::pair<MasWrapper, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
        }
    }

    correct_windings(&masMagneticsWithScoring, inputs);

    return masMagneticsWithScoring;
}


} // namespace OpenMagnetics
