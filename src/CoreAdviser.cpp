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
#include <string>
#include <vector>
#include <execution>
#include <magic_enum_utility.hpp>
#include <list>


namespace OpenMagnetics {
void normalize_scoring(std::vector<std::pair<MasWrapper, double>>* masMagneticsWithScoring, std::vector<double>* newScoring, double weight, bool invertScoring = true) {
    double maxScoring = (*std::max_element((*masMagneticsWithScoring).begin(), (*masMagneticsWithScoring).end(), [](const std::pair<MasWrapper, double>& b1, const std::pair<MasWrapper, double>& b2) {
        return b1.second < b2.second;
    })).second;
    double minScoring = (*std::min_element((*masMagneticsWithScoring).begin(), (*masMagneticsWithScoring).end(), [](const std::pair<MasWrapper, double>& b1, const std::pair<MasWrapper, double>& b2) {
        return b1.second < b2.second;
    })).second;

    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i) {
        auto magnetic = (*masMagneticsWithScoring)[i].first;
        auto score = (*newScoring)[i];
        if (maxScoring != minScoring) {
            if (invertScoring) {
                (*masMagneticsWithScoring)[i].second += weight * (1 - (score - minScoring) / (maxScoring - minScoring));
            }
            else {
                (*masMagneticsWithScoring)[i].second += weight * ((score - minScoring) / (maxScoring - minScoring));
            }
        }
        else {
            (*masMagneticsWithScoring)[i].second += 1;
        }
    }
    sort((*masMagneticsWithScoring).begin(), (*masMagneticsWithScoring).end(), [](const std::pair<MasWrapper, double>& b1, const std::pair<MasWrapper, double>& b2) {
        return b1.second > b2.second;
    });
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterAreaProduct::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, double weight) {
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
        areaProductRequiredPreCalculations.push_back(powerMean / (primaryAreaFactor * 2 * frequency * defaults.coreAdviserMaximumCurrentDensity));
    }

    std::map<std::string, double> scaledMagneticFluxDensitiesPerMaterial;
    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));
    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLossesModel();  // TODO change to factory

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
        auto core = CoreWrapper(magnetic.get_core());
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
            auto skinDepth = windingSkinEffectLossesModel.get_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
            double wireAirFillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(2 * skinDepth);
            double windingWindowUtilizationFactor = wireAirFillingFactor * bobbinFillingFactor;
            double magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
            if (!materialScaledMagneticFluxDensities.contains(core.get_material_name())) {
                auto coreLossesMethods = core.get_available_core_losses_methods();
                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
                    double referenceCoreLosses = coreLossesModelSteinmetz->get_core_losses(core, operatingPointExcitation, temperature)["totalLosses"];
                    auto aux = coreLossesModelSteinmetz->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                    magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
                }
                else {
                    double referenceCoreLosses = coreLossesModelProprietary->get_core_losses(core, operatingPointExcitation, temperature)["totalLosses"];
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
            add_scoring(core.get_name().value(), CoreAdviser::CoreAdviserFilters::AREA_PRODUCT, scoring);
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
    (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterEnergyStored::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight) {
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
            auto core = CoreWrapper(magnetic.get_core());
            add_scoring(core.get_name().value(), CoreAdviser::CoreAdviserFilters::ENERGY_STORED, scoring);
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
    (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterWindingWindowArea::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, double weight) {
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

    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLossesModel();  // TODO change to factory
    auto skinDepth = windingSkinEffectLossesModel.get_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
    double wireAirFillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(2 * skinDepth);
    double estimatedWireConductingArea = std::numbers::pi * pow(skinDepth, 2);
    double estimatedWireTotalArea = estimatedWireConductingArea / wireAirFillingFactor;
    double necessaryWireCopperArea = primaryCurrentRms / defaults.coreAdviserMaximumCurrentDensity;
    double estimatedParallels = ceil(necessaryWireCopperArea / estimatedWireConductingArea);

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = CoreWrapper(magnetic.get_core());
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
            add_scoring(core.get_name().value(), CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA, scoring);
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
    (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterCoreLosses::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight) {
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


    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));

    OpenMagnetics::MagnetizingInductance magnetizing_inductance(models);
    auto windingOhmicLosses = OpenMagnetics::WindingOhmicLosses();
    size_t numberTimeouts = 0;

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = CoreWrapper(magnetic.get_core());

        std::string shapeName = core.get_shape_name();
        if (!(shapeName.contains("PQI") || shapeName.contains("R ") || shapeName.contains("T ") || shapeName.contains("UI "))) {
            auto bobbin = OpenMagnetics::BobbinWrapper::create_quick_bobbin(core);
            magnetic.get_mutable_coil().set_bobbin(bobbin);
        }

        int64_t currentNumberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();
        OpenMagnetics::NumberTurns numberTurns(currentNumberTurns);
        std::vector<double> totalLossesPerOperatingPoint;
        std::vector<double> coreLossesPerOperatingPoint;
        std::vector<double> ohmicLossesPerOperatingPoint;
        double currentTotalLosses = DBL_MAX;
        double coreLosses = DBL_MAX;
        double ohmicLosses = DBL_MAX;
        double newTotalLosses = DBL_MAX;
        int64_t previousNumberTurnsPrimary = currentNumberTurns;

        size_t iteration = 50;

        CoilWrapper winding = CoilWrapper(magnetic.get_coil());

        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            auto operatingPoint = inputs.get_operating_point(operatingPointIndex);
            double temperature = operatingPoint.get_conditions().get_ambient_temperature();
            OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
            do {
                currentTotalLosses = newTotalLosses;
                auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
                winding.get_mutable_functional_description()[0].set_number_turns(numberTurnsCombination[0]);
                auto aux = magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint);
                auto magnetizingInductance = aux.first;
                auto magneticFluxDensity = aux.second;

                if (!check_requirement(inputs.get_design_requirements().get_magnetizing_inductance(), magnetizingInductance)) {
                    winding.get_mutable_functional_description()[0].set_number_turns(previousNumberTurnsPrimary);
                    break;
                }
                else {
                    previousNumberTurnsPrimary = numberTurnsCombination[0];
                }

                excitation.set_magnetic_flux_density(magneticFluxDensity);
                auto coreLossesMethods = core.get_available_core_losses_methods();
                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
                    coreLosses = coreLossesModelSteinmetz->get_core_losses(core, excitation, temperature)["totalLosses"];
                }
                else {
                    coreLosses = coreLossesModelProprietary->get_core_losses(core, excitation, temperature)["totalLosses"];
                }

                if (!winding.get_turns_description()) {
                    newTotalLosses = coreLosses;
                    break;
                }
                if (!(shapeName.contains("PQI") || shapeName.contains("R ") || shapeName.contains("T ") || shapeName.contains("UI "))) {
                    ohmicLosses = windingOhmicLosses.get_ohmic_losses(winding, operatingPoint, temperature);
                    newTotalLosses = coreLosses + ohmicLosses;

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


            magnetic.set_coil(winding);

            currentTotalLosses = newTotalLosses;
            totalLossesPerOperatingPoint.push_back(currentTotalLosses);
            coreLossesPerOperatingPoint.push_back(coreLosses);
            ohmicLossesPerOperatingPoint.push_back(ohmicLosses);
        }


        double meanTotalLosses = 0;
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            meanTotalLosses += totalLossesPerOperatingPoint[operatingPointIndex];
        }
        meanTotalLosses /= inputs.get_operating_points().size();
        double maximumPowerMean = *max_element(powerMeans.begin(), powerMeans.end());

        if (meanTotalLosses < maximumPowerMean * defaults.coreAdviserMaximumPercentagePowerCoreLosses / defaults.coreAdviserThresholdValidity) {
            for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
                WindingLossesOutput windingLossesOutput;
                WindingLossesPerElement windingLossesPerElement;
                OhmicLosses ohmicLossesObj;

                CoreLossesOutput coreLossesOutput;
                coreLossesOutput.set_method_used(models["coreLosses"]);
                coreLossesOutput.set_origin(ResultOrigin::SIMULATION);
                coreLossesOutput.set_core_losses(coreLossesPerOperatingPoint[operatingPointIndex]);
                ohmicLossesObj.set_losses(ohmicLossesPerOperatingPoint[operatingPointIndex]);
                ohmicLossesObj.set_method_used("Ohm");
                ohmicLossesObj.set_origin(ResultOrigin::SIMULATION);
                windingLossesPerElement.set_ohmic_losses(ohmicLossesObj);
                windingLossesOutput.set_method_used("Ohm");
                windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
                auto windingLossesPerWinding = std::vector<WindingLossesPerElement>{windingLossesPerElement};
                double windingLosses = 0;
                for (auto& loss : windingLossesPerWinding) {
                    windingLosses += loss.get_ohmic_losses().value().get_losses();
                }

                windingLossesOutput.set_winding_losses(windingLosses);
                windingLossesOutput.set_winding_losses_per_winding(windingLossesPerWinding);
                mas.get_mutable_outputs()[operatingPointIndex].set_core_losses(coreLossesOutput);
                mas.get_mutable_outputs()[operatingPointIndex].set_winding_losses(windingLossesOutput);
            }
            double scoring = meanTotalLosses;
            newScoring.push_back(scoring);
            (*unfilteredMasMagnetics)[masIndex].first = mas;
            add_scoring(core.get_name().value(), CoreAdviser::CoreAdviserFilters::CORE_LOSSES, scoring);
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
    (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterCoreTemperature::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    auto defaults = Defaults();
    std::vector<std::pair<MasWrapper, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));

    OpenMagnetics::MagnetizingInductance magnetizing_inductance(models);

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = CoreWrapper(magnetic.get_core());
        auto winding = magnetic.get_coil();
        bool validMagnetic = true;
        double meanTemperature = 0;
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            double temperature = inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature();
            auto operatingPoint = inputs.get_operating_point(operatingPointIndex);
            OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
            double coreLosses;
            if (!mas.get_outputs()[operatingPointIndex].get_core_losses()) {
                auto magneticFluxDensity = magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;
                excitation.set_magnetic_flux_density(magneticFluxDensity);
                auto coreLossesMethods = core.get_available_core_losses_methods();
                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
                    coreLosses = coreLossesModelSteinmetz->get_core_losses(core, excitation, temperature)["totalLosses"];
                }
                else {
                    coreLosses = coreLossesModelProprietary->get_core_losses(core, excitation, temperature)["totalLosses"];
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
            double calculatedTemperature = coreTemperature["maximumTemperature"];
            meanTemperature += calculatedTemperature;

            if (calculatedTemperature < defaults.coreAdviserMaximumCoreTemperature / defaults.coreAdviserThresholdValidity) {
                mas.get_mutable_outputs()[operatingPointIndex].get_core_losses().value().set_temperature(calculatedTemperature);
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
            add_scoring(core.get_name().value(), CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE, scoring);
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
    (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, &newScoring, weight);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::MagneticCoreFilterDimensions::filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, double weight) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    std::vector<double> newScoring;

    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        MasWrapper mas = (*unfilteredMasMagnetics)[masIndex].first;
        MagneticWrapper magnetic = MagneticWrapper(mas.get_magnetic());
        auto core = CoreWrapper(magnetic.get_core());
        double volume = core.get_width() * core.get_height() * core.get_depth();

        double scoring = volume;
        newScoring.push_back(scoring);
        add_scoring(core.get_name().value(), CoreAdviser::CoreAdviserFilters::DIMENSIONS, scoring);
    }

    if ((*unfilteredMasMagnetics).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string((*unfilteredMasMagnetics).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredMasMagnetics).size() > 0) {
        normalize_scoring(unfilteredMasMagnetics, &newScoring, weight);
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
    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLossesModel();  // TODO change to factory
    auto skinDepth = windingSkinEffectLossesModel.get_skin_depth("copper", frequency, temperature);  // TODO material hardcoded

    // Set round wire with diameter to two times the skin depth 
    WireWrapper wire;
    wire.set_nominal_value_conducting_diameter(skinDepth * 2);
    wire.set_nominal_value_outer_diameter(skinDepth * 2.05); // Hardcoded
    wire.set_material("copper");
    wire.set_type("round");
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
    auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/../../MAS/data/cores.ndjson");
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
                mas.set_magnetic(magnetic);
                masMagnetics.push_back(std::pair<MasWrapper, double>{mas, 0});
            }
        }
        else {
            magnetic.set_core(core);
            mas.set_magnetic(magnetic);
            masMagnetics.push_back(std::pair<MasWrapper, double>{mas, 0});
        }
    }

    return masMagnetics;
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::get_advised_core(InputsWrapper inputs, std::vector<CoreWrapper>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults, maximumNumberCores);
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, std::vector<CoreWrapper>* cores, size_t maximumNumberResults){
    auto defaults = Defaults();
    _weights = weights;
 
    CoreAdviserFilters firstFilter = CoreAdviserFilters::AREA_PRODUCT;
    double maxWeight = 0;

    for (auto const& pair : weights) {
        auto filter = pair.first;
        auto weight = pair.second;
        if (weight > maxWeight) {
            maxWeight = weight;
            firstFilter = filter;
        }
    }
    size_t maximumMagneticsAfterFiltering = defaults.coreAdviserMaximumMagneticsAfterFiltering;

    if (firstFilter == CoreAdviserFilters::CORE_LOSSES || firstFilter == CoreAdviserFilters::CORE_TEMPERATURE) {
        auto masMagneticsWithStacks = create_mas_dataset(inputs, cores, true);
        log("We start the search with " + std::to_string(masMagneticsWithStacks.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.");
        std::string firstFilterString = std::string{magic_enum::enum_name(firstFilter)};
        log("We include stacks of cores in our search because the most important selectd filter is " + firstFilterString + ".");
        auto filteredMasMagnetics = apply_filters(&masMagneticsWithStacks, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMasMagnetics.size() >= maximumNumberResults) {
            return filteredMasMagnetics;
        }

    }
    else {
        auto masMagneticsWithoutStacks = create_mas_dataset(inputs, cores, false);
        log("We start the search with " + std::to_string(masMagneticsWithoutStacks.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.");
        log("We don't include stacks of cores in our search.");
        auto filteredMasMagnetics = apply_filters(&masMagneticsWithoutStacks, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMasMagnetics.size() >= maximumNumberResults) {
            return filteredMasMagnetics;
        }
    }

    auto masMagneticsWithStacks = create_mas_dataset(inputs, cores, true);
    log("First attempt produced not enough results, so now we are searching again with " + std::to_string(masMagneticsWithStacks.size()) + " magnetics, including up to " + std::to_string(defaults.coreAdviserMaximumNumberStacks) + " cores stacked when possible.");
    weights[CoreAdviserFilters::CORE_LOSSES] = 2;
    maximumMagneticsAfterFiltering = masMagneticsWithStacks.size();
    auto filteredMasMagnetics = apply_filters(&masMagneticsWithStacks, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    return filteredMasMagnetics;
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

    sort(results.begin(), results.end(), [](const std::pair<MasWrapper, double>& b1, const std::pair<MasWrapper, double>& b2) {
        return b1.second > b2.second;
    });

    if (results.size() > maximumNumberResults) {
        results = std::vector<std::pair<MasWrapper, double>>(results.begin(), results.end() - (results.size() - maximumNumberResults));
    }
    return results;
}

void add_initial_turns(std::vector<std::pair<MasWrapper, double>> *masMagneticsWithScoring, InputsWrapper inputs) {
    MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"gapReluctance", "ZHANG"}}));
    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i){

        CoreWrapper core = CoreWrapper((*masMagneticsWithScoring)[i].first.get_magnetic().get_core());
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        double numberTurns = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs, DimensionalValues::MINIMUM);
        (*masMagneticsWithScoring)[i].first.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(numberTurns);
    }
}

std::vector<std::pair<MasWrapper, double>> CoreAdviser::apply_filters(std::vector<std::pair<MasWrapper, double>>* masMagnetics, InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
    MagneticCoreFilterAreaProduct filterAreaProduct(&_scorings);
    MagneticCoreFilterEnergyStored filterEnergyStored(&_scorings);
    MagneticCoreFilterWindingWindowArea filterWindingWindowArea(&_scorings);
    MagneticCoreFilterCoreLosses filterCoreLosses(&_scorings);
    MagneticCoreFilterCoreTemperature filterCoreTemperature(&_scorings);
    MagneticCoreFilterDimensions filterDimensions(&_scorings);

    CoreAdviserFilters firstFilter = CoreAdviserFilters::AREA_PRODUCT;
    double maxWeight = 0;

    for (auto const& pair : weights) {
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
            masMagneticsWithScoring = filterAreaProduct.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::AREA_PRODUCT]);
            break;
        case CoreAdviserFilters::ENERGY_STORED: 
            masMagneticsWithScoring = filterEnergyStored.filter_magnetics(masMagnetics, inputs, _models, weights[CoreAdviserFilters::ENERGY_STORED]);
            break;
        case CoreAdviserFilters::WINDING_WINDOW_AREA: 
            masMagneticsWithScoring = filterWindingWindowArea.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::WINDING_WINDOW_AREA]);
            break;
        case CoreAdviserFilters::CORE_LOSSES: 
            add_initial_turns(masMagnetics, inputs);
            masMagneticsWithScoring = filterCoreLosses.filter_magnetics(masMagnetics, inputs, _models, weights[CoreAdviserFilters::CORE_LOSSES]);
            break;
        case CoreAdviserFilters::CORE_TEMPERATURE: 
            add_initial_turns(masMagnetics, inputs);
            masMagneticsWithScoring = filterCoreTemperature.filter_magnetics(masMagnetics, inputs, _models, weights[CoreAdviserFilters::CORE_TEMPERATURE]);
            break;
        case CoreAdviserFilters::DIMENSIONS: 
            masMagneticsWithScoring = filterDimensions.filter_magnetics(masMagnetics, weights[CoreAdviserFilters::DIMENSIONS]);
            break;
    }

    std::string firstFilterString = std::string{magic_enum::enum_name(firstFilter)};
    log("There are " + std::to_string(masMagneticsWithScoring.size()) + " magnetics after the first filter, which was " + firstFilterString + ".");

    if (masMagneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
        masMagneticsWithScoring = std::vector<std::pair<MasWrapper, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumMagneticsAfterFiltering));
        log("There are " + std::to_string(masMagneticsWithScoring.size()) + " after culling by the score on the first filter.");
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
            log("There are " + std::to_string(masMagneticsWithScoring.size()) + " after filtering by " + filterString + ".");
        }
    });

    if (masMagneticsWithScoring.size() > maximumNumberResults) {
        masMagneticsWithScoring = std::vector<std::pair<MasWrapper, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
    }

    return masMagneticsWithScoring;
}


} // namespace OpenMagnetics
