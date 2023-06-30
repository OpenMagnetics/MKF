#include "CoreAdviser.h"
#include "MagnetizingInductance.h"
#include "WindingSkinEffectLosses.h"
#include "MagneticEnergy.h"
#include "WireWrapper.h"
#include "BobbinWrapper.h"
#include "Defaults.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <execution>


namespace OpenMagnetics {


std::map<std::string, double> normalize_scoring(std::map<std::string, double> scoring) {
    std::map<std::string, double> normalizedScoring;
    double maxScoring = (*std::max_element(scoring.begin(), scoring.end(), [](const std::pair<std::string, double>& b1, const std::pair<std::string, double>& b2) {
                        return b1.second < b2.second;
            })).second;
    double minScoring = (*std::min_element(scoring.begin(), scoring.end(), [](const std::pair<std::string, double>& b1, const std::pair<std::string, double>& b2) {
                        return b1.second < b2.second;
            })).second;
    for (auto& pair : scoring) {
        auto key = pair.first;
        auto score = pair.second;
        if (maxScoring != minScoring) {
            normalizedScoring[key] = (score - minScoring) / (maxScoring - minScoring);
        }
        else {
            normalizedScoring[key] = 1;
        }
    }
    return normalizedScoring;
}

std::pair<std::vector<Magnetic>, std::map<std::string, double>> CoreAdviser::MagneticCoreFilterAreaProduct::filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models) {
    auto defaults = Defaults();
    std::vector<Magnetic> filteredMagnetics;
    std::map<std::string, double> scoring;
    std::vector<double> voltageWaveformData = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_voltage().value().get_waveform().value().get_data();
    std::vector<double> currentWaveformData = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_current().value().get_waveform().value().get_data();

    std::map<std::string, double> materialScaledMagneticFluxDensities;
    std::map<std::string, double> bobbinFillingFactors;

    if (voltageWaveformData.size() != currentWaveformData.size()) {
        throw std::runtime_error("voltage and current waveform should have the same number of points");
    }

    double powerMean = 0;
    for (size_t i = 0; i < voltageWaveformData.size(); ++i)
    {
        powerMean += fabs(voltageWaveformData[i] * currentWaveformData[i]);
    }
    powerMean /= voltageWaveformData.size();

    double frequency = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_frequency();
    double temperature = inputs.get_operation_point(0).get_conditions().get_ambient_temperature();
    std::map<std::string, double> scaledMagneticFluxDensitiesPerMaterial;
    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));
    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLossesModel();  // TODO change to factory

    auto skinDepth = windingSkinEffectLossesModel.get_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
    double wireAirFillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(2 * skinDepth);
    double primaryAreaFactor = 1;
    if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
        primaryAreaFactor = 0.5;
    }
    double magneticFluxDensityReference = 0.18;
    double frequencyReference = 100000;
    double areaProductRequiredPreCalculation = powerMean / (primaryAreaFactor * 2 * frequency * defaults.coreAdviserMaximumCurrentDensity);
    OperationPointExcitation operationPointExcitation;
    ElectromagneticParameter magneticFluxDensity;
    Processed processed;
    operationPointExcitation.set_frequency(frequencyReference);
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);
    processed.set_peak(magneticFluxDensityReference);
    processed.set_peak_to_peak(2 * magneticFluxDensityReference);
    magneticFluxDensity.set_processed(processed);
    operationPointExcitation.set_magnetic_flux_density(magneticFluxDensity);

    std::cout << unfilteredMagnetics.size() << std::endl;

    for (auto& magnetic : unfilteredMagnetics){
        auto core = CoreWrapper(magnetic.get_core());
        double bobbinFillingFactor;
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

        double magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
        if (!materialScaledMagneticFluxDensities.contains(core.get_material_name())) {
            auto coreLossesMethods = core.get_available_core_losses_methods();
            if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
                double referenceCoreLosses = coreLossesModelSteinmetz->get_core_losses(core, operationPointExcitation, temperature)["totalLosses"];
                auto aux = coreLossesModelSteinmetz->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
            }
            else {
                double referenceCoreLosses = coreLossesModelProprietary->get_core_losses(core, operationPointExcitation, temperature)["totalLosses"];
                auto aux = coreLossesModelProprietary->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
            }
            materialScaledMagneticFluxDensities[core.get_material_name()] = magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
        }
        else {
            magneticFluxDensityPeakAtFrequencyOfReferenceLosses = materialScaledMagneticFluxDensities[core.get_material_name()];
        }

        double windingWindowUtilizationFactor = wireAirFillingFactor * bobbinFillingFactor;
        double areaProductRequired = areaProductRequiredPreCalculation / (windingWindowUtilizationFactor * magneticFluxDensityPeakAtFrequencyOfReferenceLosses);
        double areaProductCore = windingWindow.get_area().value() * windingColumn.get_area();
        if (std::isnan(magneticFluxDensityPeakAtFrequencyOfReferenceLosses) || magneticFluxDensityPeakAtFrequencyOfReferenceLosses == 0) {
            std::cout << core.get_material_name()  << std::endl;
            throw std::runtime_error("magneticFluxDensityPeakAtFrequencyOfReferenceLosses cannot be 0 or NaN");
            break;
        }
        if (std::isnan(areaProductRequired)) {
            break;
        }
        if (std::isinf(areaProductRequired) || areaProductRequired == 0) {
            throw std::runtime_error("areaProductRequired cannot be 0 or NaN");
        }

        if (areaProductCore >= areaProductRequired * defaults.coreAdviserThresholdValidity) {
            filteredMagnetics.push_back(magnetic);   
            scoring[core.get_name().value()] = fabs(areaProductCore - areaProductRequired);
        }
    }
    normalize_scoring(scoring);

    return {filteredMagnetics, scoring};
}

std::pair<std::vector<Magnetic>, std::map<std::string, double>> CoreAdviser::MagneticCoreFilterEnergyStored::filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models) {
    auto defaults = Defaults();
    OpenMagnetics::MagneticEnergy magneticEnergy(models);
    std::vector<Magnetic> filteredMagnetics;
    std::map<std::string, double> scoring;
    double requiredMagneticEnergy = magneticEnergy.required_magnetic_energy(inputs).get_nominal().value();
    auto operationPoint = inputs.get_operation_point(0);

    for (auto& magnetic : unfilteredMagnetics){  
        double totalStorableMagneticEnergy = magneticEnergy.get_core_maximum_magnetic_energy(static_cast<CoreWrapper>(magnetic.get_core()), &operationPoint);
        if (totalStorableMagneticEnergy >= requiredMagneticEnergy * defaults.coreAdviserThresholdValidity) {
            filteredMagnetics.push_back(magnetic);
            scoring[magnetic.get_core().get_name().value()] = 1.;
        }
    }

    return {filteredMagnetics, scoring};
}

std::pair<std::vector<Magnetic>, std::map<std::string, double>> CoreAdviser::MagneticCoreFilterWindingWindowArea::filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models) {
    auto defaults = Defaults();
    std::vector<Magnetic> filteredMagnetics;
    std::map<std::string, double> scoring;
    auto primaryCurrentRms = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_current().value().get_processed().value().get_rms().value();

    double frequency = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_frequency();
    double temperature = inputs.get_operation_point(0).get_conditions().get_ambient_temperature();
    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLossesModel();  // TODO change to factory
    auto skinDepth = windingSkinEffectLossesModel.get_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
    double wireAirFillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(2 * skinDepth);
    double estimatedWireConductingArea = std::numbers::pi * pow(skinDepth, 2);
    double estimatedWireTotalArea = estimatedWireConductingArea / wireAirFillingFactor;
    double necessaryWireCopperArea = primaryCurrentRms / defaults.coreAdviserMaximumCurrentDensity;
    double estimatedParallels = ceil(necessaryWireCopperArea / estimatedWireConductingArea);

    for (auto& magnetic : unfilteredMagnetics){
        auto core = CoreWrapper(magnetic.get_core());
        double primaryNumberTurns = magnetic.get_winding().get_functional_description()[0].get_number_turns();
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

            filteredMagnetics.push_back(magnetic);
            scoring[magnetic.get_core().get_name().value()] = manufacturabilityRelativeCost;
        }
    }
    normalize_scoring(scoring);

    return {filteredMagnetics, scoring};
}

std::pair<std::vector<Magnetic>, std::map<std::string, double>> CoreAdviser::MagneticCoreFilterCoreLosses::filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models) {
    auto defaults = Defaults();
    std::vector<Magnetic> filteredMagnetics;
    std::map<std::string, double> scoring;

    std::vector<double> voltageWaveformData = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_voltage().value().get_waveform().value().get_data();
    std::vector<double> currentWaveformData = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_current().value().get_waveform().value().get_data();

    double powerMean = 0;
    for (size_t i = 0; i < voltageWaveformData.size(); ++i)
    {
        powerMean += fabs(voltageWaveformData[i] * currentWaveformData[i]);
    }
    powerMean /= voltageWaveformData.size();

    double frequency = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_frequency();
    double temperature = inputs.get_operation_point(0).get_conditions().get_ambient_temperature();

    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));

    OpenMagnetics::MagnetizingInductance magnetizing_inductance(models);
    auto operationPoint = inputs.get_operation_point(0);
    OpenMagnetics::OperationPointExcitation excitation = operationPoint.get_excitations_per_winding()[0];

    for (auto& magnetic : unfilteredMagnetics){
        auto core = CoreWrapper(magnetic.get_core());
        auto winding = magnetic.get_winding();

        auto magneticFluxDensity = magnetizing_inductance.get_inductance_and_magnetic_flux_density(core, winding, &operationPoint).second;
        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double coreLosses;
        auto coreLossesMethods = core.get_available_core_losses_methods();
        if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
            coreLosses = coreLossesModelSteinmetz->get_core_losses(core, excitation, temperature)["totalLosses"];
        }
        else {
            coreLosses = coreLossesModelProprietary->get_core_losses(core, excitation, temperature)["totalLosses"];
        }

        // auto coreTemperatureModel = OpenMagnetics::CoreTemperatureModel::factory(CoreTemperatureModels::MANIKTALA);
        // auto coreTemperature = coreTemperatureModel->get_core_temperature(core, coreLosses, temperature);
        // double calculatedTemperature = coreTemperature["maximumTemperature"];

        // if (calculatedTemperature < defaults.coreAdviserMaximumCoreTemperature)
        if (coreLosses < powerMean * defaults.coreAdviserMaximumPercentagePowerCoreLosses / defaults.coreAdviserThresholdValidity) {
            filteredMagnetics.push_back(magnetic);
            scoring[core.get_name().value()] = coreLosses;
        }

    }
    normalize_scoring(scoring);

    return {filteredMagnetics, scoring};
}

std::pair<std::vector<Magnetic>, std::map<std::string, double>> CoreAdviser::MagneticCoreFilterCoreTemperature::filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, std::map<std::string, double> precalculatedCoreLosses) {
    auto defaults = Defaults();
    std::vector<Magnetic> filteredMagnetics;
    std::map<std::string, double> scoring;
    auto primaryCurrentRms = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_current().value().get_processed().value().get_rms().value();

    double frequency = InputsWrapper::get_primary_excitation(inputs.get_operation_point(0)).get_frequency();
    double temperature = inputs.get_operation_point(0).get_conditions().get_ambient_temperature();

    auto coreLossesModelSteinmetz = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = OpenMagnetics::CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));


    OpenMagnetics::MagnetizingInductance magnetizing_inductance(models);
    auto operationPoint = inputs.get_operation_point(0);
    OpenMagnetics::OperationPointExcitation excitation = operationPoint.get_excitations_per_winding()[0];

    for (auto& magnetic : unfilteredMagnetics){
        auto core = CoreWrapper(magnetic.get_core());
        auto winding = magnetic.get_winding();

        double coreLosses;
        if (!precalculatedCoreLosses.contains(core.get_name().value())) {
            std::cout << "mierda 1" << std::endl;
            auto magneticFluxDensity = magnetizing_inductance.get_inductance_and_magnetic_flux_density(core, winding, &operationPoint).second;
            excitation.set_magnetic_flux_density(magneticFluxDensity);
            auto coreLossesMethods = core.get_available_core_losses_methods();
            if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), "steinmetz") != coreLossesMethods.end()) {
                coreLosses = coreLossesModelSteinmetz->get_core_losses(core, excitation, temperature)["totalLosses"];
            }
            else {
                coreLosses = coreLossesModelProprietary->get_core_losses(core, excitation, temperature)["totalLosses"];
            }
        }
        else {
            coreLosses = precalculatedCoreLosses[core.get_name().value()];
        }

        auto coreTemperatureModel = OpenMagnetics::CoreTemperatureModel::factory(CoreTemperatureModels::MANIKTALA);
        auto coreTemperature = coreTemperatureModel->get_core_temperature(core, coreLosses, temperature);
        double calculatedTemperature = coreTemperature["maximumTemperature"];

        if (calculatedTemperature < defaults.coreAdviserMaximumCoreTemperature / defaults.coreAdviserThresholdValidity) {
            filteredMagnetics.push_back(magnetic);
            scoring[core.get_name().value()] = calculatedTemperature;
        }

    }
    normalize_scoring(scoring);

    return {filteredMagnetics, scoring};
}


WindingWrapper get_dummy_winding() {
    WindingFunctionalDescription primaryWindingFunctionalDescription;
    primaryWindingFunctionalDescription.set_isolation_side(IsolationSide::PRIMARY);
    primaryWindingFunctionalDescription.set_name("primary");
    primaryWindingFunctionalDescription.set_number_parallels(1);
    primaryWindingFunctionalDescription.set_number_turns(1);
    primaryWindingFunctionalDescription.set_wire("Dummy");

    WindingWrapper winding;
    winding.set_bobbin("Dummy");
    winding.set_functional_description({primaryWindingFunctionalDescription});
    return winding;
}

CoreWrapper CoreAdviser::get_advised_core(InputsWrapper inputs) {
    std::string file_path = __FILE__;
    auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/../../MAS/data/cores.ndjson");
    std::cout << inventory_path << std::endl;
    std::ifstream ndjsonFile(inventory_path);
    std::string jsonLine;
    std::vector<CoreWrapper> cores;
    while (std::getline(ndjsonFile, jsonLine)) {
        json jf = json::parse(jsonLine);
        CoreWrapper core(jf, false, true, false);
        cores.push_back(core);
    }

    std::vector<Magnetic> magnetics;
    MagnetizingInductance magnetizing_inductance(
        std::map<std::string, std::string>({{"gapReluctance", "ZHANG"}}));
    WindingWrapper winding = get_dummy_winding();

    Magnetic magnetic;
    magnetic.set_winding(std::move(winding));
    for (auto& core : cores){  
        magnetic.set_core(std::move(core));
        magnetics.push_back(std::move(magnetic));
    }

    auto aux = MagneticCoreFilterAreaProduct::filter_magnetics(magnetics, inputs, _models);
    magnetics = aux.first;
    aux = MagneticCoreFilterEnergyStored::filter_magnetics(magnetics, inputs, _models);
    magnetics = aux.first;

    for (size_t i = 0; i < magnetics.size(); ++i){
        double numberTurns = magnetizing_inductance.get_number_turns_from_gapping_and_inductance(cores[i], &inputs);
        magnetics[i].get_mutable_winding().get_mutable_functional_description()[0].set_number_turns(numberTurns);
    }

    aux = MagneticCoreFilterWindingWindowArea::filter_magnetics(magnetics, inputs, _models);
    magnetics = aux.first;
    aux = MagneticCoreFilterCoreLosses::filter_magnetics(magnetics, inputs, _models);
    magnetics = aux.first;
    aux = MagneticCoreFilterCoreTemperature::filter_magnetics(magnetics, inputs, _models, aux.second);
    magnetics = aux.first;

    std::cout << cores.size() << std::endl;
    std::cout << magnetics.size() << std::endl;

    return cores[0];
}

} // namespace OpenMagnetics
