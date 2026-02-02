#include "advisers/MagneticFilter.h"
#include "constructive_models/NumberTurns.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/Impedance.h"
#include <magic_enum.hpp>

#include <cfloat>
#include <cmath>
#include <algorithm>
#include "support/Exceptions.h"
#include "support/Logger.h"

namespace OpenMagnetics {


std::shared_ptr<MagneticFilter> MagneticFilter::factory(MagneticFilters filterName, std::optional<Inputs> inputs) {
    switch(filterName) {
        case MagneticFilters::AREA_PRODUCT:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter AREA_PRODUCT");
            }
            return std::make_shared<MagneticFilterAreaProduct>(inputs.value());
        case MagneticFilters::ENERGY_STORED:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter ENERGY_STORED");
            }
            return std::make_shared<MagneticFilterEnergyStored>(inputs.value());
        case MagneticFilters::ESTIMATED_COST:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter ESTIMATED_COST");
            }
            return std::make_shared<MagneticFilterEstimatedCost>(inputs.value());
        case MagneticFilters::COST:
            return std::make_shared<MagneticFilterCost>();
        case MagneticFilters::CORE_AND_DC_LOSSES:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter CORE_AND_DC_LOSSES");
            }
            return std::make_shared<MagneticFilterCoreAndDcLosses>(inputs.value());
        case MagneticFilters::CORE_DC_AND_SKIN_LOSSES:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter CORE_DC_AND_SKIN_LOSSES");
            }
            return std::make_shared<MagneticFilterCoreDcAndSkinLosses>(inputs.value());
        case MagneticFilters::LOSSES:
            return std::make_shared<MagneticFilterLosses>();
        case MagneticFilters::LOSSES_NO_PROXIMITY:
            return std::make_shared<MagneticFilterLossesNoProximity>();
        case MagneticFilters::DIMENSIONS:
            return std::make_shared<MagneticFilterDimensions>();
        case MagneticFilters::CORE_MINIMUM_IMPEDANCE:
            return std::make_shared<MagneticFilterCoreMinimumImpedance>();
        case MagneticFilters::AREA_NO_PARALLELS:
            return std::make_shared<MagneticFilterAreaNoParallels>();
        case MagneticFilters::AREA_WITH_PARALLELS:
            return std::make_shared<MagneticFilterAreaWithParallels>();
        case MagneticFilters::EFFECTIVE_RESISTANCE:
            return std::make_shared<MagneticFilterEffectiveResistance>();
        case MagneticFilters::PROXIMITY_FACTOR:
            return std::make_shared<MagneticFilterProximityFactor>();
        case MagneticFilters::SOLID_INSULATION_REQUIREMENTS:
            return std::make_shared<MagneticFilterSolidInsulationRequirements>();
        case MagneticFilters::TURNS_RATIOS:
            return std::make_shared<MagneticFilterTurnsRatios>();
        case MagneticFilters::MAXIMUM_DIMENSIONS:
            return std::make_shared<MagneticFilterMaximumDimensions>();
        case MagneticFilters::SATURATION:
            return std::make_shared<MagneticFilterSaturation>();
        case MagneticFilters::DC_CURRENT_DENSITY:
            return std::make_shared<MagneticFilterDcCurrentDensity>();
        case MagneticFilters::EFFECTIVE_CURRENT_DENSITY:
            return std::make_shared<MagneticFilterEffectiveCurrentDensity>();
        case MagneticFilters::IMPEDANCE:
            return std::make_shared<MagneticFilterImpedance>();
        case MagneticFilters::MAGNETIZING_INDUCTANCE:
            return std::make_shared<MagneticFilterMagnetizingInductance>();
        case MagneticFilters::SKIN_LOSSES_DENSITY:
            return std::make_shared<MagneticFilterSkinLossesDensity>();
        case MagneticFilters::FRINGING_FACTOR:
            return std::make_shared<MagneticFilterFringingFactor>();
        case MagneticFilters::VOLUME:
            return std::make_shared<MagneticFilterVolume>();
        case MagneticFilters::AREA:
            return std::make_shared<MagneticFilterArea>();
        case MagneticFilters::HEIGHT:
            return std::make_shared<MagneticFilterHeight>();
        case MagneticFilters::TEMPERATURE_RISE:
            return std::make_shared<MagneticFilterTemperatureRise>();
        case MagneticFilters::LOSSES_TIMES_VOLUME:
            return std::make_shared<MagneticFilterLossesTimesVolume>();
        case MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE:
            return std::make_shared<MagneticFilterVolumeTimesTemperatureRise>();
        case MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE:
            return std::make_shared<MagneticFilterLossesTimesVolumeTimesTemperatureRise>();
        case MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter LOSSES_NO_PROXIMITY_TIMES_VOLUME");
            }
            return std::make_shared<MagneticFilterLossesNoProximityTimesVolume>();
        case MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE");
            }
            return std::make_shared<MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise>();
        case MagneticFilters::MAGNETOMOTIVE_FORCE:
            return std::make_shared<MagnetomotiveForce>();
        default:
            throw ModelNotAvailableException("Unknown filter, available options are: {AREA_PRODUCT, ENERGY_STORED, ESTIMATED_COST, COST, CORE_AND_DC_LOSSES, CORE_DC_AND_SKIN_LOSSES, LOSSES, LOSSES_NO_PROXIMITY, DIMENSIONS, CORE_MINIMUM_IMPEDANCE, AREA_NO_PARALLELS, AREA_WITH_PARALLELS, EFFECTIVE_RESISTANCE, PROXIMITY_FACTOR, SOLID_INSULATION_REQUIREMENTS, TURNS_RATIOS, MAXIMUM_DIMENSIONS, SATURATION, DC_CURRENT_DENSITY, EFFECTIVE_CURRENT_DENSITY, IMPEDANCE, MAGNETIZING_INDUCTANCE, FRINGING_FACTOR, SKIN_LOSSES_DENSITY, VOLUME, AREA, HEIGHT, TEMPERATURE_RISE, LOSSES_TIMES_VOLUME, VOLUME_TIMES_TEMPERATURE_RISE, LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE, LOSSES_NO_PROXIMITY_TIMES_VOLUME, LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE}");
    }
}

MagneticFilterAreaProduct::MagneticFilterAreaProduct(Inputs inputs) {
    double frequencyReference = 100000;
    SignalDescriptor magneticFluxDensity;
    Processed processed;
    _operatingPointExcitation.set_frequency(frequencyReference);
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);
    processed.set_peak(_magneticFluxDensityReference);
    processed.set_peak_to_peak(2 * _magneticFluxDensityReference);
    magneticFluxDensity.set_processed(processed);
    _operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);
    _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));

    if (settings.get_core_adviser_include_margin() && inputs.get_design_requirements().get_insulation()) {
        auto clearanceAndCreepageDistance = InsulationCoordinator().calculate_creepage_distance(inputs, true);
        _averageMarginInWindingWindow = clearanceAndCreepageDistance;
    }

    double primaryAreaFactor = 1;
    if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
        primaryAreaFactor = 0.5;
    }

    _areaProductRequiredPreCalculations.clear();
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        auto excitation = Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex));
        auto voltageWaveform = excitation.get_voltage().value().get_waveform().value();
        auto currentWaveform = excitation.get_current().value().get_waveform().value();
        double frequency = excitation.get_frequency();
        if (voltageWaveform.get_data().size() != currentWaveform.get_data().size()) {
            size_t maximumNumberPoints = constants.numberPointsSampledWaveforms;
            maximumNumberPoints = std::max(maximumNumberPoints, voltageWaveform.get_data().size());
            maximumNumberPoints = std::max(maximumNumberPoints, currentWaveform.get_data().size());
            voltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, frequency, maximumNumberPoints);
            currentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, frequency, maximumNumberPoints);
        }

        std::vector<double> voltageWaveformData = voltageWaveform.get_data();
        std::vector<double> currentWaveformData = currentWaveform.get_data();

        double powerMean = 0;
        for (size_t i = 0; i < voltageWaveformData.size(); ++i)
        {
            powerMean += fabs(voltageWaveformData[i] * currentWaveformData[i]);
        }
        powerMean /= voltageWaveformData.size();

        double switchingFrequency = Inputs::get_switching_frequency(excitation);
        double preCalculation = 0;

        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            preCalculation = powerMean / (primaryAreaFactor * 2 * switchingFrequency * defaults.maximumCurrentDensityPlanar);
        }
        else {
            preCalculation = powerMean / (primaryAreaFactor * 2 * switchingFrequency * defaults.maximumCurrentDensity);
        }

        if (preCalculation > 1) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "maximumAreaProductRequired cannot be larger than 1 (probably)");
        }
        _areaProductRequiredPreCalculations.push_back(preCalculation);
        if (std::isinf(_areaProductRequiredPreCalculations.back()) || _areaProductRequiredPreCalculations.back() == 0) {
            std::cerr << "powerMean: " << powerMean << std::endl;
            std::cerr << "operatingPointIndex: " << operatingPointIndex << std::endl;
            std::cerr << "primaryAreaFactor: " << primaryAreaFactor << std::endl;
            std::cerr << "switchingFrequency: " << switchingFrequency << std::endl;
            std::cerr << "_areaProductRequiredPreCalculations.back(): " << _areaProductRequiredPreCalculations.back() << std::endl;
            throw NaNResultException("_areaProductRequiredPreCalculations cannot be 0 or NaN");
        }
    }
}

std::pair<bool, double> MagneticFilterAreaProduct::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    double bobbinFillingFactor;
    if (core.get_winding_windows().size() == 0)
        return {false, 0.0};
    auto windingWindow = core.get_winding_windows()[0];
    auto windingColumn = core.get_columns()[0];
    if (inputs->get_wiring_technology() == WiringTechnology::PRINTED) {
        bobbinFillingFactor = 1;
    }
    else if (!_bobbinFillingFactors.contains(core.get_shape_name())) {
        if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
            bobbinFillingFactor = Bobbin::get_filling_factor(windingWindow.get_width().value(), core.get_winding_windows()[0].get_height().value());
        }
        else {
            // For toroids: apply a realistic filling factor penalty
            // The inner circumference is smaller than outer, limiting wire packing
            // Also, manual winding of toroids is less efficient than bobbin-based winding
            // Use a conservative factor of 0.55-0.70 depending on geometry
            auto radialHeight = windingWindow.get_radial_height().value();
            double outerRadius = core.get_width() / 2;
            double innerRadius = outerRadius - radialHeight;
            // Ratio of inner to outer circumference limits packing efficiency
            double circumferenceRatio = (innerRadius > 0) ? (innerRadius / outerRadius) : 0.5;
            // Base toroid filling factor around 0.6, adjusted by geometry
            bobbinFillingFactor = 0.55 + 0.15 * circumferenceRatio;
        }
        _bobbinFillingFactors[core.get_shape_name()] = bobbinFillingFactor;
    }
    else {
        bobbinFillingFactor = _bobbinFillingFactors[core.get_shape_name()];
    }
    double windingWindowArea = windingWindow.get_area().value();
    if (_averageMarginInWindingWindow > 0) {
        if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
            if (windingWindow.get_height().value() > windingWindow.get_width().value()) {
                windingWindowArea -= windingWindow.get_width().value() * _averageMarginInWindingWindow;
            }
            else {
                windingWindowArea -= windingWindow.get_height().value() * _averageMarginInWindingWindow;
            }
        }
        else {
            auto radialHeight = windingWindow.get_radial_height().value();
            if (_averageMarginInWindingWindow > radialHeight / 2) {
                return {false, 0.0};
            }
            double wireAngle = wound_distance_to_angle(_averageMarginInWindingWindow, radialHeight / 2);
            if (std::isnan((wireAngle) / 360)) {
                throw NaNResultException("wireAngle: " + std::to_string(wireAngle));
            }
            windingWindowArea *= (wireAngle) / 360;
        }
    }
    double areaProductCore = windingWindowArea * windingColumn.get_area();
    double maximumAreaProductRequired = 0;


    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        double temperature = inputs->get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature();
        // double frequency = Inputs::get_primary_excitation(inputs->get_operating_point(operatingPointIndex)).get_frequency();
        double frequency = Inputs::get_switching_frequency(Inputs::get_primary_excitation(inputs->get_operating_point(operatingPointIndex)));
        // double switchingFrequency = Inputs::get_switching_frequency(excitation);

        auto skinDepth = _windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
        double wireAirFillingFactor = Wire::get_filling_factor_round(2 * skinDepth);
        double windingWindowUtilizationFactor = wireAirFillingFactor * bobbinFillingFactor;
        double magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
        try {
            if (!_materialScaledMagneticFluxDensities.contains(core.get_material_name())) {
                auto coreLossesMethods = core.get_available_core_losses_methods();

                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
                    double referenceCoreLosses = _coreLossesModelSteinmetz->get_core_losses(core, _operatingPointExcitation, temperature).get_core_losses();
                    auto aux = _coreLossesModelSteinmetz->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                    magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
                }
                else {
                    double referenceCoreLosses = _coreLossesModelProprietary->get_core_losses(core, _operatingPointExcitation, temperature).get_core_losses();

                    auto aux = _coreLossesModelProprietary->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                    magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
                }
                _materialScaledMagneticFluxDensities[core.get_material_name()] = magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
            }
            else {
                magneticFluxDensityPeakAtFrequencyOfReferenceLosses = _materialScaledMagneticFluxDensities[core.get_material_name()];
            }
        }
        catch (...) {
            magneticFluxDensityPeakAtFrequencyOfReferenceLosses = _magneticFluxDensityReference;
        }
        double areaProductRequired = _areaProductRequiredPreCalculations[operatingPointIndex] / (windingWindowUtilizationFactor * magneticFluxDensityPeakAtFrequencyOfReferenceLosses);
        if (std::isnan(magneticFluxDensityPeakAtFrequencyOfReferenceLosses) || magneticFluxDensityPeakAtFrequencyOfReferenceLosses == 0) {
            throw NaNResultException("magneticFluxDensityPeakAtFrequencyOfReferenceLosses cannot be 0 or NaN");
            break;
        }
        if (std::isnan(areaProductRequired)) {
            break;
        }
        if (std::isinf(areaProductRequired) || areaProductRequired == 0) {
            throw NaNResultException("areaProductRequired cannot be 0 or NaN");
        }

        maximumAreaProductRequired = std::max(maximumAreaProductRequired, areaProductRequired);
    }
    if (maximumAreaProductRequired > 1) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "maximumAreaProductRequired cannot be larger than 1 (probably)");
    }

    bool valid = areaProductCore >= maximumAreaProductRequired * defaults.coreAdviserThresholdValidity;
    double scoring = fabs(areaProductCore - maximumAreaProductRequired);

    return {valid, scoring};
}

MagneticFilterEnergyStored::MagneticFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models) {
    _models = models;
    _magneticEnergy = MagneticEnergy(models);
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
}

MagneticFilterEnergyStored::MagneticFilterEnergyStored(Inputs inputs) {
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
}

std::pair<bool, double> MagneticFilterEnergyStored::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    bool valid = true;
    double totalStorableMagneticEnergy = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_point(operatingPointIndex);
        double storableEnergy = _magneticEnergy.calculate_core_maximum_magnetic_energy(magnetic->get_core(), operatingPoint);
        totalStorableMagneticEnergy = std::max(totalStorableMagneticEnergy, storableEnergy);

        if (totalStorableMagneticEnergy >= _requiredMagneticEnergy * defaults.coreAdviserThresholdValidity) {
            if (outputs != nullptr) {
                while (outputs->size() < operatingPointIndex + 1) {
                    outputs->push_back(Outputs());
                }
                MagnetizingInductanceOutput magnetizingInductanceOutput;
                magnetizingInductanceOutput.set_maximum_magnetic_energy_core(storableEnergy);
                magnetizingInductanceOutput.set_method_used(_models["gapReluctance"]);
                magnetizingInductanceOutput.set_origin(ResultOrigin::SIMULATION);
                InductanceOutput inductanceOutput;
                if ((*outputs)[operatingPointIndex].get_inductance()) {
                    inductanceOutput = *(*outputs)[operatingPointIndex].get_inductance();
                }
                inductanceOutput.set_magnetizing_inductance(magnetizingInductanceOutput);
                (*outputs)[operatingPointIndex].set_inductance(inductanceOutput);
            }
        }
        else {
            valid = false;
            break;
        }
    }

    return {valid, totalStorableMagneticEnergy};
}

MagneticFilterEstimatedCost::MagneticFilterEstimatedCost(Inputs inputs) {
    double primaryCurrentRms = 0;
    double frequency = 0;
    double temperature = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        primaryCurrentRms = std::max(primaryCurrentRms, Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_current().value().get_processed().value().get_rms().value());
        frequency = std::max(frequency, Inputs::get_switching_frequency(Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex))));
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }

    auto windingSkinEffectLossesModel = WindingSkinEffectLosses();
    _skinDepth = windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
    _wireAirFillingFactor = Wire::get_filling_factor_round(2 * _skinDepth);
    double estimatedWireConductingArea = std::numbers::pi * pow(_skinDepth, 2);
    _estimatedWireTotalArea = estimatedWireConductingArea / _wireAirFillingFactor;
    double necessaryWireCopperArea = primaryCurrentRms / defaults.maximumCurrentDensity;
    _estimatedParallels = ceil(necessaryWireCopperArea / estimatedWireConductingArea);

    if (settings.get_core_adviser_include_margin() && inputs.get_design_requirements().get_insulation()) {
        auto clearanceAndCreepageDistance = InsulationCoordinator().calculate_creepage_distance(inputs, true);
        _averageMarginInWindingWindow = clearanceAndCreepageDistance;
    }
}

std::pair<bool, double> MagneticFilterEstimatedCost::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    double primaryNumberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
    double estimatedNeededWindingArea = primaryNumberTurns * _estimatedParallels * _estimatedWireTotalArea * (inputs->get_design_requirements().get_turns_ratios().size() + 1);
    WindingWindowElement windingWindow;

    std::string shapeName = core.get_shape_name();
    if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
        auto bobbin = Bobbin::create_quick_bobbin(core);
        windingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
    }
    else {
        windingWindow = core.get_winding_windows()[0];

    }
    double windingWindowArea = windingWindow.get_area().value();
    if (_averageMarginInWindingWindow > 0) {
        if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
            if (windingWindow.get_height().value() > windingWindow.get_width().value()) {
                windingWindowArea -= windingWindow.get_width().value() * _averageMarginInWindingWindow;
            }
            else {
                windingWindowArea -= windingWindow.get_height().value() * _averageMarginInWindingWindow;
            }
        }
        else {
            auto radialHeight = windingWindow.get_radial_height().value();
            if (_averageMarginInWindingWindow > radialHeight / 2) {
                return {false, 0.0};
            }
            double wireAngle = wound_distance_to_angle(_averageMarginInWindingWindow, radialHeight / 2);
            if (std::isnan((wireAngle) / 360)) {
                throw NaNResultException("wireAngle: " + std::to_string(wireAngle));
            }
            windingWindowArea *= (wireAngle) / 360;
        }
    }

    bool valid = windingWindowArea >= estimatedNeededWindingArea * defaults.coreAdviserThresholdValidity;


    double manufacturabilityRelativeCost;
    if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
        double estimatedNeededLayers = (primaryNumberTurns * _estimatedParallels * (2 * _skinDepth / _wireAirFillingFactor)) / windingWindow.get_height().value();
        manufacturabilityRelativeCost = estimatedNeededLayers;
    }
    else {
        if (_skinDepth >= windingWindow.get_radial_height().value()) {
            return {false, 0.0};
        }
        double layerLength = 2 * std::numbers::pi * (windingWindow.get_radial_height().value() - _skinDepth);
        double estimatedNeededLayers = (primaryNumberTurns * _estimatedParallels * (2 * _skinDepth / _wireAirFillingFactor)) / layerLength;
        if (estimatedNeededLayers < 0) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "estimatedNeededLayers cannot be negative");
        }
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

    return {valid, manufacturabilityRelativeCost};
}

std::pair<bool, double> MagneticFilterCost::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    std::vector<Wire> wires = magnetic->get_mutable_coil().get_wires();
    auto wireRelativeCosts = 0;
    for (auto wire : wires) {
        wireRelativeCosts += wire.get_relative_cost();
    }
    auto numberLayers = magnetic->get_mutable_coil().get_layers_description()->size();
    double scoring = numberLayers + wireRelativeCosts;
    return {true, scoring};
}

MagneticFilterCoreAndDcLosses::MagneticFilterCoreAndDcLosses(Inputs inputs) {
    std::map<std::string, std::string> models;
    models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
    models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
    models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);
    MagneticFilterCoreAndDcLosses(inputs, models);
}

MagneticFilterCoreAndDcLosses::MagneticFilterCoreAndDcLosses() {
    std::map<std::string, std::string> models;
    models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
    models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
    models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);

    _maximumPowerMean = 0;

    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));

    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _models = models;
}

MagneticFilterCoreAndDcLosses::MagneticFilterCoreAndDcLosses(Inputs inputs, std::map<std::string, std::string> models) {
    bool largeWaveform = false;

    std::vector<double> powerMeans(inputs.get_operating_points().size(), 0);
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        auto voltageWaveform = Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_voltage().value().get_waveform().value();
        auto currentWaveform = Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_current().value().get_waveform().value();
        double frequency = Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency();

        if (voltageWaveform.get_data().size() != currentWaveform.get_data().size()) {
            voltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, frequency, std::max(voltageWaveform.get_data().size(), currentWaveform.get_data().size()));
            currentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, frequency, std::max(voltageWaveform.get_data().size(), currentWaveform.get_data().size()));
        }
        std::vector<double> voltageWaveformData = voltageWaveform.get_data();
        std::vector<double> currentWaveformData = currentWaveform.get_data();
        if (currentWaveformData.size() > settings.get_inputs_number_points_sampled_waveforms() * 2) {
            largeWaveform = true;
        }
        for (size_t i = 0; i < voltageWaveformData.size(); ++i)
        {
            powerMeans[operatingPointIndex] += fabs(voltageWaveformData[i] * currentWaveformData[i]);
        }
        powerMeans[operatingPointIndex] /= voltageWaveformData.size();
    }


    if (largeWaveform) {
        models["coreLosses"] = to_string(CoreLossesModels::STEINMETZ);
    }

    _maximumPowerMean = *max_element(powerMeans.begin(), powerMeans.end());

    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));

    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _models = models;
}

std::pair<bool, double> MagneticFilterCoreAndDcLosses::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    std::string shapeName = core.get_shape_name();
    if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
        auto bobbin = Bobbin::create_quick_bobbin(core);
        magnetic->get_mutable_coil().set_bobbin(bobbin);
        auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();

        if (windingWindows[0].get_width()) {
            if ((windingWindows[0].get_width().value() < 0) || (windingWindows[0].get_width().value() > 1)) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened in bobbins 1:   windingWindows[0].get_width(): " + std::to_string(static_cast<int>(bool(windingWindows[0].get_width()))) +
                                         " windingWindows[0].get_width().value(): " + std::to_string(windingWindows[0].get_width().value()) + 
                                         " shapeName: " + shapeName
                                         );
            }
        }
    }

    auto currentNumberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
    NumberTurns numberTurns(currentNumberTurns);
    std::vector<double> totalLossesPerOperatingPoint;
    std::vector<CoreLossesOutput> coreLossesPerOperatingPoint;
    std::vector<WindingLossesOutput> windingLossesPerOperatingPoint;
    double currentTotalLosses = DBL_MAX;
    double coreLosses = DBL_MAX;
    CoreLossesOutput coreLossesOutput;
    double ohmicLosses = DBL_MAX;
    WindingLossesOutput windingLossesOutput;
    windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
    double newTotalLosses = DBL_MAX;
    auto previousNumberTurnsPrimary = currentNumberTurns;

    size_t iteration = 10;

    Coil coil = magnetic->get_coil();

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_point(operatingPointIndex);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();
        OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
        size_t numberTimeouts = 0;
        do {
            currentTotalLosses = newTotalLosses;
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            coil.get_mutable_functional_description()[0].set_number_turns(numberTurnsCombination[0]);
            // coil = Coil(coil);
            settings.set_coil_delimit_and_compact(false);
            coil.fast_wind();

            auto [magnetizingInductance, magneticFluxDensity] = _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint);

            if (!check_requirement(inputs->get_design_requirements().get_magnetizing_inductance(), magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                if (resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance()) < resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                    coil.get_mutable_functional_description()[0].set_number_turns(previousNumberTurnsPrimary);
                    // coil = Coil(coil);
                    settings.set_coil_delimit_and_compact(false);
                    coil.fast_wind();
                    break;
                }
            }
            else {
                previousNumberTurnsPrimary = numberTurnsCombination[0];
            }

            if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
                if (!coil.get_turns_description()) {
                    newTotalLosses = coreLosses;
                    break;
                }
            }

            excitation.set_magnetic_flux_density(magneticFluxDensity);
            auto coreLossesMethods = core.get_available_core_losses_methods(); 
            if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
                coreLossesOutput = _coreLossesModelSteinmetz->get_core_losses(core, excitation, temperature);
                coreLosses = coreLossesOutput.get_core_losses();
            }
            else {
                coreLossesOutput = _coreLossesModelProprietary->get_core_losses(core, excitation, temperature);
                coreLosses = coreLossesOutput.get_core_losses();
                if (coreLosses < 0) {
                    break;
                }
            }

            if (coreLosses < 0) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in core losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value());
            }

            if (!coil.get_turns_description()) {
                newTotalLosses = coreLosses;
                break;
            }

            if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
                windingLossesOutput = _windingOhmicLosses.calculate_ohmic_losses(coil, operatingPoint, temperature);
                ohmicLosses = windingLossesOutput.get_winding_losses();
                newTotalLosses = coreLosses + ohmicLosses;
                if (ohmicLosses < 0) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happend in ohmic losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value() + " ohmicLosses: " + std::to_string(ohmicLosses));
                }
            }
            else {
                newTotalLosses = coreLosses;
                break;
            }

            if (newTotalLosses == DBL_MAX) {
                throw CalculationException(ErrorCode::CALCULATION_DIVERGED, "Too large losses");
            }

            iteration--;
            if (iteration <=0) {
                numberTimeouts++;
                break;
            }
        }
        while(newTotalLosses < currentTotalLosses * defaults.coreAdviserThresholdValidity);


        if (coreLosses < DBL_MAX && coreLosses > 0) {
            magnetic->set_coil(coil);

            currentTotalLosses = newTotalLosses;
            totalLossesPerOperatingPoint.push_back(currentTotalLosses);
            coreLossesPerOperatingPoint.push_back(coreLossesOutput);
            windingLossesPerOperatingPoint.push_back(windingLossesOutput);
        }
    }

    bool valid = false;
    double meanTotalLosses = 0;
    if (totalLossesPerOperatingPoint.size() < inputs->get_operating_points().size()) {
        return {false, 0.0};
    }
    else {

        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            meanTotalLosses += totalLossesPerOperatingPoint[operatingPointIndex];
        }
        if (meanTotalLosses > DBL_MAX / 2) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in core losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value());
        }
        meanTotalLosses /= inputs->get_operating_points().size();

        valid = meanTotalLosses < _maximumPowerMean * defaults.coreAdviserMaximumPercentagePowerCoreLosses / defaults.coreAdviserThresholdValidity;
    }

    if (outputs != nullptr) {
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            (*outputs)[operatingPointIndex].set_core_losses(coreLossesPerOperatingPoint[operatingPointIndex]);
            (*outputs)[operatingPointIndex].set_winding_losses(windingLossesPerOperatingPoint[operatingPointIndex]);
        }
    }

    return {valid, meanTotalLosses};
}

MagneticFilterCoreDcAndSkinLosses::MagneticFilterCoreDcAndSkinLosses(Inputs inputs) {
    std::map<std::string, std::string> models;
    models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
    models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
    models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);
    MagneticFilterCoreDcAndSkinLosses(inputs, models);
}

MagneticFilterCoreDcAndSkinLosses::MagneticFilterCoreDcAndSkinLosses() {
    std::map<std::string, std::string> models;
    models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
    models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
    models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);

    _maximumPowerMean = 0;

    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));

    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _models = models;
}

MagneticFilterCoreDcAndSkinLosses::MagneticFilterCoreDcAndSkinLosses(Inputs inputs, std::map<std::string, std::string> models) {
    bool largeWaveform = false;

    std::vector<double> powerMeans(inputs.get_operating_points().size(), 0);
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        auto voltageWaveform = Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_voltage().value().get_waveform().value();
        auto currentWaveform = Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_current().value().get_waveform().value();
        double frequency = Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency();

        if (voltageWaveform.get_data().size() != currentWaveform.get_data().size()) {
            voltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, frequency, std::max(voltageWaveform.get_data().size(), currentWaveform.get_data().size()));
            currentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, frequency, std::max(voltageWaveform.get_data().size(), currentWaveform.get_data().size()));
        }
        std::vector<double> voltageWaveformData = voltageWaveform.get_data();
        std::vector<double> currentWaveformData = currentWaveform.get_data();
        if (currentWaveformData.size() > settings.get_inputs_number_points_sampled_waveforms() * 2) {
            largeWaveform = true;
        }
        for (size_t i = 0; i < voltageWaveformData.size(); ++i)
        {
            powerMeans[operatingPointIndex] += fabs(voltageWaveformData[i] * currentWaveformData[i]);
        }
        powerMeans[operatingPointIndex] /= voltageWaveformData.size();
    }


    if (largeWaveform) {
        models["coreLosses"] = to_string(CoreLossesModels::STEINMETZ);
    }

    _maximumPowerMean = *max_element(powerMeans.begin(), powerMeans.end());

    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));

    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _models = models;
}

std::pair<bool, double> MagneticFilterCoreDcAndSkinLosses::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    std::string shapeName = core.get_shape_name();
    if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
        auto bobbin = Bobbin::create_quick_bobbin(core);
        magnetic->get_mutable_coil().set_bobbin(bobbin);
        auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();

        if (windingWindows[0].get_width()) {
            if ((windingWindows[0].get_width().value() < 0) || (windingWindows[0].get_width().value() > 1)) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened in bobbins 1:   windingWindows[0].get_width(): " + std::to_string(static_cast<int>(bool(windingWindows[0].get_width()))) +
                                         " windingWindows[0].get_width().value(): " + std::to_string(windingWindows[0].get_width().value()) + 
                                         " shapeName: " + shapeName
                                         );
            }
        }
    }

    auto currentNumberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
    NumberTurns numberTurns(currentNumberTurns);
    std::vector<double> totalLossesPerOperatingPoint;
    std::vector<CoreLossesOutput> coreLossesPerOperatingPoint;
    std::vector<WindingLossesOutput> windingLossesPerOperatingPoint;
    double currentTotalLosses = DBL_MAX;
    double coreLosses = DBL_MAX;
    CoreLossesOutput coreLossesOutput;
    double ohmicAndSkinEffectLosses = DBL_MAX;
    WindingLossesOutput windingLossesOutput;
    windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
    double newTotalLosses = DBL_MAX;
    auto previousNumberTurnsPrimary = currentNumberTurns;

    size_t iteration = 10;

    Coil coil = magnetic->get_coil();

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_point(operatingPointIndex);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();
        OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
        size_t numberTimeouts = 0;
        do {
            currentTotalLosses = newTotalLosses;
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            coil.get_mutable_functional_description()[0].set_number_turns(numberTurnsCombination[0]);
            // coil = Coil(coil);
            settings.set_coil_delimit_and_compact(false);
            coil.fast_wind();

            auto [magnetizingInductance, magneticFluxDensity] = _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint);

            if (!check_requirement(inputs->get_design_requirements().get_magnetizing_inductance(), magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                if (resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance()) < resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                    coil.get_mutable_functional_description()[0].set_number_turns(previousNumberTurnsPrimary);
                    // coil = Coil(coil);
                    settings.set_coil_delimit_and_compact(false);
                    coil.fast_wind();
                    break;
                }
            }
            else {
                previousNumberTurnsPrimary = numberTurnsCombination[0];
            }

            if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
                if (!coil.get_turns_description()) {
                    newTotalLosses = coreLosses;
                    break;
                }
            }

            excitation.set_magnetic_flux_density(magneticFluxDensity);
            auto coreLossesMethods = core.get_available_core_losses_methods(); 
            if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
                coreLossesOutput = _coreLossesModelSteinmetz->get_core_losses(core, excitation, temperature);
                coreLosses = coreLossesOutput.get_core_losses();
            }
            else {
                coreLossesOutput = _coreLossesModelProprietary->get_core_losses(core, excitation, temperature);
                coreLosses = coreLossesOutput.get_core_losses();
                if (coreLosses < 0) {
                    break;
                }
            }

            if (coreLosses < 0) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in core losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value());
            }

            if (!coil.get_turns_description()) {
                newTotalLosses = coreLosses;
                break;
            }

            if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
                windingLossesOutput = _windingOhmicLosses.calculate_ohmic_losses(coil, operatingPoint, temperature);
                windingLossesOutput = _windingSkinEffectLosses.calculate_skin_effect_losses(coil, temperature, windingLossesOutput, settings.get_harmonic_amplitude_threshold());

                ohmicAndSkinEffectLosses = windingLossesOutput.get_winding_losses();
                newTotalLosses = coreLosses + ohmicAndSkinEffectLosses;
                if (ohmicAndSkinEffectLosses < 0) {
                    throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in ohmic losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value() + " ohmicAndSkinEffectLosses: " + std::to_string(ohmicAndSkinEffectLosses));
                }
            }
            else {
                newTotalLosses = coreLosses;
                break;
            }

            if (newTotalLosses == DBL_MAX) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Too large losses");
            }

            iteration--;
            if (iteration <=0) {
                numberTimeouts++;
                break;
            }
        }
        while(newTotalLosses < currentTotalLosses * defaults.coreAdviserThresholdValidity);


        if (coreLosses < DBL_MAX && coreLosses > 0) {
            magnetic->set_coil(coil);

            currentTotalLosses = newTotalLosses;
            totalLossesPerOperatingPoint.push_back(currentTotalLosses);
            coreLossesPerOperatingPoint.push_back(coreLossesOutput);
            windingLossesPerOperatingPoint.push_back(windingLossesOutput);
        }
    }

    bool valid = false;
    double meanTotalLosses = 0;
    if (totalLossesPerOperatingPoint.size() < inputs->get_operating_points().size()) {
        return {false, 0.0};
    }
    else {

        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            meanTotalLosses += totalLossesPerOperatingPoint[operatingPointIndex];
        }
        if (meanTotalLosses > DBL_MAX / 2) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in core losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value());
        }
        meanTotalLosses /= inputs->get_operating_points().size();

        valid = meanTotalLosses < _maximumPowerMean * defaults.coreAdviserMaximumPercentagePowerCoreLosses / defaults.coreAdviserThresholdValidity;
    }

    if (outputs != nullptr) {
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            (*outputs)[operatingPointIndex].set_core_losses(coreLossesPerOperatingPoint[operatingPointIndex]);
            (*outputs)[operatingPointIndex].set_winding_losses(windingLossesPerOperatingPoint[operatingPointIndex]);
        }
    }

    return {valid, meanTotalLosses};
}

MagneticFilterLosses::MagneticFilterLosses(std::map<std::string, std::string> models) {
    _models = models;
}

std::pair<bool, double> MagneticFilterLosses::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_points()[operatingPointIndex];
        auto temperature = operatingPoint.get_conditions().get_ambient_temperature();
        auto windingLosses = _magneticSimulator.calculate_winding_losses(operatingPoint, *magnetic, temperature);
        auto windingLossesValue = windingLosses.get_winding_losses();
        auto coreLosses = _magneticSimulator.calculate_core_losses(operatingPoint, *magnetic);
        auto coreLossesValue = coreLosses.get_core_losses();
        scoring += windingLossesValue + coreLossesValue;

        if (outputs != nullptr) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            (*outputs)[operatingPointIndex].set_core_losses(coreLosses);
            (*outputs)[operatingPointIndex].set_winding_losses(windingLosses);
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {true, scoring};
}

MagneticFilterLossesNoProximity::MagneticFilterLossesNoProximity(std::map<std::string, std::string> models) {
    _models = models;
}

std::pair<bool, double> MagneticFilterLossesNoProximity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_points()[operatingPointIndex];
        auto temperature = operatingPoint.get_conditions().get_ambient_temperature();
        auto windingLosses = _windingOhmicLosses.calculate_ohmic_losses(magnetic->get_coil(), operatingPoint, temperature);
        windingLosses = _windingSkinEffectLosses.calculate_skin_effect_losses(magnetic->get_coil(), temperature, windingLosses, 0.5);
        windingLosses = _windingSkinEffectLosses.calculate_skin_effect_losses(magnetic->get_coil(), temperature, windingLosses, 0.5);
        // windingLosses = WindingLosses::combine_turn_losses(windingLosses, magnetic->get_coil());
        double windingLossesValue = windingLosses.get_winding_losses();

        auto coreLosses = _magneticSimulator.calculate_core_losses(operatingPoint, *magnetic);
        auto coreLossesValue = coreLosses.get_core_losses();
        scoring += windingLossesValue + coreLossesValue;

        if (outputs != nullptr) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            (*outputs)[operatingPointIndex].set_core_losses(coreLosses);
            (*outputs)[operatingPointIndex].set_winding_losses(windingLosses);
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {true, scoring};
}

std::pair<bool, double> MagneticFilterDimensions::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    auto depth = core.get_depth();

    if (magnetic->get_coil().get_layers_description()) {
        auto coilDepth = magnetic->get_mutable_core().get_columns()[0].get_depth();
        auto layers = magnetic->get_coil().get_layers_description().value();
        for (auto layer : layers) {
            coilDepth += layer.get_dimensions()[0] * 2;
        }
        depth = std::max(depth, coilDepth);
    }

    double volume = core.get_width() * core.get_height() * depth;

    return {true, volume};
}

std::pair<bool, double> MagneticFilterCoreMinimumImpedance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    double primaryCurrentRms = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        primaryCurrentRms = std::max(primaryCurrentRms, Inputs::get_primary_excitation(inputs->get_operating_point(operatingPointIndex)).get_current().value().get_processed().value().get_rms().value());
    }

    std::string shapeName = core.get_shape_name();
    if (!((shapeName.rfind("PQI", 0) == 0) || (shapeName.rfind("UI ", 0) == 0))) {
        auto bobbin = Bobbin::create_quick_bobbin(core);
        magnetic->get_mutable_coil().set_bobbin(bobbin);
        auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();

        if (windingWindows[0].get_width()) {
            if ((windingWindows[0].get_width().value() < 0) || (windingWindows[0].get_width().value() > 1)) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened in bobbins 1:   windingWindows[0].get_width(): " + std::to_string(static_cast<int>(bool(windingWindows[0].get_width()))) +
                                         " windingWindows[0].get_width().value(): " + std::to_string(windingWindows[0].get_width().value()) + 
                                         " shapeName: " + shapeName
                                         );
            }
        }
    }

    auto currentNumberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
    NumberTurns numberTurns(currentNumberTurns);

    Coil coil = magnetic->get_coil();

    double conductingArea = primaryCurrentRms / defaults.maximumCurrentDensity;
    auto wire = Wire::get_wire_for_conducting_area(conductingArea, defaults.ambientTemperature, false);
    coil.get_mutable_functional_description()[0].set_wire(wire);
    coil.unwind();

    if (!inputs->get_design_requirements().get_minimum_impedance()) {
        throw InvalidInputException("Minimum impedance missing from requirements");
    }

    auto minimumImpedanceRequirement = inputs->get_design_requirements().get_minimum_impedance().value();
    auto windingWindowArea = magnetic->get_mutable_coil().resolve_bobbin().get_winding_window_area();

    bool validDesign = true;
    bool validMaterial = true;
    double totalImpedanceExtra;
    int timeout = 100;
    do {
        totalImpedanceExtra = 0;
        validDesign = true;
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();

        if (numberTurnsCombination[0] * std::numbers::pi * pow(resolve_dimensional_values(wire.get_outer_diameter().value()) / 2, 2) >= windingWindowArea) {
            validMaterial = false;
            break;
        }
        coil.get_mutable_functional_description()[0].set_number_turns(numberTurnsCombination[0]);
        auto selfResonantFrequency = _impedanceModel.calculate_self_resonant_frequency(core, coil);

        for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
            auto frequency = impedanceAtFrequency.get_frequency();
            if (frequency > 0.25 * selfResonantFrequency) {  // hardcoded 20% of SRF
                validDesign = false;
                break;
            }
        }

        if (!validDesign) {
            break;
        }

        for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
            auto frequency = impedanceAtFrequency.get_frequency();
            auto minimumImpedanceRequired = impedanceAtFrequency.get_impedance();
            try {
                auto impedance = abs(_impedanceModel.calculate_impedance(core, coil, frequency));
                if (impedance < minimumImpedanceRequired.get_magnitude()) {
                    validDesign = false;
                    break;
                }
                else {
                    totalImpedanceExtra += (impedance - minimumImpedanceRequired.get_magnitude());
                }

            }
            catch (const ModelNotAvailableException &exc) {
                validMaterial = false;
            }
        }

        timeout--;
    }
    while(!validDesign && validMaterial && timeout > 0);



    if (validDesign && validMaterial) {
        coil.fast_wind();

    }

    bool valid = coil.are_sections_and_layers_fitting() && coil.get_turns_description();
    
    magnetic->set_coil(std::move(coil));

    return {valid, totalImpedanceExtra};
}

MagneticFilterAreaNoParallels::MagneticFilterAreaNoParallels(int maximumNumberParallels) {
    _maximumNumberParallels = maximumNumberParallels;
}

std::pair<bool, double> MagneticFilterAreaNoParallels::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    for (auto winding : magnetic->get_coil().get_functional_description()) {
        auto section = magnetic->get_mutable_coil().get_sections_by_winding(winding.get_name())[0];
        auto [auxValid, auxScoring] = evaluate_magnetic(winding, section);
        valid &= auxValid;
        scoring += auxScoring;
    }
    scoring /= magnetic->get_coil().get_functional_description().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterAreaNoParallels::evaluate_magnetic(Winding winding, Section section) {
    auto wire = Coil::resolve_wire(winding);

    if (wire.get_type() == WireType::FOIL && winding.get_number_parallels() * winding.get_number_turns() > _maximumNumberParallels) {
        return {false, 0.0};
    }

    if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        if (wire.get_maximum_outer_width() < section.get_dimensions()[0] && wire.get_maximum_outer_height() < section.get_dimensions()[1]) {
            return {true, 0.0};
        }
        else {
            return {false, 0.0};
        }
    }
    else {
        double wireAngle = wound_distance_to_angle(wire.get_maximum_outer_height(), wire.get_maximum_outer_width());

        if (wire.get_maximum_outer_width() < section.get_dimensions()[0] && wireAngle < section.get_dimensions()[1]) {
            return {true, 0.0};
        }
        else {
            return {false, 0.0};
        }
    }
}

std::pair<bool, double> MagneticFilterAreaWithParallels::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    for (auto winding : magnetic->get_coil().get_functional_description()) {
        auto sections = magnetic->get_mutable_coil().get_sections_by_winding(winding.get_name());
        auto section = sections[0];
        double sectionArea;
        if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
            sectionArea = section.get_dimensions()[0] * section.get_dimensions()[1];
        }
        else {
            sectionArea = std::numbers::pi * pow(section.get_dimensions()[0], 2) * section.get_dimensions()[1] / 360;
        }
        auto [auxValid, auxScoring] = evaluate_magnetic(winding, section, sections.size(), sectionArea, false);
        valid &= auxValid;
        scoring += auxScoring;
    }
    scoring /= magnetic->get_coil().get_functional_description().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterAreaWithParallels::evaluate_magnetic(Winding winding, Section section, double numberSections, double sectionArea, bool allowNotFit) {
    auto wire = Coil::resolve_wire(winding);
    if (!Coil::resolve_wire(winding).get_conducting_area()) {
        throw CoilNotProcessedException("Conducting area is missing");
    }
    auto neededOuterAreaNoCompact = wire.get_maximum_outer_width() * wire.get_maximum_outer_height();

    neededOuterAreaNoCompact *= winding.get_number_parallels() * winding.get_number_turns() / numberSections;

    if (neededOuterAreaNoCompact < sectionArea) {
        // double scoring = (section.get_dimensions()[0] * section.get_dimensions()[1]) - neededOuterAreaNoCompact;
        return {true, 1.0};
    }
    else if (allowNotFit) {
        double extra = (neededOuterAreaNoCompact - sectionArea) / sectionArea;
        if (extra < 0.5) {
            return {true, extra};
        }
        else {
            return {false, 0.0};
        }
    }
    else {
        return {false, 0.0};
    }
}

std::pair<bool, double> MagneticFilterEffectiveResistance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t windingIndex = 0; windingIndex < magnetic->get_coil().get_functional_description().size(); ++windingIndex) {
        auto winding = magnetic->get_coil().get_functional_description()[windingIndex];
        auto maximumEffectiveFrequency = inputs->get_maximum_current_effective_frequency(windingIndex);
        auto temperature = inputs->get_maximum_temperature();
        auto [auxValid, auxScoring] = evaluate_magnetic(winding, maximumEffectiveFrequency, temperature);
        valid &= auxValid;
        scoring += auxScoring;
    }
    scoring /= magnetic->get_coil().get_functional_description().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterEffectiveResistance::evaluate_magnetic(Winding winding, double effectivefrequency, double temperature) {
    auto wire = Coil::resolve_wire(winding);

    double effectiveResistancePerMeter = WindingLosses::calculate_effective_resistance_per_meter(wire, effectivefrequency, temperature);

    double valid = true;
    // double valid = effectiveResistancePerMeter < defaults.maximumEffectiveCurrentDensity;

    return {valid, effectiveResistancePerMeter};
}

std::pair<bool, double> MagneticFilterProximityFactor::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t windingIndex = 0; windingIndex < magnetic->get_coil().get_functional_description().size(); ++windingIndex) {
        auto winding = magnetic->get_coil().get_functional_description()[windingIndex];
        auto maximumEffectiveFrequency = inputs->get_maximum_current_effective_frequency(windingIndex);
        auto temperature = inputs->get_maximum_temperature();
        double effectiveSkinDepth = WindingSkinEffectLosses::calculate_skin_depth("copper", maximumEffectiveFrequency, temperature);
        auto [auxValid, auxScoring] = evaluate_magnetic(winding, effectiveSkinDepth, temperature);
        valid &= auxValid;
        scoring += auxScoring;
    }
    scoring /= magnetic->get_coil().get_functional_description().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterProximityFactor::evaluate_magnetic(Winding winding, double effectiveSkinDepth, double temperature) {
    auto wire = Coil::resolve_wire(winding);

    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }
    double proximityFactor = wire.get_minimum_conducting_dimension() / effectiveSkinDepth * pow(wire.get_number_conductors().value() * winding.get_number_parallels() * winding.get_number_turns(), 2);
    // proximityFactor = wire.get_minimum_conducting_dimension() / effectiveSkinDepth * pow(winding.get_number_parallels() / (std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height())), 2);

    double valid = true;
    // double valid = effectiveResistancePerMeter < defaults.maximumEffectiveCurrentDensity;

    return {valid, proximityFactor};
}

std::pair<bool, double> MagneticFilterSolidInsulationRequirements::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = false;
    double scoring = 0;
    auto core = magnetic->get_core();
    auto coreType = core.get_functional_description().get_type();
    auto patterns = Coil::get_patterns(*inputs, coreType);
    auto repetitions = Coil::get_repetitions(*inputs, coreType);

    for (auto repetition : repetitions) {
        for (auto pattern : patterns) {
            auto aux = magnetic->get_mutable_coil().check_pattern_and_repetitions_integrity(pattern, repetition);
            pattern = aux.first;
            repetition = aux.second;
            auto combinationsSolidInsulationRequirementsForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(*inputs, pattern, repetition);
            for(size_t insulationIndex = 0; insulationIndex < combinationsSolidInsulationRequirementsForWires.size(); ++insulationIndex) {
                auto solidInsulationRequirementsForWires = combinationsSolidInsulationRequirementsForWires[insulationIndex];
                bool combinationValid = true;
                double combinationScoring = 0;
                for (size_t windingIndex = 0; windingIndex < magnetic->get_coil().get_functional_description().size(); ++windingIndex) {
                    auto winding = magnetic->get_coil().get_functional_description()[windingIndex];
                    auto [auxValid, auxScoring] = evaluate_magnetic(winding, solidInsulationRequirementsForWires[windingIndex]);
                    combinationValid &= auxValid;
                    combinationScoring += auxScoring;
                }

                valid |= combinationValid;
                if (valid) {
                    scoring = std::max(scoring, combinationScoring);
                    return {valid, scoring};
                }
            }
        }
    }
    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterSolidInsulationRequirements::evaluate_magnetic(Winding winding, WireSolidInsulationRequirements wireSolidInsulationRequirements) {
    auto wire = Coil::resolve_wire(winding);

    if (wire.get_type() == WireType::FOIL || wire.get_type() == WireType::PLANAR) {
        return {true, 0.0};
    }

    if (!wire.resolve_coating()) {
        return {false, 0.0};
    }

    auto coating = wire.resolve_coating().value();

    if (wire.get_type() == WireType::LITZ) {
        auto strand = wire.resolve_strand();
        coating = Wire::resolve_coating(strand).value();
    }

    if (!coating.get_breakdown_voltage()) {
        return {false, 0.0};
        // throw std::runtime_error("Wire " + wire.get_name().value() + " is missing breakdown voltage");
    }

    if (coating.get_breakdown_voltage().value() < wireSolidInsulationRequirements.get_minimum_breakdown_voltage()) {
        return {false, 0.0};
    }

    if (wireSolidInsulationRequirements.get_minimum_grade() && coating.get_grade()) {
        if (coating.get_grade().value() < wireSolidInsulationRequirements.get_minimum_grade().value()) {
            return {false, 0.0};
        }
    }
    else if (wireSolidInsulationRequirements.get_minimum_number_layers() && coating.get_number_layers()) {
        if (coating.get_number_layers().value() < wireSolidInsulationRequirements.get_minimum_number_layers().value()) {
            return {false, 0.0};
        }
    }
    else if (wireSolidInsulationRequirements.get_minimum_number_layers() || wireSolidInsulationRequirements.get_minimum_grade()) {
        return {false, 0.0};
    }
    
    if (wireSolidInsulationRequirements.get_maximum_grade() && coating.get_grade()) {
        if (coating.get_grade().value() > wireSolidInsulationRequirements.get_maximum_grade().value()) {
            return {false, 0.0};
        }
    }
    else if (wireSolidInsulationRequirements.get_maximum_number_layers() && coating.get_number_layers()) {
        if (coating.get_number_layers().value() > wireSolidInsulationRequirements.get_maximum_number_layers().value()) {
            return {false, 0.0};
        }
    }
    else if (wireSolidInsulationRequirements.get_maximum_number_layers() || wireSolidInsulationRequirements.get_maximum_grade()) {
        return {false, 0.0};
    }

    double scoring = 0;
    if (wireSolidInsulationRequirements.get_minimum_breakdown_voltage() > 0) {
        scoring = wireSolidInsulationRequirements.get_minimum_breakdown_voltage() - coating.get_breakdown_voltage().value();
    }

    return {true, scoring};
}

std::pair<bool, double> MagneticFilterTurnsRatios::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    if (inputs->get_design_requirements().get_turns_ratios().size() > 0) {
        auto magneticTurnsRatios = magnetic->get_turns_ratios();
        if (magneticTurnsRatios.size() != inputs->get_design_requirements().get_turns_ratios().size()) {
            return {false, 0.0};
        }
        for (size_t i = 0; i < inputs->get_design_requirements().get_turns_ratios().size(); ++i) {
            auto turnsRatioRequirement = inputs->get_design_requirements().get_turns_ratios()[i];
            // TODO: Try all combinations of turns ratios, not just the default order
            if (!check_requirement(turnsRatioRequirement, magneticTurnsRatios[i])) {
                return {false, 0.0};
            }
            scoring += abs(resolve_dimensional_values(turnsRatioRequirement) - resolve_dimensional_values(magneticTurnsRatios[i]));
        }
    }
    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterMaximumDimensions::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    if (inputs->get_design_requirements().get_maximum_dimensions()) {
        auto maximumDimensions = inputs->get_design_requirements().get_maximum_dimensions().value();
        auto magneticDimensions = magnetic->get_maximum_dimensions();
        scoring = sqrt(pow(maximumDimensions.get_width().value() - magneticDimensions[0], 2) + pow(maximumDimensions.get_height().value() - magneticDimensions[1], 2)+ pow(maximumDimensions.get_depth().value() - magneticDimensions[2], 2));
        if (!magnetic->fits(maximumDimensions, true)) {
            valid = false;
        }
    }
    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterSaturation::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    for (auto operatingPoint : inputs->get_operating_points()) {
        OpenMagnetics::MagnetizingInductance magnetizingInductanceObj;
        auto magneticFluxDensity = magnetizingInductanceObj.calculate_inductance_and_magnetic_flux_density(magnetic->get_core(), magnetic->get_coil(), &operatingPoint).second;
        auto magneticFluxDensityPeak = magneticFluxDensity.get_processed().value().get_peak().value();

        auto magneticFluxDensitySaturation = magnetic->get_mutable_core().get_magnetic_flux_density_saturation(operatingPoint.get_conditions().get_ambient_temperature());
        scoring += fabs(magneticFluxDensitySaturation - magneticFluxDensityPeak);
        if (magneticFluxDensityPeak > magneticFluxDensitySaturation) {
            return {false, 0.0};
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterDcCurrentDensity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (auto operatingPoint : inputs->get_operating_points()) {
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) { 
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
            }
            auto current = excitation.get_current().value();
            auto wire = magnetic->get_mutable_coil().resolve_wire(windingIndex);
            auto dcCurrentDensity = wire.calculate_dc_current_density(current) / magnetic->get_mutable_coil().get_functional_description()[windingIndex].get_number_parallels();

            scoring += fabs(defaults.maximumCurrentDensity - dcCurrentDensity);
            if (dcCurrentDensity > defaults.maximumCurrentDensity) {
                return {false, 0.0};
            }
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterEffectiveCurrentDensity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (auto operatingPoint : inputs->get_operating_points()) {
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
            }
            auto current = excitation.get_current().value();
            auto wire = magnetic->get_mutable_coil().resolve_wire(windingIndex);
            auto effectiveCurrentDensity = wire.calculate_effective_current_density(current, operatingPoint.get_conditions().get_ambient_temperature()) / magnetic->get_mutable_coil().get_functional_description()[windingIndex].get_number_parallels();

            scoring += fabs(defaults.maximumEffectiveCurrentDensity - effectiveCurrentDensity);
            if (effectiveCurrentDensity > defaults.maximumEffectiveCurrentDensity) {
                return {false, 0.0};
            }
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterImpedance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_design_requirements().get_minimum_impedance()) {
        auto impedanceRequirement = inputs->get_design_requirements().get_minimum_impedance().value();
        for (auto impedanceAtFrequency : impedanceRequirement) {
            auto impedance = OpenMagnetics::Impedance().calculate_impedance(*magnetic, impedanceAtFrequency.get_frequency());
            scoring += fabs(impedanceAtFrequency.get_impedance().get_magnitude() - abs(impedance));

            if (impedanceAtFrequency.get_impedance().get_magnitude() > abs(impedance)) {
                valid = false;
            }
        }
        scoring /= impedanceRequirement.size();
    }
    if (inputs->get_operating_points().size() > 0) {
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            auto operatingPoint = inputs->get_operating_points()[operatingPointIndex];
            auto impedance = OpenMagnetics::Impedance().calculate_impedance(*magnetic, operatingPoint.get_excitations_per_winding()[0].get_frequency());
            scoring += 1.0 / abs(impedance);
            std::string name = magnetic->get_coil().get_functional_description()[0].get_name();
                
            if (outputs != nullptr) {
                while (outputs->size() < operatingPointIndex + 1) {
                    outputs->push_back(Outputs());
                }
                ImpedanceOutput impedanceOutput;
                ComplexMatrixAtFrequency complexMatrixAtFrequency;
                complexMatrixAtFrequency.set_frequency(operatingPoint.get_excitations_per_winding()[0].get_frequency());
                complexMatrixAtFrequency.get_mutable_magnitude()[name][name].set_nominal(abs(impedance));
                std::vector<ComplexMatrixAtFrequency> impedanceMatrixPerFrequency;
                impedanceMatrixPerFrequency.push_back(complexMatrixAtFrequency);
                impedanceOutput.set_impedance_matrix(impedanceMatrixPerFrequency);
                impedanceOutput.set_origin(ResultOrigin::SIMULATION);
                (*outputs)[operatingPointIndex].set_impedance(impedanceOutput);
            }
        }
    }
    else {
        auto impedance = OpenMagnetics::Impedance().calculate_impedance(*magnetic, defaults.measurementFrequency);
        scoring += 1.0 / abs(impedance);
    }

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterMagnetizingInductance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_points()[operatingPointIndex];
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");

        auto aux = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic->get_mutable_core(), magnetic->get_mutable_coil(), &operatingPoint);
        double magnetizingInductance = resolve_dimensional_values(aux.get_magnetizing_inductance());
        scoring += fabs(resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance()) - magnetizingInductance);

        if (!check_requirement(inputs->get_design_requirements().get_magnetizing_inductance(), magnetizingInductance)) {
            valid = false;
        }
        else {
            if (outputs != nullptr) {
                while (outputs->size() < operatingPointIndex + 1) {
                    outputs->push_back(Outputs());
                }
                InductanceOutput inductanceOutput;
                if ((*outputs)[operatingPointIndex].get_inductance()) {
                    inductanceOutput = *(*outputs)[operatingPointIndex].get_inductance();
                }
                inductanceOutput.set_magnetizing_inductance(aux);
                (*outputs)[operatingPointIndex].set_inductance(inductanceOutput);
            }
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

MagneticFilterFringingFactor::MagneticFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models) {
    _models = models;
    _magneticEnergy = MagneticEnergy(models);
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
    ReluctanceModels reluctanceModel;
    from_json(_models["gapReluctance"], reluctanceModel);
    _reluctanceModel = ReluctanceModel::factory(reluctanceModel);
}

MagneticFilterFringingFactor::MagneticFilterFringingFactor(Inputs inputs) {
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
    _reluctanceModel = ReluctanceModel::factory();
}

std::pair<bool, double> MagneticFilterFringingFactor::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    if (core.get_shape_family() == CoreShapeFamily::T) {
        return {true, 1};
    }
    else if (core.get_gapping().size() == 0) {
        return {true, 1};
    }
    else {
        double maximumGapLength = _reluctanceModel->get_gapping_by_fringing_factor(core, _fringingFactorLitmit);
        double gapLength = core.get_gapping()[0].get_length();
        if (gapLength > maximumGapLength) {
            return {false, 1};
        }
        else {
            return {true, 1};
        }
    }
}

std::pair<bool, double> MagneticFilterSkinLossesDensity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    auto temperature = inputs->get_maximum_temperature();

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (auto operatingPoint : inputs->get_operating_points()) {
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
            }
            auto current = excitation.get_current().value();
            auto winding = magnetic->get_coil().get_functional_description()[windingIndex];
            auto [auxValid, auxScoring] = evaluate_magnetic(winding, current, temperature);
            valid &= auxValid;
            scoring += auxScoring;
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterSkinLossesDensity::evaluate_magnetic(Winding winding, SignalDescriptor current, double temperature) {
    auto wire = Coil::resolve_wire(winding);

    double skinEffectLossesPerMeter = WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature).first;
    double valid = true;
    return {valid, skinEffectLossesPerMeter};
}

std::pair<bool, double> MagneticFilterVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double volume = maximumDimensions[0] * maximumDimensions[1] * maximumDimensions[2];
    return {true, volume};
}

std::pair<bool, double> MagneticFilterArea::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double area = maximumDimensions[0] * maximumDimensions[2];
    return {true, area};
}

std::pair<bool, double> MagneticFilterHeight::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double height = maximumDimensions[1];
    return {true, height};
}

std::pair<bool, double> MagneticFilterTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto previousLosses = get_scoring(magnetic->get_reference(), MagneticFilters::LOSSES_NO_PROXIMITY);
    double losses = 0;

    if (previousLosses) {
        losses = previousLosses.value();
    }
    else {
        auto aux = _magneticFilterLossesNoProximity.evaluate_magnetic(magnetic, inputs, outputs);
        losses = aux.second;
    }

    auto coreTemperatureModel = CoreTemperatureModel::factory(defaults.coreTemperatureModelDefault);

    auto coreTemperature = coreTemperatureModel->get_core_temperature(magnetic->get_core(), losses, defaults.ambientTemperature);
    double calculatedTemperature = coreTemperature.get_maximum_temperature();

    return {true, calculatedTemperature};
}

std::pair<bool, double> MagneticFilterLossesTimesVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto previousLosses = get_scoring(magnetic->get_reference(), MagneticFilters::LOSSES);
    double losses = 0;
    if (previousLosses) {
        losses = previousLosses.value();
    }
    else {
        auto aux = MagneticFilterLosses().evaluate_magnetic(magnetic, inputs, outputs);
        losses = aux.second;
    }
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring};
}
 
std::pair<bool, double> MagneticFilterVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto previousTemperatureRise = get_scoring(magnetic->get_reference(), MagneticFilters::TEMPERATURE_RISE);
    double temperature = 0;
    if (previousTemperatureRise) {
        temperature = previousTemperatureRise.value();
    }
    else {
        auto aux = _magneticFilterTemperatureRise.evaluate_magnetic(magnetic, inputs, outputs);
        temperature = aux.second;
    }

    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, volumeScoring * temperature};
}

std::pair<bool, double> MagneticFilterLossesTimesVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto previousLosses = get_scoring(magnetic->get_reference(), MagneticFilters::LOSSES);
    double losses = 0;
    if (previousLosses) {
        losses = previousLosses.value();
    }
    else {
        auto aux = MagneticFilterLosses().evaluate_magnetic(magnetic, inputs, outputs);
        losses = aux.second;
    }
    auto previousTemperatureRise = get_scoring(magnetic->get_reference(), MagneticFilters::TEMPERATURE_RISE);
    double temperature = 0;
    if (previousTemperatureRise) {
        temperature = previousTemperatureRise.value();
    }
    else {
        auto aux = _magneticFilterTemperatureRise.evaluate_magnetic(magnetic, inputs, outputs);
        temperature = aux.second;
    }

    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring * temperature};
}

std::pair<bool, double> MagneticFilterLossesNoProximityTimesVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto previousLosses = get_scoring(magnetic->get_reference(), MagneticFilters::LOSSES_NO_PROXIMITY);
    double losses = 0;
    if (previousLosses) {
        losses = previousLosses.value();
    }
    else {
        auto aux = _magneticFilterLossesNoProximity.evaluate_magnetic(magnetic, inputs, outputs);
        losses = aux.second;
    }
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring};
}

std::pair<bool, double> MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto previousLosses = get_scoring(magnetic->get_reference(), MagneticFilters::LOSSES_NO_PROXIMITY);
    double losses = 0;
    if (previousLosses) {
        losses = previousLosses.value();
    }
    else {
        auto aux = _magneticFilterLossesNoProximity.evaluate_magnetic(magnetic, inputs, outputs);
        losses = aux.second;
    }
    auto previousTemperatureRise = get_scoring(magnetic->get_reference(), MagneticFilters::TEMPERATURE_RISE);
    double temperature = 0;
    if (previousTemperatureRise) {
        temperature = previousTemperatureRise.value();
    }
    else {
        auto aux = _magneticFilterTemperatureRise.evaluate_magnetic(magnetic, inputs, outputs);
        temperature = aux.second;
    }

    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring * temperature};
}

std::pair<bool, double> MagnetomotiveForce::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto coil = magnetic->get_coil();
    double maximumMagnetomotiveForce = 0;
    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        std::vector<double> currentRmsPerParallelPerWinding;
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = inputs->get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) {
                throw InvalidInputException("Current is missing in excitation");
            }
            if (!excitation.get_current()->get_processed()) {
                throw InvalidInputException("Current is not processed");
            }
            if (!excitation.get_current()->get_processed()->get_rms()) {
                throw InvalidInputException("Current RMS is not processed");
            }
            auto currentRms = excitation.get_current()->get_processed()->get_rms().value();
            currentRmsPerParallelPerWinding.push_back(currentRms / coil.get_functional_description()[windingIndex].get_number_parallels());
            if (!coil.get_layers_description()) {
                throw CoilNotProcessedException("Coil not wound");
            }
        }
        std::vector<double> magnetomotiveForcePerLayer;
        magnetomotiveForcePerLayer.push_back(0);
        auto layers = coil.get_layers_description().value();
        for (auto layer : layers) {
            double magnetomotiveForceThisLayer = magnetomotiveForcePerLayer.back();
            if (layer.get_type() == ElectricalType::CONDUCTION) {

                auto windingIndex = coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
                auto numberTurns = coil.get_functional_description()[windingIndex].get_number_turns();  
                auto numberPhysicalTurnsInLayer = 0;
                for (auto parallelProportion : layer.get_partial_windings()[0].get_parallels_proportion()) {
                    numberPhysicalTurnsInLayer += round(numberTurns * parallelProportion);
                }
                numberPhysicalTurnsInLayer *= layer.get_partial_windings().size();
                if (coil.get_functional_description()[windingIndex].get_isolation_side() == IsolationSide::PRIMARY) {
                    magnetomotiveForceThisLayer += numberPhysicalTurnsInLayer * currentRmsPerParallelPerWinding[windingIndex];
                }
                else {
                    magnetomotiveForceThisLayer -= numberPhysicalTurnsInLayer * currentRmsPerParallelPerWinding[windingIndex];
                }
            }
            else {
                magnetomotiveForcePerLayer.push_back(magnetomotiveForceThisLayer);
            }
        }

        double maximumMagnetomotiveForceThisOperatingPoint = *max_element(magnetomotiveForcePerLayer.begin(), magnetomotiveForcePerLayer.end());
        double minimumMagnetomotiveForceThisOperatingPoint = *min_element(magnetomotiveForcePerLayer.begin(), magnetomotiveForcePerLayer.end());
        maximumMagnetomotiveForceThisOperatingPoint = std::max(fabs(maximumMagnetomotiveForceThisOperatingPoint), fabs(minimumMagnetomotiveForceThisOperatingPoint));
        maximumMagnetomotiveForce = std::max(maximumMagnetomotiveForce, maximumMagnetomotiveForceThisOperatingPoint);

    }
    return {true, maximumMagnetomotiveForce};
}

} // namespace OpenMagnetics
