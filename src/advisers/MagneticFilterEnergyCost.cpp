#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Wire.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <string>

namespace OpenMagnetics {

MagneticFilterEnergyStored::MagneticFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models) {
    _models = models;
    _magneticEnergy = MagneticEnergy(models);
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
}

MagneticFilterEnergyStored::MagneticFilterEnergyStored(Inputs inputs) {
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
}

std::pair<bool, double> MagneticFilterEnergyStored::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double totalStorableMagneticEnergy = 0;
    // Phase 6 (perf): cache operating-points by const-ref.
    const auto& operatingPoints = inputs->get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        const auto& operatingPoint = operatingPoints[operatingPointIndex];
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
    double temperature = inputs.get_maximum_temperature();
    // Phase 6 (perf): cache operating-points by const-ref.
    const auto& operatingPoints = inputs.get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        const auto& operatingPoint = operatingPoints[operatingPointIndex];
        auto excitation = Inputs::get_primary_excitation(operatingPoint);
        primaryCurrentRms = std::max(primaryCurrentRms, excitation.get_current().value().get_processed().value().get_rms().value());
        frequency = std::max(frequency, Inputs::get_switching_frequency(excitation));
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
    const auto& core = magnetic->get_core();

    double primaryNumberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
    double estimatedNeededWindingArea = primaryNumberTurns * _estimatedParallels * _estimatedWireTotalArea * (inputs->get_design_requirements().get_turns_ratios().size() + 1);
    WindingWindowElement windingWindow;

    std::string shapeName = core.get_shape_name();
    if (!is_pqi_or_ui_shape(shapeName)) {
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
            double marginAngle = wound_distance_to_angle(_averageMarginInWindingWindow, radialHeight / 2);
            if (std::isnan((marginAngle) / 360)) {
                throw NaNResultException("marginAngle: " + std::to_string(marginAngle));
            }
            // The margin blocks marginAngle degrees of the toroidal window; the
            // remaining (360 - marginAngle) degrees are usable.
            windingWindowArea *= (360 - marginAngle) / 360;
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
    double wireRelativeCosts = 0;
    for (auto wire : wires) {
        wireRelativeCosts += wire.get_relative_cost();
    }
    if (!magnetic->get_mutable_coil().get_layers_description()) {
        throw CoilNotProcessedException("Missing layers description to evaluate cost filter");
    }
    auto numberLayers = magnetic->get_mutable_coil().get_layers_description()->size();
    double scoring = numberLayers + wireRelativeCosts;
    return {true, scoring};
}

} // namespace OpenMagnetics
