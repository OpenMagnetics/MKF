#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/Resistivity.h"
#include "Constants.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

namespace OpenMagnetics {

double WindingOhmicLosses::calculate_dc_resistance(Turn turn, const Wire& wire, double temperature) {
    double wireLength = turn.get_length();
    return calculate_dc_resistance(wireLength, wire, temperature);
}

double WindingOhmicLosses::calculate_dc_resistance(double wireLength, const Wire& wire, double temperature) {
    if (std::isnan(wireLength)) {
        throw NaNResultException("NaN found in wireLength value");
    }

    return calculate_dc_resistance_per_meter(wire, temperature) * wireLength;
}

double WindingOhmicLosses::calculate_dc_resistance_per_meter(Wire wire, double temperature) {

    WireMaterial wireMaterial = wire.resolve_material();

    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wireMaterial, temperature);

    double wireConductingArea = wire.calculate_conducting_area();

    double dcResistancePerMeter = resistivity / wireConductingArea;

    if (std::isnan(dcResistancePerMeter)) {
        throw NaNResultException("NaN found in dcResistancePerMeter value");
    }
    if (dcResistancePerMeter <= 0) {
        throw InvalidInputException(ErrorCode::CALCULATION_INVALID_RESULT, "dcResistancePerMeter must be positive");
    }
    return dcResistancePerMeter;
};

double WindingOhmicLosses::calculate_effective_resistance_per_meter(Wire wire, double frequency, double temperature) {
    WireMaterial wireMaterial = wire.resolve_material();

    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wireMaterial, temperature);

    double wireEffectiveConductingArea = wire.calculate_effective_conducting_area(frequency, temperature);

    double dcResistancePerMeter = resistivity / wireEffectiveConductingArea;
    return dcResistancePerMeter;
};

std::vector<double> WindingOhmicLosses::calculate_dc_resistance_per_winding(Coil coil, double temperature) {
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description");
    }
    auto turns = coil.get_turns_description().value();
    std::vector<std::vector<double>> seriesResistancePerWindingPerParallel;
    auto wirePerWinding = coil.get_wires();
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        seriesResistancePerWindingPerParallel.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0));
    }

    std::vector<double> dcResistancePerWinding;
    for (auto& turn : turns) {
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        double turnResistance = calculate_dc_resistance(turn, wirePerWinding[windingIndex], temperature);
        seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] += turnResistance;
    }

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        double conductance = 0;
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            conductance += 1. / seriesResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
        double parallelResistance = 1. / conductance;
        dcResistancePerWinding.push_back(parallelResistance);
    }

    return dcResistancePerWinding;
}

WindingLossesOutput WindingOhmicLosses::calculate_ohmic_losses(Coil coil, OperatingPoint operatingPoint, double temperature) {
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description");
    }
    auto turns = coil.get_turns_description().value();
    std::vector<std::vector<double>> seriesResistancePerWindingPerParallel;
    std::vector<std::vector<double>> dcCurrentPerWindingPerParallel;
    std::vector<double> dcCurrentPerWinding;
    auto wirePerWinding = coil.get_wires();
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        seriesResistancePerWindingPerParallel.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0));
        dcCurrentPerWindingPerParallel.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0));
        dcCurrentPerWinding.push_back(operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_rms().value());
    }

    std::vector<double> dcResistancePerTurn;
    std::vector<double> dcResistancePerWinding;
    for (auto& turn : turns) {
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        double turnResistance = calculate_dc_resistance(turn, wirePerWinding[windingIndex], temperature);
        dcResistancePerTurn.push_back(turnResistance);
        seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] += turnResistance;
    }

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        double conductance = 0;
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            conductance += 1. / seriesResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
        double parallelResistance = 1. / conductance;
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            dcCurrentPerWindingPerParallel[windingIndex][parallelIndex] = dcCurrentPerWinding[windingIndex] * parallelResistance / seriesResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
        dcResistancePerWinding.push_back(parallelResistance);
    }
    std::vector<WindingLossesPerElement> windingLossesPerTurn;
    std::vector<double> currentDividerPerTurn;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        Turn turn = turns[turnIndex];
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        auto currentDividerThisTurn = dcCurrentPerWinding[windingIndex] == 0? 0 : dcCurrentPerWindingPerParallel[windingIndex][parallelIndex] / dcCurrentPerWinding[windingIndex];

        double windingOhmicLossesInTurn = pow(dcCurrentPerWindingPerParallel[windingIndex][parallelIndex], 2) * dcResistancePerTurn[turnIndex];
        OhmicLosses ohmicLosses;
        WindingLossesPerElement windingLossesThisTurn;
        ohmicLosses.set_losses(windingOhmicLossesInTurn);
        ohmicLosses.set_method_used("Ohm");
        ohmicLosses.set_origin(ResultOrigin::SIMULATION);
        windingLossesThisTurn.set_ohmic_losses(ohmicLosses);
        windingLossesThisTurn.set_name(turn.get_name());
        windingLossesPerTurn.push_back(windingLossesThisTurn);

        if (std::isnan(currentDividerThisTurn)) {
            throw NaNResultException("NaN found in currentDividerThisTurn value");
        }
        currentDividerPerTurn.push_back(currentDividerThisTurn);
    }

    double windingOhmicLossesTotal = 0;

    std::vector<WindingLossesPerElement> windingLossesPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        double windingOhmicLossesInWinding = 0;

        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            windingOhmicLossesInWinding += seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] * pow(dcCurrentPerWindingPerParallel[windingIndex][parallelIndex], 2);
        }

        OhmicLosses ohmicLosses;
        WindingLossesPerElement windingLossesThisWinding;
        ohmicLosses.set_losses(windingOhmicLossesInWinding);
        ohmicLosses.set_method_used("Ohm");
        ohmicLosses.set_origin(ResultOrigin::SIMULATION);
        windingLossesThisWinding.set_ohmic_losses(ohmicLosses);
        windingLossesThisWinding.set_name(coil.get_functional_description()[windingIndex].get_name());
        windingLossesPerWinding.push_back(windingLossesThisWinding);
        windingOhmicLossesTotal += windingOhmicLossesInWinding;
    }

    WindingLossesOutput result;
    result.set_winding_losses_per_winding(windingLossesPerWinding);
    result.set_winding_losses_per_turn(windingLossesPerTurn);
    result.set_winding_losses(windingOhmicLossesTotal);
    result.set_temperature(temperature);
    result.set_temperature(temperature);
    result.set_origin(ResultOrigin::SIMULATION);
    result.set_dc_resistance_per_turn(dcResistancePerTurn);
    result.set_dc_resistance_per_winding(dcResistancePerWinding);
    result.set_current_per_winding(operatingPoint);
    result.set_current_divider_per_turn(currentDividerPerTurn);

    return result;
}

double WindingOhmicLosses::calculate_ohmic_losses_per_meter(Wire wire, SignalDescriptor current, double temperature) {

    double dcResistancePerMeter = calculate_dc_resistance_per_meter(wire, temperature);
    if (!current.get_processed()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is not processed");
    }
    if (!current.get_processed()->get_rms()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current processed is missing field RMS");
    }
    auto currentRms = current.get_processed()->get_rms().value();

    double windingOhmicLossesPerMeter = pow(currentRms, 2) * dcResistancePerMeter;

    return windingOhmicLossesPerMeter;
};
} // namespace OpenMagnetics
