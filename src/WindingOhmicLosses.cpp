#include "WindingOhmicLosses.h"
#include "Resistivity.h"
#include "Constants.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

double WindingOhmicLosses::calculate_dc_resistance(Turn turn, WireWrapper wire, double temperature) {
    double wireConductingArea;
    double wireLength = turn.get_length();
    WireMaterial wireMaterial;
    Wire realWire;
    if (wire.get_type() == WireType::LITZ) {
        realWire = WireWrapper::get_strand(wire);
    }
    else {
        realWire = wire;
    }

    auto wireMaterialDataOrString = realWire.get_material().value();
    if (std::holds_alternative<std::string>(wireMaterialDataOrString)) {
        auto wireMaterialName = std::get<std::string>(wireMaterialDataOrString);
        wireMaterial = find_wire_material_by_name(wireMaterialName);
    }
    else {
        wireMaterial = std::get<WireMaterial>(wireMaterialDataOrString);
    }
    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wireMaterial, temperature);

    switch(wire.get_type()) {
        case WireType::ROUND:
            wireConductingArea = std::numbers::pi * pow(resolve_dimensional_values(realWire.get_conducting_diameter().value()) / 2, 2);
            break;
        case WireType::LITZ:
            wireConductingArea = std::numbers::pi * pow(resolve_dimensional_values(realWire.get_conducting_diameter().value()) / 2, 2);
            break;
        case WireType::RECTANGULAR:
        case WireType::FOIL:
            wireConductingArea = resolve_dimensional_values(realWire.get_conducting_width().value()) * resolve_dimensional_values(realWire.get_conducting_height().value());
            break;
        default:
            throw std::runtime_error("Unknown wire type in WindingOhmicLosses");
    }

    if (wire.get_number_conductors()) {
        double numberConductors = wire.get_number_conductors().value();
        wireConductingArea *= numberConductors;
    }
    double dcResistance = resistivity * wireLength / wireConductingArea;
    return dcResistance;
};

WindingLossesOutput WindingOhmicLosses::calculate_ohmic_losses(CoilWrapper winding, OperatingPoint operatingPoint, double temperature) {
    if (!winding.get_turns_description()) {
        throw std::runtime_error("Missing turns description");
    }
    auto turns = winding.get_turns_description().value();
    std::vector<std::vector<double>> seriesResistancePerWindingPerParallel;
    std::vector<std::vector<double>> dcCurrentPerWindingPerParallel;
    std::vector<double> dcCurrentPerWinding;
    auto wirePerWinding = winding.get_wires();
    for (size_t windingIndex = 0; windingIndex < winding.get_functional_description().size(); ++windingIndex) {
        seriesResistancePerWindingPerParallel.push_back(std::vector<double>(winding.get_number_parallels(windingIndex), 0));
        dcCurrentPerWindingPerParallel.push_back(std::vector<double>(winding.get_number_parallels(windingIndex), 0));
        dcCurrentPerWinding.push_back(operatingPoint.get_excitations_per_winding()[windingIndex].get_current().value().get_processed().value().get_rms().value());
    }

    std::vector<double> dcResistancePerTurn;
    std::vector<double> dcResistancePerWinding;
    for (auto& turn : turns) {
        auto windingIndex = winding.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        double turnResistance = calculate_dc_resistance(turn, wirePerWinding[windingIndex], temperature);
        dcResistancePerTurn.push_back(turnResistance);
        seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] += turnResistance;
    }

    for (size_t windingIndex = 0; windingIndex < winding.get_functional_description().size(); ++windingIndex) {
        double conductance = 0;
        for (size_t parallelIndex = 0; parallelIndex < winding.get_number_parallels(windingIndex); ++parallelIndex) {
            conductance += 1. / seriesResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
        double parallelResistance = 1. / conductance;
        for (size_t parallelIndex = 0; parallelIndex < winding.get_number_parallels(windingIndex); ++parallelIndex) {
            dcCurrentPerWindingPerParallel[windingIndex][parallelIndex] = dcCurrentPerWinding[windingIndex] * parallelResistance / seriesResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
        dcResistancePerWinding.push_back(parallelResistance);
    }
    std::vector<WindingLossesPerElement> windingLossesPerTurn;
    std::vector<double> currentDividerPerTurn;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        Turn turn = turns[turnIndex];
        auto windingIndex = winding.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        auto currentDividerThisTurn = dcCurrentPerWindingPerParallel[windingIndex][parallelIndex] / dcCurrentPerWinding[windingIndex];
        double windingOhmicLossesInTurn = pow(currentDividerThisTurn, 2) * dcResistancePerTurn[turnIndex];
        OhmicLosses ohmicLosses;
        WindingLossesPerElement windingLossesThisTurn;
        ohmicLosses.set_losses(windingOhmicLossesInTurn);
        ohmicLosses.set_method_used("Ohm");
        ohmicLosses.set_origin(ResultOrigin::SIMULATION);
        windingLossesThisTurn.set_ohmic_losses(ohmicLosses);
        windingLossesPerTurn.push_back(windingLossesThisTurn);
        currentDividerPerTurn.push_back(currentDividerThisTurn);
    }

    double windingOhmicLossesTotal = 0;

    std::vector<WindingLossesPerElement> windingLossesPerWinding;
    for (size_t windingIndex = 0; windingIndex < winding.get_functional_description().size(); ++windingIndex) {
        double windingOhmicLossesInWinding = 0;

        for (size_t parallelIndex = 0; parallelIndex < winding.get_number_parallels(windingIndex); ++parallelIndex) {
            windingOhmicLossesInWinding += seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] * pow(dcCurrentPerWindingPerParallel[windingIndex][parallelIndex], 2);
        }
        OhmicLosses ohmicLosses;
        WindingLossesPerElement windingLossesThisWinding;
        ohmicLosses.set_losses(windingOhmicLossesInWinding);
        ohmicLosses.set_method_used("Ohm");
        ohmicLosses.set_origin(ResultOrigin::SIMULATION);
        windingLossesThisWinding.set_ohmic_losses(ohmicLosses);
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
};
} // namespace OpenMagnetics
