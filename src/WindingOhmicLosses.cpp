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

double WindingOhmicLosses::get_dc_resistance(Turn turn, WireWrapper wire, double temperature) {
    double wireConductingArea;
    double wireLength = turn.get_length();
    WireMaterial wireMaterial;
    WireS realWire;
    if (wire.get_type() == "litz") {
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

    if (wire.get_type() == "round" || wire.get_type() == "litz") {
        wireConductingArea = std::numbers::pi * pow(resolve_dimensional_values(realWire.get_conducting_diameter().value()) / 2, 2);
    }
    else if (wire.get_type() == "rectangular" || wire.get_type() == "foil") {
        wireConductingArea = resolve_dimensional_values(realWire.get_conducting_width().value()) * resolve_dimensional_values(realWire.get_conducting_height().value());
    }
    else {
        throw std::runtime_error("Unknown wire type in WindingOhmicLosses");
    }

    if (wire.get_number_conductors()) {
        double numberConductors = wire.get_number_conductors().value();
        wireConductingArea *= numberConductors;
    }


    double dcResistance = resistivity * wireLength / wireConductingArea;
    return dcResistance;
};

double WindingOhmicLosses::get_ohmic_losses(CoilWrapper winding, OperatingPoint operatingPoint, double temperature) {
    auto turns = winding.get_turns_description().value();
    std::vector<std::vector<double>> seriesResistancePerWindingPerParallel;
    std::vector<double> parallelResistancePerWinding;
    std::vector<std::vector<double>> dcCurrentPerWindingPerParallel;
    std::vector<double> dcCurrentPerWinding;
    auto wirePerWinding = winding.get_wires();
    for (size_t windingIndex = 0; windingIndex < winding.get_functional_description().size(); ++windingIndex) {
        seriesResistancePerWindingPerParallel.push_back(std::vector<double>(winding.get_number_parallels(windingIndex), 0));
        dcCurrentPerWindingPerParallel.push_back(std::vector<double>(winding.get_number_parallels(windingIndex), 0));
        dcCurrentPerWinding.push_back(operatingPoint.get_excitations_per_winding()[windingIndex].get_current().value().get_processed().value().get_rms().value());
    }

    for (auto& turn : turns) {
        auto windingIndex = winding.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        double turnResistance = get_dc_resistance(turn, wirePerWinding[windingIndex], temperature);
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
        parallelResistancePerWinding.push_back(parallelResistance);
    }

    std::vector<double> windingOhmicLossesPerWinding;
    double windingOhmicLossesTotal = 0;

    for (size_t windingIndex = 0; windingIndex < winding.get_functional_description().size(); ++windingIndex) {
        double windingOhmicLossesInWinding = 0;
        for (size_t parallelIndex = 0; parallelIndex < winding.get_number_parallels(windingIndex); ++parallelIndex) {
            windingOhmicLossesInWinding += seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] * pow(dcCurrentPerWindingPerParallel[windingIndex][parallelIndex], 2);
        }
        windingOhmicLossesPerWinding.push_back(windingOhmicLossesInWinding);
        windingOhmicLossesTotal += windingOhmicLossesInWinding;
    }

    return windingOhmicLossesTotal;
};
} // namespace OpenMagnetics
