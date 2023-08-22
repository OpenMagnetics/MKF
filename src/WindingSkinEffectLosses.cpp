#include "WindingSkinEffectLosses.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

double WindingSkinEffectLossesModel::get_skin_depth(WireMaterialDataOrNameUnion material, double frequency, double temperature) {
    WireMaterial wireMaterial;
    auto constants = Constants();
    if (std::holds_alternative<WireMaterial>(material)) {
        wireMaterial = std::get<WireMaterial>(material);
    }
    else if (std::holds_alternative<std::string>(material)) {
        wireMaterial = OpenMagnetics::find_wire_material_by_name(std::get<std::string>(material));
    }
    auto resistivityModel = OpenMagnetics::ResistivityModel::factory(OpenMagnetics::ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wireMaterial, temperature);

    double skinDepth = sqrt(resistivity / (std::numbers::pi * frequency * constants.vacuumPermeability * wireMaterial.get_permeability()));
    return skinDepth;
};

// double get_effective_current_density(Harmonics harmonics, WireWrapper wire, double temperature) {
//     for (size_t i=0; i<harmonics.get_amplitudes(); ++i) {
//         auto skinDepth = get_skin_depth(material, frequency, temperature);
//     }
// }

} // namespace OpenMagnetics
