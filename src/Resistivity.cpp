#include "Resistivity.h"
#include "spline.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

double ResistivityCoreMaterialModel::get_resistivity(ResistivityMaterial materialData, double temperature) {
    auto resistivityData = std::get<CoreMaterial>(materialData).get_resistivity();
    double resistivity = -1;

    if (resistivityData.size() < 2) {
        resistivity = resistivityData[0].get_value();
    }
    else {
        size_t n = resistivityData.size();
        std::vector<double> x, y;

        for (size_t i = 0; i < n; i++) {
            if (x.size() == 0 || resistivityData[i].get_temperature().value() != x.back()) {
                x.push_back(resistivityData[i].get_temperature().value());
                y.push_back(resistivityData[i].get_value());
            }
        }

        tk::spline interp(x, y, tk::spline::cspline_hermite, true);
        resistivity = interp(temperature);
    }
    return resistivity;
};

double ResistivityWireMaterialModel::get_resistivity(ResistivityMaterial materialData, double temperature) {
    auto resistivityData = std::get<WireMaterial>(materialData).get_resistivity();
    return resistivityData.get_reference_value() * (1 + resistivityData.get_temperature_coefficient() * (temperature - resistivityData.get_reference_temperature()));
};

std::shared_ptr<ResistivityModel> ResistivityModel::factory(ResistivityModels modelName) {
    if (modelName == ResistivityModels::CORE_MATERIAL) {
        std::shared_ptr<ResistivityModel> ResistivityModel(new ResistivityCoreMaterialModel);
        return ResistivityModel;
    }
    else if (modelName == ResistivityModels::WIRE_MATERIAL) {
        std::shared_ptr<ResistivityModel> ResistivityModel(new ResistivityWireMaterialModel);
        return ResistivityModel;
    }
    else
        throw std::runtime_error("Unknown Resistivity model, available options are: {CORE_MATERIAL, WIRE_MATERIAL}");
}
} // namespace OpenMagnetics
