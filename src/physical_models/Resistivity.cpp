#include "physical_models/Resistivity.h"
#include "spline.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

namespace OpenMagnetics {

double ResistivityCoreMaterialModel::get_resistivity(ResistivityMaterial materialData, double temperature) {
    auto resistivityData = std::get<CoreMaterial>(materialData).get_resistivity();
    double resistivity = -1;

    if (resistivityData.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "No resistivity data available for core material");
    }
    if (resistivityData.size() < 2) {
        resistivity = resistivityData[0].get_value();
    }
    else {
        size_t n = resistivityData.size();
        std::vector<double> x, y;

        for (size_t i = 0; i < n; i++) {
            if (x.empty() || fabs(resistivityData[i].get_temperature().value() - x.back()) > 1e-9) {
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
        return std::make_shared<ResistivityCoreMaterialModel>();
    }
    else if (modelName == ResistivityModels::WIRE_MATERIAL) {
        return std::make_shared<ResistivityWireMaterialModel>();
    }
    else
        throw ModelNotAvailableException("Unknown Resistivity model, available options are: {CORE_MATERIAL, WIRE_MATERIAL}");
}
} // namespace OpenMagnetics
