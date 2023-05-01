#include "Resistivity.h"
#include <libInterpolate/Interpolate.hpp>

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
        _1D::CubicSplineInterpolator<double> interp;
        _1D::CubicSplineInterpolator<double>::VectorType xx(n), yy(n);

        for (size_t i = 0; i < n; i++) {
            xx(i) = resistivityData[i].get_temperature().value();
            yy(i) = resistivityData[i].get_value();
        }

        interp.setData(xx, yy);
        resistivity = interp(temperature);
    }
    return resistivity;
};

double ResistivityWireMaterialModel::get_resistivity(ResistivityMaterial materialData, double temperature) {
    return 42;
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
