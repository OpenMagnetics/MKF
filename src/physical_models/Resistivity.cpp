#include "physical_models/Resistivity.h"
#include "spline.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"
#include "support/Utils.h"

namespace OpenMagnetics {

double ResistivityCoreMaterialModel::get_resistivity(ResistivityMaterial materialData, double temperature) {
    // MAS made core-material resistivity OPTIONAL (MAS b56fd13): a maker that publishes
    // permeability, saturation and a loss curve but no resistivity can now be represented
    // honestly instead of carrying an invented number. The model's contract is unchanged --
    // absent data THROWS, it does not fall back -- so a caller that needs a resistivity still
    // finds out loudly. Callers that can do without it (the cross-referencer's resistivity
    // dimension) test `get_resistivity()` for a value first, as the Curie filter already does.
    const auto& coreMaterial = std::get<CoreMaterial>(materialData);
    auto resistivityDataOptional = coreMaterial.get_resistivity();
    double resistivity = -1;

    if (!resistivityDataOptional) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "No resistivity data available for core material");
    }
    auto resistivityData = resistivityDataOptional.value();

    if (resistivityData.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "No resistivity data available for core material");
    }
    // ABT #2235: outside its tabulated temperatures the resistivity is the table extrapolated (held, for a
    // single point; the spline's end slope, for several) as before, with a WARNING whatever the flag. A point
    // without a temperature is a measurement at no stated temperature and is taken as it stands, as before.
    {
        double tableMinimum = std::numeric_limits<double>::max();
        double tableMaximum = std::numeric_limits<double>::lowest();
        bool allPointsHaveTemperature = true;
        for (const auto& datum : resistivityData) {
            if (!datum.get_temperature()) {
                allPointsHaveTemperature = false;
                break;
            }
            tableMinimum = std::min(tableMinimum, datum.get_temperature().value());
            tableMaximum = std::max(tableMaximum, datum.get_temperature().value());
        }
        if (allPointsHaveTemperature) {
            warn_if_material_temperature_outside_data(coreMaterial.get_name(), "resistivity", tableMinimum, tableMaximum, temperature);
        }
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
