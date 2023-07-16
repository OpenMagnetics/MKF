#include "Utils.h"
#include "InsulationMaterialWrapper.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

    void InsulationMaterialWrapper::extract_available_thicknesses(){

        for (auto& dielectricStrengthPoint : get_dielectric_strength()) {
            if (!dielectricStrengthPoint.get_thickness()) {
                continue;
            }
            std::pair<double, double> aux = std::pair<double, double>({dielectricStrengthPoint.get_thickness().value(), dielectricStrengthPoint.get_thickness().value() * dielectricStrengthPoint.get_value()});
            if (std::find(_available_thicknesses.begin(), _available_thicknesses.end(), aux) == _available_thicknesses.end()) {
                _available_thicknesses.push_back(aux);
            }
        }
    }

    std::vector<std::pair<double, double>> InsulationMaterialWrapper::get_available_thicknesses(){
        return _available_thicknesses;
    }

} // namespace OpenMagnetics
