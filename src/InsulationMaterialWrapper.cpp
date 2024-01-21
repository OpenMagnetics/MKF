#include "Utils.h"
#include "InsulationMaterialWrapper.h"
#include <cmath>
#include <filesystem>
#include <cfloat>
#include <numbers>
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

    std::pair<double, double> InsulationMaterialWrapper::get_thicker_tape(){
        double tapeThickness = 0;
        double tapeDielectricStrength = 0;

        for (auto dielectricStrengthElement : get_dielectric_strength()) {
            if (!dielectricStrengthElement.get_thickness()) {
                throw std::invalid_argument("Missing thickness in dielectric strength element");
            }
            if (dielectricStrengthElement.get_thickness().value() > tapeThickness) {
                tapeThickness = dielectricStrengthElement.get_thickness().value();
                tapeDielectricStrength = dielectricStrengthElement.get_value();
            }
        }
        return {tapeThickness, tapeDielectricStrength};
    }

    std::pair<double, double> InsulationMaterialWrapper::get_thinner_tape(){
        double tapeThickness = DBL_MAX;
        double tapeDielectricStrength = 0;

        for (auto dielectricStrengthElement : get_dielectric_strength()) {
            if (!dielectricStrengthElement.get_thickness()) {
                throw std::invalid_argument("Missing thickness in dielectric strength element");
            }
            if (dielectricStrengthElement.get_thickness().value() < tapeThickness) {
                tapeThickness = dielectricStrengthElement.get_thickness().value();
                tapeDielectricStrength = dielectricStrengthElement.get_value();
            }
        }
        return {tapeThickness, tapeDielectricStrength};
    }

    double InsulationMaterialWrapper::get_thicker_tape_thickness(){
        double tapeThickness = 0;

        for (auto dielectricStrengthElement : get_dielectric_strength()) {
            if (!dielectricStrengthElement.get_thickness()) {
                throw std::invalid_argument("Missing thickness in dielectric strength element");
            }
            if (dielectricStrengthElement.get_thickness().value() > tapeThickness) {
                tapeThickness = dielectricStrengthElement.get_thickness().value();
            }
        }
        return tapeThickness;
    }

    double InsulationMaterialWrapper::get_thinner_tape_thickness(){
        double tapeThickness = DBL_MAX;

        for (auto dielectricStrengthElement : get_dielectric_strength()) {
            if (!dielectricStrengthElement.get_thickness()) {
                throw std::invalid_argument("Missing thickness in dielectric strength element");
            }
            if (dielectricStrengthElement.get_thickness().value() < tapeThickness) {
                tapeThickness = dielectricStrengthElement.get_thickness().value();
            }
        }
        return tapeThickness;
    }

} // namespace OpenMagnetics
