#include "support/Utils.h"
#include "constructive_models/InsulationMaterial.h"
#include "support/Exceptions.h"
#include <cmath>
#include <filesystem>
#include <cfloat>
#include <numbers>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

    void InsulationMaterial::extract_available_thicknesses(){

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

    std::vector<std::pair<double, double>> InsulationMaterial::get_available_thicknesses(){
        return _available_thicknesses;
    }

    std::pair<double, double> InsulationMaterial::get_thicker_tape(){
        double tapeThickness = 0;
        double tapeDielectricStrength = 0;

        for (auto dielectricStrengthElement : get_dielectric_strength()) {
            if (!dielectricStrengthElement.get_thickness()) {
                continue;
                // throw std::invalid_argument("Missing thickness in dielectric strength element");
            }
            if (dielectricStrengthElement.get_thickness().value() > tapeThickness) {
                tapeThickness = dielectricStrengthElement.get_thickness().value();
                tapeDielectricStrength = dielectricStrengthElement.get_value();
            }
        }
        return {tapeThickness, tapeDielectricStrength};
    }

    std::pair<double, double> InsulationMaterial::get_thinner_tape(){
        double tapeThickness = DBL_MAX;
        double tapeDielectricStrength = 0;

        for (auto dielectricStrengthElement : get_dielectric_strength()) {
            if (!dielectricStrengthElement.get_thickness()) {
                continue;
                // throw std::invalid_argument("Missing thickness in dielectric strength element");
            }
            if (dielectricStrengthElement.get_thickness().value() < tapeThickness) {
                tapeThickness = dielectricStrengthElement.get_thickness().value();
                tapeDielectricStrength = dielectricStrengthElement.get_value();
            }
        }
        if (tapeThickness == DBL_MAX) {
            return {0, 0}; // FIX M-IMC-1: Return {0,0} when no thicknesses available
        }
        return {tapeThickness, tapeDielectricStrength};
    }

    double InsulationMaterial::get_thicker_tape_thickness(){
        double tapeThickness = 0;

        for (auto dielectricStrengthElement : get_dielectric_strength()) {
            if (!dielectricStrengthElement.get_thickness()) {
                continue;
                // throw std::invalid_argument("Missing thickness in dielectric strength element");
            }
            if (dielectricStrengthElement.get_thickness().value() > tapeThickness) {
                tapeThickness = dielectricStrengthElement.get_thickness().value();
            }
        }
        return tapeThickness;
    }

    double InsulationMaterial::get_thinner_tape_thickness(){
        double tapeThickness = DBL_MAX;

        for (auto dielectricStrengthElement : get_dielectric_strength()) {
            if (!dielectricStrengthElement.get_thickness()) {
                continue;
                // throw std::invalid_argument("Missing thickness in dielectric strength element");
            }
            if (dielectricStrengthElement.get_thickness().value() < tapeThickness) {
                tapeThickness = dielectricStrengthElement.get_thickness().value();
            }
        }
        return tapeThickness;
    }

    double InsulationMaterial::get_dielectric_strength_by_thickness(double thickness) {
        return get_dielectric_strength_by_thickness(*this, thickness);
    }
    double InsulationMaterial::get_dielectric_strength_by_thickness(InsulationMaterial materialData, double thickness) {
        auto relativePermittivityData = materialData.get_dielectric_strength();
        double relativePermittivity = -1;

        // Both reads below were unchecked (ABT #1616): [0] of an empty list, and .value() of a point
        // with no thickness when the strength has to be interpolated over thickness.
        if (relativePermittivityData.empty()) {
            throw InvalidInputException(ErrorCode::INVALID_INSULATION_DATA,
                                        "Insulation material " + materialData.get_name() + " has no dielectric strength data");
        }
        if (relativePermittivityData.size() < 2) {
            relativePermittivity = relativePermittivityData[0].get_value();
        }
        else {
            size_t n = relativePermittivityData.size();
            std::vector<double> x, y;

            for (size_t i = 0; i < n; i++) {
                if (!relativePermittivityData[i].get_thickness()) {
                    throw InvalidInputException(ErrorCode::INVALID_INSULATION_DATA,
                                                "Insulation material " + materialData.get_name() + " lists several dielectric strengths, and point " +
                                                std::to_string(i) + " has no thickness to interpolate them over");
                }
                if (x.size() == 0 || relativePermittivityData[i].get_thickness().value() != x.back()) {
                    x.push_back(relativePermittivityData[i].get_thickness().value());
                    y.push_back(relativePermittivityData[i].get_value());
                }
            }

            tk::spline interp(x, y, tk::spline::cspline_hermite, true);
            relativePermittivity = interp(thickness);
        }
        return relativePermittivity;
    };

} // namespace OpenMagnetics
