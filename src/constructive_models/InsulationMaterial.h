#pragma once

#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS;

using json = nlohmann::json;

namespace OpenMagnetics {


class InsulationMaterial : public MAS::InsulationMaterial {
private:
    std::vector<std::pair<double, double>> _available_thicknesses;
public:
    InsulationMaterial(const json& j) {
        from_json(j, *this);
        extract_available_thicknesses();
    }
    InsulationMaterial() = default;
    virtual ~InsulationMaterial() = default;

    InsulationMaterial(MAS::InsulationMaterial insulationMaterial) {
        // The WHOLE record, never a field-by-field copy: that copy silently dropped every field it
        // did not list -- `form` (ABT #1174, a sleeve then looked like a tape) and
        // `minimumBendRadius` (ABT #1296, a sleeve's rated bend was ignored).
        MAS::InsulationMaterial::operator=(std::move(insulationMaterial));
    }

    void extract_available_thicknesses();
    std::vector<std::pair<double, double>> get_available_thicknesses();
    std::pair<double, double> get_thicker_tape();
    std::pair<double, double> get_thinner_tape();
    double get_thicker_tape_thickness();
    double get_thinner_tape_thickness();

    static double get_dielectric_strength_by_thickness(InsulationMaterial materialData, double thickness);
    double get_dielectric_strength_by_thickness(double thickness);

};
} // namespace OpenMagnetics