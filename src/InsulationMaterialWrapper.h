#pragma once

#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using json = nlohmann::json;

namespace OpenMagnetics {


class InsulationMaterialWrapper : public InsulationMaterial {
  private:
    std::vector<std::pair<double, double>> _available_thicknesses;
  public:
    InsulationMaterialWrapper(const json& j) {
        from_json(j, *this);
        extract_available_thicknesses();
    }
    InsulationMaterialWrapper() = default;
    virtual ~InsulationMaterialWrapper() = default;
    void extract_available_thicknesses();
    std::vector<std::pair<double, double>> get_available_thicknesses();
    std::pair<double, double> get_thicker_tape();
    std::pair<double, double> get_thinner_tape();
    double get_thicker_tape_thickness();
    double get_thinner_tape_thickness();

};
} // namespace OpenMagnetics