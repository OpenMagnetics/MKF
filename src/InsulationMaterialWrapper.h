#pragma once

#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

namespace OpenMagnetics {
using nlohmann::json;


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

};
} // namespace OpenMagnetics