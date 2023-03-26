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

class WindingWrapper : public Winding {
  public:
    WindingWrapper(const json& j) {
        from_json(j, *this);
        wind_by_sections();
        wind_by_layers();
        wind_by_turns();
    }
    WindingWrapper() = default;
    virtual ~WindingWrapper() = default;

    void wind_by_sections();
    void wind_by_layers();
    void wind_by_turns();
};
} // namespace OpenMagnetics
