#pragma once

#include "Utils.h"
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

class WireWrapper : public WireS {
  public:
    WireWrapper() = default;
    virtual ~WireWrapper() = default;

    static InsulationWireCoating get_coating(WireS wire);
    static double get_filling_factor(double conductingDiameter, int grade = 1, WireStandard standard = WireStandard::IEC_60317, bool includeAirInCell = false);
};
} // namespace OpenMagnetics
