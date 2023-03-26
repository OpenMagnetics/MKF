#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
// #include <WindingTemplate.hpp>
#include "Constants.h"
#include "Utils.h"
#include "WindingWrapper.h"
#include "json.hpp"

#include <magic_enum.hpp>

using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

namespace OpenMagnetics {

void WindingWrapper::wind_by_sections() {}

void WindingWrapper::wind_by_layers() {}

void WindingWrapper::wind_by_turns() {}

} // namespace OpenMagnetics
