#pragma once

#include "Constants.h"
#include "json.hpp"

#include <MAS.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json-schema.hpp>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

extern std::map<std::string, json> coreMaterialDatabase;
extern std::map<std::string, json> coreShapeDatabase;
extern OpenMagnetics::Constants constants;

namespace OpenMagnetics {

template<typename T> T find_data_by_name(std::string name);

template<int decimals> double roundFloat(double value);

double roundFloat(double value, size_t decimals);

CoreShape flatten_dimensions(CoreShape shape);

double tryGetDutyCycle(Waveform waveform, double frequency);
} // namespace OpenMagnetics