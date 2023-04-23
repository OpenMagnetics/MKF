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

OpenMagnetics::CoreMaterial find_core_material_by_name(std::string name);
OpenMagnetics::CoreShape find_core_shape_by_name(std::string name);

void load_databases(bool withAliases=true);

std::vector<std::string> get_material_names();

std::vector<std::string> get_shape_names();

template<int decimals> double roundFloat(double value);

bool is_size_power_of_2(std::vector<double> data);

double roundFloat(double value, size_t decimals);

CoreShape flatten_dimensions(CoreShape shape);

double tryGetDutyCycle(Waveform waveform, double frequency);
} // namespace OpenMagnetics