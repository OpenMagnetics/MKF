#pragma once

#include <nlohmann/json-schema.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <MAS.hpp>
#include "json.hpp"
#include "Constants.h"
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

extern std::map<std::string, json> coreMaterialDatabase;
extern std::map<std::string, json> coreShapeDatabase;
extern OpenMagnetics::Constants constants;


namespace OpenMagnetics {

    template <typename T>
    T find_data_by_name(std::string name);
    template <int decimals>
    double roundFloat(double value);

}