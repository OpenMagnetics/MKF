#pragma once

#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <streambuf>
#include <nlohmann/json-schema.hpp>
#include <MAS.hpp>
#include <magic_enum.hpp>
#include "json.hpp"
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;


namespace OpenMagnetics {
    using nlohmann::json;

    class WindingWrapper:public Winding {
        public:
        WindingWrapper(const json & j) {
            from_json(j, *this);
            wind_by_sections();
            wind_by_layers();
            wind_by_turns();
        }
        virtual ~WindingWrapper() = default;

        void wind_by_sections();
        void wind_by_layers();
        void wind_by_turns();
    };
}
