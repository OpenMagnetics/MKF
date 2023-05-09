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
using nlohmann::json;

class WindingWrapper : public Winding {
    public:
        uint8_t _interleavingLevel = 1;
        WindingOrientation _windingOrientation = WindingOrientation::HORIZONTAL;
        WindingOrientation _layersOrientation = WindingOrientation::HORIZONTAL;
        WindingWrapper(const json& j, uint8_t interleavingLevel = 1,
                       WindingOrientation windingOrientation = WindingOrientation::HORIZONTAL,
                       WindingOrientation layersOrientation = WindingOrientation::HORIZONTAL) {
            _interleavingLevel = interleavingLevel;
            _windingOrientation = windingOrientation;
            _layersOrientation = layersOrientation;
            from_json(j, *this);

            std::string bobbinName = "";
            if (std::holds_alternative<std::string>(get_bobbin())) {
                bobbinName = std::get<std::string>(get_bobbin());
                if (bobbinName != "Dummy") {
                    auto bobbinData = OpenMagnetics::find_bobbin_by_name(std::get<std::string>(get_bobbin()));
                    set_bobbin(bobbinData);
                }
            }

            if (bobbinName != "Dummy") {
                wind_by_sections();
                wind_by_layers();
                wind_by_turns();
            }
        }
        WindingWrapper() = default;
        virtual ~WindingWrapper() = default;

        std::vector<bool> wind_by_consecutive_turns(std::vector<WindingFunctionalDescription> windings);
        void wind_by_sections();
        void wind_by_layers();
        void wind_by_turns();

        uint64_t get_number_turns(size_t windingIndex) {
            return get_functional_description()[windingIndex].get_number_turns();
        }

        uint64_t get_number_parallels(size_t windingIndex) {
            return get_functional_description()[windingIndex].get_number_parallels();
        }

        std::string get_name(size_t windingIndex) {
            return get_functional_description()[windingIndex].get_name();
        }

        WindingFunctionalDescription get_winding_by_name(std::string name) {
            for (auto& windingFunctionalDescription : get_functional_description()) {
                if (windingFunctionalDescription.get_name() == name) {
                    return windingFunctionalDescription;
                }
            }
            throw std::runtime_error("No such a winding name: " + name);
        }

        size_t get_winding_index_by_name(std::string name) {
            for (size_t i=0; i<get_functional_description().size(); ++i) {
                if (get_functional_description()[i].get_name() == name) {
                    return i;
                }
            }
            throw std::runtime_error("No such a winding name: " + name);
        }

};
} // namespace OpenMagnetics
