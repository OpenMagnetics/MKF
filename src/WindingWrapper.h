#pragma once

#include "WireWrapper.h"
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
        OrientationEnum _windingOrientation = OrientationEnum::HORIZONTAL;
        OrientationEnum _layersOrientation = OrientationEnum::VERTICAL;
        TurnsAlignmentEnum _turnsAlignment = TurnsAlignmentEnum::CENTERED;
        WindingWrapper(const json& j, uint8_t interleavingLevel = 1,
                       OrientationEnum orientationEnum = OrientationEnum::HORIZONTAL,
                       OrientationEnum layersOrientation = OrientationEnum::VERTICAL,
                       TurnsAlignmentEnum turnsAlignment = TurnsAlignmentEnum::CENTERED) {
            _interleavingLevel = interleavingLevel;
            _windingOrientation = orientationEnum;
            _layersOrientation = layersOrientation;
            _turnsAlignment = turnsAlignment;
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
                bool wind = true;                
                for (auto& winding : get_mutable_functional_description()) {
                    if (std::holds_alternative<std::string>(winding.get_wire())) {
                        std::string wireName = std::get<std::string>(winding.get_wire());
                        if (wireName == "Dummy") {
                            wind = false;
                            break;
                        }
                        auto wire = find_wire_by_name(wireName);
                        winding.set_wire(wire);
                    }
                }

                if (wind) {
                    double windTurns = true;
                    wind_by_sections();
                    wind_by_layers();
                    auto sections = get_sections_description().value();
                    auto layers = get_layers_description().value();
                    for (auto& section: sections) {
                        if (section.get_filling_factor().value() > 1 || horizontalFillingFactor(section) > 1 || verticalFillingFactor(section) > 1) {
                            windTurns = false;
                        }
                    }
                    for (auto& layer: layers) {
                        if (layer.get_filling_factor().value() > 1) {
                            windTurns = false;
                        }
                    }
                    if (windTurns) {
                        wind_by_turns();
                    }
                }
            }
        }
        WindingWrapper() = default;
        virtual ~WindingWrapper() = default;

        std::vector<bool> wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, uint8_t numberSlots);
        bool wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, uint8_t numberSlots);
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

        std::vector<uint64_t> get_number_turns();

        std::vector<Layer> get_layers_by_section(std::string sectionName);

        std::vector<Turn> get_turns_by_layer(std::string layerName);

        std::vector<uint64_t> get_number_parallels();

        WindingFunctionalDescription get_winding_by_name(std::string name);

        size_t get_winding_index_by_name(std::string name);

        std::vector<WireWrapper> get_wires();

        double horizontalFillingFactor(Section section);

        double verticalFillingFactor(Section section);

        double horizontalFillingFactor(Layer layer);

        double verticalFillingFactor(Layer layer);

};
} // namespace OpenMagnetics
