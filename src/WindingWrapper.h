#pragma once

#include "InputsWrapper.h"
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
    private:
        std::map<std::pair<size_t, size_t>, Section> _insulationSections;
        std::map<std::pair<size_t, size_t>, std::vector<Layer>> _insulationLayers;
        std::map<std::pair<size_t, size_t>, std::string> _insulationSectionsLog;
        std::map<std::pair<size_t, size_t>, std::string> _insulationLayersLog;
        std::string windingLog;
        bool are_sections_and_layers_fitting() {
            double windTurns = true;
            if (!get_sections_description()) {
                return false;
            }
            if (!get_layers_description()) {
                return false;
            }
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
            return windTurns;
        }
    public:
        uint8_t _interleavingLevel = 1;
        WindingOrientation _windingOrientation = WindingOrientation::HORIZONTAL;
        WindingOrientation _layersOrientation = WindingOrientation::VERTICAL;
        CoilAlignment _turnsAlignment = CoilAlignment::CENTERED;
        CoilAlignment _sectionAlignment = CoilAlignment::CENTERED;
        std::optional<InputsWrapper> _inputs;

        WindingWrapper(const json& j, uint8_t interleavingLevel = 1,
                       WindingOrientation windingOrientation = WindingOrientation::HORIZONTAL,
                       WindingOrientation layersOrientation = WindingOrientation::VERTICAL,
                       CoilAlignment turnsAlignment = CoilAlignment::CENTERED,
                       CoilAlignment sectionAlignment = CoilAlignment::CENTERED) {
            _interleavingLevel = interleavingLevel;
            _windingOrientation = windingOrientation;
            _layersOrientation = layersOrientation;
            _turnsAlignment = turnsAlignment;
            _sectionAlignment = sectionAlignment;
            from_json(j, *this);

            wind();

        }

        WindingWrapper(const Winding winding) {
            bool hasSectionsData = false;
            bool hasLayersData = false;
            bool hasTurnsData = false;

            set_functional_description(winding.get_functional_description());
            set_bobbin(winding.get_bobbin());

            if (winding.get_sections_description()) {
                hasSectionsData = true;
                set_sections_description(winding.get_sections_description());
            }
            if (winding.get_layers_description()) {
                hasLayersData = true;
                set_layers_description(winding.get_layers_description());
            }
            if (winding.get_turns_description()) {
                hasTurnsData = true;
                set_turns_description(winding.get_turns_description());
            }

            if (!hasSectionsData || !hasLayersData || (!hasTurnsData && are_sections_and_layers_fitting())) {
                wind();
            }

        }
        WindingWrapper() = default;
        virtual ~WindingWrapper() = default;

        std::vector<bool> wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, uint8_t numberSlots);
        bool wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, uint8_t numberSlots);
        void wind_by_sections();
        void wind_by_layers();
        void wind_by_turns();
        void calculate_insulation();
        void delimit_and_compact();
        void log(std::string entry) {
            windingLog += entry + "\n";
        }
        std::string read_log() {
            return windingLog;
        }

        void set_inputs(InputsWrapper inputs) {
            _inputs = inputs;
        }

        void wind() {
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
                    if (_inputs) {
                        calculate_insulation();
                    }
                    wind_by_sections();
                    wind_by_layers();
                    if (are_sections_and_layers_fitting()) {
                        wind_by_turns();
                    }
                }
            }
        }

        uint64_t get_number_turns(size_t windingIndex) {
            return get_functional_description()[windingIndex].get_number_turns();
        }

        uint64_t get_number_parallels(size_t windingIndex) {
            return get_functional_description()[windingIndex].get_number_parallels();
        }

        uint64_t get_number_turns(Section section) {
            uint64_t physicalTurnsInSection = 0;
            auto partialWinding = section.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                physicalTurnsInSection += round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
            }
            return physicalTurnsInSection;
        }

        uint64_t get_number_turns(Layer layer) {
            uint64_t physicalTurnsInLayer = 0;
            auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                physicalTurnsInLayer += round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
            }
            return physicalTurnsInLayer;
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
