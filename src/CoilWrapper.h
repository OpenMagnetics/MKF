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
#include <numbers>
#include <streambuf>
#include <vector>
using json = nlohmann::json;

namespace OpenMagnetics {

class CoilWrapper : public Coil {
    private:
        std::map<std::pair<size_t, size_t>, Section> _insulationSections;
        std::map<std::pair<size_t, size_t>, std::vector<Layer>> _insulationLayers;
        std::map<std::pair<size_t, size_t>, std::string> _insulationSectionsLog;
        std::map<std::pair<size_t, size_t>, std::string> _insulationLayersLog;
        std::string coilLog;
        bool are_sections_and_layers_fitting() {
            bool windTurns = true;
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

        CoilWrapper(const json& j, uint8_t interleavingLevel = 1,
                       WindingOrientation windingOrientation = WindingOrientation::HORIZONTAL,
                       WindingOrientation layersOrientation = WindingOrientation::VERTICAL,
                       CoilAlignment turnsAlignment = CoilAlignment::CENTERED,
                       CoilAlignment sectionAlignment = CoilAlignment::CENTERED,
                       bool delimitAndCompact = true);
        CoilWrapper(const Coil coil, bool delimitAndCompact = true);
        CoilWrapper() = default;
        virtual ~CoilWrapper() = default;

        std::vector<WindingStyle> wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, uint8_t numberSlots);
        WindingStyle wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, uint8_t numberSlots);
        bool wind_by_sections();
        bool wind_by_layers();
        bool wind_by_turns();
        bool calculate_insulation();
        bool delimit_and_compact();
        void log(std::string entry) {
            coilLog += entry + "\n";
        }
        std::string read_log() {
            return coilLog;
        }

        void set_inputs(InputsWrapper inputs) {
            _inputs = inputs;
        }

        bool wind();

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

        CoilFunctionalDescription get_winding_by_name(std::string name);

        size_t get_winding_index_by_name(std::string name);

        std::vector<WireWrapper> get_wires();

        double horizontalFillingFactor(Section section);

        double verticalFillingFactor(Section section);

        double horizontalFillingFactor(Layer layer);

        double verticalFillingFactor(Layer layer);

};
} // namespace OpenMagnetics
