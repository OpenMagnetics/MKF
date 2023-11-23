#pragma once

#include "InputsWrapper.h"
#include "WireWrapper.h"
#include "BobbinWrapper.h"
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
        bool are_sections_and_layers_fitting();
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
        bool try_wind(bool delimitAndCompact);

        std::vector<WindingStyle> wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, uint8_t numberSlots);
        WindingStyle wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, uint8_t numberSlots);
        bool wind_by_sections();
        bool wind_by_sections(std::vector<double> proportionPerWinding);
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

        void set_inputs(InputsWrapper inputs) {  // TODO: change to DesignRequirements?
            _inputs = inputs;
        }
        void set_interleaving_level(uint8_t interleavingLevel) {
            _interleavingLevel = interleavingLevel;
        }
        void set_winding_orientation(WindingOrientation windingOrientation) {
            _windingOrientation = windingOrientation;
        }
        void set_layers_orientation(WindingOrientation layersOrientation) {
            _layersOrientation = layersOrientation;
        }
        void set_turns_alignment(CoilAlignment turnsAlignment) {
            _turnsAlignment = turnsAlignment;
        }
        void set_section_alignment(CoilAlignment sectionAlignment) {
            _sectionAlignment = sectionAlignment;
        }
        uint8_t get_interleaving_level() {
            return _interleavingLevel;
        }
        WindingOrientation get_winding_orientation() {
            return _windingOrientation;
        }
        WindingOrientation get_layers_orientation() {
            return _layersOrientation;
        }
        CoilAlignment get_turns_alignment() {
            return _turnsAlignment;
        }
        CoilAlignment get_section_alignment() {
            return _sectionAlignment;
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
        void set_number_turns(std::vector<uint64_t> numberTurns);

        std::vector<Layer> get_layers_by_section(std::string sectionName);

        std::vector<Turn> get_turns_by_layer(std::string layerName);

        std::vector<uint64_t> get_number_parallels();
        void set_number_parallels(std::vector<uint64_t> numberParallels);

        CoilFunctionalDescription get_winding_by_name(std::string name);

        size_t get_winding_index_by_name(std::string name);

        std::vector<WireWrapper> get_wires();
        WireType get_wire_type(size_t windingIndex);
        static WireType get_wire_type(CoilFunctionalDescription coilFunctionalDescription);
        WireWrapper resolve_wire(size_t windingIndex);
        static WireWrapper resolve_wire(CoilFunctionalDescription coilFunctionalDescription);

        double horizontalFillingFactor(Section section);

        double verticalFillingFactor(Section section);

        double horizontalFillingFactor(Layer layer);

        double verticalFillingFactor(Layer layer);

        static BobbinWrapper resolve_bobbin(CoilWrapper coil);
        BobbinWrapper resolve_bobbin();

};
} // namespace OpenMagnetics
