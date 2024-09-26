#pragma once

#include "Insulation.h"
#include "CoreWrapper.h"
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
        std::map<std::pair<size_t, size_t>, CoilSectionInterface> _coilSectionInterfaces;
        std::map<std::pair<size_t, size_t>, std::string> _insulationSectionsLog;
        std::map<std::pair<size_t, size_t>, std::string> _insulationLayersLog;
        std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> _sectionInfoWithInsulation;
        std::vector<std::vector<double>> _marginsPerSection;
        size_t _interleavingLevel = 1;
        WindingOrientation _windingOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation _layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment _turnsAlignment = CoilAlignment::CENTERED;
        CoilAlignment _sectionAlignment = CoilAlignment::INNER_OR_TOP;
        std::optional<InputsWrapper> _inputs;
        std::map<std::string, CoilAlignment> _turnsAlignmentPerSection;
        std::map<std::string, WindingOrientation> _layersOrientationPerSection;
        std::string coilLog;
        InsulationCoordinator _standardCoordinator = InsulationCoordinator();
        std::vector<double> _currentProportionPerWinding;
        std::vector<size_t> _currentPattern;
        size_t _currentRepetitions;
        bool _strict = true;
        bool _bobbin_resolved = false;
        BobbinWrapper _bobbin;

    public:

        CoilWrapper(const json& j, size_t interleavingLevel = 1,
                       WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING,
                       WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING,
                       CoilAlignment turnsAlignment = CoilAlignment::CENTERED,
                       CoilAlignment sectionAlignment = CoilAlignment::INNER_OR_TOP);
        CoilWrapper(const Coil coil);
        CoilWrapper() = default;
        virtual ~CoilWrapper() = default;
        bool try_wind();
        bool wind();
        bool wind(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions=1);
        bool wind(std::vector<size_t> pattern, size_t repetitions=1);
        bool wind(size_t repetitions);
        void try_rewind();
        void clear();
        bool are_sections_and_layers_fitting();

        std::vector<WindingStyle> wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, std::vector<size_t> numberSlots);
        WindingStyle wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, size_t numberSlots);
        std::vector<std::pair<size_t, double>> get_ordered_sections(double spaceForSections, std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions=1);
        std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> add_insulation_to_sections(std::vector<std::pair<size_t, double>> orderedSections);
        void remove_insulation_if_margin_is_enough(std::vector<std::pair<size_t, double>> orderedSections);
        void equalize_margins(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation);

        std::vector<double> get_proportion_per_winding_based_on_wires();
        void apply_margin_tape(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation);
        std::vector<double> get_aligned_section_dimensions_rectangular_window(size_t sectionIndex);
        std::vector<double> get_aligned_section_dimensions_round_window(size_t sectionIndex);
        size_t convert_conduction_section_index_to_global(size_t conductionSectionIndex);
        std::vector<double> cartesian_to_polar(std::vector<double> value);
        static std::vector<double> cartesian_to_polar(std::vector<double> value, double radialHeight);
        std::vector<double> polar_to_cartesian(std::vector<double> value);
        static std::vector<double> polar_to_cartesian(std::vector<double> value, double radialHeight);
        void convert_turns_to_cartesian_coordinates();
        void convert_turns_to_polar_coordinates();
        std::vector<std::pair<double, std::vector<double>>> get_collision_distances(std::vector<double> turnCoordinates, std::vector<std::vector<double>> placedTurnsCoordinates, double wireHeight);

        bool wind_by_sections();
        bool wind_by_sections(size_t repetitions);
        bool wind_by_sections(std::vector<double> proportionPerWinding);
        bool wind_by_sections(std::vector<size_t> pattern, size_t repetitions);
        bool wind_by_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions);
        bool wind_by_rectangular_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions);
        bool wind_by_round_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions);
        bool wind_by_layers();
        bool wind_by_rectangular_layers();
        bool wind_by_round_layers();
        bool wind_by_turns();
        bool wind_by_rectangular_turns();
        bool wind_by_round_turns();
        bool wind_toroidal_additional_turns();
        bool calculate_insulation(bool simpleMode = false);
        bool calculate_mechanical_insulation();
        bool delimit_and_compact();
        bool delimit_and_compact_rectangular_window();
        bool delimit_and_compact_round_window();

        void log(std::string entry);
        std::string read_log();
        void set_strict(bool value);
        void set_inputs(InputsWrapper inputs);  // TODO: change to DesignRequirements?

        void set_interleaving_level(uint8_t interleavingLevel);
        void reset_margins_per_section();
        size_t get_interleaving_level();
        void set_winding_orientation(WindingOrientation windingOrientation);
        void set_layers_orientation(WindingOrientation layersOrientation, std::optional<std::string> sectionName = std::nullopt);
        void set_turns_alignment(CoilAlignment turnsAlignment, std::optional<std::string> sectionName = std::nullopt);
        void set_section_alignment(CoilAlignment sectionAlignment);

        WindingOrientation get_winding_orientation();
        WindingOrientation get_layers_orientation();
        CoilAlignment get_turns_alignment(std::optional<std::string> sectionName = std::nullopt);
        CoilAlignment get_section_alignment();

        std::vector<Section> get_sections_description_conduction();
        std::vector<Layer> get_layers_description_conduction();
        std::vector<Section> get_sections_description_insulation();
        std::vector<Layer> get_layers_description_insulation();

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

        const std::vector<Section> get_sections_by_type(ElectricalType electricalType) const;
        const Section get_section_by_name(std::string name) const;
        const Turn get_turn_by_name(std::string name) const;
        const std::vector<Section> get_sections_by_winding(std::string windingName) const;

        std::vector<Layer> get_layers_by_section(std::string sectionName);
        const std::vector<Layer> get_layers_by_type(ElectricalType electricalType) const;
        std::vector<Layer> get_layers_by_winding_index(size_t windingIndex);

        std::vector<Turn> get_turns_by_layer(std::string layerName);
        std::vector<Turn> get_turns_by_section(std::string sectionName);

        std::vector<uint64_t> get_number_parallels();
        void set_number_parallels(std::vector<uint64_t> numberParallels);

        CoilFunctionalDescription get_winding_by_name(std::string name);

        size_t get_winding_index_by_name(std::string name);
        size_t get_turn_index_by_name(std::string name);
        size_t get_layer_index_by_name(std::string name);
        size_t get_section_index_by_name(std::string name);

        std::vector<WireWrapper> get_wires();
        WireType get_wire_type(size_t windingIndex);
        static WireType get_wire_type(CoilFunctionalDescription coilFunctionalDescription);
        WireWrapper resolve_wire(size_t windingIndex);
        static WireWrapper resolve_wire(CoilFunctionalDescription coilFunctionalDescription);

        double overlapping_filling_factor(Section section);

        double contiguous_filling_factor(Section section);

        static BobbinWrapper resolve_bobbin(CoilWrapper coil);
        BobbinWrapper resolve_bobbin();

        void preload_margins(std::vector<std::vector<double>> marginPairs);
        void add_margin_to_section_by_index(size_t sectionIndex, std::vector<double> margins);
        static double calculate_external_proportion_for_wires_in_toroidal_cores(CoreWrapper core, CoilWrapper coil);

};
} // namespace OpenMagnetics
