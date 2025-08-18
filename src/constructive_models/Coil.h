#pragma once

#include "constructive_models/Insulation.h"
#include "constructive_models/Core.h"
#include "processors/Inputs.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Bobbin.h"
#include "json.hpp"

#include <MAS.hpp>
#include <vector>

using namespace MAS;
using json = nlohmann::json;

namespace OpenMagnetics {

using BobbinDataOrNameUnion = std::variant<Bobbin, std::string>;
using WireDataOrNameUnion = std::variant<Wire, std::string>;


class CoilFunctionalDescription : public MAS::CoilFunctionalDescription {
    private:
        WireDataOrNameUnion wire;

    public:
        const WireDataOrNameUnion & get_wire() const { return wire; }
        WireDataOrNameUnion & get_mutable_wire() { return wire; }
        void set_wire(const WireDataOrNameUnion & value) { this->wire = value; }
        CoilFunctionalDescription(const MAS::CoilFunctionalDescription coilFunctionalDescription) {

            if (coilFunctionalDescription.get_connections()) {
                set_connections(coilFunctionalDescription.get_connections());
            }

            set_isolation_side(coilFunctionalDescription.get_isolation_side());
            set_name(coilFunctionalDescription.get_name());
            set_number_parallels(coilFunctionalDescription.get_number_parallels());
            set_number_turns(coilFunctionalDescription.get_number_turns());
            auto wireVariant = coilFunctionalDescription.get_wire();
            if (std::holds_alternative<std::string>(wireVariant)) {
                set_wire(std::get<std::string>(wireVariant));
            }
            else {
                set_wire(std::get<MAS::Wire>(wireVariant));
            }
        };
        CoilFunctionalDescription() = default;
        virtual ~CoilFunctionalDescription() = default;

        Wire resolve_wire();
};

class Coil : public MAS::Coil {
    private:
        std::map<std::pair<size_t, size_t>, Section> _insulationSections;
        std::map<std::pair<size_t, size_t>, std::vector<Layer>> _insulationInterSectionsLayers;
        std::map<size_t, Layer> _insulationInterLayers;
        std::map<std::pair<size_t, size_t>, CoilSectionInterface> _coilSectionInterfaces;
        std::map<std::pair<size_t, size_t>, std::string> _insulationSectionsLog;
        std::map<std::pair<size_t, size_t>, std::string> _insulationInterSectionsLayersLog;
        std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> _sectionInfoWithInsulation;
        std::vector<std::vector<double>> _marginsPerSection;
        size_t _interleavingLevel = 1;
        WindingOrientation _windingOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation _layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment _turnsAlignment = CoilAlignment::CENTERED;
        CoilAlignment _sectionAlignment = CoilAlignment::INNER_OR_TOP;
        std::optional<Inputs> _inputs;
        std::map<std::string, CoilAlignment> _turnsAlignmentPerSection;
        std::map<std::string, WindingOrientation> _layersOrientationPerSection;
        std::string coilLog;
        InsulationCoordinator _standardCoordinator = InsulationCoordinator();
        std::vector<double> _currentProportionPerWinding;
        std::vector<size_t> _currentPattern;
        size_t _currentRepetitions;
        bool _strict = true;
        bool _bobbin_resolved = false;
        Bobbin _bobbin;
        BobbinDataOrNameUnion bobbin;
        std::vector<CoilFunctionalDescription> functional_description;

        bool wind_by_rectangular_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions);
        bool wind_by_round_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions);
        bool wind_by_rectangular_layers();
        bool wind_by_round_layers();
        bool wind_by_rectangular_turns();
        bool wind_by_round_turns();
        bool wind_toroidal_additional_turns();
        bool delimit_and_compact_rectangular_window();
        bool delimit_and_compact_round_window();
        bool create_default_group(Bobbin bobbin, WiringTechnology coilType = WiringTechnology::WOUND, double coreToLayerDistance = 0);

    public:
        bool wind_by_planar_sections(std::vector<size_t> stackUp, std::map<std::pair<size_t, size_t>, double> insulationThickness = {}, double coreToLayerDistance = 0);
        bool wind_by_planar_layers();
        bool wind_by_planar_turns(double borderToWireDistance, std::map<size_t, double> wireToWireDistance);

        Coil(const json& j, size_t interleavingLevel = 1,
                       WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING,
                       WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING,
                       CoilAlignment turnsAlignment = CoilAlignment::CENTERED,
                       CoilAlignment sectionAlignment = CoilAlignment::INNER_OR_TOP);
        Coil(const MAS::Coil coil);
        Coil(const json& j, bool windInConstructor);
        Coil() = default;
        virtual ~Coil() = default;
        bool fast_wind();
        bool unwind();
        bool wind();
        bool wind(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions=1);
        bool wind(std::vector<size_t> pattern, size_t repetitions=1);
        bool wind(size_t repetitions);
        bool wind_planar(std::vector<size_t> stackUp, std::optional<double> borderToWireDistance = std::nullopt, std::map<std::pair<size_t, size_t>, double> insulationThickness = {}, double coreToLayerDistance = 0);
        void try_rewind();
        void clear();
        bool are_sections_and_layers_fitting();

        const BobbinDataOrNameUnion & get_bobbin() const { return bobbin; }
        BobbinDataOrNameUnion & get_mutable_bobbin() { return bobbin; }
        void set_bobbin(const BobbinDataOrNameUnion & value) { this->bobbin = value; }

        const std::vector<CoilFunctionalDescription> & get_functional_description() const { return functional_description; }
        std::vector<CoilFunctionalDescription> & get_mutable_functional_description() { return functional_description; }
        void set_functional_description(const std::vector<CoilFunctionalDescription> & value) { this->functional_description = value; }


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
        bool wind_by_layers();
        bool wind_by_turns();
        bool calculate_insulation(bool simpleMode = false);
        bool calculate_custom_thickness_insulation(double thickness);
        bool calculate_mechanical_insulation();
        bool delimit_and_compact();

        void log(std::string entry);
        std::string read_log();
        void set_strict(bool value);
        void set_inputs(Inputs inputs);  // TODO: change to DesignRequirements?

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

        std::string get_name(size_t windingIndex);

        WiringTechnology get_coil_type(size_t groupIndex = 0);

        std::vector<uint64_t> get_number_turns();
        uint64_t get_number_turns(size_t windingIndex);
        uint64_t get_number_turns(Section section);
        uint64_t get_number_turns(Layer layer);
        void set_number_turns(std::vector<uint64_t> numberTurns);
        std::vector<IsolationSide> get_isolation_sides();
        void set_isolation_sides(std::vector<IsolationSide> isolationSides);
        std::vector<uint64_t> get_number_parallels();
        uint64_t get_number_parallels(size_t windingIndex);
        void set_number_parallels(std::vector<uint64_t> numberParallels);

        void set_interlayer_insulation(double layerThickness, std::optional<std::string> material = std::nullopt, std::optional<std::string> windingName = std::nullopt, bool autowind=true);
        void set_intersection_insulation(double layerThickness, size_t numberInsulationLayers, std::optional<std::string> material = std::nullopt, std::optional<std::pair<std::string, std::string>> windingNames = std::nullopt, bool autowind=true);

        std::vector<Section> get_sections_by_group(std::string groupName);
        const std::vector<Section> get_sections_by_type(ElectricalType electricalType) const;
        const Section get_section_by_name(std::string name) const;
        const Turn get_turn_by_name(std::string name) const;
        const std::vector<Section> get_sections_by_winding(std::string windingName) const;

        std::vector<Layer> get_layers_by_section(std::string sectionName);
        const std::vector<Layer> get_layers_by_type(ElectricalType electricalType) const;
        std::vector<Layer> get_layers_by_winding_index(size_t windingIndex);
        const Layer get_layer_by_name(std::string name) const;

        std::vector<Turn> get_turns_by_layer(std::string layerName);
        std::vector<Turn> get_turns_by_section(std::string sectionName);
        std::vector<Turn> get_turns_by_winding(std::string windingName);

        std::vector<std::string> get_layers_names_by_winding(std::string windingName);
        std::vector<std::string> get_layers_names_by_section(std::string sectionName);
        std::vector<std::string> get_turns_names_by_layer(std::string layerName);
        std::vector<std::string> get_turns_names_by_section(std::string sectionName);
        std::vector<std::string> get_turns_names_by_winding(std::string windingName);

        std::vector<size_t> get_turns_indexes_by_layer(std::string layerName);
        std::vector<size_t> get_turns_indexes_by_section(std::string sectionName);
        std::vector<size_t> get_turns_indexes_by_winding(std::string windingName);

        CoilFunctionalDescription get_winding_by_name(std::string name);

        size_t get_winding_index_by_name(std::string name);
        size_t get_turn_index_by_name(std::string name);
        size_t get_layer_index_by_name(std::string name);
        size_t get_section_index_by_name(std::string name);

        std::vector<Wire> get_wires();
        WireType get_wire_type(size_t windingIndex);
        static WireType get_wire_type(CoilFunctionalDescription coilFunctionalDescription);
        std::string get_wire_name(size_t windingIndex);
        static std::string get_wire_name(CoilFunctionalDescription coilFunctionalDescription);
        Wire resolve_wire(size_t windingIndex);
        static Wire resolve_wire(CoilFunctionalDescription coilFunctionalDescription);

        double overlapping_filling_factor(Section section);

        double contiguous_filling_factor(Section section);
        std::pair<double, std::pair<double, double>> calculate_filling_factor(size_t groupIndex = 0);

        static Bobbin resolve_bobbin(Coil coil);
        Bobbin resolve_bobbin();

        void preload_margins(std::vector<std::vector<double>> marginPairs);
        void add_margin_to_section_by_index(size_t sectionIndex, std::vector<double> margins);
        static double calculate_external_proportion_for_wires_in_toroidal_cores(Core core, Coil coil);

        void set_insulation_layers(std::map<std::pair<size_t, size_t>, std::vector<Layer>> insulationLayers);

        static InsulationMaterial resolve_insulation_layer_insulation_material(Coil coil, std::string layerName);
        InsulationMaterial resolve_insulation_layer_insulation_material(std::string layerName);
        InsulationMaterial resolve_insulation_layer_insulation_material(Layer layer);
        double get_insulation_section_thickness(std::string sectionName);
        static double get_insulation_section_thickness(Coil coil, std::string sectionName);

        double get_insulation_layer_thickness(Layer layer);
        double get_insulation_layer_thickness(std::string layerName);
        static double get_insulation_layer_thickness(Coil coil, std::string layerName);

        double get_insulation_layer_relative_permittivity(Layer layer);
        double get_insulation_layer_relative_permittivity(std::string layerName);
        static double get_insulation_layer_relative_permittivity(Coil coil, std::string layerName);
        double get_insulation_section_relative_permittivity(std::string sectionName);
        static double get_insulation_section_relative_permittivity(Coil coil, std::string sectionName);

        std::vector<double> get_turns_ratios();
        std::vector<double> get_maximum_dimensions();

        static std::vector<std::vector<size_t>> get_patterns(Inputs& inputs, CoreType coreType);
        static std::vector<size_t> get_repetitions(Inputs& inputs, CoreType coreType);
        std::pair<std::vector<size_t>, size_t> check_pattern_and_repetitions_integrity(std::vector<size_t> pattern, size_t repetitions);

        bool is_edge_wound_coil();

};
}
namespace OpenMagnetics {

void from_json(const json & j, Coil & x);
void to_json(json & j, const Coil & x);
void from_json(const json & j, CoilFunctionalDescription & x);
void to_json(json & j, const CoilFunctionalDescription & x);
} // namespace OpenMagnetics


namespace nlohmann {
template <>
struct adl_serializer<std::variant<OpenMagnetics::Bobbin, std::string>> {
    static void from_json(const json & j, std::variant<OpenMagnetics::Bobbin, std::string> & x);
    static void to_json(json & j, const std::variant<OpenMagnetics::Bobbin, std::string> & x);
};

template <>
struct adl_serializer<std::variant<OpenMagnetics::Wire, std::string>> {
    static void from_json(const json & j, std::variant<OpenMagnetics::Wire, std::string> & x);
    static void to_json(json & j, const std::variant<OpenMagnetics::Wire, std::string> & x);
};
} // namespace nlohmann


namespace OpenMagnetics {
inline void from_json(const json & j, Coil& x) {
    x.set_bobbin(j.at("bobbin").get<OpenMagnetics::BobbinDataOrNameUnion>());
    x.set_functional_description(j.at("functionalDescription").get<std::vector<CoilFunctionalDescription>>());
    x.set_layers_description(get_stack_optional<std::vector<Layer>>(j, "layersDescription"));
    x.set_sections_description(get_stack_optional<std::vector<Section>>(j, "sectionsDescription"));
    x.set_turns_description(get_stack_optional<std::vector<Turn>>(j, "turnsDescription"));
    x.set_groups_description(get_stack_optional<std::vector<Group>>(j, "groupsDescription"));
}

inline void from_json(const json & j, CoilFunctionalDescription& x) {
    x.set_connections(get_stack_optional<std::vector<ConnectionElement>>(j, "connections"));
    x.set_isolation_side(j.at("isolationSide").get<IsolationSide>());
    x.set_name(j.at("name").get<std::string>());
    x.set_number_parallels(j.at("numberParallels").get<int64_t>());
    x.set_number_turns(j.at("numberTurns").get<int64_t>());
    x.set_wire(j.at("wire").get<OpenMagnetics::WireDataOrNameUnion>());
}

inline void to_json(json & j, const Coil & x) {
    j = json::object();
    j["bobbin"] = x.get_bobbin();
    j["functionalDescription"] = x.get_functional_description();
    j["layersDescription"] = x.get_layers_description();
    j["sectionsDescription"] = x.get_sections_description();
    j["turnsDescription"] = x.get_turns_description();
    j["groupsDescription"] = x.get_groups_description();
}

inline void to_json(json & j, const CoilFunctionalDescription & x) {
    j = json::object();
    j["connections"] = x.get_connections();
    j["isolationSide"] = x.get_isolation_side();
    j["name"] = x.get_name();
    j["numberParallels"] = x.get_number_parallels();
    j["numberTurns"] = x.get_number_turns();
    j["wire"] = x.get_wire();
}
} // namespace OpenMagnetics


namespace nlohmann {
inline void adl_serializer<std::variant<OpenMagnetics::Bobbin, std::string>>::from_json(const json & j, std::variant<OpenMagnetics::Bobbin, std::string> & x) {
    if (j.is_string())
        x = j.get<std::string>();
    else if (j.is_object())
        x = j.get<OpenMagnetics::Bobbin>();
    else throw std::runtime_error("Could not deserialise!");
}

inline void adl_serializer<std::variant<OpenMagnetics::Bobbin, std::string>>::to_json(json & j, const std::variant<OpenMagnetics::Bobbin, std::string> & x) {
    switch (x.index()) {
        case 0:
            j = std::get<OpenMagnetics::Bobbin>(x);
            break;
        case 1:
            j = std::get<std::string>(x);
            break;
        default: throw std::runtime_error("Input JSON does not conform to schema!");
    }
}
inline void adl_serializer<std::variant<OpenMagnetics::Wire, std::string>>::from_json(const json & j, std::variant<OpenMagnetics::Wire, std::string> & x) {
    if (j.is_string())
        x = j.get<std::string>();
    else if (j.is_object())
        x = j.get<OpenMagnetics::Wire>();
    else throw std::runtime_error("Could not deserialise!");
}

inline void adl_serializer<std::variant<OpenMagnetics::Wire, std::string>>::to_json(json & j, const std::variant<OpenMagnetics::Wire, std::string> & x) {
    switch (x.index()) {
        case 0:
            j = std::get<OpenMagnetics::Wire>(x);
            break;
        case 1:
            j = std::get<std::string>(x);
            break;
        default: throw std::runtime_error("Input JSON does not conform to schema!");
    }
}
} // namespace nlohmann
