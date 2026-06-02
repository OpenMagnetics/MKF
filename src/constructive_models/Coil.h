#pragma once

#include "constructive_models/MasMigration.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Core.h"
#include "processors/Inputs.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Bobbin.h"
#include "json.hpp"

#include <MAS.hpp>
#include <vector>
#include "support/Exceptions.h"

using namespace MAS;
using json = nlohmann::json;

namespace OpenMagnetics {

using BobbinDataOrNameUnion = std::variant<Bobbin, std::string>;
using WireDataOrNameUnion = std::variant<Wire, std::string>;


class Winding : public MAS::CoilFunctionalDescription {
    private:
        WireDataOrNameUnion wire = std::string("Dummy");  // Initialize with valid default

    public:
        const WireDataOrNameUnion & get_wire() const { return wire; }
        WireDataOrNameUnion & get_mutable_wire() { return wire; }
        void set_wire(const WireDataOrNameUnion & value) { this->wire = value; }

        void set_isolation_side_from_index(size_t windingIndex);
        Winding(const MAS::CoilFunctionalDescription& winding) {
            if (winding.get_connections()) {
                set_connections(winding.get_connections());
            }

            set_isolation_side(winding.get_isolation_side());
            set_name(winding.get_name());
            set_number_parallels(winding.get_number_parallels() > 0 ? winding.get_number_parallels() : 1);
            set_number_turns(winding.get_number_turns());
            set_wound_with(winding.get_wound_with());
            auto wireVariant = winding.get_wire();
            if (std::holds_alternative<std::string>(wireVariant)) {
                set_wire(std::get<std::string>(wireVariant));
            }
            else {
                set_wire(std::get<MAS::Wire>(wireVariant));
            }
        };
        
        // Default constructor - initialize all members
        Winding() : MAS::CoilFunctionalDescription() {
            set_isolation_side(IsolationSide::PRIMARY);
            set_name("");
            set_number_parallels(1);
            set_number_turns(1);
            // Parent's wire is set via base class, child's wire is default-initialized above
            MAS::CoilFunctionalDescription::set_wire(std::string("Dummy"));
        }
        
        // Copy constructor
        Winding(const Winding& other) : MAS::CoilFunctionalDescription(other), wire(other.wire) {}
        
        // Move constructor
        Winding(Winding&& other) noexcept : MAS::CoilFunctionalDescription(std::move(other)), wire(std::move(other.wire)) {}
        
        // Copy assignment
        Winding& operator=(const Winding& other) {
            if (this != &other) {
                MAS::CoilFunctionalDescription::operator=(other);
                wire = other.wire;
            }
            return *this;
        }
        
        // Move assignment
        Winding& operator=(Winding&& other) noexcept {
            if (this != &other) {
                MAS::CoilFunctionalDescription::operator=(std::move(other));
                wire = std::move(other.wire);
            }
            return *this;
        }
        
        virtual ~Winding() = default;

        Wire resolve_wire();
};

// One rectangle of radial space reserved by a terminal/connection lead crossing a layer boundary,
// used when real winding geometry is enabled (Settings::get_coil_use_real_winding_geometry). It
// both feeds the section filling factor and is drawn by the Painter for debugging.
struct ConnectionReservedSpace {
    std::string section;
    std::string layer;                // the conduction layer this lead squeezes (empty for terminal leads)
    std::string winding;              // the winding whose lead reserves the space
    std::vector<double> coordinates;  // centre of the reserved rectangle (same system as turns)
    std::vector<double> dimensions;   // {width, height}
    double rotation = 0;              // degrees, for diagonal links (Z continuations); 0 = axis-aligned
    // A terminal lead routes a winding end out to the bobbin window border (entrance/exit). It is
    // drawn and its length feeds the connection loss, but it does not squeeze a conduction layer.
    // A non-terminal lead is an inter-layer transition that squeezes the crossed layer.
    bool isTerminal = false;
};

class Coil : public MAS::Coil {
    private:
        std::map<std::pair<size_t, size_t>, Section> _insulationSections;
        std::map<std::pair<size_t, size_t>, std::vector<Layer>> _insulationInterSectionsLayers;
        std::map<size_t, Layer> _insulationInterLayers;
        std::map<std::pair<size_t, size_t>, CoilSectionInterface> _coilSectionInterfaces;
        std::map<std::pair<size_t, size_t>, std::string> _insulationSectionsLog;
        std::map<std::pair<size_t, size_t>, std::string> _insulationInterSectionsLayersLog;
        std::map<std::string, size_t> _windingIndexByName;
        std::map<std::string, size_t> _turnIndexByName;
        std::map<std::string, Turn> _turnByName;
        std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> _sectionInfoWithInsulation;
        std::vector<std::vector<double>> _marginsPerSection;
        size_t _interleavingLevel = 1;
        WindingOrientation _windingOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation _layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment _turnsAlignment = CoilAlignment::CENTERED;
        CoilAlignment _sectionAlignment = CoilAlignment::INNER_OR_TOP;
        bool _sectionAlignmentExplicit = false;
        std::optional<Inputs> _inputs;
        std::map<std::string, CoilAlignment> _turnsAlignmentPerSection;
        std::map<std::string, WindingOrientation> _layersOrientationPerSection;
        // Real-winding turn blocking (global to the winding window). Filled by wind() between
        // re-wind iterations: maps a conduction layer name to the number of connection-lead slots
        // blocked at its {top, bottom}. Consumed by wind_by_rectangular_layers when
        // _applyConnectionBlocking is set. Both stay empty/false unless real winding geometry is on.
        std::map<std::string, std::pair<uint64_t, uint64_t>> _connectionBlockedSlotsPerLayer;
        bool _applyConnectionBlocking = false;
        std::string coilLog;
        InsulationCoordinator _standardCoordinator = InsulationCoordinator();
        std::vector<double> _currentProportionPerWinding;
        std::vector<size_t> _currentPattern;
        size_t _currentRepetitions = 1;
        bool _strict = true;
        bool _bobbin_resolved = false;
        Bobbin _bobbin;
        BobbinDataOrNameUnion bobbin;
        std::vector<Winding> functional_description;

        bool wind_by_rectangular_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions);
        bool wind_by_round_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions);
        bool wind_by_rectangular_layers();
        bool wind_by_round_layers();
        bool wind_by_rectangular_turns();
        bool wind_by_round_turns();
        bool wind_toroidal_additional_turns();
        // Special case: toroid with exactly one physical turn whose wire
        // outer diameter exceeds the inner-hole radius (so it cannot be
        // wound on the inner wall). Builds sections / layers / turns with
        // the wire placed at the geometric centre of the inner hole and
        // returns true on success; returns false if the case does not
        // apply or the wire physically cannot fit inside the hole.
        bool can_build_centered_single_turn_toroidal();
        bool build_centered_single_turn_toroidal();
        bool delimit_and_compact_rectangular_window();
        bool delimit_and_compact_round_window();
        bool create_default_group(Bobbin bobbin, WiringTechnology coilType = WiringTechnology::WOUND, double coreToLayerDistance = 0);
        // Multi-window dispatcher: calls create_default_group for single-window
        // bobbins, or creates one Group per winding window for multi-window
        // bobbins (all windings placed in column 0 by default — call
        // assign_windings_to_columns() to override).
        bool create_default_groups(Bobbin bobbin, WiringTechnology coilType = WiringTechnology::WOUND, double coreToLayerDistance = 0);
        // Look up which winding window a group is anchored to, by matching
        // the group's coordinates against each winding window. Returns 0 if
        // no match (safe fallback for single-column behaviour).
        size_t find_window_index_for_group(const std::string& groupName) const;

    public:
        // Distribute windings across the bobbin's winding windows. Creates one
        // Group per winding window in the bobbin and assigns windings to
        // groups per the provided indices. The outer vector size MUST equal
        // the number of winding windows; each inner vector lists the
        // functionalDescription indices to place in that column. Sets
        // groupsDescription. Must be called BEFORE wind_by_sections() if you
        // want non-default column placement.
        void assign_windings_to_columns(const std::vector<std::vector<size_t>>& windingIndicesPerColumn);
        bool wind_by_planar_sections(std::vector<size_t> stackUp, std::map<std::pair<size_t, size_t>, double> insulationThickness = {}, double coreToLayerDistance = 0);
        bool wind_by_planar_layers();
        bool wind_by_planar_turns(double borderToWireDistance, std::map<size_t, double> wireToWireDistance);

        Coil(json j, size_t interleavingLevel = 1,
                       WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING,
                       WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING,
                       CoilAlignment turnsAlignment = CoilAlignment::CENTERED,
                       CoilAlignment sectionAlignment = CoilAlignment::INNER_OR_TOP);
        Coil(const MAS::Coil coil);
        Coil(json j, bool windInConstructor);
        Coil() = default;
        virtual ~Coil() = default;
        bool fast_wind();
        bool unwind();
        bool wind();
        bool wind(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions=1);
        bool wind(std::vector<size_t> pattern, size_t repetitions=1);
        bool wind(size_t repetitions);
        bool wind_planar(std::vector<size_t> stackUp, std::optional<double> borderToWireDistance = std::nullopt, std::map<size_t, double> wireToWireDistance = {}, std::map<std::pair<size_t, size_t>, double> insulationThickness = {}, double coreToLayerDistance = 0);
        void try_rewind();
        void clear();
        bool are_sections_and_layers_fitting();

        const BobbinDataOrNameUnion & get_bobbin() const { return bobbin; }
        BobbinDataOrNameUnion & get_mutable_bobbin() { return bobbin; }
        void set_bobbin(const BobbinDataOrNameUnion & value) { this->bobbin = value; }

        const std::vector<Winding> & get_functional_description() const { return functional_description; }
        std::vector<Winding> & get_mutable_functional_description() { return functional_description; }
        void set_functional_description(const std::vector<Winding> & value) { this->functional_description = value; }


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
        void reset_insulation();
        size_t get_interleaving_level() const;
        void set_winding_orientation(WindingOrientation windingOrientation);
        void set_layers_orientation(WindingOrientation layersOrientation, std::optional<std::string> sectionName = std::nullopt);
        void set_turns_alignment(CoilAlignment turnsAlignment, std::optional<std::string> sectionName = std::nullopt);
        void set_section_alignment(CoilAlignment sectionAlignment);

        WindingOrientation get_winding_orientation();
        WindingOrientation get_layers_orientation() const;
        CoilAlignment get_turns_alignment(std::optional<std::string> sectionName = std::nullopt) const;
        CoilAlignment get_section_alignment();

        // Resolves the U/Z winding order for a section: the section's own windingOrder if set,
        // else the bobbin winding window's windingOrder, else WindingOrder::Z (current behaviour).
        WindingOrder get_winding_order(const std::string& sectionName) const;

        // Real winding geometry (only meaningful when there is more than one conduction layer in a
        // section). Returns the rectangles of space reserved by terminal/connection leads crossing
        // layer boundaries. Computed from the wound layers; independent of the real-geometry
        // setting so the Painter can also overlay it for debugging.
        std::vector<ConnectionReservedSpace> get_connection_reserved_spaces();
        // Adds the reserved-connection area into the affected section filling factors. Called at the
        // end of wind() when Settings::get_coil_use_real_winding_geometry() is true.
        void apply_connection_reserved_space();
        // Counts, per conduction layer, how many connection leads cross its {top, bottom} — derived
        // from get_connection_reserved_spaces() and the wound layer centres. This is the global
        // (window-wide) turn-blocking incidence wind() iterates on; a lead blocks any layer it
        // crosses regardless of section/winding.
        std::map<std::string, std::pair<uint64_t, uint64_t>> compute_connection_blocked_slots_per_layer();
        // Real-winding blocking makes a section's interior layers lose top/bottom slots, so an even
        // interleaving turn split leaves orphan turns in a near-empty spillover layer. Re-split each
        // winding's turns across its conduction sections (radial order) so interior sections fill
        // complete blocked layers and the remainder is pushed to the outermost section. Single
        // parallel only (bifilar needs per-parallel connections). Called in the wind() blocking pass.
        void redistribute_section_turns_for_blocking();
        // After delimit (which re-centres every layer), shift each blocked conduction layer's turns to
        // the UNblocked edge, so the slots freed by turn-blocking sit exactly where the connection
        // leads run (top edge for top-crossing leads, bottom for bottom) and no window space is wasted
        // with an empty gap on the unblocked side. Real winding geometry only.
        void align_blocked_layer_turns();

        std::vector<Section> get_sections_description_conduction() const;
        std::vector<Layer> get_layers_description_conduction() const;
        std::vector<Section> get_sections_description_insulation() const;
        std::vector<Layer> get_layers_description_insulation() const;

        std::string get_name(size_t windingIndex) const;

        WiringTechnology get_coil_type(size_t groupIndex = 0) const;

        std::vector<uint64_t> get_number_turns() const;
        uint64_t get_number_turns(size_t windingIndex) const;
        uint64_t get_number_turns(Section section);
        uint64_t get_number_turns(Layer layer);
        void set_number_turns(std::vector<uint64_t> numberTurns);
        std::vector<IsolationSide> get_isolation_sides() const;
        void set_isolation_sides(std::vector<IsolationSide> isolationSides);
        std::vector<uint64_t> get_number_parallels() const;
        uint64_t get_number_parallels(size_t windingIndex) const;
        void set_number_parallels(std::vector<uint64_t> numberParallels);

        void set_interlayer_insulation(double layerThickness, std::optional<std::string> material = std::nullopt, std::optional<std::string> windingName = std::nullopt, bool autowind=true);
        void set_intersection_insulation(double layerThickness, size_t numberInsulationLayers, std::optional<std::string> material = std::nullopt, std::optional<std::pair<std::string, std::string>> windingNames = std::nullopt, bool autowind=true);

        std::vector<Section> get_sections_by_group(std::string groupName) const;
        const std::vector<Section> get_sections_by_type(ElectricalType electricalType) const;
        const Section get_section_by_name(std::string name) const;
        Turn get_turn_by_name(std::string name);
        const std::vector<Section> get_sections_by_winding(std::string windingName) const;

        std::vector<Layer> get_layers_by_section(std::string sectionName) const;
        const std::vector<Layer> get_layers_by_type(ElectricalType electricalType) const;
        std::vector<Layer> get_layers_by_winding_index(size_t windingIndex);
        const Layer get_layer_by_name(std::string name) const;

        std::vector<Turn> get_turns_by_layer(std::string layerName) const;
        std::vector<Turn> get_turns_by_section(std::string sectionName) const;
        std::vector<Turn> get_turns_by_winding(std::string windingName) const;

        std::vector<std::string> get_layers_names_by_winding(std::string windingName) const;
        std::vector<std::string> get_layers_names_by_section(std::string sectionName) const;
        std::vector<std::string> get_turns_names_by_layer(std::string layerName) const;
        std::vector<std::string> get_turns_names_by_section(std::string sectionName) const;
        std::vector<std::string> get_turns_names_by_winding(std::string windingName) const;

        std::vector<size_t> get_turns_indexes_by_layer(std::string layerName) const;
        std::vector<size_t> get_turns_indexes_by_section(std::string sectionName) const;
        std::vector<size_t> get_turns_indexes_by_winding(std::string windingName) const;

        Winding get_winding_by_name(std::string name) const;


        size_t get_winding_index_by_name(const std::string& name);
        static size_t get_winding_index_by_name(const std::vector<Winding>& functionalDescription, const std::string& name);
        size_t get_turn_index_by_name(std::string name);
        size_t get_layer_index_by_name(std::string name) const;
        size_t get_section_index_by_name(std::string name) const;

        std::vector<Wire> get_wires();
        void set_wires(std::vector<Wire> wires);
        WireType get_wire_type(size_t windingIndex);
        std::vector<double> get_wires_length() const;
        static WireType get_wire_type(Winding winding);
        std::string get_wire_name(size_t windingIndex);
        static std::string get_wire_name(Winding winding);
        Wire resolve_wire(size_t windingIndex);
        static Wire resolve_wire(Winding winding);
        std::vector<double> resolve_margin(size_t sectionIndex);
        static std::vector<double> resolve_margin(Section section);
        static std::vector<double> resolve_margin(Margin marginVariant);
        MarginInfo resolve_margin_info(size_t sectionIndex);
        static MarginInfo resolve_margin_info(Section section);
        static MarginInfo resolve_margin_info(Margin marginVariant);

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

        std::vector<double> get_turns_ratios() const;
        std::vector<double> get_maximum_dimensions();

        static std::vector<std::vector<size_t>> get_patterns(Inputs& inputs, CoreType coreType);
        static std::vector<size_t> get_repetitions(Inputs& inputs, CoreType coreType);
        std::pair<std::vector<size_t>, size_t> check_pattern_and_repetitions_integrity(std::vector<size_t> pattern, size_t repetitions);

        bool is_edge_wound_coil();

        std::vector<Wire> guess_round_wire_from_dc_resistance(std::vector<double> dcResistances, double maxError=0.05);

        // Virtualization
        std::map<size_t, std::string> _windingNames;
        std::map<size_t, std::string> _virtualWindingNames;
        std::map<size_t, std::vector<size_t>> _virtualizationMap;
        bool needs_virtualization();
        std::map<size_t, std::vector<size_t>> create_virtualization_map();
        size_t get_winding_group_minimum_index(size_t currentWindingIndex);
        std::vector<size_t> virtualize_pattern(std::vector<size_t> pattern);
        std::vector<Winding> virtualize_functional_description();
        std::vector<size_t> compress_pattern(std::vector<size_t> pattern);
        std::vector<double> virtualize_proportion_per_winding(std::vector<double> proportionPerWinding);
        Section devirtualize_section(Section section);
        Layer devirtualize_layer(Layer layer);
        Turn devirtualize_turn(Turn turn, std::string virtualWindingName, std::string windingName, size_t newParallelIndex);
        Section virtualize_section(Section section);
        Layer virtualize_layer(Layer layer);
        Turn virtualize_turn(Turn turn, std::string virtualWindingName, std::string windingName);
        void devirtualize_sections_description();
        void devirtualize_layers_description();
        void devirtualize_turns_description();
        std::vector<Section> virtualize_sections_description();
        std::vector<Layer> virtualize_layers_description();
        std::vector<Turn> virtualize_turns_description();
        std::map<std::pair<size_t, size_t>, std::vector<Layer>> virtualize_insulation_intersections_layers();

        static Coil create_quick_coil(std::string coreShapeName,
                                      std::vector<int64_t> numberTurns,
                                      std::vector<int64_t> numberParallels = {},
                                      std::vector<OpenMagnetics::Wire> wires = {},
                                      WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING,
                                      WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING,
                                      CoilAlignment turnsAlignment = CoilAlignment::CENTERED,
                                      CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP,
                                      uint8_t interleavingLevel = 1,
                                      bool useBobbin = true,
                                      int numberStacks = 1);

        std::vector<size_t> extract_stack_up(std::vector<Section> sections);
        bool is_planar();
        std::vector<Turn> get_turns_touching_bobbin_column();
        std::vector<Turn> get_turns_touching_bobbin_walls();
        std::vector<Turn> get_turns_touching_bobbin_column(std::vector<Turn> turns);
        std::vector<Turn> get_turns_touching_bobbin_walls(std::vector<Turn> turns);
        std::vector<Turn> get_turns_touching_bobbin_column(std::vector<size_t> turnIndexes);
        std::vector<Turn> get_turns_touching_bobbin_walls(std::vector<size_t> turnIndexes);
                                               
};
}
namespace OpenMagnetics {

void from_json(const json & j, Coil & x);
void to_json(json & j, const Coil & x);
void from_json(const json & j, Winding & x);
void to_json(json & j, const Winding & x);
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
    x.set_functional_description(j.at("functionalDescription").get<std::vector<Winding>>());
    x.set_layers_description(get_stack_optional<std::vector<Layer>>(j, "layersDescription"));
    x.set_sections_description(get_stack_optional<std::vector<Section>>(j, "sectionsDescription"));
    x.set_turns_description(get_stack_optional<std::vector<Turn>>(j, "turnsDescription"));
    x.set_groups_description(get_stack_optional<std::vector<Group>>(j, "groupsDescription"));
}

inline void from_json(const json & j, Winding& x) {
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

inline void to_json(json & j, const Winding & x) {
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
        default: throw std::runtime_error("Input JSON does not conform to schema! [Coil.h:482 Bobbin variant]");
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
        default: throw std::runtime_error("Input JSON does not conform to schema! [Coil.h:501 Wire variant]");
    }
}
} // namespace nlohmann
