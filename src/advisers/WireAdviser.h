#pragma once
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

void from_json(const json & j, WireSolidInsulationRequirements & x);
void to_json(json & j, const WireSolidInsulationRequirements & x);

inline void from_json(const json & j, WireSolidInsulationRequirements& x) {
    x.set_minimum_number_layers(get_stack_optional<int64_t>(j, "minimumNumberLayers"));
    x.set_minimum_grade(get_stack_optional<int64_t>(j, "minimumGrade"));
    x.set_maximum_number_layers(get_stack_optional<int64_t>(j, "maximumNumberLayers"));
    x.set_maximum_grade(get_stack_optional<int64_t>(j, "maximumGrade"));
    x.set_minimum_breakdown_voltage(j.at("minimumBreakdownVoltage").get<double>());
}

inline void to_json(json & j, const WireSolidInsulationRequirements & x) {
    j = json::object();
    j["minimumNumberLayers"] = x.get_minimum_number_layers();
    j["minimumGrade"] = x.get_minimum_grade();
    j["maximumNumberLayers"] = x.get_maximum_number_layers();
    j["maximumGrade"] = x.get_maximum_grade();
    j["minimumBreakdownVoltage"] = x.get_minimum_breakdown_voltage();
}

/**
 * @class WireAdviser
 * @brief Recommends optimal wire types and configurations for magnetic windings.
 *
 * ## Overview
 * WireAdviser selects the best wire (round, litz, foil, rectangular, planar) and
 * parallel configuration for a given winding based on electrical requirements,
 * geometric constraints, and AC loss considerations.
 *
 * ## Scoring System
 * Each filter contributes a normalized score (0-1) that is summed to produce
 * a final ranking. Lower raw values (resistance, losses) result in higher scores.
 *
 * ## Filter Pipeline
 *
 * ### Standard Wires (round, litz, foil, rectangular)
 * Applied in order via `get_advised_wire()`:
 * 1. **filter_by_area_no_parallels**: Pre-filter eliminating wires too large for section
 * 2. **filter_by_solid_insulation_requirements**: Validates insulation grade/layers (if specified)
 * 3. **filter_by_area_with_parallels**: Validates wire fits with parallel configuration
 * 4. **filter_by_effective_resistance**: Scores by AC resistance (skin effect included)
 * 5. **filter_by_skin_losses_density**: Scores by skin effect power density
 * 6. **filter_by_proximity_factor**: Scores by proximity effect susceptibility
 *
 * ### Planar Wires
 * Applied in order via `get_advised_planar_wire()`:
 * 1. **filter_by_effective_resistance**: Scores by AC resistance
 * 2. **filter_by_skin_losses_density**: Scores by skin effect losses
 * 3. **filter_by_proximity_factor**: Scores by proximity effect
 *
 * ## Key Parameters
 * - **maximumEffectiveCurrentDensity**: Maximum allowed J (A/m²) in conductor
 * - **maximumNumberParallels**: Maximum parallel strands/wires allowed
 * - **wireSolidInsulationRequirements**: Insulation grade/layer requirements for safety
 *
 * ## Wire Type Selection
 * The adviser respects global settings to include/exclude wire types:
 * - `settings.set_wire_adviser_include_round(bool)`
 * - `settings.set_wire_adviser_include_litz(bool)`
 * - `settings.set_wire_adviser_include_foil(bool)`
 * - `settings.set_wire_adviser_include_rectangular(bool)`
 * - `settings.set_wire_adviser_include_planar(bool)`
 *
 * ## Usage Example
 * ```cpp
 * WireAdviser wireAdviser;
 * wireAdviser.set_maximum_effective_current_density(5e6);  // 5 A/mm²
 * wireAdviser.set_maximum_number_parallels(4);
 * auto results = wireAdviser.get_advised_wire(winding, section, current, temp, numSections, 5);
 * ```
 *
 * ## Industry Background
 * - Skin depth: δ = √(ρ/(π·f·μ)) determines AC current distribution
 * - Proximity effect: Increases with layer count and conductor diameter
 * - Litz wire: Reduces skin/proximity losses via transposed fine strands
 * - Foil: Low DC resistance but requires careful interleaving for AC
 */
class WireAdviser {
    protected:
        double _maximumEffectiveCurrentDensity;
        std::optional<WireSolidInsulationRequirements> _wireSolidInsulationRequirements;
        std::optional<WireStandard> _commonWireStandard;
        int _maximumNumberParallels;
        double _maximumOuterAreaProportion;
        double _wireToWireDistance = defaults.minimumWireToWireDistance;
        double _borderToWireDistance = defaults.minimumBorderToWireDistance;
        std::string _log;
    public:

        WireAdviser(double maximumEffectiveCurrentDensity, int maximumNumberParallels) {
            _maximumEffectiveCurrentDensity = maximumEffectiveCurrentDensity;
            _maximumNumberParallels = maximumNumberParallels;
        }
        WireAdviser() {
            auto defaults = Defaults();
            _maximumEffectiveCurrentDensity = defaults.maximumEffectiveCurrentDensity;
            _maximumNumberParallels = defaults.maximumNumberParallels;
        }
        virtual ~WireAdviser() = default;

        void set_wire_to_wire_distance(double wireToWireDistance) {
            _wireToWireDistance = wireToWireDistance;
        }
        double get_wire_to_wire_distance() {
            return _wireToWireDistance;
        }

        void set_border_to_wire_distance(double borderToWireDistance) {
            _borderToWireDistance = borderToWireDistance;
        }
        double get_border_to_wire_distance() {
            return _borderToWireDistance;
        }

        void set_maximum_effective_current_density(double maximumEffectiveCurrentDensity) {
            _maximumEffectiveCurrentDensity = maximumEffectiveCurrentDensity;
        }
        double get_maximum_effective_current_density() {
            return _maximumEffectiveCurrentDensity;
        }
        void set_wire_solid_insulation_requirements(WireSolidInsulationRequirements wireSolidInsulationRequirements) {
            _wireSolidInsulationRequirements = wireSolidInsulationRequirements;
        }
        void set_maximum_number_parallels(int maximumNumberParallels) {
            _maximumNumberParallels = maximumNumberParallels;
        }
        double get_maximum_number_parallels() {
            return _maximumNumberParallels;
        }
        double get_maximum_area_proportion() {
            return _maximumOuterAreaProportion;
        }
        void set_common_wire_standard(std::optional<WireStandard> commonWireStandard) {
            _commonWireStandard = commonWireStandard;
        }
        std::optional<WireStandard> get_common_wire_standard() {
            return _commonWireStandard;
        }
        std::vector<std::pair<Winding, double>> get_advised_wire(Winding winding,
                                                                                   Section section,
                                                                                   SignalDescriptor current,
                                                                                   double temperature,
                                                                                   uint8_t numberSections,
                                                                                   size_t maximumNumberResults=1);
        std::vector<std::pair<Winding, double>> get_advised_wire(std::vector<Wire>* wires,
                                                                                   Winding winding,
                                                                                   Section section,
                                                                                   SignalDescriptor current,
                                                                                   double temperature,
                                                                                   uint8_t numberSections,
                                                                                   size_t maximumNumberResults=1);

        std::vector<std::pair<Winding, double>> get_advised_planar_wire(Winding winding,
                                                                                          Section section,
                                                                                          SignalDescriptor current,
                                                                                          double temperature,
                                                                                          uint8_t numberSections,
                                                                                          size_t maximumNumberResults);

        std::vector<std::pair<Winding, double>> filter_by_area_no_parallels(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                              Section section);

        std::vector<std::pair<Winding, double>> filter_by_area_with_parallels(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                Section section,
                                                                                                double numberSections,
                                                                                                bool allowNotFit);

        std::vector<std::pair<Winding, double>> filter_by_effective_resistance(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                 SignalDescriptor current,
                                                                                                 double temperature);

        std::vector<std::pair<Winding, double>> filter_by_skin_losses_density(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                SignalDescriptor current,
                                                                                                double temperature);

        std::vector<std::pair<Winding, double>> filter_by_proximity_factor(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                 SignalDescriptor current,
                                                                                                 double temperature);

        std::vector<std::pair<Winding, double>> filter_by_solid_insulation_requirements(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                 WireSolidInsulationRequirements wireSolidInsulationRequirements);

        std::vector<std::pair<Winding, double>> create_dataset(Winding winding,
                                                                                 std::vector<Wire>* wires,
                                                                                 Section section,
                                                                                 SignalDescriptor current,
                                                                                 double temperature);
        std::vector<std::pair<Winding, double>> create_planar_dataset(Winding winding,
                                                                                        Section section,
                                                                                        uint8_t numberSections);
        void expand_wires_dataset_with_parallels(std::vector<Winding>* windings);
        void set_maximum_area_proportion(std::vector<std::pair<Winding, double>>* unfilteredCoils, Section section, uint8_t numberSections);
    
};


} // namespace OpenMagnetics