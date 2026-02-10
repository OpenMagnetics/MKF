#pragma once
#include "advisers/MagneticFilter.h"
#include "advisers/WireAdviser.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Mas.h"
#include "support/Utils.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

/**
 * @class CoilAdviser
 * @brief Recommends complete coil configurations including winding patterns, wire selection, and insulation.
 *
 * ## Overview
 * CoilAdviser extends WireAdviser to recommend complete coil designs. It handles:
 * - Winding pattern selection (interleaved vs. non-interleaved)
 * - Wire selection per winding (using WireAdviser)
 * - Section proportion calculation based on power handling
 * - Insulation coordination per safety standards
 *
 * ## Design Process
 * 1. Calculate winding window proportions based on average power per winding
 * 2. Generate candidate patterns (winding order permutations)
 * 3. For each pattern, determine insulation requirements
 * 4. Select optimal wires for each winding using WireAdviser
 * 5. Score complete coil configurations and return ranked results
 *
 * ## Scoring System
 * The default scoring filters (configurable via `load_filter_flow()`):
 * - **EFFECTIVE_RESISTANCE**: AC resistance of windings (lower = better)
 * - **EFFECTIVE_CURRENT_DENSITY**: Current density in conductors (lower = better)
 * - **MAGNETOMOTIVE_FORCE**: MMF distribution quality (lower = better)
 *
 * Each filter uses:
 * - `invert=true`: Lower raw values get higher scores
 * - `log=true`: Logarithmic normalization (compresses large differences)
 * - `weight=1.0`: Equal weight for all criteria
 *
 * ## Winding Patterns
 * For a 2-winding transformer, patterns include:
 * - `{0, 1}`: Primary then Secondary (non-interleaved)
 * - `{1, 0}`: Secondary then Primary
 * - With `repetitions > 1`: Interleaved sections (P-S-P-S, etc.)
 *
 * ## Insulation Coordination
 * Based on IEC 60664-1 and IEC 61558, determines:
 * - Required creepage/clearance distances
 * - Wire insulation grade requirements
 * - Whether margin tape or insulated wire is needed
 *
 * ## Planar vs. Wound Coils
 * - **Wound**: Traditional bobbin-wound coils with round/litz/foil wire
 * - **Planar**: PCB-based windings with copper traces
 *
 * ## Key Configuration
 * - `set_allow_margin_tape(bool)`: Allow/disallow margin tape insulation
 * - `set_allow_insulated_wire(bool)`: Allow/disallow triple-insulated wire
 * - `set_common_wire_standard(WireStandard)`: Restrict to specific wire standards
 * - `set_maximum_effective_current_density(double)`: Max current density
 * - `set_maximum_number_parallels(int)`: Max parallel conductors
 *
 * ## Usage Example
 * ```cpp
 * CoilAdviser coilAdviser;
 * coilAdviser.set_maximum_effective_current_density(5e6);
 * coilAdviser.set_allow_margin_tape(true);
 * auto results = coilAdviser.get_advised_coil(mas, 5);  // Top 5 configurations
 * ```
 *
 * ## References
 * - IEC 60664-1: Insulation coordination for low-voltage equipment
 * - IEC 61558: Safety of transformers
 */
class CoilAdviser : public WireAdviser {
    private:
        bool _allowMarginTape = true;
        bool _allowInsulatedWire = true;
        std::map<MagneticFilters, std::shared_ptr<MagneticFilter>> _filters;
        std::vector<MagneticFilterOperation> _loadedFilterFlow;
        OpenMagnetics::WireAdviser _wireAdviser;
        std::optional<WireStandard> _commonWireStandard = defaults.commonWireStandard;
        std::vector<MagneticFilterOperation> _defaultCustomMagneticFilterFlow{
            MagneticFilterOperation(MagneticFilters::EFFECTIVE_RESISTANCE, true, true, 1.0),
            MagneticFilterOperation(MagneticFilters::EFFECTIVE_CURRENT_DENSITY, true, true, 1.0),
            MagneticFilterOperation(MagneticFilters::MAGNETOMOTIVE_FORCE, true, true, 1.0),
        };

    public:

        std::vector<Mas> get_advised_coil(Mas mas, size_t maximumNumberResults=1);
        std::vector<Mas> get_advised_coil(std::vector<Wire>* wires, Mas mas, size_t maximumNumberResults=1);
        std::vector<Section> get_advised_sections(Mas mas, std::vector<size_t> pattern, size_t repetitions);
        std::vector<Section> get_advised_planar_sections(Mas mas, std::vector<size_t> pattern, size_t repetitions);
        std::vector<Mas> get_advised_coil_for_pattern(std::vector<Wire>* wires, Mas mas, std::vector<size_t> pattern, size_t repetitions, std::vector<WireSolidInsulationRequirements> solidInsulationRequirementsForWires, size_t maximumNumberResults, std::string reference);
        std::vector<Mas> get_advised_planar_coil_for_pattern(std::vector<Wire>* wires, Mas mas, std::vector<size_t> pattern, size_t repetitions, size_t maximumNumberResults, std::string reference);
        std::vector<std::pair<Mas, double>> score_magnetics(std::vector<Mas> masMagnetics, std::vector<MagneticFilterOperation> filterFlow);
        void load_filter_flow(std::vector<MagneticFilterOperation> flow, std::optional<Inputs> inputs);
        
        void set_allow_margin_tape(bool value) {
            _allowMarginTape = value;
        }
        void set_allow_insulated_wire(bool value) {
            _allowInsulatedWire = value;
        }
        void set_common_wire_standard(std::optional<WireStandard> commonWireStandard) {
            _commonWireStandard = commonWireStandard;
        }
        std::optional<WireStandard> get_common_wire_standard() {
            return _commonWireStandard;
        }

        void set_maximum_effective_current_density(double maximumEffectiveCurrentDensity) {
            _wireAdviser.set_maximum_effective_current_density(maximumEffectiveCurrentDensity);
        }
        
        double get_maximum_effective_current_density() {
            return _wireAdviser.get_maximum_effective_current_density();
        }

        void set_maximum_number_parallels(int maximumNumberParallels) {
            _wireAdviser.set_maximum_number_parallels(maximumNumberParallels);
        }

        double get_maximum_number_parallels() {
            return _wireAdviser.get_maximum_number_parallels();
        }

};


} // namespace OpenMagnetics