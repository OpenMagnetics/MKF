#pragma once
#include "support/Utils.h"
#include "advisers/MagneticFilter.h"
#include "advisers/CoreAdviser.h"
#include "advisers/CoilAdviser.h"
#include "constructive_models/Mas.h"
#include "Definitions.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

/**
 * @class MagneticAdviser
 * @brief Top-level magnetic component design optimization system.
 *
 * MagneticAdviser orchestrates the complete magnetic component design process by:
 * 1. Selecting optimal cores via CoreAdviser
 * 2. Winding coils using CoilAdviser
 * 3. Simulating complete designs with MagneticSimulator
 * 4. Scoring and ranking results using configurable filter flows
 *
 * ## Multi-Objective Optimization Approach
 *
 * This adviser implements an **a priori scalarization method** for multi-objective
 * optimization, where user-specified weights define preferences before optimization.
 * The approach follows the linear scalarization formula:
 *
 *     total_score = Σ (weight_i × normalized_score_i)
 *
 * This method is computationally efficient and guarantees Pareto-optimal solutions
 * when the Pareto front is convex. For non-convex fronts, increasing the number
 * of requested results helps explore more of the design space.
 *
 * ## Default Filter Configuration
 *
 * For custom magnetics, the default filter flow optimizes:
 * - **COST**: Minimize material and manufacturing cost (log normalization)
 * - **LOSSES**: Minimize total power losses (log normalization)
 * - **DIMENSIONS**: Minimize physical volume (linear normalization)
 *
 * For catalogue magnetics, strictly required filters ensure compatibility:
 * - TURNS_RATIOS, MAXIMUM_DIMENSIONS, SATURATION, DC_CURRENT_DENSITY,
 *   EFFECTIVE_CURRENT_DENSITY, IMPEDANCE, MAGNETIZING_INDUCTANCE
 *
 * ## Common Mode Choke (CMC) Design
 *
 * For interference suppression applications with SubApplication::COMMON_MODE_NOISE_FILTERING,
 * the adviser automatically configures for CMC optimization:
 *
 * **Core Selection:**
 * - Restricts to toroidal cores (CoreShapeFamily::T) for optimal coupling
 * - Prefers high-permeability materials (nanocrystalline, MnZn ferrite)
 * - Evaluates cores based on impedance capability, not energy storage
 *
 * **Winding Configuration:**
 * - Enables bifilar winding (repetitions = {2, 1}) for tight coupling
 * - Both windings use IsolationSide::PRIMARY (same isolation level)
 * - Turns ratio is always 1:1 for balanced common mode rejection
 *
 * **Key Filters for CMC:**
 * - CORE_MINIMUM_IMPEDANCE: Ensures minimum impedance at specified frequencies
 * - LEAKAGE_INDUCTANCE: Minimizes Lk/Lm ratio for high coupling coefficient (k ≈ 1)
 * - Operating frequency must stay below SRF margin (default 25% of self-resonant frequency)
 *
 * **CMC Design Flow:**
 * ```cpp
 * Inputs inputs;
 * inputs.get_mutable_design_requirements().set_application(Application::INTERFERENCE_SUPPRESSION);
 * inputs.get_mutable_design_requirements().set_sub_application(SubApplication::COMMON_MODE_NOISE_FILTERING);
 * inputs.get_mutable_design_requirements().set_minimum_impedance(impedanceRequirements);
 * 
 * MagneticAdviser adviser;
 * auto results = adviser.get_advised_magnetic(inputs, 5);
 * ```
 *
 * ## Weight Guidelines
 *
 * | Application          | COST | LOSSES | DIMENSIONS |
 * |---------------------|------|--------|------------|
 * | Consumer electronics | 2.0  | 0.5    | 1.5        |
 * | High-efficiency PSU  | 0.5  | 2.0    | 1.0        |
 * | Space-constrained    | 0.5  | 1.0    | 2.0        |
 * | Balanced             | 1.0  | 1.0    | 1.0        |
 *
 * ## Usage Example
 *
 * ```cpp
 * MagneticAdviser adviser;
 * adviser.set_application(Application::POWER);
 * adviser.set_core_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
 *
 * // With custom weights
 * std::map<MagneticFilters, double> weights = {
 *     {MagneticFilters::COST, 0.5},
 *     {MagneticFilters::LOSSES, 2.0},
 *     {MagneticFilters::DIMENSIONS, 1.0}
 * };
 * auto results = adviser.get_advised_magnetic(inputs, weights, 5);
 * ```
 *
 * @see CoreAdviser For core selection algorithm details
 * @see CoilAdviser For winding optimization details
 * @see MagneticFilter For available filter types
 */
class MagneticAdviser{
    public:

        std::map<MagneticFilters, std::shared_ptr<MagneticFilter>> _filters;
        std::vector<MagneticFilterOperation> _loadedFilterFlow;
        /// @brief Default filter flow for custom magnetic design.
        /// COST and LOSSES use log normalization (spans orders of magnitude).
        /// DIMENSIONS uses linear normalization (intuitive volume comparison).
        std::vector<MagneticFilterOperation> _defaultCustomMagneticFilterFlow{
            MagneticFilterOperation(MagneticFilters::COST, true, true, 1.0),
            MagneticFilterOperation(MagneticFilters::LOSSES, true, true, 1.0),
            MagneticFilterOperation(MagneticFilters::DIMENSIONS, true, false, 1.0),  // Linear for DIMENSIONS
        };

        /// @brief Default filter flow for catalogue magnetic selection.
        /// These filters are strictly required to ensure compatibility with design requirements.
        std::vector<MagneticFilterOperation> _defaultCatalogueMagneticFilterFlow{
            MagneticFilterOperation(MagneticFilters::TURNS_RATIOS, true, false, true, 1.0),
            MagneticFilterOperation(MagneticFilters::MAXIMUM_DIMENSIONS, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::SATURATION, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::DC_CURRENT_DENSITY, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::EFFECTIVE_CURRENT_DENSITY, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::IMPEDANCE, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::MAGNETIZING_INDUCTANCE, true, false, 1.0),
        };
        bool _simulateResults = true;
        bool _uniqueCoreShapes = false;
        Application _application = Application::POWER;
        CoreAdviser::CoreAdviserModes _coreAdviserMode = CoreAdviser::CoreAdviserModes::AVAILABLE_CORES;

        MagneticAdviser() {
        }
        MagneticAdviser(bool simulateResults) {
            _simulateResults = simulateResults;
        }

        void set_unique_core_shapes(bool value);
        bool get_unique_core_shapes();
        void set_application(Application value);
        Application get_application();
        void set_core_mode(CoreAdviser::CoreAdviserModes value);
        CoreAdviser::CoreAdviserModes get_core_mode();

        /**
         * @brief Get optimized magnetic designs using default weights.
         * @param inputs Design requirements and operating conditions.
         * @param maximumNumberResults Maximum number of designs to return.
         * @return Vector of (Mas, score) pairs sorted by descending score.
         */
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, size_t maximumNumberResults=1);

        /**
         * @brief Get optimized magnetic designs with custom weights.
         * @param inputs Design requirements and operating conditions.
         * @param weights Map of filter types to weight values (higher = more important).
         * @param maximumNumberResults Maximum number of designs to return.
         * @return Vector of (Mas, score) pairs sorted by descending score.
         */
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::map<MagneticFilters, double> weights, size_t maximumNumberResults);

        /**
         * @brief Get optimized magnetic designs with custom filter flow.
         * @param inputs Design requirements and operating conditions.
         * @param filterFlow Custom sequence of filter operations with full configuration.
         * @param maximumNumberResults Maximum number of designs to return.
         * @return Vector of (Mas, score) pairs sorted by descending score.
         */
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::vector<Magnetic> catalogueMagnetics, size_t maximumNumberResults=1, bool strict=true);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::vector<Magnetic> catalogueMagnetics, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults=1, bool strict=true);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::map<std::string, Magnetic> catalogueMagnetics, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults=1, bool strict=true);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(std::vector<Mas> catalogueMagneticsWithInputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults=1, bool strict=true);

        /**
         * @brief Score a collection of magnetic designs using the filter flow.
         * @param masMagneticsWithCoil Vector of complete magnetic designs to score.
         * @param filterFlow Filter operations defining scoring criteria.
         * @return Vector of (Mas, score) pairs with normalized scores.
         */
        std::vector<std::pair<Mas, double>> score_magnetics(std::vector<Mas> masMagneticsWithCoil, std::vector<MagneticFilterOperation> filterFlow);

        /// @brief Print a human-readable summary of a magnetic design.
        static void preview_magnetic(Mas mas);

        /// @brief Load and initialize filter flow for subsequent operations.
        void load_filter_flow(std::vector<MagneticFilterOperation> flow, std::optional<Inputs> inputs = std::nullopt);

        /**
         * @brief Get per-filter normalized scores for all evaluated magnetics.
         * @return Map of reference -> (filter -> normalized_score).
         */
        std::map<std::string, std::map<MagneticFilters, double>> get_scorings();
};


} // namespace OpenMagnetics