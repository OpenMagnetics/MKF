#pragma once
#include "support/Utils.h"
#include "advisers/MagneticFilter.h"
#include "advisers/CoreAdviser.h"
#include "advisers/CoilAdviser.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "Definitions.h"
#include "Defaults.h"
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
 * inputs.get_mutable_design_requirements().set_application(MAS::MagneticApplication::INTERFERENCE_SUPPRESSION);
 * inputs.get_mutable_design_requirements().set_sub_application(SubApplication::COMMON_MODE_NOISE_FILTERING);
 * inputs.get_mutable_design_requirements().set_minimum_impedance(impedanceRequirements);
 * 
 * MagneticAdviser adviser;
 * auto results = adviser.get_advised_magnetic(inputs, 5);
 * ```
 *
 *
 * ## Known Limitations & Future Improvements
 *
 * **LOGIC-1: No Thermal Coupling** — Core temperature affects winding resistance
 * which affects losses which affects temperature. Current system evaluates independently.
 * TODO: Add 2-3 iteration thermal-electromagnetic coupling loop.
 *
 * **LOGIC-2: No Loss Balance Optimization** — Industry best practice targets 50/50
 * core-to-copper loss ratio. TODO: Add loss-balance metric to scoring.
 *
 * **LOGIC-3: Linear Scalarization Only** — Cannot find solutions on non-convex
 * Pareto fronts. TODO: For small result sets, add epsilon-constraint or NSGA-II.
 *
 * **LOGIC-5: No DC Bias in Material Selection** — Permeability drops significantly
 * under DC bias for powder cores. TODO: Include DC bias in material comparison.
 *
 * **LOGIC-6: No Proximity Effect in CoreAdviser** — Multi-layer proximity losses
 * can dominate but aren't estimated during core selection, leading to undersized cores.
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
 * adviser.set_application(MAS::MagneticApplication::POWER);
 * // Uses STANDARD_CORES mode by default (calculates optimal gaps)
 * // To use existing catalog cores instead:
 * // adviser.set_core_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
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
        // ABT #164: active type constraints for the in-flight call. Set by the
        // ctx-aware overloads (via a scoped RAII guard) and consumed by the base
        // flow, which threads them to the nested CoreAdviser/CoilAdviser as
        // LOCAL filtered views. Replaces the old DatabaseFilterScope that
        // swapped the process-shared coreDatabase/wireDatabase in place (not
        // fan-out-safe: concurrent MagneticAdviser calls saw each other's
        // filtered catalogs). Default-empty => no filtering, identical to the
        // pre-#164 base-overload behaviour.
        AdviserConstraints _constraints;
        MAS::MagneticApplication _application = MAS::MagneticApplication::POWER;
        CoreAdviser::CoreAdviserModes _coreAdviserMode = CoreAdviser::CoreAdviserModes::STANDARD_CORES;

        MagneticAdviser() {
        }
        MagneticAdviser(bool simulateResults) {
            _simulateResults = simulateResults;
        }

        void set_unique_core_shapes(bool value);
        bool get_unique_core_shapes();
        void set_application(MAS::MagneticApplication value);
        MAS::MagneticApplication get_application();
        void set_core_mode(CoreAdviser::CoreAdviserModes value);
        CoreAdviser::CoreAdviserModes get_core_mode();

        /**
         * @brief Get optimized magnetic designs using default weights.
         * @param inputs Design requirements and operating conditions.
         * @param maximumNumberResults Maximum number of designs to return.
         * @return Vector of (Mas, score) pairs sorted by descending score.
         */
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, size_t maximumNumberResults=1);

        // ctx-aware variants: apply a per-call library override and type
        // constraints (shape family, core material type, wire type) before
        // delegating to the matching base flow. ctx=nullptr means "use the
        // built-in catalogs"; an empty AdviserConstraints means "no type
        // pre-filtering" — in both cases the call is equivalent to the base
        // overload.
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs,
                                                                  size_t maximumNumberResults,
                                                                  const LibraryContext* ctx,
                                                                  const AdviserConstraints& constraints);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs,
                                                                  std::map<MagneticFilters, double> weights,
                                                                  size_t maximumNumberResults,
                                                                  const LibraryContext* ctx,
                                                                  const AdviserConstraints& constraints);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs,
                                                                  std::vector<MagneticFilterOperation> filterFlow,
                                                                  size_t maximumNumberResults,
                                                                  const LibraryContext* ctx,
                                                                  const AdviserConstraints& constraints);
        std::vector<std::pair<Mas, double>> get_advised_magnetic_fast(Inputs inputs,
                                                                       size_t maximumNumberResults,
                                                                       const LibraryContext* ctx,
                                                                       const AdviserConstraints& constraints);

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
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, const std::map<std::string, Magnetic>& catalogueMagnetics, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults=1, bool strict=true);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(std::vector<Mas> catalogueMagneticsWithInputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults=1, bool strict=true);

        /**
         * @brief Fast analytical magnetic design for Pareto exploration.
         *
         * Bypasses CoilAdviser, MagneticSimulator, and iterative loss filters.
         * Uses area product filtering, analytical gap/turns, fast_wind(), and
         * DC ohmic + Steinmetz core losses. Results are approximate but physically valid.
         *
         * @param inputs Design requirements and operating conditions.
         * @param maximumNumberResults Maximum number of designs to return.
         * @return Vector of (Mas, score) pairs sorted by ascending total losses.
         */
        std::vector<std::pair<Mas, double>> get_advised_magnetic_fast(Inputs inputs, size_t maximumNumberResults = 1);
        // Fast analytical design with a caller-supplied filter flow: strictlyRequired
        // filters (e.g. DC/EFFECTIVE_CURRENT_DENSITY) DROP any wound candidate that
        // fails them, so the fast path can gate winding current density without the
        // full CoilAdviser. Empty flow == the plain fast path above.
        std::vector<std::pair<Mas, double>> get_advised_magnetic_fast(Inputs inputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults = 1);

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

        /**
         * @brief Design magnetics from a converter topology using ngspice simulation.
         * 
         * This template method accepts any converter topology (Flyback, Buck, Boost, etc.),
         * runs ngspice simulation to extract accurate voltage/current waveforms, and returns
         * optimized magnetic designs.
         * 
         * @tparam ConverterType The converter type (e.g., Flyback, AdvancedFlyback, Buck, etc.)
         * @param converter The converter instance with design parameters
         * @param maximumNumberResults Maximum number of magnetic designs to return
         * @return Vector of (Mas, score) pairs sorted by descending score
         * 
         * @note The converter must implement simulate_and_extract_operating_points() method
         * 
         * Example usage:
         * ```cpp
         * AdvancedFlyback flyback(jsonData);
         * MagneticAdviser adviser;
         * auto results = adviser.get_advised_magnetic_from_converter(flyback, 5);
         * ```
         */
        template<typename ConverterType>
        std::vector<std::pair<Mas, double>> get_advised_magnetic_from_converter(
            ConverterType& converter,
            size_t maximumNumberResults = 1,
            bool fast = false);

        /**
         * @brief Design magnetics from a converter with custom weights.
         *
         * @tparam ConverterType The converter type
         * @param converter The converter instance
         * @param weights Map of filter types to weight values
         * @param maximumNumberResults Maximum number of designs to return
         * @param fast When true, use the fast core-only adviser
         *        (get_advised_magnetic_fast) instead of the full
         *        winding+simulation adviser. The fast adviser ignores
         *        per-filter weights.
         * @return Vector of (Mas, score) pairs sorted by descending score
         */
        template<typename ConverterType>
        std::vector<std::pair<Mas, double>> get_advised_magnetic_from_converter(
            ConverterType& converter,
            std::map<MagneticFilters, double> weights,
            size_t maximumNumberResults,
            bool fast = false);
};

// Template method implementations

/**
 * @brief Helper to run ngspice simulation for transformer-based converters (with turns ratios)
 * SFINAE overload for converters with two-parameter simulate_and_extract_operating_points
 */
template<typename ConverterType>
auto run_ngspice_simulation_helper(ConverterType& converter, 
                                   const std::vector<double>& turnsRatios, 
                                   double inductance, int)
    -> decltype(converter.simulate_and_extract_operating_points(turnsRatios, inductance)) {
    return converter.simulate_and_extract_operating_points(turnsRatios, inductance);
}

/**
 * @brief Helper to run ngspice simulation for inductor-based converters (no turns ratios)
 * SFINAE overload for converters with one-parameter simulate_and_extract_operating_points
 */
template<typename ConverterType>
auto run_ngspice_simulation_helper(ConverterType& converter, 
                                   const std::vector<double>& turnsRatios, 
                                   double inductance, long)
    -> decltype(converter.simulate_and_extract_operating_points(inductance)) {
    (void)turnsRatios;  // Unused for inductors
    return converter.simulate_and_extract_operating_points(inductance);
}

/**
 * @brief Helper for converters whose simulate_and_extract_operating_points
 * takes (double scalarTurnsRatio, double magnetizingInductance) — e.g.
 * Weinberg, Flyback, IsolatedBuck, single-output forward-derived models.
 *
 * Some converter models predate the vector<double>-per-output convention used
 * by LLC/SRC/DAB. Rather than churn every model, we route them here using the
 * first element of the turnsRatios vector. process_design_requirements
 * already emits exactly one turns ratio for single-output topologies, so
 * turnsRatios.size() == 1 is the invariant we rely on. SFINAE rank-3
 * (ellipsis ".") is lower priority than the long and int overloads above so
 * the vector and one-arg paths still win when available.
 */
template<typename ConverterType>
auto run_ngspice_simulation_helper(ConverterType& converter,
                                   const std::vector<double>& turnsRatios,
                                   double inductance, ...)
    -> decltype(converter.simulate_and_extract_operating_points(turnsRatios.front(), inductance)) {
    if (turnsRatios.empty()) {
        throw std::runtime_error(
            "run_ngspice_simulation_helper: converter expects a scalar turns "
            "ratio but design requirements emitted zero turns ratios");
    }
    return converter.simulate_and_extract_operating_points(turnsRatios.front(), inductance);
}

/**
 * @brief Run ngspice simulation with proper method dispatch
 */
template<typename ConverterType>
std::vector<OperatingPoint> run_ngspice_simulation(ConverterType& converter, 
                                                   const std::vector<double>& turnsRatios, 
                                                   double inductance) {
    return run_ngspice_simulation_helper(converter, turnsRatios, inductance, 0);
}

template<typename ConverterType>
std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic_from_converter(
    ConverterType& converter,
    size_t maximumNumberResults,
    bool fast) {
    return get_advised_magnetic_from_converter(converter, std::map<MagneticFilters, double>{}, maximumNumberResults, fast);
}

template<typename ConverterType>
std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic_from_converter(
    ConverterType& converter,
    std::map<MagneticFilters, double> weights,
    size_t maximumNumberResults,
    bool fast) {
    
    // Step 1: Get design requirements from converter (calculates optimal turns ratios, inductance, etc.)
    DesignRequirements designReqs = converter.process_design_requirements();
    
    // Step 2: Extract turns ratios and inductance from calculated design requirements
    std::vector<double> turnsRatios;
    for (const auto& turnsRatio : designReqs.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(turnsRatio));
    }
    double magnetizingInductance = resolve_dimensional_values(designReqs.get_magnetizing_inductance());
    
    // Step 3: Run ngspice simulation with calculated parameters to get accurate waveforms
    // SFINAE automatically selects the correct method based on converter type
    std::vector<OperatingPoint> operatingPoints = run_ngspice_simulation(converter, turnsRatios, magnetizingInductance);

    // Step 3b: Converters that use N identical, separately-wound per-phase
    // inductors (e.g. 3-phase Vienna's La/Lb/Lc) pack those N inductors as N
    // excitations in one operating point. They are NOT one N-winding magnetic
    // (commercial designs use three independent inductors), so design ONE of
    // them: collapse each operating point to a single winding and drop the
    // multi-winding turns-ratio list. The multiplicity is recorded on the
    // result below so downstream knows N identical parts are required.
    size_t identicalInductorCount = 1;
    if (converter.uses_identical_per_phase_inductors() && !operatingPoints.empty()) {
        identicalInductorCount = operatingPoints[0].get_excitations_per_winding().size();
        if (identicalInductorCount > 1) {
            designReqs.get_mutable_turns_ratios().clear();  // single-winding inductor
            for (auto& op : operatingPoints) {
                auto& excitations = op.get_mutable_excitations_per_winding();
                if (excitations.size() > 1) {
                    excitations.resize(1);
                }
            }
        }
    }

    // Step 4: Build Inputs object with design requirements and simulated operating points
    Inputs inputs;
    inputs.set_design_requirements(designReqs);
    inputs.set_operating_points(operatingPoints);
    
    // Step 5: Process inputs to compute harmonics, magnetizing current, etc.
    // ngspice returns raw sampled waveforms; process() computes harmonics and other derived data
    inputs.process();
    
    // Step 6: Get magnetics from adviser. The fast path uses the core-only
    // adviser (get_advised_magnetic_fast), which does not take per-filter
    // weights; weights only apply to the full adviser.
    std::vector<std::pair<Mas, double>> results;
    if (fast) {
        results = get_advised_magnetic_fast(inputs, maximumNumberResults);
    } else if (weights.empty()) {
        results = get_advised_magnetic(inputs, maximumNumberResults);
    } else {
        results = get_advised_magnetic(inputs, weights, maximumNumberResults);
    }

    // Step 7: For per-phase-inductor converters, note that N identical units
    // are required (one designed, N built).
    if (identicalInductorCount > 1) {
        const std::string note = "Design represents ONE of " + std::to_string(identicalInductorCount) +
                                 " identical inductors required (one per phase).";
        for (auto& [mas, score] : results) {
            auto& magnetic = mas.get_mutable_magnetic();
            MagneticManufacturerInfo info = magnetic.get_manufacturer_info().value_or(MagneticManufacturerInfo());
            info.set_description(note);
            magnetic.set_manufacturer_info(info);
        }
    }
    return results;
}

} // namespace OpenMagnetics