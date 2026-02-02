#pragma once
#include "advisers/MagneticFilter.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Impedance.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "constructive_models/Coil.h"
#include "processors/MagneticSimulator.h"
#include "processors/Outputs.h"
#include "processors/Inputs.h"
#include "constructive_models/Core.h"
#include "constructive_models/Mas.h"
#include <cmath>
#include <MAS.hpp>
#include "support/Exceptions.h"

using namespace MAS;

namespace OpenMagnetics {

/**
 * @class CoreAdviser
 * @brief Multi-criteria magnetic core recommendation system.
 *
 * ## Overview
 * CoreAdviser selects optimal magnetic cores for power electronics applications based on
 * user-defined priorities. It evaluates cores from a database using multiple criteria and
 * returns ranked recommendations with normalized scores.
 *
 * ## Scoring System
 * The adviser uses three main filter categories, each with configurable weights (0.0-1.0):
 * - **COST**: Estimated manufacturing/purchasing cost (lower is better)
 * - **EFFICIENCY**: Power losses (core + DC winding losses, lower is better)
 * - **DIMENSIONS**: Physical size/volume (smaller is better)
 *
 * ### Score Calculation
 * 1. Each filter computes a raw score for each core candidate
 * 2. Scores are normalized using `normalize_scoring()`:
 *    - `invert=true`: Lower raw values get higher scores (used for cost, losses, size)
 *    - `log=true`: Logarithmic normalization (compresses large differences)
 *    - `log=false`: Linear normalization (preserves proportional differences)
 * 3. Final score = Σ(normalized_score × weight) for all filters
 *
 * ## Filter Pipeline
 *
 * ### Power Application (filter_available_cores_power_application)
 * The filter chain for power inductors/transformers:
 * 1. **filterAreaProduct**: Pre-filter using Area Product (AP = Aw × Ac)
 *    - Fixed weight 1.0 (binary pass/fail, not scored)
 *    - Eliminates cores too small for required energy storage
 * 2. **filterEnergyStored**: Validates energy storage capacity (L×I²/2)
 *    - Fixed weight 1.0 (binary pass/fail, not scored)
 *    - Checks if core can store required energy without saturation
 * 3. **filterCost**: Scores by estimated cost
 *    - Weight: COST user weight
 * 4. **filterDimensions**: Scores by physical volume
 *    - Weight: DIMENSIONS user weight
 *    - Uses linear normalization to preserve size differences
 * 5. **filterLosses**: Scores by total losses (core + winding DC)
 *    - Weight: EFFICIENCY user weight
 *
 * ### Suppression Application (filter_available_cores_suppression_application)
 * The filter chain for EMI/RFI suppression:
 * 1. **filterMinimumImpedance**: Pre-filter by impedance at operating frequency
 * 2. **filterCost**: Scores by estimated cost
 * 3. **filterDimensions**: Scores by physical volume
 * 4. **filterMagneticInductance**: Scores by magnetizing inductance
 * 5. **filterLosses**: Scores by total losses
 *
 * ## Operating Modes
 * - **AVAILABLE_CORES**: Uses manufacturer stock database (fastest)
 * - **STANDARD_CORES**: Uses standard core shapes with material optimization
 * - **CUSTOM_CORES**: Uses user-provided core list
 *
 * ## Toroid Handling
 * Toroidal cores use geometry-based bobbin filling factor calculation:
 *   fillingFactor = 0.55 + 0.15 × (innerRadius / outerRadius)
 * This accounts for the reduced winding area in toroid centers (range: 0.55-0.70).
 *
 * ## Usage Example
 * ```cpp
 * std::map<CoreAdviserFilters, double> weights = {
 *     {CoreAdviserFilters::COST, 0.3},
 *     {CoreAdviserFilters::EFFICIENCY, 0.5},
 *     {CoreAdviserFilters::DIMENSIONS, 0.2}
 * };
 * CoreAdviser adviser;
 * auto results = adviser.get_advised_core(inputs, weights, 5);  // Top 5 cores
 * ```
 *
 * ## References
 * - Industry practice: LI² (energy storage), Area Product method
 * - Colonel Wm. T. McLyman, "Transformer and Inductor Design Handbook"
 */
class CoreAdviser {
    public: 
        enum class CoreAdviserFilters : int {
            COST,
            EFFICIENCY,
            DIMENSIONS
        };

        enum class CoreAdviserModes : int {
            AVAILABLE_CORES,
            STANDARD_CORES,
            CUSTOM_CORES
        };
    protected:
        std::map<std::string, std::string> _models;
        bool _uniqueCoreShapes = false;
        std::map<CoreAdviserFilters, double> _weights;
        MagneticSimulator _magneticSimulator;
        WindingOhmicLosses _windingOhmicLosses;
        Application _application = Application::POWER;
        CoreAdviserModes _mode = CoreAdviserModes::AVAILABLE_CORES;


    public:

        /**
         * @brief Filter normalization configuration for each scoring category.
         *
         * Each filter category has two configuration options:
         * - **invert**: If true, lower raw values result in higher normalized scores.
         *   Used when "less is better" (cost, losses, size).
         * - **log**: If true, uses logarithmic normalization; if false, linear.
         *   Linear normalization preserves proportional differences between cores.
         *
         * Configuration rationale:
         * - COST: inverted (cheaper=better), logarithmic (diminishing returns on savings)
         * - EFFICIENCY: inverted (lower losses=better), logarithmic (diminishing returns)
         * - DIMENSIONS: inverted (smaller=better), LINEAR (size differences should be
         *   proportionally reflected; a core 2× larger should score proportionally worse)
         */
        std::map<CoreAdviserFilters, std::map<std::string, bool>> _filterConfiguration{
                { CoreAdviserFilters::COST,                  { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::EFFICIENCY,            { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::DIMENSIONS,            { {"invert", true}, {"log", false} } }
            };
        std::map<CoreAdviserFilters, std::map<std::string, double>> _scorings;
        CoreAdviser(std::map<std::string, std::string> models) {
            auto defaults = OpenMagnetics::Defaults();
            _models = models;
            if (models.find("gapReluctance") == models.end()) {
                _models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
            }
        }
        CoreAdviser() {
            _models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
            _models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
            _models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);
        }
        /**
         * @brief Get per-core scoring breakdown for debugging/analysis.
         * @param weighted If true, returns weighted scores; otherwise raw normalized scores.
         * @return Map of core names to their per-filter scores.
         */
        std::map<std::string, std::map<CoreAdviserFilters, double>> get_scorings(bool weighted = false);

        void set_unique_core_shapes(bool value);
        bool get_unique_core_shapes();
        void set_application(Application value);
        Application get_application();
        void set_mode(CoreAdviserModes value);
        CoreAdviserModes get_mode();

        /**
         * @brief Main entry point for core recommendation.
         * @param inputs Operating conditions (voltage, current, frequency, etc.)
         * @param maximumNumberResults Maximum number of cores to return.
         * @return Vector of (Mas, score) pairs, sorted by descending score.
         */
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, size_t maximumNumberResults=1);
        
        /**
         * @brief Core recommendation with custom weights.
         * @param inputs Operating conditions.
         * @param weights Map of filter weights (COST, EFFICIENCY, DIMENSIONS). Values 0.0-1.0.
         * @param maximumNumberResults Maximum number of cores to return.
         * @return Vector of (Mas, score) pairs, sorted by descending score.
         */
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<CoreShape>* shapes, size_t maximumNumberResults=1);

        Mas post_process_core(Magnetic magnetic, Inputs inputs);
        
        /**
         * @brief Filter pipeline for power inductor/transformer applications.
         *
         * Applies filters in order: AreaProduct → EnergyStored → Cost → Dimensions → Losses.
         * - AreaProduct and EnergyStored use fixed weight 1.0 (pre-filtering, not scoring)
         * - Cost, Dimensions, Losses use user-provided weights for scoring
         *
         * @param magnetics Input list of candidate magnetics with initial scores.
         * @param inputs Operating conditions.
         * @param weights User-defined weights for COST, EFFICIENCY, DIMENSIONS.
         * @param maximumMagneticsAfterFiltering Max candidates to keep after filtering.
         * @param maximumNumberResults Max results to return.
         * @return Filtered and scored list of (Mas, score) pairs.
         */
        std::vector<std::pair<Mas, double>> filter_available_cores_power_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        
        /**
         * @brief Filter pipeline for EMI/RFI suppression applications.
         *
         * Applies filters: MinimumImpedance → Cost → Dimensions → MagneticInductance → Losses.
         *
         * @param magnetics Input list of candidate magnetics with initial scores.
         * @param inputs Operating conditions (includes required impedance at frequency).
         * @param weights User-defined weights for COST, EFFICIENCY, DIMENSIONS.
         * @param maximumMagneticsAfterFiltering Max candidates to keep after filtering.
         * @param maximumNumberResults Max results to return.
         * @return Filtered and scored list of (Mas, score) pairs.
         */
        std::vector<std::pair<Mas, double>> filter_available_cores_suppression_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        
        /**
         * @brief Filter pipeline for standard core shapes (with material selection).
         *
         * Similar to filter_available_cores_power_application but works with
         * parametric core shapes rather than specific stock cores.
         */
        std::vector<std::pair<Mas, double>> filter_standard_cores_power_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        
        /**
         * @brief Filter pipeline for standard core shapes in suppression applications.
         */
        std::vector<std::pair<Mas, double>> filter_standard_cores_interference_suppression_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        std::vector<std::pair<Magnetic, double>> create_magnetic_dataset(Inputs inputs, std::vector<Core>* cores, bool includeStacks);
        std::vector<std::pair<Magnetic, double>> create_magnetic_dataset(Inputs inputs, std::vector<CoreShape>* shapes, bool includeStacks);
        void expand_magnetic_dataset_with_stacks(Inputs inputs, std::vector<Core>* cores, std::vector<std::pair<Magnetic, double>>* magnetics);
        std::vector<std::pair<Magnetic, double>> add_powder_materials(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs);
        std::vector<std::pair<Magnetic, double>> add_ferrite_materials_by_losses(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs);
        std::vector<std::pair<Magnetic, double>> add_ferrite_materials_by_impedance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs);
        bool should_include_powder(Inputs inputs);
    
    class MagneticCoreFilter {
        public:
            std::map<CoreAdviserFilters, std::map<std::string, double>>* _scorings;
            std::map<CoreAdviserFilters, std::map<std::string, bool>>* _filterConfiguration;
            bool _useCache = true;

            void add_scoring(std::string name, CoreAdviser::CoreAdviserFilters filter, double scoring) {
                if (scoring != -1) {
                    if ((*_scorings)[filter].find(name) == (*_scorings)[filter].end()) {
                        (*_scorings)[filter][name] = 0;
                    }
                    auto oldScoring = (*_scorings)[filter][name];
                    (*_scorings)[filter][name] = oldScoring + scoring;
                }
            }
            void set_scorings(std::map<CoreAdviserFilters, std::map<std::string, double>>* scorings) {
                _scorings = scorings;
            }
            void set_filter_configuration(std::map<CoreAdviserFilters, std::map<std::string, bool>>* filterConfiguration) {
                _filterConfiguration = filterConfiguration;
            }
            void set_cache_usage(bool status) {
                _useCache = status;
            }
            MagneticCoreFilter(){
            }
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, Inputs inputs, double weight=1);
    };
    
    class MagneticCoreFilterAreaProduct : public MagneticCoreFilter {
        private:
            MagneticFilterAreaProduct _filter;

        public:
            MagneticCoreFilterAreaProduct(Inputs inputs);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterEnergyStored : public MagneticCoreFilter {
        private:
            MagneticFilterEnergyStored _filter;

        public:
            MagneticCoreFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterFringingFactor : public MagneticCoreFilter {
        private:
            MagneticFilterFringingFactor _filter;

        public:
            MagneticCoreFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterCost : public MagneticCoreFilter {
        private:
            MagneticFilterEstimatedCost _filter;

        public:
            MagneticCoreFilterCost(Inputs inputs);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterLosses : public MagneticCoreFilter {
        private:
            MagneticFilterCoreAndDcLosses _filter;

        public:
            MagneticCoreFilterLosses(Inputs inputs, std::map<std::string, std::string> models);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterDimensions : public MagneticCoreFilter {
        private:
            MagneticFilterDimensions _filter;
        public:
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterMinimumImpedance : public MagneticCoreFilter {
        private:
            MagneticFilterCoreMinimumImpedance _filter;
        public:
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterMagneticInductance : public MagneticCoreFilter {
        private:
            MagneticFilterMagnetizingInductance _filter;
        public:
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };

};

void from_json(const json & j, CoreAdviser::CoreAdviserModes & x);
void to_json(json & j, const CoreAdviser::CoreAdviserModes & x);

inline void from_json(const json & j, CoreAdviser::CoreAdviserModes & x) {
    if (j == "available cores") x = CoreAdviser::CoreAdviserModes::AVAILABLE_CORES;
    else if (j == "standard cores") x = CoreAdviser::CoreAdviserModes::STANDARD_CORES;
    else if (j == "custom cores") x = CoreAdviser::CoreAdviserModes::CUSTOM_CORES;
    else { throw std::runtime_error("Input JSON does not conform to schema!"); }
}

inline void to_json(json & j, const CoreAdviser::CoreAdviserModes & x) {
    switch (x) {
        case CoreAdviser::CoreAdviserModes::AVAILABLE_CORES: j = "available cores"; break;
        case CoreAdviser::CoreAdviserModes::STANDARD_CORES: j = "standard cores"; break;
        case CoreAdviser::CoreAdviserModes::CUSTOM_CORES: j = "custom cores"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"CoreAdviser::CoreAdviserModes\": " + std::to_string(static_cast<int>(x)));
    }
}

void from_json(const json & j, CoreAdviser::CoreAdviserFilters & x);
void to_json(json & j, const CoreAdviser::CoreAdviserFilters & x);

inline void from_json(const json & j, CoreAdviser::CoreAdviserFilters & x) {
    if (j == "Cost" || j == "cost" || j == "COST") x = CoreAdviser::CoreAdviserFilters::COST;
    else if (j == "Efficiency" || j == "efficiency" || j == "EFFICIENCY") x = CoreAdviser::CoreAdviserFilters::EFFICIENCY;
    else if (j == "Dimensions" || j == "dimensions" || j == "DIMENSIONS") x = CoreAdviser::CoreAdviserFilters::DIMENSIONS;
    else { throw std::runtime_error("Input JSON does not conform to schema!"); }
}

inline void to_json(json & j, const CoreAdviser::CoreAdviserFilters & x) {
    switch (x) {
        case CoreAdviser::CoreAdviserFilters::COST: j = "Cost"; break;
        case CoreAdviser::CoreAdviserFilters::EFFICIENCY: j = "Efficiency"; break;
        case CoreAdviser::CoreAdviserFilters::DIMENSIONS: j = "Dimensions"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"CoreAdviser::CoreAdviserFilters\": " + std::to_string(static_cast<int>(x)));
    }
}

} // namespace OpenMagnetics