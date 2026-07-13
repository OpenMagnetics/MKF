#include "advisers/CoreAdviser.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"  // for is_energy_storing_topology()
#include "advisers/CoreAdviserInternal.h"      // log_stage/log_pruned, make_sinusoidal_excitation, CoreLossesModelPair
#include "constructive_models/Bobbin.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/CoilMesher.h"
#include "constructive_models/NumberTurns.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Bobbin.h"
#include <algorithm>
#include <cctype>
#include <set>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <execution>
#include <magic_enum_utility.hpp>
#include <list>
#include <cmrc/cmrc.hpp>
#include "support/Exceptions.h"
#include "support/Logger.h"
#include "support/Settings.h"

CMRC_DECLARE(data);


namespace OpenMagnetics {

namespace {
// ABT #153: marker appended to a result's manufacturer reference and core name
// when the synthetic STANDARD_CORES catalogue produced nothing and the adviser
// fell through to the AVAILABLE_CORES manufacturer catalogue. The whole point of
// the STANDARD_CORES catalogue is a broad synthesised set, yet it only pairs
// ferrite materials with concentric (E/U/PQ…) shapes and powder materials with
// toroids — so an off-the-shelf powder concentric core (e.g. a Kool Mµ E-core,
// the natural pick for a high-energy line-frequency PFC boost inductor) is never
// synthesised, even though AVAILABLE_CORES lists exactly such parts. Rather than
// returning "no core can be advised", fall through to the real manufacturer
// catalogue and TAG every result so the mode substitution is explicit in the
// returned MAS — never a silent fallback (house rule).
const std::string kAvailableFallbackTag =
    " [advised from available-cores catalogue: no standard core met this spec]";

void tag_results_as_available_fallback(std::vector<std::pair<Mas, double>>& results) {
    for (auto& [mas, scoring] : results) {
        auto& magnetic = mas.get_mutable_magnetic();
        if (magnetic.get_manufacturer_info()) {
            auto info = magnetic.get_manufacturer_info().value();
            info.set_reference(info.get_reference().value_or("unnamed") + kAvailableFallbackTag);
            magnetic.set_manufacturer_info(info);
        }
        auto& core = magnetic.get_mutable_core();
        core.set_name(core.get_name().value_or("unnamed") + kAvailableFallbackTag);
    }
}
}  // namespace


void CoreAdviser::set_unique_core_shapes(bool value) {
    _uniqueCoreShapes = value;
}

bool CoreAdviser::get_unique_core_shapes() {
    return _uniqueCoreShapes;
}

void CoreAdviser::set_application(MAS::MagneticApplication value) {
    _application = value;
}

MAS::MagneticApplication CoreAdviser::get_application() {
    return _application;
}

void CoreAdviser::set_mode(CoreAdviser::CoreAdviserModes value) {
    _mode = value;
}

CoreAdviser::CoreAdviserModes CoreAdviser::get_mode() {
    return _mode;
}

void CoreAdviser::set_weights(std::map<CoreAdviserFilters, double> weights) {
    _weights = weights;
}

std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> CoreAdviser::get_scorings(bool weighted){
    std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> swappedScorings;
    for (auto& [filter, aux] : _scorings) {
        auto filterConfiguration = _filterConfiguration[filter];

        auto normalizedScorings = OpenMagnetics::normalize_scoring(aux, weighted? _weights[filter] : 1, false, false);

        for (auto& [name, scoring] : aux) {
            swappedScorings[name][filter] = normalizedScorings[name];
        }
    }
    return swappedScorings;
}

std::vector<double> normalize_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, std::vector<double> newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    auto normalizedScorings = OpenMagnetics::normalize_scoring(newScoring, weight, filterConfiguration);

    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
        (*magneticsWithScoring)[i].second += normalizedScorings[i];
    }

    // stable_sort((*magneticsWithScoring).begin(), (*magneticsWithScoring).end(), [](std::pair<Magnetic, double>& b1, std::pair<Magnetic, double>& b2) {
    //     return b1.second > b2.second;
    // });

    return normalizedScorings;
}

void sort_magnetics_by_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring) {
    stable_sort((*magneticsWithScoring).begin(), (*magneticsWithScoring).end(), [](const std::pair<Magnetic, double>& b1, const std::pair<Magnetic, double>& b2) {
        return b1.second > b2.second;
    }); // F12 FIX: stable_sort for reproducible results
}



std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(
        Inputs inputs,
        std::map<CoreAdviserFilters, double> weights,
        size_t maximumNumberResults,
        const LibraryContext* ctx,
        const AdviserConstraints& constraints) {
    // Apply per-call library override (RAII swap on the global *Database maps).
    auto scope = ctx ? ctx->applyScoped() : LibraryContext::Scope{};

    if (constraints.empty()) {
        return get_advised_core(inputs, weights, maximumNumberResults);
    }

    // Type pre-filter at the top of the pipeline. Only AVAILABLE_CORES mode
    // operates on explicit Core instances; STANDARD_CORES / CUSTOM_CORES
    // expand from shape lists, so the shape-family filter is applied by
    // filtering coreShapeDatabase indirectly via constraints on the cores
    // dataset later in the flow (see filter inside the shapes branch).
    if (get_mode() == CoreAdviserModes::AVAILABLE_CORES) {
        if (coreDatabase.empty() && !LibraryContext::Scope::anyActive()) {
            load_cores();
        }
        auto filtered = filterCoresByConstraints(coreDatabase, constraints);
        return get_advised_core(inputs, weights, &filtered, maximumNumberResults);
    }

    if (get_mode() == CoreAdviserModes::STANDARD_CORES) {
        if (coreMaterialDatabase.empty()) load_core_materials();
        if (coreShapeDatabase.empty())    load_core_shapes();

        std::vector<MAS::CoreShape> shapes;
        std::set<std::string> seen;
        for (auto& [key, shape] : coreShapeDatabase) {
            if (!constraints.shapeFamily.empty()
                && !acceptsCoreShapeFamily(constraints.shapeFamily, shape.get_family())) {
                continue;
            }
            auto canonical = shape.get_name().value_or("");
            if (canonical.empty() || !seen.insert(canonical).second) continue;
            shapes.push_back(shape);
        }
        // Material-type constraint is applied later in the pipeline by
        // add_powder_materials / add_ferrite_materials_by_losses; for now we
        // pass the trimmed shape set in and rely on the existing material
        // selection logic. (Future refinement: thread coreMaterialType into
        // those helpers explicitly.)
        _weights = weights; // the shapes overload consumes _weights
        return get_advised_core(inputs, &shapes, maximumNumberResults);
    }

    // CUSTOM_CORES and any other mode: fall back to default flow (constraints
    // do not currently affect synthetic shape generation).
    return get_advised_core(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults){
    if (get_mode() == CoreAdviserModes::AVAILABLE_CORES) {
        if (coreDatabase.empty() && !LibraryContext::Scope::anyActive()) {
            load_cores();
        }
        return get_advised_core(inputs, weights, &coreDatabase, maximumNumberResults);
    }
    else if (get_mode() == CoreAdviserModes::STANDARD_CORES) {
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        if (coreShapeDatabase.empty()) {
            load_core_shapes();
        }
        std::vector<MAS::CoreShape> shapes;
        // Phase 1 fix: coreShapeDatabase is keyed by canonical-name AND by
        // each alias (see Utils.cpp:255-260), so a naive iteration yields the
        // same CoreShape (1 + #aliases) times. That was the root cause of
        // CoreAdviser STANDARD_CORES × POWER returning duplicate top-N
        // entries with identical names and scores (e.g. three copies of
        // "95 PQ 27/15 gapped 0.36 mm" in slots 0-2). Dedupe by canonical
        // shape.get_name() before feeding the dataset builder.
        std::set<std::string> seenShapeNames;
        for (auto& [key, shape] : coreShapeDatabase) {
            auto canonical = shape.get_name().value_or("");
            if (canonical.empty()) {
                throw std::runtime_error("CoreShape in coreShapeDatabase has no name");
            }
            if (!seenShapeNames.insert(canonical).second) {
                continue;
            }
            shapes.push_back(shape);
        }
        _weights = weights; // the shapes overload consumes _weights
        auto standardResults = get_advised_core(inputs, &shapes, maximumNumberResults);
        // ABT #153: if the synthetic standard-core catalogue served the design,
        // return it. Otherwise (a POWER design the standard synthesis cannot
        // build — see kAvailableFallbackTag) fall through to the AVAILABLE_CORES
        // manufacturer catalogue and tag the results loudly instead of returning
        // an empty set ("No core can be advised. You are on your own.").
        if (!standardResults.empty() || get_application() != MAS::MagneticApplication::POWER) {
            return standardResults;
        }
        logEntry("STANDARD_CORES synthesis returned no candidates for this power design; "
                 "falling through to the AVAILABLE_CORES manufacturer catalogue "
                 "(results are tagged to make the mode substitution explicit).", "CoreAdviser");
        if (coreDatabase.empty() && !LibraryContext::Scope::anyActive()) {
            load_cores();
        }
        auto fallbackResults = get_advised_core(inputs, weights, &coreDatabase, maximumNumberResults);
        tag_results_as_available_fallback(fallbackResults);
        return fallbackResults;
    }
    else if (get_mode() == CoreAdviserModes::CUSTOM_CORES) {
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        if (coreShapeDatabase.empty()) {
            load_core_shapes();
        }
        auto customShapes = create_custom_core_shapes(inputs);
        _weights = weights; // the shapes overload consumes _weights
        return get_advised_core(inputs, &customShapes, maximumNumberResults);
    }
    else {
        throw NotImplementedException("Unknown CoreAdviserMode");
    }
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults, maximumNumberCores);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {

    std::vector<std::pair<Mas, double>> results;


    for(size_t i = 0; i < (*cores).size(); i += maximumNumberCores) {
        auto last = std::min((*cores).size(), i + maximumNumberCores);
        std::vector<Core> partialCores;
        partialCores.reserve(last - i);

        // ABT #164: COPY (not move) the batch out of the caller's vector. The
        // batching entry point is fed &coreDatabase (the process-shared core
        // catalog) by get_advised_core(...) above, so moving elements out gutted
        // the shared database for every subsequent caller — even a sequential
        // same-thread reuse — and silently corrupted the parallel fan-out. The
        // extra memory is bounded to ONE batch (<= maximumNumberCores cores),
        // transient per iteration, not a full-database copy.
        std::copy((*cores).begin() + i, (*cores).begin() + last, std::back_inserter(partialCores));

        auto partialResult = get_advised_core(inputs, weights, &partialCores, maximumNumberResults);
        std::move(partialResult.begin(), partialResult.end(), std::back_inserter(results));
    }

    stable_sort(results.begin(), results.end(), [](const std::pair<Mas, double>& b1, const std::pair<Mas, double>& b2) {
        return b1.second > b2.second;
    });

    if (results.size() > maximumNumberResults) {
        results = std::vector<std::pair<Mas, double>>(results.begin(), results.end() - (results.size() - maximumNumberResults));
    }
    return results;
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults){
    _weights = weights;
 
    size_t maximumMagneticsAfterFiltering = settings.get_core_adviser_maximum_magnetics_after_filtering();
    std::vector<std::pair<Magnetic, double>> magnetics;

    magnetics = create_magnetic_dataset(inputs, cores, false);

    // logEntry("We start the search with " + std::to_string(magnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.", "CoreAdviser");
    // logEntry("We don't include stacks of cores in our search.", "CoreAdviser");
    std::vector<std::pair<Mas, double>> filteredMagnetics;
    
    if (get_application() == MAS::MagneticApplication::POWER) {
        filteredMagnetics = filter_available_cores_power_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMagnetics.size() >= maximumNumberResults) {
            return filteredMagnetics;
        }

        auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
        if (globalIncludeStacks) {
            expand_magnetic_dataset_with_stacks(inputs, cores, &magnetics);
        }

        logEntry("First attempt produced not enough results, so now we are searching again with " + std::to_string(magnetics.size()) + " magnetics, including up to " + std::to_string(defaults.coreAdviserMaximumNumberStacks) + " cores stacked when possible.", "CoreAdviser");
        maximumMagneticsAfterFiltering = magnetics.size();
        filteredMagnetics = filter_available_cores_power_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        return filteredMagnetics;
    }
    else {
        filteredMagnetics = filter_available_cores_suppression_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        return filteredMagnetics;
    }
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<CoreShape>* shapes, size_t maximumNumberResults) {
    auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
    auto magnetics = create_magnetic_dataset(inputs, shapes, globalIncludeStacks);

    size_t maximumMagneticsAfterFiltering = settings.get_core_adviser_maximum_magnetics_after_filtering();
    if (get_application() == MAS::MagneticApplication::POWER) {
        return filter_standard_cores_power_application(&magnetics, inputs, _weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
    else {
        return filter_standard_cores_interference_suppression_application(&magnetics, inputs, _weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
}

}