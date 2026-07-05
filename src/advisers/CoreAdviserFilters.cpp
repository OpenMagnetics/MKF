// Phase 5 (structural cleanup, nested filter classes): extracted from
// CoreAdviser.cpp. This TU owns the 12 CoreAdviser::MagneticCoreFilter*
// nested-class implementations (AreaProduct, EnergyStored, FringingFactor,
// Cost, Losses, Dimensions, TurnCount, MinimumImpedance,
// MagneticInductance, Saturation, Temperature -- 11 distinct filters,
// AreaProduct counted once).
//
// Each class is a thin adapter that delegates to the polymorphic
// MagneticFilter* family (see advisers/MagneticFilter.h) and feeds the
// resulting per-magnetic scores back through normalize_scoring +
// sort_magnetics_by_scoring (shared via advisers/CoreAdviserInternal.h).
//
// July 2026 (abt #103): the per-filter evaluate -> catch/cull -> rebuild ->
// size-check loop used to be copy-pasted verbatim 11 times. It now lives in
// one shared `evaluate_and_cull` helper (with the abt #53 per-candidate
// exception cull, now LOGGED per cull instead of silent). Each filter is a
// thin wrapper that supplies its evaluate lambda + write-back / scoring-bucket
// policy. The three real variations are preserved: write-back of the mutated
// magnetic, the TURN_COUNT custom inverted-linear scoring, and the
// Losses/Impedance/Inductance coil-compact settings restore.
//
// Sibling TU of CoreAdviser.cpp; declarations live in advisers/CoreAdviser.h.

#include "advisers/CoreAdviser.h"
#include "advisers/CoreAdviserInternal.h"
#include "advisers/MagneticFilter.h"
#include "physical_models/CoreTemperature.h"
#include "support/Logger.h"
#include "support/Settings.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace OpenMagnetics {

namespace {

using EvalFn = std::function<std::pair<bool, double>(Magnetic*, Inputs*)>;

// Shared evaluate -> cull loop (abt #103). Runs `evaluate` on each candidate; a candidate
// whose evaluation throws an OpenMagneticsException is culled (abt #53: per-candidate
// infeasibility — e.g. a degenerate gap committed for an unreachable high-Lm transformer —
// is not fatal), but the cull is now LOGGED (it used to be silent, hiding a whole pool being
// emptied by a systematic input error). Returns {survivors, their scores}; on `writeBack` the
// possibly-mutated magnetic is written back into `unfilteredMagnetics`.
std::pair<std::vector<std::pair<Magnetic, double>>, std::vector<double>>
evaluate_and_cull(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics,
                  Inputs& inputs, const EvalFn& evaluate, bool writeBack,
                  const std::string& filterName) {
    std::vector<std::pair<Magnetic, double>> filtered;
    std::vector<double> newScoring;
    std::list<size_t> listOfIndexesToErase;

    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex) {
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        bool valid = false;
        double scoring = 0.0;
        try {
            auto evalResult = evaluate(&magnetic, &inputs);
            valid = evalResult.first;
            scoring = evalResult.second;
        }
        catch (const OpenMagneticsException& e) {
            logEntry(std::string("CoreAdviser ") + filterName + ": culling infeasible candidate — " + e.what(),
                     "CoreAdviser", 2);
            valid = false;
        }

        if (valid) {
            newScoring.push_back(scoring);
            if (writeBack) {
                (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            }
        }
        else {
            listOfIndexesToErase.push_back(magneticIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filtered.push_back((*unfilteredMagnetics)[i]);
        }
    }

    if (filtered.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR,
            "Something wrong happened while filtering (" + filterName + "), size of filtered: " +
            std::to_string(filtered.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }
    return {filtered, newScoring};
}

} // anonymous namespace

// Applies the standard min-max normalized scoring for the given bucket, then sorts.
// (Member of the base filter so it can reach add_scoring / _filterConfiguration.)
void CoreAdviser::MagneticCoreFilter::apply_normalized_scoring(
        std::vector<std::pair<Magnetic, double>>* filteredMagneticsWithScoring,
        std::vector<double>& newScoring, double weight, CoreAdviser::CoreAdviserFilters bucket) {
    if (filteredMagneticsWithScoring->size() == 0) {
        return;
    }
    auto normalizedScoring = normalize_scoring(filteredMagneticsWithScoring, newScoring, weight,
                                               (*_filterConfiguration)[bucket]);
    for (size_t i = 0; i < filteredMagneticsWithScoring->size(); ++i) {
        add_scoring((*filteredMagneticsWithScoring)[i].first.get_reference(), bucket, normalizedScoring[i]);
    }
    sort_magnetics_by_scoring(filteredMagneticsWithScoring);
}

CoreAdviser::MagneticCoreFilterAreaProduct::MagneticCoreFilterAreaProduct(Inputs inputs) {
    _filter = MagneticFilterAreaProduct(inputs);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterAreaProduct::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/false, "AreaProduct");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterEnergyStored::MagneticCoreFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterEnergyStored(inputs, models);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterEnergyStored::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "EnergyStored");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterFringingFactor::MagneticCoreFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterFringingFactor(inputs, models);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterFringingFactor::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "FringingFactor");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterCost::MagneticCoreFilterCost(Inputs inputs) {
    _filter = MagneticFilterEstimatedCost(inputs);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterCost::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/false, "Cost");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::COST);
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterLosses::MagneticCoreFilterLosses(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterCoreDcAndSkinLosses(inputs, models); // F23 FIX: DC+skin effect losses
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterLosses::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    // Restore coil_delimit_and_compact even if the loss evaluation (which flips it) throws.
    auto& settings = Settings::GetInstance();
    SettingsGuard<bool> coilCompactGuard(settings, &Settings::get_coil_delimit_and_compact,
                                         &Settings::set_coil_delimit_and_compact,
                                         settings.get_coil_delimit_and_compact());
    // CAPABILITY GAP (tracked in FALLBACKS_REVIEW.md): the per-candidate loss
    // evaluators only support single-winding candidates, and the dataset stage
    // builds 1-winding dummy coils for everything except CMCs
    // (CoreAdviserDataset.cpp). For a multi-winding input every candidate
    // would be rejected by the evaluator's winding-count guard, so this filter
    // has never actually scored transformer candidates — it used to hide that
    // behind a silent all-rejected -> return-everything fallback. Keep the
    // pass-through but say so loudly until transformer-aware loss scoring
    // (secondaries populated from the requirement turns ratios) exists.
    if (!inputs.get_operating_points().empty() &&
        inputs.get_operating_points()[0].get_excitations_per_winding().size() > 1 &&
        !(*unfilteredMagnetics).empty() &&
        (*unfilteredMagnetics)[0].first.get_coil().get_functional_description().size() !=
            inputs.get_operating_points()[0].get_excitations_per_winding().size()) {
        logEntry("Losses filter cannot score multi-winding inputs against single-winding candidate coils; passing " +
                     std::to_string((*unfilteredMagnetics).size()) + " candidates through unscored",
                 "CoreAdviser");
        return *unfilteredMagnetics;
    }

    // NOTE: previously, when the loss model rejected EVERY candidate, this filter
    // silently returned the whole set unscored. That routed cores with
    // uncomputable losses into the ranking; now it returns empty like every
    // other filter (the pipeline's retry handles genuine zero-result cases).
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "Losses");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterDimensions::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/false, "Dimensions");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::DIMENSIONS);
    return filteredMagneticsWithScoring;
}

// Mirrors MagneticCoreFilterDimensions but routes scoring through the
// TURN_COUNT bucket and reads the score from MagneticFilterTurnCount
// (N_total × characteristic dimension). Used by the suppression pipelines
// to penalise low-µ candidates that need an absurd turn count to satisfy
// the inductance / impedance requirement (e.g. powder cores attempting CMC).
std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterTurnCount::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/false, "TurnCount");

    if (filteredMagneticsWithScoring.size() > 0) {
        // Compute inverted linear score directly to avoid the all-equal fallback (+1 ignoring
        // weight) in OpenMagnetics::normalize_scoring.  Fewer turns → higher score → ranked first.
        // Raw score = totalTurns (higher is worse); we invert explicitly: contribution = weight*(1 - (N-min)/(max-min)).
        double minTurns = *std::min_element(newScoring.begin(), newScoring.end());
        double maxTurns = *std::max_element(newScoring.begin(), newScoring.end());
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            double contrib;
            if (maxTurns > minTurns) {
                contrib = weight * (1.0 - (newScoring[i] - minTurns) / (maxTurns - minTurns));
            } else {
                // All same turn count — contribute weight/2 so the filter is neutral
                contrib = weight * 0.5;
            }
            filteredMagneticsWithScoring[i].second += contrib;
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::TURN_COUNT, contrib);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterMinimumImpedance::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto& settings = Settings::GetInstance();
    SettingsGuard<bool> coilCompactGuard(settings, &Settings::get_coil_delimit_and_compact,
                                         &Settings::set_coil_delimit_and_compact,
                                         settings.get_coil_delimit_and_compact());
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "MinimumImpedance");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterMagneticInductance::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto& settings = Settings::GetInstance();
    SettingsGuard<bool> coilCompactGuard(settings, &Settings::get_coil_delimit_and_compact,
                                         &Settings::set_coil_delimit_and_compact,
                                         settings.get_coil_delimit_and_compact());
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "MagneticInductance");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterSaturation::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "Saturation");
    // Saturation only contributes a headroom SCORE when it is not the first (gating) pass.
    if (!firstFilter) {
        apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterTemperature::MagneticCoreFilterTemperature(
    Inputs inputs, std::map<std::string, std::string> models, double maximumTemperature)
{
    _filter = MagneticFilterTemperature(inputs, maximumTemperature);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterTemperature::filter_magnetics(
    std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter)
{
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "Temperature");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    return filteredMagneticsWithScoring;
}

} // namespace OpenMagnetics
