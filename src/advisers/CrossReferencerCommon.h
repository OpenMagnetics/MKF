#pragma once
// Phase 3 (F2+F3): shared building blocks for cross-referencer classes.
//
// CoreCrossReferencer and CoreMaterialCrossReferencer have the same overall
// shape (gather candidates → apply a sequence of filters → normalise +
// rank scores → return top-N). Most of their bodies are genuinely
// shape-specific (different filter enums, filter-method signatures, item
// types) and don't compress well under a template. But the score-
// normalisation pass and the default-models initialisation are
// byte-identical aside from the templated item type, and were drifting
// independently (both files carry hand-applied F5/B8/XC-6/F12 fix
// comments). Centralising them here removes the duplicate-maintenance
// hazard.

#include "Defaults.h"
#include "support/Utils.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace OpenMagnetics {

// Populate the three model-name entries that every cross-referencer
// needs, if they are not already set by the caller. Idempotent.
inline void init_default_cross_referencer_models(std::map<std::string, std::string>& models) {
    auto defaults = OpenMagnetics::Defaults();
    if (models.find("gapReluctance") == models.end()) {
        models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
    }
    if (models.find("coreLosses") == models.end()) {
        models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
    }
    if (models.find("coreTemperature") == models.end()) {
        models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);
    }
}

// Update the cumulative ranking score of each item using its raw filter
// scoring. The score is rescaled to [0, 1] across the population (linear
// or log, inverted or direct, per filterConfiguration), then weighted
// and added onto the running total. Items with NaN/Inf scores are
// projected onto the population maximum (or the all-equal/all-NaN
// short-circuit is taken).
//
// Templated on the item type (Core, CoreMaterial). Behaviour and the
// embedded fix history (F5/B8/XC-6/F12) are preserved verbatim from the
// two prior independent copies.
template <typename Item>
void normalize_scoring(std::vector<std::pair<Item, double>>* rankedItems,
                       std::vector<double>* newScoring,
                       double weight,
                       std::map<std::string, bool> filterConfiguration) {
    double maximumScoring = *std::max_element(newScoring->begin(), newScoring->end());
    double minimumScoring = *std::min_element(newScoring->begin(), newScoring->end());

    // F5 FIX: Guard against all-NaN/Inf scores
    if (std::isnan(maximumScoring) || std::isinf(maximumScoring)) {
        for (size_t i = 0; i < rankedItems->size(); ++i) {
            (*rankedItems)[i].second += weight;
        }
        std::stable_sort(rankedItems->begin(), rankedItems->end(),
            [](const std::pair<Item, double>& a, const std::pair<Item, double>& b) { return a.second > b.second; });
        return;
    }

    double dataRelativeFloor = std::max(defaults.crossReferencerScoringAbsoluteFloor, maximumScoring * defaults.crossReferencerScoringDataRelativeFloorRatio); // B8 FIX: data-relative floor

    // O24 FIX: NaN/Inf protection
    for (size_t idx = 0; idx < newScoring->size(); ++idx) {
        if (std::isnan((*newScoring)[idx]) || std::isinf((*newScoring)[idx])) {
            (*newScoring)[idx] = maximumScoring;
        }
    }

    for (size_t i = 0; i < (*rankedItems).size(); ++i) {
        auto scoring = (*newScoring)[i];
        if (std::isnan(scoring)) {
            scoring = maximumScoring;
        }
        else {
            scoring = std::max(dataRelativeFloor, scoring); // B8 FIX
        }
        minimumScoring = std::max(dataRelativeFloor, minimumScoring); // B8 FIX
        if (maximumScoring != minimumScoring) {
            if (filterConfiguration["log"]) {
                if (filterConfiguration["invert"]) {
                    (*rankedItems)[i].second = (*rankedItems)[i].second + weight * (1 - (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring)));
                }
                else {
                    (*rankedItems)[i].second = (*rankedItems)[i].second + weight * (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                }
            }
            else {
                if (filterConfiguration["invert"]) {
                    (*rankedItems)[i].second = (*rankedItems)[i].second + weight * (1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring));
                }
                else {
                    (*rankedItems)[i].second = (*rankedItems)[i].second + weight * (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                }
            }
        }
        else {
            (*rankedItems)[i].second = (*rankedItems)[i].second + weight * defaults.crossReferencerNeutralScoreWhenEqual; // XC-6 FIX: neutral score when all equal (was full weight)
        }
    }
    std::stable_sort((*rankedItems).begin(), (*rankedItems).end(),
        [](const std::pair<Item, double>& b1, const std::pair<Item, double>& b2) {
            return b1.second > b2.second;
        }); // F12 FIX: stable_sort for reproducible results
}

// Build the per-name, per-filter [0, 1]-normalised scoring table that
// CoreCrossReferencer::get_scorings and CoreMaterialCrossReferencer::
// get_scorings each produced via near-identical inline code (~65 LOC
// each, drifted in trivial ways — CMCR used a hardcoded 0.5 where CCR
// used the named constant, CMCR had NaN guards CCR didn't, CCR had
// the data-relative floor CMCR didn't).
//
// Inputs: the cross-referencer's _scorings map (per filter, per item),
// the per-filter configuration (log / invert flags), the per-filter
// weights, and whether the caller wants weighted output.
//
// Behaviour preserved (union of both prior implementations):
//   * Data-relative floor on minimumScoring (was CCR-only, B8/F4 FIX).
//   * Neutral score from defaults.crossReferencerNeutralScoreWhenEqual
//     (was CCR-only, XC-6 FIX) — replaces CMCR's hardcoded 0.5.
//   * NaN check on input scoring + output (was CMCR-only).
//   * Strict min==max check from CMCR (more conservative than CCR's
//     std::abs < absoluteFloor; CCR's variant would treat values
//     differing by less than 1e-10 as "equal" and assign neutral).
//
// Templated on the cross-referencer's FilterEnum (CoreCrossReferencer-
// Filters or CoreMaterialCrossReferencerFilters).
template <typename FilterEnum>
std::map<std::string, std::map<FilterEnum, double>> compute_normalized_scorings(
    const std::map<FilterEnum, std::map<std::string, double>>& scorings,
    std::map<FilterEnum, std::map<std::string, bool>>& filterConfiguration,
    std::map<FilterEnum, double>& weights,
    bool weighted)
{
    auto defaults = OpenMagnetics::Defaults();
    std::map<std::string, std::map<FilterEnum, double>> swappedScorings;

    for (const auto& [filter, aux] : scorings) {
        if (aux.empty()) continue;
        auto& fc = filterConfiguration[filter];

        double maximumScoring = (*std::max_element(aux.begin(), aux.end(),
                                     [](const std::pair<std::string, double>& p1,
                                        const std::pair<std::string, double>& p2)
                                     { return p1.second < p2.second; })).second;
        double minimumScoring = (*std::min_element(aux.begin(), aux.end(),
                                     [](const std::pair<std::string, double>& p1,
                                        const std::pair<std::string, double>& p2)
                                     { return p1.second < p2.second; })).second;

        // CCR's F4/B8 FIX: floor minimumScoring at max(absoluteFloor,
        // maximumScoring * dataRelativeFloorRatio) so log-normalisation
        // doesn't produce wildly skewed scores when the spread is huge.
        minimumScoring = std::max(
            defaults.crossReferencerScoringAbsoluteFloor,
            std::max(minimumScoring,
                     maximumScoring * defaults.crossReferencerScoringDataRelativeFloorRatio));

        for (const auto& [name, rawScoring] : aux) {
            // CMCR-style NaN guard on input.
            if (std::isnan(rawScoring)) {
                throw std::invalid_argument(
                    "scoring cannot be nan in compute_normalized_scorings");
            }

            double out;
            if (minimumScoring == maximumScoring) {
                // XC-6 FIX: neutral score (defaults.crossReferencerNeutral
                // ScoreWhenEqual = 0.5) when the population is degenerate.
                out = weighted
                    ? weights[filter] * defaults.crossReferencerNeutralScoreWhenEqual
                    : defaults.crossReferencerNeutralScoreWhenEqual;
            }
            else if (fc["log"]) {
                double n = (std::log10(rawScoring) - std::log10(minimumScoring))
                         / (std::log10(maximumScoring) - std::log10(minimumScoring));
                if (fc["invert"]) n = 1.0 - n;
                out = weighted ? weights[filter] * n : n;
            }
            else {
                double n = (rawScoring - minimumScoring)
                         / (maximumScoring - minimumScoring);
                if (fc["invert"]) n = 1.0 - n;
                out = weighted ? weights[filter] * n : n;
            }

            // CMCR-style NaN guard on output.
            if (std::isnan(out)) {
                throw std::invalid_argument(
                    "swappedScorings[name][filter] cannot be nan in compute_normalized_scorings");
            }
            swappedScorings[name][filter] = out;
        }
    }
    return swappedScorings;
}

// Phase 7 (Option C): templated skeleton of the per-filter base class
// that both CoreCrossReferencer::MagneticCoreFilter and CoreMaterialCross
// Referencer::MagneticCoreFilter open-code. They differ only in
// (a) the FilterEnum keyed on (CoreCrossReferencerFilters vs
//     CoreMaterialCrossReferencerFilters), and
// (b) whether `_validScorings` exists (CCR has it; CMCR doesn't, and
//     its add_scoring() is a single-write to `_scorings`).
//
// We unify on the CCR shape (carries `_validScorings`); CMCR-derived
// filters can ignore it — the map stays empty, no behavioural change.
// The two existing inner MagneticCoreFilter classes become thin
// derivations so that concrete filter classes don't need their
// inheritance edited.
template <typename FilterEnum>
class CrossReferencerFilterBase {
  public:
    std::map<FilterEnum, std::map<std::string, double>>* _scorings = nullptr;
    std::map<FilterEnum, std::map<std::string, bool>>*   _validScorings = nullptr;
    std::map<FilterEnum, std::map<std::string, double>>* _scoredValues = nullptr;
    std::map<FilterEnum, std::map<std::string, bool>>*   _filterConfiguration = nullptr;

    void add_scoring(std::string name, FilterEnum filter, double scoring) {
        if (std::isnan(scoring)) {
            throw std::invalid_argument("scoring cannot be nan");
        }
        if (scoring != -1) {
            if (_validScorings) (*_validScorings)[filter][name] = true;
            (*_scorings)[filter][name] = scoring;
        }
    }
    void add_scored_value(std::string name, FilterEnum filter, double scoredValues) {
        if (scoredValues != -1) {
            (*_scoredValues)[filter][name] = scoredValues;
        }
    }
    void set_scorings(std::map<FilterEnum, std::map<std::string, double>>* scorings) {
        _scorings = scorings;
    }
    void set_valid_scorings(std::map<FilterEnum, std::map<std::string, bool>>* validScorings) {
        _validScorings = validScorings;
    }
    void set_scored_value(std::map<FilterEnum, std::map<std::string, double>>* scoredValues) {
        _scoredValues = scoredValues;
    }
    void set_filter_configuration(std::map<FilterEnum, std::map<std::string, bool>>* filterConfiguration) {
        _filterConfiguration = filterConfiguration;
    }
};

// Phase 7 (Option C, chunk ii): wire the four shared scoring maps onto
// each filter instance in one call. Replaces the 4-call boilerplate per
// filter (×5-6 filters) at the head of each apply_filters() body.
// `validScorings` is nullable — passing nullptr matches the CMCR shape
// (where add_scoring() skips the valid-scorings write).
template <typename FilterEnum>
void wire_cross_referencer_filters(
    std::initializer_list<CrossReferencerFilterBase<FilterEnum>*> filters,
    std::map<FilterEnum, std::map<std::string, double>>* scorings,
    std::map<FilterEnum, std::map<std::string, bool>>*   validScorings,
    std::map<FilterEnum, std::map<std::string, double>>* scoredValues,
    std::map<FilterEnum, std::map<std::string, bool>>*   filterConfiguration)
{
    for (auto* f : filters) {
        f->set_scorings(scorings);
        if (validScorings) f->set_valid_scorings(validScorings);
        f->set_scored_value(scoredValues);
        f->set_filter_configuration(filterConfiguration);
    }
}

// Phase 7 (Option C, dispatch unification): one declarative
// representation of a single cross-referencer step. Each step pairs the
// FilterEnum value (used solely for the "After …" log line) with a
// callable that takes the in-flight ranked list pointer and returns the
// new ranked list. The callable closes over whichever per-filter context
// the concrete cross-referencer needs (referenceCore, inputs, _models,
// temperature, weights, limit) — see CoreCrossReferencer.cpp and
// CoreMaterialCrossReferencer.cpp for the canonical use.
template <typename Item, typename FilterEnum>
struct CrossReferencerStep {
    FilterEnum filter;
    std::function<std::vector<std::pair<Item, double>>(
        std::vector<std::pair<Item, double>>*)> apply;
};

// Run a sequence of CrossReferencerStep<Item, FilterEnum> against an
// initial ranked list. After each step, the surviving list is logged as
// "There are N after filtering by <filter>." under the given log
// category — identical to the format the two open-coded
// magic_enum::enum_for_each switches used. Returns the final ranked list.
//
// Note: this intentionally does NOT iterate over the FilterEnum's value
// set. The caller controls order and may skip values, which preserves
// each cross-referencer's prior behaviour (e.g. CoreCrossReferencer
// handles CORE_LOSSES as a separate trailing step after a prune, and
// would not want it inside the iterated set).
template <typename Item, typename FilterEnum>
std::vector<std::pair<Item, double>> run_cross_referencer_pipeline(
    std::vector<std::pair<Item, double>> ranked,
    const std::vector<CrossReferencerStep<Item, FilterEnum>>& steps,
    const std::string& logCategory)
{
    for (const auto& step : steps) {
        ranked = step.apply(&ranked);
        logEntry(
            "There are " + std::to_string(ranked.size())
                + " after filtering by " + to_string(step.filter) + ".",
            logCategory, 2);
    }
    return ranked;
}

} // namespace OpenMagnetics
