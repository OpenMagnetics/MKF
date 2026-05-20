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

} // namespace OpenMagnetics
