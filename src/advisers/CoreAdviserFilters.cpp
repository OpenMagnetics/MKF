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
    // ABT #121.1: the hard VALIDITY gate runs regardless of weight; only the
    // SCORE contribution is preference-weighted. A zero user weight must drop
    // this filter's ranking influence, NOT disable its feasibility cull.
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/false, "AreaProduct");
    if (weight > 0) {
        apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterEnergyStored::MagneticCoreFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterEnergyStored(inputs, models);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterEnergyStored::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    // ABT #121.1: validity gate always; score only when weighted (see AreaProduct).
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "EnergyStored");
    if (weight > 0) {
        apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterFringingFactor::MagneticCoreFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterFringingFactor(inputs, models);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterFringingFactor::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    // ABT #121.1: the fringing-factor cull is a winding-killer feasibility gate
    // (a gap so large the coil cannot be wound / proximity losses explode) —
    // it must run even when the EFFICIENCY weight scaling it is zero. Score only
    // when weighted.
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "FringingFactor");
    if (weight > 0) {
        apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterCost::MagneticCoreFilterCost(Inputs inputs) {
    _filter = MagneticFilterEstimatedCost(inputs);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterCost::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    // ABT #121.1: the Cost filter also carries the winding-window-fit feasibility
    // gate (rejects candidates whose turns cannot physically fit the window).
    // Run it regardless of the COST weight; score only when weighted.
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/false, "Cost");
    if (weight > 0) {
        apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::COST);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterLosses::MagneticCoreFilterLosses(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterCoreDcAndSkinLosses(inputs, models); // F23 FIX: DC+skin effect losses
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterLosses::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        // ABT #121.1: unlike the cheap gate filters above, the loss filter's
        // "gate" is a per-candidate N-sweep simulation (seconds over the pool),
        // and its only hard reject — the coil could not be wound — is already
        // enforced upstream by the Cost filter's winding-window-fit gate, which
        // now runs at any weight. Losses is therefore an efficiency-optimisation
        // filter with no independent feasibility gate, so at zero EFFICIENCY
        // weight we skip it (and its expensive simulation) rather than run it
        // purely to reproduce a cull the pipeline already performed.
        logEntry("Losses filter skipped at zero weight (optimisation-only; coil-fit "
                 "feasibility is gated by the Cost filter).", "CoreAdviser", 2);
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
    // ABT #121.1: the Dimensions filter carries the planar-core-for-wound-wire
    // hard reject. Run that feasibility gate regardless of the DIMENSIONS weight;
    // score only when weighted.
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/false, "Dimensions");
    if (weight > 0) {
        apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::DIMENSIONS);
    }
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
        // ABT #121.1: the minimum-impedance gate IS a hard feasibility check, but
        // it is an expensive complex-permeability solve and the pipeline always
        // invokes it with a fixed positive gating weight (0.001) — it is never
        // user-zeroed. Guard defensively with a note rather than run the solve at
        // zero weight; if a caller ever passes weight<=0 the impedance
        // feasibility is intentionally not enforced here.
        logEntry("MinimumImpedance filter skipped at zero weight (expensive solve; "
                 "pipeline always gates it at a fixed positive weight).", "CoreAdviser", 2);
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
    // ABT #121.1: meeting the magnetizing-inductance requirement is spec
    // compliance, not an efficiency preference — the gate must run even when the
    // (small) EFFICIENCY weight scaling it is zero, or the adviser would return
    // cores that miss the required inductance. Run the gate always; score only
    // when weighted. (This per-candidate solve is not free, but a correct answer
    // outranks the runtime saved by skipping it at zero weight.)
    auto& settings = Settings::GetInstance();
    SettingsGuard<bool> coilCompactGuard(settings, &Settings::get_coil_delimit_and_compact,
                                         &Settings::set_coil_delimit_and_compact,
                                         settings.get_coil_delimit_and_compact());
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "MagneticInductance");
    if (weight > 0) {
        apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterSaturation::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    // ABT #121.1: the saturation cull (B_peak under the safe flux ceiling) is a
    // hard feasibility gate — a saturating core is never a valid design — so it
    // runs regardless of weight. Score only when weighted AND not the first
    // (pure-gating) pass.
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "Saturation");
    // Saturation only contributes a headroom SCORE when it is not the first (gating) pass.
    if (!firstFilter && weight > 0) {
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
        // ABT #121.1: the temperature-rise cull is a full thermal solve per
        // candidate (expensive) and is itself gated off by default
        // (core_adviser_enable_temperature_filter). It is an efficiency-scoped
        // preference, so at zero EFFICIENCY weight we skip it (and its solve)
        // rather than run it purely to gate.
        logEntry("Temperature filter skipped at zero weight (expensive thermal solve; "
                 "efficiency-scoped preference).", "CoreAdviser", 2);
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
