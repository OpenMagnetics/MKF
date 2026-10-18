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
#include <optional>
#include <map>
#include <numeric>
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
// survivors carry the possibly-mutated magnetic.
//
// ABT #1785: the pool is CONSUMED. Every caller replaces its pool with the result
// (`pool = filter.filter_magnetics(&pool, ...)`), so the survivors are compacted in place and the
// vector's buffer is handed back; `*unfilteredMagnetics` is left empty. This used to build the
// survivors as a second vector of copies beside the untouched input, so every filter held two
// copies of the pool at once (a Magnetic has no move operations: a "move" is a deep copy, the
// inline core material included), and the stacked second search of the available-cores adviser
// ran past 19 GB. A filter that must not consume its source uses evaluate_and_cull_from below.
std::pair<std::vector<std::pair<Magnetic, double>>, std::vector<double>>
evaluate_and_cull(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics,
                  Inputs& inputs, const EvalFn& evaluate, bool writeBack,
                  const std::string& filterName) {
    auto& pool = *unfilteredMagnetics;
    std::vector<double> newScoring;
    size_t kept = 0;

    for (size_t magneticIndex = 0; magneticIndex < pool.size(); ++magneticIndex) {
        bool valid = false;
        double scoring = 0.0;
        // A write-back filter keeps whatever its evaluation did to a valid candidate, and an
        // invalid one is dropped, so it can evaluate the pool's own element. The others must not
        // keep the evaluation's side effects: they evaluate a copy, as before.
        std::optional<Magnetic> scratch;
        if (!writeBack) {
            scratch.emplace(pool[magneticIndex].first);
        }
        Magnetic* magnetic = writeBack ? &pool[magneticIndex].first : &scratch.value();
        try {
            auto evalResult = evaluate(magnetic, &inputs);
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
            if (kept != magneticIndex) {
                pool[kept] = pool[magneticIndex];
            }
            ++kept;
        }
    }
    // Shrinking only destroys the tail: no reallocation, no copies.
    while (pool.size() > kept) {
        pool.pop_back();
    }
    std::vector<std::pair<Magnetic, double>> filtered;
    filtered.swap(pool);

    if (filtered.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR,
            "Something wrong happened while filtering (" + filterName + "), size of filtered: " +
            std::to_string(filtered.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }
    return {filtered, newScoring};
}

// The non-consuming twin of evaluate_and_cull, for a source pool that has to stay intact: each
// candidate is evaluated on a copy (side effects discarded, as a non-write-back filter does).
// Nothing is copied out here: it returns the survivors' INDICES into `source` and their scores,
// so the caller copies only the survivors it keeps.
std::pair<std::vector<size_t>, std::vector<double>>
evaluate_and_cull_from(const std::vector<std::pair<Magnetic, double>>& source,
                       Inputs& inputs, const EvalFn& evaluate, const std::string& filterName) {
    std::vector<size_t> survivors;
    std::vector<double> newScoring;
    for (size_t magneticIndex = 0; magneticIndex < source.size(); ++magneticIndex) {
        Magnetic magnetic = source[magneticIndex].first;
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
            survivors.push_back(magneticIndex);
            newScoring.push_back(scoring);
        }
    }
    return {std::move(survivors), std::move(newScoring)};
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

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterAreaProduct::filter_magnetics_from(const std::vector<std::pair<Magnetic, double>>& source, Inputs inputs, double weight, std::optional<size_t> keepBest) {
    // Same gate and scoring as filter_magnetics (writeBack=false there too), source left intact.
    auto [survivors, newScoring] = evaluate_and_cull_from(
        source, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        "AreaProduct");

    // Scores and ORDER are decided on the survivors' indices, exactly as apply_normalized_scoring
    // decides them on a copied pool: normalized over ALL survivors, every survivor's scoring
    // recorded, stable descending sort. Only then are the candidates copied, and with keepBest
    // only the first keepBest of that order (ABT #1785): the result is identical to copying every
    // survivor, sorting and resize(keepBest), without holding the survivors the cull drops (tens
    // of thousands of candidates on the available-cores adviser's stacked search).
    std::vector<double> scores(survivors.size());
    for (size_t i = 0; i < survivors.size(); ++i) {
        scores[i] = source[survivors[i]].second;
    }
    std::vector<size_t> order(survivors.size());
    std::iota(order.begin(), order.end(), size_t(0));
    if (weight > 0 && !survivors.empty()) {
        auto normalizedScoring = normalize_scoring(newScoring, weight,
                                                   (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        if (normalizedScoring.size() != survivors.size()) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR,
                "AreaProduct: " + std::to_string(normalizedScoring.size()) + " normalized scores for " +
                std::to_string(survivors.size()) + " survivors");
        }
        for (size_t i = 0; i < survivors.size(); ++i) {
            scores[i] += normalizedScoring[i];
            add_scoring(source[survivors[i]].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY,
                        normalizedScoring[i]);
        }
        std::stable_sort(order.begin(), order.end(), [&scores](size_t left, size_t right) {
            return scores[left] > scores[right];
        });
    }
    size_t kept = order.size();
    if (keepBest && kept > keepBest.value()) {
        kept = keepBest.value();
    }
    std::vector<std::pair<Magnetic, double>> filtered;
    filtered.reserve(kept);
    for (size_t position = 0; position < kept; ++position) {
        filtered.emplace_back(source[survivors[order[position]]].first, scores[order[position]]);
    }
    return filtered;
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

CoreAdviser::MagneticCoreFilterLosses::MagneticCoreFilterLosses(Inputs inputs, std::map<std::string, std::string> models, MAS::MagneticApplication application) {
    _filter = MagneticFilterCoreDcAndSkinLosses(inputs, models); // F23 FIX: DC+skin effect losses
    _filter.set_application(application);
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
        return std::exchange(*unfilteredMagnetics, {});  // consumed, not copied (ABT #1785)
    }
    // Restore coil_delimit_and_compact even if the loss evaluation (which flips it) throws.
    auto& settings = Settings::GetInstance();
    SettingsGuard<bool> coilCompactGuard(settings, &Settings::get_coil_delimit_and_compact,
                                         &Settings::set_coil_delimit_and_compact,
                                         settings.get_coil_delimit_and_compact());
    // A multi-winding input's candidates carry a one-winding stand-in coil at this stage. The
    // loss filter scores them on the coil completed from the turns ratios: core losses plus the
    // DC and skin-effect losses of every winding (MagneticFilterCoreDcAndSkinLosses::
    // evaluate_stand_in_with_every_winding). They used to pass through unscored, so the copper of
    // a 10 kW transformer never counted before the coil stage.

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
        return std::exchange(*unfilteredMagnetics, {});  // consumed, not copied (ABT #1785)
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
        return std::exchange(*unfilteredMagnetics, {});  // consumed, not copied (ABT #1785)
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
    Inputs inputs, std::map<std::string, std::string> models, double maximumTemperature, MAS::MagneticApplication application)
{
    _filter = MagneticFilterTemperature(inputs, maximumTemperature);
    _filter.set_application(application);
    _filter.set_size_stand_in_copper(true);
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
        return std::exchange(*unfilteredMagnetics, {});  // consumed, not copied (ABT #1785)
    }
    auto [filteredMagneticsWithScoring, newScoring] = evaluate_and_cull(
        unfilteredMagnetics, inputs,
        [this](Magnetic* m, Inputs* i) { return _filter.evaluate_magnetic(m, i); },
        /*writeBack=*/true, "Temperature");
    apply_normalized_scoring(&filteredMagneticsWithScoring, newScoring, weight, CoreAdviser::CoreAdviserFilters::EFFICIENCY);
    return filteredMagneticsWithScoring;
}

void CoreAdviser::filter_by_window_copper_capacity(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring,
                                                   const Inputs& inputs, const std::string& stage) {
    if (magneticsWithScoring->empty()) {
        return;
    }
    if (!MagneticFilterWindowCopperCapacity::applies_to(inputs)) {
        logEntry("Window copper capacity" + stage + ": printed winding, its copper is the board's; not screened", "CoreAdviser", 2);
        return;
    }
    MagneticFilterWindowCopperCapacity filter(inputs, defaults.maximumEffectiveCurrentDensity);
    Inputs evaluationInputs = inputs;
    std::vector<std::pair<Magnetic, double>> kept;
    kept.reserve(magneticsWithScoring->size());
    for (auto& [magnetic, scoring] : *magneticsWithScoring) {
        auto [valid, proportion] = filter.evaluate_magnetic(&magnetic, &evaluationInputs);
        if (valid) {
            kept.emplace_back(std::move(magnetic), scoring);
        }
        else {
            logEntry("Window copper capacity" + stage + ": rejected " + filter.get_last_reason(), "CoreAdviser", 2);
        }
    }
    *magneticsWithScoring = std::move(kept);
    log_stage("Window copper capacity" + stage, magneticsWithScoring->size());
}

} // namespace OpenMagnetics
