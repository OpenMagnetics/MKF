// Phase 5 (structural cleanup, pipeline orchestration): extracted from
// CoreAdviser.cpp. This translation unit owns the four big
// filter_*_application pipelines plus the pre/post-processing helpers
// (pre_process_inputs, correct_windings, cull_to_unique_core_shapes,
// CoreAdviser::post_process_core) that they call.
//
// Sibling TU of CoreAdviser.cpp; declarations live in advisers/CoreAdviser.h
// and shared helpers in advisers/CoreAdviserInternal.h.

#include "advisers/CoreAdviser.h"
#include "advisers/CoreAdviserInternal.h"
#include "advisers/CoilAdviser.h"             // ABT #4 proximity re-rank winds a representative coil
#include "Constants.h"                         // Constants::residualGap for the DMC gap pre-filter
#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"  // for is_energy_storing_topology() / is_inductor()
#include "processors/MagneticSimulator.h"     // ABT #4 proximity re-rank full-simulates the wound candidate
#include "constructive_models/Bobbin.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/NumberTurns.h"
#include "constructive_models/Wire.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/CoilMesher.h"
#include "support/Exceptions.h"
#include "support/Logger.h"
#include "support/Settings.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <execution>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <magic_enum_utility.hpp>

namespace OpenMagnetics {

void correct_windings(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, const Inputs& inputs) {
    MagnetizingInductance magnetizingInductance;

    // If the wizard supplied per-winding excitation names (e.g. LLC center-tap
    // halves are named "Secondary 0 Half 1" / "Secondary 0 Half 2"), reuse
    // those names on the Coil so the downstream UI shows the right labels,
    // and automatically group center-tap halves into the same section via
    // `wound_with`.
    std::vector<std::string> excitationNames;
    if (!inputs.get_operating_points().empty()) {
        auto& excs = inputs.get_operating_points()[0].get_excitations_per_winding();
        for (auto& exc : excs) {
            excitationNames.push_back(exc.get_name().value_or(""));
        }
        // Only adopt excitation names if they are all distinct (otherwise they
        // are generic labels like "Primary winding excitation" repeated for
        // every winding, which would produce duplicate winding names and
        // break downstream coil processing).
        std::set<std::string> uniq(excitationNames.begin(), excitationNames.end());
        if (uniq.size() != excitationNames.size() || uniq.count("") > 0) {
            excitationNames.clear();
        }
    }

    // The design requirements carry the per-winding isolation side (LLC and
    // other center-tapped topologies put BOTH halves of a center-tap on the
    // same side). Reuse those so the Coil virtualization step can merge the
    // wound-together halves; without this, correct_windings would emit
    // distinct sides from the generic index map and the merge would fail
    // with "Windings wound together must have the same isolation side".
    std::vector<IsolationSide> designRequirementIsolationSides;
    if (inputs.get_design_requirements().get_isolation_sides()) {
        designRequirementIsolationSides = inputs.get_design_requirements().get_isolation_sides().value();
    }
    auto isolation_side_for_index = [&](size_t windingIndex) -> IsolationSide {
        if (windingIndex < designRequirementIsolationSides.size()) {
            return designRequirementIsolationSides[windingIndex];
        }
        return get_isolation_side_from_index(windingIndex);
    };
    auto is_center_tap_pair = [](const std::string& a, const std::string& b) {
        // Detect the "<prefix> Half 1" / "<prefix> Half 2" naming pattern
        // (used by LLC) or the "<prefix>a" / "<prefix>b" pattern (used by
        // AHB CT, e.g. "Secondary 0a" / "Secondary 0b").
        const std::string h1 = " Half 1";
        const std::string h2 = " Half 2";
        if (a.size() > h1.size() && b.size() > h2.size() &&
            a.compare(a.size() - h1.size(), h1.size(), h1) == 0 &&
            b.compare(b.size() - h2.size(), h2.size(), h2) == 0 &&
            a.substr(0, a.size() - h1.size()) == b.substr(0, b.size() - h2.size())) {
            return true;
        }
        if (a.size() >= 2 && b.size() >= 2 &&
            a.back() == 'a' && b.back() == 'b' &&
            a.substr(0, a.size() - 1) == b.substr(0, b.size() - 1)) {
            return true;
        }
        return false;
    };

    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){

        const auto& coil = (*magneticsWithScoring)[i].first.get_coil();
        NumberTurns numberTurns(coil.get_number_turns(0), inputs.get_design_requirements());
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();

        // If a previous step (e.g. the CMC path in get_dummy_coil) already
        // populated secondaries, drop them so the loop below repopulates from
        // the design-requirements turnsRatios. Without this we'd accumulate
        // (dummy secondaries) + (turnsRatios-derived secondaries), producing
        // a coil with too many windings.
        {
            auto& fds = (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description();
            if (fds.size() > 1) {
                fds.resize(1);
            }
        }

        Bobbin bobbin;
        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            bobbin = Bobbin::create_quick_bobbin((*magneticsWithScoring)[i].first.get_core(), true);
        }
        else {
            bobbin = Bobbin::create_quick_bobbin((*magneticsWithScoring)[i].first.get_core());
        }
        (*magneticsWithScoring)[i].first.get_mutable_coil().set_bobbin(bobbin);

        // Set the PRIMARY's name from the excitation if available.
        if (!excitationNames.empty() && !excitationNames[0].empty()) {
            (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description()[0].set_name(excitationNames[0]);
        }

        // Also align the primary's isolation side with the design requirement.
        (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description()[0].set_isolation_side(isolation_side_for_index(0));

        // Every secondary below is a COPY of the primary winding, so it inherits
        // the primary's strand + parallel count. Resolve the shared (skin-depth)
        // wire and the temperature so each secondary's parallels can be re-sized
        // to ITS own current — otherwise a high-current secondary stays at the
        // primary's parallel count and is under-coppered (ABT #70).
        Wire sharedWire = (*magneticsWithScoring)[i].first.get_mutable_coil().resolve_wire(0);
        double windingTemperature = inputs.get_maximum_temperature();

        for (size_t windingIndex = 1; windingIndex < numberTurnsCombination.size(); ++windingIndex) {
            auto winding = coil.get_functional_description()[0];
            winding.set_number_turns(numberTurnsCombination[windingIndex]);
            winding.set_isolation_side(isolation_side_for_index(windingIndex));

            // Size this secondary's parallels to its own RMS current (both the DC
            // ~7 A/mm2 and effective ~12 A/mm2 limits), like get_dummy_coil does
            // for the primary. ABT #126: not for common-mode chokes — see
            // get_dummy_coil; current-sized parallels made every CMC candidate
            // unwindable and the suppression adviser returned zero results.
            if (!inputs.get_operating_points().empty() &&
                !Inputs::can_be_common_mode_choke(inputs.get_operating_points()[0])) {
                const auto& secExcitations = inputs.get_operating_points()[0].get_excitations_per_winding();
                if (windingIndex < secExcitations.size() && secExcitations[windingIndex].get_current()) {
                    auto secCurrent = secExcitations[windingIndex].get_current().value();
                    int npEffective = Wire::calculate_number_parallels_needed(
                        secCurrent, windingTemperature, sharedWire, defaults.maximumEffectiveCurrentDensity);
                    double secDcDensity = sharedWire.calculate_dc_current_density(secCurrent);
                    int npDc = static_cast<int>(std::ceil(secDcDensity / defaults.maximumCurrentDensity));
                    winding.set_number_parallels(std::max({1, npEffective, npDc}));
                }
            }

            // Prefer the excitation name over the generic isolation-side name.
            std::string name;
            if (windingIndex < excitationNames.size() && !excitationNames[windingIndex].empty()) {
                name = excitationNames[windingIndex];
            } else {
                name = get_isolation_side_name_from_index(windingIndex);
            }
            winding.set_name(name);

            // Center-tap grouping: if the previous winding's name ends in
            // "Half 1" and this one ends in "Half 2" with a matching prefix,
            // tie them together with `wound_with` so the Coil's virtualization
            // step merges them into a single virtual winding. The 2-D viewer
            // then draws both halves in the same section / layer.
            if (windingIndex >= 1 &&
                windingIndex < excitationNames.size() &&
                (windingIndex - 1) < excitationNames.size()) {
                const auto& prev = excitationNames[windingIndex - 1];
                const auto& curr = excitationNames[windingIndex];
                if (is_center_tap_pair(prev, curr)) {
                    winding.set_wound_with(std::vector<std::string>{prev});
                }
            }

            (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description().push_back(winding);
        }
    }
}

std::vector<std::pair<Magnetic, double>> cull_to_unique_core_shapes(std::vector<std::pair<Magnetic, double>> magneticsWithScoring, size_t maximumNumberResults) {
    if (magneticsWithScoring.size() > maximumNumberResults) {
        std::vector<std::pair<Magnetic, double>> magneticsWithScoringAndUniqueShapes;
        std::vector<std::string> usedShapes;
        for (size_t magneticIndex = 0; magneticIndex < magneticsWithScoring.size(); ++magneticIndex){
            Magnetic magnetic = magneticsWithScoring[magneticIndex].first;
            auto core = magnetic.get_core();
            if (std::find(usedShapes.begin(), usedShapes.end(), core.get_shape_name()) != usedShapes.end()) {
                continue;
            }
            else {
                magneticsWithScoringAndUniqueShapes.push_back(magneticsWithScoring[magneticIndex]);
                usedShapes.push_back(core.get_shape_name());
            }

            if (magneticsWithScoringAndUniqueShapes.size() == maximumNumberResults) {
                break;
            }
        }
        magneticsWithScoring.clear();

        return magneticsWithScoringAndUniqueShapes;
    }
    return magneticsWithScoring;
}

Mas CoreAdviser::post_process_core(Magnetic magnetic, Inputs inputs) {
    MagneticEnergy magneticEnergy;
    Mas mas;
    double temperature = inputs.get_maximum_temperature();

    magnetic.get_reference();
    if (magnetic.get_manufacturer_info()) {
        auto manufacturerInfo = magnetic.get_manufacturer_info().value();
        manufacturerInfo.set_reference(magnetic.get_core().get_name().value_or("unnamed"));
        magnetic.set_manufacturer_info(manufacturerInfo);
    }

    auto previousCoilDelimitAndCompact = settings.get_coil_delimit_and_compact();
    settings.set_coil_delimit_and_compact(false);
    magnetic.get_mutable_coil().fast_wind();
    settings.set_coil_delimit_and_compact(previousCoilDelimitAndCompact);

    // Store the magnetic AFTER fast_wind so the returned Mas carries the wound
    // coil (sections/layers/turns). Previously set_magnetic ran before
    // fast_wind, so the returned Mas kept the UN-wound copy while only the
    // local `magnetic` got the turns — every fast-path result came back without
    // a turnsDescription, which later threw COIL_NOT_PROCESSED on export
    // (ABT #42). The losses below are computed on the same wound `magnetic`.
    mas.set_magnetic(magnetic);

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        Outputs outputs;
        auto operatingPoint = inputs.get_operating_point(operatingPointIndex);

        MagnetizingInductanceOutput magnetizingInductanceOutput = _magneticSimulator.calculate_magnetizing_inductance(operatingPoint, magnetic);
        // calculate_magnetizing_inductance() re-evaluated the magnetizing current — the flux /
        // saturation driver, which scales with 1/Lm — at THIS core's achieved inductance (Np^2 /
        // reluctance), not the seed Lm that pre_process_inputs used. Persist that recomputed operating
        // point so the exported MAS is self-consistent with its own core: otherwise every candidate
        // ships the seed-Lm magnetizing current regardless of its actual inductance, over-stating the
        // saturating current of transformers (forward/push-pull/bridge) by several-fold and falsely
        // rejecting feasible designs. Inductors are unchanged (their magnetizing current IS the
        // winding current).
        inputs.get_mutable_operating_points()[operatingPointIndex] = operatingPoint;
        double totalStorableMagneticEnergy = magneticEnergy.calculate_core_maximum_magnetic_energy(magnetic.get_core(), operatingPoint);
        auto excitation = Inputs::get_primary_excitation(inputs.get_mutable_operating_points()[operatingPointIndex]);

        magnetizingInductanceOutput.set_maximum_magnetic_energy_core(totalStorableMagneticEnergy);
        magnetizingInductanceOutput.set_method_used(_models["gapReluctance"]);
        magnetizingInductanceOutput.set_origin(ResultOrigin::SIMULATION);


        WindingLossesOutput windingLossesOutput = _windingOhmicLosses.calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, temperature);

        CoreLossesOutput coreLossesOutput = _magneticSimulator.calculate_core_losses(operatingPoint, magnetic);

        InductanceOutput inductanceOutput;
        inductanceOutput.set_magnetizing_inductance(magnetizingInductanceOutput);
        outputs.set_inductance(inductanceOutput);
        outputs.set_winding_losses(windingLossesOutput);
        outputs.set_core_losses(coreLossesOutput);

        mas.get_mutable_outputs().push_back(outputs);
    }
    mas.set_inputs(inputs);

    return mas;
}

Inputs pre_process_inputs(Inputs inputs) {
    double magnetizingInductance = resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance());
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex){
        auto excitation = Inputs::get_primary_excitation(inputs.get_mutable_operating_points()[operatingPointIndex]);
        if (!excitation.get_voltage()) {
            if (magnetizingInductance <= 0) continue;
            inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0].set_voltage(Inputs::calculate_induced_voltage(excitation, magnetizingInductance));
            Inputs::set_current_as_magnetizing_current(&inputs.get_mutable_operating_points()[operatingPointIndex]);
        }
        else if (!excitation.get_magnetizing_current()) {
            if (magnetizingInductance <= 0) {
                // Lm=0 means "not specified" (e.g. LLC transformer). Use the primary
                // current as the magnetizing current so downstream code has it.
                Inputs::set_current_as_magnetizing_current(&inputs.get_mutable_operating_points()[operatingPointIndex]);
                continue;
            }
            auto magnetizingCurrent = Inputs::calculate_magnetizing_current(excitation, magnetizingInductance, false);
            inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
        }
    }
    return inputs;
}


std::vector<std::pair<Mas, double>> CoreAdviser::post_process_and_cut(std::vector<std::pair<Magnetic, double>>& magneticsWithScoring,
                                                                      Inputs& inputs,
                                                                      size_t maximumNumberResults,
                                                                      bool withAlternativeMaterials) {
    // See the declaration comment (ABT #126): walk in score order, keep the first
    // maximumNumberResults candidates that SURVIVE post-processing. Windings are
    // corrected (and alternative materials looked up) per candidate, so tail
    // candidates that are never reached cost nothing.
    std::vector<std::pair<Mas, double>> masWithScoring;
    std::vector<std::string> usedShapes;
    const bool uniqueShapes = get_unique_core_shapes();
    for (auto& magneticWithScoring : magneticsWithScoring) {
        if (masWithScoring.size() >= maximumNumberResults) {
            break;
        }
        std::string shapeName;
        if (uniqueShapes) {
            shapeName = magneticWithScoring.first.get_core().get_shape_name();
            if (std::find(usedShapes.begin(), usedShapes.end(), shapeName) != usedShapes.end()) {
                continue;
            }
        }
        std::vector<std::pair<Magnetic, double>> candidate{magneticWithScoring};
        correct_windings(&candidate, inputs);
        if (withAlternativeMaterials) {
            add_alternative_materials(&candidate, inputs);
        }
        try {
            auto mas = post_process_core(candidate[0].first, inputs);
            masWithScoring.push_back({mas, candidate[0].second});
            if (uniqueShapes) {
                usedShapes.push_back(shapeName);
            }
        } catch (const std::exception& e) {
            // A scored core that cannot actually host the winding (e.g. a toroid
            // whose turns don't physically fit) throws during post-processing --
            // drop it and let the next-scored candidate backfill (ABT #126: the
            // old resize-first order silently returned ZERO results when the
            // top-N were unwindable while windable survivors had been discarded).
            logEntry(std::string("CoreAdviser: dropping core that failed post-processing (backfilling): ") + e.what(), "CoreAdviser", 2);
        }
    }
    return masWithScoring;
}

std::vector<std::pair<Mas, double>> CoreAdviser::filter_available_cores_power_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
    inputs = pre_process_inputs(inputs);

    MagneticCoreFilterAreaProduct filterAreaProduct(inputs);
    MagneticCoreFilterEnergyStored filterEnergyStored(inputs, _models);
    MagneticCoreFilterCost filterCost(inputs);
    MagneticCoreFilterLosses filterLosses(inputs, _models);
    MagneticCoreFilterDimensions filterDimensions;

    filterAreaProduct.set_scorings(&_scorings);
    filterAreaProduct.set_filter_configuration(&_filterConfiguration);
    filterEnergyStored.set_scorings(&_scorings);
    filterEnergyStored.set_filter_configuration(&_filterConfiguration);
    filterCost.set_scorings(&_scorings);
    filterCost.set_filter_configuration(&_filterConfiguration);
    filterLosses.set_scorings(&_scorings);
    filterLosses.set_filter_configuration(&_filterConfiguration);
    filterDimensions.set_scorings(&_scorings);
    filterDimensions.set_filter_configuration(&_filterConfiguration);

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring = *magnetics;
    magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);  // Fixed weight: pre-filtering criterion, not efficiency scoring
    log_stage("Area Product filter", magneticsWithScoring.size());

    if (settings.get_core_adviser_enable_intermediate_pruning() && magneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
        magneticsWithScoring.resize(maximumMagneticsAfterFiltering); // F10 FIX: resize instead of copy-construct
        log_stage("first-filter score-cull", magneticsWithScoring.size());
    }

    magneticsWithScoring = filterEnergyStored.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);  // Fixed weight: pre-filtering criterion, not efficiency scoring
    log_stage("Energy Stored filter", magneticsWithScoring.size());

    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

    // Add saturation filter to reject cores that exceed magnetic flux density saturation
    MagneticCoreFilterSaturation filterSaturationAvailable;
    filterSaturationAvailable.set_scorings(&_scorings);
    filterSaturationAvailable.set_filter_configuration(&_filterConfiguration);
    magneticsWithScoring = filterSaturationAvailable.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    log_stage("Saturation filter", magneticsWithScoring.size());

    // The available-cores path has no fringing filter of its own; the gap is
    // finalized in add_initial_turns_by_inductance above and can be huge for a
    // small core forced to hit L without saturating. Reject those winding-killers.
    reject_winding_killing_gaps(&magneticsWithScoring, inputs);
    log_stage("FringingFactor (post-gap)", magneticsWithScoring.size());

    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::COST], true);
    log_stage("Cost filter", magneticsWithScoring.size());

    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::DIMENSIONS], true);
    log_stage("Dimensions filter", magneticsWithScoring.size());

    // Prune to top candidates by accumulated Cost+Dimensions score before the expensive
    // Loss filter (which sweeps N per core). Cores ranking poorly on cost+size are not
    // going to be promoted by losses alone.
    {
        const size_t preLossCap = std::max<size_t>(maximumNumberResults * 5, 50);
        if (settings.get_core_adviser_enable_intermediate_pruning() && magneticsWithScoring.size() > preLossCap) {
            magneticsWithScoring.resize(preLossCap);
            log_pruned("Core Losses filter", magneticsWithScoring.size());
        }
    }

    magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
    log_stage("Core Losses filter", magneticsWithScoring.size());

    if (settings.get_core_adviser_enable_temperature_filter()) {
        MagneticCoreFilterTemperature filterTemperature(
            inputs, _models, settings.get_core_adviser_maximum_temperature());
        filterTemperature.set_scorings(&_scorings);
        filterTemperature.set_filter_configuration(&_filterConfiguration);
        magneticsWithScoring = filterTemperature.filter_magnetics(
            &magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
        log_stage("Temperature filter", magneticsWithScoring.size());
    }

    // Retry logic: zero results most often means the score-based intermediate
    // pruning cut the larger (saturation-proof) cores BEFORE the saturation
    // filter ran — scoring favors small/cheap cores. Following the classic
    // area-product iteration ("try the next larger core", McLyman), retry with
    // pruning disabled and the saturation filter fully active. The previous
    // behavior (skip the saturation filter entirely) returned cores that can
    // saturate in hardware with no marker.
    if (magneticsWithScoring.size() == 0) {
        logEntry("No cores found with standard filters. Retrying without intermediate pruning (keeping larger cores), saturation still enforced...", "CoreAdviser");

        // Stage 1: full candidate set, saturation enforced
        magneticsWithScoring = *magnetics;
        magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);
        magneticsWithScoring = filterEnergyStored.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);
        add_initial_turns_by_inductance(&magneticsWithScoring, inputs);
        magneticsWithScoring = filterSaturationAvailable.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
        log_stage("retry Saturation (no pruning)", magneticsWithScoring.size());

        // Stage 2: still nothing — relax the saturation MARGIN to 1.0 (still
        // physically non-saturating, just without the production headroom) and
        // tag the results so the relaxation is visible downstream.
        bool reducedSaturationMargin = false;
        if (magneticsWithScoring.size() == 0 && settings.get_core_adviser_saturation_margin() > 1.0) {
            double originalSaturationMargin = settings.get_core_adviser_saturation_margin();
            logEntry("Still no cores. Relaxing saturation margin from " + std::to_string(originalSaturationMargin) + " to 1.0 (results will be tagged)...", "CoreAdviser");
            settings.set_core_adviser_saturation_margin(1.0);
            magneticsWithScoring = *magnetics;
            magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);
            magneticsWithScoring = filterEnergyStored.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);
            add_initial_turns_by_inductance(&magneticsWithScoring, inputs);
            magneticsWithScoring = filterSaturationAvailable.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
            settings.set_core_adviser_saturation_margin(originalSaturationMargin);
            reducedSaturationMargin = !magneticsWithScoring.empty();
            log_stage("retry Saturation (margin 1.0)", magneticsWithScoring.size());
        }

        // Same winding-killing-gap guard as the main flow (the retry re-seeds the
        // gap via add_initial_turns_by_inductance, so it must re-check fringing).
        reject_winding_killing_gaps(&magneticsWithScoring, inputs);
        log_stage("retry FringingFactor (post-gap)", magneticsWithScoring.size());

        magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::COST], true);
        magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::DIMENSIONS], true);
        magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);

        if (settings.get_core_adviser_enable_temperature_filter()) {
            MagneticCoreFilterTemperature filterTemperature(
                inputs, _models, settings.get_core_adviser_maximum_temperature());
            filterTemperature.set_scorings(&_scorings);
            filterTemperature.set_filter_configuration(&_filterConfiguration);
            magneticsWithScoring = filterTemperature.filter_magnetics(
                &magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
        }

        if (reducedSaturationMargin) {
            for (auto& [magnetic, scoring] : magneticsWithScoring) {
                auto coreName = magnetic.get_mutable_core().get_name().value_or("");
                magnetic.get_mutable_core().set_name(coreName + " [SATURATION MARGIN RELAXED TO 1.0]");
            }
        }

        log_stage("retry with relaxed constraints", magneticsWithScoring.size());
    }

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    // ABT #4 added a proximity-aware re-rank of the top-K here (it wound + full-
    // simulated each candidate to reorder by total loss). ABT #14 made it
    // redundant: tightening the post-gap winding-killer fringing ceiling to 1.3
    // (reject_winding_killing_gaps) prunes the catastrophic gap-fringing cores up
    // front, so the surviving pool is already proximity-sane and the default
    // inductor lands a litz design at ~3 W without the per-advise winding cost the
    // re-rank added (~1.8 s/advise, ~96% of it coil winding). Verified: ABT #4's
    // own core-adviser regression still lands ~4.3 W (was ~3.2 W) and the
    // magnetic-adviser default inductor ~3.2 W, both with the re-rank removed.
    // ABT #126: post-process before cutting (backfill; see post_process_and_cut)
    return post_process_and_cut(magneticsWithScoring, inputs, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::filter_available_cores_suppression_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
    inputs = pre_process_inputs(inputs);

    MagneticCoreFilterCost filterCost(inputs);
    MagneticCoreFilterLosses filterLosses(inputs, _models);
    MagneticCoreFilterDimensions filterDimensions;
    MagneticCoreFilterMagneticInductance filterMagneticInductance;
    MagneticCoreFilterMinimumImpedance filterMinimumImpedance;
    MagneticCoreFilterSaturation filterSaturation;
    MagneticCoreFilterTurnCount filterTurnCount;

    filterCost.set_scorings(&_scorings);
    filterCost.set_filter_configuration(&_filterConfiguration);
    filterLosses.set_scorings(&_scorings);
    filterLosses.set_filter_configuration(&_filterConfiguration);
    filterDimensions.set_scorings(&_scorings);
    filterDimensions.set_filter_configuration(&_filterConfiguration);
    filterMagneticInductance.set_scorings(&_scorings);
    filterMagneticInductance.set_filter_configuration(&_filterConfiguration);
    filterMinimumImpedance.set_scorings(&_scorings);
    filterMinimumImpedance.set_filter_configuration(&_filterConfiguration);
    filterSaturation.set_scorings(&_scorings);
    filterSaturation.set_filter_configuration(&_filterConfiguration);
    filterTurnCount.set_scorings(&_scorings);
    filterTurnCount.set_filter_configuration(&_filterConfiguration);

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring = *magnetics;

    // DMC pre-filter: differential-mode chokes carry the full line current
    // through the winding (no flux cancellation). Ungapped high-µᵢ ferrite
    // saturates immediately at typical line currents (B_sat ≈ 0.25 T vs
    // predicted 16-25 T linear B). Industry practice (Richtek AN008,
    // Magnetics CS app notes) is to use distributed-gap powder cores
    // (Sendust/KoolMu/MPP/High-Flux/iron-powder) or amorphous /
    // nanocrystalline. Gapped ferrite is acceptable too.
    auto topology = inputs.get_design_requirements().get_topology();
    if (topology == MAS::Topology::DIFFERENTIAL_MODE_CHOKE) {
        std::vector<std::pair<Magnetic, double>> dmcFiltered;
        dmcFiltered.reserve(magneticsWithScoring.size());
        for (auto& entry : magneticsWithScoring) {
            auto core = entry.first.get_core();
            auto coreMaterial = core.resolve_material();
            auto materialType = coreMaterial.get_material();
            bool isHighSatMaterial = (materialType == MAS::MaterialType::POWDER ||
                                      materialType == MAS::MaterialType::AMORPHOUS ||
                                      materialType == MAS::MaterialType::NANOCRYSTALLINE ||
                                      materialType == MAS::MaterialType::ELECTRICAL_STEEL);
            // A REAL gap, not a residual one: ~1734 catalogue cores carry a RESIDUAL gap
            // entry (Constants::residualGap = 5e-6 m > 0) representing the unavoidable
            // ungapped-core reluctance, so the old `length > 0` test treated ungapped high-µ
            // ferrite as "gapped" and let exactly the instant-saturation cores this filter
            // exists to remove pass through.
            const auto& gapping = core.get_functional_description().get_gapping();
            bool hasGap = false;
            for (const auto& gap : gapping) {
                if (gap.get_type() != MAS::GapType::RESIDUAL &&
                    gap.get_length() > Constants().residualGap) {
                    hasGap = true;
                    break;
                }
            }
            if (isHighSatMaterial || hasGap) {
                dmcFiltered.push_back(entry);
            }
        }
        magneticsWithScoring = std::move(dmcFiltered);
    }

    // (CMC powder pre-filter previously lived here. Removed in favour of the
    // generic TURN_COUNT scoring filter below: low-µ cores like powder
    // toroids in CMC duty need ~8× the turn count of high-µ ferrite to meet
    // Z_min, which the manufacturability proxy now penalises proportionally
    // — without hard-coding a material-type rule that may exclude legitimate
    // candidates in other suppression flows.)

    // NOTE: no add_initial_turns_by_inductance() here, unlike the power/energy
    // paths. For suppression the impedance requirement — not an inductance
    // target — sets the turn count, and filterMinimumImpedance below derives N
    // from scratch (NumberTurns(1) + the analytical Z∝N² jump) and overwrites
    // whatever N a candidate carries. Seeding turns by inductance first was
    // therefore pure waste: it never reached the output yet cost a full
    // per-candidate magnetizing-inductance/gapping solve (~157 ms each) on the
    // entire ~4.6k-core suppression set — minutes of work, and the hang behind
    // ABT #9. The cheap impedance filter (~0.3 ms/core) does the real culling.
    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, 0.001, true);

    magneticsWithScoring = filterTurnCount.filter_magnetics(&magneticsWithScoring, inputs, 1.0, false);

    // Cap before the two expensive per-candidate physical filters below
    // (saturation and losses each cost ~0.1 s/core: they wind the coil and run
    // the full magnetizing/loss model). The impedance survivors are already
    // sorted best-first (least over-dimensioned |Z|), so keeping a bounded
    // shortlist preserves the strongest candidates while keeping the default
    // suppression wizard inside its few-second budget. ABT #9.
    const size_t kMaxCandidatesAfterImpedance = std::max<size_t>(maximumNumberResults * 5, 20);
    if (magneticsWithScoring.size() > kMaxCandidatesAfterImpedance) {
        magneticsWithScoring.resize(kMaxCandidatesAfterImpedance);
    }

    magneticsWithScoring = filterSaturation.filter_magnetics(
        &magneticsWithScoring, inputs, 1, true);

    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::COST], true);

    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::DIMENSIONS], true);

    magneticsWithScoring = filterMagneticInductance.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);

    {
        // Interference-suppression materials are characterised by complex
        // permeability (for impedance), not by a power-loss model, so the loss
        // evaluator rejects all of them as unscorable. Losses are not the ranking
        // criterion for suppression designs (impedance/dimensions/cost/turns are),
        // so if the loss filter wipes out the entire set, pass the pre-loss
        // candidates through unscored rather than returning zero results. This is
        // scoped to the suppression pipeline and logged loudly — NOT the silent
        // global "all-rejected -> return everything" fallback removed earlier.
        auto beforeLosses = magneticsWithScoring;
        magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
        if (magneticsWithScoring.empty() && !beforeLosses.empty()) {
            logEntry("Losses filter could not score any suppression candidate (no power-loss "
                     "model for impedance-characterised materials); passing " +
                         std::to_string(beforeLosses.size()) + " candidates through unscored",
                     "CoreAdviser");
            magneticsWithScoring = std::move(beforeLosses);
        }
    }

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    // ABT #126: post-process before cutting (backfill; see post_process_and_cut)
    return post_process_and_cut(magneticsWithScoring, inputs, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::filter_standard_cores_power_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
    inputs = pre_process_inputs(inputs);

    MagneticCoreFilterAreaProduct filterAreaProduct(inputs);
    MagneticCoreFilterEnergyStored filterEnergyStored(inputs, _models);
    MagneticCoreFilterCost filterCost(inputs);
    MagneticCoreFilterLosses filterLosses(inputs, _models);
    MagneticCoreFilterDimensions filterDimensions;
    MagneticCoreFilterMagneticInductance filterMagneticInductance;
    MagneticCoreFilterFringingFactor filterFringingFactor(inputs, _models);
    MagneticCoreFilterSaturation filterSaturation;

    filterAreaProduct.set_scorings(&_scorings);
    filterAreaProduct.set_filter_configuration(&_filterConfiguration);
    filterAreaProduct.set_cache_usage(false);
    filterEnergyStored.set_scorings(&_scorings);
    filterEnergyStored.set_filter_configuration(&_filterConfiguration);
    filterEnergyStored.set_cache_usage(false);
    filterCost.set_scorings(&_scorings);
    filterCost.set_filter_configuration(&_filterConfiguration);
    filterCost.set_cache_usage(false);
    filterLosses.set_scorings(&_scorings);
    filterLosses.set_filter_configuration(&_filterConfiguration);
    filterLosses.set_cache_usage(false);
    filterDimensions.set_scorings(&_scorings);
    filterDimensions.set_filter_configuration(&_filterConfiguration);
    filterDimensions.set_cache_usage(false);
    filterMagneticInductance.set_scorings(&_scorings);
    filterMagneticInductance.set_filter_configuration(&_filterConfiguration);
    filterMagneticInductance.set_cache_usage(false);
    filterFringingFactor.set_scorings(&_scorings);
    filterFringingFactor.set_filter_configuration(&_filterConfiguration);
    filterFringingFactor.set_cache_usage(false);
    filterSaturation.set_scorings(&_scorings);
    filterSaturation.set_filter_configuration(&_filterConfiguration);
    filterSaturation.set_cache_usage(false);

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring = *magnetics;
    log_stage("Starting", magneticsWithScoring.size());

    bool usingPowderCores = should_include_powder(inputs);

    // User weights scale the tuned per-stage constants (missing entry = 1.0,
    // which reproduces the historical hardcoded behavior). Validity gates
    // (AreaProduct, EnergyStored, Saturation, TurnCount) are not user-scaled.
    auto userWeight = [&weights](CoreAdviserFilters filter) {
        auto it = weights.find(filter);
        return it == weights.end() ? 1.0 : it->second;
    };

    // ========================================================================
    // STEP 1: Area Product filter on all cores
    // ========================================================================
    magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    log_stage("AreaProduct", magneticsWithScoring.size());

    // ========================================================================
    // STEP 2: Separate ferrite (gappable) and powder/toroidal cores
    // Process them independently to avoid cross-contamination
    // ========================================================================
    std::vector<std::pair<Magnetic, double>> ferriteCores;
    std::vector<std::pair<Magnetic, double>> powderCores;
    
    for (const auto& [mag, score] : magneticsWithScoring) {
        if (mag.get_core().get_functional_description().get_type() != CoreType::TOROIDAL) {
            ferriteCores.push_back({mag, score});
        }
    }
    
    // Prune ferrite cores if needed
    size_t ferriteLimit = maximumMagneticsAfterFiltering;
    if (settings.get_core_adviser_enable_intermediate_pruning() && ferriteCores.size() > ferriteLimit) {
        ferriteCores.resize(ferriteLimit);
    }
    logEntry("Ferrite cores after pruning: " + std::to_string(ferriteCores.size()), "CoreAdviser");

    // ========================================================================
    // STEP 3: Process FERRITE cores (gapped)
    // ========================================================================
    if (!ferriteCores.empty()) {
        // Add gaps to ferrite cores
        add_gapping_standard_cores(&ferriteCores, inputs);
        log_stage("gapping ferrite", ferriteCores.size());

        // Filter by fringing factor
        ferriteCores = filterFringingFactor.filter_magnetics(&ferriteCores, inputs, 1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
        log_stage("FringingFactor", ferriteCores.size());

        // Filter by dimensions
        ferriteCores = filterDimensions.filter_magnetics(&ferriteCores, inputs, 1 * userWeight(CoreAdviserFilters::DIMENSIONS), true);
        log_stage("Dimensions (ferrite)", ferriteCores.size());

        // Assign concrete ferrite materials
        ferriteCores = add_ferrite_materials_by_losses(&ferriteCores, inputs);
        log_stage("materials (ferrite)", ferriteCores.size());

        // Calculate turns
        add_initial_turns_by_inductance(&ferriteCores, inputs);

        // Re-check fringing on the FINALIZED gap: add_initial_turns_by_inductance
        // can grow the gap (raising N + re-solving for L) to clear saturation, so
        // the gap the early fringing pass saw is stale. Reject winding-killers here.
        reject_winding_killing_gaps(&ferriteCores, inputs);
        log_stage("FringingFactor (post-gap, ferrite)", ferriteCores.size());

        // Filter by inductance
        ferriteCores = filterMagneticInductance.filter_magnetics(&ferriteCores, inputs, 0.1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
        log_stage("Inductance (ferrite)", ferriteCores.size());

        // Filter by saturation
        ferriteCores = filterSaturation.filter_magnetics(&ferriteCores, inputs, 1, true);
        log_stage("Saturation (ferrite)", ferriteCores.size());

        // Prune to top candidates by accumulated score before the expensive Loss filter
        // (which sweeps N per core). Cores already losing on cost+dimensions+inductance
        // are not going to be promoted by losses alone.
        {
            const size_t preLossCap = std::max<size_t>(maximumNumberResults * 5, 50);
            if (settings.get_core_adviser_enable_intermediate_pruning() && ferriteCores.size() > preLossCap) {
                ferriteCores.resize(preLossCap);
                log_pruned("Losses (ferrite)", ferriteCores.size());
            }
        }

        // Filter by losses
        ferriteCores = filterLosses.filter_magnetics(&ferriteCores, inputs, 1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
        log_stage("Losses (ferrite)", ferriteCores.size());

        if (settings.get_core_adviser_enable_temperature_filter()) {
            MagneticCoreFilterTemperature filterTemperature(
                inputs, _models, settings.get_core_adviser_maximum_temperature());
            filterTemperature.set_scorings(&_scorings);
            filterTemperature.set_filter_configuration(&_filterConfiguration);
            filterTemperature.set_cache_usage(false);
            ferriteCores = filterTemperature.filter_magnetics(&ferriteCores, inputs, 1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
            log_stage("Temperature (ferrite)", ferriteCores.size());
        }
    }

    // ========================================================================
    // STEP 4: Process POWDER cores (ungapped toroidal) - if appropriate or if toroidal cores explicitly requested
    // ========================================================================
    if (usingPowderCores || settings.get_use_toroidal_cores()) {
        // Get toroidal cores from original set
        for (const auto& [mag, score] : magneticsWithScoring) {
            if (mag.get_core().get_functional_description().get_type() == CoreType::TOROIDAL) {
                powderCores.push_back({mag, score});
            }
        }
        
        // Prune powder cores
        size_t powderLimit = maximumMagneticsAfterFiltering;
        if (settings.get_core_adviser_enable_intermediate_pruning() && powderCores.size() > powderLimit) {
            powderCores.resize(powderLimit);
        }
        
        if (!powderCores.empty()) {
            // Add powder materials
            powderCores = add_powder_materials(&powderCores, inputs);
            log_stage("powder materials", powderCores.size());
            
            // Filter by energy stored
            powderCores = filterEnergyStored.filter_magnetics(&powderCores, inputs, 1, true);
            log_stage("EnergyStored (powder)", powderCores.size());
            
            // Filter by dimensions
            powderCores = filterDimensions.filter_magnetics(&powderCores, inputs, 1 * userWeight(CoreAdviserFilters::DIMENSIONS), true);
            log_stage("Dimensions (powder)", powderCores.size());

            // Prune to top candidates by accumulated score before the
            // expensive `add_initial_turns_by_inductance` +
            // `filterMagneticInductance` pair. `add_powder_materials` blows
            // the candidate set up by an order of magnitude (one entry per
            // core × powder material — e.g. 400 cores × 10 materials = 4000)
            // and the inductance step takes ~70 ms per candidate. Without
            // this cap the powder path dominates PFC's 5-minute runtime
            // (≈4 min on Inductance alone, profiled). Use the same
            // top-K-by-score policy that the ferrite path already applies
            // before its Losses filter (see line 608 in this file).
            {
                const size_t preInductanceCap = std::max<size_t>(maximumNumberResults * 5, 50);
                if (settings.get_core_adviser_enable_intermediate_pruning() && powderCores.size() > preInductanceCap) {
                    powderCores.resize(preInductanceCap);
                    log_pruned("Inductance (powder)", powderCores.size());
                }
            }

            // Calculate turns
            add_initial_turns_by_inductance(&powderCores, inputs);

            // Filter by inductance
            powderCores = filterMagneticInductance.filter_magnetics(&powderCores, inputs, 0.1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
            log_stage("Inductance (powder)", powderCores.size());
            
            // Filter by saturation
            powderCores = filterSaturation.filter_magnetics(&powderCores, inputs, 1, true);
            log_stage("Saturation (powder)", powderCores.size());

            // Prune to top candidates by accumulated score before the expensive Loss filter.
            {
                const size_t preLossCap = std::max<size_t>(maximumNumberResults * 5, 50);
                if (settings.get_core_adviser_enable_intermediate_pruning() && powderCores.size() > preLossCap) {
                    powderCores.resize(preLossCap);
                    log_pruned("Losses (powder)", powderCores.size());
                }
            }

            // Filter by losses
            powderCores = filterLosses.filter_magnetics(&powderCores, inputs, 1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
            log_stage("Losses (powder)", powderCores.size());

            if (settings.get_core_adviser_enable_temperature_filter()) {
                MagneticCoreFilterTemperature filterTemperature(
                    inputs, _models, settings.get_core_adviser_maximum_temperature());
                filterTemperature.set_scorings(&_scorings);
                filterTemperature.set_filter_configuration(&_filterConfiguration);
                filterTemperature.set_cache_usage(false);
                powderCores = filterTemperature.filter_magnetics(&powderCores, inputs, 1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
                log_stage("Temperature (powder)", powderCores.size());
            }
        }
    }

    // ========================================================================
    // STEP 5: Combine ferrite and powder results
    // ========================================================================
    magneticsWithScoring.clear();
    
    // Add ferrite cores first (they typically perform better for gapped applications)
    for (const auto& core : ferriteCores) {
        magneticsWithScoring.push_back(core);
    }
    
    // Add powder cores
    for (const auto& core : powderCores) {
        magneticsWithScoring.push_back(core);
    }
    
    logEntry("Combined results: " + std::to_string(magneticsWithScoring.size()) + 
             " (ferrite: " + std::to_string(ferriteCores.size()) + 
             ", powder: " + std::to_string(powderCores.size()) + ")", "CoreAdviser");

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    // ========================================================================
    // STEP 6: Final processing
    // ========================================================================
    // ABT #4's proximity-aware re-rank was removed here in ABT #14; the post-gap
    // fringing ceiling (1.3) now prunes the catastrophic gap-fringing cores up
    // front, making the per-advise re-rank winding cost redundant. See the
    // matching note in the standard-cores pipeline above.
    // ABT #126: post-process before cutting (backfill; see post_process_and_cut)
    return post_process_and_cut(magneticsWithScoring, inputs, maximumNumberResults, /*withAlternativeMaterials=*/true);
}

std::vector<std::pair<Mas, double>> CoreAdviser::filter_standard_cores_interference_suppression_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
    inputs = pre_process_inputs(inputs);

    MagneticCoreFilterLosses filterLosses(inputs, _models);
    MagneticCoreFilterDimensions filterDimensions;
    MagneticCoreFilterMinimumImpedance filterMinimumImpedance;
    MagneticCoreFilterMagneticInductance filterMagneticInductance;
    MagneticCoreFilterCost filterCost(inputs);
    MagneticCoreFilterTurnCount filterTurnCount;


    filterLosses.set_scorings(&_scorings);
    filterLosses.set_filter_configuration(&_filterConfiguration);
    filterLosses.set_cache_usage(false);
    filterDimensions.set_scorings(&_scorings);
    filterDimensions.set_filter_configuration(&_filterConfiguration);
    filterDimensions.set_cache_usage(false);
    filterMinimumImpedance.set_scorings(&_scorings);
    filterMinimumImpedance.set_filter_configuration(&_filterConfiguration);
    filterMinimumImpedance.set_cache_usage(false);
    filterMagneticInductance.set_scorings(&_scorings);
    filterMagneticInductance.set_filter_configuration(&_filterConfiguration);
    filterMagneticInductance.set_cache_usage(false);
    filterCost.set_scorings(&_scorings);
    filterCost.set_filter_configuration(&_filterConfiguration);
    filterCost.set_cache_usage(false);
    filterTurnCount.set_scorings(&_scorings);
    filterTurnCount.set_filter_configuration(&_filterConfiguration);
    filterTurnCount.set_cache_usage(false);

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring = *magnetics;

    // if (magneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
    //     magneticsWithScoring.resize(maximumMagneticsAfterFiltering); // F10 FIX: resize instead of copy-construct
    // }


    magneticsWithScoring = add_ferrite_materials_by_impedance(&magneticsWithScoring, inputs);

    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

    // User weights scale the tuned per-stage constants (missing entry = 1.0,
    // reproducing the historical hardcoded behavior).
    auto userWeight = [&weights](CoreAdviserFilters filter) {
        auto it = weights.find(filter);
        return it == weights.end() ? 1.0 : it->second;
    };

    // Same logic as filter_available_cores_suppression_application: use a tiny
    // weight so the impedance filter is a pure gate (no score contribution).
    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, 0.001, true);
    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, 1 * userWeight(CoreAdviserFilters::COST), true);
    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, inputs, 1 * userWeight(CoreAdviserFilters::DIMENSIONS), true);
    magneticsWithScoring = filterMagneticInductance.filter_magnetics(&magneticsWithScoring, inputs, 1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
    {
        // See filter_available_cores_suppression_application: suppression materials
        // have no power-loss model, so the loss filter can reject every candidate.
        // Losses don't rank suppression designs — pass the pre-loss set through
        // unscored (logged) rather than returning zero results.
        auto beforeLosses = magneticsWithScoring;
        magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, 1 * userWeight(CoreAdviserFilters::EFFICIENCY), true);
        if (magneticsWithScoring.empty() && !beforeLosses.empty()) {
            logEntry("Losses filter could not score any suppression candidate (no power-loss "
                     "model for impedance-characterised materials); passing " +
                         std::to_string(beforeLosses.size()) + " candidates through unscored",
                     "CoreAdviser");
            magneticsWithScoring = std::move(beforeLosses);
        }
    }
    // Manufacturability proxy — see filter_available_cores_suppression_application
    // for rationale. Same internal weight (1.0) here.
    magneticsWithScoring = filterTurnCount.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    // ABT #126: post-process before cutting (backfill; see post_process_and_cut)
    return post_process_and_cut(magneticsWithScoring, inputs, maximumNumberResults, /*withAlternativeMaterials=*/true);
}

} // namespace OpenMagnetics
