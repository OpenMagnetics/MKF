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
#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"  // for is_energy_storing_topology()
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

void correct_windings(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
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

        Coil coil = Coil((*magneticsWithScoring)[i].first.get_coil());
        NumberTurns numberTurns(coil.get_number_turns(0), inputs.get_design_requirements());
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();

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

        for (size_t windingIndex = 1; windingIndex < numberTurnsCombination.size(); ++windingIndex) {
            auto winding = coil.get_functional_description()[0];
            winding.set_number_turns(numberTurnsCombination[windingIndex]);
            winding.set_isolation_side(isolation_side_for_index(windingIndex));

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
    mas.set_magnetic(magnetic);
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

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        Outputs outputs;
        auto operatingPoint = inputs.get_operating_point(operatingPointIndex);

        MagnetizingInductanceOutput magnetizingInductanceOutput = _magneticSimulator.calculate_magnetizing_inductance(operatingPoint, magnetic);
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

    // Retry logic: if no cores found, try relaxing the saturation constraint
    if (magneticsWithScoring.size() == 0) {
        logEntry("No cores found with standard filters. Retrying with relaxed saturation constraint...", "CoreAdviser");

        // Start over with the original magnetics
        magneticsWithScoring = *magnetics;
        magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);

        if (settings.get_core_adviser_enable_intermediate_pruning() && magneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
            magneticsWithScoring.resize(maximumMagneticsAfterFiltering);
        }

        magneticsWithScoring = filterEnergyStored.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);
        add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

        // Skip saturation filter - go directly to scoring filters
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

        log_stage("retry with relaxed constraints", magneticsWithScoring.size());
    }

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    if (magneticsWithScoring.size() > maximumNumberResults) {
        if (get_unique_core_shapes()) {
            magneticsWithScoring = cull_to_unique_core_shapes(magneticsWithScoring, maximumNumberResults);
        }
        else {
            magneticsWithScoring.resize(maximumNumberResults); // F10 FIX: resize instead of copy-construct
        }
    }

    correct_windings(&magneticsWithScoring, inputs);

    std::vector<std::pair<Mas, double>> masWithScoring;

    for (const auto& [magnetic, scoring] : magneticsWithScoring) {
        auto mas = post_process_core(magnetic, inputs);
        masWithScoring.push_back({mas, scoring});
    }

    return masWithScoring;
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
    if (topology == Topologies::DIFFERENTIAL_MODE_CHOKE) {
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
            bool hasGap = !core.get_functional_description().get_gapping().empty() &&
                          core.get_functional_description().get_gapping()[0].get_length() > 0;
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

    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

    // Use a tiny weight (0.001) so the impedance filter acts as a pure hard gate:
    // it eliminates cores that cannot meet |Z| spec, but contributes essentially
    // zero to the cumulative score. Without this, cores that far overshoot the
    // impedance spec (many-turn powder toroids, Z >> Z_min) score as highly as
    // 2.5 on impedance alone, overwhelming the TURN_COUNT score. A pure gate
    // lets TURN_COUNT — which correctly penalises excessive turn counts — be
    // the primary ranking signal. The tiny non-zero weight is needed because
    // weight ≤ 0 would short-circuit the filter entirely (no elimination occurs).
    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, 0.001, true);

    // Turn-count scoring: the primary ranking signal for the
    // interference-suppression pipeline. Fewer turns → less copper burden
    // → easier to wind. Weight 1.0 keeps this filter on the same scale as
    // the user's COST/DIMENSIONS/LOSSES weights (each 0–1 after normalization),
    // so high-µ ferrite (13–20 turns) naturally outranks low-µ powder (50–100
    // turns) without overriding the user's other preferences entirely.
    magneticsWithScoring = filterTurnCount.filter_magnetics(&magneticsWithScoring, inputs, 1.0, false);

    // Hard-cap candidates after the (cheap) impedance + turns-density filters to
    // keep the expensive downstream filters (saturation, losses, magnetizing-
    // inductance) bounded. The top-N here are now the ones with the best
    // combined impedance-margin + turns-density score.
    constexpr size_t kMaxCandidatesAfterImpedance = 50;
    if (magneticsWithScoring.size() > kMaxCandidatesAfterImpedance) {
        magneticsWithScoring.resize(kMaxCandidatesAfterImpedance);
    }

    // Saturation gate: the CMC path already drops the DM line-current DC bias
    // via Inputs::can_be_common_mode_choke, so this filter sees only the CM
    // ripple B. Underlying MagneticFilterSaturation compares B_peak > Bsat
    // strictly (no 0.7 headroom), so only genuinely-saturated cores drop
    // out. Without this the DIMENSIONS filter can pick a too-small toroid
    // that meets the impedance spec but saturates under normal CM noise.
    magneticsWithScoring = filterSaturation.filter_magnetics(&magneticsWithScoring, inputs, 1, true);

    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::COST], true);

    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::DIMENSIONS], true);

    magneticsWithScoring = filterMagneticInductance.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);

    magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    if (magneticsWithScoring.size() > maximumNumberResults) {
        if (get_unique_core_shapes()) {
            magneticsWithScoring = cull_to_unique_core_shapes(magneticsWithScoring, maximumNumberResults);
        }
        else {
            magneticsWithScoring.resize(maximumNumberResults); // F10 FIX: resize instead of copy-construct
        }
    }

    correct_windings(&magneticsWithScoring, inputs);

    std::vector<std::pair<Mas, double>> masWithScoring;

    for (const auto& [magnetic, scoring] : magneticsWithScoring) {
        auto mas = post_process_core(magnetic, inputs);
        masWithScoring.push_back({mas, scoring});
    }

    return masWithScoring;
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
    std::cout << "[CoreAdviser] Starting with " << magneticsWithScoring.size() << " magnetics" << std::endl;

    bool usingPowderCores = should_include_powder(inputs);

    // ========================================================================
    // STEP 1: Area Product filter on all cores
    // ========================================================================
    magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    log_stage("AreaProduct", magneticsWithScoring.size());
    std::cout << "[CoreAdviser] After AreaProduct: " << magneticsWithScoring.size() << std::endl;

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
    std::cout << "[CoreAdviser] Ferrite cores after pruning: " << ferriteCores.size() << std::endl;

    // ========================================================================
    // STEP 3: Process FERRITE cores (gapped)
    // ========================================================================
    if (!ferriteCores.empty()) {
        // Add gaps to ferrite cores
        add_gapping_standard_cores(&ferriteCores, inputs);
        log_stage("gapping ferrite", ferriteCores.size());
        std::cout << "[CoreAdviser] After gapping: " << ferriteCores.size() << std::endl;
        
        // Filter by fringing factor
        ferriteCores = filterFringingFactor.filter_magnetics(&ferriteCores, inputs, 1, true);
        log_stage("FringingFactor", ferriteCores.size());
        std::cout << "[CoreAdviser] After FringingFactor: " << ferriteCores.size() << std::endl;
        
        // Filter by dimensions
        ferriteCores = filterDimensions.filter_magnetics(&ferriteCores, inputs, 1, true);
        log_stage("Dimensions (ferrite)", ferriteCores.size());
        std::cout << "[CoreAdviser] After Dimensions: " << ferriteCores.size() << std::endl;
        
        // Assign concrete ferrite materials
        ferriteCores = add_ferrite_materials_by_losses(&ferriteCores, inputs);
        log_stage("materials (ferrite)", ferriteCores.size());
        std::cout << "[CoreAdviser] After materials: " << ferriteCores.size() << std::endl;
        
        // Calculate turns
        add_initial_turns_by_inductance(&ferriteCores, inputs);
        
        // Filter by inductance
        ferriteCores = filterMagneticInductance.filter_magnetics(&ferriteCores, inputs, 0.1, true);
        log_stage("Inductance (ferrite)", ferriteCores.size());
        std::cout << "[CoreAdviser] After Inductance: " << ferriteCores.size() << std::endl;
        
        // Filter by saturation
        ferriteCores = filterSaturation.filter_magnetics(&ferriteCores, inputs, 1, true);
        log_stage("Saturation (ferrite)", ferriteCores.size());
        std::cout << "[CoreAdviser] After Saturation: " << ferriteCores.size() << std::endl;

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
        ferriteCores = filterLosses.filter_magnetics(&ferriteCores, inputs, 1, true);
        log_stage("Losses (ferrite)", ferriteCores.size());
        std::cout << "[CoreAdviser] After Losses: " << ferriteCores.size() << std::endl;

        if (settings.get_core_adviser_enable_temperature_filter()) {
            MagneticCoreFilterTemperature filterTemperature(
                inputs, _models, settings.get_core_adviser_maximum_temperature());
            filterTemperature.set_scorings(&_scorings);
            filterTemperature.set_filter_configuration(&_filterConfiguration);
            filterTemperature.set_cache_usage(false);
            ferriteCores = filterTemperature.filter_magnetics(&ferriteCores, inputs, 1, true);
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
            powderCores = filterDimensions.filter_magnetics(&powderCores, inputs, 1, true);
            log_stage("Dimensions (powder)", powderCores.size());
            
            // Calculate turns
            add_initial_turns_by_inductance(&powderCores, inputs);
            
            // Filter by inductance
            powderCores = filterMagneticInductance.filter_magnetics(&powderCores, inputs, 0.1, true);
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
            powderCores = filterLosses.filter_magnetics(&powderCores, inputs, 1, true);
            log_stage("Losses (powder)", powderCores.size());

            if (settings.get_core_adviser_enable_temperature_filter()) {
                MagneticCoreFilterTemperature filterTemperature(
                    inputs, _models, settings.get_core_adviser_maximum_temperature());
                filterTemperature.set_scorings(&_scorings);
                filterTemperature.set_filter_configuration(&_filterConfiguration);
                filterTemperature.set_cache_usage(false);
                powderCores = filterTemperature.filter_magnetics(&powderCores, inputs, 1, true);
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
    if (magneticsWithScoring.size() > maximumNumberResults) {
        if (get_unique_core_shapes()) {
            magneticsWithScoring = cull_to_unique_core_shapes(magneticsWithScoring, maximumNumberResults);
        }
        else {
            magneticsWithScoring.resize(maximumNumberResults);
        }
    }

    correct_windings(&magneticsWithScoring, inputs);
    add_alternative_materials(&magneticsWithScoring, inputs);

    std::vector<std::pair<Mas, double>> masWithScoring;

    for (const auto& [magnetic, scoring] : magneticsWithScoring) {
        auto mas = post_process_core(magnetic, inputs);
        masWithScoring.push_back({mas, scoring});
    }

    return masWithScoring;
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
    // magneticsWithScoring = add_initial_turns_by_impedance(magneticsWithScoring, inputs);

    // Same logic as filter_available_cores_suppression_application: use a tiny
    // weight so the impedance filter is a pure gate (no score contribution).
    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, 0.001, true);
    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterMagneticInductance.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    // Manufacturability proxy — see filter_available_cores_suppression_application
    // for rationale. Same internal weight (1.0) here.
    magneticsWithScoring = filterTurnCount.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    if (magneticsWithScoring.size() > maximumNumberResults) {
        if (get_unique_core_shapes()) {
            magneticsWithScoring = cull_to_unique_core_shapes(magneticsWithScoring, maximumNumberResults);
        }
        else {
            magneticsWithScoring.resize(maximumNumberResults); // F10 FIX: resize instead of copy-construct
        }
    }

    correct_windings(&magneticsWithScoring, inputs);
    add_alternative_materials(&magneticsWithScoring, inputs);

    std::vector<std::pair<Mas, double>> masWithScoring;

    for (const auto& [magnetic, scoring] : magneticsWithScoring) {
        auto mas = post_process_core(magnetic, inputs);
        masWithScoring.push_back({mas, scoring});
    }

    return masWithScoring;
}

} // namespace OpenMagnetics
