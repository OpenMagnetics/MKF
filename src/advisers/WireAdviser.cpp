#include "advisers/MagneticFilter.h"
#include "advisers/WireAdviser.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Coil.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include <list>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <magic_enum.hpp>
#include "support/Exceptions.h"
#include "support/Logger.h"


namespace OpenMagnetics {

// Phase 1 tie-break: when two candidates end up with strictly-equal final
// scores (e.g. "Round 16.5 - Single Build" vs "Round 16.5 - Heavy Build"
// have identical conductor geometry and so identical resistance/skin/
// proximity contributions), break the tie by preferring the wire with the
// SMALLER outer dimension. That biases toward thinner insulation /
// unserved litz, which is the right default (better packing factor, less
// material cost). Strictly stable for non-tied pairs.
static double wire_outer_metric(const Winding& winding) {
    auto wire = OpenMagnetics::Coil::resolve_wire(winding);
    if (wire.get_outer_diameter()) {
        return resolve_dimensional_values(wire.get_outer_diameter().value());
    }
    if (wire.get_outer_width() && wire.get_outer_height()) {
        return resolve_dimensional_values(wire.get_outer_width().value())
             + resolve_dimensional_values(wire.get_outer_height().value());
    }
    if (wire.get_outer_width()) {
        return resolve_dimensional_values(wire.get_outer_width().value());
    }
    if (wire.get_outer_height()) {
        return resolve_dimensional_values(wire.get_outer_height().value());
    }
    // No outer dimension known — return +inf so this candidate is deprioritised
    // among ties (i.e., a candidate WITH a known outer dim wins).
    return std::numeric_limits<double>::infinity();
}

// Reorders the candidates by a permutation. A candidate is expensive to move: Winding carries its
// Wire twice, and the generated MAS classes declare a virtual destructor, so they have no move
// operations and every "move" is a deep copy. Sorting the candidates themselves cost n*log(n) such
// copies per filter; sorting their indices and placing each candidate once costs n.
static void apply_order(std::vector<std::pair<Winding, double>>* coilsWithScoring, const std::vector<size_t>& order) {
    bool identity = true;
    for (size_t i = 0; i < order.size(); ++i) {
        if (order[i] != i) {
            identity = false;
            break;
        }
    }
    if (identity) {
        return;
    }
    std::vector<std::pair<Winding, double>> reordered;
    reordered.reserve(order.size());
    for (auto index : order) {
        reordered.push_back(std::move((*coilsWithScoring)[index]));
    }
    *coilsWithScoring = std::move(reordered);
}

static void break_score_ties(std::vector<std::pair<Winding, double>>* coilsWithScoring) {
    if (coilsWithScoring->size() < 2) return;
    // Same ordering as a stable_sort of the candidates by (score desc, outer metric asc), with the
    // metric computed once per candidate instead of twice per comparison.
    std::vector<double> outerMetrics;
    outerMetrics.reserve(coilsWithScoring->size());
    for (const auto& candidate : *coilsWithScoring) {
        outerMetrics.push_back(wire_outer_metric(candidate.first));
    }
    std::vector<size_t> order(coilsWithScoring->size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const double scoreA = (*coilsWithScoring)[a].second;
        const double scoreB = (*coilsWithScoring)[b].second;
        if (scoreA != scoreB) return scoreA > scoreB;
        return outerMetrics[a] < outerMetrics[b];
    });
    apply_order(coilsWithScoring, order);
}

void normalize_scoring(std::vector<std::pair<Winding, double>>* coilsWithScoring, std::vector<double>* newScoring, bool invert=true, const std::string& filterName="") {
    // Debug: Check for NaN values before normalization
    for (size_t i = 0; i < newScoring->size(); ++i) {
        if (std::isnan((*newScoring)[i]) || std::isinf((*newScoring)[i])) {
            auto wire = OpenMagnetics::Coil::resolve_wire((*coilsWithScoring)[i].first);
            std::string wireInfo = "Wire type: " + std::string(magic_enum::enum_name(wire.get_type()));
            if (wire.get_name()) {
                wireInfo += ", name: " + wire.get_name().value();
            }
            if (wire.get_conducting_width()) {
                wireInfo += ", conducting_width: " + std::to_string(OpenMagnetics::resolve_dimensional_values(wire.get_conducting_width().value()));
            }
            else {
                wireInfo += ", conducting_width: MISSING";
            }
            if (wire.get_conducting_height()) {
                wireInfo += ", conducting_height: " + std::to_string(OpenMagnetics::resolve_dimensional_values(wire.get_conducting_height().value()));
            }
            else {
                wireInfo += ", conducting_height: MISSING";
            }
            throw std::invalid_argument("NaN/Inf scoring detected in filter '" + filterName + "' at index " + std::to_string(i) + 
                                       ". Scoring value: " + std::to_string((*newScoring)[i]) + 
                                       ". " + wireInfo);
        }
    }
    
    auto normalizedScorings = OpenMagnetics::normalize_scoring(*newScoring, 1, invert, false);

    for (size_t i = 0; i < (*coilsWithScoring).size(); ++i) {
        (*coilsWithScoring)[i].second += normalizedScorings[i];
    }
    // F12 FIX: stable ordering for reproducible results (sorted by index, see apply_order).
    std::vector<size_t> order((*coilsWithScoring).size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t b1, size_t b2) {
        return (*coilsWithScoring)[b1].second > (*coilsWithScoring)[b2].second;
    });
    apply_order(coilsWithScoring, order);
}

std::vector<std::pair<Winding, double>>  WireAdviser::filter_by_area_no_parallels(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                    Section section) {
    std::vector<std::pair<Winding, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;

    auto filter = MagneticFilterAreaNoParallels(_maximumNumberParallels);

    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        auto [valid, scoring] = filter.evaluate_magnetic((*unfilteredCoils)[coilIndex].first, section);

        if (valid) {
            newScoring.push_back(scoring);
        }
        else {
            listOfIndexesToErase.push_back(coilIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredCoils).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoilsWithScoring.push_back((*unfilteredCoils)[i]);
        }
    }
    // (*unfilteredCoils).clear();

    if (filteredCoilsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, true, "area_no_parallels");
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<Winding, double>>  WireAdviser::filter_by_area_with_parallels(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                    Section section,
                                                                                                    double numberSections,
                                                                                                    bool allowNotFit) {
    std::vector<std::pair<Winding, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    double sectionArea;
    if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        sectionArea = section.get_dimensions()[0] * section.get_dimensions()[1];
    }
    else {
        sectionArea = std::numbers::pi * pow(section.get_dimensions()[0], 2) * section.get_dimensions()[1] / 360;
    }

    auto filter = MagneticFilterAreaWithParallels();

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        auto [valid, scoring] = filter.evaluate_magnetic((*unfilteredCoils)[coilIndex].first, section, numberSections, sectionArea, allowNotFit);

        if (valid) {
            newScoring.push_back(scoring);
        }
        else {
            listOfIndexesToErase.push_back(coilIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredCoils).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoilsWithScoring.push_back((*unfilteredCoils)[i]);
        }
    }
    // (*unfilteredCoils).clear();

    if (filteredCoilsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, false, "area_with_parallels");
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<Winding, double>> WireAdviser::filter_by_effective_resistance(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                      SignalDescriptor current,
                                                                                                      double temperature) {
    std::vector<std::pair<Winding, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    if (!current.get_processed()->get_effective_frequency()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current processed is missing field effective frequency");
    }
    auto currentEffectiveFrequency = current.get_processed()->get_effective_frequency().value();

    auto filter = MagneticFilterEffectiveResistance();

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        auto [valid, scoring] = filter.evaluate_magnetic((*unfilteredCoils)[coilIndex].first, currentEffectiveFrequency, temperature);

        if (valid) {
            newScoring.push_back(scoring);
        }
        else {
            listOfIndexesToErase.push_back(coilIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredCoils).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoilsWithScoring.push_back((*unfilteredCoils)[i]);
        }
    }
    // (*unfilteredCoils).clear();

    if (filteredCoilsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, true, "effective_resistance");
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<Winding, double>> WireAdviser::filter_by_skin_losses_density(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                     SignalDescriptor current,
                                                                                                     double temperature) {
    std::vector<std::pair<Winding, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    auto filter = MagneticFilterSkinLossesDensity();

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        auto [valid, scoring] = filter.evaluate_magnetic((*unfilteredCoils)[coilIndex].first, current, temperature);

        if (valid) {
            newScoring.push_back(scoring);
        }
        else {
            listOfIndexesToErase.push_back(coilIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredCoils).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoilsWithScoring.push_back((*unfilteredCoils)[i]);
        }
    }
    // (*unfilteredCoils).clear();

    if (filteredCoilsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, true, "skin_losses_density");
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<Winding, double>> WireAdviser::filter_by_proximity_factor(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                      SignalDescriptor current,
                                                                                                      double temperature) {
    std::vector<std::pair<Winding, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    if (!current.get_processed()->get_effective_frequency()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current processed is missing field effective frequency");
    }
    auto currentEffectiveFrequency = current.get_processed()->get_effective_frequency().value();

    auto filter = MagneticFilterProximityFactor();

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        auto [valid, scoring] = filter.evaluate_magnetic((*unfilteredCoils)[coilIndex].first, currentEffectiveFrequency, temperature);

        if (valid) {
            newScoring.push_back(scoring);
        }
        else {
            listOfIndexesToErase.push_back(coilIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredCoils).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoilsWithScoring.push_back((*unfilteredCoils)[i]);
        }
    }
    // (*unfilteredCoils).clear();

    if (filteredCoilsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, true, "proximity_factor");
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<Winding, double>> WireAdviser::filter_by_layer_parity(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                            Section section) {
    std::vector<std::pair<Winding, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    auto filter = MagneticFilterLayerParity();

    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        // The parity filter never rejects a candidate: an odd layer count is manufacturable,
        // just worse, so every candidate survives and only the scoring moves.
        auto [valid, scoring] = filter.evaluate_magnetic((*unfilteredCoils)[coilIndex].first, section);
        if (!valid) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "MagneticFilterLayerParity must never reject a candidate");
        }
        newScoring.push_back(scoring);
        filteredCoilsWithScoring.push_back((*unfilteredCoils)[coilIndex]);
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, true, "layer_parity");
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<Winding, double>> WireAdviser::filter_by_solid_insulation_requirements(std::vector<std::pair<Winding, double>>* unfilteredCoils, WireSolidInsulationRequirements wireSolidInsulationRequirements) {
    std::vector<std::pair<Winding, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    auto filter = MagneticFilterSolidInsulationRequirements();

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        auto [valid, scoring] = filter.evaluate_magnetic((*unfilteredCoils)[coilIndex].first, wireSolidInsulationRequirements);

        if (valid) {
            newScoring.push_back(scoring);
        }
        else {
            listOfIndexesToErase.push_back(coilIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredCoils).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredCoilsWithScoring.push_back((*unfilteredCoils)[i]);
        }
    }
    // (*unfilteredCoils).clear();

    if (filteredCoilsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, true, "solid_insulation_requirements");
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<OpenMagnetics::Winding, double>> WireAdviser::get_advised_wire(OpenMagnetics ::Winding winding,
                                                                                        Section section,
                                                                                        SignalDescriptor current,
                                                                                        double temperature,
                                                                                        uint8_t numberSections,
                                                                                        size_t maximumNumberResults){
    std::vector<Wire> wires;
    auto& settings = OpenMagnetics::Settings::GetInstance();

    // See CoilAdviser::get_advised_coil: never lazily refill wireDatabase
    // while a LibraryContext scope is active — the scoped contents ARE the
    // library, even when empty.
    if (wireDatabase.empty() && !LibraryContext::Scope::anyActive()) {
        load_wires();
    }

    for (auto [name, wire] : wireDatabase) {
        if ((settings.get_wire_adviser_include_foil() || wire.get_type() != WireType::FOIL) &&
            (settings.get_wire_adviser_include_planar() ||  wire.get_type() != WireType::PLANAR) &&
            ((settings.get_wire_adviser_include_rectangular() && (settings.get_wire_adviser_allow_rectangular_in_toroidal_cores() || section.get_coordinate_system() == CoordinateSystem::CARTESIAN)) || wire.get_type() != WireType::RECTANGULAR) &&
            (settings.get_wire_adviser_include_litz() || wire.get_type() != WireType::LITZ) &&
            (settings.get_wire_adviser_include_round() || wire.get_type() != WireType::ROUND)) {

            if (!_commonWireStandard || !wire.get_standard()) {
                wires.push_back(wire);
            }
            else if (wire.get_standard().value() == _commonWireStandard){
                wires.push_back(wire);
            }
        }
    }
    return get_advised_wire(&wires, winding, section, current, temperature, numberSections, maximumNumberResults);
}

std::vector<std::pair<Winding, double>> WireAdviser::get_advised_wire(
        Winding winding,
        Section section,
        SignalDescriptor current,
        double temperature,
        uint8_t numberSections,
        size_t maximumNumberResults,
        const LibraryContext* ctx,
        const AdviserConstraints& constraints) {
    auto scope = ctx ? ctx->applyScoped() : LibraryContext::Scope{};

    if (constraints.wireType.empty()) {
        return get_advised_wire(winding, section, current, temperature,
                                numberSections, maximumNumberResults);
    }

    if (wireDatabase.empty() && !LibraryContext::Scope::anyActive()) load_wires();

    std::vector<Wire> wires;
    for (auto& [name, wire] : wireDatabase) {
        if (!acceptsWireType(constraints.wireType, wire.get_type())) continue;
        if (_commonWireStandard && wire.get_standard()
            && wire.get_standard().value() != _commonWireStandard) {
            continue;
        }
        wires.push_back(wire);
    }
    return get_advised_wire(&wires, winding, section, current, temperature,
                            numberSections, maximumNumberResults);
}

std::vector<std::pair<Winding, double>> WireAdviser::create_planar_dataset(Winding winding,
                                                                                       Section section,
                                                                                       SignalDescriptor current,
                                                                                       double temperature,
                                                                                       uint8_t numberSections) {
    std::vector<std::pair<Winding, double>> windings;
    auto planarWires = get_wires(WireType::PLANAR);

    // WA-LOGIC-1 NOTE: Parallels are tried first, biasing toward complex designs.
    // Consider generating no-parallels first for simpler/cheaper solutions,
    // or let scoring handle the preference between simple and parallel designs.
    // Parallels (FIRST - tried first in CoilAdviser)
    {
        auto maximumNumberTurnsPerSection = winding.get_number_turns();
        auto maximumAvailableWidthForCopper = section.get_dimensions()[0] - 2 * get_border_to_wire_distance() - (maximumNumberTurnsPerSection - 1) * get_wire_to_wire_distance();
        if (maximumAvailableWidthForCopper < 0) {
            return windings;
        }
        auto maximumAvailableWidthForTurn = maximumAvailableWidthForCopper / winding.get_number_turns();
        size_t maximumNumberParallels = _maximumNumberParallels;

        for (auto wire : planarWires) {
            // Capture catalogue display name BEFORE any per-variant mutation,
            // so each (parallels) variant gets a unique name suffix (see fix
            // for "Planar X µm appears twice in top-N" bug).
            std::string baseName = wire.get_name().value_or("Planar");
            double conductingHeight = resolve_dimensional_values(wire.get_conducting_height().value());
            if (conductingHeight < section.get_dimensions()[1]) {
                wire.set_nominal_value_outer_height(conductingHeight);
                // Set conducting_width FIRST before calculating parallels needed
                wire.set_nominal_value_conducting_width(maximumAvailableWidthForTurn);
                int minParallelsNeeded = Wire::calculate_number_parallels_needed(current, temperature, wire, _maximumEffectiveCurrentDensity);
                size_t startParallels = std::max(size_t(2), static_cast<size_t>(minParallelsNeeded));
                for (size_t numberParallels = startParallels; numberParallels <= maximumNumberParallels; ++numberParallels) {
                    wire.set_nominal_value_conducting_width(maximumAvailableWidthForTurn);
                    wire.set_nominal_value_outer_width(maximumAvailableWidthForTurn);
                    wire.set_nominal_value_conducting_area(maximumAvailableWidthForTurn * conductingHeight);
                    // Phase 1 fix: same planar catalogue entry was emitted at
                    // multiple parallel counts with identical display names
                    // and DIFFERENT downstream scores, so the top-N showed
                    // e.g. "Planar 243.59 µm" twice with no way to tell them
                    // apart. Encode parallels into the display name.
                    wire.set_name(baseName + " x" + std::to_string(numberParallels) + " parallels");
                    winding.set_wire(wire);
                    winding.set_number_parallels(numberParallels);

                    // No copper-thickness penalty: heavier copper is not ranked below its electrical merit.
                    windings.push_back(std::pair<Winding, double>{winding, 0});
                }
            }
        }
    }

    // No parallels (SECOND - only add if single parallel meets current density requirements)
    {
        auto maximumNumberTurnsPerSection = ceil(double(winding.get_number_turns()) / numberSections);
        auto maximumAvailableWidthForCopper = section.get_dimensions()[0] - 2 * get_border_to_wire_distance() - (maximumNumberTurnsPerSection - 1) * get_wire_to_wire_distance();
        if (maximumAvailableWidthForCopper < 0) {
            return windings;
        }
        auto maximumAvailableWidthForTurn = maximumAvailableWidthForCopper / maximumNumberTurnsPerSection;

        for (auto wire : planarWires) {
            double conductingHeight = resolve_dimensional_values(wire.get_conducting_height().value());
            if (conductingHeight < section.get_dimensions()[1]) {
                wire.set_nominal_value_outer_height(conductingHeight);
                // Set conducting_width FIRST before calculating parallels needed
                wire.set_nominal_value_conducting_width(maximumAvailableWidthForTurn);
                int minParallelsNeeded = Wire::calculate_number_parallels_needed(current, temperature, wire, _maximumEffectiveCurrentDensity);
                if (minParallelsNeeded <= 1) {
                    wire.set_nominal_value_conducting_width(maximumAvailableWidthForTurn);
                    wire.set_nominal_value_outer_width(maximumAvailableWidthForTurn);
                    wire.set_nominal_value_conducting_area(maximumAvailableWidthForTurn * conductingHeight);
                    winding.set_wire(wire);
                    winding.set_number_parallels(1);
                    
                    // No copper-thickness penalty: heavier copper is not ranked below its electrical merit.
                    windings.push_back(std::pair<Winding, double>{winding, 0});
                }
            }
        }
    }
    return windings;
}

// ABT #5: synthesize fine-strand litz candidates per Sullivan (strand << skin
// depth) so the WireAdviser pool contains low-proximity options for HF AC
// designs. The existing proximity/skin-aware scoring already PREFERS such litz
// when it is present — the bug was that no fine-enough litz was ever in the pool
// (catalog litz is too coarse or doesn't fit small cores). Each synthesized
// bundle carries ~the per-bundle copper; create_dataset's parallels logic sizes
// the rest. Cheap: a few wire constructions, no extra simulation. Gated to the
// regime where a solid conductor would exceed the skin depth (AC losses matter),
// so DC/LF designs are unaffected. The exact Sullivan closed-form optimum is not
// hard-coded; the strand grid spans delta/2..delta/4 and the loss-aware scoring
// (and #4's core re-rank simulate) select the empirical optimum.
static std::vector<Wire> synthesize_litz_candidates(const SignalDescriptor& current,
                                                    double temperature,
                                                    double maximumEffectiveCurrentDensity) {
    std::vector<Wire> out;
    auto& settings = Settings::GetInstance();
    if (!settings.get_wire_adviser_include_litz()) return out;
    if (!current.get_processed()) return out;
    auto processed = current.get_processed().value();
    if (!processed.get_effective_frequency() || !processed.get_rms()) return out;
    double effectiveFrequency = processed.get_effective_frequency().value();
    double rms = processed.get_rms().value();
    if (effectiveFrequency <= 0 || rms <= 0 || maximumEffectiveCurrentDensity <= 0) return out;

    // Skin depth in copper at the effective frequency (lightly temp-corrected).
    const double rhoCu = 1.724e-8 * (1.0 + 0.00393 * (temperature - 20.0));
    const double mu0 = 4e-7 * std::numbers::pi;
    double skinDepth = std::sqrt(rhoCu / (std::numbers::pi * effectiveFrequency * mu0));

    // Copper one bundle should carry (parallels multiply to meet the current).
    double targetCopperArea = rms / maximumEffectiveCurrentDensity;
    if (targetCopperArea <= 0) return out;

    // Gate: only when a solid conductor of that area would exceed the skin depth
    // (skin/proximity actually matter). Below that, solid round is already fine.
    double solidDiameter = std::sqrt(4.0 * targetCopperArea / std::numbers::pi);
    if (solidDiameter <= skinDepth) return out;

    // Sullivan-seeded grid: strand diameter = skinDepth/k (k=2..4), strand count
    // chosen to fill the bundle copper area. Strand << skin depth -> low proximity.
    for (int k : {2, 3, 4}) {
        double strandDiameter = skinDepth / k;
        if (strandDiameter < 25e-6) continue;                 // impractically fine
        double strandArea = std::numbers::pi * strandDiameter * strandDiameter / 4.0;
        int64_t numberStrands = std::max<int64_t>(3, static_cast<int64_t>(std::ceil(targetCopperArea / strandArea)));
        if (numberStrands > 2000) continue;                   // unbuildable bundle
        try {
            auto litz = Wire::create_quick_litz_wire(strandDiameter, numberStrands);
            // Name it so MAS/frontend show a meaningful label (create_quick_litz_wire
            // leaves the name unset). e.g. "Synthesized Litz 76x0.070mm".
            std::ostringstream name;
            name << "Synthesized Litz " << numberStrands << "x"
                 << std::fixed << std::setprecision(3) << (strandDiameter * 1e3) << "mm";
            litz.set_name(name.str());
            out.push_back(litz);
        }
        catch (const std::exception&) {
            // create_quick_litz_wire may reject a non-standard strand size; skip it.
        }
    }
    return out;
}

std::vector<std::pair<Winding, double>> WireAdviser::create_dataset(Winding winding,
                                                                                      std::vector<Wire>* wires,
                                                                                      Section section,
                                                                                      SignalDescriptor current,
                                                                                      double temperature){
    auto& settings = Settings::GetInstance();
    std::vector<std::pair<Winding, double>> windings;

    // The first thing get_advised_wire does with this dataset is filter_by_area_no_parallels: a
    // per-wire geometric verdict (does one wire fit the section at all) under which every survivor
    // scores 0. Applying that same verdict here, before a Winding is built, keeps exactly the same
    // candidates in the same order, and stops every call from deep-copying the whole catalogue: the
    // coil adviser calls this once per winding, pattern, insulation combination and wire
    // configuration of every core it tries (~1000 calls of ~8000 wires for a Magnetic Adviser run),
    // and on a small core nearly all of those copies were thrown away by that first filter.
    auto areaNoParallelsFilter = MagneticFilterAreaNoParallels(_maximumNumberParallels);
    // Likewise filter_by_solid_insulation_requirements, which runs next: its verdict is per wire, and
    // its scores are normalised over the wires it keeps, which are the same wires either way.
    auto solidInsulationFilter = MagneticFilterSolidInsulationRequirements();

    auto add_candidates = [&](Wire& catalogueWire) {
        if ((!settings.get_wire_adviser_include_foil() && catalogueWire.get_type() == WireType::FOIL) ||
            (!settings.get_wire_adviser_include_planar() &&  catalogueWire.get_type() == WireType::PLANAR) ||
            (!(settings.get_wire_adviser_include_rectangular() && (settings.get_wire_adviser_allow_rectangular_in_toroidal_cores() || section.get_coordinate_system() == CoordinateSystem::CARTESIAN)) && catalogueWire.get_type() == WireType::RECTANGULAR) ||
            (!settings.get_wire_adviser_include_litz() && catalogueWire.get_type() == WireType::LITZ) ||
            (!settings.get_wire_adviser_include_round() && catalogueWire.get_type() == WireType::ROUND)) {
            return;
        }
        // Round and rectangular wires are judged as they are in the catalogue; litz, foil and planar
        // first get the same preparation as before (resolved strand, cut to the section), because
        // their outer dimensions depend on it.
        const bool needsPreparation = catalogueWire.get_type() == WireType::LITZ ||
                                      catalogueWire.get_type() == WireType::FOIL ||
                                      catalogueWire.get_type() == WireType::PLANAR;
        if (!needsPreparation) {
            if (!areaNoParallelsFilter.wire_fits(catalogueWire, 1, winding.get_number_turns(), section)) {
                return;
            }
            if (_wireSolidInsulationRequirements && !solidInsulationFilter.evaluate_wire(catalogueWire, _wireSolidInsulationRequirements.value()).first) {
                return;
            }
        }
        Wire wire = catalogueWire;
        if (wire.get_type() == WireType::LITZ) {
            wire.set_strand(wire.resolve_strand());
        }
        if (wire.get_type() == WireType::FOIL) {
            wire.cut_foil_wire_to_section(section);
        }
        if (wire.get_type() == WireType::PLANAR) {
            wire.cut_planar_wire_to_section(section);
        }

        int numberParallelsNeeded;
        if (wire.get_type() == WireType::RECTANGULAR) {
            numberParallelsNeeded = 1;
        }
        else {
            numberParallelsNeeded = Wire::calculate_number_parallels_needed(current, temperature, wire, _maximumEffectiveCurrentDensity);
            if (numberParallelsNeeded > _maximumNumberParallels) {
                return;
            }
        }

        const bool fits = areaNoParallelsFilter.wire_fits(wire, numberParallelsNeeded, winding.get_number_turns(), section);
        const bool addAnotherParallel = numberParallelsNeeded < _maximumNumberParallels;
        // Only a foil's verdict depends on the number of parallels.
        const bool fitsWithAnotherParallel = addAnotherParallel &&
            areaNoParallelsFilter.wire_fits(wire, numberParallelsNeeded + 1, winding.get_number_turns(), section);
        if (!fits && !fitsWithAnotherParallel) {
            return;
        }
        if (needsPreparation && _wireSolidInsulationRequirements && !solidInsulationFilter.evaluate_wire(wire, _wireSolidInsulationRequirements.value()).first) {
            return;
        }
        winding.set_wire(std::move(wire));
        if (fits) {
            winding.set_number_parallels(numberParallelsNeeded);
            windings.push_back(std::pair<Winding, double>{winding, 0});
        }
        if (fitsWithAnotherParallel) {
            winding.set_number_parallels(numberParallelsNeeded + 1);
            windings.push_back(std::pair<Winding, double>{winding, 0});
        }
    };

    for (auto& wire : *wires) {
        add_candidates(wire);
    }
    // Extend the candidate pool with synthesized fine-strand litz (ABT #5), after the catalogue as
    // before; the caller's wire list is not mutated.
    if (_synthesizeLitz) {
        auto synthesizedLitz = synthesize_litz_candidates(current, temperature, _maximumEffectiveCurrentDensity);
        for (auto& wire : synthesizedLitz) {
            add_candidates(wire);
        }
    }

    return windings;
}

void WireAdviser::set_maximum_area_proportion(std::vector<std::pair<Winding, double>>* unfilteredCoils, Section section, uint8_t numberSections) {
    double sectionArea;

    if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        sectionArea = section.get_dimensions()[0] * section.get_dimensions()[1];
    }
    else {
        sectionArea = std::numbers::pi * pow(section.get_dimensions()[0], 2) * section.get_dimensions()[1] / 360;
    }

    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        auto wire = Coil::resolve_wire((*unfilteredCoils)[coilIndex].first);
        if (!Coil::resolve_wire((*unfilteredCoils)[coilIndex].first).get_conducting_area()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Conducting area is missing");
        }
        auto neededOuterAreaNoCompact = wire.get_maximum_outer_width() * wire.get_maximum_outer_height();

        neededOuterAreaNoCompact *= (*unfilteredCoils)[coilIndex].first.get_number_parallels() * (*unfilteredCoils)[coilIndex].first.get_number_turns() / static_cast<double>(numberSections);

        double areaProportion = neededOuterAreaNoCompact / sectionArea;
        _maximumOuterAreaProportion = std::max(_maximumOuterAreaProportion, areaProportion);

        // if (areaProportion > 1) {
        //     throw std::runtime_error("areaProportion cannot be bigger than 1");
        // }
    }
}

std::vector<std::pair<Winding, double>> WireAdviser::get_advised_planar_wire(Winding winding,
                                                                                               Section section,
                                                                                               SignalDescriptor current,
                                                                                               double temperature,
                                                                                               uint8_t numberSections,
                                                                                               size_t maximumNumberResults) {
    auto coilsWithScoring = create_planar_dataset(winding, section, current, temperature, numberSections);

    logEntry("We start the search with " + std::to_string(coilsWithScoring.size()) + " wires");

    coilsWithScoring = filter_by_effective_resistance(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " planar wires after filtering by effective resistance.");

    // WA-BUG-1 FIX: guard with harmonics check (matches get_advised_wire)
    if (current.get_harmonics()) {
        coilsWithScoring = filter_by_skin_losses_density(&coilsWithScoring, current, temperature);
        logEntry("There are " + std::to_string(coilsWithScoring.size()) + " planar wires after filtering by skin losses density.");
    }

    coilsWithScoring = filter_by_proximity_factor(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " planar wires after filtering by proximity factor.");

    break_score_ties(&coilsWithScoring);

    if (coilsWithScoring.size() > maximumNumberResults) {
        auto finalCoilsWithScoring = std::vector<std::pair<Winding, double>>(coilsWithScoring.begin(), coilsWithScoring.end() - (coilsWithScoring.size() - maximumNumberResults));
        set_maximum_area_proportion(&finalCoilsWithScoring, section, numberSections);
        return finalCoilsWithScoring;
    }
    else {
        set_maximum_area_proportion(&coilsWithScoring, section, numberSections);
    }


    return coilsWithScoring;
}

std::vector<std::pair<Winding, double>> WireAdviser::get_advised_wire(std::vector<Wire>* wires,
                                                                                        Winding winding,
                                                                                        Section section,
                                                                                        SignalDescriptor current,
                                                                                        double temperature,
                                                                                        uint8_t numberSections,
                                                                                        size_t maximumNumberResults){
    // WireAdviser instances are long-lived (CoilAdviser holds one across cores);
    // without this reset the reported maximum accumulates across advise calls.
    _maximumOuterAreaProportion = 0;
    auto coilsWithScoring = create_dataset(winding, wires, section, current, temperature);


    logEntry("We start the search with " + std::to_string(coilsWithScoring.size()) + " wires");
    coilsWithScoring = filter_by_area_no_parallels(&coilsWithScoring, section);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by area no parallels.");

    if (_wireSolidInsulationRequirements) {
        coilsWithScoring = filter_by_solid_insulation_requirements(&coilsWithScoring, _wireSolidInsulationRequirements.value());
        logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by solid insulation.");
    }

    auto tempCoilsWithScoring = filter_by_area_with_parallels(&coilsWithScoring, section, numberSections, false);
    logEntry("There are " + std::to_string(tempCoilsWithScoring.size()) + " after filtering by area with parallels.");

    if (tempCoilsWithScoring.size() == 0) {
        coilsWithScoring = filter_by_area_with_parallels(&coilsWithScoring, section, numberSections, true);
        logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by area with parallels, allowing not fitting.");
    }
    else{
        coilsWithScoring = tempCoilsWithScoring;
    }

    coilsWithScoring = filter_by_effective_resistance(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by effective resistance.");

    // Skin losses density filter requires harmonics data
    if (current.get_harmonics()) {
        coilsWithScoring = filter_by_skin_losses_density(&coilsWithScoring, current, temperature);
        logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by skin losses density.");
    }

    coilsWithScoring = filter_by_proximity_factor(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by proximity factor.");

    // ABT #1177 (WP8, DFM rule R1): opt-in, default off. Prefers a wire that fills the section
    // in an even number of layers, so the winding does not end on a drag-back.
    if (Settings::GetInstance().get_wire_adviser_penalize_odd_layer_count()) {
        coilsWithScoring = filter_by_layer_parity(&coilsWithScoring, section);
        logEntry("Scored " + std::to_string(coilsWithScoring.size()) + " wires by layer parity (DFM rule R1).");
    }

    break_score_ties(&coilsWithScoring);

    if (coilsWithScoring.size() > maximumNumberResults) {
        auto finalCoilsWithScoring = std::vector<std::pair<Winding, double>>(coilsWithScoring.begin(), coilsWithScoring.end() - (coilsWithScoring.size() - maximumNumberResults));
        set_maximum_area_proportion(&finalCoilsWithScoring, section, numberSections);
        return finalCoilsWithScoring;
    }
    else {
        set_maximum_area_proportion(&coilsWithScoring, section, numberSections);
    }


    return coilsWithScoring;
}
} // namespace OpenMagnetics
