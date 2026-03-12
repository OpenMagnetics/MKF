#include "advisers/MagneticFilter.h"
#include "advisers/WireAdviser.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Coil.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include <list>
#include <magic_enum.hpp>
#include "support/Exceptions.h"
#include "support/Logger.h"


namespace OpenMagnetics {

// F41 NOTE: This is one of 4 separate normalize_scoring implementations in the codebase.
// They should be unified into a single templated utility function.
// This version (WireAdviser) is the most complete with B7+B8+O24 fixes.
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
    std::stable_sort((*coilsWithScoring).begin(), (*coilsWithScoring).end(), [](const std::pair<Winding, double>& b1, const std::pair<Winding, double>& b2) {
        return b1.second > b2.second;
    }); // F12 FIX: stable_sort for reproducible results
}

std::vector<std::pair<Winding, double>>  WireAdviser::filter_by_area_no_parallels(std::vector<std::pair<Winding, double>>* unfilteredCoils,
                                                                                                    Section section) {
    std::vector<std::pair<Winding, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;

    auto filter = MagneticFilterAreaNoParallels(_maximumNumberParallels);

    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){
        auto winding = (*unfilteredCoils)[coilIndex].first;
        auto [valid, scoring] = filter.evaluate_magnetic(winding, section);

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

    if (wireDatabase.empty()) {
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

/**
 * @brief Calculate copper thickness penalty for planar wires (progressive).
 *
 * Penalizes thicker copper (higher oz) progressively based on real PCB industry
 * cost and manufacturing factors:
 * 
 * 1. Material cost: Scales with thickness
 * 2. Etching time: Heavy copper requires longer etch cycles
 * 3. Feature resolution: Thick copper has worse minimum trace/space
 * 4. Lead time: Heavy copper boards often have longer lead times
 * 5. Supplier availability: Fewer fabs handle 3+ oz, very few handle 6+ oz
 *
 * Industry standard copper weights and approximate cost multipliers:
 * - 0.5 oz =  17.5 µm : 1.0x (baseline)
 * - 1 oz   =  35 µm   : 1.0x (standard, still baseline)
 * - 2 oz   =  70 µm   : 1.3-1.5x cost
 * - 3 oz   = 105 µm   : 1.8-2.2x cost (heavy copper threshold)
 * - 4 oz   = 140 µm   : 2.5-3.0x cost
 * - 6 oz   = 210 µm   : 4-5x cost (specialized fabs only)
 * - 10 oz  = 350 µm   : 8-10x cost (very specialized, few suppliers)
 *
 * The penalty does NOT prevent selection of thick copper when required for
 * current density - it only makes thinner copper preferred when multiple
 * options would satisfy the electrical requirements.
 *
 * @param conductingHeight Copper thickness in meters
 * @return Penalty value (0 = best, higher = worse). Range 0 to ~1.0.
 */
double calculate_planar_thickness_penalty(double conductingHeight) {
    // Standard copper thicknesses (in meters)
    // 1 oz = 35 µm = 1.37 mil
    const double oz_1   =  35.0e-6;   // 1 oz (baseline - no penalty)
    const double oz_2   =  70.0e-6;   // 2 oz
    const double oz_3   = 105.0e-6;   // 3 oz (heavy copper threshold)
    const double oz_4   = 140.0e-6;   // 4 oz
    const double oz_6   = 210.0e-6;   // 6 oz (specialized)
    const double oz_10  = 350.0e-6;   // 10 oz (very specialized)
    
    // No penalty for standard 0.5-1 oz copper
    if (conductingHeight <= oz_1) {
        return 0.0;
    }
    
    // Progressive penalty based on PCB industry cost scaling
    // Penalty values chosen to reflect relative cost increases
    double penalty;
    
    if (conductingHeight <= oz_2) {
        // 1-2 oz: mild penalty (linear interpolation)
        // 2 oz is ~1.3-1.5x cost -> penalty ~0.08
        double ratio = (conductingHeight - oz_1) / (oz_2 - oz_1);
        penalty = 0.08 * ratio;
    }
    else if (conductingHeight <= oz_3) {
        // 2-3 oz: moderate penalty (entering heavy copper territory)
        // 3 oz is ~1.8-2.2x cost -> penalty ~0.18
        double ratio = (conductingHeight - oz_2) / (oz_3 - oz_2);
        penalty = 0.08 + 0.10 * ratio;
    }
    else if (conductingHeight <= oz_4) {
        // 3-4 oz: significant penalty (heavy copper)
        // 4 oz is ~2.5-3x cost -> penalty ~0.30
        double ratio = (conductingHeight - oz_3) / (oz_4 - oz_3);
        penalty = 0.18 + 0.12 * ratio;
    }
    else if (conductingHeight <= oz_6) {
        // 4-6 oz: steep penalty (specialized fabrication)
        // 6 oz is ~4-5x cost -> penalty ~0.50
        double ratio = (conductingHeight - oz_4) / (oz_6 - oz_4);
        penalty = 0.30 + 0.20 * ratio;
    }
    else {
        // 6-10+ oz: very steep penalty (very specialized, few suppliers)
        // 10 oz is ~8-10x cost -> penalty ~1.0
        double ratio = (conductingHeight - oz_6) / (oz_10 - oz_6);
        ratio = std::min(1.0, ratio);  // Cap at 10 oz equivalent
        penalty = 0.50 + 0.50 * ratio * ratio;  // Quadratic for extreme copper
    }
    
    return penalty;
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
                    winding.set_wire(wire);
                    winding.set_number_parallels(numberParallels);
                    
                    // Apply thickness penalty (negative because higher score = better)
                    double thicknessPenalty = -calculate_planar_thickness_penalty(conductingHeight);
                    windings.push_back(std::pair<Winding, double>{winding, thicknessPenalty});
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
                    
                    // Apply thickness penalty (negative because higher score = better)
                    double thicknessPenalty = -calculate_planar_thickness_penalty(conductingHeight);
                    windings.push_back(std::pair<Winding, double>{winding, thicknessPenalty});
                }
            }
        }
    }
    return windings;
}

std::vector<std::pair<Winding, double>> WireAdviser::create_dataset(Winding winding,
                                                                                      std::vector<Wire>* wires,
                                                                                      Section section,
                                                                                      SignalDescriptor current,
                                                                                      double temperature){
    auto& settings = Settings::GetInstance();
    std::vector<std::pair<Winding, double>> windings;

    for (auto& wire : *wires){
        if (wire.get_type() == WireType::LITZ) {
            wire.set_strand(wire.resolve_strand());
        }
    }

    for (auto& wire : *wires){
        if ((!settings.get_wire_adviser_include_foil() && wire.get_type() == WireType::FOIL) ||
            (!settings.get_wire_adviser_include_planar() &&  wire.get_type() == WireType::PLANAR) ||
            (!(settings.get_wire_adviser_include_rectangular() && (settings.get_wire_adviser_allow_rectangular_in_toroidal_cores() || section.get_coordinate_system() == CoordinateSystem::CARTESIAN)) && wire.get_type() == WireType::RECTANGULAR) ||
            (!settings.get_wire_adviser_include_litz() && wire.get_type() == WireType::LITZ) ||
            (!settings.get_wire_adviser_include_round() && wire.get_type() == WireType::ROUND)) {
            continue;
        }
        int numberParallelsNeeded;
        if (wire.get_type() == WireType::FOIL) {
            wire.cut_foil_wire_to_section(section);
        }
        if (wire.get_type() == WireType::PLANAR) {
            wire.cut_planar_wire_to_section(section);
        }


        if (wire.get_type() == WireType::RECTANGULAR) {
            numberParallelsNeeded = 1;
        }
        else {
            numberParallelsNeeded = Wire::calculate_number_parallels_needed(current, temperature, wire, _maximumEffectiveCurrentDensity);
            if (numberParallelsNeeded > _maximumNumberParallels) {
                continue;
            }
        }

        winding.set_number_parallels(numberParallelsNeeded);
        winding.set_wire(wire);
        windings.push_back(std::pair<Winding, double>{winding, 0});
        if (numberParallelsNeeded < _maximumNumberParallels) {
            winding.set_number_parallels(numberParallelsNeeded + 1);
            windings.push_back(std::pair<Winding, double>{winding, 0});
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

        neededOuterAreaNoCompact *= (*unfilteredCoils)[coilIndex].first.get_number_parallels() * (*unfilteredCoils)[coilIndex].first.get_number_turns() / numberSections;

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
