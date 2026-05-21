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
// Sibling TU of CoreAdviser.cpp; declarations live in advisers/CoreAdviser.h.

#include "advisers/CoreAdviser.h"
#include "advisers/CoreAdviserInternal.h"
#include "advisers/MagneticFilter.h"
#include "physical_models/CoreTemperature.h"
#include "support/Logger.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace OpenMagnetics {

CoreAdviser::MagneticCoreFilterAreaProduct::MagneticCoreFilterAreaProduct(Inputs inputs) {
    _filter = MagneticFilterAreaProduct(inputs);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterAreaProduct::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::map<std::string, double> scaledMagneticFluxDensitiesPerMaterial;

    std::list<size_t> listOfIndexesToErase;

    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
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
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }
    // (*unfilteredMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }


    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterEnergyStored::MagneticCoreFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterEnergyStored(inputs, models);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterEnergyStored::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){  
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
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
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }
    // (*unfilteredMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterFringingFactor::MagneticCoreFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterFringingFactor(inputs, models);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterFringingFactor::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){  
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
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
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }
    // (*unfilteredMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterCost::MagneticCoreFilterCost(Inputs inputs) {
    _filter = MagneticFilterEstimatedCost(inputs);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterCost::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
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
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }
    // (*unfilteredMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::COST]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::COST, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterLosses::MagneticCoreFilterLosses(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterCoreDcAndSkinLosses(inputs, models); // F23 FIX: DC+skin effect losses
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterLosses::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    auto coilDelimitAndCompactOld = settings.get_coil_delimit_and_compact();
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);
        if (valid) {
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            newScoring.push_back(scoring);
        }
        else {
            listOfIndexesToErase.push_back(magneticIndex);
        }
    }

    if ((*unfilteredMagnetics).size() == listOfIndexesToErase.size()) {
        settings.set_coil_delimit_and_compact(coilDelimitAndCompactOld);
        return *unfilteredMagnetics;
    }

    for (size_t i = 0; i < (*unfilteredMagnetics).size(); ++i) {

        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }
    // (*unfilteredMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }

    settings.set_coil_delimit_and_compact(coilDelimitAndCompactOld);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterDimensions::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[magneticIndex]);
            newScoring.push_back(scoring);
        }
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of filteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::DIMENSIONS]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::DIMENSIONS, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
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
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[magneticIndex]);
            newScoring.push_back(scoring);
        }
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of filteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        // Compute inverted linear score directly to avoid the all-equal fallback (+1 ignoring
        // weight) in OpenMagnetics::normalize_scoring.  Fewer turns → higher score → ranked first.
        // Raw score = totalTurns (higher is worse); we invert explicitly: contribution = weight*(1 - (N-min)/(max-min)).
        double minTurns = *std::min_element(newScoring.begin(), newScoring.end());
        double maxTurns = *std::max_element(newScoring.begin(), newScoring.end());
        std::vector<double> normalizedScoring(filteredMagneticsWithScoring.size());
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            double contrib;
            if (maxTurns > minTurns) {
                contrib = weight * (1.0 - (newScoring[i] - minTurns) / (maxTurns - minTurns));
            } else {
                // All same turn count — contribute weight/2 so the filter is neutral
                contrib = weight * 0.5;
            }
            filteredMagneticsWithScoring[i].second += contrib;
            normalizedScoring[i] = contrib;
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::TURN_COUNT, contrib);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterMinimumImpedance::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    auto coilDelimitAndCompactOld = settings.get_coil_delimit_and_compact();
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            newScoring.push_back(scoring);
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
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }


    settings.set_coil_delimit_and_compact(coilDelimitAndCompactOld);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterMagneticInductance::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    auto coilDelimitAndCompactOld = settings.get_coil_delimit_and_compact();
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            newScoring.push_back(scoring);
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
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }


    settings.set_coil_delimit_and_compact(coilDelimitAndCompactOld);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterSaturation::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            newScoring.push_back(scoring);
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
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering saturation, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0 && !firstFilter) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
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
    std::vector<std::pair<Magnetic, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex) {
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            newScoring.push_back(scoring);
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
            filteredMagneticsWithScoring.push_back((*unfilteredMagnetics)[i]);
        }
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering temperature, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }

    return filteredMagneticsWithScoring;
}

} // namespace OpenMagnetics
