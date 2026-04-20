#include "advisers/CoreAdviser.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "advisers/MagneticFilter.h"
#include "constructive_models/Bobbin.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/CoilMesher.h"
#include "constructive_models/NumberTurns.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Bobbin.h"
#include <algorithm>
#include <cctype>
#include <set>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <execution>
#include <magic_enum_utility.hpp>
#include <list>
#include <cmrc/cmrc.hpp>
#include "support/Exceptions.h"
#include "support/Logger.h"
#include "support/Settings.h"

CMRC_DECLARE(data);


namespace OpenMagnetics {

void CoreAdviser::set_unique_core_shapes(bool value) {
    _uniqueCoreShapes = value;
}

bool CoreAdviser::get_unique_core_shapes() {
    return _uniqueCoreShapes;
}

void CoreAdviser::set_application(Application value) {
    _application = value;
}

Application CoreAdviser::get_application() {
    return _application;
}

void CoreAdviser::set_mode(CoreAdviser::CoreAdviserModes value) {
    _mode = value;
}

CoreAdviser::CoreAdviserModes CoreAdviser::get_mode() {
    return _mode;
}

void CoreAdviser::set_weights(std::map<CoreAdviserFilters, double> weights) {
    _weights = weights;
}

std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> CoreAdviser::get_scorings(bool weighted){
    std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> swappedScorings;
    for (auto& [filter, aux] : _scorings) {
        auto filterConfiguration = _filterConfiguration[filter];

        auto normalizedScorings = OpenMagnetics::normalize_scoring(aux, weighted? _weights[filter] : 1, false, false);

        for (auto& [name, scoring] : aux) {
            swappedScorings[name][filter] = normalizedScorings[name];
        }
    }
    return swappedScorings;
}

std::vector<double> normalize_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, std::vector<double> newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    auto normalizedScorings = OpenMagnetics::normalize_scoring(newScoring, weight, filterConfiguration);

    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
        (*magneticsWithScoring)[i].second += normalizedScorings[i];
    }

    // stable_sort((*magneticsWithScoring).begin(), (*magneticsWithScoring).end(), [](std::pair<Magnetic, double>& b1, std::pair<Magnetic, double>& b2) {
    //     return b1.second > b2.second;
    // });

    return normalizedScorings;
}

void sort_magnetics_by_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring) {
    stable_sort((*magneticsWithScoring).begin(), (*magneticsWithScoring).end(), [](const std::pair<Magnetic, double>& b1, const std::pair<Magnetic, double>& b2) {
        return b1.second > b2.second;
    }); // F12 FIX: stable_sort for reproducible results
}

// Helper function to determine if a topology stores energy in the magnetic field
// Energy-storing topologies need gapped cores for DC bias handling
// Non-energy-storing (transformer) topologies transfer energy without storage
bool is_energy_storing_topology(std::optional<Topologies> topology) {
    if (!topology.has_value()) {
        return false; // Unknown topology, caller should use other heuristics
    }
    
    switch (topology.value()) {
        // Energy-storing topologies (inductors, flyback-derived)
        // These store energy in the magnetic field each switching cycle
        case Topologies::FLYBACK_CONVERTER:
        case Topologies::BUCK_CONVERTER:
        case Topologies::BOOST_CONVERTER:
        case Topologies::INVERTING_BUCK_BOOST_CONVERTER:
        case Topologies::ISOLATED_BUCK_BOOST_CONVERTER:
        case Topologies::CUK_CONVERTER:
        case Topologies::SEPIC:
        case Topologies::ZETA_CONVERTER:
            return true;
            
        // Transformer topologies (forward-derived)
        // These transfer energy directly without storing it
        case Topologies::SINGLE_SWITCH_FORWARD_CONVERTER:
        case Topologies::TWO_SWITCH_FORWARD_CONVERTER:
        case Topologies::ACTIVE_CLAMP_FORWARD_CONVERTER:
        case Topologies::PUSH_PULL_CONVERTER:
        case Topologies::HALF_BRIDGE_CONVERTER:
        case Topologies::FULL_BRIDGE_CONVERTER:
        case Topologies::PHASE_SHIFTED_FULL_BRIDGE_CONVERTER:
        case Topologies::ISOLATED_BUCK_CONVERTER:
        case Topologies::DUAL_ACTIVE_BRIDGE_CONVERTER:
        case Topologies::LLC_RESONANT_CONVERTER:
        case Topologies::CLLC_RESONANT_CONVERTER:
        case Topologies::WEINBERG_CONVERTER:
        case Topologies::CURRENT_TRANSFORMER:
            return false;
            
        default:
            return false; // Unknown, treat as transformer (safer - no gap)
    }
}

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
        auto core = magnetic.get_core();

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
        auto core = magnetic.get_core();

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
        auto core = magnetic.get_core();

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
        auto core = magnetic.get_core();

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

Coil get_dummy_coil(Inputs inputs) {
    double frequency = 0; 
    double temperature = 0; 
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        frequency = std::max(frequency, Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency());
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }

    // Set round wire with diameter to two times the skin depth 
    auto wire = Wire::get_wire_for_frequency(frequency, temperature, true);
    Winding primaryWinding;
    primaryWinding.set_isolation_side(IsolationSide::PRIMARY);
    primaryWinding.set_name("primary");
    primaryWinding.set_number_parallels(1);
    primaryWinding.set_number_turns(1);
    primaryWinding.set_wire(wire);

    Coil coil;
    coil.set_bobbin("Dummy");
    coil.set_functional_description({primaryWinding});
    return coil;
}

static CoreShape scale_core_shape(const CoreShape& ref, double scaleFactor) {
    CoreShape scaled = ref;
    if (!ref.get_dimensions()) return scaled;
    auto dims = ref.get_dimensions().value();
    for (auto& [key, dim] : dims) {
        if (std::holds_alternative<DimensionWithTolerance>(dim)) {
            auto& dwt = std::get<DimensionWithTolerance>(dim);
            if (dwt.get_nominal()) dwt.set_nominal(dwt.get_nominal().value() * scaleFactor);
            if (dwt.get_minimum()) dwt.set_minimum(dwt.get_minimum().value() * scaleFactor);
            if (dwt.get_maximum()) dwt.set_maximum(dwt.get_maximum().value() * scaleFactor);
        } else {
            dim = std::get<double>(dim) * scaleFactor;
        }
    }
    scaled.set_dimensions(dims);
    scaled.set_name("Custom " + ref.get_name().value_or("unknown"));
    return scaled;
}

std::vector<CoreShape> CoreAdviser::create_custom_core_shapes(Inputs inputs) {
    MagneticFilterAreaProduct apFilter(inputs);
    double apRequired = apFilter.get_estimated_area_product_required(inputs);

    if (coreShapeDatabase.empty()) {
        load_core_shapes();
    }

    std::vector<CoreShapeFamily> gappableFamilies = {
        CoreShapeFamily::E, CoreShapeFamily::ETD, CoreShapeFamily::EQ,
        CoreShapeFamily::RM, CoreShapeFamily::PQ, CoreShapeFamily::EP,
        CoreShapeFamily::EFD, CoreShapeFamily::EL, CoreShapeFamily::ER,
        CoreShapeFamily::EC, CoreShapeFamily::LP, CoreShapeFamily::EPX,
        CoreShapeFamily::U, CoreShapeFamily::T,
        CoreShapeFamily::PLANAR_E, CoreShapeFamily::PLANAR_EL, CoreShapeFamily::PLANAR_ER
    };

    std::vector<CoreShape> customShapes;
    for (auto family : gappableFamilies) {
        try {
            CoreShape refShape = find_core_shape_by_area_product(apRequired, family);
            Core refCore(refShape);
            refCore.process_data();
            auto columns = refCore.get_columns();
            if (columns.empty()) continue;
            double apRef = columns[0].get_area() * refCore.get_winding_window().get_area().value();
            if (apRef <= 0) continue;

            double s = std::pow(apRequired / apRef, 0.25);
            s = std::clamp(s, 0.5, 3.0);
            CoreShape scaled = scale_core_shape(refShape, s);
            customShapes.push_back(scaled);
        } catch (...) { continue; }
    }
    return customShapes;
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults){
    if (get_mode() == CoreAdviserModes::AVAILABLE_CORES) {
        if (coreDatabase.empty()) {
            load_cores();
        }
        return get_advised_core(inputs, weights, &coreDatabase, maximumNumberResults);
    }
    else if (get_mode() == CoreAdviserModes::STANDARD_CORES) {
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        if (coreShapeDatabase.empty()) {
            load_core_shapes();
        }
        std::vector<MAS::CoreShape> shapes;
        for (auto [name, shape] : coreShapeDatabase) {
            shapes.push_back(shape);
        }
        return get_advised_core(inputs, &shapes, maximumNumberResults);
    }
    else if (get_mode() == CoreAdviserModes::CUSTOM_CORES) {
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        if (coreShapeDatabase.empty()) {
            load_core_shapes();
        }
        auto customShapes = create_custom_core_shapes(inputs);
        return get_advised_core(inputs, &customShapes, maximumNumberResults);
    }
    else {
        throw NotImplementedException("Unknown CoreAdviserMode");
    }
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults, maximumNumberCores);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {

    std::vector<std::pair<Mas, double>> results;


    for(size_t i = 0; i < (*cores).size(); i += maximumNumberCores) {
        auto last = std::min((*cores).size(), i + maximumNumberCores);
        std::vector<Core> partialCores;
        partialCores.reserve(last - i);

        move((*cores).begin() + i, (*cores).begin() + last, std::back_inserter(partialCores));
        
        auto partialResult = get_advised_core(inputs, weights, &partialCores, maximumNumberResults);
        std::move(partialResult.begin(), partialResult.end(), std::back_inserter(results));
    }

    stable_sort(results.begin(), results.end(), [](const std::pair<Mas, double>& b1, const std::pair<Mas, double>& b2) {
        return b1.second > b2.second;
    });

    if (results.size() > maximumNumberResults) {
        results = std::vector<std::pair<Mas, double>>(results.begin(), results.end() - (results.size() - maximumNumberResults));
    }
    return results;
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults){
    _weights = weights;
 
    size_t maximumMagneticsAfterFiltering = settings.get_core_adviser_maximum_magnetics_after_filtering();
    std::vector<std::pair<Magnetic, double>> magnetics;

    magnetics = create_magnetic_dataset(inputs, cores, false);

    // logEntry("We start the search with " + std::to_string(magnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.", "CoreAdviser");
    // logEntry("We don't include stacks of cores in our search.", "CoreAdviser");
    std::vector<std::pair<Mas, double>> filteredMagnetics;
    
    if (get_application() == Application::POWER) {
        filteredMagnetics = filter_available_cores_power_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMagnetics.size() >= maximumNumberResults) {
            return filteredMagnetics;
        }

        auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
        if (globalIncludeStacks) {
            expand_magnetic_dataset_with_stacks(inputs, cores, &magnetics);
        }

        logEntry("First attempt produced not enough results, so now we are searching again with " + std::to_string(magnetics.size()) + " magnetics, including up to " + std::to_string(defaults.coreAdviserMaximumNumberStacks) + " cores stacked when possible.", "CoreAdviser");
        maximumMagneticsAfterFiltering = magnetics.size();
        filteredMagnetics = filter_available_cores_power_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        return filteredMagnetics;
    }
    else {
        filteredMagnetics = filter_available_cores_suppression_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        return filteredMagnetics;
    }
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<CoreShape>* shapes, size_t maximumNumberResults) {
    auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
    auto magnetics = create_magnetic_dataset(inputs, shapes, globalIncludeStacks);

    size_t maximumMagneticsAfterFiltering = settings.get_core_adviser_maximum_magnetics_after_filtering();
    if (get_application() == Application::POWER) {
        return filter_standard_cores_power_application(&magnetics, inputs, _weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
    else {
        return filter_standard_cores_interference_suppression_application(&magnetics, inputs, _weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::create_magnetic_dataset(Inputs inputs, std::vector<Core>* cores, bool includeStacks) {
    std::vector<std::pair<Magnetic, double>> magnetics;
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings.get_use_toroidal_cores();
    auto includeConcentricCores = settings.get_use_concentric_cores();
    auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
    auto globalIncludeDistributedGaps = settings.get_core_adviser_include_distributed_gaps();
    double maximumHeight = DBL_MAX;
    if (inputs.get_design_requirements().get_maximum_dimensions()) {
        if (inputs.get_design_requirements().get_maximum_dimensions()->get_height()) {
            maximumHeight = inputs.get_design_requirements().get_maximum_dimensions()->get_height().value();
        }
    }

    Magnetic magnetic;

    magnetic.set_coil(std::move(coil));
    for (auto& core : *cores){
        auto coreMaterial = core.resolve_material();
        if (!Core::check_material_application(coreMaterial, get_application())) {
            continue;
        }

        if (get_application() == Application::INTERFERENCE_SUPPRESSION) {
            if (core.get_type() != CoreType::TOROIDAL) {
                continue;
            }
        }
        else {
            if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
                continue;
            }

            if (!includeConcentricCores && (core.get_type() == CoreType::PIECE_AND_PLATE || core.get_type() == CoreType::TWO_PIECE_SET)) {
                continue;
            }
        }
        core.process_data();

        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            if (core.get_type() == CoreType::TOROIDAL) {
                continue;
            }
            if (core.get_columns().size() == 2) {
                continue;
            }
            auto windingWindow = core.get_winding_window();
            if (windingWindow.get_height().value() > windingWindow.get_width().value()) {
                continue;
            }
        }

        if (!core.process_gap()) {
            continue;
        }

        if (core.get_type() == CoreType::TWO_PIECE_SET) {
            if (core.get_height() > maximumHeight) {
                continue;
            }
        }

        if (!globalIncludeDistributedGaps && core.get_gapping().size() > core.get_processed_description()->get_columns().size()) {
            continue;
        }

        if (includeStacks && globalIncludeStacks && (core.get_shape_family() == CoreShapeFamily::E || core.get_shape_family() == CoreShapeFamily::PLANAR_E || core.get_shape_family() == CoreShapeFamily::T || core.get_shape_family() == CoreShapeFamily::U || core.get_shape_family() == CoreShapeFamily::C)) {
            for (size_t i = 0; i < defaults.coreAdviserMaximumNumberStacks; ++i)
            {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                // process_data() resets processed description to base values, then calls scale_to_stacks internally
                core.process_data();
                core.process_gap(); // CA-OPT-2 FIX: reprocess gap data after stacking (was commented out)
                magnetic.set_core(core);
                MagneticManufacturerInfo magneticmanufacturerinfo;
                if (i != 0) {
                    magneticmanufacturerinfo.set_reference(core.get_name().value_or("unnamed") + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticmanufacturerinfo.set_reference(core.get_name().value_or("unnamed"));
                }
                magnetic.set_manufacturer_info(magneticmanufacturerinfo);
                magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }
        else {
            magnetic.set_core(core);
            MagneticManufacturerInfo magneticmanufacturerinfo;
            magneticmanufacturerinfo.set_reference(core.get_name().value_or("unnamed"));
            magnetic.set_manufacturer_info(magneticmanufacturerinfo);
            magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
        }
    }

    return magnetics;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::create_magnetic_dataset(Inputs inputs, std::vector<CoreShape>* shapes, bool includeStacks) {
    std::vector<std::pair<Magnetic, double>> magnetics;
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings.get_use_toroidal_cores();
    auto includeConcentricCores = settings.get_use_concentric_cores();
    auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
    auto globalIncludeDistributedGaps = settings.get_core_adviser_include_distributed_gaps();
    double maximumHeight = DBL_MAX;
    if (inputs.get_design_requirements().get_maximum_dimensions()) {
        if (inputs.get_design_requirements().get_maximum_dimensions()->get_height()) {
            maximumHeight = inputs.get_design_requirements().get_maximum_dimensions()->get_height().value();
        }
    }

    Magnetic magnetic;

    magnetic.set_coil(std::move(coil));
    for (auto& shape : *shapes){
        if (shape.get_family() == CoreShapeFamily::PQI || shape.get_family() == CoreShapeFamily::UI || shape.get_family() == CoreShapeFamily::UT || shape.get_family() == CoreShapeFamily::C) {
            continue;
        }
        Core core(shape);

        if (get_application() == Application::INTERFERENCE_SUPPRESSION) {
            if (core.get_type() != CoreType::TOROIDAL) {
                continue;
            }
        }
        else {
            if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
                continue;
            }

            if (!includeConcentricCores && (core.get_type() == CoreType::PIECE_AND_PLATE || core.get_type() == CoreType::TWO_PIECE_SET)) {
                continue;
            }
        }
        core.process_data();

        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            if (core.get_type() == CoreType::TOROIDAL) {
                continue;
            }
            if (core.get_columns().size() == 2) {
                continue;
            }
            auto windingWindow = core.get_winding_window();
            if (windingWindow.get_height().value() > windingWindow.get_width().value()) {
                continue;
            }
        }


        if (!core.process_gap()) {
            continue;
        }

        if (core.get_type() == CoreType::TWO_PIECE_SET) {
            if (core.get_height() > maximumHeight) {
                continue;
            }
        }

        if (!globalIncludeDistributedGaps && core.get_gapping().size() > core.get_processed_description()->get_columns().size()) {
            continue;
        }

        if (includeStacks && globalIncludeStacks && (core.get_shape_family() == CoreShapeFamily::E || core.get_shape_family() == CoreShapeFamily::PLANAR_E || core.get_shape_family() == CoreShapeFamily::T || core.get_shape_family() == CoreShapeFamily::U || core.get_shape_family() == CoreShapeFamily::C)) {
            for (size_t i = 0; i < defaults.coreAdviserMaximumNumberStacks; ++i) {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                // process_data() resets processed description to base values, then calls scale_to_stacks internally
                core.process_data();
                core.process_gap(); // CA-OPT-2 FIX: reprocess gap data after stacking (was commented out)
                magnetic.set_core(core);
                MagneticManufacturerInfo magneticManufacturerInfo;
                if (i!=0) {
                    magneticManufacturerInfo.set_reference(core.get_name().value_or("unnamed") + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticManufacturerInfo.set_reference(core.get_name().value_or("unnamed"));
                }
                magnetic.set_manufacturer_info(magneticManufacturerInfo);
                magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }
        else {
            magnetic.set_core(core);
            MagneticManufacturerInfo magneticManufacturerInfo;
            magneticManufacturerInfo.set_reference(core.get_name().value_or("unnamed"));
            magnetic.set_manufacturer_info(magneticManufacturerInfo);
            magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
        }
    }

    return magnetics;
}

void CoreAdviser::expand_magnetic_dataset_with_stacks(Inputs inputs, std::vector<Core>* cores, std::vector<std::pair<Magnetic, double>>* magnetics) {
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings.get_use_toroidal_cores();
    double maximumHeight = DBL_MAX;
    if (inputs.get_design_requirements().get_maximum_dimensions()) {
        if (inputs.get_design_requirements().get_maximum_dimensions()->get_height()) {
            maximumHeight = inputs.get_design_requirements().get_maximum_dimensions()->get_height().value();
        }
    }

    Magnetic magnetic;

    magnetic.set_coil(std::move(coil));
    for (auto& core : *cores){
        if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
            continue;
        }

        if (core.get_type() == CoreType::TWO_PIECE_SET) {
            if (core.get_height() > maximumHeight) {
                continue;
            }
        }

        if (core.get_shape_family() == CoreShapeFamily::E || core.get_shape_family() == CoreShapeFamily::PLANAR_E || core.get_shape_family() == CoreShapeFamily::T || core.get_shape_family() == CoreShapeFamily::U || core.get_shape_family() == CoreShapeFamily::C) {

            core.process_data();
            if (!core.process_gap()) {
                continue;
            }

            for (size_t i = 1; i < defaults.coreAdviserMaximumNumberStacks; ++i)
            {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                // process_data() resets processed description to base values, then calls scale_to_stacks internally
                core.process_data();
                core.process_gap(); // CA-OPT-2 FIX: reprocess gap data after stacking (was commented out)
                MagneticManufacturerInfo magneticManufacturerInfo;
                if (i!=0) {
                    std::string name = core.get_name().value_or("unnamed");
                    name = std::regex_replace(std::string(name), std::regex(" [0-9] stacks"), "");
                    core.set_name(name + " " + std::to_string(1 + i) + " stacks" );
                    magneticManufacturerInfo.set_reference(core.get_name().value_or("unnamed"));
                }
                magnetic.set_manufacturer_info(magneticManufacturerInfo);
                magnetic.set_core(core);
                (*magnetics).push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }

    }
}

void add_initial_turns_by_inductance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    MagnetizingInductance magnetizingInductance;
    
    // Transformer vs Inductor/Energy-Storing Detection
    // =================================================
    // We need to determine if this application stores energy (needs gap) or transfers it (no gap).
    //
    // ENERGY-STORING (inductor/flyback-like):
    //   - Topologies: Flyback, Buck, Boost, Buck-Boost, SEPIC, Cuk, Zeta, etc.
    //   - Energy is stored in the magnetic field each switching cycle
    //   - Core typically needs a gap to store the required energy without saturating
    //   - Turns calculated from inductance requirement and gap
    //
    // TRANSFORMER (forward-like):
    //   - Topologies: Forward, PushPull, Half-Bridge, Full-Bridge, DAB, LLC, etc.
    //   - Energy is transferred directly, not stored in the magnetic field
    //   - Magnetizing inductance should be HIGH to minimize magnetizing current (parasitic)
    //   - Core should NOT be gapped (high permeability = high inductance)
    //   - Turns calculated from voltage/frequency to avoid saturation (Faraday's Law)
    //
    // Detection priority:
    // 1. Use topology if specified (most reliable)
    // 2. Fall back to inductance field heuristic (minimum-only = transformer)
    //
    auto topology = inputs.get_design_requirements().get_topology();
    bool isEnergyStoring = is_energy_storing_topology(topology);
    
    // If topology is not set, use the old heuristic based on inductance specification
    // minimum-only inductance suggests transformer (want high inductance, no specific value)
    // nominal or min+max suggests inductor (need specific value for energy storage)
    bool isTransformer;
    if (topology.has_value()) {
        isTransformer = !isEnergyStoring;
    } else {
        // Legacy heuristic: minimum-only = transformer
        isTransformer = inputs.get_design_requirements().get_magnetizing_inductance().get_minimum() &&
                        !inputs.get_design_requirements().get_magnetizing_inductance().get_nominal() &&
                        !inputs.get_design_requirements().get_magnetizing_inductance().get_maximum();
    }
    
    // Pre-compute shared transformer values once (not per-core)
    double transformerTemperature = 25.0;
    double maxVoltSeconds = 0.0;
    if (isTransformer) {
        for (size_t opIdx = 0; opIdx < inputs.get_operating_points().size(); ++opIdx) {
            auto op = inputs.get_operating_point(opIdx);
            transformerTemperature = std::max(transformerTemperature, op.get_conditions().get_ambient_temperature());
            auto excitation = Inputs::get_primary_excitation(op);

            if (excitation.get_voltage() && excitation.get_voltage()->get_waveform() &&
                excitation.get_voltage()->get_waveform()->get_time()) {
                auto voltageWaveform = excitation.get_voltage()->get_waveform().value();
                const auto& data = voltageWaveform.get_data();
                auto time = voltageWaveform.get_time().value();
                double integral = 0.0;
                for (size_t j = 0; j + 1 < std::min(data.size(), time.size()); ++j) {
                    integral += data[j] * (time[j + 1] - time[j]);
                    maxVoltSeconds = std::max(maxVoltSeconds, std::abs(integral));
                }
            } else if (excitation.get_voltage() && excitation.get_voltage()->get_processed() &&
                       excitation.get_voltage()->get_processed()->get_peak()) {
                double frequency = excitation.get_frequency() > 0 ? excitation.get_frequency() : 100000;
                double omega = 2 * std::numbers::pi * frequency;
                maxVoltSeconds = std::max(maxVoltSeconds,
                    excitation.get_voltage()->get_processed()->get_peak().value() / omega);
            }
        }
    }

    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){

        Core core = (*magneticsWithScoring)[i].first.get_core();
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        double initialNumberTurns = (*magneticsWithScoring)[i].first.get_coil().get_functional_description()[0].get_number_turns();

        if (initialNumberTurns == 1) {
            if (isTransformer) {
                // Feasibility seed: fewest turns that keep B_peak under the material's safe
                // operating flux. The loss filter refines N upward toward the loss optimum.
                double bMax = core.get_magnetic_flux_density_saturation(transformerTemperature, true);
                double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();

                if (maxVoltSeconds > 0 && effectiveArea > 0 && bMax > 0) {
                    initialNumberTurns = std::max(1.0, std::ceil(maxVoltSeconds / (bMax * effectiveArea)));
                } else {
                    initialNumberTurns = 5;
                }
            } else {
                // For inductors/flybacks, calculate from inductance requirement
                initialNumberTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs, DimensionalValues::MINIMUM);
            }
        }
        if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
            NumberTurns numberTurns(initialNumberTurns, inputs.get_design_requirements());
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            initialNumberTurns = numberTurnsCombination[0];
        }

        (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(initialNumberTurns);
    }
}

std::vector<std::pair<Magnetic, double>> add_initial_turns_by_impedance(std::vector<std::pair<Magnetic, double>> magneticsWithScoring, Inputs inputs) {
    Impedance impedance;
    std::vector<std::pair<Magnetic, double>> magneticsWithScoringAndTurns;
    for (size_t i = 0; i < magneticsWithScoring.size(); ++i){
        auto [magnetic, scoring] = magneticsWithScoring[i];
        Core core = magnetic.get_core();
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        Bobbin bobbin;
        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            bobbin = Bobbin::create_quick_bobbin(core, true);
        }
        else {
            bobbin = Bobbin::create_quick_bobbin(core);
        }
        magnetic.get_mutable_coil().set_bobbin(bobbin);

        double initialNumberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();

        try {
            initialNumberTurns = impedance.calculate_minimum_number_turns(magnetic, inputs);
            if (initialNumberTurns < 1) {
                continue;
            }
        }
        catch (...) { // XC-3 NOTE: consider catching std::exception for diagnostics
            continue;
        }
        if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
            NumberTurns numberTurns(initialNumberTurns, inputs.get_design_requirements());
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            initialNumberTurns = numberTurnsCombination[0];
        }

        magnetic.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(initialNumberTurns);
        magneticsWithScoringAndTurns.push_back({magnetic, scoring});
    }

    return magneticsWithScoringAndTurns;
}

void add_alternative_materials(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    CoreMaterialCrossReferencer coreMaterialCrossReferencer(std::map<std::string, std::string>{{"coreLosses", "Steinmetz"}});
    double temperature = 0; 
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }
    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){
        Core core = (*magneticsWithScoring)[i].first.get_core();
        auto coreMaterial = core.resolve_material();

        auto alternatives = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, temperature);
        std::vector<std::string> coreMaterialAlternatives;
        for (auto [alternativeCoreMaterial, scoring] : alternatives) {
            coreMaterialAlternatives.push_back(alternativeCoreMaterial.get_name());
        }
        coreMaterial.set_alternatives(coreMaterialAlternatives);
        core.set_material(coreMaterial);

        (*magneticsWithScoring)[i].first.set_core(core);
    }
}

void add_gapping(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    MagneticEnergy magneticEnergy;
    
    // Gapping Decision Logic (same as add_gapping_standard_cores):
    // - Skip gapping only if we have minimum-only inductance (want maximum L, no specific target)
    // - Energy-storing topologies or specific inductance targets need gap calculation
    auto inductanceReq = inputs.get_design_requirements().get_magnetizing_inductance();
    bool hasNominalInductance = inductanceReq.get_nominal().has_value();
    bool hasMaxInductance = inductanceReq.get_maximum().has_value();
    bool needsSpecificInductance = hasNominalInductance || hasMaxInductance;
    
    auto topology = inputs.get_design_requirements().get_topology();
    bool isEnergyStoring = is_energy_storing_topology(topology);
    bool isTransformer = topology.has_value() ? !isEnergyStoring : 
        (inductanceReq.get_minimum() && !hasNominalInductance && !hasMaxInductance);
    
    bool skipGapping = isTransformer && !needsSpecificInductance;
    
    if (skipGapping)
    {
        for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
            Core core = (*magneticsWithScoring)[i].first.get_core();
            core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            (*magneticsWithScoring)[i].first.set_core(core);
        }

        return;
    }
    auto requiredMagneticEnergy = resolve_dimensional_values(magneticEnergy.calculate_required_magnetic_energy(inputs), DimensionalValues::MAXIMUM);
    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
        Core core = (*magneticsWithScoring)[i].first.get_core();

        if (core.get_material_name() == "Dummy") {
            core.set_material_initial_permeability(defaults.ferriteInitialPermeability);
        }
        if (!core.get_processed_description()) {
            core.process_data();
        }
        if (core.get_shape_family() != CoreShapeFamily::T) {
            // Use realistic ferrite Bsat (~350 mT) instead of dummy material's 500 mT
// CA-LOGIC-3 FIX: Use actual material Bsat instead of hardcoded ferrite default
            double realisticBsat = defaults.ferriteSaturationFluxDensity; // fallback
            try {
                auto coreMat = core.resolve_material();
                double opTemp = 25.0;
                for (auto& op : inputs.get_operating_points()) {
                    opTemp = std::max(opTemp, op.get_conditions().get_ambient_temperature());
                }
                realisticBsat = Core::get_magnetic_flux_density_saturation(coreMat, opTemp);
            } catch (...) { // XC-3 NOTE: consider catching std::exception for diagnostics
                // CA-LOGIC-3: fallback to ferrite default if material data unavailable
            }
            double bSatTarget = realisticBsat * defaults.maximumProportionMagneticFluxDensitySaturation; // Use configured proportion of Bsat for safety margin
            
            // Calculate gap based on energy storage requirement
            double gapLength = roundFloat(magneticEnergy.calculate_gap_length_by_magnetic_energy(core.get_gapping()[0], realisticBsat, requiredMagneticEnergy), 5);
            double gapEnergy = gapLength;
            
            // Iterative refinement to find gap that prevents saturation
            double magnetizingCurrentPeak = 0;
            for (auto& op : inputs.get_operating_points()) {
                auto excitation = Inputs::get_primary_excitation(op);
                if (excitation.get_magnetizing_current() && excitation.get_magnetizing_current()->get_processed() && excitation.get_magnetizing_current()->get_processed()->get_peak()) {
                    magnetizingCurrentPeak = std::max(magnetizingCurrentPeak, excitation.get_magnetizing_current()->get_processed()->get_peak().value());
                }
            }
            
            // Get the actual current from primary excitation (this is what saturation filter will use)
            double actualCurrentPeak = 0;
            for (auto& op : inputs.get_operating_points()) {
                auto excitation = Inputs::get_primary_excitation(op);
                if (excitation.get_current() && excitation.get_current()->get_processed() && excitation.get_current()->get_processed()->get_peak()) {
                    actualCurrentPeak = std::max(actualCurrentPeak, excitation.get_current()->get_processed()->get_peak().value());
                }
            }
            
            // Use the larger of magnetizing current or actual current
            double currentForCalculation = std::max(magnetizingCurrentPeak, actualCurrentPeak);
            
            if (currentForCalculation > 0) {
                // Use the proper method to calculate gap from saturation constraint
                MagnetizingInductance magnetizingInductance;
                double gapFromSaturation = magnetizingInductance.calculate_gap_from_saturation_constraint(
                    core, &inputs, bSatTarget, currentForCalculation);
                
                // Use the maximum of energy gap and saturation gap
                gapLength = std::max(gapEnergy, gapFromSaturation);
            }

            core.set_ground_gapping(gapLength);
            core.process_gap();
            std::stringstream ss;
            ss << roundFloat(gapLength * 1000, 2);
            if (gapLength > 0) {
                core.set_name(core.get_name().value_or("unnamed") + " gapped " + ss.str()  + " mm");
            }
            else {
                core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            }
        }

        (*magneticsWithScoring)[i].first.set_core(core);
    }
}

// ============================================================================
// New Gapping Algorithm for STANDARD_CORES Mode
// ============================================================================

CoreAdviser::GappingConstraints CoreAdviser::calculate_gapping_constraints(Inputs inputs, Core core) {
    GappingConstraints constraints;

    // Get required magnetic energy
    MagneticEnergy magneticEnergy;
    auto requiredMagneticEnergy = resolve_dimensional_values(
        magneticEnergy.calculate_required_magnetic_energy(inputs),
        DimensionalValues::MAXIMUM);

    // Account for stacked cores - each stack shares the energy
    auto numStacksOpt = core.get_functional_description().get_number_stacks();
    int64_t numStacks = numStacksOpt ? numStacksOpt.value() : 1;
    if (numStacks > 1) {
        requiredMagneticEnergy /= numStacks;
    }

    // Get realistic Bsat for the material at operating temperature
    double opTemp = 25.0;
    for (auto& op : inputs.get_operating_points()) {
        opTemp = std::max(opTemp, op.get_conditions().get_ambient_temperature());
    }
    double realisticBsat = core.get_magnetic_flux_density_saturation(opTemp, false);

    // 1. Calculate minimum gap: energy storage requirement
    // Use the maximum allowed flux density (proportional to Bsat) for energy calculation
    double maxAllowedB = realisticBsat * defaults.maximumProportionMagneticFluxDensitySaturation;
    double minGap = magneticEnergy.calculate_gap_length_by_magnetic_energy(
        core.get_gapping()[0], maxAllowedB, requiredMagneticEnergy);
    constraints.minGap = minGap;

    // 2. Calculate maximum gap constraints
    // 2a. Maximum gap based on column width (30%)
    double columnWidth = core.get_columns()[0].get_width();
    double maxGapByColumn = 0.30 * columnWidth;

    // 2b. Maximum gap based on fringing factor (25%)
    double maxGapByFringing = calculate_gap_for_fringing_factor(0.25, core);

    // Use the more restrictive limit
    constraints.maxGap = std::min(maxGapByColumn, maxGapByFringing);

    // Ensure maxGap >= minGap (sanity check)
    if (constraints.maxGap < constraints.minGap) {
        constraints.maxGap = constraints.minGap;
    }

    // 3. Calculate saturation constraint gap (ensure we don't exceed 70% Bsat)
    MagnetizingInductance magnetizingInductance;
    double targetB = realisticBsat * defaults.maximumProportionMagneticFluxDensitySaturation;
    double peakCurrent = get_peak_current(inputs);
    
    // Account for stacked cores - current is divided among stacks
    if (numStacks > 1) {
        peakCurrent /= numStacks;
    }
    
    double saturationGap = magnetizingInductance.calculate_gap_from_saturation_constraint(
        core, &inputs, targetB, peakCurrent);

    // The minimum gap must satisfy BOTH energy storage AND saturation constraints
    double effectiveMinGap = std::max(constraints.minGap, saturationGap);
    
    // 4. Check if we can satisfy both saturation and fringing constraints
    if (effectiveMinGap > constraints.maxGap) {
        // Cannot satisfy both - prioritize saturation constraint
        // This core will have higher fringing but won't saturate
        constraints.maxGap = effectiveMinGap * 1.5;  // Allow some margin for optimization
    }
    
    // 5. Find optimal gap
    if (settings.get_gapping_strategy() == GappingOptimizationStrategy::GOLDEN_SECTION) {
        constraints.optimalGap = optimize_gap_golden_section(
            effectiveMinGap, constraints.maxGap, inputs, core);
    } else {
        // SIMPLE strategy: use effective minimum gap
        constraints.optimalGap = effectiveMinGap;
    }

    return constraints;
}

double CoreAdviser::calculate_gap_for_fringing_factor(double targetFringingFactor, Core core) {
    // Use the reluctance model to find gap length for target fringing factor
    auto reluctanceModel = ReluctanceModel::factory(_models);
    return reluctanceModel->get_gapping_by_fringing_factor(core, targetFringingFactor);
}

double CoreAdviser::get_peak_current(Inputs inputs) {
    double peakCurrent = 0.0;

    // For transformer topologies (forward converters), core saturation is driven by
    // magnetizing current only. The reflected secondary current is balanced and does
    // not contribute to net flux. Using actual current would oversize the core.
    // When topology is unset, fall back to inductance-based heuristic:
    // minimum-only inductance = transformer; nominal/max inductance = inductor.
    auto topology = inputs.get_design_requirements().get_topology();
    bool isTransformerTopology;
    if (topology.has_value()) {
        isTransformerTopology = !is_energy_storing_topology(topology);
    } else {
        auto& inductanceReq = inputs.get_design_requirements().get_magnetizing_inductance();
        isTransformerTopology = inductanceReq.get_minimum() &&
                                !inductanceReq.get_nominal() &&
                                !inductanceReq.get_maximum();
    }

    for (auto& op : inputs.get_operating_points()) {
        auto excitation = Inputs::get_primary_excitation(op);

        // Check magnetizing current
        if (excitation.get_magnetizing_current() &&
            excitation.get_magnetizing_current()->get_processed() &&
            excitation.get_magnetizing_current()->get_processed()->get_peak()) {
            peakCurrent = std::max(peakCurrent,
                excitation.get_magnetizing_current()->get_processed()->get_peak().value());
        }

        // Check actual current only for energy-storing topologies (e.g. flyback, boost)
        // and inductors (identified by having a nominal inductance).
        // For pure transformer topologies the actual current includes reflected secondary
        // load current which does not bias the core.
        if (!isTransformerTopology &&
            excitation.get_current() &&
            excitation.get_current()->get_processed() &&
            excitation.get_current()->get_processed()->get_peak()) {
            peakCurrent = std::max(peakCurrent,
                excitation.get_current()->get_processed()->get_peak().value());
        }
    }

    return peakCurrent;
}

double CoreAdviser::calculate_core_losses_for_gap(double gap, Inputs inputs, Core core) {
    try {
        // Set the gap temporarily
        core.set_ground_gapping(gap);
        core.process_gap();

        // Create Steinmetz model
        auto coreLossesModel = CoreLossesModel::factory(
            std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));

        // Calculate losses for each operating point
        double totalLosses = 0.0;
        for (auto& op : inputs.get_operating_points()) {
            auto excitation = Inputs::get_primary_excitation(op);
            double temperature = op.get_conditions().get_ambient_temperature();
            auto losses = coreLossesModel->get_core_losses(core, excitation, temperature);
            totalLosses += losses.get_core_losses();
        }

        return totalLosses;
    } catch (const std::exception& e) {
        // If we can't calculate losses (e.g., missing Steinmetz data), return a large value
        // This will make this gap less attractive during optimization
        return std::numeric_limits<double>::max();
    }
}

double CoreAdviser::optimize_gap_golden_section(double minGap, double maxGap, Inputs inputs, Core core) {
    const double PHI = 1.618033988749895;  // Golden ratio
    const double RESPHI = 2.0 - PHI;       // 1/phi^2
    const int MAX_ITERATIONS = 10;

    double a = minGap;
    double b = maxGap;

    // Initial points
    double c = a + RESPHI * (b - a);
    double d = b - RESPHI * (b - a);

    double fc = calculate_core_losses_for_gap(c, inputs, core);
    double fd = calculate_core_losses_for_gap(d, inputs, core);

    for (int i = 0; i < MAX_ITERATIONS && (b - a) > 1e-6; ++i) {
        if (fc < fd) {
            // Minimum is in [a, d]
            b = d;
            d = c;
            fd = fc;
            c = a + RESPHI * (b - a);
            fc = calculate_core_losses_for_gap(c, inputs, core);
        } else {
            // Minimum is in [c, b]
            a = c;
            c = d;
            fc = fd;
            d = b - RESPHI * (b - a);
            fd = calculate_core_losses_for_gap(d, inputs, core);
        }
    }

    // Return the midpoint of the final interval
    return (a + b) / 2.0;
}

void CoreAdviser::add_gapping_standard_cores(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring,
                                              Inputs inputs) {
    if (magneticsWithScoring->empty()) {
        return;
    }

    // Gapping Decision Logic:
    //
    // The decision to gap depends on TWO factors:
    // 1. Whether the topology stores energy (needs gap for DC bias)
    // 2. Whether a specific inductance value is required (needs gap to control AL)
    //
    // Cases:
    // - Energy-storing (Flyback, Buck, etc.): ALWAYS needs gap for energy storage
    // - Transformer with nominal/max inductance: MAY need gap to hit specific L value
    //   (e.g., LLC resonant converter needs controlled magnetizing inductance)
    // - Transformer with minimum-only inductance: NO gap needed (want maximum L)
    //
    auto topology = inputs.get_design_requirements().get_topology();
    bool isEnergyStoring = is_energy_storing_topology(topology);
    
    // Check if inductance requirement needs a specific value (nominal or bounded)
    auto inductanceReq = inputs.get_design_requirements().get_magnetizing_inductance();
    bool hasNominalInductance = inductanceReq.get_nominal().has_value();
    bool hasMaxInductance = inductanceReq.get_maximum().has_value();
    bool needsSpecificInductance = hasNominalInductance || hasMaxInductance;
    
    // Determine if we should skip gapping:
    // Only skip if it's a transformer topology AND we don't need a specific inductance value
    bool isTransformer = topology.has_value() ? !isEnergyStoring : 
        (inductanceReq.get_minimum() && !hasNominalInductance && !hasMaxInductance);
    
    bool skipGapping = isTransformer && !needsSpecificInductance;
    
    if (skipGapping) {
        // Transformer with minimum-only inductance: no gap needed (want maximum L)
        for (size_t i = 0; i < magneticsWithScoring->size(); ++i) {
            Core core = (*magneticsWithScoring)[i].first.get_core();
            core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            (*magneticsWithScoring)[i].first.set_core(core);
        }
        return;
    }

    for (size_t i = 0; i < magneticsWithScoring->size(); ++i) {
        Core core = (*magneticsWithScoring)[i].first.get_core();

        if (core.get_shape_family() == CoreShapeFamily::T) {
            // Toroidal cores don't have gaps
            core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            (*magneticsWithScoring)[i].first.set_core(core);
            continue;
        }

        // Process core data if needed
        if (!core.get_processed_description()) {
            core.process_data();
        }

        // Calculate gapping constraints
        auto constraints = calculate_gapping_constraints(inputs, core);

        // Apply the optimal gap
        core.set_ground_gapping(constraints.optimalGap);
        core.process_gap();

        // Update name with gap info (avoid duplicates)
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << (constraints.optimalGap * 1000);
        auto currentName = core.get_name().value_or("unnamed");
        // Remove existing gap info if present to avoid duplicates
        size_t gappedPos = currentName.find(" gapped ");
        if (gappedPos != std::string::npos) {
            currentName = currentName.substr(0, gappedPos);
        }
        size_t ungappedPos = currentName.find(" ungapped");
        if (ungappedPos != std::string::npos) {
            currentName = currentName.substr(0, ungappedPos);
        }
        if (constraints.optimalGap > 0) {
            core.set_name(currentName + " gapped " + ss.str() + " mm");
        } else {
            core.set_name(currentName + " ungapped");
        }

        (*magneticsWithScoring)[i].first.set_core(core);
    }
}

void CoreAdviser::refine_gaps_for_saturation(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring,
                                              Inputs inputs) {
    if (magneticsWithScoring->empty()) {
        return;
    }

    MagnetizingInductance magnetizingInductance;
    MagneticSimulator magneticSimulator;
    const int MAX_ITERATIONS = 5;  // Reduced from 10 to improve performance

    for (size_t i = 0; i < magneticsWithScoring->size();) {
        Core core = (*magneticsWithScoring)[i].first.get_core();
        
        if (core.get_shape_family() == CoreShapeFamily::T) {
            // Toroidal cores don't have gaps to refine
            ++i;
            continue;
        }

        // Skip powder cores - they have distributed gaps and are designed for high DC current
        // Powder cores are identified by their material type
        auto materialName = core.get_material_name();
        std::string mat = materialName;
        if (mat.find("Kool M") != std::string::npos || 
            mat.find("High Flux") != std::string::npos ||
            mat.find("XFlux") != std::string::npos ||
            mat.find("Edge") != std::string::npos ||
            mat.find("Fluxsan") != std::string::npos ||
            mat.find("FS ") == 0 ||  // Fluxsan prefix
            mat.find("HF ") == 0 ||  // High Flux prefix
            mat.find("GX ") == 0) {  // Another powder series
            // Powder core - skip gap refinement, keep as-is
            ++i;
            continue;
        }

        // Target safe operating flux for this core (material Bsat curve at opTemp,
        // capped by defaults.maximumProportionMagneticFluxDensitySaturation).
        double opTemp = 25.0;
        for (auto& op : inputs.get_operating_points()) {
            opTemp = std::max(opTemp, op.get_conditions().get_ambient_temperature());
        }
        double maxB = core.get_magnetic_flux_density_saturation(opTemp, /*proportion=*/true);

        // Get current gap
        double currentGap = 0;
        if (!core.get_functional_description().get_gapping().empty()) {
            currentGap = core.get_functional_description().get_gapping()[0].get_length();
        }

        // Iterative refinement to achieve target saturation
        int iteration = 0;
        double bPeak = 0;
        bool converged = false;
        bool shouldRemove = false;

        while (iteration < MAX_ITERATIONS && !converged && !shouldRemove) {
            // Apply current gap
            core.set_ground_gapping(currentGap);
            core.process_gap();


            // Simulate to get actual Bpeak
            bool gotBpeak = false;
            try {
                // Wind the coil first
                auto magnetic = (*magneticsWithScoring)[i].first;
                magnetic.get_mutable_coil().wind();
                
                auto tempMas = Mas();
                tempMas.set_magnetic(magnetic);
                tempMas.set_inputs(inputs);
                auto simulatedMas = magneticSimulator.simulate(tempMas, true);  // fastMode = true
                
                // Get Bpeak from the excitation after simulation
                if (!simulatedMas.get_inputs().get_operating_points().empty()) {
                    auto& op = simulatedMas.get_inputs().get_operating_points()[0];
                    auto excitation = Inputs::get_primary_excitation(op);
                    if (excitation.get_magnetic_flux_density().has_value()) {
                        auto bField = excitation.get_magnetic_flux_density().value();
                        if (bField.get_processed().has_value() && 
                            bField.get_processed()->get_peak().has_value()) {
                            bPeak = bField.get_processed()->get_peak().value();
                            gotBpeak = true;
                        }
                    }
                }
            } catch (const std::exception& e) {
                converged = true;
            }
            
            if (!gotBpeak) {
                converged = true;
                break;
            }

            // Check if we're at or below the safe operating flux.
            if (bPeak <= maxB * 1.01) {  // Allow 1% tolerance
                converged = true;
            } else {
                // Increase gap proportionally to the overshoot, with a 10% margin.
                double newGap = currentGap * (bPeak / maxB) * 1.1;
                
                // Check practical limits
                double columnWidth = core.get_columns()[0].get_width();
                double maxPracticalGap = columnWidth * 0.5;  // 50% of column width
                
                if (newGap > maxPracticalGap) {
                    shouldRemove = true;
                } else {
                    currentGap = newGap;
                    // Recalculate turns for the new gap to maintain inductance
                    double newTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(
                        core, &inputs, DimensionalValues::MINIMUM);
                    (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(newTurns);
                }
            }

            iteration++;
        }

        if (shouldRemove) {
            // Remove this core from the list
            magneticsWithScoring->erase(magneticsWithScoring->begin() + i);
        } else {
            // Apply final gap and update name (avoid duplicates)
            core.set_ground_gapping(currentGap);
            core.process_gap();
            
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << (currentGap * 1000);
            auto currentName = core.get_name().value_or("unnamed");
            // Remove existing gap info if present to avoid duplicates
            size_t gappedPos = currentName.find(" gapped ");
            if (gappedPos != std::string::npos) {
                currentName = currentName.substr(0, gappedPos);
            }
            size_t ungappedPos = currentName.find(" ungapped");
            if (ungappedPos != std::string::npos) {
                currentName = currentName.substr(0, ungappedPos);
            }
            if (currentGap > 0) {
                core.set_name(currentName + " gapped " + ss.str() + " mm");
            } else {
                core.set_name(currentName + " ungapped");
            }
            
            (*magneticsWithScoring)[i].first.set_core(core);
            ++i;
        }
    }
}

// ============================================================================
// Option 2: Binary Search Gap Optimization with Analytical Cost Function
// ============================================================================

double CoreAdviser::calculate_gap_cost_analytical(double gap, Inputs& inputs, Core& core, 
                                                   MagnetizingInductance& magnetizingInductance) {
    auto constants = OpenMagnetics::Constants();
    
    // Get operating conditions
    double temperature = Defaults().ambientTemperature;
    double frequency = Defaults().coreAdviserFrequencyReference;
    if (inputs.get_operating_points().size() > 0) {
        temperature = inputs.get_operating_point(0).get_conditions().get_ambient_temperature();
        frequency = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_frequency();
    }
    
    // Set gap on core
    core.set_ground_gapping(gap);
    core.process_gap();
    
    // Get target inductance
    double targetInductance = resolve_dimensional_values(
        inputs.get_design_requirements().get_magnetizing_inductance(),
        DimensionalValues::NOMINAL);
    
    // Calculate turns for this gap using MagnetizingInductance method
    double turns = magnetizingInductance.calculate_turns_for_gap(core, targetInductance, temperature, frequency);
    
    // Get peak current
    double peakCurrent = get_peak_current(inputs);
    
    // Calculate flux density using MagnetizingInductance method
    double bPeak = magnetizingInductance.calculate_flux_density_peak(core, turns, peakCurrent, temperature, frequency);
    
    // Safe operating flux cap for this material at the operating temperature.
    double maxAllowedB = core.get_magnetic_flux_density_saturation(temperature, /*proportion=*/true);

    // Check saturation constraint
    if (bPeak > maxAllowedB) {
        return std::numeric_limits<double>::max();  // Constraint violation
    }
    
    // Calculate fringing factor penalty using ReluctanceModel
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    double fringingFactor = 0.0;
    if (!core.get_gapping().empty()) {
        auto gapReluctance = reluctanceModel->get_gap_reluctance(core.get_gapping()[0]);
        fringingFactor = gapReluctance.get_fringing_factor();
    }
    
    // Estimate core losses (simplified - just use volume and frequency as proxy)
    // The actual losses depend on B, but we're already constraining B
    double coreVolume = core.get_processed_description()->get_effective_parameters().get_effective_volume();
    double lossProxy = coreVolume * frequency * bPeak * bPeak;  // Simplified Steinmetz-like estimate
    
    // Cost function: minimize losses while penalizing high fringing
    // Fringing factor > 0.25 is considered excessive
    double fringingPenalty = 0.0;
    if (fringingFactor > 0.25) {
        fringingPenalty = pow((fringingFactor - 0.25) * 10, 2);  // Quadratic penalty
    }
    
    return lossProxy + fringingPenalty;
}

std::pair<double, double> CoreAdviser::optimize_gap_and_turns_binary_search(Inputs& inputs, Core& core) {
    MagnetizingInductance magnetizingInductance(_models["gapReluctance"]);
    auto constants = OpenMagnetics::Constants();
    
    // Ensure core is processed
    if (!core.get_processed_description()) {
        core.process_data();
    }
    
    // Get operating conditions
    double temperature = Defaults().ambientTemperature;
    double frequency = Defaults().coreAdviserFrequencyReference;
    if (inputs.get_operating_points().size() > 0) {
        temperature = inputs.get_operating_point(0).get_conditions().get_ambient_temperature();
        frequency = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_frequency();
    }
    
    double targetInductance = resolve_dimensional_values(
        inputs.get_design_requirements().get_magnetizing_inductance(),
        DimensionalValues::NOMINAL);
    double peakCurrent = get_peak_current(inputs);
    double bSat = core.get_magnetic_flux_density_saturation(temperature, false);
    double maxAllowedB = bSat * defaults.maximumProportionMagneticFluxDensitySaturation;
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // Define gap search bounds
    double columnWidth = core.get_columns()[0].get_width();
    double gMin = constants.residualGap;
    double gMax = 0.4 * columnWidth;  // 40% of column width
    
    // Quick check: can this core possibly work?
    // Minimum turns from saturation: N_min = L * I / (B_max * A)
    double minTurnsFromSaturation = (targetInductance * peakCurrent) / (maxAllowedB * effectiveArea);
    if (minTurnsFromSaturation > 1000) {
        // Unreasonable number of turns - core too small
        return {-1.0, -1.0};
    }
    
    // Golden section search for optimal gap
    const double PHI = 1.618033988749895;
    const double RESPHI = 2.0 - PHI;
    const int MAX_ITERATIONS = 15;
    
    double a = gMin;
    double b = gMax;
    
    // Initial points
    double c = a + RESPHI * (b - a);
    double d = b - RESPHI * (b - a);
    
    Core coreC = core;
    Core coreD = core;
    double fc = calculate_gap_cost_analytical(c, inputs, coreC, magnetizingInductance);
    double fd = calculate_gap_cost_analytical(d, inputs, coreD, magnetizingInductance);
    
    for (int i = 0; i < MAX_ITERATIONS && (b - a) > 1e-6; ++i) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = a + RESPHI * (b - a);
            coreC = core;
            fc = calculate_gap_cost_analytical(c, inputs, coreC, magnetizingInductance);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = b - RESPHI * (b - a);
            coreD = core;
            fd = calculate_gap_cost_analytical(d, inputs, coreD, magnetizingInductance);
        }
    }
    
    // Get optimal gap (midpoint of final interval)
    double optimalGap = (a + b) / 2.0;
    
    // Verify this gap is valid
    Core finalCore = core;
    double finalCost = calculate_gap_cost_analytical(optimalGap, inputs, finalCore, magnetizingInductance);
    
    if (finalCost >= std::numeric_limits<double>::max() / 2) {
        // No valid solution found in search range
        // Try using the analytical method from Option 1 as fallback
        double maxB = bSat * defaults.maximumProportionMagneticFluxDensitySaturation;
        auto [fallbackGap, fallbackTurns] = magnetizingInductance.calculate_optimal_gap_and_turns(
            core, &inputs, maxB, peakCurrent);
        return {fallbackGap, fallbackTurns};
    }
    
    // Calculate final turns for optimal gap
    core.set_ground_gapping(optimalGap);
    core.process_gap();
    double optimalTurns = magnetizingInductance.calculate_turns_for_gap(core, targetInductance, temperature, frequency);
    
    return {optimalGap, optimalTurns};
}

void CoreAdviser::add_gapping_and_turns_analytical(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring,
                                                    Inputs inputs) {
    if (magneticsWithScoring->empty()) {
        return;
    }
    
    MagnetizingInductance magnetizingInductance(_models["gapReluctance"]);
    
    // Get operating conditions
    double temperature = Defaults().ambientTemperature;
    if (inputs.get_operating_points().size() > 0) {
        temperature = inputs.get_operating_point(0).get_conditions().get_ambient_temperature();
    }
    double peakCurrent = get_peak_current(inputs);
    
    for (size_t i = 0; i < magneticsWithScoring->size();) {
        Core core = (*magneticsWithScoring)[i].first.get_core();
        
        // Skip toroidal cores - they don't have gaps
        if (core.get_shape_family() == CoreShapeFamily::T) {
            core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            (*magneticsWithScoring)[i].first.set_core(core);
            ++i;
            continue;
        }
        
        // Ensure core is processed
        if (!core.get_processed_description()) {
            core.process_data();
        }
        
        // Get Bsat for this material
        double bSat = core.get_magnetic_flux_density_saturation(temperature, false);
        double maxAllowedB = bSat * defaults.maximumProportionMagneticFluxDensitySaturation;
        
        // Use Option 1: Direct analytical calculation (O(1) per core)
        auto [optimalGap, optimalTurns] = magnetizingInductance.calculate_optimal_gap_and_turns(
            core, &inputs, maxAllowedB, peakCurrent);
        
        if (optimalGap < 0 || optimalTurns < 0) {
            // No valid solution - remove this core
            magneticsWithScoring->erase(magneticsWithScoring->begin() + i);
            continue;
        }
        
        // Apply optimal gap
        core.set_ground_gapping(optimalGap);
        core.process_gap();
        
        // Update name with gap info
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << (optimalGap * 1000);
        auto currentName = core.get_name().value_or("unnamed");
        // Remove existing gap info if present
        size_t gappedPos = currentName.find(" gapped ");
        if (gappedPos != std::string::npos) {
            currentName = currentName.substr(0, gappedPos);
        }
        size_t ungappedPos = currentName.find(" ungapped");
        if (ungappedPos != std::string::npos) {
            currentName = currentName.substr(0, ungappedPos);
        }
        if (optimalGap > 0.0001) {
            core.set_name(currentName + " gapped " + ss.str() + " mm");
        } else {
            core.set_name(currentName + " ungapped");
        }
        
        // Apply optimal turns
        (*magneticsWithScoring)[i].first.set_core(core);
        (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(optimalTurns);
        
        ++i;
    }
}

bool CoreAdviser::should_include_powder(Inputs inputs) {
    // Check if powder cores are disabled in settings
    if (!settings.get_use_powder_cores()) {
        return false;
    }
    
    if (get_application() != Application::POWER) {
        return false;
    }
    double maximumCurrentDcBias = inputs.get_maximum_current_dc_bias();
    if (maximumCurrentDcBias > 1e-3) {
        return true;
    }
    return false;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_powder_materials(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = 10;
    double magneticFluxDensityReference = 0.18;
    std::vector<std::pair<Magnetic, double>> magneticsWithMaterials;
    std::vector<CoreMaterial> coreMaterialsToEvaluate;
    auto coreMaterials = get_core_material_names(settings.get_preferred_core_material_powder_manufacturer());
    for (auto coreMaterial : coreMaterials) {
        auto application = Core::guess_material_application(coreMaterial);
        if (application == _application) {
            coreMaterialsToEvaluate.push_back(Core::resolve_material(coreMaterial));
        }
    }
    std::vector<CoreMaterial> coreMaterialsToUse;
    std::vector<std::pair<CoreMaterial, double>> evaluations;

    double temperature = 0; 
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }
    double maximumCurrentDcBias = inputs.get_maximum_current_dc_bias();

    SignalDescriptor magneticFluxDensity;
    Processed processed;
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(maximumCurrentDcBias);
    processed.set_peak(magneticFluxDensityReference);
    processed.set_peak_to_peak(2 * magneticFluxDensityReference);
    magneticFluxDensity.set_processed(processed);
    OperatingPointExcitation operatingPointExcitation;
    operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);
    operatingPointExcitation.set_frequency(1);

    auto coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
    auto coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
    for (auto coreMaterial : coreMaterialsToEvaluate) {
        double averageVolumetricCoreLosses = 0;
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex){
            auto frequency = inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[0].get_frequency();
            Inputs::scale_time_to_frequency(operatingPointExcitation, frequency, false, false);
            auto coreLossesMethods = Core::get_available_core_losses_methods(coreMaterial);

            if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
                double volumetricCoreLosses = coreLossesModelSteinmetz->get_core_volumetric_losses(coreMaterial, operatingPointExcitation, temperature);
                averageVolumetricCoreLosses += volumetricCoreLosses;
            }
            else {
                double volumetricCoreLosses = coreLossesModelProprietary->get_core_volumetric_losses(coreMaterial, operatingPointExcitation, temperature);
                averageVolumetricCoreLosses += volumetricCoreLosses;
            }
        }
        averageVolumetricCoreLosses /= inputs.get_operating_points().size();
        evaluations.push_back({coreMaterial, averageVolumetricCoreLosses}); // NEW-BUG-FIX: was pushing 0, discarding computed losses
    }

    stable_sort((evaluations).begin(), (evaluations).end(), [](const std::pair<CoreMaterial, double>& b1, const std::pair<CoreMaterial, double>& b2) {
        return b1.second < b2.second;
    });

    for (size_t magneticIndex = 0; magneticIndex < (*magneticsWithScoring).size(); ++magneticIndex) {
        for (size_t i = 0; i < numberCoreMaterialsTouse; ++i){
            auto [magnetic, scoring] = (*magneticsWithScoring)[magneticIndex];
            magnetic.get_mutable_core().set_material(evaluations[i].first);
            magnetic.get_mutable_core().set_name(evaluations[i].first.get_name() + " " + magnetic.get_core().get_name().value_or("unnamed"));
            if (magnetic.get_manufacturer_info()) {
                auto manufacturerInfo = magnetic.get_manufacturer_info().value();
                manufacturerInfo.set_reference(evaluations[i].first.get_name() + " " + magnetic.get_reference());
                magnetic.set_manufacturer_info(manufacturerInfo);
            }
            magneticsWithMaterials.push_back({magnetic, scoring});
        }
    }
    return magneticsWithMaterials;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_ferrite_materials_by_losses(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = 2;
    double magneticFluxDensityReference = 0.18;
    std::vector<std::pair<Magnetic, double>> magneticsWithMaterials;
    std::vector<CoreMaterial> coreMaterialsToEvaluate;
    auto coreMaterials = get_core_material_names(settings.get_preferred_core_material_ferrite_manufacturer());
    for (auto coreMaterial : coreMaterials) {
        auto application = Core::guess_material_application(coreMaterial);
        if (application == _application) {
            coreMaterialsToEvaluate.push_back(Core::resolve_material(coreMaterial));
        }
    }
    std::vector<CoreMaterial> coreMaterialsToUse;
    std::vector<std::pair<CoreMaterial, double>> evaluations;

    double temperature = 0; 
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }


    SignalDescriptor magneticFluxDensity;
    Processed processed;
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);
    processed.set_peak(magneticFluxDensityReference);
    processed.set_peak_to_peak(2 * magneticFluxDensityReference);
    magneticFluxDensity.set_processed(processed);
    OperatingPointExcitation operatingPointExcitation;
    operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);
    operatingPointExcitation.set_frequency(1);

    auto coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
    auto coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
    for (auto coreMaterial : coreMaterialsToEvaluate) {
        double averageVolumetricCoreLosses = 0;
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex){
            auto frequency = inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[0].get_frequency();
            Inputs::scale_time_to_frequency(operatingPointExcitation, frequency, false, false);
            auto coreLossesMethods = Core::get_available_core_losses_methods(coreMaterial);

            if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
                double volumetricCoreLosses = coreLossesModelSteinmetz->get_core_volumetric_losses(coreMaterial, operatingPointExcitation, temperature);
                averageVolumetricCoreLosses += volumetricCoreLosses;
            }
            else {
                double volumetricCoreLosses = coreLossesModelProprietary->get_core_volumetric_losses(coreMaterial, operatingPointExcitation, temperature);
                averageVolumetricCoreLosses += volumetricCoreLosses;
            }
        }
        averageVolumetricCoreLosses /= inputs.get_operating_points().size();
        evaluations.push_back({coreMaterial, averageVolumetricCoreLosses});
    }

    stable_sort((evaluations).begin(), (evaluations).end(), [](const std::pair<CoreMaterial, double>& b1, const std::pair<CoreMaterial, double>& b2) {
        return b1.second < b2.second;
    });

    for (size_t magneticIndex = 0; magneticIndex < (*magneticsWithScoring).size(); ++magneticIndex) {
        for (size_t i = 0; i < numberCoreMaterialsTouse; ++i){
            auto [magnetic, scoring] = (*magneticsWithScoring)[magneticIndex];
            if (magnetic.get_mutable_core().get_material_name() != "Dummy") {
                magneticsWithMaterials.push_back({magnetic, scoring});
                continue;
            }
            magnetic.get_mutable_core().set_material(evaluations[i].first);
            magnetic.get_mutable_core().set_name(evaluations[i].first.get_name() + " " + (magnetic.get_core().get_name().value_or("unnamed")));
            if (magnetic.get_manufacturer_info()) {
                auto manufacturerInfo = magnetic.get_manufacturer_info().value();
                manufacturerInfo.set_reference(evaluations[i].first.get_name() + " " + magnetic.get_reference());
                magnetic.set_manufacturer_info(manufacturerInfo);
            }
            magneticsWithMaterials.push_back({magnetic, scoring});
        }
    }

    return magneticsWithMaterials;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_ferrite_materials_by_impedance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = 2;
    std::vector<std::pair<Magnetic, double>> magneticsWithMaterials;
    std::vector<CoreMaterial> coreMaterialsToEvaluate;
    auto coreMaterials = get_core_material_names(settings.get_preferred_core_material_ferrite_manufacturer());
    for (auto coreMaterial : coreMaterials) {
        auto application = Core::guess_material_application(coreMaterial);
        if (application == _application) {
            coreMaterialsToEvaluate.push_back(Core::resolve_material(coreMaterial));
        }
    }
    std::vector<CoreMaterial> coreMaterialsToUse;
    std::vector<std::pair<CoreMaterial, double>> evaluations;

    if (!inputs.get_design_requirements().get_minimum_impedance()) {
        throw std::invalid_argument("Missing impedance requirement");
    }

    auto minimumImpedanceRequirement = inputs.get_design_requirements().get_minimum_impedance().value();

    ComplexPermeability complexPermeabilityModel;
    for (auto coreMaterial : coreMaterialsToEvaluate) {
        double totalComplexPermeability = 0;
        for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
            auto frequency = impedanceAtFrequency.get_frequency();
            auto [complexPermeabilityRealPart, complexPermeabilityImaginaryPart] = complexPermeabilityModel.get_complex_permeability(coreMaterial, frequency);
            totalComplexPermeability += hypot(complexPermeabilityRealPart, complexPermeabilityImaginaryPart);
        }
        evaluations.push_back({coreMaterial, totalComplexPermeability});
    }

    stable_sort((evaluations).begin(), (evaluations).end(), [](const std::pair<CoreMaterial, double>& b1, const std::pair<CoreMaterial, double>& b2) {
        return b1.second > b2.second;
    });

    for (size_t magneticIndex = 0; magneticIndex < (*magneticsWithScoring).size(); ++magneticIndex) {
        for (size_t i = 0; i < numberCoreMaterialsTouse; ++i){
            auto [magnetic, scoring] = (*magneticsWithScoring)[magneticIndex];
            // if (magnetic.get_mutable_core().get_shape_name() != "T 90/57/23") {
            //     continue;
            // }
            // if (magnetic.get_mutable_core().get_number_stacks() > 1) {
            //     continue;
            // }
            if (magnetic.get_mutable_core().get_material_name() != "Dummy") {
                magneticsWithMaterials.push_back({magnetic, scoring});
                continue;
            }
            magnetic.get_mutable_core().set_material(evaluations[i].first);
            magnetic.get_mutable_core().set_name(evaluations[i].first.get_name() + " " + magnetic.get_core().get_name().value_or("unnamed"));
            if (magnetic.get_manufacturer_info()) {
                auto manufacturerInfo = magnetic.get_manufacturer_info().value();
                manufacturerInfo.set_reference(evaluations[i].first.get_name() + " " + magnetic.get_reference());
                magnetic.set_manufacturer_info(manufacturerInfo);
            }
            magneticsWithMaterials.push_back({magnetic, scoring});
        }
    }
    return magneticsWithMaterials;
}

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
        // Detect the "<prefix> Half 1" / "<prefix> Half 2" naming pattern.
        const std::string h1 = " Half 1";
        const std::string h2 = " Half 2";
        if (a.size() <= h1.size() || b.size() <= h2.size()) return false;
        if (a.compare(a.size() - h1.size(), h1.size(), h1) != 0) return false;
        if (b.compare(b.size() - h2.size(), h2.size(), h2) != 0) return false;
        return a.substr(0, a.size() - h1.size()) == b.substr(0, b.size() - h2.size());
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
    double temperature = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }

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
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Area Product filter.", "CoreAdviser");

    if (settings.get_core_adviser_enable_intermediate_pruning() && magneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
        magneticsWithScoring.resize(maximumMagneticsAfterFiltering); // F10 FIX: resize instead of copy-construct
        logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " after culling by the score on the first filter.", "CoreAdviser");
    }

    magneticsWithScoring = filterEnergyStored.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);  // Fixed weight: pre-filtering criterion, not efficiency scoring
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Energy Stored filter.", "CoreAdviser");

    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

    // Add saturation filter to reject cores that exceed magnetic flux density saturation
    MagneticCoreFilterSaturation filterSaturationAvailable;
    filterSaturationAvailable.set_scorings(&_scorings);
    filterSaturationAvailable.set_filter_configuration(&_filterConfiguration);
    magneticsWithScoring = filterSaturationAvailable.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Saturation filter.", "CoreAdviser");

    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::COST], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Cost filter.", "CoreAdviser");

    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::DIMENSIONS], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Dimensions filter.", "CoreAdviser");

    // Prune to top candidates by accumulated Cost+Dimensions score before the expensive
    // Loss filter (which sweeps N per core). Cores ranking poorly on cost+size are not
    // going to be promoted by losses alone.
    {
        const size_t preLossCap = std::max<size_t>(maximumNumberResults * 5, 50);
        if (settings.get_core_adviser_enable_intermediate_pruning() && magneticsWithScoring.size() > preLossCap) {
            magneticsWithScoring.resize(preLossCap);
            logEntry("Pruned to " + std::to_string(magneticsWithScoring.size()) + " magnetics before the Core Losses filter.", "CoreAdviser");
        }
    }

    magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Core Losses filter.", "CoreAdviser");

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
        
        logEntry("After retry with relaxed constraints: " + std::to_string(magneticsWithScoring.size()) + " magnetics", "CoreAdviser");
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

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring = *magnetics;

    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, 1, true);

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
    logEntry("Starting with " + std::to_string(magneticsWithScoring.size()) + " magnetics", "CoreAdviser");
    std::cout << "[CoreAdviser] Starting with " << magneticsWithScoring.size() << " magnetics" << std::endl;

    bool usingPowderCores = should_include_powder(inputs);

    // ========================================================================
    // STEP 1: Area Product filter on all cores
    // ========================================================================
    magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    logEntry("After AreaProduct: " + std::to_string(magneticsWithScoring.size()), "CoreAdviser");
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
        logEntry("After gapping ferrite: " + std::to_string(ferriteCores.size()), "CoreAdviser");
        std::cout << "[CoreAdviser] After gapping: " << ferriteCores.size() << std::endl;
        
        // Filter by fringing factor
        ferriteCores = filterFringingFactor.filter_magnetics(&ferriteCores, inputs, 1, true);
        logEntry("After FringingFactor: " + std::to_string(ferriteCores.size()), "CoreAdviser");
        std::cout << "[CoreAdviser] After FringingFactor: " << ferriteCores.size() << std::endl;
        
        // Filter by dimensions
        ferriteCores = filterDimensions.filter_magnetics(&ferriteCores, inputs, 1, true);
        logEntry("After Dimensions (ferrite): " + std::to_string(ferriteCores.size()), "CoreAdviser");
        std::cout << "[CoreAdviser] After Dimensions: " << ferriteCores.size() << std::endl;
        
        // Assign concrete ferrite materials
        ferriteCores = add_ferrite_materials_by_losses(&ferriteCores, inputs);
        logEntry("After materials (ferrite): " + std::to_string(ferriteCores.size()), "CoreAdviser");
        std::cout << "[CoreAdviser] After materials: " << ferriteCores.size() << std::endl;
        
        // Calculate turns
        add_initial_turns_by_inductance(&ferriteCores, inputs);
        
        // Filter by inductance
        ferriteCores = filterMagneticInductance.filter_magnetics(&ferriteCores, inputs, 0.1, true);
        logEntry("After Inductance (ferrite): " + std::to_string(ferriteCores.size()), "CoreAdviser");
        std::cout << "[CoreAdviser] After Inductance: " << ferriteCores.size() << std::endl;
        
        // Filter by saturation
        ferriteCores = filterSaturation.filter_magnetics(&ferriteCores, inputs, 1, true);
        logEntry("After Saturation (ferrite): " + std::to_string(ferriteCores.size()), "CoreAdviser");
        std::cout << "[CoreAdviser] After Saturation: " << ferriteCores.size() << std::endl;

        // Prune to top candidates by accumulated score before the expensive Loss filter
        // (which sweeps N per core). Cores already losing on cost+dimensions+inductance
        // are not going to be promoted by losses alone.
        {
            const size_t preLossCap = std::max<size_t>(maximumNumberResults * 5, 50);
            if (settings.get_core_adviser_enable_intermediate_pruning() && ferriteCores.size() > preLossCap) {
                ferriteCores.resize(preLossCap);
                logEntry("Pruned ferrite to " + std::to_string(ferriteCores.size()) + " before Losses", "CoreAdviser");
            }
        }

        // Filter by losses
        ferriteCores = filterLosses.filter_magnetics(&ferriteCores, inputs, 1, true);
        logEntry("After Losses (ferrite): " + std::to_string(ferriteCores.size()), "CoreAdviser");
        std::cout << "[CoreAdviser] After Losses: " << ferriteCores.size() << std::endl;
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
            logEntry("After powder materials: " + std::to_string(powderCores.size()), "CoreAdviser");
            
            // Filter by energy stored
            powderCores = filterEnergyStored.filter_magnetics(&powderCores, inputs, 1, true);
            logEntry("After EnergyStored (powder): " + std::to_string(powderCores.size()), "CoreAdviser");
            
            // Filter by dimensions
            powderCores = filterDimensions.filter_magnetics(&powderCores, inputs, 1, true);
            logEntry("After Dimensions (powder): " + std::to_string(powderCores.size()), "CoreAdviser");
            
            // Calculate turns
            add_initial_turns_by_inductance(&powderCores, inputs);
            
            // Filter by inductance
            powderCores = filterMagneticInductance.filter_magnetics(&powderCores, inputs, 0.1, true);
            logEntry("After Inductance (powder): " + std::to_string(powderCores.size()), "CoreAdviser");
            
            // Filter by saturation
            powderCores = filterSaturation.filter_magnetics(&powderCores, inputs, 1, true);
            logEntry("After Saturation (powder): " + std::to_string(powderCores.size()), "CoreAdviser");

            // Prune to top candidates by accumulated score before the expensive Loss filter.
            {
                const size_t preLossCap = std::max<size_t>(maximumNumberResults * 5, 50);
                if (settings.get_core_adviser_enable_intermediate_pruning() && powderCores.size() > preLossCap) {
                    powderCores.resize(preLossCap);
                    logEntry("Pruned powder to " + std::to_string(powderCores.size()) + " before Losses", "CoreAdviser");
                }
            }

            // Filter by losses
            powderCores = filterLosses.filter_magnetics(&powderCores, inputs, 1, true);
            logEntry("After Losses (powder): " + std::to_string(powderCores.size()), "CoreAdviser");
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

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring = *magnetics;

    // if (magneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
    //     magneticsWithScoring.resize(maximumMagneticsAfterFiltering); // F10 FIX: resize instead of copy-construct
    // }


    magneticsWithScoring = add_ferrite_materials_by_impedance(&magneticsWithScoring, inputs);

    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);
    // magneticsWithScoring = add_initial_turns_by_impedance(magneticsWithScoring, inputs);

    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterMagneticInductance.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, 1, true);

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

}