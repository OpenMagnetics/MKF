#include "advisers/CoreAdviser.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"  // for is_energy_storing_topology()
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

namespace {
// Phase 3 (F10): compact replacement for the 20+ verbatim
//   logEntry("There are " + std::to_string(c.size()) + " magnetics after the X filter.", "CoreAdviser");
// idioms scattered across the two pipeline functions. Format:
//   "After <stage>: <n>"
// (the older "There are N magnetics after the X filter." form is dropped to one
// canonical compact format — no test asserts on log content.)
inline void log_stage(const std::string& stage, size_t count) {
    logEntry("After " + stage + ": " + std::to_string(count), "CoreAdviser");
}
inline void log_pruned(const std::string& stage, size_t count) {
    logEntry("Pruned to " + std::to_string(count) + " before " + stage, "CoreAdviser");
}

// Phase 3 (F8): the sinusoidal-excitation builder + Steinmetz/Proprietary
// model-pair factory pattern appeared verbatim in add_ferrite_materials_by_losses
// and add_ferrite_materials_by_impedance. Centralised here.
inline OperatingPointExcitation make_sinusoidal_excitation(double peak, double offset, double frequency) {
    SignalDescriptor magneticFluxDensity;
    Processed processed;
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(offset);
    processed.set_peak(peak);
    processed.set_peak_to_peak(2 * peak);
    magneticFluxDensity.set_processed(processed);
    OperatingPointExcitation excitation;
    excitation.set_magnetic_flux_density(magneticFluxDensity);
    excitation.set_frequency(frequency);
    return excitation;
}

struct CoreLossesModelPair {
    std::shared_ptr<CoreLossesModel> steinmetz;
    std::shared_ptr<CoreLossesModel> proprietary;
};
inline CoreLossesModelPair make_default_core_losses_model_pair() {
    return {
        CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}})),
        CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}))
    };
}
} // namespace

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

Coil get_dummy_coil(Inputs inputs) {
    double temperature = inputs.get_maximum_temperature();
    double frequency = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        frequency = std::max(frequency, Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency());
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
    coil.set_bobbin(DUMMY_SENTINEL_NAME);
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
        // Phase 1 fix: coreShapeDatabase is keyed by canonical-name AND by
        // each alias (see Utils.cpp:255-260), so a naive iteration yields the
        // same CoreShape (1 + #aliases) times. That was the root cause of
        // CoreAdviser STANDARD_CORES × POWER returning duplicate top-N
        // entries with identical names and scores (e.g. three copies of
        // "95 PQ 27/15 gapped 0.36 mm" in slots 0-2). Dedupe by canonical
        // shape.get_name() before feeding the dataset builder.
        std::set<std::string> seenShapeNames;
        for (auto& [key, shape] : coreShapeDatabase) {
            auto canonical = shape.get_name().value_or("");
            if (canonical.empty()) {
                throw std::runtime_error("CoreShape in coreShapeDatabase has no name");
            }
            if (!seenShapeNames.insert(canonical).second) {
                continue;
            }
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
    double maximumHeight = std::numeric_limits<double>::infinity();
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
            // Common-mode chokes MUST be wound on a single closed magnetic path so the
            // common-mode flux of both windings adds in the core: only toroidal cores
            // give the tight coupling required. Differential-mode chokes have no such
            // constraint — they're regular inductors and can use any core type
            // (powdered-iron toroids, gapped ferrite E/EE/PQ, etc.). Apply the
            // toroidal-only restriction only to CMC; for DMC and other suppression
            // topologies, honour the user's includeToroidalCores / includeConcentricCores
            // settings.
            auto topology = inputs.get_design_requirements().get_topology();
            bool isCommonModeChoke = topology.has_value() && topology.value() == Topologies::COMMON_MODE_CHOKE;
            if (isCommonModeChoke) {
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
    double maximumHeight = std::numeric_limits<double>::infinity();
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
            // See comment in create_magnetic_dataset(): only CMC requires toroidal-only;
            // DMC and other suppression topologies should honour user settings.
            auto topology = inputs.get_design_requirements().get_topology();
            bool isCommonModeChoke = topology.has_value() && topology.value() == Topologies::COMMON_MODE_CHOKE;
            if (isCommonModeChoke) {
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
            // Capture the un-stacked base name once so stack variants don't
            // accumulate suffixes across loop iterations.
            std::string baseName = core.get_name().value_or("unnamed");
            for (size_t i = 0; i < defaults.coreAdviserMaximumNumberStacks; ++i) {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                // process_data() resets processed description to base values, then calls scale_to_stacks internally
                core.process_data();
                core.process_gap(); // CA-OPT-2 FIX: reprocess gap data after stacking (was commented out)
                // Phase 1 fix: previously only manufacturer_info.reference got
                // the " N stacks" suffix; core.name stayed identical across
                // stack variants. add_gapping_standard_cores then produced
                // identical display names (e.g. three copies of
                // "95 PQ 27/15 gapped 0.36 mm") for 1/2/3-stack variants.
                // Put the stack suffix in the core name too so downstream
                // names stay unique.
                std::string variantName = baseName;
                if (i != 0) {
                    variantName += " " + std::to_string(1 + i) + " stacks";
                }
                core.set_name(variantName);
                magnetic.set_core(core);
                MagneticManufacturerInfo magneticManufacturerInfo;
                magneticManufacturerInfo.set_reference(variantName);
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
    double maximumHeight = std::numeric_limits<double>::infinity();
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
    double temperature = inputs.get_maximum_temperature();
    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){
        Core core = (*magneticsWithScoring)[i].first.get_core();
        auto coreMaterial = core.resolve_material();

        // Phase 1: per-candidate alternatives lookup is annotational metadata.
        // If the cross-reference can't proceed (e.g. this candidate's material
        // has no Steinmetz coefficients at this temperature, so it can't be
        // used as a reference for ranking others), record an empty
        // alternatives list and continue. Named, logged, bounded scope.
        std::vector<std::string> coreMaterialAlternatives;
        try {
            auto alternatives = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, temperature);
            for (auto [alternativeCoreMaterial, scoring] : alternatives) {
                coreMaterialAlternatives.push_back(alternativeCoreMaterial.get_name());
            }
        }
        catch (const std::exception& e) {
            logEntry(std::string("Skipping alternative-materials lookup for candidate with material '")
                         + coreMaterial.get_name() + "': " + e.what(),
                     "CoreAdviser", 2);
        }
        coreMaterial.set_alternatives(coreMaterialAlternatives);
        core.set_material(coreMaterial);

        (*magneticsWithScoring)[i].first.set_core(core);
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
    size_t numberCoreMaterialsTouse = defaults.coreAdviserAlternativeMaterialsNumberToUse;
    double magneticFluxDensityReference = defaults.coreAdviserMagneticFluxDensityReferenceAlternative;
    std::vector<std::pair<Magnetic, double>> magneticsWithMaterials;
    std::vector<CoreMaterial> coreMaterialsToEvaluate;
    auto coreMaterials = get_core_material_names(settings.get_preferred_core_material_powder_manufacturer());
    for (auto coreMaterial : coreMaterials) {
        auto resolved = Core::resolve_material(coreMaterial);
        if (Core::check_material_application(resolved, _application)) {
            coreMaterialsToEvaluate.push_back(resolved);
        }
    }
    std::vector<CoreMaterial> coreMaterialsToUse;
    std::vector<std::pair<CoreMaterial, double>> evaluations;

    double temperature = inputs.get_maximum_temperature();
    double maximumCurrentDcBias = inputs.get_maximum_current_dc_bias();

    // Phase 3 (F8): sinusoidal excitation + loss-model pair are shared
    // with add_ferrite_materials_by_impedance.
    OperatingPointExcitation operatingPointExcitation = make_sinusoidal_excitation(
        magneticFluxDensityReference, maximumCurrentDcBias, 1);
    auto [coreLossesModelSteinmetz, coreLossesModelProprietary] = make_default_core_losses_model_pair();
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

// Phase 3 (F4): shared by add_ferrite_materials_by_losses /
// add_ferrite_materials_by_impedance — gather ferrite candidates filtered
// by current application; small but called twice.
std::vector<CoreMaterial> CoreAdviser::gather_ferrite_materials_for_application() const {
    std::vector<CoreMaterial> coreMaterialsToEvaluate;
    auto coreMaterials = get_core_material_names(settings.get_preferred_core_material_ferrite_manufacturer());
    for (auto coreMaterial : coreMaterials) {
        auto resolved = Core::resolve_material(coreMaterial);
        if (Core::check_material_application(resolved, _application)) {
            coreMaterialsToEvaluate.push_back(resolved);
        }
    }
    return coreMaterialsToEvaluate;
}

// Phase 3 (F4): the dummy-material fan-out block was duplicated verbatim
// between add_ferrite_materials_by_losses and add_ferrite_materials_by_impedance
// (24 lines × 2). Pass-through for non-Dummy candidates (the Phase 1 fix
// for the duplicate-entries bug) is preserved here.
std::vector<std::pair<Magnetic, double>> CoreAdviser::fan_out_dummy_into_top_materials(
    std::vector<std::pair<Magnetic, double>> *magneticsWithScoring,
    const std::vector<std::pair<CoreMaterial, double>>& sortedMaterialEvaluations,
    size_t numberCoreMaterialsTouse) {
    std::vector<std::pair<Magnetic, double>> magneticsWithMaterials;
    for (size_t magneticIndex = 0; magneticIndex < (*magneticsWithScoring).size(); ++magneticIndex) {
        auto [magnetic, scoring] = (*magneticsWithScoring)[magneticIndex];
        if (magnetic.get_mutable_core().get_material_name() != DUMMY_SENTINEL_NAME) {
            // Already has a concrete material — pass through ONCE.
            // Phase 1 fix: previously the inner i-loop also ran in this branch
            // (with `continue`), so non-Dummy magnetics were pushed
            // numberCoreMaterialsTouse (=2) times, producing the duplicate
            // entries in CoreAdviser STANDARD_CORES × POWER top-N output.
            magneticsWithMaterials.push_back({magnetic, scoring});
            continue;
        }
        // Dummy material — fan out into the top-N candidate ferrite materials.
        for (size_t i = 0; i < numberCoreMaterialsTouse; ++i){
            auto magneticCopy = magnetic;
            magneticCopy.get_mutable_core().set_material(sortedMaterialEvaluations[i].first);
            magneticCopy.get_mutable_core().set_name(sortedMaterialEvaluations[i].first.get_name() + " " + magneticCopy.get_core().get_name().value_or("unnamed"));
            if (magneticCopy.get_manufacturer_info()) {
                auto manufacturerInfo = magneticCopy.get_manufacturer_info().value();
                manufacturerInfo.set_reference(sortedMaterialEvaluations[i].first.get_name() + " " + magneticCopy.get_reference());
                magneticCopy.set_manufacturer_info(manufacturerInfo);
            }
            magneticsWithMaterials.push_back({magneticCopy, scoring});
        }
    }
    return magneticsWithMaterials;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_ferrite_materials_by_losses(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = defaults.coreAdviserFanOutMaterialsNumberToUse;
    double magneticFluxDensityReference = defaults.coreAdviserMagneticFluxDensityReferenceAlternative;
    auto coreMaterialsToEvaluate = gather_ferrite_materials_for_application();
    std::vector<CoreMaterial> coreMaterialsToUse;
    std::vector<std::pair<CoreMaterial, double>> evaluations;

    double temperature = inputs.get_maximum_temperature();


    // Phase 3 (F8): sinusoidal excitation + loss-model pair are shared
    // with add_ferrite_materials_by_losses.
    OperatingPointExcitation operatingPointExcitation = make_sinusoidal_excitation(
        magneticFluxDensityReference, 0, 1);
    auto [coreLossesModelSteinmetz, coreLossesModelProprietary] = make_default_core_losses_model_pair();
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

    return fan_out_dummy_into_top_materials(magneticsWithScoring, evaluations, numberCoreMaterialsTouse);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_ferrite_materials_by_impedance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = 2;
    auto coreMaterialsToEvaluate = gather_ferrite_materials_for_application();
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
        return b1.second < b2.second;
    });

    return fan_out_dummy_into_top_materials(magneticsWithScoring, evaluations, numberCoreMaterialsTouse);
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

}