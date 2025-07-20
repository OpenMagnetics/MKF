#include "advisers/CoreAdviser.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "advisers/MagneticFilter.h"
#include "constructive_models/Bobbin.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/ComplexPermeability.h"
#include "support/CoilMesher.h"
#include "constructive_models/NumberTurns.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Bobbin.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <string>
#include <vector>
#include <execution>
#include <magic_enum_utility.hpp>
#include <list>
#include <cmrc/cmrc.hpp>

CMRC_DECLARE(data);


namespace OpenMagnetics {


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

    // sort((*magneticsWithScoring).begin(), (*magneticsWithScoring).end(), [](std::pair<Magnetic, double>& b1, std::pair<Magnetic, double>& b2) {
    //     return b1.second > b2.second;
    // });

    return normalizedScorings;
}

void sort_magnetics_by_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring) {
    sort((*magneticsWithScoring).begin(), (*magneticsWithScoring).end(), [](std::pair<Magnetic, double>& b1, std::pair<Magnetic, double>& b2) {
        return b1.second > b2.second;
    });
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
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
    _filter = MagneticFilterCoreAndDcLosses(inputs, models);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterLosses::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    auto coilDelimitAndCompactOld = settings->get_coil_delimit_and_compact();
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
        settings->set_coil_delimit_and_compact(coilDelimitAndCompactOld);
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }

    settings->set_coil_delimit_and_compact(coilDelimitAndCompactOld);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterDimensions::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMagnetics;
    }
    std::vector<double> newScoring;

    for (size_t magneticIndex = 0; magneticIndex < (*unfilteredMagnetics).size(); ++magneticIndex){
        Magnetic magnetic = (*unfilteredMagnetics)[magneticIndex].first;

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, nullptr);

        newScoring.push_back(scoring);
    }

    if ((*unfilteredMagnetics).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string((*unfilteredMagnetics).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredMagnetics).size() > 0) {
        auto normalizedScoring = normalize_scoring(unfilteredMagnetics, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::DIMENSIONS]);
        for (size_t i = 0; i < (*unfilteredMagnetics).size(); ++i) {
            add_scoring((*unfilteredMagnetics)[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::DIMENSIONS, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(unfilteredMagnetics);
    }
    return (*unfilteredMagnetics);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterMinimumImpedance::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    auto coilDelimitAndCompactOld = settings->get_coil_delimit_and_compact();
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }


    settings->set_coil_delimit_and_compact(coilDelimitAndCompactOld);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::MagneticCoreFilterMagneticInductance::filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight, bool firstFilter) {
    auto coilDelimitAndCompactOld = settings->get_coil_delimit_and_compact();
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        auto normalizedScoring = normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
        for (size_t i = 0; i < filteredMagneticsWithScoring.size(); ++i) {
            add_scoring(filteredMagneticsWithScoring[i].first.get_reference(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, normalizedScoring[i]);
        }
        sort_magnetics_by_scoring(&filteredMagneticsWithScoring);
    }


    settings->set_coil_delimit_and_compact(coilDelimitAndCompactOld);
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
    CoilFunctionalDescription primaryCoilFunctionalDescription;
    primaryCoilFunctionalDescription.set_isolation_side(IsolationSide::PRIMARY);
    primaryCoilFunctionalDescription.set_name("primary");
    primaryCoilFunctionalDescription.set_number_parallels(1);
    primaryCoilFunctionalDescription.set_number_turns(1);
    primaryCoilFunctionalDescription.set_wire(wire);

    Coil coil;
    coil.set_bobbin("Dummy");
    coil.set_functional_description({primaryCoilFunctionalDescription});
    return coil;
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults){
    if (coreDatabase.empty()) {
        load_cores();
    }
    return get_advised_core(inputs, weights, &coreDatabase, maximumNumberResults);
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

    sort(results.begin(), results.end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
        return b1.second > b2.second;
    });

    if (results.size() > maximumNumberResults) {
        results = std::vector<std::pair<Mas, double>>(results.begin(), results.end() - (results.size() - maximumNumberResults));
    }
    return results;
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults){
    _weights = weights;
 
    size_t maximumMagneticsAfterFiltering = defaults.coreAdviserMaximumMagneticsAfterFiltering;
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

        auto globalIncludeStacks = settings->get_core_adviser_include_stacks();
        if (globalIncludeStacks) {
            expand_magnetic_dataset_with_stacks(inputs, cores, &magnetics);
        }

        logEntry("First attempt produced not enough results, so now we are searching again with " + std::to_string(magnetics.size()) + " magnetics, including up to " + std::to_string(defaults.coreAdviserMaximumNumberStacks) + " cores stacked when possible.", "CoreAdviser");
        maximumMagneticsAfterFiltering = magnetics.size();
        filteredMagnetics = filter_available_cores_power_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
    else {
        filteredMagnetics = filter_available_cores_suppression_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);

    }
    return filteredMagnetics;
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<CoreShape>* shapes, std::vector<CoreMaterial>* materials, size_t maximumNumberResults) {
    auto globalIncludeStacks = settings->get_core_adviser_include_stacks();
    auto magnetics = create_magnetic_dataset(inputs, shapes, false);

    size_t maximumMagneticsAfterFiltering = defaults.coreAdviserMaximumMagneticsAfterFiltering;
    if (get_application() == Application::POWER) {
        return filter_standard_cores_power_application(&magnetics, inputs, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
    else {
        return filter_standard_cores_interference_suppression_application(&magnetics, inputs, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::create_magnetic_dataset(Inputs inputs, std::vector<Core>* cores, bool includeStacks) {
    std::vector<std::pair<Magnetic, double>> magnetics;
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings->get_use_toroidal_cores();
    auto includeConcentricCores = settings->get_use_concentric_cores();
    auto globalIncludeStacks = settings->get_core_adviser_include_stacks();
    auto globalIncludeDistributedGaps = settings->get_core_adviser_include_distributed_gaps();
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
                core.scale_to_stacks(1 + i);
                // core.process_gap();
                magnetic.set_core(core);
                MagneticManufacturerInfo magneticmanufacturerinfo;
                if (i != 0) {
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticmanufacturerinfo.set_reference(core.get_name().value());
                }
                magnetic.set_manufacturer_info(magneticmanufacturerinfo);
                magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }
        else {
            magnetic.set_core(core);
            MagneticManufacturerInfo magneticmanufacturerinfo;
            magneticmanufacturerinfo.set_reference(core.get_name().value());
            magnetic.set_manufacturer_info(magneticmanufacturerinfo);
            magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
        }
    }

    return magnetics;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::create_magnetic_dataset(Inputs inputs, std::vector<CoreShape>* shapes, bool includeStacks) {
    std::vector<std::pair<Magnetic, double>> magnetics;
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings->get_use_toroidal_cores();
    auto includeConcentricCores = settings->get_use_concentric_cores();
    auto globalIncludeStacks = settings->get_core_adviser_include_stacks();
    auto globalIncludeDistributedGaps = settings->get_core_adviser_include_distributed_gaps();
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
                core.scale_to_stacks(1 + i);
                // core.process_gap();
                magnetic.set_core(core);
                MagneticManufacturerInfo magneticManufacturerInfo;
                if (i!=0) {
                    magneticManufacturerInfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticManufacturerInfo.set_reference(core.get_name().value());
                }
                magnetic.set_manufacturer_info(magneticManufacturerInfo);
                magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }
        else {
            magnetic.set_core(core);
            MagneticManufacturerInfo magneticManufacturerInfo;
            magneticManufacturerInfo.set_reference(core.get_name().value());
            magnetic.set_manufacturer_info(magneticManufacturerInfo);
            magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
        }
    }

    return magnetics;
}

void CoreAdviser::expand_magnetic_dataset_with_stacks(Inputs inputs, std::vector<Core>* cores, std::vector<std::pair<Magnetic, double>>* magnetics) {
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings->get_use_toroidal_cores();
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
                core.scale_to_stacks(1 + i);
                // core.process_gap();
                MagneticManufacturerInfo magneticManufacturerInfo;
                if (i!=0) {
                    std::string name = core.get_name().value();
                    name = std::regex_replace(std::string(name), std::regex(" [0-9] stacks"), "");
                    core.set_name(name + " " + std::to_string(1 + i) + " stacks" );
                    magneticManufacturerInfo.set_reference(core.get_name().value());
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
    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){

        Core core = (*magneticsWithScoring)[i].first.get_core();
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        double initialNumberTurns = (*magneticsWithScoring)[i].first.get_coil().get_functional_description()[0].get_number_turns();

        if (initialNumberTurns == 1) {
            initialNumberTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs, DimensionalValues::MINIMUM);
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
        magnetic.get_mutable_coil().set_bobbin(Bobbin::create_quick_bobbin(core));

        double initialNumberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();

        try {
            initialNumberTurns = impedance.calculate_minimum_number_turns(magnetic, inputs);
            if (initialNumberTurns < 1) {
                continue;
            }
        }
        catch (const std::exception &exc) {
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
    CoreMaterialCrossReferencer coreMaterialCrossReferencer(std::map<std::string, std::string>{{"coreLosses", "STEINMETZ"}});
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
    if (inputs.get_design_requirements().get_magnetizing_inductance().get_minimum() &&
        !inputs.get_design_requirements().get_magnetizing_inductance().get_nominal() &&
        !inputs.get_design_requirements().get_magnetizing_inductance().get_maximum())
    {
        for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
        Core core = (*magneticsWithScoring)[i].first.get_core();
        core.set_name(core.get_name().value() + " ungapped");
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
            double gapLength = roundFloat(magneticEnergy.calculate_gap_length_by_magnetic_energy(core.get_gapping()[0], core.get_magnetic_flux_density_saturation(), requiredMagneticEnergy), 5);
            core.set_ground_gap(gapLength);
            core.process_gap();
            std::stringstream ss;
            ss << roundFloat(gapLength * 1000, 2);
            if (gapLength > 0) {
                core.set_name(core.get_name().value() + " gapped " + ss.str()  + " mm");
            }
            else {
                core.set_name(core.get_name().value() + " ungapped");
            }
        }

        (*magneticsWithScoring)[i].first.set_core(core);
    }
}

bool CoreAdviser::should_include_powder(Inputs inputs) {
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
    auto coreMaterials = get_material_names(settings->get_preferred_core_material_powder_manufacturer());
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

    auto coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));
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
        evaluations.push_back({coreMaterial, 0});
    }

    sort((evaluations).begin(), (evaluations).end(), [](std::pair<CoreMaterial, double>& b1, std::pair<CoreMaterial, double>& b2) {
        return b1.second < b2.second;
    });

    for (size_t magneticIndex = 0; magneticIndex < (*magneticsWithScoring).size(); ++magneticIndex) {
        for (size_t i = 0; i < numberCoreMaterialsTouse; ++i){
            auto [magnetic, scoring] = (*magneticsWithScoring)[magneticIndex];
            magnetic.get_mutable_core().set_material(evaluations[i].first);
            magnetic.get_mutable_core().set_name(evaluations[i].first.get_name() + " " + magnetic.get_core().get_name().value());
            auto manufacturerInfo = magnetic.get_manufacturer_info().value();
            manufacturerInfo.set_reference(evaluations[i].first.get_name() + " " + magnetic.get_reference());
            magnetic.set_manufacturer_info(manufacturerInfo);
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
    auto coreMaterials = get_material_names(settings->get_preferred_core_material_ferrite_manufacturer());
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

    auto coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
    auto coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));
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

    sort((evaluations).begin(), (evaluations).end(), [](std::pair<CoreMaterial, double>& b1, std::pair<CoreMaterial, double>& b2) {
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
            magnetic.get_mutable_core().set_name(evaluations[i].first.get_name() + " " + magnetic.get_core().get_name().value());
            auto manufacturerInfo = magnetic.get_manufacturer_info().value();
            manufacturerInfo.set_reference(evaluations[i].first.get_name() + " " + magnetic.get_reference());
            magnetic.set_manufacturer_info(manufacturerInfo);
            magneticsWithMaterials.push_back({magnetic, scoring});
        }
    }
    return magneticsWithMaterials;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_ferrite_materials_by_impedance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = 2;
    std::vector<std::pair<Magnetic, double>> magneticsWithMaterials;
    std::vector<CoreMaterial> coreMaterialsToEvaluate;
    auto coreMaterials = get_material_names(settings->get_preferred_core_material_ferrite_manufacturer());
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

    auto temperature = inputs.get_maximum_temperature();
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

    sort((evaluations).begin(), (evaluations).end(), [](std::pair<CoreMaterial, double>& b1, std::pair<CoreMaterial, double>& b2) {
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
            magnetic.get_mutable_core().set_name(evaluations[i].first.get_name() + " " + magnetic.get_core().get_name().value());
            auto manufacturerInfo = magnetic.get_manufacturer_info().value();
            manufacturerInfo.set_reference(evaluations[i].first.get_name() + " " + magnetic.get_reference());
            magnetic.set_manufacturer_info(manufacturerInfo);
            magneticsWithMaterials.push_back({magnetic, scoring});
        }
    }
    return magneticsWithMaterials;
}

void correct_windings(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    MagnetizingInductance magnetizingInductance;
    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){

        Coil coil = Coil((*magneticsWithScoring)[i].first.get_coil());
        NumberTurns numberTurns(coil.get_number_turns(0), inputs.get_design_requirements());
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();

        (*magneticsWithScoring)[i].first.get_mutable_coil().set_bobbin(Bobbin::create_quick_bobbin((*magneticsWithScoring)[i].first.get_core()));
        for (size_t windingIndex = 1; windingIndex < numberTurnsCombination.size(); ++windingIndex) {
            auto winding = coil.get_functional_description()[0];
            winding.set_number_turns(numberTurnsCombination[windingIndex]);
            winding.set_isolation_side(get_isolation_side_from_index(windingIndex));
            winding.set_name(std::string{magic_enum::enum_name(get_isolation_side_from_index(windingIndex))});
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
    auto manufacturerInfo = magnetic.get_manufacturer_info().value();
    manufacturerInfo.set_reference(magnetic.get_core().get_name().value());
    magnetic.set_manufacturer_info(manufacturerInfo);

    auto previousCoilDelimitAndCompact = settings->get_coil_delimit_and_compact();
    settings->set_coil_delimit_and_compact(false);
    magnetic.get_mutable_coil().fast_wind();
    settings->set_coil_delimit_and_compact(previousCoilDelimitAndCompact);

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

        outputs.set_magnetizing_inductance(magnetizingInductanceOutput);
        outputs.set_winding_losses(windingLossesOutput);
        outputs.set_core_losses(coreLossesOutput);

        mas.get_mutable_outputs().push_back(outputs);
    }
    mas.set_inputs(inputs);

    return mas;
}

Inputs pre_process_inputs(Inputs inputs) {
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex){
        auto excitation = Inputs::get_primary_excitation(inputs.get_mutable_operating_points()[operatingPointIndex]);
        if (!excitation.get_voltage()) {
            inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0].set_voltage(Inputs::calculate_induced_voltage(excitation, resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance())));
            Inputs::set_current_as_magnetizing_current(&inputs.get_mutable_operating_points()[operatingPointIndex]);
        }
        else if (!excitation.get_magnetizing_current()) {
            auto magnetizingCurrent = Inputs::calculate_magnetizing_current(excitation, resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance()), false);
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
    magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Area Product filter.", "CoreAdviser");

    if (magneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
        magneticsWithScoring = std::vector<std::pair<Magnetic, double>>(magneticsWithScoring.begin(), magneticsWithScoring.end() - (magneticsWithScoring.size() - maximumMagneticsAfterFiltering));
        logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " after culling by the score on the first filter.", "CoreAdviser");
    }

    magneticsWithScoring = filterEnergyStored.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Energy Stored filter.", "CoreAdviser");

    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::COST], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Cost filter.", "CoreAdviser");

    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, weights[CoreAdviserFilters::DIMENSIONS], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Dimensions filter.", "CoreAdviser");

    magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Core Losses filter.", "CoreAdviser");

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    if (magneticsWithScoring.size() > maximumNumberResults) {
        if (get_unique_core_shapes()) {
            magneticsWithScoring = cull_to_unique_core_shapes(magneticsWithScoring, maximumNumberResults);
        }
        else {
            magneticsWithScoring = std::vector<std::pair<Magnetic, double>>(magneticsWithScoring.begin(), magneticsWithScoring.end() - (magneticsWithScoring.size() - maximumNumberResults));
        }
    }

    correct_windings(&magneticsWithScoring, inputs);

    std::vector<std::pair<Mas, double>> masWithScoring;

    for (auto [magnetic, scoring] : magneticsWithScoring) {
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

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring = *magnetics;

    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " at the beginning.", "CoreAdviser");
    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " after adding turns for inductance.", "CoreAdviser");
    // magneticsWithScoring = add_initial_turns_by_impedance(magneticsWithScoring, inputs);
    // logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " after adding turns for impedance.", "CoreAdviser");

    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Minimum Impedance filter.", "CoreAdviser");

    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::COST], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Cost filter.", "CoreAdviser");

    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, weights[CoreAdviserFilters::DIMENSIONS], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Dimensions filter.", "CoreAdviser");

    magneticsWithScoring = filterMagneticInductance.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Magnetizing Inductance.", "CoreAdviser");

    // magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
    // logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the Core Losses filter.", "CoreAdviser");

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    if (magneticsWithScoring.size() > maximumNumberResults) {
        if (get_unique_core_shapes()) {
            magneticsWithScoring = cull_to_unique_core_shapes(magneticsWithScoring, maximumNumberResults);
        }
        else {
            magneticsWithScoring = std::vector<std::pair<Magnetic, double>>(magneticsWithScoring.begin(), magneticsWithScoring.end() - (magneticsWithScoring.size() - maximumNumberResults));
        }
    }

    correct_windings(&magneticsWithScoring, inputs);

    std::vector<std::pair<Mas, double>> masWithScoring;

    for (auto [magnetic, scoring] : magneticsWithScoring) {
        auto mas = post_process_core(magnetic, inputs);
        masWithScoring.push_back({mas, scoring});
    }

    return masWithScoring;
}

std::vector<std::pair<Mas, double>> CoreAdviser::filter_standard_cores_power_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
    inputs = pre_process_inputs(inputs);

    MagneticCoreFilterAreaProduct filterAreaProduct(inputs);
    MagneticCoreFilterEnergyStored filterEnergyStored(inputs, _models);
    MagneticCoreFilterCost filterCost(inputs);
    MagneticCoreFilterLosses filterLosses(inputs, _models);
    MagneticCoreFilterDimensions filterDimensions;
    MagneticCoreFilterMagneticInductance filterMagneticInductance;
    MagneticCoreFilterFringingFactor filterFringingFactor(inputs, _models);

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

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring;
    magneticsWithScoring = filterAreaProduct.filter_magnetics(magnetics, inputs, 1, true);

    bool usingPowderCores = should_include_powder(inputs);
    if (usingPowderCores) {
        maximumMagneticsAfterFiltering /= 2;
    }

    if (magneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
        magneticsWithScoring = std::vector<std::pair<Magnetic, double>>(magneticsWithScoring.begin(), magneticsWithScoring.end() - (magneticsWithScoring.size() - maximumMagneticsAfterFiltering));
    }

    std::vector<std::pair<Magnetic, double>> ungappedmagneticsWithScoring;
    if (should_include_powder(inputs)) {
        std::copy(magneticsWithScoring.begin(), magneticsWithScoring.end(), std::back_inserter(ungappedmagneticsWithScoring));
        ungappedmagneticsWithScoring = add_powder_materials(&ungappedmagneticsWithScoring, inputs);
    }

    add_gapping(&magneticsWithScoring, inputs);
    magneticsWithScoring = filterFringingFactor.filter_magnetics(&magneticsWithScoring, inputs, 1, true);

    if (usingPowderCores) {
        ungappedmagneticsWithScoring = filterEnergyStored.filter_magnetics(&ungappedmagneticsWithScoring, inputs, 1, true);

        std::copy(ungappedmagneticsWithScoring.begin(), ungappedmagneticsWithScoring.end(), std::back_inserter(magneticsWithScoring));
    }

    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, 1, true);

    magneticsWithScoring = add_ferrite_materials_by_losses(&magneticsWithScoring, inputs);
    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

    magneticsWithScoring = filterMagneticInductance.filter_magnetics(&magneticsWithScoring, inputs, 0.1, true);
    magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, 1, true);

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    if (magneticsWithScoring.size() > maximumNumberResults) {
        if (get_unique_core_shapes()) {
            magneticsWithScoring = cull_to_unique_core_shapes(magneticsWithScoring, maximumNumberResults);
        }
        else {
            magneticsWithScoring = std::vector<std::pair<Magnetic, double>>(magneticsWithScoring.begin(), magneticsWithScoring.end() - (magneticsWithScoring.size() - maximumNumberResults));
        }
    }

    correct_windings(&magneticsWithScoring, inputs);
    add_alternative_materials(&magneticsWithScoring, inputs);

    std::vector<std::pair<Mas, double>> masWithScoring;

    for (auto [magnetic, scoring] : magneticsWithScoring) {
        auto mas = post_process_core(magnetic, inputs);
        masWithScoring.push_back({mas, scoring});
    }

    return masWithScoring;
}

std::vector<std::pair<Mas, double>> CoreAdviser::filter_standard_cores_interference_suppression_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
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
    //     magneticsWithScoring = std::vector<std::pair<Magnetic, double>>(magneticsWithScoring.begin(), magneticsWithScoring.end() - (magneticsWithScoring.size() - maximumMagneticsAfterFiltering));
    // }


    magneticsWithScoring = add_ferrite_materials_by_impedance(&magneticsWithScoring, inputs);

    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);
    // magneticsWithScoring = add_initial_turns_by_impedance(magneticsWithScoring, inputs);

    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, 1, true);
    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, 1, true);
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
            magneticsWithScoring = std::vector<std::pair<Magnetic, double>>(magneticsWithScoring.begin(), magneticsWithScoring.end() - (magneticsWithScoring.size() - maximumNumberResults));
        }
    }

    correct_windings(&magneticsWithScoring, inputs);
    add_alternative_materials(&magneticsWithScoring, inputs);

    std::vector<std::pair<Mas, double>> masWithScoring;

    for (auto [magnetic, scoring] : magneticsWithScoring) {
        auto mas = post_process_core(magnetic, inputs);
        masWithScoring.push_back({mas, scoring});
    }

    return masWithScoring;
}

} // namespace OpenMagnetics
