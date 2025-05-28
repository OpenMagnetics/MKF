#include "advisers/CoreAdviser.h"
#include "advisers/MagneticFilter.h"
#include "constructive_models/Bobbin.h"
#include "physical_models/CoreTemperature.h"
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


void CoreAdviser::logEntry(std::string entry) {
    // std::cout << entry << std::endl;
    _log += entry + "\n";
}

std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> CoreAdviser::get_scorings(bool weighted){
    std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> swappedScorings;
    for (auto& [filter, aux] : _scorings) {
        auto filterConfiguration = _filterConfiguration[filter];

        auto normalizedScorings = OpenMagnetics::normalize_scoring(aux, weighted? _weights[filter] : 1, filterConfiguration);

        for (auto& [name, scoring] : aux) {
            swappedScorings[name][filter] = normalizedScorings[name];
        }
    }
    return swappedScorings;
}

void normalize_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, std::vector<double> newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    auto normalizedScorings = OpenMagnetics::normalize_scoring(newScoring, weight, filterConfiguration);

    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
        (*magneticsWithScoring)[i].second += normalizedScorings[i];
    }

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

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::AREA_PRODUCT)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT][magnetic.get_manufacturer_info().value().get_reference().value()]) {

                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(magneticIndex);
                }
                continue;
            }
        }

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::AREA_PRODUCT, scoring, firstFilter);
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
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT]);
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

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::ENERGY_STORED)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(magneticIndex);
                }
                continue;
            }
        }
        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::ENERGY_STORED, scoring, firstFilter);
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
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED]);
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

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::COST)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::COST].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::COST][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::COST][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(magneticIndex);
                }
                continue;
            }
        }
        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::COST, scoring, firstFilter);
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
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::COST]);
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

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::EFFICIENCY)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::EFFICIENCY].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::EFFICIENCY][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::EFFICIENCY][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(magneticIndex);
                }
                continue;
            }
        }
        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, scoring, firstFilter);
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

    if (filteredMagneticsWithScoring.size() == 0) {
        settings->set_coil_delimit_and_compact(coilDelimitAndCompactOld);
        return *unfilteredMagnetics;
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
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
        add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::DIMENSIONS, scoring, firstFilter);
    }

    if ((*unfilteredMagnetics).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMagnetics: " + std::to_string((*unfilteredMagnetics).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredMagnetics).size() > 0) {
        normalize_scoring(unfilteredMagnetics, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::DIMENSIONS]);
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

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(magneticIndex);
                }
                continue;
            }
        }

        auto [valid, scoring] = _filter.evaluate_magnetic(&magnetic, &inputs);

        if (valid) {
            (*unfilteredMagnetics)[magneticIndex].first = magnetic;
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE, scoring, firstFilter);
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
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE]);
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
        if (filter == CoreAdviserFilters::MINIMUM_IMPEDANCE) {
            weights[filter] = 0;
        }
        else {
            weights[filter] = 1.0;
        }
    });
    return get_advised_core(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        if (filter == CoreAdviserFilters::MINIMUM_IMPEDANCE) {
            weights[filter] = 0;
        }
        else {
            weights[filter] = 1.0;
        }
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
        if (filter == CoreAdviserFilters::MINIMUM_IMPEDANCE) {
            weights[filter] = 0;
        }
        else {
            weights[filter] = 1.0;
        }
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
 
    CoreAdviserFilters firstFilter = CoreAdviserFilters::AREA_PRODUCT;
    double maxWeight = 0;

    for (auto& pair : weights) {
        auto filter = pair.first;
        auto weight = pair.second;
        if (weight > maxWeight) {
            maxWeight = weight;
            firstFilter = filter;
        }
    }
    size_t maximumMagneticsAfterFiltering = defaults.coreAdviserMaximumMagneticsAfterFiltering;
    std::vector<std::pair<Magnetic, double>> magnetics;
    bool needToAddStacks = false;
    bool onlyMaterialsForFilters = weights[CoreAdviserFilters::MINIMUM_IMPEDANCE] > 0;

    if (firstFilter == CoreAdviserFilters::EFFICIENCY) {
        magnetics = create_magnetic_dataset(inputs, cores, true, onlyMaterialsForFilters);
        logEntry("We start the search with " + std::to_string(magnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.");
        std::string firstFilterString = std::string{magic_enum::enum_name(firstFilter)};
        logEntry("We include stacks of cores in our search because the most important selectd filter is " + firstFilterString + ".");
        auto filteredMagnetics = apply_filters(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMagnetics.size() >= maximumNumberResults) {
            return filteredMagnetics;
        }

    }
    else {
        needToAddStacks = true;
        magnetics = create_magnetic_dataset(inputs, cores, false, onlyMaterialsForFilters);
        logEntry("We start the search with " + std::to_string(magnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.");
        logEntry("We don't include stacks of cores in our search.");
        auto filteredMagnetics = apply_filters(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMagnetics.size() >= maximumNumberResults) {
            return filteredMagnetics;
        }
    }

    auto globalIncludeStacks = settings->get_core_adviser_include_stacks();
    if (needToAddStacks && globalIncludeStacks) {
        expand_magnetic_dataset_with_stacks(inputs, cores, &magnetics);
    }

    logEntry("First attempt produced not enough results, so now we are searching again with " + std::to_string(magnetics.size()) + " magnetics, including up to " + std::to_string(defaults.coreAdviserMaximumNumberStacks) + " cores stacked when possible.");
    maximumMagneticsAfterFiltering = magnetics.size();
    auto filteredMagnetics = apply_filters(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    return filteredMagnetics;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::create_magnetic_dataset(Inputs inputs, std::vector<Core>* cores, bool includeStacks, bool onlyMaterialsForFilters) {
    std::vector<std::pair<Magnetic, double>> magnetics;
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings->get_use_toroidal_cores();
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
        if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
            continue;
        }

        if (onlyMaterialsForFilters && !core.can_be_used_for_filtering()) {
            continue;
        }

        core.process_data();

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
                if (i!=0) {
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stack" );
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
                MagneticManufacturerInfo magneticmanufacturerinfo;
                if (i!=0) {
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticmanufacturerinfo.set_reference(core.get_name().value() + " " + std::to_string(1 + i) + " stack" );
                }
                magnetic.set_manufacturer_info(magneticmanufacturerinfo);
                magnetic.set_core(core);
                (*magnetics).push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }

    }
}

void add_initial_turns(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
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

Mas CoreAdviser::post_process_core(Magnetic magnetic, Inputs inputs) {
    MagneticEnergy magneticEnergy;
    Mas mas;
    mas.set_magnetic(magnetic);
    double temperature = 0; 
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        temperature = std::max(temperature, inputs.get_operating_point(operatingPointIndex).get_conditions().get_ambient_temperature());
    }

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

std::vector<std::pair<Mas, double>> CoreAdviser::apply_filters(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex){
        auto excitation = Inputs::get_primary_excitation(inputs.get_mutable_operating_points()[operatingPointIndex]);
        if (!excitation.get_voltage()) {
            inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0].set_voltage(Inputs::calculate_induced_voltage(excitation, resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance())));
            Inputs::set_current_as_magnetizing_current(&inputs.get_mutable_operating_points()[operatingPointIndex]);
        }
        else if (!excitation.get_magnetizing_current()) {

            auto magnetizingCurrent = Inputs::calculate_magnetizing_current(excitation,
                                                                                   resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance()),
                                                                                   false);

            inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
        }

    }

    MagneticCoreFilterAreaProduct filterAreaProduct(inputs);
    MagneticCoreFilterEnergyStored filterEnergyStored(inputs, _models);
    MagneticCoreFilterCost filterCost(inputs);
    MagneticCoreFilterLosses filterLosses(inputs, _models);
    MagneticCoreFilterDimensions filterDimensions;
    MagneticCoreFilterMinimumImpedance filterMinimumImpedance;

    filterAreaProduct.set_scorings(&_scorings);
    filterAreaProduct.set_valid_scorings(&_validScorings);
    filterAreaProduct.set_filter_configuration(&_filterConfiguration);
    filterEnergyStored.set_scorings(&_scorings);
    filterEnergyStored.set_valid_scorings(&_validScorings);
    filterEnergyStored.set_filter_configuration(&_filterConfiguration);
    filterCost.set_scorings(&_scorings);
    filterCost.set_valid_scorings(&_validScorings);
    filterCost.set_filter_configuration(&_filterConfiguration);
    filterLosses.set_scorings(&_scorings);
    filterLosses.set_valid_scorings(&_validScorings);
    filterLosses.set_filter_configuration(&_filterConfiguration);
    filterDimensions.set_scorings(&_scorings);
    filterDimensions.set_valid_scorings(&_validScorings);
    filterDimensions.set_filter_configuration(&_filterConfiguration);
    filterMinimumImpedance.set_scorings(&_scorings);
    filterMinimumImpedance.set_valid_scorings(&_validScorings);
    filterMinimumImpedance.set_filter_configuration(&_filterConfiguration);

    CoreAdviserFilters firstFilter = CoreAdviserFilters::AREA_PRODUCT;

    if (weights[CoreAdviserFilters::EFFICIENCY] != 0 && weights[CoreAdviserFilters::MINIMUM_IMPEDANCE] != 0) {
        throw std::runtime_error("EFFICIENCY and MINIMUM_IMPEDANCE filters cannot be used together in the core Adviser");
    }

    double maxWeight = 0;

    for (auto& pair : weights) {
        auto filter = pair.first;
        auto weight = pair.second;
        if (weight > maxWeight) {
            maxWeight = weight;
            firstFilter = filter;
        }
    }

    std::vector<std::pair<Magnetic, double>> magneticsWithScoring;

    switch (firstFilter) {
        case CoreAdviserFilters::AREA_PRODUCT: 
            magneticsWithScoring = filterAreaProduct.filter_magnetics(magnetics, inputs, weights[CoreAdviserFilters::AREA_PRODUCT], true);
            break;
        case CoreAdviserFilters::ENERGY_STORED: 
            magneticsWithScoring = filterEnergyStored.filter_magnetics(magnetics, inputs, weights[CoreAdviserFilters::ENERGY_STORED], true);
            break;
        case CoreAdviserFilters::COST: 
            magneticsWithScoring = filterCost.filter_magnetics(magnetics, inputs, weights[CoreAdviserFilters::COST], true);
            break;
        case CoreAdviserFilters::EFFICIENCY: 
            add_initial_turns(magnetics, inputs);
            magneticsWithScoring = filterLosses.filter_magnetics(magnetics, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
            break;
        case CoreAdviserFilters::DIMENSIONS: 
            magneticsWithScoring = filterDimensions.filter_magnetics(magnetics, weights[CoreAdviserFilters::DIMENSIONS], true);
            break;
        case CoreAdviserFilters::MINIMUM_IMPEDANCE: 
            magneticsWithScoring = filterMinimumImpedance.filter_magnetics(magnetics, inputs, weights[CoreAdviserFilters::MINIMUM_IMPEDANCE], true);
            break;
    }

    std::string firstFilterString = std::string{magic_enum::enum_name(firstFilter)};
    logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " magnetics after the first filter, which was " + firstFilterString + ".");

    if (magneticsWithScoring.size() == 0) {
        return {};
    }

    if (magneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
        magneticsWithScoring = std::vector<std::pair<Magnetic, double>>(magneticsWithScoring.begin(), magneticsWithScoring.end() - (magneticsWithScoring.size() - maximumMagneticsAfterFiltering));
        logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " after culling by the score on the first filter.");
    }

    if (firstFilter != CoreAdviserFilters::EFFICIENCY) {
        logEntry("Adding initial number of turns to " + std::to_string(magneticsWithScoring.size()) + " magnetics.");

        add_initial_turns(&magneticsWithScoring, inputs);
        logEntry("Added initial number of turns to " + std::to_string(magneticsWithScoring.size()) + " magnetics.");
    }
    
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        if (filter != firstFilter) {
            std::string filterString = std::string{magic_enum::enum_name(filter)};
            logEntry("Filtering by " + filterString + ".");
            switch (filter) {
                case CoreAdviserFilters::AREA_PRODUCT: 
                    magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::AREA_PRODUCT]);
                    break;
                case CoreAdviserFilters::ENERGY_STORED: 
                    magneticsWithScoring = filterEnergyStored.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::ENERGY_STORED]);
                    break;
                case CoreAdviserFilters::COST: 
                    magneticsWithScoring = filterCost.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::COST]);
                    break;
                case CoreAdviserFilters::EFFICIENCY: 
                    magneticsWithScoring = filterLosses.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY]);
                    break;
                case CoreAdviserFilters::DIMENSIONS: 
                    magneticsWithScoring = filterDimensions.filter_magnetics(&magneticsWithScoring, weights[CoreAdviserFilters::DIMENSIONS]);
                    break;
                case CoreAdviserFilters::MINIMUM_IMPEDANCE: 
                    magneticsWithScoring = filterMinimumImpedance.filter_magnetics(&magneticsWithScoring, inputs, weights[CoreAdviserFilters::MINIMUM_IMPEDANCE]);
                    break;
            }    
            logEntry("There are " + std::to_string(magneticsWithScoring.size()) + " after filtering by " + filterString + ".");
        }
    });

    if (magneticsWithScoring.size() > maximumNumberResults) {
        if (_uniqueCoreShapes) {
            std::vector<std::pair<Magnetic, double>> magneticsWithScoringAndUniqueShapes;
            std::vector<std::string> useShapes;
            for (size_t magneticIndex = 0; magneticIndex < magneticsWithScoring.size(); ++magneticIndex){
                Magnetic magnetic = magneticsWithScoring[magneticIndex].first;
                auto core = magnetic.get_core();
                if (std::find(useShapes.begin(), useShapes.end(), core.get_shape_name()) != useShapes.end()) {
                    continue;
                }
                else {
                    magneticsWithScoringAndUniqueShapes.push_back(magneticsWithScoring[magneticIndex]);
                    useShapes.push_back(core.get_shape_name());
                }

                if (magneticsWithScoringAndUniqueShapes.size() == maximumNumberResults) {
                    break;
                }
            }
            magneticsWithScoring.clear();
            magneticsWithScoring = magneticsWithScoringAndUniqueShapes;
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

} // namespace OpenMagnetics
