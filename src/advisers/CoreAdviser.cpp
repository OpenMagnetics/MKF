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

void normalize_scoring(std::vector<std::pair<Mas, double>>* masMagneticsWithScoring, std::vector<double> newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    auto normalizedScorings = OpenMagnetics::normalize_scoring(newScoring, weight, filterConfiguration);

    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i) {
        (*masMagneticsWithScoring)[i].second += normalizedScorings[i];
    }

    sort((*masMagneticsWithScoring).begin(), (*masMagneticsWithScoring).end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
        return b1.second > b2.second;
    });
}

CoreAdviser::MagneticCoreFilterAreaProduct::MagneticCoreFilterAreaProduct(Inputs inputs) {
    _filter = MagneticFilterAreaProduct(inputs);
}

std::vector<std::pair<Mas, double>> CoreAdviser::MagneticCoreFilterAreaProduct::filter_magnetics(std::vector<std::pair<Mas, double>>* unfilteredMasMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    std::vector<std::pair<Mas, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::map<std::string, double> scaledMagneticFluxDensitiesPerMaterial;

    std::list<size_t> listOfIndexesToErase;

    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        Mas mas = (*unfilteredMasMagnetics)[masIndex].first;
        Magnetic magnetic = mas.get_magnetic();

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::AREA_PRODUCT)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT][magnetic.get_manufacturer_info().value().get_reference().value()]) {

                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }

        auto [valid, scoring] = _filter.evaluate_magnetic(&mas, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::AREA_PRODUCT, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }


    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT]);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterEnergyStored::MagneticCoreFilterEnergyStored(std::map<std::string, std::string> models) {
    _filter = MagneticFilterEnergyStored(models);
}

std::vector<std::pair<Mas, double>> CoreAdviser::MagneticCoreFilterEnergyStored::filter_magnetics(std::vector<std::pair<Mas, double>>* unfilteredMasMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    std::vector<std::pair<Mas, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){  
        Mas mas = (*unfilteredMasMagnetics)[masIndex].first;
        Magnetic magnetic = Magnetic(mas.get_magnetic());

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::ENERGY_STORED)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }
        auto [valid, scoring] = _filter.evaluate_magnetic(&mas, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
            (*unfilteredMasMagnetics)[masIndex].first = mas;
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::ENERGY_STORED, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }
    }


    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::ENERGY_STORED]);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterCost::MagneticCoreFilterCost(Inputs inputs) {
    _filter = MagneticFilterCost(inputs);
}

std::vector<std::pair<Mas, double>> CoreAdviser::MagneticCoreFilterCost::filter_magnetics(std::vector<std::pair<Mas, double>>* unfilteredMasMagnetics, Inputs inputs, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    std::vector<std::pair<Mas, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        Mas mas = (*unfilteredMasMagnetics)[masIndex].first;
        Magnetic magnetic = Magnetic(mas.get_magnetic());
        auto core = magnetic.get_core();

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::COST)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::COST].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::COST][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::COST][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }
        auto [valid, scoring] = _filter.evaluate_magnetic(&mas, &inputs);

        if (valid) {
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::COST, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }
    }


    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::COST]);
    }
    return filteredMagneticsWithScoring;
}

CoreAdviser::MagneticCoreFilterLosses::MagneticCoreFilterLosses(Inputs inputs, std::map<std::string, std::string> models) {
    _filter = MagneticFilterLosses(inputs, models);
}

std::vector<std::pair<Mas, double>> CoreAdviser::MagneticCoreFilterLosses::filter_magnetics(std::vector<std::pair<Mas, double>>* unfilteredMasMagnetics, Inputs inputs, double weight, bool firstFilter) {
    auto coilDelimitAndCompactOld = settings->get_coil_delimit_and_compact();
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    std::vector<std::pair<Mas, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        Mas mas = (*unfilteredMasMagnetics)[masIndex].first;
        Magnetic magnetic = Magnetic(mas.get_magnetic());
        auto core = magnetic.get_core();

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::EFFICIENCY)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::EFFICIENCY].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::EFFICIENCY][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::EFFICIENCY][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }
        auto [valid, scoring] = _filter.evaluate_magnetic(&mas, &inputs);

        if (valid) {
            (*unfilteredMasMagnetics)[masIndex].first = mas;
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::EFFICIENCY, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }
    }

    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }
    // (*unfilteredMasMagnetics).clear();

    if (filteredMagneticsWithScoring.size() == 0) {
        settings->set_coil_delimit_and_compact(coilDelimitAndCompactOld);
        return *unfilteredMasMagnetics;
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredMagneticsWithScoring.size() > 0) {
        normalize_scoring(&filteredMagneticsWithScoring, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::EFFICIENCY]);
    }
    settings->set_coil_delimit_and_compact(coilDelimitAndCompactOld);
    return filteredMagneticsWithScoring;
}

std::vector<std::pair<Mas, double>> CoreAdviser::MagneticCoreFilterDimensions::filter_magnetics(std::vector<std::pair<Mas, double>>* unfilteredMasMagnetics, double weight, bool firstFilter) {
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    std::vector<double> newScoring;

    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        Mas mas = (*unfilteredMasMagnetics)[masIndex].first;
        Magnetic magnetic = Magnetic(mas.get_magnetic());

        auto [valid, scoring] = _filter.evaluate_magnetic(&mas, nullptr);

        newScoring.push_back(scoring);
        add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::DIMENSIONS, scoring, firstFilter);
    }

    if ((*unfilteredMasMagnetics).size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string((*unfilteredMasMagnetics).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredMasMagnetics).size() > 0) {
        normalize_scoring(unfilteredMasMagnetics, newScoring, weight, (*_filterConfiguration)[CoreAdviser::CoreAdviserFilters::DIMENSIONS]);
    }
    return (*unfilteredMasMagnetics);
}

std::vector<std::pair<Mas, double>> CoreAdviser::MagneticCoreFilterMinimumImpedance::filter_magnetics(std::vector<std::pair<Mas, double>>* unfilteredMasMagnetics, Inputs inputs, double weight, bool firstFilter) {
    auto coilDelimitAndCompactOld = settings->get_coil_delimit_and_compact();
    if (weight <= 0) {
        return *unfilteredMasMagnetics;
    }
    std::vector<std::pair<Mas, double>> filteredMagneticsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t masIndex = 0; masIndex < (*unfilteredMasMagnetics).size(); ++masIndex){
        Mas mas = (*unfilteredMasMagnetics)[masIndex].first;
        Magnetic magnetic = Magnetic(mas.get_magnetic());
        auto core = magnetic.get_core();

        if ((*_validScorings).contains(CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE)) {
            if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE].contains(magnetic.get_manufacturer_info().value().get_reference().value())) {
                if ((*_validScorings)[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE][magnetic.get_manufacturer_info().value().get_reference().value()]) {
                    newScoring.push_back((*_scorings)[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE][magnetic.get_manufacturer_info().value().get_reference().value()]);
                }
                else {
                    listOfIndexesToErase.push_back(masIndex);
                }
                continue;
            }
        }

        auto [valid, scoring] = _filter.evaluate_magnetic(&mas, &inputs);

        if (valid) {
            (*unfilteredMasMagnetics)[masIndex].first = mas;
            newScoring.push_back(scoring);
            add_scoring(magnetic.get_manufacturer_info().value().get_reference().value(), CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE, scoring, firstFilter);
        }
        else {
            listOfIndexesToErase.push_back(masIndex);
        }

    }

    for (size_t i = 0; i < (*unfilteredMasMagnetics).size(); ++i) {
        if (listOfIndexesToErase.size() > 0 && i == listOfIndexesToErase.front()) {
            listOfIndexesToErase.pop_front();
        }
        else {
            filteredMagneticsWithScoring.push_back((*unfilteredMasMagnetics)[i]);
        }
    }

    if (filteredMagneticsWithScoring.size() != newScoring.size()) {
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredMasMagnetics: " + std::to_string(filteredMagneticsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
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
    std::vector<std::pair<Mas, double>> masMagnetics;
    bool needToAddStacks = false;
    bool onlyMaterialsForFilters = weights[CoreAdviserFilters::MINIMUM_IMPEDANCE] > 0;

    if (firstFilter == CoreAdviserFilters::EFFICIENCY) {
        masMagnetics = create_mas_dataset(inputs, cores, true, onlyMaterialsForFilters);
        logEntry("We start the search with " + std::to_string(masMagnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.");
        std::string firstFilterString = std::string{magic_enum::enum_name(firstFilter)};
        logEntry("We include stacks of cores in our search because the most important selectd filter is " + firstFilterString + ".");
        auto filteredMasMagnetics = apply_filters(&masMagnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMasMagnetics.size() >= maximumNumberResults) {
            return filteredMasMagnetics;
        }

    }
    else {
        needToAddStacks = true;
        masMagnetics = create_mas_dataset(inputs, cores, false, onlyMaterialsForFilters);
        logEntry("We start the search with " + std::to_string(masMagnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.");
        logEntry("We don't include stacks of cores in our search.");
        auto filteredMasMagnetics = apply_filters(&masMagnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMasMagnetics.size() >= maximumNumberResults) {
            return filteredMasMagnetics;
        }
    }

    auto globalIncludeStacks = settings->get_core_adviser_include_stacks();
    if (needToAddStacks && globalIncludeStacks) {
        expand_mas_dataset_with_stacks(inputs, cores, &masMagnetics);
    }

    logEntry("First attempt produced not enough results, so now we are searching again with " + std::to_string(masMagnetics.size()) + " magnetics, including up to " + std::to_string(defaults.coreAdviserMaximumNumberStacks) + " cores stacked when possible.");
    maximumMagneticsAfterFiltering = masMagnetics.size();
    auto filteredMasMagnetics = apply_filters(&masMagnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    return filteredMasMagnetics;
}

std::vector<std::pair<Mas, double>> CoreAdviser::create_mas_dataset(Inputs inputs, std::vector<Core>* cores, bool includeStacks, bool onlyMaterialsForFilters) {
    std::vector<std::pair<Mas, double>> masMagnetics;
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
    Mas mas;
    Outputs outputs;
    for (size_t i = 0; i < inputs.get_operating_points().size(); ++i)
    {
        mas.get_mutable_outputs().push_back(outputs);
    }
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
                mas.set_magnetic(magnetic);
                masMagnetics.push_back(std::pair<Mas, double>{mas, 0});
            }
        }
        else {
            magnetic.set_core(core);
            MagneticManufacturerInfo magneticmanufacturerinfo;
            magneticmanufacturerinfo.set_reference(core.get_name().value());
            magnetic.set_manufacturer_info(magneticmanufacturerinfo);
            mas.set_magnetic(magnetic);
            masMagnetics.push_back(std::pair<Mas, double>{mas, 0});
        }
    }

    return masMagnetics;
}

void CoreAdviser::expand_mas_dataset_with_stacks(Inputs inputs, std::vector<Core>* cores, std::vector<std::pair<Mas, double>>* masMagnetics) {
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings->get_use_toroidal_cores();
    double maximumHeight = DBL_MAX;
    if (inputs.get_design_requirements().get_maximum_dimensions()) {
        if (inputs.get_design_requirements().get_maximum_dimensions()->get_height()) {
            maximumHeight = inputs.get_design_requirements().get_maximum_dimensions()->get_height().value();
        }
    }


    Magnetic magnetic;
    Mas mas;
    Outputs outputs;
    for (size_t i = 0; i < inputs.get_operating_points().size(); ++i)
    {
        mas.get_mutable_outputs().push_back(outputs);
    }
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
                mas.set_magnetic(magnetic);
                (*masMagnetics).push_back(std::pair<Mas, double>{mas, 0});
            }
        }

    }
}

void add_initial_turns(std::vector<std::pair<Mas, double>> *masMagneticsWithScoring, Inputs inputs) {
    MagnetizingInductance magnetizingInductance;
    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i){

        Core core = (*masMagneticsWithScoring)[i].first.get_magnetic().get_core();
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        double initialNumberTurns = (*masMagneticsWithScoring)[i].first.get_magnetic().get_coil().get_functional_description()[0].get_number_turns();

        if (initialNumberTurns == 1) {
            initialNumberTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs, DimensionalValues::MINIMUM);
        }
        if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
            NumberTurns numberTurns(initialNumberTurns, inputs.get_design_requirements());
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            initialNumberTurns = numberTurnsCombination[0];
        }

        (*masMagneticsWithScoring)[i].first.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(initialNumberTurns);
    }
}

void correct_windings(std::vector<std::pair<Mas, double>> *masMagneticsWithScoring, Inputs inputs) {
    MagnetizingInductance magnetizingInductance;
    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i){

        Coil coil = Coil((*masMagneticsWithScoring)[i].first.get_magnetic().get_coil());
        NumberTurns numberTurns(coil.get_number_turns(0), inputs.get_design_requirements());
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();

        (*masMagneticsWithScoring)[i].first.set_inputs(inputs);
        (*masMagneticsWithScoring)[i].first.get_mutable_magnetic().get_mutable_coil().set_bobbin(Bobbin::create_quick_bobbin((*masMagneticsWithScoring)[i].first.get_magnetic().get_core()));
        for (size_t windingIndex = 1; windingIndex < numberTurnsCombination.size(); ++windingIndex) {
            auto winding = coil.get_functional_description()[0];
            winding.set_number_turns(numberTurnsCombination[windingIndex]);
            winding.set_isolation_side(get_isolation_side_from_index(windingIndex));
            winding.set_name(std::string{magic_enum::enum_name(get_isolation_side_from_index(windingIndex))});
            (*masMagneticsWithScoring)[i].first.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description().push_back(winding);
        }
    }
}

std::vector<std::pair<Mas, double>> CoreAdviser::apply_filters(std::vector<std::pair<Mas, double>>* masMagnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults){
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
    MagneticCoreFilterEnergyStored filterEnergyStored(_models);
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

    std::vector<std::pair<Mas, double>> masMagneticsWithScoring;

    switch (firstFilter) {
        case CoreAdviserFilters::AREA_PRODUCT: 
            masMagneticsWithScoring = filterAreaProduct.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::AREA_PRODUCT], true);
            break;
        case CoreAdviserFilters::ENERGY_STORED: 
            masMagneticsWithScoring = filterEnergyStored.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::ENERGY_STORED], true);
            break;
        case CoreAdviserFilters::COST: 
            masMagneticsWithScoring = filterCost.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::COST], true);
            break;
        case CoreAdviserFilters::EFFICIENCY: 
            add_initial_turns(masMagnetics, inputs);
            masMagneticsWithScoring = filterLosses.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::EFFICIENCY], true);
            break;
        case CoreAdviserFilters::DIMENSIONS: 
            masMagneticsWithScoring = filterDimensions.filter_magnetics(masMagnetics, weights[CoreAdviserFilters::DIMENSIONS], true);
            break;
        case CoreAdviserFilters::MINIMUM_IMPEDANCE: 
            masMagneticsWithScoring = filterMinimumImpedance.filter_magnetics(masMagnetics, inputs, weights[CoreAdviserFilters::MINIMUM_IMPEDANCE], true);
            break;
    }

    std::string firstFilterString = std::string{magic_enum::enum_name(firstFilter)};
    logEntry("There are " + std::to_string(masMagneticsWithScoring.size()) + " magnetics after the first filter, which was " + firstFilterString + ".");

    if (masMagneticsWithScoring.size() == 0) {
        return masMagneticsWithScoring;
    }

    if (masMagneticsWithScoring.size() > maximumMagneticsAfterFiltering) {
        masMagneticsWithScoring = std::vector<std::pair<Mas, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumMagneticsAfterFiltering));
        logEntry("There are " + std::to_string(masMagneticsWithScoring.size()) + " after culling by the score on the first filter.");
    }

    if (firstFilter != CoreAdviserFilters::EFFICIENCY) {
        logEntry("Adding initial number of turns to " + std::to_string(masMagneticsWithScoring.size()) + " magnetics.");

        add_initial_turns(&masMagneticsWithScoring, inputs);
        logEntry("Added initial number of turns to " + std::to_string(masMagneticsWithScoring.size()) + " magnetics.");
    }
    
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        if (filter != firstFilter) {
            std::string filterString = std::string{magic_enum::enum_name(filter)};
            logEntry("Filtering by " + filterString + ".");
            switch (filter) {
                case CoreAdviserFilters::AREA_PRODUCT: 
                    masMagneticsWithScoring = filterAreaProduct.filter_magnetics(&masMagneticsWithScoring, inputs, weights[CoreAdviserFilters::AREA_PRODUCT]);
                    break;
                case CoreAdviserFilters::ENERGY_STORED: 
                    masMagneticsWithScoring = filterEnergyStored.filter_magnetics(&masMagneticsWithScoring, inputs, weights[CoreAdviserFilters::ENERGY_STORED]);
                    break;
                case CoreAdviserFilters::COST: 
                    masMagneticsWithScoring = filterCost.filter_magnetics(&masMagneticsWithScoring, inputs, weights[CoreAdviserFilters::COST]);
                    break;
                case CoreAdviserFilters::EFFICIENCY: 
                    masMagneticsWithScoring = filterLosses.filter_magnetics(&masMagneticsWithScoring, inputs, weights[CoreAdviserFilters::EFFICIENCY]);
                    break;
                case CoreAdviserFilters::DIMENSIONS: 
                    masMagneticsWithScoring = filterDimensions.filter_magnetics(&masMagneticsWithScoring, weights[CoreAdviserFilters::DIMENSIONS]);
                    break;
                case CoreAdviserFilters::MINIMUM_IMPEDANCE: 
                    masMagneticsWithScoring = filterMinimumImpedance.filter_magnetics(&masMagneticsWithScoring, inputs, weights[CoreAdviserFilters::MINIMUM_IMPEDANCE]);
                    break;
            }    
            logEntry("There are " + std::to_string(masMagneticsWithScoring.size()) + " after filtering by " + filterString + ".");
        }
    });

    if (masMagneticsWithScoring.size() > maximumNumberResults) {
        if (_uniqueCoreShapes) {
            std::vector<std::pair<Mas, double>> masMagneticsWithScoringAndUniqueShapes;
            std::vector<std::string> useShapes;
            for (size_t masIndex = 0; masIndex < masMagneticsWithScoring.size(); ++masIndex){
                Mas mas = masMagneticsWithScoring[masIndex].first;
                auto core = mas.get_magnetic().get_core();
                if (std::find(useShapes.begin(), useShapes.end(), core.get_shape_name()) != useShapes.end()) {
                    continue;
                }
                else {
                    masMagneticsWithScoringAndUniqueShapes.push_back(masMagneticsWithScoring[masIndex]);
                    useShapes.push_back(core.get_shape_name());
                }

                if (masMagneticsWithScoringAndUniqueShapes.size() == maximumNumberResults) {
                    break;
                }
            }
            masMagneticsWithScoring.clear();
            masMagneticsWithScoring = masMagneticsWithScoringAndUniqueShapes;
        }
        else {
            masMagneticsWithScoring = std::vector<std::pair<Mas, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
        }
    }

    correct_windings(&masMagneticsWithScoring, inputs);

    return masMagneticsWithScoring;
}

} // namespace OpenMagnetics
