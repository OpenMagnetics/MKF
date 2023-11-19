#include "WireAdviser.h"
#include "WireWrapper.h"
#include "CoilWrapper.h"
#include "WindingLosses.h"
#include "WindingSkinEffectLosses.h"
#include <list>


namespace OpenMagnetics {

void normalize_scoring(std::vector<std::pair<CoilFunctionalDescription, double>>* coilsWithScoring, std::vector<double>* newScoring, bool invert=true) {
    double maximumScoring = *std::max_element(newScoring->begin(), newScoring->end());
    double minimumScoring = *std::min_element(newScoring->begin(), newScoring->end());

    for (size_t i = 0; i < (*coilsWithScoring).size(); ++i) {
        auto mas = (*coilsWithScoring)[i].first;
        auto scoring = (*newScoring)[i];
        if (maximumScoring != minimumScoring) {
            if (invert) {
                (*coilsWithScoring)[i].second += 1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring);
            }
            else {
                (*coilsWithScoring)[i].second += (scoring - minimumScoring) / (maximumScoring - minimumScoring);
            }
        }
        else {
            (*coilsWithScoring)[i].second += 1;
        }
    }
    sort((*coilsWithScoring).begin(), (*coilsWithScoring).end(), [](std::pair<CoilFunctionalDescription, double>& b1, std::pair<CoilFunctionalDescription, double>& b2) {
        return b1.second > b2.second;
    });
}

std::vector<std::pair<CoilFunctionalDescription, double>>  WireAdviser::filter_by_area_no_parallels(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                    Section section) {
    std::vector<std::pair<CoilFunctionalDescription, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){  
        auto wire = CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first);
        if (wire.get_maximum_outer_width() < section.get_dimensions()[0] && wire.get_maximum_outer_height() < section.get_dimensions()[1]) {
            double scoring = 0;
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring);
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<CoilFunctionalDescription, double>>  WireAdviser::filter_by_area_with_parallels(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                    Section section,
                                                                                                    double numberSections) {
    std::vector<std::pair<CoilFunctionalDescription, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){  
        auto wire = CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first);
        if (!CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first).get_conducting_area()) {
            throw std::runtime_error("Conducting area is missing");
        }
        auto neededOuterAreaNoCompact = wire.get_maximum_outer_width() * wire.get_maximum_outer_height();

        neededOuterAreaNoCompact *= (*unfilteredCoils)[coilIndex].first.get_number_parallels() * (*unfilteredCoils)[coilIndex].first.get_number_turns() / numberSections;

        if (neededOuterAreaNoCompact < (section.get_dimensions()[0] * section.get_dimensions()[1])) {
            // double scoring = (section.get_dimensions()[0] * section.get_dimensions()[1]) - neededOuterAreaNoCompact;
            double scoring = 1;
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, false);
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::filter_by_skin_depth_resistance(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                      SignalDescriptor current,
                                                                                                      double temperature) {
    std::vector<std::pair<CoilFunctionalDescription, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){  
        auto wire = CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first);

        double effective_resistance_per_meter = WindingLosses::calculate_skin_effect_resistance_per_meter(wire, current, temperature);

        double scoring = effective_resistance_per_meter;
        newScoring.push_back(scoring);
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring);
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::filter_by_effective_resistance(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                      SignalDescriptor current,
                                                                                                      double temperature) {
    std::vector<std::pair<CoilFunctionalDescription, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    if (!current.get_processed()->get_effective_frequency()) {
        throw std::runtime_error("Current processed is missing field effective frequency");
    }
    auto currentEffectiveFrequency = current.get_processed()->get_effective_frequency().value();

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){  
        auto wire = CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first);

        double effective_resistance_per_meter = WindingLosses::calculate_effective_resistance_per_meter(wire, currentEffectiveFrequency, temperature);

        double scoring = effective_resistance_per_meter;
        newScoring.push_back(scoring);
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring);
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::filter_by_proximity_factor(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                      SignalDescriptor current,
                                                                                                      double temperature) {
    std::vector<std::pair<CoilFunctionalDescription, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    if (!current.get_processed()->get_effective_frequency()) {
        throw std::runtime_error("Current processed is missing field effective frequency");
    }
    auto currentEffectiveFrequency = current.get_processed()->get_effective_frequency().value();

    double effectiveSkinDepth = WindingSkinEffectLosses::calculate_skin_depth("copper", currentEffectiveFrequency, temperature);

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){  
        auto wire = CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first);

        if (!wire.get_number_conductors()) {
            wire.set_number_conductors(1);
        }
        double proximityFactor = wire.get_minimum_conducting_dimension() / effectiveSkinDepth * pow(wire.get_number_conductors().value(), 2);
        proximityFactor = wire.get_minimum_conducting_dimension() / effectiveSkinDepth * pow((*unfilteredCoils)[coilIndex].first.get_number_parallels() / (std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height())), 2);

        double scoring = proximityFactor;
        newScoring.push_back(scoring);
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring);
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::get_advised_wire(CoilFunctionalDescription coilFunctionalDescription,
                                                                                        Section section,
                                                                                        SignalDescriptor current,
                                                                                        double temperature,
                                                                                        uint8_t numberSections,
                                                                                        size_t maximumNumberResults){
    std::string file_path = __FILE__;
    auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/../../MAS/data/wires.ndjson");
    std::ifstream ndjsonFile(inventory_path);
    std::string jsonLine;
    std::vector<WireWrapper> wires;
    while (std::getline(ndjsonFile, jsonLine)) {
        json jf = json::parse(jsonLine);
        WireWrapper wire(jf);
        if ((_includeFoil || wire.get_type() != WireType::FOIL) &&
            (_includeRectangular || wire.get_type() != WireType::RECTANGULAR) &&
            (_includeLitz || wire.get_type() != WireType::LITZ) &&
            (_includeRound || wire.get_type() != WireType::ROUND)) {
            wires.push_back(wire);
        }
    }
    return get_advised_wire(&wires, coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
}


std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::create_dataset(CoilFunctionalDescription coilFunctionalDescription,
                                                                                      std::vector<WireWrapper>* wires,
                                                                                      Section section,
                                                                                      SignalDescriptor current,
                                                                                      double temperature){
    std::vector<std::pair<CoilFunctionalDescription, double>> coilFunctionalDescriptions;

    for (auto& wire : *wires){
        if ((!_includeFoil && wire.get_type() == WireType::FOIL) ||
            (!_includeRectangular && wire.get_type() == WireType::RECTANGULAR) ||
            (!_includeLitz && wire.get_type() == WireType::LITZ) ||
            (!_includeRound && wire.get_type() == WireType::ROUND)) {
            continue;
        }
        int numberParallelsNeeded;
        if (wire.get_type() == WireType::FOIL) {
            wire.cut_foil_wire_to_section(section);
        }

        if (wire.get_type() == WireType::RECTANGULAR) {
            numberParallelsNeeded = 1;
        }
        else {
            numberParallelsNeeded = WireWrapper::calculate_number_parallels_needed(current, temperature, wire, _maximumEffectiveCurrentDensity);
            if (numberParallelsNeeded > _maximumNumberParallels) {
                continue;
            }
        }

        coilFunctionalDescription.set_number_parallels(numberParallelsNeeded);
        coilFunctionalDescription.set_wire(wire);
        coilFunctionalDescriptions.push_back(std::pair<CoilFunctionalDescription, double>{coilFunctionalDescription, 0});
        if (numberParallelsNeeded < _maximumNumberParallels) {
            coilFunctionalDescription.set_number_parallels(numberParallelsNeeded);
            coilFunctionalDescriptions.push_back(std::pair<CoilFunctionalDescription, double>{coilFunctionalDescription, 0});
        }
    }

    return coilFunctionalDescriptions;
}

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::get_advised_wire(std::vector<WireWrapper>* wires,
                                                                                        CoilFunctionalDescription coilFunctionalDescription,
                                                                                        Section section,
                                                                                        SignalDescriptor current,
                                                                                        double temperature,
                                                                                        uint8_t numberSections,
                                                                                        size_t maximumNumberResults){
    auto coilsWithScoring = create_dataset(coilFunctionalDescription, wires, section, current, temperature);
    coilsWithScoring = filter_by_area_no_parallels(&coilsWithScoring, section);
    coilsWithScoring = filter_by_area_with_parallels(&coilsWithScoring, section, numberSections);
    coilsWithScoring = filter_by_effective_resistance(&coilsWithScoring, current, temperature);
    coilsWithScoring = filter_by_proximity_factor(&coilsWithScoring, current, temperature);


    if (coilsWithScoring.size() > maximumNumberResults) {
        return std::vector<std::pair<CoilFunctionalDescription, double>>(coilsWithScoring.begin(), coilsWithScoring.end() - (coilsWithScoring.size() - maximumNumberResults));
    }

    return coilsWithScoring;
}
} // namespace OpenMagnetics
