#include "WireAdviser.h"
#include "WireWrapper.h"
#include "CoilWrapper.h"
#include "WindingLosses.h"
#include "WindingSkinEffectLosses.h"
#include "Settings.h"
#include <list>


namespace OpenMagnetics {

void WireAdviser::logEntry(std::string entry) {
    // std::cout << entry << std::endl;
    _log += entry + "\n";
}

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

        if (wire.get_type() == WireType::FOIL && (*unfilteredCoils)[coilIndex].first.get_number_parallels() * (*unfilteredCoils)[coilIndex].first.get_number_turns() > _maximumNumberParallels) {
            listOfIndexesToErase.push_back(coilIndex);
            continue;
        }

        if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
            if (wire.get_maximum_outer_width() < section.get_dimensions()[0] && wire.get_maximum_outer_height() < section.get_dimensions()[1]) {
                double scoring = 0;
                newScoring.push_back(scoring);
            }
            else {
                listOfIndexesToErase.push_back(coilIndex);
            }
        }
        else {
            double wireAngle = wound_distance_to_angle(wire.get_maximum_outer_height(), wire.get_maximum_outer_width());

            if (wire.get_maximum_outer_width() < section.get_dimensions()[0] && wireAngle < section.get_dimensions()[1]) {
                double scoring = 0;
                newScoring.push_back(scoring);
            }
            else {
                listOfIndexesToErase.push_back(coilIndex);
            }

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
                                                                                                    double numberSections,
                                                                                                    bool allowNotFit) {
    std::vector<std::pair<CoilFunctionalDescription, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    double sectionArea;
    if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        sectionArea = section.get_dimensions()[0] * section.get_dimensions()[1];
    }
    else {
        sectionArea = std::numbers::pi * pow(section.get_dimensions()[0], 2) * section.get_dimensions()[1] / 360;
    }


    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){  
        auto wire = CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first);
        if (!CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first).get_conducting_area()) {
            throw std::runtime_error("Conducting area is missing");
        }
        auto neededOuterAreaNoCompact = wire.get_maximum_outer_width() * wire.get_maximum_outer_height();

        neededOuterAreaNoCompact *= (*unfilteredCoils)[coilIndex].first.get_number_parallels() * (*unfilteredCoils)[coilIndex].first.get_number_turns() / numberSections;

        if (neededOuterAreaNoCompact < sectionArea) {
            // double scoring = (section.get_dimensions()[0] * section.get_dimensions()[1]) - neededOuterAreaNoCompact;
            double scoring = 1;
            newScoring.push_back(scoring);
        }
        else if (allowNotFit) {
            double extra = (neededOuterAreaNoCompact - sectionArea) / sectionArea;
            if (extra < 0.5) {
                double scoring = extra;
                newScoring.push_back(scoring);
            }
            else {
                listOfIndexesToErase.push_back(coilIndex);
            }
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
        double proximityFactor = wire.get_minimum_conducting_dimension() / effectiveSkinDepth * pow(wire.get_number_conductors().value() * (*unfilteredCoils)[coilIndex].first.get_number_parallels() * (*unfilteredCoils)[coilIndex].first.get_number_turns(), 2);
        // proximityFactor = wire.get_minimum_conducting_dimension() / effectiveSkinDepth * pow((*unfilteredCoils)[coilIndex].first.get_number_parallels() / (std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height())), 2);

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

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::filter_by_solid_insulation_requirements(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils, WireSolidInsulationRequirements wireSolidInsulationRequirements) {
    std::vector<std::pair<CoilFunctionalDescription, double>> filteredCoilsWithScoring;
    std::vector<double> newScoring;

    std::list<size_t> listOfIndexesToErase;
    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){  
        auto wire = CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first);
        bool isValid = true;

        if (wire.get_type() == WireType::FOIL || wire.get_type() == WireType::PLANAR) {
            newScoring.push_back(0);
            continue;
        }

        if (!wire.resolve_coating()) {
            listOfIndexesToErase.push_back(coilIndex);
            continue;
        }

        auto coating = wire.resolve_coating().value();

        if (wire.get_type() == WireType::LITZ) {
            auto strand = wire.resolve_strand();
            coating = WireWrapper::resolve_coating(strand).value();
        }

        if (!coating.get_breakdown_voltage()) {
            throw std::runtime_error("Wire " + wire.get_name().value() + " is missing breakdown voltage");
        }


        if (coating.get_breakdown_voltage().value() < wireSolidInsulationRequirements.get_minimum_breakdown_voltage()) {
            isValid = false;
        }

        if (wireSolidInsulationRequirements.get_minimum_grade() && coating.get_grade()) {
            if (coating.get_grade().value() < wireSolidInsulationRequirements.get_minimum_grade().value()) {
                isValid = false;
            }
        }
        else if (wireSolidInsulationRequirements.get_minimum_number_layers() && coating.get_number_layers()) {
            if (coating.get_number_layers().value() < wireSolidInsulationRequirements.get_minimum_number_layers().value()) {
                isValid = false;
            }
        }
        else if (wireSolidInsulationRequirements.get_minimum_number_layers() || wireSolidInsulationRequirements.get_minimum_grade()) {
            isValid = false;
        }
        
        if (wireSolidInsulationRequirements.get_maximum_grade() && coating.get_grade()) {
            if (coating.get_grade().value() > wireSolidInsulationRequirements.get_maximum_grade().value()) {
                isValid = false;
            }
        }
        else if (wireSolidInsulationRequirements.get_maximum_number_layers() && coating.get_number_layers()) {
            if (coating.get_number_layers().value() > wireSolidInsulationRequirements.get_maximum_number_layers().value()) {
                isValid = false;
            }
        }
        else if (wireSolidInsulationRequirements.get_maximum_number_layers() || wireSolidInsulationRequirements.get_maximum_grade()) {
            isValid = false;
        }

        if (isValid) {
            double scoring = 0;
            if (wireSolidInsulationRequirements.get_minimum_breakdown_voltage() > 0) {
                scoring = wireSolidInsulationRequirements.get_minimum_breakdown_voltage() - coating.get_breakdown_voltage().value();
            }
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

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::get_advised_wire(CoilFunctionalDescription coilFunctionalDescription,
                                                                                        Section section,
                                                                                        SignalDescriptor current,
                                                                                        double temperature,
                                                                                        uint8_t numberSections,
                                                                                        size_t maximumNumberResults){
    std::vector<WireWrapper> wires;
    auto settings = OpenMagnetics::Settings::GetInstance();

    if (wireDatabase.empty()) {
        load_wires();
    }

    for (auto [name, wire] : wireDatabase) {
        if ((settings->get_wire_adviser_include_foil() || wire.get_type() != WireType::FOIL) &&
            (settings->get_wire_adviser_include_planar() ||  wire.get_type() != WireType::PLANAR) &&
            ((settings->get_wire_adviser_include_rectangular() && (settings->get_wire_adviser_allow_rectangular_in_toroidal_cores() || section.get_coordinate_system() == CoordinateSystem::CARTESIAN)) || wire.get_type() != WireType::RECTANGULAR) &&
            (settings->get_wire_adviser_include_litz() || wire.get_type() != WireType::LITZ) &&
            (settings->get_wire_adviser_include_round() || wire.get_type() != WireType::ROUND)) {

            if (!_commonWireStandard || !wire.get_standard()) {
                wires.push_back(wire);
            }
            else if (wire.get_standard().value() == _commonWireStandard){
                wires.push_back(wire);
            }
        }
    }
    return get_advised_wire(&wires, coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
}


std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::create_dataset(CoilFunctionalDescription coilFunctionalDescription,
                                                                                      std::vector<WireWrapper>* wires,
                                                                                      Section section,
                                                                                      SignalDescriptor current,
                                                                                      double temperature){
    auto settings = Settings::GetInstance();
    std::vector<std::pair<CoilFunctionalDescription, double>> coilFunctionalDescriptions;

    for (auto& wire : *wires){
        if (wire.get_type() == WireType::LITZ) {
            wire.set_strand(wire.resolve_strand());
        }
    }

    for (auto& wire : *wires){
        if ((!settings->get_wire_adviser_include_foil() && wire.get_type() == WireType::FOIL) ||
            (!settings->get_wire_adviser_include_planar() &&  wire.get_type() == WireType::PLANAR) ||
            (!(settings->get_wire_adviser_include_rectangular() && (settings->get_wire_adviser_allow_rectangular_in_toroidal_cores() || section.get_coordinate_system() == CoordinateSystem::CARTESIAN)) && wire.get_type() == WireType::RECTANGULAR) ||
            (!settings->get_wire_adviser_include_litz() && wire.get_type() == WireType::LITZ) ||
            (!settings->get_wire_adviser_include_round() && wire.get_type() == WireType::ROUND)) {
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
            numberParallelsNeeded = WireWrapper::calculate_number_parallels_needed(current, temperature, wire, _maximumEffectiveCurrentDensity);
            if (numberParallelsNeeded > _maximumNumberParallels) {
                continue;
            }
        }

        coilFunctionalDescription.set_number_parallels(numberParallelsNeeded);
        coilFunctionalDescription.set_wire(wire);
        coilFunctionalDescriptions.push_back(std::pair<CoilFunctionalDescription, double>{coilFunctionalDescription, 0});
        if (numberParallelsNeeded < _maximumNumberParallels) {
            coilFunctionalDescription.set_number_parallels(numberParallelsNeeded + 1);
            coilFunctionalDescriptions.push_back(std::pair<CoilFunctionalDescription, double>{coilFunctionalDescription, 0});
        }
    }

    return coilFunctionalDescriptions;
}

void WireAdviser::set_maximum_area_proportion(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils, Section section, uint8_t numberSections) {
    double sectionArea;

    if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        sectionArea = section.get_dimensions()[0] * section.get_dimensions()[1];
    }
    else {
        sectionArea = std::numbers::pi * pow(section.get_dimensions()[0], 2) * section.get_dimensions()[1] / 360;
    }

    for (size_t coilIndex = 0; coilIndex < (*unfilteredCoils).size(); ++coilIndex){  
        auto wire = CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first);
        if (!CoilWrapper::resolve_wire((*unfilteredCoils)[coilIndex].first).get_conducting_area()) {
            throw std::runtime_error("Conducting area is missing");
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
std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::get_advised_wire(std::vector<WireWrapper>* wires,
                                                                                        CoilFunctionalDescription coilFunctionalDescription,
                                                                                        Section section,
                                                                                        SignalDescriptor current,
                                                                                        double temperature,
                                                                                        uint8_t numberSections,
                                                                                        size_t maximumNumberResults){
    auto coilsWithScoring = create_dataset(coilFunctionalDescription, wires, section, current, temperature);

    logEntry("We start the search with " + std::to_string(coilsWithScoring.size()) + " wires");
    coilsWithScoring = filter_by_area_no_parallels(&coilsWithScoring, section);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by area no parallels.");

    if (_wireSolidInsulationRequirements) {
        coilsWithScoring = filter_by_solid_insulation_requirements(&coilsWithScoring, _wireSolidInsulationRequirements.value());
        logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by solid insulation.");
    }

    auto tempCoilsWithScoring = filter_by_area_with_parallels(&coilsWithScoring, section, numberSections, false);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by area with parallels.");

    if (tempCoilsWithScoring.size() == 0) {
        coilsWithScoring = filter_by_area_with_parallels(&coilsWithScoring, section, numberSections, true);
        logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by area with parallels, allowing not fitting.");
    }
    else{
        coilsWithScoring = tempCoilsWithScoring;
    }

    coilsWithScoring = filter_by_effective_resistance(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by effective resistance.");

    coilsWithScoring = filter_by_proximity_factor(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " after filtering by proximity factor.");

    if (coilsWithScoring.size() > maximumNumberResults) {
        auto finalCoilsWithScoring = std::vector<std::pair<CoilFunctionalDescription, double>>(coilsWithScoring.begin(), coilsWithScoring.end() - (coilsWithScoring.size() - maximumNumberResults));
        set_maximum_area_proportion(&finalCoilsWithScoring, section, numberSections);
        return finalCoilsWithScoring;
    }
    else {
        set_maximum_area_proportion(&coilsWithScoring, section, numberSections);
    }


    return coilsWithScoring;
}
} // namespace OpenMagnetics
