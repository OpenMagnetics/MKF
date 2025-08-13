#include "advisers/MagneticFilter.h"
#include "advisers/WireAdviser.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Coil.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Settings.h"
#include <list>


namespace OpenMagnetics {

void normalize_scoring(std::vector<std::pair<CoilFunctionalDescription, double>>* coilsWithScoring, std::vector<double>* newScoring, bool invert=true) {
    auto normalizedScorings = OpenMagnetics::normalize_scoring(*newScoring, 1, invert, false);

    for (size_t i = 0; i < (*coilsWithScoring).size(); ++i) {
        (*coilsWithScoring)[i].second += normalizedScorings[i];
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring, false);
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring);
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::filter_by_skin_losses_density(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                     SignalDescriptor current,
                                                                                                     double temperature) {
    std::vector<std::pair<CoilFunctionalDescription, double>> filteredCoilsWithScoring;
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
        throw std::runtime_error("Something wrong happened while filtering, size of unfilteredCoils: " + std::to_string(filteredCoilsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if (filteredCoilsWithScoring.size() > 0) {
        normalize_scoring(&filteredCoilsWithScoring, &newScoring);
    }
    return filteredCoilsWithScoring;
}

std::vector<std::pair<OpenMagnetics::CoilFunctionalDescription, double>> WireAdviser::get_advised_wire(OpenMagnetics ::CoilFunctionalDescription coilFunctionalDescription,
                                                                                        Section section,
                                                                                        SignalDescriptor current,
                                                                                        double temperature,
                                                                                        uint8_t numberSections,
                                                                                        size_t maximumNumberResults){
    std::vector<Wire> wires;
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

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::create_planar_dataset(CoilFunctionalDescription coilFunctionalDescription,
                                                                                             Section section,
                                                                                             uint8_t numberSections) {
    std::vector<std::pair<CoilFunctionalDescription, double>> coilFunctionalDescriptions;
    auto planarWires = get_wires(WireType::PLANAR);

    // No paralells
    {
        auto maximumNumberTurnsPerSection = ceil(double(coilFunctionalDescription.get_number_turns()) / numberSections);
        auto maximumAvailableWidthForCopper = section.get_dimensions()[0] - 2 * get_border_to_wire_distance() - (maximumNumberTurnsPerSection - 1) * get_wire_to_wire_distance();
        if (maximumAvailableWidthForCopper < 0) {
            return coilFunctionalDescriptions;
            // throw std::runtime_error("maximumAvailableWidthForCopper cannot be negative");
        }
        auto maximumAvailableWidthForTurn = maximumAvailableWidthForCopper / maximumNumberTurnsPerSection;

        for (auto wire : planarWires) {
            if (resolve_dimensional_values(wire.get_conducting_height().value()) < section.get_dimensions()[1]) {
                wire.set_nominal_value_outer_height(resolve_dimensional_values(wire.get_conducting_height().value()));
                wire.set_nominal_value_conducting_width(maximumAvailableWidthForTurn);
                wire.set_nominal_value_outer_width(maximumAvailableWidthForTurn);
                wire.set_nominal_value_conducting_area(maximumAvailableWidthForTurn * resolve_dimensional_values(wire.get_conducting_height().value()));
                coilFunctionalDescription.set_wire(wire);
                coilFunctionalDescription.set_number_parallels(1);
                coilFunctionalDescriptions.push_back(std::pair<CoilFunctionalDescription, double>{coilFunctionalDescription, 0});
            }
        }

    }

    // Paralells
    {
        auto maximumNumberTurnsPerSection = coilFunctionalDescription.get_number_turns();
        auto maximumAvailableWidthForCopper = section.get_dimensions()[0] - 2 * get_border_to_wire_distance() - (maximumNumberTurnsPerSection - 1) * get_wire_to_wire_distance();
        if (maximumAvailableWidthForCopper < 0) {
            return coilFunctionalDescriptions;
            // throw std::runtime_error("maximumAvailableWidthForCopper cannot be negative");
        }
        auto maximumAvailableWidthForTurn = maximumAvailableWidthForCopper / coilFunctionalDescription.get_number_turns();
        size_t maximumNumberParallels = numberSections;

        for (auto wire : planarWires) {
            if (resolve_dimensional_values(wire.get_conducting_height().value()) < section.get_dimensions()[1]) {
                wire.set_nominal_value_outer_height(resolve_dimensional_values(wire.get_conducting_height().value()));
                for (size_t numberParallels = 2; numberParallels <= maximumNumberParallels; ++numberParallels) {
                    wire.set_nominal_value_conducting_width(maximumAvailableWidthForTurn);
                    wire.set_nominal_value_outer_width(maximumAvailableWidthForTurn);
                    wire.set_nominal_value_conducting_area(maximumAvailableWidthForTurn * resolve_dimensional_values(wire.get_conducting_height().value()));
                    coilFunctionalDescription.set_wire(wire);
                    coilFunctionalDescription.set_number_parallels(numberParallels);
                    coilFunctionalDescriptions.push_back(std::pair<CoilFunctionalDescription, double>{coilFunctionalDescription, 0});
                }
            }
        }


    }
    return coilFunctionalDescriptions;
}

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::create_dataset(CoilFunctionalDescription coilFunctionalDescription,
                                                                                      std::vector<Wire>* wires,
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
            numberParallelsNeeded = Wire::calculate_number_parallels_needed(current, temperature, wire, _maximumEffectiveCurrentDensity);
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
        auto wire = Coil::resolve_wire((*unfilteredCoils)[coilIndex].first);
        if (!Coil::resolve_wire((*unfilteredCoils)[coilIndex].first).get_conducting_area()) {
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

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::get_advised_planar_wire(CoilFunctionalDescription coilFunctionalDescription,
                                                                                               Section section,
                                                                                               SignalDescriptor current,
                                                                                               double temperature,
                                                                                               uint8_t numberSections,
                                                                                               size_t maximumNumberResults) {
    auto coilsWithScoring = create_planar_dataset(coilFunctionalDescription, section, numberSections);

    logEntry("We start the search with " + std::to_string(coilsWithScoring.size()) + " wires");

    coilsWithScoring = filter_by_effective_resistance(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " planar wires after filtering by effective resistance.");

    coilsWithScoring = filter_by_skin_losses_density(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " planar wires after filtering by skin losses density.");

    coilsWithScoring = filter_by_proximity_factor(&coilsWithScoring, current, temperature);
    logEntry("There are " + std::to_string(coilsWithScoring.size()) + " planar wires after filtering by proximity factor.");

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

std::vector<std::pair<CoilFunctionalDescription, double>> WireAdviser::get_advised_wire(std::vector<Wire>* wires,
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
