#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <cfloat>
#include <cmath>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Utils.h"
#include "constructive_models/Coil.h"
#include "json.hpp"
#include "constructive_models/InsulationMaterial.h"

#include <magic_enum.hpp>

using json = nlohmann::json;

namespace OpenMagnetics {


std::vector<double> Coil::cartesian_to_polar(std::vector<double> value) {
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        throw std::runtime_error("Not supposed to convert for these cores");
    }
    else {
        std::vector<double> convertedValue;
        double angle = atan2(value[1], value[0]) * 180 / std::numbers::pi;
        if (angle < 0) {
            angle += 360;
        }
        double radius = hypot(value[0], value[1]);
        double radialHeight = windingWindows[0].get_radial_height().value() - radius;
        return {radialHeight, angle};
    }
}

std::vector<double> Coil::cartesian_to_polar(std::vector<double> value, double radialHeight) {
    std::vector<double> convertedValue;
    double angle = atan2(value[1], value[0]) * 180 / std::numbers::pi;
    if (angle < 0) {
        angle += 360;
    }
    double radius = hypot(value[0], value[1]);
    double turnRadialHeight = radialHeight - radius;
    return {turnRadialHeight, angle};
}

std::vector<double> Coil::polar_to_cartesian(std::vector<double> value) {
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        throw std::runtime_error("Not supposed to convert for these cores");
    }
    else {
        std::vector<double> convertedValue;
        double radius = windingWindows[0].get_radial_height().value() - value[0];
        double angleRadians = value[1] / 180 * std::numbers::pi;
        double x = radius * cos(angleRadians);
        double y = radius * sin(angleRadians);
        return {x, y};
    }
}

std::vector<double> Coil::polar_to_cartesian(std::vector<double> value, double radialHeight) {
    std::vector<double> convertedValue;
    double radius = radialHeight - value[0];
    double angleRadians = value[1] / 180 * std::numbers::pi;
    double x = radius * cos(angleRadians);
    double y = radius * sin(angleRadians);
    return {x, y};
}

void Coil::convert_turns_to_cartesian_coordinates() {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    if (!get_turns_description()) {
        throw std::runtime_error("Missing turns");
    }

    auto turns = get_turns_description().value();
    if (turns[0].get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        return;
    }

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto cartesianCoordinates = polar_to_cartesian(turns[turnIndex].get_coordinates(), windingWindowRadialHeight);
        turns[turnIndex].set_coordinates(cartesianCoordinates);
        turns[turnIndex].set_coordinate_system(CoordinateSystem::CARTESIAN);
        if (turns[turnIndex].get_additional_coordinates()) {
            auto additionalCoordinates = turns[turnIndex].get_additional_coordinates().value();
            for (size_t additionalTurnIndex = 0; additionalTurnIndex < additionalCoordinates.size(); ++additionalTurnIndex) {
                additionalCoordinates[additionalTurnIndex] = polar_to_cartesian(additionalCoordinates[additionalTurnIndex], windingWindowRadialHeight);
            }
            turns[turnIndex].set_additional_coordinates(additionalCoordinates);
        }
    }

    set_turns_description(turns);
}

void Coil::convert_turns_to_polar_coordinates() {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    if (!get_turns_description()) {
        throw std::runtime_error("Missing turns");
    }

    auto turns = get_turns_description().value();
    if (turns[0].get_coordinate_system().value() == CoordinateSystem::POLAR) {
        return;
    }

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto cartesianCoordinates = cartesian_to_polar(turns[turnIndex].get_coordinates(), windingWindowRadialHeight);
        turns[turnIndex].set_coordinates(cartesianCoordinates);
        turns[turnIndex].set_coordinate_system(CoordinateSystem::POLAR);
        if (turns[turnIndex].get_additional_coordinates()) {
            auto additionalCoordinates = turns[turnIndex].get_additional_coordinates().value();
            for (size_t additionalTurnIndex = 0; additionalTurnIndex < additionalCoordinates.size(); ++additionalTurnIndex) {
                additionalCoordinates[additionalTurnIndex] = cartesian_to_polar(additionalCoordinates[additionalTurnIndex], windingWindowRadialHeight);
            }
            turns[turnIndex].set_additional_coordinates(additionalCoordinates);
        }
    }

    set_turns_description(turns);
}

Coil::Coil(const json& j, size_t interleavingLevel,
                               WindingOrientation windingOrientation,
                               WindingOrientation layersOrientation,
                               CoilAlignment turnsAlignment,
                               CoilAlignment sectionAlignment) {
    _interleavingLevel = interleavingLevel;
    _windingOrientation = windingOrientation;
    _layersOrientation = layersOrientation;
    _turnsAlignment = turnsAlignment;
    _sectionAlignment = sectionAlignment;
    from_json(j, *this);

    wind();
}

Coil::Coil(const json& j, bool windInConstructor) {
    from_json(j, *this);

    if (windInConstructor) {
        wind();
    }
}

Coil::Coil(const MAS::Coil coil) {
    bool hasSectionsData = false;
    bool hasLayersData = false;
    bool hasTurnsData = false;

    set_functional_description({});
    for (auto winding : coil.get_functional_description()) {
        get_mutable_functional_description().push_back(winding);
    }

    auto bobbinVariant = coil.get_bobbin();
    if (std::holds_alternative<std::string>(bobbinVariant)) {
        auto bobbinName = std::get<std::string>(bobbinVariant);
        set_bobbin(bobbinName);
    }
    else {
        auto bobbin = OpenMagnetics::Bobbin(bobbinVariant);
        set_bobbin(bobbin);
    }

    if (coil.get_sections_description()) {
        hasSectionsData = true;
        set_sections_description(coil.get_sections_description());
    }
    if (coil.get_layers_description()) {
        hasLayersData = true;
        set_layers_description(coil.get_layers_description());
    }
    if (coil.get_turns_description()) {
        hasTurnsData = true;
        set_turns_description(coil.get_turns_description());
    }
    auto delimitAndCompact = settings->get_coil_delimit_and_compact();

    if (!hasSectionsData || !hasLayersData || (!hasTurnsData && are_sections_and_layers_fitting())) {
        if (wind() && delimitAndCompact) {
            delimit_and_compact();
        }
    }

}

void Coil::log(std::string entry) {
    coilLog += entry + "\n";
}

std::string Coil::read_log() {
    return coilLog;
}

void Coil::set_strict(bool value) {
    _strict = value;
}

void Coil::set_inputs(Inputs inputs) {
    _inputs = inputs;
}

void Coil::set_interleaving_level(uint8_t interleavingLevel) {
    _interleavingLevel = interleavingLevel;
    _marginsPerSection = std::vector<std::vector<double>>(interleavingLevel, {0, 0});
}

void Coil::reset_margins_per_section() {
    _marginsPerSection.clear();
}

size_t Coil::get_interleaving_level() {
    return _currentRepetitions;
}

void Coil::set_winding_orientation(WindingOrientation windingOrientation) {
    _windingOrientation = windingOrientation;
    auto bobbin = resolve_bobbin();
    if (bobbin.get_processed_description()) {
        bobbin.set_winding_orientation(windingOrientation);
        set_bobbin(bobbin);
    }
}

void Coil::set_layers_orientation(WindingOrientation layersOrientation, std::optional<std::string> sectionName) {
    if (sectionName) {
        _layersOrientationPerSection[sectionName.value()] = layersOrientation;
    }
    else {
        _layersOrientation = layersOrientation;
    }
}

void Coil::set_turns_alignment(CoilAlignment turnsAlignment, std::optional<std::string> sectionName) {
    if (sectionName) {
        _turnsAlignmentPerSection[sectionName.value()] = turnsAlignment;
    }
    else {
        _turnsAlignment = turnsAlignment;
    }
}

void Coil::set_section_alignment(CoilAlignment sectionAlignment) {
    _sectionAlignment = sectionAlignment;
}

WindingOrientation Coil::get_winding_orientation() {
    auto bobbin = resolve_bobbin();
    auto windingOrientationFromBobbin = bobbin.get_winding_orientation();

    if (!windingOrientationFromBobbin) {
        return _windingOrientation;
    }
    else {
        return windingOrientationFromBobbin.value();
    }
}

WindingOrientation Coil::get_layers_orientation() {
    return _layersOrientation;
}

CoilAlignment Coil::get_turns_alignment(std::optional<std::string> sectionName) {
    if (sectionName) {
        if (_turnsAlignmentPerSection.find(sectionName.value()) != _turnsAlignmentPerSection.end()) {
            return _turnsAlignmentPerSection[sectionName.value()];
        }
        else {
            return _turnsAlignment;
        }
    }
    else {
        return _turnsAlignment;
    }
}

CoilAlignment Coil::get_section_alignment() {
    auto bobbin = resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        return _sectionAlignment;
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windingWindows.size() > 1) {
        throw std::runtime_error("Bobbins with more than winding window not implemented yet");
    }
    if (windingWindows[0].get_sections_alignment()) {
        return windingWindows[0].get_sections_alignment().value();
    }
    return _sectionAlignment;
}

bool Coil::fast_wind() {
    _strict = false;

    wind_by_sections();
    if (!get_sections_description()) {
        return false;
    }
    wind_by_layers();
    if (!get_layers_description()) {
        return false;
    }
    auto previousIncludeAdditionalCoordinates = settings->get_coil_include_additional_coordinates();
    settings->set_coil_include_additional_coordinates(false);
    wind_by_turns();
    settings->set_coil_include_additional_coordinates(previousIncludeAdditionalCoordinates);

    if (!get_turns_description()) {
        return false;
    }
    return true;
}

bool Coil::unwind() {
    set_sections_description(std::nullopt);
    set_layers_description(std::nullopt);
    set_turns_description(std::nullopt);
    return true;
}

bool Coil::wind() {
    std::vector<double> proportionPerWinding;

    proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    return wind(proportionPerWinding, pattern, _interleavingLevel);
}

bool Coil::wind(size_t repetitions){
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind(proportionPerWinding, pattern, repetitions);
}

bool Coil::wind(std::vector<size_t> pattern, size_t repetitions){
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind(proportionPerWinding, pattern, repetitions);
}

bool Coil::wind(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    bool windEvenIfNotFit = settings->get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings->get_coil_delimit_and_compact();
    bool tryRewind = settings->get_coil_try_rewind();
    std::string bobbinName = "";
    if (std::holds_alternative<std::string>(get_bobbin())) {
        bobbinName = std::get<std::string>(get_bobbin());
        if (bobbinName != "Dummy") {
            auto bobbinData = find_bobbin_by_name(std::get<std::string>(get_bobbin()));
            set_bobbin(bobbinData);
        }
    }
    _currentProportionPerWinding = proportionPerWinding;
    _currentPattern = pattern;
    _currentRepetitions = repetitions;

    if (bobbinName != "Dummy") {
        bool wind = true;                
        for (auto& winding : get_mutable_functional_description()) {
            if (std::holds_alternative<std::string>(winding.get_wire())) {
                std::string wireName = std::get<std::string>(winding.get_wire());
                if (wireName == "Dummy") {
                    wind = false;
                    break;
                }
                auto wire = find_wire_by_name(wireName);
                winding.set_wire(wire);
            }
        }

        if (wind) {
            set_sections_description(std::nullopt);
            set_layers_description(std::nullopt);
            set_turns_description(std::nullopt);

            if (_inputs) {
                if (_inputs->get_design_requirements().get_insulation()) {
                    logEntry("Calculating Required Insulation", "Coil", 2);
                    calculate_insulation();
                }
                else {
                    logEntry("Calculating Mechanical Insulation", "Coil", 2);
                    calculate_mechanical_insulation();
                }
            }
            else {
                logEntry("Calculating Mechanical Insulation", "Coil", 2);
                calculate_mechanical_insulation();
            }
            logEntry("Winding by sections", "Coil", 2);
            wind_by_sections(proportionPerWinding, pattern, repetitions);
            logEntry("Winding by layers", "Coil", 2);
            wind_by_layers();

            if (!get_layers_description()) {
                return false;
            }

            auto sections = get_sections_description().value();

            if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
                logEntry("Winding by turns", "Coil", 2);
                wind_by_turns();
                if (delimitAndCompact) {
                    logEntry("Delimiting and compacting", "Coil", 2);
                    delimit_and_compact();
                }
            }
            if (tryRewind && (!are_sections_and_layers_fitting() || !get_turns_description())) {
                logEntry("Trying to rewind", "Coil", 2);
                try_rewind();
            }
        }
    }
    return are_sections_and_layers_fitting() && get_turns_description();
}

bool Coil::wind_planar(std::vector<size_t> stackUp, std::optional<double> borderToWireDistance, std::optional<double> wireToWireDistance, std::optional<double> insulationThickness, double coreToLayerDistance) {
    bool windEvenIfNotFit = settings->get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings->get_coil_delimit_and_compact();
    std::string bobbinName = "";
    if (std::holds_alternative<std::string>(get_bobbin())) {
        bobbinName = std::get<std::string>(get_bobbin());
        if (bobbinName != "Dummy") {
            auto bobbinData = find_bobbin_by_name(std::get<std::string>(get_bobbin()));
            set_bobbin(bobbinData);
        }
    }

    if (bobbinName != "Dummy") {
        bool wind = true;                
        for (auto& winding : get_mutable_functional_description()) {
            if (std::holds_alternative<std::string>(winding.get_wire())) {
                std::string wireName = std::get<std::string>(winding.get_wire());
                if (wireName == "Dummy") {
                    wind = false;
                    break;
                }
                auto wire = find_wire_by_name(wireName);
                winding.set_wire(wire);
            }
        }

        if (wind) {
            set_groups_description(std::nullopt);
            set_sections_description(std::nullopt);
            set_layers_description(std::nullopt);
            set_turns_description(std::nullopt);

            if (_inputs) {
                if (_inputs->get_design_requirements().get_insulation()) {
                    logEntry("Calculating Required Insulation", "Coil", 2);
                    auto standardCoordinator = InsulationCoordinator();
                    auto clearance = standardCoordinator.calculate_clearance(_inputs.value());
                    if (!borderToWireDistance) {
                        borderToWireDistance = std::max(defaults.minimumBorderToWireDistance, clearance);
                    }
                    if (!wireToWireDistance) {
                        wireToWireDistance = std::max(defaults.minimumBorderToWireDistance, clearance);
                    }
                }
            }

            if (!borderToWireDistance) {
                borderToWireDistance = defaults.minimumBorderToWireDistance;
            }
            if (!wireToWireDistance) {
                wireToWireDistance = defaults.minimumWireToWireDistance;
            }

            logEntry("Winding by sections", "Coil", 2);
            auto result = wind_by_planar_sections(stackUp, insulationThickness, coreToLayerDistance);
            logEntry("Winding by layers", "Coil", 2);
            result = wind_by_planar_layers();

            if (!get_layers_description()) {
                return false;
            }

            if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
                logEntry("Winding by turns", "Coil", 2);
                result = wind_by_planar_turns(borderToWireDistance.value(), wireToWireDistance.value());
                if (delimitAndCompact) {
                    logEntry("Delimiting and compacting", "Coil", 2);
                    delimit_and_compact();
                }
            }
        }
    }
    return are_sections_and_layers_fitting() && get_turns_description();
}

std::vector<WindingStyle> Coil::wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, std::vector<size_t> numberSlots) {
    std::vector<WindingStyle> windByConsecutiveTurns;
    for (size_t i = 0; i < numberTurns.size(); ++i) {
        if (numberSlots[i] <= 0) {
            throw std::runtime_error("Number of slots cannot be less than 1, please verify your isolation sides requirement");
        }
        if (numberTurns[i] == numberSlots[i]) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
            log("Winding winding " + std::to_string(i) + " by putting together parallels of the same turn, as the number of turns is equal to the number of sections.");
            continue;
        }
        if (numberParallels[i] == numberSlots[i]) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            log("Winding winding " + std::to_string(i) + " by putting together turns of the same parallel, as the number of parallels is equal to the number of sections.");
            continue;
        }
        if (numberParallels[i] % numberSlots[i] == 0) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            log("Winding winding " + std::to_string(i) + " by putting together turns of the same parallel, as the number of parallels is divisible by the number of sections.");
            continue;
        }
        if (numberTurns[i] % numberSlots[i] == 0) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
            log("Winding winding " + std::to_string(i) + " by putting together parallels of the same turn, as the number of turns is divisible by the number of sections.");
            continue;
        }
        windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
        log("Winding winding " + std::to_string(i) + " by putting together turns of the same parallel, as the number of parallels is smaller than the number of turns.");

    }
    return windByConsecutiveTurns;
}

WindingStyle Coil::wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, size_t numberSlots) {
    if (numberTurns == numberSlots) {
        log("Winding layer by putting together parallels of the same turn, as the number of turns is equal to the number of layers.");
        return WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS;
    }
    if (numberParallels == numberSlots) {
        log("Winding layer by putting together turns of the same parallel, as the number of parallels is equal to the number of layers.");
        return WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
    }
    if (numberParallels % numberSlots == 0) {
        log("Winding layer by putting together turns of the same parallel, as the number of parallels is divisible by the number of layers.");
        return WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
    }
    if (numberTurns % numberSlots == 0) {
        log("Winding layer by putting together parallels of the same turn, as the number of turns is divisible by the number of layers.");
        return WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS;
    }
    log("Winding layer by putting together turns of the same parallel, as neither the number of parallels nor the number of turns is divisible by the number of turns.");
    return WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
}


uint64_t Coil::get_number_turns(size_t windingIndex) {
    return get_functional_description()[windingIndex].get_number_turns();
}

uint64_t Coil::get_number_parallels(size_t windingIndex) {
    return get_functional_description()[windingIndex].get_number_parallels();
}

uint64_t Coil::get_number_turns(Section section) {
    uint64_t physicalTurnsInSection = 0;
    auto partialWinding = section.get_partial_windings()[0];  // TODO: Support multiwinding in layers
    auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
        physicalTurnsInSection += round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
    }
    return physicalTurnsInSection;
}

uint64_t Coil::get_number_turns(Layer layer) {
    uint64_t physicalTurnsInLayer = 0;
    auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
    auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
        physicalTurnsInLayer += round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
    }
    return physicalTurnsInLayer;
}

std::string Coil::get_name(size_t windingIndex) {
    return get_functional_description()[windingIndex].get_name();
}

std::vector<uint64_t> Coil::get_number_turns() {
    std::vector<uint64_t> numberTurns;
    for (auto & winding : get_functional_description()) {
        numberTurns.push_back(winding.get_number_turns());
    }
    return numberTurns;
}

void Coil::set_number_turns(std::vector<uint64_t> numberTurns) {
    for (size_t i=0; i< get_functional_description().size(); ++i) {
        get_mutable_functional_description()[i].set_number_turns(numberTurns[i]);
    }
}

std::vector<IsolationSide> Coil::get_isolation_sides() {
    std::vector<IsolationSide> isolationSides;
    for (auto & winding : get_functional_description()) {
        isolationSides.push_back(winding.get_isolation_side());
    }
    return isolationSides;
}

void Coil::set_isolation_sides(std::vector<IsolationSide> isolationSides) {
    for (size_t i=0; i< get_functional_description().size(); ++i) {
        get_mutable_functional_description()[i].set_isolation_side(isolationSides[i]);
    }
}

std::vector<Layer> Coil::get_layers_by_section(std::string sectionName) {
    auto layers = get_layers_description().value();
    std::vector<Layer> foundLayers;
    for (auto & layer : layers) {
        auto layerSectionName = layer.get_section().value();
        if (layerSectionName == sectionName) {
            foundLayers.push_back(layer);
        }
    }
    return foundLayers;
}

std::vector<Turn> Coil::get_turns_by_layer(std::string layerName) {
    auto turns = get_turns_description().value();
    std::vector<Turn> foundTurns;
    for (auto & turn : turns) {
        auto turnLayerName = turn.get_layer().value();
        if (turnLayerName == layerName) {
            foundTurns.push_back(turn);
        }
    }
    return foundTurns;
}

std::vector<Turn> Coil::get_turns_by_winding(std::string windingName) {
    auto turns = get_turns_description().value();
    std::vector<Turn> foundTurns;
    for (auto & turn : turns) {
        auto turnSectionName = turn.get_winding();
        if (turnSectionName == windingName) {
            foundTurns.push_back(turn);
        }
    }
    return foundTurns;
}

std::vector<Turn> Coil::get_turns_by_section(std::string sectionName) {
    auto turns = get_turns_description().value();
    std::vector<Turn> foundTurns;
    for (auto & turn : turns) {
        auto turnSectionName = turn.get_section().value();
        if (turnSectionName == sectionName) {
            foundTurns.push_back(turn);
        }
    }
    return foundTurns;
}

std::vector<std::string> Coil::get_layers_names_by_winding(std::string windingName) {
    auto layers = get_layers_description().value();
    std::vector<std::string> foundLayers;
    for (auto & layer : layers) {
        auto layerWindings = layer.get_partial_windings();
        for (auto & winding : layerWindings) {
            if (winding.get_winding() == windingName) {
                foundLayers.push_back(layer.get_name());
                break;
            }
        }
    }
    return foundLayers;
}

std::vector<std::string> Coil::get_layers_names_by_section(std::string sectionName) {
    auto layers = get_layers_description().value();
    std::vector<std::string> foundLayers;
    for (auto & layer : layers) {
        auto layerSectionName = layer.get_section().value();
        if (layerSectionName == sectionName) {
            foundLayers.push_back(layer.get_name());
        }
    }
    return foundLayers;
}

std::vector<std::string> Coil::get_turns_names_by_layer(std::string layerName) {
    auto turns = get_turns_description().value();
    std::vector<std::string> foundTurns;
    for (auto & turn : turns) {
        auto turnLayerName = turn.get_layer().value();
        if (turnLayerName == layerName) {
            foundTurns.push_back(turn.get_name());
        }
    }
    return foundTurns;
}

std::vector<std::string> Coil::get_turns_names_by_winding(std::string windingName) {
    auto turns = get_turns_description().value();
    std::vector<std::string> foundTurns;
    for (auto & turn : turns) {
        auto turnWindingName = turn.get_winding();
        if (turnWindingName == windingName) {
            foundTurns.push_back(turn.get_name());
        }
    }
    return foundTurns;
}

std::vector<std::string> Coil::get_turns_names_by_section(std::string sectionName) {
    auto turns = get_turns_description().value();
    std::vector<std::string> foundTurns;
    for (auto & turn : turns) {
        auto turnSectionName = turn.get_section().value();
        if (turnSectionName == sectionName) {
            foundTurns.push_back(turn.get_name());
        }
    }
    return foundTurns;
}
    
std::vector<size_t> Coil::get_turns_indexes_by_layer(std::string layerName) {
    auto turns = get_turns_description().value();
    std::vector<size_t> foundTurns;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turnLayerName = turns[turnIndex].get_layer().value();
        if (turnLayerName == layerName) {
            foundTurns.push_back(turnIndex);
        }
    }
    return foundTurns;
}

std::vector<size_t> Coil::get_turns_indexes_by_winding(std::string windingName) {
    auto turns = get_turns_description().value();
    std::vector<size_t> foundTurns;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turnWindingName = turns[turnIndex].get_winding();
        if (turnWindingName == windingName) {
            foundTurns.push_back(turnIndex);
        }
    }
    return foundTurns;
}

std::vector<size_t> Coil::get_turns_indexes_by_section(std::string sectionName) {
    auto turns = get_turns_description().value();
    std::vector<size_t> foundTurns;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turnSectionName = turns[turnIndex].get_section().value();
        if (turnSectionName == sectionName) {
            foundTurns.push_back(turnIndex);
        }
    }
    return foundTurns;
}

std::vector<Section> Coil::get_sections_by_group(std::string groupName) {
    auto sections = get_sections_description().value();
    std::vector<Section> foundSections;
    for (auto & section : sections) {
        if (section.get_group()) {
            auto sectionSectionGroup = section.get_group().value();
            if (sectionSectionGroup == groupName) {
                foundSections.push_back(section);
            }
        }
    }
    return foundSections;
}

const std::vector<Section> Coil::get_sections_by_type(ElectricalType electricalType) const {
    auto sections = get_sections_description().value();
    std::vector<Section> foundSections;
    for (auto & section : sections) {
        auto sectionSectionType = section.get_type();
        if (sectionSectionType == electricalType) {
            foundSections.push_back(section);
        }
    }
    return foundSections;
}

const std::vector<Section> Coil::get_sections_by_winding(std::string windingName) const {
    auto sections = get_sections_description().value();
    std::vector<Section> foundSections;
    for (auto & section : sections) {
        for (auto & winding : section.get_partial_windings()) {
            if (winding.get_winding() == windingName) {
                foundSections.push_back(section);
            }
        }
    }
    return foundSections;
}

const Section Coil::get_section_by_name(std::string name) const {
    auto sections = get_sections_description().value();
    for (auto & section : sections) {
        if (section.get_name() == name) {
            return section;
        }
    }
    throw std::runtime_error("Not found section with name:" + name);
}

const Layer Coil::get_layer_by_name(std::string name) const {
    if (!get_layers_description()) {
        throw std::runtime_error("Coil is missing layers description");
    }

    auto layers = get_layers_description().value();
    for (auto & layer : layers) {
        if (layer.get_name() == name) {
            return layer;
        }
    }
    throw std::runtime_error("Not found layer with name:" + name);
}

const Turn Coil::get_turn_by_name(std::string name) const {
    if (!get_turns_description()) {
        throw std::runtime_error("Turns description not set, did you forget to wind?");
    }
    auto turns = get_turns_description().value();
    for (auto & turn : turns) {
        if (turn.get_name() == name) {
            return turn;
        }
    }
    throw std::runtime_error("Not found turn with name:" + name);
}

const std::vector<Layer> Coil::get_layers_by_type(ElectricalType electricalType) const {
    auto layers = get_layers_description().value();
    std::vector<Layer> foundLayers;
    for (auto & layer : layers) {
        auto layerSectionType = layer.get_type();
        if (layerSectionType == electricalType) {
            foundLayers.push_back(layer);
        }
    }
    return foundLayers;
}

std::vector<Layer> Coil::get_layers_by_winding_index(size_t windingIndex) {
    auto layers = get_layers_by_type(ElectricalType::CONDUCTION);
    std::vector<Layer> foundLayers;
    for (auto & layer : layers) {
        auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
        auto winding = get_winding_by_name(partialWinding.get_winding());
        auto layerWindingIndex = get_winding_index_by_name(partialWinding.get_winding());
        if (layerWindingIndex == windingIndex) {
            foundLayers.push_back(layer);
        }
    }
    return foundLayers;
}

std::vector<uint64_t> Coil::get_number_parallels() {
    std::vector<uint64_t> numberParallels;
    for (auto & winding : get_functional_description()) {
        numberParallels.push_back(winding.get_number_parallels());
    }
    return numberParallels;
}

void Coil::set_number_parallels(std::vector<uint64_t> numberParallels){
    for (size_t i=0; i< get_functional_description().size(); ++i) {
        get_mutable_functional_description()[i].set_number_parallels(numberParallels[i]);
    }
}

CoilFunctionalDescription Coil::get_winding_by_name(std::string name) {
    for (auto& CoilFunctionalDescription : get_functional_description()) {
        if (CoilFunctionalDescription.get_name() == name) {
            return CoilFunctionalDescription;
        }
    }
    throw std::runtime_error("No such a winding name: " + name);
}

size_t Coil::get_winding_index_by_name(std::string name) {
    for (size_t i=0; i<get_functional_description().size(); ++i) {
        if (get_functional_description()[i].get_name() == name) {
            return i;
        }
    }
    throw std::runtime_error("No such a winding name: " + name);
}

size_t Coil::get_turn_index_by_name(std::string name) {
    if (!get_turns_description()) {
        throw std::runtime_error("Turns description not set, did you forget to wind?");
    }
    auto turns = get_turns_description().value();
    for (size_t i=0; i<turns.size(); ++i) {
        if (turns[i].get_name() == name) {
            return i;
        }
    }
    throw std::runtime_error("No such a turn name: " + name);
}

size_t Coil::get_layer_index_by_name(std::string name) {
    if (!get_layers_description()) {
        throw std::runtime_error("Layers description not set, did you forget to wind?");
    }
    auto layers = get_layers_description().value();
    for (size_t i=0; i<layers.size(); ++i) {
        if (layers[i].get_name() == name) {
            return i;
        }
    }
    throw std::runtime_error("No such a layer name: " + name);
}

size_t Coil::get_section_index_by_name(std::string name) {
    if (!get_sections_description()) {
        throw std::runtime_error("Sections description not set, did you forget to wind?");
    }
    auto sections = get_sections_description().value();
    for (size_t i=0; i<sections.size(); ++i) {
        if (sections[i].get_name() == name) {
            return i;
        }
    }
    throw std::runtime_error("No such a section name: " + name);
}

bool Coil::are_sections_and_layers_fitting() {
    bool windTurns = true;
    if (!get_sections_description()) {
        return false;
    }
    if (!get_layers_description()) {
        return false;
    }
    auto sections = get_sections_description().value();
    auto layers = get_layers_description().value();

    for (auto& section: sections) {

        if (roundFloat(section.get_filling_factor().value(), 6) > 1 || roundFloat(overlapping_filling_factor(section), 6) > 1 || roundFloat(contiguous_filling_factor(section), 6) > 1) {
            windTurns = false;
        }
    }
    for (auto& layer: layers) {
        if (roundFloat(layer.get_filling_factor().value(), 6) > 1) {
            windTurns = false;
        }
    }
    return windTurns;
}

double Coil::overlapping_filling_factor(Section section) {
    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto layers = get_layers_by_section(section.get_name());

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        double sectionWidth = section.get_dimensions()[0];
        double layersWidth = 0;
        for (auto& layer : layers) {
            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                layersWidth += layer.get_dimensions()[0];
            }
            else {
                layersWidth = std::max(layersWidth, layer.get_dimensions()[0]);
            }
        }
        return layersWidth / sectionWidth;
    }
    else {
        double sectionRadialHeight = section.get_dimensions()[0];
        double layersRadialHeight = 0;
        for (auto& layer : layers) {
            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                layersRadialHeight += layer.get_dimensions()[0];
            }
            else {
                layersRadialHeight = std::max(layersRadialHeight, layer.get_dimensions()[0]);
            }
        }
        return layersRadialHeight / sectionRadialHeight;
    }
}

double Coil::contiguous_filling_factor(Section section) {
    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto layers = get_layers_by_section(section.get_name());

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        double sectionHeight = section.get_dimensions()[1];
        double layersHeight = 0;
        for (auto& layer : layers) {
            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                layersHeight = std::max(layersHeight, layer.get_dimensions()[1]);
            }
            else {
                layersHeight += layer.get_dimensions()[1];
            }
        }
        return layersHeight / sectionHeight;
    }
    else {
        double sectionAngle = section.get_dimensions()[1];
        double layersAngle = 0;
        for (auto& layer : layers) {
            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                layersAngle = std::max(layersAngle, layer.get_dimensions()[1]);
            }
            else {
                layersAngle += layer.get_dimensions()[1];
            }
        }
        return layersAngle / sectionAngle;

    }
}

std::pair<uint64_t, std::vector<double>> get_parallels_proportions(size_t slotIndex, size_t slots, uint64_t numberTurns, uint64_t numberParallels, 
                                                                   std::vector<double> remainingParallelsProportion, WindingStyle windByConsecutiveTurns,
                                                                   std::vector<double> totalParallelsProportion, double slotRelativeProportion=1,
                                                                   std::optional<double> slotAbsolutePhysicalTurns = std::nullopt) {
    uint64_t physicalTurnsThisSlot = 0;
    std::vector<double> slotParallelsProportion(numberParallels, 0);
    if (windByConsecutiveTurns == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
        size_t remainingPhysicalTurns = 0;
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            remainingPhysicalTurns += round(remainingParallelsProportion[parallelIndex] * numberTurns);
        }
        if (slotAbsolutePhysicalTurns)
            physicalTurnsThisSlot = slotAbsolutePhysicalTurns.value();
        else
            physicalTurnsThisSlot = std::min(uint64_t(remainingPhysicalTurns), uint64_t(ceil(double(remainingPhysicalTurns) / (slots - slotIndex) * slotRelativeProportion)));
        uint64_t remainingPhysicalTurnsThisSection = physicalTurnsThisSlot;

        size_t currentParallel = 0;
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            if (remainingParallelsProportion[parallelIndex] > 0) {
                currentParallel = parallelIndex;
                break;
            }
        }

        while (remainingPhysicalTurnsThisSection > 0) {
            uint64_t numberTurnsToFitInCurrentParallel = round(remainingParallelsProportion[currentParallel] * numberTurns);
            if (remainingPhysicalTurnsThisSection >= numberTurnsToFitInCurrentParallel) {
                remainingPhysicalTurnsThisSection -= numberTurnsToFitInCurrentParallel;
                slotParallelsProportion[currentParallel] = double(numberTurnsToFitInCurrentParallel) / numberTurns;
                currentParallel++;
            }
            else {
                double proportionParallelsThisSection = double(remainingPhysicalTurnsThisSection) / numberTurns;
                slotParallelsProportion[currentParallel] += proportionParallelsThisSection;
                remainingPhysicalTurnsThisSection = 0;
            }
        }
    }
    else {
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            double remainingSlots = slots - slotIndex;
            double remainingTurnsBeforeThisParallel = numberTurns * remainingParallelsProportion[parallelIndex];
            double numberTurnsToAddToCurrentParallel = ceil(roundFloat(remainingTurnsBeforeThisParallel / remainingSlots * slotRelativeProportion, 10));
            // double numberTurnsToAddToCurrentParallel = ceil(numberTurns * totalParallelsProportion[parallelIndex] / slots);
            double remainingTurnsAfterThisParallel = remainingTurnsBeforeThisParallel - numberTurnsToAddToCurrentParallel;
            double remainingSlotsAfterThisOne = remainingSlots - 1;
            if (remainingTurnsAfterThisParallel < remainingSlotsAfterThisOne) {
                numberTurnsToAddToCurrentParallel = ceil(roundFloat(remainingTurnsBeforeThisParallel / remainingSlots, 10));
            }
            double proportionParallelsThisSection = std::min(remainingParallelsProportion[parallelIndex], numberTurnsToAddToCurrentParallel / numberTurns);
            physicalTurnsThisSlot += numberTurnsToAddToCurrentParallel;
            slotParallelsProportion[parallelIndex] = proportionParallelsThisSection;
        }
    }

    return {physicalTurnsThisSlot, slotParallelsProportion};
}

double get_area_used_in_wires(OpenMagnetics::Wire wire, uint64_t physicalTurns) {
    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        double wireDiameter = wire.get_maximum_outer_width();
        return physicalTurns * pow(wireDiameter, 2);
    }
    else {
        double wireWidth = wire.get_maximum_outer_width();
        double wireHeight = wire.get_maximum_outer_height();
        return physicalTurns * wireWidth * wireHeight;
    }
}

void Coil::set_insulation_layers(std::map<std::pair<size_t, size_t>, std::vector<Layer>> insulationLayers) {
    _insulationLayers = insulationLayers;
}

bool Coil::calculate_custom_thickness_insulation(double thickness) {
    // Insulation layers just for mechanical reasons, one layer between sections at least
    auto wirePerWinding = get_wires();

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto layersOrientation = _layersOrientation;

    // TODO: Properly think about insulation layers with weird windings
    auto windingOrientation = get_winding_orientation();

    if (windingOrientation == WindingOrientation::CONTIGUOUS && _layersOrientation == WindingOrientation::OVERLAPPING) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::CONTIGUOUS;
        }
    }
    if (windingOrientation == WindingOrientation::OVERLAPPING && _layersOrientation == WindingOrientation::CONTIGUOUS) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::OVERLAPPING;
        }
    }

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            // if (leftTopWindingIndex == rightBottomWindingIndex) {
            //     continue;
            // }
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            CoilSectionInterface coilSectionInterface;
            coilSectionInterface.set_number_layers_insulation(1);
            InsulationMaterial defaultInsulationMaterial = find_insulation_material_by_name(defaults.defaultInsulationMaterial);
            coilSectionInterface.set_solid_insulation_thickness(defaultInsulationMaterial.get_thinner_tape_thickness());
            coilSectionInterface.set_total_margin_tape_distance(0);
            coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);

            _insulationLayers[windingsMapKey] = std::vector<Layer>();
            _coilSectionInterfaces[windingsMapKey] = coilSectionInterface;

            Layer layer;
            layer.set_partial_windings(std::vector<PartialWinding>{});
            // layer.set_section(section.get_name());
            layer.set_type(ElectricalType::INSULATION);
            layer.set_name("temp");
            layer.set_orientation(layersOrientation);
            layer.set_turns_alignment(CoilAlignment::SPREAD); // HARDCODED, maybe in the future configure for shields made of turns?

            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                double windingWindowHeight = windingWindows[0].get_height().value();
                double windingWindowWidth = windingWindows[0].get_width().value();
                if (layersOrientation == WindingOrientation::OVERLAPPING) {
                    layer.set_dimensions(std::vector<double>{thickness, windingWindowHeight});
                }
                else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
                    layer.set_dimensions(std::vector<double>{windingWindowWidth, thickness});
                }
            }
            else {
                layer.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    layer.set_dimensions(std::vector<double>{thickness, windingWindowAngle});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    double tapeThicknessInAngle = wound_distance_to_angle(thickness, windingWindowRadialHeight);
                    layer.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                }
            }
            // layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
            layer.set_filling_factor(1);
            _insulationLayers[windingsMapKey].push_back(layer);

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(layersOrientation);
            section.set_type(ElectricalType::INSULATION);

            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                section.set_coordinate_system(CoordinateSystem::CARTESIAN);
                double windingWindowHeight = windingWindows[0].get_height().value();
                double windingWindowWidth = windingWindows[0].get_width().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{thickness, windingWindowHeight});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    section.set_dimensions(std::vector<double>{windingWindowWidth, thickness});
                }
            }
            else {
                section.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{thickness, windingWindowAngle});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    double tapeThicknessInAngle = wound_distance_to_angle(thickness, windingWindowRadialHeight);
                    section.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                }
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

bool Coil::calculate_mechanical_insulation() {
    // Insulation layers just for mechanical reasons, one layer between sections at least
    auto wirePerWinding = get_wires();

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto layersOrientation = _layersOrientation;

    // TODO: Properly think about insulation layers with weird windings
    auto windingOrientation = get_winding_orientation();

    if (windingOrientation == WindingOrientation::CONTIGUOUS && _layersOrientation == WindingOrientation::OVERLAPPING) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::CONTIGUOUS;
        }
    }
    if (windingOrientation == WindingOrientation::OVERLAPPING && _layersOrientation == WindingOrientation::CONTIGUOUS) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::OVERLAPPING;
        }
    }

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            if (leftTopWindingIndex == rightBottomWindingIndex) {
                continue;
            }
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            CoilSectionInterface coilSectionInterface;
            coilSectionInterface.set_number_layers_insulation(1);
            InsulationMaterial defaultInsulationMaterial = find_insulation_material_by_name(defaults.defaultInsulationMaterial);
            coilSectionInterface.set_solid_insulation_thickness(defaultInsulationMaterial.get_thinner_tape_thickness());
            coilSectionInterface.set_total_margin_tape_distance(0);
            coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);

            _insulationLayers[windingsMapKey] = std::vector<Layer>();
            _coilSectionInterfaces[windingsMapKey] = coilSectionInterface;

            for (size_t layerIndex = 0; layerIndex < coilSectionInterface.get_number_layers_insulation(); ++layerIndex) {
                Layer layer;
                layer.set_partial_windings(std::vector<PartialWinding>{});
                // layer.set_section(section.get_name());
                layer.set_type(ElectricalType::INSULATION);
                layer.set_name("temp");
                layer.set_orientation(layersOrientation);
                layer.set_turns_alignment(CoilAlignment::SPREAD); // HARDCODED, maybe in the future configure for shields made of turns?

                if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                    layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                    double windingWindowHeight = windingWindows[0].get_height().value();
                    double windingWindowWidth = windingWindows[0].get_width().value();
                    if (layersOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{defaultInsulationMaterial.get_thinner_tape_thickness(), windingWindowHeight});
                    }
                    else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
                        layer.set_dimensions(std::vector<double>{windingWindowWidth, defaultInsulationMaterial.get_thinner_tape_thickness()});
                    }
                }
                else {
                    layer.set_coordinate_system(CoordinateSystem::POLAR);
                    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                    double windingWindowAngle = windingWindows[0].get_angle().value();
                    if (windingOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{defaultInsulationMaterial.get_thinner_tape_thickness(), windingWindowAngle});
                    }
                    else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                        double tapeThicknessInAngle = wound_distance_to_angle(defaultInsulationMaterial.get_thinner_tape_thickness(), windingWindowRadialHeight);
                        layer.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                    }
                }
                // layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layer.set_filling_factor(1);
                _insulationLayers[windingsMapKey].push_back(layer);
            }
            // _insulationLayersLog[windingsMapKey] = "Adding " + std::to_string(coilSectionInterface.get_number_layers_insulation()) + " insulation layers, as we need a thickness of " + std::to_string(smallestInsulationThicknessCoveringRemaining * 1000) + " mm to achieve " + neededInsulationTypeString + " insulation";

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(layersOrientation);
            section.set_type(ElectricalType::INSULATION);

            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                section.set_coordinate_system(CoordinateSystem::CARTESIAN);
                double windingWindowHeight = windingWindows[0].get_height().value();
                double windingWindowWidth = windingWindows[0].get_width().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowHeight});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    section.set_dimensions(std::vector<double>{windingWindowWidth, coilSectionInterface.get_solid_insulation_thickness()});
                }
            }
            else {
                section.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowAngle});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    double tapeThicknessInAngle = wound_distance_to_angle(coilSectionInterface.get_solid_insulation_thickness(), windingWindowRadialHeight);
                    section.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                }
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

bool Coil::calculate_insulation(bool simpleMode) {
    auto inputs = _inputs.value();

    if (!inputs.get_design_requirements().get_insulation()) {
        return false;
    }

    auto wirePerWinding = get_wires();

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto layersOrientation = _layersOrientation;
    auto windingOrientation = get_winding_orientation();

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            if (leftTopWindingIndex == rightBottomWindingIndex) {
                continue;
            }
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            CoilSectionInterface coilSectionInterface;
            coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::INSULATING);
            InsulationMaterial chosenInsulationMaterial;

            if (simpleMode) {
                InsulationMaterial defaultInsulationMaterial = find_insulation_material_by_name(defaults.defaultInsulationMaterial);
                chosenInsulationMaterial = defaultInsulationMaterial;
                coilSectionInterface.set_solid_insulation_thickness(defaultInsulationMaterial.get_thinner_tape_thickness());
                if (settings->get_coil_allow_margin_tape()) {
                    coilSectionInterface.set_number_layers_insulation(1);
                    coilSectionInterface.set_total_margin_tape_distance(_standardCoordinator.calculate_creepage_distance(inputs, true));
                }
                else {
                    coilSectionInterface.set_number_layers_insulation(3);
                    coilSectionInterface.set_total_margin_tape_distance(0);
                }
            }
            else {
                coilSectionInterface.set_solid_insulation_thickness(DBL_MAX);
                coilSectionInterface.set_number_layers_insulation(ULONG_MAX);

                if (insulationMaterialDatabase.empty()) {
                    load_insulation_materials();
                }

                for (auto& insulationMaterial : insulationMaterialDatabase) {
                    auto auxCoilSectionInterface = _standardCoordinator.calculate_coil_section_interface_layers(inputs, wireLeftTopWinding, wireRightBottomWinding, insulationMaterial.second);
                    if (auxCoilSectionInterface) {
                        if (auxCoilSectionInterface.value().get_solid_insulation_thickness() < coilSectionInterface.get_solid_insulation_thickness()) {
                            coilSectionInterface = auxCoilSectionInterface.value();
                            chosenInsulationMaterial = insulationMaterial.second;
                        }
                    }
                }

                if (coilSectionInterface.get_number_layers_insulation() == ULONG_MAX) {
                    throw std::runtime_error("No insulation found with current requirements");
                }
            }

            _insulationLayers[windingsMapKey] = std::vector<Layer>();
            _coilSectionInterfaces[windingsMapKey] = coilSectionInterface;

            for (size_t layerIndex = 0; layerIndex < coilSectionInterface.get_number_layers_insulation(); ++layerIndex) {
                Layer layer;
                layer.set_partial_windings(std::vector<PartialWinding>{});
                // layer.set_section(section.get_name());
                layer.set_type(ElectricalType::INSULATION);
                layer.set_name("temp");
                layer.set_orientation(_layersOrientation);
                layer.set_turns_alignment(CoilAlignment::SPREAD); // HARDCODED, maybe in the future configure for shields made of turns?

                if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                    layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                    double windingWindowHeight = windingWindows[0].get_height().value();
                    double windingWindowWidth = windingWindows[0].get_width().value();
                    if (layersOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{chosenInsulationMaterial.get_thinner_tape_thickness(), windingWindowHeight});
                    }
                    else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
                        layer.set_dimensions(std::vector<double>{windingWindowWidth, chosenInsulationMaterial.get_thinner_tape_thickness()});
                    }
                }
                else {
                    layer.set_coordinate_system(CoordinateSystem::POLAR);
                    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                    double windingWindowAngle = windingWindows[0].get_angle().value();
                    if (windingOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{chosenInsulationMaterial.get_thinner_tape_thickness(), windingWindowAngle});
                    }
                    else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                        double tapeThicknessInAngle = wound_distance_to_angle(chosenInsulationMaterial.get_thinner_tape_thickness(), windingWindowRadialHeight);
                        layer.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                    }
                }

                layer.set_filling_factor(1);
                _insulationLayers[windingsMapKey].push_back(layer);
            }
            // _insulationLayersLog[windingsMapKey] = "Adding " + std::to_string(coilSectionInterface.get_number_layers_insulation()) + " insulation layers, as we need a thickness of " + std::to_string(smallestInsulationThicknessCoveringRemaining * 1000) + " mm to achieve " + neededInsulationTypeString + " insulation";

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(_layersOrientation);
            section.set_type(ElectricalType::INSULATION);

            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                section.set_coordinate_system(CoordinateSystem::CARTESIAN);
                double windingWindowHeight = windingWindows[0].get_height().value();
                double windingWindowWidth = windingWindows[0].get_width().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowHeight});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    section.set_dimensions(std::vector<double>{windingWindowWidth, coilSectionInterface.get_solid_insulation_thickness()});
                }
            }
            else {
                section.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowAngle});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    double tapeThicknessInAngle = wound_distance_to_angle(coilSectionInterface.get_solid_insulation_thickness(), windingWindowRadialHeight);
                    section.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                }
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

std::vector<std::pair<size_t, double>> Coil::get_ordered_sections(double spaceForSections, std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    std::vector<std::pair<size_t, double>> orderedSections;
    double numberWindings = get_functional_description().size();
    auto numberSectionsPerWinding = std::vector<size_t>(numberWindings, 0);
    for (auto windingIndex : pattern) {
        if (windingIndex >= numberWindings) {
            throw std::invalid_argument("Winding index does not exist in winding");
        }
        numberSectionsPerWinding[windingIndex] += repetitions;
    }

    for (size_t repetitionIndex = 0; repetitionIndex < repetitions; ++repetitionIndex) {
        for (auto windingIndex : pattern) {
            if (roundFloat(proportionPerWinding[windingIndex], 6) > 1) {
                throw std::invalid_argument("proportionPerWinding[windingIndex] cannot be greater than 1: " + std::to_string(proportionPerWinding[windingIndex]));
            }
            double spaceForSection = roundFloat(spaceForSections * proportionPerWinding[windingIndex] / numberSectionsPerWinding[windingIndex], 9);
            orderedSections.push_back({windingIndex, spaceForSection});
        }
    }

    return orderedSections;
}

std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> Coil::add_insulation_to_sections(std::vector<std::pair<size_t, double>> orderedSections){
    std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation;
    auto windingOrientation = get_winding_orientation();
    for (size_t sectionIndex = 1; sectionIndex < orderedSections.size(); ++sectionIndex) {
        auto leftWindingIndex = orderedSections[sectionIndex - 1].first;
        auto rightWindingIndex = orderedSections[sectionIndex].first;
        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex}; 
        if (!_insulationSections.contains(windingsMapKey)) {
            continue;
        }
        auto currentSpaceForLeftSection = orderedSections[sectionIndex - 1].second;
        auto currentSpaceForRightSection = orderedSections[sectionIndex].second;

        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            std::pair<size_t, double> leftSectionInfo = {leftWindingIndex, currentSpaceForLeftSection - _insulationSections[windingsMapKey].get_dimensions()[0] / 2};
            orderedSections[sectionIndex - 1] = leftSectionInfo;
            std::pair<size_t, double> rightSectionInfo = {rightWindingIndex, currentSpaceForRightSection - _insulationSections[windingsMapKey].get_dimensions()[0] / 2};
            orderedSections[sectionIndex] = rightSectionInfo;
        }
        else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            std::pair<size_t, double> leftSectionInfo = {leftWindingIndex, currentSpaceForLeftSection - _insulationSections[windingsMapKey].get_dimensions()[1] / 2};
            orderedSections[sectionIndex - 1] = leftSectionInfo;
            std::pair<size_t, double> rightSectionInfo = {rightWindingIndex, currentSpaceForRightSection - _insulationSections[windingsMapKey].get_dimensions()[1] / 2};
            orderedSections[sectionIndex] = rightSectionInfo;
        }
    }

    orderedSectionsWithInsulation.push_back({ElectricalType::CONDUCTION, orderedSections[0]});
    for (size_t sectionIndex = 1; sectionIndex < orderedSections.size(); ++sectionIndex) {
        auto leftWindingIndex = orderedSections[sectionIndex - 1].first;
        auto rightWindingIndex = orderedSections[sectionIndex].first;
        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex}; 
        if (_insulationSections.contains(windingsMapKey)) {
            std::pair<size_t, double> insulationSectionInfo;
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[0]};
            }
            else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[1]};
            }

            orderedSectionsWithInsulation.push_back({ElectricalType::INSULATION, insulationSectionInfo});
        }
        orderedSectionsWithInsulation.push_back({ElectricalType::CONDUCTION, orderedSections[sectionIndex]});
    }

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();


    // last insulation layer we compare between last and first
    if (windingOrientation != WindingOrientation::CONTIGUOUS || bobbinWindingWindowShape != WindingWindowShape::RECTANGULAR) {
        // We don't add one in the sections are contiguous, as they end in the bobbin
        auto leftWindingIndex = orderedSections[orderedSections.size() - 1].first;
        auto rightWindingIndex = orderedSections[0].first;
        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex}; 

        if (_insulationSections.contains(windingsMapKey)) {
            std::pair<size_t, double> insulationSectionInfo;
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[0]};
            }
            else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[1]};
            }

            orderedSectionsWithInsulation.push_back({ElectricalType::INSULATION, insulationSectionInfo});
        }
    }

    return orderedSectionsWithInsulation;
}

std::vector<double> Coil::get_proportion_per_winding_based_on_wires() {
    std::vector<double> physicalTurnsAreaPerWinding;
    double totalPhysicalTurnsArea = 0;
    auto wirePerWinding = get_wires();
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex){
        double wireOuterRectangularArea = wirePerWinding[windingIndex].get_maximum_outer_width() * wirePerWinding[windingIndex].get_maximum_outer_height();
        double totalAreaThisWinding = wireOuterRectangularArea * get_functional_description()[windingIndex].get_number_turns() * get_functional_description()[windingIndex].get_number_parallels();
        physicalTurnsAreaPerWinding.push_back(totalAreaThisWinding);
        totalPhysicalTurnsArea += totalAreaThisWinding;
    }
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex){
        physicalTurnsAreaPerWinding[windingIndex] /= totalPhysicalTurnsArea;
    }

    return physicalTurnsAreaPerWinding;
}

std::pair<double, double> get_section_round_dimensions(std::pair<ElectricalType, std::pair<size_t, double>> sectionWithInsulationScaledWithArea,
                                                 WindingOrientation windingOrientation, double windingWindowRadialHeight, double windingWindowAngle) {

    auto sectionInfo = sectionWithInsulationScaledWithArea.second;
    auto spaceForSection = sectionInfo.second;

    double currentSectionRadialHeight = 0;
    double currentSectionAngle = 0;
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        currentSectionRadialHeight = spaceForSection;
        currentSectionAngle = windingWindowAngle;
    }
    else {
        currentSectionRadialHeight = windingWindowRadialHeight;
        currentSectionAngle = spaceForSection;
    }

    return {currentSectionRadialHeight, currentSectionAngle};
}

std::vector<double> get_physical_turns_proportions(std::vector<int64_t> physicalTurns) {
    std::vector<double> physicalTurnsProportions;
    double average = 0;
    for (size_t index = 0; index < physicalTurns.size(); ++index) {
        average += double(physicalTurns[index]);
    }
    average /= physicalTurns.size();

    for (size_t index = 0; index < physicalTurns.size(); ++index) {
        if (index < physicalTurns.size() - 1)
            physicalTurnsProportions.push_back(double(physicalTurns[index]) / average);
        else
            physicalTurnsProportions.push_back(1 + double(physicalTurns[index]) / average);
    }

    return physicalTurnsProportions;
}

std::vector<double> get_length_proportions(std::vector<double> lengths, std::vector<size_t> windingIndexes) {
    std::vector<size_t> uniqueIndexes;
    for (size_t windingIndex = 0; windingIndex < windingIndexes.size(); ++windingIndex) {
        if(std::find(uniqueIndexes.begin(), uniqueIndexes.end(), windingIndexes[windingIndex]) == uniqueIndexes.end()) {
            uniqueIndexes.push_back(windingIndexes[windingIndex]);
        }
    }

    std::vector<double> lengthProportions;
    std::vector<double> averages(uniqueIndexes.size(), 0);
    std::vector<double> numberSectionsPerWinding(uniqueIndexes.size(), 0);

    for (size_t index = 0; index < lengths.size(); ++index) {
        averages[windingIndexes[index]] += lengths[index];
        numberSectionsPerWinding[windingIndexes[index]]++;
    }

    for (size_t windingIndex = 0; windingIndex < averages.size(); ++windingIndex) {
        averages[windingIndex] /= numberSectionsPerWinding[windingIndex];
    }

    for (size_t index = 0; index < lengths.size(); ++index) {
        if (index < lengths.size() - 1)
            lengthProportions.push_back(lengths[index] / averages[windingIndexes[index]]);
        else
            lengthProportions.push_back(1 + lengths[index] / averages[windingIndexes[index]]);
    }

    return lengthProportions;
}

std::vector<double> get_section_lengths(std::vector<double> currentSectionRadialHeights, std::vector<double> currentSectionAngles, double windingWindowRadialHeight) {
    std::vector<double> sectionLengths;
    double radialHeightIncrease = windingWindowRadialHeight / currentSectionRadialHeights.size();
    for (size_t sectionIndex = 0; sectionIndex < currentSectionRadialHeights.size(); ++sectionIndex) {
        double radius = windingWindowRadialHeight - radialHeightIncrease * sectionIndex - radialHeightIncrease;
        sectionLengths.push_back(2 * std::numbers::pi * radius * currentSectionAngles[sectionIndex] / 360);
    }
    return sectionLengths;
}

std::vector<double> get_section_areas(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulationScaledWithArea, std::vector<double> currentSectionAngles, double windingWindowRadialHeight) {
    std::vector<double> sectionAreas;
    double currentRadius = windingWindowRadialHeight;
    size_t currentConductionSectionIndex = 0;
    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulationScaledWithArea.size(); ++sectionIndex) {
        if (orderedSectionsWithInsulationScaledWithArea[sectionIndex].first == ElectricalType::CONDUCTION) {
            auto sectionInfo = orderedSectionsWithInsulationScaledWithArea[sectionIndex].second;
            auto spaceForSection = sectionInfo.second;
            double outerRadius = currentRadius;
            double innerRadius = currentRadius - spaceForSection;
            currentRadius -= spaceForSection;
            sectionAreas.push_back(std::numbers::pi * (pow(outerRadius, 2) - pow(innerRadius, 2)) * currentSectionAngles[currentConductionSectionIndex] / 360);
            currentConductionSectionIndex++;
        }

    }
    return sectionAreas;
}

std::pair<size_t, std::vector<int64_t>> get_number_layers_needed_and_number_physical_turns(double radialHeight, double angle, Wire wire, int64_t physicalTurnsInSection, double windingWindowRadius) {
    int64_t reaminingPhysicalTurnsInSection = physicalTurnsInSection;
    double wireWidth = resolve_dimensional_values(wire.get_maximum_outer_width());
    double wireHeight = resolve_dimensional_values(wire.get_maximum_outer_height());
    double currentRadialHeight = radialHeight;
    double currentRadius;
    if (wire.get_type() == WireType::FOIL){
        throw std::runtime_error("Foil is not supported in toroids");
    }
    if (wire.get_type() == WireType::PLANAR){
        throw std::runtime_error("Planar is not supported in toroids");
    }
    if (wire.get_type() == WireType::RECTANGULAR){
        currentRadius = windingWindowRadius - wireWidth - currentRadialHeight;
    }
    else {
        currentRadius = windingWindowRadius - wireWidth / 2 - currentRadialHeight;
    }
    double sectionAvailableAngle = angle;
    std::vector<int64_t> layerPhysicalTurns;
    size_t numberLayers = 0;
    while (reaminingPhysicalTurnsInSection > 0) {
        double wireAngle = wound_distance_to_angle(wireHeight, std::max(wireWidth, currentRadius));
        int64_t numberTurnsFittingThisLayer = std::max(1.0, floor(sectionAvailableAngle / wireAngle));
        reaminingPhysicalTurnsInSection -= numberTurnsFittingThisLayer;

        layerPhysicalTurns.push_back(numberTurnsFittingThisLayer);
        numberLayers++;
        if (currentRadius > wireWidth) {
            currentRadius -= wireWidth;
        }
    }

    int64_t numberTurnsToCorrect = -reaminingPhysicalTurnsInSection;
    size_t currentIndex = numberLayers - 1;
    while (numberTurnsToCorrect > 0) {
        layerPhysicalTurns[currentIndex]--;
        numberTurnsToCorrect--;
        if (currentIndex == 0)
            currentIndex = numberLayers - 1;
        else
            currentIndex--;
    }

    return {numberLayers, layerPhysicalTurns};
}

std::pair<size_t, std::vector<int64_t>> get_number_layers_needed_and_number_physical_turns(Section section, Wire wire, int64_t physicalTurnsInSection, double windingWindowRadius) {
    return get_number_layers_needed_and_number_physical_turns(section.get_coordinates()[0] - section.get_dimensions()[0] / 2, section.get_dimensions()[1], wire, physicalTurnsInSection, windingWindowRadius);
}

void Coil::apply_margin_tape(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation) {
    if (_marginsPerSection.size() < orderedSectionsWithInsulation.size()) {
        _marginsPerSection = std::vector<std::vector<double>>(orderedSectionsWithInsulation.size(), {0, 0});
    }

    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
        if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
            if (sectionIndex > 0 && !_coilSectionInterfaces.empty()) {

                if (orderedSectionsWithInsulation[sectionIndex - 1].first != ElectricalType::INSULATION) {
                    throw std::runtime_error("There cannot be two sections without insulation in between");
                }
                auto windingIndex = orderedSectionsWithInsulation[sectionIndex].second.first;
                auto previousWindingIndex = orderedSectionsWithInsulation[sectionIndex - 2].second.first;
                auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, windingIndex};
                auto coilSectionInterface = _coilSectionInterfaces[windingsMapKey];
                _marginsPerSection[sectionIndex][0] =  std::max(_marginsPerSection[sectionIndex][0], coilSectionInterface.get_total_margin_tape_distance() / 2);
                _marginsPerSection[sectionIndex][1] =  std::max(_marginsPerSection[sectionIndex][1], coilSectionInterface.get_total_margin_tape_distance() / 2);
                _marginsPerSection[sectionIndex - 2][0] =  std::max(_marginsPerSection[sectionIndex - 2][0], coilSectionInterface.get_total_margin_tape_distance() / 2);
                _marginsPerSection[sectionIndex - 2][1] =  std::max(_marginsPerSection[sectionIndex - 2][1], coilSectionInterface.get_total_margin_tape_distance() / 2);
            }
        }
    }
}

void Coil::equalize_margins(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation) {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
        if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
            if (!_coilSectionInterfaces.empty()) {

                size_t indexForMarginLeftSection = sectionIndex;
                size_t indexForMarginRightSection;
                if (sectionIndex != (orderedSectionsWithInsulation.size() - 2)) {
                    indexForMarginRightSection = sectionIndex + 2;
                }
                else {
                    indexForMarginRightSection = 0;
                }

                auto windingIndex = orderedSectionsWithInsulation[indexForMarginLeftSection].second.first;
                auto previousWindingIndex = orderedSectionsWithInsulation[indexForMarginRightSection].second.first;
                auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, windingIndex};
                auto coilSectionInterface = _coilSectionInterfaces[windingsMapKey];
                double totalMargin = _marginsPerSection[indexForMarginLeftSection][1] + _marginsPerSection[indexForMarginRightSection][0];
                double leftAvailableSpace = orderedSectionsWithInsulation[indexForMarginLeftSection].second.second;
                double rightAvailableSpace = orderedSectionsWithInsulation[indexForMarginRightSection].second.second;
                double totalAvailableSpace = leftAvailableSpace + rightAvailableSpace;
                _marginsPerSection[indexForMarginLeftSection][1] = leftAvailableSpace / totalAvailableSpace * totalMargin;
                _marginsPerSection[indexForMarginRightSection][0] = rightAvailableSpace / totalAvailableSpace * totalMargin;
            }
        }
    }
}

bool Coil::wind_by_sections() {
    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    std::vector<double> proportionPerWinding;
    auto windingOrientation = get_winding_orientation();
    auto sectionAlignment = get_section_alignment();

    if (bobbinWindingWindowShape == WindingWindowShape::ROUND && windingOrientation == WindingOrientation::CONTIGUOUS && sectionAlignment == CoilAlignment::SPREAD ) {
        proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    }
    else {
        proportionPerWinding = get_proportion_per_winding_based_on_wires();
    }
    return wind_by_sections(proportionPerWinding);
}

bool Coil::wind_by_sections(size_t repetitions){
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind_by_sections(proportionPerWinding, pattern, repetitions);
}

bool Coil::wind_by_sections(std::vector<size_t> pattern, size_t repetitions) {
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind_by_sections(proportionPerWinding, pattern, repetitions);
}

bool Coil::wind_by_sections(std::vector<double> proportionPerWinding) {
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    return wind_by_sections(proportionPerWinding, pattern, _interleavingLevel);
}

bool Coil::create_default_group(Bobbin bobbin, WiringTechnology coilType) {
    Group group;
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    group.set_coordinates(windingWindows[0].get_coordinates().value());
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        group.set_dimensions(std::vector<double>{windingWindows[0].get_width().value(), windingWindows[0].get_height().value()});
        group.set_coordinate_system(CoordinateSystem::CARTESIAN);
    }
    else {
        group.set_dimensions(std::vector<double>{windingWindows[0].get_radial_height().value(), windingWindows[0].get_angle().value()});
        group.set_coordinate_system(CoordinateSystem::POLAR);
    }
    group.set_name("Default group");
    std::vector<PartialWinding> partialWindings;

    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        PartialWinding partialWinding;
        partialWinding.set_winding(get_name(windingIndex));
        partialWinding.set_parallels_proportion(std::vector<double>(get_number_parallels(windingIndex), 1));
        partialWindings.push_back(partialWinding);
    }
    group.set_partial_windings(partialWindings);
    group.set_sections_orientation(get_winding_orientation());
    group.set_type(coilType);
    set_groups_description(std::vector<Group>{group});

    return true;
}

bool Coil::wind_by_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    _currentProportionPerWinding = proportionPerWinding;
    _currentPattern = pattern;
    _currentRepetitions = repetitions;

    if (repetitions <= 0) {
        throw std::runtime_error("Interleaving levels must be greater than 0");
    }

    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (!bobbin.get_processed_description()) {
        throw std::runtime_error("Bobbin not processed");
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windingWindows.size() > 1) {
        throw std::runtime_error("Bobbins with more than winding window not implemented yet");
    }
    if (!windingWindows[0].get_sections_orientation()) {
        windingWindows[0].set_sections_orientation(_windingOrientation);
    }
    if (!windingWindows[0].get_sections_alignment()) {
        windingWindows[0].set_sections_alignment(_sectionAlignment);
    }
    bobbinProcessedDescription.set_winding_windows(windingWindows);
    bobbin.set_processed_description(bobbinProcessedDescription);
    set_bobbin(bobbin);

    if (!get_groups_description()) {
        create_default_group(bobbin);
    }

    set_sections_description(std::nullopt);
    set_layers_description(std::nullopt);
    set_turns_description(std::nullopt);

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        return wind_by_rectangular_sections(proportionPerWinding, pattern, repetitions);
    }
    else {
        return wind_by_round_sections(proportionPerWinding, pattern, repetitions);
    }
}

bool Coil::wind_by_rectangular_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    set_sections_description(std::nullopt);
    std::vector<Section> sectionsDescription;

    if (!get_groups_description()) {
        throw std::runtime_error("At least default group must be defined at this point.");
    }

    auto groups = get_groups_description().value();
    std::vector<std::vector<double>> remainingParallelsProportion;

    for (auto group : groups) {
        double availableWidth = group.get_dimensions()[0];
        double availableHeight = group.get_dimensions()[1];
        double spaceForSections = 0;
        auto windingOrientation = group.get_sections_orientation();

        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            spaceForSections = availableWidth;
        }
        else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            spaceForSections = availableHeight;
        }

        auto orderedSections = get_ordered_sections(spaceForSections, proportionPerWinding, pattern, repetitions);

        auto orderedSectionsWithInsulation = add_insulation_to_sections(orderedSections);

        double numberWindings = get_functional_description().size();
        auto numberSectionsPerWinding = std::vector<size_t>(numberWindings, 0);
        auto currentSectionPerWinding = std::vector<size_t>(numberWindings, 0);
        for (auto orderedSection : orderedSectionsWithInsulation) {
            if (orderedSection.first == ElectricalType::CONDUCTION) {
                auto windingIndex = orderedSection.second.first;
                numberSectionsPerWinding[windingIndex]++;
            }
        }

        auto windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(), get_number_parallels(), numberSectionsPerWinding);
       
        auto wirePerWinding = get_wires();
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            remainingParallelsProportion.push_back(std::vector<double>(get_number_parallels(windingIndex), 1));
        }
        double currentSectionCenterWidth = DBL_MAX;
        double currentSectionCenterHeight = DBL_MAX;

        apply_margin_tape(orderedSectionsWithInsulation);

        for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
            if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
                auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
                auto windingIndex = sectionInfo.first;
                auto spaceForSection = sectionInfo.second;

                double currentSectionHeight = 0;
                double currentSectionWidth = 0;
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionWidth = spaceForSection;
                    // currentSectionHeight = availableHeight - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1];
                    currentSectionHeight = availableHeight;

                    if (currentSectionCenterWidth == DBL_MAX) {
                        currentSectionCenterWidth = group.get_coordinates()[0] - availableWidth / 2;
                    }
                    if (currentSectionCenterHeight == DBL_MAX) {
                        currentSectionCenterHeight = group.get_coordinates()[1];
                    }
                }
                else {
                    currentSectionWidth = availableWidth;
                    currentSectionHeight = spaceForSection;
                    if (currentSectionCenterWidth == DBL_MAX) {
                        currentSectionCenterWidth = group.get_coordinates()[0];
                    }
                    if (currentSectionCenterHeight == DBL_MAX) {
                        currentSectionCenterHeight = group.get_coordinates()[1] + availableHeight / 2;
                    }
                }

                PartialWinding partialWinding;
                Section section;
                partialWinding.set_winding(get_name(windingIndex));

                auto parallelsProportions = get_parallels_proportions(currentSectionPerWinding[windingIndex],
                                                                       numberSectionsPerWinding[windingIndex],
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportion[windingIndex],
                                                                       windByConsecutiveTurns[windingIndex],
                                                                       std::vector<double>(get_number_parallels(windingIndex), 1));

                std::vector<double> sectionParallelsProportion = parallelsProportions.second;
                uint64_t physicalTurnsThisSection = parallelsProportions.first;

                partialWinding.set_parallels_proportion(sectionParallelsProportion);
                section.set_name(get_name(windingIndex) +  " section " + std::to_string(currentSectionPerWinding[windingIndex]));
                section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
                section.set_group(group.get_name());
                section.set_type(ElectricalType::CONDUCTION);
                section.set_margin(_marginsPerSection[sectionIndex]);
                section.set_layers_orientation(_layersOrientation);
                section.set_coordinate_system(CoordinateSystem::CARTESIAN);
                
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{currentSectionWidth, currentSectionHeight - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1]});
                }
                else {
                    section.set_dimensions(std::vector<double>{currentSectionWidth - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1], currentSectionHeight});
                }

                if (wirePerWinding[windingIndex].get_type() == WireType::FOIL && !wirePerWinding[windingIndex].get_conducting_height()) {
                    wirePerWinding[windingIndex].cut_foil_wire_to_section(section);
                    get_mutable_functional_description()[windingIndex].set_wire(wirePerWinding[windingIndex]);
                }

                if (wirePerWinding[windingIndex].get_type() == WireType::PLANAR && !wirePerWinding[windingIndex].get_conducting_width()) {
                    wirePerWinding[windingIndex].cut_planar_wire_to_section(section);
                    get_mutable_functional_description()[windingIndex].set_wire(wirePerWinding[windingIndex]);
                }

                if (windingOrientation == WindingOrientation::OVERLAPPING) {

                    if ((section.get_margin().value()[0] + section.get_margin().value()[1] + resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_height())) > currentSectionHeight) {
                        std::string wireType = std::string(magic_enum::enum_name(wirePerWinding[windingIndex].get_type()));
                        return false;
                        // throw std::runtime_error("Margin plus a turn cannot larger than winding window" + 
                        //                          std::string{", margin:"} + std::to_string(section.get_margin().value()[0] + section.get_margin().value()[1]) + 
                        //                          ", wire type: " + wireType + 
                        //                          ", wire height: " + std::to_string(resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_height())) + 
                        //                          ", section height: " + std::to_string(currentSectionHeight)
                        //                          );
                    }
                }
                else {
                    if ((section.get_margin().value()[0] + section.get_margin().value()[1] + resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_width())) > currentSectionWidth) {
                        return false;
                        // throw std::runtime_error("Margin plus a turn cannot larger than winding window" + 
                        //                          std::string{", margin:"} + std::to_string(section.get_margin().value()[0] + section.get_margin().value()[1]) + 
                        //                          ", wire width:" + std::to_string(resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_width())) + 
                        //                          ", section width:" + std::to_string(currentSectionWidth)
                        //                          );
                    }
                }

                if (section.get_dimensions()[0] < 0) {
                    throw std::runtime_error("Something wrong happened in section dimensions 0: " + std::to_string(section.get_dimensions()[0]) +
                                             " availableWidth: " + std::to_string(availableWidth) +
                                             " currentSectionWidth: " + std::to_string(currentSectionWidth) +
                                             " currentSectionHeight: " + std::to_string(currentSectionHeight) + 
                                             " _marginsPerSection[sectionIndex][0]: " + std::to_string(_marginsPerSection[sectionIndex][0])
                                             );
                }
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_coordinates(std::vector<double>{currentSectionCenterWidth + currentSectionWidth / 2, currentSectionCenterHeight, 0});
                }
                else {
                    section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight - currentSectionHeight / 2, 0});
                }

                if (section.get_coordinates()[0] < -1) {
                    throw std::runtime_error("Something wrong happened in section coordiantes 0: " + std::to_string(section.get_coordinates()[0]) +
                                             " currentSectionCenterWidth: " + std::to_string(currentSectionCenterWidth) +
                                             " group.get_coordinates()[0]: " + std::to_string(group.get_coordinates()[0]) +
                                             " group.get_dimensions()[0]: " + std::to_string(group.get_dimensions()[0]) +
                                             " availableWidth: " + std::to_string(availableWidth) +
                                             " currentSectionWidth: " + std::to_string(currentSectionWidth) +
                                             " currentSectionCenterHeight: " + std::to_string(currentSectionCenterHeight)
                                             );
                }

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / (currentSectionWidth * (currentSectionHeight - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1])));
                } 
                else {
                    section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / ((currentSectionWidth - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1]) * currentSectionHeight));
                }
                section.set_winding_style(windByConsecutiveTurns[windingIndex]);
                sectionsDescription.push_back(section);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportion[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
                }

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterWidth += currentSectionWidth;
                }
                else {
                    currentSectionCenterHeight -= currentSectionHeight;
                }

                currentSectionPerWinding[windingIndex]++;
            }
            else {
                if (sectionIndex == 0) {
                    throw std::runtime_error("Insulation layers cannot be the first one (for now)");
                }

                auto previousWindingIndex = orderedSectionsWithInsulation[sectionIndex - 1].second.first;
                size_t nextWindingIndex;
                if (sectionIndex != (orderedSectionsWithInsulation.size() - 1)) {
                    nextWindingIndex = orderedSectionsWithInsulation[sectionIndex + 1].second.first;
                }
                else {
                    nextWindingIndex = orderedSectionsWithInsulation[0].second.first;
                }

                auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, nextWindingIndex};
                if (!_insulationSections.contains(windingsMapKey)) {
                    log(_insulationSectionsLog[windingsMapKey]);
                    continue;
                }

                auto insulationSection = _insulationSections[windingsMapKey];

                insulationSection.set_group(group.get_name());
                insulationSection.set_name("Insulation between " + get_name(previousWindingIndex) + " and " + get_name(nextWindingIndex) + " section " + std::to_string(sectionIndex));
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth + insulationSection.get_dimensions()[0] / 2,
                                                                          currentSectionCenterHeight,
                                                                          0});
                }
                else {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth,
                                                                          currentSectionCenterHeight - insulationSection.get_dimensions()[1] / 2,
                                                                          0});
                }

                sectionsDescription.push_back(insulationSection);
                log(_insulationSectionsLog[windingsMapKey]);

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterWidth += insulationSection.get_dimensions()[0];
                }
                else {
                    currentSectionCenterHeight -= insulationSection.get_dimensions()[1];
                }

            }
        }
    }

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            if (roundFloat(remainingParallelsProportion[windingIndex][parallelIndex], 9) != 0) {
                throw std::runtime_error("There are unassigned parallel proportion in a rectangular section, something went wrong");
            }
        }
    }

    set_sections_description(sectionsDescription);
    return true;
}

void Coil::remove_insulation_if_margin_is_enough(std::vector<std::pair<size_t, double>> orderedSections) {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    size_t multiplier;
    if (_marginsPerSection.size() > orderedSections.size()) {
        multiplier = 2;
    }
    else {
        multiplier = 1;
    }

    for (size_t sectionIndex = 0; sectionIndex < orderedSections.size(); ++sectionIndex) {
        size_t indexForMarginLeftSection = sectionIndex * multiplier;
        size_t indexForMarginRightSection;
        if (sectionIndex != (orderedSections.size() - 1)) {
            indexForMarginRightSection = (sectionIndex + 1) * multiplier;
        }
        else {
            indexForMarginRightSection = 0;
        }
        while (indexForMarginLeftSection >= _marginsPerSection.size() || indexForMarginRightSection >= _marginsPerSection.size()) {
            _marginsPerSection.push_back({0, 0});
        }
    }    

    for (size_t sectionIndex = 0; sectionIndex < orderedSections.size(); ++sectionIndex) {
        size_t indexForMarginLeftSection = sectionIndex * multiplier;
        // size_t indexForMarginLeftSection = sectionIndex;
        size_t indexForMarginRightSection;
        size_t leftWindingIndex = orderedSections[sectionIndex].first;
        size_t rightWindingIndex;
        if (sectionIndex != (orderedSections.size() - 1)) {
            indexForMarginRightSection = (sectionIndex + 1) * multiplier;
            // indexForMarginRightSection = sectionIndex + 1;
            rightWindingIndex = orderedSections[sectionIndex + 1].first;
        }
        else {
            indexForMarginRightSection = 0;
            rightWindingIndex = orderedSections[0].first;
        }

        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex};
        double totalMargin = 0;
        if (_insulationSections.contains(windingsMapKey)) {
            auto coilSectionInterface = _coilSectionInterfaces[windingsMapKey];
            totalMargin = coilSectionInterface.get_total_margin_tape_distance();
        }

        if (_marginsPerSection.size() != 0) {
            double leftMargin = _marginsPerSection[indexForMarginLeftSection][1];
            double rightMargin = _marginsPerSection[indexForMarginRightSection][0];
            totalMargin = std::max(totalMargin, leftMargin + rightMargin);
        }

        double totalMarginAngle = wound_distance_to_angle(totalMargin, windingWindowRadialHeight);

        double totalInsulationDimension = 0;
        if (_insulationSections.contains(windingsMapKey)) {
            totalInsulationDimension = _insulationSections[windingsMapKey].get_dimensions()[1];

            if (totalMarginAngle > totalInsulationDimension) {
                _insulationSections[windingsMapKey].set_dimensions({_insulationSections[windingsMapKey].get_dimensions()[0], 0});
            }
        }
    }
}

bool Coil::wind_by_round_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    set_sections_description(std::nullopt);
    std::vector<Section> sectionsDescription;

    if (!get_groups_description()) {
        throw std::runtime_error("At least default group must be defined at this point.");
    }

    auto groups = get_groups_description().value();
    std::vector<std::vector<double>> remainingParallelsProportion;

    for (auto group : groups) {

        auto bobbin = resolve_bobbin();
        auto bobbinProcessedDescription = bobbin.get_processed_description().value();
        auto windingWindows = bobbinProcessedDescription.get_winding_windows();
        double availableRadialHeight = group.get_dimensions()[0];
        double availableAngle = group.get_dimensions()[1];

        double spaceForSections = 0;
        auto windingOrientation = get_winding_orientation();

        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            spaceForSections = availableRadialHeight;
        }
        else {
            spaceForSections = availableAngle;
        }

        auto orderedSections = get_ordered_sections(spaceForSections, proportionPerWinding, pattern, repetitions);

        if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            remove_insulation_if_margin_is_enough(orderedSections);
        }
        auto orderedSectionsWithInsulation = add_insulation_to_sections(orderedSections);

        double numberWindings = get_functional_description().size();
        auto numberSectionsPerWinding = std::vector<size_t>(numberWindings, 0);
        auto currentSectionPerWinding = std::vector<size_t>(numberWindings, 0);
        for (auto orderedSection : orderedSectionsWithInsulation) {
            if (orderedSection.first == ElectricalType::CONDUCTION) {
                auto windingIndex = orderedSection.second.first;
                numberSectionsPerWinding[windingIndex]++;
            }
        }
        auto windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(), get_number_parallels(), numberSectionsPerWinding);
       
        auto wirePerWinding = get_wires();
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            remainingParallelsProportion.push_back(std::vector<double>(get_number_parallels(windingIndex), 1));
        }
        double currentSectionCenterAngle = DBL_MAX;
        double currentSectionCenterRadialHeight = DBL_MAX;

        apply_margin_tape(orderedSectionsWithInsulation);
        if (settings->get_coil_equalize_margins()) {
            equalize_margins(orderedSectionsWithInsulation);
        }

        std::vector<double> currentSectionRadialHeights;
        std::vector<double> currentSectionAngles;
        std::vector<size_t> windingIndexes;

        for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
            if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
                auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
                auto windingIndex = sectionInfo.first;
                auto aux = get_section_round_dimensions(orderedSectionsWithInsulation[sectionIndex], windingOrientation, availableRadialHeight, availableAngle);
                currentSectionRadialHeights.push_back(aux.first);
                currentSectionAngles.push_back(aux.second);
                windingIndexes.push_back(windingIndex);
            }
        }
        
        // std::vector<double> sectionLengths = get_section_lengths(currentSectionRadialHeights, currentSectionAngles, availableRadialHeight);

        std::vector<double> sectionPhysicalTurnsProportions;
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            std::vector<double> sectionAreas = get_section_areas(orderedSectionsWithInsulation, currentSectionAngles, availableRadialHeight);
            sectionPhysicalTurnsProportions = get_length_proportions(sectionAreas, windingIndexes);
        }
        else {
            sectionPhysicalTurnsProportions = std::vector<double>(orderedSectionsWithInsulation.size(), 1);
        }

        size_t conductingSectionIndex = 0;
        for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
            if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
                auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
                auto windingIndex = sectionInfo.first;

                double currentSectionRadialHeight = currentSectionRadialHeights[conductingSectionIndex];
                double currentSectionAngle = currentSectionAngles[conductingSectionIndex];

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    if (currentSectionCenterRadialHeight == DBL_MAX) {
                        currentSectionCenterRadialHeight = 0;
                    }
                    if (currentSectionCenterAngle == DBL_MAX) {
                        currentSectionCenterAngle = 180;
                    }
                }
                else {
                    if (currentSectionCenterRadialHeight == DBL_MAX) {
                        currentSectionCenterRadialHeight = 0;
                    }
                    if (currentSectionCenterAngle == DBL_MAX) {
                        currentSectionCenterAngle = 0;
                    }
                }
                
                PartialWinding partialWinding;
                Section section;
                partialWinding.set_winding(get_name(windingIndex));

                auto parallelsProportions = get_parallels_proportions(currentSectionPerWinding[windingIndex],
                                                                       numberSectionsPerWinding[windingIndex],
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportion[windingIndex],
                                                                       windByConsecutiveTurns[windingIndex],
                                                                       std::vector<double>(get_number_parallels(windingIndex), 1),
                                                                       sectionPhysicalTurnsProportions[windingIndex]);

                std::vector<double> sectionParallelsProportion = parallelsProportions.second;
                uint64_t physicalTurnsThisSection = parallelsProportions.first;

                double marginAngle0 = 0;
                double marginAngle1 = 0;
                size_t numberLayers = ULONG_MAX;
                size_t prevNumberLayers = 0;

                // We correct the radial height to exactly what we need, so afterwards we can calculate exactly how many turns we need
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    auto aux = get_number_layers_needed_and_number_physical_turns(currentSectionCenterRadialHeight + _marginsPerSection[sectionIndex][0], currentSectionAngle, wirePerWinding[windingIndex], physicalTurnsThisSection, availableRadialHeight);
                    numberLayers = aux.first;
                    currentSectionRadialHeight = numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width();
                }
                else {
                    while (numberLayers != prevNumberLayers) {
                        prevNumberLayers = numberLayers;
                        double currentSectionAngleMinusMargin = currentSectionAngle - marginAngle0 - marginAngle1;
                        auto aux = get_number_layers_needed_and_number_physical_turns(currentSectionCenterRadialHeight, currentSectionAngleMinusMargin, wirePerWinding[windingIndex], physicalTurnsThisSection, availableRadialHeight);
                        numberLayers = aux.first;
                        if (_strict) {
                            currentSectionRadialHeight = numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width();
                        }
                        double lastLayerMaximumRadius = availableRadialHeight - (currentSectionCenterRadialHeight + numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width());
                        if (lastLayerMaximumRadius < 0) {
                            break;
                        }
                        marginAngle0 = wound_distance_to_angle(_marginsPerSection[sectionIndex][0], lastLayerMaximumRadius);
                        marginAngle1 = wound_distance_to_angle(_marginsPerSection[sectionIndex][1], lastLayerMaximumRadius);
                    }                
                    currentSectionAngle -= marginAngle0 + marginAngle1;
                }


                // if (currentSectionRadialHeight > availableRadialHeight) {
                //     return false;
                // }
                if (currentSectionAngle < 0) {
                    return false;
                }

                partialWinding.set_parallels_proportion(sectionParallelsProportion);
                section.set_name(get_name(windingIndex) +  " section " + std::to_string(currentSectionPerWinding[windingIndex]));
                section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
                section.set_type(ElectricalType::CONDUCTION);
                section.set_group(group.get_name());
                section.set_margin(_marginsPerSection[sectionIndex]);
                section.set_layers_orientation(_layersOrientation);
                section.set_coordinate_system(CoordinateSystem::POLAR);
                
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{currentSectionRadialHeight, currentSectionAngle});
                    section.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + currentSectionRadialHeight / 2 + _marginsPerSection[sectionIndex][0], currentSectionCenterAngle, 0});
                }
                else {
                    section.set_dimensions(std::vector<double>{currentSectionRadialHeight, currentSectionAngle});
                    section.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + currentSectionRadialHeight / 2, currentSectionCenterAngle  + currentSectionAngle / 2 + marginAngle0, 0});
                }

                if (section.get_dimensions()[0] < 0) {
                    throw std::runtime_error("Something wrong happened in section dimensions 0: " + std::to_string(section.get_dimensions()[0]) +
                                             " currentSectionRadialHeight: " + std::to_string(currentSectionRadialHeight) +
                                             " currentSectionAngle: " + std::to_string(currentSectionAngle)
                                             );
                }


                if (section.get_dimensions()[1] < 0) {
                    throw std::runtime_error("Something wrong happened in section dimensions 1: " + std::to_string(section.get_dimensions()[1]) +
                                             " currentSectionRadialHeight: " + std::to_string(currentSectionRadialHeight) +
                                             " currentSectionAngle: " + std::to_string(currentSectionAngle)
                                             );
                }
                
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    double ringArea = std::numbers::pi * pow(availableRadialHeight - currentSectionCenterRadialHeight, 2) - std::numbers::pi * pow(availableRadialHeight - (currentSectionCenterRadialHeight + currentSectionRadialHeight), 2);

                    section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / ringArea);
                } 
                else {
                    double ringArea = std::numbers::pi * pow(availableRadialHeight - currentSectionCenterRadialHeight, 2) * currentSectionAngle / 360;
                    section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / ringArea);
                }
                section.set_winding_style(windByConsecutiveTurns[windingIndex]);
                sectionsDescription.push_back(section);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportion[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
                }

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterRadialHeight += currentSectionRadialHeight + _marginsPerSection[sectionIndex][0] + _marginsPerSection[sectionIndex][1];
                }
                else {
                    currentSectionCenterAngle += currentSectionAngle + marginAngle0 + marginAngle1;
                }

                currentSectionPerWinding[windingIndex]++;
                conductingSectionIndex++;
            }
            else {
                if (sectionIndex == 0) {
                    throw std::runtime_error("Insulation layers cannot be the first one (for now)");
                }
     
                auto previousWindingIndex = orderedSectionsWithInsulation[sectionIndex - 1].second.first;
                size_t nextWindingIndex;
                if (sectionIndex != (orderedSectionsWithInsulation.size() - 1)) {
                    nextWindingIndex = orderedSectionsWithInsulation[sectionIndex + 1].second.first;
                }
                else {
                    nextWindingIndex = orderedSectionsWithInsulation[0].second.first;
                }

                auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, nextWindingIndex};
                if (!_insulationSections.contains(windingsMapKey)) {
                    log(_insulationSectionsLog[windingsMapKey]);
                    continue;
                }

                auto insulationSection = _insulationSections[windingsMapKey];

                insulationSection.set_group(group.get_name());
                insulationSection.set_name("Insulation between " + get_name(previousWindingIndex) + " and " + get_name(previousWindingIndex) + " section " + std::to_string(sectionIndex));
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + insulationSection.get_dimensions()[0] / 2,
                                                                          currentSectionCenterAngle,
                                                                          0});
                }
                else {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + insulationSection.get_dimensions()[0] / 2,
                                                                          currentSectionCenterAngle - insulationSection.get_dimensions()[1] / 2,
                                                                          0});
                }

                sectionsDescription.push_back(insulationSection);
                log(_insulationSectionsLog[windingsMapKey]);

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterRadialHeight += insulationSection.get_dimensions()[0];
                }
                else {
                    currentSectionCenterAngle += insulationSection.get_dimensions()[1];
                }
            }
        }
    }


    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            if (roundFloat(remainingParallelsProportion[windingIndex][parallelIndex], 9) != 0) {
                throw std::runtime_error("There are unassigned parallel proportion in a round section, something went wrong");
            }
        }
    }

    set_sections_description(sectionsDescription);
    return true;
}

bool Coil::wind_by_planar_sections(std::vector<size_t> stackUpForThisGroup, std::optional<double> insulationThickness, double coreToLayerDistance) {
    // In planar coils each section will have only one layer
    set_layers_description(std::nullopt);
    std::vector<Section> sections;

    if (!insulationThickness) {
        insulationThickness = defaults.pcbInsulationThickness;
    }

    auto bobbin = resolve_bobbin();
    if (!get_groups_description()) {
        create_default_group(bobbin, WiringTechnology::PRINTED);
    }

    if (!get_groups_description()) {
        throw std::runtime_error("At least default group must be defined at this point.");
    }

    auto groups = get_groups_description().value();
    if (groups.size() > 1) {
        throw std::runtime_error("Only one group supported for now.");
    }
    auto group = groups[0];

    auto wirePerWinding = get_wires();
    if (wirePerWinding.size() == 0) {
        throw std::runtime_error("Wires missing");
    }

    std::vector<size_t> numberSectionsPerWinding = std::vector<size_t>(wirePerWinding.size(), 0);
    std::vector<std::vector<double>> totalParallelsProportionPerWinding;
    std::vector<std::vector<double>> remainingParallelsProportionPerWinding;
    std::vector<size_t> currentSectionIndexPerwinding = std::vector<size_t>(wirePerWinding.size(), 0);

    for (auto windingIndex : stackUpForThisGroup) {
        numberSectionsPerWinding[windingIndex]++;
    }

    for (auto winding : group.get_partial_windings()) {
        totalParallelsProportionPerWinding.push_back(winding.get_parallels_proportion());
        remainingParallelsProportionPerWinding.push_back(winding.get_parallels_proportion());
    }
    for (auto partialWinding : group.get_partial_windings()) {
        auto parallelsProportion = partialWinding.get_parallels_proportion();
        totalParallelsProportionPerWinding.push_back(parallelsProportion);
        remainingParallelsProportionPerWinding.push_back(parallelsProportion);
    }

    std::vector<double> sectionWidthPerWinding;
    std::vector<double> sectionHeightPerWinding;
    double totalSectionHeight = 0;

    for (size_t stackUpIndex = 0; stackUpIndex < stackUpForThisGroup.size(); ++stackUpIndex) {
        auto windingIndex = stackUpForThisGroup[stackUpIndex];
        sectionWidthPerWinding.push_back(group.get_dimensions()[0]);
        double sectionHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
        sectionHeightPerWinding.push_back(sectionHeight);
        totalSectionHeight += sectionHeight;
        if (stackUpIndex < stackUpForThisGroup.size() - 1) {
            totalSectionHeight += insulationThickness.value();
        }
    }
    double currentSectionCenterWidth = roundFloat(group.get_coordinates()[0], 9);
    double currentSectionCenterHeight = roundFloat(group.get_coordinates()[1] + totalSectionHeight / 2, 9);

    for (size_t stackUpIndex = 0; stackUpIndex < stackUpForThisGroup.size(); ++stackUpIndex) {
        Section section;
        auto windingIndex = stackUpForThisGroup[stackUpIndex];
        auto remainingParallelsProportionInWinding = remainingParallelsProportionPerWinding[windingIndex];
        auto totalParallelsProportionInWinding = totalParallelsProportionPerWinding[windingIndex];
        auto numberSections = numberSectionsPerWinding[windingIndex];
        auto winding = get_functional_description()[windingIndex];
        auto sectionIndex = currentSectionIndexPerwinding[windingIndex];
        double sectionWidth = sectionWidthPerWinding[windingIndex] - coreToLayerDistance * 2;
        double sectionHeight = sectionHeightPerWinding[windingIndex];
        currentSectionCenterHeight -= sectionHeight / 2;

        WindingStyle windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberSections);

        double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
 
        auto parallelsProportions = get_parallels_proportions(sectionIndex,
                                                               numberSections,
                                                               get_number_turns(windingIndex),
                                                               get_number_parallels(windingIndex),
                                                               remainingParallelsProportionInWinding,
                                                               windByConsecutiveTurns,
                                                               totalParallelsProportionInWinding);

        std::vector<double> sectionParallelsProportion = parallelsProportions.second;

        size_t numberParallelsProportionsToZero = 0;
        for (auto parallelProportion : sectionParallelsProportion) {
            if (parallelProportion == 0) {
                numberParallelsProportionsToZero++;
            }
        }

        if (numberParallelsProportionsToZero == sectionParallelsProportion.size()) {
            throw std::runtime_error("Parallel proportion in section cannot be all be 0");
        }

        uint64_t physicalTurnsThisSection = parallelsProportions.first;

        auto partialWinding = group.get_partial_windings()[windingIndex];
        partialWinding.set_parallels_proportion(sectionParallelsProportion);
        section.set_partial_windings(std::vector<PartialWinding>{partialWinding});
        section.set_group(group.get_name());
        section.set_type(ElectricalType::CONDUCTION);
        section.set_name(winding.get_name() + " section " + std::to_string(sectionIndex));
        section.set_layers_orientation(WindingOrientation::CONTIGUOUS);
        section.set_dimensions(std::vector<double>{sectionWidth, sectionHeight});
        section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
        section.set_coordinate_system(CoordinateSystem::CARTESIAN);

        section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / (sectionWidth * sectionHeight));
        section.set_winding_style(windByConsecutiveTurns);
        sections.push_back(section);

        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            remainingParallelsProportionPerWinding[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
        }

        currentSectionCenterHeight = roundFloat(currentSectionCenterHeight -= sectionHeight / 2, 9);
        currentSectionIndexPerwinding[windingIndex]++;

        if (stackUpIndex < stackUpForThisGroup.size() - 1 && insulationThickness.value() > 0) {
            currentSectionCenterHeight -= insulationThickness.value() / 2;

            Section insulationSection;
            insulationSection.set_type(ElectricalType::INSULATION);
            insulationSection.set_name("Insulation section between stack index" + std::to_string(stackUpIndex) + " and " + std::to_string(stackUpIndex + 1));
            insulationSection.set_dimensions(std::vector<double>{sectionWidth, insulationThickness.value()});
            insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            insulationSection.set_coordinate_system(CoordinateSystem::CARTESIAN);
            insulationSection.set_layers_orientation(WindingOrientation::CONTIGUOUS);
            insulationSection.set_filling_factor(1);
            sections.push_back(insulationSection);
            currentSectionCenterHeight -= insulationThickness.value() / 2;
        }

    }
    set_sections_description(sections);
    return true;
}

bool Coil::wind_by_layers() {
    set_layers_description(std::nullopt);
    if (!get_sections_description()) {
        return false;
    }
    auto bobbin = resolve_bobbin();

    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        return wind_by_rectangular_layers();
    }
    else {
        return wind_by_round_layers();
    }
}

bool Coil::wind_by_rectangular_layers() {
    set_layers_description(std::nullopt);
    if (!get_sections_description()) {
        return false;
    }

    auto wirePerWinding = get_wires();

    auto sections = get_sections_description().value();

    std::vector<Layer> layers;
    for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
        if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {

            uint64_t maximumNumberLayersFittingInSection;
            uint64_t maximumNumberPhysicalTurnsPerLayer;
            uint64_t minimumNumberLayerNeeded;
            uint64_t numberLayers;
            uint64_t physicalTurnsInSection = 0;
            double layerWidth = 0;
            double layerHeight = 0;
            auto remainingParallelsProportionInSection = sections[sectionIndex].get_partial_windings()[0].get_parallels_proportion();
            auto totalParallelsProportionInSection = sections[sectionIndex].get_partial_windings()[0].get_parallels_proportion();
            if (sections[sectionIndex].get_partial_windings().size() > 1) {
                throw std::runtime_error("More than one winding per layer not supported yet");
            }
            auto partialWinding = sections[sectionIndex].get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                physicalTurnsInSection += round(remainingParallelsProportionInSection[parallelIndex] * get_number_turns(windingIndex));
            }

            if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                if (!wirePerWinding[windingIndex].get_outer_diameter()) {
                    throw std::runtime_error("Missing wire outer diameter");
                }
                double wireDiameter = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_diameter().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[0] / wireDiameter;
                    maximumNumberPhysicalTurnsPerLayer = floor(sections[sectionIndex].get_dimensions()[1] / wireDiameter);
                    layerWidth = wireDiameter;
                    layerHeight = sections[sectionIndex].get_dimensions()[1];
                } else {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[1] / wireDiameter;
                    maximumNumberPhysicalTurnsPerLayer = floor(sections[sectionIndex].get_dimensions()[0] / wireDiameter);
                    layerWidth = sections[sectionIndex].get_dimensions()[0];
                    layerHeight = wireDiameter;
                }
            }
            else {
                if (!wirePerWinding[windingIndex].get_outer_width()) {
                    throw std::runtime_error("Missing wire outer width");
                }
                if (!wirePerWinding[windingIndex].get_outer_height()) {
                    throw std::runtime_error("Missing wire outer height");
                }
                double wireWidth = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_width().value());
                double wireHeight = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_height().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[0] / wireWidth;
                    if (wirePerWinding[windingIndex].get_type() == WireType::FOIL) {
                        maximumNumberPhysicalTurnsPerLayer = 1;
                    }
                    else {
                        maximumNumberPhysicalTurnsPerLayer = floor(sections[sectionIndex].get_dimensions()[1] / wireHeight);
                    }
                    layerWidth = wireWidth;
                    layerHeight = sections[sectionIndex].get_dimensions()[1];
                } else {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[1] / wireHeight;
                    if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR && settings->get_coil_only_one_turn_per_layer_in_contiguous_rectangular()) {
                        maximumNumberPhysicalTurnsPerLayer = 1; 
                    }
                    else {
                        maximumNumberPhysicalTurnsPerLayer = floor(sections[sectionIndex].get_dimensions()[0] / wireWidth);
                    }
                    layerWidth = sections[sectionIndex].get_dimensions()[0];
                    layerHeight = wireHeight;
                }
            }


            if (maximumNumberLayersFittingInSection == 0) {
                numberLayers = ceil(double(physicalTurnsInSection) / maximumNumberPhysicalTurnsPerLayer);
            }
            else if (maximumNumberPhysicalTurnsPerLayer == 0) {
                numberLayers = maximumNumberLayersFittingInSection;
            }
            else {
                minimumNumberLayerNeeded = ceil(double(physicalTurnsInSection) / maximumNumberPhysicalTurnsPerLayer);
                numberLayers = std::min(minimumNumberLayerNeeded, maximumNumberLayersFittingInSection);
            }

            // We cannot have more layers than physical turns
            numberLayers = std::min(numberLayers, physicalTurnsInSection);
            auto turnsAlignment = get_turns_alignment(sections[sectionIndex].get_name());

            double currentLayerCenterWidth;
            double currentLayerCenterHeight;
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerWidth / 2, 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0], 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);

                if (turnsAlignment == CoilAlignment::CENTERED || turnsAlignment == CoilAlignment::SPREAD) {
                    currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + (numberLayers * layerHeight) / 2 - layerHeight / 2, 9);
                }
                else if (turnsAlignment == CoilAlignment::INNER_OR_TOP) {
                    currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);
                }
                else if (turnsAlignment == CoilAlignment::OUTER_OR_BOTTOM) {
                    currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] - sections[sectionIndex].get_dimensions()[1] / 2 + (numberLayers * layerHeight) - layerHeight / 2, 9);
                }
                else {
                    throw std::invalid_argument("Unknown turns alignment");
                }
            }

            WindingStyle windByConsecutiveTurns;
            if (sections[sectionIndex].get_winding_style()) {
                windByConsecutiveTurns = sections[sectionIndex].get_winding_style().value();
            }
            else {
                windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberLayers);
            }

            if (sections[sectionIndex].get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS && maximumNumberPhysicalTurnsPerLayer < get_number_parallels(windingIndex)) {
                windByConsecutiveTurns = WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
            }

            for (size_t layerIndex = 0; layerIndex < numberLayers; ++layerIndex) {
                Layer layer;

                auto parallelsProportions = get_parallels_proportions(layerIndex,
                                                                       numberLayers,
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportionInSection,
                                                                       windByConsecutiveTurns,
                                                                       totalParallelsProportionInSection);

                std::vector<double> layerParallelsProportion = parallelsProportions.second;

                size_t numberParallelsProportionsToZero = 0;
                for (auto parallelProportion : layerParallelsProportion) {
                    if (parallelProportion == 0) {
                        numberParallelsProportionsToZero++;
                    }
                }

                if (numberParallelsProportionsToZero == layerParallelsProportion.size()) {
                    throw std::runtime_error("Parallel proportion in layer cannot be all be 0");
                }

                uint64_t physicalTurnsThisLayer = parallelsProportions.first;

                partialWinding.set_parallels_proportion(layerParallelsProportion);
                layer.set_partial_windings(std::vector<PartialWinding>{partialWinding});
                layer.set_section(sections[sectionIndex].get_name());
                layer.set_type(ElectricalType::CONDUCTION);
                layer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                layer.set_orientation(sections[sectionIndex].get_layers_orientation());
                layer.set_turns_alignment(turnsAlignment);
                layer.set_dimensions(std::vector<double>{layerWidth, layerHeight});
                layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layer.set_coordinate_system(CoordinateSystem::CARTESIAN);

                layer.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisLayer) / (layerWidth * layerHeight));
                layer.set_winding_style(windByConsecutiveTurns);
                layers.push_back(layer);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportionInSection[parallelIndex] -= layerParallelsProportion[parallelIndex];
                }

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                    currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight, 9);
                }
                else {
                    currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth, 9);
                }
            }
        }
        else {
            if (sectionIndex == 0) {
                throw std::runtime_error("outer insulation layers not implemented");
            }

            auto partialWinding = sections[sectionIndex - 1].get_partial_windings()[0];
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            Section nextSection;
            if (sectionIndex != (sections.size() - 1)) {
                if (sections[sectionIndex - 1].get_type() != ElectricalType::CONDUCTION || sections[sectionIndex + 1].get_type() != ElectricalType::CONDUCTION) {
                    throw std::runtime_error("Previous and next sections must be conductive");
                }
                nextSection = sections[sectionIndex + 1];
            }
            else {
                nextSection = sections[0];
            }
            // auto nextSection = sections[sectionIndex + 1];
            auto nextPartialWinding = nextSection.get_partial_windings()[0];
            auto nextWindingIndex = get_winding_index_by_name(nextPartialWinding.get_winding());

            auto windingsMapKey = std::pair<size_t, size_t>{windingIndex, nextWindingIndex};
            if (!_insulationLayers.contains(windingsMapKey)) {
                log(_insulationLayersLog[windingsMapKey]);
                continue;
            }

            auto insulationLayers = _insulationLayers[windingsMapKey];
            if (insulationLayers.size() == 0) {
                throw std::runtime_error("There must be at least one insulation layer between layers");
            }

            double layerWidth = insulationLayers[0].get_dimensions()[0];
            double layerHeight = insulationLayers[0].get_dimensions()[1];

            double currentLayerCenterWidth;
            double currentLayerCenterHeight;
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerWidth / 2, 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0], 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);
            }

            for (size_t layerIndex = 0; layerIndex < insulationLayers.size(); ++layerIndex) {
                auto insulationLayer = insulationLayers[layerIndex];
                insulationLayer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                insulationLayer.set_section(sections[sectionIndex].get_name());
                insulationLayer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layers.push_back(insulationLayer);

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                    currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight, 9);
                }
                else {
                    currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth, 9);
                }
            }
        }
    }
    set_layers_description(layers);
    return true;
}

bool Coil::wind_by_round_layers() {
    set_layers_description(std::nullopt);
    if (!get_sections_description()) {
        return false;
    }
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    auto wirePerWinding = get_wires();

    auto sections = get_sections_description().value();

    std::vector<Layer> layers;
    for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
        if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
            uint64_t maximumNumberLayersFittingInSection;
            uint64_t maximumNumberPhysicalTurnsPerLayer;
            uint64_t minimumNumberLayerNeeded = 0;
            uint64_t numberLayers;
            std::vector<int64_t> layerPhysicalTurns;
            uint64_t physicalTurnsInSection = 0;
            double layerRadialHeight = 0;
            double layerAngle = 0;
            auto remainingParallelsProportionInSection = sections[sectionIndex].get_partial_windings()[0].get_parallels_proportion();
            auto totalParallelsProportionInSection = sections[sectionIndex].get_partial_windings()[0].get_parallels_proportion();
            if (sections[sectionIndex].get_partial_windings().size() > 1) {
                throw std::runtime_error("More than one winding per layer not supported yet");
            }
            auto partialWinding = sections[sectionIndex].get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                physicalTurnsInSection += round(remainingParallelsProportionInSection[parallelIndex] * get_number_turns(windingIndex));
            }

            if (std::isnan(sections[sectionIndex].get_coordinates()[1]) ||sections[sectionIndex].get_coordinates()[1] < 0) {
                return false;
                // throw std::runtime_error("sections[sectionIndex].get_coordinates()[1] cannot be negative: " + std::to_string(sections[sectionIndex].get_coordinates()[1]));
            }


            if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                double wireDiameter = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_diameter().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = roundFloat(sections[sectionIndex].get_dimensions()[0] / wireDiameter, 9);
                    double averageLayerPerimeter = 2 * std::numbers::pi * (sections[sectionIndex].get_dimensions()[1] / 360) * (windingWindowRadialHeight - sections[sectionIndex].get_coordinates()[0]);
                    maximumNumberPhysicalTurnsPerLayer = floor(averageLayerPerimeter / wireDiameter);
                    layerRadialHeight = wireDiameter;
                    layerAngle = sections[sectionIndex].get_dimensions()[1];
                } else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }
            }
            else {
                double wireWidth = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_width().value());
                double wireHeight = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_height().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = roundFloat(sections[sectionIndex].get_dimensions()[0] / wireWidth, 9);
                    double averageLayerPerimeter = 2 * std::numbers::pi * (sections[sectionIndex].get_dimensions()[1] / 360) * (windingWindowRadialHeight - sections[sectionIndex].get_coordinates()[0]);
                    if (wirePerWinding[windingIndex].get_type() == WireType::FOIL) {
                        throw std::invalid_argument("Cannot have foil in toroids");
                    }
                    else {
                        maximumNumberPhysicalTurnsPerLayer = floor(averageLayerPerimeter / wireHeight);
                    }
                    layerRadialHeight = wireWidth;
                    layerAngle = sections[sectionIndex].get_dimensions()[1];
                } else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }
            }

            if (maximumNumberLayersFittingInSection == 0) {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight);
                numberLayers = aux.first;
                layerPhysicalTurns = aux.second;
            }
            else if (maximumNumberPhysicalTurnsPerLayer == 0) {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight);
                layerPhysicalTurns = aux.second;
                numberLayers = maximumNumberLayersFittingInSection;
            }
            else {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight);
                minimumNumberLayerNeeded = aux.first;
                layerPhysicalTurns = aux.second;
                numberLayers = std::min(minimumNumberLayerNeeded, maximumNumberLayersFittingInSection);
            }

            // We cannot have more layers than physical turns
            numberLayers = std::min(numberLayers, physicalTurnsInSection);

            if (minimumNumberLayerNeeded > numberLayers) {
                return false;
            }

            double currentLayerCenterRadialHeight;
            double currentLayerCenterAngle;
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                currentLayerCenterRadialHeight = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerRadialHeight / 2, 9);
                currentLayerCenterAngle = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                throw std::invalid_argument("Only overlapping layers allowed in toroids");
            }

            WindingStyle windByConsecutiveTurns;
            if (sections[sectionIndex].get_winding_style()) {
                windByConsecutiveTurns = sections[sectionIndex].get_winding_style().value();
            }
            else {
                windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberLayers);
            }

            if (sections[sectionIndex].get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS && maximumNumberPhysicalTurnsPerLayer < get_number_parallels(windingIndex)) {
                windByConsecutiveTurns = WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
            }

            for (size_t layerIndex = 0; layerIndex < numberLayers; ++layerIndex) {
                Layer layer;

                // TODO: probably will have to add a factor to have more or less turns per layer than the section average

                auto parallelsProportions = get_parallels_proportions(layerIndex,
                                                                       numberLayers,
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportionInSection,
                                                                       windByConsecutiveTurns,
                                                                       totalParallelsProportionInSection,
                                                                       1,
                                                                       layerPhysicalTurns[layerIndex]);

                std::vector<double> layerParallelsProportion = parallelsProportions.second;

                size_t numberParallelsProportionsToZero = 0;
                for (auto parallelProportion : layerParallelsProportion) {
                    if (parallelProportion == 0) {
                        numberParallelsProportionsToZero++;
                    }
                }

                if (numberParallelsProportionsToZero == layerParallelsProportion.size()) {
                    throw std::runtime_error("Parallel proportion in layer cannot be all be 0");
                }

                uint64_t physicalTurnsThisLayer = parallelsProportions.first;
                auto turnsAlignment = get_turns_alignment(sections[sectionIndex].get_name());

                partialWinding.set_parallels_proportion(layerParallelsProportion);
                layer.set_partial_windings(std::vector<PartialWinding>{partialWinding});
                layer.set_section(sections[sectionIndex].get_name());
                layer.set_type(ElectricalType::CONDUCTION);
                layer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                layer.set_orientation(sections[sectionIndex].get_layers_orientation());
                layer.set_turns_alignment(turnsAlignment);
                layer.set_dimensions(std::vector<double>{layerRadialHeight, layerAngle});
                layer.set_coordinates(std::vector<double>{currentLayerCenterRadialHeight, currentLayerCenterAngle, 0});
                layer.set_coordinate_system(CoordinateSystem::POLAR);

                double layerPerimeter = 2 * std::numbers::pi * (layerAngle / 360) * (windingWindowRadialHeight - layerRadialHeight);
                layer.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisLayer) / (layerPerimeter * layerRadialHeight));
                layer.set_winding_style(windByConsecutiveTurns);
                layers.push_back(layer);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportionInSection[parallelIndex] -= layerParallelsProportion[parallelIndex];
                }

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    currentLayerCenterRadialHeight = roundFloat(currentLayerCenterRadialHeight + layerRadialHeight, 9);
                }
                else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }
            }
        }
        else {
            if (sectionIndex == 0) {
                throw std::runtime_error("Inner insulation layers not implemented");
            }

            auto partialWinding = sections[sectionIndex - 1].get_partial_windings()[0];
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

            Section nextSection;
            if (sectionIndex != (sections.size() - 1)) {
                if (sections[sectionIndex - 1].get_type() != ElectricalType::CONDUCTION || sections[sectionIndex + 1].get_type() != ElectricalType::CONDUCTION) {
                    throw std::runtime_error("Previous and next sections must be conductive");
                }
                nextSection = sections[sectionIndex + 1];
            }
            else {
                nextSection = sections[0];
            }
            // auto nextSection = sections[sectionIndex + 1];
            auto nextPartialWinding = nextSection.get_partial_windings()[0];
            auto nextWindingIndex = get_winding_index_by_name(nextPartialWinding.get_winding());

            // If the angle of the section is 0 it means that the margin is enoguh and we don't need to add insulation layers
            if (sections[sectionIndex].get_dimensions()[1] == 0) {
                continue;
            }

            auto windingsMapKey = std::pair<size_t, size_t>{windingIndex, nextWindingIndex};
            if (!_insulationLayers.contains(windingsMapKey)) {
                log(_insulationLayersLog[windingsMapKey]);
                continue;
            }

            auto insulationLayers = _insulationLayers[windingsMapKey];
            if (insulationLayers.size() == 0) {
                throw std::runtime_error("There must be at least one insulation layer between layers");
            }

            double layerRadialHeight = insulationLayers[0].get_dimensions()[0];

            double currentLayerCenterRadialHeight;
            double currentLayerCenterAngle;

            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                currentLayerCenterRadialHeight = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerRadialHeight / 2, 9);
                currentLayerCenterAngle = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                throw std::invalid_argument("Only overlapping layers allowed in toroids");
            }

            for (size_t layerIndex = 0; layerIndex < insulationLayers.size(); ++layerIndex) {
                auto insulationLayer = insulationLayers[layerIndex];
                insulationLayer.set_section(sections[sectionIndex].get_name());
                insulationLayer.set_coordinate_system(CoordinateSystem::POLAR);
                insulationLayer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterRadialHeight, currentLayerCenterAngle, 0});
                layers.push_back(insulationLayer);


                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    currentLayerCenterRadialHeight = roundFloat(currentLayerCenterRadialHeight + layerRadialHeight, 9);
                }
                else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }
            }
        }
    }
    set_layers_description(layers);
    return true;
}

bool Coil::wind_by_planar_layers() {
    set_layers_description(std::nullopt);
    std::vector<Layer> layers;
    if (!get_sections_description()) {
        return false;
    }

    auto sections = get_sections_description().value();

    for (auto section : sections) {
        Layer layer;
        layer.set_partial_windings(section.get_partial_windings());
        layer.set_section(section.get_name());
        layer.set_type(section.get_type());
        layer.set_orientation(section.get_layers_orientation());
        layer.set_dimensions(section.get_dimensions());
        layer.set_coordinates(section.get_coordinates());
        layer.set_coordinate_system(section.get_coordinate_system());
        layer.set_winding_style(section.get_winding_style());
        layer.set_filling_factor(section.get_filling_factor());
        layer.set_name(std::regex_replace(std::string(section.get_name()), std::regex("section"), "layer"));
        if (section.get_type() == ElectricalType::CONDUCTION) {
            layer.set_turns_alignment(CoilAlignment::SPREAD);
        }
        else {
            layer.set_insulation_material(defaults.defaultPcbInsulationMaterial);
        }

        layers.push_back(layer);

    }
    set_layers_description(layers);
    return true;
}

bool Coil::wind_by_turns() {
    set_turns_description(std::nullopt);
    if (!get_layers_description()) {
        return false;
    }
    auto bobbin = resolve_bobbin();

    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        return wind_by_rectangular_turns();
    }
    else {
        return wind_by_round_turns();
    }
}

bool Coil::wind_by_rectangular_turns() {
    set_turns_description(std::nullopt);
    if (!get_layers_description()) {
        return false;
    }
    auto wirePerWinding = get_wires();
    std::vector<std::vector<int64_t>> currentTurnIndex;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        currentTurnIndex.push_back(std::vector<int64_t>(get_number_parallels(windingIndex), 0));
    }
    auto bobbin = resolve_bobbin();
    auto bobbinColumnShape = bobbin.get_processed_description().value().get_column_shape();
    auto bobbinColumnDepth = bobbin.get_processed_description().value().get_column_depth();
    double bobbinColumnWidth;
    if (bobbin.get_processed_description().value().get_column_width()) {
        bobbinColumnWidth = bobbin.get_processed_description().value().get_column_width().value();
    }
    else {
        auto bobbinWindingWindow = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows()[0];
        double bobbinWindingWindowWidth = bobbinWindingWindow.get_width().value();
        double bobbinWindingWindowCenterWidth = bobbinWindingWindow.get_coordinates().value()[0];
        bobbinColumnWidth = bobbinWindingWindowCenterWidth - bobbinWindingWindowWidth / 2;
    }

    auto layers = get_layers_description().value();

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        if (wirePerWinding[windingIndex].get_type() == WireType::PLANAR) {
            auto conductionLayers = get_layers_by_type(ElectricalType::CONDUCTION);
            if (conductionLayers.size() > settings->get_coil_maximum_layers_planar()) {
                return false;
            }
        }

    }

    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            double currentTurnCenterWidth = 0;
            double currentTurnCenterHeight = 0;
            double currentTurnWidthIncrement = 0;
            double currentTurnHeightIncrement = 0;
            double totalLayerHeight;
            double totalLayerWidth;
            if (layer.get_partial_windings().size() > 1) {
                throw std::runtime_error("More than one winding per layer not supported yet");
            }
            auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
            double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
            auto physicalTurnsInLayer = get_number_turns(layer);
            auto alignment = layer.get_turns_alignment().value();

            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                totalLayerWidth = layer.get_dimensions()[0];
                totalLayerHeight = roundFloat(physicalTurnsInLayer * wireHeight, 9);

                currentTurnWidthIncrement = 0;
                currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0], 9);
                switch (alignment) {
                    case CoilAlignment::CENTERED:
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + totalLayerHeight / 2 - wireHeight / 2, 9);
                        currentTurnHeightIncrement = wireHeight;
                        break;

                    case CoilAlignment::INNER_OR_TOP:
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - wireHeight / 2, 9);
                        currentTurnHeightIncrement = wireHeight;
                        break;

                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 + totalLayerHeight - wireHeight / 2, 9);
                        currentTurnHeightIncrement = wireHeight;
                        break;

                    case CoilAlignment::SPREAD:
                        currentTurnHeightIncrement = roundFloat(layer.get_dimensions()[1] / physicalTurnsInLayer, 9);
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - currentTurnHeightIncrement / 2, 9);
                        break;
                }

            } 
            else {
                totalLayerWidth = roundFloat(physicalTurnsInLayer * wireWidth, 9);
                if (totalLayerWidth > std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows()[0].get_width().value()) {
                    return false;
                }
                totalLayerHeight = layer.get_dimensions()[1];
                currentTurnHeightIncrement = 0;
                currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1], 9);
                switch (alignment) {
                    case CoilAlignment::CENTERED:
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - totalLayerWidth / 2 + wireWidth / 2, 9);
                        currentTurnWidthIncrement = wireWidth;
                        break;

                    case CoilAlignment::INNER_OR_TOP:
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + wireWidth / 2, 9);
                        currentTurnWidthIncrement = wireWidth;
                        break;

                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + (layer.get_dimensions()[0] - totalLayerWidth) + wireWidth / 2, 9);
                        currentTurnWidthIncrement = wireWidth;
                        break;

                    case CoilAlignment::SPREAD:
                        currentTurnWidthIncrement = roundFloat(layer.get_dimensions()[0] / physicalTurnsInLayer, 9);
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + wireWidth / 2 + currentTurnWidthIncrement / 2, 9);
                        break;
                }
            }

            if (!layer.get_winding_style()) {
                layer.set_winding_style(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            }


            if (layer.get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    int64_t numberTurns = round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
                    for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                        Turn turn;
                        turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                        turn.set_layer(layer.get_name());
                        if (bobbinColumnShape == ColumnShape::ROUND) {
                            turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth);
                            if (turn.get_length() < 0) {
                                return false;
                                // throw std::runtime_error("Something wrong happened in turn length 2: " + std::to_string(turn.get_length()) +
                                //                          " currentTurnCenterWidth: " + std::to_string(currentTurnCenterWidth) +
                                //                          " layer.get_coordinates()[0]: " + std::to_string(layer.get_coordinates()[0]));
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::OBLONG) {
                            turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                            double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                            turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);

                            if (turn.get_length() < 0) {
                                return false;
                                // throw std::runtime_error("Something wrong happened in turn length 2: " + std::to_string(turn.get_length()) + " bobbinColumnDepth: " + std::to_string(bobbinColumnDepth)  + " bobbinColumnWidth: " + std::to_string(bobbinColumnWidth)  + " currentTurnCornerRadius: " + std::to_string(currentTurnCornerRadius));
                            }
                        }
                        else {
                            throw std::runtime_error("only round or rectangular columns supported for bobbins");
                        }
                        turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                        turn.set_orientation(TurnOrientation::CLOCKWISE);
                        turn.set_parallel(parallelIndex);
                        turn.set_section(layer.get_section().value());
                        turn.set_winding(partialWinding.get_winding());
                        turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                        turn.set_rotation(0);
                        turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

                        turns.push_back(turn);
                        currentTurnCenterWidth += currentTurnWidthIncrement;
                        currentTurnCenterHeight -= currentTurnHeightIncrement;
                        currentTurnIndex[windingIndex][parallelIndex]++; 
                    }
                }
            }
            else {
                int64_t firstParallelIndex = 0;
                while (roundFloat(partialWinding.get_parallels_proportion()[firstParallelIndex], 10) == 0) {
                    firstParallelIndex++;
                }
                int64_t numberTurns = round(partialWinding.get_parallels_proportion()[firstParallelIndex] * get_number_turns(windingIndex));
                for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                        if (roundFloat(partialWinding.get_parallels_proportion()[parallelIndex], 10) > 0) {
                            Turn turn;
                            turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                            turn.set_layer(layer.get_name());
                            if (bobbinColumnShape == ColumnShape::ROUND) {
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth);
                                    if (turn.get_length() < 0) {
                                        return false;
                                        // throw std::runtime_error("Something wrong happened in turn length 3: " + std::to_string(turn.get_length()) + " currentTurnCenterWidth: " + std::to_string(currentTurnCenterWidth));
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::OBLONG) {
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                                    if (turn.get_length() < 0) {
                                        return false;
                                        // throw std::runtime_error("Something wrong happened in turn length 3: " + std::to_string(turn.get_length()) + " currentTurnCenterWidth: " + std::to_string(currentTurnCenterWidth));
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                                turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                                if (turn.get_length() < 0) {
                                    return false;
                                    // throw std::runtime_error("Something wrong happened in turn length 3: " + std::to_string(turn.get_length()) + " bobbinColumnDepth: " + std::to_string(bobbinColumnDepth)  + " bobbinColumnWidth: " + std::to_string(bobbinColumnWidth)  + " currentTurnCornerRadius: " + std::to_string(currentTurnCornerRadius));
                                }
                            }
                            else {
                                throw std::runtime_error("only round or rectangular columns supported for bobbins");
                            }
                            turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                            turn.set_orientation(TurnOrientation::CLOCKWISE);
                            turn.set_parallel(parallelIndex);
                            turn.set_section(layer.get_section().value());
                            turn.set_winding(partialWinding.get_winding());
                            turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                            turn.set_rotation(0);
                            turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

                            turns.push_back(turn);
                            currentTurnCenterWidth += currentTurnWidthIncrement;
                            currentTurnCenterHeight -= currentTurnHeightIncrement;
                            currentTurnIndex[windingIndex][parallelIndex]++; 
                        }
                    }
                }
            }
        }
    }

    set_turns_description(turns);
    return true;
}

bool Coil::wind_by_round_turns() {
    set_turns_description(std::nullopt);
    if (!get_layers_description()) {
        return false;
    }
    auto wirePerWinding = get_wires();
    std::vector<std::vector<int64_t>> currentTurnIndex;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        currentTurnIndex.push_back(std::vector<int64_t>(get_number_parallels(windingIndex), 0));
    }
    auto bobbinColumnShape = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_shape();
    auto bobbinColumnDepth = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_depth();
    double bobbinColumnWidth;
    if (std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_width()) {
        bobbinColumnWidth = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_width().value();
    }
    else {
        throw std::runtime_error("Toroids must have their bobbin column set");
    }

    auto layers = get_layers_description().value();

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR) {
            auto layersInWinding = get_layers_by_winding_index(windingIndex);
            if (layersInWinding.size() > 1) {
                return false;
            }
        }
    }

    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            double currentTurnCenterRadialHeight = 0;
            double currentTurnCenterAngle = 0;
            double currentTurnRadialHeightIncrement = 0;
            double currentTurnAngleIncrement = 0;
            double totalLayerAngle;
            if (layer.get_partial_windings().size() > 1) {
                throw std::runtime_error("More than one winding per layer not supported yet");
            }
            auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
            double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
            auto physicalTurnsInLayer = get_number_turns(layer);
            auto alignment = layer.get_turns_alignment().value();

            auto bobbin = resolve_bobbin();
            auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
            double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

            double wireRadius;
            if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR) {
                wireRadius = windingWindowRadialHeight - layer.get_coordinates()[0] - wireWidth / 2;
            }
            else {
                wireRadius = windingWindowRadialHeight - layer.get_coordinates()[0];
            }
            double wireAngle = wound_distance_to_angle(wireHeight, wireRadius);
            if (_strict && (wireRadius <= 0 || wireAngle > 180 || std::isnan(wireAngle))) {
                // Turns won't fit
                return false;
            }

            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                totalLayerAngle = physicalTurnsInLayer * wireAngle;

                currentTurnRadialHeightIncrement = 0;
                currentTurnCenterRadialHeight = roundFloat(layer.get_coordinates()[0], 9);
                switch (alignment) {
                    case CoilAlignment::CENTERED:
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] - totalLayerAngle / 2 + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::INNER_OR_TOP:
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - totalLayerAngle + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::SPREAD:
                        currentTurnAngleIncrement = roundFloat(layer.get_dimensions()[1] / physicalTurnsInLayer, 9);
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 + currentTurnAngleIncrement / 2, 9);
                        break;
                }

            } 
            else {
                throw std::invalid_argument("Only overlapping layers allowed in toroids");
            }

            if (!layer.get_winding_style()) {
                layer.set_winding_style(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            }

            if (layer.get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    int64_t numberTurns = round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
                    for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                        Turn turn;
                        turn.set_coordinates(std::vector<double>{currentTurnCenterRadialHeight, currentTurnCenterAngle});
                        turn.set_layer(layer.get_name());
                        if (bobbinColumnShape == ColumnShape::ROUND) {
                            turn.set_length(2 * std::numbers::pi * (currentTurnCenterRadialHeight + bobbinColumnWidth));
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::OBLONG) {
                            turn.set_length(2 * std::numbers::pi * (currentTurnCenterRadialHeight + bobbinColumnWidth) + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                            double currentTurnCornerRadius = turn.get_coordinates()[0];
                            turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else {
                            throw std::runtime_error("only round or rectangular columns supported for bobbins");
                        }
                        turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                        turn.set_orientation(TurnOrientation::CLOCKWISE);
                        turn.set_parallel(parallelIndex);
                        turn.set_section(layer.get_section().value());
                        turn.set_winding(partialWinding.get_winding());
                        turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                        turn.set_rotation(currentTurnCenterAngle);
                        turn.set_coordinate_system(CoordinateSystem::POLAR);

                        turns.push_back(turn);
                        currentTurnCenterRadialHeight += currentTurnRadialHeightIncrement;
                        currentTurnCenterAngle += currentTurnAngleIncrement;
                        currentTurnIndex[windingIndex][parallelIndex]++; 
                    }
                }
            }
            else {
                int64_t firstParallelIndex = 0;
                while (roundFloat(partialWinding.get_parallels_proportion()[firstParallelIndex], 10) == 0) {
                    firstParallelIndex++;
                }
                int64_t numberTurns = round(partialWinding.get_parallels_proportion()[firstParallelIndex] * get_number_turns(windingIndex));
                for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                        if (roundFloat(partialWinding.get_parallels_proportion()[parallelIndex], 10) > 0) {
                            Turn turn;
                            turn.set_coordinates(std::vector<double>{currentTurnCenterRadialHeight, currentTurnCenterAngle});
                            turn.set_layer(layer.get_name());
                            if (bobbinColumnShape == ColumnShape::ROUND) {
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterRadialHeight);
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::OBLONG) {
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterRadialHeight + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                double currentTurnCornerRadius = currentTurnCenterRadialHeight;
                                turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                                if (turn.get_length() < 0) {
                                    return false;
                                }
                            }
                            else {
                                throw std::runtime_error("only round or rectangular columns supported for bobbins");
                            }
                            turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                            turn.set_orientation(TurnOrientation::CLOCKWISE);
                            turn.set_parallel(parallelIndex);
                            turn.set_section(layer.get_section().value());
                            turn.set_winding(partialWinding.get_winding());
                            turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                            turn.set_rotation(currentTurnCenterAngle);
                            turn.set_coordinate_system(CoordinateSystem::POLAR);


                            turns.push_back(turn);
                            currentTurnCenterRadialHeight += currentTurnRadialHeightIncrement;
                            currentTurnCenterAngle += currentTurnAngleIncrement;
                            currentTurnIndex[windingIndex][parallelIndex]++; 
                        }
                    }
                }
            }
        }
    }

    set_turns_description(turns);

    convert_turns_to_cartesian_coordinates();
    return true;
}

bool Coil::wind_by_planar_turns(double borderToWireDistance, double wireToWireDistance) {
    set_turns_description(std::nullopt);
    if (!get_layers_description()) {
        return false;
    }
    auto wirePerWinding = get_wires();

    std::vector<std::vector<int64_t>> currentTurnIndex;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        currentTurnIndex.push_back(std::vector<int64_t>(get_number_parallels(windingIndex), 0));
    }
    auto bobbin = resolve_bobbin();
    auto bobbinColumnShape = bobbin.get_processed_description().value().get_column_shape();
    auto bobbinColumnDepth = bobbin.get_processed_description().value().get_column_depth();
    double bobbinColumnWidth;
    if (bobbin.get_processed_description().value().get_column_width()) {
        bobbinColumnWidth = bobbin.get_processed_description().value().get_column_width().value();
    }
    else {
        auto bobbinWindingWindow = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows()[0];
        double bobbinWindingWindowWidth = bobbinWindingWindow.get_width().value();
        double bobbinWindingWindowCenterWidth = bobbinWindingWindow.get_coordinates().value()[0];
        bobbinColumnWidth = bobbinWindingWindowCenterWidth - bobbinWindingWindowWidth / 2;
    }

    auto layers = get_layers_description().value();

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        if (wirePerWinding[windingIndex].get_type() == WireType::PLANAR) {
            auto conductionLayers = get_layers_by_type(ElectricalType::CONDUCTION);
            if (conductionLayers.size() > settings->get_coil_maximum_layers_planar()) {
                return false;
            }
        }

    }

    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            double totalLayerHeight;
            double totalLayerWidth;
            if (layer.get_partial_windings().size() > 1) {
                throw std::runtime_error("More than one winding per layer not supported yet");
            }
            auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
            double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
            auto physicalTurnsInLayer = get_number_turns(layer);
            double currentTurnWidthIncrement = wireWidth + wireToWireDistance;
            double currentTurnHeightIncrement = 0;
            double currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + wireWidth / 2, 9) + borderToWireDistance;
            double currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1], 9);

            auto alignment = layer.get_turns_alignment().value();

            if (!layer.get_winding_style()) {
                layer.set_winding_style(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            }

            if (layer.get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    int64_t numberTurns = round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
                    double totalWidthNeeded = borderToWireDistance * 2 + numberTurns * wireWidth + (numberTurns - 1) * wireToWireDistance;
                    if(_strict && totalWidthNeeded > layer.get_dimensions()[0]) {
                        return false;
                    }

                    for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                        Turn turn;
                        turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                        turn.set_layer(layer.get_name());
                        if (bobbinColumnShape == ColumnShape::ROUND) {
                            turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth);
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::OBLONG) {
                            turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                            double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                            turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);

                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else {
                            throw std::runtime_error("only round or rectangular columns supported for bobbins");
                        }
                        turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                        turn.set_orientation(TurnOrientation::CLOCKWISE);
                        turn.set_parallel(parallelIndex);
                        turn.set_section(layer.get_section().value());
                        turn.set_winding(partialWinding.get_winding());
                        turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                        turn.set_rotation(0);
                        turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

                        turns.push_back(turn);
                        currentTurnCenterWidth += currentTurnWidthIncrement;
                        currentTurnCenterHeight -= currentTurnHeightIncrement;
                        currentTurnIndex[windingIndex][parallelIndex]++; 
                    }
                }
            }
            else {
                int64_t firstParallelIndex = 0;
                while (roundFloat(partialWinding.get_parallels_proportion()[firstParallelIndex], 10) == 0) {
                    firstParallelIndex++;
                }
                int64_t numberTurns = round(partialWinding.get_parallels_proportion()[firstParallelIndex] * get_number_turns(windingIndex));
                for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                        if (roundFloat(partialWinding.get_parallels_proportion()[parallelIndex], 10) > 0) {
                            Turn turn;
                            turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                            turn.set_layer(layer.get_name());
                            if (bobbinColumnShape == ColumnShape::ROUND) {
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth);
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::OBLONG) {
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                                turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                                if (turn.get_length() < 0) {
                                    return false;
                                }
                            }
                            else {
                                throw std::runtime_error("only round or rectangular columns supported for bobbins");
                            }
                            turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                            turn.set_orientation(TurnOrientation::CLOCKWISE);
                            turn.set_parallel(parallelIndex);
                            turn.set_section(layer.get_section().value());
                            turn.set_winding(partialWinding.get_winding());
                            turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                            turn.set_rotation(0);
                            turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

                            turns.push_back(turn);
                            currentTurnCenterWidth += currentTurnWidthIncrement;
                            currentTurnCenterHeight -= currentTurnHeightIncrement;
                            currentTurnIndex[windingIndex][parallelIndex]++; 
                        }
                    }
                }
            }
        }
    }

    set_turns_description(turns);
    return true;
}

std::vector<std::pair<double, std::vector<double>>> Coil::get_collision_distances(std::vector<double> turnCoordinates, std::vector<std::vector<double>> placedTurnsCoordinates, double wireHeight) {
    std::vector<std::pair<double, std::vector<double>>> collisions;
    auto turnCartesianCoordinates = polar_to_cartesian(turnCoordinates);
    for (auto placedTurnCoordinates : placedTurnsCoordinates) {
        auto placedTurnCartesianCoordinates = polar_to_cartesian(placedTurnCoordinates);
        double distance = sqrt(pow(turnCartesianCoordinates[0] - placedTurnCartesianCoordinates[0], 2) + pow(turnCartesianCoordinates[1] - placedTurnCartesianCoordinates[1], 2));
        if (distance - wireHeight < 0) {
            double collisionDistance = wireHeight - distance;
            auto placedCoordinates = placedTurnCoordinates;
            collisions.push_back({collisionDistance, placedCoordinates});
        }

        if (collisions.size() == 2) {
            break;
        }
    }

    return collisions;
}

bool Coil::wind_toroidal_additional_turns() {
    if (!get_layers_description()) {
        return false;
    }
    if (!get_turns_description()) {
        return false;
    }
    auto wirePerWinding = get_wires();
    std::vector<std::vector<int64_t>> currentTurnIndex;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        currentTurnIndex.push_back(std::vector<int64_t>(get_number_parallels(windingIndex), 0));
    }
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    double bobbinColumnWidth;
    if (bobbin.get_processed_description().value().get_column_width()) {
        bobbinColumnWidth = bobbin.get_processed_description().value().get_column_width().value();
    }
    else {
        throw std::runtime_error("Toroids must have their bobbin column set");
    }
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
    auto bobbinColumnShape = bobbin.get_processed_description().value().get_column_shape();
    auto bobbinColumnDepth = bobbin.get_processed_description().value().get_column_depth();

    auto sections = get_sections_description().value();
    auto layers = get_layers_description().value();
    auto turns = get_turns_description().value();
    double currentBaseRadialHeight = -bobbinColumnWidth * 2;
    std::map<size_t, double> maximumAdditionalRadialHeightPerSectionByIndex;
    auto windingOrientation = get_winding_orientation();

    for (auto section : sections) { 
        if (section.get_type() == ElectricalType::CONDUCTION) {
            std::vector<std::vector<double>> placedTurnsCoordinates;
            auto turnsInSection = get_turns_by_section(section.get_name());
            auto partialWinding = section.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            // double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
            double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentBaseRadialHeight -= turnsInSection[0].get_dimensions().value()[0] / 2;
            }
            else {
                currentBaseRadialHeight = -bobbinColumnWidth * 2 - turnsInSection[0].get_dimensions().value()[0] / 2;
            }
            double currentSectionMaximumAdditionalRadialHeight = 0;
            for (auto turn : turnsInSection) {
                auto turnIndex = get_turn_index_by_name(turn.get_name());
                std::vector<double> additionalCoordinates = {-bobbinColumnWidth * 2 - turn.get_coordinates()[0], turn.get_coordinates()[1]};

                // If the turn is not on the first layer of the section, we try to place it in a slot there
                if (roundFloat(turn.get_coordinates()[0] - turn.get_dimensions().value()[0] / 2, 9) > 0) {
                    std::vector<double> newCoordinates = {additionalCoordinates[0], additionalCoordinates[1]};
                    newCoordinates[0] = currentBaseRadialHeight;
                    auto collisions = get_collision_distances(newCoordinates, placedTurnsCoordinates, wireHeight);

                    if (collisions.size() > 0) {
                        bool tryAngularMove = collisions.size() > 0;
                        bool tryReversedAngularMove = collisions.size() > 0;
                        bool previouslyAdditionAngularMovement = false;
                        bool try0Degrees = true;
                        bool tryMinus0Degrees = true;
                        bool try30Degrees = true;
                        bool tryMinus30Degrees = true;
                        bool try45Degrees = true;
                        bool tryMinus45Degrees = true;
                        bool try60Degrees = true;
                        bool tryMinus60Degrees = true;
                        bool tryAvoidingCollisionDistance = true;
                        double previousCollisionDistance = 0;
                        std::vector<double> originalCollidedCoordinate;
                        double restoredHeightAfter60Degrees = 0;

                        double collisionDistance = collisions[0].first;
                        auto collidedCoordinate = collisions[0].second;

                        uint64_t timeout = 1000;
                        while (newCoordinates[0] > additionalCoordinates[0]) {
                            timeout--;
                            if (timeout == 0) {
                                throw std::runtime_error("timeout in wind_toroidal_additional_turns");
                            }
                            if (tryAvoidingCollisionDistance && collisionDistance < 1e-6) {
                                tryAvoidingCollisionDistance = false;
                                if (collidedCoordinate[1] > newCoordinates[1]) {
                                    newCoordinates[1] -= ceilFloat(collisionDistance / std::numbers::pi * 180, 3);
                                }
                                else {
                                    newCoordinates[1] += ceilFloat(collisionDistance / std::numbers::pi * 180, 3);
                                }
                            }
                            else if (tryAngularMove) {
                                tryAngularMove = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                double increment = ceilFloat(wound_distance_to_angle(collisionDistance, currentRadius), 3);
                                if (collidedCoordinate[1] > newCoordinates[1]) {
                                    previouslyAdditionAngularMovement = false;
                                    if (newCoordinates[1] - increment > (section.get_coordinates()[1] - section.get_dimensions()[1] / 2))
                                        newCoordinates[1] -= increment;
                                }
                                else {
                                    previouslyAdditionAngularMovement = true;
                                    if (newCoordinates[1] + increment < (section.get_coordinates()[1] - section.get_dimensions()[1] / 2))
                                        newCoordinates[1] += increment;
                                }
                            }
                            else if (tryReversedAngularMove) {
                                tryReversedAngularMove = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                double currentAngleCollision = ceilFloat(wound_distance_to_angle(previousCollisionDistance, currentRadius), 3);
                                double currentWireAngle = ceilFloat(wound_distance_to_angle(wireHeight, currentRadius), 3);
                                double currentAngleMovement = currentWireAngle + (currentWireAngle - currentAngleCollision);

                                if (previouslyAdditionAngularMovement) {
                                    if (newCoordinates[1] - currentAngleMovement > (section.get_coordinates()[1] - section.get_dimensions()[1] / 2))
                                        newCoordinates[1] -= currentAngleMovement;
                                }
                                else {
                                    if (newCoordinates[1] + currentAngleMovement < (section.get_coordinates()[1] - section.get_dimensions()[1] / 2))
                                        newCoordinates[1] += currentAngleMovement;
                                }
                            }
                            else if (try0Degrees) {
                                try0Degrees = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                restoredHeightAfter60Degrees = newCoordinates[0];
                                newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(0);
                                newCoordinates[1] = originalCollidedCoordinate[1] + ceilFloat(wound_distance_to_angle(wireHeight * cos(0), currentRadius), 3);
                            }
                            else if (tryMinus0Degrees) {
                                tryMinus0Degrees = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(0);
                                newCoordinates[1] = originalCollidedCoordinate[1] - ceilFloat(wound_distance_to_angle(wireHeight * cos(0), currentRadius), 3);
                            }
                            else if (try30Degrees) {
                                try30Degrees = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 6);
                                newCoordinates[1] = originalCollidedCoordinate[1] + ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 6), currentRadius), 3);
                            }
                            else if (tryMinus30Degrees) {
                                tryMinus30Degrees = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 6);
                                newCoordinates[1] = originalCollidedCoordinate[1] - ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 6), currentRadius), 3);
                            }
                            else if (try45Degrees) {
                                try45Degrees = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 4);
                                newCoordinates[1] = originalCollidedCoordinate[1] + ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 4), currentRadius), 3);
                            }
                            else if (tryMinus45Degrees) {
                                tryMinus45Degrees = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 4);
                                newCoordinates[1] = originalCollidedCoordinate[1] - ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 4), currentRadius), 3);
                            }
                            else if (try60Degrees) {
                                try60Degrees = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 3);
                                newCoordinates[1] = originalCollidedCoordinate[1] + ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 3), currentRadius), 3);
                            }
                            else if (tryMinus60Degrees) {
                                tryMinus60Degrees = false;
                                double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 3);
                                newCoordinates[1] = originalCollidedCoordinate[1] - ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 3), currentRadius), 3);
                            }
                            else {
                                try0Degrees = true;
                                tryMinus0Degrees = true;
                                try30Degrees = true;
                                tryMinus30Degrees = true;
                                try45Degrees = true;
                                tryMinus45Degrees = true;
                                try60Degrees = true;
                                tryMinus60Degrees = true;
                                tryAngularMove = true;
                                tryAvoidingCollisionDistance = true;
                                previousCollisionDistance = 0;
                                if (restoredHeightAfter60Degrees != 0) {
                                    newCoordinates[0] = restoredHeightAfter60Degrees;
                                    restoredHeightAfter60Degrees = 0;
                                }
                                newCoordinates[0] -= turn.get_dimensions().value()[0] / 2;
                                newCoordinates[1] = additionalCoordinates[1];
                            }
                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                            double currentWireAngle = ceilFloat(wound_distance_to_angle(wireHeight, currentRadius), 3);

                            if (newCoordinates[1] - currentWireAngle / 2 < (section.get_coordinates()[1] - section.get_dimensions()[1] / 2)) {
                                newCoordinates[1] = additionalCoordinates[1];
                            }
                            if (newCoordinates[1] + currentWireAngle / 2 > (section.get_coordinates()[1] + section.get_dimensions()[1] / 2)) {
                                newCoordinates[1] = additionalCoordinates[1];
                            }

                            collisions = get_collision_distances(newCoordinates, placedTurnsCoordinates, wireHeight);
                            if (collisions.size() == 0) {
                                break;
                            }
                            collidedCoordinate = collisions[0].second;
                            if (previousCollisionDistance == 0) {
                                originalCollidedCoordinate = collidedCoordinate;
                            }
                            previousCollisionDistance = collisionDistance;
                            collisionDistance = collisions[0].first;
                        }
                    }
                    additionalCoordinates = newCoordinates;
                }
                currentSectionMaximumAdditionalRadialHeight = std::min(currentSectionMaximumAdditionalRadialHeight, additionalCoordinates[0]);
                turn.set_additional_coordinates(std::vector<std::vector<double>>{additionalCoordinates});

                if (bobbinColumnShape == ColumnShape::ROUND) {
                    double b = (turn.get_coordinates()[0] - turn.get_additional_coordinates().value()[0][0]) / 2;
                    double a = turn.get_coordinates()[0];
                    // Ramanujan  approximation for ellipse perimeter
                    double perimeter = std::numbers::pi * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b)));
                    turns[turnIndex].set_length(perimeter);
                    if (turns[turnIndex].get_length() < 0) {
                        throw std::runtime_error("Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " turns[turnIndex].get_coordinates()[0]: " + std::to_string(turns[turnIndex].get_coordinates()[0]));
                    }
                }
                else if (bobbinColumnShape == ColumnShape::OBLONG) {
                    double b = (turn.get_coordinates()[0] - turn.get_additional_coordinates().value()[0][0]) / 2;
                    double a = turn.get_coordinates()[0];
                    // Ramanujan  approximation for ellipse perimeter
                    double perimeter = std::numbers::pi * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b))) + 4 * (bobbinColumnDepth - bobbinColumnWidth);
                    turns[turnIndex].set_length(perimeter);
                    if (turns[turnIndex].get_length() < 0) {
                        throw std::runtime_error("Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " turns[turnIndex].get_coordinates()[0]: " + std::to_string(turns[turnIndex].get_coordinates()[0]));
                    }
                }
                else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                    double currentInternalTurnCornerRadius = turns[turnIndex].get_coordinates()[0];
                    double currentExternalTurnCornerRadius = -turn.get_additional_coordinates().value()[0][0] - 2 * bobbinColumnWidth;
                    double maximumVerticalDistance = currentInternalTurnCornerRadius * 2 + 2 * bobbinColumnDepth;
                    double externalVerticalStraightDistance = maximumVerticalDistance - 2 * currentExternalTurnCornerRadius;
                    turns[turnIndex].set_length(2 * bobbinColumnDepth + 4 * bobbinColumnWidth + externalVerticalStraightDistance + std::numbers::pi * currentInternalTurnCornerRadius + std::numbers::pi * currentExternalTurnCornerRadius);

                    if (turns[turnIndex].get_length() < 0) {
                        throw std::runtime_error("Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " bobbinColumnDepth: " + std::to_string(bobbinColumnDepth)  + " bobbinColumnWidth: " + std::to_string(bobbinColumnWidth)  + " currentExternalTurnCornerRadius: " + std::to_string(currentExternalTurnCornerRadius));
                    }
                }
                else {
                    throw std::runtime_error("only round or rectangular columns supported for bobbins");
                }

                turns[turnIndex] = turn;
                placedTurnsCoordinates.push_back(additionalCoordinates);
            }

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentSectionMaximumAdditionalRadialHeight -= turnsInSection[0].get_dimensions().value()[0] / 2;
                auto sectionIndex = get_section_index_by_name(section.get_name());
                currentBaseRadialHeight = currentSectionMaximumAdditionalRadialHeight;
                if (sectionIndex < sections.size() - 1) {
                    maximumAdditionalRadialHeightPerSectionByIndex[sectionIndex + 1] = currentBaseRadialHeight;
                }
            }
        }
        else {
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentBaseRadialHeight -= section.get_dimensions()[0];
            }
        }

    }
    set_turns_description(turns);

    for (auto [sectionIndex, radialHeight] : maximumAdditionalRadialHeightPerSectionByIndex) {
        auto layersInSection = get_layers_by_section(sections[sectionIndex].get_name());

        double currentRadialHeight = radialHeight;
        for (auto layer : layersInSection) {
            if (layer.get_type() == ElectricalType::INSULATION) {
                auto layerIndex = get_layer_index_by_name(layer.get_name());
                currentRadialHeight -= layer.get_dimensions()[0] / 2;
                std::vector<double> additionalCoordinates = {currentRadialHeight, layer.get_coordinates()[1]};
                layers[layerIndex].set_additional_coordinates(std::vector<std::vector<double>>{additionalCoordinates});
                currentRadialHeight -= layer.get_dimensions()[0] / 2;
            }
        }
    }
    set_layers_description(layers);

    return true;
}

std::vector<double> Coil::get_aligned_section_dimensions_rectangular_window(size_t sectionIndex) {
    auto sections = get_sections_description().value();
    if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
        sections[sectionIndex].set_margin(std::vector<double>{0, 0});
    }

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();
    auto windingOrientation = get_winding_orientation();

    if (sections.size() == 0) {
        throw std::runtime_error("No sections in coil");
    }
    double totalSectionsWidth = 0;
    double totalSectionsHeight = 0;
    for (size_t auxSectionIndex = 0; auxSectionIndex < sections.size(); ++auxSectionIndex) {
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            totalSectionsWidth += sections[auxSectionIndex].get_dimensions()[0];
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsHeight = std::max(totalSectionsHeight, sections[auxSectionIndex].get_dimensions()[1]);
            }
        }
        else {
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsWidth = std::max(totalSectionsWidth, sections[auxSectionIndex].get_dimensions()[0]);
            }
            totalSectionsHeight += sections[auxSectionIndex].get_dimensions()[1];
        }
    }

    double currentCoilWidth;
    double currentCoilHeight;
    double paddingAmongSectionWidth = 0;
    double paddingAmongSectionHeight = 0;
    auto turnsAlignment = get_turns_alignment(sections[sectionIndex].get_name());

    auto sectionAlignment = get_section_alignment();
    switch (sectionAlignment) {
        case CoilAlignment::INNER_OR_TOP:

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - sections[sectionIndex].get_margin().value()[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + sections[sectionIndex].get_margin().value()[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilHeight = 0;
                            double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - sections[sectionIndex].get_margin().value()[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                            double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + sections[sectionIndex].get_margin().value()[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                            currentCoilHeight = std::min(currentCoilHeight, currentCoilHeightTop);
                            currentCoilHeight = std::max(currentCoilHeight, currentCoilHeightBottom);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilHeight = -sections[sectionIndex].get_margin().value()[0] / 2 + sections[sectionIndex].get_margin().value()[1] / 2;
                        break;
                    default:
                        throw std::runtime_error("No such section alignment");
                }
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - sections[sectionIndex].get_margin().value()[1] - sections[sectionIndex].get_dimensions()[0];
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                            double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                            double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - sections[sectionIndex].get_margin().value()[1] - sections[sectionIndex].get_dimensions()[0];
                            currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                            currentCoilWidth = std::min(currentCoilWidth, currentCoilWidthRight);
                            break;
                        }
                    case CoilAlignment::SPREAD:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                        break;
                    default:
                        throw std::runtime_error("No such section alignment");
                }
            }
            break;
        case CoilAlignment::OUTER_OR_BOTTOM:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - totalSectionsWidth;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - sections[sectionIndex].get_margin().value()[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + sections[sectionIndex].get_margin().value()[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilHeight = 0;
                            double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - sections[sectionIndex].get_margin().value()[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                            double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + sections[sectionIndex].get_margin().value()[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                            currentCoilHeight = std::min(currentCoilHeight, currentCoilHeightTop);
                            currentCoilHeight = std::max(currentCoilHeight, currentCoilHeightBottom);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilHeight = -sections[sectionIndex].get_margin().value()[0] / 2 + sections[sectionIndex].get_margin().value()[1] / 2;
                        break;
                    default:
                        throw std::runtime_error("No such section alignment");
                }
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + totalSectionsHeight;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - sections[sectionIndex].get_margin().value()[1] - sections[sectionIndex].get_dimensions()[0];
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                            double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                            double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - sections[sectionIndex].get_margin().value()[1] - sections[sectionIndex].get_dimensions()[0];
                            currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                            currentCoilWidth = std::min(currentCoilWidth, currentCoilWidthRight);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                        break;
                    default:
                        throw std::runtime_error("No such section alignment");
                }
            }
            break;
        case CoilAlignment::SPREAD:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - sections[sectionIndex].get_margin().value()[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + sections[sectionIndex].get_margin().value()[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilHeight = 0;
                            double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - sections[sectionIndex].get_margin().value()[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                            double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + sections[sectionIndex].get_margin().value()[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                            currentCoilHeight = std::min(currentCoilHeight, currentCoilHeightTop);
                            currentCoilHeight = std::max(currentCoilHeight, currentCoilHeightBottom);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilHeight = -sections[sectionIndex].get_margin().value()[0] / 2 + sections[sectionIndex].get_margin().value()[1] / 2;
                        break;
                    default:
                        throw std::runtime_error("No such section alignment");
                }
                paddingAmongSectionWidth = windingWindows[0].get_width().value() - totalSectionsWidth;
                if (sections.size() > 1) {
                    paddingAmongSectionWidth /= sections.size() - 1;
                }
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2;
                paddingAmongSectionHeight = windingWindows[0].get_height().value() - totalSectionsHeight;
                if (sections.size() > 1) {
                    paddingAmongSectionHeight /= sections.size() - 1;
                }
                else {
                    currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + totalSectionsHeight / 2;
                }

                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - sections[sectionIndex].get_margin().value()[1] - sections[sectionIndex].get_dimensions()[0];
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                            double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                            double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - sections[sectionIndex].get_margin().value()[1] - sections[sectionIndex].get_dimensions()[0];
                            currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                            currentCoilWidth = std::min(currentCoilWidth, currentCoilWidthRight);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                        break;
                    default:
                        throw std::runtime_error("No such section alignment");
                }
            }
            break;
        case CoilAlignment::CENTERED:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - sections[sectionIndex].get_margin().value()[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + sections[sectionIndex].get_margin().value()[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilHeight = 0;
                            double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - sections[sectionIndex].get_margin().value()[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                            double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + sections[sectionIndex].get_margin().value()[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                            currentCoilHeight = std::min(currentCoilHeight, currentCoilHeightTop);
                            currentCoilHeight = std::max(currentCoilHeight, currentCoilHeightBottom);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilHeight = -sections[sectionIndex].get_margin().value()[0] / 2 + sections[sectionIndex].get_margin().value()[1] / 2;
                        break;
                    default:
                        throw std::runtime_error("No such section alignment");
                }
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + totalSectionsHeight / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - sections[sectionIndex].get_margin().value()[1] - sections[sectionIndex].get_dimensions()[0];
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                            double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                            double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - sections[sectionIndex].get_margin().value()[1] - sections[sectionIndex].get_dimensions()[0];
                            if (currentCoilWidthLeft < 0) {
                                throw std::invalid_argument("currentCoilWidthLeft cannot be less than 0: " + std::to_string(currentCoilWidthLeft));
                            }
                            if (currentCoilWidthRight < 0) {
                                throw std::invalid_argument("currentCoilWidthRight cannot be less than 0: " + std::to_string(currentCoilWidthRight));
                            }
                            currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                            if (currentCoilWidthRight >= 0) {
                                currentCoilWidth = std::min(currentCoilWidth, currentCoilWidthRight);
                            }
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + sections[sectionIndex].get_margin().value()[0];
                        break;
                    default:
                        throw std::runtime_error("No such section alignment");
                }
            }
            break;
        default:
            throw std::runtime_error("No such section alignment");

    }

    return {currentCoilWidth, currentCoilHeight, paddingAmongSectionWidth, paddingAmongSectionHeight};
}

std::vector<double> Coil::get_aligned_section_dimensions_round_window(size_t sectionIndex) {
    auto sections = get_sections_description().value();
    if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
        sections[sectionIndex].set_margin(std::vector<double>{0, 0});
    }

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowAngle = windingWindows[0].get_angle().value();
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
    auto windingOrientation = get_winding_orientation();

    if (sections.size() == 0) {
        throw std::runtime_error("No sections in coil");
    }
    double totalSectionsRadialHeight = 0;
    double totalSectionsAngle = 0;
    for (size_t auxSectionIndex = 0; auxSectionIndex < sections.size(); ++auxSectionIndex) {
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            totalSectionsRadialHeight += sections[auxSectionIndex].get_dimensions()[0];
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsAngle = std::max(totalSectionsAngle, sections[auxSectionIndex].get_dimensions()[1]);
            }
        }
        else {
            double marginAngle0 = 0;
            double marginAngle1 = 0;
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsRadialHeight = std::max(totalSectionsRadialHeight, sections[auxSectionIndex].get_dimensions()[0]);
                double lastLayerMaximumRadius = windingWindowRadialHeight - (sections[auxSectionIndex].get_coordinates()[0] + sections[auxSectionIndex].get_dimensions()[0] / 2);
                marginAngle0 = wound_distance_to_angle(sections[auxSectionIndex].get_margin().value()[0], lastLayerMaximumRadius);
                marginAngle1 = wound_distance_to_angle(sections[auxSectionIndex].get_margin().value()[1], lastLayerMaximumRadius);
            }
            totalSectionsAngle += sections[auxSectionIndex].get_dimensions()[1] + marginAngle0 + marginAngle1;
        }
    }

    double currentCoilRadialHeight;
    double currentCoilAngle;
    double paddingAmongSectionRadialHeight = 0;
    double paddingAmongSectionAngle = 0;
    double marginAngle0 = 0;

    if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
        double lastLayerMaximumRadius = windingWindowRadialHeight - (sections[sectionIndex].get_coordinates()[0] + sections[sectionIndex].get_dimensions()[0] / 2);
        marginAngle0 = wound_distance_to_angle(sections[sectionIndex].get_margin().value()[0], lastLayerMaximumRadius);
    }
    auto turnsAlignment = get_turns_alignment(sections[sectionIndex].get_name());

    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        currentCoilRadialHeight = 0;
        switch (turnsAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                currentCoilAngle = sections[sectionIndex].get_dimensions()[1] / 2;
                break;
            case CoilAlignment::OUTER_OR_BOTTOM:
                currentCoilAngle = windingWindowAngle - sections[sectionIndex].get_dimensions()[1] / 2;
                break;
            case CoilAlignment::CENTERED:
                {
                    currentCoilAngle = 180;
                    break;
                }
                break;
            case CoilAlignment::SPREAD:
                currentCoilAngle = sections[sectionIndex].get_dimensions()[1] / 2;
                break;
            default:
                throw std::runtime_error("No such section alignment");
        }
    }
    else {
        auto sectionAlignment = get_section_alignment();
        switch (sectionAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                currentCoilAngle = sections[sectionIndex].get_coordinates()[1] - sections[sectionIndex].get_dimensions()[1] / 2 - marginAngle0;
                break;
            case CoilAlignment::OUTER_OR_BOTTOM:
                currentCoilAngle = windingWindowAngle - totalSectionsAngle;
                break;
            case CoilAlignment::SPREAD:
                currentCoilAngle = sections[sectionIndex].get_coordinates()[1];
                paddingAmongSectionAngle = windingWindows[0].get_angle().value() - totalSectionsAngle;
                if (sections.size() > 1) {
                    paddingAmongSectionAngle /= sections.size() - 1;
                }
                else {
                    currentCoilAngle = windingWindowAngle / 2 + totalSectionsAngle / 2;
                }
                break;
            case CoilAlignment::CENTERED:
                currentCoilAngle = windingWindowAngle / 2 - totalSectionsAngle / 2;
                break;
            default:
                throw std::runtime_error("No such section alignment");
        }
    }

    return {currentCoilRadialHeight, currentCoilAngle, paddingAmongSectionRadialHeight, paddingAmongSectionAngle};
}

bool Coil::delimit_and_compact() {
    auto bobbin = resolve_bobbin();

    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        return delimit_and_compact_rectangular_window();
    }
    else {
        return delimit_and_compact_round_window();
    }
}

WiringTechnology Coil::get_coil_type(size_t groupIndex) {
    if (!get_groups_description()) {
        return WiringTechnology::WOUND;
    }
    auto groups = get_groups_description().value();
    if (groupIndex >= groups.size()) {
        throw std::runtime_error("Non existing group index");
    }
    auto group = get_groups_description().value()[groupIndex];
    return group.get_type();
}


bool Coil::delimit_and_compact_rectangular_window() {
    // Delimit
        
    auto groupType = get_coil_type();

    if (!get_sections_description()) {
        throw std::runtime_error("No sections to delimit");
    }

    if (get_layers_description()) {
        auto layers = get_layers_description().value();
        if (get_turns_description() && groupType == WiringTechnology::WOUND) {
            for (size_t i = 0; i < layers.size(); ++i) {
                if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                    auto turnsInLayer = get_turns_by_layer(layers[i].get_name());
                    auto layerCoordinates = layers[i].get_coordinates();
                    double currentLayerMaximumWidth = (turnsInLayer[0].get_coordinates()[0] - layerCoordinates[0]) + turnsInLayer[0].get_dimensions().value()[0] / 2;
                    double currentLayerMinimumWidth = (turnsInLayer[0].get_coordinates()[0] - layerCoordinates[0]) - turnsInLayer[0].get_dimensions().value()[0] / 2;
                    double currentLayerMaximumHeight = (turnsInLayer[0].get_coordinates()[1] - layerCoordinates[1]) + turnsInLayer[0].get_dimensions().value()[1] / 2;
                    double currentLayerMinimumHeight = (turnsInLayer[0].get_coordinates()[1] - layerCoordinates[1]) - turnsInLayer[0].get_dimensions().value()[1] / 2;
                    for (auto& turn : turnsInLayer) {
                        currentLayerMaximumWidth = std::max(currentLayerMaximumWidth, (turn.get_coordinates()[0] - layerCoordinates[0]) + turn.get_dimensions().value()[0] / 2);
                        currentLayerMinimumWidth = std::min(currentLayerMinimumWidth, (turn.get_coordinates()[0] - layerCoordinates[0]) - turn.get_dimensions().value()[0] / 2);
                        currentLayerMaximumHeight = std::max(currentLayerMaximumHeight, (turn.get_coordinates()[1] - layerCoordinates[1]) + turn.get_dimensions().value()[1] / 2);
                        currentLayerMinimumHeight = std::min(currentLayerMinimumHeight, (turn.get_coordinates()[1] - layerCoordinates[1]) - turn.get_dimensions().value()[1] / 2);
                    }
                    layers[i].set_coordinates(std::vector<double>({layerCoordinates[0] + (currentLayerMaximumWidth + currentLayerMinimumWidth) / 2,
                                                               layerCoordinates[1] + (currentLayerMaximumHeight + currentLayerMinimumHeight) / 2}));
                    layers[i].set_dimensions(std::vector<double>({currentLayerMaximumWidth - currentLayerMinimumWidth,
                                                               currentLayerMaximumHeight - currentLayerMinimumHeight}));

                }
                set_layers_description(layers);
            }
        }

        auto sections = get_sections_description().value();
        for (size_t i = 0; i < sections.size(); ++i) {
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                auto layersInSection = get_layers_by_section(sections[i].get_name());
                if (layersInSection.size() == 0) {
                    throw std::runtime_error("No layers in section: " + sections[i].get_name());
                }
                auto sectionCoordinates = sections[i].get_coordinates();
                double currentSectionMaximumWidth = (layersInSection[0].get_coordinates()[0] - sectionCoordinates[0]) + layersInSection[0].get_dimensions()[0] / 2;
                double currentSectionMinimumWidth = (layersInSection[0].get_coordinates()[0] - sectionCoordinates[0]) - layersInSection[0].get_dimensions()[0] / 2;
                double currentSectionMaximumHeight = (layersInSection[0].get_coordinates()[1] - sectionCoordinates[1]) + layersInSection[0].get_dimensions()[1] / 2;
                double currentSectionMinimumHeight = (layersInSection[0].get_coordinates()[1] - sectionCoordinates[1]) - layersInSection[0].get_dimensions()[1] / 2;

                for (auto& layer : layersInSection) {
                    currentSectionMaximumWidth = std::max(currentSectionMaximumWidth, (layer.get_coordinates()[0] - sectionCoordinates[0]) + layer.get_dimensions()[0] / 2);
                    currentSectionMinimumWidth = std::min(currentSectionMinimumWidth, (layer.get_coordinates()[0] - sectionCoordinates[0]) - layer.get_dimensions()[0] / 2);
                    currentSectionMaximumHeight = std::max(currentSectionMaximumHeight, (layer.get_coordinates()[1] - sectionCoordinates[1]) + layer.get_dimensions()[1] / 2);
                    currentSectionMinimumHeight = std::min(currentSectionMinimumHeight, (layer.get_coordinates()[1] - sectionCoordinates[1]) - layer.get_dimensions()[1] / 2);
                }
                sections[i].set_coordinates(std::vector<double>({sectionCoordinates[0] + (currentSectionMaximumWidth + currentSectionMinimumWidth) / 2,
                                                           sectionCoordinates[1] + (currentSectionMaximumHeight + currentSectionMinimumHeight) / 2}));
                sections[i].set_dimensions(std::vector<double>({currentSectionMaximumWidth - currentSectionMinimumWidth,
                                                       currentSectionMaximumHeight - currentSectionMinimumHeight}));
            }
        }
        set_sections_description(sections);

        if (get_groups_description() && groupType == WiringTechnology::PRINTED) {
            auto groups = get_groups_description().value();
            for (size_t i = 0; i < groups.size(); ++i) {
                auto sectionsInGroup = get_sections_by_group(groups[i].get_name());
                if (sectionsInGroup.size() == 0) {
                    throw std::runtime_error("No sections in group: " + groups[i].get_name());
                }
                auto groupCoordinates = groups[i].get_coordinates();
                double currentGroupMaximumHeight = (sectionsInGroup[0].get_coordinates()[1] - groupCoordinates[1]) + sectionsInGroup[0].get_dimensions()[1] / 2;
                double currentGroupMinimumHeight = (sectionsInGroup[0].get_coordinates()[1] - groupCoordinates[1]) - sectionsInGroup[0].get_dimensions()[1] / 2;

                for (auto& section : sectionsInGroup) {
                    currentGroupMaximumHeight = std::max(currentGroupMaximumHeight, (section.get_coordinates()[1] - groupCoordinates[1]) + section.get_dimensions()[1] / 2);
                    currentGroupMinimumHeight = std::min(currentGroupMinimumHeight, (section.get_coordinates()[1] - groupCoordinates[1]) - section.get_dimensions()[1] / 2);
                }
                groups[i].set_coordinates(std::vector<double>({groupCoordinates[0], groupCoordinates[1] + (currentGroupMaximumHeight + currentGroupMinimumHeight) / 2}));
                groups[i].set_dimensions(std::vector<double>({groups[i].get_dimensions()[0], currentGroupMaximumHeight - currentGroupMinimumHeight}));
            }
            set_groups_description(groups);
        }
    }

     // Compact
    if (get_sections_description() && groupType == WiringTechnology::WOUND) {
        auto sections = get_sections_description().value();

        std::vector<std::vector<double>> alignedSectionDimensionsPerSection;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            alignedSectionDimensionsPerSection.push_back(get_aligned_section_dimensions_rectangular_window(sectionIndex));
        }

        double currentCoilWidth = alignedSectionDimensionsPerSection[0][0];
        double currentCoilHeight = alignedSectionDimensionsPerSection[0][1];
        double paddingAmongSectionWidth = alignedSectionDimensionsPerSection[0][2];
        double paddingAmongSectionHeight = alignedSectionDimensionsPerSection[0][3];

        std::vector<Turn> turns;
        if (get_turns_description()) {
            turns = get_turns_description().value();
        }

        std::vector<Layer> layers;
        if (get_layers_description()) {
            layers = get_layers_description().value();
        }

        auto bobbinColumnShape = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_shape();
        auto bobbinColumnDepth = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_depth();
        double bobbinColumnWidth;
        if (std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_width()) {
            bobbinColumnWidth = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_width().value();
        }
        else {
            auto bobbinWindingWindow = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows()[0];
            double bobbinWindingWindowWidth = bobbinWindingWindow.get_width().value();
            double bobbinWindingWindowCenterWidth = bobbinWindingWindow.get_coordinates().value()[0];
            bobbinColumnWidth = bobbinWindingWindowCenterWidth - bobbinWindingWindowWidth / 2;
        }

        auto windingOrientation = get_winding_orientation();

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilHeight = alignedSectionDimensionsPerSection[sectionIndex][1];
                currentCoilWidth += sections[sectionIndex].get_dimensions()[0] / 2;
            }
            else {
                currentCoilHeight -= sections[sectionIndex].get_dimensions()[1] / 2;
                currentCoilWidth = alignedSectionDimensionsPerSection[sectionIndex][0];
            }

            double compactingShiftWidth = sections[sectionIndex].get_coordinates()[0] - currentCoilWidth;
            double compactingShiftHeight = sections[sectionIndex].get_coordinates()[1] - currentCoilHeight;

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                // compactingShiftWidth += sections[sectionIndex].get_dimensions()[0] / 2;
                if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
                    compactingShiftHeight = 0;
                }

            }
            else {
                compactingShiftWidth -= sections[sectionIndex].get_dimensions()[0] / 2;
                if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
                    compactingShiftWidth = 0;
                }
            }

            if (compactingShiftWidth != 0 || compactingShiftHeight != 0) {
                sections[sectionIndex].set_coordinates(std::vector<double>({
                    sections[sectionIndex].get_coordinates()[0] - compactingShiftWidth,
                    sections[sectionIndex].get_coordinates()[1] - compactingShiftHeight
                }));

                for (size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
                    if (layers[layerIndex].get_section().value() == sections[sectionIndex].get_name()){
                        layers[layerIndex].set_coordinates(std::vector<double>({
                            layers[layerIndex].get_coordinates()[0] - compactingShiftWidth,
                            layers[layerIndex].get_coordinates()[1] - compactingShiftHeight
                        }));
                        for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
                            if (turns[turnIndex].get_layer().value() == layers[layerIndex].get_name()){

                                if (bobbinColumnShape == ColumnShape::ROUND || bobbinColumnShape == ColumnShape::OBLONG || bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                    if (turns[turnIndex].get_coordinates()[0] < compactingShiftWidth) {
                                        throw std::runtime_error("Something wrong happened with compactingShiftWidth: " + std::to_string(compactingShiftWidth) +
                                                                 "\nsections[sectionIndex].get_coordinates()[0]: " + std::to_string(sections[sectionIndex].get_coordinates()[0]) +
                                                                 "\ncurrentCoilWidth: " + std::to_string(currentCoilWidth) +
                                                                 "\nturns[turnIndex].get_coordinates()[0]: " + std::to_string(turns[turnIndex].get_coordinates()[0])
                                                                 );
                                    }
                                }
                                else {
                                    throw std::runtime_error("only round or rectangular columns supported for bobbins");
                                }

                                turns[turnIndex].set_coordinates(std::vector<double>({
                                    turns[turnIndex].get_coordinates()[0] - compactingShiftWidth,
                                    turns[turnIndex].get_coordinates()[1] - compactingShiftHeight
                                }));

                                if (bobbinColumnShape == ColumnShape::ROUND) {
                                    turns[turnIndex].set_length(2 * std::numbers::pi * turns[turnIndex].get_coordinates()[0]);
                                    if (turns[turnIndex].get_length() < 0) {
                                        throw std::runtime_error("Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " turns[turnIndex].get_coordinates()[0]: " + std::to_string(turns[turnIndex].get_coordinates()[0]));
                                    }
                                }
                                else if (bobbinColumnShape == ColumnShape::OBLONG) {
                                    turns[turnIndex].set_length(2 * std::numbers::pi * turns[turnIndex].get_coordinates()[0] + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                                    if (turns[turnIndex].get_length() < 0) {
                                        throw std::runtime_error("Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " turns[turnIndex].get_coordinates()[0]: " + std::to_string(turns[turnIndex].get_coordinates()[0]));
                                    }
                                }
                                else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                    double currentTurnCornerRadius = turns[turnIndex].get_coordinates()[0] - bobbinColumnWidth;
                                    turns[turnIndex].set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);

                                    if (turns[turnIndex].get_length() < 0) {
                                        throw std::runtime_error("Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " bobbinColumnDepth: " + std::to_string(bobbinColumnDepth)  + " bobbinColumnWidth: " + std::to_string(bobbinColumnWidth)  + " currentTurnCornerRadius: " + std::to_string(currentTurnCornerRadius));
                                    }
                                }
                                else {
                                    throw std::runtime_error("only round or rectangular columns supported for bobbins");
                                }
                            }
                        }
                    }
                }
            }
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth += sections[sectionIndex].get_dimensions()[0] / 2 + paddingAmongSectionWidth;
            }
            else {
                currentCoilHeight -= sections[sectionIndex].get_dimensions()[1] / 2 + paddingAmongSectionHeight;
            }
        }
        if (get_turns_description()) {
            set_turns_description(turns);
        }
        if (get_layers_description()) {
            set_layers_description(layers);
        }
        set_sections_description(sections);
    }

    // Add extra margin for support if required
    bool fillCoilSectionsWithMarginTape = settings->get_coil_fill_sections_with_margin_tape();


    if (fillCoilSectionsWithMarginTape) {
        auto bobbin = resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        double windingWindowHeight = windingWindowDimensions[1];
        double windingWindowWidth = windingWindowDimensions[0];
        auto sections = get_sections_description().value();
        for (size_t i = 0; i < sections.size(); ++i) {
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                auto sectionOrientation = bobbin.get_winding_window_sections_orientation(0);
                if (sectionOrientation == WindingOrientation::OVERLAPPING) {
                    auto topSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[1] + windingWindowHeight / 2) - (sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2));
                    auto bottomSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[1] - windingWindowHeight / 2) - (sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2));
                    sections[i].set_margin(std::vector<double>{topSpaceBetweenSectionAndBobbin, bottomSpaceBetweenSectionAndBobbin});
                }
                else if (sectionOrientation == WindingOrientation::CONTIGUOUS) {
                    auto innerSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[0] - windingWindowWidth / 2) - (sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2));
                    auto outerSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[0] + windingWindowWidth / 2) - (sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2));
                    sections[i].set_margin(std::vector<double>{innerSpaceBetweenSectionAndBobbin, outerSpaceBetweenSectionAndBobbin});
                }
            }
        }
        set_sections_description(sections);
    }

    return true;
}

bool Coil::delimit_and_compact_round_window() {
    if (get_turns_description()) {
        convert_turns_to_polar_coordinates();
    }

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto windingWindowsRadius = windingWindows[0].get_radial_height().value();


    // Radial Delimit
    if (get_layers_description()) {
        auto layers = get_layers_description().value();
        if (get_turns_description()) {
            for (size_t i = 0; i < layers.size(); ++i) {
                if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                    auto turnsInLayer = get_turns_by_layer(layers[i].get_name());
                    auto layerCoordinates = layers[i].get_coordinates();
                    auto section = get_section_by_name(layers[i].get_section().value());

                    double currentLayerMaximumRadialHeight = (turnsInLayer[0].get_coordinates()[0] - layerCoordinates[0]) + turnsInLayer[0].get_dimensions().value()[0] / 2;
                    double currentLayerMinimumRadialHeight = (turnsInLayer[0].get_coordinates()[0] - layerCoordinates[0]) - turnsInLayer[0].get_dimensions().value()[0] / 2;

                    for (auto& turn : turnsInLayer) {
                        currentLayerMaximumRadialHeight = std::max(currentLayerMaximumRadialHeight, (turn.get_coordinates()[0] - layerCoordinates[0]) + turn.get_dimensions().value()[0] / 2);
                        currentLayerMinimumRadialHeight = std::min(currentLayerMinimumRadialHeight, (turn.get_coordinates()[0] - layerCoordinates[0]) - turn.get_dimensions().value()[0] / 2);
                    }

                    layers[i].set_coordinates(std::vector<double>({layerCoordinates[0] + (currentLayerMaximumRadialHeight + currentLayerMinimumRadialHeight) / 2, layers[i].get_coordinates()[1]}));
                    layers[i].set_dimensions(std::vector<double>({currentLayerMaximumRadialHeight - currentLayerMinimumRadialHeight, layers[i].get_dimensions()[1]}));
                }
                set_layers_description(layers);
            }
        }

        auto sections = get_sections_description().value();
        for (size_t i = 0; i < sections.size(); ++i) {
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                auto layersInSection = get_layers_by_section(sections[i].get_name());
                if (layersInSection.size() == 0) {
                    throw std::runtime_error("No layers in section: " + sections[i].get_name());
                }
                auto sectionCoordinates = sections[i].get_coordinates();
                double currentSectionMaximumRadialHeight = (layersInSection[0].get_coordinates()[0] - sectionCoordinates[0]) + layersInSection[0].get_dimensions()[0] / 2;
                double currentSectionMinimumRadialHeight = (layersInSection[0].get_coordinates()[0] - sectionCoordinates[0]) - layersInSection[0].get_dimensions()[0] / 2;

                for (auto& layer : layersInSection) {
                    currentSectionMaximumRadialHeight = std::max(currentSectionMaximumRadialHeight, (layer.get_coordinates()[0] - sectionCoordinates[0]) + layer.get_dimensions()[0] / 2);
                    currentSectionMinimumRadialHeight = std::min(currentSectionMinimumRadialHeight, (layer.get_coordinates()[0] - sectionCoordinates[0]) - layer.get_dimensions()[0] / 2);
                }
                sections[i].set_coordinates(std::vector<double>({sectionCoordinates[0] + (currentSectionMaximumRadialHeight + currentSectionMinimumRadialHeight) / 2,
                                                           sections[i].get_coordinates()[1]}));
                sections[i].set_dimensions(std::vector<double>({currentSectionMaximumRadialHeight - currentSectionMinimumRadialHeight,
                                                       sections[i].get_dimensions()[1]}));
            }
        }
        set_sections_description(sections);
    }

    // Angular Delimit
    if (get_layers_description()) {
        auto wirePerWinding = get_wires();
        auto layers = get_layers_description().value();
        if (get_turns_description()) {
            for (size_t i = 0; i < layers.size(); ++i) {
                if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                    auto turnsInLayer = get_turns_by_layer(layers[i].get_name());
                    auto layerCoordinates = layers[i].get_coordinates();
                    auto section = get_section_by_name(layers[i].get_section().value());

                    auto windingIndex = get_winding_index_by_name(turnsInLayer[0].get_winding());
                    double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();

                    double wireRadius;
                    if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR) {
                        wireRadius = windingWindows[0].get_radial_height().value() - turnsInLayer[0].get_coordinates()[0] - wireWidth / 2;
                    }
                    else {
                        wireRadius = windingWindows[0].get_radial_height().value() - turnsInLayer[0].get_coordinates()[0];
                    }

                    double turnDimensionAngle = wound_distance_to_angle(turnsInLayer[0].get_dimensions().value()[1], wireRadius);

                    // for (auto& turn : turnsInLayer) {
                        // turnDimensionAngle = wound_distance_to_angle(turn.get_dimensions().value()[1], windingWindows[0].get_radial_height().value() - turn.get_coordinates()[0]);
                    // }
                    double layerAngle = turnDimensionAngle * turnsInLayer.size();
                    double layerCenterAngle = 0;
                    // double marginAngle0 = 0;
                    // double marginAngle1 = 0;

                    // if (section.get_type() == ElectricalType::CONDUCTION) {
                    //     double lastLayerMaximumRadius = windingWindowsRadius - (section.get_coordinates()[0] + section.get_dimensions()[0] / 2);
                    //     marginAngle0 = wound_distance_to_angle(section.get_margin().value()[0], lastLayerMaximumRadius);
                    //     marginAngle1 = wound_distance_to_angle(section.get_margin().value()[1], lastLayerMaximumRadius);
                    // }

                    switch (layers[i].get_turns_alignment().value()) {
                        case CoilAlignment::INNER_OR_TOP:
                            layerCenterAngle = section.get_coordinates()[1] - section.get_dimensions()[1] / 2 + layerAngle / 2;
                            break;
                        case CoilAlignment::OUTER_OR_BOTTOM:
                            layerCenterAngle = section.get_coordinates()[1] + section.get_dimensions()[1] / 2 - layerAngle / 2;
                            break;
                        case CoilAlignment::CENTERED:
                            layerCenterAngle = section.get_coordinates()[1];
                            break;
                        case CoilAlignment::SPREAD:
                            layerCenterAngle = section.get_coordinates()[1];
                            layerAngle = section.get_dimensions()[1];
                            break;
                        default:
                            throw std::runtime_error("No such section alignment");
                    }
                    layers[i].set_coordinates(std::vector<double>({layers[i].get_coordinates()[0], layerCenterAngle}));
                    layers[i].set_dimensions(std::vector<double>({layers[i].get_dimensions()[0], layerAngle}));
                }
                set_layers_description(layers);
            }
        }

        auto sections = get_sections_description().value();
        for (size_t i = 0; i < sections.size(); ++i) {
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                auto layersInSection = get_layers_by_section(sections[i].get_name());
                if (layersInSection.size() == 0) {
                    throw std::runtime_error("No layers in section: " + sections[i].get_name());
                }
                auto sectionCoordinates = sections[i].get_coordinates();
                double currentSectionMaximumAngle = (layersInSection[0].get_coordinates()[1] - sectionCoordinates[1]) + layersInSection[0].get_dimensions()[1] / 2;
                double currentSectionMinimumAngle = (layersInSection[0].get_coordinates()[1] - sectionCoordinates[1]) - layersInSection[0].get_dimensions()[1] / 2;

                for (auto& layer : layersInSection) {
                    currentSectionMaximumAngle = std::max(currentSectionMaximumAngle, (layer.get_coordinates()[1] - sectionCoordinates[1]) + layer.get_dimensions()[1] / 2);
                    currentSectionMinimumAngle = std::min(currentSectionMinimumAngle, (layer.get_coordinates()[1] - sectionCoordinates[1]) - layer.get_dimensions()[1] / 2);
                }
                sections[i].set_coordinates(std::vector<double>({sections[i].get_coordinates()[0], sectionCoordinates[1] + (currentSectionMaximumAngle + currentSectionMinimumAngle) / 2}));
                sections[i].set_dimensions(std::vector<double>({sections[i].get_dimensions()[0], currentSectionMaximumAngle - currentSectionMinimumAngle}));
            }
        }
        set_sections_description(sections);
    }

     // Angular Compact
    if (get_sections_description()) {
        auto sections = get_sections_description().value();

        std::vector<std::vector<double>> alignedSectionDimensionsPerSection;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            alignedSectionDimensionsPerSection.push_back(get_aligned_section_dimensions_round_window(sectionIndex));
        }

        double currentCoilAngle = alignedSectionDimensionsPerSection[0][1];
        double paddingAmongSectionAngle = alignedSectionDimensionsPerSection[0][3];
        std::vector<Turn> turns;
        if (get_turns_description()) {
            turns = get_turns_description().value();
        }

        std::vector<Layer> layers;
        if (get_layers_description()) {
            layers = get_layers_description().value();
        }

        auto bobbinColumnShape = bobbin.get_processed_description().value().get_column_shape();
        auto bobbinColumnDepth = bobbin.get_processed_description().value().get_column_depth();
        double bobbinColumnWidth;
        if (bobbin.get_processed_description().value().get_column_width()) {
            bobbinColumnWidth = bobbin.get_processed_description().value().get_column_width().value();
        }
        else {
            auto bobbinWindingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
            double bobbinWindingWindowWidth = bobbinWindingWindow.get_width().value();
            double bobbinWindingWindowCenterWidth = bobbinWindingWindow.get_coordinates().value()[0];
            bobbinColumnWidth = bobbinWindingWindowCenterWidth - bobbinWindingWindowWidth / 2;
        }

        auto windingOrientation = get_winding_orientation();

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {

            double marginAngle0 = 0;
            double marginAngle1 = 0;

            if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
                double lastLayerMaximumRadius = windingWindowsRadius - (sections[sectionIndex].get_coordinates()[0] + sections[sectionIndex].get_dimensions()[0] / 2);
                marginAngle0 = wound_distance_to_angle(sections[sectionIndex].get_margin().value()[0], lastLayerMaximumRadius);
                marginAngle1 = wound_distance_to_angle(sections[sectionIndex].get_margin().value()[1], lastLayerMaximumRadius);
            }


            auto sectionAlignment = get_section_alignment();
            if (windingOrientation == WindingOrientation::OVERLAPPING || sectionAlignment == CoilAlignment::SPREAD) {
                currentCoilAngle = alignedSectionDimensionsPerSection[sectionIndex][1];
            }
            else {
                currentCoilAngle += sections[sectionIndex].get_dimensions()[1] / 2 + marginAngle0;
            }

            double compactingShiftAngle = sections[sectionIndex].get_coordinates()[1] - currentCoilAngle;

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
                    compactingShiftAngle = 0;
                }
            }

            sections[sectionIndex].set_coordinates(std::vector<double>({
                sections[sectionIndex].get_coordinates()[0],
                sections[sectionIndex].get_coordinates()[1] - compactingShiftAngle
            }));

            for (size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
                if (layers[layerIndex].get_section().value() == sections[sectionIndex].get_name()){
                    layers[layerIndex].set_coordinates(std::vector<double>({
                        layers[layerIndex].get_coordinates()[0],
                        layers[layerIndex].get_coordinates()[1] - compactingShiftAngle
                    }));
                    size_t turnInThisLayerIndex = 0;
                    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
                        if (turns[turnIndex].get_layer().value() == layers[layerIndex].get_name()){

                            turns[turnIndex].set_coordinates(std::vector<double>({
                                turns[turnIndex].get_coordinates()[0],
                                turns[turnIndex].get_coordinates()[1] - compactingShiftAngle
                            }));


                            if (bobbinColumnShape == ColumnShape::ROUND) {
                                turns[turnIndex].set_length(2 * std::numbers::pi * (turns[turnIndex].get_coordinates()[0] + bobbinColumnWidth));
                                if (turns[turnIndex].get_length() < 0) {
                                    return false;
                                }
                            }
                            else if (bobbinColumnShape == ColumnShape::OBLONG) {
                                turns[turnIndex].set_length(2 * std::numbers::pi * (turns[turnIndex].get_coordinates()[0] + bobbinColumnWidth) + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                                if (turns[turnIndex].get_length() < 0) {
                                    return false;
                                }
                            }
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                double currentTurnCornerRadius = turns[turnIndex].get_coordinates()[0];
                                turns[turnIndex].set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                                if (turns[turnIndex].get_length() < 0) {
                                    return false;
                                }
                            }
                            else {
                                throw std::runtime_error("only round or rectangular columns supported for bobbins");
                            }

                            turnInThisLayerIndex++;
                        }
                    }
                }
                }
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
            }
            else {
                currentCoilAngle += sections[sectionIndex].get_dimensions()[1] / 2 + paddingAmongSectionAngle + marginAngle1;
            }
        }
        if (get_turns_description()) {
            set_turns_description(turns);
        }
        if (get_layers_description()) {
            set_layers_description(layers);
        }
        set_sections_description(sections);

        if (settings->get_coil_include_additional_coordinates()) {
            wind_toroidal_additional_turns();
        }
    }

    if (get_turns_description()) {
        convert_turns_to_cartesian_coordinates();
    }
    return true;
}

std::vector<Wire> Coil::get_wires() {
    std::vector<Wire> wirePerWinding;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        Wire wire = resolve_wire(get_functional_description()[windingIndex]);
        wirePerWinding.push_back(wire);
    }
    return wirePerWinding;
}

Wire Coil::resolve_wire(size_t windingIndex) {
    return resolve_wire(get_functional_description()[windingIndex]);
}

Wire CoilFunctionalDescription::resolve_wire() {
    auto wireOrString = get_wire();
    Wire wire;
    if (std::holds_alternative<std::string>(wireOrString)) {
        try {
            wire = find_wire_by_name(std::get<std::string>(wireOrString));
        }
        catch (const std::exception &exc) {
            // If wire is not found because it is "Dummy", we return a small Round, as it should only happening when get an advised wire
            if (std::get<std::string>(wireOrString) == "Dummy") {
                wire = find_wire_by_name("Round 0.01 - Grade 1");
            }
            else {
                throw std::runtime_error("wire not found: " + std::get<std::string>(wireOrString));
            }
        }
    }
    else {
        wire = std::get<OpenMagnetics::Wire>(wireOrString);
    }
    return wire;
}

Wire Coil::resolve_wire(CoilFunctionalDescription coilFunctionalDescription) {
    return coilFunctionalDescription.resolve_wire();
}

WireType Coil::get_wire_type(CoilFunctionalDescription coilFunctionalDescription) {
    return resolve_wire(coilFunctionalDescription).get_type();
}

WireType Coil::get_wire_type(size_t windingIndex) {
    return get_wire_type(get_functional_description()[windingIndex]);
}

std::string Coil::get_wire_name(CoilFunctionalDescription coilFunctionalDescription) {
    auto name = resolve_wire(coilFunctionalDescription).get_name();
    if (name) {
        return name.value();
    }
    else {
        return "Custom";
    }
}

std::string Coil::get_wire_name(size_t windingIndex) {
    return get_wire_name(get_functional_description()[windingIndex]);
}

Bobbin Coil::resolve_bobbin(Coil coil) { 
    return coil.resolve_bobbin();
}

Bobbin Coil::resolve_bobbin() {
    if (_bobbin_resolved) {
        return _bobbin;
    }

    auto bobbinDataOrNameUnion = get_bobbin();
    if (std::holds_alternative<std::string>(bobbinDataOrNameUnion)) {
        if (std::get<std::string>(bobbinDataOrNameUnion) == "Dummy")
            throw std::runtime_error("Bobbin is dummy");

        auto bobbin = find_bobbin_by_name(std::get<std::string>(bobbinDataOrNameUnion));
        _bobbin = bobbin;
        return bobbin;
    }
    else {
        _bobbin = Bobbin(std::get<Bobbin>(bobbinDataOrNameUnion));
        return _bobbin;
    }
}

size_t Coil::convert_conduction_section_index_to_global(size_t conductionSectionIndex) {
    size_t currentConductionSectionIndex = 0;
    if (!get_sections_description()) {
        throw std::runtime_error("In Convert Conduction Sections: Section description empty, wind coil first");
    }
    auto sections = get_sections_description().value();
    for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
        if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
            if (currentConductionSectionIndex == conductionSectionIndex) {
                return sectionIndex;
            }
            currentConductionSectionIndex++;
        }
    }
    throw std::runtime_error("Index not found");
}

void Coil::clear() {
    set_groups_description(std::nullopt);
    set_sections_description(std::nullopt);
    set_layers_description(std::nullopt);
    set_turns_description(std::nullopt);
}

void Coil::try_rewind() {
    if (!get_sections_description()) {
        return;
    }
    if (!get_layers_description()) {
        return;
    }
 
    if (!get_turns_description()) {
        wind_by_turns();
        delimit_and_compact();
    }
    auto electricalSections = get_sections_by_type(ElectricalType::CONDUCTION);

    if (electricalSections.size() == 1 || get_functional_description().size() == 1) {
        return;
    }

    bool windEvenIfNotFit = settings->get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings->get_coil_delimit_and_compact();

    auto sections = get_sections_description().value();
    std::vector<double> extraSpaceNeededPerSection;
    double totalExtraSpaceNeeded = 0;
    auto bobbin = resolve_bobbin();
    auto sectionOrientation = bobbin.get_winding_window_sections_orientation(0);
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    double windingWindowRemainingRestrictiveDimension;
    double windingWindowRemainingRestrictiveDimensionAccordingToSections;
    double windingWindowRestrictiveDimension;
    if (sectionOrientation == WindingOrientation::OVERLAPPING) {
        windingWindowRemainingRestrictiveDimensionAccordingToSections = windingWindowDimensions[0];
        windingWindowRemainingRestrictiveDimension = windingWindowDimensions[0];
        windingWindowRestrictiveDimension = windingWindowDimensions[0];
    }
    else {
        windingWindowRemainingRestrictiveDimensionAccordingToSections = windingWindowDimensions[1];
        windingWindowRemainingRestrictiveDimension = windingWindowDimensions[1];
        windingWindowRestrictiveDimension = windingWindowDimensions[1];
    }


    for (auto& section : sections) {
        if (section.get_type() == ElectricalType::INSULATION) {
            if (sectionOrientation == WindingOrientation::OVERLAPPING) {
                windingWindowRestrictiveDimension -= section.get_dimensions()[0];
            }
            else {
                windingWindowRestrictiveDimension -= section.get_dimensions()[1];
            }
        }
    }

    for (auto& section : sections) {
        double sectionRestrictiveDimension;
        double layersRestrictiveDimension = 0;
        double sectionFillingFactor;
        double extraSpaceNeededThisSection = 0;

        auto layers = get_layers_by_section(section.get_name());
        if (sectionOrientation == WindingOrientation::OVERLAPPING) {
            if (section.get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                for (auto layer : layers) {
                    double layerRestrictiveDimension = layer.get_dimensions()[0];
                    double layerFillingFactor = layer.get_filling_factor().value();
                    layersRestrictiveDimension += layerRestrictiveDimension;

                    extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                    windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
                }
            }
            if (section.get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                double layerRestrictiveDimension = 0;
                double layerFillingFactor = 0;
                for (auto layer : layers) {
                    layerRestrictiveDimension = std::max(layerRestrictiveDimension, layer.get_dimensions()[0]);
                    layerFillingFactor = std::max(layerFillingFactor, layer.get_filling_factor().value());
                }
                layersRestrictiveDimension = layerRestrictiveDimension;
                extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
            }
        }
        else if (sectionOrientation == WindingOrientation::CONTIGUOUS) {
            if (section.get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                double layerRestrictiveDimension = 0;
                double layerFillingFactor = 0;
                for (auto layer : layers) {
                    layerRestrictiveDimension = std::max(layerRestrictiveDimension, layer.get_dimensions()[1]);
                    layerFillingFactor = std::max(layerFillingFactor, layer.get_filling_factor().value());
                }

                layersRestrictiveDimension = layerRestrictiveDimension;
                extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
            }
            if (section.get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                for (auto layer : layers) {
                    double layerRestrictiveDimension = layer.get_dimensions()[1];
                    double layerFillingFactor = layer.get_filling_factor().value();
                    layersRestrictiveDimension += layerRestrictiveDimension;

                    extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                    windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
                }
            }
        }

        if (sectionOrientation == WindingOrientation::OVERLAPPING) {
            sectionRestrictiveDimension = section.get_dimensions()[0];
            sectionFillingFactor = overlapping_filling_factor(section);
        }
        else {
            sectionRestrictiveDimension = section.get_dimensions()[1];
            sectionFillingFactor = contiguous_filling_factor(section);
        }
        windingWindowRemainingRestrictiveDimensionAccordingToSections -= sectionRestrictiveDimension;


        extraSpaceNeededThisSection = std::max(extraSpaceNeededThisSection, (sectionFillingFactor - 1) * sectionRestrictiveDimension);
        if (extraSpaceNeededThisSection < 0 || std::isnan(extraSpaceNeededThisSection)) {
            throw std::runtime_error("extraSpaceNeededThisSection cannot be negative or nan: " + std::to_string(extraSpaceNeededThisSection));
        }
        if (section.get_type() == ElectricalType::CONDUCTION) {
        }
        extraSpaceNeededPerSection.push_back(extraSpaceNeededThisSection);
        totalExtraSpaceNeeded += extraSpaceNeededThisSection;
    }

    if (windingWindowRemainingRestrictiveDimensionAccordingToSections <= 0 || totalExtraSpaceNeeded <= 0) {
        return;
    }

    std::vector<double> newProportions;
    double numberWindings = get_functional_description().size();

    if (totalExtraSpaceNeeded < 0 || std::isnan(totalExtraSpaceNeeded)) {
        throw std::runtime_error("totalExtraSpaceNeeded cannot be negative or nan: " + std::to_string(totalExtraSpaceNeeded));
    }

    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        // double currentProportion = _currentProportionPerWinding[windingIndex];
        double currentSpace = 0;
        double extraSpaceNeededThisWinding = 0;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            for (auto & winding : sections[sectionIndex].get_partial_windings()) {
                if (winding.get_winding() == get_functional_description()[windingIndex].get_name()) {
                    if (sectionOrientation == WindingOrientation::OVERLAPPING) {
                        currentSpace += sections[sectionIndex].get_dimensions()[0];

                        // We need to add half the insulation space after it, in case there is
                        if (sectionIndex < sections.size() - 1) {
                            if (sections[sectionIndex + 1].get_type() == ElectricalType::INSULATION) {
                                // throw std::runtime_error("Consecutive layer to CONDUCTION must always be INSULATION");
                                if (sectionIndex == 0) {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[0] / 2;
                                }
                                else if (sectionIndex == sections.size() - 2) {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[0] * 3 / 2;
                                }
                                else {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[0];
                                }
                            }
                        }
                    }
                    else {
                        currentSpace += sections[sectionIndex].get_dimensions()[1];

                        // We need to add half the insulation space after it, in case there is
                        if (sectionIndex < sections.size() - 1) {
                            if (sections[sectionIndex + 1].get_type() != ElectricalType::INSULATION) {
                                // throw std::runtime_error("Consecutive layer to CONDUCTION must always be INSULATION");
                                if (sectionIndex == 0 || sectionIndex == sections.size() - 2) {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[1] / 2;
                                }
                                else {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[1];
                                }
                            }
                        }
                    }

                    extraSpaceNeededThisWinding += extraSpaceNeededPerSection[sectionIndex];
                    continue;
                }
            }
        }
        if (extraSpaceNeededThisWinding < 0 || std::isnan(extraSpaceNeededThisWinding)) {
            throw std::runtime_error("extraSpaceNeededThisWinding cannot be negative or nan: " + std::to_string(extraSpaceNeededThisWinding));
        }
        // double proportionOfNeededForThisWinding = extraSpaceNeededThisWinding / totalExtraSpaceNeeded;
        double extraSpaceGottenByThisWinding = windingWindowRemainingRestrictiveDimensionAccordingToSections * extraSpaceNeededThisWinding / totalExtraSpaceNeeded;
        double newSpaceGottenByThisWinding = currentSpace + extraSpaceGottenByThisWinding;
        double newProportionGottenByThisWinding = newSpaceGottenByThisWinding / windingWindowRestrictiveDimension;

        if (extraSpaceGottenByThisWinding < 0 || std::isnan(extraSpaceGottenByThisWinding)) {
            throw std::runtime_error("extraSpaceGottenByThisWinding cannot be negative or nan: " + std::to_string(extraSpaceGottenByThisWinding));
        }
        if (newProportionGottenByThisWinding < 0 || std::isnan(newProportionGottenByThisWinding)) {
            throw std::runtime_error("newProportionGottenByThisWinding cannot be negative or nan: " + std::to_string(newProportionGottenByThisWinding));
        }
        if (roundFloat(newProportionGottenByThisWinding, 6) > 1 || std::isnan(newProportionGottenByThisWinding)) {
            throw std::runtime_error("newProportionGottenByThisWinding cannot be greater than 1 or nan: " + std::to_string(newProportionGottenByThisWinding));
        }

        newProportions.push_back(newProportionGottenByThisWinding);
    }

    wind_by_sections(newProportions, _currentPattern, _currentRepetitions);



    wind_by_layers();

    if (!get_layers_description()) {
        return;
    }
    // set_turns_description(std::nullopt);
    if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
        wind_by_turns();
        if (delimitAndCompact) {
            delimit_and_compact();
        }
    }
}

void Coil::preload_margins(std::vector<std::vector<double>> marginPairs) {
    for (auto margins : marginPairs) {
        _marginsPerSection.push_back(margins);
        // Add an extra one for the insulation layer
        _marginsPerSection.push_back(margins);
    }
}

void Coil::add_margin_to_section_by_index(size_t sectionIndex, std::vector<double> margins) {
    if (!get_sections_description()) {
        throw std::runtime_error("In Add Margin to Section: Section description empty, wind coil first");
    }
    if (margins.size() != 2) {
        throw std::runtime_error("Margin vector must have two elements");
    }
    auto sections = get_sections_description().value();
    _marginsPerSection[convert_conduction_section_index_to_global(sectionIndex)] = margins;
    sections[convert_conduction_section_index_to_global(sectionIndex)].set_margin(margins);

    set_sections_description(sections);


    bool windEvenIfNotFit = settings->get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings->get_coil_delimit_and_compact();
    bool tryRewind = settings->get_coil_try_rewind();

    wind_by_sections();
    wind_by_layers();
    if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
        wind_by_turns();
        if (delimitAndCompact) {
            delimit_and_compact();
        }
    }
    if (tryRewind && !are_sections_and_layers_fitting()) {
        try_rewind();
    }
}

std::vector<Section> Coil::get_sections_description_conduction() {
    std::vector<Section> sectionsConduction;
    if (!get_sections_description()) {
        throw std::runtime_error("Not wound by sections");
    }
    std::vector<Section> sections = get_sections_description().value();
    for (auto section : sections) {
        if (section.get_type() == ElectricalType::CONDUCTION) {
            sectionsConduction.push_back(section);
        }
    }

    return sectionsConduction;
}

std::vector<Layer> Coil::get_layers_description_conduction() {
    std::vector<Layer> layersConduction;
    if (!get_layers_description()) {
        throw std::runtime_error("Not wound by layers");
    }
    std::vector<Layer> layers = get_layers_description().value();
    for (auto layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            layersConduction.push_back(layer);
        }
    }

    return layersConduction;
}

std::vector<Section> Coil::get_sections_description_insulation() {
    std::vector<Section> sectionsInsulation;
    if (!get_sections_description()) {
        throw std::runtime_error("Not wound by sections");
    }
    std::vector<Section> sections = get_sections_description().value();
    for (auto section : sections) {
        if (section.get_type() == ElectricalType::INSULATION) {
            sectionsInsulation.push_back(section);
        }
    }

    return sectionsInsulation;
}

std::vector<Layer> Coil::get_layers_description_insulation() {
    std::vector<Layer> layersInsulation;
    if (!get_layers_description()) {
        throw std::runtime_error("Not wound by layers");
    }
    std::vector<Layer> layers = get_layers_description().value();
    for (auto layer : layers) {
        if (layer.get_type() == ElectricalType::INSULATION) {
            layersInsulation.push_back(layer);
        }
    }

    return layersInsulation;
}

double Coil::calculate_external_proportion_for_wires_in_toroidal_cores(Core core, Coil coil) {
    CoreShape shape = std::get<CoreShape>(core.get_functional_description().get_shape());
    auto processedDescription = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});

    double coreWidth = processedDescription.get_width();

    if (!coil.get_turns_description()) {
        return 1;
    }

    auto turns = coil.get_turns_description().value();
    double maximumAdditionalRadialCoordinate = 0;
    for (size_t i = 0; i < turns.size(); ++i){
        if (turns[i].get_additional_coordinates()) {
            auto additionalCoordinates = turns[i].get_additional_coordinates().value();
            for (auto additionalCoordinate : additionalCoordinates){
                maximumAdditionalRadialCoordinate = std::max(maximumAdditionalRadialCoordinate, hypot(additionalCoordinate[0], additionalCoordinate[1]) + turns[i].get_dimensions().value()[0] / 2);
            }
        }
    }
    auto bobbin = coil.resolve_bobbin();

    auto sectionsOrientation = bobbin.get_winding_window_sections_orientation();

    if (maximumAdditionalRadialCoordinate > 0 && sectionsOrientation == WindingOrientation::OVERLAPPING) {
        auto sections = coil.get_sections_by_type(ElectricalType::INSULATION);
        for (auto section : sections){
            maximumAdditionalRadialCoordinate += section.get_dimensions()[0];
        }
    }

    if (maximumAdditionalRadialCoordinate == 0) {
        return 1;
    }

    return (2 * maximumAdditionalRadialCoordinate) / coreWidth;
}


double Coil::get_insulation_section_thickness(std::string sectionName) {
    return get_insulation_section_thickness(*this, sectionName);
}

double Coil::get_insulation_section_thickness(Coil coil, std::string sectionName) {
    if (!coil.get_sections_description()) {
        throw std::runtime_error("Coil is missing sections description");
    }
    if (!coil.get_layers_description()) {
        throw std::runtime_error("Coil is missing layers description");
    }

    auto layers = coil.get_layers_by_section(sectionName);

    double thickness = 0;

    for (auto layer : layers) {
        thickness += coil.get_insulation_layer_thickness(layer);
    }

    return thickness;
}

double Coil::get_insulation_layer_thickness(Coil coil, std::string layerName) {
    return coil.get_insulation_layer_thickness(layerName);
}

double Coil::get_insulation_layer_thickness(std::string layerName) {
    if (!get_layers_description()) {
        throw std::runtime_error("Coil is missing layers description");
    }
    auto layer = get_layer_by_name(layerName);
    return get_insulation_layer_thickness(layer);
}

double Coil::get_insulation_layer_thickness(Layer layer) {
    if (!layer.get_coordinate_system()) {
        layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
    }
    if (layer.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        if (layer.get_orientation() == WindingOrientation::CONTIGUOUS) {
            return layer.get_dimensions()[1];
        }
        else {
            return layer.get_dimensions()[0];
        }
    }
    else {
        if (layer.get_orientation() == WindingOrientation::CONTIGUOUS) {
            auto bobbin = resolve_bobbin();
            auto bobbinProcessedDescription = bobbin.get_processed_description().value();
            auto windingWindows = bobbinProcessedDescription.get_winding_windows();

            double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
            double layerRadialHeight = layer.get_dimensions()[0];
            double radius = windingWindowRadialHeight - layerRadialHeight;
            double layerAngle = layer.get_dimensions()[1];
            double layerPerimeter = std::numbers::pi * (layerAngle / 180) * radius;
            return layerPerimeter;
        }
        else {
            return layer.get_dimensions()[0];
        }
    }
}

InsulationMaterial Coil::resolve_insulation_layer_insulation_material(std::string layerName){
    auto layer = get_layer_by_name(layerName);
    return resolve_insulation_layer_insulation_material(layer);
}

InsulationMaterial Coil::resolve_insulation_layer_insulation_material(Coil coil, std::string layerName) {
    auto layer = coil.get_layer_by_name(layerName);
    return coil.resolve_insulation_layer_insulation_material(layer);
}

InsulationMaterial Coil::resolve_insulation_layer_insulation_material(Layer layer) {
    if (!layer.get_insulation_material()) {
        layer.set_insulation_material(defaults.defaultLayerInsulationMaterial);
        // throw std::runtime_error("Layer is missing material information");
    }

    auto insulationMaterial = layer.get_insulation_material().value();
    // If the material is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(insulationMaterial)) {
        auto insulationMaterialData = find_insulation_material_by_name(std::get<std::string>(insulationMaterial));

        return insulationMaterialData;
    }
    else {
        return InsulationMaterial(std::get<MAS::InsulationMaterial>(insulationMaterial));
    }
}

double Coil::get_insulation_layer_relative_permittivity(std::string layerName) {
    auto layer = get_layer_by_name(layerName);
    return get_insulation_layer_relative_permittivity(layer);
}
double Coil::get_insulation_layer_relative_permittivity(Coil coil, std::string layerName) {
    return coil.get_insulation_layer_relative_permittivity(layerName);
}
double Coil::get_insulation_layer_relative_permittivity(Layer layer) {
    auto coatingInsulationMaterial = resolve_insulation_layer_insulation_material(layer);
    if (!coatingInsulationMaterial.get_relative_permittivity())
        throw std::runtime_error("Coating insulation material is missing dielectric constant");
    return coatingInsulationMaterial.get_relative_permittivity().value();
}

double Coil::get_insulation_section_relative_permittivity(std::string sectionName) {
    return get_insulation_section_relative_permittivity(*this, sectionName);
}
double Coil::get_insulation_section_relative_permittivity(Coil coil, std::string sectionName) {
    auto layers = coil.get_layers_by_section(sectionName);
    if (layers.size() == 0)
        throw std::runtime_error("No layers in this section");

    double averagerelativePermittivity = 0;
    for (auto layer : layers) {
        averagerelativePermittivity += coil.get_insulation_layer_relative_permittivity(layer);
    }
    return averagerelativePermittivity / layers.size();
}

std::vector<double> Coil::get_turns_ratios() {
    std::vector<double>  turnsRatios;
    for (size_t windingIndex = 1; windingIndex < get_functional_description().size(); ++windingIndex) {
        turnsRatios.push_back(double(get_functional_description()[0].get_number_turns()) / get_functional_description()[windingIndex].get_number_turns());
    }

    return turnsRatios;
}

std::vector<double> Coil::get_maximum_dimensions() {
    std::vector<double> bobbinMaximumDimensions = resolve_bobbin().get_maximum_dimensions();

    if (!get_turns_description()) {
        throw std::runtime_error("Missing turns");
    }
    auto turns = get_turns_description().value();

    double width = 0;
    double height = 0;
    double depth = 0;

    for (auto turn : turns) {
        double turnMaxWidthPosition = 0;
        double turnMaxHeightPosition = 0;
        if (turn.get_additional_coordinates()) {
            auto additionalCoordinates = turn.get_additional_coordinates().value();
            turnMaxWidthPosition = fabs(additionalCoordinates[0][0]) + turn.get_dimensions().value()[0] / 2;
            turnMaxHeightPosition = fabs(additionalCoordinates[0][1]) + turn.get_dimensions().value()[1] / 2;
        }
        else {
            turnMaxWidthPosition = fabs(turn.get_coordinates()[0]) + turn.get_dimensions().value()[0] / 2;
            turnMaxHeightPosition = fabs(turn.get_coordinates()[1]) + turn.get_dimensions().value()[1] / 2;
        }

        width = std::max(width, turnMaxWidthPosition);
        height = std::max(height, turnMaxHeightPosition);
    }

    double bobbinExtraDepthDimension = bobbinMaximumDimensions[0] - bobbinMaximumDimensions[2];
    depth = width + bobbinExtraDepthDimension;

    width = std::max(width, bobbinMaximumDimensions[0]);
    height = std::max(height, bobbinMaximumDimensions[1]);
    depth = std::max(depth, bobbinMaximumDimensions[2]);

    return std::vector<double>{width, height, depth};
}


std::vector<std::vector<size_t>> Coil::get_patterns(Inputs& inputs, CoreType coreType) {
    auto isolationSidesRequired = inputs.get_isolation_sides_used();

    if (!inputs.get_design_requirements().get_isolation_sides()) {
        throw std::runtime_error("Missing isolation sides requirement");
    }

    auto isolationSidesRequirement = inputs.get_design_requirements().get_isolation_sides().value();

    std::vector<std::vector<size_t>> sectionPatterns;
    for(size_t i = 0; i < tgamma(isolationSidesRequired.size() + 1) / 2; ++i) {
        std::vector<size_t> sectionPattern;
        for (auto isolationSide : isolationSidesRequired) {
            for (size_t windingIndex = 0; windingIndex < inputs.get_mutable_design_requirements().get_turns_ratios().size() + 1; ++windingIndex) {
                if (isolationSidesRequirement[windingIndex] == isolationSide) {
                    sectionPattern.push_back(windingIndex);
                }
            }
        }
        sectionPatterns.push_back(sectionPattern);
        if (sectionPatterns.size() > defaults.maximumCoilPattern) {
            break;
        }

        std::next_permutation(isolationSidesRequired.begin(), isolationSidesRequired.end());
    }

    if (coreType == CoreType::TOROIDAL) {
        // We remove the last combination as in toroids they go around
        size_t elementsToKeep = std::max(size_t(1), isolationSidesRequired.size() - 1);
        sectionPatterns = std::vector<std::vector<size_t>>(sectionPatterns.begin(), sectionPatterns.end() - (sectionPatterns.size() - elementsToKeep));
    }

    return sectionPatterns;
}


std::vector<size_t> Coil::get_repetitions(Inputs& inputs, CoreType coreType) {
    if (inputs.get_design_requirements().get_turns_ratios().size() == 0 || coreType == CoreType::TOROIDAL) {
        return {1};  // hardcoded
    }
    if (inputs.get_design_requirements().get_wiring_technology()) {
        if (inputs.get_design_requirements().get_wiring_technology().value() == WiringTechnology::PRINTED) {
            std::vector<size_t> repetitions;
            for (size_t repetition = 1; repetition <= (settings->get_coil_maximum_layers_planar() / (inputs.get_design_requirements().get_turns_ratios().size() + 1)); ++repetition) {
                repetitions.push_back(repetition);
            }
            return repetitions;
        }
    }
    if (inputs.get_design_requirements().get_leakage_inductance()) {
        return {2, 1};  // hardcoded
    }
    else{
        return {1, 2};  // hardcoded
    }
}

std::pair<std::vector<size_t>, size_t> Coil::check_pattern_and_repetitions_integrity(std::vector<size_t> pattern, size_t repetitions) {
    bool needsMerge = false;
    for (auto winding : get_functional_description()) {
        // TODO expand for more than one winding per layer
        size_t numberPhysicalTurns = winding.get_number_turns() * winding.get_number_parallels();
        if (numberPhysicalTurns < repetitions) {
            needsMerge = true;
        }
    }

    std::vector<size_t> newPattern;
    if (needsMerge) {
        for (size_t repetition = 1; repetition <= repetitions; ++repetition) {
            for (auto windingIndex : pattern) {
                auto winding = get_functional_description()[windingIndex];
                size_t numberPhysicalTurns = winding.get_number_turns() * winding.get_number_parallels();
                if (numberPhysicalTurns >= repetition) {
                    newPattern.push_back(windingIndex);
                }
            }
        }
        return {newPattern, 1};
    }
    return {pattern, repetitions};
}

bool Coil::is_edge_wound_coil() {
    auto wires = get_wires();
    for (auto wire : wires) {
        if (wire.get_type() != WireType::RECTANGULAR) {
            return false;
        }
    }

    return true;
}


} // namespace OpenMagnetics
 
