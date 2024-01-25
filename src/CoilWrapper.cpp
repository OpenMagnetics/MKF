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
#include "Constants.h"
#include "Utils.h"
#include "CoilWrapper.h"
#include "json.hpp"
#include "InsulationMaterialWrapper.h"
#include "Settings.h"

#include <magic_enum.hpp>

using json = nlohmann::json;

namespace OpenMagnetics {


CoilWrapper::CoilWrapper(const json& j, uint8_t interleavingLevel,
                               WindingOrientation windingOrientation,
                               WindingOrientation layersOrientation,
                               CoilAlignment turnsAlignment,
                               CoilAlignment sectionAlignment,
                               bool delimitAndCompact) {
    _interleavingLevel = interleavingLevel;
    _windingOrientation = windingOrientation;
    _layersOrientation = layersOrientation;
    _turnsAlignment = turnsAlignment;
    _sectionAlignment = sectionAlignment;
    from_json(j, *this);

    if (wind() && delimitAndCompact) {
        delimit_and_compact();
    }

}

CoilWrapper::CoilWrapper(const Coil coil, bool delimitAndCompact) {
    bool hasSectionsData = false;
    bool hasLayersData = false;
    bool hasTurnsData = false;

    set_functional_description(coil.get_functional_description());
    set_bobbin(coil.get_bobbin());

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

    if (!hasSectionsData || !hasLayersData || (!hasTurnsData && are_sections_and_layers_fitting())) {
        if (wind() && delimitAndCompact) {
            delimit_and_compact();
        }
    }

}

bool CoilWrapper::try_wind(bool delimitAndCompact) {
    auto settings = Settings::GetInstance();
    bool windEvenIfNotFit = settings->get_wind_even_if_not_fit();
    bool hasSectionsData = false;
    bool hasLayersData = false;
    bool hasTurnsData = false;

    if (get_sections_description()) {
        hasSectionsData = true;
    }
    if (get_layers_description()) {
        hasLayersData = true;
    }
    if (get_turns_description()) {
        hasTurnsData = true;
    }

    if (!hasSectionsData || !hasLayersData || (!hasTurnsData && (windEvenIfNotFit || are_sections_and_layers_fitting()))) {
        if (wind() && delimitAndCompact) {
            return delimit_and_compact();
        }
    }
    return false;
}

bool CoilWrapper::wind() {
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    return wind(proportionPerWinding, pattern, _interleavingLevel);
}

bool CoilWrapper::wind(std::vector<size_t> pattern, size_t repetitions){
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind(proportionPerWinding, pattern, repetitions);
}

bool CoilWrapper::wind(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    auto settings = Settings::GetInstance();
    bool windEvenIfNotFit = settings->get_wind_even_if_not_fit();
    std::string bobbinName = "";
    if (std::holds_alternative<std::string>(get_bobbin())) {
        bobbinName = std::get<std::string>(get_bobbin());
        if (bobbinName != "Dummy") {
            auto bobbinData = OpenMagnetics::find_bobbin_by_name(std::get<std::string>(get_bobbin()));
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
            if (_inputs) {
                calculate_insulation();
            }
            else {
                calculate_mechanical_insulation();
            }
            wind_by_sections(proportionPerWinding, pattern, repetitions);
            wind_by_layers(); 
    // if (!are_sections_and_layers_fitting()) {
    //     try_rewind();
    // }
            if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
                wind_by_turns();
                return true;
            }

        }
    }
    return false;
}


std::vector<WindingStyle> CoilWrapper::wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, std::vector<size_t> numberSlots) {
    std::vector<WindingStyle> windByConsecutiveTurns;
    for (size_t i = 0; i < numberTurns.size(); ++i) {
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

WindingStyle CoilWrapper::wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, uint8_t numberSlots) {
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

std::vector<uint64_t> CoilWrapper::get_number_turns() {
    std::vector<uint64_t> numberTurns;
    for (auto & winding : get_functional_description()) {
        numberTurns.push_back(winding.get_number_turns());
    }
    return numberTurns;
}

void CoilWrapper::set_number_turns(std::vector<uint64_t> numberTurns) {
    for (size_t i=0; i< get_functional_description().size(); ++i) {
        get_mutable_functional_description()[i].set_number_turns(numberTurns[i]);
    }
}


std::vector<Layer> CoilWrapper::get_layers_by_section(std::string sectionName) {
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

std::vector<Turn> CoilWrapper::get_turns_by_layer(std::string layerName) {
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

const std::vector<Section> CoilWrapper::get_sections_by_type(ElectricalType electricalType) const {
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

const std::vector<Section> CoilWrapper::get_sections_by_winding(std::string windingName) const {
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

const Section CoilWrapper::get_section_by_name(std::string name) const {
    auto sections = get_sections_description().value();
    for (auto & section : sections) {
        if (section.get_name() == name) {
            return section;
        }
    }
    throw std::runtime_error("Not found section with name:" + name);
}

const std::vector<Layer> CoilWrapper::get_layers_by_type(ElectricalType electricalType) const {
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

std::vector<uint64_t> CoilWrapper::get_number_parallels() {
    std::vector<uint64_t> numberParallels;
    for (auto & winding : get_functional_description()) {
        numberParallels.push_back(winding.get_number_parallels());
    }
    return numberParallels;
}

void CoilWrapper::set_number_parallels(std::vector<uint64_t> numberParallels){
    for (size_t i=0; i< get_functional_description().size(); ++i) {
        get_mutable_functional_description()[i].set_number_parallels(numberParallels[i]);
    }
}

CoilFunctionalDescription CoilWrapper::get_winding_by_name(std::string name) {
    for (auto& CoilFunctionalDescription : get_functional_description()) {
        if (CoilFunctionalDescription.get_name() == name) {
            return CoilFunctionalDescription;
        }
    }
    throw std::runtime_error("No such a winding name: " + name);
}

size_t CoilWrapper::get_winding_index_by_name(std::string name) {
    for (size_t i=0; i<get_functional_description().size(); ++i) {
        if (get_functional_description()[i].get_name() == name) {
            return i;
        }
    }
    throw std::runtime_error("No such a winding name: " + name);
}

bool CoilWrapper::are_sections_and_layers_fitting() {
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
        if (section.get_filling_factor().value() > 1 || horizontal_filling_factor(section) > 1 || vertical_filling_factor(section) > 1) {
            windTurns = false;
        }
    }
    for (auto& layer: layers) {
        if (layer.get_filling_factor().value() > 1) {
            windTurns = false;
        }
    }
    return windTurns;
}

double CoilWrapper::horizontal_filling_factor(Section section) {
    auto layers = get_layers_by_section(section.get_name());
    double sectionWidth = section.get_dimensions()[0];
    double layersWidth = 0;
    for (auto& layer : layers) {
        if (layer.get_orientation() == WindingOrientation::VERTICAL) {
            layersWidth += layer.get_dimensions()[0];
        }
        else {
            layersWidth = std::max(layersWidth, layer.get_dimensions()[0]);
        }
    }
    return layersWidth / sectionWidth;
}

double CoilWrapper::vertical_filling_factor(Section section) {
    auto layers = get_layers_by_section(section.get_name());
    double sectionHeight = section.get_dimensions()[1];
    double layersHeight = 0;
    for (auto& layer : layers) {
        if (layer.get_orientation() == WindingOrientation::VERTICAL) {
            layersHeight = std::max(layersHeight, layer.get_dimensions()[1]);
        }
        else {
            layersHeight += layer.get_dimensions()[1];
        }
    }
    return layersHeight / sectionHeight;
}

double CoilWrapper::horizontal_filling_factor(Layer layer) {
    auto turns = get_turns_by_layer(layer.get_name());
    double layerWidth = layer.get_dimensions()[0];
    double turnsWidth = 0;
    for (auto& turn : turns) {
        turnsWidth += turn.get_dimensions().value()[0];
    }
    return turnsWidth / layerWidth;
}

double CoilWrapper::vertical_filling_factor(Layer layer) {
    auto turns = get_turns_by_layer(layer.get_name());
    double layerHeight = layer.get_dimensions()[1];
    double turnsHeight = 0;
    for (auto& turn : turns) {
        turnsHeight += turn.get_dimensions().value()[1];
    }
    return turnsHeight / layerHeight;
}

std::pair<uint64_t, std::vector<double>> get_parallels_proportions(size_t slotIndex, uint8_t slots, uint64_t numberTurns, uint64_t numberParallels, 
                                                                   std::vector<double> remainingParallelsProportion, WindingStyle windByConsecutiveTurns,
                                                                   std::vector<double> totalParallelsProportion) {
    uint64_t physicalTurnsThisSlot = 0;
    std::vector<double> slotParallelsProportion(numberParallels, 0);
    if (windByConsecutiveTurns == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
        size_t remainingPhysicalTurns = 0;
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            remainingPhysicalTurns += round(remainingParallelsProportion[parallelIndex] * numberTurns);
        }
        physicalTurnsThisSlot = std::min(uint64_t(remainingPhysicalTurns), uint64_t(ceil(double(remainingPhysicalTurns) / (slots - slotIndex))));
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
            double numberTurnsToAddToCurrentParallel = ceil(numberTurns * totalParallelsProportion[parallelIndex] / slots);
            double remainingTurnsBeforeThisParallel = numberTurns * remainingParallelsProportion[parallelIndex];
            double remainingTurnsAfterThisParallel = remainingTurnsBeforeThisParallel - numberTurnsToAddToCurrentParallel;
            double remainingSlots = slots - slotIndex;
            double remainingSlotsAfterThisOne = remainingSlots - 1;
            if (remainingTurnsAfterThisParallel < remainingSlotsAfterThisOne) {
                numberTurnsToAddToCurrentParallel = roundFloat(remainingTurnsBeforeThisParallel / remainingSlots, 10);
            }

            double proportionParallelsThisSection = std::min(remainingParallelsProportion[parallelIndex], numberTurnsToAddToCurrentParallel / numberTurns);
            physicalTurnsThisSlot += numberTurnsToAddToCurrentParallel;
            slotParallelsProportion[parallelIndex] = proportionParallelsThisSection;
        }
    }

    return {physicalTurnsThisSlot, slotParallelsProportion};
}

double get_area_used_in_wires(WireWrapper wire, uint64_t physicalTurns) {
    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        double wireDiameter = resolve_dimensional_values(wire.get_outer_diameter().value());
        return physicalTurns * pow(wireDiameter, 2);
    }
    else {
        double wireWidth = resolve_dimensional_values(wire.get_outer_width().value());
        double wireHeight = resolve_dimensional_values(wire.get_outer_height().value());
        return physicalTurns * wireWidth * wireHeight;
    }
}


bool CoilWrapper::calculate_mechanical_insulation() {
    // Insulation layers just for mechanical reasons, one layer between sections at least
    auto wirePerWinding = get_wires();

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();

    auto layersOrientation = _layersOrientation;

    // TODO: Properly think about insualtion layers with weird windings
    if (_windingOrientation == WindingOrientation::VERTICAL && _layersOrientation == WindingOrientation::VERTICAL) {
        layersOrientation = WindingOrientation::HORIZONTAL;
    }
    if (_windingOrientation == WindingOrientation::HORIZONTAL && _layersOrientation == WindingOrientation::HORIZONTAL) {
        layersOrientation = WindingOrientation::VERTICAL;
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
            InsulationMaterialWrapper defaultInsulationMaterial = find_insulation_material_by_name(Defaults().defaultInsulationMaterial);
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
                layer.set_turns_alignment(_turnsAlignment);
                if (layersOrientation == WindingOrientation::VERTICAL) {
                    layer.set_dimensions(std::vector<double>{defaultInsulationMaterial.get_thinner_tape_thickness(), windingWindowHeight});
                }
                else if (layersOrientation == WindingOrientation::HORIZONTAL) {
                    layer.set_dimensions(std::vector<double>{windingWindowWidth, defaultInsulationMaterial.get_thinner_tape_thickness()});
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
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowHeight});
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                section.set_dimensions(std::vector<double>{windingWindowWidth, coilSectionInterface.get_solid_insulation_thickness()});
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

bool CoilWrapper::calculate_insulation() {
    auto wirePerWinding = get_wires();

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            if (leftTopWindingIndex == rightBottomWindingIndex) {
                continue;
            }
            auto inputs = _inputs.value();
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            CoilSectionInterface coilSectionInterface;
            coilSectionInterface.set_solid_insulation_thickness(DBL_MAX);
            coilSectionInterface.set_number_layers_insulation(ULONG_MAX);
            InsulationMaterialWrapper chosenInsulationMaterial;

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

            _insulationLayers[windingsMapKey] = std::vector<Layer>();
            _coilSectionInterfaces[windingsMapKey] = coilSectionInterface;

            for (size_t layerIndex = 0; layerIndex < coilSectionInterface.get_number_layers_insulation(); ++layerIndex) {
                Layer layer;
                layer.set_partial_windings(std::vector<PartialWinding>{});
                // layer.set_section(section.get_name());
                layer.set_type(ElectricalType::INSULATION);
                layer.set_name("temp");
                layer.set_orientation(_layersOrientation);
                layer.set_turns_alignment(_turnsAlignment);
                if (_layersOrientation == WindingOrientation::VERTICAL) {
                    layer.set_dimensions(std::vector<double>{chosenInsulationMaterial.get_thicker_tape_thickness(), windingWindowHeight});
                }
                else if (_layersOrientation == WindingOrientation::HORIZONTAL) {
                    layer.set_dimensions(std::vector<double>{windingWindowWidth, chosenInsulationMaterial.get_thicker_tape_thickness()});
                }
                // layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layer.set_filling_factor(1);
                _insulationLayers[windingsMapKey].push_back(layer);
            }
            // _insulationLayersLog[windingsMapKey] = "Adding " + std::to_string(coilSectionInterface.get_number_layers_insulation()) + " insulation layers, as we need a thickness of " + std::to_string(smallestInsulationThicknessCoveringRemaining * 1000) + " mm to achieve " + neededInsulationTypeString + " insulation";

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(_layersOrientation);
            section.set_type(ElectricalType::INSULATION);
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowHeight});
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                section.set_dimensions(std::vector<double>{windingWindowWidth, coilSectionInterface.get_solid_insulation_thickness()});
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

std::vector<std::pair<size_t, double>> CoilWrapper::get_ordered_sections(double spaceForSections, std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
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
            double spaceForSection = roundFloat(spaceForSections * proportionPerWinding[windingIndex] / numberSectionsPerWinding[windingIndex], 9);
            orderedSections.push_back({windingIndex, spaceForSection});
        }
    }

    return orderedSections;
}

std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> CoilWrapper::add_insulation_to_sections(std::vector<std::pair<size_t, double>> orderedSections){
    std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation;
    for (size_t sectionIndex = 1; sectionIndex < orderedSections.size(); ++sectionIndex) {
        auto leftWindingIndex = orderedSections[sectionIndex - 1].first;
        auto rightWindingIndex = orderedSections[sectionIndex].first;
        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex}; 
        if (!_insulationSections.contains(windingsMapKey)) {
            continue;
        }
        auto currentSpaceForLeftSection = orderedSections[sectionIndex - 1].second;
        auto currentSpaceForRightSection = orderedSections[sectionIndex].second;
        if (_windingOrientation == WindingOrientation::HORIZONTAL) {
            std::pair<size_t, double> leftSectionInfo = {leftWindingIndex, currentSpaceForLeftSection - _insulationSections[windingsMapKey].get_dimensions()[0] / 2};
            orderedSections[sectionIndex - 1] = leftSectionInfo;
            std::pair<size_t, double> rightSectionInfo = {rightWindingIndex, currentSpaceForRightSection - _insulationSections[windingsMapKey].get_dimensions()[0] / 2};
            orderedSections[sectionIndex] = rightSectionInfo;
        }
        else if (_windingOrientation == WindingOrientation::VERTICAL) {
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
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[0]};
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[1]};
            }

            orderedSectionsWithInsulation.push_back({ElectricalType::INSULATION, insulationSectionInfo});
        }
        orderedSectionsWithInsulation.push_back({ElectricalType::CONDUCTION, orderedSections[sectionIndex]});
    }

    return orderedSectionsWithInsulation;
}

std::vector<double> CoilWrapper::get_proportion_per_winding_based_on_wires() {
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

bool CoilWrapper::wind_by_sections() {
    // auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    auto proportionPerWinding = get_proportion_per_winding_based_on_wires();
    return wind_by_sections(proportionPerWinding);
}
bool CoilWrapper::wind_by_sections(std::vector<size_t> pattern, size_t repetitions) {
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind_by_sections(proportionPerWinding, pattern, repetitions);
}

bool CoilWrapper::wind_by_sections(std::vector<double> proportionPerWinding) {
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    return wind_by_sections(proportionPerWinding, pattern, _interleavingLevel);
}

void CoilWrapper::apply_margin_tape(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation) {
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

bool CoilWrapper::wind_by_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    if (_interleavingLevel <= 0) {
        throw std::runtime_error("Interleaving levels must be greater than 0");
    }
    set_sections_description(std::nullopt);
    std::vector<Section> sectionsDescription;
    auto bobbin = resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        throw std::runtime_error("Bobbin not processed");
    }
    auto bobbinProcessedDescription = std::get<Bobbin>(get_bobbin()).get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windingWindows.size() > 1) {
        throw std::runtime_error("Bobbins with more than winding window not implemented yet");
    }
    windingWindows[0].set_sections_orientation(_windingOrientation);
    bobbinProcessedDescription.set_winding_windows(windingWindows);
    bobbin.set_processed_description(bobbinProcessedDescription);
    set_bobbin(bobbin);

    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();

    double spaceForSections = 0;
    if (_windingOrientation == WindingOrientation::HORIZONTAL) {
        spaceForSections = windingWindowWidth;
    }
    else if (_windingOrientation == WindingOrientation::VERTICAL) {
        spaceForSections = windingWindowHeight;
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
   
    std::vector<std::vector<double>> remainingParallelsProportion;
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
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentSectionWidth = spaceForSection;
                // currentSectionHeight = windingWindowHeight - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1];
                currentSectionHeight = windingWindowHeight;

                if (currentSectionCenterWidth == DBL_MAX) {
                    currentSectionCenterWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                }
                if (currentSectionCenterHeight == DBL_MAX) {
                    currentSectionCenterHeight = windingWindows[0].get_coordinates().value()[1];
                }
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                currentSectionWidth = windingWindowWidth;
                currentSectionHeight = spaceForSection;
                if (currentSectionCenterWidth == DBL_MAX) {
                    currentSectionCenterWidth = windingWindows[0].get_coordinates().value()[0];
                }
                if (currentSectionCenterHeight == DBL_MAX) {
                    currentSectionCenterHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2;
                }
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
            }

            
            PartialWinding partialWinding;
            Section section;
            partialWinding.set_winding(get_name(windingIndex));

            auto parallels_proportions = get_parallels_proportions(currentSectionPerWinding[windingIndex],
                                                                   numberSectionsPerWinding[windingIndex],
                                                                   get_number_turns(windingIndex),
                                                                   get_number_parallels(windingIndex),
                                                                   remainingParallelsProportion[windingIndex],
                                                                   windByConsecutiveTurns[windingIndex],
                                                                   std::vector<double>(get_number_parallels(windingIndex), 1));

            std::vector<double> sectionParallelsProportion = parallels_proportions.second;
            uint64_t physicalTurnsThisSection = parallels_proportions.first;

            partialWinding.set_parallels_proportion(sectionParallelsProportion);
            section.set_name(get_name(windingIndex) +  " section " + std::to_string(currentSectionPerWinding[windingIndex]));
            section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
            section.set_type(ElectricalType::CONDUCTION);
            section.set_margin(_marginsPerSection[sectionIndex]);
            section.set_layers_orientation(_layersOrientation);
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                section.set_dimensions(std::vector<double>{currentSectionWidth, currentSectionHeight - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1]});
            }
            else {
                section.set_dimensions(std::vector<double>{currentSectionWidth - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1], currentSectionHeight});
            }

            if (section.get_dimensions()[0] < 0) {
                throw std::runtime_error("Something wrong happened in section dimn 1: " + std::to_string(section.get_dimensions()[0]) +
                                         " currentSectionWidth: " + std::to_string(currentSectionWidth) +
                                         " currentSectionHeight: " + std::to_string(currentSectionHeight)
                                         );
            }
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                section.set_coordinates(std::vector<double>{currentSectionCenterWidth + currentSectionWidth / 2, currentSectionCenterHeight, 0});
            }
            else {
                section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight - currentSectionHeight / 2, 0});
            }

            if (section.get_coordinates()[0] < -1) {
                throw std::runtime_error("Something wrong happened in section dimn 1: " + std::to_string(section.get_coordinates()[0]) +
                                         " currentSectionCenterWidth: " + std::to_string(currentSectionCenterWidth) +
                                         " windingWindows[0].get_coordinates().value()[0]: " + std::to_string(windingWindows[0].get_coordinates().value()[0]) +
                                         " windingWindows[0].get_width(): " + std::to_string(bool(windingWindows[0].get_width())) +
                                         " windingWindows[0].get_width().value(): " + std::to_string(windingWindows[0].get_width().value()) +
                                         " windingWindowWidth: " + std::to_string(windingWindowWidth) +
                                         " currentSectionWidth: " + std::to_string(currentSectionWidth) +
                                         " currentSectionCenterHeight: " + std::to_string(currentSectionCenterHeight)
                                         );
            }

            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
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

            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentSectionCenterWidth += currentSectionWidth;
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                currentSectionCenterHeight -= currentSectionHeight;
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
            }

            currentSectionPerWinding[windingIndex]++;
        }
        else {
            if (sectionIndex == 0 || sectionIndex == (orderedSectionsWithInsulation.size() - 1)) {
                throw std::runtime_error("Insulation layers cannot be the first or the last one (for now)");
            }

            auto previousWindingIndex = orderedSectionsWithInsulation[sectionIndex - 1].second.first;
            auto nextWindingIndex = orderedSectionsWithInsulation[sectionIndex + 1].second.first;

            auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, nextWindingIndex};
            if (!_insulationSections.contains(windingsMapKey)) {
                log(_insulationSectionsLog[windingsMapKey]);
                continue;
            }

            auto insulationSection = _insulationSections[windingsMapKey];

            insulationSection.set_name("Insulation between " + get_name(previousWindingIndex) + " and " + get_name(previousWindingIndex) + " section " + std::to_string(sectionIndex));
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth + insulationSection.get_dimensions()[0] / 2,
                                                                      currentSectionCenterHeight,
                                                                      0});
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth,
                                                                      currentSectionCenterHeight - insulationSection.get_dimensions()[1] / 2,
                                                                      0});
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
            }

            sectionsDescription.push_back(insulationSection);
            log(_insulationSectionsLog[windingsMapKey]);

            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentSectionCenterWidth += insulationSection.get_dimensions()[0];
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                currentSectionCenterHeight -= insulationSection.get_dimensions()[1];
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
            }

        }
    }

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            if (roundFloat(remainingParallelsProportion[windingIndex][parallelIndex], 9) != 0) {
                throw std::runtime_error("There are unassigned parallel proportion, something went wrong");
            }
        }
    }

    set_sections_description(sectionsDescription);
    return true;
}

bool CoilWrapper::wind_by_layers() {
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
                double wireDiameter = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_diameter().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::VERTICAL) {
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
                double wireWidth = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_width().value());
                double wireHeight = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_height().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::VERTICAL) {
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
                    if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR) {  // TODO: Fix for future PLANAR case
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

            double currentLayerCenterWidth;
            double currentLayerCenterHeight;
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::VERTICAL) {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerWidth / 2, 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0], 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);
                if (_turnsAlignment == CoilAlignment::CENTERED || _turnsAlignment == CoilAlignment::SPREAD) {
                    currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + (numberLayers * layerHeight) / 2 - layerHeight / 2, 9);
                }
                else if (_turnsAlignment == CoilAlignment::INNER_OR_TOP) {
                    currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);
                }
                else if (_turnsAlignment == CoilAlignment::OUTER_OR_BOTTOM) {
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

                auto parallels_proportions = get_parallels_proportions(layerIndex,
                                                                       numberLayers,
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportionInSection,
                                                                       windByConsecutiveTurns,
                                                                       totalParallelsProportionInSection);

                std::vector<double> layerParallelsProportion = parallels_proportions.second;

                size_t numberParallelsProportionsToZero = 0;
                for (auto parallelProportion : layerParallelsProportion) {
                    if (parallelProportion == 0) {
                        numberParallelsProportionsToZero++;
                    }
                }

                if (numberParallelsProportionsToZero == layerParallelsProportion.size()) {
                    throw std::runtime_error("Parallel proportion in layer cannot be all be 0");
                }

                uint64_t physicalTurnsThisLayer = parallels_proportions.first;

                partialWinding.set_parallels_proportion(layerParallelsProportion);
                layer.set_partial_windings(std::vector<PartialWinding>{partialWinding});
                layer.set_section(sections[sectionIndex].get_name());
                layer.set_type(ElectricalType::CONDUCTION);
                layer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                layer.set_orientation(sections[sectionIndex].get_layers_orientation());
                layer.set_turns_alignment(_turnsAlignment);
                layer.set_dimensions(std::vector<double>{layerWidth, layerHeight});
                layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});

                layer.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisLayer) / (layerWidth * layerHeight));
                layer.set_winding_style(windByConsecutiveTurns);
                layers.push_back(layer);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportionInSection[parallelIndex] -= layerParallelsProportion[parallelIndex];
                }

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::HORIZONTAL) {
                    currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight, 9);
                }
                else if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::VERTICAL) {
                    currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth, 9);
                }
                else {
                    throw std::runtime_error("Toroidal windings not implemented");
                }
            }
        }
        else {
            if ((sectionIndex == sections.size() - 1) || (sectionIndex == 0)) {
                throw std::runtime_error("Inner or outer insulation layers not implemented");
            }

            if (sections[sectionIndex - 1].get_type() != ElectricalType::CONDUCTION || sections[sectionIndex + 1].get_type() != ElectricalType::CONDUCTION) {
                throw std::runtime_error("Previous and next sections must be conductive");
            }
            auto partialWinding = sections[sectionIndex - 1].get_partial_windings()[0];
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            auto nextSection = sections[sectionIndex + 1];
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
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::VERTICAL) {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerWidth / 2, 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0], 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);
            }

            for (size_t layerIndex = 0; layerIndex < insulationLayers.size(); ++layerIndex) {
                auto insulationLayer = insulationLayers[layerIndex];
                insulationLayer.set_section(sections[sectionIndex].get_name());
                insulationLayer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layers.push_back(insulationLayer);

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::HORIZONTAL) {
                    currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight, 9);
                }
                else if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::VERTICAL) {
                    currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth, 9);
                }
                else {
                    throw std::runtime_error("Toroidal windings not implemented");
                }
            }
        }
    }
    set_layers_description(layers);
    return true;
}

bool CoilWrapper::wind_by_turns() {
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
        auto bobbinWindingWindow = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows()[0];
        double bobbinWindingWindowWidth = bobbinWindingWindow.get_width().value();
        double bobbinWindingWindowCenterWidth = bobbinWindingWindow.get_coordinates().value()[0];
        bobbinColumnWidth = bobbinWindingWindowCenterWidth - bobbinWindingWindowWidth / 2;
    }

    auto layers = get_layers_description().value();

    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            double currentTurnCenterWidth = 0;
            double currentTurnCenterHeight = 0;
            double currentTurnWidthIncrement = 0;
            double currentTurnHeightIncrement = 0;
            double totalLayerHeight;
            double totalLayerWidth;
            double wireHeight = 0;
            double wireWidth = 0;
            if (layer.get_partial_windings().size() > 1) {
                throw std::runtime_error("More than one winding per layer not supported yet");
            }
            auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            auto physicalTurnsInLayer = get_number_turns(layer);
            auto alignment = layer.get_turns_alignment().value();
            if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                wireWidth = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_diameter().value());
                wireHeight = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_diameter().value());
            }
            else {
                wireWidth = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_width().value());
                wireHeight = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_height().value());
            }
            if (layer.get_orientation() == WindingOrientation::VERTICAL) {
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
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - wireHeight / 2 - currentTurnHeightIncrement / 2, 9);
                        break;
                }

            } 
            else {
                totalLayerWidth = roundFloat(physicalTurnsInLayer * wireWidth, 9);
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
                        else if (bobbinColumnShape == ColumnShape::RECTANGULAR) {
                            double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                            turn.set_length(2 * bobbinColumnDepth + 2 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);

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
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR) {
                                double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                                turn.set_length(2 * bobbinColumnDepth + 2 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
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


std::vector<double> CoilWrapper::get_aligned_section_dimensions(size_t sectionIndex) {
    auto sections = get_sections_description().value();
    if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
        sections[sectionIndex].set_margin(std::vector<double>{0, 0});
    }

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();

    if (sections.size() == 0) {
        throw std::runtime_error("No sections in coil");
    }
    double totalSectionsWidth = 0;
    double totalSectionsHeight = 0;
    for (size_t auxSectionIndex = 0; auxSectionIndex < sections.size(); ++auxSectionIndex) {
        if (_windingOrientation == WindingOrientation::HORIZONTAL) {
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

    switch (_sectionAlignment) {
        case CoilAlignment::INNER_OR_TOP:
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (_turnsAlignment) {
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
                switch (_turnsAlignment) {
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
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - totalSectionsWidth;
                switch (_turnsAlignment) {
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
                switch (_turnsAlignment) {
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
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (_turnsAlignment) {
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

                switch (_turnsAlignment) {
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
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (_turnsAlignment) {
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
                switch (_turnsAlignment) {
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
        default:
            throw std::runtime_error("No such section alignment");

    }

    return {currentCoilWidth, currentCoilHeight, paddingAmongSectionWidth, paddingAmongSectionHeight};
}

bool CoilWrapper::delimit_and_compact() {
    // Delimit
    if (get_layers_description()) {
        auto layers = get_layers_description().value();
        if (get_turns_description()) {
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

                for (auto& turn : layersInSection) {
                    currentSectionMaximumWidth = std::max(currentSectionMaximumWidth, (turn.get_coordinates()[0] - sectionCoordinates[0]) + turn.get_dimensions()[0] / 2);
                    currentSectionMinimumWidth = std::min(currentSectionMinimumWidth, (turn.get_coordinates()[0] - sectionCoordinates[0]) - turn.get_dimensions()[0] / 2);
                    currentSectionMaximumHeight = std::max(currentSectionMaximumHeight, (turn.get_coordinates()[1] - sectionCoordinates[1]) + turn.get_dimensions()[1] / 2);
                    currentSectionMinimumHeight = std::min(currentSectionMinimumHeight, (turn.get_coordinates()[1] - sectionCoordinates[1]) - turn.get_dimensions()[1] / 2);
                }
                sections[i].set_coordinates(std::vector<double>({sectionCoordinates[0] + (currentSectionMaximumWidth + currentSectionMinimumWidth) / 2,
                                                           sectionCoordinates[1] + (currentSectionMaximumHeight + currentSectionMinimumHeight) / 2}));
                sections[i].set_dimensions(std::vector<double>({currentSectionMaximumWidth - currentSectionMinimumWidth,
                                                       currentSectionMaximumHeight - currentSectionMinimumHeight}));
            }
        }
        set_sections_description(sections);
    }

     // Compact
    if (get_sections_description()) {
        auto sections = get_sections_description().value();

        std::vector<std::vector<double>> alignedSectionDimensionsPerSection;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            alignedSectionDimensionsPerSection.push_back(get_aligned_section_dimensions(sectionIndex));
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

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentCoilHeight = alignedSectionDimensionsPerSection[sectionIndex][1];
                currentCoilWidth += sections[sectionIndex].get_dimensions()[0] / 2;
            }
            else {
                currentCoilHeight -= sections[sectionIndex].get_dimensions()[1] / 2;
                currentCoilWidth = alignedSectionDimensionsPerSection[sectionIndex][0];
            }

            double compactingShiftWidth = sections[sectionIndex].get_coordinates()[0] - currentCoilWidth;
            double compactingShiftHeight = sections[sectionIndex].get_coordinates()[1] - currentCoilHeight;

            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
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
                                else if (bobbinColumnShape == ColumnShape::RECTANGULAR) {
                                    double currentTurnCornerRadius = turns[turnIndex].get_coordinates()[0] - bobbinColumnWidth;
                                    turns[turnIndex].set_length(2 * bobbinColumnDepth + 2 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);

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
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
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
    auto settings = Settings::GetInstance();
    bool fillCoilSectionsWithMarginTape = settings->get_fill_coil_sections_with_margin_tape();


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
                if (sectionOrientation == WindingOrientation::HORIZONTAL) {
                    auto topSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[1] + windingWindowHeight / 2) - (sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2));
                    auto bottomSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[1] - windingWindowHeight / 2) - (sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2));
                    // if (sections[i].get_margin() && topSpaceBetweenSectionAndBobbin < sections[i].get_margin().value()[0]) {
                    //     throw std::runtime_error("topSpaceBetweenSectionAndBobbin cannot be smaller than top margin");
                    // }
                    // if (sections[i].get_margin() && bottomSpaceBetweenSectionAndBobbin < sections[i].get_margin().value()[1]) {
                    //     throw std::runtime_error("bottomSpaceBetweenSectionAndBobbin cannot be smaller than bottom margin");
                    // }
                    sections[i].set_margin(std::vector<double>{topSpaceBetweenSectionAndBobbin, bottomSpaceBetweenSectionAndBobbin});
                }
                else if (sectionOrientation == WindingOrientation::VERTICAL) {
                    auto innerSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[0] - windingWindowWidth / 2) - (sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2));
                    auto outerSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[0] + windingWindowWidth / 2) - (sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2));
                    // if (sections[i].get_margin() && innerSpaceBetweenSectionAndBobbin < sections[i].get_margin().value()[0]) {
                    //     throw std::runtime_error("innerSpaceBetweenSectionAndBobbin cannot be smaller than inner margin");
                    // }
                    // if (sections[i].get_margin() && outerSpaceBetweenSectionAndBobbin < sections[i].get_margin().value()[1]) {
                    //     throw std::runtime_error("outerSpaceBetweenSectionAndBobbin cannot be smaller than outer margin");
                    // }
                    sections[i].set_margin(std::vector<double>{innerSpaceBetweenSectionAndBobbin, outerSpaceBetweenSectionAndBobbin});
                }
            }
        }
        set_sections_description(sections);
    }

    return true;
}

std::vector<WireWrapper> CoilWrapper::get_wires() {
    std::vector<WireWrapper> wirePerWinding;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        WireWrapper wire = resolve_wire(get_functional_description()[windingIndex]);
        wirePerWinding.push_back(wire);
    }
    return wirePerWinding;
}

WireWrapper CoilWrapper::resolve_wire(size_t windingIndex) {
    return resolve_wire(get_functional_description()[windingIndex]);
}

WireWrapper CoilWrapper::resolve_wire(CoilFunctionalDescription coilFunctionalDescription) {
    auto wireOrString = coilFunctionalDescription.get_wire();
    WireWrapper wire;
    if (std::holds_alternative<std::string>(wireOrString)) {
        wire = find_wire_by_name(std::get<std::string>(wireOrString));
    }
    else {
        wire = std::get<Wire>(wireOrString);
    }
    return wire;
}

WireType CoilWrapper::get_wire_type(CoilFunctionalDescription coilFunctionalDescription) {
    return resolve_wire(coilFunctionalDescription).get_type();
}

WireType CoilWrapper::get_wire_type(size_t windingIndex) {
    return get_wire_type(get_functional_description()[windingIndex]);
}

BobbinWrapper CoilWrapper::resolve_bobbin(CoilWrapper coil) { 
    return coil.resolve_bobbin();
}

BobbinWrapper CoilWrapper::resolve_bobbin() { 
    auto bobbinDataOrNameUnion = get_bobbin();
    if (std::holds_alternative<std::string>(bobbinDataOrNameUnion)) {
        if (std::get<std::string>(bobbinDataOrNameUnion) == "Dummy")
            throw std::runtime_error("Bobbin is dummy");

        auto bobbin = find_bobbin_by_name(std::get<std::string>(bobbinDataOrNameUnion));
        return BobbinWrapper(bobbin);
    }
    else {
        return BobbinWrapper(std::get<Bobbin>(bobbinDataOrNameUnion));
    }
}

size_t CoilWrapper::convert_conduction_section_index_to_global(size_t conductionSectionIndex) {
    size_t currentConductionSectionIndex = 0;
    if (!get_sections_description()) {
        throw std::runtime_error("Section description empty, wind coil first");
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

void CoilWrapper::try_rewind() {
    if (!get_sections_description()) {
        throw std::runtime_error("Section description empty, wind coil first");
    }
    auto sections = get_sections_description().value();
    std::vector<double> extraSpaceNeededPerSection;
    double totalExtraSpaceNeeded = 0;
    auto bobbin = resolve_bobbin();
    auto sectionOrientation = bobbin.get_winding_window_sections_orientation(0);
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    double windingWindowRemainingRestrictiveDimension;
    double windingWindowRestrictiveDimension;
    if (sectionOrientation == WindingOrientation::HORIZONTAL) {
        windingWindowRemainingRestrictiveDimension = windingWindowDimensions[0];
        windingWindowRestrictiveDimension = windingWindowDimensions[0];
    }
    else {
        windingWindowRemainingRestrictiveDimension = windingWindowDimensions[1];
        windingWindowRestrictiveDimension = windingWindowDimensions[1];
    }

    for (auto& section : sections) {
        if (section.get_type() == ElectricalType::INSULATION) {
            if (sectionOrientation == WindingOrientation::HORIZONTAL) {
                windingWindowRestrictiveDimension -= section.get_dimensions()[0];
            }
            else {
                windingWindowRestrictiveDimension -= section.get_dimensions()[1];
            }

        }

        double sectionRestrictiveDimension;
        double sectionFillingFactor;
        double extraSpaceNeededThisSection = 0;

        auto layers = get_layers_by_section(section.get_name());
        for (auto layer : layers) {
            double layerRestrictiveDimension;
            double sectionFillingFactor = layer.get_filling_factor().value();
            if (section.get_layers_orientation() == WindingOrientation::VERTICAL) {
                layerRestrictiveDimension = layer.get_dimensions()[0];
            }
            else {
                layerRestrictiveDimension = layer.get_dimensions()[1];
            }

            extraSpaceNeededThisSection += std::max(0.0, (sectionFillingFactor - 1) * layerRestrictiveDimension);
            windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
        }

        if (sectionOrientation == WindingOrientation::HORIZONTAL) {
            sectionRestrictiveDimension = section.get_dimensions()[0];
            sectionFillingFactor = horizontal_filling_factor(section);
        }
        else {
            sectionRestrictiveDimension = section.get_dimensions()[1];
            sectionFillingFactor = vertical_filling_factor(section);
        }


        extraSpaceNeededThisSection = std::max(extraSpaceNeededThisSection, (sectionFillingFactor - 1) * sectionRestrictiveDimension);
        extraSpaceNeededPerSection.push_back(extraSpaceNeededThisSection);
        totalExtraSpaceNeeded += extraSpaceNeededThisSection;
    }

    std::vector<double> newProportions;

    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        double currentProportion = _currentProportionPerWinding[windingIndex];
        double currentSpace = 0;
        double extraSpaceNeededThisWinding = 0;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            for (auto & winding : sections[sectionIndex].get_partial_windings()) {
                if (winding.get_winding() == get_functional_description()[windingIndex].get_name()) {
                    if (sectionOrientation == WindingOrientation::HORIZONTAL) {
                        currentSpace += sections[sectionIndex].get_dimensions()[0];
                    }
                    else {
                        currentSpace += sections[sectionIndex].get_dimensions()[1];
                    }
                    extraSpaceNeededThisWinding += extraSpaceNeededPerSection[sectionIndex];
                    continue;
                }
            }
        }
        double proportionOfNeededForThisWinding = extraSpaceNeededThisWinding / totalExtraSpaceNeeded;
        double extraSpaceGottenByThisWinding = windingWindowRemainingRestrictiveDimension * extraSpaceNeededThisWinding / totalExtraSpaceNeeded;
        double newSpaceGottenByThisWinding = currentSpace + extraSpaceGottenByThisWinding;
        double newProportionGottenByThisWinding = newSpaceGottenByThisWinding / windingWindowRestrictiveDimension;
        newProportions.push_back(newProportionGottenByThisWinding);
    }

    wind_by_sections(newProportions, _currentPattern, _currentRepetitions);
    wind_by_layers();
    wind_by_turns();
    delimit_and_compact();
}


void CoilWrapper::add_margin_to_section_by_index(size_t sectionIndex, std::vector<double> margins) {
    if (!get_sections_description()) {
        throw std::runtime_error("Section description empty, wind coil first");
    }
    if (margins.size() != 2) {
        throw std::runtime_error("Margin vector must have two elements");
    }
    auto sections = get_sections_description().value();
    _marginsPerSection[convert_conduction_section_index_to_global(sectionIndex)] = margins;
    sections[convert_conduction_section_index_to_global(sectionIndex)].set_margin(margins);

    set_sections_description(sections);


    auto settings = Settings::GetInstance();
    bool windEvenIfNotFit = settings->get_wind_even_if_not_fit();
    wind_by_sections();
    wind_by_layers();
    wind_by_turns();
    delimit_and_compact();
    // std::cout << are_sections_and_layers_fitting() << std::endl;
    if (!are_sections_and_layers_fitting()) {
        try_rewind();
    }

}

std::vector<Section> CoilWrapper::get_sections_description_conduction() {
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

std::vector<Layer> CoilWrapper::get_layers_description_conduction() {
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

std::vector<Section> CoilWrapper::get_sections_description_insulation() {
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

std::vector<Layer> CoilWrapper::get_layers_description_insulation() {
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

} // namespace OpenMagnetics
 
