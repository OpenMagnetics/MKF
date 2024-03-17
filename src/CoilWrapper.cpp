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


std::vector<double> CoilWrapper::cartesian_to_polar(std::vector<double> value) {
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


CoilWrapper::CoilWrapper(const json& j, size_t interleavingLevel,
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

    auto settings = Settings::GetInstance();
    auto delimitAndCompact = settings->get_coil_delimit_and_compact();

    wind();
}

CoilWrapper::CoilWrapper(const Coil coil) {
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
    auto settings = Settings::GetInstance();
    auto delimitAndCompact = settings->get_coil_delimit_and_compact();

    if (!hasSectionsData || !hasLayersData || (!hasTurnsData && are_sections_and_layers_fitting())) {
        if (wind() && delimitAndCompact) {
            delimit_and_compact();
        }
    }

}

bool CoilWrapper::try_wind() {
    auto settings = Settings::GetInstance();
    bool windEvenIfNotFit = settings->get_coil_wind_even_if_not_fit();
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
    auto delimitAndCompact = settings->get_coil_delimit_and_compact();

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

bool CoilWrapper::wind(size_t repetitions){
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind(proportionPerWinding, pattern, repetitions);
}

bool CoilWrapper::wind(std::vector<size_t> pattern, size_t repetitions){
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind(proportionPerWinding, pattern, repetitions);
}

bool CoilWrapper::wind(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    auto settings = Settings::GetInstance();
    bool windEvenIfNotFit = settings->get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings->get_coil_delimit_and_compact();
    bool tryRewind = settings->get_coil_try_rewind();
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
    }
    return are_sections_and_layers_fitting();
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

WindingStyle CoilWrapper::wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, size_t numberSlots) {
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

double CoilWrapper::overlapping_filling_factor(Section section) {
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

double CoilWrapper::contiguous_filling_factor(Section section) {
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
            double numberTurnsToAddToCurrentParallel = ceil(remainingTurnsBeforeThisParallel / remainingSlots * slotRelativeProportion);
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

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto layersOrientation = _layersOrientation;

    // TODO: Properly think about insulation layers with weird windings

    if (_windingOrientation == WindingOrientation::CONTIGUOUS && _layersOrientation == WindingOrientation::OVERLAPPING) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::CONTIGUOUS;
        }
    }
    if (_windingOrientation == WindingOrientation::OVERLAPPING && _layersOrientation == WindingOrientation::CONTIGUOUS) {
        layersOrientation = WindingOrientation::OVERLAPPING;
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
                    if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{defaultInsulationMaterial.get_thinner_tape_thickness(), windingWindowAngle});
                    }
                    else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
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
                if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowHeight});
                }
                else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
                    section.set_dimensions(std::vector<double>{windingWindowWidth, coilSectionInterface.get_solid_insulation_thickness()});
                }
            }
            else {
                section.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowAngle});
                }
                else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
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

bool CoilWrapper::calculate_insulation(bool simpleMode) {
    auto inputs = _inputs.value();

    if (!inputs.get_design_requirements().get_insulation()) {
        return false;
    }

    auto settings = Settings::GetInstance();
    auto wirePerWinding = get_wires();

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            if (leftTopWindingIndex == rightBottomWindingIndex) {
                continue;
            }
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            CoilSectionInterface coilSectionInterface;
            InsulationMaterialWrapper chosenInsulationMaterial;

            if (simpleMode) {
                InsulationMaterialWrapper defaultInsulationMaterial = find_insulation_material_by_name(Defaults().defaultInsulationMaterial);
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
                layer.set_turns_alignment(_turnsAlignment);
                if (_layersOrientation == WindingOrientation::OVERLAPPING) {
                    layer.set_dimensions(std::vector<double>{chosenInsulationMaterial.get_thicker_tape_thickness(), windingWindowHeight});
                }
                else if (_layersOrientation == WindingOrientation::CONTIGUOUS) {
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
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowHeight});
            }
            else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
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
            if (roundFloat(proportionPerWinding[windingIndex], 6) > 1) {
                throw std::invalid_argument("proportionPerWinding[windingIndex] cannot be greater than 1: " + std::to_string(proportionPerWinding[windingIndex]));
            }
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
        if (_windingOrientation == WindingOrientation::OVERLAPPING) {
            std::pair<size_t, double> leftSectionInfo = {leftWindingIndex, currentSpaceForLeftSection - _insulationSections[windingsMapKey].get_dimensions()[0] / 2};
            orderedSections[sectionIndex - 1] = leftSectionInfo;
            std::pair<size_t, double> rightSectionInfo = {rightWindingIndex, currentSpaceForRightSection - _insulationSections[windingsMapKey].get_dimensions()[0] / 2};
            orderedSections[sectionIndex] = rightSectionInfo;
        }
        else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
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
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[0]};
            }
            else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[1]};
            }

            orderedSectionsWithInsulation.push_back({ElectricalType::INSULATION, insulationSectionInfo});
        }
        orderedSectionsWithInsulation.push_back({ElectricalType::CONDUCTION, orderedSections[sectionIndex]});
    }


    // last insulation layer we compare between last and first
    auto leftWindingIndex = orderedSections[orderedSections.size() - 1].first;
    auto rightWindingIndex = orderedSections[0].first;
    auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex}; 

    if (_insulationSections.contains(windingsMapKey)) {
        std::pair<size_t, double> insulationSectionInfo;
        if (_windingOrientation == WindingOrientation::OVERLAPPING) {
            insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[0]};
        }
        else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
            insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[1]};
        }

        orderedSectionsWithInsulation.push_back({ElectricalType::INSULATION, insulationSectionInfo});
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

std::pair<double, double> get_section_round_dimensions(std::pair<ElectricalType, std::pair<size_t, double>> sectionWithInsulationScaledWithArea,
                                                 WindingOrientation windingOrientation, WindingWindowElement windingWindow) {

    double windingWindowRadialHeight = windingWindow.get_radial_height().value();
    double windingWindowAngle = windingWindow.get_angle().value();

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

std::pair<size_t, std::vector<int64_t>> get_number_layers_needed_and_number_physical_turns(double radialHeight, double angle, WireWrapper wire, int64_t physicalTurnsInSection, double windingWindowRadius) {
    int64_t reaminingPhysicalTurnsInSection = physicalTurnsInSection;
    double wireWidth = resolve_dimensional_values(wire.get_maximum_outer_width());
    double wireHeight = resolve_dimensional_values(wire.get_maximum_outer_height());
    double currentRadialHeight = radialHeight;
    double currentRadius = windingWindowRadius - wireWidth / 2 - currentRadialHeight;
    double sectionAvailableAngle = angle;
    std::vector<int64_t> layerPhysicalTurns;
    size_t numberLayers = 0;
    while (reaminingPhysicalTurnsInSection > 0) {
        double wireAngle = wound_distance_to_angle(wireHeight, currentRadius);
        int64_t numberTurnsFittingThisLayer = floor(sectionAvailableAngle / wireAngle);
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

std::pair<size_t, std::vector<int64_t>> get_number_layers_needed_and_number_physical_turns(Section section, WireWrapper wire, int64_t physicalTurnsInSection, double windingWindowRadius) {
    return get_number_layers_needed_and_number_physical_turns(section.get_coordinates()[0] - section.get_dimensions()[0] / 2, section.get_dimensions()[1], wire, physicalTurnsInSection, windingWindowRadius);
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

bool CoilWrapper::wind_by_sections() {
    // auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    auto proportionPerWinding = get_proportion_per_winding_based_on_wires();
    return wind_by_sections(proportionPerWinding);
}

bool CoilWrapper::wind_by_sections(size_t repetitions){
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind_by_sections(proportionPerWinding, pattern, repetitions);
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

bool CoilWrapper::wind_by_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    _currentProportionPerWinding = proportionPerWinding;
    _currentPattern = pattern;
    _currentRepetitions = repetitions;

    if (repetitions <= 0) {
        throw std::runtime_error("Interleaving levels must be greater than 0");
    }

    auto bobbin = resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        throw std::runtime_error("Bobbin not processed");
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windingWindows.size() > 1) {
        throw std::runtime_error("Bobbins with more than winding window not implemented yet");
    }
    windingWindows[0].set_sections_orientation(_windingOrientation);
    bobbinProcessedDescription.set_winding_windows(windingWindows);
    bobbin.set_processed_description(bobbinProcessedDescription);
    set_bobbin(bobbin);

    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        return wind_by_rectangular_sections(proportionPerWinding, pattern, repetitions);
    }
    else {
        return wind_by_round_sections(proportionPerWinding, pattern, repetitions);
    }
}

bool CoilWrapper::wind_by_rectangular_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    set_sections_description(std::nullopt);
    std::vector<Section> sectionsDescription;
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();
    double spaceForSections = 0;
    if (_windingOrientation == WindingOrientation::OVERLAPPING) {
        spaceForSections = windingWindowWidth;
    }
    else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
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
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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
            else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
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
            section.set_coordinate_system(CoordinateSystem::CARTESIAN);


            
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                section.set_dimensions(std::vector<double>{currentSectionWidth, currentSectionHeight - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1]});
            }
            else {
                section.set_dimensions(std::vector<double>{currentSectionWidth - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1], currentSectionHeight});
            }

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                if ((section.get_margin().value()[0] + section.get_margin().value()[1] + resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_height())) > currentSectionHeight) {
                    throw std::runtime_error("Margin plus a turn cannot larger than winding window" + 
                                             std::string{", margin:"} + std::to_string(section.get_margin().value()[0] + section.get_margin().value()[1]) + 
                                             ", wire height:" + std::to_string(resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_height())) + 
                                             ", section height:" + std::to_string(currentSectionHeight)
                                             );
                }
            }
            else {
                if ((section.get_margin().value()[0] + section.get_margin().value()[1] + resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_width())) > currentSectionWidth) {
                    throw std::runtime_error("Margin plus a turn cannot larger than winding window");
                }
            }

            if (section.get_dimensions()[0] < 0) {
                throw std::runtime_error("Something wrong happened in section dimensions 0: " + std::to_string(section.get_dimensions()[0]) +
                                         " currentSectionWidth: " + std::to_string(currentSectionWidth) +
                                         " currentSectionHeight: " + std::to_string(currentSectionHeight)
                                         );
            }
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                section.set_coordinates(std::vector<double>{currentSectionCenterWidth + currentSectionWidth / 2, currentSectionCenterHeight, 0});
            }
            else {
                section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight - currentSectionHeight / 2, 0});
            }

            if (section.get_coordinates()[0] < -1) {
                throw std::runtime_error("Something wrong happened in section coordiantes 0: " + std::to_string(section.get_coordinates()[0]) +
                                         " currentSectionCenterWidth: " + std::to_string(currentSectionCenterWidth) +
                                         " windingWindows[0].get_coordinates().value()[0]: " + std::to_string(windingWindows[0].get_coordinates().value()[0]) +
                                         " windingWindows[0].get_width(): " + std::to_string(bool(windingWindows[0].get_width())) +
                                         " windingWindows[0].get_width().value(): " + std::to_string(windingWindows[0].get_width().value()) +
                                         " windingWindowWidth: " + std::to_string(windingWindowWidth) +
                                         " currentSectionWidth: " + std::to_string(currentSectionWidth) +
                                         " currentSectionCenterHeight: " + std::to_string(currentSectionCenterHeight)
                                         );
            }

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                currentSectionCenterWidth += currentSectionWidth;
            }
            else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
                currentSectionCenterHeight -= currentSectionHeight;
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
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

            insulationSection.set_name("Insulation between " + get_name(previousWindingIndex) + " and " + get_name(previousWindingIndex) + " section " + std::to_string(sectionIndex));
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth + insulationSection.get_dimensions()[0] / 2,
                                                                      currentSectionCenterHeight,
                                                                      0});
            }
            else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
                insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth,
                                                                      currentSectionCenterHeight - insulationSection.get_dimensions()[1] / 2,
                                                                      0});
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
            }

            sectionsDescription.push_back(insulationSection);
            log(_insulationSectionsLog[windingsMapKey]);

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                currentSectionCenterWidth += insulationSection.get_dimensions()[0];
            }
            else if (_windingOrientation == WindingOrientation::CONTIGUOUS) {
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
                throw std::runtime_error("There are unassigned parallel proportion in a rectangular section, something went wrong");
            }
        }
    }

    set_sections_description(sectionsDescription);
    return true;
}

bool CoilWrapper::wind_by_round_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    set_sections_description(std::nullopt);
    std::vector<Section> sectionsDescription;
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
    double windingWindowAngle = windingWindows[0].get_angle().value();

    double spaceForSections = 0;
    if (_windingOrientation == WindingOrientation::OVERLAPPING) {
        spaceForSections = windingWindowRadialHeight;
    }
    else {
        spaceForSections = windingWindowAngle;
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
    double currentSectionCenterAngle = DBL_MAX;
    double currentSectionCenterRadialHeight = DBL_MAX;

    apply_margin_tape(orderedSectionsWithInsulation);

    std::vector<double> currentSectionRadialHeights;
    std::vector<double> currentSectionAngles;
    std::vector<size_t> windingIndexes;

    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
        if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
            auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
            auto windingIndex = sectionInfo.first;
            auto aux = get_section_round_dimensions(orderedSectionsWithInsulation[sectionIndex], _windingOrientation, windingWindows[0]);
            currentSectionRadialHeights.push_back(aux.first);
            currentSectionAngles.push_back(aux.second);
            windingIndexes.push_back(windingIndex);
        }
    }
    
    // std::vector<double> sectionLengths = get_section_lengths(currentSectionRadialHeights, currentSectionAngles, windingWindowRadialHeight);

    std::vector<double> sectionPhysicalTurnsProportions;
    if (_windingOrientation == WindingOrientation::OVERLAPPING) {
        std::vector<double> sectionAreas = get_section_areas(orderedSectionsWithInsulation, currentSectionAngles, windingWindowRadialHeight);
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



            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                if (currentSectionCenterRadialHeight == DBL_MAX) {
                    currentSectionCenterRadialHeight = 0;
                }
                if (currentSectionCenterAngle == DBL_MAX) {
                    currentSectionCenterAngle = 180;
                }
            }
            else {
                if (currentSectionCenterRadialHeight == DBL_MAX) {
                    currentSectionCenterRadialHeight = windingWindowRadialHeight / 2;
                }
                if (currentSectionCenterAngle == DBL_MAX) {
                    currentSectionCenterAngle = 0;
                }
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
                                                                   std::vector<double>(get_number_parallels(windingIndex), 1),
                                                                   sectionPhysicalTurnsProportions[windingIndex]);

            std::vector<double> sectionParallelsProportion = parallels_proportions.second;
            uint64_t physicalTurnsThisSection = parallels_proportions.first;

            std::cout << "sectionIndex: " <<  sectionIndex << std::endl;
            std::cout << "sectionIndex: " <<  sectionIndex << std::endl;
            std::cout << "physicalTurnsThisSection: " <<  physicalTurnsThisSection << std::endl;
            std::cout << "get_number_turns(windingIndex) * get_number_parallels(windingIndex): " <<  (get_number_turns(windingIndex) * get_number_parallels(windingIndex)) << std::endl;

            // If OVERLAPPING, we correct the radial height to exactly what we need, so afterwards we can calculate exactly how many turns we need
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                auto aux = get_number_layers_needed_and_number_physical_turns(currentSectionCenterRadialHeight + _marginsPerSection[sectionIndex][0], currentSectionAngle, wirePerWinding[windingIndex], physicalTurnsThisSection, windingWindowRadialHeight);
                auto numberLayers = aux.first;
                currentSectionRadialHeight = numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width();
            }

            partialWinding.set_parallels_proportion(sectionParallelsProportion);
            section.set_name(get_name(windingIndex) +  " section " + std::to_string(currentSectionPerWinding[windingIndex]));
            section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
            section.set_type(ElectricalType::CONDUCTION);
            section.set_margin(_marginsPerSection[sectionIndex]);
            section.set_layers_orientation(_layersOrientation);
            section.set_coordinate_system(CoordinateSystem::POLAR);
            
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                section.set_dimensions(std::vector<double>{currentSectionRadialHeight, currentSectionAngle});
                section.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + currentSectionRadialHeight / 2 + _marginsPerSection[sectionIndex][0], currentSectionCenterAngle, 0});
            }
            else {
                double firstLayerMinimumRadius = windingWindowRadialHeight - (currentSectionCenterRadialHeight - currentSectionRadialHeight / 2) - wirePerWinding[windingIndex].get_maximum_outer_dimension();
                double marginAngle0 = wound_distance_to_angle(_marginsPerSection[sectionIndex][0], firstLayerMinimumRadius);
                double marginAngle1 = wound_distance_to_angle(_marginsPerSection[sectionIndex][1], firstLayerMinimumRadius);
                section.set_dimensions(std::vector<double>{currentSectionRadialHeight, currentSectionAngle - marginAngle0 - marginAngle1});
                section.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight, currentSectionCenterAngle  + currentSectionAngle / 2 - marginAngle0 / 2 + marginAngle1 / 2, 0});
            }

            if (section.get_dimensions()[0] < 0) {
                throw std::runtime_error("Something wrong happened in section dimensions 0: " + std::to_string(section.get_dimensions()[0]) +
                                         " currentSectionRadialHeight: " + std::to_string(currentSectionRadialHeight) +
                                         " currentSectionAngle: " + std::to_string(currentSectionAngle)
                                         );
            }

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                double ringArea = std::numbers::pi * pow(windingWindowRadialHeight - currentSectionCenterRadialHeight, 2) - std::numbers::pi * pow(windingWindowRadialHeight - (currentSectionCenterRadialHeight - currentSectionRadialHeight), 2);
                section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / ringArea);
            } 
            else {
                double ringArea = std::numbers::pi * pow(windingWindowRadialHeight - currentSectionCenterRadialHeight, 2) * currentSectionAngle / 360;
                section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / ringArea);
            }
            section.set_winding_style(windByConsecutiveTurns[windingIndex]);
            sectionsDescription.push_back(section);

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                remainingParallelsProportion[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
            }

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                currentSectionCenterRadialHeight += currentSectionRadialHeight + _marginsPerSection[sectionIndex][0] + _marginsPerSection[sectionIndex][1];
            }
            else {
                currentSectionCenterAngle += currentSectionAngle;
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

            insulationSection.set_name("Insulation between " + get_name(previousWindingIndex) + " and " + get_name(previousWindingIndex) + " section " + std::to_string(sectionIndex));
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                insulationSection.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + insulationSection.get_dimensions()[0] / 2,
                                                                      currentSectionCenterAngle,
                                                                      0});
            }
            else {
                insulationSection.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight,
                                                                      currentSectionCenterAngle - insulationSection.get_dimensions()[1] / 2,
                                                                      0});
            }

            sectionsDescription.push_back(insulationSection);
            log(_insulationSectionsLog[windingsMapKey]);

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                currentSectionCenterRadialHeight += insulationSection.get_dimensions()[0];
            }
            else {
                currentSectionCenterAngle += insulationSection.get_dimensions()[1];
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

bool CoilWrapper::wind_by_layers() {
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

bool CoilWrapper::wind_by_rectangular_layers() {
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
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
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
                else if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth, 9);
                }
                else {
                    throw std::runtime_error("Toroidal windings not implemented");
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
                insulationLayer.set_section(sections[sectionIndex].get_name());
                insulationLayer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layers.push_back(insulationLayer);

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                    currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight, 9);
                }
                else if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
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

bool CoilWrapper::wind_by_round_layers() {
    set_layers_description(std::nullopt);
    if (!get_sections_description()) {
        return false;
    }
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();


    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
    double windingWindowAngle = windingWindows[0].get_angle().value();

    auto wirePerWinding = get_wires();

    auto sections = get_sections_description().value();

    std::vector<Layer> layers;
    for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
        if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
            uint64_t maximumNumberLayersFittingInSection;
            uint64_t maximumNumberPhysicalTurnsPerLayer;
            uint64_t minimumNumberLayerNeeded;
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

            if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                double wireDiameter = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_diameter().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[0] / wireDiameter;
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
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[0] / wireWidth;
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

                auto parallels_proportions = get_parallels_proportions(layerIndex,
                                                                       numberLayers,
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportionInSection,
                                                                       windByConsecutiveTurns,
                                                                       totalParallelsProportionInSection,
                                                                       1,
                                                                       layerPhysicalTurns[layerIndex]);

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
            double layerHeight = insulationLayers[0].get_dimensions()[1];

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

bool CoilWrapper::wind_by_turns() {
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

bool CoilWrapper::wind_by_rectangular_turns() {
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
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - wireHeight / 2 - currentTurnHeightIncrement / 2, 9);
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

bool CoilWrapper::wind_by_round_turns() {
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

    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            double currentTurnCenterRadialHeight = 0;
            double currentTurnCenterAngle = 0;
            double currentTurnRadialHeightIncrement = 0;
            double currentTurnAngleIncrement = 0;
            double totalLayerRadialHeight;
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

            double currentLayerRadius = windingWindowRadialHeight - layer.get_coordinates()[0];
            double wireAngle = wound_distance_to_angle(wireHeight, currentLayerRadius);

            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                totalLayerRadialHeight = layer.get_dimensions()[0];
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
                            // TODO: Calculate proper turn length for overlapping layers
                            turn.set_length(2 * std::numbers::pi * currentTurnCenterRadialHeight);
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::RECTANGULAR) {
                            // TODO: Calculate proper turn length for overlapping layers
                            double currentTurnCornerRadius = currentTurnCenterRadialHeight - bobbinColumnWidth;
                            turn.set_length(2 * bobbinColumnDepth + 2 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);

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
                        turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

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
                            // TODO: Calculate proper turn length for overlapping layers
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterRadialHeight);
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR) {
                            // TODO: Calculate proper turn length for overlapping layers
                                double currentTurnCornerRadius = currentTurnCenterRadialHeight - bobbinColumnWidth;
                                turn.set_length(2 * bobbinColumnDepth + 2 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
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
                            turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

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
    return true;
}

std::vector<double> CoilWrapper::get_aligned_section_dimensions_rectangular_window(size_t sectionIndex) {
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
        if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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

std::vector<double> CoilWrapper::get_aligned_section_dimensions_round_window(size_t sectionIndex) {
    auto sections = get_sections_description().value();
    if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
        sections[sectionIndex].set_margin(std::vector<double>{0, 0});
    }

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowAngle = windingWindows[0].get_angle().value();
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    if (sections.size() == 0) {
        throw std::runtime_error("No sections in coil");
    }
    double totalSectionsRadialHeight = 0;
    double totalSectionsAngle = 0;
    for (size_t auxSectionIndex = 0; auxSectionIndex < sections.size(); ++auxSectionIndex) {
        if (_windingOrientation == WindingOrientation::OVERLAPPING) {
            totalSectionsRadialHeight += sections[auxSectionIndex].get_dimensions()[0];
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsAngle = std::max(totalSectionsAngle, sections[auxSectionIndex].get_dimensions()[1]);
            }
        }
        else {
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsRadialHeight = std::max(totalSectionsRadialHeight, sections[auxSectionIndex].get_dimensions()[0]);
            }
            totalSectionsAngle += sections[auxSectionIndex].get_dimensions()[1];
        }
    }

    double currentCoilRadialHeight;
    double currentCoilAngle;
    double paddingAmongSectionRadialHeight = 0;
    double paddingAmongSectionAngle = 0;
    double marginAngle0 = 0;
    double marginAngle1 = 0;

    if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {

        auto partialWinding = sections[sectionIndex].get_partial_windings()[0];  // TODO: Support multiwinding in layers
        auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

        auto wirePerWinding = get_wires();
        double firstLayerMinimumRadius = windingWindowRadialHeight - (sections[sectionIndex].get_coordinates()[0] - totalSectionsRadialHeight / 2) - wirePerWinding[windingIndex].get_maximum_outer_dimension();
        marginAngle0 = wound_distance_to_angle(sections[sectionIndex].get_margin().value()[0], firstLayerMinimumRadius);
        marginAngle1 = wound_distance_to_angle(sections[sectionIndex].get_margin().value()[1], firstLayerMinimumRadius);
    }

    if (_windingOrientation == WindingOrientation::OVERLAPPING) {
        currentCoilRadialHeight = 0;
        switch (_turnsAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                currentCoilAngle = marginAngle0 + sections[sectionIndex].get_dimensions()[1] / 2;
                break;
            case CoilAlignment::OUTER_OR_BOTTOM:
                currentCoilAngle = windingWindowAngle - sections[sectionIndex].get_dimensions()[1] / 2 - marginAngle1;
                break;
            case CoilAlignment::CENTERED:
                {
                    currentCoilAngle = 180;
                    if (sections[sectionIndex].get_coordinates()[1] - sections[sectionIndex].get_dimensions()[1] / 2 < marginAngle0) {
                        currentCoilAngle = 180 - marginAngle0 / 2;
                    }
                    if (sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 > windingWindowAngle - marginAngle1) {
                        currentCoilAngle = 180 + marginAngle1 / 2;
                    }
                    break;
                }
                break;
            case CoilAlignment::SPREAD:
                currentCoilAngle = marginAngle0 + sections[sectionIndex].get_dimensions()[1] / 2;
                break;
            default:
                throw std::runtime_error("No such section alignment");
        }
    }

    else {
        switch (_sectionAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                currentCoilAngle = sections[sectionIndex].get_coordinates()[1] - sections[sectionIndex].get_dimensions()[1] / 2 + marginAngle1;
                break;
            case CoilAlignment::OUTER_OR_BOTTOM:
                currentCoilAngle = windingWindowAngle - totalSectionsAngle - marginAngle1;
                break;
            case CoilAlignment::SPREAD:
                currentCoilAngle = sections[sectionIndex].get_coordinates()[1] + marginAngle1;
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

bool CoilWrapper::delimit_and_compact() {
    auto bobbin = resolve_bobbin();

    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        return delimit_and_compact_rectangular_window();
    }
    else {
        return delimit_and_compact_round_window();
    }
}

bool CoilWrapper::delimit_and_compact_rectangular_window() {
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


        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilHeight = alignedSectionDimensionsPerSection[sectionIndex][1];
                currentCoilWidth += sections[sectionIndex].get_dimensions()[0] / 2;
            }
            else {
                currentCoilHeight -= sections[sectionIndex].get_dimensions()[1] / 2;
                currentCoilWidth = alignedSectionDimensionsPerSection[sectionIndex][0];
            }

            double compactingShiftWidth = sections[sectionIndex].get_coordinates()[0] - currentCoilWidth;
            double compactingShiftHeight = sections[sectionIndex].get_coordinates()[1] - currentCoilHeight;

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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

                                if (bobbinColumnShape == ColumnShape::ROUND || bobbinColumnShape == ColumnShape::RECTANGULAR) {
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
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
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

bool CoilWrapper::delimit_and_compact_round_window() {
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();

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
        auto layers = get_layers_description().value();
        if (get_turns_description()) {
            for (size_t i = 0; i < layers.size(); ++i) {
                if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                    auto turnsInLayer = get_turns_by_layer(layers[i].get_name());
                    auto layerCoordinates = layers[i].get_coordinates();
                    auto section = get_section_by_name(layers[i].get_section().value());

                    double turnDimensionAngle = wound_distance_to_angle(turnsInLayer[0].get_dimensions().value()[1], windingWindows[0].get_radial_height().value() - turnsInLayer[0].get_coordinates()[0]);

                    for (auto& turn : turnsInLayer) {
                        turnDimensionAngle = wound_distance_to_angle(turn.get_dimensions().value()[1], windingWindows[0].get_radial_height().value() - turn.get_coordinates()[0]);
                    }
                    double layerAngle = turnDimensionAngle * turnsInLayer.size();
                    double layerCenterAngle = 0;

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

        auto bobbin = resolve_bobbin();
        auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
        auto windingWindowsRadius = windingWindows[0].get_radial_height().value();
        double windingWindowAngle = windingWindows[0].get_angle().value();
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


        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            if (_windingOrientation == WindingOrientation::OVERLAPPING || _sectionAlignment == CoilAlignment::SPREAD) {
                currentCoilAngle = alignedSectionDimensionsPerSection[sectionIndex][1];
            }
            else {
                currentCoilAngle += sections[sectionIndex].get_dimensions()[1] / 2;
            }

            double compactingShiftAngle = sections[sectionIndex].get_coordinates()[1] - currentCoilAngle;

            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
                    compactingShiftAngle = 0;
                }

            }

            double marginAngle0 = 0;
            double marginAngle1 = 0;

            if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {

                auto partialWinding = sections[sectionIndex].get_partial_windings()[0];  // TODO: Support multiwinding in layers
                auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

                auto wirePerWinding = get_wires();
                double firstLayerMinimumRadius = windingWindowsRadius - (sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2) - wirePerWinding[windingIndex].get_maximum_outer_dimension();
                marginAngle0 = wound_distance_to_angle(sections[sectionIndex].get_margin().value()[0], firstLayerMinimumRadius);
                marginAngle1 = wound_distance_to_angle(sections[sectionIndex].get_margin().value()[1], firstLayerMinimumRadius);
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

                            // angleShiftDueToChangeOfRadius corrects the distance between turns after expanding the radius
                            double angleShiftDueToChangeOfRadius = 0;
                            double angleBetweenTurns = 0;
                            double anglePaddingBetweenTurns = 0;
                            double wireAngle = wound_distance_to_angle(turns[turnIndex].get_dimensions().value()[1], windingWindowsRadius - (turns[turnIndex].get_coordinates()[0]));
                            if (turnInThisLayerIndex == 0) {
                                if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                                    switch (_turnsAlignment) {
                                        case CoilAlignment::INNER_OR_TOP:
                                            angleShiftDueToChangeOfRadius = -(turns[turnIndex].get_coordinates()[1] - wireAngle / 2 - marginAngle0);
                                            break;

                                        case CoilAlignment::OUTER_OR_BOTTOM:
                                            angleShiftDueToChangeOfRadius = (layers[layerIndex].get_coordinates()[1] - layers[layerIndex].get_dimensions()[1] / 2 + wireAngle / 2) - turns[turnIndex].get_coordinates()[1];
                                            break;

                                        case CoilAlignment::CENTERED:
                                            angleShiftDueToChangeOfRadius = (layers[layerIndex].get_coordinates()[1] - layers[layerIndex].get_dimensions()[1] / 2 + wireAngle / 2) - turns[turnIndex].get_coordinates()[1];
                                            break;

                                        case CoilAlignment::SPREAD:
                                            angleShiftDueToChangeOfRadius = 0;
                                            break;

                                        default:

                                            throw std::runtime_error("No such section alignment");
                                    }
                                }
                                else {
                                }
                            }
                            else {
                                if (_windingOrientation == WindingOrientation::OVERLAPPING) {
                                    angleBetweenTurns = turns[turnIndex].get_coordinates()[1] - turns[turnIndex - 1].get_coordinates()[1];
                                    anglePaddingBetweenTurns = turns[turnIndex].get_coordinates()[1] - turns[turnIndex - 1].get_coordinates()[1];
                                    angleShiftDueToChangeOfRadius = -angleBetweenTurns + wireAngle - compactingShiftAngle;
                                    switch (_turnsAlignment) {
                                        case CoilAlignment::INNER_OR_TOP:
                                        case CoilAlignment::OUTER_OR_BOTTOM:
                                        case CoilAlignment::CENTERED:
                                            anglePaddingBetweenTurns = 0;
                                            break;

                                        case CoilAlignment::SPREAD:
                                            anglePaddingBetweenTurns = turns[turnIndex].get_coordinates()[1] - turns[turnIndex - 1].get_coordinates()[1];
                                            angleShiftDueToChangeOfRadius = 0;
                                            break;

                                        default:

                                            throw std::runtime_error("No such section alignment");
                                    }
                                }
                                else {
                                }
                            }

                            turns[turnIndex].set_coordinates(std::vector<double>({
                                turns[turnIndex].get_coordinates()[0],
                                turns[turnIndex].get_coordinates()[1] - compactingShiftAngle
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

                            turnInThisLayerIndex++;
                        }
                    }
                }
                }
            if (_windingOrientation == WindingOrientation::OVERLAPPING) {
            }
            else {
                currentCoilAngle += sections[sectionIndex].get_dimensions()[1] / 2 + paddingAmongSectionAngle;
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

    // // Add extra margin for support if required
    // auto settings = Settings::GetInstance();
    // bool fillCoilSectionsWithMarginTape = settings->get_coil_fill_sections_with_margin_tape();


    // if (fillCoilSectionsWithMarginTape) {
    //     auto bobbin = resolve_bobbin();
    //     auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    //     auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    //     double windingWindowHeight = windingWindowDimensions[1];
    //     double windingWindowWidth = windingWindowDimensions[0];
    //     auto sections = get_sections_description().value();
    //     for (size_t i = 0; i < sections.size(); ++i) {
    //         if (sections[i].get_type() == ElectricalType::CONDUCTION) {
    //             auto sectionOrientation = bobbin.get_winding_window_sections_orientation(0);
    //             if (sectionOrientation == WindingOrientation::OVERLAPPING) {
    //                 auto topSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[1] + windingWindowHeight / 2) - (sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2));
    //                 auto bottomSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[1] - windingWindowHeight / 2) - (sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2));
    //                 sections[i].set_margin(std::vector<double>{topSpaceBetweenSectionAndBobbin, bottomSpaceBetweenSectionAndBobbin});
    //             }
    //             else if (sectionOrientation == WindingOrientation::CONTIGUOUS) {
    //                 auto innerSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[0] - windingWindowWidth / 2) - (sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2));
    //                 auto outerSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[0] + windingWindowWidth / 2) - (sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2));
    //                 sections[i].set_margin(std::vector<double>{innerSpaceBetweenSectionAndBobbin, outerSpaceBetweenSectionAndBobbin});
    //             }
    //         }
    //     }
    //     set_sections_description(sections);
    // }

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
    auto settings = Settings::GetInstance();

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
                    }
                    else {
                        currentSpace += sections[sectionIndex].get_dimensions()[1];
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
    // set_turns_description(std::nullopt);
    if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
        wind_by_turns();
        if (delimitAndCompact) {
            delimit_and_compact();
        }
    }
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
 
