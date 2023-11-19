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

    if (!hasSectionsData || !hasLayersData || (!hasTurnsData && are_sections_and_layers_fitting())) {
        if (wind() && delimitAndCompact) {
            return delimit_and_compact();
        }
    }
    return false;
}

bool CoilWrapper::wind() {
    std::string bobbinName = "";
    if (std::holds_alternative<std::string>(get_bobbin())) {
        bobbinName = std::get<std::string>(get_bobbin());
        if (bobbinName != "Dummy") {
            auto bobbinData = OpenMagnetics::find_bobbin_by_name(std::get<std::string>(get_bobbin()));
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
            if (_inputs) {
                calculate_insulation();
            }
            wind_by_sections();
            wind_by_layers();
            if (are_sections_and_layers_fitting()) {
                wind_by_turns();
                return true;
            }

        }
    }
    return false;
}


std::vector<WindingStyle> CoilWrapper::wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, uint8_t numberSlots) {
    std::vector<WindingStyle> windByConsecutiveTurns;
    for (size_t i = 0; i < numberTurns.size(); ++i) {
        if (numberTurns[i] == numberSlots) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
            log("Winding winding " + std::to_string(i) + " by putting together parallels of the same turn, as the number of turns is equal to the number of sections.");
            continue;
        }
        if (numberParallels[i] == numberSlots) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            log("Winding winding " + std::to_string(i) + " by putting together turns of the same parallel, as the number of parallels is equal to the number of sections.");
            continue;
        }
        if (numberParallels[i] % numberSlots == 0) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            log("Winding winding " + std::to_string(i) + " by putting together turns of the same parallel, as the number of parallels is divisible by the number of sections.");
            continue;
        }
        if (numberTurns[i] % numberSlots == 0) {
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

double CoilWrapper::horizontalFillingFactor(Section section) {
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

double CoilWrapper::verticalFillingFactor(Section section) {
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

double CoilWrapper::horizontalFillingFactor(Layer layer) {
    auto turns = get_turns_by_layer(layer.get_name());
    double layerWidth = layer.get_dimensions()[0];
    double turnsWidth = 0;
    for (auto& turn : turns) {
        turnsWidth += turn.get_dimensions().value()[0];
    }
    return turnsWidth / layerWidth;
}

double CoilWrapper::verticalFillingFactor(Layer layer) {
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
            double numberTurnsToAddToCurrentParallel = ceil(numberTurns * totalParallelsProportion[parallelIndex] / slots) ;
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


bool CoilWrapper::calculate_insulation() {
    auto wirePerWinding = get_wires();

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            auto inputs = _inputs.value();
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            double totalVoltageToInsulate = 0;
            double dielectricVoltageToInsulate = 0;

            for (auto& operatingPoint : inputs.get_operating_points()) {
                auto excitationLeftTopWinding = operatingPoint.get_excitations_per_winding()[leftTopWindingIndex];
                auto excitationRightBottomWinding = operatingPoint.get_excitations_per_winding()[rightBottomWindingIndex];
                totalVoltageToInsulate = std::max(totalVoltageToInsulate, excitationLeftTopWinding.get_voltage().value().get_processed().value().get_rms().value() + excitationRightBottomWinding.get_voltage().value().get_processed().value().get_rms().value());
            }

            int timesVoltageNeedsToBeCovered;
            InsulationType neededInsulationType;
            bool insulationTypeSpecified = false;
            if (inputs.get_design_requirements().get_insulation()) {
                if (inputs.get_design_requirements().get_insulation().value().get_insulation_type()) {
                    insulationTypeSpecified = true;
                }
            }

            if (insulationTypeSpecified) {
                neededInsulationType = inputs.get_design_requirements().get_insulation().value().get_insulation_type().value();
                switch (neededInsulationType) {
                    case InsulationType::BASIC:
                    case InsulationType::FUNCTIONAL:
                        timesVoltageNeedsToBeCovered = 1;
                        dielectricVoltageToInsulate = totalVoltageToInsulate;
                        break;
                    case InsulationType::SUPPLEMENTARY:
                        timesVoltageNeedsToBeCovered = 2;
                        dielectricVoltageToInsulate = 2 * totalVoltageToInsulate + 1000;
                        break;
                    case InsulationType::DOUBLE:
                        timesVoltageNeedsToBeCovered = 3;
                        dielectricVoltageToInsulate = 2 * totalVoltageToInsulate + 1000;
                        break;
                    case InsulationType::REINFORCED:
                        timesVoltageNeedsToBeCovered = 3;
                        dielectricVoltageToInsulate = 2 * (2 * totalVoltageToInsulate + 1000);
                        break;
                    default:
                        timesVoltageNeedsToBeCovered = 1;
                        dielectricVoltageToInsulate = 2 * totalVoltageToInsulate + 1000;
                        break;
                }
            }
            else {
                neededInsulationType = InsulationType::FUNCTIONAL;
                timesVoltageNeedsToBeCovered = 1;
                dielectricVoltageToInsulate = totalVoltageToInsulate;
            }

            auto timesVoltageIsCoveredByWireInsulationLeftTopWinding = wireLeftTopWinding.get_equivalent_insulation_layers(dielectricVoltageToInsulate);
            auto timesVoltageIsCoveredByWireInsulationRightBottomWinding = wireRightBottomWinding.get_equivalent_insulation_layers(dielectricVoltageToInsulate);
            std::string neededInsulationTypeString = std::string{magic_enum::enum_name(neededInsulationType)};

            int timesVoltageNotCoveredByWires = timesVoltageNeedsToBeCovered - timesVoltageIsCoveredByWireInsulationLeftTopWinding - timesVoltageIsCoveredByWireInsulationRightBottomWinding;
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            if (timesVoltageNotCoveredByWires > 0) {
                _insulationSectionsLog[windingsMapKey] = "Adding an insulation section, because wires are counting for " + std::to_string(timesVoltageNotCoveredByWires) + " full isolation, and " + neededInsulationTypeString + " needs " + std::to_string(timesVoltageNeedsToBeCovered) + " times.";
            }
            else {
                dielectricVoltageToInsulate = 1;  // Just to have minimum mechanical layer
                _insulationSectionsLog[windingsMapKey] = "No insulation section needed, because wires are enough for covering " + neededInsulationTypeString + " Insulation. Just adding minimum mechanical layer";
            }


            if (insulationMaterialDatabase.empty()) {
                load_databases(true);
            }

            double maxAmbientTemperature = 0;
            for (auto& operatingPoint : inputs.get_operating_points()) {
                 maxAmbientTemperature = std::max(maxAmbientTemperature, operatingPoint.get_conditions().get_ambient_temperature());
            }

            double smallestInsulationThicknessCoveringRemaining = DBL_MAX;
            double chosenMaterialThickness = 0;
            size_t chosenNumberLayers = 0;
            InsulationMaterialWrapper chosenMaterial;


            for (auto& insulationMaterial : insulationMaterialDatabase) {

                if (insulationMaterial.second.get_melting_point()) {
                    if (insulationMaterial.second.get_melting_point().value() < maxAmbientTemperature) {
                        continue;
                    }
                }

                for (auto& aux : insulationMaterial.second.get_available_thicknesses()) {
                    int layersNeeded = ceil(dielectricVoltageToInsulate / aux.second);
                    double totalThicknessNeeded = layersNeeded * aux.first;
                    if (totalThicknessNeeded < smallestInsulationThicknessCoveringRemaining) {
                        smallestInsulationThicknessCoveringRemaining = totalThicknessNeeded;
                        chosenMaterial = insulationMaterial.second;
                        chosenMaterialThickness = aux.first;
                        chosenNumberLayers = layersNeeded;
                    }
                }
            }

            _insulationLayers[windingsMapKey] = std::vector<Layer>();

            for (size_t layerIndex = 0; layerIndex < chosenNumberLayers; ++layerIndex) {
                Layer layer;
                layer.set_partial_windings(std::vector<PartialWinding>{});
                // layer.set_section(section.get_name());
                layer.set_type(ElectricalType::INSULATION);
                layer.set_name("temp");
                layer.set_orientation(_layersOrientation);
                layer.set_turns_alignment(_turnsAlignment);
                if (_layersOrientation == WindingOrientation::VERTICAL) {
                    layer.set_dimensions(std::vector<double>{chosenMaterialThickness, windingWindowHeight});
                }
                else if (_layersOrientation == WindingOrientation::HORIZONTAL) {
                    layer.set_dimensions(std::vector<double>{windingWindowWidth, chosenMaterialThickness});
                }
                // layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layer.set_filling_factor(1);
                _insulationLayers[windingsMapKey].push_back(layer);
            }
            _insulationLayersLog[windingsMapKey] = "Adding " + std::to_string(chosenNumberLayers) + " insulation layers, as we need a thickness of " + std::to_string(smallestInsulationThicknessCoveringRemaining * 1000) + " mm to achieve " + neededInsulationTypeString + " insulation";

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(_layersOrientation);
            section.set_type(ElectricalType::INSULATION);
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                section.set_dimensions(std::vector<double>{smallestInsulationThicknessCoveringRemaining, windingWindowHeight});
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                section.set_dimensions(std::vector<double>{windingWindowWidth, smallestInsulationThicknessCoveringRemaining});
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

bool CoilWrapper::wind_by_sections() {
    auto proportionPerWinding = std::vector<double>(get_functional_description().size(), 1.0 / get_functional_description().size());
    return wind_by_sections(proportionPerWinding);
}

bool CoilWrapper::wind_by_sections(std::vector<double> proportionPerWinding) {
    if (_interleavingLevel <= 0) {
        throw std::runtime_error("Interleaving levels must be greater than 0");
    }
    set_sections_description(std::nullopt);
    std::vector<Section> sectionsDescription;
    auto windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(), get_number_parallels(), _interleavingLevel);
    auto bobbin = resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        throw std::runtime_error("Bobbin not processed");
    }
    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    if (windingWindows.size() > 1) {
        throw std::runtime_error("Bobbins with more than winding window not implemented yet");
    }
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();
    std::vector<std::vector<double>> remainingParallelsProportion;
    auto wirePerWinding = get_wires();
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        remainingParallelsProportion.push_back(std::vector<double>(get_number_parallels(windingIndex), 1));
    }
    double interleavedHeight = 0;
    double interleavedWidth = 0;
    double currentSectionCenterWidth = DBL_MAX;
    double currentSectionCenterHeight = DBL_MAX;

    double totalInsulationWidth = 0;
    double totalInsulationHeight = 0;
    for (size_t sectionIndex = 0; sectionIndex < _interleavingLevel; ++sectionIndex) {
        for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
            if (_inputs && !(sectionIndex == (_interleavingLevel - 1U) && windingIndex == (get_functional_description().size() - 1))) {
                auto nextWindingIndex = (windingIndex + 1) % get_functional_description().size();
                auto windingsMapKey = std::pair<size_t, size_t>{windingIndex, nextWindingIndex}; 
                if (!_insulationSections.contains(windingsMapKey)) {
                    continue;
                }
                if (_windingOrientation == WindingOrientation::HORIZONTAL) {

                    totalInsulationWidth += _insulationSections[windingsMapKey].get_dimensions()[0];
                    totalInsulationHeight = std::max(totalInsulationHeight, _insulationSections[windingsMapKey].get_dimensions()[1]);
                }
                else if (_windingOrientation == WindingOrientation::VERTICAL) {
                    totalInsulationWidth = std::max(totalInsulationWidth, _insulationSections[windingsMapKey].get_dimensions()[0]);
                    totalInsulationHeight += _insulationSections[windingsMapKey].get_dimensions()[1];
                }
            }
        }
    }
    for (size_t sectionIndex = 0; sectionIndex < _interleavingLevel; ++sectionIndex) {
        for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {

            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                interleavedWidth = roundFloat((windingWindowWidth - totalInsulationWidth) * proportionPerWinding[windingIndex] / _interleavingLevel, 9);
                interleavedHeight = windingWindowHeight;

                if (currentSectionCenterWidth == DBL_MAX) {
                    currentSectionCenterWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + interleavedWidth / 2;
                }
                if (currentSectionCenterHeight == DBL_MAX) {
                    currentSectionCenterHeight = windingWindows[0].get_coordinates().value()[1];
                }
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                interleavedWidth = windingWindowWidth;
                interleavedHeight = roundFloat((windingWindowHeight - totalInsulationHeight) * proportionPerWinding[windingIndex] / _interleavingLevel, 9);
                if (currentSectionCenterWidth == DBL_MAX) {
                    currentSectionCenterWidth = windingWindows[0].get_coordinates().value()[0];
                }
                if (currentSectionCenterHeight == DBL_MAX) {
                    currentSectionCenterHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - interleavedHeight / 2;
                }
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
            }

            
            PartialWinding partialWinding;
            Section section;
            partialWinding.set_winding(get_name(windingIndex));

            auto parallels_proportions = get_parallels_proportions(sectionIndex,
                                                                   _interleavingLevel,
                                                                   get_number_turns(windingIndex),
                                                                   get_number_parallels(windingIndex),
                                                                   remainingParallelsProportion[windingIndex],
                                                                   windByConsecutiveTurns[windingIndex],
                                                                   std::vector<double>(get_number_parallels(windingIndex), 1));

            std::vector<double> sectionParallelsProportion = parallels_proportions.second;
            uint64_t physicalTurnsThisSection = parallels_proportions.first;

            partialWinding.set_parallels_proportion(sectionParallelsProportion);
            section.set_name(get_name(windingIndex) +  " section " + std::to_string(sectionIndex));
            section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
            section.set_type(ElectricalType::CONDUCTION);
            section.set_layers_orientation(_layersOrientation);
            section.set_dimensions(std::vector<double>{interleavedWidth, interleavedHeight});
            if (section.get_dimensions()[0] < 0) {
                return false;
                // throw std::runtime_error("Something wrong happened in section dimn 1: " + std::to_string(section.get_dimensions()[0]) +
                //                          " interleavedWidth: " + std::to_string(interleavedWidth) +
                //                          " interleavedHeight: " + std::to_string(interleavedHeight)
                //                          );
            }
            section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});

            if (section.get_coordinates()[0] < -1) {
                return false;
                // throw std::runtime_error("Something wrong happened in section dimn 1: " + std::to_string(section.get_coordinates()[0]) +
                //                          " currentSectionCenterWidth: " + std::to_string(currentSectionCenterWidth) +
                //                          " windingWindows[0].get_coordinates().value()[0]: " + std::to_string(windingWindows[0].get_coordinates().value()[0]) +
                //                          " windingWindows[0].get_width(): " + std::to_string(bool(windingWindows[0].get_width())) +
                //                          " windingWindows[0].get_width().value(): " + std::to_string(windingWindows[0].get_width().value()) +
                //                          " windingWindowWidth: " + std::to_string(windingWindowWidth) +
                //                          " interleavedWidth: " + std::to_string(interleavedWidth) +
                //                          " currentSectionCenterHeight: " + std::to_string(currentSectionCenterHeight)
                //                          );
            }

            section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / (interleavedWidth * interleavedHeight));
            section.set_winding_style(windByConsecutiveTurns[windingIndex]);
            sectionsDescription.push_back(section);

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                remainingParallelsProportion[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
            }

            if (_inputs && !(sectionIndex == (_interleavingLevel - 1U) && windingIndex == (get_functional_description().size() - 1))) {
                auto nextWindingIndex = (windingIndex + 1) % get_functional_description().size();

                auto windingsMapKey = std::pair<size_t, size_t>{windingIndex, nextWindingIndex};
                if (!_insulationSections.contains(windingsMapKey)) {
                    log(_insulationSectionsLog[windingsMapKey]);
                    continue;
                }

                auto insulationSection = _insulationSections[windingsMapKey];


                insulationSection.set_name("Insulation " + get_name(windingIndex) +  " section " + std::to_string(sectionIndex));
                if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth + interleavedWidth / 2 + insulationSection.get_dimensions()[0] / 2,
                                                                          currentSectionCenterHeight,
                                                                          0});
                }
                else if (_windingOrientation == WindingOrientation::VERTICAL) {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth,
                                                                          currentSectionCenterHeight - interleavedHeight / 2 - insulationSection.get_dimensions()[1] / 2,
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

            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentSectionCenterWidth += interleavedWidth;
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                currentSectionCenterHeight -= interleavedHeight;
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

            double currentLayerCenterWidth;
            double currentLayerCenterHeight;
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::VERTICAL) {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerWidth / 2, 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0], 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);
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
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - wireHeight / 2, 9);
                        currentTurnHeightIncrement = roundFloat((layer.get_dimensions()[1] - wireHeight) / (physicalTurnsInLayer - 1), 9);
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
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + wireWidth / 2, 9);
                        currentTurnWidthIncrement = roundFloat((layer.get_dimensions()[0] - wireWidth) / (physicalTurnsInLayer - 1), 9);
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

bool CoilWrapper::delimit_and_compact() {
    // Delimit
    if (get_layers_description()) {
        auto layers = get_layers_description().value();
        if (get_turns_description()) {
            for (size_t i = 0; i < layers.size(); ++i) {
                if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                    auto turnsInLayer = get_turns_by_layer(layers[i].get_name());
                    if (turnsInLayer.size() == 0) {
                        throw std::runtime_error("No turns in layer: " + layers[i].get_name());
                    }
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

        auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
        double windingWindowHeight = windingWindows[0].get_height().value();
        double windingWindowWidth = windingWindows[0].get_width().value();

        if (sections.size() == 0) {
            throw std::runtime_error("No sections in coil");
        }
        double totalSectionsWidth = 0;
        double totalSectionsHeight = 0;
        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                totalSectionsWidth += sections[sectionIndex].get_dimensions()[0];
                if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
                    totalSectionsHeight = std::max(totalSectionsHeight, sections[sectionIndex].get_dimensions()[1]);
                }
            }
            else {
                if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
                    totalSectionsWidth = std::max(totalSectionsWidth, sections[sectionIndex].get_dimensions()[0]);
                }
                totalSectionsHeight += sections[sectionIndex].get_dimensions()[1];
            }
        }

        double currentCoilWidth;
        double currentCoilHeight;
        double paddingAmongSectionWidth = 0;
        double paddingAmongSectionHeight = 0;

        switch (_sectionAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                    switch (_turnsAlignment) {
                        case CoilAlignment::INNER_OR_TOP:
                            currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - totalSectionsHeight / 2;
                            break;
                        case CoilAlignment::OUTER_OR_BOTTOM:
                            currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + totalSectionsHeight / 2;
                            break;
                        case CoilAlignment::CENTERED:
                        case CoilAlignment::SPREAD:
                            currentCoilHeight = 0;
                            break;
                        default:
                            throw std::runtime_error("No such section alignment");
                    }
                }
                else {
                    currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2;
                }
                break;
            case CoilAlignment::OUTER_OR_BOTTOM:
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - totalSectionsWidth;
                if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                    switch (_turnsAlignment) {
                        case CoilAlignment::INNER_OR_TOP:
                            currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - totalSectionsHeight / 2;
                            break;
                        case CoilAlignment::OUTER_OR_BOTTOM:
                            currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + totalSectionsHeight / 2;
                            break;
                        case CoilAlignment::CENTERED:
                        case CoilAlignment::SPREAD:
                            currentCoilHeight = 0;
                            break;
                        default:
                            throw std::runtime_error("No such section alignment");
                    }
                }
                else {
                    currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + totalSectionsHeight;
                }
                break;
            case CoilAlignment::SPREAD:
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                    switch (_turnsAlignment) {
                        case CoilAlignment::INNER_OR_TOP:
                            currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - totalSectionsHeight / 2;
                            break;
                        case CoilAlignment::OUTER_OR_BOTTOM:
                            currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + totalSectionsHeight / 2;
                            break;
                        case CoilAlignment::CENTERED:
                        case CoilAlignment::SPREAD:
                            currentCoilHeight = 0;
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
                }
                break;
            case CoilAlignment::CENTERED:
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 ;
                if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                    switch (_turnsAlignment) {
                        case CoilAlignment::INNER_OR_TOP:
                            currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - totalSectionsHeight / 2;
                            break;
                        case CoilAlignment::OUTER_OR_BOTTOM:
                            currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + totalSectionsHeight / 2;
                            break;
                        case CoilAlignment::CENTERED:
                        case CoilAlignment::SPREAD:
                            currentCoilHeight = 0;
                            break;
                        default:
                            throw std::runtime_error("No such section alignment");
                    }
                }
                else {
                    currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + totalSectionsHeight / 2;
                }
                break;
            default:
                throw std::runtime_error("No such section alignment");

        }

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
                currentCoilWidth += sections[sectionIndex].get_dimensions()[0] / 2;
            }
            else {
                currentCoilHeight -= sections[sectionIndex].get_dimensions()[1] / 2;
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
} // namespace OpenMagnetics
