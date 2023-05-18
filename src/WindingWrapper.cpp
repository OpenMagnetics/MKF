#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include "Constants.h"
#include "Utils.h"
#include "WindingWrapper.h"
#include "json.hpp"

#include <magic_enum.hpp>
#include "../tests/TestingUtils.h"

using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

namespace OpenMagnetics {
std::vector<bool> WindingWrapper::wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, uint8_t numberSlots) {
    std::vector<bool> windByConsecutiveTurns;
    for (size_t i = 0; i < numberTurns.size(); ++i) {
        if (numberTurns[i] == numberSlots) {
            windByConsecutiveTurns.push_back(false);
            continue;
        }
        if (numberParallels[i] == numberSlots) {
            windByConsecutiveTurns.push_back(true);
            continue;
        }
        if (numberParallels[i] % numberSlots == 0) {
            windByConsecutiveTurns.push_back(true);
            continue;
        }
        if (numberTurns[i] % numberSlots == 0) {
            windByConsecutiveTurns.push_back(false);
            continue;
        }
        windByConsecutiveTurns.push_back(true);
    }
    return windByConsecutiveTurns;
}

bool WindingWrapper::wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, uint8_t numberSlots) {
    if (numberTurns == numberSlots) {
        return false;
    }
    if (numberParallels == numberSlots) {
        return true;
    }
    if (numberParallels % numberSlots == 0) {
        return true;
    }
    if (numberTurns % numberSlots == 0) {
        return false;
    }
    return true;
}

std::vector<uint64_t> WindingWrapper::get_number_turns() {
    std::vector<uint64_t> numberTurns;
    for (auto & winding : get_functional_description()) {
        numberTurns.push_back(winding.get_number_turns());
    }
    return numberTurns;
}

std::vector<Layer> WindingWrapper::get_layers_by_section(std::string sectionName) {
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

std::vector<Turn> WindingWrapper::get_turns_by_layer(std::string layerName) {
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

std::vector<uint64_t> WindingWrapper::get_number_parallels() {
    std::vector<uint64_t> numberParallels;
    for (auto & winding : get_functional_description()) {
        numberParallels.push_back(winding.get_number_parallels());
    }
    return numberParallels;
}


WindingFunctionalDescription WindingWrapper::get_winding_by_name(std::string name) {
    for (auto& windingFunctionalDescription : get_functional_description()) {
        if (windingFunctionalDescription.get_name() == name) {
            return windingFunctionalDescription;
        }
    }
    throw std::runtime_error("No such a winding name: " + name);
}

size_t WindingWrapper::get_winding_index_by_name(std::string name) {
    for (size_t i=0; i<get_functional_description().size(); ++i) {
        if (get_functional_description()[i].get_name() == name) {
            return i;
        }
    }
    throw std::runtime_error("No such a winding name: " + name);
}

std::vector<WireWrapper> WindingWrapper::get_wires() {
    std::vector<WireWrapper> wirePerWinding;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        WireWrapper wire(std::get<WireS>(get_functional_description()[windingIndex].get_wire()));
        wirePerWinding.push_back(wire);
    }
    return wirePerWinding;
}

double WindingWrapper::horizontalFillingFactor(Section section) {
    auto layers = get_layers_by_section(section.get_name());
    double sectionWidth = section.get_dimensions()[0];
    double layersWidth = 0;
    for (auto& layer : layers) {
        layersWidth += layer.get_dimensions()[0];
    }
    return layersWidth / sectionWidth;
}

double WindingWrapper::verticalFillingFactor(Section section) {
    auto layers = get_layers_by_section(section.get_name());
    double sectionHeight = section.get_dimensions()[1];
    double layersHeight = 0;
    for (auto& layer : layers) {
        layersHeight += layer.get_dimensions()[1];
    }
    return layersHeight / sectionHeight;
}

double WindingWrapper::horizontalFillingFactor(Layer layer) {
    auto turns = get_turns_by_layer(layer.get_name());
    double layerWidth = layer.get_dimensions()[0];
    double turnsWidth = 0;
    for (auto& turn : turns) {
        turnsWidth += turn.get_dimensions().value()[0];
    }
    return turnsWidth / layerWidth;
}

double WindingWrapper::verticalFillingFactor(Layer layer) {
    auto turns = get_turns_by_layer(layer.get_name());
    double layerHeight = layer.get_dimensions()[1];
    double turnsHeight = 0;
    for (auto& turn : turns) {
        turnsHeight += turn.get_dimensions().value()[1];
    }
    return turnsHeight / layerHeight;
}

std::pair<uint64_t, std::vector<double>> get_parallels_proportions(size_t slotIndex, uint8_t slots, uint64_t numberTurns, uint64_t numberParallels, 
                                                                   std::vector<double> remainingParallelsProportion, bool windByConsecutiveTurns,
                                                                   std::vector<double> totalParallelsProportion) {
    uint64_t physicalTurnsThisSlot = 0;
    std::vector<double> slotParallelsProportion(numberParallels, 0);
    if (windByConsecutiveTurns) {
        size_t remainingPhysicalTurns = 0;
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            remainingPhysicalTurns += round(remainingParallelsProportion[parallelIndex] * numberTurns);
        }
        physicalTurnsThisSlot = std::min(remainingPhysicalTurns, uint64_t(ceil(double(remainingPhysicalTurns) / (slots - slotIndex))));
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

double get_area_used_in_wires(WireS wire, uint64_t physicalTurns) {
        if (wire.get_type() == "round") {
            double wireDiameter = resolve_dimensional_values<DimensionalValues::NOMINAL>(wire.get_outer_diameter().value());
            return physicalTurns * pow(wireDiameter, 2);
        }
        else {
            double wireWidth = resolve_dimensional_values<DimensionalValues::NOMINAL>(wire.get_outer_width().value());
            double wireHeight = resolve_dimensional_values<DimensionalValues::NOMINAL>(wire.get_outer_height().value());
            return physicalTurns * wireWidth * wireHeight;
        }
}


void WindingWrapper::wind_by_sections() {
    std::vector<Section> sectionsDescription;
    auto windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(), get_number_parallels(), _interleavingLevel);
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
    double currentSectionCenterWidth = 0;
    double currentSectionCenterHeight = 0;
    if (_windingOrientation == OrientationEnum::HORIZONTAL) {
        interleavedWidth = roundFloat(windingWindowWidth / _interleavingLevel / get_functional_description().size(), 9);
        interleavedHeight = windingWindowHeight;
        currentSectionCenterWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + interleavedWidth / 2;
        currentSectionCenterHeight = windingWindows[0].get_coordinates().value()[1];
    }
    else if (_windingOrientation == OrientationEnum::VERTICAL) {
        interleavedWidth = windingWindowWidth;
        interleavedHeight = roundFloat(windingWindowHeight / _interleavingLevel / get_functional_description().size(), 9);
        currentSectionCenterWidth = windingWindows[0].get_coordinates().value()[0];
        currentSectionCenterHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - interleavedHeight / 2;
    }
    else {
        throw std::runtime_error("Toroidal windings not implemented");
    }

    for (size_t sectionIndex = 0; sectionIndex < _interleavingLevel; ++sectionIndex) {

        for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
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
            section.set_layers_orientation(_layersOrientation);
            section.set_dimensions(std::vector<double>{interleavedWidth, interleavedHeight});
            section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / (interleavedWidth * interleavedHeight));
            sectionsDescription.push_back(section);

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                remainingParallelsProportion[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
            }
            if (_windingOrientation == OrientationEnum::HORIZONTAL) {
                currentSectionCenterWidth += interleavedWidth;
            }
            else if (_windingOrientation == OrientationEnum::VERTICAL) {
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

}

void WindingWrapper::wind_by_layers() {
    if (!get_sections_description()) {
        return;
    }
    auto wirePerWinding = get_wires();

    auto sections = get_sections_description().value();

    std::vector<Layer> layers;
    for (auto& section : sections) {
        uint64_t maximumNumberLayersFittingInSection;
        uint64_t maximumNumberPhysicalTurnsPerLayer;
        uint64_t minimumNumberLayerNeeded;
        uint64_t numberLayers;
        uint64_t physicalTurnsInSection = 0;
        double currentLayerCenterWidth;
        double currentLayerCenterHeight;
        if (_layersOrientation == OrientationEnum::VERTICAL) {
            currentLayerCenterWidth = roundFloat(section.get_coordinates()[0] - section.get_dimensions()[0] / 2, 9);
            currentLayerCenterHeight = roundFloat(section.get_coordinates()[1], 9);
        } else {
            currentLayerCenterWidth = roundFloat(section.get_coordinates()[0], 9);
            currentLayerCenterHeight = roundFloat(section.get_coordinates()[1] + section.get_dimensions()[1] / 2, 9);
        }
        double layerWidth = 0;
        double layerHeight = 0;
        auto remainingParallelsProportionInSection = section.get_partial_windings()[0].get_parallels_proportion();
        auto totalParallelsProportionInSection = section.get_partial_windings()[0].get_parallels_proportion();
        if (section.get_partial_windings().size() > 1) {
            throw std::runtime_error("More than one winding per layer not supported yet");
        }
        auto partialWinding = section.get_partial_windings()[0];  // TODO: Support multiwinding in layers
        auto winding = get_winding_by_name(partialWinding.get_winding());
        auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            physicalTurnsInSection += round(remainingParallelsProportionInSection[parallelIndex] * get_number_turns(windingIndex));
        }

        if (wirePerWinding[windingIndex].get_type() == "round") {
            double wireDiameter = resolve_dimensional_values<DimensionalValues::NOMINAL>(wirePerWinding[windingIndex].get_outer_diameter().value());
            maximumNumberLayersFittingInSection = section.get_dimensions()[0] / wireDiameter;
            maximumNumberPhysicalTurnsPerLayer = floor(section.get_dimensions()[1] / wireDiameter);
            if (_layersOrientation == OrientationEnum::VERTICAL) {
                layerWidth = wireDiameter;
                layerHeight = section.get_dimensions()[1];
            } else {
                layerWidth = section.get_dimensions()[0];
                layerHeight = wireDiameter;
            }
        }
        else {
            double wireWidth = resolve_dimensional_values<DimensionalValues::NOMINAL>(wirePerWinding[windingIndex].get_outer_width().value());
            double wireHeight = resolve_dimensional_values<DimensionalValues::NOMINAL>(wirePerWinding[windingIndex].get_outer_height().value());
            if (_layersOrientation == OrientationEnum::VERTICAL) {
                maximumNumberLayersFittingInSection = section.get_dimensions()[0] / wireWidth;
                maximumNumberPhysicalTurnsPerLayer = floor(section.get_dimensions()[1] / wireHeight);
                layerWidth = wireWidth;
                layerHeight = section.get_dimensions()[1];
            } else {
                maximumNumberLayersFittingInSection = section.get_dimensions()[1] / wireHeight;
                maximumNumberPhysicalTurnsPerLayer = floor(section.get_dimensions()[0] / wireWidth);
                layerWidth = section.get_dimensions()[0];
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


        bool windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberLayers);
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
            layer.set_section(section.get_name());
            layer.set_type(LayersDescriptionType::WIRING);
            layer.set_name(section.get_name() +  " layer " + std::to_string(layerIndex));
            layer.set_orientation(_layersOrientation);
            layer.set_turns_alignment(_turnsAlignment);
            layer.set_dimensions(std::vector<double>{layerWidth, layerHeight});
            layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
            layer.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisLayer) / (layerWidth * layerHeight));
            layers.push_back(layer);

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                remainingParallelsProportionInSection[parallelIndex] -= layerParallelsProportion[parallelIndex];
            }

            if (_layersOrientation == OrientationEnum::HORIZONTAL) {
                currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight, 9);
            }
            else if (_layersOrientation == OrientationEnum::VERTICAL) {
                currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth, 9);
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
            }
        }
    }
    set_layers_description(layers);
}

void WindingWrapper::wind_by_turns() {
    if (!get_layers_description()) {
        return;
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
        auto physicalTurns = get_number_turns(windingIndex) * get_number_parallels(windingIndex);
        auto alignment = layer.get_turns_alignment().value();
        if (wirePerWinding[windingIndex].get_type() == "round") {
            wireWidth = resolve_dimensional_values<DimensionalValues::NOMINAL>(wirePerWinding[windingIndex].get_outer_diameter().value());
            wireHeight = resolve_dimensional_values<DimensionalValues::NOMINAL>(wirePerWinding[windingIndex].get_outer_diameter().value());
        }
        else {
            wireWidth = resolve_dimensional_values<DimensionalValues::NOMINAL>(wirePerWinding[windingIndex].get_outer_width().value());
            wireHeight = resolve_dimensional_values<DimensionalValues::NOMINAL>(wirePerWinding[windingIndex].get_outer_height().value());
        }
        if (_layersOrientation == OrientationEnum::VERTICAL) {
            totalLayerWidth = layer.get_dimensions()[0];
            totalLayerHeight = roundFloat(physicalTurns * wireHeight, 9);
            currentTurnWidthIncrement = 0;
            currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0], 9);
            switch (alignment) {
                case TurnsAlignmentEnum::CENTERED:
                    currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + totalLayerHeight / 2 - wireHeight / 2, 9);
                    currentTurnHeightIncrement = wireHeight;
                    break;

                case TurnsAlignmentEnum::INNER_OR_TOP:
                    currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - wireHeight / 2, 9);
                    currentTurnHeightIncrement = wireHeight;
                    break;

                case TurnsAlignmentEnum::OUTER_OR_BOTTOM:
                    currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 + totalLayerHeight - wireHeight / 2, 9);
                    currentTurnHeightIncrement = wireHeight;
                    break;

                case TurnsAlignmentEnum::SPREAD:
                    currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - wireHeight / 2, 9);
                    currentTurnHeightIncrement = roundFloat(layer.get_dimensions()[1] / (physicalTurns - 1), 9);
                    break;
            }

        } 
        else {
            totalLayerWidth = roundFloat(physicalTurns * wireWidth, 9);
            totalLayerHeight = layer.get_dimensions()[1];
            currentTurnHeightIncrement = 0;
            currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1], 9);
            switch (alignment) {
                case TurnsAlignmentEnum::CENTERED:
                    currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] + totalLayerWidth / 2 - wireWidth / 2, 9);
                    currentTurnWidthIncrement = wireWidth;
                    break;

                case TurnsAlignmentEnum::INNER_OR_TOP:
                    currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2 - wireWidth / 2, 9);
                    currentTurnWidthIncrement = wireWidth;
                    break;

                case TurnsAlignmentEnum::OUTER_OR_BOTTOM:
                    currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + totalLayerWidth - wireWidth / 2, 9);
                    currentTurnWidthIncrement = wireWidth;
                    break;

                case TurnsAlignmentEnum::SPREAD:
                    currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2 - wireWidth / 2, 9);
                    currentTurnWidthIncrement = roundFloat(layer.get_dimensions()[0] / (physicalTurns - 1), 9);
                    break;
            }
        }

        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            int64_t numberTurns = round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
            for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                Turn turn;
                turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                turn.set_layer(layer.get_name());
                if (bobbinColumnShape == ColumnShape::ROUND) {
                    turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth);
                }
                else if (bobbinColumnShape == ColumnShape::RECTANGULAR) {
                    double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                    turn.set_length(2 * bobbinColumnDepth + 2 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                }
                else {
                    throw std::runtime_error("only round or rectangular columns supported for bobbins");
                }
                turn.set_name("Parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(turnIndex));
                turn.set_orientation(TurnOrientation::CLOCKWISE);
                turn.set_parallel(parallelIndex);
                turn.set_section(layer.get_section().value());
                turn.set_winding(partialWinding.get_winding());
                turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});

                turns.push_back(turn);
                currentTurnCenterWidth += currentTurnWidthIncrement;
                currentTurnCenterHeight += currentTurnHeightIncrement;
            }
        }

    }

    set_turns_description(turns);
}

} // namespace OpenMagnetics
