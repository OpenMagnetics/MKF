#include "StrayCapacitance.h"
#include "Settings.h"
#include "Defaults.h"
#include "Constants.h"
#include "Models.h"
#include "MAS.hpp"
#include "Utils.h"
#include "json.hpp"
#include <matplot/matplot.h>
#include <cfloat>
#include <set>

namespace OpenMagnetics {

    std::vector<Turn> StrayCapacitance::get_surrounding_turns(Turn currentTurn, std::vector<Turn> turnsDescription) {
        std::vector<Turn> surroundingTurns;
        for (auto potentiallySurroundingTurn : turnsDescription) {
            auto factor = Defaults().overlappingFactorSurroundingTurns;
            auto x1 = currentTurn.get_coordinates()[0];
            auto y1 = currentTurn.get_coordinates()[1];
            auto x2 = potentiallySurroundingTurn.get_coordinates()[0];
            auto y2 = potentiallySurroundingTurn.get_coordinates()[1];
            if (x1 == x2 && y1 == y2) {
                continue;
            }

            double maximumDimensionOf12 = (std::max(potentiallySurroundingTurn.get_dimensions().value()[0], potentiallySurroundingTurn.get_dimensions().value()[1]) + 
                                                   std::max(currentTurn.get_dimensions().value()[0], currentTurn.get_dimensions().value()[1])) / 2;
            bool thereIsTurnBetween12 = false;
            for (auto potentiallyCollidingTurn : turnsDescription) {
                auto x0 = potentiallyCollidingTurn.get_coordinates()[0];
                auto y0 = potentiallyCollidingTurn.get_coordinates()[1];
                auto dx0 = potentiallyCollidingTurn.get_dimensions().value()[0];
                auto dy0 = potentiallyCollidingTurn.get_dimensions().value()[1];
                if ((x1 == x0 && y1 == y0) || (x2 == x0 && y2 == y0)) {
                    continue;
                }

                if ((x0 + dx0 / 2 * factor) < std::min(x1, x2)) {
                    continue;
                }
                if ((x0 - dx0 / 2 * factor) > std::max(x1, x2)) {
                    continue;
                }
                if ((y0 + dy0 / 2 * factor) < std::min(y1, y2)) {
                    continue;
                }
                if ((y0 - dy0 / 2 * factor) > std::max(y1, y2)) {
                    continue;
                }

                double maximumDimensionOf0 = std::max(potentiallyCollidingTurn.get_dimensions().value()[0], potentiallyCollidingTurn.get_dimensions().value()[1]);
                auto distanceFrom0toLine12 = fabs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt(pow(y2 - y1, 2) + pow(x2 - x1, 2));
                if (maximumDimensionOf12 / 2 + maximumDimensionOf0 / 2 * factor > distanceFrom0toLine12) {
                    thereIsTurnBetween12 = true;
                    break;
                }
            }

            if (thereIsTurnBetween12) {
                continue;
            }
            else {
                surroundingTurns.push_back(potentiallySurroundingTurn);
            }
        }
        return surroundingTurns;
    }

    std::vector<Layer> StrayCapacitance::get_insulation_layers_between_two_turns(Turn firstTurn, Turn secondTurn, CoilWrapper coil) {
        if (!coil.get_sections_description()) {
            throw std::invalid_argument("Missing sections description");
        }
        if (!coil.get_layers_description()) {
            throw std::invalid_argument("Missing layers description");
        }
        std::vector<Layer> layersInBetween;
        auto layers = coil.get_layers_description_insulation();

        auto bobbin = coil.resolve_bobbin();
        auto layerOrientation = bobbin.get_winding_window_sections_orientation(0);
        if (layerOrientation == WindingOrientation::OVERLAPPING) {
            double x1 = firstTurn.get_coordinates()[0];
            double x2 = secondTurn.get_coordinates()[0];

            for (auto layer : layers) {
                if (layer.get_coordinates()[0] > std::min(x1, x2) && layer.get_coordinates()[0] < std::max(x1, x2)) {
                    layersInBetween.push_back(layer);
                }
            }
        }
        else {
            double y1 = firstTurn.get_coordinates()[1];
            double y2 = secondTurn.get_coordinates()[1];

            // TODO: take into account round winding windows?
            for (auto layer : layers) {
                if (layer.get_coordinates()[1] > std::min(y1, y2) && layer.get_coordinates()[1] < std::max(y1, y2)) {
                    layersInBetween.push_back(layer);
                }
            }
        }

        return layersInBetween;
    }

    StrayCapacitanceOutput StrayCapacitance::calculate_voltages_per_turn(CoilWrapper coil, OperatingPoint operatingPoint) {
        if (!coil.get_turns_description()) {
            throw std::invalid_argument("Missing turns description");
        }
        auto turns = coil.get_turns_description().value();
        std::map<std::string, double> voltageRmsPerWinding;
        std::map<std::string, int64_t> numberTurnsPerWinding;
        for (size_t windingIndex = 0; windingIndex <  coil.get_functional_description().size(); ++windingIndex) {
            if (windingIndex >= operatingPoint.get_excitations_per_winding().size()) {
                throw std::runtime_error("Missing excitation for windingIndex: " + std::to_string(windingIndex));
            }
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_voltage()) {
                throw std::invalid_argument("Missing voltage");
            }
            if (!excitation.get_voltage()->get_processed()) {
                throw std::invalid_argument("Voltage is not processed");
            }
            voltageRmsPerWinding[coil.get_functional_description()[windingIndex].get_name()] = excitation.get_voltage()->get_processed()->get_rms().value();
            numberTurnsPerWinding[coil.get_functional_description()[windingIndex].get_name()] = coil.get_functional_description()[windingIndex].get_number_turns();
        }

        std::map<std::string, std::map<size_t, size_t>> turnIndexPerWindingPerParallel;
            for (auto winding : coil.get_functional_description()) {
                for (size_t parallelIndex = 0; parallelIndex < winding.get_number_parallels(); ++parallelIndex) {
                    turnIndexPerWindingPerParallel[winding.get_name()][parallelIndex] = 0;
                }
            }

        std::vector<double> voltageDividerStartPerTurn;
        std::vector<double> voltageDividerEndPerTurn;
        std::vector<double> voltagePerTurn;
        StrayCapacitanceOutput strayCapacitanceOutput;
        for (auto turn : turns) {
            auto turn_winding = turn.get_winding();
            auto turn_parallel = turn.get_parallel();
            double voltageDividerStartThisTurn = (double(numberTurnsPerWinding[turn_winding]) - turnIndexPerWindingPerParallel[turn_winding][turn_parallel]) / numberTurnsPerWinding[turn_winding];
            double voltageDividerEndThisTurn = (double(numberTurnsPerWinding[turn_winding]) - turnIndexPerWindingPerParallel[turn_winding][turn_parallel] - 1) / numberTurnsPerWinding[turn_winding];
            voltageDividerStartPerTurn.push_back(voltageDividerStartThisTurn);
            voltageDividerEndPerTurn.push_back(voltageDividerEndThisTurn);
            voltagePerTurn.push_back(voltageRmsPerWinding[turn_winding] * voltageDividerStartThisTurn);
            turnIndexPerWindingPerParallel[turn_winding][turn_parallel]++;
        }

        strayCapacitanceOutput.set_voltage_divider_start_per_turn(voltageDividerStartPerTurn);
        strayCapacitanceOutput.set_voltage_divider_end_per_turn(voltageDividerEndPerTurn);
        strayCapacitanceOutput.set_voltage_per_turn(voltagePerTurn);

        return strayCapacitanceOutput;
    }

    double get_effective_relative_permittivity(double firstThickness, double firstrelativePermittivity, double secondThickness, double secondrelativePermittivity) {
        return firstrelativePermittivity * secondrelativePermittivity * (firstThickness + secondThickness) / (firstThickness * secondrelativePermittivity + secondThickness * firstrelativePermittivity);
    }

    double get_wire_insulation_relative_permittivity(WireWrapper wire) {
        if (wire.get_type() != WireType::ROUND) {
            throw std::runtime_error("Other wires not implemented yet, check Albach's book");
        }
        auto conductingRadius = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
        return 2.5 + 0.7 / sqrt(2 * conductingRadius * 1000);
    }


    std::shared_ptr<StrayCapacitanceModel> StrayCapacitanceModel::factory(StrayCapacitanceModels modelName) {
        if (modelName == StrayCapacitanceModels::KOCH) {
            return std::make_shared<StrayCapacitanceKochModel>();
        }
        else if (modelName == StrayCapacitanceModels::ALBACH) {
            return std::make_shared<StrayCapacitanceAlbachModel>();
        }

        else
            throw std::runtime_error("Unknown Stray capacitance model, available options are: {KOCH, ALBACH}");
    }

    std::vector<double> StrayCapacitanceModel::preprocess_data(Turn firstTurn, WireWrapper firstWire, Turn secondTurn, WireWrapper secondWire, CoilWrapper coil) {
        if (firstWire.get_type() != WireType::ROUND || secondWire.get_type() != WireType::ROUND) {
            throw std::runtime_error("Other wires not implemented yet, check Albach's book");
        }

        InsulationMaterialWrapper insulationMaterial = find_insulation_material_by_name(Defaults().defaultInsulationMaterial);

        auto epsilonDFirstWire = get_wire_insulation_relative_permittivity(firstWire);
        auto epsilonDSecondWire = get_wire_insulation_relative_permittivity(secondWire);
        auto epsilonD = (epsilonDFirstWire + epsilonDSecondWire) / 2;

        auto insulationThicknessFirstWire = firstWire.get_coating_thickness();
        auto insulationThicknessSecondWire = secondWire.get_coating_thickness();
        auto insulationThickness = (insulationThicknessFirstWire + insulationThicknessSecondWire) / 2;

        auto conductingDiameterFirstWire = firstWire.get_maximum_conducting_width();
        auto conductingDiameterSecondWire = secondWire.get_maximum_conducting_width();
        auto outerDiameterFirstWire = firstWire.get_maximum_outer_width();
        auto outerDiameterSecondWire = secondWire.get_maximum_outer_width();
        auto conductingRadius = (conductingDiameterFirstWire + conductingDiameterSecondWire) / 2;

        double distanceBetweenTurns = hypot(firstTurn.get_coordinates()[0] - secondTurn.get_coordinates()[0], firstTurn.get_coordinates()[1] - secondTurn.get_coordinates()[1]);
        distanceBetweenTurns -= outerDiameterFirstWire / 2 + outerDiameterSecondWire / 2;
        distanceBetweenTurns = roundFloat(distanceBetweenTurns, 9);
        double distanceThroughLayers = 0;
        std::vector<double> distancesThroughLayers;
        std::vector<double> relativePermittivityLayers;
        double effectiveRelativePermittivityLayers = 1;
                
        std::vector<Layer> insulationLayersInBetween = StrayCapacitance::get_insulation_layers_between_two_turns(firstTurn, secondTurn, coil);

        for (auto layer : insulationLayersInBetween) {
            auto distance = coil.get_insulation_layer_thickness(layer);
            auto relativePermittivity = coil.get_insulation_layer_relative_permittivity(layer);
            distanceThroughLayers += distance;
            distancesThroughLayers.push_back(distance);
            relativePermittivityLayers.push_back(relativePermittivity);
        }

        for (size_t index = 0; index < distancesThroughLayers.size(); ++index) {
            if (index == 0) {
                effectiveRelativePermittivityLayers = relativePermittivityLayers[index];
            }
            else {
                effectiveRelativePermittivityLayers = get_effective_relative_permittivity(distancesThroughLayers[index - 1], effectiveRelativePermittivityLayers, distancesThroughLayers[index], relativePermittivityLayers[index]);
            }
        }
        double distanceThroughAir = distanceBetweenTurns - distanceThroughLayers;

        if (distanceBetweenTurns < 0) {
            throw std::invalid_argument("distanceBetweenTurns cannot be negative");
        }
        auto averageTurnLength = (firstTurn.get_length() + secondTurn.get_length()) / 2;

        return {insulationThickness, averageTurnLength, conductingRadius, distanceThroughLayers, distanceThroughAir, epsilonD, effectiveRelativePermittivityLayers};
    }

    double StrayCapacitanceAlbachModel::calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF) {
        auto vacuumPermittivity = Constants().vacuumPermittivity;

        double effectiveRelativePermittivity = epsilonF;
        double distanceThroughLayersAndAir = distanceThroughAir + distanceThroughLayers;
        if (distanceThroughAir > 0 && distanceThroughLayers > 0) {
            effectiveRelativePermittivity = get_effective_relative_permittivity(distanceThroughLayers, epsilonF, distanceThroughAir, vacuumPermittivity);
        }
        else if (distanceThroughAir > 0 && distanceThroughLayers == 0) {
            effectiveRelativePermittivity = vacuumPermittivity;
        }

        double zeta = 1 - insulationThickness / (epsilonD * conductingRadius);
        double beta = 1.0 / zeta * (1 + distanceThroughLayers / (2 * effectiveRelativePermittivity * conductingRadius));
        double V = beta / sqrt(pow(beta, 2) - 1) * atan(sqrt((beta + 1) / (beta - 1)));
        double Z = 1.0 / (pow(beta, 2) - 1)*((pow(beta, 2) - 2) * V - beta / 2) - std::numbers::pi / 4;
        double Y1 = 1.0 / zeta * (V - std::numbers::pi / 4 + 1.0 / (2 * epsilonD) * pow(distanceThroughLayers / conductingRadius, 2) * Z / zeta);
        double C0 = 2.0 / 3 * vacuumPermittivity * averageTurnLength * Y1;

        if (std::isnan(beta)) {
            throw std::invalid_argument("beta is NAN");
        }

        return C0;
    }

    double StrayCapacitanceKochModel::calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF) {
        auto vacuumPermittivity = Constants().vacuumPermittivity;

        double alpha = 1 - insulationThickness / (epsilonD * conductingRadius);
        double beta = 1.0 / alpha * (1 + distanceThroughLayers / (2 * epsilonF * conductingRadius));
        double V = beta / sqrt(pow(beta, 2) - 1) * atan(sqrt((beta + 1) / (beta - 1))) - std::numbers::pi / 4;
        double Z = beta * (pow(beta, 2) - 2) / pow(pow(beta, 2) - 1, 1.5) * atan(sqrt((beta + 1) / (beta - 1))) - beta / (2 * (pow(beta, 2) - 1))  - std::numbers::pi / 4;
        double C0 = vacuumPermittivity * averageTurnLength / (1 - insulationThickness / (epsilonD * conductingRadius)) * (V + 1.0 / (8 * epsilonD) * pow(2 * insulationThickness / conductingRadius, 2) * Z / (1 - insulationThickness / (epsilonD * conductingRadius)));

        if (std::isnan(beta)) {
            throw std::invalid_argument("beta is NAN");
        }

        return C0;
    }

    double StrayCapacitance::calculate_static_capacitance_between_two_turns(Turn firstTurn, WireWrapper firstWire, Turn secondTurn, WireWrapper secondWire, CoilWrapper coil) {
        auto aux = _model->preprocess_data(firstTurn, firstWire, secondTurn, secondWire, coil);
        double insulationThickness = aux[0];
        double averageTurnLength = aux[1];
        double conductingRadius = aux[2];
        double distanceThroughLayers = aux[3];
        double distanceThroughAir = aux[4];
        double epsilonD = aux[5];
        double epsilonF = aux[6];
        return _model->calculate_static_capacitance_between_two_turns(insulationThickness, averageTurnLength, conductingRadius, distanceThroughLayers, distanceThroughAir, epsilonD, epsilonF);
    }

    std::map<std::pair<std::string, std::string>, double> StrayCapacitance::calculate_capacitance_among_turns(CoilWrapper coil) {
        if (!coil.get_turns_description()) {
            throw std::invalid_argument("Missing turns description");
        }

        std::map<std::pair<std::string, std::string>, double> capacitanceAmongTurns;

        auto turns = coil.get_turns_description().value();
        auto vacuumPermittivity = Constants().vacuumPermittivity;
        auto wirePerWinding = coil.get_wires();

        std::set<std::pair<std::string, std::string>> turnsCombinations;

        for (size_t turnIndex = 0; turnIndex <  turns.size(); ++turnIndex) {
            auto turnWindingIndex = coil.get_winding_index_by_name(turns[turnIndex].get_winding());
            auto turnWire = wirePerWinding[turnWindingIndex];
            auto firstTurnName = turns[turnIndex].get_name();
            auto surroundingTurns = OpenMagnetics::StrayCapacitance::get_surrounding_turns(turns[turnIndex], turns);
            for (auto surroundingTurn : surroundingTurns) {
                auto secondTurnName = surroundingTurn.get_name();
                auto key = std::make_pair(firstTurnName, secondTurnName);
                auto surroundingTurnIndex = coil.get_turn_index_by_name(surroundingTurn.get_name());
                if (turnsCombinations.contains(key) || turnsCombinations.contains(std::make_pair(secondTurnName, firstTurnName))) {
                    continue;
                }
                auto surroundingTurnWindingIndex = coil.get_winding_index_by_name(surroundingTurn.get_winding());
                auto surroundingTurnWire = wirePerWinding[surroundingTurnWindingIndex];
                double capacitance = calculate_static_capacitance_between_two_turns(turns[turnIndex], turnWire, surroundingTurn, surroundingTurnWire, coil);
                capacitanceAmongTurns[key] = capacitance;
                turnsCombinations.insert(key);
            }
        }

        return capacitanceAmongTurns;
    }

    std::map<std::pair<std::string, std::string>, double> StrayCapacitance::calculate_capacitance_among_layers(CoilWrapper coil) {
        if (!coil.get_layers_description()) {
            throw std::invalid_argument("Missing layers description");
        }
        auto capacitanceAmongTurns = calculate_capacitance_among_turns(coil);

        auto layers = coil.get_layers_description_conduction();
        std::map<std::pair<std::string, std::string>, double> capacitanceMapPerLayers;
        for (auto firstLayer : layers) {
            auto turnsInFirstLayer = coil.get_turns_names_by_layer(firstLayer.get_name());
            for (auto secondLayer : layers) {
                auto layersKey = std::make_pair(firstLayer.get_name(), secondLayer.get_name());
                if (capacitanceMapPerLayers.contains(layersKey) || capacitanceMapPerLayers.contains(std::make_pair(secondLayer.get_name(), firstLayer.get_name()))) {
                    continue;
                }

                double capacitanceInBetweenTheseLayers = 0;
                auto turnsInSecondLayer = coil.get_turns_names_by_layer(secondLayer.get_name());
                for (auto turnInFirstLayer : turnsInFirstLayer) {
                    for (auto turnInSecondLayer : turnsInSecondLayer) {
                        auto turnsKey = std::make_pair(turnInFirstLayer, turnInSecondLayer);
                        if (capacitanceAmongTurns.contains(turnsKey)) {
                            capacitanceInBetweenTheseLayers += capacitanceAmongTurns[turnsKey];
                        }
                    }
                }
                capacitanceMapPerLayers[layersKey] = capacitanceInBetweenTheseLayers;
            }

        }

        return capacitanceMapPerLayers;
    }

    std::map<std::pair<std::string, std::string>, double> StrayCapacitance::calculate_capacitance_among_windings(CoilWrapper coil) {
        auto capacitanceAmongLayers = calculate_capacitance_among_layers(coil);

        std::map<std::pair<std::string, std::string>, double> capacitanceMapPerWindings;
        for (auto firstWinding : coil.get_functional_description()) {
            auto layersInFirstWinding = coil.get_layers_names_by_winding(firstWinding.get_name());
            for (auto secondwinding : coil.get_functional_description()) {
                auto windingsKey = std::make_pair(firstWinding.get_name(), secondwinding.get_name());
                if (capacitanceMapPerWindings.contains(windingsKey) || capacitanceMapPerWindings.contains(std::make_pair(secondwinding.get_name(), firstWinding.get_name()))) {
                    continue;
                }

                double capacitanceInBetweenTheseWindings = 0;
                auto layersInSecondWinding = coil.get_layers_names_by_winding(secondwinding.get_name());
                for (auto layerInFirstWinding : layersInFirstWinding) {
                    for (auto layerInSecondWinding : layersInSecondWinding) {
                        auto layersKey = std::make_pair(layerInFirstWinding, layerInSecondWinding);
                        if (capacitanceAmongLayers.contains(layersKey)) {
                            capacitanceInBetweenTheseWindings += capacitanceAmongLayers[layersKey];
                        }
                    }
                }
                capacitanceMapPerWindings[windingsKey] = capacitanceInBetweenTheseWindings;
            }

        }

        return capacitanceMapPerWindings;
    }



} // namespace OpenMagnetics
