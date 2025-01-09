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

        auto dx1 = currentTurn.get_dimensions().value()[0];
        auto dy1 = currentTurn.get_dimensions().value()[1];
        auto dx2 = potentiallySurroundingTurn.get_dimensions().value()[0];
        auto dy2 = potentiallySurroundingTurn.get_dimensions().value()[1];

        double distance = hypot(x2 - x1, y2 - y1) - std::max(dx1, dy1) / 2 - std::max(dx2, dy2) / 2;

        if (distance > 4e-5) {
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
    auto windingWindowShape = bobbin.get_winding_window_shape(0);
    if (layerOrientation == WindingOrientation::OVERLAPPING) {
        double x1 = firstTurn.get_coordinates()[0];
        double x2 = secondTurn.get_coordinates()[0];

        if (windingWindowShape == WindingWindowShape::ROUND) {
            x1 = coil.cartesian_to_polar(firstTurn.get_coordinates())[0];
            x2 = coil.cartesian_to_polar(secondTurn.get_coordinates())[0];
        }


        for (auto layer : layers) {
            if (layer.get_coordinates()[0] > std::min(x1, x2) && layer.get_coordinates()[0] < std::max(x1, x2)) {
                layersInBetween.push_back(layer);
            }
        }
    }
    else {

        if (windingWindowShape == WindingWindowShape::ROUND) {
            double y1 = coil.cartesian_to_polar(firstTurn.get_coordinates())[1];
            double y2 = coil.cartesian_to_polar(secondTurn.get_coordinates())[1];
            if (y1 < 90 && y2 > 270) {
                for (auto layer : layers) {
                    if (layer.get_coordinates()[1] > y2) {
                        layersInBetween.push_back(layer);
                    }
                }
            }
            else if (y2 < 90 && y1 > 270) {
                for (auto layer : layers) {
                    if (layer.get_coordinates()[1] > y1) {
                        layersInBetween.push_back(layer);
                    }
                }
            }
            else {
                for (auto layer : layers) {
                    if (layer.get_coordinates()[1] > std::min(y1, y2) && layer.get_coordinates()[1] < std::max(y1, y2)) {
                        layersInBetween.push_back(layer);
                    }
                }
            }
        }
        else {
            double y1 = firstTurn.get_coordinates()[1];
            double y2 = secondTurn.get_coordinates()[1];
            for (auto layer : layers) {
                if (layer.get_coordinates()[1] > std::min(y1, y2) && layer.get_coordinates()[1] < std::max(y1, y2)) {
                    layersInBetween.push_back(layer);
                }
            }
        }

    }

    return layersInBetween;
}

StrayCapacitanceOutput StrayCapacitance::calculate_voltages_per_turn(CoilWrapper coil, OperatingPoint operatingPoint) {
    std::map<std::string, double> voltageRmsPerWinding;
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
    }
    return StrayCapacitance::calculate_voltages_per_turn(coil, voltageRmsPerWinding);
}

StrayCapacitanceOutput StrayCapacitance::calculate_voltages_per_turn(CoilWrapper coil, std::map<std::string, double> voltageRmsPerWinding) {
    if (!coil.get_turns_description()) {
        throw std::invalid_argument("Missing turns description");
    }
    auto turns = coil.get_turns_description().value();
    std::map<std::string, int64_t> numberTurnsPerWinding;
    for (size_t windingIndex = 0; windingIndex <  coil.get_functional_description().size(); ++windingIndex) {
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
        auto turnWinding = turn.get_winding();
        auto turnParallel = turn.get_parallel();

        double voltageDividerCenterThisTurn;
        double voltageDividerStartThisTurn;
        double voltageDividerEndThisTurn;
        if (numberTurnsPerWinding[turnWinding] > 1) {
            voltageDividerCenterThisTurn = (double(numberTurnsPerWinding[turnWinding] - 1) - turnIndexPerWindingPerParallel[turnWinding][turnParallel]) / (numberTurnsPerWinding[turnWinding] - 1);
            voltageDividerStartThisTurn = (double(numberTurnsPerWinding[turnWinding]) - turnIndexPerWindingPerParallel[turnWinding][turnParallel]) / numberTurnsPerWinding[turnWinding];
            voltageDividerEndThisTurn = (double(numberTurnsPerWinding[turnWinding]) - turnIndexPerWindingPerParallel[turnWinding][turnParallel] - 1) / numberTurnsPerWinding[turnWinding];
        }
        else  {
            voltageDividerCenterThisTurn = 0.5;
            voltageDividerStartThisTurn = 1;
            voltageDividerEndThisTurn = 0;
        }
        voltageDividerStartPerTurn.push_back(voltageDividerStartThisTurn);
        voltageDividerEndPerTurn.push_back(voltageDividerEndThisTurn);
        double voltage = voltageRmsPerWinding[turnWinding] * voltageDividerCenterThisTurn;
        if (std::isnan(voltage)) {
            throw std::runtime_error("voltage cannot be nan");
        }
        if (std::isinf(voltage)) {
            throw std::runtime_error("voltage cannot be inf");
        }
        voltagePerTurn.push_back(voltageRmsPerWinding[turnWinding] * voltageDividerCenterThisTurn);
        turnIndexPerWindingPerParallel[turnWinding][turnParallel]++;
    }

    strayCapacitanceOutput.set_voltage_divider_start_per_turn(voltageDividerStartPerTurn);
    strayCapacitanceOutput.set_voltage_divider_end_per_turn(voltageDividerEndPerTurn);
    strayCapacitanceOutput.set_voltage_per_turn(voltagePerTurn);

    return strayCapacitanceOutput;
}

double get_effective_relative_permittivity(double firstThickness, double firstRelativePermittivity, double secondThickness, double secondRelativePermittivity) {
    return firstRelativePermittivity * secondRelativePermittivity * (firstThickness + secondThickness) / (firstThickness * secondRelativePermittivity + secondThickness * firstRelativePermittivity);
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
    else if (modelName == StrayCapacitanceModels::DUERDOTH) {
        return std::make_shared<StrayCapacitanceDuerdothModel>();
    }
    else if (modelName == StrayCapacitanceModels::MASSARINI) {
        return std::make_shared<StrayCapacitanceMassariniModel>();
    }

    else
        throw std::runtime_error("Unknown Stray capacitance model, available options are: {KOCH, ALBACH, DUERDOTH, MASSARINI}");
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
    distanceBetweenTurns = roundFloat(distanceBetweenTurns, 6);
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

double StrayCapacitanceMassariniModel::calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    double Dc = conductingRadius * 2;
    double D0 = (conductingRadius + insulationThickness) * 2;
    double epsilonR = get_effective_relative_permittivity(insulationThickness, epsilonD, distanceThroughAir + distanceThroughLayers, epsilonF);
    double aux0 = 2 * epsilonR + log(D0 / Dc);
    double aux1 = sqrt(log(D0 / Dc) * (2 * epsilonR + log(D0 / Dc)));
    double aux2 = sqrt(2 * epsilonR * log(D0 / Dc) + pow(log(D0 / Dc), 2));

    double C0 = vacuumPermittivity * averageTurnLength * 2 * epsilonR * atan(((-1 + sqrt(3)) * aux0) / ((1 + sqrt(3)) * aux1)) / aux2;

    return C0;
}

double StrayCapacitanceDuerdothModel::calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    double h = distanceThroughAir + distanceThroughLayers;
    double delta = insulationThickness;
    double r0 = conductingRadius;
    double dtt = h + 2 * r0 + 2 * insulationThickness;
    double dPrima = 2 * r0 + h;
    double dEff = dPrima - 2.3 * (r0 + delta) + 0.26 * dtt;
    double epsilonEff = get_effective_relative_permittivity(delta, epsilonD, h, epsilonF);

    double C0 = vacuumPermittivity * epsilonEff * averageTurnLength * 2 * r0 / dEff;

    return C0;
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
    double beta = 1.0 / zeta * (1 + distanceThroughLayersAndAir / (2 * effectiveRelativePermittivity * conductingRadius));
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


std::map<std::string, double> StrayCapacitance::calculate_capacitance_matrix(double energy, double voltageDrop, double relativeTurnsRatio) {
    std::map<std::string, double> result;

    double C0 = energy * 2 / pow(voltageDrop, 2);

    auto gamma1 = -C0 / 6;
    auto gamma2 = -C0 / 6;
    auto gamma3 = C0 / 3;
    auto gamma4 = C0 / 3;
    auto gamma5 = C0 / 6;
    auto gamma6 = C0 / 6;

    result["C11"] = gamma1 + relativeTurnsRatio * (gamma4 + gamma5);
    result["C12"] = -2 * gamma4;
    result["C13"] = 2 * relativeTurnsRatio * gamma5;
    result["C22"] = gamma2 + gamma4 + gamma6;
    result["C23"] = 2 * gamma6;
    result["C33"] = gamma3 + gamma5 + gamma6;


    result["C1"] = gamma1 + relativeTurnsRatio * gamma2;
    result["C2"] = gamma5 + gamma6;
    result["C3"] = gamma3;

    return result;
}


std::map<std::pair<std::string, std::string>, double> StrayCapacitance::calculate_capacitance_among_windings(CoilWrapper coil) {
    auto capacitanceAmongTurns = calculate_capacitance_among_turns(coil);

    std::map<std::string, double> voltageRmsPerWinding;
    double primaryNumberTurns = coil.get_functional_description()[0].get_number_turns();
    for (auto winding : coil.get_functional_description()) {
        double turnsRatio = primaryNumberTurns /  winding.get_number_turns();
        voltageRmsPerWinding[winding.get_name()] = 10.0 / turnsRatio;
    }


    auto voltagesPerTurn = StrayCapacitance::calculate_voltages_per_turn(coil, voltageRmsPerWinding).get_voltage_per_turn().value();
    auto windings = coil.get_functional_description();
    std::map<std::pair<std::string, std::string>, double> capacitanceMapPerWindings;
    for (auto firstWinding : windings) {
        auto turnsInFirstWinding = coil.get_turns_names_by_winding(firstWinding.get_name());
        auto firstWindingName = firstWinding.get_name();
        double minVoltageInFirstWinding = 1;
        double maxVoltageInFirstWinding = 0;
        double minVoltageInSecondWinding = 1;
        double maxVoltageInSecondWinding = 0;
        for (auto secondWinding : windings) {
            auto windingsKey = std::make_pair(firstWinding.get_name(), secondWinding.get_name());
            if (capacitanceMapPerWindings.contains(windingsKey) || capacitanceMapPerWindings.contains(std::make_pair(secondWinding.get_name(), firstWinding.get_name()))) {
                continue;
            }
            auto secondWindingName = secondWinding.get_name();

            double V3 = 42;
            double V3calculated = 0;

            if (firstWindingName == secondWindingName) {
                V3calculated = 0;
                // capacitanceMapPerWindings[windingsKey] = 0;
                // continue;
            }

            while (V3 != V3calculated) {
                V3 = V3calculated;
                double energyInBetweenTheseWindings = 0;
                double voltageDropInBetweenTheseWindings = 0;
                // double C0 = 0;
                auto turnsInSecondWinding = coil.get_turns_names_by_winding(secondWinding.get_name());
                bool windingAreNotAdjacent = true;
                for (auto turnInFirstWinding : turnsInFirstWinding) {
                    auto firstTurnVoltage = voltagesPerTurn[coil.get_turn_index_by_name(turnInFirstWinding)];
                    minVoltageInFirstWinding = std::min(minVoltageInFirstWinding, firstTurnVoltage);
                    maxVoltageInFirstWinding = std::max(maxVoltageInFirstWinding, firstTurnVoltage);
                    for (auto turnInSecondWinding : turnsInSecondWinding) {
                        auto secondTurnVoltage = voltagesPerTurn[coil.get_turn_index_by_name(turnInSecondWinding)];
                        if (firstWindingName != secondWindingName) {
                            secondTurnVoltage = -secondTurnVoltage;
                        }
                        minVoltageInSecondWinding = std::min(minVoltageInSecondWinding, secondTurnVoltage);
                        maxVoltageInSecondWinding = std::max(maxVoltageInSecondWinding, secondTurnVoltage);
                        auto turnsKey = std::make_pair(turnInFirstWinding, turnInSecondWinding);
                        if (capacitanceAmongTurns.contains(turnsKey)) {
                            windingAreNotAdjacent = false;
                            energyInBetweenTheseWindings += 0.5 * capacitanceAmongTurns[turnsKey] * pow(V3 + firstTurnVoltage - secondTurnVoltage, 2);
                            if (std::isnan(energyInBetweenTheseWindings)) {
                                throw std::runtime_error("Energy cannot be nan");
                            }
                        }
                    }
                }
                if (windingAreNotAdjacent) {
                    capacitanceMapPerWindings[windingsKey] = 0;
                    continue;

                }
                double voltageDropBetweenWindings = maxVoltageInFirstWinding - minVoltageInSecondWinding + V3;
                double relativeTurnsRatio = static_cast<double>(firstWinding.get_number_turns()) / secondWinding.get_number_turns(); 
                auto matrix = calculate_capacitance_matrix(energyInBetweenTheseWindings, voltageDropBetweenWindings, relativeTurnsRatio);
                if (firstWindingName != secondWindingName) {
                    V3calculated = fabs(-(matrix["C13"] * maxVoltageInFirstWinding + matrix["C23"] * fabs(minVoltageInSecondWinding)) / matrix["C33"]);
                }

                capacitanceMapPerWindings[windingsKey] = energyInBetweenTheseWindings * 2 / pow(voltageDropBetweenWindings, 2);
            }
        }
    }
    return capacitanceMapPerWindings;
}


std::map<std::pair<std::string, std::string>, double> StrayCapacitance::calculate_maxwell_capacitance_matrix(CoilWrapper coil) {
    auto capacitanceMapPerWindings = calculate_capacitance_among_windings(coil);
    std::map<std::pair<std::string, std::string>, double> result;
    auto windings = coil.get_functional_description();

    for (auto firstWinding : windings) {
        double capacitanceSum = 0;
        for (auto secondWinding : windings) {
            auto key = std::make_pair(firstWinding.get_name(), secondWinding.get_name());
            auto mirroredKey = std::make_pair(secondWinding.get_name(), firstWinding.get_name());
            if (capacitanceMapPerWindings.contains(key)) {
                capacitanceSum += capacitanceMapPerWindings[key];
                if (firstWinding.get_name() != secondWinding.get_name()) {
                    result[key] = -capacitanceMapPerWindings[key];
                    result[mirroredKey] = -capacitanceMapPerWindings[key];
                }
            }
            else {
                capacitanceSum += capacitanceMapPerWindings[mirroredKey];
                if (firstWinding.get_name() != secondWinding.get_name()) {
                    result[key] = -capacitanceMapPerWindings[mirroredKey];
                    result[mirroredKey] = -capacitanceMapPerWindings[mirroredKey];
                }
            }
        }
        result[std::make_pair(firstWinding.get_name(), firstWinding.get_name())] = capacitanceSum;
    }

    return result;
}


double capacitance_turn_to_turn(double turnDiameter, double wireRadius, double centerSeparation) {
    // According to https://sci-hub.st/https://ieeexplore.ieee.org/document/793378
    double epsilon0 = Constants().vacuumPermittivity;
    double ctt = pow(std::numbers::pi, 2) * turnDiameter * epsilon0 / log(centerSeparation / (2 * wireRadius) + sqrt(pow(centerSeparation / (2 * wireRadius), 2) - 1));
    return ctt;
}

double capacitance_turn_to_shield(double turnDiameter, double wireRadius, double distance) {
    // According to https://sci-hub.st/https://ieeexplore.ieee.org/document/793378
    double epsilon0 = Constants().vacuumPermittivity;
    double cts = 2 * pow(std::numbers::pi, 2) * turnDiameter * epsilon0 / log(distance / wireRadius + sqrt(pow(distance / wireRadius, 2) - 1));
    return cts;
}

double cab(double n, double ctt, double cts) {
    if (n == 2) {
        return ctt + cts / 2;
    }
    else if (n == 3) {
        return ctt / 2 + cts / 2;
    }
    else {
        double cabValue = cab(n - 2, ctt, cts);
        return (cabValue * ctt / 2) / (cabValue  + ctt / 2) + cts / 2;
    }
}

double cas(double n, double ctt, double cts) {
    if (n == 1) {
        return cts;
    }
    else {
        double casValue = cas(n - 1, ctt, cts);
        return (casValue * ctt) / (casValue  + ctt) + cts;
    }
}

double StrayCapacitanceOneLayer::calculate_capacitance(CoilWrapper coil) {
    // Baed on https://sci-hub.st/https://ieeexplore.ieee.org/document/793378
    double numberTurns = coil.get_functional_description()[0].get_number_turns();
    auto wireRadius = coil.resolve_wire(0).get_maximum_conducting_width() / 2;
    double distanceTurnsToCore = coil.resolve_bobbin().get_processed_description()->get_column_thickness() + coil.resolve_wire(0).get_maximum_outer_width() / 2;
    double turnDiameter = 2 * std::numbers::pi * (coil.resolve_bobbin().get_processed_description()->get_column_width().value() + wireRadius);
    double centerSeparation = coil.resolve_wire(0).get_maximum_outer_width();
    if (coil.get_turns_description()) {
        if (coil.get_turns_description().value().size() > 1) {
            auto firstTurn = coil.get_turns_description().value()[0];
            auto secondTurn = coil.get_turns_description().value()[1];
            centerSeparation = sqrt(pow(firstTurn.get_coordinates()[0] - secondTurn.get_coordinates()[0], 2) + pow(firstTurn.get_coordinates()[1] - secondTurn.get_coordinates()[1], 2));
        }
    }

    double ctt = capacitance_turn_to_turn(turnDiameter, wireRadius, centerSeparation);
    double cts = capacitance_turn_to_shield(turnDiameter, wireRadius, distanceTurnsToCore);
    double C2;
    double casValue = cas(numberTurns, ctt, cts);

    if (std::isnan(casValue)) {
        throw std::runtime_error("capacitance cannot be NaN");
    }

    if (numberTurns > 1) {
        double cabValue = cab(numberTurns, ctt, cts);
        C2 = 2 * cabValue * casValue / (4 * cabValue - casValue);
        double C1 = cabValue - cabValue * casValue / (4 * cabValue - casValue);
        // std::cout << "C2: " << C2 << std::endl;
        // std::cout << "C1: " << C1 << std::endl;

        C2 = C2 * 2;

        if (C1 > 1e-13) {
            C2 =  1.0 / (1.0 / C2 + 1.0 / C1);
        }

    }
    else {
        C2 = casValue;
    }


    auto capacitance = C2;
    if (coil.get_layers_description()) {
        capacitance *= coil.get_layers_by_winding_index(0).size();
    }

    return capacitance;
}

} // namespace OpenMagnetics
