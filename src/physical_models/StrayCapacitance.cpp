#include "physical_models/StrayCapacitance.h"
#include "support/Settings.h"
#include "Defaults.h"
#include "Constants.h"
#include "Models.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#include <matplot/matplot.h>
#include <cfloat>
#include <set>

namespace OpenMagnetics {

std::vector<std::pair<Turn, size_t>> StrayCapacitance::get_surrounding_turns(Turn currentTurn, std::vector<Turn> turnsDescription) {
    std::vector<std::pair<Turn, size_t>> surroundingTurns;
    for (size_t turnIndex = 0; turnIndex < turnsDescription.size(); ++turnIndex) {
        auto potentiallySurroundingTurn = turnsDescription[turnIndex];
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

        double minimumDimension = std::min(std::max(dx1, dy1), std::max(dx2, dy2));
        double distance = hypot(x2 - x1, y2 - y1) - std::max(dx1, dy1) / 2 - std::max(dx2, dy2) / 2;

        if (distance > minimumDimension / 2) {
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
            surroundingTurns.push_back({potentiallySurroundingTurn, turnIndex});
        }
    }
    return surroundingTurns;
}

double StrayCapacitance::calculate_area_between_two_turns(Turn firstTurn, Turn secondTurn) {
    return std::max(calculate_area_between_two_turns_using_diagonals(firstTurn, secondTurn), calculate_area_between_two_turns_using_vecticals_and_horizontals(firstTurn, secondTurn));
}
double StrayCapacitance::calculate_area_between_two_turns_using_diagonals(Turn firstTurn, Turn secondTurn) {

    double x1 = 0;
    double y1 = 0;
    double x2 = 0;
    double y2 = 0;
    double x3 = 0;
    double y3 = 0;
    double x4 = 0;
    double y4 = 0;
    double sqrt2 = sqrt(2);

    double angleBetweenTurns = atan2(secondTurn.get_coordinates()[1] - firstTurn.get_coordinates()[1], secondTurn.get_coordinates()[0] - firstTurn.get_coordinates()[0]);
    if ((angleBetweenTurns > 0 && angleBetweenTurns < std::numbers::pi / 2) || (angleBetweenTurns < 0 && angleBetweenTurns < -std::numbers::pi / 2)) {
        // Top right quadrant

        if (firstTurn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::RECTANGULAR) {
            x1 = firstTurn.get_coordinates()[0] - firstTurn.get_dimensions().value()[0] / 2;
            y1 = firstTurn.get_coordinates()[1] + firstTurn.get_dimensions().value()[1] / 2;
            x2 = firstTurn.get_coordinates()[0] + firstTurn.get_dimensions().value()[0] / 2;
            y2 = firstTurn.get_coordinates()[1] - firstTurn.get_dimensions().value()[1] / 2;
        }
        else {
            x1 = firstTurn.get_coordinates()[0] - firstTurn.get_dimensions().value()[0] / 2 / sqrt2;
            y1 = firstTurn.get_coordinates()[1] + firstTurn.get_dimensions().value()[1] / 2 / sqrt2;
            x2 = firstTurn.get_coordinates()[0] + firstTurn.get_dimensions().value()[0] / 2 / sqrt2;
            y2 = firstTurn.get_coordinates()[1] - firstTurn.get_dimensions().value()[1] / 2 / sqrt2;
        }
        if (secondTurn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::RECTANGULAR) {
            x3 = secondTurn.get_coordinates()[0] - secondTurn.get_dimensions().value()[0] / 2;
            y3 = secondTurn.get_coordinates()[1] + secondTurn.get_dimensions().value()[1] / 2;
            x4 = secondTurn.get_coordinates()[0] + secondTurn.get_dimensions().value()[0] / 2;
            y4 = secondTurn.get_coordinates()[1] - secondTurn.get_dimensions().value()[1] / 2;
        }
        else {
            x3 = secondTurn.get_coordinates()[0] - secondTurn.get_dimensions().value()[0] / 2 / sqrt2;
            y3 = secondTurn.get_coordinates()[1] + secondTurn.get_dimensions().value()[1] / 2 / sqrt2;
            x4 = secondTurn.get_coordinates()[0] + secondTurn.get_dimensions().value()[0] / 2 / sqrt2;
            y4 = secondTurn.get_coordinates()[1] - secondTurn.get_dimensions().value()[1] / 2 / sqrt2;
        }
    }
    else {
        // Top left quadrant
        if (firstTurn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::RECTANGULAR) {
            x1 = firstTurn.get_coordinates()[0] + firstTurn.get_dimensions().value()[0] / 2;
            y1 = firstTurn.get_coordinates()[1] + firstTurn.get_dimensions().value()[1] / 2;
            x2 = firstTurn.get_coordinates()[0] - firstTurn.get_dimensions().value()[0] / 2;
            y2 = firstTurn.get_coordinates()[1] - firstTurn.get_dimensions().value()[1] / 2;
        }
        else {
            x1 = firstTurn.get_coordinates()[0] + firstTurn.get_dimensions().value()[0] / 2 / sqrt2;
            y1 = firstTurn.get_coordinates()[1] + firstTurn.get_dimensions().value()[1] / 2 / sqrt2;
            x2 = firstTurn.get_coordinates()[0] - firstTurn.get_dimensions().value()[0] / 2 / sqrt2;
            y2 = firstTurn.get_coordinates()[1] - firstTurn.get_dimensions().value()[1] / 2 / sqrt2;
        }
        if (secondTurn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::RECTANGULAR) {
            x3 = secondTurn.get_coordinates()[0] + secondTurn.get_dimensions().value()[0] / 2;
            y3 = secondTurn.get_coordinates()[1] + secondTurn.get_dimensions().value()[1] / 2;
            x4 = secondTurn.get_coordinates()[0] - secondTurn.get_dimensions().value()[0] / 2;
            y4 = secondTurn.get_coordinates()[1] - secondTurn.get_dimensions().value()[1] / 2;
        }
        else {
            x3 = secondTurn.get_coordinates()[0] + secondTurn.get_dimensions().value()[0] / 2 / sqrt2;
            y3 = secondTurn.get_coordinates()[1] + secondTurn.get_dimensions().value()[1] / 2 / sqrt2;
            x4 = secondTurn.get_coordinates()[0] - secondTurn.get_dimensions().value()[0] / 2 / sqrt2;
            y4 = secondTurn.get_coordinates()[1] - secondTurn.get_dimensions().value()[1] / 2 / sqrt2;
        }
    }


    // Helper: Cross product (p2-p1) × (p3-p1)
    auto cross = [](double x1, double y1, double x2, double y2, 
                    double x3, double y3) -> double {
        return (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
    };
    
    // Store points
    double x[4] = {x1, x2, x3, x4};
    double y[4] = {y1, y2, y3, y4};
    
    // Find leftmost point (for convex hull start)
    int start = 0;
    for (int i = 1; i < 4; i++) {
        if (x[i] < x[start] || (x[i] == x[start] && y[i] < y[start])) {
            start = i;
        }
    }
    
    // Build convex hull (Graham scan for 4 points)
    int hull[5];
    int hullSize = 0;
    int p = start, q;
    
    do {
        hull[hullSize++] = p;
        
        // Find next point that makes the largest left turn
        q = (p + 1) % 4;
        for (int i = 0; i < 4; i++) {
            if (i == p || i == q) continue;
            
            double orient = cross(x[p], y[p], x[q], y[q], x[i], y[i]);
            if (orient < 0) {
                q = i;
            }
        }
        p = q;
    } while (p != start && hullSize < 5);
    
    // Close the polygon
    hull[hullSize] = hull[0];
    
    // Shoelace formula on convex hull
    double area = 0.0;
    for (int i = 0; i < hullSize; i++) {
        area += x[hull[i]] * y[hull[i+1]];
        area -= y[hull[i]] * x[hull[i+1]];
    }

    double areaWithoutWires = 0.5 * std::abs(area);


    if (firstTurn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
        double halfCircleArea = std::numbers::pi * pow(firstTurn.get_dimensions().value()[0], 2) / 8;
        areaWithoutWires -= halfCircleArea;
    }
    else {
        double halfRectangularArea = firstTurn.get_dimensions().value()[0] * firstTurn.get_dimensions().value()[1] / 2;
        areaWithoutWires -= halfRectangularArea;
    }

    if (secondTurn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
        double halfCircleArea = std::numbers::pi * pow(secondTurn.get_dimensions().value()[0], 2) / 8;
        areaWithoutWires -= halfCircleArea;
    }
    else {
        double halfRectangularArea = secondTurn.get_dimensions().value()[0] * secondTurn.get_dimensions().value()[1] / 2;
        areaWithoutWires -= halfRectangularArea;
    }

    return areaWithoutWires;
}
double StrayCapacitance::calculate_area_between_two_turns_using_vecticals_and_horizontals(Turn firstTurn, Turn secondTurn) {

    double x1 = 0;
    double y1 = 0;
    double x2 = 0;
    double y2 = 0;
    double x3 = 0;
    double y3 = 0;
    double x4 = 0;
    double y4 = 0;

    double angleBetweenTurns = atan2(secondTurn.get_coordinates()[1] - firstTurn.get_coordinates()[1], secondTurn.get_coordinates()[0] - firstTurn.get_coordinates()[0]);
    if ((angleBetweenTurns > std::numbers::pi / 4 && angleBetweenTurns < 3 * std::numbers::pi / 4) || (angleBetweenTurns < -std::numbers::pi / 4 && angleBetweenTurns > -3 * std::numbers::pi / 4)) {
        // Top right quadrant

        x1 = firstTurn.get_coordinates()[0] - firstTurn.get_dimensions().value()[0] / 2;
        y1 = firstTurn.get_coordinates()[1];
        x2 = firstTurn.get_coordinates()[0] + firstTurn.get_dimensions().value()[0] / 2;
        y2 = firstTurn.get_coordinates()[1];

        x3 = secondTurn.get_coordinates()[0] - secondTurn.get_dimensions().value()[0] / 2;
        y3 = secondTurn.get_coordinates()[1];
        x4 = secondTurn.get_coordinates()[0] + secondTurn.get_dimensions().value()[0] / 2;
        y4 = secondTurn.get_coordinates()[1];
    }
    else {
        x1 = firstTurn.get_coordinates()[0];
        y1 = firstTurn.get_coordinates()[1] + firstTurn.get_dimensions().value()[1] / 2;
        x2 = firstTurn.get_coordinates()[0];
        y2 = firstTurn.get_coordinates()[1] - firstTurn.get_dimensions().value()[1] / 2;

        x3 = secondTurn.get_coordinates()[0];
        y3 = secondTurn.get_coordinates()[1] + secondTurn.get_dimensions().value()[1] / 2;
        x4 = secondTurn.get_coordinates()[0];
        y4 = secondTurn.get_coordinates()[1] - secondTurn.get_dimensions().value()[1] / 2;
    }

    // Helper: Cross product (p2-p1) × (p3-p1)
    auto cross = [](double x1, double y1, double x2, double y2, 
                    double x3, double y3) -> double {
        return (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
    };
    
    // Store points
    double x[4] = {x1, x2, x3, x4};
    double y[4] = {y1, y2, y3, y4};
    
    // Find leftmost point (for convex hull start)
    int start = 0;
    for (int i = 1; i < 4; i++) {
        if (x[i] < x[start] || (x[i] == x[start] && y[i] < y[start])) {
            start = i;
        }
    }
    
    // Build convex hull (Graham scan for 4 points)
    int hull[5];
    int hullSize = 0;
    int p = start, q;
    
    do {
        hull[hullSize++] = p;
        
        // Find next point that makes the largest left turn
        q = (p + 1) % 4;
        for (int i = 0; i < 4; i++) {
            if (i == p || i == q) continue;
            
            double orient = cross(x[p], y[p], x[q], y[q], x[i], y[i]);
            if (orient < 0) {
                q = i;
            }
        }
        p = q;
    } while (p != start && hullSize < 5);
    
    // Close the polygon
    hull[hullSize] = hull[0];
    
    // Shoelace formula on convex hull
    double area = 0.0;
    for (int i = 0; i < hullSize; i++) {
        area += x[hull[i]] * y[hull[i+1]];
        area -= y[hull[i]] * x[hull[i+1]];
    }

    double areaWithoutWires = 0.5 * std::abs(area);


    if (firstTurn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
        double halfCircleArea = std::numbers::pi * pow(firstTurn.get_dimensions().value()[0], 2) / 8;
        areaWithoutWires -= halfCircleArea;
    }
    else {
        double halfRectangularArea = firstTurn.get_dimensions().value()[0] * firstTurn.get_dimensions().value()[1] / 2;
        areaWithoutWires -= halfRectangularArea;
    }

    if (secondTurn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
        double halfCircleArea = std::numbers::pi * pow(secondTurn.get_dimensions().value()[0], 2) / 8;
        areaWithoutWires -= halfCircleArea;
    }
    else {
        double halfRectangularArea = secondTurn.get_dimensions().value()[0] * secondTurn.get_dimensions().value()[1] / 2;
        areaWithoutWires -= halfRectangularArea;
    }

    return areaWithoutWires;
}

std::vector<Layer> StrayCapacitance::get_insulation_layers_between_two_turns(Turn firstTurn, Turn secondTurn, Coil coil) {
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
            if (layer.get_dimensions()[0] > 0 && layer.get_coordinates()[0] > std::min(x1, x2) && layer.get_coordinates()[0] < std::max(x1, x2)) {
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
                    if (layer.get_dimensions()[1] > 0 && layer.get_coordinates()[1] > y2) {
                        layersInBetween.push_back(layer);
                    }
                }
            }
            else if (y2 < 90 && y1 > 270) {
                for (auto layer : layers) {
                    if (layer.get_dimensions()[1] > 0 && layer.get_coordinates()[1] > y1) {
                        layersInBetween.push_back(layer);
                    }
                }
            }
            else {
                for (auto layer : layers) {
                    if (layer.get_dimensions()[1] > 0 && layer.get_coordinates()[1] > std::min(y1, y2) && layer.get_coordinates()[1] < std::max(y1, y2)) {
                        layersInBetween.push_back(layer);
                    }
                }
            }
        }
        else {
            double y1 = firstTurn.get_coordinates()[1];
            double y2 = secondTurn.get_coordinates()[1];
            for (auto layer : layers) {
                if (layer.get_dimensions()[1] > 0 && layer.get_coordinates()[1] > std::min(y1, y2) && layer.get_coordinates()[1] < std::max(y1, y2)) {
                    layersInBetween.push_back(layer);
                }
            }
        }

    }

    return layersInBetween;
}

StrayCapacitanceOutput StrayCapacitance::calculate_voltages_per_turn(Coil coil, OperatingPoint operatingPoint) {
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

StrayCapacitanceOutput StrayCapacitance::calculate_voltages_per_turn(Coil coil, std::map<std::string, double> voltageRmsPerWinding) {
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
            for (int64_t parallelIndex = 0; parallelIndex < winding.get_number_parallels(); ++parallelIndex) {
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

double get_wire_insulation_relative_permittivity(Wire wire) {
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

std::vector<double> StrayCapacitanceModel::preprocess_data_for_round_wires(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, std::optional<Coil> coil) {
    if (firstWire.get_type() != WireType::ROUND || secondWire.get_type() != WireType::ROUND) {
        throw std::runtime_error("Other wires not implemented yet, check Albach's book");
    }

    InsulationMaterial insulationMaterial = find_insulation_material_by_name(Defaults().defaultInsulationMaterial);

    auto relativePermittivityWireCoatingFirstWire = get_wire_insulation_relative_permittivity(firstWire);
    auto relativePermittivityWireCoatingSecondWire = get_wire_insulation_relative_permittivity(secondWire);
    auto relativePermittivityWireCoating = (relativePermittivityWireCoatingFirstWire + relativePermittivityWireCoatingSecondWire) / 2;

    auto wireCoatingThicknessFirstWire = firstWire.get_coating_thickness();
    auto wireCoatingThicknessSecondWire = secondWire.get_coating_thickness();
    auto wireCoatingThickness = (wireCoatingThicknessFirstWire + wireCoatingThicknessSecondWire) / 2;

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
            
    if (coil) {
        std::vector<Layer> insulationLayersInBetween = StrayCapacitance::get_insulation_layers_between_two_turns(firstTurn, secondTurn, coil.value());

        for (auto layer : insulationLayersInBetween) {
            auto distance = coil->get_insulation_layer_thickness(layer);
            auto relativePermittivity = coil->get_insulation_layer_relative_permittivity(layer);
            distanceThroughLayers += distance;
            distancesThroughLayers.push_back(distance);
            relativePermittivityLayers.push_back(relativePermittivity);
        }
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
        distanceBetweenTurns = DBL_MAX;
        distanceThroughAir = DBL_MAX;
        distanceThroughLayers = DBL_MAX;
        // throw std::invalid_argument("distanceBetweenTurns cannot be negative");
    }
    auto averageTurnLength = (firstTurn.get_length() + secondTurn.get_length()) / 2;

    return {wireCoatingThickness, averageTurnLength, conductingRadius, distanceThroughLayers, distanceThroughAir, relativePermittivityWireCoating, effectiveRelativePermittivityLayers};
}

std::vector<double> StrayCapacitanceParallelPlateModel::preprocess_data_for_planar_wires(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire) {

    InsulationMaterial insulationMaterial = find_insulation_material_by_name(Defaults().defaultInsulationMaterial);

    double conductingDimensionFirstWire = 0;
    double conductingDimensionSecondWire = 0;
    double conductingDimension = 0;

    double distanceBetweenTurns = 0;
    double overlappingDimension = 0;
    if (firstTurn.get_coordinates()[1] == secondTurn.get_coordinates()[1]) {
        distanceBetweenTurns = fabs(firstTurn.get_coordinates()[0] - secondTurn.get_coordinates()[0]);
        distanceBetweenTurns -= firstWire.get_maximum_conducting_width() / 2 + secondWire.get_maximum_conducting_width() / 2;
        conductingDimensionFirstWire = firstWire.get_maximum_conducting_height();
        conductingDimensionSecondWire = secondWire.get_maximum_conducting_height();
        conductingDimension = (conductingDimensionFirstWire + conductingDimensionSecondWire) / 2;
        overlappingDimension = conductingDimension;
    }
    else {
        distanceBetweenTurns = fabs(firstTurn.get_coordinates()[1] - secondTurn.get_coordinates()[1]);
        distanceBetweenTurns -= firstWire.get_maximum_conducting_height() / 2 + secondWire.get_maximum_conducting_height() / 2;
        conductingDimensionFirstWire = firstWire.get_maximum_conducting_width();
        conductingDimensionSecondWire = secondWire.get_maximum_conducting_width();
        double firstWireLeftCoordinate = firstTurn.get_coordinates()[0] - firstTurn.get_dimensions().value()[0] / 2;
        double firstWireRightCoordinate = firstTurn.get_coordinates()[0] + firstTurn.get_dimensions().value()[0] / 2;
        double secondWireLeftCoordinate = secondTurn.get_coordinates()[0] - secondTurn.get_dimensions().value()[0] / 2;
        double secondWireRightCoordinate = secondTurn.get_coordinates()[0] + secondTurn.get_dimensions().value()[0] / 2;
        conductingDimension = std::min(conductingDimensionFirstWire, conductingDimensionSecondWire);
        overlappingDimension = fabs(std::min(firstWireRightCoordinate, secondWireRightCoordinate) - std::max(firstWireLeftCoordinate, secondWireLeftCoordinate));
    }

    distanceBetweenTurns = roundFloat(distanceBetweenTurns, 6);

    double distanceThroughLayers = distanceBetweenTurns;

    std::vector<double> distancesThroughLayers;
    std::vector<double> relativePermittivityLayers;

    auto coatingInsulationMaterial = find_insulation_material_by_name(defaults.defaultPcbInsulationMaterial);
    if (!coatingInsulationMaterial.get_relative_permittivity())
        throw std::runtime_error("FR4 insulation material is missing dielectric constant");
    double effectiveRelativePermittivityLayers = coatingInsulationMaterial.get_relative_permittivity().value();

    auto averageTurnLength = (firstTurn.get_length() + secondTurn.get_length()) / 2;

    return {averageTurnLength, overlappingDimension, distanceThroughLayers, effectiveRelativePermittivityLayers};
}

double StrayCapacitanceMassariniModel::calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    double Dc = conductingRadius * 2;
    double D0;
    double epsilonR;
    if (wireCoatingThickness > 0) {
        D0 = (conductingRadius + wireCoatingThickness) * 2;
        epsilonR = get_effective_relative_permittivity(wireCoatingThickness, relativePermittivityWireCoating, distanceThroughAir + distanceThroughLayers, relativePermittivityInsulationLayers);
    }
    else {
        D0 = (conductingRadius + distanceThroughAir / 2) * 2;
        epsilonR = get_effective_relative_permittivity(distanceThroughAir / 2, vacuumPermittivity, distanceThroughAir + distanceThroughLayers, relativePermittivityInsulationLayers);
    }
    double aux0 = 2 * epsilonR + log(D0 / Dc);
    double aux1 = sqrt(log(D0 / Dc) * (2 * epsilonR + log(D0 / Dc)));
    double aux2 = sqrt(2 * epsilonR * log(D0 / Dc) + pow(log(D0 / Dc), 2));
    double C0 = vacuumPermittivity * averageTurnLength * 2 * epsilonR * atan(((-1 + sqrt(3)) * aux0) / ((1 + sqrt(3)) * aux1)) / aux2;

    return C0;
}

double StrayCapacitanceDuerdothModel::calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    double h = distanceThroughAir + distanceThroughLayers;
    double delta = wireCoatingThickness;
    double r0 = conductingRadius;
    double dtt = 2 * r0 + 2 * wireCoatingThickness;
    double dPrima = 2 * (r0 + delta) + h;
    // double dEff = dPrima - 2.3 * (r0 + delta) + 0.26 * dtt;
    double dEff = dPrima - 0.15 * 2 * (r0 + delta) + 0.26 * (dtt);
    double epsilonEff = get_effective_relative_permittivity(delta, relativePermittivityWireCoating, h, relativePermittivityInsulationLayers);
    if (std::isnan(epsilonEff)) {
        throw std::invalid_argument("epsilonEff is NAN");
    }

    double C0 = vacuumPermittivity * epsilonEff * averageTurnLength * 2 * r0 / dEff;

    return C0;
}

double StrayCapacitanceAlbachModel::calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    double effectiveRelativePermittivity = relativePermittivityInsulationLayers;
    double distanceThroughLayersAndAir = distanceThroughAir + distanceThroughLayers;
    if (distanceThroughAir > 0 && distanceThroughLayers > 0) {
        effectiveRelativePermittivity = get_effective_relative_permittivity(distanceThroughLayers, relativePermittivityInsulationLayers, distanceThroughAir, 1);
    }
    else if (distanceThroughAir > 0 && distanceThroughLayers == 0) {
        effectiveRelativePermittivity = 1;
    }

    double zeta = 1 - wireCoatingThickness / (relativePermittivityWireCoating * (conductingRadius + wireCoatingThickness));
    double beta = 1.0 / zeta * (1 + distanceThroughLayersAndAir / (2 * effectiveRelativePermittivity * (conductingRadius + wireCoatingThickness)));
    double V = beta / sqrt(pow(beta, 2) - 1) * atan(sqrt((beta + 1) / (beta - 1)));
    double Z = 1.0 / (pow(beta, 2) - 1)*((pow(beta, 2) - 2) * V - beta / 2) - std::numbers::pi / 4;
    double Y1 = 1.0 / zeta * (V - std::numbers::pi / 4 + 1.0 / (2 * relativePermittivityWireCoating) * pow(distanceThroughLayers / (conductingRadius + wireCoatingThickness), 2) * Z / zeta);
    double C0 = 2.0 / 3 * vacuumPermittivity * averageTurnLength * Y1;

    if (std::isnan(beta)) {
        throw std::invalid_argument("beta is NAN");
    }

    return C0;
}

double StrayCapacitanceKochModel::calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    double alpha = 1 - wireCoatingThickness / (relativePermittivityWireCoating * conductingRadius);
    double beta;
    if (distanceThroughLayers > 0) {
        beta = 1.0 / alpha * (1 + distanceThroughLayers / (2 * relativePermittivityInsulationLayers * conductingRadius));
    }
    else {
        beta = 1.0 / alpha * (1 + distanceThroughAir / (2 * vacuumPermittivity * conductingRadius));
    }
    double V = beta / sqrt(pow(beta, 2) - 1) * atan(sqrt((beta + 1) / (beta - 1))) - std::numbers::pi / 4;
    double Z = beta * (pow(beta, 2) - 2) / pow(pow(beta, 2) - 1, 1.5) * atan(sqrt((beta + 1) / (beta - 1))) - beta / (2 * (pow(beta, 2) - 1))  - std::numbers::pi / 4;
    double C0 = vacuumPermittivity * averageTurnLength / (1 - wireCoatingThickness / (relativePermittivityWireCoating * conductingRadius)) * (V + 1.0 / (8 * relativePermittivityWireCoating) * pow(2 * wireCoatingThickness / conductingRadius, 2) * Z / (1 - wireCoatingThickness / (relativePermittivityWireCoating * conductingRadius)));

    if (std::isnan(beta)) {
        throw std::invalid_argument("beta is NAN");
    }

    return C0;
}

double StrayCapacitanceParallelPlateModel::calculate_static_capacitance_between_two_turns(double overlappingDimension, double averageTurnLength, double distanceThroughLayers, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;
    return vacuumPermittivity * relativePermittivityInsulationLayers * overlappingDimension * averageTurnLength / distanceThroughLayers;
}

double StrayCapacitance::calculate_static_capacitance_between_two_turns(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, std::optional<Coil> coil) {
    if (firstWire.get_type() == WireType::PLANAR && secondWire.get_type() == WireType::PLANAR) {
        StrayCapacitanceParallelPlateModel model;
        auto aux = model.preprocess_data_for_planar_wires(firstTurn, firstWire, secondTurn, secondWire);
        double averageTurnLength = aux[0];
        double overlappingDimension = aux[1];
        double distanceThroughLayers = aux[2];
        double relativePermittivityInsulationLayers = aux[3];
        double capacitance = model.calculate_static_capacitance_between_two_turns(overlappingDimension, averageTurnLength, distanceThroughLayers, relativePermittivityInsulationLayers);
        return capacitance;
    }
    else {

        auto aux = _model->preprocess_data_for_round_wires(firstTurn, firstWire, secondTurn, secondWire, coil);
        double wireCoatingThickness = aux[0];
        double averageTurnLength = aux[1];
        double conductingRadius = aux[2];
        double distanceThroughLayers = aux[3];
        double distanceThroughAir = aux[4];
        double relativePermittivityWireCoating = aux[5];
        double relativePermittivityInsulationLayers = aux[6];
        double capacitance = _model->calculate_static_capacitance_between_two_turns(wireCoatingThickness, averageTurnLength, conductingRadius, distanceThroughLayers, distanceThroughAir, relativePermittivityWireCoating, relativePermittivityInsulationLayers);

        return capacitance;

    }
}

double StrayCapacitance::calculate_energy_between_two_turns(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, double voltageDrop, std::optional<Coil> coil) {
    double capacitance = calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire, coil);
    double energy = 0.5 * capacitance * pow(voltageDrop, 2);
    return energy;
}

double StrayCapacitance::calculate_energy_density_between_two_turns(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, double voltageDrop, std::optional<Coil> coil) {
    double energy = calculate_energy_between_two_turns(firstTurn, firstWire, secondTurn, secondWire, voltageDrop, coil);
    double area = calculate_area_between_two_turns(firstTurn, secondTurn);
    return energy / area;
}

std::map<std::pair<size_t, size_t>, double> StrayCapacitance::calculate_capacitance_among_turns(Coil coil) {
    if (!coil.get_turns_description()) {
        throw std::invalid_argument("Missing turns description");
    }

    std::map<std::pair<size_t, size_t>, double> capacitanceAmongTurns;

    auto turns = coil.get_turns_description().value();
    auto wirePerWinding = coil.get_wires();

    std::set<std::pair<size_t, size_t>> turnsCombinations;

    for (size_t turnIndex = 0; turnIndex <  turns.size(); ++turnIndex) {
        auto turnWindingIndex = coil.get_winding_index_by_name(turns[turnIndex].get_winding());
        auto turnWire = wirePerWinding[turnWindingIndex];
        auto surroundingTurns = OpenMagnetics::StrayCapacitance::get_surrounding_turns(turns[turnIndex], turns);
        for (auto [surroundingTurn, surroundingTurnIndex] : surroundingTurns) {
            auto key = std::make_pair(turnIndex, surroundingTurnIndex);
            auto inverseKey = std::make_pair(surroundingTurnIndex, turnIndex);
            if (turnsCombinations.contains(key) || turnsCombinations.contains(inverseKey)) {
                continue;
            }
            auto surroundingTurnWindingIndex = coil.get_winding_index_by_name(surroundingTurn.get_winding());
            auto surroundingTurnWire = wirePerWinding[surroundingTurnWindingIndex];
            double capacitance = calculate_static_capacitance_between_two_turns(turns[turnIndex], turnWire, surroundingTurn, surroundingTurnWire, coil);
            capacitanceAmongTurns[key] = capacitance;
            // capacitanceAmongTurns[inverseKey] = capacitance;
            turnsCombinations.insert(key);
        }
    }

    return capacitanceAmongTurns;
}


ScalarMatrixAtFrequency StrayCapacitance::calculate_capacitance_matrix_between_windings(double energy, double voltageDrop, double relativeTurnsRatio) {
    ScalarMatrixAtFrequency scalarMatrixAtFrequency;
    double C0 = energy * 2 / pow(voltageDrop, 2);
    scalarMatrixAtFrequency.set_frequency(0);  // Static result

    auto gamma1 = -C0 / 6;
    auto gamma2 = -C0 / 6;
    auto gamma3 = C0 / 3;
    auto gamma4 = C0 / 3;
    auto gamma5 = C0 / 6;
    auto gamma6 = C0 / 6;

    scalarMatrixAtFrequency.get_mutable_magnitude()["1"]["1"].set_nominal(gamma1 + relativeTurnsRatio * (gamma4 + gamma5));
    scalarMatrixAtFrequency.get_mutable_magnitude()["1"]["2"].set_nominal(-2 * gamma4);
    scalarMatrixAtFrequency.get_mutable_magnitude()["1"]["3"].set_nominal(2 * relativeTurnsRatio * gamma5);
    scalarMatrixAtFrequency.get_mutable_magnitude()["2"]["2"].set_nominal(gamma2 + gamma4 + gamma6);
    scalarMatrixAtFrequency.get_mutable_magnitude()["2"]["3"].set_nominal(2 * gamma6);
    scalarMatrixAtFrequency.get_mutable_magnitude()["3"]["3"].set_nominal(gamma3 + gamma5 + gamma6);

    return scalarMatrixAtFrequency;

}

std::pair<SixCapacitorNetworkPerWinding, TripoleCapacitancePerWinding> StrayCapacitance::calculate_capacitance_models_between_windings(double energy, double voltageDrop, double relativeTurnsRatio) {
    std::map<std::string, double> result;

    double C0 = energy * 2 / pow(voltageDrop, 2);

    auto gamma1 = -C0 / 6;
    auto gamma2 = -C0 / 6;
    auto gamma3 = C0 / 3;
    auto gamma4 = C0 / 3;
    auto gamma5 = C0 / 6;
    auto gamma6 = C0 / 6;

    SixCapacitorNetworkPerWinding sixCapacitorNetworkPerWinding;
    sixCapacitorNetworkPerWinding.set_c1(gamma1);
    sixCapacitorNetworkPerWinding.set_c2(gamma2);
    sixCapacitorNetworkPerWinding.set_c3(gamma3);
    sixCapacitorNetworkPerWinding.set_c4(gamma4);
    sixCapacitorNetworkPerWinding.set_c5(gamma5);
    sixCapacitorNetworkPerWinding.set_c6(gamma6);

    auto C1 = gamma1 + relativeTurnsRatio * gamma2;
    auto C2 = gamma5 + gamma6;
    auto C3 = gamma3;

    TripoleCapacitancePerWinding tripoleCapacitancePerWinding;
    tripoleCapacitancePerWinding.set_c1(C1);
    tripoleCapacitancePerWinding.set_c2(C2);
    tripoleCapacitancePerWinding.set_c3(C3);

    return {sixCapacitorNetworkPerWinding, tripoleCapacitancePerWinding};
}


StrayCapacitanceOutput StrayCapacitance::calculate_capacitance(Coil coil) {
    std::map<std::pair<size_t, size_t>, double> electricEnergyBetweenTurnsMap;
    std::map<std::pair<size_t, size_t>, double> voltageDropBetweenTurnsMap;
    std::map<std::string, std::map<std::string, ScalarMatrixAtFrequency>> capacitanceMatrix;
    std::map<std::string, std::map<std::string, SixCapacitorNetworkPerWinding>> sixCapacitorNetworkPerWinding;
    std::map<std::string, std::map<std::string, TripoleCapacitancePerWinding>> tripoleCapacitancePerWinding;

    auto capacitanceAmongTurns = calculate_capacitance_among_turns(coil);

    std::map<std::string, double> voltageRmsPerWinding;
    double primaryNumberTurns = coil.get_functional_description()[0].get_number_turns();
    for (auto winding : coil.get_functional_description()) {
        double turnsRatio = primaryNumberTurns /  winding.get_number_turns();
        voltageRmsPerWinding[winding.get_name()] = 10.0 / turnsRatio;
    }

    auto strayCapacitanceOutput = StrayCapacitance::calculate_voltages_per_turn(coil, voltageRmsPerWinding);
    auto voltagesPerTurn = strayCapacitanceOutput.get_voltage_per_turn().value();
    auto windings = coil.get_functional_description();
    std::map<std::pair<std::string, std::string>, double> capacitanceMapPerWindings;
    for (auto firstWinding : windings) {
        auto firstWindingName = firstWinding.get_name();
        auto turnsInFirstWinding = coil.get_turns_indexes_by_winding(firstWindingName);
        double minVoltageInFirstWinding = 1;
        double maxVoltageInFirstWinding = 0;
        double minVoltageInSecondWinding = 1;
        double maxVoltageInSecondWinding = 0;
        for (auto secondWinding : windings) {
            auto secondWindingName = secondWinding.get_name();
            auto windingsKey = std::make_pair(firstWindingName, secondWindingName);
            if (capacitanceMapPerWindings.contains(windingsKey) || capacitanceMapPerWindings.contains(std::make_pair(secondWindingName, firstWindingName))) {
                continue;
            }

            double V3 = 42;
            double V3calculated = 0;

            if (firstWindingName == secondWindingName) {
                V3calculated = 0;
                // capacitanceMapPerWindings[windingsKey] = 0;
                // continue;
            }

            double energyInBetweenTheseWindings = 0;
            double voltageDropBetweenWindings = 0;
            double relativeTurnsRatio = 0;
            ScalarMatrixAtFrequency capacitanceMatrixBetweenWindings;
            // std::cout << "firstWindingName: " << firstWindingName << std::endl;
            while (fabs(V3 - V3calculated) / V3 > 0.001) {
                energyInBetweenTheseWindings = 0;
                V3 = V3calculated;
                // double C0 = 0;
                auto turnsInSecondWinding = coil.get_turns_indexes_by_winding(secondWindingName);
                bool windingAreNotAdjacent = true;
                for (auto turnInFirstWinding : turnsInFirstWinding) {
                    auto firstTurnVoltage = voltagesPerTurn[turnInFirstWinding];
                    minVoltageInFirstWinding = std::min(minVoltageInFirstWinding, firstTurnVoltage);
                    maxVoltageInFirstWinding = std::max(maxVoltageInFirstWinding, firstTurnVoltage);
                    for (auto turnInSecondWinding : turnsInSecondWinding) {
                        auto secondTurnVoltage = voltagesPerTurn[turnInSecondWinding];
                        if (firstWindingName != secondWindingName) {
                            secondTurnVoltage = -secondTurnVoltage;
                        }
                        minVoltageInSecondWinding = std::min(minVoltageInSecondWinding, secondTurnVoltage);
                        maxVoltageInSecondWinding = std::max(maxVoltageInSecondWinding, secondTurnVoltage);
                        auto turnsKey = std::make_pair(turnInFirstWinding, turnInSecondWinding);
                    // std::cout << "turnInFirstWinding: " << turnInFirstWinding << std::endl;
                    // std::cout << "turnInSecondWinding: " << turnInSecondWinding << std::endl;
                        if (capacitanceAmongTurns.contains(turnsKey)) {
                            windingAreNotAdjacent = false;
                            double voltageDropAmongTurns = V3 + firstTurnVoltage - secondTurnVoltage;
                            double energyBetweenTurns = 0.5 * capacitanceAmongTurns[turnsKey] * pow(voltageDropAmongTurns, 2);
                            energyInBetweenTheseWindings += energyBetweenTurns;
                            electricEnergyBetweenTurnsMap[turnsKey] = energyBetweenTurns;
                            voltageDropBetweenTurnsMap[turnsKey] = voltageDropAmongTurns;
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

                voltageDropBetweenWindings = maxVoltageInFirstWinding - minVoltageInSecondWinding + V3;
                relativeTurnsRatio = static_cast<double>(firstWinding.get_number_turns()) / secondWinding.get_number_turns(); 
                capacitanceMatrixBetweenWindings = calculate_capacitance_matrix_between_windings(energyInBetweenTheseWindings, voltageDropBetweenWindings, relativeTurnsRatio);
                if (firstWindingName != secondWindingName) {
                    V3calculated = fabs(-(resolve_dimensional_values(capacitanceMatrixBetweenWindings.get_mutable_magnitude()["1"]["3"]) * maxVoltageInFirstWinding + resolve_dimensional_values(capacitanceMatrixBetweenWindings.get_mutable_magnitude()["2"]["3"]) * fabs(minVoltageInSecondWinding)) / resolve_dimensional_values(capacitanceMatrixBetweenWindings.get_mutable_magnitude()["3"]["3"]));
                }
                capacitanceMapPerWindings[windingsKey] = energyInBetweenTheseWindings * 2 / pow(voltageDropBetweenWindings, 2);
            }
            capacitanceMatrix[firstWindingName][secondWindingName] = capacitanceMatrixBetweenWindings;
            capacitanceMatrix[secondWindingName][firstWindingName] = capacitanceMatrixBetweenWindings;
            auto [sixCapacitorNetworkPerWindingBetweenWindings, tripoleCapacitancePerWindingBetweenWindings] = calculate_capacitance_models_between_windings(energyInBetweenTheseWindings, voltageDropBetweenWindings, relativeTurnsRatio);
            sixCapacitorNetworkPerWinding[firstWindingName][secondWindingName] = sixCapacitorNetworkPerWindingBetweenWindings;
            sixCapacitorNetworkPerWinding[secondWindingName][firstWindingName] = sixCapacitorNetworkPerWindingBetweenWindings;
            tripoleCapacitancePerWinding[firstWindingName][secondWindingName] = tripoleCapacitancePerWindingBetweenWindings;
            tripoleCapacitancePerWinding[secondWindingName][firstWindingName] = tripoleCapacitancePerWindingBetweenWindings;
        }
    }

    std::map<std::string, std::map<std::string, double>> staticCapacitanceMapPerWindings;
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        auto winding = windings[windingIndex];
        staticCapacitanceMapPerWindings[winding.get_name()] = std::map<std::string, double>();
    }

    for (size_t firstWindingIndex = 0; firstWindingIndex < windings.size(); ++firstWindingIndex) {
        auto firstWinding = windings[firstWindingIndex];
        for (size_t secondWindingIndex = 0; secondWindingIndex < windings.size(); ++secondWindingIndex) {
            auto secondWinding = windings[secondWindingIndex];
            auto windingsKey = std::make_pair(firstWinding.get_name(), secondWinding.get_name());
            if (!capacitanceMapPerWindings.contains(windingsKey)) {
                continue;
            }
            auto capacitance = capacitanceMapPerWindings[windingsKey];
            staticCapacitanceMapPerWindings[firstWinding.get_name()][secondWinding.get_name()] = capacitance;
            staticCapacitanceMapPerWindings[secondWinding.get_name()][firstWinding.get_name()] = capacitance;
        }
    }


    std::map<std::string, std::map<std::string, double>> electricEnergyAmongTurns;
    std::map<std::string, std::map<std::string, double>> voltageDropAmongTurns;
    std::map<std::string, std::map<std::string, double>> capacitanceAmongTurnsOutput;
    auto turns = coil.get_turns_description().value();
    for (size_t firstTurnIndex = 0; firstTurnIndex < turns.size(); ++firstTurnIndex) {
        auto firstTurnName = turns[firstTurnIndex].get_name();
        for (size_t secondTurnIndex = firstTurnIndex + 1; secondTurnIndex < turns.size(); ++secondTurnIndex) {
            auto secondTurnName = turns[secondTurnIndex].get_name();
            auto turnsKey = std::make_pair(firstTurnIndex, secondTurnIndex);
            electricEnergyAmongTurns[firstTurnName][secondTurnName] = electricEnergyBetweenTurnsMap[turnsKey];
            electricEnergyAmongTurns[secondTurnName][firstTurnName] = electricEnergyBetweenTurnsMap[turnsKey];
            voltageDropAmongTurns[firstTurnName][secondTurnName] = voltageDropBetweenTurnsMap[turnsKey];
            voltageDropAmongTurns[secondTurnName][firstTurnName] = -voltageDropBetweenTurnsMap[turnsKey];
            capacitanceAmongTurnsOutput[firstTurnName][secondTurnName] = capacitanceAmongTurns[turnsKey];
            capacitanceAmongTurnsOutput[secondTurnName][firstTurnName] = capacitanceAmongTurns[turnsKey];
        }
    }

    strayCapacitanceOutput.set_capacitance_among_turns(capacitanceAmongTurnsOutput);
    strayCapacitanceOutput.set_capacitance_among_windings(staticCapacitanceMapPerWindings);
    strayCapacitanceOutput.set_electric_energy_among_turns(electricEnergyAmongTurns);
    strayCapacitanceOutput.set_voltage_drop_among_turns(voltageDropAmongTurns);
    strayCapacitanceOutput.set_maxwell_capacitance_matrix(calculate_maxwell_capacitance_matrix(coil, staticCapacitanceMapPerWindings));
    strayCapacitanceOutput.set_capacitance_matrix(capacitanceMatrix);
    strayCapacitanceOutput.set_six_capacitor_network_per_winding(sixCapacitorNetworkPerWinding);
    strayCapacitanceOutput.set_tripole_capacitance_per_winding(tripoleCapacitancePerWinding);

    return strayCapacitanceOutput;
}


std::vector<ScalarMatrixAtFrequency> StrayCapacitance::calculate_maxwell_capacitance_matrix(Coil coil, std::map<std::string, std::map<std::string, double>> capacitanceAmongWindings) {
    ScalarMatrixAtFrequency scalarMatrixAtFrequency;
    scalarMatrixAtFrequency.set_frequency(0);  // Static result

    auto windings = coil.get_functional_description();

    for (auto firstWinding : windings) {
        std::string firstWindingName = firstWinding.get_name();
        double capacitanceSum = 0;
        for (auto secondWinding : windings) {
            std::string secondWindingName = secondWinding.get_name();
            if (capacitanceAmongWindings.contains(firstWindingName)) {
                if (capacitanceAmongWindings[firstWindingName].contains(secondWindingName)) {
                    auto capacitance = capacitanceAmongWindings[firstWindingName][secondWindingName];
                    capacitanceSum += capacitance;
                    if (firstWindingName != secondWindingName) {
                        scalarMatrixAtFrequency.get_mutable_magnitude()[firstWindingName][secondWindingName].set_nominal(-capacitance);
                        scalarMatrixAtFrequency.get_mutable_magnitude()[secondWindingName][firstWindingName].set_nominal(-capacitance);
                    }
                }
            }
            else {
                capacitanceSum += capacitanceAmongWindings[secondWindingName][firstWindingName];
                throw std::runtime_error("Old code, this should not happen");
                if (firstWindingName != secondWindingName) {
                    auto capacitance = capacitanceAmongWindings[secondWindingName][firstWindingName];
                    scalarMatrixAtFrequency.get_mutable_magnitude()[firstWindingName][secondWindingName].set_nominal(-capacitance);
                    scalarMatrixAtFrequency.get_mutable_magnitude()[secondWindingName][firstWindingName].set_nominal(-capacitance);
                }
            }
        }
        scalarMatrixAtFrequency.get_mutable_magnitude()[firstWindingName][firstWindingName].set_nominal(capacitanceSum);
    }

    return {scalarMatrixAtFrequency};
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

double StrayCapacitanceOneLayer::calculate_capacitance(Coil coil) {
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
