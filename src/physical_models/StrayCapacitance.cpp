#include "physical_models/StrayCapacitance.h"
#include "physical_models/Impedance.h"
#include <algorithm>
#include "support/Settings.h"
#include "Defaults.h"
#include "Constants.h"
#include "Models.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#include <cfloat>
#include <cmath>
#include <limits>
#include <numbers>
#include <set>
#include "support/Exceptions.h"
#include <magic_enum.hpp>

namespace OpenMagnetics {

// Helper function to get all coordinate positions for a turn (main + additional for toroidal)
static std::vector<std::vector<double>> get_all_turn_coordinates(const Turn& turn) {
    std::vector<std::vector<double>> coords;
    coords.push_back(turn.get_coordinates());
    if (turn.get_additional_coordinates()) {
        auto additionalCoords = turn.get_additional_coordinates().value();
        for (const auto& addCoord : additionalCoords) {
            if (addCoord.size() >= 2) {
                coords.push_back(addCoord);
            }
        }
    }
    return coords;
}

// Helper function to compute the global minimum surface-to-surface gap between any two turns
static double compute_global_minimum_gap(const std::vector<Turn>& turnsDescription) {
    double globalMinGap = DBL_MAX;
    
    for (size_t i = 0; i < turnsDescription.size(); ++i) {
        auto coords1 = get_all_turn_coordinates(turnsDescription[i]);
        double maxDim1 = std::max(turnsDescription[i].get_dimensions().value()[0], 
                                   turnsDescription[i].get_dimensions().value()[1]);
        
        for (size_t j = i + 1; j < turnsDescription.size(); ++j) {
            auto coords2 = get_all_turn_coordinates(turnsDescription[j]);
            double maxDim2 = std::max(turnsDescription[j].get_dimensions().value()[0], 
                                       turnsDescription[j].get_dimensions().value()[1]);
            
            // Find minimum gap between any coordinate pair
            for (const auto& c1 : coords1) {
                for (const auto& c2 : coords2) {
                    double centerDist = hypot(c1[0] - c2[0], c1[1] - c2[1]);
                    double surfaceGap = centerDist - maxDim1 / 2 - maxDim2 / 2;
                    // Only consider positive gaps (non-overlapping turns)
                    if (surfaceGap >= 0 && surfaceGap < globalMinGap) {
                        globalMinGap = surfaceGap;
                    }
                }
            }
        }
    }
    
    // If no valid gap found (all turns overlap), return a small default
    // This ensures at least one turn pair will be found
    if (globalMinGap == DBL_MAX) {
        return 1e-6;
    }
    
    return globalMinGap;
}

// The core "coating" that sits in the winding-to-core dielectric path is the insulating JACKET
// over the ferrite surface the turns rest on (parylene / epoxy / nylon / glass). A
// MAGNETIC_EPOXY coating is a different thing: the powder-loaded shield cap moulded over the
// winding of a semishielded drum. It is not between a turn and the core surface, and its
// material is a powder core material, not an insulation (looking it up as one threw
// MISSING_DATA "Insulation material not found: Kool Mµ 26" and made every semishielded drum
// uncomputable under the full model). So it contributes no jacket layer here: thickness 0,
// relative permittivity 1; the wire enamel alone bounds the element, as on a bare core.
static std::pair<double, double> resolve_core_jacket(const Core& core) {
    auto coating = core.get_functional_description().get_coating();
    if (coating && std::holds_alternative<CoreCoating>(coating.value())) {
        auto type = std::get<CoreCoating>(coating.value()).get_type();
        if (type && type.value() == CoatingType::MAGNETIC_EPOXY) {
            return {0.0, 1.0};
        }
    }
    double thickness = core.get_coating_thickness();
    // A bare (uncoated) ferrite has no coating layer; the wire enamel alone then bounds the
    // element. get_coating_relative_permittivity() throws when there is no coating, so only
    // resolve it when a coating actually exists.
    double relativePermittivity = thickness > 0 ? core.get_coating_relative_permittivity() : 1.0;
    return {thickness, relativePermittivity};
}

std::vector<std::pair<Turn, size_t>> StrayCapacitance::get_surrounding_turns(Turn currentTurn, std::vector<Turn> turnsDescription, double globalMinimumGap) {
    std::vector<std::pair<Turn, size_t>> surroundingTurns;
    auto factor = Defaults().overlappingFactorSurroundingTurns;

    auto dx1 = currentTurn.get_dimensions().value()[0];
    auto dy1 = currentTurn.get_dimensions().value()[1];

    // Get all coordinate positions for the current turn (including additional_coordinates for toroidal)
    auto currentCoords = get_all_turn_coordinates(currentTurn);

    for (size_t turnIndex = 0; turnIndex < turnsDescription.size(); ++turnIndex) {
        auto potentiallySurroundingTurn = turnsDescription[turnIndex];

        // Skip if this is the same turn as currentTurn
        if (potentiallySurroundingTurn.get_coordinates()[0] == currentTurn.get_coordinates()[0] &&
            potentiallySurroundingTurn.get_coordinates()[1] == currentTurn.get_coordinates()[1]) {
            continue;
        }

        auto dx2 = potentiallySurroundingTurn.get_dimensions().value()[0];
        auto dy2 = potentiallySurroundingTurn.get_dimensions().value()[1];

        double minimumDimension = std::min(std::max(dx1, dy1), std::max(dx2, dy2));
        double maximumDim1 = std::max(dx1, dy1);
        double maximumDim2 = std::max(dx2, dy2);

        // Get all coordinate positions for the potentially surrounding turn
        auto potentialCoords = get_all_turn_coordinates(potentiallySurroundingTurn);

        // Find minimum distance across all coordinate combinations for threshold check
        double minDistance = DBL_MAX;
        double bestX1 = 0, bestY1 = 0, bestX2 = 0, bestY2 = 0;
        bool foundValidPair = false;

        for (const auto& c1 : currentCoords) {
            for (const auto& c2 : potentialCoords) {
                double x1 = c1[0], y1 = c1[1];
                double x2 = c2[0], y2 = c2[1];

                // Skip if same position
                if (x1 == x2 && y1 == y2) {
                    continue;
                }

                double distance = hypot(x2 - x1, y2 - y1) - maximumDim1 / 2 - maximumDim2 / 2;
                if (distance < minDistance) {
                    minDistance = distance;
                    bestX1 = x1; bestY1 = y1;
                    bestX2 = x2; bestY2 = y2;
                    foundValidPair = true;
                }
            }
        }

        // Skip if no valid coordinate pair found
        if (!foundValidPair) {
            continue;
        }

        // Dual threshold logic:
        // threshold1: 1x minimumDimension (wire diameter) - for tightly wound coils
        // threshold2: 1.5x globalMinimumGap - adaptive based on actual geometry
        // Use the maximum of the two thresholds to ensure we capture adjacent turns
        // while not including too many distant turns
        double threshold1 = minimumDimension;  // 1x wire diameter
        double threshold2 = (globalMinimumGap > 0) ? globalMinimumGap * 1.5 : 0;
        double distanceThreshold = std::max(threshold1, threshold2);
        
        if (minDistance > distanceThreshold) {
            continue;
        }

        double maximumDimensionOf12 = (maximumDim1 + maximumDim2) / 2;

        // Use the closest coordinate pair for the "turn between" check
        double x1 = bestX1, y1 = bestY1;
        double x2 = bestX2, y2 = bestY2;

        bool thereIsTurnBetween12 = false;
        for (auto potentiallyCollidingTurn : turnsDescription) {
            auto dx0 = potentiallyCollidingTurn.get_dimensions().value()[0];
            auto dy0 = potentiallyCollidingTurn.get_dimensions().value()[1];

            // Check all coordinates of the potentially colliding turn
            auto collidingCoords = get_all_turn_coordinates(potentiallyCollidingTurn);
            for (const auto& c0 : collidingCoords) {
                double x0 = c0[0], y0 = c0[1];

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

                double maximumDimensionOf0 = std::max(dx0, dy0);
                auto distanceFrom0toLine12 = fabs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt(pow(y2 - y1, 2) + pow(x2 - x1, 2));
                if (maximumDimensionOf12 / 2 + maximumDimensionOf0 / 2 * factor > distanceFrom0toLine12) {
                    thereIsTurnBetween12 = true;
                    break;
                }
            }
            if (thereIsTurnBetween12) break;
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
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing excitation for windingIndex: " + std::to_string(windingIndex));
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
            throw NaNResultException("voltage cannot be nan");
        }
        if (std::isinf(voltage)) {
            throw InvalidInputException(ErrorCode::CALCULATION_INVALID_RESULT, "voltage cannot be inf");
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
    // Handle edge case where distances are DBL_MAX (turns overlapping)
    // In this case, return the larger permittivity as a safe fallback
    if (std::isinf(firstThickness) || std::isinf(secondThickness)) {
        return std::max(firstRelativePermittivity, secondRelativePermittivity);
    }
    
    // Handle edge case where one thickness is zero
    if (firstThickness == 0) {
        return secondRelativePermittivity;
    }
    if (secondThickness == 0) {
        return firstRelativePermittivity;
    }
    
    return firstRelativePermittivity * secondRelativePermittivity * (firstThickness + secondThickness) / (firstThickness * secondRelativePermittivity + secondThickness * firstRelativePermittivity);
}

/**
 * @brief Returns the effective relative permittivity of a wire's insulation coating.
 *
 * For ROUND wires: Uses empirical formula from Albach Chapter 3.
 * For LITZ wires: Uses serving material permittivity (fallback to strand enamel formula).
 * For FOIL/RECTANGULAR: Uses coating material permittivity (fallback to typical film value).
 * For PLANAR: Returns typical FR4 permittivity.
 *
 * @param wire The wire to get the insulation permittivity for
 * @return Effective relative permittivity of the wire insulation
 */
double get_wire_insulation_relative_permittivity(Wire wire) {
    if (wire.get_type() == WireType::ROUND) {
        // Empirical formula for enamelled round wire coating permittivity
        // Based on Albach, "Induktivitaeten in der Leistungselektronik", Chapter 3
        auto conductingRadius = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
        return 2.5 + 0.7 / sqrt(2 * conductingRadius * 1000);
    }
    else if (wire.get_type() == WireType::LITZ) {
        // For litz wire, the relevant insulation for turn-to-turn capacitance is the serving
        // (outer insulation of the bundle). Use the serving material's permittivity if available,
        // otherwise fall back to the strand enamel empirical formula.
        // Reference: Albach "Induktivitaeten in der Leistungselektronik", Chapter 3
        // Reference: Biela/Kolar "Using Transformer Parasitics for Resonant Converters"
        try {
            return wire.get_coating_relative_permittivity();
        }
        catch (...) {
            // Fallback: use strand enamel empirical formula with strand radius
            auto strand = wire.resolve_strand();
            auto strandConductingRadius = resolve_dimensional_values(strand.get_conducting_diameter()) / 2;
            return 2.5 + 0.7 / sqrt(2 * strandConductingRadius * 1000);
        }
    }
    else if (wire.get_type() == WireType::FOIL || wire.get_type() == WireType::RECTANGULAR) {
        // For foil and rectangular wires, use the coating material's permittivity.
        // These wires normally use parallel plate model, but permittivity may still be needed
        // for mixed-type turn pairs.
        try {
            return wire.get_coating_relative_permittivity();
        }
        catch (const InvalidInputException&) {
            // Only missing coating/permittivity DATA falls back to the typical
            // film value; any other failure must propagate.
            return 3.5; // Default: typical polyester/polyimide film
        }
    }
    else {
        // PLANAR or unknown: FR4 or similar substrate material
        return 4.5;
    }
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
        throw ModelNotAvailableException("Unknown Stray capacitance model, available options are: {KOCH, ALBACH, DUERDOTH, MASSARINI}");
}

std::vector<double> StrayCapacitanceModel::preprocess_data_for_round_wires(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, std::optional<Coil> coil) {
    // Accept ROUND and LITZ wire types. For LITZ, we treat the bundle as an equivalent
    // round conductor with the outer bundle diameter as the "conducting" surface and
    // the serving/insulation as the "coating".
    // Reference: Albach "Induktivitaeten in der Leistungselektronik", Chapter 3
    // Reference: Biela/Kolar "Using Transformer Parasitics for Resonant Converters"
    
    auto isRoundLike = [](WireType type) {
        return type == WireType::ROUND || type == WireType::LITZ;
    };
    
    if (!isRoundLike(firstWire.get_type()) || !isRoundLike(secondWire.get_type())) {
        throw NotImplementedException("preprocess_data_for_round_wires only supports ROUND and LITZ wires, got: " + std::string(magic_enum::enum_name(firstWire.get_type())) + " and " + std::string(magic_enum::enum_name(secondWire.get_type())));
    }

    InsulationMaterial insulationMaterial = find_insulation_material_by_name(Defaults().defaultInsulationMaterial);

    auto relativePermittivityWireCoatingFirstWire = get_wire_insulation_relative_permittivity(firstWire);
    auto relativePermittivityWireCoatingSecondWire = get_wire_insulation_relative_permittivity(secondWire);
    auto relativePermittivityWireCoating = (relativePermittivityWireCoatingFirstWire + relativePermittivityWireCoatingSecondWire) / 2;

    // For LITZ: The "conducting diameter" for capacitance purposes is the bare bundle diameter
    //           (without serving). The bare bundle surface acts as the effective conductor
    //           boundary for inter-turn capacitance. We must compute this correctly using
    //           get_outer_diameter_bare_litz() rather than relying on get_coating_thickness()
    //           which includes strand packing space and gives incorrect results.
    // For ROUND: The conducting diameter is simply the bare wire diameter.
    double conductingDiameterFirstWire;
    double conductingDiameterSecondWire;
    double outerDiameterFirstWire;
    double outerDiameterSecondWire;
    double wireCoatingThicknessFirstWire;
    double wireCoatingThicknessSecondWire;
    
    if (firstWire.get_type() == WireType::LITZ) {
        // Overall outer diameter (includes serving/insulation)
        outerDiameterFirstWire = firstWire.calculate_outer_diameter();
        
        // Compute bare bundle diameter directly from strand geometry
        // This is the effective "conductor surface" for capacitance
        auto strand = firstWire.resolve_strand();
        auto strandCoating = Wire::resolve_coating(strand);
        auto standard = WireStandard::IEC_60317;
        if (firstWire.get_standard()) {
            standard = firstWire.get_standard().value();
        }
        int grade = 1;
        if (strandCoating && strandCoating->get_grade()) {
            grade = strandCoating->get_grade().value();
        }
        int64_t numConductors = firstWire.get_number_conductors().value();
        double strandConductingDiameter = resolve_dimensional_values(strand.get_conducting_diameter());
        
        double bareBundleDiameter = Wire::get_outer_diameter_bare_litz(
            strandConductingDiameter, numConductors, grade, standard);
        
        // The "conducting diameter" for capacitance is the bare bundle surface
        conductingDiameterFirstWire = bareBundleDiameter;
        
        // The "coating thickness" is only the serving/outer insulation
        wireCoatingThicknessFirstWire = (outerDiameterFirstWire - bareBundleDiameter) / 2;
        if (wireCoatingThicknessFirstWire < 0) {
            wireCoatingThicknessFirstWire = 0;
        }
    } else {
        // ROUND wire
        conductingDiameterFirstWire = firstWire.get_maximum_conducting_width();
        outerDiameterFirstWire = firstWire.get_maximum_outer_width();
        wireCoatingThicknessFirstWire = firstWire.get_coating_thickness();
    }
    
    if (secondWire.get_type() == WireType::LITZ) {
        outerDiameterSecondWire = secondWire.calculate_outer_diameter();
        
        auto strand = secondWire.resolve_strand();
        auto strandCoating = Wire::resolve_coating(strand);
        auto standard = WireStandard::IEC_60317;
        if (secondWire.get_standard()) {
            standard = secondWire.get_standard().value();
        }
        int grade = 1;
        if (strandCoating && strandCoating->get_grade()) {
            grade = strandCoating->get_grade().value();
        }
        int64_t numConductors = secondWire.get_number_conductors().value();
        double strandConductingDiameter = resolve_dimensional_values(strand.get_conducting_diameter());
        
        double bareBundleDiameter = Wire::get_outer_diameter_bare_litz(
            strandConductingDiameter, numConductors, grade, standard);
        
        conductingDiameterSecondWire = bareBundleDiameter;
        
        wireCoatingThicknessSecondWire = (outerDiameterSecondWire - bareBundleDiameter) / 2;
        if (wireCoatingThicknessSecondWire < 0) {
            wireCoatingThicknessSecondWire = 0;
        }
    } else {
        // ROUND wire
        conductingDiameterSecondWire = secondWire.get_maximum_conducting_width();
        outerDiameterSecondWire = secondWire.get_maximum_outer_width();
        wireCoatingThicknessSecondWire = secondWire.get_coating_thickness();
    }
    
    // Now average the coating thickness (after per-wire computation)
    auto wireCoatingThickness = (wireCoatingThicknessFirstWire + wireCoatingThicknessSecondWire) / 2;
    // ABT #851: this is the average conducting RADIUS. It used to be written as
    // (diameter1 + diameter2) / 2 — the average DIAMETER — and every round-wire pair model
    // downstream consumed that as the radius, i.e. 2x the real value for equal wires.
    auto conductingRadius = (conductingDiameterFirstWire + conductingDiameterSecondWire) / 4;

    // For toroidal cores, check distance between all coordinate combinations (inner/outer halves)
    auto coords1 = get_all_turn_coordinates(firstTurn);
    auto coords2 = get_all_turn_coordinates(secondTurn);
    
    // ABT #851: the separation that sets a turn pair's capacitance is their CLOSEST
    // APPROACH. A toroidal turn carries inner and outer crossing coordinates, so the
    // combinations of two ADJACENT turns mix the true near-zero gap (inner-inner) with
    // across-the-core distances of 10-20 mm (inner-outer). This used to AVERAGE them,
    // telling every pair model that touching turns sit millimetres apart: the plate-family
    // statics (Albach/Koch/Duerdoth) collapsed to ~0.07 pF per pair where the near-contact
    // value is ~3-4 pF, and only Massarini's logarithmic form survived — which is why the
    // choice of static model appeared to matter enormously on toroids. Adjacency detection
    // (get_surrounding_turns) already uses the minimum over the same combinations; the
    // capacitance must see the same geometry it was selected on. Off toroids each turn has
    // a single coordinate, so min == mean and nothing changes there.
    double minDistance = std::numeric_limits<double>::max();
    int validDistancesCount = 0;
    for (const auto& c1 : coords1) {
        for (const auto& c2 : coords2) {
            // Skip if same position (same turn)
            if (c1[0] == c2[0] && c1[1] == c2[1]) {
                continue;
            }
            double dist = hypot(c1[0] - c2[0], c1[1] - c2[1]);
            dist -= outerDiameterFirstWire / 2 + outerDiameterSecondWire / 2;
            minDistance = std::min(minDistance, dist);
            validDistancesCount++;
        }
    }
    
    double distanceBetweenTurns;
    if (validDistancesCount > 0) {
        distanceBetweenTurns = minDistance;
    } else {
        // Fallback to single coordinate calculation
        double x1 = firstTurn.get_coordinates()[0];
        double y1 = firstTurn.get_coordinates()[1];
        double x2 = secondTurn.get_coordinates()[0];
        double y2 = secondTurn.get_coordinates()[1];
        distanceBetweenTurns = hypot(x2 - x1, y2 - y1) - outerDiameterFirstWire / 2 - outerDiameterSecondWire / 2;
    }
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
    // The insulation layers are collected along the winding path, not the turns'
    // closest approach, so their summed thickness can exceed the actual
    // surface-to-surface separation (seen on the boundary turns of a
    // separated-winding toroidal CMC: 4.7 mm of layers across a 2.8 mm gap).
    // The dielectric between two turns cannot be thicker than their gap — cap it,
    // or distanceThroughAir goes negative and corrupts every downstream formula
    // (ABT #173: negative static capacitance).
    if (distanceThroughLayers > std::max(distanceBetweenTurns, 0.0)) {
        distanceThroughLayers = std::max(distanceBetweenTurns, 0.0);
    }
    double distanceThroughAir = distanceBetweenTurns - distanceThroughLayers;

    if (distanceBetweenTurns < 0) {
        // For toroidal cores, negative distance can occur when turns overlap in 2D projection
        // but are actually at different radial positions. Use a small positive distance instead.
        distanceBetweenTurns = 1e-9;  // 1nm minimum gap
        distanceThroughAir = distanceBetweenTurns - distanceThroughLayers;
        if (distanceThroughAir < 0) {
            distanceThroughAir = 0;
        }
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
        throw InvalidInputException(ErrorCode::INVALID_INSULATION_DATA, "FR4 insulation material is missing dielectric constant");
    double effectiveRelativePermittivityLayers = coatingInsulationMaterial.get_relative_permittivity().value();

    auto averageTurnLength = (firstTurn.get_length() + secondTurn.get_length()) / 2;

    return {averageTurnLength, overlappingDimension, distanceThroughLayers, effectiveRelativePermittivityLayers};
}

/**
 * @brief Calculates static capacitance between two turns using the Massarini model.
 * 
 * Based on: "Self-Capacitance of Inductors" by Antonio Massarini and Marian K. Kazimierczuk
 * IEEE Transactions on Power Electronics, Vol. 12, No. 4, July 1997
 * https://ieeexplore.ieee.org/document/602562
 * 
 * The model derives turn-to-turn capacitance for round conductors with insulation coating.
 * 
 * Key equations implemented:
 * - Equation (9): Effective permittivity calculation for multi-layer dielectrics
 *   ε_eff = (d1 + d2) / (d1/ε1 + d2/ε2)
 * 
 * - Equation (14): Capacitance per unit length between adjacent conductors
 *   The formula accounts for the geometry of round wires using the conductor diameter (Dc)
 *   and outer diameter (D0) including insulation.
 * 
 * - Equation (15): Simplified form using arctangent function
 *   C0 = ε0 * lt * 2εr * arctan(((-1+√3)(2εr + ln(D0/Dc))) / ((1+√3)√(ln(D0/Dc)(2εr + ln(D0/Dc))))) 
 *        / √(2εr*ln(D0/Dc) + ln²(D0/Dc))
 *   where:
 *   - lt = average turn length
 *   - Dc = conductor diameter (conducting radius × 2)
 *   - D0 = outer diameter including insulation
 *   - εr = effective relative permittivity
 * 
 * - Equation (20): Self-capacitance of single-layer coil from sum of turn-to-turn capacitances
 *   Cs = (1/n) * Σ Ctt
 * 
 * @param wireCoatingThickness Thickness of wire insulation coating [m]
 * @param averageTurnLength Average length of a single turn [m]
 * @param conductingRadius Radius of the conducting wire core [m]
 * @param distanceThroughLayers Distance through insulation layers between turns [m]
 * @param distanceThroughAir Distance through air gap between turns [m]
 * @param relativePermittivityWireCoating Relative permittivity of wire insulation
 * @param relativePermittivityInsulationLayers Relative permittivity of layer insulation
 * @return Static capacitance between turns [F]
 */
double StrayCapacitanceMassariniModel::calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    // Dc: Conductor diameter (bare wire)
    double Dc = conductingRadius * 2;
    // D0: Outer diameter including insulation
    double D0;
    // εr: Effective relative permittivity
    double epsilonR;
    if (wireCoatingThickness > 0) {
        // Wire has insulation coating - use coating as outer boundary
        D0 = (conductingRadius + wireCoatingThickness) * 2;
        // Calculate effective permittivity using Eq. (9) for series dielectric layers
        epsilonR = get_effective_relative_permittivity(wireCoatingThickness, relativePermittivityWireCoating, distanceThroughAir + distanceThroughLayers, relativePermittivityInsulationLayers);
    }
    else {
        // Bare wire - use half the air gap as effective outer boundary
        D0 = (conductingRadius + distanceThroughAir / 2) * 2;
        epsilonR = get_effective_relative_permittivity(distanceThroughAir / 2, 1.0, distanceThroughAir + distanceThroughLayers, relativePermittivityInsulationLayers);
    }
    
    // Implementation of Eq. (15) from Massarini paper
    // Intermediate terms for the arctangent-based capacitance formula
    double aux0 = 2 * epsilonR + log(D0 / Dc);
    double aux1 = sqrt(log(D0 / Dc) * (2 * epsilonR + log(D0 / Dc)));
    double aux2 = sqrt(2 * epsilonR * log(D0 / Dc) + pow(log(D0 / Dc), 2));
    
    // C0 = ε0 * lt * 2εr * arctan(...) / √(...)
    double C0 = vacuumPermittivity * averageTurnLength * 2 * epsilonR * atan(((-1 + sqrt(3)) * aux0) / ((1 + sqrt(3)) * aux1)) / aux2;

    return C0;
}

/**
 * @brief Calculates static capacitance between two turns using the Duerdoth model.
 * 
 * Based on: "Equivalent Capacitances of Transformer Windings" by W. T. Duerdoth
 * Electronic Engineering, Vol. 18, 1946
 * 
 * The Duerdoth model is a simplified parallel-plate approximation for round wires
 * with empirical correction factors for the effective plate separation distance.
 * 
 * This model was one of the earliest analytical approaches and uses geometric
 * simplifications to approximate the capacitance between adjacent round conductors.
 * 
 * Key parameters:
 * - h: Total gap distance (air + insulation layers)
 * - δ: Wire insulation coating thickness
 * - r0: Conductor radius
 * - dtt: Total wire outer diameter = 2(r0 + δ)
 * - d': Center-to-center distance = 2(r0 + δ) + h
 * 
 * - deff: Effective parallel-plate distance (empirically corrected)
 *   deff = d' - 0.15 * 2(r0 + δ) + 0.26 * dtt
 *   The coefficients 0.15 and 0.26 are empirical corrections for round wire geometry
 * 
 * - εeff: Effective permittivity for series dielectric layers
 * 
 * - C0: Parallel-plate capacitance with effective distance
 *   C0 = ε0 * εeff * lt * 2r0 / deff
 * 
 * @note This model is less accurate than Koch/Massarini but useful for quick estimates
 * 
 * @param wireCoatingThickness Thickness of wire insulation coating [m]
 * @param averageTurnLength Average length of a single turn [m]
 * @param conductingRadius Radius of the conducting wire core [m]
 * @param distanceThroughLayers Distance through insulation layers between turns [m]
 * @param distanceThroughAir Distance through air gap between turns [m]
 * @param relativePermittivityWireCoating Relative permittivity of wire insulation
 * @param relativePermittivityInsulationLayers Relative permittivity of layer insulation
 * @return Static capacitance between turns [F]
 */
double StrayCapacitanceDuerdothModel::calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    // h: Total gap distance between wire surfaces
    double h = distanceThroughAir + distanceThroughLayers;
    // δ: Insulation coating thickness
    double delta = wireCoatingThickness;
    // r0: Bare conductor radius
    double r0 = conductingRadius;
    // dtt: Total outer diameter of insulated wire
    double dtt = 2 * r0 + 2 * wireCoatingThickness;
    // d': Center-to-center distance between conductors
    double dPrima = 2 * (r0 + delta) + h;
    
    // deff: Effective parallel-plate distance with empirical corrections
    // Coefficients 0.15 and 0.26 correct for round wire geometry
    // double dEff = dPrima - 2.3 * (r0 + delta) + 0.26 * dtt;  // Alternative formulation
    double dEff = dPrima - 0.15 * 2 * (r0 + delta) + 0.26 * (dtt);
    
    // εeff: Effective relative permittivity for series dielectric layers
    double epsilonEff = get_effective_relative_permittivity(delta, relativePermittivityWireCoating, h, relativePermittivityInsulationLayers);
    if (std::isnan(epsilonEff)) {
        throw std::invalid_argument("epsilonEff is NAN");
    }

    // C0: Parallel-plate capacitance with effective distance
    // C = ε0 * εeff * A / d, where A = lt * 2r0 (projected area)
    double C0 = vacuumPermittivity * epsilonEff * averageTurnLength * 2 * r0 / dEff;

    return C0;
}

/**
 * @brief Calculates static capacitance between two turns using the Albach model.
 * 
 * Based on: "Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre 
 * parasitären Eigenschaften" by Manfred Albach
 * Springer Vieweg, 2017, ISBN 978-3-658-15081-5
 * Chapter 3: "Die Kapazität von Luftspulen", pages 45-68
 * 
 * This model calculates turn-to-turn capacitance for round conductors with insulation
 * by computing the electric field energy stored between wires.
 * 
 * Key equations from Chapter 3:
 * 
 * - Eq. (3.9): Parameter definitions
 *   θ = 1/(δ/(εD*r0))  where δ = coating thickness, εD = coating permittivity
 *   β = (1/θ) * (1 + h/(2*εF*(r0+δ)))  where h = gap distance, εF = layer permittivity
 * 
 * - Eq. (3.16): Y1 geometric factor for energy calculation
 *   Y1 = V/4 + (1/(2*εD)) * (δ/r0)² * Z
 *   where:
 *   V = (β/√(β²-1)) * arctan(√((β+1)/(β-1)))
 *   Z = ((β²-2)*V - β/2) / (β²-1) - π/4
 * 
 * - Final capacitance from Eq. (3.14), (3.19):
 *   Wges = ε0 * lw * U² * Y1  (energy in rectangular region between turns)
 *   C0 = (2/3) * ε0 * lt * Y1  (capacitance with 2/3 factor for voltage distribution)
 * 
 * The 2/3 factor comes from the quadratic voltage distribution integration
 * in Eq. (3.21): ∫(l/lmax * U0)² dl = lmax/3 * U0²
 * 
 * @note The Albach model uses the OUTER conductor radius (r0 + δ) in the ζ and β
 *       parameters, unlike the Koch model which uses bare conductor radius r0.
 *       This provides improved accuracy for thick insulation coatings.
 * 
 * @param wireCoatingThickness Thickness of wire insulation coating δ [m]
 * @param averageTurnLength Average length of a single turn lt [m]
 * @param conductingRadius Radius of the conducting wire core r0 [m]
 * @param distanceThroughLayers Distance through insulation layers between turns [m]
 * @param distanceThroughAir Distance through air gap between turns [m]
 * @param relativePermittivityWireCoating Relative permittivity of wire insulation εD
 * @param relativePermittivityInsulationLayers Relative permittivity of layer insulation εF
 * @return Static capacitance between turns [F]
 */
double StrayCapacitanceAlbachModel::calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    // Calculate effective relative permittivity for combined air and insulation layers
    double effectiveRelativePermittivity = relativePermittivityInsulationLayers;
    double distanceThroughLayersAndAir = distanceThroughAir + distanceThroughLayers;
    if (distanceThroughAir > 0 && distanceThroughLayers > 0) {
        effectiveRelativePermittivity = get_effective_relative_permittivity(distanceThroughLayers, relativePermittivityInsulationLayers, distanceThroughAir, 1);
    }
    else if (distanceThroughAir > 0 && distanceThroughLayers == 0) {
        effectiveRelativePermittivity = 1;
    }

    // Handle edge case: when turns overlap (distance is DBL_MAX), return a large capacitance
    if (std::isinf(distanceThroughLayersAndAir)) {
        // When turns overlap or are extremely close, capacitance approaches infinity
        // Return a very large value instead of trying to calculate
        return 1e-6;  // 1 µF - effectively infinite for practical purposes
    }

    // ζ: Modified insulation parameter using OUTER radius (r0 + δ)
    // This differs from Koch model which uses bare conductor radius
    double zeta = 1 - wireCoatingThickness / (relativePermittivityWireCoating * (conductingRadius + wireCoatingThickness));
    
    // β: Gap geometry parameter
    double beta = 1.0 / zeta * (1 + distanceThroughLayersAndAir / (2 * effectiveRelativePermittivity * (conductingRadius + wireCoatingThickness)));
    
    // Handle edge case where beta is too close to 1 - use parallel plate approximation
    // When beta <= 1, sqrt(beta^2 - 1) returns NaN
    if (beta <= 1.001) {
        // Use parallel plate approximation as fallback
        double totalDistance = wireCoatingThickness + distanceThroughLayersAndAir;
        double epsilonEff = get_effective_relative_permittivity(wireCoatingThickness, relativePermittivityWireCoating, distanceThroughLayersAndAir, relativePermittivityInsulationLayers);
        double projectedWidth = conductingRadius * 2;  // Projected width of round wire
        return vacuumPermittivity * epsilonEff * projectedWidth * averageTurnLength / totalDistance;
    }
    
    // V: Arctangent auxiliary function
    double V = beta / sqrt(pow(beta, 2) - 1) * atan(sqrt((beta + 1) / (beta - 1)));
    
    // Z: Second-order correction term
    double Z = 1.0 / (pow(beta, 2) - 1)*((pow(beta, 2) - 2) * V - beta / 2) - std::numbers::pi / 4;
    
    // Y1: Combined auxiliary function with insulation thickness correction.
    // The second-order term is a perturbation in the WIRE COATING thickness δ
    // (cf. Koch eq. 7: (1/(8εr))·(2δ/r0)² — identical at leading order to
    // (1/(2εr))·(δ/(r0+δ))² used here). It previously plugged in
    // distanceThroughLayers (the inter-turn insulation stack): harmless for
    // adjacent same-layer turns (0), but with millimetre layer stacks the
    // "correction" dominated (Z < 0) and drove the capacitance NEGATIVE
    // (ABT #173: -0.57 pF on the boundary turns of a separated-winding CMC).
    double Y1 = 1.0 / zeta * (V - std::numbers::pi / 4 + 1.0 / (2 * relativePermittivityWireCoating) * pow(wireCoatingThickness / (conductingRadius + wireCoatingThickness), 2) * Z / zeta);
    
    // C0: Final capacitance with 2/3 geometric factor from Albach
    double C0 = 2.0 / 3 * vacuumPermittivity * averageTurnLength * Y1;

    if (std::isnan(beta)) {
        throw std::invalid_argument("beta is NAN");
    }

    return C0;
}

/**
 * @brief Calculates static capacitance between two turns using the Koch model.
 * 
 * Based on: "Berechnung der Kapazität von Spulen, insbesondere in Schalenkernen" by K. Koch
 * Reproduced in: "Using Transformer Parasitics for Resonant Converters—A Review of the 
 * Calculation of the Stray Capacitance of Transformers" by Juergen Biela and Johann W. Kolar
 * IEEE Transactions on Industry Applications, Vol. 44, No. 1, January/February 2008
 * https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/biela_IEEETrans_ReviewStrayCap.pdf
 * 
 * The Koch model provides capacitance calculation for round conductors with insulation,
 * using an elliptic integral approach approximated with arctangent functions.
 * 
 * Key equations from Biela/Kolar review (Section III-A):
 * 
 * - Equation (3): Parameter α accounting for insulation coating
 *   α = 1 - δ/(εr_ins * r0)
 *   where δ = insulation thickness, r0 = conductor radius
 * 
 * - Equation (4): Parameter β for air/insulation gap geometry
 *   β = (1/α) * (1 + h/(2 * εr_layer * r0))
 *   where h = gap distance between wires
 * 
 * - Equation (5): Auxiliary function V (arctangent form of elliptic integral)
 *   V = (β/√(β²-1)) * arctan(√((β+1)/(β-1))) - π/4
 * 
 * - Equation (6): Auxiliary function Z (second-order correction)
 *   Z = β(β²-2)/((β²-1)^(3/2)) * arctan(√((β+1)/(β-1))) - β/(2(β²-1)) - π/4
 * 
 * - Equation (7): Final capacitance per unit length
 *   C = ε0 * lt / α * (V + (1/(8*εr_ins)) * (2δ/r0)² * Z/α)
 * 
 * @note This model is particularly accurate for tightly wound coils with thin insulation
 *       and is commonly used in shell-type core calculations.
 * 
 * @param wireCoatingThickness Thickness of wire insulation coating [m]
 * @param averageTurnLength Average length of a single turn [m]
 * @param conductingRadius Radius of the conducting wire core [m]
 * @param distanceThroughLayers Distance through insulation layers between turns [m]
 * @param distanceThroughAir Distance through air gap between turns [m]
 * @param relativePermittivityWireCoating Relative permittivity of wire insulation
 * @param relativePermittivityInsulationLayers Relative permittivity of layer insulation
 * @return Static capacitance between turns [F]
 */
double StrayCapacitanceKochModel::calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) {
    auto vacuumPermittivity = Constants().vacuumPermittivity;

    // Handle edge case: when turns overlap (distance is DBL_MAX), return a large capacitance
    double totalDistance = distanceThroughLayers + distanceThroughAir;
    if (std::isinf(totalDistance)) {
        // When turns overlap or are extremely close, capacitance approaches infinity
        // Return a very large value instead of trying to calculate
        return 1e-6;  // 1 µF - effectively infinite for practical purposes
    }

    // α: Parameter accounting for insulation coating effect - Eq. (3)
    // α = 1 - δ/(εr * r0), where δ = coating thickness
    double alpha = 1 - wireCoatingThickness / (relativePermittivityWireCoating * conductingRadius);
    
    // β: Parameter for gap geometry - Eq. (4)
    // β = (1/α) * (1 + h/(2 * εr * r0))
    double beta;
    if (distanceThroughLayers > 0) {
        // Gap contains insulation layers
        beta = 1.0 / alpha * (1 + distanceThroughLayers / (2 * relativePermittivityInsulationLayers * conductingRadius));
    }
    else {
        // Gap is pure air (relative permittivity 1)
        beta = 1.0 / alpha * (1 + distanceThroughAir / (2 * 1.0 * conductingRadius));
    }
    
    // Handle edge case where beta is too close to 1 - use parallel plate approximation
    // When beta <= 1, sqrt(beta^2 - 1) returns NaN
    if (beta <= 1.001) {
        // Use parallel plate approximation as fallback
        double totalDistance = wireCoatingThickness + distanceThroughLayers + distanceThroughAir;
        double epsilonEff = get_effective_relative_permittivity(wireCoatingThickness, relativePermittivityWireCoating, distanceThroughLayers + distanceThroughAir, relativePermittivityInsulationLayers);
        double projectedWidth = conductingRadius * 2;  // Projected width of round wire
        return vacuumPermittivity * epsilonEff * projectedWidth * averageTurnLength / totalDistance;
    }
    
    // V: Auxiliary function (arctangent approximation of elliptic integral) - Eq. (5)
    // V = (β/√(β²-1)) * arctan(√((β+1)/(β-1))) - π/4
    double V = beta / sqrt(pow(beta, 2) - 1) * atan(sqrt((beta + 1) / (beta - 1))) - std::numbers::pi / 4;
    
    // Z: Second-order correction term - Eq. (6)
    // Z = β(β²-2)/((β²-1)^(3/2)) * arctan(...) - β/(2(β²-1)) - π/4
    double Z = beta * (pow(beta, 2) - 2) / pow(pow(beta, 2) - 1, 1.5) * atan(sqrt((beta + 1) / (beta - 1))) - beta / (2 * (pow(beta, 2) - 1))  - std::numbers::pi / 4;
    
    // Final capacitance - Eq. (7)
    // C = ε0 * lt / α * (V + (1/(8*εr)) * (2δ/r0)² * Z/α)
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

// ABT #395: a non-finite turn-to-turn capacitance used to escape this function silently and only
// surface much later as "Energy cannot be nan", naming neither the turns nor the cause. The
// arithmetic that converted it: the energy of a pair is 0.5*C*dV^2, and dV is EXACTLY zero between
// corresponding turns of two parallels (they are the same electrical node, so the per-turn voltage
// divider gives them identical potentials) -- so an infinite C became 0.5*inf*0 = NaN. That is why
// the failure separated perfectly on numberParallels: with a single parallel no turn pair is ever
// at identical potential, so the infinity stayed an infinity and never turned into a NaN.
//
// The infinity itself is real, not a rounding artefact: two BARE conductors whose surfaces are in
// contact have no dielectric and no separation between them, and every model here diverges. That
// is inconsistent input rather than something to approximate around, so report it where the
// geometry is still in hand instead of fabricating a minimum gap.
static void throw_if_capacitance_is_not_finite(double capacitance, Turn firstTurn, Wire firstWire,
                                               Turn secondTurn, Wire secondWire) {
    if (std::isfinite(capacitance)) {
        return;
    }
    auto describeTurn = [](Turn turn) {
        std::string description = turn.get_name();
        auto coordinates = turn.get_coordinates();
        return description + " at (" + std::to_string(coordinates[0]) + ", " + std::to_string(coordinates[1]) + ")";
    };
    double firstCoatingThickness = firstWire.get_coating_thickness();
    double secondCoatingThickness = secondWire.get_coating_thickness();
    double centerToCenterDistance = hypot(firstTurn.get_coordinates()[0] - secondTurn.get_coordinates()[0],
                                          firstTurn.get_coordinates()[1] - secondTurn.get_coordinates()[1]);
    double surfaceSeparation = centerToCenterDistance - firstWire.get_maximum_outer_width() / 2 -
                               secondWire.get_maximum_outer_width() / 2;
    std::string reason;
    if (firstCoatingThickness <= 0 || secondCoatingThickness <= 0) {
        reason = " The wires report zero coating thickness (" + std::to_string(firstCoatingThickness) + " m and " +
                 std::to_string(secondCoatingThickness) +
                 " m), so there is no dielectric between the conductors; give the wire its insulation"
                 " (a coating thickness, or an outer diameter larger than the conducting diameter).";
    }
    throw NaNResultException("Turn-to-turn capacitance is " + std::to_string(capacitance) + " between " +
                             describeTurn(firstTurn) + " and " + describeTurn(secondTurn) +
                             ", whose conducting surfaces are " + std::to_string(surfaceSeparation) +
                             " m apart." + reason);
}

// ABT #406: report the ELECTRICAL fault, from the geometry, before any capacitance model runs.
//
// The condition is that the two CONDUCTORS touch. Coating thickness is what normally keeps them
// apart: in a close-wound coil the OUTER surfaces are in contact by design, and the copper is
// still separated by t1 + t2, which is perfectly valid. Take the coating away and the outer
// surface IS the copper, so the same close-wound geometry puts bare metal against bare metal —
// shorted turns.
//
// Testing the geometry rather than the resulting number matters for one case in particular: with
// a single parallel the divergent capacitance never becomes a NaN (it takes two turns at
// identical potential to form 0.5*inf*0), so it used to travel all the way into the exported
// netlist without raising anything — a silently wrong answer rather than a crash. A geometric
// test fires there too, and does not depend on which capacitance formula happens to diverge.
//
// Deliberately conservative: it also requires a missing coating. Turns whose coated conductors
// merely overlap slightly through layout approximation are left to the existing non-finite guard,
// so this cannot start rejecting coils that model fine today.
static void throw_if_turns_are_shorted(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire) {
    // Wire::get_coating_thickness() now answers this correctly for litz (it reads the
    // strand enamel, not the bundle serving), so the local special-case that used to live
    // here is gone — one definition of "what insulates this conductor", not two.
    double firstCoatingThickness = firstWire.get_coating_thickness();
    double secondCoatingThickness = secondWire.get_coating_thickness();
    if (firstCoatingThickness > 0 && secondCoatingThickness > 0) {
        return;
    }

    auto isRoundLike = [](WireType type) {
        return type == WireType::ROUND || type == WireType::LITZ;
    };
    double deltaX = firstTurn.get_coordinates()[0] - secondTurn.get_coordinates()[0];
    double deltaY = firstTurn.get_coordinates()[1] - secondTurn.get_coordinates()[1];

    // Gap between the CONDUCTING surfaces. Round-like pairs get the true centre-to-centre test;
    // anything involving a flat wire gets the bounding-box test, whose axes are the ones its
    // width and height are defined along (a circle test would be wrong for a strip).
    double conductingSurfaceGap;
    if (isRoundLike(firstWire.get_type()) && isRoundLike(secondWire.get_type())) {
        conductingSurfaceGap = hypot(deltaX, deltaY) - firstWire.get_maximum_conducting_width() / 2 -
                               secondWire.get_maximum_conducting_width() / 2;
    }
    else {
        double gapAlongWidth = fabs(deltaX) - firstWire.get_maximum_conducting_width() / 2 -
                               secondWire.get_maximum_conducting_width() / 2;
        double gapAlongHeight = fabs(deltaY) - firstWire.get_maximum_conducting_height() / 2 -
                                secondWire.get_maximum_conducting_height() / 2;
        // Separated along EITHER axis is enough to keep them apart.
        conductingSurfaceGap = std::max(gapAlongWidth, gapAlongHeight);
    }

    // Contact, not "gap > 0": turns laid exactly one conductor width apart come out of the
    // winder with a gap of a few 1e-19 m (0.010997 - 0.010797 - 0.0002 in doubles), and a
    // strict > 0 let that pair through to the models, which then diverged and surfaced as an
    // anonymous non-finite capacitance instead of this message (seen on a 3-turn bare foil
    // shield, ABT #853). One nanometre is far below any real dielectric and far above the
    // rounding of millimetre-scale coordinates, so anything closer than that is touching.
    constexpr double contactTolerance = 1e-9;
    if (conductingSurfaceGap > contactTolerance) {
        return;
    }

    auto describeTurn = [](Turn turn) {
        auto coordinates = turn.get_coordinates();
        return turn.get_name() + " at (" + std::to_string(coordinates[0]) + ", " +
               std::to_string(coordinates[1]) + ")";
    };
    throw ShortedTurnsException(
        "Turns " + describeTurn(firstTurn) + " and " + describeTurn(secondTurn) +
        " are in electrical contact: their conducting surfaces are " +
        std::to_string(conductingSurfaceGap) + " m apart and the wires report coating thicknesses of " +
        std::to_string(firstCoatingThickness) + " m and " + std::to_string(secondCoatingThickness) +
        " m, so there is no insulation between the conductors. This describes shorted turns, not a"
        " capacitance to compute; give the wire its insulation (a coating thickness, or an outer"
        " diameter larger than the conducting diameter).");
}

double StrayCapacitance::calculate_static_capacitance_between_two_turns(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, std::optional<Coil> coil) {
    throw_if_turns_are_shorted(firstTurn, firstWire, secondTurn, secondWire);

    auto isFlatWire = [](WireType type) {
        return type == WireType::PLANAR || type == WireType::FOIL || type == WireType::RECTANGULAR;
    };
    auto isRoundLike = [](WireType type) {
        return type == WireType::ROUND || type == WireType::LITZ;
    };

    if (isFlatWire(firstWire.get_type()) && isFlatWire(secondWire.get_type())) {
        // Both wires are flat: use parallel plate model
        StrayCapacitanceParallelPlateModel model;
        auto aux = model.preprocess_data_for_planar_wires(firstTurn, firstWire, secondTurn, secondWire);
        double averageTurnLength = aux[0];
        double overlappingDimension = aux[1];
        double distanceThroughLayers = aux[2];
        double relativePermittivityInsulationLayers = aux[3];
        double capacitance = model.calculate_static_capacitance_between_two_turns(overlappingDimension, averageTurnLength, distanceThroughLayers, relativePermittivityInsulationLayers);
        throw_if_capacitance_is_not_finite(capacitance, firstTurn, firstWire, secondTurn, secondWire);
        return capacitance;
    }
    else if (isRoundLike(firstWire.get_type()) && isRoundLike(secondWire.get_type())) {
        // Both wires are round-like (ROUND or LITZ): use cylindrical wire model
        auto aux = _model->preprocess_data_for_round_wires(firstTurn, firstWire, secondTurn, secondWire, coil);
        double wireCoatingThickness = aux[0];
        double averageTurnLength = aux[1];
        double conductingRadius = aux[2];
        double distanceThroughLayers = aux[3];
        double distanceThroughAir = aux[4];
        double relativePermittivityWireCoating = aux[5];
        double relativePermittivityInsulationLayers = aux[6];
        double capacitance = _model->calculate_static_capacitance_between_two_turns(wireCoatingThickness, averageTurnLength, conductingRadius, distanceThroughLayers, distanceThroughAir, relativePermittivityWireCoating, relativePermittivityInsulationLayers);
        throw_if_capacitance_is_not_finite(capacitance, firstTurn, firstWire, secondTurn, secondWire);
        return capacitance;
    }
    else {
        // Mixed types (e.g., ROUND next to FOIL, LITZ next to RECTANGULAR)
        // Use parallel plate model as a conservative approximation
        StrayCapacitanceParallelPlateModel model;
        auto aux = model.preprocess_data_for_planar_wires(firstTurn, firstWire, secondTurn, secondWire);
        double averageTurnLength = aux[0];
        double overlappingDimension = aux[1];
        double distanceThroughLayers = aux[2];
        double relativePermittivityInsulationLayers = aux[3];
        double capacitance = model.calculate_static_capacitance_between_two_turns(overlappingDimension, averageTurnLength, distanceThroughLayers, relativePermittivityInsulationLayers);
        throw_if_capacitance_is_not_finite(capacitance, firstTurn, firstWire, secondTurn, secondWire);
        return capacitance;
    }
}


// Effective bare-conductor radius for the turn-to-core field-spreading integral.
// Round/litz wires pass their true radius. Flat conductors (planar / foil /
// rectangular) are modelled as an AREA-EQUIVALENT cylinder, r_eq = sqrt(w*h/pi):
// the turn-to-core element is a fringing/curvature field that leaves the conductor
// surface and spreads to the core, and calculate_turn_to_core_capacitance stays
// finite at contact only because of that curvature. A true flat-face parallel
// plate would diverge at zero separation (touching plates), and these planar parts
// carry no modelled trace-to-core dielectric to bound it — so reusing the curvature
// integral with an equivalent radius is the physically consistent, finite choice
// (rather than fabricating a dielectric thickness absent from the data).
static double turn_to_core_equivalent_radius(Wire wire) {
    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        return wire.get_maximum_conducting_width() / 2.0;
    }
    double conductingWidth = wire.get_maximum_conducting_width();
    double conductingHeight = wire.get_maximum_conducting_height();
    return std::sqrt(conductingWidth * conductingHeight / std::numbers::pi);
}

double StrayCapacitance::calculate_turn_to_core_capacitance(double conductingRadius, double turnLength,
                                                            double wireCoatingThickness, double wireCoatingRelativePermittivity,
                                                            double airGapToCore,
                                                            double coreCoatingThickness, double coreCoatingRelativePermittivity) {
    // Capacitance of a turn to the equipotential ferrite core, after Kovacic et al.,
    // "Analytical Wideband Model of a Common-Mode Choke" (IEEE TPEL 2012), eqs (36)-(38).
    // The displacement field leaves the bare conductor surface and reaches the ferrite
    // along field lines that spread with the angle phi; each line crosses, in series, the
    // dielectric stack
    //   wire enamel (t_we / eps_we) | air gap | core coating (t_co / eps_co),
    // where the air-gap length grows with phi as (phi/sin phi)*(w_s + r*(1-cos phi)).
    // Integrating dC = eps0 * dS / sum(x_i/eps_ri), with dS = L*(d_w/2) dphi, over phi in
    // [-pi/2, pi/2] (symmetric, so 0..pi/2 doubled). The ferrite's own permittivity does
    // not enter — it is the electrode. Unlike the old Massarini reuse (bounded by the wire
    // enamel and almost insensitive to the core gap/coating), this is properly sensitive to
    // the wire-to-core air gap and the core coating thickness — the terms that set the real,
    // sub-pF turn-to-core value. The wire enamel is kept in the stack so the element stays
    // finite even at zero air gap and no core coating (close-wound, uncoated).
    //
    // ABT #848 (2026-08-23): the Kovacic construction above — radial field lines leaving the
    // lower half of the conductor, each stretched by phi/sin(phi) — is an engineering
    // approximation, and against the exact solution it is low by 1.3-2x for the geometries the
    // WE chokes actually have (checked on 230 measured parts: 0.5-0.8x). The exact electrostatic
    // solution for a conducting cylinder of radius r whose axis sits a height h above a
    // conducting plane (image method, Smythe §4.13) is
    //     C' = 2 pi eps0 / acosh(h / r)   per unit length,
    // and it counts the whole circumference, not a half. The series dielectric layers between
    // the copper and the ferrite — wire enamel, the wire-to-core air gap and the core coating or
    // case — enter as their air-equivalent thickness t_i / eps_ri (exact in the parallel-plate
    // limit t << r; for a thick case wall the result lies between the all-air and
    // all-dielectric solutions, which is the right order). Same inputs, same bounded behaviour:
    // finite at zero air gap (the enamel keeps h > r), monotone in the gap and the coating.
    if (conductingRadius <= 0 || turnLength <= 0) {
        return 0;
    }
    const double vacuumPermittivity = Constants().vacuumPermittivity;
    const double wireRadius = conductingRadius;  // contract: conductingRadius is the radius (Dc = 2*r)
    const double enamelTerm = wireCoatingRelativePermittivity > 0 ? wireCoatingThickness / wireCoatingRelativePermittivity : 0.0;
    const double coatingTerm = coreCoatingRelativePermittivity > 0 ? coreCoatingThickness / coreCoatingRelativePermittivity : 0.0;
    const double airEquivalentGap = enamelTerm + std::max(0.0, airGapToCore) + coatingTerm;
    if (airEquivalentGap <= 0) {
        // Bare copper on a bare conductor: the image solution diverges, and so does the physics
        // — that is a short, which the callers' geometry guards report; never a number here.
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
            "Turn-to-core capacitance asked for a bare conductor in contact with the core (no wire"
            " coating, no air gap, no core coating): that is a short circuit, not a capacitance");
    }
    const double heightOverRadius = 1.0 + airEquivalentGap / wireRadius;
    return 2.0 * std::numbers::pi * vacuumPermittivity * turnLength / std::acosh(heightOverRadius);
}

double StrayCapacitance::core_image_factor(const Core& core, double frequency) {
    // ABT #848: the floating-core network (turn -> core -> turn) treats the core as an
    // EQUIPOTENTIAL ELECTRODE — a perfect image plane for every turn. That is exact for a
    // conductor and, for AC fields, holds for any body whose complex permittivity dwarfs the
    // dielectric between it and the turns: the field cannot penetrate it, and its surface
    // polarisation charge follows the field like free charge would. A body whose permittivity
    // is only a few times the surroundings images a line charge with the dielectric half-space
    // fraction
    //     beta = (|eps2| - eps1) / (|eps2| + eps1),
    // with eps2 the core's COMPLEX relative permittivity at the frequency of interest,
    //     eps2 = eps' - j eps'',   eps'' = eps''_dielectric + 1 / (omega eps0 rho),
    // and eps1 the dielectric on the core surface (jacket, else air). MnZn at 1-100 MHz:
    // |eps2| ~ 1e4-1e5 -> beta = 1.000; nanocrystalline ribbon (rho ~ 1e-6 Ohm.m): |eps2| ~ 1e10
    // -> beta = 1; NiZn K07 at 10 MHz: eps' 15, conduction negligible -> beta ~ 0.9 bare,
    // ~0.6 under an epoxy jacket. Every quantity is a MAS material property; nothing chosen.
    //
    // Measured (WE NiZn chokes, geometry right): the four K07 parts resonate at x1.20 of
    // RedExpert with beta = 1, x1.08 with this fraction, and x3.3 with the network switched
    // OFF — so a dielectric core IS an image plane to its turns (the displacement-based
    // picture above), and gating the network on conduction alone, which had seemed the safer
    // physics, is contradicted by the data. The fraction is what the measurements support.
    //
    // The external medium eps1 is the dielectric the core surface touches: its jacket when it
    // has one (parylene/epoxy/nylon case), else air.
    auto material = core.resolve_material();
    double epsExternal = 1.0;
    {
        auto coating = core.get_functional_description().get_coating();
        bool magneticShield = coating && std::holds_alternative<CoreCoating>(coating.value()) &&
                              std::get<CoreCoating>(coating.value()).get_type() &&
                              std::get<CoreCoating>(coating.value()).get_type().value() == CoatingType::MAGNETIC_EPOXY;
        if (!magneticShield && core.get_coating_thickness() > 0) {
            epsExternal = core.get_coating_relative_permittivity();
        }
    }
    double omega = 2.0 * std::numbers::pi * frequency;
    const double vacuumPermittivity = Constants().vacuumPermittivity;
    double epsReal = 0.0, epsImag = 0.0;
    bool haveData = false;
    if (material.get_permittivity() && material.get_permittivity()->get_complex()) {
        auto complexPermittivity = material.get_permittivity()->get_complex().value();
        epsReal = interpolate_permittivity_points(complexPermittivity.get_real(), frequency);
        epsImag = interpolate_permittivity_points(complexPermittivity.get_imaginary(), frequency);
        haveData = epsReal > 0;
    }
    else if (!material.get_resistivity().empty() && frequency > 0) {
        // Conduction only (metallic ribbon cores, or a ferrite whose permittivity the database
        // does not carry but whose resistivity it does): eps'' = 1 / (omega eps0 rho).
        auto resistivity = material.get_resistivity()[0];
        double rho = resistivity.get_value();
        if (rho > 0) {
            epsImag = 1.0 / (omega * vacuumPermittivity * rho);
            haveData = true;
        }
    }
    if (!haveData) {
        // No permittivity and no resistivity in the database: the model keeps its pre-#848
        // meaning — the core is the electrode. (A ferrite with no data is not a reason to
        // invent a permittivity; it is a reason to add it to MAS.)
        return 1.0;
    }
    double epsCoreMagnitude = std::hypot(epsReal, epsImag);
    double beta = (epsCoreMagnitude - epsExternal) / (epsCoreMagnitude + epsExternal);
    return std::clamp(beta, 0.0, 1.0);
}

// ABT #848: the air gap between a turn's conductor surface and the core surface it faces, from
// the coil's real geometry instead of "close-wound, 0". Toroid (round window): the first layer
// sits on the (jacketed) core, a deeper layer sits on the layer below it, so the gap is how far
// the turn's radial position lies inside the bore surface. Bobbin-wound (rectangular window):
// the turn is separated from the central ferrite column by the bobbin wall and any inner
// layers, i.e. by its distance from the ferrite column surface; the bobbin plastic in that gap
// is counted as air (air-equivalent thickness would be t/eps_r — this is the conservative side,
// and documented rather than invented). A bare conductor at zero gap is then only ever asked
// for when it really is in contact, and calculate_turn_to_core_capacitance refuses it as the
// short it is.
static double turn_to_core_air_gap(Coil& coil, const Turn& turn, double conductingRadius) {
    auto bobbin = coil.resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        return 0.0;
    }
    auto processed = bobbin.get_processed_description().value();
    auto windows = processed.get_winding_windows();
    if (windows.empty()) {
        return 0.0;
    }
    auto coordinates = turn.get_coordinates();
    if (bobbin.get_winding_window_shape() == WindingWindowShape::ROUND) {
        // Bore radius of the surface the first layer rests on; the turn's radius from the axis.
        if (!windows[0].get_radial_height()) {
            return 0.0;
        }
        double boreRadius = windows[0].get_radial_height().value();
        double turnRadius = std::hypot(coordinates[0], coordinates.size() > 1 ? coordinates[1] : 0.0);
        return std::max(0.0, (boreRadius - conductingRadius) - turnRadius);
    }
    // Rectangular window: the ferrite column surface sits at column_width - column_thickness
    // from the axis (column_width includes the bobbin wall), the turn's conductor surface at
    // |x| - r.
    if (!processed.get_column_width()) {
        return 0.0;
    }
    double ferriteColumnSurface = processed.get_column_width().value() - processed.get_column_thickness();
    return std::max(0.0, std::abs(coordinates[0]) - conductingRadius - ferriteColumnSurface);
}

double StrayCapacitance::calculate_winding_to_core_capacitance(Coil coil, Core core, std::string windingName, std::optional<double> frequency) {
    // Total capacitance from one whole winding to the (equipotential) ferrite core,
    // = the parallel sum of every turn's turn-to-core element (all turns of the
    // winding share the single core node). This is the building block for the
    // through-core inter-winding coupling that sets the differential-mode resonance
    // of a separated-winding CMC (the two windings have no adjacent turns, so the
    // only capacitive path between them is turn -> core -> turn).
    //
    // The core coating is the dielectric floor that keeps the turn-to-core element
    // finite at contact; the winding is taken as close-wound onto the coated core
    // (air gap = 0), matching the literature shortcut C_winding-to-core ~ 2*C_turn-to-turn
    // for bobbin-less toroids (Tu et al., CPSS 2025; Massarini & Kazimierczuk 1997).
    if (!coil.get_turns_description()) {
        coil.wind();
    }

    auto [coreCoatingThickness, coreCoatingRelativePermittivity] = resolve_core_jacket(core);
    double imageFactor = frequency ? core_image_factor(core, frequency.value()) : 1.0;

    auto turns = coil.get_turns_description().value();
    auto wirePerWinding = coil.get_wires();
    auto windingIndex = coil.get_winding_index_by_name(windingName);
    auto wire = wirePerWinding[windingIndex];

    // calculate_turn_to_core_capacitance takes the bare-conductor RADIUS. Round/litz
    // wires use their true radius; flat conductors use an area-equivalent radius
    // (see turn_to_core_equivalent_radius), so planar/foil/rectangular windings get a
    // finite through-core element instead of throwing.
    double conductingRadius = turn_to_core_equivalent_radius(wire);
    double wireCoatingThickness = wire.get_coating_thickness();
    double wireCoatingRelativePermittivity = get_wire_insulation_relative_permittivity(wire);

    double windingToCore = 0;
    for (auto& turn : turns) {
        if (turn.get_winding() != windingName) {
            continue;
        }
        windingToCore += imageFactor * calculate_turn_to_core_capacitance(
            conductingRadius, turn.get_length(),
            wireCoatingThickness, wireCoatingRelativePermittivity,
            turn_to_core_air_gap(coil, turn, conductingRadius),
            coreCoatingThickness, coreCoatingRelativePermittivity);
    }
    return windingToCore;
}

double StrayCapacitance::calculate_through_core_capacitance(Coil coil, Core core,
        const std::string& firstWindingName, const std::string& secondWindingName,
        const std::vector<double>& voltagesPerTurn, std::optional<double> frequency) {
    // Inter-winding capacitance between two separated windings through the floating,
    // equipotential ferrite core. Energy method (CPSS 2025 core-potential approach):
    //   1. each turn i has a turn-to-core element C_i and sits at potential V_i;
    //   2. the two windings carry opposing DM currents, so the second winding's turns
    //      enter with the opposite sign (V_i -> -V_i);
    //   3. the core is a single floating node — its potential balances the charge:
    //      V_core = sum(C_i V_i) / sum(C_i);
    //   4. the energy stored in all the turn-to-core capacitors is
    //      W = sum 0.5 C_i (V_i - V_core)^2, and the effective inter-winding
    //      capacitance referenced to the differential terminal voltage V_dm is
    //      C = 2W / V_dm^2.
    // Weighting by the actual per-turn potential (instead of summing every turn at the
    // full winding voltage) is what removes the gross overestimate of the naive sum.
    if (!coil.get_turns_description()) {
        coil.wind();
    }
    auto turns = coil.get_turns_description().value();
    auto wirePerWinding = coil.get_wires();

    auto [coreCoatingThickness, coreCoatingRelativePermittivity] = resolve_core_jacket(core);
    double imageFactor = frequency ? core_image_factor(core, frequency.value()) : 1.0;

    std::vector<double> turnCoreCapacitance;
    std::vector<double> turnPotential;
    double sumCV = 0;
    double sumC = 0;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto windingName = turns[turnIndex].get_winding();
        double sign;
        if (windingName == firstWindingName) {
            sign = 1.0;
        }
        else if (windingName == secondWindingName) {
            sign = -1.0;  // opposing DM current
        }
        else {
            continue;
        }

        auto wire = wirePerWinding[coil.get_winding_index_by_name(windingName)];
        // Bare-conductor radius for the image solution; flat conductors
        // use an area-equivalent radius (see turn_to_core_equivalent_radius).
        double capacitance = imageFactor * calculate_turn_to_core_capacitance(
            turn_to_core_equivalent_radius(wire), turns[turnIndex].get_length(),
            wire.get_coating_thickness(), get_wire_insulation_relative_permittivity(wire),
            turn_to_core_air_gap(coil, turns[turnIndex], turn_to_core_equivalent_radius(wire)),
            coreCoatingThickness, coreCoatingRelativePermittivity);
        double potential = sign * voltagesPerTurn[turnIndex];

        turnCoreCapacitance.push_back(capacitance);
        turnPotential.push_back(potential);
        sumCV += capacitance * potential;
        sumC += capacitance;
    }

    if (sumC <= 0) {
        return 0;
    }
    double corePotential = sumCV / sumC;

    double energy = 0;
    double maxFirstPotential = -std::numeric_limits<double>::max();
    double minSecondPotential = std::numeric_limits<double>::max();
    for (size_t k = 0; k < turnCoreCapacitance.size(); ++k) {
        energy += 0.5 * turnCoreCapacitance[k] * std::pow(turnPotential[k] - corePotential, 2);
        if (turnPotential[k] >= 0) {
            maxFirstPotential = std::max(maxFirstPotential, turnPotential[k]);
        }
        else {
            minSecondPotential = std::min(minSecondPotential, turnPotential[k]);
        }
    }
    // Differential terminal voltage across the two windings (first winding spans
    // [0, maxFirst]; second, negated, spans [minSecond, 0]).
    double differentialVoltage = maxFirstPotential - minSecondPotential;
    if (differentialVoltage <= 0) {
        return 0;
    }
    return 2.0 * energy / (differentialVoltage * differentialVoltage);
}

double StrayCapacitance::calculate_winding_to_core_self_energy(Coil coil, Core core,
        const std::string& windingName,
        const std::vector<double>& voltagesPerTurn, std::optional<double> frequency) {
    // One winding against the floating core (ABT #848). Mirrors
    // calculate_through_core_capacitance exactly — same per-turn elements, same
    // charge-balanced core node — but over a single winding's turns, all with
    // positive sign, and returning the ENERGY rather than a terminal-referenced
    // capacitance: the caller folds it into the self-pair energy sum, which is
    // already referenced to the winding's own voltage span.
    if (!coil.get_turns_description()) {
        coil.wind();
    }
    auto turns = coil.get_turns_description().value();
    auto wirePerWinding = coil.get_wires();

    auto [coreCoatingThickness, coreCoatingRelativePermittivity] = resolve_core_jacket(core);
    double imageFactor = frequency ? core_image_factor(core, frequency.value()) : 1.0;

    std::vector<double> turnCoreCapacitance;
    std::vector<double> turnPotential;
    double sumCV = 0;
    double sumC = 0;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        if (turns[turnIndex].get_winding() != windingName) {
            continue;
        }
        auto wire = wirePerWinding[coil.get_winding_index_by_name(windingName)];
        double capacitance = imageFactor * calculate_turn_to_core_capacitance(
            turn_to_core_equivalent_radius(wire), turns[turnIndex].get_length(),
            wire.get_coating_thickness(), get_wire_insulation_relative_permittivity(wire),
            turn_to_core_air_gap(coil, turns[turnIndex], turn_to_core_equivalent_radius(wire)),
            coreCoatingThickness, coreCoatingRelativePermittivity);
        turnCoreCapacitance.push_back(capacitance);
        turnPotential.push_back(voltagesPerTurn[turnIndex]);
        sumCV += capacitance * voltagesPerTurn[turnIndex];
        sumC += capacitance;
    }
    if (sumC <= 0) {
        return 0;
    }
    double corePotential = sumCV / sumC;   // floating node: charge balance
    double energy = 0;
    for (size_t k = 0; k < turnCoreCapacitance.size(); ++k) {
        energy += 0.5 * turnCoreCapacitance[k] * std::pow(turnPotential[k] - corePotential, 2);
    }
    return energy;
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

    // Compute global minimum gap once for all turns
    double globalMinimumGap = compute_global_minimum_gap(turns);

    std::set<std::pair<size_t, size_t>> turnsCombinations;

    // ABT #853: corresponding turns of two parallels of the SAME winding are one electrical
    // node -- the parallels join at the terminals, and the per-turn voltage divider
    // (calculate_voltages_per_turn) hands turn k of every parallel the identical potential.
    // There is no voltage across such a pair, hence no capacitor to compute and nothing that
    // can short, however tightly a bifilar/trifilar pair is wound. get_surrounding_turns()
    // filters by coordinates only, so without this the pair reached the short-circuit guard
    // as an ordinary neighbour and a bare (or zero-derived) coating declared the winding
    // shorted. The ordinal is counted exactly as the voltage divider counts it (per winding,
    // per parallel, in turns order) so "same node" here means "same potential" there.
    std::map<std::string, std::map<int64_t, size_t>> nextOrdinalPerWindingPerParallel;
    std::vector<size_t> ordinalWithinParallel;
    ordinalWithinParallel.reserve(turns.size());
    for (auto& turn : turns) {
        ordinalWithinParallel.push_back(nextOrdinalPerWindingPerParallel[turn.get_winding()][turn.get_parallel()]++);
    }
    auto areSameElectricalNode = [&](size_t firstIndex, size_t secondIndex) {
        return turns[firstIndex].get_winding() == turns[secondIndex].get_winding() &&
               turns[firstIndex].get_parallel() != turns[secondIndex].get_parallel() &&
               ordinalWithinParallel[firstIndex] == ordinalWithinParallel[secondIndex];
    };

    for (size_t turnIndex = 0; turnIndex <  turns.size(); ++turnIndex) {
        auto turnWindingIndex = coil.get_winding_index_by_name(turns[turnIndex].get_winding());
        auto turnWire = wirePerWinding[turnWindingIndex];
        auto surroundingTurns = OpenMagnetics::StrayCapacitance::get_surrounding_turns(turns[turnIndex], turns, globalMinimumGap);

        for (auto [surroundingTurn, surroundingTurnIndex] : surroundingTurns) {
            auto key = std::make_pair(turnIndex, surroundingTurnIndex);
            auto inverseKey = std::make_pair(surroundingTurnIndex, turnIndex);
            if (turnsCombinations.contains(key) || turnsCombinations.contains(inverseKey)) {
                continue;
            }
            if (areSameElectricalNode(turnIndex, surroundingTurnIndex)) {
                turnsCombinations.insert(key);
                continue;
            }
            auto surroundingTurnWindingIndex = coil.get_winding_index_by_name(surroundingTurn.get_winding());
            auto surroundingTurnWire = wirePerWinding[surroundingTurnWindingIndex];
            double capacitance = calculate_static_capacitance_between_two_turns(turns[turnIndex], turnWire, surroundingTurn, surroundingTurnWire, coil);
            capacitanceAmongTurns[key] = capacitance;
            capacitanceAmongTurns[inverseKey] = capacitance;
            turnsCombinations.insert(key);
        }
    }

    return capacitanceAmongTurns;
}


/**
 * @brief Calculates the 3x3 capacitance matrix between two windings.
 * 
 * Based on: "Using Transformer Parasitics for Resonant Converters—A Review of the 
 * Calculation of the Stray Capacitance of Transformers" by Juergen Biela and Johann W. Kolar
 * IEEE Transactions on Industry Applications, Vol. 44, No. 1, January/February 2008
 * Section V: "Calculation of the Equivalent Capacitance Network"
 * 
 * This method converts stored electric field energy between windings into a six-capacitor
 * network model, then transforms it to a 3x3 Maxwell capacitance matrix representation.
 * 
 * The six-capacitor network (γ1-γ6) represents capacitances between:
 * - Primary start/end terminals
 * - Secondary start/end terminals  
 * - Cross-winding capacitances
 * 
 * The transformation to 3-terminal representation uses the **turns ratio** (n = N1/N2)
 * to properly scale cross-winding terms, as capacitance depends on voltage ratios which
 * scale with turns ratio (unlike resistance matrix which uses inductance ratio √(L1/L2)).
 * 
 * Key equations from Biela/Kolar (Section V-A):
 * - C0 = 2 * E / V² (total equivalent capacitance from stored energy)
 * - γ1 = γ2 = -C0/6 (self-capacitance negative terms)
 * - γ3 = γ4 = C0/3 (mutual capacitance terms)
 * - γ5 = γ6 = C0/6 (cross-coupling terms)
 * 
 * Matrix elements with turns ratio (n) scaling:
 * - C[1][1] = γ1 + n*(γ4 + γ5)
 * - C[1][2] = -2*γ4
 * - C[1][3] = 2*n*γ5
 * - C[2][2] = γ2 + γ4 + γ6
 * - C[2][3] = 2*γ6
 * - C[3][3] = γ3 + γ5 + γ6
 * 
 * @param energy Total electric field energy stored between windings [J]
 * @param voltageDrop Voltage difference between the windings [V]
 * @param relativeTurnsRatio Turns ratio N_primary/N_secondary
 * @return 3x3 capacitance matrix representation at DC (frequency = 0)
 */
ScalarMatrixAtFrequency StrayCapacitance::calculate_capacitance_matrix_between_windings(double energy, double voltageDrop, double relativeTurnsRatio) {
    ScalarMatrixAtFrequency scalarMatrixAtFrequency;
    // C0: Total equivalent capacitance from energy method, C = 2E/V² (from E = ½CV²).
    // This assembles the internal three-input-multipole coefficients whose ["*"]["3"] entries drive
    // the floating-node (V3) convergence in calculate_capacitance_with_voltages. It is an internal
    // solver representation, NOT the terminal Maxwell matrix (see calculate_maxwell_capacitance_matrix).
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

// Canonical Biela & Kolar six-capacitor two-port network (IEEE Trans. Ind. Appl. 2008, eq. 30) for
// a winding pair, from the static inter-winding capacitance C0 by the electrostatic energy method.
// Node topology (their Fig. 16), terminals P1+ P1- P2+ P2-:
//   C1: P1+ - P1-  (primary self)     = -C0/6
//   C2: P2+ - P2-  (secondary self)   = -C0/6
//   C3: P1+ - P2+  (top-to-top)       =  C0/3
//   C4: P1- - P2-  (bottom-to-bottom) =  C0/3
//   C5: P1+ - P2-  (cross)            =  C0/6
//   C6: P1- - P2+  (cross)            =  C0/6
// The -C0/6 self terms are INTRINSIC to this complete energy-equivalent representation (a valid
// mathematical two-port, not a passive netlist). Do NOT emit this to transient SPICE — use the
// positive three-capacitor pi/tripole model for simulation.
SixCapacitorNetworkPerWinding StrayCapacitance::calculate_six_capacitor_network(double staticInterwindingCapacitance) {
    const double C0 = staticInterwindingCapacitance;
    SixCapacitorNetworkPerWinding n;
    n.set_c1(-C0 / 6.0);
    n.set_c2(-C0 / 6.0);
    n.set_c3(C0 / 3.0);
    n.set_c4(C0 / 3.0);
    n.set_c5(C0 / 6.0);
    n.set_c6(C0 / 6.0);
    return n;
}


StrayCapacitanceOutput StrayCapacitance::calculate_capacitance(Coil coil, std::optional<Core> core, std::optional<double> frequency) {
    std::map<std::string, double> voltageRmsPerWinding;
    double primaryNumberTurns = coil.get_functional_description()[0].get_number_turns();
    for (auto winding : coil.get_functional_description()) {
        double turnsRatio = primaryNumberTurns /  winding.get_number_turns();
        voltageRmsPerWinding[winding.get_name()] = 10.0 / turnsRatio;
    }
    return calculate_capacitance_with_voltages(coil, voltageRmsPerWinding, core, frequency);
}

StrayCapacitanceOutput StrayCapacitance::calculate_capacitance(Coil coil, OperatingPoint operatingPoint, std::optional<Core> core, std::optional<double> frequency) {
    // Extract actual RMS voltages from operating point
    std::map<std::string, double> voltageRmsPerWinding;
    
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        if (windingIndex >= operatingPoint.get_excitations_per_winding().size()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing excitation for winding index: " + std::to_string(windingIndex));
        }
        
        auto& excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        
        if (!excitation.get_voltage()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing voltage in excitation for winding: " + coil.get_functional_description()[windingIndex].get_name());
        }
        
        if (!excitation.get_voltage()->get_processed()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Voltage is not processed for winding: " + coil.get_functional_description()[windingIndex].get_name());
        }
        
        auto rmsVoltage = excitation.get_voltage()->get_processed()->get_rms();
        if (!rmsVoltage) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing RMS voltage value for winding: " + coil.get_functional_description()[windingIndex].get_name());
        }
        
        voltageRmsPerWinding[coil.get_functional_description()[windingIndex].get_name()] = rmsVoltage.value();
    }

    return calculate_capacitance_with_voltages(coil, voltageRmsPerWinding, core, frequency);
}

StrayCapacitanceOutput StrayCapacitance::calculate_capacitance_with_voltages(Coil coil, std::map<std::string, double> voltageRmsPerWinding, std::optional<Core> core, std::optional<double> frequency) {
    std::map<std::pair<size_t, size_t>, double> electricEnergyBetweenTurnsMap;
    std::map<std::pair<size_t, size_t>, double> voltageDropBetweenTurnsMap;
    std::map<std::string, std::map<std::string, ScalarMatrixAtFrequency>> capacitanceMatrix;
    std::map<std::string, std::map<std::string, SixCapacitorNetworkPerWinding>> sixCapacitorNetworkPerWinding;
    std::map<std::string, std::map<std::string, TripoleCapacitancePerWinding>> tripoleCapacitancePerWinding;

    auto capacitanceAmongTurns = calculate_capacitance_among_turns(coil);

    auto strayCapacitanceOutput = StrayCapacitance::calculate_voltages_per_turn(coil, voltageRmsPerWinding);
    auto voltagesPerTurn = strayCapacitanceOutput.get_voltage_per_turn().value();
    auto windings = coil.get_functional_description();
    std::map<std::pair<std::string, std::string>, double> capacitanceMapPerWindings;
    for (auto firstWinding : windings) {
        auto firstWindingName = firstWinding.get_name();
        auto turnsInFirstWinding = coil.get_turns_indexes_by_winding(firstWindingName);
        for (auto secondWinding : windings) {
            auto secondWindingName = secondWinding.get_name();
            auto windingsKey = std::make_pair(firstWindingName, secondWindingName);
            if (capacitanceMapPerWindings.contains(windingsKey) || capacitanceMapPerWindings.contains(std::make_pair(secondWindingName, firstWindingName))) {
                continue;
            }

            double minVoltageInFirstWinding = std::numeric_limits<double>::max();
            double maxVoltageInFirstWinding = std::numeric_limits<double>::lowest();
            double minVoltageInSecondWinding = std::numeric_limits<double>::max();
            double maxVoltageInSecondWinding = std::numeric_limits<double>::lowest();

            double V3 = 0;
            double V3calculated = 0;
            bool firstConvergenceIteration = true;

            double energyInBetweenTheseWindings = 0;
            double voltageDropBetweenWindings = 0;
            double relativeTurnsRatio = 0;
            ScalarMatrixAtFrequency capacitanceMatrixBetweenWindings;
            while (firstConvergenceIteration || fabs(V3 - V3calculated) > 0.001 * std::max(fabs(V3), std::numeric_limits<double>::min())) {
                firstConvergenceIteration = false;
                energyInBetweenTheseWindings = 0;
                V3 = V3calculated;
                // double C0 = 0;
                auto turnsInSecondWinding = coil.get_turns_indexes_by_winding(secondWindingName);
                bool windingAreNotAdjacent = true;
                std::set<std::pair<size_t, size_t>> turnsCombinations;
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
                        auto inverseTurnsKey = std::make_pair(turnInSecondWinding, turnInFirstWinding);

                        if (turnsCombinations.contains(turnsKey) || turnsCombinations.contains(inverseTurnsKey)) {
                            continue;
                        }
                        turnsCombinations.insert(turnsKey);
                        if (capacitanceAmongTurns.contains(turnsKey)) {
                            windingAreNotAdjacent = false;
                            double voltageDropAmongTurns = V3 + firstTurnVoltage - secondTurnVoltage;
                            double energyBetweenTurns = 0.5 * capacitanceAmongTurns[turnsKey] * pow(voltageDropAmongTurns, 2);
                            energyInBetweenTheseWindings += energyBetweenTurns;
                            electricEnergyBetweenTurnsMap[turnsKey] = energyBetweenTurns;
                            voltageDropBetweenTurnsMap[turnsKey] = voltageDropAmongTurns;

                            if (std::isnan(energyInBetweenTheseWindings)) {
                                // ABT #395: say WHICH pair and WHICH factor, so a NaN arriving by a
                                // route the capacitance guard does not cover is still localisable
                                // instead of being an anonymous failure of the whole export.
                                throw NaNResultException(
                                    "Energy cannot be nan: turns " + std::to_string(turnInFirstWinding) + " (" +
                                    firstWindingName + ", " + std::to_string(firstTurnVoltage) + " V) and " +
                                    std::to_string(turnInSecondWinding) + " (" + secondWindingName + ", " +
                                    std::to_string(secondTurnVoltage) + " V) have capacitance " +
                                    std::to_string(capacitanceAmongTurns[turnsKey]) + " F across a voltage drop of " +
                                    std::to_string(voltageDropAmongTurns) + " V");
                            }
                        }
                    }
                }
                // ABT #848: a winding's self-capacitance is chain PLUS core. The pair sum
                // above only carries the turn-to-turn chain, which SHRINKS with turn count
                // (adjacent turns differ by V/(N-1)); the turn-to-core elements against the
                // floating core GROW with it, and on toroids they dominate — without them
                // the self term came out ~50-100x low against 107 measured WE chokes
                // (near-zero for single-layer windings). Same energy method and per-turn
                // elements as the through-core inter-winding path below.
                if (firstWindingName == secondWindingName && core) {
                    energyInBetweenTheseWindings += calculate_winding_to_core_self_energy(
                        coil, core.value(), firstWindingName, voltagesPerTurn, frequency);
                }

                if (windingAreNotAdjacent) {
                    // No turn of the first winding is adjacent to any turn of the second,
                    // so there is no direct turn-to-turn capacitive path between them (the
                    // separated-winding CMC case: the two windings sit on opposite arcs of
                    // the toroid). The only path is through the core: turn -> core -> turn.
                    // The energy method (potential-weighted, floating-core node) builds it
                    // from the per-turn turn-to-core elements. Requires the core (its coating
                    // sets the dielectric floor); without it we fall back to 0 (old behaviour).
                    double throughCore = 0;
                    if (core && firstWindingName != secondWindingName) {
                        throughCore = calculate_through_core_capacitance(coil, core.value(), firstWindingName, secondWindingName, voltagesPerTurn, frequency);
                    }
                    capacitanceMapPerWindings[windingsKey] = throughCore;
                    continue;
                }

                voltageDropBetweenWindings = maxVoltageInFirstWinding - minVoltageInSecondWinding + V3;
                relativeTurnsRatio = static_cast<double>(firstWinding.get_number_turns()) / secondWinding.get_number_turns(); 
                capacitanceMatrixBetweenWindings = calculate_capacitance_matrix_between_windings(energyInBetweenTheseWindings, voltageDropBetweenWindings, relativeTurnsRatio);
                if (firstWindingName != secondWindingName) {
                    V3calculated = fabs(-(resolve_dimensional_values(capacitanceMatrixBetweenWindings.get_mutable_magnitude()["1"]["3"]) * maxVoltageInFirstWinding + resolve_dimensional_values(capacitanceMatrixBetweenWindings.get_mutable_magnitude()["2"]["3"]) * fabs(minVoltageInSecondWinding)) / resolve_dimensional_values(capacitanceMatrixBetweenWindings.get_mutable_magnitude()["3"]["3"]));
                }
                // ABT #278: a single-turn winding has no intra-winding turn pairs, so its
                // self-capacitance term has zero energy AND zero voltage drop — the formula's 0/0
                // produced NaN (found by the example battery on the 1-turn current-sense primary).
                // Zero energy means zero capacitance exactly, whatever the voltage drop; a nonzero
                // energy with no voltage drop is inconsistent and must fail loudly, not divide.
                if (energyInBetweenTheseWindings == 0) {
                    capacitanceMapPerWindings[windingsKey] = 0;
                }
                else if (voltageDropBetweenWindings == 0) {
                    throw NaNResultException("Nonzero inter-turn energy with zero voltage drop between windings '"
                                             + firstWindingName + "' and '" + secondWindingName + "'");
                }
                else {
                    capacitanceMapPerWindings[windingsKey] = energyInBetweenTheseWindings * 2 / pow(voltageDropBetweenWindings, 2);
                }
            }
            capacitanceMatrix[firstWindingName][secondWindingName] = capacitanceMatrixBetweenWindings;
            capacitanceMatrix[secondWindingName][firstWindingName] = capacitanceMatrixBetweenWindings;
            // 6-capacitor network + positive 3-capacitor (tripole/pi) model are derived below, once
            // the full positive self+inter capacitance map is available (the tripole needs BOTH
            // windings' self-capacitances, which are only complete after this loop).
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

    // Positive 3-capacitor (pi / tripole) model + canonical 6-capacitor network per winding PAIR,
    // now that the full positive self+inter capacitance map is known. The tripole is the measurable,
    // SPICE-friendly model (Lu/Zhu/Hui 2003; Coilcraft/OMICRON vendor models): C1 = self-capacitance
    // of the first winding, C2 = self of the second, C3 = the inter-winding capacitance — all
    // positive, straight from the energy method (the same values the SPICE export uses). The 6C is
    // the complete Biela/Kolar energy-equivalent (with its intrinsic negative self terms) built from
    // the inter-winding static capacitance.
    for (size_t i = 0; i < windings.size(); ++i) {
        std::string wi = windings[i].get_name();
        for (size_t j = i + 1; j < windings.size(); ++j) {
            std::string wj = windings[j].get_name();
            if (!staticCapacitanceMapPerWindings.contains(wi) || !staticCapacitanceMapPerWindings.at(wi).contains(wj)) {
                continue;
            }
            double cinter = staticCapacitanceMapPerWindings.at(wi).at(wj);
            double cselfI = staticCapacitanceMapPerWindings.at(wi).contains(wi) ? staticCapacitanceMapPerWindings.at(wi).at(wi) : 0.0;
            double cselfJ = staticCapacitanceMapPerWindings.at(wj).contains(wj) ? staticCapacitanceMapPerWindings.at(wj).at(wj) : 0.0;

            TripoleCapacitancePerWinding tripole;
            tripole.set_c1(cselfI);
            tripole.set_c2(cselfJ);
            tripole.set_c3(cinter);
            tripoleCapacitancePerWinding[wi][wj] = tripole;
            tripoleCapacitancePerWinding[wj][wi] = tripole;

            auto sixCap = calculate_six_capacitor_network(cinter);
            sixCapacitorNetworkPerWinding[wi][wj] = sixCap;
            sixCapacitorNetworkPerWinding[wj][wi] = sixCap;
        }
    }


    std::map<std::string, std::map<std::string, double>> electricEnergyAmongTurns;
    std::map<std::string, std::map<std::string, double>> voltageDropAmongTurns;
    std::map<std::string, std::map<std::string, double>> capacitanceAmongTurnsOutput;
    auto turns = coil.get_turns_description().value();
    for (size_t firstTurnIndex = 0; firstTurnIndex < turns.size(); ++firstTurnIndex) {
        auto firstTurnName = turns[firstTurnIndex].get_name();
        for (size_t secondTurnIndex = 0; secondTurnIndex < turns.size(); ++secondTurnIndex) {
        // for (size_t secondTurnIndex = firstTurnIndex + 1; secondTurnIndex < turns.size(); ++secondTurnIndex) {
            auto secondTurnName = turns[secondTurnIndex].get_name();
            auto turnsKey = std::make_pair(firstTurnIndex, secondTurnIndex);
            auto inverseTurnsKey = std::make_pair(secondTurnIndex, firstTurnIndex);
            
            // Electric energy - check both key directions
            if (electricEnergyBetweenTurnsMap.contains(turnsKey)) {
                electricEnergyAmongTurns[firstTurnName][secondTurnName] = electricEnergyBetweenTurnsMap[turnsKey];
            } else if (electricEnergyBetweenTurnsMap.contains(inverseTurnsKey)) {
                electricEnergyAmongTurns[firstTurnName][secondTurnName] = electricEnergyBetweenTurnsMap[inverseTurnsKey];
            } else {
                electricEnergyAmongTurns[firstTurnName][secondTurnName] = 0;
            }
            
            // Voltage drop - check both key directions (with sign flip for inverse)
            if (voltageDropBetweenTurnsMap.contains(turnsKey)) {
                voltageDropAmongTurns[firstTurnName][secondTurnName] = voltageDropBetweenTurnsMap[turnsKey];
            } else if (voltageDropBetweenTurnsMap.contains(inverseTurnsKey)) {
                // For inverse key, negate the voltage drop (V_ab = -V_ba)
                voltageDropAmongTurns[firstTurnName][secondTurnName] = -voltageDropBetweenTurnsMap[inverseTurnsKey];
            } else {
                voltageDropAmongTurns[firstTurnName][secondTurnName] = 0;
            }
            
            // Capacitance - check both key directions (capacitance is symmetric)
            if (capacitanceAmongTurns.contains(turnsKey)) {
                capacitanceAmongTurnsOutput[firstTurnName][secondTurnName] = capacitanceAmongTurns[turnsKey];
            } else if (capacitanceAmongTurns.contains(inverseTurnsKey)) {
                capacitanceAmongTurnsOutput[firstTurnName][secondTurnName] = capacitanceAmongTurns[inverseTurnsKey];
            } else {
                capacitanceAmongTurnsOutput[firstTurnName][secondTurnName] = 0;
            }
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
    strayCapacitanceOutput.set_origin(ResultOrigin::SIMULATION);
    strayCapacitanceOutput.set_method_used(_model->methodName);

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
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Old code, this should not happen");
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
    // Iterative ladder (ABT #346). The former n-2 recursion with exact
    // n==2 / n==3 base cases never terminated for n < 2 or non-integer n
    // (0-turn placeholder coils), and legitimately-huge turn counts (open
    // magnetic circuits needing tens of thousands of turns) exceeded the
    // WASM call-depth budget and killed the whole Magnetic Adviser.
    if (n < 2 || n != std::floor(n)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Stray capacitance ladder cab() needs an integer number of turns >= 2, got " + std::to_string(n));
    }
    bool even = std::fmod(n, 2) == 0;
    double value = even ? (ctt + cts / 2) : (ctt / 2 + cts / 2);
    for (double k = even ? 2 : 3; k < n; k += 2) {
        value = (value * ctt / 2) / (value + ctt / 2) + cts / 2;
    }
    return value;
}

double cas(double n, double ctt, double cts) {
    // Iterative ladder (ABT #346), same reasoning as cab().
    if (n < 1 || n != std::floor(n)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Stray capacitance ladder cas() needs an integer number of turns >= 1, got " + std::to_string(n));
    }
    double value = cts;
    for (double k = 1; k < n; k += 1) {
        value = (value * ctt) / (value + ctt) + cts;
    }
    return value;
}

double StrayCapacitanceOneLayer::calculate_capacitance(Coil coil, std::optional<Core> core) {
    // Baed on https://sci-hub.st/https://ieeexplore.ieee.org/document/793378
    double numberTurns = coil.get_functional_description()[0].get_number_turns();
    auto wireRadius = coil.resolve_wire(0).get_maximum_conducting_width() / 2;
    double columnThickness = coil.resolve_bobbin().get_processed_description()->get_column_thickness();
    double distanceTurnsToCore = columnThickness + coil.resolve_wire(0).get_maximum_outer_width() / 2;

    // ABT #846: capacitance_turn_to_turn / capacitance_turn_to_shield are the textbook
    // per-unit-length forms (pi*e0/acosh, 2*pi*e0/acosh) times pi*turnDiameter, so the
    // parameter must be a DIAMETER — the turn length is then pi*D. column_width is the
    // column RADIUS (Bobbin sets it to coreColumn.width/2 + columnThickness), so the turn
    // wraps at radius (column_width + wireRadius) and its diameter is twice that. The old
    // code passed the CIRCUMFERENCE 2*pi*(column_width + wireRadius), making every
    // capacitance pi times too large (verified against Medhurst's 1947 solenoid data).
    double turnDiameter = 2 * (coil.resolve_bobbin().get_processed_description()->get_column_width().value() + wireRadius);

    // A toroid's winding window is ROUND; its turn does not circle a round column but wraps
    // the core CROSS-SECTION — length 2*(radialThickness + height) plus the rounding of the
    // four corners at radius (columnThickness + wireRadius). The core height never enters
    // the bobbin's column data, so this needs the core itself (ABT #845/#846).
    bool toroidal = coil.resolve_bobbin().get_winding_window_shape(0) == WindingWindowShape::ROUND;
    if (toroidal && numberTurns > 1) {
        // The toroid model below is only correct with the core's cross-section: the bobbin
        // column carries no height, so without the core the turn length would silently fall
        // back to the column-based circle — neither the old behaviour nor the new one.
        // Require it loudly instead of guessing (the optional default serves non-toroids).
        if (!core || !core->get_processed_description()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "StrayCapacitanceOneLayer needs the processed core for a toroidal coil: the "
                "turn wraps the core cross-section (width and height), which the bobbin's "
                "column data cannot provide. Pass the core, or process its description first.");
        }
        auto column = core->get_processed_description()->get_columns().at(0);
        double turnLength = 2 * (column.get_width() + column.get_depth()) + 2 * std::numbers::pi * (columnThickness + wireRadius);
        turnDiameter = turnLength / std::numbers::pi;
    }
    double centerSeparation = coil.resolve_wire(0).get_maximum_outer_width();
    if (coil.get_turns_description()) {
        if (coil.get_turns_description().value().size() > 1) {
            auto firstTurn = coil.get_turns_description().value()[0];
            auto secondTurn = coil.get_turns_description().value()[1];
            centerSeparation = sqrt(pow(firstTurn.get_coordinates()[0] - secondTurn.get_coordinates()[0], 2) + pow(firstTurn.get_coordinates()[1] - secondTurn.get_coordinates()[1], 2));
        }
    }

    // The two building blocks are cylinder-to-cylinder / cylinder-to-plane
    // capacitances, C ~ 1/acosh(x), valid only for x > 1 (conductors not touching):
    //   turn-to-turn:   x = centerSeparation / (2*wireRadius)
    //   turn-to-shield: x = distanceTurnsToCore / wireRadius
    // x -> 1 is the contact singularity (acosh(1)=0 => C=inf => the cas/cab
    // recursion produces NaN). This happens for a flat/planar conductor (whose
    // half-WIDTH is used as "wireRadius" while the stacking pitch is the thin
    // dimension) and for thick round Grade-1 wire flush to the core (proportionally
    // thin coating + zero bobbin column thickness => distanceTurnsToCore == wireRadius).
    // In those degenerate cases this fast single-layer round-wire approximation is
    // outside its validity domain, so defer to the full StrayCapacitance model
    // (parallel-plate for planar, Albach for round) which is finite and wire-type aware.
    // The 1.001 margin matches the acosh-domain guard the Albach/Koch models use and
    // also avoids the near-singular huge-but-finite value just above contact.
    constexpr double acoshDomainFloor = 1.001;
    auto wireType = coil.resolve_wire(0).get_type();
    bool roundLike = (wireType == WireType::ROUND || wireType == WireType::LITZ);
    double turnToTurnArg = centerSeparation / (2 * wireRadius);
    double turnToShieldArg = distanceTurnsToCore / wireRadius;
    if (!roundLike || turnToTurnArg <= acoshDomainFloor || turnToShieldArg <= acoshDomainFloor) {
        auto capacitanceAmongWindings = StrayCapacitance().calculate_capacitance(coil).get_capacitance_among_windings().value();
        auto firstWindingName = coil.get_functional_description()[0].get_name();
        return capacitanceAmongWindings[firstWindingName][firstWindingName];
    }

    double ctt = capacitance_turn_to_turn(turnDiameter, wireRadius, centerSeparation);
    double cts = capacitance_turn_to_shield(turnDiameter, wireRadius, distanceTurnsToCore);

    // ABT #845: the cas/cab ladder models a GROUNDED shield, which screens distant turns and
    // makes the result converge to the fixed point sqrt(cts*ctt), independent of turn count
    // above ~15 turns. A magnetic core has no terminal: it is a FLOATING conductor that
    // settles at the potential where its displacement currents balance (Pasko, Kazimierczuk,
    // Grzesik, IEEE Trans. EMC 57(2), 2015). For a single-layer toroid with the standard
    // linear voltage distribution V_k = V*k/(N-1), charge balance puts the core at V/2 and
    // the stored energies reduce to two closed forms:
    //   turn-to-turn chain: (N-1) series gaps of V/(N-1) -> C_chain = ctt / (N-1)
    //   turn-to-core:       (1/2)*cts*Sum((k/(N-1) - 1/2)^2)*V^2 -> C_core ~ N*cts/12
    // The core term GROWS with N — measured WE common-mode chokes scale as C ~ N^+1.3 while
    // the grounded ladder scales as N^-0.05, resonating 3x high on many-turn parts and
    // parking a spurious low resonance on few-turn parts (ABT #845 has the data).
    //
    // N=1 deliberately falls through to the ladder below (which reduces to cts): the chain
    // term ctt/(N-1) is singular at one turn, and a single turn has no turn-to-turn chain
    // for the floating-core reduction to act on. The model discontinuity at N=2 is known
    // and intended — do not "fix" the guard to >= 1.
    if (toroidal && numberTurns > 1) {
        double capacitance = ctt / (numberTurns - 1) + numberTurns * cts / 12;
        if (std::isnan(capacitance)) {
            throw NaNResultException("capacitance cannot be NaN");
        }
        // Deliberately NO layers multiplier here (unlike the ladder path below): the
        // floating-core energy sum already counts every turn once, and on a multi-layer
        // toroid the outer layers sit FARTHER from the core, weakening — not multiplying —
        // their coupling. Applying the legacy x-layers heuristic doubled the capacitance of
        // the 110-turn T12.5/7.5/5 anchor (2 layers) and moved its self-resonance from the
        // measured ~180 kHz to 130 kHz; without it the single-layer form lands on the pin.
        return capacitance;
    }

    double C2;
    double casValue = cas(numberTurns, ctt, cts);

    if (std::isnan(casValue)) {
        throw NaNResultException("capacitance cannot be NaN");
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


// ==================== Bipolar Coordinate System for Round-Round Pairs ====================

StrayCapacitance::BipolarParams StrayCapacitance::compute_bipolar_params(const Turn& t1, const Turn& t2) {
    BipolarParams params;
    auto& coords1 = t1.get_coordinates();
    auto& coords2 = t2.get_coordinates();
    if (coords1.size() < 2 || coords2.size() < 2) {
        params.focalHalfDistance = 0; params.tau1 = 0; params.tau2 = 0;
        params.cosAngle = 1; params.sinAngle = 0;
        params.midX = 0; params.midY = 0;
        return params;
    }
    double cx1 = coords1[0], cy1 = coords1[1];
    double cx2 = coords2[0], cy2 = coords2[1];
    double dx = cx2 - cx1, dy = cy2 - cy1;
    double c = hypot(dx, dy);
    params.midX = (cx1 + cx2) / 2.0;
    params.midY = (cy1 + cy2) / 2.0;
    if (c < 1e-15) {
        params.focalHalfDistance = 0; params.tau1 = 0; params.tau2 = 0;
        params.cosAngle = 1; params.sinAngle = 0;
        return params;
    }
    params.cosAngle = dx / c;
    params.sinAngle = dy / c;
    if (!t1.get_dimensions() || !t2.get_dimensions() ||
        t1.get_dimensions().value().empty() || t2.get_dimensions().value().empty()) {
        params.focalHalfDistance = 0; params.tau1 = 0; params.tau2 = 0;
        return params;
    }
    double r1 = t1.get_dimensions().value()[0] / 2.0;
    double r2 = t2.get_dimensions().value()[0] / 2.0;
    double r1mr2 = r1 - r2, r1pr2 = r1 + r2, c2 = c * c;
    double term1 = c2 - r1mr2 * r1mr2;
    double term2 = c2 - r1pr2 * r1pr2;
    double a2 = term1 * term2 / (4.0 * c2);
    if (a2 <= 0) {
        params.focalHalfDistance = 0; params.tau1 = 0; params.tau2 = 0;
        return params;
    }
    params.focalHalfDistance = sqrt(a2);
    double coshTau1 = std::max(1.0, (c2 + r1 * r1 - r2 * r2) / (2.0 * r1 * c));
    double coshTau2 = std::max(1.0, (c2 + r2 * r2 - r1 * r1) / (2.0 * r2 * c));
    params.tau1 = -acosh(coshTau1);
    params.tau2 =  acosh(coshTau2);
    return params;
}

double StrayCapacitance::bipolar_tau_at_point(double lx, double ly, double a) {
    if (a < 1e-15) return 0;
    double dPlus2  = (lx + a) * (lx + a) + ly * ly;
    double dMinus2 = (lx - a) * (lx - a) + ly * ly;
    if (dPlus2 < 1e-30 || dMinus2 < 1e-30) return 0;
    return 0.5 * log(dPlus2 / dMinus2);
}

double StrayCapacitance::bipolar_energy_density_at_point(
    double px, double py, const BipolarParams& params,
    double voltageDrop, double epsilonEff)
{
    double a = params.focalHalfDistance;
    if (a < 1e-15) return 0;
    double tauDiff = params.tau2 - params.tau1;
    if (fabs(tauDiff) < 1e-15) return 0;
    double lx =  params.cosAngle * (px - params.midX) + params.sinAngle * (py - params.midY);
    double ly = -params.sinAngle * (px - params.midX) + params.cosAngle * (py - params.midY);
    double tau = bipolar_tau_at_point(lx, ly, a);
    double denomSigma = lx * lx + ly * ly - a * a;
    double sigma = atan2(2.0 * a * ly, denomSigma);
    double hDenom = cosh(tau) - cos(sigma);
    if (fabs(hDenom) < 1e-15) return 0;
    double h = a / hDenom;
    double gradV2 = (voltageDrop * voltageDrop) / (tauDiff * tauDiff * h * h);
    return 0.5 * epsilonEff * gradV2;
}


} // namespace OpenMagnetics
