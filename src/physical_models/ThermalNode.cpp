#include "ThermalNode.h"
#include "ThermalResistance.h"
#include "Definitions.h"
#include "Wire.h"
#include <numbers>
#include <algorithm>

namespace OpenMagnetics {

// ============================================================================
// WireCoatingUtils Implementation
// ============================================================================

double WireCoatingUtils::calculateCoatingResistance(const InsulationWireCoating& coating, double area) {
    if (area <= 0.0) return 0.0;
    
    // Use Wire::get_coating_thickness() and Wire::get_coating_thermal_conductivity()
    double thickness = Wire::get_coating_thickness(coating);
    double k = Wire::get_coating_thermal_conductivity(coating);
    
    if (thickness <= 0.0 || k <= 0.0) {
        return 0.0;  // No coating resistance
    }
    
    return thickness / (k * area);
}

// ============================================================================
// ThermalNodeQuadrant Implementation
// ============================================================================

double ThermalNodeQuadrant::calculateConductionResistance(const ThermalNodeQuadrant& other, double contactArea) const {
    if (contactArea <= 0.0) {
        return 1e9;  // Very high resistance for no contact
    }
    
    // =========================================================================
    // IMP-NEW-01: Proper series conduction resistance R = R1 + R2
    // =========================================================================
    // WHAT: Calculate individual resistances for each material in series.
    // WHY:  Previous code used avgLength and min(k1,k2) which artificially
    //       increases resistance when materials differ significantly.
    //       Correct: R_total = L1/(k1*A) + L2/(k2*A)
    // PREVIOUS: avgLength/(minK*A) approximation
    // =========================================================================
    double R1 = ThermalResistance::calculateConductionResistance(
        this->length, this->thermalConductivity, contactArea);
    double R2 = ThermalResistance::calculateConductionResistance(
        other.length, other.thermalConductivity, contactArea);
    return R1 + R2;
}

double ThermalNodeQuadrant::calculateConvectionResistance(double heatTransferCoefficient) const {
    // Effective surface area accounts for coverage by turns
    // surfaceCoverage = 1.0: fully exposed, surfaceCoverage = 0.0: fully covered
    double effectiveArea = surfaceArea * surfaceCoverage;
    
    if (heatTransferCoefficient <= 0.0 || effectiveArea <= 0.0) {
        return 1e9;  // Very high resistance (no convection if fully covered)
    }
    return 1.0 / (heatTransferCoefficient * effectiveArea);
}

// ============================================================================
// ThermalNodeQuadrant Implementation
// ============================================================================

json ThermalNodeQuadrant::toJson() const {
    json j;
    j["face"] = magic_enum::enum_name(face);
    j["surfaceArea"] = surfaceArea;
    j["length"] = length;
    j["thermalConductivity"] = thermalConductivity;
    j["surfaceCoverage"] = surfaceCoverage;
    if (coating.has_value()) {
        // Use Wire::get_coating_thickness() and Wire::get_coating_thermal_conductivity()
        j["coating"]["thickness"] = Wire::get_coating_thickness(coating.value());
        j["coating"]["thermalConductivity"] = Wire::get_coating_thermal_conductivity(coating.value());
    }
    if (connectedNodeId.has_value()) {
        j["connectedNodeId"] = connectedNodeId.value();
        j["connectedQuadrant"] = magic_enum::enum_name(connectedQuadrant.value());
    }
    j["connectionType"] = magic_enum::enum_name(connectionType);
    return j;
}

// ============================================================================
// ThermalNetworkNode Implementation
// ============================================================================

ThermalNodeQuadrant* ThermalNetworkNode::getQuadrant(ThermalNodeFace face) {
    for (auto& q : quadrants) {
        if (q.face == face) {
            return &q;
        }
    }
    return nullptr;
}

const ThermalNodeQuadrant* ThermalNetworkNode::getQuadrant(ThermalNodeFace face) const {
    for (const auto& q : quadrants) {
        if (q.face == face) {
            return &q;
        }
    }
    return nullptr;
}

void ThermalNetworkNode::initializeToroidalQuadrants(double wireWidth, double wireHeight, 
                                                      double turnLength,
                                                      double thermalCond, 
                                                      bool isInner,
                                                      double centerRadius,
                                                      const std::optional<InsulationWireCoating>& wireCoating,
                                                      TurnCrossSectionalShape shape) {
    isInnerTurn = isInner;
    
    // Store geometry
    crossSectionalShape = shape;
    if (shape == TurnCrossSectionalShape::ROUND) {
        dimensions = NodeDimensions::cylindrical((wireWidth + wireHeight) / 2.0, turnLength / 2.0);
    } else {
        dimensions = NodeDimensions(wireWidth, wireHeight, turnLength / 2.0);
    }
    
    // For toroidal: each node represents half the turn (inner or outer)
    // Calculate radii for inner and outer faces
    // If centerRadius not provided, estimate from turnLength = 2*pi*R
    double rCenter = (centerRadius > 0) ? centerRadius : turnLength / (2.0 * std::numbers::pi);
    double rInner = rCenter - wireWidth / 2.0;   // Radius to inner face
    double rOuter = rCenter + wireWidth / 2.0;   // Radius to outer face
    
    // Ensure radii are positive
    rInner = std::max(rInner, 1e-6);
    rOuter = std::max(rOuter, 1e-6);
    
    // Length at each face = 2 * pi * radius (full turn) / 2 (half turn)
    double lengthInner = std::numbers::pi * rInner;   // Length of inner face
    double lengthOuter = std::numbers::pi * rOuter;   // Length of outer face
    double lengthTangential = std::numbers::pi * rCenter;  // Length at center for tangential
    
    // Surface areas:
    // - RADIAL_INNER face: height x lengthInner
    // - RADIAL_OUTER face: height x lengthOuter  
    // - TANGENTIAL faces: width x lengthTangential
    double radialInnerArea = wireHeight * lengthInner;
    double radialOuterArea = wireHeight * lengthOuter;
    double tangentialFaceArea = wireWidth * lengthTangential;
    
    // Get node position for limit coordinate calculations
    // IMP-9
    double nodeX = safeCoord(physicalCoordinates, 0);
    double nodeY = safeCoord(physicalCoordinates, 1);
    double nodeAngle = std::atan2(nodeY, nodeX);
    
    // Quadrant 0: RADIAL_INNER (facing toward center) - shorter length
    quadrants[0].face = ThermalNodeFace::RADIAL_INNER;
    quadrants[0].surfaceArea = radialInnerArea;
    quadrants[0].length = lengthInner;
    quadrants[0].thermalConductivity = thermalCond;
    quadrants[0].coating = wireCoating;
    // Limit coordinate: toward center by wireWidth/2 from node position
    {
        double xOffset = (wireWidth / 2.0) * std::cos(nodeAngle + std::numbers::pi);  // Toward center
        double yOffset = (wireWidth / 2.0) * std::sin(nodeAngle + std::numbers::pi);
        quadrants[0].limitCoordinates = {nodeX + xOffset, nodeY + yOffset, 0.0};
    }
    
    // Quadrant 1: RADIAL_OUTER (facing away from center) - longer length
    quadrants[1].face = ThermalNodeFace::RADIAL_OUTER;
    quadrants[1].surfaceArea = radialOuterArea;
    quadrants[1].length = lengthOuter;
    quadrants[1].thermalConductivity = thermalCond;
    quadrants[1].coating = wireCoating;
    // Limit coordinate: away from center by wireWidth/2 from node position
    {
        double xOffset = (wireWidth / 2.0) * std::cos(nodeAngle);  // Away from center
        double yOffset = (wireWidth / 2.0) * std::sin(nodeAngle);
        quadrants[1].limitCoordinates = {nodeX + xOffset, nodeY + yOffset, 0.0};
    }
    
    // Quadrant 2: TANGENTIAL_LEFT (facing CCW along winding)
    quadrants[2].face = ThermalNodeFace::TANGENTIAL_LEFT;
    quadrants[2].surfaceArea = tangentialFaceArea;
    quadrants[2].length = lengthTangential;
    quadrants[2].thermalConductivity = thermalCond;
    quadrants[2].coating = wireCoating;
    // Limit coordinate: at the tangential border of the wire (wireHeight/2 from center in tangential direction)
    {
        double angleLeft = nodeAngle + std::numbers::pi / 2.0;  // 90° CCW
        double xOffset = (wireHeight / 2.0) * std::cos(angleLeft);
        double yOffset = (wireHeight / 2.0) * std::sin(angleLeft);
        quadrants[2].limitCoordinates = {nodeX + xOffset, nodeY + yOffset, 0.0};
    }
    
    // Quadrant 3: TANGENTIAL_RIGHT (facing CW along winding)
    quadrants[3].face = ThermalNodeFace::TANGENTIAL_RIGHT;
    quadrants[3].surfaceArea = tangentialFaceArea;
    quadrants[3].length = lengthTangential;
    quadrants[3].thermalConductivity = thermalCond;
    quadrants[3].coating = wireCoating;
    // Limit coordinate: at the tangential border of the wire (wireHeight/2 from center in tangential direction)
    {
        double angleRight = nodeAngle - std::numbers::pi / 2.0;  // 90° CW
        double xOffset = (wireHeight / 2.0) * std::cos(angleRight);
        double yOffset = (wireHeight / 2.0) * std::sin(angleRight);
        quadrants[3].limitCoordinates = {nodeX + xOffset, nodeY + yOffset, 0.0};
    }
}

void ThermalNetworkNode::initializeConcentricQuadrant(double wireWidth, double wireHeight,
                                                       double turnLength,
                                                       double thermalCond,
                                                       const std::optional<InsulationWireCoating>& wireCoating) {
    isInnerTurn = false;
    
    // Store geometry
    crossSectionalShape = TurnCrossSectionalShape::RECTANGULAR;
    dimensions = NodeDimensions(wireWidth, wireHeight, turnLength);
    
    // For concentric: single node represents the full turn
    // Only use quadrants[0] with face = NONE (indicates non-toroidal)
    
    // Calculate total surface area of the wire
    // For rectangular wire: 2*(width*height) + 2*(width+height)*length
    // But we approximate as the lateral surface area: perimeter * length
    double perimeter = 2.0 * (wireWidth + wireHeight);
    double totalSurfaceArea = perimeter * turnLength;
    
    quadrants[0].face = ThermalNodeFace::NONE;
    quadrants[0].surfaceArea = totalSurfaceArea;
    quadrants[0].length = turnLength;
    quadrants[0].thermalConductivity = thermalCond;
    quadrants[0].coating = wireCoating;
    
    // Other quadrants are unused (face = NONE, area = 0)
    for (size_t i = 1; i < 4; ++i) {
        quadrants[i].face = ThermalNodeFace::NONE;
        quadrants[i].surfaceArea = 0.0;
        quadrants[i].length = 0.0;
        quadrants[i].thermalConductivity = thermalCond;
        quadrants[i].coating = wireCoating;
    }
}

void ThermalNetworkNode::initializeConcentricTurnQuadrants(double wireWidth, double wireHeight, 
                                                            double turnLength,
                                                            double thermalCond,
                                                            const std::optional<InsulationWireCoating>& wireCoating,
                                                            TurnCrossSectionalShape shape) {
    isInnerTurn = false;
    
    // Store geometry
    crossSectionalShape = shape;
    dimensions = NodeDimensions(wireWidth, wireHeight, turnLength);
    
    // For concentric turns: quadrants map to cardinal directions (Cartesian-aligned)
    // This is different from toroidal where directions are angle-dependent
    //
    // Quadrant mapping for concentric turns:
    // - RADIAL_INNER  (index 0): LEFT face (-X) - facing toward bobbin/center
    // - RADIAL_OUTER  (index 1): RIGHT face (+X) - facing away from center
    // - TANGENTIAL_LEFT  (index 2): TOP face (+Y) - facing up
    // - TANGENTIAL_RIGHT (index 3): BOTTOM face (-Y) - facing down
    
    double sideArea = wireHeight * turnLength;    // LEFT/RIGHT faces: height * length
    double topBottomArea = wireWidth * turnLength; // TOP/BOTTOM faces: width * length
    
    // Get node center position
    // IMP-9
    double nodeX = safeCoord(physicalCoordinates, 0);
    double nodeY = safeCoord(physicalCoordinates, 1);
    
    // LEFT face (-X direction) - toward center/bobbin
    quadrants[0].face = ThermalNodeFace::RADIAL_INNER;
    quadrants[0].surfaceArea = sideArea;
    quadrants[0].length = turnLength;
    quadrants[0].thermalConductivity = thermalCond;
    quadrants[0].coating = wireCoating;
    quadrants[0].limitCoordinates = {nodeX - wireWidth/2, nodeY, 0.0};
    
    // RIGHT face (+X direction) - away from center
    quadrants[1].face = ThermalNodeFace::RADIAL_OUTER;
    quadrants[1].surfaceArea = sideArea;
    quadrants[1].length = turnLength;
    quadrants[1].thermalConductivity = thermalCond;
    quadrants[1].coating = wireCoating;
    quadrants[1].limitCoordinates = {nodeX + wireWidth/2, nodeY, 0.0};
    
    // TOP face (+Y direction)
    quadrants[2].face = ThermalNodeFace::TANGENTIAL_LEFT;
    quadrants[2].surfaceArea = topBottomArea;
    quadrants[2].length = turnLength;
    quadrants[2].thermalConductivity = thermalCond;
    quadrants[2].coating = wireCoating;
    quadrants[2].limitCoordinates = {nodeX, nodeY + wireHeight/2, 0.0};
    
    // BOTTOM face (-Y direction)
    quadrants[3].face = ThermalNodeFace::TANGENTIAL_RIGHT;
    quadrants[3].surfaceArea = topBottomArea;
    quadrants[3].length = turnLength;
    quadrants[3].thermalConductivity = thermalCond;
    quadrants[3].coating = wireCoating;
    quadrants[3].limitCoordinates = {nodeX, nodeY - wireHeight/2, 0.0};
}

void ThermalNetworkNode::initializeConcentricCoreQuadrants(double width, double height, double depth, 
                                                           double thermalCond) {
    // Store geometry
    dimensions = NodeDimensions(width, height, depth);
    
    // For concentric cores, map quadrants to cardinal directions:
    // TANGENTIAL_LEFT  (index 2) ↔ TOP    (+Y direction)
    // TANGENTIAL_RIGHT (index 3) ↔ BOTTOM (-Y direction)
    // RADIAL_INNER     (index 1) ↔ LEFT   (-X direction, toward center)
    // RADIAL_OUTER     (index 0) ↔ RIGHT  (+X direction, away from center)
    
    // Calculate surface areas for each face
    // RIGHT/LEFT faces: height * depth
    double sideArea = height * depth;
    // TOP/BOTTOM faces: width * depth
    double topBottomArea = width * depth;
    
    // RIGHT face (+X direction)
    quadrants[0].face = ThermalNodeFace::RADIAL_OUTER;
    quadrants[0].surfaceArea = sideArea;
    quadrants[0].length = depth;
    quadrants[0].thermalConductivity = thermalCond;
    
    // LEFT face (-X direction, toward center)
    quadrants[1].face = ThermalNodeFace::RADIAL_INNER;
    quadrants[1].surfaceArea = sideArea;
    quadrants[1].length = depth;
    quadrants[1].thermalConductivity = thermalCond;
    
    // TOP face (+Y direction)
    quadrants[2].face = ThermalNodeFace::TANGENTIAL_LEFT;
    quadrants[2].surfaceArea = topBottomArea;
    quadrants[2].length = width;
    quadrants[2].thermalConductivity = thermalCond;
    
    // BOTTOM face (-Y direction)
    quadrants[3].face = ThermalNodeFace::TANGENTIAL_RIGHT;
    quadrants[3].surfaceArea = topBottomArea;
    quadrants[3].length = width;
    quadrants[3].thermalConductivity = thermalCond;
    
    // Set limit coordinates for each quadrant (relative to node center)
    // IMP-9
    double nodeX = safeCoord(physicalCoordinates, 0);
    double nodeY = safeCoord(physicalCoordinates, 1);
    
    // For toroidal insulation layers: limit coordinates should be angled based on node position
    // Calculate angle from origin to node center
    double nodeAngle = std::atan2(nodeY, nodeX);
    double cosA = std::cos(nodeAngle);
    double sinA = std::sin(nodeAngle);
    
    // width = radial thickness, height = arc length (tangential)
    // RADIAL_OUTER: away from center (in direction of nodeAngle)
    quadrants[0].limitCoordinates = {nodeX + cosA * width/2, nodeY + sinA * width/2, 0.0};
    // RADIAL_INNER: toward center (opposite to nodeAngle)
    quadrants[1].limitCoordinates = {nodeX - cosA * width/2, nodeY - sinA * width/2, 0.0};
    // TANGENTIAL_LEFT: perpendicular, counterclockwise from radial (+90°)
    quadrants[2].limitCoordinates = {nodeX - sinA * height/2, nodeY + cosA * height/2, 0.0};
    // TANGENTIAL_RIGHT: perpendicular, clockwise from radial (-90°)
    quadrants[3].limitCoordinates = {nodeX + sinA * height/2, nodeY - cosA * height/2, 0.0};
}

void ThermalNetworkNode::initializeInsulationLayerQuadrants(double width, double height, double depth,
                                                             double thermalCond) {
    // Store geometry
    dimensions = NodeDimensions(width, height, depth);
    crossSectionalShape = TurnCrossSectionalShape::RECTANGULAR;
    
    // Insulation layers are rectangular nodes with 4 faces for heat conduction.
    // Unlike core/bobbin nodes, insulation layers can have turns connected to ANY face.
    // 
    // Quadrant mapping:
    // Index 0 (RADIAL_OUTER): RIGHT face (+X) - facing outer turns
    // Index 1 (RADIAL_INNER): LEFT face (-X) - facing inner turns/core
    // Index 2 (TANGENTIAL_LEFT): TOP face (+Y) - facing top turns/yoke
    // Index 3 (TANGENTIAL_RIGHT): BOTTOM face (-Y) - facing bottom turns/yoke
    
    // Calculate surface areas for each face
    // RIGHT/LEFT faces: height * depth
    double sideArea = height * depth;
    // TOP/BOTTOM faces: width * depth
    double topBottomArea = width * depth;
    
    // RIGHT face (+X direction) - facing outward toward outer turns
    quadrants[0].face = ThermalNodeFace::RADIAL_OUTER;
    quadrants[0].surfaceArea = sideArea;
    quadrants[0].length = depth;
    quadrants[0].thermalConductivity = thermalCond;
    
    // LEFT face (-X direction) - facing inward toward inner turns/core
    quadrants[1].face = ThermalNodeFace::RADIAL_INNER;
    quadrants[1].surfaceArea = sideArea;
    quadrants[1].length = depth;
    quadrants[1].thermalConductivity = thermalCond;
    
    // TOP face (+Y direction) - facing toward top turns/yoke
    quadrants[2].face = ThermalNodeFace::TANGENTIAL_LEFT;
    quadrants[2].surfaceArea = topBottomArea;
    quadrants[2].length = width;
    quadrants[2].thermalConductivity = thermalCond;
    
    // BOTTOM face (-Y direction) - facing toward bottom turns/yoke
    quadrants[3].face = ThermalNodeFace::TANGENTIAL_RIGHT;
    quadrants[3].surfaceArea = topBottomArea;
    quadrants[3].length = width;
    quadrants[3].thermalConductivity = thermalCond;
    
    // Set limit coordinates for each quadrant (relative to node center)
    // IMP-9
    double nodeX = safeCoord(physicalCoordinates, 0);
    double nodeY = safeCoord(physicalCoordinates, 1);
    
    // For toroidal insulation layers: limit coordinates should be angled based on node position
    // Calculate angle from origin to node center
    double nodeAngle = std::atan2(nodeY, nodeX);
    double cosA = std::cos(nodeAngle);
    double sinA = std::sin(nodeAngle);
    
    // width = radial thickness, height = arc length (tangential)
    // RADIAL_OUTER: away from center (in direction of nodeAngle)
    quadrants[0].limitCoordinates = {nodeX + cosA * width/2, nodeY + sinA * width/2, 0.0};
    // RADIAL_INNER: toward center (opposite to nodeAngle)
    quadrants[1].limitCoordinates = {nodeX - cosA * width/2, nodeY - sinA * width/2, 0.0};
    // TANGENTIAL_LEFT: perpendicular, counterclockwise from radial (+90 deg)
    quadrants[2].limitCoordinates = {nodeX - sinA * height/2, nodeY + cosA * height/2, 0.0};
    // TANGENTIAL_RIGHT: perpendicular, clockwise from radial (-90 deg)
    quadrants[3].limitCoordinates = {nodeX + sinA * height/2, nodeY - cosA * height/2, 0.0};
}

// ============================================================================
// Surface Coverage Calculations
// ============================================================================

double ThermalNetworkNode::calculateToroidalSurfaceCoverage(double coreRadius, double segmentAngle,
                                                             const std::vector<double>& turnWidths) {
    // Calculate total arc length for this core segment
    // Arc length = radius * angle
    double totalArcLength = coreRadius * segmentAngle;
    
    if (totalArcLength <= 0.0) {
        return 1.0;  // No surface, fully exposed (degenerate case)
    }
    
    // Calculate total width covered by turns
    double coveredLength = 0.0;
    for (double turnWidth : turnWidths) {
        // Each turn covers its width along the arc
        coveredLength += turnWidth;
    }
    
    // Coverage ratio = exposed / total
    // If turns cover more than the arc (overlap), coverage = 0
    double exposedLength = std::max(0.0, totalArcLength - coveredLength);
    double coverage = exposedLength / totalArcLength;
    
    // Clamp to [0, 1]
    return std::clamp(coverage, 0.0, 1.0);
}

double ThermalNetworkNode::calculateConcentricSurfaceCoverage(double bobbinHeight,
                                                               const std::vector<double>& turnHeights) {
    // For concentric bobbin, the RIGHT face height is covered by turns
    // Coverage ratio = exposed / total
    
    if (bobbinHeight <= 0.0) {
        return 1.0;  // No surface, fully exposed (degenerate case)
    }
    
    // Calculate total height covered by turns
    double coveredHeight = 0.0;
    for (double turnHeight : turnHeights) {
        coveredHeight += turnHeight;
    }
    
    // Coverage ratio = exposed / total
    double exposedHeight = std::max(0.0, bobbinHeight - coveredHeight);
    double coverage = exposedHeight / bobbinHeight;
    
    // Clamp to [0, 1]
    return std::clamp(coverage, 0.0, 1.0);
}

} // namespace OpenMagnetics
