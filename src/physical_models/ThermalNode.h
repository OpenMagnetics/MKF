#pragma once

/**
 * @file ThermalNode.h
 * @brief Thermal node data structures and utilities for magnetic component thermal modeling
 * 
 * @author OpenMagnetics Team
 * @date 2024-2025
 * 
 * @section node_architecture Thermal Node Architecture
 * 
 * The thermal model represents magnetic components as a network of discrete thermal nodes.
 * Each node has:
 * - Geometric properties (position, dimensions)
 * - Thermal properties (temperature, power dissipation)
 * - Surface subdivision into quadrants for directional heat transfer
 * 
 * @section quadrant_system Quadrant System
 * 
 * Turn nodes are subdivided into quadrants representing different surface directions:
 * 
 * TOROIDAL CORES (polar coordinate system):
 * @code
 *                    TL (Tangential Left)
 *                         ↑
 *                         |
 *    RI (Radial Inner) ←──●──→ RO (Radial Outer)
 *                         |
 *                         ↓
 *                    TR (Tangential Right)
 * @endcode
 * 
 * CONCENTRIC CORES (Cartesian coordinate system):
 * @code
 *                         T (Top, +Y)
 *                         ↑
 *                         |
 *         L (Left, -X) ←──●──→ R (Right, +X)
 *                         |
 *                         ↓
 *                         B (Bottom, -Y)
 * @endcode
 * 
 * Each quadrant stores:
 * - surfaceArea: Physical surface area [m²]
 * - surfaceCoverage: Fraction exposed to air (0-1)
 * - limitCoordinates: Point on the quadrant surface for connection routing
 * - connectionType: Primary heat transfer mechanism
 * 
 * @section resistance_calc Resistance Calculation Methods
 * 
 * 1. CONDUCTION between touching quadrants:
 *    R = (t₁/k₁ + t₂/k₂ + gap/k_air) / A_contact
 *    where t = coating thickness, k = conductivity, gap = air gap between turns
 * 
 * 2. CONVECTION to ambient:
 *    R = 1 / (h_conv × A_exposed)
 *    h_conv from Churchill-Chu correlation for natural convection
 * 
 * 3. RADIATION to ambient:
 *    R = 1 / (h_rad × A_exposed)
 *    h_rad = εσ(T_s² + T_a²)(T_s + T_a)
 * 
 * @section references References
 * 
 * - IEC 60317: Specifications for particular types of winding wires
 * - ASTM D2214: Standard Test Method for Estimating Thermal Conductivity
 * - JEDEC JESD51: Methodology for Thermal Measurement of Component Packages
 */

#include "json.hpp"
#include <magic_enum.hpp>
#include <string>
#include <vector>
#include <optional>
#include <array>
#include <MAS.hpp>

using json = nlohmann::json;
using namespace MAS;  // For TurnCrossSectionalShape, InsulationWireCoating, etc.

namespace OpenMagnetics {

// IMP-9: Bounds-checked coordinate access (inside namespace to be visible to members)
inline double safeCoord(const std::vector<double>& coords, size_t idx, double def = 0.0) {
    return idx < coords.size() ? coords[idx] : def;
}


/**
 * @brief Types of heat transfer mechanisms
 * 
 * Used to classify thermal resistances and specify how heat flows
 * between nodes or from nodes to ambient.
 * 
 * @var CONDUCTION Heat transfer through solid materials
 *      - Dominant in metal components (core, windings)
 *      - Calculated as R = L / (k × A)
 * 
 * @var NATURAL_CONVECTION Heat transfer to surrounding fluid (air) via buoyancy
 *      - Occurs when fluid motion is driven by temperature differences
 *      - Calculated using Churchill-Chu correlation
 * 
 * @var FORCED_CONVECTION Heat transfer with externally-driven fluid flow
 *      - Occurs with fans, blowers, or forced air cooling
 *      - Calculated using flat plate correlations with Reynolds number
 * 
 * @var RADIATION Heat transfer via electromagnetic waves
 *      - Significant at high temperatures (>100°C above ambient)
 *      - Calculated using Stefan-Boltzmann law
 */
enum class HeatTransferType {
    CONDUCTION,
    NATURAL_CONVECTION,
    FORCED_CONVECTION,
    RADIATION,
    HEATSINK_CONVECTION  // Connection through heatsink to ambient
};

/**
 * @brief Surface orientation for convection calculations
 */
enum class SurfaceOrientation {
    VERTICAL,
    HORIZONTAL_TOP,     // Hot surface facing up
    HORIZONTAL_BOTTOM   // Hot surface facing down
};

/**
 * @brief Face direction for quadrant-based thermal nodes in polar coordinates
 * For toroidal cores: each physical object is divided into 4 quadrants
 * representing heat transfer in radial and tangential directions
 */
enum class ThermalNodeFace {
    NONE,                 // Not a quadrant node (used for concentric cores)
    RADIAL_INNER,        // Pointing toward toroidal center (radially inward)
    RADIAL_OUTER,        // Pointing away from toroidal center (radially outward)
    TANGENTIAL_LEFT,     // Tangential direction (left side when going inner->outer)
    TANGENTIAL_RIGHT     // Tangential direction (right side when going inner->outer)
};

/**
 * @brief Geometric shape of a thermal node
 * 
 * Uses MAS::TurnCrossSectionalShape for consistency with turn/wire models:
 * - ROUND: Circular cross-section (round wires, cylindrical cores)
 * - RECTANGULAR: Rectangular cross-section (foil wires, core columns)
 */

/**
 * @brief Types of physical parts that thermal nodes can belong to
 */
enum class ThermalNodePartType {
    CORE_CENTRAL_COLUMN,
    CORE_LATERAL_COLUMN,
    CORE_TOP_YOKE,
    CORE_BOTTOM_YOKE,
    CORE_TOROIDAL_SEGMENT,
    BOBBIN_CENTRAL_COLUMN,
    BOBBIN_TOP_YOKE,
    BOBBIN_BOTTOM_YOKE,
    TURN,
    INSULATION_LAYER,  // Solid insulation layer between turns/sections
    AMBIENT
};

/**
 * @brief Helper functions for MAS::InsulationWireCoating thermal calculations
 * 
 * These functions operate directly on MAS::InsulationWireCoating objects
 * to calculate thermal properties and resistances.
 */
namespace WireCoatingUtils {
    
    /**
     * @brief Calculate thermal resistance of coating layer
     * R_coating = thickness / (k * area)
     * Uses Wire::get_coating_thickness() and Wire::get_coating_thermal_conductivity()
     */
    double calculateCoatingResistance(const InsulationWireCoating& coating, double area);
    
} // namespace WireCoatingUtils

/**
 * @brief Physical properties for a single quadrant of a thermal node
 * 
 * Each node (for toroidal cores) has 4 quadrants representing different
 * faces/surfaces for heat exchange. Each quadrant stores its own
 * physical properties: surface area, length, and material conductivity.
 */
struct ThermalNodeQuadrant {
    ThermalNodeFace face;           // Which quadrant/face (cannot be NONE for toroidal)
    double surfaceArea;             // Exposed surface area for this quadrant (m²)
    double length;                  // Length of wire in this quadrant (m) - half turn length for toroidal
    double thermalConductivity;     // Thermal conductivity of wire material (W/(m·K))
    
    // Surface coverage factor: 1.0 = fully covered/exposed, 0.0 = fully covered by turns
    // For turns: always 1.0 (the turn itself covers the surface)
    // For cores: 0.0 to 1.0 representing exposed surface after turn coverage
    // Used to calculate effective convection area for core surfaces
    double surfaceCoverage = 1.0;
    
    // Coating on this quadrant surface (wire enamel, insulation, etc.)
    // Uses MAS::InsulationWireCoating directly - all thermal properties extracted from it
    std::optional<InsulationWireCoating> coating;
    
    // Limit coordinates: actual surface position for this quadrant face
    // Used for accurate distance calculations in conduction detection
    // For round wires: physicalCoords + radius in quadrant direction
    // For rectangular wires: physicalCoords + width/2 or height/2 in quadrant direction
    std::array<double, 3> limitCoordinates = {0.0, 0.0, 0.0};
    
    // Connection info - which neighboring node/quadrant this connects to
    std::optional<size_t> connectedNodeId;      // Index of connected node in _nodes vector
    std::optional<ThermalNodeFace> connectedQuadrant;  // Which quadrant of connected node
    HeatTransferType connectionType;            // How heat transfers (CONDUCTION, CONVECTION)
    
    ThermalNodeQuadrant()
        : face(ThermalNodeFace::NONE),
          surfaceArea(0.0),
          length(0.0),
          thermalConductivity(385.0),  // Default to copper
          surfaceCoverage(1.0),
          limitCoordinates({0.0, 0.0, 0.0}),
          connectionType(HeatTransferType::CONDUCTION) {}
    
    ThermalNodeQuadrant(ThermalNodeFace f, double area, double len, double k, double coverage = 1.0)
        : face(f),
          surfaceArea(area),
          length(len),
          thermalConductivity(k),
          surfaceCoverage(coverage),
          limitCoordinates({0.0, 0.0, 0.0}),
          connectionType(HeatTransferType::CONDUCTION) {}
    
    /**
     * @brief Calculate conduction resistance to another quadrant
     * R = L / (k * A) where L is average length, k is min conductivity, A is contact area
     * Does NOT include coating resistance - that is added separately
     */
    double calculateConductionResistance(const ThermalNodeQuadrant& other, double contactArea) const;
    
    /**
     * @brief Calculate convection resistance to ambient
     * R = 1 / (h * surfaceArea)
     */
    double calculateConvectionResistance(double heatTransferCoefficient) const;
    
    json toJson() const;
};

/**
 * @brief Geometric dimensions of a thermal node
 * 
 * Self-contained shape information so no need to reference external objects
 */
struct NodeDimensions {
    double width = 0.0;      // X dimension (m) - radial for toroidal
    double height = 0.0;     // Y dimension (m) - axial for toroidal  
    double depth = 0.0;      // Z dimension (m) - into page for 2D view
    double diameter = 0.0;   // For cylindrical shapes (alternative to width/height)
    
    NodeDimensions() = default;
    
    NodeDimensions(double w, double h, double d)
        : width(w), height(h), depth(d), diameter(0.0) {}
    
    // Factory method for cylindrical shape
    static NodeDimensions cylindrical(double diameter, double depth) {
        NodeDimensions dim;
        dim.width = diameter;
        dim.height = diameter;
        dim.depth = depth;
        dim.diameter = diameter;
        return dim;
    }
    
    /**
     * @brief Get cross-sectional area perpendicular to specified axis
     */
    double getCrossSectionalArea(char axis) const {
        switch (axis) {
            case 'x': return height * depth;
            case 'y': return width * depth;
            case 'z': return width * height;
            default: return width * height;
        }
    }
    
    /**
     * @brief Get surface area for specified faces
     */
    double getSurfaceArea(const std::vector<std::string>& faces) const {
        double area = 0.0;
        for (const auto& face : faces) {
            if (face == "top" || face == "bottom") {
                area += width * depth;
            } else if (face == "left" || face == "right") {
                area += height * depth;
            } else if (face == "front" || face == "back") {
                area += width * height;
            }
        }
        return area;
    }
    
    json toJson() const {
        json j;
        j["width"] = width;
        j["height"] = height;
        j["depth"] = depth;
        j["diameter"] = diameter;
        return j;
    }
};

/**
 * @brief Represents a thermal node in the thermal network
 * 
 * A thermal node represents a discrete element in the thermal equivalent circuit.
 * For toroidal cores, each node has 4 quadrants with independent physical properties.
 * For concentric cores, the node has a single set of properties.
 */
class ThermalNetworkNode {
public:
    // Identification
    std::string name;
    ThermalNodePartType part;
    
    // Physical properties
    std::vector<double> physicalCoordinates;  // [x, y, z] in meters - center of node
    double temperature;                        // Temperature in Celsius (uniform across node)
    double powerDissipation;                   // Total power dissipated by this node in Watts
    
    // Geometric properties (self-contained)
    TurnCrossSectionalShape crossSectionalShape = TurnCrossSectionalShape::ROUND;  // From MAS
    NodeDimensions dimensions;                                // Width, height, depth
    double angle = 0.0;                                      // Angular position (for toroidal, radians)
    
    // For toroidal cores: 4 quadrants with independent properties
    // For concentric cores: only index 0 is used (face = NONE)
    std::array<ThermalNodeQuadrant, 4> quadrants;
    
    // Schematic layout
    std::vector<double> schematicCoordinates;  // [x, y] in SVG units
    
    // Optional reference to the actual magnetic component
    std::optional<size_t> windingIndex;
    std::optional<size_t> turnIndex;
    std::optional<size_t> coreSegmentIndex;
    std::optional<size_t> insulationLayerIndex;  // For INSULATION_LAYER nodes
    
    // True if this is an inner half of a toroidal turn
    bool isInnerTurn = false;
    
    // NEW: Fixed temperature flag for cold plate / potting with fixed temperature
    bool isFixedTemperature = false;
    
    // IMP-NEW-10: Thermal capacitance for future transient support
    // C_th = rho * c_p * V [J/K]. Not used by steady-state solver.
    double thermalCapacitance = 0.0;
    
    /**
     * @brief Default constructor
     */
    ThermalNetworkNode() 
        : temperature(25.0), 
          powerDissipation(0.0),
          angle(0.0),
          isInnerTurn(false),
          isFixedTemperature(false) {}
    
    /**
     * @brief Constructor with name and part
     */
    ThermalNetworkNode(const std::string& nodeName, ThermalNodePartType nodePart)
        : name(nodeName),
          part(nodePart),
          temperature(25.0),
          powerDissipation(0.0),
          angle(0.0),
          isInnerTurn(false),
          isFixedTemperature(false) {}
    
    /**
     * @brief Check if this is an ambient node
     */
    bool isAmbient() const {
        return part == ThermalNodePartType::AMBIENT;
    }
    
    /**
     * @brief Check if this is a toroidal node (uses quadrants)
     */
    bool isToroidal() const {
        return quadrants[0].face != ThermalNodeFace::NONE;
    }
    
    /**
     * @brief Get quadrant by face type
     * @return Pointer to quadrant or nullptr if not found
     */
    ThermalNodeQuadrant* getQuadrant(ThermalNodeFace face);
    const ThermalNodeQuadrant* getQuadrant(ThermalNodeFace face) const;
    
    /**
     * @brief Get total surface area across all quadrants
     */
    double getTotalSurfaceArea() const {
        double total = 0.0;
        for (const auto& q : quadrants) {
            total += q.surfaceArea;
        }
        return total;
    }
    
    /**
     * @brief Set geometric properties from wire/turn dimensions
     */
    void setGeometry(double width, double height, double depth, double nodeAngle = 0.0) {
        dimensions = NodeDimensions(width, height, depth);
        angle = nodeAngle;
        crossSectionalShape = TurnCrossSectionalShape::RECTANGULAR;
    }
    
    void setGeometry(double diameter, double depth, double nodeAngle, bool isCylindrical) {
        dimensions = NodeDimensions::cylindrical(diameter, depth);
        angle = nodeAngle;
        crossSectionalShape = TurnCrossSectionalShape::ROUND;
    }
    
    /**
     * @brief Set geometry from MAS turn cross-sectional shape
     */
    void setGeometryFromTurnShape(TurnCrossSectionalShape shape, double dim1, double dim2, double depth, double nodeAngle = 0.0) {
        crossSectionalShape = shape;
        angle = nodeAngle;
        
        if (shape == TurnCrossSectionalShape::ROUND) {
            dimensions = NodeDimensions::cylindrical(dim1, depth);  // dim1 = diameter
        } else if (shape == TurnCrossSectionalShape::RECTANGULAR) {
            dimensions = NodeDimensions(dim1, dim2, depth);  // dim1 = width, dim2 = height
        } else {
            // Default to rectangular
            dimensions = NodeDimensions(dim1, dim2, depth);
        }
    }
    
    /**
     * @brief Initialize quadrants for toroidal turn
     * @param wireWidth Width of rectangular wire (radial dimension)
     * @param wireHeight Height of rectangular wire (axial dimension)
     * @param turnLength Total length of the turn at center radius
     * @param thermalCond Thermal conductivity of wire material
     * @param isInner true if this is inner half of turn
     * @param centerRadius Radius to center of wire (for length calculation)
     * @param wireCoating Optional MAS wire coating object (contains thickness, material, etc.)
     */
    void initializeToroidalQuadrants(double wireWidth, double wireHeight, double turnLength,
                                      double thermalCond, bool isInner, double centerRadius = 0.0,
                                      const std::optional<InsulationWireCoating>& wireCoating = std::nullopt,
                                      TurnCrossSectionalShape shape = TurnCrossSectionalShape::RECTANGULAR);
    
    /**
     * @brief Initialize single quadrant for concentric turn
     * @param wireWidth Width of wire
     * @param wireHeight Height of wire  
     * @param turnLength Total length of the turn
     * @param thermalCond Thermal conductivity of wire material
     * @param wireCoating Optional MAS wire coating object (contains thickness, material, etc.)
     */
    void initializeConcentricQuadrant(double wireWidth, double wireHeight, double turnLength,
                                       double thermalCond,
                                       const std::optional<InsulationWireCoating>& wireCoating = std::nullopt);
    
    /**
     * @brief Initialize quadrants for concentric turn nodes with cardinal directions
     * 
     * For concentric turns arranged in horizontal layers, maps quadrants to
     * cardinal directions (aligned with Cartesian axes):
     * - RADIAL_INNER     (index 0) ↔ LEFT   (-X direction, toward center/bobbin)
     * - RADIAL_OUTER     (index 1) ↔ RIGHT  (+X direction, away from center)
     * - TANGENTIAL_LEFT  (index 2) ↔ TOP    (+Y direction)
     * - TANGENTIAL_RIGHT (index 3) ↔ BOTTOM (-Y direction)
     * 
     * @param wireWidth Wire width in X direction (radial/horizontal)
     * @param wireHeight Wire height in Y direction (axial/vertical)
     * @param turnLength Length of the turn (circumferential)
     * @param thermalCond Thermal conductivity of wire material
     * @param wireCoating Optional wire coating information
     */
    void initializeConcentricTurnQuadrants(double wireWidth, double wireHeight, double turnLength,
                                            double thermalCond,
                                            const std::optional<InsulationWireCoating>& wireCoating = std::nullopt,
                                            TurnCrossSectionalShape shape = TurnCrossSectionalShape::ROUND);
    
    /**
     * @brief Initialize quadrants for concentric core/bobbin nodes with cardinal directions
     * 
     * For concentric cores, maps quadrants to cardinal directions:
     * - TANGENTIAL_LEFT  ↔ TOP    (+Y direction)
     * - TANGENTIAL_RIGHT ↔ BOTTOM (-Y direction)  
     * - RADIAL_INNER     ↔ LEFT   (-X direction, toward center)
     * - RADIAL_OUTER     ↔ RIGHT  (+X direction, away from center)
     * 
     * @param width Width in X direction (horizontal size)
     * @param height Height in Y direction (vertical size)
     * @param depth Depth in Z direction (into the page)
     * @param thermalCond Thermal conductivity of material
     */
    void initializeConcentricCoreQuadrants(double width, double height, double depth, double thermalCond);
    
    /**
     * @brief Initialize quadrants for insulation layer nodes (always rectangular)
     * 
     * Insulation layers are rectangular nodes that sit between turns or sections.
     * They have 4 quadrants representing their 4 faces, allowing heat conduction
     * from any adjacent turn.
     * 
     * Quadrant mapping for insulation layers:
     * - Index 0 (RADIAL_OUTER): RIGHT face (+X direction) - facing outer turns
     * - Index 1 (RADIAL_INNER): LEFT face (-X direction) - facing inner turns/core
     * - Index 2 (TANGENTIAL_LEFT): TOP face (+Y direction) - facing top turns/yoke
     * - Index 3 (TANGENTIAL_RIGHT): BOTTOM face (-Y direction) - facing bottom turns/yoke
     * 
     * @param width Layer width in X direction (m) - typically thickness of insulation
     * @param height Layer height in Y direction (m) - typically span of winding window
     * @param depth Layer depth in Z direction (m) - typically depth of coil
     * @param thermalCond Thermal conductivity of insulation material (W/m·K)
     */
    void initializeInsulationLayerQuadrants(double width, double height, double depth, double thermalCond);
    
    /**
     * @brief Calculate surface coverage for toroidal core quadrants
     * 
     * For toroidal cores, calculates what fraction of the inner/outer arc surface
     * is covered by turns vs exposed to air. This is used to determine effective
     * convection area for core cooling.
     * 
     * @param coreRadius Inner or outer radius of the toroidal core (m)
     * @param segmentAngle Angular span of this core segment (radians)
     * @param turnWidths Vector of widths (radial dimension) of turns touching this surface
     * @return Coverage ratio 0.0 to 1.0 (1.0 = fully exposed, 0.0 = fully covered)
     */
    static double calculateToroidalSurfaceCoverage(double coreRadius, double segmentAngle,
                                                    const std::vector<double>& turnWidths);
    
    /**
     * @brief Calculate surface coverage for concentric bobbin nodes
     * 
     * For concentric bobbins, calculates what fraction of the RIGHT face (facing turns)
     * is covered by turns vs exposed to air. The bobbin height is divided among turns,
     * so coverage depends on how many turns touch the bobbin surface.
     * 
     * @param bobbinHeight Height of the bobbin central column (m)
     * @param turnHeights Vector of heights (Y dimension) of turns touching this surface
     * @return Coverage ratio 0.0 to 1.0 (1.0 = fully exposed, 0.0 = fully covered)
     */
    static double calculateConcentricSurfaceCoverage(double bobbinHeight,
                                                      const std::vector<double>& turnHeights);
    
    /**
     * @brief Convert to JSON for serialization
     */
    json toJson() const {
        json j;
        j["name"] = name;
        j["part"] = magic_enum::enum_name(part);
        j["physicalCoordinates"] = physicalCoordinates;
        j["schematicCoordinates"] = schematicCoordinates;
        j["temperature"] = temperature;
        j["powerDissipation"] = powerDissipation;
        j["isInnerTurn"] = isInnerTurn;
        j["crossSectionalShape"] = magic_enum::enum_name(crossSectionalShape);
        j["dimensions"] = dimensions.toJson();
        j["angle"] = angle;
        
        json quadrantsJson = json::array();
        for (const auto& q : quadrants) {
            quadrantsJson.push_back(q.toJson());
        }
        j["quadrants"] = quadrantsJson;
        
        if (windingIndex.has_value()) j["windingIndex"] = windingIndex.value();
        if (turnIndex.has_value()) j["turnIndex"] = turnIndex.value();
        if (coreSegmentIndex.has_value()) j["coreSegmentIndex"] = coreSegmentIndex.value();
        
        return j;
    }
};

/**
 * @brief Electrical/thermal insulation layer between components
 * 
 * Represents insulation tape, sheets, or barriers that are NOT part of
 * the component itself but placed between components for electrical isolation.
 * This adds series thermal resistance to the connection.
 */
struct InsulationLayer {
    double thickness = 0.0;              // Insulation thickness in meters
    double thermalConductivity = 0.2;     // Insulation thermal conductivity W/(m·K)
    std::string materialName;            // Material identifier (e.g., "kapton", "mylar", "nomex")
    std::string description;             // Human-readable description (e.g., "layer tape between turns")
    
    InsulationLayer() = default;
    
    InsulationLayer(double t, double k, const std::string& name, const std::string& desc = "")
        : thickness(t), thermalConductivity(k), materialName(name), description(desc) {}
    
    /**
     * @brief Calculate thermal resistance of this insulation layer
     * R_insulation = thickness / (k * area)
     */
    double calculateInsulationResistance(double area) const {
        if (thickness <= 0.0 || thermalConductivity <= 0.0 || area <= 0.0) {
            return 0.0;  // No insulation resistance
        }
        return thickness / (thermalConductivity * area);
    }
    
    json toJson() const {
        json j;
        j["thickness"] = thickness;
        j["thermalConductivity"] = thermalConductivity;
        j["materialName"] = materialName;
        j["description"] = description;
        return j;
    }
};

/**
 * @brief Thermal resistance element connecting two node quadrants
 * 
 * This represents a thermal connection between specific quadrants of two nodes.
 * The resistance is calculated from the quadrant properties.
 * Can include additional insulation layers between the components.
 */
struct ThermalResistanceElement {
    size_t nodeFromId;              // Index of first node in _nodes vector
    ThermalNodeFace quadrantFrom;   // Quadrant of nodeFrom that this resistor connects to
    size_t nodeToId;                // Index of second node in _nodes vector
    ThermalNodeFace quadrantTo;     // Quadrant of nodeTo that this resistor connects to
    HeatTransferType type;          // CONDUCTION, NATURAL_CONVECTION, FORCED_CONVECTION, RADIATION
    double resistance = 10.0;       // Thermal resistance in K/W (calculated)
    
    // Electrical/thermal insulation layers between components (NOT coating on the components themselves)
    // These are added in series with the base resistance
    std::vector<InsulationLayer> insulationLayers;
    
    // Backward compatibility fields
    double area = 0.0;              // Contact/exposed area in m²
    double length = 0.0;            // Conduction path length in m
    double thermalConductivity = 0.0; // Material thermal conductivity W/(m·K)
    SurfaceOrientation orientation = SurfaceOrientation::VERTICAL;
    
    // Backward compatibility: transferType maps to type
    HeatTransferType getTransferType() const { return type; }
    void setTransferType(HeatTransferType t) { type = t; }
    
    ThermalResistanceElement()
        : nodeFromId(0),
          quadrantFrom(ThermalNodeFace::NONE),
          nodeToId(0),
          quadrantTo(ThermalNodeFace::NONE),
          type(HeatTransferType::CONDUCTION),
          resistance(10.0),
          area(0.0),
          length(0.0),
          thermalConductivity(0.0),
          orientation(SurfaceOrientation::VERTICAL) {}
    
    // Backward compatibility: access type as transferType
    HeatTransferType transferType() const { return type; }
    void transferType(HeatTransferType t) { type = t; }
    
    /**
     * @brief Add an insulation layer to this connection
     */
    void addInsulationLayer(double thickness, double conductivity, 
                           const std::string& material, const std::string& desc = "") {
        insulationLayers.emplace_back(thickness, conductivity, material, desc);
    }
    
    /**
     * @brief Calculate total insulation resistance from all layers
     * R_total = sum(R_layer_i) for all layers
     */
    double calculateTotalInsulationResistance(double area) const {
        double totalR = 0.0;
        for (const auto& layer : insulationLayers) {
            totalR += layer.calculateInsulationResistance(area);
        }
        return totalR;
    }
    
    /**
     * @brief Get total resistance including insulation layers
     * R_total = R_base + R_insulation_layers
     */
    double getTotalResistance(double area) const {
        return resistance + calculateTotalInsulationResistance(area);
    }
    
    json toJson() const {
        json j;
        j["nodeFromId"] = nodeFromId;
        j["quadrantFrom"] = magic_enum::enum_name(quadrantFrom);
        j["nodeToId"] = nodeToId;
        j["quadrantTo"] = magic_enum::enum_name(quadrantTo);
        j["type"] = magic_enum::enum_name(type);
        j["resistance"] = resistance;
        j["area"] = area;
        j["length"] = length;
        j["thermalConductivity"] = thermalConductivity;
        
        json layersJson = json::array();
        for (const auto& layer : insulationLayers) {
            layersJson.push_back(layer.toJson());
        }
        j["insulationLayers"] = layersJson;
        
        return j;
    }
};

} // namespace OpenMagnetics
