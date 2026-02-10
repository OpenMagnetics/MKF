#pragma once

/**
 * @file Temperature.h
 * @brief Thermal analysis system for magnetic components
 * 
 * @author OpenMagnetics Team
 * @date 2024-2025
 * 
 * @section overview Overview
 * 
 * This module provides comprehensive thermal analysis for magnetic components
 * (inductors, transformers) by modeling heat generation, conduction, convection,
 * and radiation. It creates a thermal equivalent circuit and solves for
 * temperature distribution.
 * 
 * @section theory Physical Theory
 * 
 * The thermal model is based on the thermal-electrical analogy:
 * - Temperature (T) ↔ Voltage (V)
 * - Heat flow (Q) ↔ Current (I)
 * - Thermal resistance (Rth) ↔ Electrical resistance (R)
 * 
 * Key equations:
 * - Conduction: R_th = L / (k × A)  [K/W]
 *   where L = length, k = thermal conductivity, A = cross-sectional area
 * 
 * - Convection: R_th = 1 / (h × A)  [K/W]
 *   where h = convection coefficient [W/(m²·K)]
 * 
 * - Natural convection coefficient (Churchill & Chu correlation):
 *   Nu_L = {0.825 + 0.387 Ra_L^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/9)}²
 *   h = Nu_L × k_air / L
 * 
 * - Radiation: R_th = 1 / (h_rad × A)
 *   where h_rad = ε × σ × (T_s² + T_a²)(T_s + T_a)  [Stefan-Boltzmann law]
 * 
 * @section architecture Architecture
 * 
 * The system uses a node-based thermal network:
 * 
 * 1. **ThermalNetworkNode**: Represents discrete temperature points
 *    - Core nodes: Central column, yokes, toroidal segments
 *    - Turn nodes: Individual windings with quadrant subdivision
 *    - Ambient node: Reference temperature boundary
 * 
 * 2. **ThermalNodeQuadrant**: Subdivision of turn surfaces for directional heat transfer
 *    - TOROIDAL: RI (Radial Inner), RO (Radial Outer), TL (Tangential Left), TR (Tangential Right)
 *    - CONCENTRIC: L (Left), R (Right), T (Top), B (Bottom)
 * 
 * 3. **ThermalResistanceElement**: Connections between node quadrants
 *    - Type: CONDUCTION, NATURAL_CONVECTION, FORCED_CONVECTION, RADIATION
 *    - Calculated based on geometry and material properties
 * 
 * @section algorithm Connection Algorithm
 * 
 * Turn-to-Turn Connections (Purely Geometric):
 * - Condition: surface_distance < min(wire_width, wire_height) / 4
 * - Surface distance = center_distance - extent1 - extent2
 * - Extent = min(wire_dimensions) / 2 for each wire
 * - Facing quadrants determined by dot product of direction vectors
 * 
 * Convection to Air (Purely Geometric):
 * - Blocked if object exists in quadrant direction within wire_max_dim distance
 * - Minimum radial difference > wire_min_dim/4 to avoid tangential neighbor blocking
 * - Angular tolerance: 0.3 rad (~17°)
 * 
 * @section references References
 * 
 * 1. Churchill & Chu (1975) - Correlating equations for laminar and turbulent free
 *    convection from a vertical plate. Int. J. Heat Mass Transfer.
 * 
 * 2. Incropera & DeWitt (2002) - Fundamentals of Heat and Mass Transfer, 6th Ed.
 *    Wiley. Chapter 7 (External Convection).
 * 
 * 3. Kyaw et al. (2018) - Thermal modeling of litz wire windings. IEEE Trans. Magn.
 * 
 * 4. López et al. (2019) - Anisotropic thermal conductivity in magnetic component
 *    windings. IEEE Trans. Power Electron.
 * 
 * 5. Mukosiej et al. (2022) - Thermal modeling of high-frequency magnetic components
 *    using homogenization techniques. IEEE Trans. Power Electron.
 * 
 * @section ferrite_materials Ferrite Thermal Properties
 * 
 * Material         Thermal Conductivity (W/m·K)  Source
 * ------------------------------------------------------------------------
 * 3C90 (Ferroxcube)        4.0-4.5            Ferroxcube Datasheet
 * 3C95 (Ferroxcube)        4.0-4.5            Ferroxcube Datasheet
 * 3C97 (Ferroxcube)        4.5-5.0            Ferroxcube Datasheet
 * N87 (TDK)                4.0-4.5            TDK Datasheet
 * N97 (TDK)                4.0-4.5            TDK Datasheet
 * PC40 (TDK)               4.0                TDK Datasheet
 * PC95 (TDK)               4.0                TDK Datasheet
 * 
 * @section copper_wire Copper Wire Properties
 * 
 * Property                    Value                      Source
 * ------------------------------------------------------------------------
 * Thermal Conductivity        385-400 W/(m·K)           CRC Handbook
 * Electrical Resistivity      1.68×10⁻⁸ Ω·m            @ 20°C
 * Temperature Coefficient     0.00393 /°C              α (resistance)
 * Enamel Insulation Thickness 20-50 μm                 IEC 60317
 * Enamel k                    0.2-0.4 W/(m·K)          Polyimide
 * 
 * @section air_properties Air Properties (@ 25°C, 1 atm)
 * 
 * Property                    Value
 * ------------------------------------------------------------------------
 * Density                     1.184 kg/m³
 * Thermal Conductivity        0.0261 W/(m·K)
 * Kinematic Viscosity         15.7×10⁻⁶ m²/s
 * Prandtl Number              0.71
 * Thermal Expansion Coeff.    3.35×10⁻³ /K
 */

#include "ThermalNode.h"
#include "ThermalResistance.h"
#include "Magnetic.h"
#include "json.hpp"
#include <MAS.hpp>
#include <vector>
#include <memory>
#include <map>

using json = nlohmann::json;
using namespace MAS;

namespace OpenMagnetics {

// Old thermal types for compatibility with BasicPainter
enum class ThermalNodeType {
    CORE_CENTRAL_COLUMN,
    CORE_LATERAL_COLUMN,
    CORE_TOP_YOKE,
    CORE_BOTTOM_YOKE,
    CORE_TOROIDAL_SEGMENT,
    COIL_TURN,
    COIL_TURN_INNER,
    COIL_TURN_OUTER,
    BOBBIN_INNER,
    BOBBIN_OUTER,
    BOBBIN_COLUMN,
    BOBBIN_TOP_WALL,
    BOBBIN_BOTTOM_WALL,
    AMBIENT
};

/**
 * @brief Cooling types for thermal analysis
 */
enum class CoolingType {
    NATURAL_CONVECTION,   // Default when no cooling specified
    FORCED_CONVECTION,    // velocity is set
    HEATSINK,             // thermal_resistance is set, no maximum_temperature
    COLD_PLATE,           // maximum_temperature is set
    UNKNOWN
};

/**
 * @brief Utility class for cooling configuration
 */
class CoolingUtils {
public:
    /**
     * @brief Detect cooling type from MAS::Cooling object
     */
    static CoolingType detectCoolingType(const MAS::Cooling& cooling);
    
    /**
     * @brief Check if cooling is natural convection
     */
    static bool isNaturalConvection(const MAS::Cooling& cooling);
    
    /**
     * @brief Check if cooling is forced convection
     */
    static bool isForcedConvection(const MAS::Cooling& cooling);
    
    /**
     * @brief Check if cooling is heatsink
     */
    static bool isHeatsink(const MAS::Cooling& cooling);
    
    /**
     * @brief Check if cooling is cold plate
     */
    static bool isColdPlate(const MAS::Cooling& cooling);
    
    /**
     * @brief Calculate forced convection coefficient
     */
    static double calculateForcedConvectionCoefficient(
        double surfaceTemp,
        double ambientTemp,
        double velocity,
        double characteristicLength,
        double fluidConductivity = 0.0261,  // W/m·K for air at 25°C
        double kinematicViscosity = 15.7e-6,  // m²/s for air at 25°C
        double prandtlNumber = 0.71);
    
    /**
     * @brief Calculate mixed convection coefficient (natural + forced)
     */
    static double calculateMixedConvectionCoefficient(
        double h_natural,
        double h_forced);
};

struct ThermalNode {
    size_t id;
    ThermalNodeType type;
    std::string name;
    double temperature;
    double powerDissipation;
    std::vector<double> coordinates;
    double volume = 0.0;
    double exposedSurfaceArea = 0.0;
    double emissivity = 0.9;
    
    bool isAmbient() const { return type == ThermalNodeType::AMBIENT; }
};

/**
 * @brief Output structure for thermal analysis results
 */
struct ThermalResult {
    std::string methodUsed;
    double maximumTemperature;
    double averageCoreTemperature;
    double averageCoilTemperature;
    std::map<std::string, double> nodeTemperatures;
    double totalThermalResistance;
    size_t iterationsToConverge;
    bool converged;
    std::vector<ThermalResistanceElement> thermalResistances;
};

/**
 * @brief Configuration for thermal analysis
 * 
 * This structure controls all aspects of the thermal simulation, including:
 * - Model fidelity (per-turn vs per-layer vs single coil node)
 * - Physics inclusion (convection types, radiation)
 * - Material properties (thermal conductivities)
 * - Insulation modeling (inter-turn, turn-to-core)
 * - Solver settings (convergence tolerance, max iterations)
 * - Output options (schematic generation)
 * 
 * @note Default values are set for typical ferrite core + copper winding
 *       applications at natural convection conditions.
 */
struct TemperatureConfig {
    // Core settings
    int toroidalSegments = 8;           // Number of segments for toroidal cores
    double ambientTemperature = 25.0;   // Ambient temperature in Celsius
    
    // Loss inputs
    double coreLosses = 0.0;            // Total core losses in W
    double windingLosses = 0.0;         // Total winding losses in W (if no per-turn distribution)
    std::optional<WindingLossesOutput> windingLossesOutput;  // Per-turn losses from MKF simulation (preferred)
    
    // Thermal model configuration
    double convergenceTolerance = 0.1;  // Temperature convergence criterion (°C)
    size_t maxIterations = 100;         // Maximum solver iterations
    double coreThermalConductivity = 4.0;  // Ferrite thermal conductivity (W/m·K)
    
    // Inter-turn insulation (electrical insulation tape between turns)
    bool useInterTurnInsulation = false;        // Enable inter-turn insulation layers
    double interTurnInsulationThickness = 0.0;  // Thickness of insulation between turns (m)
    double interTurnInsulationConductivity = 0.2; // Thermal conductivity of insulation (W/m·K)
    
    // Turn-to-core insulation
    bool useTurnToCoreInsulation = false;       // Enable insulation between turns and core
    double turnToCoreInsulationThickness = 0.0; // Thickness of insulation between turn and core (m)
    double turnToCoreInsulationConductivity = 0.2; // Thermal conductivity (W/m·K)
    
    // Convection settings
    bool includeForcedConvection = false;
    double airVelocity = 0.0;           // Air velocity for forced convection (m/s)
    bool includeRadiation = true;
    double surfaceEmissivity = 0.9;
    
    // Model resolution
    bool nodePerCoilLayer = false;      // Create separate nodes for coil layers
    bool nodePerCoilTurn = false;       // Create separate nodes for each turn (per-turn model)
    
    // Output settings
    bool plotSchematic = true;          // Whether to generate schematic visualization
    std::string schematicOutputPath = "output/thermal_schematic.svg";
    
    // ===== NEW: MAS Cooling Configuration =====
    // Optional MAS cooling object (from OperatingConditions)
    std::optional<MAS::Cooling> masCooling;
    
    // Factory method to create config from MAS inputs
    static TemperatureConfig fromMasOperatingConditions(
        const MAS::OperatingConditions& conditions);
};

/**
 * @brief Main class for thermal analysis of magnetic components
 * 
 * This class acts as a wrapper that:
 * 1. Analyzes a Magnetic component and breaks it into thermal nodes
 * 2. Determines thermal connections between nodes (conduction/convection)
 * 3. Calculates thermal resistances
 * 4. Creates and solves the thermal equivalent circuit
 * 5. Generates schematic visualization
 * 
 * Architecture:
 * - For toroidal cores: Each turn has 2 nodes (inner/outer), each with 4 quadrants
 * - For concentric cores: Each turn has 1 node with a single quadrant
 * - Core nodes: Fixed number of segments for toroidal, 4 nodes for concentric
 * - Thermal resistances connect specific quadrants of specific nodes
 */
class Temperature {
private:
    Magnetic _magnetic;
    TemperatureConfig _config;
    
    std::vector<ThermalNetworkNode> _nodes;
    std::vector<ThermalResistanceElement> _resistances;
    
    double _scaleFactor = 1.0;  // Scale factor for schematic coordinates
    
    // Wire properties cache
    double _wireWidth = 0.001;      // Wire width (radial dimension) in meters
    double _wireHeight = 0.001;     // Wire height (axial dimension) in meters
    double _wireThermalCond = 385.0; // Copper thermal conductivity
    bool _isToroidal = false;
    bool _isPlanar = false;         // True for planar windings (PCB traces)
    bool _isRoundWire = false;      // True for round wires and litz
    std::optional<InsulationWireCoating> _wireCoating;  // Wire coating for thermal calculations
    
    double getMinimumDistanceForConduction() const {
        // Threshold for conduction: accounts for actual gap between surfaces
        // For turns on a toroidal core, the gap is typically wireWidth/2 (half wire extends inward/outward from center)
        // Using wireWidth ensures we capture conduction across reasonable gaps including bobbin/padding
        if (_isRoundWire) {
            return std::max(_wireWidth, _wireHeight) * 0.25;  // 25% of wire diameter
        } else {
            return std::min(_wireWidth, _wireHeight) * 0.25;  // 25% of wire thickness
        }
    }
    
    double getMaximumDistanceForConvection() const {
        if (_isRoundWire) {
            return std::max(_wireWidth, _wireHeight);  // wireDiameter
        } else {
            return std::min(_wireWidth, _wireHeight);  // min(width,height)
        }
    }
    
public:
    /**
     * @brief Constructor with magnetic component only
     */
    explicit Temperature(const Magnetic& magnetic, 
                        const TemperatureConfig& config = TemperatureConfig())
        : _magnetic(magnetic), _config(config) {}
    
    /**
     * @brief Main method to calculate temperatures
     * 
     * Process flow:
     * 1. Create core nodes (toroidal segments or concentric structure)
     * 2. Create turn nodes with quadrant initialization
     * 3. Calculate distances between node surfaces
     * 4. Create thermal resistances between connecting quadrants
     * 5. Solve thermal circuit
     * 
     * @return Temperature calculation result
     */
    ThermalResult calculateTemperatures();
    
    /**
     * @brief Get the list of thermal nodes
     */
    const std::vector<ThermalNetworkNode>& getNodes() const { return _nodes; }
    
    /**
     * @brief Get the list of thermal resistances
     */
    const std::vector<ThermalResistanceElement>& getResistances() const { return _resistances; }
    
    /**
     * @brief Get configuration
     */
    const TemperatureConfig& getConfig() const { return _config; }
    
    /**
     * @brief Get temperature at a specific 3D point using nearest node interpolation
     * @param point 3D coordinates [x, y, z] in meters
     * @return Temperature at that point in °C
     */
    double getTemperatureAtPoint(const std::vector<double>& point) const;
    
    /**
     * @brief Set configuration
     */
    void setConfig(const TemperatureConfig& config) { _config = config; }
    
    /**
     * @brief Get bulk thermal resistance (Tmax - Tambient) / Ptotal
     */
    double getBulkThermalResistance() const;
    
    // Legacy API functions (maintained for backward compatibility)
    static double calculate_temperature_from_core_thermal_resistance(Core core, double totalLosses);
    static double calculate_temperature_from_core_thermal_resistance(double thermalResistance, double totalLosses);

private:
    // =========================================================================
    // Node Creation
    // =========================================================================
    
    /**
     * @brief Extract wire properties from coil and cache them
     */
    void extractWireProperties();
    
    /**
     * @brief Extract wire from first winding in coil
     * @return Optional wire object (null if no wire found)
     */
    std::optional<Wire> extractWire() const;
    
    /**
     * @brief Create all thermal nodes (core and turns)
     */
    void createThermalNodes();
    
    /**
     * @brief Create core nodes for concentric core types (E, ETD, PQ, etc.)
     * TODO: Implement proper quadrant support for concentric cores
     */
    void createConcentricCoreNodes();
    
    /**
     * @brief Create core nodes for toroidal cores (segmented ring)
     */
    void createToroidalCoreNodes();
    
    /**
     * @brief Create bobbin nodes if present
     */
    void createBobbinNodes();
    
    /**
     * @brief Create turn nodes with quadrant initialization
     * 
     * For toroidal: 2 nodes per turn (inner/outer), each with 4 quadrants
     * For concentric: 1 node per turn, single quadrant
     */
    void createTurnNodes();
    
    // =========================================================================
    // Resistance Creation
    // =========================================================================
    
    /**
     * @brief Create all thermal resistances between node quadrants
     * 
     * Checks distances between node surfaces and creates resistances
     * between connecting quadrants.
     * 
     * Connection types created:
     * - Core-to-core conduction (through core material)
     * - Turn-to-turn conduction (when surfaces are close)
     * - Turn-to-bobbin conduction (when in contact)
     * - Convection/radiation to ambient (exposed surfaces)
     */
    void createThermalResistances();
    
    /**
     * @brief Create core-to-core conduction resistances
     * 
     * For concentric cores: connects central column to yokes, yokes to lateral columns
     * Resistance formula: R = L / (k_core × A)
     * where L is path length, A is cross-sectional area
     */
    void createConcentricCoreConnections();
    
    /**
     * @brief Create core-to-core conduction resistances for toroidal cores
     * 
     * Connects adjacent segments in the toroidal ring.
     * Also creates connection from innermost segments to ambient (inner window).
     */
    void createToroidalCoreConnections();
    
    /**
     * @brief Create bobbin-related conduction resistances
     */
    void createBobbinConnections();
    
    /**
     * @brief Create connections from bobbin yokes to nearby turns
     * @param bobbinTopYokeIdx Index of top bobbin yoke node
     * @param bobbinBottomYokeIdx Index of bottom bobbin yoke node  
     * @param turnNodeIndices Vector of turn node indices
     */
    void createBobbinYokeToTurnConnections(size_t bobbinTopYokeIdx, size_t bobbinBottomYokeIdx,
                                           const std::vector<size_t>& turnNodeIndices);
    
    /**
     * @brief Create turn-to-turn conduction resistances based on distance
     * 
     * For toroidal: connects quadrants of adjacent turns if surfaces are close
     * For concentric: standard distance-based connections
     */
    void createTurnToTurnConnections();
    
    /**
     * @brief Create turn-to-turn connections for concentric cores
     * 
     * Connects adjacent turns through appropriate quadrants based on relative position.
     * 
     * Connection criteria:
     * - surface_distance < min(wire_dims) / 4
     * - surface_distance = center_distance - extent1 - extent2
     * - extent = min(wire_width, wire_height) / 2
     * 
     * @param turnNodeIndices List of turn node indices to connect
     * @param minConductionDist Minimum distance for conduction threshold
     */
    void createConcentricTurnToTurnConnections(const std::vector<size_t>& turnNodeIndices, double minConductionDist);
    
    /**
     * @brief Create turn-to-turn connections for toroidal cores
     * 
     * Connects turns through radial and tangential faces based on geometric proximity.
     * 
     * Connection Logic:
     * 1. Calculate center distance between all turn pairs
     * 2. Calculate surface distance using wire extents
     * 3. If surface_distance < threshold, find facing quadrants using dot product
     * 4. Create conduction resistance between facing quadrants
     * 
     * Facing quadrant determination:
     * - Calculate turn angles using atan2(y, x)
     * - Quadrant directions derived from turn angle
     * - Best facing pair found via maximum dot product with connection direction
     * 
     * @param turnNodeIndices List of turn node indices to connect
     * @param minConductionDist Minimum distance for conduction threshold
     */
    void createToroidalTurnToTurnConnections(const std::vector<size_t>& turnNodeIndices, double minConductionDist);
    
    /**
     * @brief Create turn-to-bobbin conduction resistances for concentric cores
     * 
     * Checks turn quadrant limit coordinates against bobbin surfaces
     * Creates connections when turns are close to bobbin walls
     */
    void createTurnToBobbinConnections();
    
    /**
     * @brief Create turn-to-core conduction resistances based on distance
     * 
     * For toroidal: inner turn quadrants connect to inner core, outer to outer
     * For concentric: turns connect to core if no bobbin present
     */
    void createTurnToSolidConnections();
    
    /**
     * @brief Create thermal nodes for insulation layers
     * 
     * Creates nodes for solid insulation layers between turns or sections.
     * Each insulation layer becomes a thermal node with 4 quadrants (faces)
     * that can conduct heat to/from adjacent turns.
     * 
     * Insulation layers have no heat generation but provide conduction paths
     * between turns that are separated by insulation material.
     */
    void createInsulationLayerNodes();
    
    /**
     * @brief Create connections between turns and insulation layers
     * 
     * Geometrically determines which turns are adjacent to each insulation
     * layer quadrant and creates conduction resistances between them.
     * 
     * Connections are based on spatial proximity - if a turn is close enough
     * to an insulation layer face, a thermal connection is created.
     */
    void createTurnToInsulationConnections();
    
    /**
     * @brief Create convection resistances to ambient
     * 
     * Creates convection connections for exposed quadrants
     */
    void createConvectionConnections();
    
    // =========================================================================
    // Helper Methods
    // =========================================================================
    
    /**
     * @brief Calculate surface-to-surface distance between two nodes
     * @return Distance in meters, negative if overlapping
     */
    double calculateSurfaceDistance(const ThermalNetworkNode& node1, 
                                    const ThermalNetworkNode& node2) const;
    
    /**
     * @brief Check if two turn quadrants should be connected
     * @return true if surfaces are within contact threshold
     */
    bool shouldConnectQuadrants(const ThermalNetworkNode& node1, ThermalNodeFace face1,
                                const ThermalNetworkNode& node2, ThermalNodeFace face2) const;
    
    /**
     * @brief Get insulation layer thermal resistance between turns
     */
    double getInsulationLayerThermalResistance(int turnIdx1, int turnIdx2, double contactArea);
    
    /**
     * @brief Calculate contact area between two quadrants
     */
    double calculateContactArea(const ThermalNodeQuadrant& q1, const ThermalNodeQuadrant& q2) const;
    
    void calculateSchematicScaling();
    void plotSchematic();
    ThermalResult solveThermalCircuit();
    
    bool hasBobbinNodes() const;
    
    // =========================================================================
    // Cooling Application Methods
    // =========================================================================
    
    /**
     * @brief Apply MAS cooling configuration
     */
    void applyMasCooling(const MAS::Cooling& cooling);
    
    /**
     * @brief Apply forced convection cooling
     */
    void applyForcedConvection(const MAS::Cooling& cooling);
    
    /**
     * @brief Apply heatsink cooling (top mount for concentric cores)
     */
    void applyHeatsinkCooling(const MAS::Cooling& cooling);
    
    /**
     * @brief Apply cold plate cooling (bottom mount)
     */
    void applyColdPlateCooling(const MAS::Cooling& cooling);
};

} // namespace OpenMagnetics
