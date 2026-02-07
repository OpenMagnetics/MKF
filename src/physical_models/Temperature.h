#pragma once

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
    bool _isRoundWire = false;      // True for round wires and litz
    std::optional<InsulationWireCoating> _wireCoating;  // Wire coating for thermal calculations
    
    double getMinimumDistanceForConduction() const {
        // Threshold for conduction: accounts for actual gap between surfaces
        // For turns on a toroidal core, the gap is typically wireWidth/2 (half wire extends inward/outward from center)
        // Using wireWidth ensures we capture conduction across reasonable gaps including bobbin/padding
        if (_isRoundWire) {
            return std::max(_wireWidth, _wireHeight) * 0.75;  // 75% of wire diameter
        } else {
            return std::min(_wireWidth, _wireHeight) * 0.75;  // 75% of wire thickness
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
     */
    void createThermalResistances();
    
    /**
     * @brief Create core-to-core conduction resistances
     */
    void createConcentricCoreConnections();
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
     * Connects adjacent turns through appropriate quadrants based on relative position
     */
    void createConcentricTurnToTurnConnections(const std::vector<size_t>& turnNodeIndices, double minConductionDist);
    
    /**
     * @brief Create turn-to-turn connections for toroidal cores
     * 
     * Connects inner/outer turns through radial faces and tangentially adjacent turns
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
};

} // namespace OpenMagnetics
