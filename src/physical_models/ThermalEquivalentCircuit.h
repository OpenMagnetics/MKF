#pragma once

#include "Constants.h"
#include "Defaults.h"
#include "constructive_models/Magnetic.h"
#include <MAS.hpp>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <Eigen/Dense>

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief Types of thermal nodes in the equivalent circuit
 */
enum class ThermalNodeType {
    CORE_CENTRAL_COLUMN,
    CORE_LATERAL_COLUMN,
    CORE_TOP_YOKE,
    CORE_BOTTOM_YOKE,
    COIL_SECTION,
    COIL_LAYER,
    COIL_TURN,
    BOBBIN_INNER,
    BOBBIN_OUTER,
    AMBIENT
};

/**
 * @brief Types of heat transfer mechanisms
 */
enum class HeatTransferType {
    CONDUCTION,
    CONVECTION_NATURAL,
    CONVECTION_FORCED,
    RADIATION
};

/**
 * @brief Surface orientation for convection calculations
 */
enum class SurfaceOrientation {
    HORIZONTAL_TOP,
    HORIZONTAL_BOTTOM,
    VERTICAL
};

/**
 * @brief Fluid properties for convection calculations
 */
struct FluidProperties {
    double density;                    // kg/m³
    double thermalConductivity;        // W/(m·K)
    double dynamicViscosity;           // Pa·s
    double specificHeat;               // J/(kg·K)
    double thermalExpansionCoefficient; // 1/K
    double kinematicViscosity;         // m²/s
    double prandtlNumber;              // dimensionless
    
    /**
     * @brief Get air properties at a given temperature
     * @param temperature Temperature in Celsius
     * @return FluidProperties for air at that temperature
     */
    static FluidProperties getAirProperties(double temperature);
};

/**
 * @brief Represents a node in the thermal equivalent circuit
 */
struct ThermalNode {
    size_t id;
    ThermalNodeType type;
    std::string name;
    double temperature;           // Current temperature (°C)
    double powerDissipation;      // Heat generated at node (W)
    std::vector<double> coordinates; // 3D coordinates (m)
    double volume;                // Volume of the element (m³)
    double exposedSurfaceArea;    // Surface exposed to fluid (m²)
    double emissivity;            // Surface emissivity for radiation
    
    bool isAmbient() const { return type == ThermalNodeType::AMBIENT; }
};

/**
 * @brief Represents a thermal resistance between two nodes
 */
struct ThermalResistanceElement {
    size_t nodeFromId;
    size_t nodeToId;
    HeatTransferType transferType;
    double resistance;            // K/W
    double area;                  // Cross-sectional or surface area (m²)
    double length;                // Conduction length (m)
    double thermalConductivity;   // Material thermal conductivity (W/(m·K))
    SurfaceOrientation orientation; // For convection calculations
};

/**
 * @brief Configuration for the thermal model
 */
struct ThermalModelConfiguration {
    double ambientTemperature = 25.0;           // °C
    double convergenceTolerance = 0.1;          // °C
    size_t maxIterations = 100;
    bool includeForcedConvection = false;
    double airVelocity = 0.0;                   // m/s (for forced convection)
    double coreThermalConductivity = 4.0;       // W/(m·K) - typical ferrite
    double bobbinThermalConductivity = 0.2;     // W/(m·K) - typical plastic
    bool includeRadiation = true;
    double defaultEmissivity = 0.9;             // For painted/dark surfaces
    
    // Granularity options
    bool nodePerCoreColumn = true;
    bool nodePerCoilLayer = true;               // If false, one node per section
    bool nodePerCoilTurn = false;               // Maximum granularity (expensive)
};

/**
 * @brief Output structure for thermal analysis results
 */
struct ThermalAnalysisOutput {
    std::string methodUsed;
    double maximumTemperature;
    double averageCoreTemperature;
    double averageCoilTemperature;
    std::map<std::string, double> nodeTemperatures;
    double totalThermalResistance;              // Junction-to-ambient
    size_t iterationsToConverge;
    bool converged;
    
    // Per-element thermal resistances
    std::vector<ThermalResistanceElement> thermalResistances;
};

/**
 * @brief Main class for thermal equivalent circuit modeling
 * 
 * This class creates and solves a thermal equivalent circuit for magnetic components.
 * It models:
 * - Conduction through core, bobbin, and windings
 * - Convection (natural or forced) from exposed surfaces
 * - Radiation from exposed surfaces
 * 
 * The circuit is solved iteratively since convection/radiation coefficients
 * depend on temperature.
 * 
 * Based on:
 * - Van den Bossche & Valchev: "Thermal Modeling of E-type Magnetic Components"
 * - Dey et al.: "Lumped Parameter Thermal Network Modelling of Power Transformers"
 * - Salinas: PhD Thesis on Thermal Modelling of High-Frequency Magnetic Components
 */
class ThermalEquivalentCircuit {
private:
    ThermalModelConfiguration config;
    std::vector<ThermalNode> nodes;
    std::vector<ThermalResistanceElement> resistances;
    size_t ambientNodeId;
    
    // Conductance matrix and power vector for solving
    Eigen::MatrixXd conductanceMatrix;
    Eigen::VectorXd powerVector;
    Eigen::VectorXd temperatureVector;
    
    /**
     * @brief Build the thermal network from magnetic geometry
     */
    void buildNetwork(Magnetic magnetic, 
                      const std::map<std::string, double>& coreLossesPerElement,
                      const std::map<std::string, double>& windingLossesPerElement);
    
    /**
     * @brief Create nodes for the core
     */
    void createCoreNodes(Core core);
    
    /**
     * @brief Create nodes for the coil
     */
    void createCoilNodes(const Coil& coil);
    
    /**
     * @brief Create bobbin nodes
     */
    void createBobbinNodes(Bobbin bobbin);
    
    /**
     * @brief Create thermal resistances between nodes
     */
    void createThermalResistances(Magnetic magnetic);
    
    /**
     * @brief Assemble the conductance matrix
     */
    void assembleMatrix();
    
    /**
     * @brief Solve the thermal circuit for temperatures
     */
    void solveCircuit();
    
    /**
     * @brief Update temperature-dependent resistances
     */
    void updateResistances();
    
    /**
     * @brief Check convergence of temperature values
     */
    bool checkConvergence(const Eigen::VectorXd& oldTemperatures);

public:
    ThermalEquivalentCircuit() = default;
    ThermalEquivalentCircuit(const ThermalModelConfiguration& config);
    ~ThermalEquivalentCircuit() = default;
    
    /**
     * @brief Set the model configuration
     */
    void setConfiguration(const ThermalModelConfiguration& config);
    
    /**
     * @brief Calculate temperatures for a magnetic component
     * 
     * @param magnetic The magnetic component (core + coil)
     * @param coreLosses Total core losses (W)
     * @param windingLosses Total winding losses (W)
     * @return ThermalAnalysisOutput Results of the analysis
     */
    ThermalAnalysisOutput calculateTemperatures(
        Magnetic magnetic,
        double coreLosses,
        double windingLosses);
    
    /**
     * @brief Calculate temperatures with detailed loss distribution
     * 
     * @param magnetic The magnetic component
     * @param coreLossesPerElement Map of element name to loss (W)
     * @param windingLossesPerElement Map of layer/section name to loss (W)
     * @return ThermalAnalysisOutput Results of the analysis
     */
    ThermalAnalysisOutput calculateTemperaturesDetailed(
        Magnetic magnetic,
        const std::map<std::string, double>& coreLossesPerElement,
        const std::map<std::string, double>& windingLossesPerElement);
    
    /**
     * @brief Get the temperature at a specific point
     * 
     * @param coordinates 3D coordinates (m)
     * @return Temperature at that point (°C)
     */
    double getTemperatureAtPoint(const std::vector<double>& coordinates);
    
    /**
     * @brief Get the equivalent bulk thermal resistance
     * 
     * @return Thermal resistance from hotspot to ambient (K/W)
     */
    double getBulkThermalResistance() const;
    
    // Static utility methods for thermal resistance calculations
    
    /**
     * @brief Calculate conduction thermal resistance
     * R = L / (k * A)
     */
    static double calculateConductionResistance(double length, double thermalConductivity, double area);
    
    /**
     * @brief Calculate natural convection heat transfer coefficient
     * Uses Churchill-Chu correlation for vertical surfaces
     * Uses McAdams correlation for horizontal surfaces
     * 
     * @param surfaceTemperature Surface temperature (°C)
     * @param ambientTemperature Ambient temperature (°C)
     * @param characteristicLength Height for vertical, diameter for horizontal (m)
     * @param orientation Surface orientation
     * @return Heat transfer coefficient h (W/(m²·K))
     */
    static double calculateNaturalConvectionCoefficient(
        double surfaceTemperature,
        double ambientTemperature,
        double characteristicLength,
        SurfaceOrientation orientation);
    
    /**
     * @brief Calculate forced convection heat transfer coefficient
     * 
     * @param airVelocity Air velocity (m/s)
     * @param characteristicLength Characteristic length (m)
     * @param temperature Air temperature (°C)
     * @return Heat transfer coefficient h (W/(m²·K))
     */
    static double calculateForcedConvectionCoefficient(
        double airVelocity,
        double characteristicLength,
        double temperature);
    
    /**
     * @brief Calculate convection thermal resistance
     * R = 1 / (h * A)
     */
    static double calculateConvectionResistance(double heatTransferCoefficient, double area);
    
    /**
     * @brief Calculate radiation heat transfer coefficient
     * h_rad = ε * σ * (Ts² + T∞²) * (Ts + T∞)
     */
    static double calculateRadiationCoefficient(
        double surfaceTemperature,
        double ambientTemperature,
        double emissivity);
    
    /**
     * @brief Calculate radiation thermal resistance
     * R = 1 / (h_rad * A)
     */
    static double calculateRadiationResistance(
        double surfaceTemperature,
        double ambientTemperature,
        double emissivity,
        double area);
    
    /**
     * @brief Get thermal conductivity for common materials
     * Tries to look up from MAS databases first, then falls back to defaults
     */
    static double getMaterialThermalConductivity(const std::string& materialName);
    
    /**
     * @brief Get thermal conductivity for wire material with temperature interpolation
     * @param wireMaterial The wire material to get conductivity for
     * @param temperature Temperature in Celsius
     * @return Thermal conductivity in W/(m·K)
     */
    static double getWireMaterialThermalConductivity(const WireMaterial& wireMaterial, double temperature);
    
    /**
     * @brief Get thermal conductivity from core material (from MAS data)
     * @param coreMaterial The core material 
     * @return Thermal conductivity in W/(m·K), or default ferrite value if not available
     */
    static double getCoreMaterialThermalConductivity(const CoreMaterial& coreMaterial);
    
    /**
     * @brief Calculate gap thermal resistance
     * Gaps are mostly air with low thermal conductivity
     * @param gap The core gap geometry
     * @return Thermal resistance in K/W
     */
    static double calculateGapThermalResistance(const CoreGap& gap);
    
    // Accessors
    const std::vector<ThermalNode>& getNodes() const { return nodes; }
    const std::vector<ThermalResistanceElement>& getResistances() const { return resistances; }
};

/**
 * @brief Factory class for thermal models
 */
class ThermalModel {
public:
    enum class ModelType {
        BULK_MANIKTALA,           // Simple empirical model
        EQUIVALENT_CIRCUIT,       // Full thermal equivalent circuit
        EQUIVALENT_CIRCUIT_SIMPLE // Simplified equivalent circuit
    };
    
    /**
     * @brief Create a thermal model
     * @param type Type of model to create
     * @return Shared pointer to the model
     */
    static std::shared_ptr<ThermalEquivalentCircuit> factory(ModelType type = ModelType::EQUIVALENT_CIRCUIT);
};

} // namespace OpenMagnetics
