#include "physical_models/Temperature.h"
#include "physical_models/ThermalResistance.h"
#include "StrayCapacitance.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "Constants.h"
#include <filesystem>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <set>
#include <numbers>
#include <stdexcept>

namespace OpenMagnetics {

// Debug flag for thermal model - set to true for verbose output
constexpr bool THERMAL_DEBUG = false;

// Contact threshold: surfaces must be within this distance to conduct
constexpr double CONTACT_THRESHOLD_FACTOR = 0.25;  // wireDiameter / 4

// ============================================================================
// Wire Property Helper Functions (extract properties from Wire object locally)
// ============================================================================

/**
 * @brief Extract wire dimensions from Wire object
 * @return Pair of {width, height} in meters
 */
static std::pair<double, double> getWireDimensions(const Wire& wire) {
    double width = 0.001;   // Default 1mm
    double height = 0.001;  // Default 1mm
    
    bool isRound = (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ);
    
    // Get wire dimensions
    if (wire.get_conducting_diameter()) {
        auto condDiam = wire.get_conducting_diameter().value();
        if (condDiam.get_nominal()) {
            width = condDiam.get_nominal().value();
            height = width;  // Round wire
        }
    }
    
    if (wire.get_outer_diameter()) {
        auto outerDiam = wire.get_outer_diameter().value();
        if (outerDiam.get_nominal()) {
            double outer = outerDiam.get_nominal().value();
            if (!isRound) {
                width = outer;  // For rectangular, width = radial
            }
        }
    }
    
    // For rectangular wires, dimensions are estimated from conducting diameter
    if (wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        if (wire.get_outer_diameter()) {
            auto outerDiam = wire.get_outer_diameter().value();
            if (outerDiam.get_nominal()) {
                width = outerDiam.get_nominal().value();
                height = width * 0.5;  // Assume 2:1 aspect ratio
            }
        }
    }
    
    return {width, height};
}

/**
 * @brief Check if wire is round (including litz)
 */
static bool isRoundWire(const Wire& wire) {
    return (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ);
}

/**
 * @brief Get wire thermal conductivity from Wire object
 */
static double getWireThermalConductivity(const Wire& wire) {
    double thermalCond = 385.0;  // Default copper
    
    if (wire.get_material()) {
        auto materialVariant = wire.get_material().value();
        if (std::holds_alternative<std::string>(materialVariant)) {
            std::string materialName = std::get<std::string>(materialVariant);
            try {
                auto wireMaterial = find_wire_material_by_name(materialName);
                auto tc = wireMaterial.get_thermal_conductivity();
                if (tc && !tc->empty()) {
                    thermalCond = (*tc)[0].get_value();
                }
            } catch (...) {
                thermalCond = 385.0;
            }
        }
    }
    
    return thermalCond;
}

/**
 * @brief Calculate minimum distance for conduction detection
 */
static double getMinimumConductionDistance(double wireWidth, double wireHeight, bool round) {
    if (round) {
        return std::max(wireWidth, wireHeight) * 0.75;  // 75% of wire diameter
    } else {
        return std::min(wireWidth, wireHeight) * 0.75;  // 75% of wire thickness
    }
}

/**
 * @brief Calculate maximum distance for convection detection
 */
static double getMaximumConvectionDistance(double wireWidth, double wireHeight, bool round) {
    if (round) {
        return std::max(wireWidth, wireHeight);  // wire diameter
    } else {
        return std::min(wireWidth, wireHeight);  // min(width, height)
    }
}

/**
 * @brief Simple matrix class for thermal circuit solver
 */
class SimpleMatrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows_, cols_;
    
public:
    SimpleMatrix() : rows_(0), cols_(0) {}
    SimpleMatrix(size_t rows, size_t cols, double val = 0.0) 
        : data(rows, std::vector<double>(cols, val)), rows_(rows), cols_(cols) {}
    
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    
    double& operator()(size_t i, size_t j) { return data[i][j]; }
    const double& operator()(size_t i, size_t j) const { return data[i][j]; }
    
    void setZero() {
        for (auto& row : data) {
            std::fill(row.begin(), row.end(), 0.0);
        }
    }
    
    void setRowZero(size_t row) {
        std::fill(data[row].begin(), data[row].end(), 0.0);
    }
    
    void setColZero(size_t col) {
        for (size_t i = 0; i < rows_; ++i) {
            data[i][col] = 0.0;
        }
    }
    
    // Solve Ax = b using Gauss-Jordan elimination with partial pivoting
    static std::vector<double> solve(SimpleMatrix A, std::vector<double> b) {
        size_t n = A.rows();
        if (n == 0 || b.size() != n) {
            throw std::invalid_argument("Matrix dimensions mismatch");
        }
        
        // Create augmented matrix [A|b]
        std::vector<std::vector<double>> aug(n, std::vector<double>(n + 1));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                aug[i][j] = A(i, j);
            }
            aug[i][n] = b[i];
        }
        
        // Forward elimination with partial pivoting
        for (size_t col = 0; col < n; ++col) {
            // Find pivot
            size_t maxRow = col;
            double maxVal = std::abs(aug[col][col]);
            for (size_t row = col + 1; row < n; ++row) {
                if (std::abs(aug[row][col]) > maxVal) {
                    maxVal = std::abs(aug[row][col]);
                    maxRow = row;
                }
            }
            
            // Swap rows
            if (maxRow != col) {
                std::swap(aug[col], aug[maxRow]);
            }
            
            // Check for singular matrix
            if (std::abs(aug[col][col]) < 1e-15) {
                throw std::runtime_error("Matrix is singular or nearly singular");
            }
            
            // Eliminate column entries below pivot
            for (size_t row = col + 1; row < n; ++row) {
                double factor = aug[row][col] / aug[col][col];
                for (size_t j = col; j <= n; ++j) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }
        
        // Back substitution
        std::vector<double> x(n);
        for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
            x[i] = aug[i][n];
            for (size_t j = i + 1; j < n; ++j) {
                x[i] -= aug[i][j] * x[j];
            }
            x[i] /= aug[i][i];
        }
        
        return x;
    }
};

// ============================================================================
// Legacy API Functions
// ============================================================================

double Temperature::calculate_temperature_from_core_thermal_resistance(Core core, double totalLosses) {
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Core is missing processed description");
    }

    double thermalResistance;
    if (!core.get_processed_description()->get_thermal_resistance()) {
        auto thermalResistanceModel = CoreThermalResistanceModel::factory();
        thermalResistance = thermalResistanceModel->get_core_thermal_resistance_reluctance(core);
    }
    else {
        thermalResistance = core.get_processed_description()->get_thermal_resistance().value();
    }

    return thermalResistance * totalLosses;
}

double Temperature::calculate_temperature_from_core_thermal_resistance(double thermalResistance, double totalLosses) {
    return thermalResistance * totalLosses;
}

// ============================================================================
// Main Analysis Entry Point
// ============================================================================

ThermalResult Temperature::calculateTemperatures() {
    if (THERMAL_DEBUG) {
        std::cout << "Temperature::calculateTemperatures() started" << std::endl;
    }
    
    // Step 1: Extract wire properties from coil
    extractWireProperties();
    
    // Step 2: Create thermal nodes (core, turns, bobbin)
    createThermalNodes();
    
    // Step 3: Create thermal resistances between nodes
    createThermalResistances();
    
    // Step 4: Calculate schematic scaling
    calculateSchematicScaling();
    
    // Step 5: Create and solve thermal equivalent circuit
    ThermalResult result = solveThermalCircuit();
    
    // Step 6: Plot schematic if requested
    if (_config.plotSchematic) {
        plotSchematic();
    }
    
    if (THERMAL_DEBUG) {
        std::cout << "Temperature::calculateTemperatures() completed" << std::endl;
    }
    
    return result;
}

// ============================================================================
// Wire Property Extraction
// ============================================================================

void Temperature::extractWireProperties() {
    auto coil = _magnetic.get_coil();
    auto windings = coil.get_functional_description();
    
    if (windings.empty()) {
        return;
    }
    
    auto wireVariant = windings[0].get_wire();
    if (!std::holds_alternative<Wire>(wireVariant)) {
        return;
    }
    
    auto wire = std::get<Wire>(wireVariant);
    _isRoundWire = (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ);
    
    // Get wire dimensions
    if (wire.get_conducting_diameter()) {
        auto condDiam = wire.get_conducting_diameter().value();
        if (condDiam.get_nominal()) {
            _wireWidth = condDiam.get_nominal().value();
            _wireHeight = _wireWidth;  // Round wire
        }
    }
    
    if (wire.get_outer_diameter()) {
        auto outerDiam = wire.get_outer_diameter().value();
        if (outerDiam.get_nominal()) {
            double outer = outerDiam.get_nominal().value();
            if (!_isRoundWire) {
                _wireWidth = outer;  // For rectangular, width = radial
            }
        }
    }
    
    // For rectangular wires, dimensions are estimated from conducting diameter
    if (wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        if (wire.get_outer_diameter()) {
            auto outerDiam = wire.get_outer_diameter().value();
            if (outerDiam.get_nominal()) {
                _wireWidth = outerDiam.get_nominal().value();
                _wireHeight = _wireWidth * 0.5;  // Assume 2:1 aspect ratio
            }
        }
    }
    
    // Get wire material thermal conductivity
    if (wire.get_material()) {
        auto materialVariant = wire.get_material().value();
        if (std::holds_alternative<std::string>(materialVariant)) {
            std::string materialName = std::get<std::string>(materialVariant);
            try {
                auto wireMaterial = find_wire_material_by_name(materialName);
                auto thermalCond = wireMaterial.get_thermal_conductivity();
                if (thermalCond && !thermalCond->empty()) {
                    _wireThermalCond = (*thermalCond)[0].get_value();
                }
            } catch (...) {
                _wireThermalCond = 385.0;
            }
        }
    }
    
    // Get wire coating for thermal calculations
    _wireCoating = wire.resolve_coating();
    
    if (THERMAL_DEBUG) {
        std::cout << "Wire properties: width=" << _wireWidth * 1000 << "mm, height=" << _wireHeight * 1000 
                  << "mm, k=" << _wireThermalCond << " W/(m·K), round=" << _isRoundWire 
                  << ", minCondDist=" << getMinimumDistanceForConduction() * 1000 << "mm"
                  << ", maxConvDist=" << getMaximumDistanceForConvection() * 1000 << "mm" << std::endl;
    }
}

std::optional<Wire> Temperature::extractWire() const {
    auto coil = _magnetic.get_coil();
    auto windings = coil.get_functional_description();
    
    if (windings.empty()) {
        return std::nullopt;
    }
    
    auto wireVariant = windings[0].get_wire();
    if (!std::holds_alternative<Wire>(wireVariant)) {
        return std::nullopt;
    }
    
    return std::get<Wire>(wireVariant);
}

// ============================================================================
// Node Creation
// ============================================================================

void Temperature::createThermalNodes() {
    _nodes.clear();
    
    auto core = _magnetic.get_core();
    _isToroidal = (core.get_shape_family() == CoreShapeFamily::T);
    
    // Ensure coil is wound for proper turn coordinates
    if (!_magnetic.get_coil().get_turns_description()) {
        try {
            _magnetic.get_mutable_coil().wind();
            if (THERMAL_DEBUG) {
                std::cout << "Coil wound successfully for thermal analysis" << std::endl;
            }
        } catch (const std::exception& e) {
            if (THERMAL_DEBUG) {
                std::cout << "Coil winding failed: " << e.what() << std::endl;
            }
        }
    }
    
    if (_isToroidal) {
        createToroidalCoreNodes();
    } else {
        createConcentricCoreNodes();
    }
    
    createBobbinNodes();
    createInsulationLayerNodes();
    createTurnNodes();
    
    // Add ambient node (always last)
    ThermalNetworkNode ambientNode;
    ambientNode.part = ThermalNodePartType::AMBIENT;
    ambientNode.name = "Ambient";
    ambientNode.temperature = _config.ambientTemperature;
    ambientNode.powerDissipation = 0.0;
    ambientNode.physicalCoordinates = {0, 0, 0};
    _nodes.push_back(ambientNode);
    
    if (THERMAL_DEBUG) {
        std::cout << "Created " << _nodes.size() << " thermal nodes (including ambient)" << std::endl;
    }
}

void Temperature::createToroidalCoreNodes() {
    auto core = _magnetic.get_core();
    auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    
    double outerDiameter = dimensions["A"];
    double innerDiameter = dimensions["B"];
    double height = dimensions["C"];
    
    double windingWindowInnerRadius = innerDiameter / 2.0;
    double windingWindowOuterRadius = outerDiameter / 2.0;
    double meanRadius = (windingWindowInnerRadius + windingWindowOuterRadius) / 2.0;
    
    size_t numSegments = _config.toroidalSegments;
    double angularStep = 2.0 * M_PI / numSegments;
    
    double coreK = _config.coreThermalConductivity;
    
    // Pre-calculate turn positions and widths for segment-specific coverage
    // Map each turn to its angular position and width
    struct TurnCoverageInfo {
        double angle;           // Angular position (radians)
        double innerWidth;      // Width at inner radius
        double outerWidth;      // Width at outer radius
        bool isInner;           // Is this turn on inner or outer half?
    };
    std::vector<TurnCoverageInfo> turnCoverageInfo;
    
    auto coil = _magnetic.get_coil();
    auto turns = coil.get_turns_description();
    if (turns) {
        for (const auto& turn : *turns) {
            auto coords = turn.get_coordinates();
            if (coords.size() >= 2) {
                // Calculate angular position from turn center coordinates
                double turnAngle = std::atan2(coords[1], coords[0]);
                
                // Determine if turn is on inner or outer half based on radius
                double turnRadius = std::sqrt(coords[0]*coords[0] + coords[1]*coords[1]);
                bool isInner = turnRadius < meanRadius;
                
                // Calculate turn width at inner and outer radii
                // Use turn's actual dimensions if available, otherwise use defaults
                double innerWidth = 0.001;  // Default 1mm
                double outerWidth = 0.001;
                if (turn.get_dimensions() && turn.get_dimensions()->size() >= 1) {
                    innerWidth = (*turn.get_dimensions())[0];
                    outerWidth = innerWidth;
                }
                
                turnCoverageInfo.push_back({turnAngle, innerWidth, outerWidth, isInner});
            }
        }
    }
    
    for (size_t i = 0; i < numSegments; ++i) {
        double angle = i * angularStep;
        double midAngle = angle + angularStep / 2.0;
        
        ThermalNetworkNode node;
        node.part = ThermalNodePartType::CORE_TOROIDAL_SEGMENT;
        node.name = "Core_Segment_" + std::to_string(i);
        node.temperature = _config.ambientTemperature;
        node.powerDissipation = _config.coreLosses / numSegments;
        
        node.physicalCoordinates = {
            meanRadius * std::cos(midAngle),
            meanRadius * std::sin(midAngle),
            0
        };
        
        // Initialize quadrants for core segment
        double segmentArcLength = meanRadius * angularStep;
        double radialDepth = windingWindowOuterRadius - windingWindowInnerRadius;
        double radialLength = radialDepth / 2.0;
        
        double innerCircumference = 2.0 * M_PI * windingWindowInnerRadius;
        double radialInnerArea = (innerCircumference / numSegments) * height;
        
        double outerCircumference = 2.0 * M_PI * windingWindowOuterRadius;
        double radialOuterArea = (outerCircumference / numSegments) * height;
        
        double tangentialArea = radialDepth * height;
        
        // Calculate segment-specific surface coverage
        // Find turns that are adjacent to this core segment's angular span
        std::vector<double> innerTurnWidths;
        std::vector<double> outerTurnWidths;
        
        for (const auto& turnInfo : turnCoverageInfo) {
            // Normalize turn angle to [0, 2π)
            double turnAngle = turnInfo.angle;
            while (turnAngle < 0) turnAngle += 2.0 * M_PI;
            
            // Check if turn is within this segment's angular span [angle, angle + angularStep)
            // Add small tolerance for edge cases
            double tol = angularStep * 0.1;
            double segStart = angle - tol;
            double segEnd = angle + angularStep + tol;
            
            bool inSegment = false;
            if (segStart <= turnAngle && turnAngle < segEnd) {
                inSegment = true;
            } else if (segStart < 0 && turnAngle >= (2.0 * M_PI + segStart)) {
                // Handle wrap-around at 0/2π boundary
                inSegment = true;
            } else if (segEnd > 2.0 * M_PI && turnAngle < (segEnd - 2.0 * M_PI)) {
                // Handle wrap-around at 0/2π boundary
                inSegment = true;
            }
            
            if (inSegment) {
                innerTurnWidths.push_back(turnInfo.innerWidth);
                outerTurnWidths.push_back(turnInfo.outerWidth);
            }
        }
        
        // Calculate surface coverage for inner and outer radial surfaces
        // Inner surface: turns cover from the inside
        double innerCoverage = ThermalNetworkNode::calculateToroidalSurfaceCoverage(
            windingWindowInnerRadius, angularStep, innerTurnWidths);
        
        // Outer surface: turns cover from the outside
        double outerCoverage = ThermalNetworkNode::calculateToroidalSurfaceCoverage(
            windingWindowOuterRadius, angularStep, outerTurnWidths);
        
        node.quadrants[0].face = ThermalNodeFace::RADIAL_INNER;
        node.quadrants[0].surfaceArea = radialInnerArea;
        node.quadrants[0].length = radialLength;
        node.quadrants[0].thermalConductivity = coreK;
        node.quadrants[0].surfaceCoverage = innerCoverage;
        
        node.quadrants[1].face = ThermalNodeFace::RADIAL_OUTER;
        node.quadrants[1].surfaceArea = radialOuterArea;
        node.quadrants[1].length = radialLength;
        node.quadrants[1].thermalConductivity = coreK;
        node.quadrants[1].surfaceCoverage = outerCoverage;
        
        node.quadrants[2].face = ThermalNodeFace::TANGENTIAL_LEFT;
        node.quadrants[2].surfaceArea = tangentialArea;
        node.quadrants[2].length = segmentArcLength / 2.0;
        node.quadrants[2].thermalConductivity = coreK;
        node.quadrants[2].surfaceCoverage = 1.0;  // Tangential faces fully exposed
        
        node.quadrants[3].face = ThermalNodeFace::TANGENTIAL_RIGHT;
        node.quadrants[3].surfaceArea = tangentialArea;
        node.quadrants[3].length = segmentArcLength / 2.0;
        node.quadrants[3].thermalConductivity = coreK;
        node.quadrants[3].surfaceCoverage = 1.0;  // Tangential faces fully exposed
        
        // Set limit coordinates for core segment quadrants (surface positions for conduction detection)
        // RADIAL_INNER: at inner radius facing toward winding window
        node.quadrants[0].limitCoordinates = {
            windingWindowInnerRadius * std::cos(midAngle),
            windingWindowInnerRadius * std::sin(midAngle),
            0.0
        };
        // RADIAL_OUTER: at outer radius facing away from center
        node.quadrants[1].limitCoordinates = {
            windingWindowOuterRadius * std::cos(midAngle),
            windingWindowOuterRadius * std::sin(midAngle),
            0.0
        };
        // TANGENTIAL_LEFT: at segment start angle
        node.quadrants[2].limitCoordinates = {
            meanRadius * std::cos(angle),
            meanRadius * std::sin(angle),
            0.0
        };
        // TANGENTIAL_RIGHT: at segment end angle
        node.quadrants[3].limitCoordinates = {
            meanRadius * std::cos(angle + angularStep),
            meanRadius * std::sin(angle + angularStep),
            0.0
        };
        
        _nodes.push_back(node);
    }
}

void Temperature::createConcentricCoreNodes() {
    auto core = _magnetic.get_core();
    auto columns = core.get_columns();
    auto processedDesc = core.get_processed_description();
    
    if (columns.empty() || !processedDesc) {
        return;
    }
    
    double coreK = _config.coreThermalConductivity;
    double coreWidth = processedDesc->get_width();
    double coreHeight = processedDesc->get_height();
    double coreDepth = processedDesc->get_depth();
    
    // Get gaps for each column
    auto gaps = core.get_functional_description().get_gapping();
    
    // Find main (central) column and right lateral column
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    auto rightColumn = core.find_closest_column_by_coordinates(
        std::vector<double>({coreWidth / 2, 0, -coreDepth / 2}));
    
    // Helper to count gaps in a column
    auto countGapsInColumn = [&](const ColumnElement& col) -> int {
        int count = 0;
        for (const auto& gap : gaps) {
            if (!gap.get_coordinates()) continue;
            auto gapCoords = gap.get_coordinates().value();
            // Check if gap is in this column (by x-coordinate proximity)
            double colX = (col.get_type() == ColumnType::CENTRAL) ? 0 : 
                         (col.get_coordinates()[0] > 0 ? coreWidth/2 : -coreWidth/2);
            if (std::abs(gapCoords[0] - colX) < coreWidth * 0.25) {
                count++;
            }
        }
        return count;
    };
    
    // Create central column node(s) - using HALF depth for symmetry (left half only)
    int mainColGaps = countGapsInColumn(mainColumn);
    double halfDepth = mainColumn.get_depth() / 2.0;  // Model half the core
    if (mainColGaps == 0) {
        // Single central column node
        ThermalNetworkNode node;
        node.part = ThermalNodePartType::CORE_CENTRAL_COLUMN;
        node.name = "Core_CentralColumn";
        node.temperature = _config.ambientTemperature;
        node.powerDissipation = _config.coreLosses * 0.4;  // 40% of core losses
        node.physicalCoordinates = {0, 0, 0};
        node.initializeConcentricCoreQuadrants(mainColumn.get_width(), 
                                                mainColumn.get_height(), halfDepth, coreK);
        _nodes.push_back(node);
    } else {
        // Split central column by gaps
        double chunkHeight = mainColumn.get_height() / (mainColGaps + 1);
        for (int i = 0; i <= mainColGaps; ++i) {
            ThermalNetworkNode node;
            node.part = ThermalNodePartType::CORE_CENTRAL_COLUMN;
            node.name = "Core_CentralColumn_Chunk_" + std::to_string(i);
            node.temperature = _config.ambientTemperature;
            node.powerDissipation = (_config.coreLosses * 0.4) / (mainColGaps + 1);
            // Position each chunk along Y axis
            double yPos = -mainColumn.get_height()/2 + chunkHeight * (i + 0.5);
            node.physicalCoordinates = {0, yPos, 0};
            node.initializeConcentricCoreQuadrants(mainColumn.get_width(), 
                                                    chunkHeight, halfDepth, coreK);
            _nodes.push_back(node);
        }
    }
    
    // Create top yoke node - using HALF depth for symmetry
    // Position between center (x=0) and right lateral column (x=coreWidth/2)
    ThermalNetworkNode topYoke;
    topYoke.part = ThermalNodePartType::CORE_TOP_YOKE;
    topYoke.name = "Core_TopYoke";
    topYoke.temperature = _config.ambientTemperature;
    topYoke.powerDissipation = _config.coreLosses * 0.1;  // 10% of core losses
    topYoke.physicalCoordinates = {coreWidth / 4, coreHeight / 2 - mainColumn.get_width()/4, 0};
    topYoke.initializeConcentricCoreQuadrants(coreWidth / 2, mainColumn.get_width()/2, coreDepth / 2, coreK);
    _nodes.push_back(topYoke);
    
    // Create bottom yoke node - using HALF depth for symmetry
    // Position between center (x=0) and right lateral column (x=coreWidth/2)
    ThermalNetworkNode bottomYoke;
    bottomYoke.part = ThermalNodePartType::CORE_BOTTOM_YOKE;
    bottomYoke.name = "Core_BottomYoke";
    bottomYoke.temperature = _config.ambientTemperature;
    bottomYoke.powerDissipation = _config.coreLosses * 0.1;  // 10% of core losses
    bottomYoke.physicalCoordinates = {coreWidth / 4, -coreHeight / 2 + mainColumn.get_width()/4, 0};
    bottomYoke.initializeConcentricCoreQuadrants(coreWidth / 2, mainColumn.get_width()/2, coreDepth / 2, coreK);
    _nodes.push_back(bottomYoke);
    
    // Create lateral column node (RIGHT side only for symmetry - half core model)
    if (columns.size() > 1) {
        int latColGaps = countGapsInColumn(rightColumn);
        double halfDepth = rightColumn.get_depth() / 2.0;  // Model half the core
        
        double offset = coreWidth / 2 - rightColumn.get_width() / 2;
        
        if (latColGaps == 0) {
            // Single lateral column node (RIGHT side only)
            ThermalNetworkNode rightNode;
            rightNode.part = ThermalNodePartType::CORE_LATERAL_COLUMN;
            rightNode.name = "Core_LateralColumn_Right";
            rightNode.temperature = _config.ambientTemperature;
            rightNode.powerDissipation = _config.coreLosses * 0.2;  // 20% of core losses
            rightNode.physicalCoordinates = {offset, 0, 0};
            rightNode.initializeConcentricCoreQuadrants(rightColumn.get_width(), 
                                                         rightColumn.get_height(), halfDepth, coreK);
            _nodes.push_back(rightNode);
        } else {
            // Split lateral column by gaps
            double chunkHeight = rightColumn.get_height() / (latColGaps + 1);
            for (int i = 0; i <= latColGaps; ++i) {
                double yPos = -rightColumn.get_height()/2 + chunkHeight * (i + 0.5);
                
                // Right lateral chunk only
                ThermalNetworkNode rightNode;
                rightNode.part = ThermalNodePartType::CORE_LATERAL_COLUMN;
                rightNode.name = "Core_LateralColumn_Right_Chunk_" + std::to_string(i);
                rightNode.temperature = _config.ambientTemperature;
                rightNode.powerDissipation = (_config.coreLosses * 0.2) / (latColGaps + 1);
                rightNode.physicalCoordinates = {offset, yPos, 0};
                rightNode.initializeConcentricCoreQuadrants(rightColumn.get_width(), 
                                                             chunkHeight, halfDepth, coreK);
                _nodes.push_back(rightNode);
            }
        }
    }
}

void Temperature::createBobbinNodes() {
    auto bobbinOpt = _magnetic.get_mutable_coil().get_bobbin();
    if (!std::holds_alternative<Bobbin>(bobbinOpt)) {
        return;
    }
    
    if (_isToroidal) {
        return;
    }
    
    auto bobbin = std::get<Bobbin>(bobbinOpt);
    auto processedDesc = bobbin.get_processed_description();
    if (!processedDesc) {
        return;
    }
    
    double bobbinK = 0.2;  // Thermal conductivity of typical bobbin material (W/m·K)
    double wallThickness = processedDesc->get_wall_thickness();
    
    // Only create bobbin nodes if there's actual wall thickness
    if (wallThickness <= 0) {
        return;
    }
    
    auto core = _magnetic.get_core();
    auto processedCore = core.get_processed_description();
    if (!processedCore) {
        return;
    }
    
    double coreWidth = processedCore->get_width();
    double coreHeight = processedCore->get_height();
    double coreDepth = processedCore->get_depth();
    
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    double columnWidth = mainColumn.get_width();
    
    // Get winding window for proper bobbin positioning
    auto windingWindows = core.get_winding_windows();
    double windingWindowWidth = 0.01;  // Default 10mm
    double windingWindowHeight = coreHeight * 0.8;
    if (!windingWindows.empty()) {
        windingWindowWidth = windingWindows[0].get_width().value_or(0.01);
        windingWindowHeight = windingWindows[0].get_height().value_or(coreHeight * 0.8);
    }
    
    // Calculate bobbin positions
    // Bobbin column is at the INNER edge of the winding window (closer to center)
    double bobbinColumnX = columnWidth / 2 + wallThickness / 2;
    
    // Bobbin yokes are aligned with core yokes (same X coordinate)
    // Core yokes are at coreWidth/4, so bobbin yokes should be at the same X
    double bobbinYokeX = coreWidth / 4;
    
    // 1. Bobbin Central Column Wall
    // Positioned at the inner surface of the central column winding window
    ThermalNetworkNode columnWall;
    columnWall.part = ThermalNodePartType::BOBBIN_CENTRAL_COLUMN;
    columnWall.name = "Bobbin_CentralColumn_Wall";
    columnWall.temperature = _config.ambientTemperature;
    columnWall.powerDissipation = 0.0;
    // Position at the winding window inner edge (closer to center)
    columnWall.physicalCoordinates = {bobbinColumnX, 0, 0};
    columnWall.initializeConcentricCoreQuadrants(wallThickness, windingWindowHeight, coreDepth / 2, bobbinK);
    
    // Calculate surface coverage for the RIGHT face (facing turns)
    // Find turns that are adjacent to the bobbin central column
    std::vector<double> turnHeights;
    auto coil = _magnetic.get_mutable_coil();
    if (coil.get_turns_description()) {
        auto turns = coil.get_turns_description().value();
        for (const auto& turn : turns) {
            auto coords = turn.get_coordinates();
            if (coords.size() >= 2) {
                double turnX = coords[0];
                double turnY = coords[1];
                double turnWidth = 0;
                double turnHeight = 0;
                if (turn.get_dimensions().has_value() && turn.get_dimensions().value().size() >= 2) {
                    turnWidth = turn.get_dimensions().value()[0];   // X dimension
                    turnHeight = turn.get_dimensions().value()[1];  // Y dimension
                }
                // Check if turn is close to bobbin (within turn width + bobbin thickness)
                double turnLeftEdge = turnX - turnWidth / 2;
                double bobbinRightEdge = bobbinColumnX + wallThickness / 2;
                if (std::abs(turnLeftEdge - bobbinRightEdge) < (turnWidth + wallThickness) && 
                    turnY >= -windingWindowHeight/2 && turnY <= windingWindowHeight/2) {
                    turnHeights.push_back(turnHeight);
                }
            }
        }
    }
    double rightFaceCoverage = ThermalNetworkNode::calculateConcentricSurfaceCoverage(
        windingWindowHeight, turnHeights);
    columnWall.quadrants[0].surfaceCoverage = rightFaceCoverage;  // RIGHT face
    
    _nodes.push_back(columnWall);
    
    // 2. Bobbin Top Yoke Wall
    // Positioned at the top of the winding window
    ThermalNetworkNode topYokeWall;
    topYokeWall.part = ThermalNodePartType::BOBBIN_TOP_YOKE;
    topYokeWall.name = "Bobbin_TopYoke_Wall";
    topYokeWall.temperature = _config.ambientTemperature;
    topYokeWall.powerDissipation = 0.0;
    // X: outer edge of winding window, Y: top of winding window
    double windingWindowTop = windingWindowHeight / 2;
    topYokeWall.physicalCoordinates = {bobbinYokeX, windingWindowTop - wallThickness / 2, 0};
    topYokeWall.initializeConcentricCoreQuadrants(windingWindowWidth / 2, wallThickness, coreDepth / 2, bobbinK);
    
    // Calculate surface coverage for the RIGHT face (facing turns)
    std::vector<double> topYokeTurnWidths;
    if (coil.get_turns_description()) {
        auto turns = coil.get_turns_description().value();
        for (const auto& turn : turns) {
            auto coords = turn.get_coordinates();
            if (coords.size() >= 2) {
                double turnX = coords[0];
                double turnY = coords[1];
                double turnWidth = 0;
                double turnHeight = 0;
                if (turn.get_dimensions().has_value() && turn.get_dimensions().value().size() >= 2) {
                    turnWidth = turn.get_dimensions().value()[0];   // X dimension
                    turnHeight = turn.get_dimensions().value()[1];  // Y dimension
                }
                // Check if turn is near the top yoke's right face
                double turnBottom = turnY - turnHeight / 2;
                double yokeTop = windingWindowTop;
                if (std::abs(turnBottom - yokeTop) < (turnHeight + wallThickness) &&
                    turnX >= bobbinYokeX - windingWindowWidth/4 && turnX <= bobbinYokeX + windingWindowWidth/2) {
                    topYokeTurnWidths.push_back(turnWidth);
                }
            }
        }
    }
    double topYokeRightCoverage = ThermalNetworkNode::calculateConcentricSurfaceCoverage(
        windingWindowWidth / 2, topYokeTurnWidths);
    topYokeWall.quadrants[0].surfaceCoverage = topYokeRightCoverage;  // RIGHT face
    
    _nodes.push_back(topYokeWall);
    
    // 3. Bobbin Bottom Yoke Wall
    // Positioned at the bottom of the winding window
    ThermalNetworkNode bottomYokeWall;
    bottomYokeWall.part = ThermalNodePartType::BOBBIN_BOTTOM_YOKE;
    bottomYokeWall.name = "Bobbin_BottomYoke_Wall";
    bottomYokeWall.temperature = _config.ambientTemperature;
    bottomYokeWall.powerDissipation = 0.0;
    double windingWindowBottom = -windingWindowHeight / 2;
    bottomYokeWall.physicalCoordinates = {bobbinYokeX, windingWindowBottom + wallThickness / 2, 0};
    bottomYokeWall.initializeConcentricCoreQuadrants(windingWindowWidth / 2, wallThickness, coreDepth / 2, bobbinK);
    
    // Calculate surface coverage for the RIGHT face (facing turns)
    std::vector<double> bottomYokeTurnWidths;
    if (coil.get_turns_description()) {
        auto turns = coil.get_turns_description().value();
        for (const auto& turn : turns) {
            auto coords = turn.get_coordinates();
            if (coords.size() >= 2) {
                double turnX = coords[0];
                double turnY = coords[1];
                double turnWidth = 0;
                double turnHeight = 0;
                if (turn.get_dimensions().has_value() && turn.get_dimensions().value().size() >= 2) {
                    turnWidth = turn.get_dimensions().value()[0];   // X dimension
                    turnHeight = turn.get_dimensions().value()[1];  // Y dimension
                }
                // Check if turn is near the bottom yoke's right face
                double turnTop = turnY + turnHeight / 2;
                double yokeBottom = windingWindowBottom;
                if (std::abs(turnTop - yokeBottom) < (turnHeight + wallThickness) &&
                    turnX >= bobbinYokeX - windingWindowWidth/4 && turnX <= bobbinYokeX + windingWindowWidth/2) {
                    bottomYokeTurnWidths.push_back(turnWidth);
                }
            }
        }
    }
    double bottomYokeRightCoverage = ThermalNetworkNode::calculateConcentricSurfaceCoverage(
        windingWindowWidth / 2, bottomYokeTurnWidths);
    bottomYokeWall.quadrants[0].surfaceCoverage = bottomYokeRightCoverage;  // RIGHT face
    
    _nodes.push_back(bottomYokeWall);
}

void Temperature::createInsulationLayerNodes() {
    // Skip insulation layers for toroidal cores - coordinate system is complex
    // and requires proper handling of polar coordinates relative to section positions
    if (_isToroidal) {
        if (THERMAL_DEBUG) {
            std::cout << "Skipping insulation layer nodes for toroidal core" << std::endl;
        }
        return;
    }
    
    auto coil = _magnetic.get_mutable_coil();
    
    // Check if we have layers description
    if (!coil.get_layers_description()) {
        if (THERMAL_DEBUG) {
            std::cout << "No layers description, skipping insulation layer nodes" << std::endl;
        }
        return;
    }
    
    // Get all insulation layers
    auto insulationLayers = coil.get_layers_description_insulation();
    if (insulationLayers.empty()) {
        if (THERMAL_DEBUG) {
            std::cout << "No insulation layers found" << std::endl;
        }
        return;
    }
    
    if (THERMAL_DEBUG) {
        std::cout << "Found " << insulationLayers.size() << " insulation layers" << std::endl;
    }
    
    // Get core depth for insulation layer depth
    auto core = _magnetic.get_core();
    auto processedCore = core.get_processed_description();
    double coreDepth = 0.01;  // Default 10mm
    if (processedCore) {
        coreDepth = processedCore->get_depth();
    }
    
    size_t layerIdx = 0;
    size_t createdCount = 0;
    for (const auto& layer : insulationLayers) {
        // Skip layers without proper dimensions or coordinates
        if (layer.get_dimensions().size() < 2) {
            continue;
        }
        if (layer.get_coordinates().size() < 2) {
            continue;
        }
        
        // Get layer center coordinates
        double layerX = layer.get_coordinates()[0];
        double layerY = layer.get_coordinates()[1];
        
        // Get layer geometry
        double layerWidth = layer.get_dimensions()[0];   // X dimension (thickness)
        double layerHeight = layer.get_dimensions()[1];  // Y dimension (span)
        double layerDepth = coreDepth / 2.0;  // Half depth for single side modeling
        
        // Check coordinate system and convert if necessary
        auto coordSystem = layer.get_coordinate_system();
        if (coordSystem.has_value() && coordSystem.value() == MAS::CoordinateSystem::POLAR) {
            // Polar coordinates: layerX = radius, layerY = angle in degrees
            // Dimensions: width = radial thickness, height = angular span (in degrees * 1000)
            double radius = layerX;
            double angleDeg = layerY;
            double angleRad = angleDeg * M_PI / 180.0;
            
            // Convert center position to Cartesian (relative to section center)
            layerX = radius * std::cos(angleRad);
            layerY = radius * std::sin(angleRad);
            
            // Convert dimensions: width is radial (already in meters), height is arc length
            double radialThickness = layerWidth;  // Already in meters (radial direction)
            double angularSpanDeg = layerHeight / 1000.0;  // Convert from millidegrees to degrees
            double angularSpanRad = angularSpanDeg * M_PI / 180.0;
            double arcLength = radius * angularSpanRad;  // Arc length = radius * angle
            
            // For insulation layers in polar coordinates:
            // - layerWidth becomes the radial thickness
            // - layerHeight becomes the arc length (circumferential span)
            layerWidth = radialThickness;
            layerHeight = arcLength;
        }
        
        // Get layer thermal conductivity from material
        double layerK = 0.2;  // Default for typical insulation material
        try {
            auto insulationMaterial = coil.resolve_insulation_layer_insulation_material(layer);
            if (insulationMaterial.get_thermal_conductivity()) {
                layerK = insulationMaterial.get_thermal_conductivity().value();
            }
        } catch (...) {
            // Use default if material cannot be resolved
        }
        
        // Create insulation layer node
        ThermalNetworkNode insulationNode;
        insulationNode.part = ThermalNodePartType::INSULATION_LAYER;
        insulationNode.name = "L_" + std::to_string(layerIdx);
        insulationNode.insulationLayerIndex = layerIdx;
        insulationNode.temperature = _config.ambientTemperature;
        insulationNode.powerDissipation = 0.0;  // No heat generation in insulation
        insulationNode.physicalCoordinates = {layerX, layerY, 0};
        
        // Initialize quadrants for rectangular insulation layer
        insulationNode.initializeInsulationLayerQuadrants(layerWidth, layerHeight, layerDepth, layerK);
        
        _nodes.push_back(insulationNode);
        layerIdx++;
        createdCount++;
        
        if (THERMAL_DEBUG) {
            std::cout << "Created insulation layer node: " << insulationNode.name 
                      << " at (" << layerX * 1000 << "mm, " << layerY * 1000 << "mm)"
                      << " size (" << layerWidth * 1000 << "mm x " << layerHeight * 1000 << "mm)"
                      << std::endl;
        }
    }
    
    if (THERMAL_DEBUG) {
        std::cout << "Created " << createdCount << " insulation layer nodes" << std::endl;
    }
}

void Temperature::createTurnNodes() {
    auto coil = _magnetic.get_coil();
    auto turns = coil.get_turns_description();
    
    if (!turns || turns->empty()) {
        return;
    }
    
    // Extract wire properties locally from the coil
    auto wireOpt = extractWire();
    double defaultWireWidth = 0.001;
    double defaultWireHeight = 0.001;
    double defaultWireThermalCond = 385.0;
    bool defaultIsRoundWire = false;    std::optional<InsulationWireCoating> defaultWireCoating;
    
    if (wireOpt) {
        auto [w, h] = getWireDimensions(wireOpt.value());
        defaultWireWidth = w;
        defaultWireHeight = h;
        defaultWireThermalCond = getWireThermalConductivity(wireOpt.value());
        defaultIsRoundWire = isRoundWire(wireOpt.value());
        defaultWireCoating = wireOpt.value().resolve_coating();
    }
    
    if (_isToroidal) {
        auto core = _magnetic.get_core();
        auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
        
        double outerDiameter = dimensions["A"];
        double innerDiameter = dimensions["B"];
        double windingWindowInnerRadius = innerDiameter / 2.0;
        double windingWindowOuterRadius = outerDiameter / 2.0;
        double meanRadius = (windingWindowInnerRadius + windingWindowOuterRadius) / 2.0;
        
        // Get per-turn losses from simulation output
        // Real losses per turn are REQUIRED - no mock/equal distribution allowed
        std::vector<double> turnLosses;
        
        if (!_config.windingLossesOutput) {
            throw std::runtime_error("WindingLossesOutput is required for thermal analysis. "
                                     "Use MagneticSimulator to calculate real losses per turn.");
        }
        
        auto lossesPerTurn = _config.windingLossesOutput->get_winding_losses_per_turn();
        if (!lossesPerTurn || lossesPerTurn->empty()) {
            throw std::runtime_error("WindingLossesOutput must contain per-turn losses (winding_losses_per_turn). "
                                     "Use MagneticSimulator to calculate real losses per turn.");
        }
        
        for (const auto& elem : *lossesPerTurn) {
            double loss = 0.0;
            if (elem.get_ohmic_losses()) {
                loss += elem.get_ohmic_losses()->get_losses();
            }
            if (elem.get_skin_effect_losses()) {
                auto harmonics = elem.get_skin_effect_losses()->get_losses_per_harmonic();
                for (double h : harmonics) {
                    loss += h;
                }
            }
            if (elem.get_proximity_effect_losses()) {
                auto harmonics = elem.get_proximity_effect_losses()->get_losses_per_harmonic();
                for (double h : harmonics) {
                    loss += h;
                }
            }
            turnLosses.push_back(loss);
        }
        
        for (size_t t = 0; t < turns->size(); ++t) {
            const auto& turn = (*turns)[t];
            double turnLoss = (t < turnLosses.size()) ? turnLosses[t] : 0.0;
            
            auto coords = turn.get_coordinates();
            double angle = 0;
            double turnCenterRadius = meanRadius;  // fallback
            if (coords.size() >= 2) {
                angle = std::atan2(coords[1], coords[0]);
                // Use actual turn radius from coordinates
                turnCenterRadius = std::sqrt(coords[0]*coords[0] + coords[1]*coords[1]);
            }
            
            // Get wire dimensions from turn (use turn-specific if available, otherwise use defaults)
            double wireWidth = defaultWireWidth;
            double wireHeight = defaultWireHeight;
            if (turn.get_dimensions() && turn.get_dimensions()->size() >= 2) {
                wireWidth = (*turn.get_dimensions())[0];   // radial dimension
                wireHeight = (*turn.get_dimensions())[1];  // axial dimension
            }
            
            double totalTurnLength = turn.get_length();
            double halfLength = totalTurnLength / 2.0;
            
            // Check if we have proper additionalCoordinates for this turn
            bool hasAdditionalCoords = turn.get_additional_coordinates().has_value() && 
                                       !turn.get_additional_coordinates().value().empty() &&
                                       turn.get_additional_coordinates().value()[0].size() >= 2;
            
            // INNER NODE - at turn.get_coordinates() (inner surface, facing toward center/winding window)
            ThermalNetworkNode innerNode;
            innerNode.part = ThermalNodePartType::TURN;
            innerNode.name = turn.get_name() + "_Inner";
            innerNode.temperature = _config.ambientTemperature;
            innerNode.turnIndex = static_cast<int>(t);
            innerNode.isInnerTurn = true;
            
            // Use actual turn coordinates for inner surface
            double innerSurfaceX = coords[0];
            double innerSurfaceY = coords[1];
            double innerSurfaceRadius = turnCenterRadius;
            
            innerNode.physicalCoordinates = {innerSurfaceX, innerSurfaceY, 0};
            
            // If no additionalCoordinates, inner node gets ALL the loss and full length
            // If additionalCoordinates exist, split loss and length between inner and outer nodes
            if (hasAdditionalCoords) {
                innerNode.powerDissipation = turnLoss / 2.0;
                innerNode.initializeToroidalQuadrants(wireWidth, wireHeight, halfLength, 
                                                       defaultWireThermalCond, true, innerSurfaceRadius,
                                                       defaultWireCoating,
                                                       defaultIsRoundWire ? TurnCrossSectionalShape::ROUND : TurnCrossSectionalShape::RECTANGULAR);
            } else {
                innerNode.powerDissipation = turnLoss;  // All loss to inner node
                innerNode.initializeToroidalQuadrants(wireWidth, wireHeight, totalTurnLength, 
                                                       defaultWireThermalCond, true, innerSurfaceRadius,
                                                       defaultWireCoating,
                                                       defaultIsRoundWire ? TurnCrossSectionalShape::ROUND : TurnCrossSectionalShape::RECTANGULAR);
            }
            _nodes.push_back(innerNode);
            
            // OUTER NODE - only created if additionalCoordinates exist (proper turn geometry)
            if (hasAdditionalCoords) {
                ThermalNetworkNode outerNode;
                outerNode.part = ThermalNodePartType::TURN;
                outerNode.name = turn.get_name() + "_Outer";
                outerNode.temperature = _config.ambientTemperature;
                outerNode.powerDissipation = turnLoss / 2.0;
                outerNode.turnIndex = static_cast<int>(t);
                outerNode.isInnerTurn = false;
                
                auto addCoords = turn.get_additional_coordinates().value()[0];
                double outerSurfaceX = addCoords[0];
                double outerSurfaceY = addCoords[1];
                double outerSurfaceRadius = std::sqrt(outerSurfaceX*outerSurfaceX + outerSurfaceY*outerSurfaceY);
                
                outerNode.physicalCoordinates = {outerSurfaceX, outerSurfaceY, 0};
                
                outerNode.initializeToroidalQuadrants(wireWidth, wireHeight, halfLength,
                                                       defaultWireThermalCond, false, outerSurfaceRadius,
                                                       defaultWireCoating,
                                                       defaultIsRoundWire ? TurnCrossSectionalShape::ROUND : TurnCrossSectionalShape::RECTANGULAR);
                _nodes.push_back(outerNode);
            }
        }
    } else {
        // Concentric cores - single node per turn
        for (size_t t = 0; t < turns->size(); ++t) {
            const auto& turn = (*turns)[t];
            
            ThermalNetworkNode node;
            node.part = ThermalNodePartType::TURN;
            node.name = turn.get_name();
            node.temperature = _config.ambientTemperature;
            node.powerDissipation = _config.windingLosses / turns->size();
            node.turnIndex = static_cast<int>(t);
            
            auto coords = turn.get_coordinates();
            if (coords.size() >= 3) {
                node.physicalCoordinates = {coords[0], coords[1], coords[2]};
            } else if (coords.size() >= 2) {
                node.physicalCoordinates = {coords[0], coords[1], 0};
            } else {
                node.physicalCoordinates = {0, 0, 0};
            }
            
            // Get wire dimensions from turn (each turn may have different wire size)
            double wireWidth = defaultWireWidth;
            double wireHeight = defaultWireHeight;
            bool isRoundWire = defaultIsRoundWire;
            if (turn.get_dimensions() && turn.get_dimensions()->size() >= 2) {
                wireWidth = (*turn.get_dimensions())[0];
                wireHeight = (*turn.get_dimensions())[1];
            }
            // Check turn's cross-sectional shape
            if (turn.get_cross_sectional_shape().has_value()) {
                isRoundWire = (turn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND);
            }
            
            double turnLength = turn.get_length();
            // For concentric cores, use toroidal quadrant initialization with 0 angle
            // This gives proper RADIAL_INNER/OUTER and TANGENTIAL_LEFT/RIGHT quadrants
            node.initializeToroidalQuadrants(wireWidth, wireHeight, turnLength, 
                                             defaultWireThermalCond, true, 0.0,
                                             defaultWireCoating,
                                             isRoundWire ? TurnCrossSectionalShape::ROUND : TurnCrossSectionalShape::RECTANGULAR);
            _nodes.push_back(node);
        }
    }
}

// ============================================================================
// Thermal Resistance Creation
// ============================================================================

void Temperature::createThermalResistances() {
    _resistances.clear();
    
    if (_isToroidal) {
        createToroidalCoreConnections();
    } else {
        createConcentricCoreConnections();
    }
    
    createBobbinConnections();
    createTurnToTurnConnections();
    if (!_isToroidal) {
        createTurnToBobbinConnections();
    }
    createTurnToInsulationConnections();
    createTurnToSolidConnections();
    createConvectionConnections();
    
    if (THERMAL_DEBUG) {
        std::cout << "Created " << _resistances.size() << " thermal resistances" << std::endl;
    }
}

void Temperature::createToroidalCoreConnections() {
    size_t numSegments = _config.toroidalSegments;
    
    for (size_t i = 0; i < numSegments; ++i) {
        size_t nextIdx = (i + 1) % numSegments;
        
        ThermalResistanceElement r;
        r.nodeFromId = i;
        r.quadrantFrom = ThermalNodeFace::TANGENTIAL_RIGHT;
        r.nodeToId = nextIdx;
        r.quadrantTo = ThermalNodeFace::TANGENTIAL_LEFT;
        r.type = HeatTransferType::CONDUCTION;
        
        auto* q1 = _nodes[i].getQuadrant(ThermalNodeFace::TANGENTIAL_RIGHT);
        auto* q2 = _nodes[nextIdx].getQuadrant(ThermalNodeFace::TANGENTIAL_LEFT);
        
        if (q1 && q2) {
            double contactArea = std::min(q1->surfaceArea, q2->surfaceArea);
            r.resistance = q1->calculateConductionResistance(*q2, contactArea);
            
            // Add coating resistances if coatings exist
            if (q1->coating.has_value()) {
                r.resistance += WireCoatingUtils::calculateCoatingResistance(q1->coating.value(), contactArea);
            }
            if (q2->coating.has_value()) {
                r.resistance += WireCoatingUtils::calculateCoatingResistance(q2->coating.value(), contactArea);
            }
            
            _resistances.push_back(r);
        }
    }
}

void Temperature::createConcentricCoreConnections() {
    // Collect all core node indices by type
    std::vector<size_t> centralColumnIndices;
    std::vector<size_t> lateralColumnIndices;
    std::vector<size_t> topYokeIndices;
    std::vector<size_t> bottomYokeIndices;
    
    for (size_t i = 0; i < _nodes.size(); ++i) {
        switch (_nodes[i].part) {
            case ThermalNodePartType::CORE_CENTRAL_COLUMN:
                centralColumnIndices.push_back(i);
                break;
            case ThermalNodePartType::CORE_LATERAL_COLUMN:
                lateralColumnIndices.push_back(i);
                break;
            case ThermalNodePartType::CORE_TOP_YOKE:
                topYokeIndices.push_back(i);
                break;
            case ThermalNodePartType::CORE_BOTTOM_YOKE:
                bottomYokeIndices.push_back(i);
                break;
            default:
                break;
        }
    }
    
    auto core = _magnetic.get_core();
    auto processedDesc = core.get_processed_description();
    double coreK = _config.coreThermalConductivity;
    
    // Helper to create conduction connection between core nodes with specific quadrants
    auto createCoreConnection = [&](size_t fromIdx, ThermalNodeFace fromFace,
                                    size_t toIdx, ThermalNodeFace toFace,
                                    double contactArea) {
        ThermalResistanceElement r;
        r.nodeFromId = fromIdx;
        r.quadrantFrom = fromFace;
        r.nodeToId = toIdx;
        r.quadrantTo = toFace;
        r.type = HeatTransferType::CONDUCTION;
        
        // Distance is the distance between node centers
        double dx = _nodes[fromIdx].physicalCoordinates[0] - _nodes[toIdx].physicalCoordinates[0];
        double dy = _nodes[fromIdx].physicalCoordinates[1] - _nodes[toIdx].physicalCoordinates[1];
        double dist = std::sqrt(dx*dx + dy*dy);
        
        r.resistance = ThermalResistance::calculateConductionResistance(dist, coreK, contactArea);
        _resistances.push_back(r);
    };
    
    // Helper to get appropriate cross-sectional area based on connection direction
    auto getColumnCrossSection = [&](size_t idx) -> double {
        // For columns, use the area perpendicular to the column axis (x-axis for central, y for lateral)
        return _nodes[idx].dimensions.height * _nodes[idx].dimensions.depth;
    };
    auto getYokeCrossSection = [&](size_t idx) -> double {
        // For yokes, use width * depth
        return _nodes[idx].dimensions.width * _nodes[idx].dimensions.depth;
    };
    
    // Helper to find closest node in a list
    auto findClosestNode = [&](size_t fromIdx, const std::vector<size_t>& candidates) -> size_t {
        size_t closestIdx = candidates[0];
        double minDist = std::numeric_limits<double>::max();
        for (size_t candidateIdx : candidates) {
            if (candidateIdx == fromIdx) continue;
            double dx = _nodes[fromIdx].physicalCoordinates[0] - _nodes[candidateIdx].physicalCoordinates[0];
            double dy = _nodes[fromIdx].physicalCoordinates[1] - _nodes[candidateIdx].physicalCoordinates[1];
            double dist = std::sqrt(dx*dx + dy*dy);
            if (dist < minDist) {
                minDist = dist;
                closestIdx = candidateIdx;
            }
        }
        return closestIdx;
    };
    
    // 1. Connect central column chunks to closest neighbors only (vertical chain)
    if (centralColumnIndices.size() > 1) {
        // Sort by Y position
        std::sort(centralColumnIndices.begin(), centralColumnIndices.end(), 
                  [&](size_t a, size_t b) {
                      return _nodes[a].physicalCoordinates[1] < _nodes[b].physicalCoordinates[1];
                  });
        // Connect each chunk to its vertical neighbor
        // Lower chunk uses TANGENTIAL_LEFT (up), higher chunk uses TANGENTIAL_RIGHT (down)
        for (size_t i = 0; i < centralColumnIndices.size() - 1; ++i) {
            double area = getColumnCrossSection(centralColumnIndices[i]) * 0.5;
            createCoreConnection(centralColumnIndices[i], ThermalNodeFace::TANGENTIAL_LEFT,
                                centralColumnIndices[i+1], ThermalNodeFace::TANGENTIAL_RIGHT,
                                area);
        }
    }
    
    // 2. Connect lateral column chunks to closest neighbors (vertical chain)
    if (lateralColumnIndices.size() > 1) {
        std::sort(lateralColumnIndices.begin(), lateralColumnIndices.end(),
                  [&](size_t a, size_t b) {
                      return _nodes[a].physicalCoordinates[1] < _nodes[b].physicalCoordinates[1];
                  });
        // Connect each chunk to its vertical neighbor
        // Lower chunk uses TANGENTIAL_LEFT (up), higher chunk uses TANGENTIAL_RIGHT (down)
        for (size_t i = 0; i < lateralColumnIndices.size() - 1; ++i) {
            double area = getColumnCrossSection(lateralColumnIndices[i]) * 0.5;
            createCoreConnection(lateralColumnIndices[i], ThermalNodeFace::TANGENTIAL_LEFT,
                                lateralColumnIndices[i+1], ThermalNodeFace::TANGENTIAL_RIGHT,
                                area);
        }
    }
    
    // 3. Connect TOP central column chunk to closest top yoke only
    if (!centralColumnIndices.empty() && !topYokeIndices.empty()) {
        // Find top-most central column chunk (highest Y)
        size_t topColIdx = centralColumnIndices[0];
        double maxY = _nodes[topColIdx].physicalCoordinates[1];
        for (size_t colIdx : centralColumnIndices) {
            double y = _nodes[colIdx].physicalCoordinates[1];
            if (y > maxY) {
                maxY = y;
                topColIdx = colIdx;
            }
        }
        size_t closestYoke = findClosestNode(topColIdx, topYokeIndices);
        double area = getColumnCrossSection(topColIdx);
        // With cardinal mapping: Column TOP ↔ Yoke BOTTOM
        createCoreConnection(topColIdx, ThermalNodeFace::TANGENTIAL_LEFT,  // Column TOP
                            closestYoke, ThermalNodeFace::TANGENTIAL_RIGHT, // Yoke BOTTOM
                            area);
    }
    
    // 4. Connect BOTTOM central column chunk to closest bottom yoke only
    if (!centralColumnIndices.empty() && !bottomYokeIndices.empty()) {
        // Find bottom-most central column chunk (lowest Y)
        size_t botColIdx = centralColumnIndices[0];
        double minY = _nodes[botColIdx].physicalCoordinates[1];
        for (size_t colIdx : centralColumnIndices) {
            double y = _nodes[colIdx].physicalCoordinates[1];
            if (y < minY) {
                minY = y;
                botColIdx = colIdx;
            }
        }
        size_t closestYoke = findClosestNode(botColIdx, bottomYokeIndices);
        double area = getColumnCrossSection(botColIdx);
        // With cardinal mapping: Column BOTTOM ↔ Yoke TOP
        createCoreConnection(botColIdx, ThermalNodeFace::TANGENTIAL_RIGHT, // Column BOTTOM
                            closestYoke, ThermalNodeFace::TANGENTIAL_LEFT,  // Yoke TOP
                            area);
    }
    
    // 5. Connect TOP lateral column chunk to closest top yoke only
    if (!lateralColumnIndices.empty() && !topYokeIndices.empty()) {
        // Find top-most lateral column chunk
        size_t topColIdx = lateralColumnIndices[0];
        double maxY = _nodes[topColIdx].physicalCoordinates[1];
        for (size_t colIdx : lateralColumnIndices) {
            double y = _nodes[colIdx].physicalCoordinates[1];
            if (y > maxY) {
                maxY = y;
                topColIdx = colIdx;
            }
        }
        size_t closestYoke = findClosestNode(topColIdx, topYokeIndices);
        double area = getColumnCrossSection(topColIdx);
        // With cardinal mapping: Column TOP ↔ Yoke BOTTOM
        createCoreConnection(topColIdx, ThermalNodeFace::TANGENTIAL_LEFT,  // Column TOP
                            closestYoke, ThermalNodeFace::TANGENTIAL_RIGHT, // Yoke BOTTOM
                            area);
    }
    
    // 6. Connect BOTTOM lateral column chunk to closest bottom yoke only
    if (!lateralColumnIndices.empty() && !bottomYokeIndices.empty()) {
        // Find bottom-most lateral column chunk
        size_t botColIdx = lateralColumnIndices[0];
        double minY = _nodes[botColIdx].physicalCoordinates[1];
        for (size_t colIdx : lateralColumnIndices) {
            double y = _nodes[colIdx].physicalCoordinates[1];
            if (y < minY) {
                minY = y;
                botColIdx = colIdx;
            }
        }
        size_t closestYoke = findClosestNode(botColIdx, bottomYokeIndices);
        double area = getColumnCrossSection(botColIdx);
        // With cardinal mapping: Column BOTTOM ↔ Yoke TOP
        createCoreConnection(botColIdx, ThermalNodeFace::TANGENTIAL_RIGHT, // Column BOTTOM
                            closestYoke, ThermalNodeFace::TANGENTIAL_LEFT,  // Yoke TOP
                            area);
    }
}

void Temperature::createBobbinConnections() {
    if (!hasBobbinNodes()) {
        return;
    }
    
    // Collect bobbin node indices
    size_t bobbinColumnIdx = std::numeric_limits<size_t>::max();
    size_t bobbinTopYokeIdx = std::numeric_limits<size_t>::max();
    size_t bobbinBottomYokeIdx = std::numeric_limits<size_t>::max();
    
    std::vector<size_t> coreColumnIndices;
    std::vector<size_t> coreTopYokeIndices;
    std::vector<size_t> coreBottomYokeIndices;
    
    for (size_t i = 0; i < _nodes.size(); ++i) {
        switch (_nodes[i].part) {
            case ThermalNodePartType::BOBBIN_CENTRAL_COLUMN:
                bobbinColumnIdx = i;
                break;
            case ThermalNodePartType::BOBBIN_TOP_YOKE:
                bobbinTopYokeIdx = i;
                break;
            case ThermalNodePartType::BOBBIN_BOTTOM_YOKE:
                bobbinBottomYokeIdx = i;
                break;
            case ThermalNodePartType::CORE_CENTRAL_COLUMN:
                coreColumnIndices.push_back(i);
                break;
            case ThermalNodePartType::CORE_TOP_YOKE:
                coreTopYokeIndices.push_back(i);
                break;
            case ThermalNodePartType::CORE_BOTTOM_YOKE:
                coreBottomYokeIndices.push_back(i);
                break;
            default:
                break;
        }
    }
    
    // Helper to create conduction connection with specific quadrants
    auto createConnectionWithQuadrants = [&](size_t fromIdx, ThermalNodeFace fromFace,
                                              size_t toIdx, ThermalNodeFace toFace,
                                              double dist, double area, double k = 0.2) {
        ThermalResistanceElement r;
        r.nodeFromId = fromIdx;
        r.quadrantFrom = fromFace;
        r.nodeToId = toIdx;
        r.quadrantTo = toFace;
        r.type = HeatTransferType::CONDUCTION;
        
        r.resistance = ThermalResistance::calculateConductionResistance(dist, k, area);
        _resistances.push_back(r);
        
        if (THERMAL_DEBUG) {
            std::cout << "Connection: " << _nodes[fromIdx].name << "[" << magic_enum::enum_name(fromFace) 
                      << "] -> " << _nodes[toIdx].name << "[" << magic_enum::enum_name(toFace)
                      << "] (R=" << r.resistance << " K/W)" << std::endl;
        }
    };
    
    // Helper to calculate contact area between bobbin and core surfaces
    auto calculateBobbinColumnContactArea = [&](size_t bobbinIdx, size_t coreIdx) -> double {
        double bobbinHeight = _nodes[bobbinIdx].dimensions.height;
        double bobbinDepth = _nodes[bobbinIdx].dimensions.depth;
        double coreDepth = _nodes[coreIdx].dimensions.depth;
        return bobbinHeight * std::min(bobbinDepth, coreDepth);
    };
    
    auto calculateBobbinYokeContactArea = [&](size_t bobbinIdx, size_t coreIdx) -> double {
        double wallThickness = _nodes[bobbinIdx].dimensions.height;
        double bobbinDepth = _nodes[bobbinIdx].dimensions.depth;
        double coreDepth = _nodes[coreIdx].dimensions.depth;
        return wallThickness * std::min(bobbinDepth, coreDepth);
    };
    
    // 1. Bobbin column connects to ALL core column nodes
    // Bobbin uses RADIAL_INNER (facing toward core), Core uses RADIAL_OUTER (facing toward bobbin)
    if (bobbinColumnIdx != std::numeric_limits<size_t>::max()) {
        for (size_t coreIdx : coreColumnIndices) {
            double dx = _nodes[bobbinColumnIdx].physicalCoordinates[0] - _nodes[coreIdx].physicalCoordinates[0];
            double dy = _nodes[bobbinColumnIdx].physicalCoordinates[1] - _nodes[coreIdx].physicalCoordinates[1];
            double dist = std::sqrt(dx*dx + dy*dy);
            double area = calculateBobbinColumnContactArea(bobbinColumnIdx, coreIdx);
            createConnectionWithQuadrants(bobbinColumnIdx, ThermalNodeFace::RADIAL_INNER,
                                          coreIdx, ThermalNodeFace::RADIAL_OUTER,
                                          dist, area, 0.2);
        }
    }
    
    // 2. Bobbin top yoke connects ONLY to closest top core yoke
    // Bobbin uses RADIAL_INNER (facing down/toward core), Core uses TANGENTIAL_LEFT (facing up/toward bobbin)
    if (bobbinTopYokeIdx != std::numeric_limits<size_t>::max() && !coreTopYokeIndices.empty()) {
        size_t closestIdx = coreTopYokeIndices[0];
        double minDist = std::numeric_limits<double>::max();
        for (size_t coreIdx : coreTopYokeIndices) {
            double dx = _nodes[bobbinTopYokeIdx].physicalCoordinates[0] - _nodes[coreIdx].physicalCoordinates[0];
            double dy = _nodes[bobbinTopYokeIdx].physicalCoordinates[1] - _nodes[coreIdx].physicalCoordinates[1];
            double dist = std::sqrt(dx*dx + dy*dy);
            if (dist < minDist) {
                minDist = dist;
                closestIdx = coreIdx;
            }
        }
        double area = calculateBobbinYokeContactArea(bobbinTopYokeIdx, closestIdx);
        // With cardinal mapping: Bobbin yoke TOP ↔ Core yoke BOTTOM
        createConnectionWithQuadrants(bobbinTopYokeIdx, ThermalNodeFace::TANGENTIAL_LEFT,  // Bobbin TOP
                                      closestIdx, ThermalNodeFace::TANGENTIAL_RIGHT,       // Core yoke BOTTOM
                                      minDist, area, 0.2);
    }
    
    // 3. Bobbin bottom yoke connects ONLY to closest bottom core yoke
    // With cardinal mapping: Bobbin yoke BOTTOM ↔ Core yoke TOP
    if (bobbinBottomYokeIdx != std::numeric_limits<size_t>::max() && !coreBottomYokeIndices.empty()) {
        size_t closestIdx = coreBottomYokeIndices[0];
        double minDist = std::numeric_limits<double>::max();
        for (size_t coreIdx : coreBottomYokeIndices) {
            double dx = _nodes[bobbinBottomYokeIdx].physicalCoordinates[0] - _nodes[coreIdx].physicalCoordinates[0];
            double dy = _nodes[bobbinBottomYokeIdx].physicalCoordinates[1] - _nodes[coreIdx].physicalCoordinates[1];
            double dist = std::sqrt(dx*dx + dy*dy);
            if (dist < minDist) {
                minDist = dist;
                closestIdx = coreIdx;
            }
        }
        double area = calculateBobbinYokeContactArea(bobbinBottomYokeIdx, closestIdx);
        // With cardinal mapping: Bobbin yoke BOTTOM ↔ Core yoke TOP
        createConnectionWithQuadrants(bobbinBottomYokeIdx, ThermalNodeFace::TANGENTIAL_RIGHT, // Bobbin BOTTOM
                                      closestIdx, ThermalNodeFace::TANGENTIAL_LEFT,           // Core yoke TOP
                                      minDist, area, 0.2);
    }
    
    // 4. Bobbin column connects to bobbin yokes (internal bobbin conduction)
    // With cardinal mapping:
    // - Column TOP (TANGENTIAL_LEFT) connects to Top yoke LEFT (RADIAL_INNER)
    // - Column BOTTOM (TANGENTIAL_RIGHT) connects to Bottom yoke LEFT (RADIAL_INNER)
    if (bobbinColumnIdx != std::numeric_limits<size_t>::max()) {
        if (bobbinTopYokeIdx != std::numeric_limits<size_t>::max()) {
            double dx = _nodes[bobbinColumnIdx].physicalCoordinates[0] - _nodes[bobbinTopYokeIdx].physicalCoordinates[0];
            double dy = _nodes[bobbinColumnIdx].physicalCoordinates[1] - _nodes[bobbinTopYokeIdx].physicalCoordinates[1];
            double dist = std::sqrt(dx*dx + dy*dy);
            double area = _nodes[bobbinColumnIdx].dimensions.depth * _nodes[bobbinColumnIdx].dimensions.width * 0.5;
            // Column TOP ↔ Top yoke LEFT
            createConnectionWithQuadrants(bobbinColumnIdx, ThermalNodeFace::TANGENTIAL_LEFT,  // Column TOP
                                          bobbinTopYokeIdx, ThermalNodeFace::RADIAL_INNER,     // Yoke LEFT
                                          dist, area, 0.2);
        }
        if (bobbinBottomYokeIdx != std::numeric_limits<size_t>::max()) {
            double dx = _nodes[bobbinColumnIdx].physicalCoordinates[0] - _nodes[bobbinBottomYokeIdx].physicalCoordinates[0];
            double dy = _nodes[bobbinColumnIdx].physicalCoordinates[1] - _nodes[bobbinBottomYokeIdx].physicalCoordinates[1];
            double dist = std::sqrt(dx*dx + dy*dy);
            double area = _nodes[bobbinColumnIdx].dimensions.depth * _nodes[bobbinColumnIdx].dimensions.width * 0.5;
            // Column BOTTOM ↔ Bottom yoke LEFT
            createConnectionWithQuadrants(bobbinColumnIdx, ThermalNodeFace::TANGENTIAL_RIGHT, // Column BOTTOM
                                          bobbinBottomYokeIdx, ThermalNodeFace::RADIAL_INNER,  // Yoke LEFT
                                          dist, area, 0.2);
        }
    }
    
    // 5. Bobbin yokes connect to nearby turns on their tangential faces
    // Top yoke: TANGENTIAL_RIGHT faces toward turns (downward)
    // Bottom yoke: TANGENTIAL_LEFT faces toward turns (upward)
    std::vector<size_t> turnNodeIndices;
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part == ThermalNodePartType::TURN) {
            turnNodeIndices.push_back(i);
        }
    }
    createBobbinYokeToTurnConnections(bobbinTopYokeIdx, bobbinBottomYokeIdx, turnNodeIndices);
}

// ============================================================================
// Connection Creation - Per Quadrant Logic
// ============================================================================

void Temperature::createTurnToTurnConnections() {
    double minConductionDist = getMinimumDistanceForConduction();
    
    std::vector<size_t> turnNodeIndices;
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part == ThermalNodePartType::TURN) {
            turnNodeIndices.push_back(i);
        }
    }
    
    if (_isToroidal) {
        createToroidalTurnToTurnConnections(turnNodeIndices, minConductionDist);
    } else {
        createConcentricTurnToTurnConnections(turnNodeIndices, minConductionDist);
    }
}

void Temperature::createConcentricTurnToTurnConnections(const std::vector<size_t>& turnNodeIndices, double minConductionDist) {
    // Pure geometry-based connections - no layer or consecutiveness logic
    for (size_t i = 0; i < turnNodeIndices.size(); i++) {
        size_t node1Idx = turnNodeIndices[i];
        const auto& node1 = _nodes[node1Idx];
        
        for (size_t j = i + 1; j < turnNodeIndices.size(); j++) {
            size_t node2Idx = turnNodeIndices[j];
            const auto& node2 = _nodes[node2Idx];
            
            double dx = node1.physicalCoordinates[0] - node2.physicalCoordinates[0];
            double dy = node1.physicalCoordinates[1] - node2.physicalCoordinates[1];
            double centerDistance = std::sqrt(dx*dx + dy*dy);
            
            // Connection condition: surface distance < min(min_outer_dim1, min_outer_dim2) / 4
            double minOuterDim1 = std::min(node1.dimensions.width, node1.dimensions.height);
            double minOuterDim2 = std::min(node2.dimensions.width, node2.dimensions.height);
            double thresholdDist = std::min(minOuterDim1, minOuterDim2) / 4.0;
            
            // For rectangular wires, calculate surface distance based on relative orientation
            // The surface distance is center distance minus the projected size in the connection direction
            double surfaceDistance;
            if (centerDistance < 1e-9) {
                surfaceDistance = 0;
            } else {
                // For each node, use half the min dimension as the conservative extent estimate
                // This works for both rectangular and round wires
                double extent1 = minOuterDim1 / 2.0;
                double extent2 = minOuterDim2 / 2.0;
                
                surfaceDistance = centerDistance - extent1 - extent2;
            }
            
            if (surfaceDistance > thresholdDist) continue;
            
            // Find which quadrants are facing each other
            ThermalNodeFace face1, face2;
            
            // Determine primary direction of connection
            bool moreHorizontal = std::abs(dx) > std::abs(dy);
            
            if (moreHorizontal) {
                // Horizontal connection - use RADIAL faces
                if (dx > 0) {
                    // node1 is to the right of node2
                    face1 = ThermalNodeFace::RADIAL_INNER;  // LEFT face of node1
                    face2 = ThermalNodeFace::RADIAL_OUTER;  // RIGHT face of node2
                } else {
                    // node1 is to the left of node2
                    face1 = ThermalNodeFace::RADIAL_OUTER;  // RIGHT face of node1
                    face2 = ThermalNodeFace::RADIAL_INNER;  // LEFT face of node2
                }
            } else {
                // Vertical connection - use TANGENTIAL faces
                if (dy > 0) {
                    // node1 is above node2
                    face1 = ThermalNodeFace::TANGENTIAL_RIGHT; // BOTTOM face of node1
                    face2 = ThermalNodeFace::TANGENTIAL_LEFT;  // TOP face of node2
                } else {
                    // node1 is below node2
                    face1 = ThermalNodeFace::TANGENTIAL_LEFT;  // TOP face of node1
                    face2 = ThermalNodeFace::TANGENTIAL_RIGHT; // BOTTOM face of node2
                }
            }
            
            ThermalResistanceElement r;
            r.nodeFromId = node1Idx;
            r.quadrantFrom = face1;
            r.nodeToId = node2Idx;
            r.quadrantTo = face2;
            r.type = HeatTransferType::CONDUCTION;
            
            auto* q1 = _nodes[node1Idx].getQuadrant(face1);
            auto* q2 = _nodes[node2Idx].getQuadrant(face2);
            
            if (q1 && q2) {
                double contactArea = std::min(q1->surfaceArea, q2->surfaceArea);
                int turn1Idx = _nodes[node1Idx].turnIndex.value_or(0);
                int turn2Idx = _nodes[node2Idx].turnIndex.value_or(0);
                
                // Check if there's a solid insulation layer between these turns
                // If so, skip direct connection - turns will connect through insulation layer node
                auto coil = _magnetic.get_coil();
                auto turnsDescription = coil.get_turns_description();
                bool hasSolidInsulationBetween = false;
                
                if (turnsDescription && turn1Idx >= 0 && turn1Idx < static_cast<int>(turnsDescription->size()) &&
                    turn2Idx >= 0 && turn2Idx < static_cast<int>(turnsDescription->size())) {
                    try {
                        auto& turn1 = (*turnsDescription)[turn1Idx];
                        auto& turn2 = (*turnsDescription)[turn2Idx];
                        auto layersBetween = StrayCapacitance::get_insulation_layers_between_two_turns(turn1, turn2, coil);
                        hasSolidInsulationBetween = !layersBetween.empty();
                    } catch (...) {
                        // If we can't determine insulation layers (e.g., missing layer description),
                        // assume no solid insulation and create direct connection
                        hasSolidInsulationBetween = false;
                    }
                }
                
                // Also check geometrically if an insulation layer node is between these turns
                // This is important for concentric cores where insulation layers are explicit nodes
                bool insulationNodeBetween = false;
                for (size_t k = 0; k < _nodes.size(); k++) {
                    if (_nodes[k].part != ThermalNodePartType::INSULATION_LAYER) continue;
                    
                    const auto& insulationNode = _nodes[k];
                    double insX = insulationNode.physicalCoordinates[0];
                    double insY = insulationNode.physicalCoordinates[1];
                    double insWidth = insulationNode.dimensions.width;
                    double insHeight = insulationNode.dimensions.height;
                    
                    // Check if insulation layer is between the two turns
                    // For horizontal (radial) connections: insulation X should be between turn1 X and turn2 X
                    // For vertical (tangential) connections: insulation Y should be between turn1 Y and turn2 Y
                    double turn1X = node1.physicalCoordinates[0];
                    double turn1Y = node1.physicalCoordinates[1];
                    double turn2X = node2.physicalCoordinates[0];
                    double turn2Y = node2.physicalCoordinates[1];
                    
                    if (moreHorizontal) {
                        // Horizontal connection - check if insulation is between turns in X direction
                        double minTurnX = std::min(turn1X, turn2X);
                        double maxTurnX = std::max(turn1X, turn2X);
                        double insLeft = insX - insWidth / 2;
                        double insRight = insX + insWidth / 2;
                        
                        // Check if insulation overlaps the line between turns
                        bool xOverlap = (insLeft < maxTurnX) && (insRight > minTurnX);
                        bool yOverlap = std::abs(insY - turn1Y) < (insHeight / 2 + node1.dimensions.height / 2);
                        
                        if (xOverlap && yOverlap) {
                            insulationNodeBetween = true;
                            break;
                        }
                    } else {
                        // Vertical connection - check if insulation is between turns in Y direction
                        double minTurnY = std::min(turn1Y, turn2Y);
                        double maxTurnY = std::max(turn1Y, turn2Y);
                        double insBottom = insY - insHeight / 2;
                        double insTop = insY + insHeight / 2;
                        
                        // Check if insulation overlaps the line between turns
                        bool yOverlap = (insBottom < maxTurnY) && (insTop > minTurnY);
                        bool xOverlap = std::abs(insX - turn1X) < (insWidth / 2 + node1.dimensions.width / 2);
                        
                        if (xOverlap && yOverlap) {
                            insulationNodeBetween = true;
                            break;
                        }
                    }
                }
                
                // Only create direct turn-to-turn connection if there's NO solid insulation layer
                // Solid insulation layers are now modeled as separate thermal nodes
                if (!hasSolidInsulationBetween && !insulationNodeBetween) {
                    r.resistance = getInsulationLayerThermalResistance(turn1Idx, turn2Idx, contactArea);
                    
                    if (q1->coating.has_value()) {
                        r.resistance += WireCoatingUtils::calculateCoatingResistance(q1->coating.value(), contactArea);
                    }
                    if (q2->coating.has_value()) {
                        r.resistance += WireCoatingUtils::calculateCoatingResistance(q2->coating.value(), contactArea);
                    }
                    
                    _resistances.push_back(r);
                }
            }
        }
    }
}

void Temperature::createToroidalTurnToTurnConnections(const std::vector<size_t>& turnNodeIndices, double minConductionDist) {
    // Pure geometry-based connections using turn rotation to determine facing quadrants
    for (size_t i = 0; i < turnNodeIndices.size(); i++) {
        size_t node1Idx = turnNodeIndices[i];
        const auto& node1 = _nodes[node1Idx];
        
        for (size_t j = i + 1; j < turnNodeIndices.size(); j++) {
            size_t node2Idx = turnNodeIndices[j];
            const auto& node2 = _nodes[node2Idx];
            
            double dx = node1.physicalCoordinates[0] - node2.physicalCoordinates[0];
            double dy = node1.physicalCoordinates[1] - node2.physicalCoordinates[1];
            double centerDistance = std::sqrt(dx*dx + dy*dy);
            
            // Get minimum dimensions from both wires (used for both threshold and surface distance)
            double minDim1 = std::min(node1.dimensions.width, node1.dimensions.height);
            double minDim2 = std::min(node2.dimensions.width, node2.dimensions.height);
            double thresholdDist = std::min(minDim1, minDim2) / 4.0;
            
            // Calculate surface distance: center distance minus half of min dimensions
            // This is more accurate than using max dimension for rectangular wires
            double extent1 = minDim1 / 2.0;
            double extent2 = minDim2 / 2.0;
            double surfaceDistance = centerDistance - extent1 - extent2;
            
            if (surfaceDistance > thresholdDist) continue;
            
            // Get the angle of each turn (rotation around the toroid)
            double angle1 = std::atan2(node1.physicalCoordinates[1], node1.physicalCoordinates[0]);
            double angle2 = std::atan2(node2.physicalCoordinates[1], node2.physicalCoordinates[0]);
            
            // Direction from turn1 to turn2
            double dirX = -dx / centerDistance;  // unit vector from 1 to 2
            double dirY = -dy / centerDistance;
            
            // For each turn, determine which quadrant faces the other turn
            // Quadrant directions (in global coordinates):
            // RADIAL_OUTER:  direction = angle
            // RADIAL_INNER:  direction = angle + pi
            // TANGENTIAL_LEFT:  direction = angle + pi/2 (CCW tangent)
            // TANGENTIAL_RIGHT: direction = angle - pi/2 (CW tangent)
            
            auto getQuadrantDirection = [](double turnAngle, ThermalNodeFace face) -> std::pair<double, double> {
                switch (face) {
                    case ThermalNodeFace::RADIAL_OUTER:
                        return {std::cos(turnAngle), std::sin(turnAngle)};
                    case ThermalNodeFace::RADIAL_INNER:
                        return {std::cos(turnAngle + std::numbers::pi), std::sin(turnAngle + std::numbers::pi)};
                    case ThermalNodeFace::TANGENTIAL_LEFT:
                        return {std::cos(turnAngle + std::numbers::pi/2), std::sin(turnAngle + std::numbers::pi/2)};
                    case ThermalNodeFace::TANGENTIAL_RIGHT:
                        return {std::cos(turnAngle - std::numbers::pi/2), std::sin(turnAngle - std::numbers::pi/2)};
                    default:
                        return {0, 0};
                }
            };
            
            // Find which quadrant of turn1 best faces turn2
            // (quadrant whose direction has highest dot product with dirX,dirY)
            std::array<ThermalNodeFace, 4> faces = {
                ThermalNodeFace::RADIAL_OUTER, ThermalNodeFace::RADIAL_INNER,
                ThermalNodeFace::TANGENTIAL_LEFT, ThermalNodeFace::TANGENTIAL_RIGHT
            };
            
            double bestDot1 = -1.0;
            ThermalNodeFace face1 = ThermalNodeFace::NONE;
            for (auto face : faces) {
                auto [qx, qy] = getQuadrantDirection(angle1, face);
                double dot = qx * dirX + qy * dirY;
                if (dot > bestDot1) {
                    bestDot1 = dot;
                    face1 = face;
                }
            }
            
            // Find which quadrant of turn2 best faces turn1
            // (direction opposite to dirX, dirY)
            double bestDot2 = -1.0;
            ThermalNodeFace face2 = ThermalNodeFace::NONE;
            for (auto face : faces) {
                auto [qx, qy] = getQuadrantDirection(angle2, face);
                double dot = qx * (-dirX) + qy * (-dirY);  // opposite direction
                if (dot > bestDot2) {
                    bestDot2 = dot;
                    face2 = face;
                }
            }
            
            // Only connect if we found valid facing quadrants
            if (face1 == ThermalNodeFace::NONE || face2 == ThermalNodeFace::NONE) continue;
            
            ThermalResistanceElement r;
            r.nodeFromId = node1Idx;
            r.quadrantFrom = face1;
            r.nodeToId = node2Idx;
            r.quadrantTo = face2;
            r.type = HeatTransferType::CONDUCTION;
            
            auto* q1 = _nodes[node1Idx].getQuadrant(face1);
            auto* q2 = _nodes[node2Idx].getQuadrant(face2);
            
            if (q1 && q2) {
                double contactArea = std::min(q1->surfaceArea, q2->surfaceArea);
                int turn1Idx = _nodes[node1Idx].turnIndex.value_or(0);
                int turn2Idx = _nodes[node2Idx].turnIndex.value_or(0);
                
                // Check if there's a solid insulation layer between these turns
                // If so, skip direct connection - turns will connect through insulation layer node
                auto coil = _magnetic.get_coil();
                auto turnsDescription = coil.get_turns_description();
                bool hasSolidInsulationBetween = false;
                
                if (turnsDescription && turn1Idx >= 0 && turn1Idx < static_cast<int>(turnsDescription->size()) &&
                    turn2Idx >= 0 && turn2Idx < static_cast<int>(turnsDescription->size())) {
                    try {
                        auto& turn1 = (*turnsDescription)[turn1Idx];
                        auto& turn2 = (*turnsDescription)[turn2Idx];
                        auto layersBetween = StrayCapacitance::get_insulation_layers_between_two_turns(turn1, turn2, coil);
                        hasSolidInsulationBetween = !layersBetween.empty();
                    } catch (...) {
                        // If we can't determine insulation layers (e.g., missing layer description),
                        // assume no solid insulation and create direct connection
                        hasSolidInsulationBetween = false;
                    }
                }
                
                // Only create direct turn-to-turn connection if there's NO solid insulation layer
                // Solid insulation layers are now modeled as separate thermal nodes
                if (!hasSolidInsulationBetween) {
                    double baseResistance = getInsulationLayerThermalResistance(turn1Idx, turn2Idx, contactArea);
                    
                    if (q1->coating.has_value()) {
                        baseResistance += WireCoatingUtils::calculateCoatingResistance(q1->coating.value(), contactArea);
                    }
                    if (q2->coating.has_value()) {
                        baseResistance += WireCoatingUtils::calculateCoatingResistance(q2->coating.value(), contactArea);
                    }
                    
                    r.resistance = baseResistance;
                    
                    if (_config.useInterTurnInsulation && _config.interTurnInsulationThickness > 0) {
                        r.addInsulationLayer(
                            _config.interTurnInsulationThickness,
                            _config.interTurnInsulationConductivity,
                            "inter_turn_insulation",
                            "Additional insulation between turns from config"
                        );
                        r.resistance += r.calculateTotalInsulationResistance(contactArea);
                    }
                    
                    _resistances.push_back(r);
                }
            }
        }
    }
}

void Temperature::createBobbinYokeToTurnConnections(size_t bobbinTopYokeIdx, size_t bobbinBottomYokeIdx, 
                                                     const std::vector<size_t>& turnNodeIndices) {
    auto findNearbyTurns = [&](size_t yokeIdx, ThermalNodeFace yokeFace) -> std::vector<std::pair<size_t, double>> {
        std::vector<std::pair<size_t, double>> nearbyTurns;
        if (yokeIdx == std::numeric_limits<size_t>::max()) return nearbyTurns;
        
        const auto& yokeNode = _nodes[yokeIdx];
        double yokeX = yokeNode.physicalCoordinates[0];
        double yokeY = yokeNode.physicalCoordinates[1];
        
        for (size_t turnIdx : turnNodeIndices) {
            const auto& turnNode = _nodes[turnIdx];
            double turnX = turnNode.physicalCoordinates[0];
            double turnY = turnNode.physicalCoordinates[1];
            double turnWidth = turnNode.dimensions.width;
            double turnHeight = turnNode.dimensions.height;
            
            // Threshold distance for turn-to-bobbin contact based on this turn's dimensions
            double contactThreshold = std::max(turnWidth, turnHeight) / 4.0;
            
            // Calculate distance from yoke to turn
            double dx = turnX - yokeX;
            double dy = turnY - yokeY;
            double dist = std::sqrt(dx*dx + dy*dy);
            
            // Check if turn is within contact threshold
            if (dist < contactThreshold + std::max(turnWidth, turnHeight)) {
                // Check if turn is on the correct side of the yoke
                bool correctSide = false;
                if (yokeFace == ThermalNodeFace::TANGENTIAL_RIGHT) {
                    // Top yoke: turn should be below yoke (dy < 0)
                    correctSide = dy < 0;
                } else if (yokeFace == ThermalNodeFace::TANGENTIAL_LEFT) {
                    // Bottom yoke: turn should be above yoke (dy > 0)
                    correctSide = dy > 0;
                }
                
                if (correctSide) {
                    nearbyTurns.push_back({turnIdx, dist});
                }
            }
        }
        
        return nearbyTurns;
    };
    
    // Helper to create bobbin-to-turn connection
    auto createBobbinToTurnConnection = [&](size_t bobbinIdx, ThermalNodeFace bobbinFace, 
                                            size_t turnIdx, double dist) {
        const auto& turnNode = _nodes[turnIdx];
        double turnWidth = turnNode.dimensions.width;
        double turnHeight = turnNode.dimensions.height;
        
        ThermalResistanceElement r;
        r.nodeFromId = bobbinIdx;
        r.quadrantFrom = bobbinFace;
        r.nodeToId = turnIdx;
        r.quadrantTo = (bobbinFace == ThermalNodeFace::TANGENTIAL_RIGHT) ? 
                       ThermalNodeFace::TANGENTIAL_LEFT : ThermalNodeFace::TANGENTIAL_RIGHT;
        r.type = HeatTransferType::CONDUCTION;
        
        // Calculate resistance through air/bobbin gap - use turn's actual dimensions
        double contactArea = turnWidth * turnHeight * 0.5;
        double k_air = 0.025;  // W/m·K
        r.resistance = dist / (k_air * contactArea);
        
        _resistances.push_back(r);
    };
    
    // Connect top bobbin yoke to nearby turns on TANGENTIAL_RIGHT face
    if (bobbinTopYokeIdx != std::numeric_limits<size_t>::max()) {
        auto nearbyTurns = findNearbyTurns(bobbinTopYokeIdx, ThermalNodeFace::TANGENTIAL_RIGHT);
        
        if (!nearbyTurns.empty()) {
            // Connect to all nearby turns (like toroidal core chunk approach)
            double totalContactArea = 0;
            for (const auto& [turnIdx, dist] : nearbyTurns) {
                // Estimate contact area for this turn using its stored dimensions
                const auto& turnNode = _nodes[turnIdx];
                double turnContactArea = turnNode.dimensions.width * turnNode.dimensions.height * 0.25;
                totalContactArea += turnContactArea;
            }
            
            // Create connections with proportional area
            for (const auto& [turnIdx, dist] : nearbyTurns) {
                createBobbinToTurnConnection(bobbinTopYokeIdx, ThermalNodeFace::TANGENTIAL_RIGHT, 
                                            turnIdx, dist);
            }
        }
        // If no nearby turns, convection to air will be handled by createConvectionConnections
    }
    
    // Connect bottom bobbin yoke to nearby turns on TANGENTIAL_LEFT face
    if (bobbinBottomYokeIdx != std::numeric_limits<size_t>::max()) {
        auto nearbyTurns = findNearbyTurns(bobbinBottomYokeIdx, ThermalNodeFace::TANGENTIAL_LEFT);
        
        if (!nearbyTurns.empty()) {
            double totalContactArea = 0;
            for (const auto& [turnIdx, dist] : nearbyTurns) {
                const auto& turnNode = _nodes[turnIdx];
                double turnContactArea = turnNode.dimensions.width * turnNode.dimensions.height * 0.25;
                totalContactArea += turnContactArea;
            }
            
            for (const auto& [turnIdx, dist] : nearbyTurns) {
                createBobbinToTurnConnection(bobbinBottomYokeIdx, ThermalNodeFace::TANGENTIAL_LEFT, 
                                            turnIdx, dist);
            }
        }
    }
}

void Temperature::createTurnToBobbinConnections() {
    // Find bobbin nodes
    std::vector<size_t> bobbinNodeIndices;
    std::vector<size_t> turnNodeIndices;
    
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part == ThermalNodePartType::BOBBIN_CENTRAL_COLUMN ||
            _nodes[i].part == ThermalNodePartType::BOBBIN_TOP_YOKE ||
            _nodes[i].part == ThermalNodePartType::BOBBIN_BOTTOM_YOKE) {
            bobbinNodeIndices.push_back(i);
        } else if (_nodes[i].part == ThermalNodePartType::TURN) {
            turnNodeIndices.push_back(i);
        }
    }
    
    if (bobbinNodeIndices.empty()) {
        return;  // No bobbin, turns will connect directly to core
    }
    
    double minConductionDist = getMinimumDistanceForConduction();
    
    // For each turn, create conduction connections for ALL quadrants touching any bobbin
    for (size_t turnIdx : turnNodeIndices) {
        auto& turnNode = _nodes[turnIdx];
        
        double turnX = turnNode.physicalCoordinates[0];
        double turnY = turnNode.physicalCoordinates[1];
        
        // Track which (face, bobbin) pairs we've already connected to avoid duplicates
        std::set<std::pair<ThermalNodeFace, size_t>> connectedPairs;
        
        for (const auto& quadrant : turnNode.quadrants) {
            if (quadrant.face == ThermalNodeFace::NONE) continue;
            
            // For concentric turns, calculate limit position based on horizontal/vertical offset
            // (not angled like toroidal turns) - use turn node's stored dimensions
            double turnWidth = turnNode.dimensions.width;
            double turnHeight = turnNode.dimensions.height;
            double limitX, limitY;
            switch (quadrant.face) {
                case ThermalNodeFace::RADIAL_INNER:   // Left face (-X)
                    limitX = turnX - turnWidth / 2.0;
                    limitY = turnY;
                    break;
                case ThermalNodeFace::RADIAL_OUTER:   // Right face (+X)
                    limitX = turnX + turnWidth / 2.0;
                    limitY = turnY;
                    break;
                case ThermalNodeFace::TANGENTIAL_LEFT:  // Top face (+Y)
                    limitX = turnX;
                    limitY = turnY + turnHeight / 2.0;
                    break;
                case ThermalNodeFace::TANGENTIAL_RIGHT: // Bottom face (-Y)
                    limitX = turnX;
                    limitY = turnY - turnHeight / 2.0;
                    break;
                default:
                    limitX = quadrant.limitCoordinates[0];
                    limitY = quadrant.limitCoordinates[1];
                    break;
            }
            
            // Check all bobbin nodes for this quadrant
            for (size_t bobbinIdx : bobbinNodeIndices) {
                // Skip if we already connected this face to this bobbin
                if (connectedPairs.count({quadrant.face, bobbinIdx})) continue;
                
                const auto& bobbinNode = _nodes[bobbinIdx];
                bool isCentralColumn = (bobbinNode.part == ThermalNodePartType::BOBBIN_CENTRAL_COLUMN);
                bool isTopYoke = (bobbinNode.part == ThermalNodePartType::BOBBIN_TOP_YOKE);
                bool isBottomYoke = (bobbinNode.part == ThermalNodePartType::BOBBIN_BOTTOM_YOKE);
                
                // For concentric cores, only connect specific faces to specific bobbin parts:
                // - RADIAL_INNER (LEFT face) -> Central column only
                // - TANGENTIAL_LEFT (TOP face) -> Top yoke only (if turn is at top)
                // - TANGENTIAL_RIGHT (BOTTOM face) -> Bottom yoke only (if turn is at bottom)
                // - RADIAL_OUTER (RIGHT face) -> No bobbin connection
                if (quadrant.face == ThermalNodeFace::RADIAL_INNER && !isCentralColumn) continue;
                if (quadrant.face == ThermalNodeFace::TANGENTIAL_LEFT && !isTopYoke) continue;
                if (quadrant.face == ThermalNodeFace::TANGENTIAL_RIGHT && !isBottomYoke) continue;
                if (quadrant.face == ThermalNodeFace::RADIAL_OUTER) continue;
                
                double bobbinWidth = bobbinNode.dimensions.width;
                double bobbinHeight = bobbinNode.dimensions.height;
                
                // Calculate distance to bobbin surface (not center)
                double bobbinLeft = bobbinNode.physicalCoordinates[0] - bobbinWidth / 2;
                double bobbinRight = bobbinNode.physicalCoordinates[0] + bobbinWidth / 2;
                double bobbinBottom = bobbinNode.physicalCoordinates[1] - bobbinHeight / 2;
                double bobbinTop = bobbinNode.physicalCoordinates[1] + bobbinHeight / 2;
                
                // Closest point on bobbin to the turn's limit point
                double closestX = std::max(bobbinLeft, std::min(limitX, bobbinRight));
                double closestY = std::max(bobbinBottom, std::min(limitY, bobbinTop));
                
                double dx = limitX - closestX;
                double dy = limitY - closestY;
                double dist = std::sqrt(dx*dx + dy*dy);
                
                // Check if within the expanded bounding box for conduction
                bool withinX = std::abs(limitX - bobbinNode.physicalCoordinates[0]) < (bobbinWidth / 2 + minConductionDist);
                bool withinY = std::abs(limitY - bobbinNode.physicalCoordinates[1]) < (bobbinHeight / 2 + minConductionDist);
                
                // Create conduction connection if within range
                if (withinX && withinY && dist < (minConductionDist * 2)) {
                    ThermalResistanceElement r;
                    r.nodeFromId = turnIdx;
                    r.quadrantFrom = quadrant.face;
                    r.nodeToId = bobbinIdx;
                    r.quadrantTo = ThermalNodeFace::NONE;
                    r.type = HeatTransferType::CONDUCTION;
                    
                    auto* q = turnNode.getQuadrant(quadrant.face);
                    double contactArea = q ? q->surfaceArea * 0.5 : (turnWidth * turnHeight * 0.25);
                    double distance = std::max(dist, 1e-6);
                    
                    double k_air = 0.025;  // W/m·K
                    r.resistance = distance / (k_air * contactArea);
                    
                    _resistances.push_back(r);
                    connectedPairs.insert({quadrant.face, bobbinIdx});
                    
                    if (THERMAL_DEBUG) {
                        std::cout << "Turn " << turnNode.name << " face=" << static_cast<int>(quadrant.face) 
                                  << " -> Bobbin " << bobbinNode.name 
                                  << " (distance=" << dist * 1000 << "mm)" << std::endl;
                    }
                }
            }
        }
    }
}

void Temperature::createTurnToInsulationConnections() {
    // Find all insulation layer nodes
    std::vector<size_t> insulationNodeIndices;
    std::vector<size_t> turnNodeIndices;
    
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part == ThermalNodePartType::INSULATION_LAYER) {
            insulationNodeIndices.push_back(i);
        } else if (_nodes[i].part == ThermalNodePartType::TURN) {
            turnNodeIndices.push_back(i);
        }
    }
    
    if (insulationNodeIndices.empty()) {
        return;  // No insulation layers to connect
    }
    
    if (THERMAL_DEBUG) {
        std::cout << "Creating turn-to-insulation connections: " << insulationNodeIndices.size() 
                  << " insulation layers, " << turnNodeIndices.size() << " turns" << std::endl;
    }
    
    double minConductionDist = getMinimumDistanceForConduction();
    
    // For each insulation layer, find adjacent turns on each face
    for (size_t insulationIdx : insulationNodeIndices) {
        const auto& insulationNode = _nodes[insulationIdx];
        
        // Get insulation layer geometry
        double insulationX = insulationNode.physicalCoordinates[0];
        double insulationY = insulationNode.physicalCoordinates[1];
        double insulationWidth = insulationNode.dimensions.width;
        double insulationHeight = insulationNode.dimensions.height;
        
        // Get insulation material properties for resistance calculation
        double insulationK = 0.2;  // Default
        if (insulationNode.quadrants[0].thermalConductivity > 0) {
            insulationK = insulationNode.quadrants[0].thermalConductivity;
        }
        
        // Define the four faces of the insulation layer
        struct FaceCheck {
            ThermalNodeFace face;
            double edgeX;  // X coordinate of the face edge
            double edgeY;  // Y coordinate of the face edge
            bool isVertical;  // true for LEFT/RIGHT faces, false for TOP/BOTTOM
            int direction;  // +1 or -1 indicating outward normal direction
        };
        
        std::vector<FaceCheck> faces = {
            {ThermalNodeFace::RADIAL_OUTER,  insulationX + insulationWidth/2, insulationY, true, +1},   // RIGHT
            {ThermalNodeFace::RADIAL_INNER,  insulationX - insulationWidth/2, insulationY, true, -1},   // LEFT
            {ThermalNodeFace::TANGENTIAL_LEFT,  insulationX, insulationY + insulationHeight/2, false, +1}, // TOP
            {ThermalNodeFace::TANGENTIAL_RIGHT, insulationX, insulationY - insulationHeight/2, false, -1}  // BOTTOM
        };
        
        // Check each face for adjacent turns
        // Connection threshold: surface-to-surface distance must be < 5% of turn diameter
        for (const auto& face : faces) {
            std::vector<std::pair<size_t, double>> adjacentTurns;  // (turnIdx, distance)
            
            for (size_t turnIdx : turnNodeIndices) {
                const auto& turnNode = _nodes[turnIdx];
                double turnX = turnNode.physicalCoordinates[0];
                double turnY = turnNode.physicalCoordinates[1];
                double turnWidth = turnNode.dimensions.width;
                double turnHeight = turnNode.dimensions.height;
                
                // Calculate turn "diameter" (max dimension for round, min for rectangular)
                double turnDiameter = std::max(turnWidth, turnHeight);
                // Threshold: 5% of turn diameter for surface-to-surface connection
                double connectionThreshold = turnDiameter * 0.05;
                
                // Calculate distance from turn surface to insulation face
                double dist;
                if (face.isVertical) {
                    // For vertical faces (LEFT/RIGHT), check X distance
                    double turnEdgeX = (face.direction > 0) ? turnX - turnWidth/2 : turnX + turnWidth/2;
                    dist = std::abs(face.edgeX - turnEdgeX);
                    
                    // Check if turn is within Y range of the face
                    double yOverlap = std::min(
                        std::abs(turnY - (insulationY + insulationHeight/2)),
                        std::abs(turnY - (insulationY - insulationHeight/2))
                    );
                    
                    // Connect only if surface distance < 5% of turn diameter
                    if (dist < connectionThreshold && yOverlap < (insulationHeight/2 + turnHeight/2)) {
                        adjacentTurns.push_back({turnIdx, dist});
                    }
                } else {
                    // For horizontal faces (TOP/BOTTOM), check Y distance
                    double turnEdgeY = (face.direction > 0) ? turnY - turnHeight/2 : turnY + turnHeight/2;
                    dist = std::abs(face.edgeY - turnEdgeY);
                    
                    // Check if turn is within X range of the face
                    double xOverlap = std::min(
                        std::abs(turnX - (insulationX + insulationWidth/2)),
                        std::abs(turnX - (insulationX - insulationWidth/2))
                    );
                    
                    // Connect only if surface distance < 5% of turn diameter
                    if (dist < connectionThreshold && xOverlap < (insulationWidth/2 + turnWidth/2)) {
                        adjacentTurns.push_back({turnIdx, dist});
                    }
                }
            }
            
            // Create connections to all adjacent turns
            for (const auto& [turnIdx, dist] : adjacentTurns) {
                const auto& turnNode = _nodes[turnIdx];
                double turnWidth = turnNode.dimensions.width;
                double turnHeight = turnNode.dimensions.height;
                
                ThermalResistanceElement r;
                r.nodeFromId = insulationIdx;
                r.quadrantFrom = face.face;
                r.nodeToId = turnIdx;
                
                // Map insulation face to corresponding turn face
                // Turn face should be opposite of insulation face (they touch)
                switch (face.face) {
                    case ThermalNodeFace::RADIAL_OUTER:
                        r.quadrantTo = ThermalNodeFace::RADIAL_INNER;  // Turn's LEFT touches insulation's RIGHT
                        break;
                    case ThermalNodeFace::RADIAL_INNER:
                        r.quadrantTo = ThermalNodeFace::RADIAL_OUTER;  // Turn's RIGHT touches insulation's LEFT
                        break;
                    case ThermalNodeFace::TANGENTIAL_LEFT:
                        r.quadrantTo = ThermalNodeFace::TANGENTIAL_RIGHT;  // Turn's BOTTOM touches insulation's TOP
                        break;
                    case ThermalNodeFace::TANGENTIAL_RIGHT:
                        r.quadrantTo = ThermalNodeFace::TANGENTIAL_LEFT;  // Turn's TOP touches insulation's BOTTOM
                        break;
                    default:
                        r.quadrantTo = ThermalNodeFace::NONE;
                        break;
                }
                
                // For concentric cores, only allow connections to RADIAL faces of turns
                // TANGENTIAL (top/bottom) faces should only connect to ambient via convection
                if (!_isToroidal && 
                    (r.quadrantTo == ThermalNodeFace::TANGENTIAL_LEFT || 
                     r.quadrantTo == ThermalNodeFace::TANGENTIAL_RIGHT)) {
                    continue;
                }
                
                r.type = HeatTransferType::CONDUCTION;
                
                // Calculate contact area using turn's actual dimensions
                double contactArea = turnWidth * turnHeight * 0.5;
                
                // Resistance through insulation: R = thickness / (k * A)
                // Use half the insulation thickness as the conduction path
                double conductionDistance = insulationWidth / 2.0;
                if (!face.isVertical) {
                    conductionDistance = insulationHeight / 2.0;
                }
                
                r.resistance = ThermalResistance::calculateConductionResistance(
                    conductionDistance, insulationK, contactArea);
                
                _resistances.push_back(r);
                
                if (THERMAL_DEBUG) {
                    std::cout << "Insulation " << insulationNode.name << " face=" << static_cast<int>(face.face) 
                              << " -> Turn " << _nodes[turnIdx].name 
                              << " R=" << r.resistance << " K/W" << std::endl;
                }
            }
        }
    }
}

void Temperature::createTurnToSolidConnections() {
    if (!_isToroidal) return;
    
    auto core = _magnetic.get_core();
    auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    
    double coreOuterDiameter = dimensions["A"];
    double coreInnerDiameter = dimensions["B"];
    double windingWindowInnerRadius = coreInnerDiameter / 2.0;
    double windingWindowOuterRadius = coreOuterDiameter / 2.0;
    
    double minConductionDist = getMinimumDistanceForConduction();
    
    if (THERMAL_DEBUG) {
        std::cout << "Turn-to-core conduction: innerRadius=" << windingWindowInnerRadius 
                  << ", outerRadius=" << windingWindowOuterRadius 
                  << ", minConductionDist=" << minConductionDist << std::endl;
    }
    
    std::vector<size_t> coreNodeIndices;
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part == ThermalNodePartType::CORE_TOROIDAL_SEGMENT) {
            coreNodeIndices.push_back(i);
        }
    }
    
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part != ThermalNodePartType::TURN) continue;
        
        const auto& turnNode = _nodes[i];
        
        // Determine if this is inner or outer surface node based on name
        bool isInnerNode = turnNode.name.find("_Inner") != std::string::npos;
        bool isOuterNode = turnNode.name.find("_Outer") != std::string::npos;
        
        // Inner surface nodes: check RADIAL_OUTER quadrant for contact with core's inner surface
        if (isInnerNode) {
            auto* qOuter = turnNode.getQuadrant(ThermalNodeFace::RADIAL_OUTER);
            if (THERMAL_DEBUG && qOuter) {
                double limitX = qOuter->limitCoordinates[0];
                double limitY = qOuter->limitCoordinates[1];
                double limitRadialPos = std::sqrt(limitX*limitX + limitY*limitY);
                double distToInnerCore = limitRadialPos - windingWindowInnerRadius;
                std::cout << "Inner turn node " << turnNode.name << ": limitPos=[" << limitX << "," << limitY 
                          << "], radialPos=" << limitRadialPos << ", distToCore=" << distToInnerCore << std::endl;
            }
            if (qOuter) {
                // Use limitCoordinates for accurate surface position
                double limitX = qOuter->limitCoordinates[0];
                double limitY = qOuter->limitCoordinates[1];
                double limitRadialPos = std::sqrt(limitX*limitX + limitY*limitY);
                
                // Inner node conducts to core's inner surface (windingWindowInnerRadius)
                double distToInnerCore = limitRadialPos - windingWindowInnerRadius;
                
                if (std::abs(distToInnerCore) <= minConductionDist) {
                    size_t closestCoreIdx = 0;
                    double minDist = 1e9;
                    
                    for (size_t coreIdx : coreNodeIndices) {
                        // Use core's RADIAL_INNER limitCoordinates (facing toward winding window)
                        auto* coreQuadrant = _nodes[coreIdx].getQuadrant(ThermalNodeFace::RADIAL_INNER);
                        double coreX = coreQuadrant ? coreQuadrant->limitCoordinates[0] : _nodes[coreIdx].physicalCoordinates[0];
                        double coreY = coreQuadrant ? coreQuadrant->limitCoordinates[1] : _nodes[coreIdx].physicalCoordinates[1];
                        
                        double dx = limitX - coreX;
                        double dy = limitY - coreY;
                        double dist = std::sqrt(dx*dx + dy*dy);
                        if (dist < minDist) {
                            minDist = dist;
                            closestCoreIdx = coreIdx;
                        }
                    }
                    
                    ThermalResistanceElement r;
                    r.nodeFromId = i;
                    r.quadrantFrom = ThermalNodeFace::RADIAL_OUTER;  // Turn's outer face
                    r.nodeToId = closestCoreIdx;
                    r.quadrantTo = ThermalNodeFace::RADIAL_INNER;    // Core's inner face
                    r.type = HeatTransferType::CONDUCTION;
                    

                    
                    int turnIdx = turnNode.turnIndex.value_or(0);
                    
                    // Copper conduction: from turn node center to surface
                    // Node is at wire surface, so conduction length is half the turn's width
                    double turnWidth = turnNode.dimensions.width;
                    double copperLength = turnWidth / 2.0;
                    double copperResistance = ThermalResistance::calculateConductionResistance(
                        copperLength, qOuter->thermalConductivity, qOuter->surfaceArea);
                    
                    // Insulation/enamel resistance
                    double enamelResistance = getInsulationLayerThermalResistance(turnIdx, -1, qOuter->surfaceArea);
                    
                    // Coating resistance
                    double coatingResistance = 0.0;
                    if (qOuter->coating.has_value()) {
                        coatingResistance = WireCoatingUtils::calculateCoatingResistance(qOuter->coating.value(), qOuter->surfaceArea);
                    }
                    
                    double totalResistance = copperResistance + enamelResistance + coatingResistance;
                    r.resistance = totalResistance;
                    
                    // Add turn-to-core insulation from config if enabled
                    if (_config.useTurnToCoreInsulation && _config.turnToCoreInsulationThickness > 0) {
                        r.addInsulationLayer(
                            _config.turnToCoreInsulationThickness,
                            _config.turnToCoreInsulationConductivity,
                            "turn_to_core_insulation",
                            "Insulation between turn and core from config"
                        );
                        r.resistance += r.calculateTotalInsulationResistance(qOuter->surfaceArea);
                    }
                    
                    _resistances.push_back(r);
                }
            }
        }
        
        // Outer surface nodes: check RADIAL_INNER quadrant for contact with core's outer surface
        if (isOuterNode) {
            auto* qInner = turnNode.getQuadrant(ThermalNodeFace::RADIAL_INNER);
            if (THERMAL_DEBUG && qInner) {
                double limitX = qInner->limitCoordinates[0];
                double limitY = qInner->limitCoordinates[1];
                double limitRadialPos = std::sqrt(limitX*limitX + limitY*limitY);
                double distToOuterCore = windingWindowOuterRadius - limitRadialPos;
                std::cout << "Outer turn node " << turnNode.name << ": limitPos=[" << limitX << "," << limitY 
                          << "], radialPos=" << limitRadialPos << ", distToCore=" << distToOuterCore << std::endl;
            }
            if (qInner) {
                // Use limitCoordinates for accurate surface position
                double limitX = qInner->limitCoordinates[0];
                double limitY = qInner->limitCoordinates[1];
                double limitRadialPos = std::sqrt(limitX*limitX + limitY*limitY);
                
                // Outer node conducts to core's outer surface (windingWindowOuterRadius)
                double distToOuterCore = windingWindowOuterRadius - limitRadialPos;
                
                if (std::abs(distToOuterCore) <= minConductionDist) {
                    size_t closestCoreIdx = 0;
                    double minDist = 1e9;
                    
                    for (size_t coreIdx : coreNodeIndices) {
                        // Use core's RADIAL_OUTER limitCoordinates (facing away from center)
                        auto* coreQuadrant = _nodes[coreIdx].getQuadrant(ThermalNodeFace::RADIAL_OUTER);
                        double coreX = coreQuadrant ? coreQuadrant->limitCoordinates[0] : _nodes[coreIdx].physicalCoordinates[0];
                        double coreY = coreQuadrant ? coreQuadrant->limitCoordinates[1] : _nodes[coreIdx].physicalCoordinates[1];
                        
                        double dx = limitX - coreX;
                        double dy = limitY - coreY;
                        double dist = std::sqrt(dx*dx + dy*dy);
                        if (dist < minDist) {
                            minDist = dist;
                            closestCoreIdx = coreIdx;
                        }
                    }
                    
                    ThermalResistanceElement r;
                    r.nodeFromId = i;
                    r.quadrantFrom = ThermalNodeFace::RADIAL_INNER;  // Turn's inner face (toward center)
                    r.nodeToId = closestCoreIdx;
                    r.quadrantTo = ThermalNodeFace::RADIAL_OUTER;     // Core's outer face (away from center)
                    r.type = HeatTransferType::CONDUCTION;
                    

                    
                    int turnIdx = turnNode.turnIndex.value_or(0);
                    
                    // Copper conduction: from turn node center to surface
                    double turnWidth = turnNode.dimensions.width;
                    double copperLength = turnWidth / 2.0;
                    double copperResistance = ThermalResistance::calculateConductionResistance(
                        copperLength, qInner->thermalConductivity, qInner->surfaceArea);
                    
                    // Insulation/enamel resistance
                    double enamelResistance = getInsulationLayerThermalResistance(turnIdx, -1, qInner->surfaceArea);
                    
                    // Coating resistance
                    double coatingResistance = 0.0;
                    if (qInner->coating.has_value()) {
                        coatingResistance = WireCoatingUtils::calculateCoatingResistance(qInner->coating.value(), qInner->surfaceArea);
                    }
                    
                    double totalResistance = copperResistance + enamelResistance + coatingResistance;
                    r.resistance = totalResistance;
                    
                    // Add turn-to-core insulation from config if enabled
                    if (_config.useTurnToCoreInsulation && _config.turnToCoreInsulationThickness > 0) {
                        r.addInsulationLayer(
                            _config.turnToCoreInsulationThickness,
                            _config.turnToCoreInsulationConductivity,
                            "turn_to_core_insulation",
                            "Insulation between turn and core from config"
                        );
                        r.resistance += r.calculateTotalInsulationResistance(qInner->surfaceArea);
                    }
                    
                    _resistances.push_back(r);
                }
            }
        }
    }
}

// ============================================================================
// Convection Connections - Exposed Quadrant Detection
// ============================================================================

void Temperature::createConvectionConnections() {
    size_t ambientIdx = _nodes.size() - 1;
    
    // Extract wire properties locally for convection calculations
    auto wireOpt = extractWire();
    double wireWidth = 0.001;  // Default 1mm
    double wireHeight = 0.001;
    bool isRound = false;
    if (wireOpt) {
        auto [w, h] = getWireDimensions(wireOpt.value());
        wireWidth = w;
        wireHeight = h;
        isRound = isRoundWire(wireOpt.value());
    }
    
    double maxConvectionDist = getMaximumConvectionDistance(wireWidth, wireHeight, isRound);
    double minConductionDist = getMinimumConductionDistance(wireWidth, wireHeight, isRound);
    
    // Get convection coefficient
    double surfaceTemp = _config.ambientTemperature + 30.0;
    double h_conv;
    
    if (_config.includeForcedConvection) {
        h_conv = ThermalResistance::calculateForcedConvectionCoefficient(
            _config.airVelocity, wireWidth, _config.ambientTemperature);
    } else {
        h_conv = ThermalResistance::calculateNaturalConvectionCoefficient(
            surfaceTemp, _config.ambientTemperature, wireWidth, SurfaceOrientation::VERTICAL);
    }
    
    if (_config.includeRadiation) {
        double h_rad = ThermalResistance::calculateRadiationCoefficient(
            surfaceTemp, _config.ambientTemperature, _config.surfaceEmissivity);
        h_conv += h_rad;
    }
    
    if (_isToroidal) {
        // Build a map of which quadrants are already connected by conduction
        // Key: "nodeId_quadrantFace" -> true if connected
        std::set<std::string> connectedQuadrants;
        for (const auto& res : _resistances) {
            if (res.type == HeatTransferType::CONDUCTION) {
                std::string key1 = std::to_string(res.nodeFromId) + "_" + 
                                   std::string(magic_enum::enum_name(res.quadrantFrom));
                std::string key2 = std::to_string(res.nodeToId) + "_" + 
                                   std::string(magic_enum::enum_name(res.quadrantTo));
                connectedQuadrants.insert(key1);
                connectedQuadrants.insert(key2);
            }
        }
        
        // For each turn node, check each quadrant for exposure to air
        for (size_t i = 0; i < _nodes.size(); i++) {
            if (_nodes[i].part != ThermalNodePartType::TURN) continue;
            
            const auto& node = _nodes[i];
            double nodeX = node.physicalCoordinates[0];
            double nodeY = node.physicalCoordinates[1];
            double nodeR = std::sqrt(nodeX*nodeX + nodeY*nodeY);
            double nodeAngle = std::atan2(nodeY, nodeX);
            double nodeWidth = node.dimensions.width;
            double nodeHeight = node.dimensions.height;
            
            // Check each quadrant
            for (int qIdx = 0; qIdx < 4; ++qIdx) {
                ThermalNodeFace face = node.quadrants[qIdx].face;
                if (face == ThermalNodeFace::NONE) continue;
                
                // Skip if this quadrant is already connected by conduction
                std::string qKey = std::to_string(i) + "_" + std::string(magic_enum::enum_name(face));
                if (connectedQuadrants.count(qKey) > 0) continue;
                
                bool isExposed = true;
                
                // Check for blocking objects in the quadrant's direction
                // Purely geometric: block if there's a turn in the quadrant direction
                // that is significantly offset radially (not a tangential neighbor)
                // AND within this node's max dimension distance
                double maxBlockingDist = std::max(nodeWidth, nodeHeight);
                // Minimum radial difference to distinguish radial from tangential neighbors
                double minRadialDiff = std::min(nodeWidth, nodeHeight) / 4.0;
                
                if (face == ThermalNodeFace::RADIAL_INNER) {
                    // Check for any object significantly closer to center in this direction
                    for (size_t j = 0; j < _nodes.size(); j++) {
                        if (i == j) continue;
                        
                        double otherX = _nodes[j].physicalCoordinates[0];
                        double otherY = _nodes[j].physicalCoordinates[1];
                        double otherR = std::sqrt(otherX*otherX + otherY*otherY);
                        double otherAngle = std::atan2(otherY, otherX);
                        
                        double radialDiff = nodeR - otherR;
                        // Must be: significantly closer (not tangential) AND within max blocking distance
                        if (radialDiff > minRadialDiff && radialDiff < maxBlockingDist) {
                            double angleDiff = std::abs(nodeAngle - otherAngle);
                            while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                            
                            // Block if there's an object in roughly same angular direction
                            if (std::abs(angleDiff) < 0.3) {
                                isExposed = false;
                                break;
                            }
                        }
                    }
                }
                else if (face == ThermalNodeFace::RADIAL_OUTER) {
                    // Check for any object significantly farther from center in this direction
                    for (size_t j = 0; j < _nodes.size(); j++) {
                        if (i == j) continue;
                        
                        double otherX = _nodes[j].physicalCoordinates[0];
                        double otherY = _nodes[j].physicalCoordinates[1];
                        double otherR = std::sqrt(otherX*otherX + otherY*otherY);
                        double otherAngle = std::atan2(otherY, otherX);
                        
                        double radialDiff = otherR - nodeR;
                        // Must be: significantly farther (not tangential) AND within max blocking distance
                        if (radialDiff > minRadialDiff && radialDiff < maxBlockingDist) {
                            double angleDiff = std::abs(nodeAngle - otherAngle);
                            while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                            
                            // Block if there's an object in roughly same angular direction
                            if (std::abs(angleDiff) < 0.3) {
                                isExposed = false;
                                break;
                            }
                        }
                    }
                }
                else if (face == ThermalNodeFace::TANGENTIAL_LEFT || face == ThermalNodeFace::TANGENTIAL_RIGHT) {
                    // For tangential faces, check adjacent turns
                    // If there's a conduction connection, this face is not exposed
                    // Otherwise, check if there's a blocking turn within convection distance
                    for (size_t j = 0; j < _nodes.size(); j++) {
                        if (i == j) continue;
                        if (_nodes[j].part != ThermalNodePartType::TURN) continue;
                        
                        double otherX = _nodes[j].physicalCoordinates[0];
                        double otherY = _nodes[j].physicalCoordinates[1];
                        double otherR = std::sqrt(otherX*otherX + otherY*otherY);
                        double otherAngle = std::atan2(otherY, otherX);
                        
                        // Similar radius (same "layer" inner or outer)
                        if (std::abs(otherR - nodeR) < nodeWidth) {
                            double angleDiff = otherAngle - nodeAngle;
                            while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                            while (angleDiff < -M_PI) angleDiff += 2 * M_PI;
                            
                            bool isLeft = (face == ThermalNodeFace::TANGENTIAL_LEFT);
                            double distAlongTangent = std::abs(angleDiff) * nodeR;
                            
                            // Check if there's a turn blocking this direction
                            if (isLeft && angleDiff < 0 && distAlongTangent < maxConvectionDist) {
                                isExposed = false;
                                break;
                            }
                            if (!isLeft && angleDiff > 0 && distAlongTangent < maxConvectionDist) {
                                isExposed = false;
                                break;
                            }
                        }
                    }
                }
                
                // If exposed, connect to ambient
                if (isExposed) {
                    auto* q = node.getQuadrant(face);
                    if (q && q->surfaceArea > 0) {
                        ThermalResistanceElement r;
                        r.nodeFromId = i;
                        r.quadrantFrom = face;
                        r.nodeToId = ambientIdx;
                        r.quadrantTo = ThermalNodeFace::NONE;
                        r.type = _config.includeForcedConvection ? 
                                 HeatTransferType::FORCED_CONVECTION : 
                                 HeatTransferType::NATURAL_CONVECTION;
                        // Use coating-aware calculation if coating exists
                        r.resistance = q->calculateConvectionResistance(h_conv);
                        if (q->coating.has_value()) {
                            r.resistance += WireCoatingUtils::calculateCoatingResistance(q->coating.value(), q->surfaceArea);
                        }
                        _resistances.push_back(r);
                    } else if (THERMAL_DEBUG) {
                        std::cout << "    NOT CREATED: q==null or surfaceArea=0" << std::endl;
                    }
                }
            }
        }
        
        // Core convection - connect exposed core quadrants
        auto core = _magnetic.get_core();
        auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
        double coreInnerR = dimensions["B"] / 2.0;
        double coreOuterR = dimensions["A"] / 2.0;
        
        for (size_t i = 0; i < _nodes.size(); i++) {
            if (_nodes[i].part != ThermalNodePartType::CORE_TOROIDAL_SEGMENT) continue;
            
            // Check RADIAL_OUTER - exposed if no turn is close to outer core
            bool outerBlocked = false;
            for (size_t j = 0; j < _nodes.size(); j++) {
                if (_nodes[j].part != ThermalNodePartType::TURN) continue;
                
                double turnR = std::sqrt(
                    _nodes[j].physicalCoordinates[0]*_nodes[j].physicalCoordinates[0] +
                    _nodes[j].physicalCoordinates[1]*_nodes[j].physicalCoordinates[1]
                );
                double turnWidth = _nodes[j].dimensions.width;
                
                if (std::abs(turnR - coreOuterR) < turnWidth) {
                    outerBlocked = true;
                    break;
                }
            }
            
            if (!outerBlocked) {
                std::string qKey = std::to_string(i) + "_" + 
                                  std::string(magic_enum::enum_name(ThermalNodeFace::RADIAL_OUTER));
                if (connectedQuadrants.count(qKey) == 0) {
                    auto* q = _nodes[i].getQuadrant(ThermalNodeFace::RADIAL_OUTER);
                    if (q) {
                        ThermalResistanceElement r;
                        r.nodeFromId = i;
                        r.quadrantFrom = ThermalNodeFace::RADIAL_OUTER;
                        r.nodeToId = ambientIdx;
                        r.quadrantTo = ThermalNodeFace::NONE;
                        r.type = _config.includeForcedConvection ? 
                                 HeatTransferType::FORCED_CONVECTION : 
                                 HeatTransferType::NATURAL_CONVECTION;
                        r.resistance = q->calculateConvectionResistance(h_conv);
                        _resistances.push_back(r);
                    }
                }
            }
            
            // Check RADIAL_INNER - exposed if no turn is close to inner core
            bool innerBlocked = false;
            for (size_t j = 0; j < _nodes.size(); j++) {
                if (_nodes[j].part != ThermalNodePartType::TURN) continue;
                
                double turnR = std::sqrt(
                    _nodes[j].physicalCoordinates[0]*_nodes[j].physicalCoordinates[0] +
                    _nodes[j].physicalCoordinates[1]*_nodes[j].physicalCoordinates[1]
                );
                double turnWidth = _nodes[j].dimensions.width;
                
                if (std::abs(turnR - coreInnerR) < turnWidth) {
                    innerBlocked = true;
                    break;
                }
            }
            
            if (!innerBlocked) {
                std::string qKey = std::to_string(i) + "_" + 
                                  std::string(magic_enum::enum_name(ThermalNodeFace::RADIAL_INNER));
                if (connectedQuadrants.count(qKey) == 0) {
                    auto* q = _nodes[i].getQuadrant(ThermalNodeFace::RADIAL_INNER);
                    if (q) {
                        ThermalResistanceElement r;
                        r.nodeFromId = i;
                        r.quadrantFrom = ThermalNodeFace::RADIAL_INNER;
                        r.nodeToId = ambientIdx;
                        r.quadrantTo = ThermalNodeFace::NONE;
                        r.type = _config.includeForcedConvection ? 
                                 HeatTransferType::FORCED_CONVECTION : 
                                 HeatTransferType::NATURAL_CONVECTION;
                        r.resistance = q->calculateConvectionResistance(h_conv);
                        _resistances.push_back(r);
                    }
                }
            }
        }
    } else {
        // Check if we have concentric core nodes
        bool hasConcentricCoreNodes = false;
        for (const auto& node : _nodes) {
            if (node.part == ThermalNodePartType::CORE_CENTRAL_COLUMN ||
                node.part == ThermalNodePartType::CORE_TOP_YOKE ||
                node.part == ThermalNodePartType::CORE_BOTTOM_YOKE) {
                hasConcentricCoreNodes = true;
                break;
            }
        }
        
        if (hasConcentricCoreNodes) {
            // Concentric core - quadrant-specific convection with symmetry boundary
            // Build a map of which quadrants are already connected by conduction
            std::set<std::string> connectedQuadrants;
            for (const auto& res : _resistances) {
                if (res.type == HeatTransferType::CONDUCTION) {
                    std::string key1 = std::to_string(res.nodeFromId) + "_" + 
                                       std::string(magic_enum::enum_name(res.quadrantFrom));
                    std::string key2 = std::to_string(res.nodeToId) + "_" + 
                                       std::string(magic_enum::enum_name(res.quadrantTo));
                    connectedQuadrants.insert(key1);
                    connectedQuadrants.insert(key2);
                }
            }
            
            for (size_t i = 0; i < _nodes.size(); i++) {
                if (_nodes[i].part == ThermalNodePartType::AMBIENT) continue;
                
                // Check if this is a central column node or yoke (symmetry boundary on LEFT side)
                bool isCentralColumn = (_nodes[i].part == ThermalNodePartType::CORE_CENTRAL_COLUMN);
                bool isYoke = (_nodes[i].part == ThermalNodePartType::CORE_TOP_YOKE ||
                              _nodes[i].part == ThermalNodePartType::CORE_BOTTOM_YOKE);
                bool isTurn = (_nodes[i].part == ThermalNodePartType::TURN);
                
                // Check each quadrant for convection exposure
                for (int qIdx = 0; qIdx < 4; ++qIdx) {
                    ThermalNodeFace face = _nodes[i].quadrants[qIdx].face;
                    if (face == ThermalNodeFace::NONE) continue;
                    
                    // Skip if this quadrant is already connected by conduction
                    std::string qKey = std::to_string(i) + "_" + std::string(magic_enum::enum_name(face));
                    if (connectedQuadrants.count(qKey) > 0) continue;
                    
                    // Symmetry boundary: Central column's LEFT face (RADIAL_INNER) is adiabatic
                    if (isCentralColumn && face == ThermalNodeFace::RADIAL_INNER) continue;
                    
                    // Symmetry boundary: Yoke's LEFT face (RADIAL_INNER) is adiabatic
                    if (isYoke && face == ThermalNodeFace::RADIAL_INNER) continue;
                    
                    // For turns: check if quadrant is blocked by another turn
                    if (isTurn) {
                        bool isBlocked = false;
                        double turnX = _nodes[i].physicalCoordinates[0];
                        double turnY = _nodes[i].physicalCoordinates[1];
                        
                        // Determine direction based on face
                        double dirX = 0.0, dirY = 0.0;
                        if (face == ThermalNodeFace::RADIAL_INNER) dirX = -1.0;
                        else if (face == ThermalNodeFace::RADIAL_OUTER) dirX = 1.0;
                        else if (face == ThermalNodeFace::TANGENTIAL_LEFT) dirY = 1.0;
                        else if (face == ThermalNodeFace::TANGENTIAL_RIGHT) dirY = -1.0;
                        
                        // Check for blocking turns in that direction
                        double turnWidth = _nodes[i].dimensions.width;
                        double turnHeight = _nodes[i].dimensions.height;
                        
                        for (size_t j = 0; j < _nodes.size(); j++) {
                            if (i == j) continue;
                            if (_nodes[j].part != ThermalNodePartType::TURN) continue;
                            
                            double otherX = _nodes[j].physicalCoordinates[0];
                            double otherY = _nodes[j].physicalCoordinates[1];
                            double dx = otherX - turnX;
                            double dy = otherY - turnY;
                            
                            // Check if other turn is in the blocking direction
                            bool inDirection = false;
                            if (dirX != 0.0 && std::abs(dy) < turnHeight && dx * dirX > 0 && std::abs(dx) < turnWidth * 1.5) {
                                inDirection = true;
                            }
                            if (dirY != 0.0 && std::abs(dx) < turnWidth && dy * dirY > 0 && std::abs(dy) < turnHeight * 1.5) {
                                inDirection = true;
                            }
                            
                            if (inDirection) {
                                isBlocked = true;
                                break;
                            }
                        }
                        
                        // Check if there's an insulation layer to the right of this turn
                        // If so, block RIGHT, TOP, and BOTTOM faces from convection
                        bool hasInsulationToRight = false;
                        double turnRightEdge = turnX + turnWidth / 2;
                        double turnTop = turnY + turnHeight / 2;
                        double turnBottom = turnY - turnHeight / 2;
                        
                        for (size_t j = 0; j < _nodes.size(); j++) {
                            if (_nodes[j].part != ThermalNodePartType::INSULATION_LAYER) continue;
                            
                            double insX = _nodes[j].physicalCoordinates[0];
                            double insY = _nodes[j].physicalCoordinates[1];
                            double insWidth = _nodes[j].dimensions.width;
                            double insHeight = _nodes[j].dimensions.height;
                            
                            double insLeftEdge = insX - insWidth / 2;
                            double insTop = insY + insHeight / 2;
                            double insBottom = insY - insHeight / 2;
                            
                            // Check if insulation layer is to the right of this turn
                            // AND overlaps in Y (vertically)
                            bool isToTheRight = insLeftEdge >= turnRightEdge;
                            bool overlapsVertically = !(insBottom > turnTop || insTop < turnBottom);
                            
                            if (isToTheRight && overlapsVertically) {
                                hasInsulationToRight = true;
                                if (THERMAL_DEBUG) {
                                    std::cout << "Turn " << _nodes[i].name << " has insulation layer to right: " 
                                              << _nodes[j].name << std::endl;
                                }
                                break;
                            }
                        }
                        
                        // Block RIGHT, TOP, and BOTTOM faces if there's insulation to the right
                        if (hasInsulationToRight && 
                            (face == ThermalNodeFace::RADIAL_OUTER || 
                             face == ThermalNodeFace::TANGENTIAL_LEFT || 
                             face == ThermalNodeFace::TANGENTIAL_RIGHT)) {
                            isBlocked = true;
                            // Update surface coverage to show as covered (for schematic visualization)
                            auto* q = _nodes[i].getQuadrant(face);
                            if (q) {
                                q->surfaceCoverage = 0.0;
                            }
                            if (THERMAL_DEBUG) {
                                std::cout << "Turn " << _nodes[i].name << " " 
                                          << magic_enum::enum_name(face) 
                                          << " face blocked by insulation to right" << std::endl;
                            }
                        }
                        
                        if (isBlocked) continue;
                    }
                    
                    // All other exposed quadrants get convection
                    auto* q = _nodes[i].getQuadrant(face);
                    if (q && q->surfaceArea > 0) {
                        ThermalResistanceElement r;
                        r.nodeFromId = i;
                        r.quadrantFrom = face;
                        r.nodeToId = ambientIdx;
                        r.quadrantTo = ThermalNodeFace::NONE;
                        r.type = _config.includeForcedConvection ? 
                                 HeatTransferType::FORCED_CONVECTION : 
                                 HeatTransferType::NATURAL_CONVECTION;
                        r.resistance = q->calculateConvectionResistance(h_conv);
                        _resistances.push_back(r);
                    }
                }
            }
        } else {
            // Non-toroidal, non-concentric core (e.g., E/ETD shapes) - simplified convection
            for (size_t i = 0; i < _nodes.size(); i++) {
                if (_nodes[i].part == ThermalNodePartType::AMBIENT) continue;
                
                ThermalResistanceElement r;
                r.nodeFromId = i;
                r.quadrantFrom = ThermalNodeFace::NONE;
                r.nodeToId = ambientIdx;
                r.quadrantTo = ThermalNodeFace::NONE;
                r.type = _config.includeForcedConvection ? 
                         HeatTransferType::FORCED_CONVECTION : 
                         HeatTransferType::NATURAL_CONVECTION;
                r.resistance = ThermalResistance::calculateConvectionResistance(
                    h_conv, _nodes[i].getTotalSurfaceArea());
                _resistances.push_back(r);
            }
        }
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

double Temperature::calculateSurfaceDistance(const ThermalNetworkNode& node1, 
                                              const ThermalNetworkNode& node2) const {
    double dx = node1.physicalCoordinates[0] - node2.physicalCoordinates[0];
    double dy = node1.physicalCoordinates[1] - node2.physicalCoordinates[1];
    double dz = node1.physicalCoordinates[2] - node2.physicalCoordinates[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

bool Temperature::shouldConnectQuadrants(const ThermalNetworkNode& node1, ThermalNodeFace face1,
                                          const ThermalNetworkNode& node2, ThermalNodeFace face2) const {
    double dist = calculateSurfaceDistance(node1, node2);
    // Use average dimensions for conduction threshold
    double avgWidth = (node1.dimensions.width + node2.dimensions.width) / 2.0;
    double avgHeight = (node1.dimensions.height + node2.dimensions.height) / 2.0;
    // Assume rectangular for threshold calculation (conservative)
    double minConductionDist = std::min(avgWidth, avgHeight) * 0.75;
    return dist <= minConductionDist;
}

double Temperature::calculateContactArea(const ThermalNodeQuadrant& q1, const ThermalNodeQuadrant& q2) const {
    double minLength = std::min(q1.length, q2.length);
    // Use quadrant dimensions if available, otherwise use a default
    double height = q1.surfaceArea / (q1.length > 1e-9 ? q1.length : 0.001);
    return height * minLength;
}

double Temperature::getInsulationLayerThermalResistance(int turnIdx1, int turnIdx2, double contactArea) {
    if (contactArea <= 0.0) {
        return 1e9;
    }
    
    try {
        auto coil = _magnetic.get_coil();
        auto turnsDescription = coil.get_turns_description();
        
        // For turn-to-solid connections (turnIdx2 = -1)
        if (turnIdx2 < 0) {
            double insulationThickness = 0.00005;  // 50 microns enamel
            double insulationK = 0.2;  // Polyurethane
            
            double resistance = ThermalResistance::calculateConductionResistance(
                insulationThickness, insulationK, contactArea);
            
            return std::max(resistance, 0.001);
        }
        
        if (!turnsDescription || turnIdx1 < 0 || turnIdx1 >= static_cast<int>(turnsDescription->size()) || 
            turnIdx2 >= static_cast<int>(turnsDescription->size())) {
            return 0.001;
        }
        
        auto& turn1 = (*turnsDescription)[turnIdx1];
        auto& turn2 = (*turnsDescription)[turnIdx2];
        
        auto layersBetween = StrayCapacitance::get_insulation_layers_between_two_turns(turn1, turn2, coil);
        
        double totalLayerResistance = 0.0;
        for (const auto& layer : layersBetween) {
            double layerThickness = coil.get_insulation_layer_thickness(layer);
            double layerK = 0.2;
            
            auto insulationMaterialOpt = layer.get_insulation_material();
            if (insulationMaterialOpt.has_value()) {
                auto& insulationMaterial = insulationMaterialOpt.value();
                if (std::holds_alternative<MAS::InsulationMaterial>(insulationMaterial)) {
                    auto material = std::get<MAS::InsulationMaterial>(insulationMaterial);
                    if (material.get_thermal_conductivity().has_value()) {
                        layerK = material.get_thermal_conductivity().value();
                    }
                } else {
                    std::string materialName = std::get<std::string>(insulationMaterial);
                    layerK = ThermalResistance::getMaterialThermalConductivity(materialName);
                }
            }
            
            totalLayerResistance += ThermalResistance::calculateConductionResistance(
                layerThickness, layerK, contactArea);
        }
        
        return totalLayerResistance;
    } catch (...) {
        return 0.001;
    }
}

bool Temperature::hasBobbinNodes() const {
    for (const auto& node : _nodes) {
        if (node.part == ThermalNodePartType::BOBBIN_CENTRAL_COLUMN) {
            return true;
        }
    }
    return false;
}

void Temperature::calculateSchematicScaling() {
    if (_nodes.empty()) {
        _scaleFactor = 1.0;
        return;
    }
    
    if (THERMAL_DEBUG) {
        std::cout << "calculateSchematicScaling: " << _nodes.size() << " nodes" << std::endl;
    }
    
    // Find bounding box of all physical coordinates
    double minX = 1e9, maxX = -1e9;
    double minY = 1e9, maxY = -1e9;
    
    for (const auto& node : _nodes) {
        if (node.physicalCoordinates.size() >= 2) {
            minX = std::min(minX, node.physicalCoordinates[0]);
            maxX = std::max(maxX, node.physicalCoordinates[0]);
            minY = std::min(minY, node.physicalCoordinates[1]);
            maxY = std::max(maxY, node.physicalCoordinates[1]);
        }
    }
    
    // Calculate dimensions
    double width = maxX - minX;
    double height = maxY - minY;
    
    // Target schematic size (pixels)
    const double targetWidth = 400.0;
    const double targetHeight = 400.0;
    const double margin = 50.0;
    
    // Calculate scale factor to fit within target size
    double scaleX = (targetWidth - 2 * margin) / (width > 1e-6 ? width : 1.0);
    double scaleY = (targetHeight - 2 * margin) / (height > 1e-6 ? height : 1.0);
    _scaleFactor = std::min(scaleX, scaleY);
    
    // Center offset
    double centerX = (minX + maxX) / 2.0;
    double centerY = (minY + maxY) / 2.0;
    double schematicCenterX = targetWidth / 2.0;
    double schematicCenterY = targetHeight / 2.0;
    
    // Set schematic coordinates for each node based on physical coordinates
    for (auto& node : _nodes) {
        if (node.physicalCoordinates.size() >= 2) {
            double x = node.physicalCoordinates[0];
            double y = node.physicalCoordinates[1];
            
            // Transform: scale and center
            double schematicX = schematicCenterX + (x - centerX) * _scaleFactor;
            double schematicY = schematicCenterY + (y - centerY) * _scaleFactor;
            
            node.schematicCoordinates = {schematicX, schematicY};
        } else {
            // Default position for nodes without coordinates
            node.schematicCoordinates = {schematicCenterX, schematicCenterY};
        }
    }
    
    if (THERMAL_DEBUG) {
        std::cout << "Schematic scaling: scale=" << _scaleFactor 
                  << ", bounds=[" << minX*1000 << "," << minY*1000 << " to " 
                  << maxX*1000 << "," << maxY*1000 << "]mm" << std::endl;
    }
}

void Temperature::plotSchematic() {
    if (_config.schematicOutputPath.empty()) return;
    
    std::filesystem::path outFile = _config.schematicOutputPath;
    std::filesystem::create_directories(outFile.parent_path());
    
    json schematic;
    schematic["nodes"] = json::array();
    for (const auto& node : _nodes) {
        schematic["nodes"].push_back(node.toJson());
    }
    
    schematic["resistances"] = json::array();
    for (const auto& res : _resistances) {
        json r;
        r["nodeFromId"] = res.nodeFromId;
        r["quadrantFrom"] = magic_enum::enum_name(res.quadrantFrom);
        r["nodeToId"] = res.nodeToId;
        r["quadrantTo"] = magic_enum::enum_name(res.quadrantTo);
        r["type"] = magic_enum::enum_name(res.type);
        r["resistance"] = res.resistance;
        schematic["resistances"].push_back(r);
    }
    
    std::ofstream file(outFile.string() + ".json");
    file << schematic.dump(2);
    file.close();
}

// ============================================================================
// Solver
// ============================================================================

ThermalResult Temperature::solveThermalCircuit() {
    size_t n = _nodes.size();
    if (n == 0) {
        ThermalResult result;
        result.converged = false;
        result.maximumTemperature = _config.ambientTemperature;
        result.totalThermalResistance = 0;
        return result;
    }
    
    size_t ambientIdx = n - 1;
    
    std::vector<double> temperatures(n, _config.ambientTemperature);
    std::vector<double> powerInputs(n, 0.0);
    
    for (size_t i = 0; i < n; ++i) {
        powerInputs[i] = _nodes[i].powerDissipation;
    }
    
    size_t iteration = 0;
    bool converged = false;
    std::vector<double> oldTemperatures = temperatures;
    
    while (iteration < _config.maxIterations && !converged) {
        SimpleMatrix G(n, n, 0.0);
        
        for (const auto& res : _resistances) {
            double g = 1.0 / std::max(res.resistance, 1e-9);
            
            size_t i = res.nodeFromId;
            size_t j = res.nodeToId;
            
            G(i, i) += g;
            if (j < n) {
                G(j, j) += g;
                G(i, j) -= g;
                G(j, i) -= g;
            }
        }
        
        G.setRowZero(ambientIdx);
        G(ambientIdx, ambientIdx) = 1.0;
        powerInputs[ambientIdx] = _config.ambientTemperature;
        
        try {
            temperatures = SimpleMatrix::solve(G, powerInputs);
        } catch (const std::exception& e) {
            if (THERMAL_DEBUG) {
                std::cerr << "Solver error: " << e.what() << std::endl;
            }
            break;
        }
        
        converged = true;
        for (size_t i = 0; i < n; ++i) {
            if (std::abs(temperatures[i] - oldTemperatures[i]) > _config.convergenceTolerance) {
                converged = false;
                break;
            }
        }
        
        oldTemperatures = temperatures;
        iteration++;
    }
    
    for (size_t i = 0; i < n; ++i) {
        _nodes[i].temperature = temperatures[i];
    }
    
    ThermalResult result;
    result.converged = converged;
    result.iterationsToConverge = iteration;
    result.thermalResistances = _resistances;
    
    result.maximumTemperature = _config.ambientTemperature;
    for (size_t i = 0; i < n - 1; ++i) {
        result.nodeTemperatures[_nodes[i].name] = temperatures[i];
        if (temperatures[i] > result.maximumTemperature) {
            result.maximumTemperature = temperatures[i];
        }
    }
    
    double totalPower = 0.0;
    for (size_t i = 0; i < n; ++i) {
        totalPower += _nodes[i].powerDissipation;
    }
    
    if (totalPower > 0) {
        result.totalThermalResistance = (result.maximumTemperature - _config.ambientTemperature) / totalPower;
    } else {
        result.totalThermalResistance = 0;
    }
    
    double coreTempSum = 0.0;
    size_t coreCount = 0;
    double coilTempSum = 0.0;
    size_t coilCount = 0;
    
    for (size_t i = 0; i < n; ++i) {
        if (_nodes[i].part == ThermalNodePartType::CORE_TOROIDAL_SEGMENT ||
            _nodes[i].part == ThermalNodePartType::CORE_CENTRAL_COLUMN ||
            _nodes[i].part == ThermalNodePartType::CORE_LATERAL_COLUMN) {
            coreTempSum += temperatures[i];
            coreCount++;
        } else if (_nodes[i].part == ThermalNodePartType::TURN) {
            coilTempSum += temperatures[i];
            coilCount++;
        }
    }
    
    result.averageCoreTemperature = coreCount > 0 ? coreTempSum / coreCount : _config.ambientTemperature;
    result.averageCoilTemperature = coilCount > 0 ? coilTempSum / coilCount : _config.ambientTemperature;
    result.methodUsed = "Quadrant-based Thermal Equivalent Circuit";
    
    return result;
}

double Temperature::getBulkThermalResistance() const {
    double totalPower = 0.0;
    for (const auto& node : _nodes) {
        totalPower += node.powerDissipation;
    }
    
    if (totalPower <= 0) return 0.0;
    
    double maxTemp = _config.ambientTemperature;
    for (const auto& node : _nodes) {
        if (node.temperature > maxTemp) {
            maxTemp = node.temperature;
        }
    }
    
    return (maxTemp - _config.ambientTemperature) / totalPower;
}

double Temperature::getTemperatureAtPoint(const std::vector<double>& point) const {
    if (point.size() < 2) return _config.ambientTemperature;
    
    double minDist = 1e9;
    double nearestTemp = _config.ambientTemperature;
    
    for (const auto& node : _nodes) {
        if (node.part == ThermalNodePartType::AMBIENT) continue;
        
        double dx = node.physicalCoordinates[0] - point[0];
        double dy = node.physicalCoordinates[1] - point[1];
        double dist = std::sqrt(dx*dx + dy*dy);
        
        if (dist < minDist) {
            minDist = dist;
            nearestTemp = node.temperature;
        }
    }
    
    return nearestTemp;
}

} // namespace OpenMagnetics
