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
                // For rectangular wires, width is constant
                // For round wires, projected width varies with radius
                double innerWidth = _wireWidth;
                double outerWidth = _wireWidth;
                
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
    
    if (columns.empty()) {
        return;
    }
    
    double coreK = _config.coreThermalConductivity;
    
    for (size_t i = 0; i < columns.size(); ++i) {
        ThermalNetworkNode node;
        node.part = (columns[i].get_type() == ColumnType::CENTRAL) 
                    ? ThermalNodePartType::CORE_CENTRAL_COLUMN 
                    : ThermalNodePartType::CORE_LATERAL_COLUMN;
        node.name = "Core_Column_" + std::to_string(i);
        node.temperature = _config.ambientTemperature;
        node.powerDissipation = _config.coreLosses / columns.size();
        
        double width = columns[i].get_width();
        double depth = columns[i].get_depth();
        double height = columns[i].get_height();
        
        if (columns[i].get_type() == ColumnType::CENTRAL) {
            node.physicalCoordinates = {0, 0, 0};
        } else {
            double offset = width * 2.0;
            node.physicalCoordinates = {(i == 1) ? offset : -offset, 0, 0};
        }
        
        node.initializeConcentricQuadrant(width, depth, height, coreK);
        _nodes.push_back(node);
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
    
    ThermalNetworkNode columnNode;
    columnNode.part = ThermalNodePartType::BOBBIN_CENTRAL_COLUMN;
    columnNode.name = "Bobbin_CentralColumn";
    columnNode.temperature = _config.ambientTemperature;
    columnNode.powerDissipation = 0.0;
    columnNode.physicalCoordinates = {0, 0, 0};
    
    double bobbinK = 0.2;
    columnNode.initializeConcentricQuadrant(0.01, 0.01, 0.01, bobbinK);
    _nodes.push_back(columnNode);
}

void Temperature::createTurnNodes() {
    auto coil = _magnetic.get_coil();
    auto turns = coil.get_turns_description();
    
    if (!turns || turns->empty()) {
        return;
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
            
            // Get wire dimensions from turn
            double wireWidth = _wireWidth;
            double wireHeight = _wireHeight;
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
                                                       _wireThermalCond, true, innerSurfaceRadius,
                                                       _wireCoating,
                                                       _isRoundWire ? TurnCrossSectionalShape::ROUND : TurnCrossSectionalShape::RECTANGULAR);
            } else {
                innerNode.powerDissipation = turnLoss;  // All loss to inner node
                innerNode.initializeToroidalQuadrants(wireWidth, wireHeight, totalTurnLength, 
                                                       _wireThermalCond, true, innerSurfaceRadius,
                                                       _wireCoating,
                                                       _isRoundWire ? TurnCrossSectionalShape::ROUND : TurnCrossSectionalShape::RECTANGULAR);
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
                                                       _wireThermalCond, false, outerSurfaceRadius,
                                                       _wireCoating,
                                                       _isRoundWire ? TurnCrossSectionalShape::ROUND : TurnCrossSectionalShape::RECTANGULAR);
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
            
            double turnLength = turn.get_length();
            node.initializeConcentricQuadrant(_wireWidth, _wireHeight, turnLength, _wireThermalCond,
                                               _wireCoating);
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
    std::vector<size_t> coreNodeIndices;
    for (size_t i = 0; i < _nodes.size(); ++i) {
        if (_nodes[i].part == ThermalNodePartType::CORE_CENTRAL_COLUMN ||
            _nodes[i].part == ThermalNodePartType::CORE_LATERAL_COLUMN) {
            coreNodeIndices.push_back(i);
        }
    }
    
    for (size_t i = 0; i < coreNodeIndices.size(); ++i) {
        for (size_t j = i + 1; j < coreNodeIndices.size(); ++j) {
            size_t node1Idx = coreNodeIndices[i];
            size_t node2Idx = coreNodeIndices[j];
            
            double dist = calculateSurfaceDistance(_nodes[node1Idx], _nodes[node2Idx]);
            
            if (dist < 0.05) {
                ThermalResistanceElement r;
                r.nodeFromId = node1Idx;
                r.quadrantFrom = ThermalNodeFace::NONE;
                r.nodeToId = node2Idx;
                r.quadrantTo = ThermalNodeFace::NONE;
                r.type = HeatTransferType::CONDUCTION;
                
                double avgK = (_nodes[node1Idx].quadrants[0].thermalConductivity + 
                              _nodes[node2Idx].quadrants[0].thermalConductivity) / 2.0;
                double area = std::min(_nodes[node1Idx].getTotalSurfaceArea(),
                                      _nodes[node2Idx].getTotalSurfaceArea()) * 0.1;
                
                r.resistance = ThermalResistance::calculateConductionResistance(dist, avgK, area);
                _resistances.push_back(r);
            }
        }
    }
}

void Temperature::createBobbinConnections() {
    if (!hasBobbinNodes()) {
        return;
    }
    
    size_t bobbinIdx = 0;
    for (size_t i = 0; i < _nodes.size(); ++i) {
        if (_nodes[i].part == ThermalNodePartType::BOBBIN_CENTRAL_COLUMN) {
            bobbinIdx = i;
            break;
        }
    }
    
    for (size_t i = 0; i < _nodes.size(); ++i) {
        if (_nodes[i].part == ThermalNodePartType::CORE_CENTRAL_COLUMN ||
            _nodes[i].part == ThermalNodePartType::CORE_LATERAL_COLUMN) {
            
            ThermalResistanceElement r;
            r.nodeFromId = bobbinIdx;
            r.quadrantFrom = ThermalNodeFace::NONE;
            r.nodeToId = i;
            r.quadrantTo = ThermalNodeFace::NONE;
            r.type = HeatTransferType::CONDUCTION;
            
            double dist = 0.001;
            double area = _nodes[bobbinIdx].getTotalSurfaceArea() * 0.5;
            double bobbinK = 0.2;
            
            r.resistance = ThermalResistance::calculateConductionResistance(dist, bobbinK, area);
            _resistances.push_back(r);
        }
    }
}

// ============================================================================
// Connection Creation - Per Quadrant Logic
// ============================================================================

void Temperature::createTurnToTurnConnections() {
    if (!_isToroidal) {
        return;
    }
    
    double minConductionDist = getMinimumDistanceForConduction();
    
    std::vector<size_t> turnNodeIndices;
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part == ThermalNodePartType::TURN) {
            turnNodeIndices.push_back(i);
        }
    }
    
    for (size_t i = 0; i < turnNodeIndices.size(); i++) {
        size_t node1Idx = turnNodeIndices[i];
        const auto& node1 = _nodes[node1Idx];
        
        for (size_t j = i + 1; j < turnNodeIndices.size(); j++) {
            size_t node2Idx = turnNodeIndices[j];
            const auto& node2 = _nodes[node2Idx];
            
            double dx = node1.physicalCoordinates[0] - node2.physicalCoordinates[0];
            double dy = node1.physicalCoordinates[1] - node2.physicalCoordinates[1];
            double centerDistance = std::sqrt(dx*dx + dy*dy);
            
            double surfaceDistance = centerDistance - _wireWidth;
            
            if (surfaceDistance <= minConductionDist) {
                double angle1 = std::atan2(node1.physicalCoordinates[1], node1.physicalCoordinates[0]);
                double angle2 = std::atan2(node2.physicalCoordinates[1], node2.physicalCoordinates[0]);
                
                double angleDiff = angle1 - angle2;
                while (angleDiff > std::numbers::pi) angleDiff -= 2 * std::numbers::pi;
                while (angleDiff < -std::numbers::pi) angleDiff += 2 * std::numbers::pi;
                
                double r1 = std::sqrt(node1.physicalCoordinates[0]*node1.physicalCoordinates[0] + 
                                      node1.physicalCoordinates[1]*node1.physicalCoordinates[1]);
                double r2 = std::sqrt(node2.physicalCoordinates[0]*node2.physicalCoordinates[0] + 
                                      node2.physicalCoordinates[1]*node2.physicalCoordinates[1]);
                double radialDiff = std::abs(r1 - r2);
                
                // Tangential connection
                if (radialDiff < _wireWidth && std::abs(angleDiff) > 0.1) {
                    size_t leftNodeIdx, rightNodeIdx;
                    
                    if (angle1 > angle2) {
                        leftNodeIdx = node1Idx;
                        rightNodeIdx = node2Idx;
                    } else {
                        leftNodeIdx = node2Idx;
                        rightNodeIdx = node1Idx;
                    }
                    
                    ThermalResistanceElement r;
                    r.nodeFromId = leftNodeIdx;
                    r.quadrantFrom = ThermalNodeFace::TANGENTIAL_RIGHT;
                    r.nodeToId = rightNodeIdx;
                    r.quadrantTo = ThermalNodeFace::TANGENTIAL_LEFT;
                    r.type = HeatTransferType::CONDUCTION;
                    
                    auto* qLeft = _nodes[leftNodeIdx].getQuadrant(ThermalNodeFace::TANGENTIAL_RIGHT);
                    auto* qRight = _nodes[rightNodeIdx].getQuadrant(ThermalNodeFace::TANGENTIAL_LEFT);
                    
                    if (qLeft && qRight) {
                        double contactArea = std::min(qLeft->surfaceArea, qRight->surfaceArea) * 0.5;
                        int turn1Idx = _nodes[leftNodeIdx].turnIndex.value_or(0);
                        int turn2Idx = _nodes[rightNodeIdx].turnIndex.value_or(0);
                        
                        // Base resistance from insulation layers
                        double baseResistance = getInsulationLayerThermalResistance(turn1Idx, turn2Idx, contactArea);
                        
                        // Add coating resistances from both quadrants if they exist
                        if (qLeft->coating.has_value()) {
                            baseResistance += WireCoatingUtils::calculateCoatingResistance(qLeft->coating.value(), contactArea);
                        }
                        if (qRight->coating.has_value()) {
                            baseResistance += WireCoatingUtils::calculateCoatingResistance(qRight->coating.value(), contactArea);
                        }
                        
                        r.resistance = baseResistance;
                        
                        // Add additional inter-turn insulation from config if enabled
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
                
                // Radial connection
                if (radialDiff < _wireWidth * 1.5 && std::abs(angleDiff) < 0.3) {
                    size_t innerNodeIdx, outerNodeIdx;
                    
                    if (r1 < r2) {
                        innerNodeIdx = node1Idx;
                        outerNodeIdx = node2Idx;
                    } else {
                        innerNodeIdx = node2Idx;
                        outerNodeIdx = node1Idx;
                    }
                    
                    ThermalResistanceElement r;
                    r.nodeFromId = innerNodeIdx;
                    r.quadrantFrom = ThermalNodeFace::RADIAL_OUTER;
                    r.nodeToId = outerNodeIdx;
                    r.quadrantTo = ThermalNodeFace::RADIAL_INNER;
                    r.type = HeatTransferType::CONDUCTION;
                    
                    auto* qInner = _nodes[innerNodeIdx].getQuadrant(ThermalNodeFace::RADIAL_OUTER);
                    auto* qOuter = _nodes[outerNodeIdx].getQuadrant(ThermalNodeFace::RADIAL_INNER);
                    
                    if (qInner && qOuter) {
                        double contactArea = std::min(qInner->surfaceArea, qOuter->surfaceArea) * 0.5;
                        int turn1Idx = _nodes[innerNodeIdx].turnIndex.value_or(0);
                        int turn2Idx = _nodes[outerNodeIdx].turnIndex.value_or(0);
                        
                        // Base resistance from insulation layers
                        double baseResistance = getInsulationLayerThermalResistance(turn1Idx, turn2Idx, contactArea);
                        
                        // Add coating resistances from both quadrants
                        if (qInner->coating.has_value()) {
                            baseResistance += WireCoatingUtils::calculateCoatingResistance(qInner->coating.value(), contactArea);
                        }
                        if (qOuter->coating.has_value()) {
                            baseResistance += WireCoatingUtils::calculateCoatingResistance(qOuter->coating.value(), contactArea);
                        }
                        
                        r.resistance = baseResistance;
                        
                        // Add additional inter-turn insulation from config if enabled
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
                    // Node is at wire surface, so conduction length is half the distance to opposite surface
                    double copperLength = _wireWidth / 2.0;
                    double copperResistance = ThermalResistance::calculateConductionResistance(
                        copperLength, _wireThermalCond, qOuter->surfaceArea);
                    
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
                    double copperLength = _wireWidth / 2.0;
                    double copperResistance = ThermalResistance::calculateConductionResistance(
                        copperLength, _wireThermalCond, qInner->surfaceArea);
                    
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
    
    double maxConvectionDist = getMaximumDistanceForConvection();
    double minConductionDist = getMinimumDistanceForConduction();
    
    // Get convection coefficient
    double surfaceTemp = _config.ambientTemperature + 30.0;
    double h_conv;
    
    if (_config.includeForcedConvection) {
        h_conv = ThermalResistance::calculateForcedConvectionCoefficient(
            _config.airVelocity, _wireWidth, _config.ambientTemperature);
    } else {
        h_conv = ThermalResistance::calculateNaturalConvectionCoefficient(
            surfaceTemp, _config.ambientTemperature, _wireWidth, SurfaceOrientation::VERTICAL);
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
            
            // Check each quadrant
            for (int qIdx = 0; qIdx < 4; ++qIdx) {
                ThermalNodeFace face = node.quadrants[qIdx].face;
                if (face == ThermalNodeFace::NONE) continue;
                
                // Skip if this quadrant is already connected by conduction
                std::string qKey = std::to_string(i) + "_" + std::string(magic_enum::enum_name(face));
                if (connectedQuadrants.count(qKey) > 0) continue;
                
                bool isExposed = true;
                
                // Check for blocking objects in the quadrant's direction
                if (face == ThermalNodeFace::RADIAL_INNER) {
                    // Check for turns or core closer to center
                    for (size_t j = 0; j < _nodes.size(); j++) {
                        if (i == j) continue;
                        
                        double otherX = _nodes[j].physicalCoordinates[0];
                        double otherY = _nodes[j].physicalCoordinates[1];
                        double otherR = std::sqrt(otherX*otherX + otherY*otherY);
                        double otherAngle = std::atan2(otherY, otherX);
                        
                        // Check if other object is in the radial inward direction
                        if (otherR < nodeR) {
                            double angleDiff = std::abs(nodeAngle - otherAngle);
                            while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                            
                            // Roughly same angular position and within max convection distance
                            if (std::abs(angleDiff) < 0.3 && (nodeR - otherR) < maxConvectionDist) {
                                isExposed = false;
                                break;
                            }
                        }
                    }
                }
                else if (face == ThermalNodeFace::RADIAL_OUTER) {
                    // Check for turns or core farther from center
                    for (size_t j = 0; j < _nodes.size(); j++) {
                        if (i == j) continue;
                        
                        double otherX = _nodes[j].physicalCoordinates[0];
                        double otherY = _nodes[j].physicalCoordinates[1];
                        double otherR = std::sqrt(otherX*otherX + otherY*otherY);
                        double otherAngle = std::atan2(otherY, otherX);
                        
                        if (otherR > nodeR) {
                            double angleDiff = std::abs(nodeAngle - otherAngle);
                            while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                            
                            if (std::abs(angleDiff) < 0.3 && (otherR - nodeR) < maxConvectionDist) {
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
                        if (std::abs(otherR - nodeR) < _wireWidth) {
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
                
                if (std::abs(turnR - coreOuterR) < _wireWidth) {
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
                
                if (std::abs(turnR - coreInnerR) < _wireWidth) {
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
        // Concentric core - simplified convection for all exposed surfaces
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
    return dist <= getMinimumDistanceForConduction();
}

double Temperature::calculateContactArea(const ThermalNodeQuadrant& q1, const ThermalNodeQuadrant& q2) const {
    double minLength = std::min(q1.length, q2.length);
    return _wireHeight * minLength;
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
