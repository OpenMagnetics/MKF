#include "physical_models/Temperature.h"
#include "constructive_models/CorePiece.h"
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
using namespace ThermalDefaults; // IMP-2


// Debug flag for thermal model - set to true for verbose output
constexpr bool THERMAL_DEBUG = false;

// Contact threshold: surfaces must be within this distance to conduct
constexpr double CONTACT_THRESHOLD_FACTOR = 0.25;  // wireDiameter / 4

// ============================================================================
// Wire Property Helper Functions (extract properties from Wire object locally)
// ============================================================================

// Dense matrix solver for the thermal circuit. File-local — not part of the
// public Temperature API. For networks >100 nodes, consider a sparse solver.
namespace {
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
    void setZero() { for (auto& row : data) std::fill(row.begin(), row.end(), 0.0); }
    void setRowZero(size_t row) { std::fill(data[row].begin(), data[row].end(), 0.0); }
    void setColZero(size_t col) { for (size_t i = 0; i < rows_; ++i) data[i][col] = 0.0; }

    // Non-throwing solve with success flag
    static std::vector<double> solve(SimpleMatrix A, std::vector<double> b, bool& success) {
        success = true;
        size_t n = A.rows();
        if (n == 0 || b.size() != n) { success = false; return std::vector<double>(n, 0.0); }
        std::vector<std::vector<double>> aug(n, std::vector<double>(n + 1));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) aug[i][j] = A(i, j);
            aug[i][n] = b[i];
        }
        for (size_t col = 0; col < n; ++col) {
            size_t maxRow = col;
            double maxVal = std::abs(aug[col][col]);
            for (size_t row = col + 1; row < n; ++row) {
                if (std::abs(aug[row][col]) > maxVal) { maxVal = std::abs(aug[row][col]); maxRow = row; }
            }
            if (maxRow != col) std::swap(aug[col], aug[maxRow]);
            if (std::abs(aug[col][col]) < 1e-15) { success = false; return std::vector<double>(n, 0.0); }
            for (size_t row = col + 1; row < n; ++row) {
                double factor = aug[row][col] / aug[col][col];
                for (size_t j = col; j <= n; ++j) aug[row][j] -= factor * aug[col][j];
            }
        }
        std::vector<double> x(n);
        for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
            x[i] = aug[i][n];
            for (size_t j = i + 1; j < n; ++j) x[i] -= aug[i][j] * x[j];
            x[i] /= aug[i][i];
        }
        return x;
    }
};
} // anonymous namespace

// ============================================================================
// CoolingUtils Implementation
// ============================================================================

CoolingType CoolingUtils::detectCoolingType(const MAS::Cooling& cooling) {
    // Cold plate: has maximum_temperature
    if (cooling.get_maximum_temperature().has_value()) {
        return CoolingType::COLD_PLATE;
    }
    
    // Forced convection: has velocity
    if (cooling.get_velocity().has_value() && 
        !cooling.get_velocity().value().empty()) {
        return CoolingType::FORCED_CONVECTION;
    }
    
    // Heatsink: has thermal_resistance (and no maximum_temperature)
    if (cooling.get_thermal_resistance().has_value()) {
        return CoolingType::HEATSINK;
    }
    
    // Natural convection: has temperature but no velocity/thermal_resistance/max_temp
    if (cooling.get_temperature().has_value()) {
        return CoolingType::NATURAL_CONVECTION;
    }
    
    return CoolingType::UNKNOWN;
}

bool CoolingUtils::isNaturalConvection(const MAS::Cooling& cooling) {
    return detectCoolingType(cooling) == CoolingType::NATURAL_CONVECTION;
}

bool CoolingUtils::isForcedConvection(const MAS::Cooling& cooling) {
    return detectCoolingType(cooling) == CoolingType::FORCED_CONVECTION;
}

bool CoolingUtils::isHeatsink(const MAS::Cooling& cooling) {
    return detectCoolingType(cooling) == CoolingType::HEATSINK;
}

bool CoolingUtils::isColdPlate(const MAS::Cooling& cooling) {
    return detectCoolingType(cooling) == CoolingType::COLD_PLATE;
}

double CoolingUtils::calculateForcedConvectionCoefficient(
    double surfaceTemp,
    double ambientTemp,
    double velocity,
    double characteristicLength,
    double fluidConductivity,
    double kinematicViscosity,
    double prandtlNumber) {
    
    // Calculate Reynolds number: Re = V * L / nu
    double reynolds = velocity * characteristicLength / kinematicViscosity;
    
    // For flat plate laminar flow (Re < 5e5): Nu = 0.664 * Re^0.5 * Pr^(1/3)
    // For turbulent flow (Re >= 5e5): Nu = 0.037 * Re^0.8 * Pr^(1/3)
    double nusselt;
    if (reynolds < 5e5) {
        nusselt = 0.664 * std::sqrt(reynolds) * std::cbrt(prandtlNumber);
    } else {
        nusselt = 0.037 * std::pow(reynolds, 0.8) * std::cbrt(prandtlNumber);
    }
    
    // h = Nu * k / L
    return nusselt * fluidConductivity / characteristicLength;
}

double CoolingUtils::calculateMixedConvectionCoefficient(
    double h_natural,
    double h_forced) {
    
    // Mixed convection: h_total = (h_natural^3 + h_forced^3)^(1/3)
    return std::cbrt(std::pow(h_natural, 3) + std::pow(h_forced, 3));
}

TemperatureConfig TemperatureConfig::fromMasOperatingConditions(
    const MAS::OperatingConditions& conditions) {
    
    TemperatureConfig config;
    config.ambientTemperature = conditions.get_ambient_temperature();
    
    if (conditions.get_cooling().has_value()) {
        config.masCooling = conditions.get_cooling().value();
    }
    
    return config;
}

// ============================================================================
// Material Properties
// ============================================================================

double Temperature::getCoreThermalConductivity() const {
    try {
        auto coreMaterial = Core(_magnetic.get_core()).resolve_material();
        return ThermalResistance::getCoreMaterialThermalConductivity(coreMaterial);
    } catch (...) {
        return _config.coreThermalConductivity;
    }
}

double Temperature::getBobbinThermalConductivity() const {
    try {
        auto bobbinOpt = _magnetic.get_coil().get_bobbin();
        if (std::holds_alternative<Bobbin>(bobbinOpt)) {
            auto bobbin = std::get<Bobbin>(bobbinOpt);
            auto funcDesc = bobbin.get_functional_description();
            if (funcDesc) {
                auto materialOpt = funcDesc->get_material();
                if (materialOpt) {
                    auto& materialVariant = *materialOpt;
                    if (std::holds_alternative<MAS::InsulationMaterial>(materialVariant)) {
                        OpenMagnetics::InsulationMaterial resolved(
                            std::get<MAS::InsulationMaterial>(materialVariant));
                        return ThermalResistance::getInsulationMaterialThermalConductivity(resolved);
                    } else if (std::holds_alternative<std::string>(materialVariant)) {
                        auto material = find_insulation_material_by_name(
                            std::get<std::string>(materialVariant));
                        return ThermalResistance::getInsulationMaterialThermalConductivity(material);
                    }
                }
            }
        }
        // No material on bobbin: use default bobbin material (PET)
        auto defaultMaterial = find_insulation_material_by_name(Defaults().defaultBobbinMaterial);
        return ThermalResistance::getInsulationMaterialThermalConductivity(defaultMaterial);
    } catch (...) {
        return kBobbin_ThermalConductivity;
    }
}

// ============================================================================
// Main Analysis Entry Point
// ============================================================================

ThermalResult Temperature::calculateTemperatures() {
    if (THERMAL_DEBUG) {
    }
    
    // Step 1: Extract wire properties from coil
    extractWireProperties();
    
    // Step 2: Create thermal nodes (core, turns, bobbin)
    createThermalNodes();
    
    // Step 3: Create thermal resistances between nodes
    createThermalResistances();
    
    // Step 4: Calculate schematic scaling
    calculateSchematicScaling();
    
    // Step 5: Plot schematic if requested (before solving, so it's available even if solver fails)
    if (_config.plotSchematic) {
        plotSchematic();
    }
    
    // Step 6: Create and solve thermal equivalent circuit
    ThermalResult result = solveThermalCircuit();
    
    if (THERMAL_DEBUG) {
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
        throw std::runtime_error("Temperature::extractWireProperties: No windings found in coil. "
                                 "Cannot extract wire properties for thermal analysis.");
    }
    
    auto wireVariant = windings[0].get_wire();
    if (!std::holds_alternative<Wire>(wireVariant)) {
        throw std::runtime_error("Temperature::extractWireProperties: Wire in winding[0] is not a resolved Wire object. "
                                 "Wire must be resolved before thermal analysis.");
    }
    
    auto wire = std::get<Wire>(wireVariant);
    _isRoundWire = (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ);
    _isPlanar = (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::FOIL);
    
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
    
    // For rectangular/foil wires, query actual width and height
    if (wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL || wire.get_type() == WireType::PLANAR) {
        // Prefer conducting_width / conducting_height if available
        if (wire.get_conducting_width() && wire.get_conducting_height()) {
            _wireWidth = resolve_dimensional_values(wire.get_conducting_width().value());
            _wireHeight = resolve_dimensional_values(wire.get_conducting_height().value());
        } else if (wire.get_outer_width() && wire.get_outer_height()) {
            _wireWidth = resolve_dimensional_values(wire.get_outer_width().value());
            _wireHeight = resolve_dimensional_values(wire.get_outer_height().value());
        } else if (wire.get_outer_diameter()) {
            auto outerDiam = wire.get_outer_diameter().value();
            if (outerDiam.get_nominal()) {
                _wireWidth = outerDiam.get_nominal().value();
                // Fallback: use outer_diameter for both if no separate width/height
                _wireHeight = _wireWidth;
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
            } catch (const std::exception& e) {
                throw std::runtime_error("Temperature::extractWireProperties: Failed to lookup wire material '" + 
                                         materialName + "' for thermal conductivity: " + e.what());
            }
        }
    }
    
    // Get wire coating for thermal calculations
    _wireCoating = wire.resolve_coating();
    
    // =========================================================================
    // IMP-NEW-06: Populate per-winding wire properties from EACH winding's wire
    // =========================================================================
    for (size_t wIdx = 0; wIdx < windings.size(); ++wIdx) {
        WindingWireProperties wProps;
        auto wv = windings[wIdx].get_wire();
        if (std::holds_alternative<Wire>(wv)) {
            auto w = std::get<Wire>(wv);
            wProps.isRoundWire = (w.get_type() == WireType::ROUND || w.get_type() == WireType::LITZ);
            wProps.isPlanar = (w.get_type() == WireType::PLANAR || w.get_type() == WireType::FOIL);
            
            // Extract dimensions from THIS winding's wire (not winding[0])
            wProps.wireWidth = kWire_DefaultWidth;
            wProps.wireHeight = kWire_DefaultHeight;
            
            if (wProps.isRoundWire) {
                if (w.get_conducting_diameter()) {
                    auto condDiam = w.get_conducting_diameter().value();
                    if (condDiam.get_nominal()) {
                        wProps.wireWidth = condDiam.get_nominal().value();
                        wProps.wireHeight = wProps.wireWidth;
                    }
                } else if (w.get_outer_diameter()) {
                    auto outerDiam = w.get_outer_diameter().value();
                    if (outerDiam.get_nominal()) {
                        wProps.wireWidth = outerDiam.get_nominal().value();
                        wProps.wireHeight = wProps.wireWidth;
                    }
                }
            } else if (w.get_type() == WireType::RECTANGULAR || w.get_type() == WireType::FOIL || w.get_type() == WireType::PLANAR) {
                if (w.get_conducting_width() && w.get_conducting_height()) {
                    wProps.wireWidth = resolve_dimensional_values(w.get_conducting_width().value());
                    wProps.wireHeight = resolve_dimensional_values(w.get_conducting_height().value());
                } else if (w.get_outer_width() && w.get_outer_height()) {
                    wProps.wireWidth = resolve_dimensional_values(w.get_outer_width().value());
                    wProps.wireHeight = resolve_dimensional_values(w.get_outer_height().value());
                }
            }
            
            // Validate that wire dimensions were actually extracted
            if (wProps.wireWidth == kWire_DefaultWidth && wProps.wireHeight == kWire_DefaultHeight) {
                throw std::runtime_error("Temperature::extractWireProperties: Could not extract wire dimensions for winding " +
                                         std::to_string(wIdx) + ". Wire must have conducting_diameter (round/litz) or " +
                                         "conducting_width/height (rectangular/foil/planar) defined.");
            }
            
            // Extract thermal conductivity from THIS winding's wire material
            wProps.wireThermalCond = kWire_CopperThermalConductivity;
            if (w.get_material()) {
                auto materialVariant = w.get_material().value();
                if (std::holds_alternative<std::string>(materialVariant)) {
                    try {
                        auto wireMaterial = find_wire_material_by_name(std::get<std::string>(materialVariant));
                        auto thermalCond = wireMaterial.get_thermal_conductivity();
                        if (thermalCond && !thermalCond->empty()) {
                            wProps.wireThermalCond = (*thermalCond)[0].get_value();
                        }
                    } catch (const std::exception& e) {
                        throw std::runtime_error("Temperature::extractWireProperties: Failed to lookup wire material for winding " + 
                                                 std::to_string(wIdx) + ": " + e.what());
                    }
                }
            }
            
            wProps.wireCoating = w.resolve_coating();
            _perWindingWireProps[wIdx] = wProps;
        }
    }
    
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
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("Temperature::createThermalNodes: Failed to wind coil: " + std::string(e.what()));
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
}

void Temperature::createToroidalCoreNodes() {
    auto core = _magnetic.get_core();
    auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    
    if (dimensions.find("A") == dimensions.end() || 
        dimensions.find("B") == dimensions.end() || 
        dimensions.find("C") == dimensions.end()) {
        throw std::runtime_error("Toroidal core dimensions must include keys A (outer diameter), "
                                 "B (inner diameter), and C (height)");
    }
    double outerDiameter = dimensions["A"];
    double innerDiameter = dimensions["B"];
    double height = dimensions["C"];
    
    double windingWindowInnerRadius = innerDiameter / 2.0;
    double windingWindowOuterRadius = outerDiameter / 2.0;
    double meanRadius = (windingWindowInnerRadius + windingWindowOuterRadius) / 2.0;
    
    size_t numSegments = _config.toroidalSegments;
    double angularStep = 2.0 * M_PI / numSegments;
    
    double coreK = getCoreThermalConductivity();
    
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
                double innerWidth = _wireWidth;  // Use cached wire width
                double outerWidth = _wireWidth;
                if (turn.get_dimensions() && turn.get_dimensions()->size() >= 1) {
                    innerWidth = (*turn.get_dimensions())[0];
                    outerWidth = innerWidth;
                }
                
                turnCoverageInfo.push_back({turnAngle, innerWidth, outerWidth, isInner});
            }
        }
    }
    
    // Pre-calculate loss weighting proportional to 1/r^2 (B ~ 1/r, losses ~ B^2)
    // Each segment's mean radius determines its share of total core losses
    double lossWeightSum = 0.0;
    std::vector<double> segmentLossWeights(numSegments);
    for (size_t i = 0; i < numSegments; ++i) {
        // All segments have same mean radius in toroidal azimuthal segmentation,
        // so losses are equal per segment (radial variation is within each segment)
        // For azimuthal segmentation: the loss distribution is uniform around the ring
        segmentLossWeights[i] = 1.0;
        lossWeightSum += segmentLossWeights[i];
    }
    
    for (size_t i = 0; i < numSegments; ++i) {
        double angle = i * angularStep;
        double midAngle = angle + angularStep / 2.0;
        
        ThermalNetworkNode node;
        node.part = ThermalNodePartType::CORE_TOROIDAL_SEGMENT;
        node.name = "Core_Segment_" + std::to_string(i);
        node.temperature = _config.ambientTemperature;
        // Distribute losses proportionally (currently uniform for azimuthal segments,
        // but framework supports non-uniform distribution if segments are radially split)
        node.powerDissipation = _config.coreLosses * (segmentLossWeights[i] / lossWeightSum);
        
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
    
    if (columns.empty()) {
        throw std::runtime_error("Temperature::createConcentricCoreNodes: Core has no columns. "
                                 "Cannot build thermal model without core column geometry.");
    }
    if (!processedDesc) {
        throw std::runtime_error("Temperature::createConcentricCoreNodes: Core is missing processed description. "
                                 "Cannot build thermal model without core dimensions.");
    }
    
    double coreK = getCoreThermalConductivity();
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
    

    // =========================================================================
    // Volume-proportional loss distribution (replaces hardcoded 40/20/10/10)
    // =========================================================================
    auto corePiece = CorePiece::factory(core.resolve_shape());
    auto lossFractions = corePiece->calculate_core_loss_fractions();

    double centralColumnLosses = _config.coreLosses * lossFractions.centralColumn;
    double lateralColumnLosses = _config.coreLosses * lossFractions.lateralColumn;
    double topYokeLosses       = _config.coreLosses * lossFractions.topYoke;
    double bottomYokeLosses    = _config.coreLosses * lossFractions.bottomYoke;

    // Core loss distribution calculated
    // =========================================================================

    // =========================================================================
    // IMP-NEW-12: HALF-CORE SYMMETRY MODEL
    // =========================================================================
    // This model uses HALF the core depth for concentric cores (left-right symmetry).
    // - Node volumes are halved -> power dissipation per node is correct (same W/m3)
    // - Surface areas are halved -> convection area represents one half only
    // - The symmetry correction is applied by doubling all convection surface areas
    //   to account for the other symmetric half
    // =========================================================================
    // Create central column node(s) - using HALF depth for symmetry (left half only)
    double halfDepth = mainColumn.get_depth() / 2.0;  // Model half the core
    // Always create a single central column node with standard name for painter compatibility
    // (Chunked nodes would require painter support for matching)
    ThermalNetworkNode node;
    node.part = ThermalNodePartType::CORE_CENTRAL_COLUMN;
    node.name = "Core_Column_0";
    node.temperature = _config.ambientTemperature;
    node.powerDissipation = centralColumnLosses; // Volume-proportional
    node.physicalCoordinates = {0, 0, 0};
    node.initializeConcentricCoreQuadrants(mainColumn.get_width(), 
                                            mainColumn.get_height(), halfDepth, coreK);
    _nodes.push_back(node);
    
    // Create top yoke node - using HALF depth for symmetry
    // Position between center (x=0) and right lateral column (x=coreWidth/2)
    ThermalNetworkNode topYoke;
    topYoke.part = ThermalNodePartType::CORE_TOP_YOKE;
    topYoke.name = "Core_Top_Yoke";
    topYoke.temperature = _config.ambientTemperature;
    topYoke.powerDissipation = topYokeLosses; // Volume-proportional
    topYoke.physicalCoordinates = {coreWidth / 4, coreHeight / 2 - mainColumn.get_width()/4, 0};
    topYoke.initializeConcentricCoreQuadrants(coreWidth / 2, mainColumn.get_width()/2, coreDepth / 2, coreK);
    _nodes.push_back(topYoke);
    
    // Create bottom yoke node - using HALF depth for symmetry
    // Position between center (x=0) and right lateral column (x=coreWidth/2)
    ThermalNetworkNode bottomYoke;
    bottomYoke.part = ThermalNodePartType::CORE_BOTTOM_YOKE;
    bottomYoke.name = "Core_Bottom_Yoke";
    bottomYoke.temperature = _config.ambientTemperature;
    bottomYoke.powerDissipation = bottomYokeLosses; // Volume-proportional
    bottomYoke.physicalCoordinates = {coreWidth / 4, -coreHeight / 2 + mainColumn.get_width()/4, 0};
    bottomYoke.initializeConcentricCoreQuadrants(coreWidth / 2, mainColumn.get_width()/2, coreDepth / 2, coreK);
    _nodes.push_back(bottomYoke);
    
    // Create lateral column node (RIGHT side only for symmetry - half core model)
    if (columns.size() > 1) {
        double latHalfDepth = rightColumn.get_depth() / 2.0; // Model half the core
        
        double offset = coreWidth / 2 - rightColumn.get_width() / 2;
        
        // Always create a single lateral column node with standard name for painter compatibility
        ThermalNetworkNode rightNode;
        rightNode.part = ThermalNodePartType::CORE_LATERAL_COLUMN;
        rightNode.name = "Core_Column_1";
        rightNode.temperature = _config.ambientTemperature;
        rightNode.powerDissipation = lateralColumnLosses; // Volume-proportional
        rightNode.physicalCoordinates = {offset, 0, 0};
        rightNode.initializeConcentricCoreQuadrants(rightColumn.get_width(), rightColumn.get_height(), latHalfDepth, coreK);
        _nodes.push_back(rightNode);
    }
}

void Temperature::createBobbinNodes() {
    auto bobbinOpt = _magnetic.get_mutable_coil().get_bobbin();
    if (!std::holds_alternative<Bobbin>(bobbinOpt)) {
        throw std::runtime_error("Temperature::createBobbinNodes: Bobbin is not resolved to a Bobbin object. "
                                 "Bobbin reference must be resolved before thermal analysis.");
    }
    
    if (_isToroidal) {
        return;
    }
    
    auto bobbin = std::get<Bobbin>(bobbinOpt);
    auto processedDesc = bobbin.get_processed_description();
    if (!processedDesc) {
        throw std::runtime_error("Temperature::createBobbinNodes: Bobbin is missing processed description. "
                                 "Cannot build thermal model without bobbin geometry.");
    }
    
    double bobbinK = getBobbinThermalConductivity();
    double wallThickness = processedDesc->get_wall_thickness();
    
    // Only create bobbin nodes if there's actual wall thickness
    if (wallThickness <= 0) {
        return;
    }
    
    auto core = _magnetic.get_core();
    auto processedCore = core.get_processed_description();
    if (!processedCore) {
        throw std::runtime_error("Temperature::createBobbinNodes: Core is missing processed description. "
                                 "Cannot build thermal model without core dimensions.");
    }
    
    double coreWidth = processedCore->get_width();
    double coreDepth = processedCore->get_depth();
    
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    double columnWidth = mainColumn.get_width();
    
    // Get winding window for proper bobbin positioning
    auto windingWindows = core.get_winding_windows();
    if (windingWindows.empty()) {
        throw std::runtime_error("Temperature::createBobbinNodes: Core has no winding windows. "
                                 "Cannot create bobbin nodes without winding window geometry.");
    }
    if (!windingWindows[0].get_width().has_value() || !windingWindows[0].get_height().has_value()) {
        throw std::runtime_error("Temperature::createBobbinNodes: Winding window is missing width or height. "
                                 "Both dimensions must be specified for thermal analysis.");
    }
    double windingWindowWidth = windingWindows[0].get_width().value();
    double windingWindowHeight = windingWindows[0].get_height().value();
    if (windingWindowWidth <= 0 || windingWindowHeight <= 0) {
        throw std::runtime_error("Temperature::createBobbinNodes: Winding window has zero or negative dimensions. "
                                 "Width=" + std::to_string(windingWindowWidth) + ", Height=" + std::to_string(windingWindowHeight));
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
    auto coil = _magnetic.get_mutable_coil();
    
    // Check if we have layers description
    if (!coil.get_layers_description()) {
        if (THERMAL_DEBUG) {
        }
        return;
    }
    
    // Get all insulation layers
    auto insulationLayers = coil.get_layers_description_insulation();
    if (insulationLayers.empty()) {
        if (THERMAL_DEBUG) {
        }
        return;
    }
    
    if (THERMAL_DEBUG) {
    }
    
    // Get core depth for insulation layer depth
    auto core = _magnetic.get_core();
    auto processedCore = core.get_processed_description();
    if (!processedCore) {
        throw std::runtime_error("Temperature::createInsulationLayerNodes: Core is missing processed description. "
                                 "Cannot determine core depth for insulation layer geometry.");
    }
    double coreDepth = processedCore->get_depth();
    
    size_t layerIdx = 0;
    size_t createdCount = 0;
    
    if (THERMAL_DEBUG) {
    }
    
    for (const auto& layer : insulationLayers) {
        if (THERMAL_DEBUG) {
        }
        // Skip layers without proper dimensions or coordinates
        if (layer.get_dimensions().size() < 2) {
            continue;
        }
        if (layer.get_coordinates().size() < 2) {
            continue;
        }
        
        // Note: dimensions[1] is angular span, which can be in degrees or millidegrees depending on source
        double angularSpan = layer.get_dimensions()[1];
        // If it's > 1000, it's probably millidegrees (360000 = 360°), otherwise it's degrees (360 = 360°)
        double angularSpanDegrees = (angularSpan > 1000.0) ? (angularSpan / 1000.0) : angularSpan;
        bool isFullCircle = (angularSpanDegrees >= 359.0);
        
        // Get layer center coordinates
        double layerX = layer.get_coordinates()[0];
        double layerY = layer.get_coordinates()[1];
        
        // Get layer geometry
        double layerWidth = layer.get_dimensions()[0];   // X dimension (thickness)
        double layerHeight = layer.get_dimensions()[1];  // Y dimension (span)
        double layerDepth = coreDepth / 2.0;  // Half depth for single side modeling
        
        // For concentric cores, if insulation layer has zero thickness, use a default
        // based on typical inter-layer insulation or coil's insulation specification
        if (!_isToroidal && layerWidth < 1e-9) {
            // Try to get thickness from coil's insulation layer specification
            try {
                double specifiedThickness = coil.get_insulation_layer_thickness(layer);
                if (specifiedThickness > 1e-9) {
                    layerWidth = specifiedThickness;
                } else {
                    layerWidth = kInsulation_DefaultThickness;
                }
            } catch (const std::exception& e) { // IMP-8
                layerWidth = kInsulation_DefaultThickness;
                if (THERMAL_DEBUG) {
                    std::cout << "Temperature: insulation-layer thickness query failed, using default: "
                              << e.what() << std::endl;
                }
            }
        }
        
        // Get layer thermal conductivity from material
        double layerK = kInsulation_DefaultConductivity;
        try {
            auto insulationMaterial = coil.resolve_insulation_layer_insulation_material(layer);
            if (insulationMaterial.get_thermal_conductivity()) {
                layerK = insulationMaterial.get_thermal_conductivity().value();
            }
        } catch (const std::exception& e) { // IMP-8
            // Use default if material cannot be resolved
            if (THERMAL_DEBUG) {
                std::cout << "Temperature: insulation material resolution failed, using default k: "
                          << e.what() << std::endl;
            }
        }
        
        // Check coordinate system
        auto coordSystem = layer.get_coordinate_system();
        bool isPolar = coordSystem.has_value() && 
                       coordSystem.value() == MAS::CoordinateSystem::POLAR;
        
        if (_isToroidal && isPolar) {
            // For toroidal cores with polar coordinates, chunk the insulation layer
            // into segments matching the core's angular segmentation
            // Each segment gets INNER and OUTER nodes (like turns)
            
            // layerX = radial height (distance from winding window inner surface to layer center)
            // layerY = angle in degrees (center angle of the layer)
            // layerWidth = radial thickness  
            // layerHeight = angular span in millidegrees
            
            double radialHeight = layerX;  // Distance from winding window inner surface to layer center
            double centerAngleDeg = layerY;
            double radialThickness = layerWidth;
            double angularSpanDeg = layerHeight / 1000.0;  // Convert from millidegrees

            if (radialThickness < 1e-9) {
                // Zero-thickness = geometric boundary marker from coil winder, not physical insulation.
                // Skip node creation; turns will connect directly across this boundary.
                layerIdx++;
                continue;
            }

            // Get core dimensions
            auto coreDims = flatten_dimensions(core.resolve_shape().get_dimensions().value());
            double innerDiameter = coreDims["B"];  // Inner hole diameter
            double outerDiameter = coreDims["A"];  // Outer core diameter
            
            // For insulation layers wrapping around the core:
            // - Inner nodes: inside the hole, toward center (B/2 - radialHeight)
            // - Outer nodes: outside the core, away from center (A/2 + radialHeight)
            double innerSurfaceRadius = innerDiameter / 2.0 - radialHeight;
            double outerSurfaceRadius = outerDiameter / 2.0 + radialHeight;
            
            // Use same number of segments as toroidal core
            size_t numSegments = _config.toroidalSegments;
            double angularStep = 2.0 * M_PI / numSegments;
            
            // Calculate which segments this insulation layer covers
            double layerStartAngle = (centerAngleDeg - angularSpanDeg / 2.0) * M_PI / 180.0;
            double layerEndAngle = (centerAngleDeg + angularSpanDeg / 2.0) * M_PI / 180.0;
            
            // Normalize angles to [0, 2π)
            auto normalizeAngle = [](double a) {
                while (a < 0) a += 2.0 * M_PI;
                while (a >= 2.0 * M_PI) a -= 2.0 * M_PI;
                return a;
            };
            layerStartAngle = normalizeAngle(layerStartAngle);
            layerEndAngle = normalizeAngle(layerEndAngle);
            
            // Create INNER and OUTER nodes for each overlapping segment
            for (size_t segIdx = 0; segIdx < numSegments; ++segIdx) {
                double segAngle = segIdx * angularStep;
                double segMidAngle = segAngle + angularStep / 2.0;
                double segEnd = segAngle + angularStep;
                
                // Check overlap (full-circle layers overlap with all segments)
                bool overlaps = isFullCircle;
                if (!overlaps) {
                    if (layerStartAngle <= layerEndAngle) {
                        overlaps = (segAngle < layerEndAngle && segEnd > layerStartAngle);
                    } else {
                        overlaps = (segAngle < layerEndAngle || segEnd > layerStartAngle);
                    }
                }
                
                if (!overlaps) continue;
                
                // Calculate segment arc length (at average radius)
                double avgRadius = (innerSurfaceRadius + outerSurfaceRadius) / 2.0;
                double segmentArcLength = avgRadius * angularStep;
                
                // Create INNER node (at inner surface of insulation layer)
                ThermalNetworkNode innerNode;
                innerNode.part = ThermalNodePartType::INSULATION_LAYER;
                innerNode.name = "IL_" + std::to_string(layerIdx) + "_" + std::to_string(segIdx) + "_i";
                innerNode.insulationLayerIndex = layerIdx;
                innerNode.temperature = _config.ambientTemperature;
                innerNode.powerDissipation = 0.0;
                innerNode.isInnerTurn = true;  // Mark as inner surface
                
                double innerX = innerSurfaceRadius * std::cos(segMidAngle);
                double innerY = innerSurfaceRadius * std::sin(segMidAngle);
                innerNode.physicalCoordinates = {innerX, innerY, 0};
                innerNode.initializeInsulationLayerQuadrants(radialThickness, segmentArcLength, layerDepth, layerK);
                
                // Override tangential limit coordinates to match segment boundaries
                // (generic ±90° perpendicular from initializeInsulationLayerQuadrants doesn't
                //  align with angular segment edges for toroidal layout)
                innerNode.quadrants[2].limitCoordinates = {
                    innerSurfaceRadius * std::cos(segAngle),
                    innerSurfaceRadius * std::sin(segAngle),
                    0.0
                };
                innerNode.quadrants[3].limitCoordinates = {
                    innerSurfaceRadius * std::cos(segEnd),
                    innerSurfaceRadius * std::sin(segEnd),
                    0.0
                };
                
                _nodes.push_back(innerNode);
                createdCount++;
                
                // Create OUTER node (at outer surface of insulation layer)
                ThermalNetworkNode outerNode;
                outerNode.part = ThermalNodePartType::INSULATION_LAYER;
                outerNode.name = "IL_" + std::to_string(layerIdx) + "_" + std::to_string(segIdx) + "_o";
                outerNode.insulationLayerIndex = layerIdx;
                outerNode.temperature = _config.ambientTemperature;
                outerNode.powerDissipation = 0.0;
                outerNode.isInnerTurn = false;  // Mark as outer surface
                
                double outerX = outerSurfaceRadius * std::cos(segMidAngle);
                double outerY = outerSurfaceRadius * std::sin(segMidAngle);
                outerNode.physicalCoordinates = {outerX, outerY, 0};
                outerNode.initializeInsulationLayerQuadrants(radialThickness, segmentArcLength, layerDepth, layerK);
                
                // Override tangential limit coordinates to match segment boundaries
                outerNode.quadrants[2].limitCoordinates = {
                    outerSurfaceRadius * std::cos(segAngle),
                    outerSurfaceRadius * std::sin(segAngle),
                    0.0
                };
                outerNode.quadrants[3].limitCoordinates = {
                    outerSurfaceRadius * std::cos(segEnd),
                    outerSurfaceRadius * std::sin(segEnd),
                    0.0
                };
                
                _nodes.push_back(outerNode);
                createdCount++;
                
                if (THERMAL_DEBUG) {
                    std::cout << "Created insulation nodes: " << innerNode.name << " at radius=" 
                              << innerSurfaceRadius * 1000 << "mm, " << outerNode.name << " at radius=" 
                              << outerSurfaceRadius * 1000 << "mm" << std::endl;
                }
            }
            layerIdx++;
        } else {
            // For concentric cores or Cartesian coordinates: single node
            if (isPolar) {
                // Convert polar to Cartesian
                double radius = layerX;
                double angleDeg = layerY;
                double angleRad = angleDeg * M_PI / 180.0;
                
                layerX = radius * std::cos(angleRad);
                layerY = radius * std::sin(angleRad);
                
                double radialThickness = layerWidth;
                double angularSpanDeg = layerHeight / 1000.0;
                double angularSpanRad = angularSpanDeg * M_PI / 180.0;
                double arcLength = radius * angularSpanRad;
                
                layerWidth = radialThickness;
                layerHeight = arcLength;
            }
            
            // Create single insulation layer node
            ThermalNetworkNode insulationNode;
            insulationNode.part = ThermalNodePartType::INSULATION_LAYER;
            insulationNode.name = "L_" + std::to_string(layerIdx);
            insulationNode.insulationLayerIndex = layerIdx;
            insulationNode.temperature = _config.ambientTemperature;
            insulationNode.powerDissipation = 0.0;
            insulationNode.physicalCoordinates = {layerX, layerY, 0};
            
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
    }
    
    if (THERMAL_DEBUG) {
    }
}

void Temperature::createTurnNodes() {
    auto coil = _magnetic.get_coil();
    auto turns = coil.get_turns_description();
    
    if (!turns) {
        throw std::runtime_error("Temperature::createTurnNodes: Coil has no turns description. "
                                 "Cannot build thermal model without turn geometry.");
    }
    if (turns->empty()) {
        throw std::runtime_error("Temperature::createTurnNodes: Coil has zero turns. "
                                 "Cannot build thermal model without turns to assign losses to.");
    }
    
    // Extract wire properties locally from the coil
    auto wireOpt = extractWire();
    if (!wireOpt) {
        throw std::runtime_error("Temperature::createTurnNodes: Cannot extract wire from coil. "
                                 "Wire properties are required for thermal analysis.");
    }
    
    // IMP-3: Use cached member variables
    double defaultWireWidth = _wireWidth;
    double defaultWireHeight = _wireHeight;
    double defaultWireThermalCond = _wireThermalCond;
    bool defaultIsRoundWire = _isRoundWire;
    std::optional<InsulationWireCoating> defaultWireCoating = wireOpt.value().resolve_coating();
    
    // Pre-compute turn counts per winding and per-turn index within winding
    // Uses incrementing counter per winding to handle interleaved windings correctly
    std::map<size_t, size_t> turnsPerWinding;
    std::map<size_t, size_t> windingTurnCounter;  // Incrementing counter per winding
    std::vector<size_t> turnIdxInWindingVec(turns->size(), 0);  // Per-turn index within its winding
    for (size_t i = 0; i < turns->size(); ++i) {
        size_t wIdx = coil.get_winding_index_by_name((*turns)[i].get_winding());
        turnIdxInWindingVec[i] = windingTurnCounter[wIdx];
        windingTurnCounter[wIdx]++;
        turnsPerWinding[wIdx]++;
    }
    
    // Get per-turn losses from simulation output if available
    // Per-turn losses are REQUIRED for accurate thermal analysis
    std::vector<double> turnLosses;
    bool hasPerTurnLosses = false;
    
    if (_config.windingLossesOutput) {
        auto lossesPerTurn = _config.windingLossesOutput->get_winding_losses_per_turn();
        if (lossesPerTurn && !lossesPerTurn->empty()) {
            hasPerTurnLosses = true;
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
        }
    } else if (_config.windingLosses <= 0.0) {
        // When total winding losses are zero, all turns have zero losses
        // No windingLossesOutput needed - this is an explicit zero, not a fallback
        hasPerTurnLosses = true;
        turnLosses.assign(turns->size(), 0.0);
    }
    
    if (_isToroidal) {
        
        
        // Per-turn losses are already calculated above in turnLosses vector
        // For toroidal cores, per-turn losses are strictly required
        if (!hasPerTurnLosses) {
            throw std::runtime_error("WindingLossesOutput with per-turn losses is required for thermal analysis. "
                                     "Use MagneticSimulator to calculate real losses per turn.");
        }
        
        // Validate that per-turn losses count matches turn count
        if (turnLosses.size() != turns->size()) {
            throw std::runtime_error("Per-turn losses count (" + std::to_string(turnLosses.size()) + 
                                     ") does not match turn count (" + std::to_string(turns->size()) + 
                                     "). The winding losses simulation and coil description are inconsistent.");
        }
        
        for (size_t t = 0; t < turns->size(); ++t) {
            const auto& turn = (*turns)[t];
            double turnLoss = turnLosses[t];
            
            // Calculate winding index and turn index within winding
            size_t windingIdx = coil.get_winding_index_by_name(turn.get_winding());
            size_t turnIdxInWinding = turnIdxInWindingVec[t];
            
            auto coords = turn.get_coordinates();
            if (coords.size() < 2) {
                throw std::runtime_error("Temperature::createToroidalTurnNodes: Turn " + turn.get_name() + 
                                         " has fewer than 2 coordinates. Valid coordinates are required for thermal analysis.");
            }
            double turnCenterRadius = std::sqrt(coords[0]*coords[0] + coords[1]*coords[1]);
            
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
            
            if (THERMAL_DEBUG && t < 3) {
                std::cout << "Turn " << t << " coordinates: [" << coords[0] << ", " << coords[1] 
                          << "], radius=" << turnCenterRadius*1000 << "mm" << std::endl;
            }
            innerNode.temperature = _config.ambientTemperature;
            innerNode.windingIndex = static_cast<int>(windingIdx);
            innerNode.turnIndex = static_cast<int>(turnIdxInWinding);
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
                outerNode.windingIndex = static_cast<int>(windingIdx);
                outerNode.turnIndex = static_cast<int>(turnIdxInWinding);
                outerNode.isInnerTurn = false;
                
                auto addCoords = turn.get_additional_coordinates().value()[0];
                double outerSurfaceX = addCoords[0];
                double outerSurfaceY = addCoords[1];
                double outerSurfaceRadius = std::sqrt(outerSurfaceX*outerSurfaceX + outerSurfaceY*outerSurfaceY);
                
                if (THERMAL_DEBUG && t < 3) {
                    std::cout << "Turn " << t << " additionalCoords: [" << outerSurfaceX << ", " << outerSurfaceY 
                              << "], radius=" << outerSurfaceRadius*1000 << "mm" << std::endl;
                }
                
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
        // Get column geometry for turn length calculations at each face's radial position
        auto bobbinOpt = _magnetic.get_mutable_coil().get_bobbin();
        ColumnShape colShape = ColumnShape::ROUND;  // default
        double colWidth = 0.0;
        double colDepth = 0.0;
        if (std::holds_alternative<Bobbin>(bobbinOpt)) {
            auto bobbin = std::get<Bobbin>(bobbinOpt);
            auto bobbinDesc = bobbin.get_processed_description();
            if (bobbinDesc) {
                colShape = bobbinDesc->get_column_shape();
                colWidth = bobbinDesc->get_column_width().value_or(0.0);
                colDepth = bobbinDesc->get_column_depth();
            }
        }
        
        for (size_t t = 0; t < turns->size(); ++t) {
            const auto& turn = (*turns)[t];
            
            // Calculate winding index and turn index within winding
            size_t windingIdx = coil.get_winding_index_by_name(turn.get_winding());
            size_t turnIdxInWinding = turnIdxInWindingVec[t];
            
            ThermalNetworkNode node;
            node.part = ThermalNodePartType::TURN;
            node.name = turn.get_name();
            node.temperature = _config.ambientTemperature;
            // Per-turn losses are REQUIRED - no fallback allowed
            if (!hasPerTurnLosses || t >= turnLosses.size()) {
                throw std::runtime_error("Per-turn winding losses are required for thermal analysis. "
                                         "Use MagneticSimulator to calculate real losses per turn.");
            }
            node.powerDissipation = turnLosses[t];
            node.windingIndex = static_cast<int>(windingIdx);
            node.turnIndex = static_cast<int>(turnIdxInWinding);
            
            auto coords = turn.get_coordinates();
            if (coords.size() >= 3) {
                node.physicalCoordinates = {coords[0], coords[1], coords[2]};
            } else if (coords.size() >= 2) {
                node.physicalCoordinates = {coords[0], coords[1], 0};
            } else {
                throw std::runtime_error("Temperature::createConcentricTurnNodes: Turn " + turn.get_name() + 
                                         " has fewer than 2 coordinates. Valid coordinates are required for thermal analysis.");
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
            
            // Use turn's x-coordinate (radial position) as centerRadius
            double centerRadius = coords[0];
            node.initializeConcentricTurnQuadrants(wireWidth, wireHeight, centerRadius,
                                                   colShape, colWidth, colDepth,
                                                   defaultWireThermalCond,
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
    createInsulationLayerConnections();
    createTurnToSolidConnections();
    createConvectionConnections();
    
    // NEW: Apply MAS cooling configuration if specified
    if (_config.masCooling.has_value()) {
        applyMasCooling(_config.masCooling.value());
    }
    
    if (THERMAL_DEBUG) {
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
    double coreK = getCoreThermalConductivity();
    
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
    
    double bobbinK = getBobbinThermalConductivity();
    
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
                                               double dist, double area, double k) {
        ThermalResistanceElement r;
        r.nodeFromId = fromIdx;
        r.quadrantFrom = fromFace;
        r.nodeToId = toIdx;
        r.quadrantTo = toFace;
        r.type = HeatTransferType::CONDUCTION;
        
        r.resistance = ThermalResistance::calculateConductionResistance(dist, k, area);
        _resistances.push_back(r);
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
                                          dist, area, bobbinK);
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
                                      minDist, area, bobbinK);
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
                                      minDist, area, bobbinK);
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
                                          dist, area, bobbinK);
        }
        if (bobbinBottomYokeIdx != std::numeric_limits<size_t>::max()) {
            double dx = _nodes[bobbinColumnIdx].physicalCoordinates[0] - _nodes[bobbinBottomYokeIdx].physicalCoordinates[0];
            double dy = _nodes[bobbinColumnIdx].physicalCoordinates[1] - _nodes[bobbinBottomYokeIdx].physicalCoordinates[1];
            double dist = std::sqrt(dx*dx + dy*dy);
            double area = _nodes[bobbinColumnIdx].dimensions.depth * _nodes[bobbinColumnIdx].dimensions.width * 0.5;
            // Column BOTTOM ↔ Bottom yoke LEFT
            createConnectionWithQuadrants(bobbinColumnIdx, ThermalNodeFace::TANGENTIAL_RIGHT, // Column BOTTOM
                                          bobbinBottomYokeIdx, ThermalNodeFace::RADIAL_INNER,  // Yoke LEFT
                                          dist, area, bobbinK);
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
    // NOTE: O(N^2) complexity over turn nodes. For designs with >1000 turns,
    // consider implementing spatial indexing (grid hash or k-d tree) for neighbor search.
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
                if (!_nodes[node1Idx].turnIndex.has_value()) {
                    throw std::runtime_error("Temperature::createConcentricTurnToTurnConnections: Turn node " +
                                             _nodes[node1Idx].name + " is missing turnIndex.");
                }
                if (!_nodes[node2Idx].turnIndex.has_value()) {
                    throw std::runtime_error("Temperature::createConcentricTurnToTurnConnections: Turn node " +
                                             _nodes[node2Idx].name + " is missing turnIndex.");
                }
                int turn1Idx = _nodes[node1Idx].turnIndex.value();
                int turn2Idx = _nodes[node2Idx].turnIndex.value();

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
                    } catch (const std::exception& e) { // IMP-8
                        // If we can't determine insulation layers (e.g., missing layer description),
                        // assume no solid insulation and create direct connection
                        hasSolidInsulationBetween = false;
                        if (THERMAL_DEBUG) {
                            std::cout << "Temperature: insulation-layer query failed for concentric turn pair, assuming none: "
                                      << e.what() << std::endl;
                        }
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
    // NOTE: O(N^2) complexity over turn nodes. For designs with >1000 turns,
    // consider implementing spatial indexing (grid hash or k-d tree) for neighbor search.
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

            // Guard against coincident turn centers (degenerate geometry).
            // Without this, the unit-vector divisions below produce NaN.
            if (centerDistance < 1e-12) continue;

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
                if (!_nodes[node1Idx].turnIndex.has_value()) {
                    throw std::runtime_error("Temperature::createToroidalTurnToTurnConnections: Turn node " +
                                             _nodes[node1Idx].name + " is missing turnIndex.");
                }
                if (!_nodes[node2Idx].turnIndex.has_value()) {
                    throw std::runtime_error("Temperature::createToroidalTurnToTurnConnections: Turn node " +
                                             _nodes[node2Idx].name + " is missing turnIndex.");
                }
                int turn1Idx = _nodes[node1Idx].turnIndex.value();
                int turn2Idx = _nodes[node2Idx].turnIndex.value();

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
                        // Zero-thickness layers are geometric boundary markers, not physical insulation.
                        // Filter them out so they don't block direct turn-to-turn connections.
                        layersBetween.erase(
                            std::remove_if(layersBetween.begin(), layersBetween.end(),
                                [](const auto& layer) {
                                    return layer.get_dimensions().empty() || layer.get_dimensions()[0] < 1e-9;
                                }),
                            layersBetween.end()
                        );
                        hasSolidInsulationBetween = !layersBetween.empty();
                    } catch (const std::exception& e) { // IMP-8
                        // If we can't determine insulation layers (e.g., missing layer description),
                        // assume no solid insulation and create direct connection
                        hasSolidInsulationBetween = false;
                        if (THERMAL_DEBUG) {
                            std::cout << "Temperature: insulation-layer query failed for toroidal turn pair, assuming none: "
                                      << e.what() << std::endl;
                        }
                    }
                }

                // Also check geometrically if an insulation layer node is between these turns
                // For toroidal cores: check if insulation layer is radially between the two turns
                bool insulationNodeBetween = false;
                for (size_t k = 0; k < _nodes.size(); k++) {
                    if (_nodes[k].part != ThermalNodePartType::INSULATION_LAYER) continue;

                    const auto& insulationNode = _nodes[k];
                    double insX = insulationNode.physicalCoordinates[0];
                    double insY = insulationNode.physicalCoordinates[1];

                    // Calculate radial distance of insulation layer from core center
                    double insRadius = std::sqrt(insX * insX + insY * insY);

                    // Calculate radial distances of the two turns from core center
                    double turn1X = node1.physicalCoordinates[0];
                    double turn1Y = node1.physicalCoordinates[1];
                    double turn2X = node2.physicalCoordinates[0];
                    double turn2Y = node2.physicalCoordinates[1];
                    double turn1Radius = std::sqrt(turn1X * turn1X + turn1Y * turn1Y);
                    double turn2Radius = std::sqrt(turn2X * turn2X + turn2Y * turn2Y);

                    // Calculate angular positions
                    double turn1Angle = std::atan2(turn1Y, turn1X);
                    double turn2Angle = std::atan2(turn2Y, turn2X);
                    double insAngle = std::atan2(insY, insX);

                    // Normalize angle differences to [-pi, pi]
                    auto normalizeAngleDiff = [](double diff) {
                        while (diff > std::numbers::pi) diff -= 2.0 * std::numbers::pi;
                        while (diff < -std::numbers::pi) diff += 2.0 * std::numbers::pi;
                        return diff;
                    };

                    double angleDiff1 = normalizeAngleDiff(insAngle - turn1Angle);
                    double angleDiff2 = normalizeAngleDiff(insAngle - turn2Angle);

                    // Check if insulation is angularly near both turns (within reasonable tolerance)
                    double angularTolerance = 0.5;  // ~30 degrees in radians
                    bool angularlyNearBoth = (std::abs(angleDiff1) < angularTolerance &&
                                             std::abs(angleDiff2) < angularTolerance);

                    // Check if insulation is radially between the two turns
                    double minTurnRadius = std::min(turn1Radius, turn2Radius);
                    double maxTurnRadius = std::max(turn1Radius, turn2Radius);
                    double radiusTolerance = std::max(node1.dimensions.width, node1.dimensions.height);

                    bool radiallyBetween = (insRadius > minTurnRadius - radiusTolerance) &&
                                          (insRadius < maxTurnRadius + radiusTolerance);

                    if (angularlyNearBoth && radiallyBetween) {
                        insulationNodeBetween = true;
                        break;
                    }
                }

                // Only create direct turn-to-turn connection if there's NO solid insulation layer
                // Solid insulation layers are now modeled as separate thermal nodes
                if (!hasSolidInsulationBetween && !insulationNodeBetween) {
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
        
        // Calculate resistance through air/bobbin gap - use quadrant surface area
        auto* turnQuadrant = _nodes[turnIdx].getQuadrant(r.quadrantTo);
        double contactArea = turnQuadrant ? turnQuadrant->surfaceArea : (turnWidth * turnHeight * 0.5);
        // Issue 15: Use configurable interface conductivity for turn-to-bobbin contact
        double k_interface = _config.turnToBobbinInterfaceConductivity;
        r.resistance = dist / (k_interface * contactArea);
        
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
                auto* tq = turnNode.getQuadrant(ThermalNodeFace::TANGENTIAL_LEFT);
                double turnContactArea = tq ? tq->surfaceArea : (turnNode.dimensions.width * turnNode.dimensions.height * 0.25);
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
                auto* tq = turnNode.getQuadrant(ThermalNodeFace::TANGENTIAL_RIGHT);
                double turnContactArea = tq ? tq->surfaceArea : (turnNode.dimensions.width * turnNode.dimensions.height * 0.25);
                totalContactArea += turnContactArea;
            }
            
            for (const auto& [turnIdx, dist] : nearbyTurns) {
                createBobbinToTurnConnection(bobbinBottomYokeIdx, ThermalNodeFace::TANGENTIAL_LEFT, 
                                            turnIdx, dist);
            }
        }
    }
}

std::set<std::pair<size_t, ThermalNodeFace>> Temperature::getConnectedQuadrants() const {
    std::set<std::pair<size_t, ThermalNodeFace>> connectedQuadrants;
    
    for (const auto& res : _resistances) {
        if (res.type == HeatTransferType::CONDUCTION) {
            // Add both from and to quadrants
            connectedQuadrants.insert({res.nodeFromId, res.quadrantFrom});
            connectedQuadrants.insert({res.nodeToId, res.quadrantTo});
        }
    }
    
    return connectedQuadrants;
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
    
    // Get already-connected quadrants to prevent multiple connections on same face
    auto connectedQuadrants = getConnectedQuadrants();
    
    // For each turn, create conduction connections for ALL quadrants touching any bobbin
    for (size_t turnIdx : turnNodeIndices) {
        auto& turnNode = _nodes[turnIdx];
        
        double turnX = turnNode.physicalCoordinates[0];
        double turnY = turnNode.physicalCoordinates[1];
        
        // Track which (face, bobbin) pairs we've already connected to avoid duplicates
        std::set<std::pair<ThermalNodeFace, size_t>> connectedPairs;
        
        for (const auto& quadrant : turnNode.quadrants) {
            if (quadrant.face == ThermalNodeFace::NONE) continue;
            
            // Skip if this quadrant is already connected to something else
            if (connectedQuadrants.count({turnIdx, quadrant.face})) continue;
            
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
                    double contactArea = q ? q->surfaceArea : (turnWidth * turnHeight * 0.5);
                    double distance = std::max(dist, 1e-6);
                    
                    // IMP-NEW-05: Configurable interface conductivity (was k_air=0.025)
                    double k_interface = _config.turnToBobbinInterfaceConductivity;
                    r.resistance = distance / (k_interface * contactArea);
                    
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
    

    // Get already-connected quadrants to prevent multiple connections on same face
    auto connectedQuadrants = getConnectedQuadrants();

    // Skip toroidal handling here - connections will be created after processing all insulation nodes
    // This ensures each turn connects only to the closest angular chunk of each layer
    if (!_isToroidal) {
        // For concentric cores: iterate through turns and find closest insulation layer on each side
        // This ensures turns connect to the correct insulation layer based on perpendicular distance
        
        for (size_t turnIdx : turnNodeIndices) {
            const auto& turnNode = _nodes[turnIdx];
            double turnX = turnNode.physicalCoordinates[0];
            double turnY = turnNode.physicalCoordinates[1];
            double turnWidth = turnNode.dimensions.width;
            double turnHeight = turnNode.dimensions.height;
            bool turnIsRound = (turnNode.crossSectionalShape == TurnCrossSectionalShape::ROUND);
            double turnRadius = turnIsRound ? (turnWidth / 2.0) : 0.0;

            // Connection threshold: For insulation connections, use generous distance
            // Must be larger than convection blocking distance to ensure blocked faces have conduction path
            double connectionThreshold = 10.0 * std::max(turnWidth, turnHeight);  // 10x turn size

            // Find closest insulation layer on LEFT and RIGHT sides of the turn
            struct InsulationCandidate {
                size_t insulationIdx = 0;
                double distance = 0.0;  // Surface-to-surface distance
                double insulationK = 0.2;
                double insulationWidth = 0.0;
                double insulationHeight = 0.0;
            };

            InsulationCandidate leftCandidate;
            InsulationCandidate rightCandidate;
            bool hasLeft = false;
            bool hasRight = false;

            for (size_t insulationIdx : insulationNodeIndices) {
                const auto& insulationNode = _nodes[insulationIdx];
                double insulationX = insulationNode.physicalCoordinates[0];
                double insulationY = insulationNode.physicalCoordinates[1];
                double insulationW = insulationNode.dimensions.width;
                double insulationH = insulationNode.dimensions.height;

                // Get insulation thermal conductivity
                double insulationK = kInsulation_DefaultConductivity;
                if (insulationNode.quadrants[0].thermalConductivity > 0) {
                    insulationK = insulationNode.quadrants[0].thermalConductivity;
                }

                // Calculate insulation layer edges (insulation is always rectangular)
                double insulationLeftEdge = insulationX - insulationW / 2.0;
                double insulationRightEdge = insulationX + insulationW / 2.0;

                // Check if turn is within the insulation layer's Y span (with tolerance)
                double insulationBottom = insulationY - insulationH / 2.0;
                double insulationTop = insulationY + insulationH / 2.0;
                double turnBottom = turnY - turnHeight / 2.0;
                double turnTop = turnY + turnHeight / 2.0;

                // Use generous Y-overlap check: allow partial overlap or nearby turns
                double yTolerance = std::max(turnHeight, insulationH);
                bool yOverlaps = !(turnTop + yTolerance < insulationBottom || turnBottom - yTolerance > insulationTop);
                if (!yOverlaps) continue;

                // Calculate perpendicular distance from turn to insulation surface
                // For round turns: distance from center to nearest surface of rectangular insulation
                // For rectangular turns: edge-to-edge distance

                double distToLeft, distToRight;

                if (turnIsRound) {
                    // Round wire: calculate perpendicular distance from center to insulation surface
                    // Distance to LEFT insulation (insulation's right surface)
                    distToLeft = turnX - turnRadius - insulationRightEdge;
                    // Distance to RIGHT insulation (insulation's left surface)
                    distToRight = insulationLeftEdge - (turnX + turnRadius);
                } else {
                    // Rectangular wire: edge-to-edge distance
                    double turnLeftEdge = turnX - turnWidth / 2.0;
                    double turnRightEdge = turnX + turnWidth / 2.0;
                    distToLeft = turnLeftEdge - insulationRightEdge;
                    distToRight = insulationLeftEdge - turnRightEdge;
                }

                // Check if insulation is to the LEFT of the turn
                if (distToLeft >= -connectionThreshold && distToLeft <= connectionThreshold) {
                    if (!hasLeft || std::abs(distToLeft) < std::abs(leftCandidate.distance)) {
                        leftCandidate = {insulationIdx, distToLeft, insulationK, insulationW, insulationH};
                        hasLeft = true;
                    }
                }

                // Check if insulation is to the RIGHT of the turn
                if (distToRight >= -connectionThreshold && distToRight <= connectionThreshold) {
                    if (!hasRight || std::abs(distToRight) < std::abs(rightCandidate.distance)) {
                        rightCandidate = {insulationIdx, distToRight, insulationK, insulationW, insulationH};
                        hasRight = true;
                    }
                }
            }
            
            // Create connection to LEFT insulation layer (connects to turn's LEFT face)
            // Skip if turn's LEFT face is already connected
            if (hasLeft && !connectedQuadrants.count({turnIdx, ThermalNodeFace::RADIAL_INNER})) {
                const auto& insNode = _nodes[leftCandidate.insulationIdx];
                
                ThermalResistanceElement r;
                r.nodeFromId = leftCandidate.insulationIdx;
                r.quadrantFrom = ThermalNodeFace::RADIAL_OUTER;  // Insulation's RIGHT face touches turn's LEFT
                r.nodeToId = turnIdx;
                r.quadrantTo = ThermalNodeFace::RADIAL_INNER;    // Turn's LEFT face
                r.type = HeatTransferType::CONDUCTION;
                
                auto* leftQ = turnNode.getQuadrant(ThermalNodeFace::RADIAL_INNER);
                double contactArea = leftQ ? leftQ->surfaceArea : (turnWidth * turnHeight * 0.5);
                double conductionDistance = leftCandidate.insulationWidth / 2.0;
                
                r.resistance = ThermalResistance::calculateConductionResistance(
                    conductionDistance, leftCandidate.insulationK, contactArea);
                
                _resistances.push_back(r);
                
                if (THERMAL_DEBUG) {
                    std::cout << "Turn " << turnNode.name << " (winding=" << turnNode.windingIndex.value_or(-1) 
                              << ") LEFT -> Insulation " << insNode.name 
                              << " RIGHT (dist=" << leftCandidate.distance * 1000 << "mm, R=" 
                              << r.resistance << " K/W)" << std::endl;
                }
            }
            
            // Create connection to RIGHT insulation layer (connects to turn's RIGHT face)
            // Skip if turn's RIGHT face is already connected
            if (hasRight && !connectedQuadrants.count({turnIdx, ThermalNodeFace::RADIAL_OUTER})) {
                const auto& insNode = _nodes[rightCandidate.insulationIdx];
                
                ThermalResistanceElement r;
                r.nodeFromId = rightCandidate.insulationIdx;
                r.quadrantFrom = ThermalNodeFace::RADIAL_INNER;  // Insulation's LEFT face touches turn's RIGHT
                r.nodeToId = turnIdx;
                r.quadrantTo = ThermalNodeFace::RADIAL_OUTER;    // Turn's RIGHT face
                r.type = HeatTransferType::CONDUCTION;
                
                auto* rightQ = turnNode.getQuadrant(ThermalNodeFace::RADIAL_OUTER);
                double contactArea = rightQ ? rightQ->surfaceArea : (turnWidth * turnHeight * 0.5);
                double conductionDistance = rightCandidate.insulationWidth / 2.0;
                
                r.resistance = ThermalResistance::calculateConductionResistance(
                    conductionDistance, rightCandidate.insulationK, contactArea);
                
                _resistances.push_back(r);
                
                if (THERMAL_DEBUG) {
                    std::cout << "Turn " << turnNode.name << " (winding=" << turnNode.windingIndex.value_or(-1) 
                              << ") RIGHT -> Insulation " << insNode.name 
                              << " LEFT (dist=" << rightCandidate.distance * 1000 << "mm, R=" 
                              << r.resistance << " K/W)" << std::endl;
                }
            }
        }
    }
    
    // ============================================================================
    // Toroidal Turn-to-Insulation Connections (Radial Proximity + Angular Matching)
    // ============================================================================
    // For each turn node:
    // 1. Check both radial sides (inner and outer surface) against all insulation layers
    // 2. If turn surface is close to layer surface (within 15% of wire radius), mark for connection
    // 3. Find the closest insulation node by angle for each matching layer
    // 4. Choose the quadrant combination (RI-RO or RO-RI) with closest limit coordinates
    size_t connectionCount = 0;
    if (_isToroidal && !insulationNodeIndices.empty()) {
        // Get core dimensions for inner/outer region detection
        auto core = _magnetic.get_core();
        auto coreDims = flatten_dimensions(core.resolve_shape().get_dimensions().value());
        
        // Group insulation nodes by layer index
        std::map<int, std::vector<size_t>> insNodesByLayer;
        for (size_t insIdx : insulationNodeIndices) {
            if (!_nodes[insIdx].insulationLayerIndex.has_value()) {
                throw std::runtime_error("Temperature::createTurnToInsulationConnections: Insulation layer node " +
                                         _nodes[insIdx].name + " is missing insulationLayerIndex.");
            }
            int layerIdx = _nodes[insIdx].insulationLayerIndex.value();
            insNodesByLayer[layerIdx].push_back(insIdx);
        }
        
        // Pre-calculate layer data for efficiency
        // Each layer has both INNER and OUTER nodes - handle them separately
        struct LayerData {
            int layerIdx;
            double layerRadius;      // Center radius from core center
            double semiThickness;    // Half of radial thickness
            bool isInner;            // true = inner side (inside core hole), false = outer side (outside core)
        };
        std::vector<LayerData> allLayers;
        for (const auto& [layerIdx, layerInsIndices] : insNodesByLayer) {
            if (layerInsIndices.empty()) continue;
            
            // Find one inner node and one outer node for this layer
            const ThermalNetworkNode* innerNode = nullptr;
            const ThermalNetworkNode* outerNode = nullptr;
            
            for (size_t insIdx : layerInsIndices) {
                const auto& node = _nodes[insIdx];
                if (node.isInnerTurn && !innerNode) {
                    innerNode = &node;
                } else if (!node.isInnerTurn && !outerNode) {
                    outerNode = &node;
                }
                if (innerNode && outerNode) break;
            }
            
            // Add inner side layer data
            if (innerNode) {
                double insX = innerNode->physicalCoordinates[0];
                double insY = innerNode->physicalCoordinates[1];
                double layerRadius = std::sqrt(insX*insX + insY*insY);
                double semiThickness = innerNode->dimensions.width / 2.0;
                allLayers.push_back({layerIdx, layerRadius, semiThickness, true});
            }
            
            // Add outer side layer data
            if (outerNode) {
                double insX = outerNode->physicalCoordinates[0];
                double insY = outerNode->physicalCoordinates[1];
                double layerRadius = std::sqrt(insX*insX + insY*insY);
                double semiThickness = outerNode->dimensions.width / 2.0;
                allLayers.push_back({layerIdx, layerRadius, semiThickness, false});
            }
        }
        
        for (size_t turnIdx : turnNodeIndices) {
            const auto& turnNode = _nodes[turnIdx];
            double turnX = turnNode.physicalCoordinates[0];
            double turnY = turnNode.physicalCoordinates[1];
            double turnRadiusPos = std::sqrt(turnX*turnX + turnY*turnY);
            double turnAngle = std::atan2(turnY, turnX);
            
            // Get turn physical dimensions
            double turnDiameter = std::max(turnNode.dimensions.width, turnNode.dimensions.height);
            double turnRadiusActual = turnDiameter / 2.0;  // Physical radius of the turn wire
            double proximityThreshold = turnRadiusActual * 0.15;  // 15% of wire radius
            
            // Check inner/outer consistency
            bool isTurnInner = turnNode.isInnerTurn;
            
            // Calculate turn's inner and outer surface radii
            double turnSurfaceInnerRadius = turnRadiusPos - turnRadiusActual;
            double turnSurfaceOuterRadius = turnRadiusPos + turnRadiusActual;
            
            // Find layers that are close on either radial side of the turn
            std::vector<std::pair<int, double>> matchingLayers;  // (layerIdx, whichSide: -1=inner, +1=outer)
            
            for (const auto& layer : allLayers) {
                // Skip layers on opposite side of core (inner vs outer)
                if (layer.isInner != isTurnInner) continue;
                
                // Calculate layer's inner and outer surface radii
                double layerSurfaceInnerRadius = layer.layerRadius - layer.semiThickness;
                double layerSurfaceOuterRadius = layer.layerRadius + layer.semiThickness;
                
                // Check if turn's INNER surface is close to layer's OUTER surface
                // (layer is closer to center than turn)
                double gapInner = turnSurfaceInnerRadius - layerSurfaceOuterRadius;
                if (std::abs(gapInner) < proximityThreshold) {
                    matchingLayers.push_back({layer.layerIdx, -1.0});  // -1 = turn's inner side
                }
                
                // Check if turn's OUTER surface is close to layer's INNER surface
                // (layer is farther from center than turn)
                double gapOuter = layerSurfaceInnerRadius - turnSurfaceOuterRadius;
                if (std::abs(gapOuter) < proximityThreshold) {
                    matchingLayers.push_back({layer.layerIdx, +1.0});  // +1 = turn's outer side
                }
            }
            
            // For each matching layer, find closest node by angle and create connection
            for (const auto& [layerIdx, whichSide] : matchingLayers) {
                const auto& layerInsIndices = insNodesByLayer[layerIdx];
                
                // Find closest insulation node in this layer by angular proximity
                // Must match inner/outer side
                size_t bestInsIdx = 0;
                double bestAngleDiff = std::numeric_limits<double>::max();
                bool foundNode = false;
                
                for (size_t insIdx : layerInsIndices) {
                    const auto& insNode = _nodes[insIdx];
                    
                    // Skip nodes on opposite side (inner vs outer)
                    if (insNode.isInnerTurn != isTurnInner) continue;
                    
                    double insX = insNode.physicalCoordinates[0];
                    double insY = insNode.physicalCoordinates[1];
                    double insAngle = std::atan2(insY, insX);
                    
                    double angleDiff = std::abs(insAngle - turnAngle);
                    if (angleDiff > M_PI) angleDiff = 2.0 * M_PI - angleDiff;
                    
                    if (angleDiff < bestAngleDiff) {
                        bestAngleDiff = angleDiff;
                        bestInsIdx = insIdx;
                        foundNode = true;
                    }
                }
                
                if (!foundNode) continue;
                
                const auto& insNode = _nodes[bestInsIdx];
                
                // Get quadrants for distance checking
                auto* turnQInner = turnNode.getQuadrant(ThermalNodeFace::RADIAL_INNER);
                auto* turnQOuter = turnNode.getQuadrant(ThermalNodeFace::RADIAL_OUTER);
                auto* insQInner = insNode.getQuadrant(ThermalNodeFace::RADIAL_INNER);
                auto* insQOuter = insNode.getQuadrant(ThermalNodeFace::RADIAL_OUTER);
                
                // Determine which quadrant combination to use based on whichSide and physical proximity
                ThermalNodeFace bestTurnFace = ThermalNodeFace::NONE;
                ThermalNodeFace bestInsFace = ThermalNodeFace::NONE;
                double minQuadrantDistance = 1e9;
                
                if (whichSide < 0) {
                    // Turn's INNER surface faces the layer
                    // Valid combination: Turn RI <-> Insulation RO
                    if (turnQInner && insQOuter) {
                        double dx = turnQInner->limitCoordinates[0] - insQOuter->limitCoordinates[0];
                        double dy = turnQInner->limitCoordinates[1] - insQOuter->limitCoordinates[1];
                        minQuadrantDistance = std::sqrt(dx*dx + dy*dy);
                        bestTurnFace = ThermalNodeFace::RADIAL_INNER;
                        bestInsFace = ThermalNodeFace::RADIAL_OUTER;
                    }
                } else {
                    // Turn's OUTER surface faces the layer
                    // Valid combination: Turn RO <-> Insulation RI
                    if (turnQOuter && insQInner) {
                        double dx = turnQOuter->limitCoordinates[0] - insQInner->limitCoordinates[0];
                        double dy = turnQOuter->limitCoordinates[1] - insQInner->limitCoordinates[1];
                        minQuadrantDistance = std::sqrt(dx*dx + dy*dy);
                        bestTurnFace = ThermalNodeFace::RADIAL_OUTER;
                        bestInsFace = ThermalNodeFace::RADIAL_INNER;
                    }
                }
                
                // Only create connection if we found valid facing quadrants
                if (bestTurnFace == ThermalNodeFace::NONE || bestInsFace == ThermalNodeFace::NONE) continue;
                
                // Create the thermal resistance connection
                ThermalResistanceElement r;
                r.nodeFromId = bestInsIdx;
                r.quadrantFrom = bestInsFace;
                r.nodeToId = turnIdx;
                r.quadrantTo = bestTurnFace;
                r.type = HeatTransferType::CONDUCTION;
                
                // Calculate contact area (minimum of the two quadrant surface areas)
                auto* turnQ = turnNode.getQuadrant(bestTurnFace);
                auto* insQ = insNode.getQuadrant(bestInsFace);
                double contactArea = std::min(
                    turnQ ? turnQ->surfaceArea : 0.0,
                    insQ ? insQ->surfaceArea : 0.0
                );
                if (contactArea <= 0) {
                    // Fallback: use quadrant surface area from the turn
                    auto* fallbackQ = turnNode.getQuadrant(bestTurnFace);
                    contactArea = fallbackQ ? fallbackQ->surfaceArea : 
                                  (turnNode.dimensions.width * turnNode.dimensions.height * 0.5);
                }
                
                // Conduction distance is the surface-to-surface distance
                double conductionDistance = std::max(minQuadrantDistance, 1e-6);
                
                // Get insulation thermal conductivity
                double insulationK = kInsulation_DefaultConductivity;
                if (insNode.quadrants[0].thermalConductivity > 0) {
                    insulationK = insNode.quadrants[0].thermalConductivity;
                }
                
                r.resistance = ThermalResistance::calculateConductionResistance(
                    conductionDistance, insulationK, contactArea);
                
                _resistances.push_back(r);
                connectionCount++;
                
                if (THERMAL_DEBUG) {
                    std::cout << "Toroidal: " << insNode.name << "[" << magic_enum::enum_name(bestInsFace) 
                              << "] -> " << turnNode.name << "[" << magic_enum::enum_name(bestTurnFace)
                              << "] side=" << (whichSide < 0 ? "inner" : "outer") 
                              << ", angleDiff=" << bestAngleDiff * 180.0 / M_PI << "deg, R=" 
                              << r.resistance << " K/W" << std::endl;
                }
            }
        }
        
        if (THERMAL_DEBUG) {
        }
    }
}

void Temperature::createInsulationLayerConnections() {
    if (!_isToroidal) return;
    
    // Group insulation layer nodes by (layerIdx, isInner) and segment index
    // Node naming convention: IL_{layerIdx}_{segIdx}_{i|o}
    // We need to connect adjacent segments tangentially and inner↔outer radially
    
    // Structure: layerGroups[layerIdx][isInner] = vector of (segIdx, nodeIdx) sorted by segIdx
    struct SegNode {
        size_t segIdx;
        size_t nodeIdx;
    };
    // Map from layerIdx -> inner nodes and outer nodes
    std::map<size_t, std::vector<SegNode>> innerNodes;  // layerIdx -> sorted by segIdx
    std::map<size_t, std::vector<SegNode>> outerNodes;
    
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part != ThermalNodePartType::INSULATION_LAYER) continue;
        
        const std::string& name = _nodes[i].name;
        // Parse IL_{layerIdx}_{segIdx}_{i|o}
        if (name.substr(0, 3) != "IL_") continue;  // Only toroidal segmented insulation nodes
        
        // Extract layerIdx and segIdx from name
        size_t firstUnderscore = name.find('_', 0);   // After "IL"
        size_t secondUnderscore = name.find('_', firstUnderscore + 1);
        size_t thirdUnderscore = name.find('_', secondUnderscore + 1);
        
        if (secondUnderscore == std::string::npos || thirdUnderscore == std::string::npos) continue;
        
        size_t layerIdx = std::stoul(name.substr(firstUnderscore + 1, secondUnderscore - firstUnderscore - 1));
        size_t segIdx = std::stoul(name.substr(secondUnderscore + 1, thirdUnderscore - secondUnderscore - 1));
        std::string suffix = name.substr(thirdUnderscore + 1);
        
        if (suffix == "i") {
            innerNodes[layerIdx].push_back({segIdx, i});
        } else if (suffix == "o") {
            outerNodes[layerIdx].push_back({segIdx, i});
        }
    }
    
    // Sort each group by segment index
    for (auto& [layerIdx, nodes] : innerNodes) {
        std::sort(nodes.begin(), nodes.end(), [](const SegNode& a, const SegNode& b) {
            return a.segIdx < b.segIdx;
        });
    }
    for (auto& [layerIdx, nodes] : outerNodes) {
        std::sort(nodes.begin(), nodes.end(), [](const SegNode& a, const SegNode& b) {
            return a.segIdx < b.segIdx;
        });
    }
    
    // Helper: create tangential connections in circular fashion for a set of segment nodes
    auto connectTangentially = [&](const std::vector<SegNode>& segNodes) {
        if (segNodes.size() < 2) return;
        
        for (size_t i = 0; i < segNodes.size(); i++) {
            size_t nextI = (i + 1) % segNodes.size();
            
            size_t fromIdx = segNodes[i].nodeIdx;
            size_t toIdx = segNodes[nextI].nodeIdx;
            
            ThermalResistanceElement r;
            r.nodeFromId = fromIdx;
            r.quadrantFrom = ThermalNodeFace::TANGENTIAL_RIGHT;
            r.nodeToId = toIdx;
            r.quadrantTo = ThermalNodeFace::TANGENTIAL_LEFT;
            r.type = HeatTransferType::CONDUCTION;
            
            auto* q1 = _nodes[fromIdx].getQuadrant(ThermalNodeFace::TANGENTIAL_RIGHT);
            auto* q2 = _nodes[toIdx].getQuadrant(ThermalNodeFace::TANGENTIAL_LEFT);
            
            if (q1 && q2) {
                double contactArea = std::min(q1->surfaceArea, q2->surfaceArea);
                r.resistance = q1->calculateConductionResistance(*q2, contactArea);
                _resistances.push_back(r);
            }
        }
    };
    
    // Connect tangentially for each layer's inner and outer node rings
    for (auto& [layerIdx, nodes] : innerNodes) {
        connectTangentially(nodes);
    }
    for (auto& [layerIdx, nodes] : outerNodes) {
        connectTangentially(nodes);
    }
    
    if (THERMAL_DEBUG) {
        std::cout << "Created insulation layer tangential connections" << std::endl;
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

                    if (!turnNode.turnIndex.has_value()) {
                        throw std::runtime_error("Temperature::createTurnToCoreConnections: Turn node " +
                                                 turnNode.name + " is missing turnIndex.");
                    }
                    int turnIdx = turnNode.turnIndex.value();

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

                    if (!turnNode.turnIndex.has_value()) {
                        throw std::runtime_error("Temperature::createTurnToCoreConnections: Turn node " +
                                                 turnNode.name + " is missing turnIndex.");
                    }
                    int turnIdx = turnNode.turnIndex.value();
                    
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
    // IMP-3: Use cached member variables
    double wireWidth = _wireWidth;
    
    
    // Get convection coefficient
    double surfaceTemp = _config.ambientTemperature + kConvection_InitialDeltaT;
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
    
    // IMP-4: Dispatch to geometry-specific method
    if (_isToroidal) {
        createToroidalConvectionConnections(ambientIdx, h_conv);
    } else if (_isPlanar) {
        createPlanarConvectionConnections(ambientIdx, h_conv);
    } else {
        createConcentricConvectionConnections(ambientIdx, h_conv);
    }
}

// IMP-4: Toroidal convection connections
void Temperature::createToroidalConvectionConnections(size_t ambientIdx, double h_conv) {
    double maxConvectionDist = getMaximumDistanceForConvection();

    auto core = _magnetic.get_core();
    auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    double coreInnerR = dimensions["B"] / 2.0;
    double coreOuterR = dimensions["A"] / 2.0;

    std::set<std::string> connectedQuadrants;
    for (const auto& res : _resistances) {
        if (res.type == HeatTransferType::CONDUCTION) {
            connectedQuadrants.insert(std::to_string(res.nodeFromId) + "_" + std::string(magic_enum::enum_name(res.quadrantFrom)));
            connectedQuadrants.insert(std::to_string(res.nodeToId) + "_" + std::string(magic_enum::enum_name(res.quadrantTo)));
        }
    }

    bool hasInsulationLayers = false;
    for (size_t i = 0; i < _nodes.size(); i++) {
        if (_nodes[i].part == ThermalNodePartType::INSULATION_LAYER) { hasInsulationLayers = true; break; }
    }

    // For each turn node, check each quadrant for exposure to air
        for (size_t i = 0; i < _nodes.size(); i++) {
            if (_nodes[i].part != ThermalNodePartType::TURN) continue;
            
            // If insulation layers exist, turns should NOT have convection to ambient
            // Only insulation layer nodes should have convection (on radial faces)
            if (hasInsulationLayers) continue;
            
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
                            // Proper angle normalization to [-pi, pi]
                            double angleDiff = nodeAngle - otherAngle;
                            while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                            while (angleDiff < -M_PI) angleDiff += 2 * M_PI;
                            
                            // Block if there's an object in roughly same angular direction
                            if (std::abs(angleDiff) < kConvection_AngularBlockingTolerance) {
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
                            // Proper angle normalization to [-pi, pi]
                            double angleDiff = nodeAngle - otherAngle;
                            while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                            while (angleDiff < -M_PI) angleDiff += 2 * M_PI;
                            
                            // Block if there's an object in roughly same angular direction
                            if (std::abs(angleDiff) < kConvection_AngularBlockingTolerance) {
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
                    
                    // Also check if tangential face is covered by an insulation layer
                    // This happens when the turn is at the same radial position as an insulation layer
                    if (isExposed) {
                        for (size_t j = 0; j < _nodes.size(); j++) {
                            if (_nodes[j].part != ThermalNodePartType::INSULATION_LAYER) continue;
                            
                            double insX = _nodes[j].physicalCoordinates[0];
                            double insY = _nodes[j].physicalCoordinates[1];
                            double insR = std::sqrt(insX*insX + insY*insY);
                            double insAngle = std::atan2(insY, insX);
                            
                            // Check if insulation is at similar radial position (within wire width)
                            // and in the same inner/outer region
                            bool isInsInner = (insR < coreInnerR);
                            if (node.isInnerTurn != isInsInner) continue;
                            
                            // Check if insulation is at similar radial position
                            // Use full wire width as threshold to account for wire radius + insulation thickness
                            double radialThreshold = nodeWidth;
                            if (std::abs(insR - nodeR) < radialThreshold) {
                                // Check angular proximity - is the insulation in the direction of this face?
                                double angleDiff = insAngle - nodeAngle;
                                while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                                while (angleDiff < -M_PI) angleDiff += 2 * M_PI;
                                
                                bool isLeft = (face == ThermalNodeFace::TANGENTIAL_LEFT);
                                double distAlongTangent = std::abs(angleDiff) * nodeR;
                                
                                // If insulation is in this tangential direction, block convection
                                // TANGENTIAL_LEFT: insulation at larger angle (positive angleDiff)
                                // TANGENTIAL_RIGHT: insulation at smaller angle (negative angleDiff)
                                if (isLeft && angleDiff > 0 && distAlongTangent < maxConvectionDist) {
                                    isExposed = false;
                                    break;
                                }
                                if (!isLeft && angleDiff < 0 && distAlongTangent < maxConvectionDist) {
                                    isExposed = false;
                                    break;
                                }
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
                        r.area = q->surfaceArea;  // Store area for forced convection calculation
                        if (q->coating.has_value()) {
                            r.resistance += WireCoatingUtils::calculateCoatingResistance(q->coating.value(), q->surfaceArea);
                        }
                        _resistances.push_back(r);
                    } else if (THERMAL_DEBUG) {
                    }
                }
            }
        }
        
        // Core convection - connect exposed core quadrants
        // When insulation layers are present, core is covered by insulation and has no convection
        if (!hasInsulationLayers) {
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
                        r.area = q->surfaceArea;  // Store area for forced convection calculation
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
        }
        
        // Insulation layer convection - only the outermost insulation layer has convection to ambient
        // Find the outermost insulation layer (highest layer index)
        int maxLayerIdx = -1;
        for (size_t i = 0; i < _nodes.size(); i++) {
            if (_nodes[i].part == ThermalNodePartType::INSULATION_LAYER) {
                if (!_nodes[i].insulationLayerIndex.has_value()) {
                    throw std::runtime_error("Temperature::createToroidalConvectionConnections: Insulation layer node " +
                                             _nodes[i].name + " is missing insulationLayerIndex.");
                }
                int layerIdx = _nodes[i].insulationLayerIndex.value();
                if (layerIdx > maxLayerIdx) {
                    maxLayerIdx = layerIdx;
                }
            }
        }

        // Only create convection for the outermost insulation layer
        for (size_t i = 0; i < _nodes.size(); i++) {
            if (_nodes[i].part != ThermalNodePartType::INSULATION_LAYER) continue;

            const auto& node = _nodes[i];
            if (!node.insulationLayerIndex.has_value()) {
                throw std::runtime_error("Temperature::createToroidalConvectionConnections: Insulation layer node " +
                                         node.name + " is missing insulationLayerIndex.");
            }
            int layerIdx = node.insulationLayerIndex.value();

            // Skip if not the outermost layer
            if (layerIdx != maxLayerIdx) continue;
            
            // Determine if this is inner or outer insulation node by name
            bool isOuterNode = (node.name.find("_o") != std::string::npos);
            bool isInnerNode = (node.name.find("_i") != std::string::npos);
            
            // Outer nodes: connect RADIAL_OUTER face to ambient (exposed to air)
            if (isOuterNode) {
                auto* q = node.getQuadrant(ThermalNodeFace::RADIAL_OUTER);
                if (q && q->surfaceArea > 0) {
                    ThermalResistanceElement r;
                    r.nodeFromId = i;
                    r.quadrantFrom = ThermalNodeFace::RADIAL_OUTER;
                    r.nodeToId = ambientIdx;
                    r.quadrantTo = ThermalNodeFace::NONE;
                    r.type = _config.includeForcedConvection ? 
                             HeatTransferType::FORCED_CONVECTION : 
                             HeatTransferType::NATURAL_CONVECTION;
                    r.resistance = q->calculateConvectionResistance(h_conv);
                    r.area = q->surfaceArea;  // Store area for forced convection calculation
                    _resistances.push_back(r);
                    if (THERMAL_DEBUG) {
                    }
                } else if (THERMAL_DEBUG) {
                }
            }
            
            // Inner nodes: connect RADIAL_INNER face to ambient (through core hole)
            if (isInnerNode) {
                auto* q = node.getQuadrant(ThermalNodeFace::RADIAL_INNER);
                if (q && q->surfaceArea > 0) {
                    ThermalResistanceElement r;
                    r.nodeFromId = i;
                    r.quadrantFrom = ThermalNodeFace::RADIAL_INNER;
                    r.nodeToId = ambientIdx;
                    r.quadrantTo = ThermalNodeFace::NONE;
                    r.type = _config.includeForcedConvection ? 
                             HeatTransferType::FORCED_CONVECTION : 
                             HeatTransferType::NATURAL_CONVECTION;
                    r.resistance = q->calculateConvectionResistance(h_conv);
                    r.area = q->surfaceArea;  // Store area for forced convection calculation
                    _resistances.push_back(r);
                    if (THERMAL_DEBUG) {
                    }
                } else if (THERMAL_DEBUG) {
                }
            }
        }
}

// IMP-4: Planar convection connections
void Temperature::createPlanarConvectionConnections(size_t ambientIdx, double h_conv) {
        // ============================================================================
        // Planar case: Turn quadrants connect to closest FR4 layer (conduction),
        // except top-most quadrants of top-most turns and bottom-most quadrants
        // of bottom-most turns, which connect to Tamb (convection).
        // Rule: at most ONE connection per quadrant.
        // ============================================================================

        // 1. Find FR4 insulation layer nodes
        // For planar configurations, all insulation layers are PCB substrates (FR4)
        std::vector<size_t> fr4LayerIndices;
        for (size_t i = 0; i < _nodes.size(); i++) {
            if (_nodes[i].part == ThermalNodePartType::INSULATION_LAYER) {
                fr4LayerIndices.push_back(i);
            }
        }

        // 2. Helper lambda: Find turns at LOWEST Y-coordinate (top PCB layer)
        //    Note: SVG Y-axis goes top-to-bottom, so lower Y = top of image
        auto findTopLayerTurns = [&]() -> std::vector<size_t> {
            std::vector<size_t> topTurns;
            double minY = 1e9;

            // First pass: find LOWEST Y-coordinate (top in SVG display)
            for (size_t i = 0; i < _nodes.size(); i++) {
                if (_nodes[i].part == ThermalNodePartType::TURN) {
                    double y = _nodes[i].physicalCoordinates[1];
                    if (y < minY) minY = y;
                }
            }

            // Second pass: collect all turns at this Y-coordinate (within tolerance)
            double yTolerance = 1e-5;  // Very small tolerance to detect distinct layers
            for (size_t i = 0; i < _nodes.size(); i++) {
                if (_nodes[i].part == ThermalNodePartType::TURN) {
                    double y = _nodes[i].physicalCoordinates[1];
                    if (std::abs(y - minY) < yTolerance) {
                        topTurns.push_back(i);
                    }
                }
            }

            return topTurns;
        };

        // 3. Helper lambda: Find turns at HIGHEST Y-coordinate (bottom PCB layer)
        //    Note: SVG Y-axis goes top-to-bottom, so higher Y = bottom of image
        auto findBottomLayerTurns = [&]() -> std::vector<size_t> {
            std::vector<size_t> bottomTurns;
            double maxY = -1e9;

            // First pass: find HIGHEST Y-coordinate (bottom in SVG display)
            for (size_t i = 0; i < _nodes.size(); i++) {
                if (_nodes[i].part == ThermalNodePartType::TURN) {
                    double y = _nodes[i].physicalCoordinates[1];
                    if (y > maxY) maxY = y;
                }
            }

            // Second pass: collect all turns at this Y-coordinate (within tolerance)
            double yTolerance = 1e-5;
            for (size_t i = 0; i < _nodes.size(); i++) {
                if (_nodes[i].part == ThermalNodePartType::TURN) {
                    double y = _nodes[i].physicalCoordinates[1];
                    if (std::abs(y - maxY) < yTolerance) {
                        bottomTurns.push_back(i);
                    }
                }
            }
            return bottomTurns;
        };

        // Get top and bottom layer turns
        std::vector<size_t> topLayerTurnIndices = findTopLayerTurns();
        std::vector<size_t> bottomLayerTurnIndices = findBottomLayerTurns();

        // Create a set for fast lookup
        std::set<size_t> topLayerSet(topLayerTurnIndices.begin(), topLayerTurnIndices.end());
        std::set<size_t> bottomLayerSet(bottomLayerTurnIndices.begin(), bottomLayerTurnIndices.end());

        // 4. For each turn, connect quadrants:
        //    - Top quadrant (TANGENTIAL_LEFT) of top layer turns → Ambient (convection)
        //    - Bottom quadrant (TANGENTIAL_RIGHT) of bottom layer turns → Ambient (convection)
        //    - All other quadrants → nearest FR4 quadrant (conduction)
        for (size_t i = 0; i < _nodes.size(); i++) {
            if (_nodes[i].part != ThermalNodePartType::TURN) continue;

            bool isTopLayerTurn = topLayerSet.count(i) > 0;
            bool isBottomLayerTurn = bottomLayerSet.count(i) > 0;

            for (int qIdx = 0; qIdx < 4; ++qIdx) {
                ThermalNodeFace face = _nodes[i].quadrants[qIdx].face;
                if (face == ThermalNodeFace::NONE) continue;

                auto* turnQuad = _nodes[i].getQuadrant(face);
                if (!turnQuad || turnQuad->surfaceArea <= 0) continue;

                // Exception 1: Top quadrant of top layer turns connects to ambient
                if (isTopLayerTurn && face == ThermalNodeFace::TANGENTIAL_LEFT) {
                    ThermalResistanceElement r;
                    r.nodeFromId = i;
                    r.quadrantFrom = face;
                    r.nodeToId = ambientIdx;
                    r.quadrantTo = ThermalNodeFace::NONE;
                    r.type = _config.includeForcedConvection ?
                        HeatTransferType::FORCED_CONVECTION :
                        HeatTransferType::NATURAL_CONVECTION;
                    r.resistance = turnQuad->calculateConvectionResistance(h_conv);
                    r.area = turnQuad->surfaceArea;
                    _resistances.push_back(r);

                    continue;  // Skip FR4 connection for this quadrant
                }

                // Exception 2: Bottom quadrant of bottom layer turns connects to ambient
                if (isBottomLayerTurn && face == ThermalNodeFace::TANGENTIAL_RIGHT) {
                    ThermalResistanceElement r;
                    r.nodeFromId = i;
                    r.quadrantFrom = face;
                    r.nodeToId = ambientIdx;
                    r.quadrantTo = ThermalNodeFace::NONE;
                    r.type = _config.includeForcedConvection ?
                        HeatTransferType::FORCED_CONVECTION :
                        HeatTransferType::NATURAL_CONVECTION;
                    r.resistance = turnQuad->calculateConvectionResistance(h_conv);
                    r.area = turnQuad->surfaceArea;
                    _resistances.push_back(r);

                    continue;  // Skip FR4 connection for this quadrant
                }

                // Default: Connect this quadrant to the closest FR4 layer quadrant via conduction
                if (!fr4LayerIndices.empty()) {
                    // No distance limit — find the nearest FR4 quadrant
                    size_t closestFr4Node = SIZE_MAX;
                    ThermalNodeFace closestFr4Quadrant = ThermalNodeFace::NONE;
                    double minDist = 1e9;

                    // Get turn quadrant limit coordinates for distance calculation
                    auto turnLimitCoords = turnQuad->limitCoordinates;

                    // Search all FR4 nodes and their quadrants
                    for (size_t fr4Idx : fr4LayerIndices) {
                        const auto& fr4Node = _nodes[fr4Idx];

                        // Check all 4 quadrants of this FR4 layer node
                        for (int qIdx = 0; qIdx < 4; ++qIdx) {
                            ThermalNodeFace fr4Face = fr4Node.quadrants[qIdx].face;
                            if (fr4Face == ThermalNodeFace::NONE) continue;

                            auto* fr4Quad = fr4Node.getQuadrant(fr4Face);
                            if (!fr4Quad || fr4Quad->surfaceArea <= 0) continue;

                            // Calculate distance between turn quadrant and FR4 quadrant
                            auto fr4LimitCoords = fr4Quad->limitCoordinates;
                            double dx = turnLimitCoords[0] - fr4LimitCoords[0];
                            double dy = turnLimitCoords[1] - fr4LimitCoords[1];
                            double dz = turnLimitCoords[2] - fr4LimitCoords[2];
                            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                            if (dist < minDist) {
                                minDist = dist;
                                closestFr4Node = fr4Idx;
                                closestFr4Quadrant = fr4Face;
                            }
                        }
                    }

                    if (closestFr4Node != SIZE_MAX && closestFr4Quadrant != ThermalNodeFace::NONE) {
                        ThermalResistanceElement r;
                        r.nodeFromId = i;
                        r.quadrantFrom = face;
                        r.nodeToId = closestFr4Node;
                        r.quadrantTo = closestFr4Quadrant;  // Connect to specific FR4 quadrant!
                        r.type = HeatTransferType::CONDUCTION;
                        double contactArea = turnQuad->surfaceArea;
                        double thickness = std::max(minDist, 1e-6);
                        double k = kFR4_ThermalConductivity; // FR4 through-plane thermal conductivity W/(m·K)
                        r.resistance = ThermalResistance::calculateConductionResistance(
                            thickness, k, contactArea);
                        r.area = contactArea;
                        _resistances.push_back(r);

                    }
                }
            } // end quadrant loop
        } // end turn loop

        // 5. Connect FR4 insulation layer quadrants to ambient (for heat dissipation)
        // Only topmost and bottommost FR4 layers are exposed to ambient
        // Surface area is reduced by copper coverage from top/bottom turns

        // Find topmost and bottommost FR4 layers by Y-coordinate
        double fr4MaxY = -1e9;
        double fr4MinY = 1e9;
        for (size_t fr4Idx : fr4LayerIndices) {
            double fr4Y = _nodes[fr4Idx].physicalCoordinates[1];
            if (fr4Y > fr4MaxY) fr4MaxY = fr4Y;
            if (fr4Y < fr4MinY) fr4MinY = fr4Y;
        }

        // Calculate copper coverage area on top and bottom surfaces
        // Top surface: covered by top layer turns
        // Bottom surface: covered by bottom layer turns
        double topCopperArea = 0.0;
        double bottomCopperArea = 0.0;

        for (size_t turnIdx : topLayerTurnIndices) {
            // Top layer turns cover the FR4 bottom surface (they're on top, blocking bottom FR4 surface)
            // Actually, in planar, top layer is at minY (top in SVG), so they sit on the topmost FR4
            // The copper area is the turn's cross-sectional area in the XZ plane
            if (_nodes[turnIdx].dimensions.width > 0 && _nodes[turnIdx].dimensions.height > 0) {
                topCopperArea += _nodes[turnIdx].dimensions.width * _nodes[turnIdx].dimensions.height;
            }
        }

        for (size_t turnIdx : bottomLayerTurnIndices) {
            // Bottom layer turns cover the FR4 top surface
            if (_nodes[turnIdx].dimensions.width > 0 && _nodes[turnIdx].dimensions.height > 0) {
                bottomCopperArea += _nodes[turnIdx].dimensions.width * _nodes[turnIdx].dimensions.height;
            }
        }

        double fr4Tolerance = 1e-5;  // Same small tolerance to detect distinct layers
        for (size_t fr4Idx : fr4LayerIndices) {
            const auto& fr4Node = _nodes[fr4Idx];
            double fr4Y = fr4Node.physicalCoordinates[1];

            bool isTopMostFR4 = (std::abs(fr4Y - fr4MaxY) < fr4Tolerance);
            bool isBottomMostFR4 = (std::abs(fr4Y - fr4MinY) < fr4Tolerance);

            // Connect top surface of topmost FR4 layer to ambient
            // This surface is covered by bottom layer turns
            if (isTopMostFR4) {
                auto* topQuad = fr4Node.getQuadrant(ThermalNodeFace::TANGENTIAL_LEFT);
                if (topQuad && topQuad->surfaceArea > 0) {
                    // Calculate effective convection area (FR4 area - copper coverage)
                    double effectiveArea = std::max(topQuad->surfaceArea - bottomCopperArea, topQuad->surfaceArea * 0.1);
                    double areaProportion = effectiveArea / topQuad->surfaceArea;

                    ThermalResistanceElement r;
                    r.nodeFromId = fr4Idx;
                    r.quadrantFrom = ThermalNodeFace::TANGENTIAL_LEFT;
                    r.nodeToId = ambientIdx;
                    r.quadrantTo = ThermalNodeFace::NONE;
                    r.type = _config.includeForcedConvection ?
                        HeatTransferType::FORCED_CONVECTION :
                        HeatTransferType::NATURAL_CONVECTION;
                    // Resistance scales inversely with area
                    r.resistance = topQuad->calculateConvectionResistance(h_conv) / areaProportion;
                    r.area = effectiveArea;
                    _resistances.push_back(r);

                }
            }

            // Connect bottom surface of bottommost FR4 layer to ambient
            // This surface is covered by top layer turns (remember: minY = top in SVG = topmost turns)
            if (isBottomMostFR4) {
                auto* bottomQuad = fr4Node.getQuadrant(ThermalNodeFace::TANGENTIAL_RIGHT);
                if (bottomQuad && bottomQuad->surfaceArea > 0) {
                    // Calculate effective convection area (FR4 area - copper coverage)
                    double effectiveArea = std::max(bottomQuad->surfaceArea - topCopperArea, bottomQuad->surfaceArea * 0.1);
                    double areaProportion = effectiveArea / bottomQuad->surfaceArea;

                    ThermalResistanceElement r;
                    r.nodeFromId = fr4Idx;
                    r.quadrantFrom = ThermalNodeFace::TANGENTIAL_RIGHT;
                    r.nodeToId = ambientIdx;
                    r.quadrantTo = ThermalNodeFace::NONE;
                    r.type = _config.includeForcedConvection ?
                        HeatTransferType::FORCED_CONVECTION :
                        HeatTransferType::NATURAL_CONVECTION;
                    // Resistance scales inversely with area
                    r.resistance = bottomQuad->calculateConvectionResistance(h_conv) / areaProportion;
                    r.area = effectiveArea;
                    _resistances.push_back(r);

                }
            }
        }

        // 6. Connect core nodes to ambient via convection
        // In planar designs, core surfaces are exposed to air
        for (size_t i = 0; i < _nodes.size(); i++) {
            if (_nodes[i].part == ThermalNodePartType::AMBIENT) continue;

            bool isCoreNode = (_nodes[i].part == ThermalNodePartType::CORE_CENTRAL_COLUMN ||
                              _nodes[i].part == ThermalNodePartType::CORE_LATERAL_COLUMN ||
                              _nodes[i].part == ThermalNodePartType::CORE_TOP_YOKE ||
                              _nodes[i].part == ThermalNodePartType::CORE_BOTTOM_YOKE);

            if (!isCoreNode) continue;

            // Build connected quadrants set for planar core convection check
            std::set<std::string> planarConnectedQuadrants;
            for (const auto& res : _resistances) {
                if (res.type == HeatTransferType::CONDUCTION) {
                    planarConnectedQuadrants.insert(std::to_string(res.nodeFromId) + "_" + std::string(magic_enum::enum_name(res.quadrantFrom)));
                    planarConnectedQuadrants.insert(std::to_string(res.nodeToId) + "_" + std::string(magic_enum::enum_name(res.quadrantTo)));
                }
            }

            // Connect only exposed (non-conduction-connected) core quadrants to ambient
            for (int qIdx = 0; qIdx < 4; ++qIdx) {
                ThermalNodeFace face = _nodes[i].quadrants[qIdx].face;
                if (face == ThermalNodeFace::NONE) continue;

                // Skip quadrants already connected by conduction (e.g., facing bobbin/winding)
                std::string qKey = std::to_string(i) + "_" + std::string(magic_enum::enum_name(face));
                if (planarConnectedQuadrants.count(qKey) > 0) continue;

                auto* coreQuad = _nodes[i].getQuadrant(face);
                if (!coreQuad || coreQuad->surfaceArea <= 0) continue;

                ThermalResistanceElement r;
                r.nodeFromId = i;
                r.quadrantFrom = face;
                r.nodeToId = ambientIdx;
                r.quadrantTo = ThermalNodeFace::NONE;
                r.type = _config.includeForcedConvection ?
                    HeatTransferType::FORCED_CONVECTION :
                    HeatTransferType::NATURAL_CONVECTION;
                r.resistance = coreQuad->calculateConvectionResistance(h_conv);
                r.area = coreQuad->surfaceArea;
                _resistances.push_back(r);

            }
        }
}

// IMP-4: Concentric convection connections
void Temperature::createConcentricConvectionConnections(size_t ambientIdx, double h_conv) {
    // Track initial resistance count for symmetry correction at the end
    size_t initialResistanceCount = _resistances.size();
    
    // Build connected quadrants set ONCE (used throughout this method)
    std::set<std::string> connectedQuadrants;
    for (const auto& res : _resistances) {
        if (res.type == HeatTransferType::CONDUCTION) {
            connectedQuadrants.insert(std::to_string(res.nodeFromId) + "_" + std::string(magic_enum::enum_name(res.quadrantFrom)));
            connectedQuadrants.insert(std::to_string(res.nodeToId) + "_" + std::string(magic_enum::enum_name(res.quadrantTo)));
        }
    }

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
            // connectedQuadrants already built above - reuse it
            
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
                        r.area = q->surfaceArea;  // Store area for forced convection calculation
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
                double surfaceArea = _nodes[i].getTotalSurfaceArea();
                r.resistance = ThermalResistance::calculateConvectionResistance(h_conv, surfaceArea);
                r.area = surfaceArea;  // Store area for forced convection calculation
                _resistances.push_back(r);
            }
        }
    
    // Half-core symmetry correction: the model uses half the core depth,
    // so convection areas represent only one half. The real core has double
    // the convection area (both halves), so halve all convection resistances
    // created in this method (R = 1/(h*A), double A => R/2)
    for (size_t i = initialResistanceCount; i < _resistances.size(); ++i) {
        if (_resistances[i].type == HeatTransferType::NATURAL_CONVECTION ||
            _resistances[i].type == HeatTransferType::FORCED_CONVECTION) {
            _resistances[i].resistance /= 2.0;
            _resistances[i].area *= 2.0;
        }
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

double Temperature::calculateSurfaceDistance(const ThermalNetworkNode& node1, 
                                              const ThermalNetworkNode& node2) const {
    // Calculate true surface-to-surface distance by subtracting wire extents
    double dx = node1.physicalCoordinates[0] - node2.physicalCoordinates[0];
    double dy = node1.physicalCoordinates[1] - node2.physicalCoordinates[1];
    double dz = (node1.physicalCoordinates.size() > 2 && node2.physicalCoordinates.size() > 2)
                ? (node1.physicalCoordinates[2] - node2.physicalCoordinates[2]) : 0.0;
    double centerDistance = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // Subtract extents (half of minimum dimension) for each node
    double extent1 = std::min(node1.dimensions.width, node1.dimensions.height) / 2.0;
    double extent2 = std::min(node2.dimensions.width, node2.dimensions.height) / 2.0;
    double surfaceDistance = centerDistance - extent1 - extent2;
    
    return std::max(surfaceDistance, 0.0);
}

bool Temperature::shouldConnectQuadrants(const ThermalNetworkNode& node1, ThermalNodeFace face1,
                                          const ThermalNetworkNode& node2, ThermalNodeFace face2) const {
    // Use evaluateTurnPairDistance for consistent threshold logic
    auto geom = evaluateTurnPairDistance(node1, node2);
    return geom.shouldConnect;
}

double Temperature::calculateContactArea(const ThermalNodeQuadrant& q1, const ThermalNodeQuadrant& q2) const {
    if (q1.length <= 1e-9) {
        throw std::runtime_error("Temperature::calculateContactArea: Quadrant length is effectively zero (" +
                                 std::to_string(q1.length) + "). Cannot calculate contact area.");
    }
    if (q2.length <= 1e-9) {
        throw std::runtime_error("Temperature::calculateContactArea: Quadrant length is effectively zero (" +
                                 std::to_string(q2.length) + "). Cannot calculate contact area.");
    }
    double minLength = std::min(q1.length, q2.length);
    double height = q1.surfaceArea / q1.length;
    return height * minLength;
}

double Temperature::getInsulationLayerThermalResistance(int turnIdx1, int turnIdx2, double contactArea) {
    if (contactArea <= 0.0) {
        throw std::runtime_error("Temperature::getInsulationLayerThermalResistance: Contact area is zero or negative (" +
                                 std::to_string(contactArea) + "). Cannot calculate thermal resistance.");
    }
    
    try {
        auto coil = _magnetic.get_coil();
        auto turnsDescription = coil.get_turns_description();
        
        // For turn-to-solid connections (turnIdx2 = -1)
        if (turnIdx2 < 0) {
            double insulationThickness = kWire_DefaultEnamelThickness;
            double insulationK = kWire_DefaultEnamelConductivity;
            
            double resistance = ThermalResistance::calculateConductionResistance(
                insulationThickness, insulationK, contactArea);
            
            return std::max(resistance, 1e-6);  // Physically plausible minimum
        }
        
        if (!turnsDescription || turnIdx1 < 0 || turnIdx1 >= static_cast<int>(turnsDescription->size()) || 
            turnIdx2 >= static_cast<int>(turnsDescription->size())) {
            // Return resistance based on typical enamel insulation rather than near-zero
            double defaultResistance = ThermalResistance::calculateConductionResistance(
                kWire_DefaultEnamelThickness, kWire_DefaultEnamelConductivity, contactArea);
            return std::max(defaultResistance, 1e-6);
        }
        
        auto& turn1 = (*turnsDescription)[turnIdx1];
        auto& turn2 = (*turnsDescription)[turnIdx2];
        
        auto layersBetween = StrayCapacitance::get_insulation_layers_between_two_turns(turn1, turn2, coil);
        
        double totalLayerResistance = 0.0;
        for (const auto& layer : layersBetween) {
            double layerThickness = coil.get_insulation_layer_thickness(layer);
            double layerK = kInsulation_DefaultConductivity;
            
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
    } catch (const std::exception& e) { // IMP-8
        // Return resistance based on typical enamel insulation rather than near-zero
        if (THERMAL_DEBUG) {
            std::cout << "Temperature: wire insulation resistance calc failed, using enamel default: "
                      << e.what() << std::endl;
        }
        double defaultResistance = ThermalResistance::calculateConductionResistance(
            kWire_DefaultEnamelThickness, kWire_DefaultEnamelConductivity, contactArea);
        return std::max(defaultResistance, 1e-6);
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
    }
    
    // Find bounding box of all physical coordinates
    // Skip insulation layer nodes to avoid excessive scaling from widely spaced layers
    double minX = 1e9, maxX = -1e9;
    double minY = 1e9, maxY = -1e9;
    
    for (const auto& node : _nodes) {
        if (node.physicalCoordinates.size() >= 2) {
            // Skip insulation layer nodes for scaling calculation
            if (node.part == ThermalNodePartType::INSULATION_LAYER) {
                continue;
            }
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


// IMP-5: Recalculate temperature-dependent convection/radiation resistances
void Temperature::recalculateConvectionResistances(const std::vector<double>& temperatures) {
    for (auto& res : _resistances) {
        if (res.nodeFromId >= temperatures.size()) continue;
        double surfaceTemp = temperatures[res.nodeFromId];
        double ambientTemp = _config.ambientTemperature;
        if (res.type == HeatTransferType::NATURAL_CONVECTION && res.area > 0) {
            double charLength = std::sqrt(std::max(res.area, 1e-9));
            double h_new = ThermalResistance::calculateNaturalConvectionCoefficient(
                // IMP-NEW-04: Use stored orientation instead of always VERTICAL
                // WHY: Horizontal surfaces have very different h values.
                // PREVIOUS: SurfaceOrientation::VERTICAL (hardcoded)
                surfaceTemp, ambientTemp, charLength, res.orientation);
            if (h_new > 0) res.resistance = 1.0 / (h_new * res.area);
        } else if (res.type == HeatTransferType::FORCED_CONVECTION && res.area > 0) {
            double charLength = std::sqrt(std::max(res.area, 1e-9));
            double h_new = ThermalResistance::calculateForcedConvectionCoefficient(
                _config.airVelocity, charLength, surfaceTemp);
            if (h_new > 0) res.resistance = 1.0 / (h_new * res.area);
        } else if (res.type == HeatTransferType::RADIATION && res.area > 0) {
            double h_rad = ThermalResistance::calculateRadiationCoefficient(
                surfaceTemp, ambientTemp, _config.surfaceEmissivity);
            if (h_rad > 0) res.resistance = 1.0 / (h_rad * res.area);
        }
    }
}

ThermalResult Temperature::solveThermalCircuit() {
    size_t n = _nodes.size();
    if (n == 0) {
        throw std::runtime_error("Temperature::solveThermalCircuit: Thermal network has zero nodes. "
                                 "The thermal model was not properly built.");
    }

    // Find the ambient node (look for AMBIENT part type, not just last node)
    size_t ambientIdx = n - 1;  // Default to last node
    for (size_t i = 0; i < n; ++i) {
        if (_nodes[i].isAmbient()) {
            ambientIdx = i;
            break;
        }
    }
    
    std::vector<double> temperatures(n, _config.ambientTemperature);
    std::vector<double> powerInputs(n, 0.0);
    
    for (size_t i = 0; i < n; ++i) {
        powerInputs[i] = _nodes[i].powerDissipation;
    }
    
    size_t iteration = 0;
    bool converged = false;
    std::vector<double> oldTemperatures = temperatures;
    
    while (iteration < _config.maxIterations && !converged) {
        // IMP-5
        if (iteration > 0) { recalculateConvectionResistances(temperatures); }

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

        // Set ambient node as fixed temperature
        G.setRowZero(ambientIdx);
        G(ambientIdx, ambientIdx) = 1.0;
        powerInputs[ambientIdx] = _config.ambientTemperature;
        
        // Set any fixed temperature nodes (cold plate, etc.)
        for (size_t i = 0; i < n; ++i) {
            if (_nodes[i].isFixedTemperature && i != ambientIdx) {
                G.setRowZero(i);
                G(i, i) = 1.0;
                powerInputs[i] = _nodes[i].temperature;
            }
        }
        
        // Diagnose disconnected nodes before attempting solve
        {
            std::vector<size_t> disconnectedNodes;
            for (size_t i = 0; i < n; ++i) {
                if (i == ambientIdx) continue;
                if (_nodes[i].isFixedTemperature) continue;
                if (G(i, i) < 1e-12) {
                    disconnectedNodes.push_back(i);
                }
            }
            if (!disconnectedNodes.empty()) {
                std::string msg = "Temperature::solveThermalCircuit: " + std::to_string(disconnectedNodes.size()) + " disconnected node(s) found (no thermal path to any other node):\n";
                for (size_t idx : disconnectedNodes) {
                    msg += "  Node " + std::to_string(idx) + ": name='" + _nodes[idx].name + "' part=";
                    switch (_nodes[idx].part) {
                        case ThermalNodePartType::TURN:                  msg += "TURN"; break;
                        case ThermalNodePartType::CORE_CENTRAL_COLUMN:   msg += "CORE_CENTRAL_COLUMN"; break;
                        case ThermalNodePartType::CORE_LATERAL_COLUMN:   msg += "CORE_LATERAL_COLUMN"; break;
                        case ThermalNodePartType::CORE_TOP_YOKE:         msg += "CORE_TOP_YOKE"; break;
                        case ThermalNodePartType::CORE_BOTTOM_YOKE:      msg += "CORE_BOTTOM_YOKE"; break;
                        case ThermalNodePartType::CORE_TOROIDAL_SEGMENT: msg += "CORE_TOROIDAL_SEGMENT"; break;
                        case ThermalNodePartType::AMBIENT:               msg += "AMBIENT"; break;
                        default:                                          msg += "OTHER"; break;
                    }
                    msg += " power=" + std::to_string(_nodes[idx].powerDissipation) + "W\n";
                }
                throw std::runtime_error(msg);
            }
        }

        try {
            // IMP-10: Non-throwing solve
            bool solveSuccess = true;
            temperatures = SimpleMatrix::solve(G, powerInputs, solveSuccess);
            if (!solveSuccess) {
                throw std::runtime_error("Temperature::solveThermalCircuit: Thermal matrix is singular (not solvable). "
                                         "This indicates disconnected thermal nodes or missing boundary conditions.");
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("Temperature::solveThermalCircuit: Solver failed with exception: " +
                                     std::string(e.what()));
        }
        
        // Check for NaN or infinite temperatures
        bool hasInvalidTemps = false;
        for (size_t i = 0; i < n; ++i) {
            if (!std::isfinite(temperatures[i])) {
                hasInvalidTemps = true;
            }
        }
        if (hasInvalidTemps) {
            throw std::runtime_error("Temperature::solveThermalCircuit: Solver produced NaN or infinite temperatures at iteration " +
                                     std::to_string(iteration) + ". This indicates a numerical instability in the thermal network.");
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
        if (THERMAL_DEBUG && _nodes[i].part == ThermalNodePartType::TURN) {
            std::cout << "  " << _nodes[i].name << ": " << temperatures[i] << "°C" << std::endl;
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
            _nodes[i].part == ThermalNodePartType::CORE_LATERAL_COLUMN ||
            _nodes[i].part == ThermalNodePartType::CORE_TOP_YOKE ||
            _nodes[i].part == ThermalNodePartType::CORE_BOTTOM_YOKE) {
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

std::map<std::string, double> Temperature::getTemperaturesByComponentType() const {
    std::map<std::string, double> result;
    
    for (const auto& node : _nodes) {
        if (node.part == ThermalNodePartType::AMBIENT) continue;
        
        std::string category;
        switch (node.part) {
            case ThermalNodePartType::CORE_CENTRAL_COLUMN:
            case ThermalNodePartType::CORE_LATERAL_COLUMN:
            case ThermalNodePartType::CORE_TOP_YOKE:
            case ThermalNodePartType::CORE_BOTTOM_YOKE:
            case ThermalNodePartType::CORE_TOROIDAL_SEGMENT:
                category = "core";
                break;
            case ThermalNodePartType::BOBBIN_CENTRAL_COLUMN:
            case ThermalNodePartType::BOBBIN_TOP_YOKE:
            case ThermalNodePartType::BOBBIN_BOTTOM_YOKE:
                category = "bobbin";
                break;
            case ThermalNodePartType::TURN:
                if (node.windingIndex.has_value()) {
                    category = "winding " + std::to_string(node.windingIndex.value());
                } else {
                    category = "winding";
                }
                break;
            case ThermalNodePartType::INSULATION_LAYER:
                category = "insulation";
                break;
            case ThermalNodePartType::HEATSINK:
                category = "heatsink";
                break;
            case ThermalNodePartType::COLD_PLATE:
                category = "cold plate";
                break;
            default:
                continue;
        }
        
        auto it = result.find(category);
        if (it == result.end()) {
            result[category] = node.temperature;
        } else {
            it->second = std::max(it->second, node.temperature);
        }
    }
    
    return result;
}

std::map<std::string, double> Temperature::getTemperaturePerTurn() const {
    // For toroidal cores, a single turn is split into inner/outer nodes
    // (e.g. "Turn_W0_T3_Inner", "Turn_W0_T3_Outer"). We merge them by
    // stripping the suffix and taking the maximum.
    std::map<std::string, double> result;
    
    for (const auto& node : _nodes) {
        if (node.part != ThermalNodePartType::TURN) continue;
        
        std::string key = node.name;
        // Strip "_Inner" or "_Outer" suffix to group toroidal half-nodes
        const std::string innerSuffix = "_Inner";
        const std::string outerSuffix = "_Outer";
        if (key.size() > innerSuffix.size() &&
            key.compare(key.size() - innerSuffix.size(), innerSuffix.size(), innerSuffix) == 0) {
            key = key.substr(0, key.size() - innerSuffix.size());
        } else if (key.size() > outerSuffix.size() &&
                   key.compare(key.size() - outerSuffix.size(), outerSuffix.size(), outerSuffix) == 0) {
            key = key.substr(0, key.size() - outerSuffix.size());
        }
        
        auto it = result.find(key);
        if (it == result.end()) {
            result[key] = node.temperature;
        } else {
            it->second = std::max(it->second, node.temperature);
        }
    }
    
    return result;
}

double Temperature::getTemperatureAtPoint(const std::vector<double>& point) const {
    if (point.size() < 2) {
        throw std::runtime_error("Temperature::getTemperatureAtPoint: Point must have at least 2 coordinates (x, y).");
    }
    
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

// ============================================================================
// Cooling Application Methods
// ============================================================================

void Temperature::applyMasCooling(const MAS::Cooling& cooling) {
    auto type = CoolingUtils::detectCoolingType(cooling);
    
    if (THERMAL_DEBUG) {
    }
    
    switch (type) {
        case CoolingType::FORCED_CONVECTION:
            applyForcedConvection(cooling);
            break;
        case CoolingType::HEATSINK:
            applyHeatsinkCooling(cooling);
            break;
        case CoolingType::COLD_PLATE:
            applyColdPlateCooling(cooling);
            break;
        case CoolingType::NATURAL_CONVECTION:
            // Default behavior, no special handling needed
            break;
        case CoolingType::UNKNOWN:
            std::cerr << "Warning: Unknown cooling type in MAS::Cooling" << std::endl;
            break;
    }
}

void Temperature::applyForcedConvection(const MAS::Cooling& cooling) {
    if (!cooling.get_velocity().has_value() || cooling.get_velocity().value().empty()) {
        throw std::runtime_error("Temperature::applyForcedConvection: Forced convection requested but velocity is missing or empty.");
    }
    
    double velocity = cooling.get_velocity().value()[0]; // m/s
    
    // Store velocity in config so recalculateConvectionResistances can use it
    _config.airVelocity = velocity;
    _config.includeForcedConvection = true;
    
    if (THERMAL_DEBUG) {
    }
    
    for (auto& resistance : _resistances) {
        if (resistance.type == HeatTransferType::NATURAL_CONVECTION) {
            // Skip resistances with zero or invalid area
            if (resistance.area <= 1e-15) {
                if (THERMAL_DEBUG) {
                    std::cerr << "Warning: Skipping forced convection for resistance with invalid area: " << resistance.area << std::endl;
                }
                continue;
            }
            
            // Get current surface temperature and ambient
            double surfaceTemp = _nodes[resistance.nodeFromId].temperature;
            double ambientTemp = _config.ambientTemperature;
            
            // Calculate forced convection coefficient
            double charLength = std::sqrt(resistance.area);
            double h_forced = CoolingUtils::calculateForcedConvectionCoefficient(
                surfaceTemp, ambientTemp, velocity, charLength);
            
            // Calculate natural convection coefficient (from existing resistance)
            double h_natural = 1.0 / (std::max(resistance.resistance, 1e-9) * resistance.area);
            
            // Mixed convection formula
            double h_total = CoolingUtils::calculateMixedConvectionCoefficient(h_natural, h_forced);
            
            // Update resistance (ensure we don't get infinity)
            if (h_total > 1e-9) {
                resistance.resistance = 1.0 / (h_total * resistance.area);
                resistance.type = HeatTransferType::FORCED_CONVECTION;
            }
        }
    }
}

void Temperature::applyHeatsinkCooling(const MAS::Cooling& cooling) {
    // NOTE: Heatsink currently always mounts on top yoke (for concentric cores).
    // Future improvement: add heatsinkPlacement config option (TOP, BOTTOM, SIDE)
    // to support different mounting orientations.
    // Find top yoke node (for concentric cores)
    size_t topYokeIdx = static_cast<size_t>(-1);
    for (size_t i = 0; i < _nodes.size(); ++i) {
        if (_nodes[i].part == ThermalNodePartType::CORE_TOP_YOKE ||
            _nodes[i].part == ThermalNodePartType::BOBBIN_TOP_YOKE) {
            topYokeIdx = i;
            break;
        }
    }
    
    if (topYokeIdx == static_cast<size_t>(-1)) {
        throw std::runtime_error("Temperature::applyHeatsinkCooling: No top yoke node found for heatsink mounting. "
                                 "Heatsink cooling cannot be applied without a top yoke or bobbin top yoke node.");
    }
    
    if (THERMAL_DEBUG) {
    }
    
    // Create heatsink node
    ThermalNetworkNode heatsinkNode;
    heatsinkNode.name = "Heatsink";
    heatsinkNode.part = ThermalNodePartType::HEATSINK;
    heatsinkNode.temperature = _config.ambientTemperature;
    heatsinkNode.physicalCoordinates = _nodes[topYokeIdx].physicalCoordinates;
    heatsinkNode.physicalCoordinates[1] += 0.02;  // 20mm above top yoke
    
    size_t heatsinkIdx = _nodes.size();
    _nodes.push_back(heatsinkNode);
    
    // Create TIM resistance if interface properties provided
    double timResistance = kTIM_DefaultResistance;
    if (cooling.get_interface_thickness().has_value() && 
        cooling.get_interface_thermal_resistance().has_value()) {
        // Issue 25: interface_thermal_resistance is area-specific resistance (K*m^2/W)
        // R_total = R_specific / area
        double areaSpecificResistance = cooling.get_interface_thermal_resistance().value();
        double area = _nodes[topYokeIdx].getTotalSurfaceArea();
        if (area > 0) {
            timResistance = areaSpecificResistance / area;
        }
    }
    
    ThermalResistanceElement timR;
    timR.nodeFromId = topYokeIdx;
    timR.quadrantFrom = ThermalNodeFace::NONE;
    timR.nodeToId = heatsinkIdx;
    timR.quadrantTo = ThermalNodeFace::NONE;
    timR.type = HeatTransferType::CONDUCTION;
    timR.resistance = timResistance;
    timR.area = _nodes[topYokeIdx].getTotalSurfaceArea();
    _resistances.push_back(timR);
    
    // Create heatsink-to-ambient resistance
    double r_heatsink = cooling.get_thermal_resistance().value();  // K/W
    
    size_t ambientIdx = static_cast<size_t>(-1);
    for (size_t i = 0; i < _nodes.size(); ++i) {
        if (_nodes[i].isAmbient()) {
            ambientIdx = i;
            break;
        }
    }
    
    if (ambientIdx != static_cast<size_t>(-1)) {
        ThermalResistanceElement hsR;
        hsR.nodeFromId = heatsinkIdx;
        hsR.quadrantFrom = ThermalNodeFace::NONE;
        hsR.nodeToId = ambientIdx;
        hsR.quadrantTo = ThermalNodeFace::NONE;
        hsR.type = HeatTransferType::HEATSINK_CONVECTION;
        hsR.resistance = r_heatsink;
        hsR.area = 1.0;  // Heatsink thermal resistance is already total
        _resistances.push_back(hsR);
    }
}

void Temperature::applyColdPlateCooling(const MAS::Cooling& cooling) {
    if (!cooling.get_maximum_temperature().has_value()) {
        throw std::runtime_error("Temperature::applyColdPlateCooling: Cold plate cooling requested but maximum_temperature is missing.");
    }

    double coldPlateTemp = cooling.get_maximum_temperature().value();
    
    if (THERMAL_DEBUG) {
    }
    
    // Find bottom surface nodes (core bottom or bobbin bottom for concentric, 
    // bottom segments for toroidal)
    std::vector<size_t> surfaceNodes;
    double minY = 1e9;
    
    // First pass: find minimum Y coordinate of core nodes
    for (size_t i = 0; i < _nodes.size(); ++i) {
        if (_nodes[i].part == ThermalNodePartType::CORE_BOTTOM_YOKE ||
            _nodes[i].part == ThermalNodePartType::BOBBIN_BOTTOM_YOKE ||
            _nodes[i].part == ThermalNodePartType::CORE_TOROIDAL_SEGMENT) {
            minY = std::min(minY, _nodes[i].physicalCoordinates[1]);
        }
    }
    
    // Calculate proximity threshold based on core size
    // For toroidal cores, use 10% of core outer radius as proximity threshold
    double coldPlateProximityThreshold = 0.005;  // Default 5mm
    auto core = _magnetic.get_core();
    if (core.get_shape_family() == CoreShapeFamily::T) {
        auto coreDims = flatten_dimensions(core.resolve_shape().get_dimensions().value());
        if (coreDims.find("A") != coreDims.end()) {
            coldPlateProximityThreshold = 0.1 * coreDims["A"] / 2.0;  // 10% of outer radius
        }
    }
    
    // Second pass: collect nodes near the bottom (within proximity threshold of minY)
    for (size_t i = 0; i < _nodes.size(); ++i) {
        if (_nodes[i].part == ThermalNodePartType::CORE_BOTTOM_YOKE ||
            _nodes[i].part == ThermalNodePartType::BOBBIN_BOTTOM_YOKE) {
            surfaceNodes.push_back(i);
        } else if (_nodes[i].part == ThermalNodePartType::CORE_TOROIDAL_SEGMENT) {
            // For toroidal, include segments near the bottom
            if (std::abs(_nodes[i].physicalCoordinates[1] - minY) < coldPlateProximityThreshold) {
                surfaceNodes.push_back(i);
            }
        }
    }
    
    if (surfaceNodes.empty()) {
        throw std::runtime_error("Temperature::applyColdPlateCooling: No surface nodes found for cold plate attachment. "
                                 "Cold plate cooling requires core bottom yoke, bobbin bottom yoke, or toroidal segments.");
    }
    
    // Create cold plate node with fixed temperature
    // Note: We add this as a regular node (not AMBIENT) with fixed temperature flag
    // The solver will treat it as a fixed temperature boundary
    ThermalNetworkNode coldPlateNode;
    coldPlateNode.name = "ColdPlate";
    coldPlateNode.part = ThermalNodePartType::COLD_PLATE;
    coldPlateNode.temperature = coldPlateTemp;
    coldPlateNode.isFixedTemperature = true;
    
    // Calculate center position of surface nodes
    double avgX = 0.0, avgY = 0.0;
    for (size_t idx : surfaceNodes) {
        avgX += _nodes[idx].physicalCoordinates[0];
        avgY += _nodes[idx].physicalCoordinates[1];
    }
    avgX /= surfaceNodes.size();
    avgY /= surfaceNodes.size();
    
    coldPlateNode.physicalCoordinates = {avgX, avgY - 0.01, 0.0};  // 10mm below
    
    size_t coldPlateIdx = _nodes.size();
    _nodes.push_back(coldPlateNode);
    
    // Connect surface nodes to cold plate
    for (size_t nodeIdx : surfaceNodes) {
        double timResistance = kTIM_DefaultResistance;
        
        if (cooling.get_interface_thickness().has_value() && 
            cooling.get_interface_thermal_resistance().has_value()) {
            // Issue 25: interface_thermal_resistance is area-specific resistance (K*m^2/W)
            // R_total = R_specific / area
            double areaSpecificResistance = cooling.get_interface_thermal_resistance().value();
            double area = _nodes[nodeIdx].getTotalSurfaceArea();
            if (area > 0) {
                timResistance = areaSpecificResistance / area;
            }
        }
        
        ThermalResistanceElement r;
        r.nodeFromId = nodeIdx;
        r.quadrantFrom = ThermalNodeFace::NONE;
        r.nodeToId = coldPlateIdx;
        r.quadrantTo = ThermalNodeFace::NONE;
        r.type = HeatTransferType::CONDUCTION;
        r.resistance = timResistance;
        r.area = _nodes[nodeIdx].getTotalSurfaceArea();
        _resistances.push_back(r);
    }
}

} // namespace OpenMagnetics
