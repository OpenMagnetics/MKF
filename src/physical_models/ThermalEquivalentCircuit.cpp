#include "ThermalEquivalentCircuit.h"
#include "support/Utils.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <set>

namespace OpenMagnetics {

// Physical constants
constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;  // W/(m²·K⁴)
constexpr double GRAVITY = 9.81;                      // m/s²
constexpr double KELVIN_OFFSET = 273.15;

//=============================================================================
// SimpleMatrix Implementation
//=============================================================================

std::vector<double> SimpleMatrix::solve(SimpleMatrix A, std::vector<double> b) {
    // Gauss-Jordan elimination with partial pivoting
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

//=============================================================================
// FluidProperties Implementation
//=============================================================================

FluidProperties FluidProperties::getAirProperties(double temperature) {
    // Air properties at 1 atm, interpolated for temperature range 0-200°C
    // Reference: Engineering Toolbox, various heat transfer textbooks
    FluidProperties props;
    
    double T = temperature + KELVIN_OFFSET;  // Convert to Kelvin for calculations
    
    // Density (ideal gas approximation)
    props.density = 101325.0 / (287.0 * T);  // ρ = P/(R*T)
    
    // Thermal conductivity (linear approximation 20-100°C: k ≈ 0.024 - 0.032 W/(m·K))
    props.thermalConductivity = 0.0241 + 7.7e-5 * temperature;
    
    // Dynamic viscosity (Sutherland's formula approximation)
    props.dynamicViscosity = 1.716e-5 * std::pow(T / 273.0, 1.5) * (273.0 + 111.0) / (T + 111.0);
    
    // Specific heat (relatively constant for air)
    props.specificHeat = 1006.0;  // J/(kg·K)
    
    // Thermal expansion coefficient (for ideal gas: β = 1/T)
    props.thermalExpansionCoefficient = 1.0 / T;
    
    // Kinematic viscosity
    props.kinematicViscosity = props.dynamicViscosity / props.density;
    
    // Prandtl number (approximately constant for air ≈ 0.71)
    props.prandtlNumber = props.specificHeat * props.dynamicViscosity / props.thermalConductivity;
    
    return props;
}

//=============================================================================
// ThermalEquivalentCircuit Implementation
//=============================================================================

ThermalEquivalentCircuit::ThermalEquivalentCircuit(const ThermalModelConfiguration& config)
    : config(config), ambientNodeId(0) {
}

void ThermalEquivalentCircuit::setConfiguration(const ThermalModelConfiguration& config) {
    this->config = config;
}

double ThermalEquivalentCircuit::calculateConductionResistance(double length, double thermalConductivity, double area) {
    if (thermalConductivity <= 0 || area <= 0) {
        throw std::invalid_argument("Thermal conductivity and area must be positive");
    }
    if (length <= 0) {
        return 0.0;  // No resistance for zero length
    }
    return length / (thermalConductivity * area);
}

double ThermalEquivalentCircuit::calculateNaturalConvectionCoefficient(
    double surfaceTemperature,
    double ambientTemperature,
    double characteristicLength,
    SurfaceOrientation orientation) {
    
    // Get air properties at film temperature
    double filmTemp = (surfaceTemperature + ambientTemperature) / 2.0;
    FluidProperties air = FluidProperties::getAirProperties(filmTemp);
    
    // Temperature difference
    double deltaT = std::abs(surfaceTemperature - ambientTemperature);
    if (deltaT < 0.1) {
        deltaT = 0.1;  // Avoid division by zero, assume small gradient
    }
    
    // Grashof number: Gr = g * β * ΔT * L³ / ν²
    double Gr = GRAVITY * air.thermalExpansionCoefficient * deltaT * 
                std::pow(characteristicLength, 3) / std::pow(air.kinematicViscosity, 2);
    
    // Rayleigh number: Ra = Gr * Pr
    double Ra = Gr * air.prandtlNumber;
    
    double Nu = 0.0;  // Nusselt number
    
    switch (orientation) {
        case SurfaceOrientation::VERTICAL: {
            // Churchill-Chu correlation for vertical surfaces (all Ra)
            // Nu = {0.825 + 0.387*Ra^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27)}²
            double term = std::pow(1.0 + std::pow(0.492 / air.prandtlNumber, 9.0/16.0), 8.0/27.0);
            Nu = std::pow(0.825 + 0.387 * std::pow(Ra, 1.0/6.0) / term, 2);
            break;
        }
        case SurfaceOrientation::HORIZONTAL_TOP: {
            // Hot surface facing up (or cold facing down)
            if (Ra < 1e4) {
                Nu = 0.54 * std::pow(Ra, 0.25);
            } else if (Ra < 1e7) {
                Nu = 0.54 * std::pow(Ra, 0.25);
            } else {
                Nu = 0.15 * std::pow(Ra, 1.0/3.0);
            }
            break;
        }
        case SurfaceOrientation::HORIZONTAL_BOTTOM: {
            // Hot surface facing down (or cold facing up)
            if (Ra < 1e10) {
                Nu = 0.27 * std::pow(Ra, 0.25);
            } else {
                Nu = 0.27 * std::pow(Ra, 0.25);
            }
            break;
        }
    }
    
    // Ensure minimum Nusselt number
    Nu = std::max(Nu, 0.5);
    
    // Heat transfer coefficient: h = Nu * k / L
    double h = Nu * air.thermalConductivity / characteristicLength;
    
    // Typical range for natural convection: 5-25 W/(m²·K)
    return std::max(h, 5.0);  // Minimum practical value
}

double ThermalEquivalentCircuit::calculateForcedConvectionCoefficient(
    double airVelocity,
    double characteristicLength,
    double temperature) {
    
    if (airVelocity <= 0) {
        return calculateNaturalConvectionCoefficient(temperature + 20, temperature, 
                                                     characteristicLength, SurfaceOrientation::VERTICAL);
    }
    
    FluidProperties air = FluidProperties::getAirProperties(temperature);
    
    // Reynolds number: Re = V * L / ν
    double Re = airVelocity * characteristicLength / air.kinematicViscosity;
    
    double Nu = 0.0;
    
    if (Re < 5e5) {
        // Laminar flow: Nu = 0.664 * Re^0.5 * Pr^(1/3)
        Nu = 0.664 * std::pow(Re, 0.5) * std::pow(air.prandtlNumber, 1.0/3.0);
    } else {
        // Turbulent flow: Nu = 0.037 * Re^0.8 * Pr^(1/3)
        Nu = 0.037 * std::pow(Re, 0.8) * std::pow(air.prandtlNumber, 1.0/3.0);
    }
    
    // Heat transfer coefficient
    double h = Nu * air.thermalConductivity / characteristicLength;
    
    // Typical range for forced convection: 25-250 W/(m²·K)
    return std::max(h, 10.0);
}

double ThermalEquivalentCircuit::calculateConvectionResistance(double heatTransferCoefficient, double area) {
    if (heatTransferCoefficient <= 0 || area <= 0) {
        throw std::invalid_argument("Heat transfer coefficient and area must be positive");
    }
    return 1.0 / (heatTransferCoefficient * area);
}

double ThermalEquivalentCircuit::calculateRadiationCoefficient(
    double surfaceTemperature,
    double ambientTemperature,
    double emissivity) {
    
    // Convert to Kelvin
    double Ts = surfaceTemperature + KELVIN_OFFSET;
    double Ta = ambientTemperature + KELVIN_OFFSET;
    
    // Linearized radiation coefficient: h_rad = ε * σ * (Ts² + Ta²) * (Ts + Ta)
    double h_rad = emissivity * STEFAN_BOLTZMANN * (Ts*Ts + Ta*Ta) * (Ts + Ta);
    
    return h_rad;
}

double ThermalEquivalentCircuit::calculateRadiationResistance(
    double surfaceTemperature,
    double ambientTemperature,
    double emissivity,
    double area) {
    
    double h_rad = calculateRadiationCoefficient(surfaceTemperature, ambientTemperature, emissivity);
    
    if (h_rad <= 0 || area <= 0) {
        return 1e9;  // Very high resistance (effectively no radiation)
    }
    
    return 1.0 / (h_rad * area);
}

double ThermalEquivalentCircuit::getMaterialThermalConductivity(const std::string& materialName) {
    // First try to look up from MAS databases
    // Convert to lowercase for database lookup
    std::string lowerName = materialName;
    std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);
    
    // Try wire material database (copper, aluminium have thermalConductivity as array)
    try {
        auto wireMaterial = find_wire_material_by_name(lowerName);
        // Use the interpolation function at room temperature (25°C)
        return getWireMaterialThermalConductivity(wireMaterial, 25.0);
    } catch (...) {
        // Not a wire material, continue
    }
    
    // Try insulation material database
    try {
        auto insulationMaterial = find_insulation_material_by_name(lowerName);
        auto thermalCond = insulationMaterial.get_thermal_conductivity();
        if (thermalCond) {
            return *thermalCond;
        }
    } catch (...) {
        // Not an insulation material, continue
    }
    
    // Try core material database (uses heatConductivity field)
    try {
        auto coreMaterial = find_core_material_by_name(materialName);  // Core materials may be case-sensitive (e.g., "N87")
        auto heatCond = coreMaterial.get_heat_conductivity();
        if (heatCond) {
            // DimensionWithTolerance has nominal value
            auto nominal = heatCond->get_nominal();
            if (nominal) {
                return *nominal;
            }
        }
    } catch (...) {
        // Not a core material, continue
    }
    
    // Fallback to common thermal conductivities in W/(m·K)
    static const std::map<std::string, double> conductivities = {
        {"copper", 385.0},
        {"aluminium", 237.0},
        {"aluminum", 237.0},
        {"ferrite", 4.0},
        {"iron_powder", 20.0},
        {"air", 0.026},
        {"epoxy", 0.25},
        {"polyamide", 0.25},
        {"nylon", 0.25},
        {"pet", 0.15},
        {"pbt", 0.2},
        {"lcp", 0.3},
        {"mica", 0.5},
        {"kapton", 0.12},
        {"mylar", 0.15},
        {"teflon", 0.25},
        {"silicone", 0.2},
        {"thermal_compound", 1.0},
        {"solder", 50.0}
    };
    
    std::string name = materialName;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    auto it = conductivities.find(name);
    if (it != conductivities.end()) {
        return it->second;
    }
    
    // Default for unknown materials (assume plastic-like)
    return 0.25;
}

double ThermalEquivalentCircuit::getWireMaterialThermalConductivity(const WireMaterial& wireMaterial, double temperature) {
    auto thermalCond = wireMaterial.get_thermal_conductivity();
    if (!thermalCond || thermalCond->empty()) {
        // Default to copper if no data
        return 385.0;
    }
    
    const auto& data = *thermalCond;
    
    // Find surrounding points for interpolation
    if (data.size() == 1) {
        return data[0].get_value();
    }
    
    // Sort by temperature
    std::vector<std::pair<double, double>> points;
    for (const auto& elem : data) {
        points.emplace_back(elem.get_temperature(), elem.get_value());
    }
    std::sort(points.begin(), points.end());
    
    // Linear interpolation
    if (temperature <= points.front().first) {
        return points.front().second;
    }
    if (temperature >= points.back().first) {
        return points.back().second;
    }
    
    for (size_t i = 0; i < points.size() - 1; ++i) {
        if (temperature >= points[i].first && temperature <= points[i + 1].first) {
            double t = (temperature - points[i].first) / (points[i + 1].first - points[i].first);
            return points[i].second + t * (points[i + 1].second - points[i].second);
        }
    }
    
    return 385.0;  // Default copper
}

double ThermalEquivalentCircuit::getCoreMaterialThermalConductivity(const CoreMaterial& coreMaterial) {
    auto heatCond = coreMaterial.get_heat_conductivity();
    if (heatCond) {
        auto nominal = heatCond->get_nominal();
        if (nominal) {
            return *nominal;
        }
    }
    
    // Default based on material type
    auto materialType = coreMaterial.get_material();
    switch (materialType) {
        case MaterialType::FERRITE:
            return 4.0;     // Typical ferrite k = 3.5-5 W/(m·K)
        case MaterialType::POWDER:
            return 20.0;    // Iron powder composites
        case MaterialType::NANOCRYSTALLINE:
            return 10.0;    // Nanocrystalline
        case MaterialType::AMORPHOUS:
            return 10.0;    // Amorphous
        case MaterialType::ELECTRICAL_STEEL:
            return 30.0;    // Electrical steel
        default:
            return 4.0;     // Default ferrite
    }
}

double ThermalEquivalentCircuit::getInsulationMaterialThermalConductivity(const InsulationMaterial& material) {
    auto thermalCond = material.get_thermal_conductivity();
    if (thermalCond) {
        return *thermalCond;
    }
    
    // Default for insulation materials (typical plastics)
    return 0.2;  // W/(m·K)
}

double ThermalEquivalentCircuit::calculateGapThermalResistance(const CoreGap& gap) {
    double gapLength = gap.get_length();
    if (gapLength <= 0) {
        return 0.0;  // No resistance for zero-length gap
    }
    
    // Gap area - use provided area or estimate from section dimensions
    double gapArea;
    auto areaOpt = gap.get_area();
    if (areaOpt) {
        gapArea = *areaOpt;
    } else {
        auto sectionDims = gap.get_section_dimensions();
        if (sectionDims && sectionDims->size() >= 2) {
            // Assuming rectangular/round section
            auto shape = gap.get_shape();
            if (shape && *shape == ColumnShape::ROUND) {
                double d = (*sectionDims)[0];
                gapArea = M_PI * d * d / 4.0;
            } else {
                gapArea = (*sectionDims)[0] * (*sectionDims)[1];
            }
        } else {
            // Default small gap area (1 cm²)
            gapArea = 1e-4;
        }
    }
    
    // Air thermal conductivity at room temperature
    // Could be enhanced with spacer material properties for additive gaps
    double gapK = 0.026;  // W/(m·K) for air
    
    auto gapType = gap.get_type();
    if (gapType == GapType::ADDITIVE) {
        // Additive gaps may have a spacer material
        // Could look up spacer material here if available
        gapK = 0.2;  // Typical plastic spacer
    }
    
    // Thermal resistance: R = L / (k * A)
    return gapLength / (gapK * gapArea);
}

void ThermalEquivalentCircuit::createCoreNodes(Core core) {
    // Get core geometry
    auto geoDesc = core.create_geometrical_description();
    if (!geoDesc) {
        throw std::runtime_error("Failed to create core geometrical description");
    }
    
    auto columns = core.get_columns();
    
    // Create node for each column
    for (size_t i = 0; i < columns.size(); ++i) {
        ThermalNode node;
        node.id = nodes.size();
        node.type = (columns[i].get_type() == ColumnType::CENTRAL) ? 
                    ThermalNodeType::CORE_CENTRAL_COLUMN : ThermalNodeType::CORE_LATERAL_COLUMN;
        node.name = "Core_Column_" + std::to_string(i);
        node.temperature = config.ambientTemperature;
        node.powerDissipation = 0.0;
        
        // Get column dimensions
        auto colShape = columns[i].get_shape();
        double height = columns[i].get_height();
        double depth = columns[i].get_depth();
        double width = columns[i].get_width();
        
        node.volume = height * depth * width;  // Approximate
        
        // Exposed surface area (vertical surfaces of column)
        node.exposedSurfaceArea = 2 * (height * depth + height * width);
        
        node.coordinates = columns[i].get_coordinates();
        node.emissivity = config.defaultEmissivity;
        
        nodes.push_back(node);
    }
    
    // Create nodes for top and bottom yokes (simplified as single nodes)
    double coreHeight = core.get_height();
    double coreDepth = core.get_depth();
    double coreWidth = core.get_width();
    
    // Top yoke
    ThermalNode topYoke;
    topYoke.id = nodes.size();
    topYoke.type = ThermalNodeType::CORE_TOP_YOKE;
    topYoke.name = "Core_Top_Yoke";
    topYoke.temperature = config.ambientTemperature;
    topYoke.powerDissipation = 0.0;
    topYoke.volume = coreWidth * coreDepth * (coreHeight * 0.2);  // Estimate yoke as 20% of height
    topYoke.exposedSurfaceArea = coreWidth * coreDepth;  // Top surface
    topYoke.coordinates = {0, coreHeight / 2.0, 0};
    topYoke.emissivity = config.defaultEmissivity;
    nodes.push_back(topYoke);
    
    // Bottom yoke
    ThermalNode bottomYoke;
    bottomYoke.id = nodes.size();
    bottomYoke.type = ThermalNodeType::CORE_BOTTOM_YOKE;
    bottomYoke.name = "Core_Bottom_Yoke";
    bottomYoke.temperature = config.ambientTemperature;
    bottomYoke.powerDissipation = 0.0;
    bottomYoke.volume = coreWidth * coreDepth * (coreHeight * 0.2);
    bottomYoke.exposedSurfaceArea = coreWidth * coreDepth;  // Bottom surface
    bottomYoke.coordinates = {0, -coreHeight / 2.0, 0};
    bottomYoke.emissivity = config.defaultEmissivity;
    nodes.push_back(bottomYoke);
}

void ThermalEquivalentCircuit::createCoilNodes(const Coil& coil, const Core& core) {
    auto sections = coil.get_sections_description();
    auto layers = coil.get_layers_description();
    auto turns = coil.get_turns_description();
    
    // If per-turn nodes requested and turns are available, use createCoilTurnNodes
    if (config.nodePerCoilTurn && turns && !turns->empty()) {
        // Note: createCoilTurnNodes takes Coil by copy, so we need to pass a copy
        Coil coilCopy = coil;
        createCoilTurnNodes(coilCopy, core);
        return;
    }
    
    // Try layer-based nodes first if configured and layers exist
    if (config.nodePerCoilLayer && layers && !layers->empty()) {
        // Create one node per conduction layer (skip insulation layers)
        size_t conductionLayerIndex = 0;
        for (size_t l = 0; l < layers->size(); ++l) {
            auto& layer = (*layers)[l];
            
            // Skip insulation layers - they don't generate heat
            if (layer.get_type() == ElectricalType::INSULATION) {
                continue;
            }
            
            ThermalNode node;
            node.id = nodes.size();
            node.type = ThermalNodeType::COIL_LAYER;
            node.name = layer.get_name();
            node.temperature = config.ambientTemperature;
            node.powerDissipation = 0.0;  // Will be set by distributeLayerLosses
            
            auto layerDims = layer.get_dimensions();
            auto layerCoords = layer.get_coordinates();
            
            // Estimate volume and surface area
            if (layerDims.size() >= 2 && layerCoords.size() >= 2) {
                // Cylindrical layer approximation
                double width = layerDims[0];
                double height = layerDims[1];
                double radius = layerCoords[0] + width / 2;
                node.volume = M_PI * width * radius * 2 * height;
                node.exposedSurfaceArea = 2 * M_PI * radius * height;
                node.coordinates = {layerCoords[0], layerCoords[1], 0};
            } else {
                node.volume = 1e-6;  // 1 mm³ default
                node.exposedSurfaceArea = 1e-4;  // 1 cm² default
                node.coordinates = {0, 0, 0};
            }
            
            node.emissivity = config.defaultEmissivity;
            nodes.push_back(node);
            conductionLayerIndex++;
        }
    } else if (sections && !sections->empty()) {
        // Fallback to section-based nodes if no layers available
        for (size_t s = 0; s < sections->size(); ++s) {
            auto& section = (*sections)[s];
            auto dims = section.get_dimensions();
            auto coords = section.get_coordinates();
            
            ThermalNode node;
            node.id = nodes.size();
            node.type = ThermalNodeType::COIL_SECTION;
            node.name = section.get_name();
            node.temperature = config.ambientTemperature;
            node.powerDissipation = 0.0;  // Will be set by distributeSectionLosses
            
            if (dims.size() >= 2 && coords.size() >= 2) {
                double width = dims[0];
                double height = dims[1];
                double radius = coords[0] + width / 2;
                node.volume = M_PI * width * radius * 2 * height;
                node.exposedSurfaceArea = 2 * M_PI * (coords[0] + width) * height;
                node.coordinates = {coords[0], coords[1], 0};
            } else {
                node.volume = 1e-6;
                node.exposedSurfaceArea = 1e-4;
                node.coordinates = {0, 0, 0};
            }
            
            node.emissivity = config.defaultEmissivity;
            nodes.push_back(node);
        }
    } else {
        // Last resort: create a generic coil node based on winding window dimensions
        auto windings = coil.get_functional_description();
        if (!windings.empty()) {
            ThermalNode node;
            node.id = nodes.size();
            node.type = ThermalNodeType::COIL_SECTION;
            node.name = "Coil";
            node.temperature = config.ambientTemperature;
            node.powerDissipation = 0.0;
            
            // Estimate coil dimensions from core's winding window
            auto processedCore = core.get_processed_description();
            double coilVolume = 1e-6;  // Default 1 mm³
            double coilSurfaceArea = 1e-4;  // Default 1 cm²
            
            if (processedCore) {
                auto windingWindows = processedCore->get_winding_windows();
                if (!windingWindows.empty()) {
                    auto& window = windingWindows[0];
                    double windowWidth = 0.010;   // 10mm default
                    double windowHeight = 0.020;  // 20mm default
                    
                    // Get actual dimensions - WindingWindowElement has get_width(), get_height(), get_radial_height()
                    auto radialHeight = window.get_radial_height();
                    auto height = window.get_height();
                    if (radialHeight) windowWidth = *radialHeight;
                    if (height) windowHeight = *height;
                    
                    // Estimate fill factor around 50%
                    double fillFactor = 0.5;
                    
                    // For ETD/E type cores, coil wraps around central column
                    auto columns = processedCore->get_columns();
                    double centralRadius = 0.005;  // 5mm default
                    if (!columns.empty()) {
                        auto& col = columns[0];
                        // ColumnElement has get_width() directly
                        centralRadius = col.get_width() / 2.0;
                    }
                    
                    // Coil is an annular cylinder
                    double innerRadius = centralRadius;
                    double outerRadius = centralRadius + windowWidth * fillFactor;
                    double coilHeight = windowHeight * 0.9;  // 90% fill
                    
                    coilVolume = M_PI * (outerRadius * outerRadius - innerRadius * innerRadius) * coilHeight;
                    // Surface area: inner + outer cylindrical surfaces
                    coilSurfaceArea = 2 * M_PI * outerRadius * coilHeight + 
                                     2 * M_PI * innerRadius * coilHeight;
                }
            }
            
            node.volume = coilVolume;
            node.exposedSurfaceArea = coilSurfaceArea;
            node.coordinates = {0, 0, 0};
            node.emissivity = config.defaultEmissivity;
            nodes.push_back(node);
        }
    }
}

void ThermalEquivalentCircuit::createBobbinNodes(Bobbin bobbin) {
    // Create inner and outer bobbin surface nodes
    auto processed = bobbin.get_processed_description();
    if (!processed) return;
    
    double wallThickness = processed->get_wall_thickness();
    
    // Get winding window dimensions [width, height]
    auto windingWindowDims = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoords = bobbin.get_winding_window_coordinates(0);
    
    double width = windingWindowDims.size() >= 1 ? windingWindowDims[0] : 0.01;
    double height = windingWindowDims.size() >= 2 ? windingWindowDims[1] : 0.01;
    double radius = windingWindowCoords.size() >= 1 ? windingWindowCoords[0] : 0.01;
    
    // Inner bobbin surface (in contact with windings)
    ThermalNode innerNode;
    innerNode.id = nodes.size();
    innerNode.type = ThermalNodeType::BOBBIN_INNER;
    innerNode.name = "Bobbin_Inner";
    innerNode.temperature = config.ambientTemperature;
    innerNode.powerDissipation = 0.0;
    innerNode.volume = 2 * M_PI * radius * height * wallThickness;
    innerNode.exposedSurfaceArea = 2 * M_PI * radius * height;
    innerNode.coordinates = {radius, 0, 0};
    innerNode.emissivity = config.defaultEmissivity;
    nodes.push_back(innerNode);
    
    // Outer bobbin surface (in contact with core or air)
    ThermalNode outerNode;
    outerNode.id = nodes.size();
    outerNode.type = ThermalNodeType::BOBBIN_OUTER;
    outerNode.name = "Bobbin_Outer";
    outerNode.temperature = config.ambientTemperature;
    outerNode.powerDissipation = 0.0;
    double outerRadius = radius + width;
    outerNode.volume = 2 * M_PI * outerRadius * height * wallThickness;
    outerNode.exposedSurfaceArea = 2 * M_PI * outerRadius * height;
    outerNode.coordinates = {outerRadius, 0, 0};
    outerNode.emissivity = config.defaultEmissivity;
    nodes.push_back(outerNode);
}

void ThermalEquivalentCircuit::createCoilTurnNodes(Coil coil, const Core& core) {
    auto turns = coil.get_turns_description();
    if (!turns || turns->empty()) {
        return;
    }
    
    for (size_t t = 0; t < turns->size(); ++t) {
        const auto& turn = (*turns)[t];
        
        ThermalNode node;
        node.id = nodes.size();
        node.type = ThermalNodeType::COIL_TURN;
        node.name = turn.get_name();
        node.temperature = config.ambientTemperature;
        node.powerDissipation = 0.0;  // Will be set by distributeTurnLosses
        
        // Get turn dimensions and coordinates
        auto coords = turn.get_coordinates();
        auto dims = turn.get_dimensions();
        double turnLength = turn.get_length();
        
        // Estimate turn cross-sectional area
        double turnArea = 0;
        if (dims && dims->size() >= 2) {
            // Rectangular or defined dimensions
            turnArea = (*dims)[0] * (*dims)[1];
        } else if (dims && dims->size() >= 1) {
            // Round wire - diameter given
            double d = (*dims)[0];
            turnArea = M_PI * d * d / 4.0;
        } else {
            // Fallback: try to get from wire
            turnArea = 1e-6;  // 1 mm² default
        }
        
        node.volume = turnArea * turnLength;
        node.exposedSurfaceArea = std::sqrt(turnArea * M_PI) * turnLength;  // Approximate circumference * length
        
        if (coords.size() >= 2) {
            node.coordinates = {coords[0], coords[1], coords.size() >= 3 ? coords[2] : 0};
        } else {
            node.coordinates = {0, 0, 0};
        }
        
        node.emissivity = config.defaultEmissivity;
        nodes.push_back(node);
    }
}

double ThermalEquivalentCircuit::getTotalLossFromElement(const WindingLossesPerElement& element) {
    double total = 0.0;
    
    // Add ohmic losses
    if (element.get_ohmic_losses()) {
        total += element.get_ohmic_losses()->get_losses();
    }
    
    // Add skin effect losses (sum of harmonics)
    if (element.get_skin_effect_losses()) {
        auto skinLosses = element.get_skin_effect_losses()->get_losses_per_harmonic();
        for (double loss : skinLosses) {
            total += loss;
        }
    }
    
    // Add proximity effect losses (sum of harmonics)
    if (element.get_proximity_effect_losses()) {
        auto proxLosses = element.get_proximity_effect_losses()->get_losses_per_harmonic();
        for (double loss : proxLosses) {
            total += loss;
        }
    }
    
    return total;
}

void ThermalEquivalentCircuit::distributeTurnLosses(const Coil& coil, const WindingLossesOutput& windingLosses) {
    auto lossesPerTurnOpt = windingLosses.get_winding_losses_per_turn();
    if (!lossesPerTurnOpt) {
        throw std::runtime_error("Per-turn losses are required but not available in WindingLossesOutput");
    }
    
    auto& lossesPerTurn = *lossesPerTurnOpt;
    
    // Match turn nodes to their losses by index
    size_t turnIndex = 0;
    for (auto& node : nodes) {
        if (node.type == ThermalNodeType::COIL_TURN) {
            if (turnIndex < lossesPerTurn.size()) {
                node.powerDissipation = getTotalLossFromElement(lossesPerTurn[turnIndex]);
            }
            turnIndex++;
        }
    }
}

double ThermalEquivalentCircuit::calculateNodeDistance(const ThermalNode& node1, const ThermalNode& node2) {
    double dist = 0;
    size_t n = std::min(node1.coordinates.size(), node2.coordinates.size());
    for (size_t i = 0; i < n; ++i) {
        dist += std::pow(node1.coordinates[i] - node2.coordinates[i], 2);
    }
    return std::sqrt(dist);
}

bool ThermalEquivalentCircuit::isTurnAdjacentToCore(const Turn& turn, Core core) const {
    // A turn is adjacent to core if it's in the innermost layer (closest to central column)
    // or if it's at the top/bottom of the winding window (close to yokes)
    auto coords = turn.get_coordinates();
    if (coords.empty()) return false;
    
    // Get core column radius
    auto columns = core.get_columns();
    double centralColumnRadius = 0.005;  // Default 5mm
    for (const auto& col : columns) {
        if (col.get_type() == ColumnType::CENTRAL) {
            centralColumnRadius = col.get_width() / 2.0;
            break;
        }
    }
    
    // Check radial proximity to central column
    double turnRadius = coords[0];  // x-coordinate is radial distance
    double gapToCore = turnRadius - centralColumnRadius;
    
    // Consider adjacent if within 2mm of core surface
    return gapToCore < 0.002;
}

bool ThermalEquivalentCircuit::isTurnOnOuterLayer(const Turn& turn, const Coil& coil) const {
    // Check if turn is on the outermost layer (furthest from core)
    auto turns = coil.get_turns_description();
    if (!turns || turns->empty()) return false;
    
    auto coords = turn.get_coordinates();
    if (coords.empty()) return false;
    
    double turnRadius = coords[0];
    
    // Find maximum radius among all turns
    double maxRadius = 0;
    for (const auto& t : *turns) {
        auto tc = t.get_coordinates();
        if (!tc.empty()) {
            maxRadius = std::max(maxRadius, tc[0]);
        }
    }
    
    // Consider on outer layer if within 1mm of max radius
    return (maxRadius - turnRadius) < 0.001;
}

bool ThermalEquivalentCircuit::isTurnAdjacentToBobbin(const Turn& turn, const Coil& coil) const {
    // A turn is adjacent to bobbin if it's in the innermost layer
    // (bobbin is between core and winding)
    auto turns = coil.get_turns_description();
    if (!turns || turns->empty()) return false;
    
    auto coords = turn.get_coordinates();
    if (coords.empty()) return false;
    
    double turnRadius = coords[0];
    
    // Find minimum radius among all turns
    double minRadius = std::numeric_limits<double>::max();
    for (const auto& t : *turns) {
        auto tc = t.get_coordinates();
        if (!tc.empty()) {
            minRadius = std::min(minRadius, tc[0]);
        }
    }
    
    // Consider adjacent to bobbin if within 0.5mm of min radius
    return (turnRadius - minRadius) < 0.0005;
}

std::vector<size_t> ThermalEquivalentCircuit::findNeighboringTurnNodes(size_t turnNodeId, const Coil& coil) {
    std::vector<size_t> neighbors;
    
    const ThermalNode& turnNode = nodes[turnNodeId];
    
    // Find other turn nodes that are close (within 3x turn diameter)
    double searchRadius = 0.003;  // 3mm default
    if (!turnNode.coordinates.empty()) {
        // Estimate turn size from volume
        double turnDiameter = std::pow(turnNode.volume, 1.0/3.0) * 2;
        searchRadius = turnDiameter * 3;
    }
    
    for (const auto& node : nodes) {
        if (node.id == turnNodeId) continue;
        if (node.type != ThermalNodeType::COIL_TURN) continue;
        
        double dist = calculateNodeDistance(turnNode, node);
        if (dist < searchRadius && dist > 0) {
            neighbors.push_back(node.id);
        }
    }
    
    return neighbors;
}

void ThermalEquivalentCircuit::createTurnThermalResistances(Magnetic magnetic) {
    Core core = magnetic.get_core();
    Coil coil = magnetic.get_coil();
    auto turns = coil.get_turns_description();
    
    if (!turns || turns->empty()) {
        return;
    }
    
    // Find node IDs by type
    std::vector<size_t> turnNodes;
    std::vector<size_t> coreColumnNodes;
    size_t topYokeNode = 0, bottomYokeNode = 0;
    size_t bobbinInnerNode = 0;
    
    for (const auto& node : nodes) {
        switch (node.type) {
            case ThermalNodeType::COIL_TURN:
                turnNodes.push_back(node.id);
                break;
            case ThermalNodeType::CORE_CENTRAL_COLUMN:
            case ThermalNodeType::CORE_LATERAL_COLUMN:
                coreColumnNodes.push_back(node.id);
                break;
            case ThermalNodeType::CORE_TOP_YOKE:
                topYokeNode = node.id;
                break;
            case ThermalNodeType::CORE_BOTTOM_YOKE:
                bottomYokeNode = node.id;
                break;
            case ThermalNodeType::BOBBIN_INNER:
                bobbinInnerNode = node.id;
                break;
            default:
                break;
        }
    }
    
    // Get thermal conductivities
    double windingK = config.windingEffectiveThermalConductivity;
    // double bobbinK = config.bobbinThermalConductivity;  // For future use
    
    // Create resistances between neighboring turns (turn-to-turn conduction)
    std::set<std::pair<size_t, size_t>> connectedPairs;
    
    for (size_t i = 0; i < turnNodes.size(); ++i) {
        size_t turnNodeId = turnNodes[i];
        
        // Find the corresponding Turn by matching node name
        const Turn* turnPtr = nullptr;
        for (const auto& t : *turns) {
            if (t.get_name() == nodes[turnNodeId].name) {
                turnPtr = &t;
                break;
            }
        }
        if (!turnPtr) continue;  // Skip if we can't find the turn
        const Turn& turn = *turnPtr;
        
        // Connect to neighboring turns
        auto neighbors = findNeighboringTurnNodes(turnNodeId, coil);
        for (size_t neighborId : neighbors) {
            // Avoid duplicate connections
            auto pairKey = std::make_pair(std::min(turnNodeId, neighborId), std::max(turnNodeId, neighborId));
            if (connectedPairs.count(pairKey) > 0) continue;
            connectedPairs.insert(pairKey);
            
            ThermalResistanceElement res;
            res.nodeFromId = turnNodeId;
            res.nodeToId = neighborId;
            res.transferType = HeatTransferType::CONDUCTION;
            
            double dist = calculateNodeDistance(nodes[turnNodeId], nodes[neighborId]);
            res.length = dist;
            res.thermalConductivity = windingK;
            
            // Contact area - approximate as turn cross-section
            double area = std::min(nodes[turnNodeId].volume, nodes[neighborId].volume) / 
                         std::max(res.length, 0.0001);
            res.area = area;
            
            res.resistance = calculateConductionResistance(res.length, res.thermalConductivity, res.area);
            if (res.resistance > 0 && res.resistance < 1e6) {
                resistances.push_back(res);
            }
        }
        
        // Connect turns adjacent to bobbin
        if (isTurnAdjacentToBobbin(turn, coil) && bobbinInnerNode > 0) {
            ThermalResistanceElement res;
            res.nodeFromId = turnNodeId;
            res.nodeToId = bobbinInnerNode;
            res.transferType = HeatTransferType::CONDUCTION;
            res.length = 0.0005;  // 0.5mm insulation gap
            res.thermalConductivity = 0.2;  // Insulation/air
            res.area = nodes[turnNodeId].exposedSurfaceArea * 0.3;
            res.resistance = calculateConductionResistance(res.length, res.thermalConductivity, res.area);
            if (res.resistance > 0) {
                resistances.push_back(res);
            }
        }
        
        // Connect turns on outer layer to ambient (convection)
        if (isTurnOnOuterLayer(turn, coil)) {
            ThermalResistanceElement convRes;
            convRes.nodeFromId = turnNodeId;
            convRes.nodeToId = ambientNodeId;
            convRes.transferType = config.includeForcedConvection ? 
                                   HeatTransferType::CONVECTION_FORCED : HeatTransferType::CONVECTION_NATURAL;
            convRes.area = nodes[turnNodeId].exposedSurfaceArea * 0.5;  // Partial exposure
            convRes.orientation = SurfaceOrientation::VERTICAL;
            
            double charLength = std::pow(nodes[turnNodeId].volume, 1.0/3.0);
            double h = config.includeForcedConvection ?
                      calculateForcedConvectionCoefficient(config.airVelocity, charLength, config.ambientTemperature) :
                      calculateNaturalConvectionCoefficient(nodes[turnNodeId].temperature, config.ambientTemperature,
                                                           charLength, convRes.orientation);
            convRes.resistance = calculateConvectionResistance(h, convRes.area);
            resistances.push_back(convRes);
        }
        
        // Connect turns adjacent to core
        if (isTurnAdjacentToCore(turn, core) && !coreColumnNodes.empty()) {
            // Connect to nearest core column
            size_t nearestCoreNode = coreColumnNodes[0];
            double minDist = std::numeric_limits<double>::max();
            for (size_t colId : coreColumnNodes) {
                double dist = calculateNodeDistance(nodes[turnNodeId], nodes[colId]);
                if (dist < minDist) {
                    minDist = dist;
                    nearestCoreNode = colId;
                }
            }
            
            ThermalResistanceElement res;
            res.nodeFromId = turnNodeId;
            res.nodeToId = nearestCoreNode;
            res.transferType = HeatTransferType::CONDUCTION;
            res.length = minDist;
            res.thermalConductivity = 0.05;  // Air/insulation between turn and core
            res.area = nodes[turnNodeId].exposedSurfaceArea * 0.2;
            res.resistance = calculateConductionResistance(res.length, res.thermalConductivity, res.area);
            if (res.resistance > 0 && res.resistance < 1e6) {
                resistances.push_back(res);
            }
        }
    }
}

void ThermalEquivalentCircuit::createThermalResistances(Magnetic magnetic) {
    Core core = magnetic.get_core();
    Coil coil = magnetic.get_coil();
    
    // Check if we have turn nodes - if so, use turn-based resistances
    bool hasTurnNodes = false;
    for (const auto& node : nodes) {
        if (node.type == ThermalNodeType::COIL_TURN) {
            hasTurnNodes = true;
            break;
        }
    }
    
    if (hasTurnNodes) {
        // Use turn-based thermal network
        createTurnThermalResistances(magnetic);
        
        // Still need to create core resistances
        // (createTurnThermalResistances handles turn-to-core connections)
    }
    
    // Find node IDs by type
    std::vector<size_t> coreColumnNodes;
    std::vector<size_t> coilNodes;
    size_t topYokeNode = 0, bottomYokeNode = 0;
    size_t bobbinInnerNode = 0, bobbinOuterNode = 0;
    
    for (const auto& node : nodes) {
        switch (node.type) {
            case ThermalNodeType::CORE_CENTRAL_COLUMN:
            case ThermalNodeType::CORE_LATERAL_COLUMN:
                coreColumnNodes.push_back(node.id);
                break;
            case ThermalNodeType::CORE_TOP_YOKE:
                topYokeNode = node.id;
                break;
            case ThermalNodeType::CORE_BOTTOM_YOKE:
                bottomYokeNode = node.id;
                break;
            case ThermalNodeType::COIL_SECTION:
            case ThermalNodeType::COIL_LAYER:
                coilNodes.push_back(node.id);
                break;
            case ThermalNodeType::BOBBIN_INNER:
                bobbinInnerNode = node.id;
                break;
            case ThermalNodeType::BOBBIN_OUTER:
                bobbinOuterNode = node.id;
                break;
            default:
                break;
        }
    }
    
    // Get core thermal conductivity from MAS material data
    double coreK = config.coreThermalConductivity;  // Default fallback
    try {
        auto coreMaterial = find_core_material_by_name(core.get_material_name());
        coreK = getCoreMaterialThermalConductivity(coreMaterial);
    } catch (...) {
        // Use config default if material lookup fails
    }
    
    double bobbinK = config.bobbinThermalConductivity;
    
    // Add gap thermal resistances
    auto gapping = core.get_gapping();
    double totalGapResistance = 0.0;
    for (const auto& gap : gapping) {
        double gapR = calculateGapThermalResistance(gap);
        totalGapResistance += gapR;
    }
    
    // Distribute gap resistance among central column nodes (gaps are typically in central column)
    double gapResistancePerColumn = (coreColumnNodes.empty() || totalGapResistance == 0) ? 
                                     0 : totalGapResistance / coreColumnNodes.size();
    
    // Core internal conduction: columns to yokes (including gap resistance in series)
    for (size_t colId : coreColumnNodes) {
        // Column to top yoke
        ThermalResistanceElement res;
        res.nodeFromId = colId;
        res.nodeToId = topYokeNode;
        res.transferType = HeatTransferType::CONDUCTION;
        
        // Estimate path length and area
        double colHeight = nodes[colId].coordinates.size() > 1 ? 
                          std::abs(nodes[colId].coordinates[1] - nodes[topYokeNode].coordinates[1]) : 0.01;
        double area = std::pow(nodes[colId].volume / colHeight, 0.5) * colHeight * 0.5;  // Rough estimate
        
        res.length = colHeight / 2;
        res.thermalConductivity = coreK;
        res.area = area;
        
        // Core conduction resistance plus gap resistance in series
        double coreCondR = calculateConductionResistance(res.length, res.thermalConductivity, res.area);
        
        // Add gap resistance if this is a central column (gaps typically in center leg)
        bool isCentralColumn = (nodes[colId].type == ThermalNodeType::CORE_CENTRAL_COLUMN);
        if (isCentralColumn && gapResistancePerColumn > 0) {
            res.resistance = coreCondR + gapResistancePerColumn / 2;  // Half to each yoke direction
        } else {
            res.resistance = coreCondR;
        }
        
        resistances.push_back(res);
        
        // Column to bottom yoke
        res.nodeToId = bottomYokeNode;
        if (isCentralColumn && gapResistancePerColumn > 0) {
            res.resistance = coreCondR + gapResistancePerColumn / 2;
        }
        resistances.push_back(res);
    }
    
    // Core surface convection/radiation to ambient
    for (const auto& node : nodes) {
        if (node.type == ThermalNodeType::CORE_CENTRAL_COLUMN ||
            node.type == ThermalNodeType::CORE_LATERAL_COLUMN ||
            node.type == ThermalNodeType::CORE_TOP_YOKE ||
            node.type == ThermalNodeType::CORE_BOTTOM_YOKE) {
            
            // Natural convection
            ThermalResistanceElement convRes;
            convRes.nodeFromId = node.id;
            convRes.nodeToId = ambientNodeId;
            convRes.transferType = config.includeForcedConvection ? 
                                   HeatTransferType::CONVECTION_FORCED : HeatTransferType::CONVECTION_NATURAL;
            convRes.area = node.exposedSurfaceArea;
            
            // Determine orientation
            if (node.type == ThermalNodeType::CORE_TOP_YOKE) {
                convRes.orientation = SurfaceOrientation::HORIZONTAL_TOP;
            } else if (node.type == ThermalNodeType::CORE_BOTTOM_YOKE) {
                convRes.orientation = SurfaceOrientation::HORIZONTAL_BOTTOM;
            } else {
                convRes.orientation = SurfaceOrientation::VERTICAL;
            }
            
            // Initial estimate of h
            double charLength = std::pow(node.volume, 1.0/3.0);
            double h = config.includeForcedConvection ?
                      calculateForcedConvectionCoefficient(config.airVelocity, charLength, config.ambientTemperature) :
                      calculateNaturalConvectionCoefficient(node.temperature, config.ambientTemperature, 
                                                           charLength, convRes.orientation);
            convRes.resistance = calculateConvectionResistance(h, convRes.area);
            resistances.push_back(convRes);
            
            // Radiation
            if (config.includeRadiation) {
                ThermalResistanceElement radRes;
                radRes.nodeFromId = node.id;
                radRes.nodeToId = ambientNodeId;
                radRes.transferType = HeatTransferType::RADIATION;
                radRes.area = node.exposedSurfaceArea;
                radRes.resistance = calculateRadiationResistance(node.temperature, config.ambientTemperature,
                                                                 node.emissivity, radRes.area);
                resistances.push_back(radRes);
            }
        }
    }
    
    // Bobbin conduction (inner to outer)
    if (bobbinInnerNode > 0 && bobbinOuterNode > 0) {
        ThermalResistanceElement bobbinRes;
        bobbinRes.nodeFromId = bobbinInnerNode;
        bobbinRes.nodeToId = bobbinOuterNode;
        bobbinRes.transferType = HeatTransferType::CONDUCTION;
        
        // Use wall thickness as length
        Bobbin bobbin = magnetic.get_bobbin();
        auto processed = bobbin.get_processed_description();
        if (processed) {
            bobbinRes.length = processed->get_wall_thickness();
        } else {
            bobbinRes.length = 0.001;  // Default 1mm
        }
        
        bobbinRes.thermalConductivity = bobbinK;
        bobbinRes.area = nodes[bobbinInnerNode].exposedSurfaceArea;
        bobbinRes.resistance = calculateConductionResistance(bobbinRes.length, bobbinRes.thermalConductivity, 
                                                             bobbinRes.area);
        resistances.push_back(bobbinRes);
    }
    
    // Get layer descriptions for insulation data
    auto layersOpt = coil.get_layers_description();
    std::vector<Layer> layers;
    if (layersOpt) {
        layers = *layersOpt;
    }
    
    // Find insulation layers (they sit between conduction layers)
    std::vector<InsulationMaterial> interLayerInsulation;
    for (const auto& layer : layers) {
        if (layer.get_type() == ElectricalType::INSULATION) {
            try {
                auto insulationMaterial = coil.resolve_insulation_layer_insulation_material(layer);
                interLayerInsulation.push_back(insulationMaterial);
            } catch (...) {
                // Use default if resolution fails
                try {
                    auto defaultMaterial = find_insulation_material_by_name("mylar");
                    interLayerInsulation.push_back(defaultMaterial);
                } catch (...) {
                    // No insulation material data available
                }
            }
        }
    }
    
    // Coil to bobbin inner surface (conduction through first insulation layer)
    double innerInsulationK = 0.2;  // Default
    double innerInsulationThickness = 0.0005;  // Default 0.5mm
    
    // Try to get the innermost insulation layer properties
    for (const auto& layer : layers) {
        if (layer.get_type() == ElectricalType::INSULATION) {
            // First insulation layer is typically between bobbin and first winding layer
            auto dims = layer.get_dimensions();
            if (!dims.empty() && dims[0] > 0) {
                innerInsulationThickness = dims[0];  // Radial thickness
            }
            try {
                auto material = coil.resolve_insulation_layer_insulation_material(layer);
                innerInsulationK = getInsulationMaterialThermalConductivity(material);
            } catch (...) {
                // Use default
            }
            break;  // Only need first insulation layer
        }
    }
    
    for (size_t coilId : coilNodes) {
        if (bobbinInnerNode > 0) {
            // Connect coil to bobbin inner surface through insulation
            ThermalResistanceElement coilBobbinRes;
            coilBobbinRes.nodeFromId = coilId;
            coilBobbinRes.nodeToId = bobbinInnerNode;
            coilBobbinRes.transferType = HeatTransferType::CONDUCTION;
            coilBobbinRes.length = innerInsulationThickness;
            coilBobbinRes.thermalConductivity = innerInsulationK;
            coilBobbinRes.area = nodes[coilId].exposedSurfaceArea * 0.5;  // Partial contact
            coilBobbinRes.resistance = calculateConductionResistance(coilBobbinRes.length, 
                                                                     coilBobbinRes.thermalConductivity,
                                                                     coilBobbinRes.area);
            resistances.push_back(coilBobbinRes);
        } else {
            // No bobbin - connect coil directly to core columns through air/insulation
            for (size_t colId : coreColumnNodes) {
                ThermalResistanceElement coilCoreRes;
                coilCoreRes.nodeFromId = coilId;
                coilCoreRes.nodeToId = colId;
                coilCoreRes.transferType = HeatTransferType::CONDUCTION;
                coilCoreRes.length = 0.002;  // 2mm gap estimate
                coilCoreRes.thermalConductivity = 0.05;  // Air/insulation composite
                coilCoreRes.area = nodes[coilId].exposedSurfaceArea * 0.3;
                coilCoreRes.resistance = calculateConductionResistance(coilCoreRes.length, 
                                                                       coilCoreRes.thermalConductivity,
                                                                       coilCoreRes.area);
                resistances.push_back(coilCoreRes);
            }
        }
    }
    
    // Inter-layer conduction in coil (through insulation layers)
    size_t insulationIdx = 0;
    for (size_t i = 0; i + 1 < coilNodes.size(); ++i) {
        ThermalResistanceElement layerRes;
        layerRes.nodeFromId = coilNodes[i];
        layerRes.nodeToId = coilNodes[i + 1];
        layerRes.transferType = HeatTransferType::CONDUCTION;
        
        // Try to use actual insulation layer thickness and material
        double layerK = 0.2;  // Default insulation k
        double layerThickness = 0.0002;  // Default 0.2mm
        
        // Get insulation properties from corresponding layer
        if (insulationIdx < interLayerInsulation.size()) {
            layerK = getInsulationMaterialThermalConductivity(interLayerInsulation[insulationIdx]);
        }
        
        // Try to find insulation layer dimensions between these conduction layers
        size_t insulationLayerCount = 0;
        for (const auto& layer : layers) {
            if (layer.get_type() == ElectricalType::INSULATION) {
                if (insulationLayerCount == insulationIdx + 1) {  // +1 to skip first (bobbin) insulation
                    auto dims = layer.get_dimensions();
                    if (!dims.empty() && dims[0] > 0) {
                        layerThickness = dims[0];
                    }
                    break;
                }
                insulationLayerCount++;
            }
        }
        
        layerRes.length = layerThickness;
        layerRes.thermalConductivity = layerK;
        layerRes.area = std::min(nodes[coilNodes[i]].exposedSurfaceArea, 
                                 nodes[coilNodes[i + 1]].exposedSurfaceArea);
        layerRes.resistance = calculateConductionResistance(layerRes.length, layerRes.thermalConductivity,
                                                           layerRes.area);
        resistances.push_back(layerRes);
        insulationIdx++;
    }
    
    // Outer coil layer to ambient (if exposed)
    if (!coilNodes.empty()) {
        size_t outerCoilId = coilNodes.back();
        
        ThermalResistanceElement convRes;
        convRes.nodeFromId = outerCoilId;
        convRes.nodeToId = ambientNodeId;
        convRes.transferType = config.includeForcedConvection ?
                               HeatTransferType::CONVECTION_FORCED : HeatTransferType::CONVECTION_NATURAL;
        convRes.area = nodes[outerCoilId].exposedSurfaceArea * 0.3;  // Partial exposure
        convRes.orientation = SurfaceOrientation::VERTICAL;
        
        double charLength = std::pow(nodes[outerCoilId].volume, 1.0/3.0);
        double h = config.includeForcedConvection ?
                  calculateForcedConvectionCoefficient(config.airVelocity, charLength, config.ambientTemperature) :
                  calculateNaturalConvectionCoefficient(nodes[outerCoilId].temperature, config.ambientTemperature,
                                                       charLength, convRes.orientation);
        convRes.resistance = calculateConvectionResistance(h, convRes.area);
        resistances.push_back(convRes);
    }
}

void ThermalEquivalentCircuit::buildNetworkWithWindingLosses(Magnetic magnetic,
                                                             double coreLosses,
                                                             const WindingLossesOutput& windingLosses) {
    nodes.clear();
    resistances.clear();
    
    // Create ambient node first
    ThermalNode ambientNode;
    ambientNode.id = 0;
    ambientNode.type = ThermalNodeType::AMBIENT;
    ambientNode.name = "Ambient";
    ambientNode.temperature = config.ambientTemperature;
    ambientNode.powerDissipation = 0.0;
    ambientNode.volume = 0;
    ambientNode.exposedSurfaceArea = 0;
    ambientNode.emissivity = 0;
    nodes.push_back(ambientNode);
    ambientNodeId = 0;
    
    // Create nodes for each part
    Core core = magnetic.get_core();
    createCoreNodes(core);
    
    Coil coil = magnetic.get_coil();
    createCoilNodes(coil, core);  // This will use createCoilTurnNodes if nodePerCoilTurn is enabled
    
    // Get bobbin and create nodes
    Bobbin bobbin = magnetic.get_bobbin();
    if (bobbin.get_functional_description()) {
        createBobbinNodes(bobbin);
    }
    
    // Distribute core losses proportionally to volume
    double totalCoreVolume = 0;
    for (const auto& node : nodes) {
        if (node.type == ThermalNodeType::CORE_CENTRAL_COLUMN ||
            node.type == ThermalNodeType::CORE_LATERAL_COLUMN ||
            node.type == ThermalNodeType::CORE_TOP_YOKE ||
            node.type == ThermalNodeType::CORE_BOTTOM_YOKE) {
            totalCoreVolume += node.volume;
        }
    }
    
    for (auto& node : nodes) {
        if (node.type == ThermalNodeType::CORE_CENTRAL_COLUMN ||
            node.type == ThermalNodeType::CORE_LATERAL_COLUMN ||
            node.type == ThermalNodeType::CORE_TOP_YOKE ||
            node.type == ThermalNodeType::CORE_BOTTOM_YOKE) {
            if (totalCoreVolume > 0) {
                node.powerDissipation = coreLosses * (node.volume / totalCoreVolume);
            }
        }
    }
    
    // Distribute winding losses using per-turn losses (required)
    distributeTurnLosses(coil, windingLosses);
    
    // Create thermal resistances
    createThermalResistances(magnetic);
}

void ThermalEquivalentCircuit::assembleMatrix() {
    size_t n = nodes.size();
    
    conductanceMatrix = SimpleMatrix(n, n, 0.0);
    powerVector.assign(n, 0.0);
    temperatureVector.assign(n, 0.0);
    
    // Initialize temperatures
    for (size_t i = 0; i < n; ++i) {
        temperatureVector[i] = nodes[i].temperature;
        powerVector[i] = nodes[i].powerDissipation;
    }
    
    // Fill conductance matrix
    for (const auto& res : resistances) {
        if (res.resistance <= 0) continue;
        
        double G = 1.0 / res.resistance;
        
        size_t i = res.nodeFromId;
        size_t j = res.nodeToId;
        
        // Add to diagonal elements
        conductanceMatrix(i, i) += G;
        conductanceMatrix(j, j) += G;
        
        // Add to off-diagonal elements
        conductanceMatrix(i, j) -= G;
        conductanceMatrix(j, i) -= G;
    }
    
    // Handle ambient node (fixed temperature boundary condition)
    // Method: Modify the power vector to account for known ambient temperature
    // For each node i connected to ambient: P_i -= G_i,amb * T_amb
    // Then set the ambient equation to T_amb = config.ambientTemperature
    
    double T_amb = config.ambientTemperature;
    
    // Move known T_amb terms to the right-hand side
    for (size_t i = 0; i < n; ++i) {
        if (i != ambientNodeId) {
            powerVector[i] -= conductanceMatrix(i, ambientNodeId) * T_amb;
        }
    }
    
    // Now zero the ambient row/column and set diagonal to 1
    conductanceMatrix.setRowZero(ambientNodeId);
    conductanceMatrix.setColZero(ambientNodeId);
    conductanceMatrix(ambientNodeId, ambientNodeId) = 1.0;
    powerVector[ambientNodeId] = T_amb;
}

void ThermalEquivalentCircuit::solveCircuit() {
    // Solve [G][T] = [P] using Gauss-Jordan elimination
    temperatureVector = SimpleMatrix::solve(conductanceMatrix, powerVector);
    
    // Update node temperatures
    for (size_t i = 0; i < nodes.size(); ++i) {
        nodes[i].temperature = temperatureVector[i];
    }
}

void ThermalEquivalentCircuit::updateResistances() {
    // Update temperature-dependent resistances (convection and radiation)
    for (auto& res : resistances) {
        if (res.transferType == HeatTransferType::CONVECTION_NATURAL ||
            res.transferType == HeatTransferType::CONVECTION_FORCED) {
            
            double surfaceTemp = nodes[res.nodeFromId].temperature;
            double charLength = std::pow(nodes[res.nodeFromId].volume, 1.0/3.0);
            if (charLength <= 0) charLength = 0.01;
            
            double h;
            if (res.transferType == HeatTransferType::CONVECTION_FORCED) {
                h = calculateForcedConvectionCoefficient(config.airVelocity, charLength, 
                                                        config.ambientTemperature);
            } else {
                h = calculateNaturalConvectionCoefficient(surfaceTemp, config.ambientTemperature,
                                                         charLength, res.orientation);
            }
            
            res.resistance = calculateConvectionResistance(h, res.area);
        }
        else if (res.transferType == HeatTransferType::RADIATION) {
            double surfaceTemp = nodes[res.nodeFromId].temperature;
            double emissivity = nodes[res.nodeFromId].emissivity;
            res.resistance = calculateRadiationResistance(surfaceTemp, config.ambientTemperature,
                                                          emissivity, res.area);
        }
    }
}

bool ThermalEquivalentCircuit::checkConvergence(const std::vector<double>& oldTemperatures) {
    double maxDiff = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (!nodes[i].isAmbient()) {
            maxDiff = std::max(maxDiff, std::abs(temperatureVector[i] - oldTemperatures[i]));
        }
    }
    return maxDiff < config.convergenceTolerance;
}

ThermalAnalysisOutput ThermalEquivalentCircuit::calculateTemperatures(
    Magnetic magnetic,
    double coreLosses,
    const WindingLossesOutput& windingLosses) {
    
    // Always use per-turn model
    config.nodePerCoilTurn = true;
    
    ThermalAnalysisOutput output;
    output.methodUsed = "ThermalEquivalentCircuit";
    output.converged = false;
    output.iterationsToConverge = 0;
    
    // Build the thermal network using actual per-element losses
    buildNetworkWithWindingLosses(magnetic, coreLosses, windingLosses);
    
    if (nodes.size() < 2) {
        throw std::runtime_error("Insufficient nodes in thermal network");
    }
    
    // Iterative solution (for temperature-dependent properties)
    std::vector<double> oldTemperatures(nodes.size(), 0.0);
    
    for (size_t iter = 0; iter < config.maxIterations; ++iter) {
        // Store old temperatures
        for (size_t i = 0; i < nodes.size(); ++i) {
            oldTemperatures[i] = nodes[i].temperature;
        }
        
        // Assemble and solve
        assembleMatrix();
        solveCircuit();
        
        // Check convergence
        if (checkConvergence(oldTemperatures)) {
            output.converged = true;
            output.iterationsToConverge = iter + 1;
            break;
        }
        
        // Update temperature-dependent resistances
        updateResistances();
    }
    
    // Extract results
    output.maximumTemperature = config.ambientTemperature;
    double sumCoreTemp = 0, numCoreNodes = 0;
    double sumCoilTemp = 0, numCoilNodes = 0;
    
    for (const auto& node : nodes) {
        output.nodeTemperatures[node.name] = node.temperature;
        output.maximumTemperature = std::max(output.maximumTemperature, node.temperature);
        
        if (node.type == ThermalNodeType::CORE_CENTRAL_COLUMN ||
            node.type == ThermalNodeType::CORE_LATERAL_COLUMN ||
            node.type == ThermalNodeType::CORE_TOP_YOKE ||
            node.type == ThermalNodeType::CORE_BOTTOM_YOKE) {
            sumCoreTemp += node.temperature;
            numCoreNodes++;
        }
        
        if (node.type == ThermalNodeType::COIL_SECTION ||
            node.type == ThermalNodeType::COIL_LAYER ||
            node.type == ThermalNodeType::COIL_TURN) {
            sumCoilTemp += node.temperature;
            numCoilNodes++;
        }
    }
    
    output.averageCoreTemperature = numCoreNodes > 0 ? sumCoreTemp / numCoreNodes : config.ambientTemperature;
    output.averageCoilTemperature = numCoilNodes > 0 ? sumCoilTemp / numCoilNodes : config.ambientTemperature;
    
    // Calculate bulk thermal resistance
    double totalLosses = coreLosses + windingLosses.get_winding_losses();
    
    if (totalLosses > 0) {
        output.totalThermalResistance = (output.maximumTemperature - config.ambientTemperature) / totalLosses;
    } else {
        output.totalThermalResistance = 0;
    }
    
    output.thermalResistances = resistances;
    
    return output;
}

ThermalAnalysisOutput ThermalEquivalentCircuit::calculateTemperatures(
    Magnetic magnetic,
    double coreLosses,
    double windingLosses) {
    
    // Create WindingLossesOutput with losses distributed proportionally by number of turns
    WindingLossesOutput windingLossesOutput;
    windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
    windingLossesOutput.set_method_used("ProportionalDistribution");
    windingLossesOutput.set_winding_losses(windingLosses);
    
    // Ensure coil has turns description
    Coil coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        magnetic.get_mutable_coil().wind();
        coil = magnetic.get_coil();
    }
    
    auto turnsOpt = coil.get_turns_description();
    
    if (turnsOpt && !turnsOpt->empty()) {
        size_t totalTurns = turnsOpt->size();
        double lossPerTurn = windingLosses / totalTurns;
        
        std::vector<WindingLossesPerElement> perTurnLosses;
        for (size_t i = 0; i < totalTurns; ++i) {
            WindingLossesPerElement turnLoss;
            
            // Set ohmic losses for this turn
            OhmicLosses ohmic;
            ohmic.set_losses(lossPerTurn);
            ohmic.set_origin(ResultOrigin::SIMULATION);
            turnLoss.set_ohmic_losses(ohmic);
            
            perTurnLosses.push_back(turnLoss);
        }
        
        windingLossesOutput.set_winding_losses_per_turn(perTurnLosses);
    }
    
    return calculateTemperatures(magnetic, coreLosses, windingLossesOutput);
}

double ThermalEquivalentCircuit::getTemperatureAtPoint(const std::vector<double>& coordinates) {
    if (nodes.empty()) {
        return config.ambientTemperature;
    }
    
    // Find nearest node (simple nearest-neighbor interpolation)
    double minDist = std::numeric_limits<double>::max();
    double nearestTemp = config.ambientTemperature;
    
    for (const auto& node : nodes) {
        if (node.isAmbient() || node.coordinates.empty()) continue;
        
        double dist = 0;
        for (size_t i = 0; i < std::min(coordinates.size(), node.coordinates.size()); ++i) {
            dist += std::pow(coordinates[i] - node.coordinates[i], 2);
        }
        dist = std::sqrt(dist);
        
        if (dist < minDist) {
            minDist = dist;
            nearestTemp = node.temperature;
        }
    }
    
    return nearestTemp;
}

double ThermalEquivalentCircuit::getBulkThermalResistance() const {
    // Calculate total power dissipation
    double totalPower = 0;
    double maxTemp = config.ambientTemperature;
    
    for (const auto& node : nodes) {
        totalPower += node.powerDissipation;
        if (!node.isAmbient()) {
            maxTemp = std::max(maxTemp, node.temperature);
        }
    }
    
    if (totalPower > 0) {
        return (maxTemp - config.ambientTemperature) / totalPower;
    }
    return 0;
}

//=============================================================================
// ThermalModel Factory Implementation
//=============================================================================

std::shared_ptr<ThermalEquivalentCircuit> ThermalModel::factory(ModelType type) {
    switch (type) {
        case ModelType::EQUIVALENT_CIRCUIT:
        case ModelType::EQUIVALENT_CIRCUIT_SIMPLE:
            return std::make_shared<ThermalEquivalentCircuit>();
        case ModelType::BULK_MANIKTALA:
            // Fall through to default - could create specialized model here
        default:
            return std::make_shared<ThermalEquivalentCircuit>();
    }
}

} // namespace OpenMagnetics
