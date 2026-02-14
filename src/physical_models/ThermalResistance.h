#pragma once
#include "Constants.h"
#include "constructive_models/Core.h"
#include "ThermalNode.h"  // For HeatTransferType
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

// SurfaceOrientation is now defined in ThermalNode.h

/**
 * @brief Air and fluid properties for convection calculations
 */
struct FluidProperties {
    double density;                  // kg/m³
    double thermalConductivity;      // W/(m·K)
    double dynamicViscosity;         // Pa·s
    double kinematicViscosity;       // m²/s
    double prandtlNumber;
    double thermalExpansionCoefficient;  // 1/K
    
    /**
     * @brief Get air properties at a given temperature
     * @param temperature Temperature in °C
     * @return FluidProperties for air
     */
    static FluidProperties getAirProperties(double temperature) {
        FluidProperties air;
        double T = temperature + 273.15;  // Convert to Kelvin
        
        // Ideal gas approximation for density
        air.density = 101325.0 / (287.05 * T);
        
        // Sutherland's law for viscosity
        double mu0 = 1.716e-5;
        double T0 = 273.15;
        double S = 110.4;
        air.dynamicViscosity = mu0 * std::pow(T / T0, 1.5) * (T0 + S) / (T + S);
        
        air.kinematicViscosity = air.dynamicViscosity / air.density;
        
        // Thermal conductivity (linear approximation)
        air.thermalConductivity = 0.0241 + 7.5e-5 * temperature;
        
        // IMP-NEW-02: Temperature-dependent Prandtl number
        // Pr varies from ~0.715 (0C) to ~0.700 (200C) per Incropera Table A.4
        // PREVIOUS: air.prandtlNumber = 0.71;
        air.prandtlNumber = 0.7150 - 7.5e-5 * temperature;
        air.prandtlNumber = std::clamp(air.prandtlNumber, 0.680, 0.720);
        
        // Thermal expansion coefficient (ideal gas)
        air.thermalExpansionCoefficient = 1.0 / T;
        
        return air;
    }
};

class ThermalResistance {
private:
protected:
public:
    ThermalResistance() = default;
    ~ThermalResistance() = default;
    
    // Static utility methods for thermal resistance calculations
    
    /**
     * @brief Calculate conduction thermal resistance
     * R = L / (k * A)
     */
    static double calculateConductionResistance(double length, double thermalConductivity, double area) {
        if (length <= 0) return 0.0;
        if (thermalConductivity <= 0 || area <= 0) return 1e9;  // Return high resistance for invalid params
        return length / (thermalConductivity * area);
    }
    
    /**
     * @brief Calculate natural convection heat transfer coefficient
     * Uses Churchill-Chu correlation for vertical surfaces
     * Uses McAdams correlation for horizontal surfaces
     */
    static double calculateNaturalConvectionCoefficient(
        double surfaceTemperature,
        double ambientTemperature,
        double characteristicLength,
        SurfaceOrientation orientation);
    
    /**
     * @brief Calculate forced convection heat transfer coefficient
     */
    static double calculateForcedConvectionCoefficient(
        double airVelocity,
        double characteristicLength,
        double temperature);
    
    /**
     * @brief Calculate convection thermal resistance
     * R = 1 / (h * A)
     */
    static double calculateConvectionResistance(double heatTransferCoefficient, double area) {
        if (heatTransferCoefficient <= 0 || area <= 0) {
            throw std::runtime_error("Invalid parameters for convection resistance");
        }
        return 1.0 / (heatTransferCoefficient * area);
    }
    
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
        double area) {
        double h_rad = calculateRadiationCoefficient(surfaceTemperature, ambientTemperature, emissivity);
        return 1.0 / (h_rad * area);
    }
    
    /**
     * @brief Get thermal conductivity for common materials
     */
    static double getMaterialThermalConductivity(const std::string& materialName);
    
    /**
     * @brief Get thermal conductivity for wire material with temperature interpolation
     */
    static double getWireMaterialThermalConductivity(const WireMaterial& wireMaterial, double temperature);
    
    /**
     * @brief Get thermal conductivity from core material
     */
    static double getCoreMaterialThermalConductivity(const CoreMaterial& coreMaterial);
    
    /**
     * @brief Get thermal conductivity from insulation material
     */
    static double getInsulationMaterialThermalConductivity(const InsulationMaterial& material);
    
    // IMP-NEW-07: Proper characteristic length for convection
    static double calculateCharacteristicLength(
        double surfaceArea, double height, double width, SurfaceOrientation orientation) {
        switch (orientation) {
            case SurfaceOrientation::VERTICAL: return std::max(height, 1e-4);
            case SurfaceOrientation::HORIZONTAL_TOP:
            case SurfaceOrientation::HORIZONTAL_BOTTOM: {
                if (width > 0 && surfaceArea > 0) {
                    double depth = surfaceArea / std::max(width, 1e-6);
                    return std::max(surfaceArea / (2.0*(width+depth)), 1e-4);
                }
                return std::max(std::sqrt(std::max(surfaceArea,1e-12))/4.0, 1e-4);
            }
            default: return std::max(height, 1e-4);
        }
    }
    
    // IMP-NEW-08: View-factor-aware radiation coefficient
    static double calculateRadiationCoefficientWithViewFactor(
        double surfaceTemperature, double ambientTemperature,
        double emissivity, double viewFactor);
};

class CoreThermalResistanceModel {
private:
protected:
public:
    CoreThermalResistanceModel() = default;
    ~CoreThermalResistanceModel() = default;
    virtual double get_core_thermal_resistance_reluctance(Core core) = 0;
    static std::shared_ptr<CoreThermalResistanceModel> factory(CoreThermalResistanceModels modelName);
    static std::shared_ptr<CoreThermalResistanceModel> factory();
};


// Based on Maniktala fomula for thermal resistance
// https://www.e-magnetica.pl/doku.php/thermal_resistance_of_ferrite_cores
class CoreThermalResistanceManiktalaModel : public CoreThermalResistanceModel {
  public:
    std::string methodName = "Maniktala";
    double get_core_thermal_resistance_reluctance(Core core);
};


} // namespace OpenMagnetics