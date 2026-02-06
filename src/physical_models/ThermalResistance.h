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
        
        // Prandtl number (relatively constant for air)
        air.prandtlNumber = 0.71;
        
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