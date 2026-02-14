#include "support/Utils.h"
#include "physical_models/ThermalResistance.h"
#include "support/Exceptions.h"
#include <magic_enum.hpp>
#include <algorithm>
#include <cmath>

namespace OpenMagnetics {
using namespace ThermalDefaults; // IMP-2


std::shared_ptr<CoreThermalResistanceModel> CoreThermalResistanceModel::factory(CoreThermalResistanceModels modelName) {
    if (modelName == CoreThermalResistanceModels::MANIKTALA) {
        return std::make_shared<CoreThermalResistanceManiktalaModel>();
    }
    else
        throw ModelNotAvailableException("Unknown Thermal Resistance mode, available options are: {MANIKTALA}");
}

std::shared_ptr<CoreThermalResistanceModel> CoreThermalResistanceModel::factory() {
    return factory(defaults.coreThermalResistanceModelDefault);
}

/**
 * @brief Calculates core thermal resistance using Maniktala's empirical formula.
 */
double CoreThermalResistanceManiktalaModel::get_core_thermal_resistance_reluctance(Core core) {
    double effectiveVolume = core.get_effective_volume();
    return 53 * pow(effectiveVolume, -0.54);
}

//=============================================================================
// Static Thermal Resistance Calculation Methods
//=============================================================================

double ThermalResistance::calculateNaturalConvectionCoefficient(
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
    Constants constants;
    double Gr = constants.gravityAcceleration * air.thermalExpansionCoefficient * deltaT * 
                std::pow(characteristicLength, 3) / std::pow(air.kinematicViscosity, 2);
    
    // Rayleigh number: Ra = Gr * Pr
    double Ra = Gr * air.prandtlNumber;
    
    double Nu = 0.0;  // Nusselt number
    
    switch (orientation) {
        case SurfaceOrientation::VERTICAL: {
            // Churchill-Chu correlation for vertical surfaces (all Ra)
            double term = std::pow(1.0 + std::pow(0.492 / air.prandtlNumber, 9.0/16.0), 8.0/27.0);
            Nu = std::pow(0.825 + 0.387 * std::pow(Ra, 1.0/6.0) / term, 2);
            break;
        }
        case SurfaceOrientation::HORIZONTAL_TOP: {
            // Hot surface facing up (or cold facing down)
            if (Ra < 1e7) {
                Nu = 0.54 * std::pow(Ra, 0.25);
            } else {
                Nu = 0.15 * std::pow(Ra, 1.0/3.0);
            }
            break;
        }
        case SurfaceOrientation::HORIZONTAL_BOTTOM: {
            // Hot surface facing down (or cold facing up)
            Nu = 0.27 * std::pow(Ra, 0.25);
            break;
        }
    }
    
    // Ensure minimum Nusselt number
    Nu = std::max(Nu, 0.5);
    
    // Heat transfer coefficient: h = Nu * k / L
    double h = Nu * air.thermalConductivity / characteristicLength;
    
    // IMP-NEW-03: Lowered minimum from 5.0 to 2.0 W/(m2 K)
    // WHY: For small dT or small components, real h can be < 5.
    //      Clamping at 5.0 overestimates cooling -> under-predicts temperatures.
    // PREVIOUS: return std::max(h, 5.0);
    return std::max(h, 2.0); // Lowered minimum practical value
}

double ThermalResistance::calculateForcedConvectionCoefficient(
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
    return std::max(h, kConvection_MinForcedH); // IMP-7
}

double ThermalResistance::calculateRadiationCoefficient(
    double surfaceTemperature,
    double ambientTemperature,
    double emissivity) {
    
    // Convert to Kelvin
    Constants constants;
    double Ts = surfaceTemperature + constants.kelvinOffset;
    double Ta = ambientTemperature + constants.kelvinOffset;
    
    // Linearized radiation coefficient: h_rad = ε * σ * (Ts² + Ta²) * (Ts + Ta)
    double h_rad = emissivity * constants.stefanBoltzmannConstant * (Ts*Ts + Ta*Ta) * (Ts + Ta);
    
    return h_rad;
}

// IMP-NEW-08: View-factor-aware radiation coefficient
double ThermalResistance::calculateRadiationCoefficientWithViewFactor(
    double surfaceTemperature, double ambientTemperature,
    double emissivity, double viewFactor) {
    double h_rad = calculateRadiationCoefficient(surfaceTemperature, ambientTemperature, emissivity);
    return h_rad * std::clamp(viewFactor, 0.0, 1.0);
}

double ThermalResistance::getMaterialThermalConductivity(const std::string& materialName) {
    // Convert to lowercase for database lookup
    std::string lowerName = materialName;
    std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);
    
    // Try wire material database (copper, aluminium have thermalConductivity as array)
    try {
        auto wireMaterial = find_wire_material_by_name(lowerName);
        // Use the interpolation function at room temperature (25°C)
        return getWireMaterialThermalConductivity(wireMaterial, 25.0);
    } catch (const std::exception& /* e */) { // IMP-8
        // Not a wire material, continue
    }
    
    // Try insulation material database
    try {
        auto insulationMaterial = find_insulation_material_by_name(lowerName);
        auto thermalCond = insulationMaterial.get_thermal_conductivity();
        if (thermalCond) {
            return *thermalCond;
        }
    } catch (const std::exception& /* e */) { // IMP-8
        // Not an insulation material, continue
    }
    
    // Try core material database (uses heatConductivity field)
    try {
        auto coreMaterial = find_core_material_by_name(materialName);  // Core materials may be case-sensitive
        auto heatCond = coreMaterial.get_heat_conductivity();
        if (heatCond) {
            auto nominal = heatCond->get_nominal();
            if (nominal) {
                return *nominal;
            }
        }
    } catch (const std::exception& /* e */) { // IMP-8
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
    
    // Default fallback (generic material)
    return 1.0;
}

double ThermalResistance::getWireMaterialThermalConductivity(const WireMaterial& wireMaterial, double temperature) {
    // Get thermal conductivity from MAS database with temperature interpolation
    auto thermalConductivityOpt = wireMaterial.get_thermal_conductivity();
    
    if (!thermalConductivityOpt) {
        // Fallback to material name lookup  
        std::string materialName = wireMaterial.get_name();
        return getMaterialThermalConductivity(materialName);
    }
    
    auto& thermalCond = *thermalConductivityOpt;
    
    // Check if it's a constant value (just one entry)
    if (thermalCond.size() == 1) {
        return thermalCond[0].get_value();
    }
    
    // Interpolate based on temperature
    // thermalCond is array of {temperature, value} pairs
    if (thermalCond.size() == 0) {
        return getMaterialThermalConductivity("copper");  // Reasonable default
    }
    
    // Find bracketing temperatures
    for (size_t i = 0; i < thermalCond.size() - 1; ++i) {
        double T1 = thermalCond[i].get_temperature();
        double T2 = thermalCond[i+1].get_temperature();
        double k1 = thermalCond[i].get_value();
        double k2 = thermalCond[i+1].get_value();
        
        if (temperature >= T1 && temperature <= T2) {
            // Linear interpolation
            double alpha = (temperature - T1) / (T2 - T1);
            return k1 + alpha * (k2 - k1);
        }
    }
    
    // If temperature is out of range, use nearest endpoint
    if (temperature < thermalCond[0].get_temperature()) {
        return thermalCond[0].get_value();
    } else {
        return thermalCond.back().get_value();
    }
}

double ThermalResistance::getCoreMaterialThermalConductivity(const CoreMaterial& coreMaterial) {
    auto heatCond = coreMaterial.get_heat_conductivity();
    if (heatCond) {
        auto nominal = heatCond->get_nominal();
        if (nominal) {
            return *nominal;
        }
    }
    
    // Fallback for ferrite
    return 4.0;
}

double ThermalResistance::getInsulationMaterialThermalConductivity(const InsulationMaterial& material) {
    auto thermalCond = material.get_thermal_conductivity();
    if (thermalCond) {
        return *thermalCond;
    }
    
    // Default for insulation
    return 0.2;
}

} // namespace OpenMagnetics

