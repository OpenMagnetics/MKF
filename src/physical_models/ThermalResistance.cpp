#include "support/Utils.h"
#include "physical_models/ThermalResistance.h"
#include "support/Exceptions.h"

namespace OpenMagnetics {

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
 * 
 * Based on: "Switching Power Supplies A-Z" by Sanjaya Maniktala
 * 2nd Edition, Newnes/Elsevier, 2012, ISBN 978-0-12-386533-5
 * Chapter 3: Off-Line Converter Design and Magnetics, page 154
 * 
 * This empirical formula provides the thermal resistance for EE, EI, ETD, and EC
 * type ferrite cores based on their effective volume:
 * 
 *   Rth = 53 × Ve^(-0.54)  [°C/W]
 * 
 * where:
 * - Ve = effective core volume [cm³]
 * - Rth = thermal resistance from core surface to ambient [°C/W]
 * 
 * Example from book: For ETD-34 with Ve = 7.64 cm³:
 *   Rth = 53 × 7.64^(-0.54) = 17.67 °C/W
 * 
 * This formula assumes natural convection cooling and is suitable for
 * estimating temperature rise: ΔT = Rth × (Pcore + Pcu)
 * 
 * @param core Core object containing effective volume
 * @return Thermal resistance [°C/W]
 */
double CoreThermalResistanceManiktalaModel::get_core_thermal_resistance_reluctance(Core core) {
    double effectiveVolume = core.get_effective_volume();
    return 53 * pow(effectiveVolume, -0.54);
}

} // namespace OpenMagnetics
