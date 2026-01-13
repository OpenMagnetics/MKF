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

double CoreThermalResistanceManiktalaModel::get_core_thermal_resistance_reluctance(Core core) {
    double effectiveVolume = core.get_effective_volume();
    return 53 * pow(effectiveVolume, -0.54);
}

} // namespace OpenMagnetics
