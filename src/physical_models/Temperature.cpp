#include "physical_models/ThermalResistance.h"
#include "physical_models/Temperature.h"
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double Temperature::calculate_temperature_from_core_thermal_resistance(Core core, double totalLosses){
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
    double Temperature::calculate_temperature_from_core_thermal_resistance(double thermalResistance, double totalLosses){
    	return thermalResistance * totalLosses;
    }
} // namespace OpenMagnetics
