#include "Temperature.h"

namespace OpenMagnetics {

    double Temperature::calculate_temperature_from_core_thermal_resistance(CoreWrapper core, double totalLosses){
    	if (!core.get_processed_description()) {
    		throw std::runtime_error("Core is missing processed description");
    	}
    	if (!core.get_processed_description()->get_thermal_resistance()) {
    		throw std::runtime_error("Core is missing thermal resistance");
    	}

    	return core.get_processed_description()->get_thermal_resistance().value() * totalLosses;
    }
    double Temperature::calculate_temperature_from_core_thermal_resistance(double thermalResistance, double totalLosses){
    	return thermalResistance * totalLosses;
    }
} // namespace OpenMagnetics
