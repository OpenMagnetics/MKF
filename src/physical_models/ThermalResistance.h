#pragma once
#include "Constants.h"
#include "constructive_models/Core.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class ThermalResistance {
private:
protected:
public:
    ThermalResistance() = default;
    ~ThermalResistance() = default;
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