#pragma once
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/WindingLosses.h"
#include "advisers/CoreAdviser.h"
#include "advisers/CoilAdviser.h"
#include "support/Utils.h"
#include "constructive_models/Mas.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class MagneticSimulator {
    private:
        bool _enableTemperatureConvergence = false;

        CoreTemperatureModels _coreTemperatureModelName;
        ReluctanceModels _reluctanceModelName;

        std::shared_ptr<CoreTemperatureModel> _coreTemperatureModel;
        MagnetizingInductance _magnetizingInductanceModel;
        CoreLosses _coreLossesModel;

    public:

        MagneticSimulator() {
            _coreTemperatureModelName = defaults.coreTemperatureModelDefault;
            
            _reluctanceModelName = defaults.reluctanceModelDefault;

            _coreTemperatureModel = CoreTemperatureModel::factory(_coreTemperatureModelName);
            _magnetizingInductanceModel = MagnetizingInductance(std::string(magic_enum::enum_name(_reluctanceModelName)));
        }

        void set_core_temperature_model_name(CoreTemperatureModels model) {
            _coreTemperatureModelName = model;
        }
        void set_reluctance_model_name(ReluctanceModels model) {
            _reluctanceModelName = model;
        }
        void set_core_losses_model_name(CoreLossesModels model) {
            _coreLossesModel.set_core_losses_model_name(model);
        }

        Mas simulate(Mas mas, bool fastMode=false);
        Mas simulate(const Inputs& inputs, const Magnetic& magnetic, bool fastMode=false);
        CoreLossesOutput calculate_core_losses(OperatingPoint& operatingPoint, Magnetic magnetic);
        LeakageInductanceOutput calculate_leakage_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        MagnetizingInductanceOutput calculate_magnetizing_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        WindingLossesOutput calculate_winding_losses(OperatingPoint& operatingPoint, Magnetic magnetic, std::optional<double> temperature = std::nullopt);

};


} // namespace OpenMagnetics