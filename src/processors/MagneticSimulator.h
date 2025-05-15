#pragma once
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/WindingLosses.h"
#include "advisers/CoreAdviser.h"
#include "advisers/CoilAdviser.h"
#include "support/Settings.h"
#include "constructive_models/MasWrapper.h"
#include <MAS.hpp>


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
            auto settings = Settings::GetInstance();
            _coreTemperatureModelName = Defaults().coreTemperatureModelDefault;
            
            _reluctanceModelName = Defaults().reluctanceModelDefault;

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

        MasWrapper simulate(MasWrapper mas, bool fastMode=false);
        MasWrapper simulate(const InputsWrapper& inputs, const MagneticWrapper& magnetic, bool fastMode=false);
        CoreLossesOutput calculate_core_losses(OperatingPoint& operatingPoint, MagneticWrapper magnetic);
        LeakageInductanceOutput calculate_leakage_inductance(OperatingPoint& operatingPoint, MagneticWrapper magnetic);
        MagnetizingInductanceOutput calculate_magnetizing_inductance(OperatingPoint& operatingPoint, MagneticWrapper magnetic);
        WindingLossesOutput calculate_winding_losses(OperatingPoint& operatingPoint, MagneticWrapper magnetic, std::optional<double> temperature = std::nullopt);

};


} // namespace OpenMagnetics