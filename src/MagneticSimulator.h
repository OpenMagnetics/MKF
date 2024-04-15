#pragma once
#include "MagnetizingInductance.h"
#include "CoreLosses.h"
#include "CoreTemperature.h"
#include "WindingLosses.h"
#include "CoreAdviser.h"
#include "CoilAdviser.h"
#include "MasWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class MagneticSimulator {
    private:
        bool _enableTemperatureConvergence = false;

        std::vector<CoreLossesModels> _coreLossesModelNames;
        CoreTemperatureModels _coreTemperatureModelName;
        ReluctanceModels _reluctanceModelName;

        std::vector<std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>> _coreLossesModels;
        std::shared_ptr<CoreTemperatureModel> _coreTemperatureModel;
        MagnetizingInductance _magnetizingInductanceModel;

    public:

        MagneticSimulator() {
            _coreLossesModelNames = {Defaults().coreLossesModelDefault, CoreLossesModels::PROPRIETARY, CoreLossesModels::STEINMETZ, CoreLossesModels::ROSHEN};
            _coreTemperatureModelName = Defaults().coreTemperatureModelDefault;
            
            _reluctanceModelName = Defaults().reluctanceModelDefault;

            for (auto modelName : _coreLossesModelNames) {
                _coreLossesModels.push_back(std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>{modelName, CoreLossesModel::factory(modelName)});
            }
            _coreTemperatureModel = CoreTemperatureModel::factory(_coreTemperatureModelName);
            _magnetizingInductanceModel = MagnetizingInductance(std::string(magic_enum::enum_name(_reluctanceModelName)));
        }

        void set_core_losses_model_name(CoreLossesModels model) {
            _coreLossesModelNames = {model, CoreLossesModels::PROPRIETARY, CoreLossesModels::STEINMETZ, CoreLossesModels::ROSHEN};
            _coreLossesModels.clear();
            for (auto modelName : _coreLossesModelNames) {
                _coreLossesModels.push_back(std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>{modelName, CoreLossesModel::factory(modelName)});
            }
        }
        void set_core_temperature_model_name(CoreTemperatureModels model) {
            _coreTemperatureModelName = model;
        }
        void set_reluctance_model_name(ReluctanceModels model) {
            _reluctanceModelName = model;
        }

        MasWrapper simulate(MasWrapper mas);
        MasWrapper simulate(const InputsWrapper& inputs, const MagneticWrapper& magnetic);
        CoreLossesOutput calculate_core_losses(OperatingPoint& operatingPoint, MagneticWrapper magnetic);
        MagnetizingInductanceOutput calculate_magnetizing_inductance(OperatingPoint& operatingPoint, MagneticWrapper magnetic);
        WindingLossesOutput calculate_winding_losses(OperatingPoint& operatingPoint, MagneticWrapper magnetic, std::optional<double> temperature = std::nullopt);

};


} // namespace OpenMagnetics