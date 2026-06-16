#pragma once
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/WindingLosses.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include "constructive_models/Mas.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class MagneticSimulator {
    private:
        bool _enableTemperatureConvergence = false;

        ReluctanceModels _reluctanceModelName;

        MagnetizingInductance _magnetizingInductanceModel;
        CoreLosses _coreLossesModel;

    public:

        MagneticSimulator() {
            // Use Settings instead of defaults
            auto& settings = Settings::GetInstance();
            _reluctanceModelName = settings.get_reluctance_model();

            _magnetizingInductanceModel = MagnetizingInductance(to_string(_reluctanceModelName));
        }

        void set_reluctance_model_name(ReluctanceModels model) {
            _reluctanceModelName = model;
        }
        void set_core_losses_model_name(CoreLossesModels model) {
            _coreLossesModel.set_core_losses_model_name(model);
        }

        Mas simulate(Mas mas, bool fastMode=false);
        Mas simulate(const Inputs& inputs, const Magnetic& magnetic, bool fastMode=false);

        // Builds a manufacturer-style datasheet (MagneticManufacturerInfo.datasheetInfo) from an
        // already-simulated MAS, reading the per-operating-point Outputs and the magnetic geometry.
        // Existing manufacturerInfo fields (name, reference, description, family, ...) are preserved;
        // only the datasheetInfo block is populated/overwritten. The result is also attached to
        // mas.magnetic.manufacturerInfo. Requires mas to carry simulation outputs (run simulate first).
        MagneticManufacturerInfo build_datasheet(Mas& mas);
        CoreLossesOutput calculate_core_losses(OperatingPoint& operatingPoint, Magnetic magnetic);
        LeakageInductanceOutput calculate_leakage_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        static LeakageInductanceOutput calculate_leakage_inductance(Magnetic magnetic, double frequency);
        MagnetizingInductanceOutput calculate_magnetizing_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        WindingLossesOutput calculate_winding_losses(OperatingPoint& operatingPoint, Magnetic magnetic, std::optional<double> temperature = std::nullopt);

};


} // namespace OpenMagnetics