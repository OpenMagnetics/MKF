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
        // Full thermal-network hot-spot for one simulated operating point, configured through
        // TemperatureConfig::fromSimulatedOutput (ABT #906) — the same config the temperature-map
        // plot wrappers use, so the exported MAS and the UI can never disagree.
        TemperatureOutput calculate_temperature(OperatingPoint& operatingPoint, Magnetic magnetic, const Outputs& output);
        // ABT #838: `knownWindingLosses`, when given, is the winding-loss result of the pass just
        // run for this operating point. The core-loss loop needs the CORE's temperature, and the
        // winding heats the core: without it the loop solves a core-only network with core losses
        // alone, so on a winding-dominated design the core sits cooler than it really does and its
        // losses are read off the material curve at the wrong temperature. With it, the loop uses
        // the full network (turn nodes carry the winding power) exactly as the exported
        // temperature does. Absent = the historical core-only behaviour, which is also the first
        // pass, when no winding-loss result exists yet.
        CoreLossesOutput calculate_core_losses(OperatingPoint& operatingPoint, Magnetic magnetic,
                                               std::optional<WindingLossesOutput> knownWindingLosses = std::nullopt);
        LeakageInductanceOutput calculate_leakage_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        static LeakageInductanceOutput calculate_leakage_inductance(Magnetic magnetic, double frequency);
        MagnetizingInductanceOutput calculate_magnetizing_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        WindingLossesOutput calculate_winding_losses(OperatingPoint& operatingPoint, Magnetic magnetic, std::optional<double> temperature = std::nullopt);

};


} // namespace OpenMagnetics