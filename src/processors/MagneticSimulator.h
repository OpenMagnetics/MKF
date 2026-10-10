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

// ABT #838: which winding-loss model the loss <-> temperature solve evaluates for a wound coil.
// The advisers' loss filters each publish their own winding-loss model (MagneticFilterLosses the
// full one, MagneticFilterLossesNoProximity DC + skin effect, the core adviser's post-processing
// DC only); the solve heats the network with exactly the winding losses the caller publishes, so
// the core and winding losses it returns are one self-consistent pair and are computed once.
enum class WindingLossesScope {
    FULL,                       // WindingLosses: DC + skin + proximity (simulate())
    OHMIC_AND_SKIN_EFFECT,      // DC + skin effect, no proximity
    OHMIC                       // DC only
};

class MagneticSimulator {
    private:
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
        // Default preference: heads the core-loss model cascade (see CoreLosses::set_core_losses_model_name).
        void set_core_losses_model_name(CoreLossesModels model) {
            _coreLossesModel.set_core_losses_model_name(model);
        }
        // Explicit request: exactly this model, throwing for a material that cannot run it (ABT #1497).
        void set_core_losses_requested_model_name(CoreLossesModels model) {
            _coreLossesModel.set_core_losses_requested_model_name(model);
        }

        // The design's topology for the magnetizing-inductance model (see
        // MagnetizingInductance::set_topology). simulate() sets it from its Inputs; callers of the
        // per-operating-point methods below set it from theirs.
        void set_topology(std::optional<Topology> topology) {
            _magnetizingInductanceModel.set_topology(topology);
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
        // ABT #838: core losses at the temperature the part's dissipation produces, by
        // fixed-point iteration from ambient on the thermal network (the core-loss <-> temperature
        // loop this function has always held, now always run to convergence):
        //   - coil WOUND (turns description present): core AND winding heat the same full network,
        //     so neither temperature can be found without the other's losses. Core losses are
        //     evaluated at T_core = hottest CORE node, winding losses at T_winding = mean TURN-node
        //     temperature (WindingLosses takes one temperature for the whole coil). simulate() uses
        //     this, and so do the advisers' loss filters and the PyOM/WASM wrappers when they pass a
        //     wound coil.
        //   - no turns (the core advisers' sentinel-coil candidates): no winding losses exist to
        //     compute; the core-only network carries the core losses alone.
        // Converged when the network returns both temperatures within
        // kOperatingTemperatureTolerance of the ones the losses were evaluated at. Each
        // temperature's step is halved whenever its residual changes sign (losses falling steeply
        // with temperature make plain substitution oscillate), and the core is kept below its
        // material's Curie point, where its loss and permeability data end. Not settling within
        // kOperatingTemperatureMaximumIterations, or being driven to the Curie point, THROWS
        // (thermal runaway: there is no operating point to report, and the last iterate is not one).
        // `windingLossesOutput` / `temperatureOutput`, when given, receive the winding losses and
        // the full-network temperature solved together with the returned core losses (asking for
        // winding losses of an unwound coil throws). `windingLossesScope` selects the winding-loss
        // model the solve evaluates (and returns) for a wound coil.
        static constexpr double kOperatingTemperatureTolerance = 0.5;      // K, both temperatures
        static constexpr size_t kOperatingTemperatureMaximumIterations = 100;
        CoreLossesOutput calculate_core_losses(OperatingPoint& operatingPoint, Magnetic magnetic,
                                               std::optional<WindingLossesOutput>* windingLossesOutput = nullptr,
                                               std::optional<TemperatureOutput>* temperatureOutput = nullptr,
                                               WindingLossesScope windingLossesScope = WindingLossesScope::FULL);
        LeakageInductanceOutput calculate_leakage_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        static LeakageInductanceOutput calculate_leakage_inductance(Magnetic magnetic, double frequency);
        MagnetizingInductanceOutput calculate_magnetizing_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        WindingLossesOutput calculate_winding_losses(OperatingPoint& operatingPoint, Magnetic magnetic, std::optional<double> temperature = std::nullopt);

};


} // namespace OpenMagnetics