#pragma once
#include "MAS.hpp"
#include "physical_models/MagneticField.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/WindingProximityEffectLosses.h"

using namespace MAS; // QUAL-001 TODO: qualify types and remove

namespace OpenMagnetics {

// Struct to hold all optional model choices for winding losses calculation
struct WindingLossesModelSet;  // physical_models/WindingLossesModelSets.h (ABT #1304)

struct WindingLossesModels {
    std::optional<MagneticFieldStrengthModels> magneticFieldStrengthModel = std::nullopt;
    std::optional<MagneticFieldStrengthFringingEffectModels> magneticFieldStrengthFringingEffectModel = std::nullopt;
    std::optional<WindingSkinEffectLossesModels> skinEffectModel = std::nullopt;
    std::optional<WindingProximityEffectLossesModels> proximityEffectModel = std::nullopt;
    // ABT #1304: a named set fixes every model above; combining it with an individual model throws. nullopt: the
    // set in Settings (default CUSTOM, the individually selected models).
    std::optional<WindingLossesModelSets> modelSet = std::nullopt;
};

class WindingLosses
{
    private:
        int64_t _quickModeForManyTurnsThreshold = 1000;
        WindingLossesModels _models;
        bool use_quick_harmonic_threshold(const Magnetic& magnetic) const;
        // The model set the losses of `magnetic` are computed with (nullptr = CUSTOM); throws when
        // a named set is combined with individually chosen models or does not model the magnetic.
        const WindingLossesModelSet* resolve_model_set(Magnetic& magnetic) const;

    public:
        // The field step of the winding losses (ABT #838): the operating point the losses are
        // computed for (normalised as calculate_losses normalises it), the magnetic field strength
        // H at every turn, and what the field was solved with: the winding-loss model set
        // (nullptr = CUSTOM, the individually selected models; ABT #1304 / #1759) and the split of
        // the harmonic currents between a winding's parallels (ABT #1314; empty without parallels).
        // Without parallels H follows from the currents and the geometry and does not depend on
        // the temperature: a caller that evaluates the same coil and operating point at several
        // temperatures (the loss <-> temperature solve) computes it once and runs only the
        // per-turn kernels at each temperature. With parallels the split is solved at the
        // conductors' temperature and feeds H, so the field belongs to that temperature
        // (depends_on_temperature) and calculate_losses refuses it at another one.
        struct WindingLossesField {
            OperatingPoint operatingPoint;
            WindingWindowMagneticStrengthFieldPhasorOutput field;
            std::optional<ParallelCurrentSplit> parallelCurrentSplit;
            const WindingLossesModelSet* modelSet;
            double temperature;
            bool depends_on_temperature() const { return parallelCurrentSplit.has_value(); }
        };

        WindingLosses() {}
        
        WindingLosses(WindingLossesModels models) : _models(models) {}

        // Legacy constructor for backwards compatibility
        WindingLosses(std::optional<MagneticFieldStrengthModels> magneticFieldStrengthModel,
                      MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModel = Defaults().magneticFieldStrengthFringingEffectModelDefault) {
            _models.magneticFieldStrengthModel = magneticFieldStrengthModel;
            _models.magneticFieldStrengthFringingEffectModel = magneticFieldStrengthFringingEffectModel;
        }
        
        WindingLossesModels get_models() const { return _models; }
        void set_models(WindingLossesModels models) { _models = models; }

        // Field step, then the per-turn kernels at `temperature`.
        WindingLossesOutput calculate_losses(Magnetic magnetic, OperatingPoint operatingPoint, double temperature);
        // `temperature`: the conductors' temperature the parallel current split is solved at.
        WindingLossesField calculate_magnetic_field(Magnetic magnetic, OperatingPoint operatingPoint, double temperature);
        // The per-turn kernels (ohmic, skin effect, proximity from the field step's H) at
        // `temperature`. `field` must come from calculate_magnetic_field of this magnetic, and
        // from that same temperature when it depends on it.
        WindingLossesOutput calculate_losses(Magnetic magnetic, const WindingLossesField& field, double temperature);
        static WindingLossesOutput combine_turn_losses(WindingLossesOutput windingLossesOutput, Coil coil);
        static double get_total_winding_losses(std::vector<WindingLossesPerElement> windingLossesPerElement);
        static double get_total_winding_losses(WindingLossesPerElement windingLossesThisElement);
        static double calculate_effective_resistance_of_winding(Magnetic magnetic, size_t windingIndex, double frequency, double temperature);
        ScalarMatrixAtFrequency calculate_resistance_matrix(Magnetic magnetic, double temperature, double frequency);
        static double calculate_losses_per_meter(Wire wire, SignalDescriptor current, double temperature);
        static double calculate_effective_resistance_per_meter(const Wire& wire, double effectiveFrequency, double temperature);
        static double calculate_skin_effect_resistance_per_meter(Wire wire, SignalDescriptor current, double temperature);

};

} // namespace OpenMagnetics