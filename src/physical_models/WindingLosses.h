#pragma once
#include "MAS.hpp"
#include "physical_models/MagneticField.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/WindingProximityEffectLosses.h"

using namespace MAS; // QUAL-001 TODO: qualify types and remove

namespace OpenMagnetics {

// Struct to hold all optional model choices for winding losses calculation
struct WindingLossesModels {
    std::optional<MagneticFieldStrengthModels> magneticFieldStrengthModel = std::nullopt;
    std::optional<MagneticFieldStrengthFringingEffectModels> magneticFieldStrengthFringingEffectModel = std::nullopt;
    std::optional<WindingSkinEffectLossesModels> skinEffectModel = std::nullopt;
    std::optional<WindingProximityEffectLossesModels> proximityEffectModel = std::nullopt;
};

class WindingLosses
{
    private:
        int64_t _quickModeForManyTurnsThreshold = 1000;
        WindingLossesModels _models;

    public:

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

        WindingLossesOutput calculate_losses(Magnetic magnetic, OperatingPoint operatingPoint, double temperature);
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