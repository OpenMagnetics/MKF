#pragma once
#include "MAS.hpp"
#include "physical_models/MagneticField.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/WindingProximityEffectLosses.h"

namespace OpenMagnetics {

class WindingLosses
{
    private:
        int64_t _quickModeForManyTurnsThreshold = 1000;
        MagneticFieldStrengthModels _magneticFieldStrengthModel;
        MagneticFieldStrengthFringingEffectModels _magneticFieldStrengthFringingEffectModel;

    public:

        WindingLosses(MagneticFieldStrengthModels magneticFieldStrengthModel = Defaults().magneticFieldStrengthModelDefault,
                      MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModel = Defaults().magneticFieldStrengthFringingEffectModelDefault) {
            _magneticFieldStrengthModel = magneticFieldStrengthModel;
            _magneticFieldStrengthFringingEffectModel = magneticFieldStrengthFringingEffectModel;
        }

        WindingLossesOutput calculate_losses(MagneticWrapper magnetic, OperatingPoint operatingPoint, double temperature);
        static double calculate_effective_resistance_of_winding(MagneticWrapper magnetic, size_t windingIndex, double frequency, double temperature);
        ResistanceMatrixAtFrequency calculate_resistance_matrix(MagneticWrapper magnetic, double temperature, double frequency);
        static double calculate_losses_per_meter(WireWrapper wire, SignalDescriptor current, double temperature);
        static double calculate_effective_resistance_per_meter(const WireWrapper& wire, double effectiveFrequency, double temperature);
        static double calculate_skin_effect_resistance_per_meter(WireWrapper wire, SignalDescriptor current, double temperature);

};

} // namespace OpenMagnetics