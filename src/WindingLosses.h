#pragma once
#include "MAS.hpp"
#include "MagneticField.h"
#include "WindingOhmicLosses.h"
#include "WindingSkinEffectLosses.h"
#include "WindingProximityEffectLosses.h"

namespace OpenMagnetics {

class WindingLosses
{
    private:
        bool _includeFringing = true;
        bool _quickModeForManyTurns = true;
        int64_t _quickModeForManyTurnsThreshold = 1000;
        int _mirroringDimension = Defaults().magneticFieldMirroringDimension;
        double _windingLossesHarmonicAmplitudeThreshold = Defaults().windingLossesHarmonicAmplitudeThreshold;
        MagneticFieldStrengthModels _magneticFieldStrengthModel;
        MagneticFieldStrengthFringingEffectModels _magneticFieldStrengthFringingEffectModel;

    public:

        WindingLosses(MagneticFieldStrengthModels magneticFieldStrengthModel = Defaults().magneticFieldStrengthModelDefault,
                      MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModel = Defaults().magneticFieldStrengthFringingEffectModelDefault) {
            _magneticFieldStrengthModel = magneticFieldStrengthModel;
            _magneticFieldStrengthFringingEffectModel = magneticFieldStrengthFringingEffectModel;
        }

        void set_quick_mode_for_many_turns(bool value) {
            _quickModeForManyTurns = value;
        }

        void set_fringing_effect(bool value) {
            _includeFringing = value;
        }

        void set_mirroring_dimension(int mirroringDimension) {
            _mirroringDimension = mirroringDimension;
        }

        void set_winding_losses_harmonic_amplitude_threshold(double windingLossesHarmonicAmplitudeThreshold) {
            _windingLossesHarmonicAmplitudeThreshold = windingLossesHarmonicAmplitudeThreshold;
        }

        WindingLossesOutput calculate_losses(MagneticWrapper magnetic, OperatingPoint operatingPoint, double temperature);
        static double calculate_losses_per_meter(WireWrapper wire, SignalDescriptor current, double temperature);
        static double calculate_effective_resistance_per_meter(const WireWrapper& wire, double effectiveFrequency, double temperature);
        static double calculate_skin_effect_resistance_per_meter(WireWrapper wire, SignalDescriptor current, double temperature);

};

} // namespace OpenMagnetics