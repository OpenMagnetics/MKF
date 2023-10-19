#pragma once
#include "MAS.hpp"
#include "WindingOhmicLosses.h"
#include "WindingSkinEffectLosses.h"

namespace OpenMagnetics {

class WindingLosses
{
    public:
    static WindingLossesOutput calculate_losses(CoilWrapper coil, OperatingPoint operatingPoint, double temperature)
    {
        auto windingLossesOutput = WindingOhmicLosses::calculate_ohmic_losses(coil, operatingPoint, temperature);
        windingLossesOutput = WindingSkinEffectLosses::calculate_skin_effect_losses(coil, temperature, windingLossesOutput);

        return windingLossesOutput;
    }

};

} // namespace OpenMagnetics