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
        std::cout << "Mierda 1" << std::endl;
        auto windingLossesOutput = WindingOhmicLosses::calculate_ohmic_losses(coil, operatingPoint, temperature);
        std::cout << "Mierda 2" << std::endl;
        windingLossesOutput = WindingSkinEffectLosses::calculate_skin_effect_losses(coil, temperature, windingLossesOutput);
        std::cout << "Mierda 3" << std::endl;

        return windingLossesOutput;
    }

};

} // namespace OpenMagnetics