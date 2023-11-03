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

    static double calculate_losses_per_meter(WireWrapper wire, SignalDescriptor current, double temperature)
    {
        return WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature);
    }

    static double calculate_effective_resistance_per_meter(WireWrapper wire, double effectiveFrequency, double temperature)
    {
        return WindingOhmicLosses::calculate_effective_resistance_per_meter(wire, effectiveFrequency, temperature);
    }


    static double calculate_skin_effect_resistance_per_meter(WireWrapper wire, SignalDescriptor current, double temperature)
    {
        if (!current.get_processed()->get_rms()) {
            throw std::runtime_error("Current processed is missing field RMS");
        }
        auto currentRms = current.get_processed()->get_rms().value();
        return WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature) / pow(currentRms, 2);
    }

};

} // namespace OpenMagnetics