#include "Settings.h"
#include "WindingLosses.h"
#include "MAS.hpp"
#include "MagneticField.h"
#include "WindingOhmicLosses.h"
#include "WindingSkinEffectLosses.h"
#include "WindingProximityEffectLosses.h"

namespace OpenMagnetics {

WindingLossesOutput WindingLosses::calculate_losses(MagneticWrapper magnetic, OperatingPoint operatingPoint, double temperature) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto windingLossesOutput = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, temperature);
    windingLossesOutput = WindingSkinEffectLosses::calculate_skin_effect_losses(magnetic.get_coil(), temperature, windingLossesOutput, settings->get_harmonic_amplitude_threshold());

    MagneticField magneticField(_magneticFieldStrengthModel, OpenMagnetics::MagneticFieldStrengthFringingEffectModels::ROSHEN);

    int64_t totalNumberPhysicalTurns = 0;
    for (auto winding : magnetic.get_coil().get_functional_description()) {
        totalNumberPhysicalTurns += winding.get_number_turns() * winding.get_number_parallels();
    }

    auto previoustHarmonicAmplitudeThreshold = settings->get_harmonic_amplitude_threshold();
    if (settings->get_harmonic_amplitude_threshold_quick_mode() && totalNumberPhysicalTurns > _quickModeForManyTurnsThreshold) {
        settings->set_harmonic_amplitude_threshold(settings->get_harmonic_amplitude_threshold() * 2);
    }
    auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic);

    windingLossesOutput = WindingProximityEffectLosses::calculate_proximity_effect_losses(magnetic.get_coil(), temperature, windingLossesOutput, windingWindowMagneticStrengthFieldOutput);

    if (settings->get_harmonic_amplitude_threshold_quick_mode() && totalNumberPhysicalTurns > _quickModeForManyTurnsThreshold) {
        settings->set_harmonic_amplitude_threshold(previoustHarmonicAmplitudeThreshold);
    }
    return windingLossesOutput;
}

double WindingLosses::calculate_losses_per_meter(WireWrapper wire, SignalDescriptor current, double temperature)
{
    return WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature).first;
}

double WindingLosses::calculate_effective_resistance_per_meter(const WireWrapper& wire, double effectiveFrequency, double temperature)
{
    return WindingOhmicLosses::calculate_effective_resistance_per_meter(wire, effectiveFrequency, temperature);
}


double WindingLosses::calculate_skin_effect_resistance_per_meter(WireWrapper wire, SignalDescriptor current, double temperature)
{
    if (!current.get_processed()->get_rms()) {
        throw std::runtime_error("Current processed is missing field RMS");
    }
    auto currentRms = current.get_processed()->get_rms().value();
    return WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature).first / pow(currentRms, 2);
}


} // namespace OpenMagnetics