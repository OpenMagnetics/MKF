#if !defined(DEFAULTS_H)
#    define DEFAULTS_H 1
#    pragma once

#include "Models.h"

namespace OpenMagnetics {
struct Defaults {
    Defaults() {};
    const OpenMagnetics::CoreLossesModels coreLossesModelDefault = OpenMagnetics::CoreLossesModels::IGSE;
    const OpenMagnetics::CoreTemperatureModels coreTemperatureModelDefault =
        OpenMagnetics::CoreTemperatureModels::MANIKTALA;
    const OpenMagnetics::ReluctanceModels reluctanceModelDefault = OpenMagnetics::ReluctanceModels::ZHANG;
    const OpenMagnetics::MagneticFieldStrengthModels magneticFieldStrengthModelDefault = OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON;
    const OpenMagnetics::MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModelDefault = OpenMagnetics::MagneticFieldStrengthFringingEffectModels::ROSHEN;
    const double maximumProportionMagneticFluxDensitySaturation = 0.7;
    const double coreAdviserFrequencyReference = 100000;
    const double coreAdviserMagneticFluxDensityReference = 0.5;
    const double coreAdviserThresholdValidity = 0.9;
    const double coreAdviserMaximumCoreTemperature = 150;
    const double coreAdviserMaximumPercentagePowerCoreLosses = 0.05;
    const uint64_t coreAdviserMaximumMagneticsAfterFiltering = 500;
    const uint64_t coreAdviserMaximumNumberStacks = 4;
    const double maximumCurrentDensity = 7000000;
    const double maximumEffectiveCurrentDensity = 12000000;
    const double maximumNumberParallels = 5;
    const double magneticFluxDensitySaturation = 0.5;
    const double magnetizingInductanceThresholdValidity = 0.25;
    const double windingLossesHarmonicAmplitudeThreshold = 0.05;
    const double ambientTemperature = 25;
    const double magneticFieldMirroringDimension = 1;
    const double maximumCoilPattern = 6;
    const WindingOrientation defaultSectionsOrientation = WindingOrientation::HORIZONTAL;
    const std::string defaultInsulationMaterial = "ETFE";
};
} // namespace OpenMagnetics

#endif