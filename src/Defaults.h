#if !defined(DEFAULTS_H)
#    define DEFAULTS_H 1
#    pragma once

#include "Models.h"
#include "cstdint"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {
struct Defaults {
    Defaults() {};
    const OpenMagnetics::CoreLossesModels coreLossesModelDefault = OpenMagnetics::CoreLossesModels::IGSE;
    const OpenMagnetics::CoreTemperatureModels coreTemperatureModelDefault = OpenMagnetics::CoreTemperatureModels::MANIKTALA;
    const OpenMagnetics::ReluctanceModels reluctanceModelDefault = OpenMagnetics::ReluctanceModels::ZHANG;
    // ALBACH + ROSHEN selected based on comprehensive validation (Jan 2026):
    // Lowest overall error: 9.6% ± 8.4% across 87 tests covering ROUND, LITZ, RECTANGULAR, and FOIL wires
    // See MKF/docs/WindingLossesModelValidation.md for full methodology and results
    const OpenMagnetics::MagneticFieldStrengthModels magneticFieldStrengthModelDefault = OpenMagnetics::MagneticFieldStrengthModels::ALBACH;
    const OpenMagnetics::MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModelDefault = OpenMagnetics::MagneticFieldStrengthFringingEffectModels::ROSHEN;
    const OpenMagnetics::CoreThermalResistanceModels coreThermalResistanceModelDefault = OpenMagnetics::CoreThermalResistanceModels::MANIKTALA;
    const double maximumProportionMagneticFluxDensitySaturation = 0.7;
    const double coreAdviserFrequencyReference = 100000;
    const double coreAdviserMagneticFluxDensityReference = 0.5;
    const double coreAdviserThresholdValidity = 0.9;
    const double coreAdviserMaximumCoreTemperature = 150;
    const double coreAdviserMaximumPercentagePowerCoreLosses = 0.05;
    const uint64_t coreAdviserMaximumMagneticsAfterFiltering = 500;
    const uint64_t coreAdviserMaximumNumberStacks = 4;
    const double maximumCurrentDensity = 7000000;
    const double maximumCurrentDensityPlanar = 2000000;
    const double maximumEffectiveCurrentDensity = 12000000;
    const int maximumNumberParallels = 5;
    const double magneticFluxDensitySaturation = 0.5;
    const double ferriteInitialPermeability = 2000;
    const double magnetizingInductanceThresholdValidity = 0.25;
    const double selfResonantFrequencyMargin = 0.25;  // Maximum operating frequency as fraction of SRF (e.g., 0.25 = f_op < 0.25 * SRF)
    const double harmonicAmplitudeThreshold = 0.05;
    const double ambientTemperature = 25;
    const double measurementFrequency = 10000;
    const double maximumFrequency = 100e6;
    const int magneticFieldMirroringDimension = 1;
    const double maximumCoilPattern = 6;
    const WindingOrientation defaultRoundWindowSectionsOrientation = WindingOrientation::CONTIGUOUS;
    const WindingOrientation defaultRectangularWindowSectionsOrientation = WindingOrientation::OVERLAPPING;
    const CoilAlignment defaultRoundWindowSectionsAlignment = CoilAlignment::SPREAD;
    const CoilAlignment defaultRectangularWindowSectionsAlignment = CoilAlignment::INNER_OR_TOP;
    const std::string defaultEnamelledInsulationMaterial = "Polyurethane 155";
    const std::string defaultInsulationMaterial = "ETFE";
    const std::string defaultLayerInsulationMaterial = "Kapton HN";
    const std::string defaultConductorMaterial = "copper";
    const std::string defaultPcbInsulationMaterial = "FR4";
    const double overlappingFactorSurroundingTurns = 0.7;
    const WireStandard commonWireStandard = WireStandard::NEMA_MW_1000_C;
    const WiringTechnology wiringTechnology = WiringTechnology::WOUND;
    const double pcbInsulationThickness = 100e-6;
    const double minimumWireToWireDistance = 90e-6;
    const double minimumBorderToWireDistance = 90e-6;
    const double coreToLayerDistance = 250e-6;
};
} // namespace OpenMagnetics

#endif