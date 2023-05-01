#if !defined(DEFAULTS_H)
#    define DEFAULTS_H 1
#    pragma once

#    include "CoreLosses.h"
#    include "CoreTemperature.h"
#    include "Reluctance.h"

namespace OpenMagnetics {
struct Defaults {
    Defaults() {};
    const OpenMagnetics::CoreLossesModels coreLossesModelDefault = OpenMagnetics::CoreLossesModels::IGSE;
    const OpenMagnetics::CoreTemperatureModels coreTemperatureModelDefault =
        OpenMagnetics::CoreTemperatureModels::MANIKTALA;
    const OpenMagnetics::ReluctanceModels reluctanceModelDefault = OpenMagnetics::ReluctanceModels::ZHANG;
    const double maximumProportionMagneticFluxDensitySaturation = 0.7;
    const double coreAdviserReferenceFrequency = 100000;
    const double coreAdviserReferenceMagneticFluxDensityProportion = 0.5;
};
} // namespace OpenMagnetics

#endif