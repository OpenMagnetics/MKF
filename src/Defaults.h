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
};
} // namespace OpenMagnetics

#endif