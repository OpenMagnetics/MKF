#pragma once
#include <MAS.hpp>

namespace OpenMagnetics {

enum class MagneticFieldStrengthModels : int {
    BINNS_LAWRENSON,
    LAMMERANER,
    DOWELL
};

enum class MagneticFieldStrengthFringingEffectModels : int {
    ROSHEN,
    ALBACH
};


enum class ReluctanceModels : int {
    ZHANG,
    PARTRIDGE,
    EFFECTIVE_AREA,
    EFFECTIVE_LENGTH,
    MUEHLETHALER,
    STENGLEIN,
    BALAKRISHNAN,
    CLASSIC
};


enum class CoreLossesModels : int {
    PROPRIETARY,
    STEINMETZ,
    IGSE,
    BARG,
    ROSHEN,
    ALBACH,
    NSE,
    MSE
};


enum class CoreTemperatureModels : int {
    KAZIMIERCZUK,
    MANIKTALA,
    TDK,
    DIXON,
    AMIDON
};


enum class WindingSkinEffectLossesModels : int {
    DOWELL,
    WOJDA,
    ALBACH,
    PAYNE,
    LOTFI,
    KAZIMIERCZUK,
    KUTKUT,
    FERREIRA,
    DIMITRAKAKIS,
    WANG,
    HOLGUIN,
    PERRY
};


enum class WindingProximityEffectLossesModels : int {
    ROSSMANITH,
    WANG,
    FERREIRA,
    LAMMERANER,
    ALBACH,
    DOWELL
};

} // namespace OpenMagnetics