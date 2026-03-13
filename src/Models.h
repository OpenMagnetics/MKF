#pragma once
#include "json.hpp"

using json = nlohmann::json;

namespace OpenMagnetics {

enum class MagneticFieldStrengthModels : int {
    BINNS_LAWRENSON,
    LAMMERANER,
    DOWELL,
    WANG,
    ALBACH
};

enum class MagneticFieldStrengthFringingEffectModels : int {
    ROSHEN,
    ALBACH,
    SULLIVAN
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
    CIGSE,
    BARG,
    ROSHEN,
    ALBACH,
    NSE,
    MSE,
    LOSS_FACTOR
};


enum class CoreThermalResistanceModels : int {
    MANIKTALA,
    DIXON
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
    XI_NAN,
    KAZIMIERCZUK,
    KUTKUT,
    FERREIRA,
    DIMITRAKAKIS,
    MUEHLETHALER,
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
    DOWELL,
    XI_NAN,
    WOJDA,
    SULLIVAN,
    BARTOLI,
    VANDELAC
};


enum class StrayCapacitanceModels : int {
    KOCH,
    ALBACH,
    DUERDOTH,
    MASSARINI
};

enum class ElectricFieldOutputUnit : int {
    JOULES_PER_CUBIC_METER,  // J/m³ - Energy density (natural output of SDF_PHYSICS)
    VOLTS_PER_METER          // V/m - Electric field magnitude (requires conversion)
};

} // namespace OpenMagnetics
