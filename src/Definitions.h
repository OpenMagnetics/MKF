#pragma once
#include <MAS.hpp>
#include "json.hpp"
#include "Models.h"

using json = nlohmann::json;

namespace OpenMagnetics {

enum class DimensionalValues : int {
    MAXIMUM,
    NOMINAL,
    MINIMUM
};

void from_json(const json & j, DimensionalValues & x);
void to_json(json & j, const DimensionalValues & x);

inline void from_json(const json & j, DimensionalValues & x) {
    if (j == "Maximum") x = DimensionalValues::MAXIMUM;
    else if (j == "Nominal") x = DimensionalValues::NOMINAL;
    else if (j == "Minimum") x = DimensionalValues::MINIMUM;
    else { throw std::runtime_error("Input JSON does not conform to DimensionalValues schema: " + to_string(j)); }
}

inline void to_json(json & j, const DimensionalValues & x) {
    switch (x) {
        case DimensionalValues::MAXIMUM: j = "Maximum"; break;
        case DimensionalValues::NOMINAL: j = "Nominal"; break;
        case DimensionalValues::MINIMUM: j = "Minimum"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"DimensionalValues\": " + std::to_string(static_cast<int>(x)));
    }
}

enum class GappingType : int {
    GROUND,
    SPACER,
    RESIDUAL,
    DISTRIBUTED
};

void from_json(const json & j, GappingType & x);
void to_json(json & j, const GappingType & x);

inline void from_json(const json & j, GappingType & x) {
    if (j == "Ground" || j == "ground" || j == "GROUND") x = GappingType::GROUND;
    else if (j == "Spacer" || j == "spacer" || j == "SPACER") x = GappingType::SPACER;
    else if (j == "Residual" || j == "residual" || j == "RESIDUAL") x = GappingType::RESIDUAL;
    else if (j == "Distributed" || j == "distributed" || j == "DISTRIBUTED") x = GappingType::DISTRIBUTED;
    else { throw std::runtime_error("Input JSON does not conform to GappingType schema: " + to_string(j)); }
}

inline void to_json(json & j, const GappingType & x) {
    switch (x) {
        case GappingType::GROUND: j = "Ground"; break;
        case GappingType::SPACER: j = "Spacer"; break;
        case GappingType::RESIDUAL: j = "Residual"; break;
        case GappingType::DISTRIBUTED: j = "Distributed"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"GappingType\": " + std::to_string(static_cast<int>(x)));
    }
}

enum class CoreCrossReferencerFilters : int {
    PERMEANCE,
    CORE_LOSSES, 
    SATURATION, 
    WINDING_WINDOW_AREA, 
    ENVELOPING_VOLUME, 
    EFFECTIVE_AREA
};

void from_json(const json & j, CoreCrossReferencerFilters & x);
void to_json(json & j, const CoreCrossReferencerFilters & x);

inline void from_json(const json & j, CoreCrossReferencerFilters & x) {
    if (j == "Permeance" || j == "permeance" || j == "PERMEANCE") x = CoreCrossReferencerFilters::PERMEANCE;
    else if (j == "CoreLosses" || j == "corelosses" || j == "CORELOSSES") x = CoreCrossReferencerFilters::CORE_LOSSES;
    else if (j == "Saturation" || j == "saturation" || j == "SATURATION") x = CoreCrossReferencerFilters::SATURATION;
    if (j == "WindingWindowArea" || j == "windingwindowarea" || j == "WINDINGWINDOWAREA") x = CoreCrossReferencerFilters::WINDING_WINDOW_AREA;
    else if (j == "EnvelopingVolume" || j == "envelopingvolume" || j == "ENVELOPINGVOLUME") x = CoreCrossReferencerFilters::ENVELOPING_VOLUME;
    else if (j == "EffectiveArea" || j == "effectivearea" || j == "EFFECTIVEAREA") x = CoreCrossReferencerFilters::EFFECTIVE_AREA;
    else { throw std::runtime_error("Input JSON does not conform to CoreCrossReferencerFilters schema: " + to_string(j)); }
}

inline void to_json(json & j, const CoreCrossReferencerFilters & x) {
    switch (x) {
        case CoreCrossReferencerFilters::PERMEANCE: j = "Permeance"; break;
        case CoreCrossReferencerFilters::CORE_LOSSES: j = "CoreLosses"; break;
        case CoreCrossReferencerFilters::SATURATION: j = "Saturation"; break;
        case CoreCrossReferencerFilters::WINDING_WINDOW_AREA: j = "WindingWindowArea"; break;
        case CoreCrossReferencerFilters::ENVELOPING_VOLUME: j = "EnvelopingVolume"; break;
        case CoreCrossReferencerFilters::EFFECTIVE_AREA: j = "EffectiveArea"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"CoreCrossReferencerFilters\": " + std::to_string(static_cast<int>(x)));
    }
}

enum class CoreMaterialCrossReferencerFilters : int {
    INITIAL_PERMEABILITY,
    REMANENCE, 
    COERCIVE_FORCE, 
    SATURATION,
    CURIE_TEMPERATURE,
    VOLUMETRIC_LOSSES,
    RESISTIVITY
};

void from_json(const json & j, CoreMaterialCrossReferencerFilters & x);
void to_json(json & j, const CoreMaterialCrossReferencerFilters & x);

inline void from_json(const json & j, CoreMaterialCrossReferencerFilters & x) {
    if (j == "InitialPermeability" || j == "initialpermeability" || j == "INITIALPERMEABILITY") x = CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY;
    else if (j == "Remanence" || j == "remanence" || j == "REMANENCE") x = CoreMaterialCrossReferencerFilters::REMANENCE;
    else if (j == "CoerciveForce" || j == "coerciveforce" || j == "COERCIVEFORCE") x = CoreMaterialCrossReferencerFilters::COERCIVE_FORCE;
    else if (j == "Saturation" || j == "saturation" || j == "SATURATION") x = CoreMaterialCrossReferencerFilters::SATURATION;
    if (j == "CurieTemperature" || j == "curietemperature" || j == "CURIETEMPERATURE") x = CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE;
    else if (j == "VolumetricLosses" || j == "volumetriclosses" || j == "VOLUMETRICLOSSES") x = CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES;
    else if (j == "Resistivity" || j == "resistivity" || j == "RESISTIVITY") x = CoreMaterialCrossReferencerFilters::RESISTIVITY;
    else { throw std::runtime_error("Input JSON does not conform to CoreMaterialCrossReferencerFilters schema: " + to_string(j)); }
}

inline void to_json(json & j, const CoreMaterialCrossReferencerFilters & x) {
    switch (x) {
        case CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY: j = "InitialPermeability"; break;
        case CoreMaterialCrossReferencerFilters::REMANENCE: j = "Remanence"; break;
        case CoreMaterialCrossReferencerFilters::COERCIVE_FORCE: j = "CoerciveForce"; break;
        case CoreMaterialCrossReferencerFilters::SATURATION: j = "Saturation"; break;
        case CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE: j = "CurieTemperature"; break;
        case CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES: j = "VolumetricLosses"; break;
        case CoreMaterialCrossReferencerFilters::RESISTIVITY: j = "Resistivity"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"CoreMaterialCrossReferencerFilters\": " + std::to_string(static_cast<int>(x)));
    }
}

enum class OrderedIsolationSide : int { 
    PRIMARY,
    SECONDARY,
    TERTIARY,
    QUATERNARY,
    QUINARY,
    SENARY,
    SEPTENARY,
    OCTONARY,
    NONARY,
    DENARY,
    UNDENARY,
    DUODENARY
};

void from_json(const json & j, OrderedIsolationSide & x);
void to_json(json & j, const OrderedIsolationSide & x);

inline void from_json(const json & j, OrderedIsolationSide & x) {
    if (j == "primary" || j == "primary" || j == "PRIMARY") x = OrderedIsolationSide::PRIMARY;
    else if (j == "secondary" || j == "secondary" || j == "SECONDARY") x = OrderedIsolationSide::SECONDARY;
    else if (j == "tertiary" || j == "tertiary" || j == "TERTIARY") x = OrderedIsolationSide::TERTIARY;
    if (j == "quaternary" || j == "quaternary" || j == "QUATERNARY") x = OrderedIsolationSide::QUATERNARY;
    else if (j == "quinary" || j == "quinary" || j == "QUINARY") x = OrderedIsolationSide::QUINARY;
    else if (j == "senary" || j == "senary" || j == "SENARY") x = OrderedIsolationSide::SENARY;
    if (j == "septenary" || j == "septenary" || j == "SEPTENARY") x = OrderedIsolationSide::SEPTENARY;
    else if (j == "octonary" || j == "octonary" || j == "OCTONARY") x = OrderedIsolationSide::OCTONARY;
    else if (j == "nonary" || j == "nonary" || j == "NONARY") x = OrderedIsolationSide::NONARY;
    if (j == "denary" || j == "denary" || j == "DENARY") x = OrderedIsolationSide::DENARY;
    else if (j == "undenary" || j == "undenary" || j == "UNDENARY") x = OrderedIsolationSide::UNDENARY;
    else if (j == "duodenary" || j == "duodenary" || j == "DUODENARY") x = OrderedIsolationSide::DUODENARY;
    else { throw std::runtime_error("Input JSON does not conform to OrderedIsolationSide schema: " + to_string(j)); }
}

inline void to_json(json & j, const OrderedIsolationSide & x) {
    switch (x) {
        case OrderedIsolationSide::PRIMARY: j = "primary"; break;
        case OrderedIsolationSide::SECONDARY: j = "secondary"; break;
        case OrderedIsolationSide::TERTIARY: j = "tertiary"; break;
        case OrderedIsolationSide::QUATERNARY: j = "quaternary"; break;
        case OrderedIsolationSide::QUINARY: j = "quinary"; break;
        case OrderedIsolationSide::SENARY: j = "senary"; break;
        case OrderedIsolationSide::SEPTENARY: j = "septenary"; break;
        case OrderedIsolationSide::OCTONARY: j = "octonary"; break;
        case OrderedIsolationSide::NONARY: j = "nonary"; break;
        case OrderedIsolationSide::DENARY: j = "denary"; break;
        case OrderedIsolationSide::UNDENARY: j = "undenary"; break;
        case OrderedIsolationSide::DUODENARY: j = "duodenary"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"OrderedIsolationSide\": " + std::to_string(static_cast<int>(x)));
    }
}


enum class PainterModes : int {
    CONTOUR,
    QUIVER,
    SCATTER,
};

void from_json(const json & j, PainterModes & x);
void to_json(json & j, const PainterModes & x);

inline void from_json(const json & j, PainterModes & x) {
    if (j == "Contour" || j == "contour" || j == "CONTOUR") x = PainterModes::CONTOUR;
    else if (j == "Quiver" || j == "quiver" || j == "QUIVER") x = PainterModes::QUIVER;
    else if (j == "Scatter" || j == "scatter" || j == "SCATTER") x = PainterModes::SCATTER;
    else { throw std::runtime_error("Input JSON does not conform to PainterModes schema: " + to_string(j)); }
}

inline void to_json(json & j, const PainterModes & x) {
    switch (x) {
        case PainterModes::CONTOUR: j = "Contour"; break;
        case PainterModes::QUIVER: j = "Quiver"; break;
        case PainterModes::SCATTER: j = "Scatter"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"PainterModes\": " + std::to_string(static_cast<int>(x)));
    }
}

void from_json(const json & j, MagneticFieldStrengthModels & x);
void to_json(json & j, const MagneticFieldStrengthModels & x);

inline void to_json(json & j, const MagneticFieldStrengthModels & x) {
    switch (x) {
        case MagneticFieldStrengthModels::BINNS_LAWRENSON: j = "BinnsLawrenson"; break;
        case MagneticFieldStrengthModels::LAMMERANER: j = "Lammeraner"; break;
        case MagneticFieldStrengthModels::DOWELL: j = "Dowell"; break;
        case MagneticFieldStrengthModels::WANG: j = "Wang"; break;
        case MagneticFieldStrengthModels::ALBACH: j = "Albach"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"MagneticFieldStrengthModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, MagneticFieldStrengthModels & x) {
    if (j == "BinnsLawrenson" || j == "binnslawrenson" || j == "BINNSLAWRENSON") x = MagneticFieldStrengthModels::BINNS_LAWRENSON;
    else if (j == "Lammeraner" || j == "lammeraner" || j == "LAMMERANER") x = MagneticFieldStrengthModels::LAMMERANER;
    else if (j == "Dowell" || j == "dowell" || j == "DOWELL") x = MagneticFieldStrengthModels::DOWELL;
    else if (j == "Wang" || j == "wang" || j == "WANG") x = MagneticFieldStrengthModels::WANG;
    else if (j == "Albach" || j == "albach" || j == "ALBACH" || j == "Albach2D" || j == "albach2d" || j == "ALBACH_2D") x = MagneticFieldStrengthModels::ALBACH;
    else { throw std::runtime_error("Input JSON does not conform to MagneticFieldStrengthModels schema: " + to_string(j)); }
}

void from_json(const json & j, MagneticFieldStrengthFringingEffectModels & x);
void to_json(json & j, const MagneticFieldStrengthFringingEffectModels & x);

inline void to_json(json & j, const MagneticFieldStrengthFringingEffectModels & x) {
    switch (x) {
        case MagneticFieldStrengthFringingEffectModels::ROSHEN: j = "Roshen"; break;
        case MagneticFieldStrengthFringingEffectModels::ALBACH: j = "Albach"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"MagneticFieldStrengthFringingEffectModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, MagneticFieldStrengthFringingEffectModels & x) {
    if (j == "Roshen" || j == "roshen" || j == "ROSHEN") x = MagneticFieldStrengthFringingEffectModels::ROSHEN;
    else if (j == "Albach" || j == "albach" || j == "ALBACH") x = MagneticFieldStrengthFringingEffectModels::ALBACH;
    else { throw std::runtime_error("Input JSON does not conform to MagneticFieldStrengthFringingEffectModels schema: " + to_string(j)); }
}

void from_json(const json & j, ReluctanceModels & x);
void to_json(json & j, const ReluctanceModels & x);

inline void to_json(json & j, const ReluctanceModels & x) {
    switch (x) {
        case ReluctanceModels::ZHANG: j = "Zhang"; break;
        case ReluctanceModels::PARTRIDGE: j = "Partridge"; break;
        case ReluctanceModels::EFFECTIVE_AREA: j = "EffectiveArea"; break;
        case ReluctanceModels::EFFECTIVE_LENGTH: j = "EffectiveLength"; break;
        case ReluctanceModels::MUEHLETHALER: j = "Muehlethaler"; break;
        case ReluctanceModels::STENGLEIN: j = "Stenglein"; break;
        case ReluctanceModels::BALAKRISHNAN: j = "Balakrishnan"; break;
        case ReluctanceModels::CLASSIC: j = "Classic"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"ReluctanceModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, ReluctanceModels & x) {
    if (j == "Zhang" || j == "zhang" || j == "ZHANG") x = ReluctanceModels::ZHANG;
    else if (j == "Partridge" || j == "partridge" || j == "PARTRIDGE") x = ReluctanceModels::PARTRIDGE;
    else if (j == "EffectiveArea" || j == "effectivearea" || j == "EFFECTIVEAREA") x = ReluctanceModels::EFFECTIVE_AREA;
    else if (j == "EffectiveLength" || j == "effectivelength" || j == "EFFECTIVELENGTH") x = ReluctanceModels::EFFECTIVE_LENGTH;
    else if (j == "Muehlethaler" || j == "muehlethaler" || j == "MUEHLETHALER") x = ReluctanceModels::MUEHLETHALER;
    else if (j == "Stenglein" || j == "stenglein" || j == "STENGLEIN") x = ReluctanceModels::STENGLEIN;
    else if (j == "Balakrishnan" || j == "balakrishnan" || j == "BALAKRISHNAN") x = ReluctanceModels::BALAKRISHNAN;
    else if (j == "Classic" || j == "classic" || j == "CLASSIC") x = ReluctanceModels::CLASSIC;
    else { throw std::runtime_error("Input JSON does not conform to ReluctanceModels schema: " + to_string(j)); }
}

void from_json(const json & j, CoreLossesModels & x);
void to_json(json & j, const CoreLossesModels & x);

inline void to_json(json & j, const CoreLossesModels & x) {
    switch (x) {
        case CoreLossesModels::PROPRIETARY: j = "Proprietary"; break;
        case CoreLossesModels::STEINMETZ: j = "Steinmetz"; break;
        case CoreLossesModels::IGSE: j = "IGSE"; break;
        case CoreLossesModels::BARG: j = "Barg"; break;
        case CoreLossesModels::ROSHEN: j = "Roshen"; break;
        case CoreLossesModels::ALBACH: j = "Albach"; break;
        case CoreLossesModels::NSE: j = "NSE"; break;
        case CoreLossesModels::MSE: j = "MSE"; break;
        case CoreLossesModels::LOSS_FACTOR: j = "LossFactor"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"CoreLossesModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, CoreLossesModels & x) {
    if (j == "Proprietary" || j == "proprietary" || j == "PROPRIETARY") x = CoreLossesModels::PROPRIETARY;
    else if (j == "Steinmetz" || j == "steinmetz" || j == "STEINMETZ") x = CoreLossesModels::STEINMETZ;
    else if (j == "Igse" || j == "igse" || j == "IGSE") x = CoreLossesModels::IGSE;
    else if (j == "Barg" || j == "barg" || j == "BARG") x = CoreLossesModels::BARG;
    else if (j == "Roshen" || j == "roshen" || j == "ROSHEN") x = CoreLossesModels::ROSHEN;
    else if (j == "Albach" || j == "albach" || j == "ALBACH") x = CoreLossesModels::ALBACH;
    else if (j == "Nse" || j == "nse" || j == "NSE") x = CoreLossesModels::NSE;
    else if (j == "Mse" || j == "mse" || j == "MSE") x = CoreLossesModels::MSE;
    else if (j == "LossFactor" || j == "lossfactor" || j == "LOSSFACTOR") x = CoreLossesModels::LOSS_FACTOR;
    else { throw std::runtime_error("Input JSON does not conform to CoreLossesModels schema: " + to_string(j)); }
}

void from_json(const json & j, CoreThermalResistanceModels & x);
void to_json(json & j, const CoreThermalResistanceModels & x);

inline void to_json(json & j, const CoreThermalResistanceModels & x) {
    switch (x) {
        case CoreThermalResistanceModels::MANIKTALA: j = "Maniktala"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"CoreThermalResistanceModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, CoreThermalResistanceModels & x) {
    if (j == "Maniktala" || j == "maniktala" || j == "MANIKTALA") x = CoreThermalResistanceModels::MANIKTALA;
    else { throw std::runtime_error("Input JSON does not conform to CoreThermalResistanceModels schema: " + to_string(j)); }
}

void from_json(const json & j, CoreTemperatureModels & x);
void to_json(json & j, const CoreTemperatureModels & x);

inline void to_json(json & j, const CoreTemperatureModels & x) {
    switch (x) {
        case CoreTemperatureModels::KAZIMIERCZUK: j = "Kazimierczuk"; break;
        case CoreTemperatureModels::MANIKTALA: j = "Maniktala"; break;
        case CoreTemperatureModels::TDK: j = "Tdk"; break;
        case CoreTemperatureModels::DIXON: j = "Dixon"; break;
        case CoreTemperatureModels::AMIDON: j = "Amidon"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"CoreTemperatureModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, CoreTemperatureModels & x) {
    if (j == "Kazimierczuk" || j == "kazimierczuk" || j == "KAZIMIERCZUK") x = CoreTemperatureModels::KAZIMIERCZUK;
    else if (j == "Maniktala" || j == "maniktala" || j == "MANIKTALA") x = CoreTemperatureModels::MANIKTALA;
    else if (j == "Tdk" || j == "tdk" || j == "TDK") x = CoreTemperatureModels::TDK;
    else if (j == "Dixon" || j == "dixon" || j == "DIXON") x = CoreTemperatureModels::DIXON;
    else if (j == "Amidon" || j == "amidon" || j == "AMIDON") x = CoreTemperatureModels::AMIDON;
    else { throw std::runtime_error("Input JSON does not conform to CoreTemperatureModels schema: " + to_string(j)); }
}

void from_json(const json & j, WindingSkinEffectLossesModels & x);
void to_json(json & j, const WindingSkinEffectLossesModels & x);

inline void to_json(json & j, const WindingSkinEffectLossesModels & x) {
    switch (x) {
        case WindingSkinEffectLossesModels::DOWELL: j = "Dowell"; break;
        case WindingSkinEffectLossesModels::WOJDA: j = "Wojda"; break;
        case WindingSkinEffectLossesModels::ALBACH: j = "Albach"; break;
        case WindingSkinEffectLossesModels::PAYNE: j = "Payne"; break;
        case WindingSkinEffectLossesModels::LOTFI: j = "Lotfi"; break;
        case WindingSkinEffectLossesModels::KAZIMIERCZUK: j = "Kazimierczuk"; break;
        case WindingSkinEffectLossesModels::KUTKUT: j = "Kutkut"; break;
        case WindingSkinEffectLossesModels::FERREIRA: j = "Ferreira"; break;
        case WindingSkinEffectLossesModels::DIMITRAKAKIS: j = "Dimitrakakis"; break;
        case WindingSkinEffectLossesModels::WANG: j = "Wang"; break;
        case WindingSkinEffectLossesModels::HOLGUIN: j = "Holguin"; break;
        case WindingSkinEffectLossesModels::PERRY: j = "Perry"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"WindingSkinEffectLossesModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, WindingSkinEffectLossesModels & x) {
    if (j == "Dowell" || j == "dowell" || j == "DOWELL") x = WindingSkinEffectLossesModels::DOWELL;
    else if (j == "Wojda" || j == "wojda" || j == "WOJDA") x = WindingSkinEffectLossesModels::WOJDA;
    else if (j == "Albach" || j == "albach" || j == "ALBACH") x = WindingSkinEffectLossesModels::ALBACH;
    else if (j == "Payne" || j == "payne" || j == "PAYNE") x = WindingSkinEffectLossesModels::PAYNE;
    else if (j == "Lotfi" || j == "lotfi" || j == "LOTFI") x = WindingSkinEffectLossesModels::LOTFI;
    else if (j == "Kazimierczuk" || j == "kazimierczuk" || j == "KAZIMIERCZUK") x = WindingSkinEffectLossesModels::KAZIMIERCZUK;
    else if (j == "Kutkut" || j == "kutkut" || j == "KUTKUT") x = WindingSkinEffectLossesModels::KUTKUT;
    else if (j == "Ferreira" || j == "ferreira" || j == "FERREIRA") x = WindingSkinEffectLossesModels::FERREIRA;
    else if (j == "Dimitrakakis" || j == "dimitrakakis" || j == "DIMITRAKAKIS") x = WindingSkinEffectLossesModels::DIMITRAKAKIS;
    else if (j == "Wang" || j == "wang" || j == "WANG") x = WindingSkinEffectLossesModels::WANG;
    else if (j == "Holguin" || j == "holguin" || j == "HOLGUIN") x = WindingSkinEffectLossesModels::HOLGUIN;
    else if (j == "Perry" || j == "perry" || j == "PERRY") x = WindingSkinEffectLossesModels::PERRY;
    else { throw std::runtime_error("Input JSON does not conform to WindingSkinEffectLossesModels schema: " + to_string(j)); }
}

void from_json(const json & j, WindingProximityEffectLossesModels & x);
void to_json(json & j, const WindingProximityEffectLossesModels & x);

inline void to_json(json & j, const WindingProximityEffectLossesModels & x) {
    switch (x) {
        case WindingProximityEffectLossesModels::ROSSMANITH: j = "Rossmanith"; break;
        case WindingProximityEffectLossesModels::WANG: j = "Wang"; break;
        case WindingProximityEffectLossesModels::FERREIRA: j = "Ferreira"; break;
        case WindingProximityEffectLossesModels::LAMMERANER: j = "Lammeraner"; break;
        case WindingProximityEffectLossesModels::ALBACH: j = "Albach"; break;
        case WindingProximityEffectLossesModels::DOWELL: j = "Dowell"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"WindingProximityEffectLossesModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, WindingProximityEffectLossesModels & x) {
    if (j == "Rossmanith" || j == "rossmanith" || j == "ROSSMANITH") x = WindingProximityEffectLossesModels::ROSSMANITH;
    else if (j == "Wang" || j == "wang" || j == "WANG") x = WindingProximityEffectLossesModels::WANG;
    else if (j == "Ferreira" || j == "ferreira" || j == "FERREIRA") x = WindingProximityEffectLossesModels::FERREIRA;
    else if (j == "Lammeraner" || j == "lammeraner" || j == "LAMMERANER") x = WindingProximityEffectLossesModels::LAMMERANER;
    else if (j == "Albach" || j == "albach" || j == "ALBACH") x = WindingProximityEffectLossesModels::ALBACH;
    else if (j == "Dowell" || j == "dowell" || j == "DOWELL") x = WindingProximityEffectLossesModels::DOWELL;
    else { throw std::runtime_error("Input JSON does not conform to WindingProximityEffectLossesModels schema: " + to_string(j)); }
}

void from_json(const json & j, StrayCapacitanceModels & x);
void to_json(json & j, const StrayCapacitanceModels & x);

inline void to_json(json & j, const StrayCapacitanceModels & x) {
    switch (x) {
        case StrayCapacitanceModels::KOCH: j = "Koch"; break;
        case StrayCapacitanceModels::ALBACH: j = "Albach"; break;
        case StrayCapacitanceModels::DUERDOTH: j = "Duerdoth"; break;
        case StrayCapacitanceModels::MASSARINI: j = "Massarini"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"StrayCapacitanceModels\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json & j, StrayCapacitanceModels & x) {
    if (j == "Koch" || j == "koch" || j == "KOCH") x = StrayCapacitanceModels::KOCH;
    else if (j == "Albach" || j == "albach" || j == "ALBACH") x = StrayCapacitanceModels::ALBACH;
    else if (j == "Duerdoth" || j == "duerdoth" || j == "DUERDOTH") x = StrayCapacitanceModels::DUERDOTH;
    else if (j == "Massarini" || j == "massarini" || j == "MASSARINI") x = StrayCapacitanceModels::MASSARINI;
    else { throw std::runtime_error("Input JSON does not conform to StrayCapacitanceModels schema: " + to_string(j)); }
}


class Curve2D {

    std::vector<double> xPoints;
    std::vector<double> yPoints;
    std::string title;

    public: 
        Curve2D() = default;
        virtual ~Curve2D() = default;
        Curve2D(std::vector<double> xPoints, std::vector<double> yPoints, std::string title) : xPoints(xPoints), yPoints(yPoints), title(title){}

        const std::vector<double> & get_x_points() const { return xPoints; }
        std::vector<double> & get_mutable_x_points() { return xPoints; }
        void set_x_points(const std::vector<double> & value) { this->xPoints = value; }

        const std::vector<double> & get_y_points() const { return yPoints; }
        std::vector<double> & get_mutable_y_points() { return yPoints; }
        void set_y_points(const std::vector<double> & value) { this->yPoints = value; }

        const std::string & get_title() const { return title; }
        std::string & get_mutable_title() { return title; }
        void set_title(const std::string & value) { this->title = value; }
};
                
void from_json(const json & j, Curve2D & x);
void to_json(json & j, const Curve2D & x);

inline void from_json(const json & j, Curve2D& x) {
    x.set_x_points(j.at("xPoints").get<std::vector<double>>());
    x.set_y_points(j.at("yPoints").get<std::vector<double>>());
    x.set_title(j.at("title").get<std::string>());
}

inline void to_json(json & j, const Curve2D & x) {
    j = json::object();
    j["xPoints"] = x.get_x_points();
    j["yPoints"] = x.get_y_points();
    j["title"] = x.get_title();
}

enum class MagneticFilters : int {
    AREA_PRODUCT,
    ENERGY_STORED,
    ESTIMATED_COST,
    COST,
    CORE_AND_DC_LOSSES,
    CORE_DC_AND_SKIN_LOSSES,
    LOSSES,
    LOSSES_NO_PROXIMITY,
    DIMENSIONS,
    CORE_MINIMUM_IMPEDANCE,
    AREA_NO_PARALLELS,
    AREA_WITH_PARALLELS,
    EFFECTIVE_RESISTANCE,
    PROXIMITY_FACTOR,
    SOLID_INSULATION_REQUIREMENTS,
    TURNS_RATIOS,
    MAXIMUM_DIMENSIONS,
    SATURATION,
    DC_CURRENT_DENSITY,
    EFFECTIVE_CURRENT_DENSITY,
    IMPEDANCE,
    MAGNETIZING_INDUCTANCE,
    FRINGING_FACTOR,
    SKIN_LOSSES_DENSITY,
    VOLUME,
    AREA,
    HEIGHT,
    TEMPERATURE_RISE,
    VOLUME_TIMES_TEMPERATURE_RISE,
    LOSSES_TIMES_VOLUME,
    LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE,
    LOSSES_NO_PROXIMITY_TIMES_VOLUME,
    LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE,
    MAGNETOMOTIVE_FORCE,
    LEAKAGE_INDUCTANCE  // For CMC optimization - minimize leakage for tight coupling
};

class MagneticFilterOperation {
    public:
    MagneticFilterOperation(MagneticFilters filter, bool invert, bool log, double weight) : filter(filter), invert(invert), log(log), weight(weight) {};
    MagneticFilterOperation(MagneticFilters filter, bool invert, bool log, bool strictlyRequired, double weight) : filter(filter), invert(invert), log(log), strictlyRequired(strictlyRequired), weight(weight) {};
    MagneticFilterOperation() {};
    virtual ~MagneticFilterOperation() = default;

    private:
    MagneticFilters filter = MagneticFilters::DIMENSIONS;
    bool invert = true;
    bool log = false;
    bool strictlyRequired = false;
    double weight = 1;

    public:

    MagneticFilters get_filter() const { return filter; }
    void set_filter(const MagneticFilters & value) { this->filter = value; }

    bool get_invert() const { return invert; }
    void set_invert(const bool & value) { this->invert = value; }

    bool get_log() const { return log; }
    void set_log(const bool & value) { this->log = value; }

    bool get_strictly_required() const { return strictlyRequired; }
    void set_strictly_required(const bool & value) { this->strictlyRequired = value; }

    double get_weight() const { return weight; }
    void set_weight(const double & value) { this->weight = value; }
};


void from_json(const json & j, MagneticFilters & x);
void to_json(json & j, const MagneticFilters & x);
void from_json(const json& j, MagneticFilterOperation& x);
void to_json(json& j, const MagneticFilterOperation& x);
void from_json(const json& j, std::vector<MagneticFilterOperation>& v);
void to_json(json& j, const std::vector<MagneticFilterOperation>& v);


inline void from_json(const json & j, MagneticFilters & x) {
    if (j == "Area Product") x = MagneticFilters::AREA_PRODUCT;
    else if (j == "Energy Stored") x = MagneticFilters::ENERGY_STORED;
    else if (j == "Cost") x = MagneticFilters::COST;
    else if (j == "Estimated Cost") x = MagneticFilters::ESTIMATED_COST;
    else if (j == "Core And DC Losses") x = MagneticFilters::CORE_AND_DC_LOSSES;
    else if (j == "Core DC And Skin Losses") x = MagneticFilters::CORE_DC_AND_SKIN_LOSSES;
    else if (j == "Losses") x = MagneticFilters::LOSSES;
    else if (j == "Losses No Proximity") x = MagneticFilters::LOSSES_NO_PROXIMITY;
    else if (j == "Dimensions") x = MagneticFilters::DIMENSIONS;
    else if (j == "Core Minimum Impedance") x = MagneticFilters::CORE_MINIMUM_IMPEDANCE;
    else if (j == "Area No Parallels") x = MagneticFilters::AREA_NO_PARALLELS;
    else if (j == "Area With Parallels") x = MagneticFilters::AREA_WITH_PARALLELS;
    else if (j == "Effective Resistance") x = MagneticFilters::EFFECTIVE_RESISTANCE;
    else if (j == "Proximity Factor") x = MagneticFilters::PROXIMITY_FACTOR;
    else if (j == "Solid Insulation Requirements") x = MagneticFilters::SOLID_INSULATION_REQUIREMENTS;
    else if (j == "Turns Ratios") x = MagneticFilters::TURNS_RATIOS;
    else if (j == "Maximum Dimensions") x = MagneticFilters::MAXIMUM_DIMENSIONS;
    else if (j == "Saturation") x = MagneticFilters::SATURATION;
    else if (j == "Dc Current Density") x = MagneticFilters::DC_CURRENT_DENSITY;
    else if (j == "Effective Current Density") x = MagneticFilters::EFFECTIVE_CURRENT_DENSITY;
    else if (j == "Impedance") x = MagneticFilters::IMPEDANCE;
    else if (j == "Magnetizing Inductance") x = MagneticFilters::MAGNETIZING_INDUCTANCE;
    else if (j == "Fringing Factor") x = MagneticFilters::FRINGING_FACTOR;
    else if (j == "Skin Losses Density") x = MagneticFilters::SKIN_LOSSES_DENSITY;
    else if (j == "Volume") x = MagneticFilters::VOLUME;
    else if (j == "Area") x = MagneticFilters::AREA;
    else if (j == "Height") x = MagneticFilters::HEIGHT;
    else if (j == "Temperature Rise") x = MagneticFilters::TEMPERATURE_RISE;
    else if (j == "Losses Times Volume") x = MagneticFilters::LOSSES_TIMES_VOLUME;
    else if (j == "Volume Times Temperature Rise") x = MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE;
    else if (j == "Losses Times Volume Times Temperature Rise") x = MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE;
    else if (j == "Losses No Proximity Times Volume") x = MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME;
    else if (j == "Losses No Proximity Times Volume Times Temperature Rise") x = MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE;
    else if (j == "MagnetomotiveForce") x = MagneticFilters::MAGNETOMOTIVE_FORCE;
    else { throw std::runtime_error("Input JSON does not conform to MagneticFilters schema!"); }
}

inline void to_json(json & j, const MagneticFilters & x) {
    switch (x) {
        case MagneticFilters::AREA_PRODUCT: j = "Area Product"; break;
        case MagneticFilters::ENERGY_STORED: j = "Energy Stored"; break;
        case MagneticFilters::COST: j = "Cost"; break;
        case MagneticFilters::ESTIMATED_COST: j = "Estimated Cost"; break;
        case MagneticFilters::CORE_AND_DC_LOSSES: j = "Core And DC Losses"; break;
        case MagneticFilters::CORE_DC_AND_SKIN_LOSSES: j = "Core DC And Skin Losses"; break;
        case MagneticFilters::LOSSES: j = "Losses"; break;
        case MagneticFilters::LOSSES_NO_PROXIMITY: j = "Losses No Proximity"; break;
        case MagneticFilters::DIMENSIONS: j = "Dimensions"; break;
        case MagneticFilters::CORE_MINIMUM_IMPEDANCE: j = "Core Minimum Impedance"; break;
        case MagneticFilters::AREA_NO_PARALLELS: j = "Area No Parallels"; break;
        case MagneticFilters::AREA_WITH_PARALLELS: j = "Area With Parallels"; break;
        case MagneticFilters::EFFECTIVE_RESISTANCE: j = "Effective Resistance"; break;
        case MagneticFilters::PROXIMITY_FACTOR: j = "Proximity Factor"; break;
        case MagneticFilters::SOLID_INSULATION_REQUIREMENTS: j = "Solid Insulation Requirements"; break;
        case MagneticFilters::TURNS_RATIOS: j = "Turns Ratios"; break;
        case MagneticFilters::MAXIMUM_DIMENSIONS: j = "Maximum Dimensions"; break;
        case MagneticFilters::SATURATION: j = "Saturation"; break;
        case MagneticFilters::DC_CURRENT_DENSITY: j = "Dc Current Density"; break;
        case MagneticFilters::EFFECTIVE_CURRENT_DENSITY: j = "Effective Current Density"; break;
        case MagneticFilters::IMPEDANCE: j = "Impedance"; break;
        case MagneticFilters::MAGNETIZING_INDUCTANCE: j = "Magnetizing Inductance"; break;
        case MagneticFilters::FRINGING_FACTOR: j = "Fringing Factor"; break;
        case MagneticFilters::SKIN_LOSSES_DENSITY: j = "Skin Losses Density"; break;
        case MagneticFilters::VOLUME: j = "Volume"; break;
        case MagneticFilters::AREA: j = "Area"; break;
        case MagneticFilters::HEIGHT: j = "Height"; break;
        case MagneticFilters::TEMPERATURE_RISE: j = "Temperature Rise"; break;
        case MagneticFilters::LOSSES_TIMES_VOLUME: j = "Losses Times Volume"; break;
        case MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE: j = "Volume Times Temperature Rise"; break;
        case MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE: j = "Losses Times Volume Times Temperature Rise"; break;
        case MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME: j = "Losses No Proximity Times Volume"; break;
        case MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE: j = "Losses No Proximity Times Volume Times Temperature Rise"; break;
        case MagneticFilters::MAGNETOMOTIVE_FORCE: j = "MagnetomotiveForce"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"[object Object]\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json& j, MagneticFilterOperation& x) {
    x.set_filter(j.at("filter").get<MagneticFilters>());
    x.set_invert(j.at("invert").get<bool>());
    x.set_log(j.at("log").get<bool>());
    x.set_strictly_required(j.at("strictlyRequired").get<bool>());
    x.set_weight(j.at("weight").get<double>());
}

inline void to_json(json& j, const MagneticFilterOperation& x) {
    j = json::object();
    j["filter"] = x.get_filter();
    j["invert"] = x.get_invert();
    j["log"] = x.get_log();
    j["strictlyRequired"] = x.get_strictly_required();
    j["weight"] = x.get_weight();
}

inline void from_json(const json& j, std::vector<MagneticFilterOperation>& v) {
    for (auto e : j) {
        MagneticFilterOperation x(
            e.at("filter").get<MagneticFilters>(),
            e.at("invert").get<bool>(),
            e.at("log").get<bool>(),
            e.at("strictlyRequired").get<bool>(),
            e.at("weight").get<double>()
        );
        v.push_back(x);
    }
}

inline void to_json(json& j, const std::vector<MagneticFilterOperation>& v) {
    j = json::array();
    for (auto x : v) {
        json e;
        e["filter"] = x.get_filter();
        e["invert"] = x.get_invert();
        e["log"] = x.get_log();
        e["strictlyRequired"] = x.get_strictly_required();
        e["weight"] = x.get_weight();
        j.push_back(e);
    }
}

double resolve_dimensional_values(MAS::Dimension dimensionValue, DimensionalValues preferredValue = DimensionalValues::NOMINAL);
inline double resolve_dimensional_values(MAS::Dimension dimensionValue, DimensionalValues preferredValue) {
    double doubleValue = 0;
    if (std::holds_alternative<MAS::DimensionWithTolerance>(dimensionValue)) {
        switch (preferredValue) {
            case DimensionalValues::MAXIMUM:
                if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                break;
            case DimensionalValues::NOMINAL:
                if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().has_value() &&
                         std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue =
                        (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().value() +
                         std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().value()) /
                        2;
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                break;
            case DimensionalValues::MINIMUM:
                if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                break;
            default:
                throw std::runtime_error("Unknown type of dimension, options are {MAXIMUM, NOMINAL, MINIMUM}");
        }
    }
    else if (std::holds_alternative<double>(dimensionValue)) {
        doubleValue = std::get<double>(dimensionValue);
    }
    else {
        throw std::runtime_error("Unknown variant in dimensionValue, holding variant");
    }
    return doubleValue;
}


class MagneticCoreSearchElement {
    public:
        MagneticCoreSearchElement() = default;
        virtual ~MagneticCoreSearchElement() = default;
    private:
        std::string name;
        std::optional<double> effective_length;
        std::optional<double> effective_area;
        std::optional<double> effective_volume;
        std::optional<double> winding_window_area;
        std::optional<double> area_product;
        std::optional<double> winding_width;
        std::optional<double> winding_height;

    public:

        const std::string & get_name() const { return name; }
        std::string & get_mutable_name() { return name; }
        void set_name(const std::string & value) { this->name = value; }

        std::optional<double> get_effective_length() const { return effective_length; }
        void set_effective_length(std::optional<double> value) { this->effective_length = value; }
        
        std::optional<double> get_effective_area() const { return effective_area; }
        void set_effective_area(std::optional<double> value) { this->effective_area = value; }
        
        std::optional<double> get_effective_volume() const { return effective_volume; }
        void set_effective_volume(std::optional<double> value) { this->effective_volume = value; }
        
        std::optional<double> get_winding_window_area() const { return winding_window_area; }
        void set_winding_window_area(std::optional<double> value) { this->winding_window_area = value; }
        
        std::optional<double> get_area_product() const { return area_product; }
        void set_area_product(std::optional<double> value) { this->area_product = value; }
        
        std::optional<double> get_winding_width() const { return winding_width; }
        void set_winding_width(std::optional<double> value) { this->winding_width = value; }
        
        std::optional<double> get_winding_height() const { return winding_height; }
        void set_winding_height(std::optional<double> value) { this->winding_height = value; }
};

}