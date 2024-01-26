#pragma once
#include "Defaults.h"
#include "MagneticWrapper.h"
#include "Utils.h"
#include <MAS.hpp>

namespace OpenMagnetics {

class LeakageInductance{
    public:

        LeakageInductance(){
        };
        virtual ~LeakageInductance() = default;


    LeakageInductanceOutput calculate_leakage_inductance(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);
    ComplexField calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex = 1);

};
} // namespace OpenMagnetics
