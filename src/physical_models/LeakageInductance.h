#pragma once
#include "Defaults.h"
#include "constructive_models/MagneticWrapper.h"
#include "support/Utils.h"
#include <MAS.hpp>

namespace OpenMagnetics {

class LeakageInductance{
    public:

        LeakageInductance(){
        };
        virtual ~LeakageInductance() = default;


    LeakageInductanceOutput calculate_leakage_inductance(MagneticWrapper magnetic, double frequency, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1);
    ComplexField calculate_leakage_magnetic_field(MagneticWrapper magnetic, double frequency, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1);
    LeakageInductanceOutput calculate_leakage_inductance_all_windings(MagneticWrapper magnetic, double frequency, size_t sourceIndex = 0, size_t harmonicIndex = 1);
    std::pair<ComplexField, double> calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1);

};
} // namespace OpenMagnetics
