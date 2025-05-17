#pragma once
#include "Defaults.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class LeakageInductance{
    public:

        LeakageInductance(){
        };
        virtual ~LeakageInductance() = default;


    LeakageInductanceOutput calculate_leakage_inductance(Magnetic magnetic, double frequency, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1);
    ComplexField calculate_leakage_magnetic_field(Magnetic magnetic, double frequency, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1);
    LeakageInductanceOutput calculate_leakage_inductance_all_windings(Magnetic magnetic, double frequency, size_t sourceIndex = 0, size_t harmonicIndex = 1);
    std::pair<ComplexField, double> calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1);

};
} // namespace OpenMagnetics
