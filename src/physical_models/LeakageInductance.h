#pragma once
#include "Defaults.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include "support/CoilMesher.h"
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
    std::pair<size_t, size_t> calculate_number_points_needed_for_leakage(Coil coil);

    private:
        static constexpr double NEGLIGIBLE_CURRENT = 1e-9;
        static constexpr double PLANAR_THICKNESS_RATIO_THRESHOLD = 10.0;
        static constexpr double SINUSOIDAL_PEAK_TO_PEAK = 2.0;
        static constexpr double SINUSOIDAL_DUTY_CYCLE = 0.5;
        static constexpr double SINUSOIDAL_OFFSET = 0.0;
        static constexpr double DEGREES_IN_CIRCLE = 360.0;
        static constexpr double MINIMUM_PRECISION_LEVEL = 1.0;

        OperatingPoint create_leakage_operating_point(Magnetic& magnetic, size_t sourceIndex, size_t destinationIndex, double frequency);
        CoilMesherModels select_mesh_model(Magnetic& magnetic);
        std::pair<size_t, size_t> calculate_grid_points(Magnetic& magnetic, double frequency);

};
} // namespace OpenMagnetics
