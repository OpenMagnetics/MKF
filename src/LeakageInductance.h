#pragma once
#include "Defaults.h"
#include "MagneticWrapper.h"
#include "Utils.h"
#include <MAS.hpp>

namespace OpenMagnetics {

class LeakageInductance{
    protected:
        size_t _numberPointsX = 25;
        size_t _numberPointsY = 50;
        bool _includeFringing = false;
        int _mirroringDimension = Defaults().magneticFieldMirroringDimension;

    public:

        LeakageInductance(){
        };
        virtual ~LeakageInductance() = default;

    void set_number_points_x(double numberPoints) {
        _numberPointsX = numberPoints;
    }
    void set_number_points_y(double numberPoints) {
        _numberPointsY = numberPoints;
    }
    void set_fringing_effect(bool value) {
        _includeFringing = value;
    }
    void set_mirroring_dimension(int mirroringDimension) {
        _mirroringDimension = mirroringDimension;
    }
    
    LeakageInductanceOutput calculate_leakage_inductance(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);
    ComplexField calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex = 1);

};
} // namespace OpenMagnetics
