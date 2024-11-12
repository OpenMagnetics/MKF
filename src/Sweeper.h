#pragma once
#include <MAS.hpp>
#include "Utils.h"
#include "Defaults.h"
#include "Impedance.h"
#include "MagneticWrapper.h"
#include <complex>

namespace OpenMagnetics {

class Sweeper {
    private:
    protected:
    public:
    Curve2D sweep_impedance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, std::string title = "Impedance over frequency");

    Curve2D sweep_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, size_t windingIndex, double temperature = Defaults().ambientTemperature, std::string title = "Resistance over frequency");

};

} // namespace OpenMagnetics