#pragma once
#include <MAS.hpp>
#include "Utils.h"
#include "Impedance.h"
#include "MagneticWrapper.h"
#include <complex>

namespace OpenMagnetics {

class Sweeper {
    private:
    protected:
    public:
    Curve2D sweep_impedance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, std::string title = "Impedance over frequency");

};

} // namespace OpenMagnetics