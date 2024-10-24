#pragma once
#include "Constants.h"
#include "MagneticWrapper.h"
#include "CoreWrapper.h"
#include "CoilWrapper.h"
#include <MAS.hpp>
#include <complex>

namespace OpenMagnetics {

class Impedance {
    private:
    protected:
    public:
    std::complex<double> calculate_impedance(MagneticWrapper magnetic, double frequency);
    std::complex<double> calculate_impedance(CoreWrapper core, CoilWrapper coil, double frequency);

};

} // namespace OpenMagnetics