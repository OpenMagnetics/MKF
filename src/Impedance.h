#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "MagnetizingInductance.h"
#include "MagneticWrapper.h"
#include "CoreLosses.h"
#include "CoreWrapper.h"
#include "CoilWrapper.h"
#include <MAS.hpp>
#include <complex>

namespace OpenMagnetics {

class Impedance {
    private:
        MagnetizingInductance _magnetizingInductanceModel;
        CoreLosses _coreLossesModel;
    protected:
    public:
    std::complex<double> calculate_impedance(MagneticWrapper magnetic, double frequency, double temperature = Defaults().ambientTemperature);
    std::complex<double> calculate_impedance(CoreWrapper core, CoilWrapper coil, double frequency, double temperature = Defaults().ambientTemperature);
    double calculate_self_resonant_frequency(MagneticWrapper magnetic, double temperature = Defaults().ambientTemperature);
    double calculate_self_resonant_frequency(CoreWrapper core, CoilWrapper coil, double temperature = Defaults().ambientTemperature);

};

} // namespace OpenMagnetics