#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "physical_models/MagnetizingInductance.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/CoreLosses.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include <MAS.hpp>
#include <complex>

using namespace MAS;

namespace OpenMagnetics {

class Impedance {
    private:
        bool _fastCapacitance;
        size_t _maximumNumberTurns = 200;
    protected:
    public:
    Impedance(bool fastCapacitance=true) {
        _fastCapacitance = fastCapacitance;
    } 

    virtual ~Impedance() = default;

    std::complex<double> calculate_impedance(Magnetic magnetic, double frequency, double temperature = Defaults().ambientTemperature);
    std::complex<double> calculate_impedance(Core core, Coil coil, double frequency, double temperature = Defaults().ambientTemperature);
    double calculate_q_factor(Magnetic magnetic, double frequency, double temperature = Defaults().ambientTemperature);
    double calculate_q_factor(Core core, Coil coil, double frequency, double temperature = Defaults().ambientTemperature);
    double calculate_self_resonant_frequency(Magnetic magnetic, double temperature = Defaults().ambientTemperature);
    double calculate_self_resonant_frequency(Core core, Coil coil, double temperature = Defaults().ambientTemperature);
    int64_t calculate_minimum_number_turns(Magnetic magnetic, Inputs inputs);

};

} // namespace OpenMagnetics