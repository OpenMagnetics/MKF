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

// Frequency-independent building blocks of the differential-mode impedance:
// the leakage inductance, the winding DC resistance and the inter-winding
// (through-core) capacitance. Computed once per magnetic so a frequency sweep
// does not recompute the (expensive) stray-capacitance model at every point.
struct DifferentialModeParameters {
    double leakageInductance;
    double windingResistance;
    double interWindingCapacitance;
};

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
    std::complex<double> calculate_impedance(Core core, Coil coil, double frequency, double magneticFieldDcBias, double temperature);
    // Differential-mode impedance of a coupled (e.g. common-mode-choke) magnetic:
    // the opposing winding currents cancel the core flux, so the inductance seen
    // is the leakage inductance, which parallel-resonates with the inter-winding
    // capacitance (the second self-resonance).
    std::complex<double> calculate_differential_mode_impedance(Magnetic magnetic, double frequency, double temperature = Defaults().ambientTemperature);
    std::complex<double> calculate_differential_mode_impedance(Core core, Coil coil, double frequency, double temperature = Defaults().ambientTemperature);
    // Compute the frequency-independent DM parameters once (the leakage is taken at
    // referenceFrequency — it is essentially air-cored and flat). Then evaluate the
    // impedance cheaply at each frequency. A sweep should call the first once and the
    // second per point rather than calling calculate_differential_mode_impedance N times.
    DifferentialModeParameters calculate_differential_mode_parameters(Core core, Coil coil, double referenceFrequency, double temperature = Defaults().ambientTemperature);
    std::complex<double> differential_mode_impedance_from_parameters(const DifferentialModeParameters& parameters, double frequency);
    double calculate_q_factor(Magnetic magnetic, double frequency, double temperature = Defaults().ambientTemperature);
    double calculate_q_factor(Core core, Coil coil, double frequency, double temperature = Defaults().ambientTemperature);
    double calculate_self_resonant_frequency(Magnetic magnetic, double temperature = Defaults().ambientTemperature);
    double calculate_self_resonant_frequency(Core core, Coil coil, double temperature = Defaults().ambientTemperature);
    int64_t calculate_minimum_number_turns(Magnetic magnetic, Inputs inputs);

};

} // namespace OpenMagnetics