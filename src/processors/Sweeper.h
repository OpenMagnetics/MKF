#pragma once
#include <MAS.hpp>
#include "support/Utils.h"
#include "Defaults.h"
#include "physical_models/Impedance.h"
#include "constructive_models/Magnetic.h"
#include <complex>

using namespace MAS;

namespace OpenMagnetics {

class Sweeper {
    private:
    protected:
    public:
    static Curve2D sweep_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode="log", std::string title = "Impedance over frequency");
    static Curve2D sweep_q_factor_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode="log", std::string title = "Impedance over frequency");
    static Curve2D sweep_magnetizing_inductance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Magnetizing Inductance over frequency");
    static Curve2D sweep_magnetizing_inductance_over_temperature(Magnetic magnetic, double start, double stop, size_t numberElements, double frequency = defaults.measurementFrequency, std::string mode="linear", std::string title = "Magnetizing Inductance over temperature");
    static Curve2D sweep_magnetizing_inductance_over_dc_bias(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="linear", std::string title = "Magnetizing Inductance over DC bias");
    static Curve2D sweep_winding_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, size_t windingIndex, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Winding Resistance over frequency");
    static Curve2D sweep_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Resistance over frequency");
    static Curve2D sweep_core_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Core Resistance over frequency");
    static Curve2D sweep_core_losses_over_frequency(Magnetic magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Core Losses over frequency");
    static Curve2D sweep_winding_losses_over_frequency(Magnetic magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Winding Losses over frequency");

};

} // namespace OpenMagnetics