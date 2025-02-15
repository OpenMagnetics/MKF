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
    static Curve2D sweep_impedance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, std::string mode="log", std::string title = "Impedance over frequency");
    static Curve2D sweep_winding_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, size_t windingIndex, double temperature = Defaults().ambientTemperature, std::string mode="log", std::string title = "Winding Resistance over frequency");
    static Curve2D sweep_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, double temperature = Defaults().ambientTemperature, std::string mode="log", std::string title = "Resistance over frequency");
    static Curve2D sweep_core_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, double temperature = Defaults().ambientTemperature, std::string mode="log", std::string title = "Core Resistance over frequency");
    static Curve2D sweep_core_losses_over_frequency(MagneticWrapper magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature = Defaults().ambientTemperature, std::string mode="log", std::string title = "Core Losses over frequency");
    static Curve2D sweep_winding_losses_over_frequency(MagneticWrapper magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature = Defaults().ambientTemperature, std::string mode="log", std::string title = "Winding Losses over frequency");

};

} // namespace OpenMagnetics