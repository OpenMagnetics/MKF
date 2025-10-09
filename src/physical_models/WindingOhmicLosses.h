#pragma once
#include "MAS.hpp"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

class WindingOhmicLosses {
    private:
    protected:
    public:
        WindingOhmicLosses() {
        }
        static double calculate_effective_resistance_per_meter(Wire wire, double frequency, double temperature);
        static double calculate_dc_resistance_per_meter(Wire wire, double temperature);
        static double calculate_dc_resistance(Turn turn, const Wire& wire, double temperature);
        static double calculate_dc_resistance(double wireLength, const Wire& wire, double temperature);
        static std::vector<double> calculate_dc_resistance_per_winding(Coil coil, double temperature);
        static WindingLossesOutput calculate_ohmic_losses(Coil coil, OperatingPoint operatingPoint, double temperature);
        static double calculate_ohmic_losses_per_meter(Wire wire, SignalDescriptor current, double temperature);
};

} // namespace OpenMagnetics