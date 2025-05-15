#pragma once
#include "MAS.hpp"
#include "constructive_models/CoilWrapper.h"
#include "constructive_models/WireWrapper.h"
#include "support/Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

class WindingOhmicLosses {
    private:
    protected:
    public:
        WindingOhmicLosses() {
        }
        static double calculate_effective_resistance_per_meter(WireWrapper wire, double frequency, double temperature);
        static double calculate_dc_resistance_per_meter(WireWrapper wire, double temperature);
        static double calculate_dc_resistance(Turn turn, const WireWrapper& wire, double temperature);
        static double calculate_dc_resistance(double wireLength, const WireWrapper& wire, double temperature);
        static std::vector<double> calculate_dc_resistance_per_winding(CoilWrapper coil, double temperature);
        static WindingLossesOutput calculate_ohmic_losses(CoilWrapper coil, OperatingPoint operatingPoint, double temperature);
        static double calculate_ohmic_losses_per_meter(WireWrapper wire, SignalDescriptor current, double temperature);
};

} // namespace OpenMagnetics