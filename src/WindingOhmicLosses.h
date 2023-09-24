#pragma once
#include "MAS.hpp"
#include "CoilWrapper.h"
#include "WireWrapper.h"
#include "Utils.h"

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
        double get_dc_resistance(Turn turn, WireWrapper wire, double temperature);
        WindingLossesOutput get_ohmic_losses(CoilWrapper winding, OperatingPoint operatingPoint, double temperature);
};

} // namespace OpenMagnetics