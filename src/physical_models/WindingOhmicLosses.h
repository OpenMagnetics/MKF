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

using namespace MAS; // QUAL-001 TODO: qualify types and remove

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
        // Centerline copper length of each winding PARALLEL's terminal/connection routing (provided
        // lead lengths, or the geometric route to the window border, never both, plus the in-winding
        // links across layer boundaries: U stubs/edge runs, Z diagonals, YZ-face dragbacks). The ONE
        // length source for both the DC stage (below) and the skin stage
        // (WindingSkinEffectLosses::calculate_skin_effect_losses), so they can never disagree about
        // what copper exists. Indexed [winding][parallel]. Returns zeros unless real winding
        // geometry is enabled; throws for planar (PCB) constructions (ABT #492 owner ruling).
        static std::vector<std::vector<double>> calculate_connection_length_per_winding_per_parallel(Coil coil);
        // DC resistance contributed by each winding PARALLEL's terminal/connection leads (provided
        // lead lengths plus the derived in-winding routing across layer boundaries). Each parallel of a
        // bifilar/N-filar group is its own conductor with its own leads, in series with that parallel's
        // turns. Indexed [winding][parallel]. Returns zeros unless real winding geometry is enabled.
        static std::vector<std::vector<double>> calculate_connection_resistance_per_winding_per_parallel(Coil coil, double temperature);
        static WindingLossesOutput calculate_ohmic_losses(Coil coil, OperatingPoint operatingPoint, double temperature);
        static double calculate_ohmic_losses_per_meter(Wire wire, SignalDescriptor current, double temperature);
};

} // namespace OpenMagnetics