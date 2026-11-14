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

struct ParallelCurrentSplit;  // physical_models/MagneticField.h (ABT #1314)

class WindingOhmicLosses {
    private:
    protected:
    public:
        WindingOhmicLosses() {
        }
        static double calculate_effective_resistance_per_meter(Wire wire, double frequency, double temperature);
        static double calculate_dc_resistance_per_meter(Wire wire, double temperature);
        // A turn's DC resistance, over its conductor's own length (calculate_turn_conductor_length). `toroidal`: the coil is
        // wound on a toroid (its bobbin's winding window is round).
        static double calculate_dc_resistance(Turn turn, const Wire& wire, double temperature, bool toroidal);
        // A turn is a closed loop round its column, and its conductor is a family of parallel loops: the one at a distance u
        // outwards of the turn's centre line is L + 2 pi u long (Steiner's formula for the parallel curves of a convex loop),
        // whatever the column's shape, as the corners of one conductor bend round one centre. A band of the conductor between
        // the offsets u0 < u1 (along that normal, full height h) therefore conducts as a ring,
        //     G = sigma h ln((L + 2 pi u1) / (L + 2 pi u0)) / (2 pi),
        // the same as a straight band of its width over the length 2 pi (u1 - u0) / ln((L + 2 pi u1) / (L + 2 pi u0)), which
        // this returns (L + 2 pi u0 when u1 = u0). It is the length every loss of the band is booked over: its DC resistance
        // and the losses of its currents. The lengths of adjacent bands combine as their conductances do (the logarithms add),
        // so a conductor cut into bands has exactly its own length. A turn of zero length (the real winding can leave one) is
        // no loop: every band of it is 0 long. Throws when the band's inner loop would not be longer than 0.
        static double calculate_band_length(double turnLength, double lowerOffset, double upperOffset);
        // The length a turn's DC resistance (rho length / conducting area) is booked over: its conductor's own length.
        //  - Solid round wire of radius a: integrating 1 / (L + 2 pi u) over the disc gives (L + sqrt(L^2 - (2 pi a)^2)) / 2.
        //  - Rectangular, foil and planar wire in a window: calculate_band_length over its extent along the window's x (the
        //    normal of its loop; a quarter-turn rotation swaps width and height).
        //  - Litz and multi-strand wire: the turn's length. Their strands are transposed, so each runs at every offset in turn.
        //  - Rectangular, foil and planar wire on a toroid: the turn's length (ABT #2384 left them: the cells of a flat
        //    conductor on a toroid are not weighed by their offsets either, so the two stay consistent).
        static double calculate_turn_conductor_length(const Turn& turn, const Wire& wire, bool toroidal);
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
        // ABT #1314: the copper loss that the parallel current split adds over the DC (resistive) split, R_dc sum_h(|I_strand,h|^2 -
        // (divider |I_winding,h|)^2)/2 for every harmonic it solved, charged to the AC (skin) losses of the turns and of the winding's
        // connection leads at that harmonic. Ohmic stays the DC split (Alf, 2026-10-04). Call after the skin effect losses.
        static WindingLossesOutput apply_parallel_current_split(Coil coil, WindingLossesOutput windingLossesOutput, const ParallelCurrentSplit& parallelCurrentSplit, double temperature);
        static double calculate_ohmic_losses_per_meter(Wire wire, SignalDescriptor current, double temperature);
};

} // namespace OpenMagnetics