#pragma once

#include "processors/CircuitSimulatorInterface.h"

#include <string>
#include <vector>

namespace OpenMagnetics {

// Saturation model parameters for circuit simulation
struct SaturationParameters {
    double Bsat;           // Saturation flux density [T]
    double mu_r;           // Relative permeability (initial)
    double Ae;             // Effective core area [m²]
    double le;             // Effective magnetic path length [m]
    double primaryTurns;   // Number of primary turns
    double Lmag;           // Magnetizing inductance [H]
    bool valid;            // Whether parameters are valid

    double Isat() const {
        // I_sat = lambda_sat / L = B_sat·N·A_e / L_actual — the SAME identity the authoritative
        // Magnetic::calculate_saturation_current uses. The old H-field form (Hsat·le/N) back-solved through
        // the INITIAL-permeability BARE-CORE inductance (mu0·mu_r·N²·Ae/le), which disagrees with the real
        // gapped/effective L_actual by the turns factor — exporting an I_sat ~N× too LOW. ngspice's
        // saturating-L source L0/(1+(I/Isat)²) then collapsed and aborted ('timestep too small') on
        // perfectly-valid cores (e.g. a 5.2× FALSE saturation on a core whose true I_sat is 1.75× the
        // current). Lmag here is the very calculate_inductance_from_number_turns_and_gapping value that
        // calculate_saturation_current uses, so the exported subcircuit now agrees with it. (ABT #33.)
        return fluxLinkageSat() / Lmag;
    }

    double fluxLinkageSat() const {
        return primaryTurns * Ae * Bsat;
    }
};

// Extract saturation parameters from a magnetic component (shared across vendors).
SaturationParameters get_saturation_parameters(const Magnetic& magnetic, double temperature);

// Shared SPICE-like circuit fragment emitters used by both Ngspice and Ltspice exporters.
std::string emit_fracpole_winding_spice(
    const FractionalPoleNetwork& net, size_t windingIndex, const std::string& is,
    double dcResistance);

std::string emit_core_ladder_spice(
    const std::vector<double>& coeffs, size_t numWindings);

std::string emit_core_rosano_spice(
    const std::vector<double>& coeffs, size_t numWindings);

std::string emit_winding_rosano_spice(
    const std::vector<double>& coeffs,
    const std::string& is,
    double dcResistance);

std::string emit_mutual_resistance_network_spice(
    const std::vector<CircuitSimulatorExporter::MutualResistanceCoefficients>& mutualCoeffs,
    double magnetizingInductance,
    size_t numWindings);

} // namespace OpenMagnetics
