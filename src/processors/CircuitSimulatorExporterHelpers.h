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
        // Saturation current of the as-designed (wound + gapped) magnetic:
        //
        //   I_sat = B_sat · N · A_e / L_actual   (== fluxLinkageSat / Lmag)
        //
        // This is the SAME identity Magnetic::calculate_saturation_current uses,
        // with L_actual the gapped-magnetic inductance already resolved into
        // `Lmag` (calculate_inductance_from_number_turns_and_gapping).
        //
        // The previous form `Bsat/(mu0*mu_r) * le / N` used the core path `le`
        // with the *material* initial permeability `mu_r`, i.e. the saturation
        // current of an UNGAPPED core. For a gapped E/RM core the gap dominates
        // reluctance, so that form under-reported I_sat by the gap factor
        // (mu_r/mu_e, ~10×) — collapsing the exported saturable inductor for any
        // converter whose current exceeded the bogus low value (ABT #40).
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
