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
    double saturationCurrent = 0;  // I_sat [A] — from Magnetic::calculate_saturation_current (authoritative)
    bool valid = false;    // Whether parameters are valid

    // The saturation current that parameterises the exported saturable-L model. DELEGATES to the
    // authoritative Magnetic::calculate_saturation_current (computed once in get_saturation_parameters) — it
    // does NOT re-derive the physics here. The previous H-field form (Hsat·le/N, the initial-permeability
    // bare-core inductance) was a PARALLEL formula that drifted ~N× below calculate_saturation_current after
    // the latter was fixed to use the real gapped/effective L — so ngspice's L0/(1+(I/Isat)²) collapsed and
    // aborted on perfectly-valid cores (a false 5.2× saturation on a core whose true margin is 1.75×). One
    // source of truth, so it cannot drift again. (ABT #33.)
    double Isat() const { return saturationCurrent; }

    double fluxLinkageSat() const {
        return primaryTurns * Ae * Bsat;
    }
};

// Extract saturation parameters from a magnetic component (shared across vendors). Takes the magnetic BY
// VALUE because it delegates I_sat to the non-const Magnetic::calculate_saturation_current.
SaturationParameters get_saturation_parameters(Magnetic magnetic, double temperature);

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
