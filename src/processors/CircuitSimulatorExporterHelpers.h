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
    double Isat;           // Saturation current [A] — from Magnetic::calculate_saturation_current
                           // (gapped B_sat·N·A_e/L). NOT recomputed here: the previous local
                           // formula Bsat/(mu0*mu_r)*le/N used the UNGAPPED core path and
                           // under-reported I_sat by the gap factor (~10×), collapsing the
                           // exported saturable inductor above the bogus value (ABT #40).
    bool valid;            // Whether parameters are valid

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

// PD-safe behavioural realization of the mutual (cross-coupling) resistance loss for
// n>=3 windings (abt #50/#72): grounded uncoupled ladders + sense/behavioural sources,
// adding NO magnetically coupled inductors, so the Lmag coupled-L matrix stays
// positive-definite. Requires the winding series path to be routed through the
// Node_Wtop_<k> nodes (see export_magnetic_as_subcircuit, LADDER mode).
std::string emit_mutual_resistance_behavioural_spice(
    const std::vector<CircuitSimulatorExporter::MutualResistanceCoefficients>& mutualCoeffs,
    size_t numWindings,
    double frequency);

} // namespace OpenMagnetics
