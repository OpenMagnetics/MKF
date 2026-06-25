#pragma once
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"

namespace OpenMagnetics {

using namespace MAS;

// 2D Partial-Element-Equivalent-Circuit (PEEC) winding-loss model.
//
// Computes skin + proximity + ohmic losses TOGETHER by discretising every turn
// cross-section into filaments, building a (R + jwL) impedance matrix and
// solving the constrained current distribution per harmonic. Unlike the
// analytical 1D-per-direction models (Wang/Dowell/Ferreira), this captures the
// true 2D current crowding (edges, corners, and the strongly non-uniform
// air-gap fringing field).
//
// Reuses the existing infrastructure:
//   - CoilMesher's method-of-images (rectangular winding window) to impose the
//     high-permeability ferrite-core boundary condition on the partial
//     inductances (each filament's images contribute to L with the same
//     (mu-n)/(mu+n) attenuation weights used by the analytical field models).
//   - ResistivityModel / WindingSkinEffectLosses for material data.
//
// Status: evaluation prototype (feature/peec-winding-losses). Compiles, runs,
// reuses CoilMesher's core-image boundary condition, and is a drop-in alternate
// of WindingLosses::calculate_losses.
//
// KNOWN LIMITATION (do not use for production yet): the current formulation
// under-captures the frequency-dependent proximity/skin redistribution — on the
// planar FEM cases it returns ~the ohmic baseline, roughly flat vs frequency,
// instead of the expected f-dependent rise. The partial-inductance scaling /
// constrained-solve needs validation against a 2D FEM reference (and the
// log-kernel reference-length handling re-checked) before this is trustworthy.
//
// Scope: non-toroidal winding windows only for now; litz must be routed to the
// analytical model (meshing every strand is impractical — homogenise the bundle
// to enable PEEC).
//
// Reference: I. Kovacevic-Badstuebner et al., "A Fast Method for the
// Calculation of Foil Winding Losses", EPE 2015 (ETH Zurich).
class WindingLossesPEEC {
  public:
    struct Filament {
        double x = 0;          // centroid [m]
        double y = 0;
        double area = 0;       // cross-section area [m^2]
        double width = 0;      // filament cross-section dimensions [m]
        double height = 0;
        size_t turnIndex = 0;
        size_t windingIndex = 0;
        double turnLength = 0; // [m]
    };

    // Number of filaments across each cross-section dimension (so a turn is
    // discretised into _filamentsPerDimension^2 filaments). Higher = more
    // accurate, more expensive (dense solve scales ~ (N_filaments)^3).
    size_t _filamentsPerDimension = 6;

    // Highest harmonic loss must exceed this fraction of the fundamental's loss
    // contribution to be included (mirrors the analytical amplitude threshold).
    double _harmonicAmplitudeThreshold = 0.005;

    WindingLossesPEEC() = default;

    // Drop-in alternate of WindingLosses::calculate_losses(): returns the same
    // WindingLossesOutput populated with total winding losses (and per-winding).
    WindingLossesOutput calculate_losses(Magnetic magnetic, OperatingPoint operatingPoint, double temperature);
};

} // namespace OpenMagnetics
