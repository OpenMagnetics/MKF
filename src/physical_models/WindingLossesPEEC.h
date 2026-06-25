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
// VALIDATION STATUS (vs the trusted analytical models on rectangular/foil
// configs, PEEC/analytical ratio):
//   - DC / low frequency (<=~10 kHz): 0.96-1.06  -> matches (scale is correct).
//   - mid frequency (~10-100 kHz):    ~1.0-1.36  -> within test tolerance.
//   - high frequency (>=~500 kHz):    ~0.34-0.65 -> UNDER-predicts.
// Root cause of the high-f gap: a UNIFORM filament mesh cannot resolve the
// sub-skin-depth surface crowding without a very fine mesh, and refining the
// uniform mesh makes the dense 2D log-kernel system ill-conditioned (the
// solution collapses). So the model is physically correct at low/mid frequency
// but plateaus below the true loss at high frequency.
//
// TO FINISH (not a constant tweak — real solver work): graded/non-uniform
// meshing refined toward conductor surfaces (resolve skin depth with few
// filaments), and a better-conditioned PEEC formulation (mesh/loop-current or
// regularised partial inductances). Until then this is NOT production-ready and
// must stay opt-in behind the analytical default.
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
