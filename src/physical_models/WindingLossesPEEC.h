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
//   - DC / low frequency (<=~10 kHz): ~0.96-1.06  -> matches (scale is correct).
//   - mid/high frequency: MESH-DEPENDENT and NOT YET CONVERGED:
//       * uniform mesh:  ~1.04 @100 kHz, falling to ~0.34-0.65 @>=500 kHz (under);
//       * graded mesh:   ~2.09 @100 kHz (over).
// The true answer is bracketed but the result still depends strongly on the
// mesh -> the partial-inductance formulation is not yet convergent. The graded
// (surface-refined) mesh fixed the ill-conditioning of the fine uniform mesh
// and captures more of the high-f crowding, but over/under-shoots depending on
// resolution.
//
// TO FINISH (real EM-solver work, NOT a constant tweak): make the partial
// inductances convergent — correct self/mutual 2D partial-inductance integrals
// (exact filament GMD / area integrals instead of the point-log approximation),
// a consistent reference handling, and a mesh-convergence study. Consider a
// vetted PEEC formulation/library. Until it converges to the analytical models
// (and an independent 2D FEM) across frequency, this is NOT production-ready and
// stays opt-in behind the analytical default.
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

    // Mesh controls (exposed for convergence studies).
    double _cellMinFactor = 3.0;     // surface cell size = skinDepth / factor
    double _gradingRatio = 1.5;      // cell growth ratio toward the centre
    size_t _maxCellsPerSide = 12;    // per-dimension cell cap

    // Highest harmonic loss must exceed this fraction of the fundamental's loss
    // contribution to be included (mirrors the analytical amplitude threshold).
    double _harmonicAmplitudeThreshold = 0.005;

    WindingLossesPEEC() = default;

    // Drop-in alternate of WindingLosses::calculate_losses(): returns the same
    // WindingLossesOutput populated with total winding losses (and per-winding).
    WindingLossesOutput calculate_losses(Magnetic magnetic, OperatingPoint operatingPoint, double temperature);
};

} // namespace OpenMagnetics
