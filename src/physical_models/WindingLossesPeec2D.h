#pragma once
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "processors/Inputs.h"

using namespace MAS;

namespace OpenMagnetics {

// =============================================================================
// ABT #836 / #837 / #139: 2D PEEC winding-loss engine.
//
// The analytical pipeline (CoilMesher -> MagneticField -> skin/proximity kernels)
// superposes the fields of IMPOSED turn currents and samples them at passive
// points: no unknown lives on any conductor, so neighbouring turns neither
// screen each other nor redistribute their own current. That single structural
// fact is the foil-next-to-a-gap failure (OMFEM R_ac/R_dc 11-21 where the
// kernels read 1.3-2.0, ABT #836), the packing slope (MKF/FEM 0.66x at 2 turns
// -> 1.33x at 12, ABT #837) and the stacked-parallel over-prediction (ABT #139).
//
// This engine solves the current DISTRIBUTION instead, following
// I. Kovacevic-Badstuebner, R. Burkart, C. Dittli, A. Muesing, J.W. Kolar,
// "A Fast Method for the Calculation of Foil Winding Losses", ECCE Europe 2015
// (2D PEEC; validated -0.12% vs ANSYS Maxwell2D and -6%/+2.7%/+16% vs
// calorimetric measurement on a gapped foil E-core inductor):
//   - every turn's conducting cross-section is subdivided into rectangular
//     cells, each carrying an unknown complex current;
//   - cells couple through closed-form 2D partial inductances (per-metre. The
//     divergent (mu0/2pi)*log(l) term of the 2D limit is a constant added to
//     EVERY matrix entry; under the per-conductor current constraints it is
//     absorbed by the constraint multipliers, so it is simply never computed —
//     the paper's Eq. (8) argument, asserted by a unit test);
//   - the core is the same (2M+1)x(2N+1) method-of-images lattice the
//     CoilMesher uses, with the same (mu - max(|m|,|n|))/(mu + max(|m|,|n|))
//     weights and the same window frames — images add MUTUAL TERMS ONLY, no
//     new unknowns;
//   - each FUNCTIONAL (subtractive/additive) gap is a fictitious conductor at
//     the gap position carrying the per-harmonic MMF across that gap,
//     phi(h) * R_gap, participating in the imaging (Van den Bossche's device;
//     residual mating-surface gaps are excluded, the ABT #832 doctrine);
//   - per retained harmonic, the bordered system [Z, C^T; C, 0] is solved for
//     the cell currents with one constraint row per conductor (turn); loss is
//     sum of 0.5 * R_cell * |I_cell|^2 * turnLength (peak-amplitude
//     convention, consistent with the rest of the loss pipeline).
//
// The MAS output schema wants ohmic + skin + proximity. PEEC produces one
// self-consistent number per turn, so the classical decomposition is recovered
// BY ITS DEFINITION: each turn is also solved ISOLATED (its own cells only, no
// neighbours, no images, no gap) — skin(h) = isolated(h) - dc-equivalent(h),
// proximity(h) = full(h) - isolated(h). Ohmic (including connection copper,
// which this engine does not model) stays with WindingOhmicLosses.
//
// Scope of this first implementation (throws loudly outside it — no fallback):
//   - concentric windows only (toroids keep the analytical path: a round
//     window has no rectangular image frame);
//   - ROUND / RECTANGULAR / FOIL / PLANAR wires (litz keeps the analytical
//     path: solid-equivalent PEEC would misrepresent strand-level behaviour);
//   - single parallel per winding (the constraint machinery generalises to
//     solving the parallel split — a TODO that will REPLACE the DC-only
//     currentDividerPerTurn — but none of the FEM-arbitrated fixtures needs
//     it yet).
// =============================================================================

class WindingLossesPeec2D {
  public:
    struct TurnHarmonicLoss {
        double frequency = 0;
        double fullLoss = 0;       // W, self-consistent solve (skin + proximity together)
        double isolatedLoss = 0;   // W, the turn alone in free space (classical "skin")
        double dcEquivalentLoss = 0; // W, same harmonic current at uniform density
    };
    struct Diagnostics {
        size_t totalCells = 0;
        size_t totalConductors = 0;
        std::vector<std::vector<TurnHarmonicLoss>> perTurnPerHarmonic; // [turn][harmonic]
        std::vector<double> gapMmfPerHarmonicFundamental;              // A, per functional gap
    };

    WindingLossesPeec2D() = default;

    // Same contract as WindingLosses::calculate_losses. Ohmic comes from
    // WindingOhmicLosses; skin/proximity per turn per harmonic come from the
    // PEEC solves, method_used = "Peec2D".
    WindingLossesOutput calculate_losses(Magnetic magnetic, OperatingPoint operatingPoint, double temperature);

    const Diagnostics& get_diagnostics() const { return _diagnostics; }

    // Knobs (deliberately few; every one has a physical rationale).
    size_t maximumCells = 6000;      // hard cap; LOGS when it binds, never silently truncates
    size_t maximumCellsWasm = 1500;  // documented lower cap for the browser build (not wired yet)
    double cellsPerSkinDepth = 3.0;  // across the penetrated (thin) dimension
    size_t minimumCellsThin = 2;
    size_t maximumCellsThin = 12;
    size_t minimumCellsWide = 4;
    size_t maximumCellsWide = 48;

  private:
    Diagnostics _diagnostics;
};

} // namespace OpenMagnetics
