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

    // ABT #837: 2.5D ANGULAR COVERAGE. The image lattice models a window bounded by
    // high-permeability material on all four sides — true for a pot core, but an E/PQ/ETD
    // window is OPEN over most of its azimuth: the winding is a full loop and the core
    // legs sit behind it only for a fraction f of the revolution. Solving the fully-imaged
    // problem everywhere over-confines the field for those families. The engine therefore
    // solves TWICE — with images and in open air — and blends the LOSSES,
    //     P = f * P_covered + (1 - f) * P_open,
    // which is the same angular-sector average OMFEM applies for its thermal model (its
    // [coverage] line reports the same quantity, measured by ray-marching the real 3D
    // solids: 0.522 for PQ 28/20, 1.000 for the P 3.3/2.6 pot core).
    //
    // DEFAULT IS 1.0 (fully imaged) — the blend is available but NOT applied unless a
    // caller asks for it, because measuring it did not support switching it on:
    //   rect5 (PQ 20/16, 5 wide rect turns)  0.65 -> 0.73 of FEM   better
    //   foil beside a 1 mm gap (ETD 34)      0.94 -> 0.80 of FEM   WORSE
    //   12t round (P 3.3/2.6, coverage 1.0)  unchanged at 1.000
    // The foil result is the tell. OMFEM routes a round-post core WITH a functional gap to
    // an AXISYMMETRIC solve (MasMesher.cpp: !hasFunctionalGap && nLateral>=2 -> planar),
    // i.e. a full solid of revolution — its reference has no open sector at all, so
    // blending toward open air can only move away from it. And the cases where blending
    // helped (rect5) are routed to PLANAR CARTESIAN, the very frame PEEC already solves in,
    // where an "open azimuthal sector" has no meaning either. So the rect5 gain is loss
    // added in the right direction for the wrong reason. Left in, defaulted off, because
    // the effect is real for a 3D winding and will matter once a reference exists that
    // resolves it; it must not be tuned against references that cannot see it.
    //
    // std::nullopt = 1.0. Set it explicitly (e.g. to estimate_angular_coverage(core), or to
    // a measured OMFEM [coverage] value) to enable the blend.
    std::optional<double> angularCoverage;

    // The geometric estimator, exposed for testing. Sums each lateral column's angular
    // subtense as seen from the central axis; true pot cores (closed shells) are 1.0 by
    // definition. It is APPROXIMATE — it reads 0.44 where OMFEM's ray march gives 0.52 on
    // PQ 28/20 — because MKF's column model records a bounding depth, not the real curved
    // segment. Prefer an explicit angularCoverage where a measured value exists.
    static double estimate_angular_coverage(const Core& core);

  private:
    Diagnostics _diagnostics;
};

} // namespace OpenMagnetics
