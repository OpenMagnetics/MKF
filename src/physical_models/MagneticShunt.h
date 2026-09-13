#pragma once
#include "Defaults.h"
#include "constructive_models/Magnetic.h"
#include <MAS.hpp>
#include <optional>
#include <string>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

// Magnetic shunts (MAS-RFC 0015, ABT #1176): the reluctance ladder that carries the leakage flux
// through a permeable sheet laid across the winding window, after
//   J. Zhang, Z. Ouyang, M. C. Duffy, M. A. E. Andersen, W. G. Hurley, "Leakage Inductance
//   Calculation for Planar Transformers with a Magnetic Shunt", IEEE Trans. Ind. Appl. 50(6), 2014,
//   eqs. (10)-(16), and
//   M. Li, Z. Ouyang, M. A. E. Andersen, "High Frequency LLC Resonant Converter with Magnetic Shunt
//   Integrated Planar Transformer", IEEE Trans. Power Electron. 34(3), 2019 (DTU Orbit 2018),
//   eqs. (13), (17), (27).
//
// Topology, per winding window the shunt crosses (node A = inner column at the shunt level, node B =
// outer column at the shunt level):
//   - the shunt branch A-B: the in-window span of the sheet, R_s2 = span / (mu0 mus t d), in series
//     with the air gaps to the columns when the sheet stops short of them (gapToColumns,
//     R_g = g / (mu0 (t + g)(d + g)), Li eq. (13) fringing), and with the gaps of a segmented sheet;
//   - the upper core loop A-B, driven by the ampere-turns above the sheet: the core path
//     (Zhang eq. (10)-(11): R_c1/2 + R_c2 with R_c1 = R_c2 = l_e / (2 mu0 mur A_e)) plus, where the
//     sheet runs into a column through that column's gap, the air left above it
//     (R = l_a / (mu0 (w + l_a)(d + l_a)), l_a = (g_core - t) / 2) and half the sheet thickness
//     (Zhang R_s1 = (t/2) / (mu0 mus w d));
//   - the lower core loop, identical, driven by the (equal and opposite) ampere-turns below.
// phi_s = F / (R_branch + R_up || R_low); the network stores E = F^2 / (2 R_w), so the window
// contributes L = F^2 / R_w per ampere of source current. With no gaps and a stiff core this reduces
// to Zhang's dominant term L = 2 n^2 mu0 mus t l_w / b_w for the two windows of an E core.
//
// Only rectangular-window cores with a central column and two rectangular lateral columns (E, planar
// E, EFD-type) are modelled; anything else throws.

struct MagneticShuntWindowCrossing {
    // +1 for the window on the +x side of the main column, -1 for the mirrored one.
    int side = 0;
    // Length of sheet lying inside the window along x (m).
    double spanInWindow = 0;
    // Depth of sheet material across the window (m).
    double depth = 0;
    // Width of the window the sheet crosses (m).
    double windowWidth = 0;
    double shuntBranchReluctance = 0;
    double upperLoopReluctance = 0;
    double lowerLoopReluctance = 0;
    // R_w = R_branch + R_up || R_low.
    double totalReluctance = 0;
    // Column-overlap areas per side (m^2), zero when the sheet stops short of that column.
    double innerOverlapArea = 0;
    double outerOverlapArea = 0;
};

// Geometry of one shunt as the network sees it: everything except the sheet permeability.
struct MagneticShuntNetworkInputs {
    std::string name;
    MagneticShuntPlacement placement = MagneticShuntPlacement::BETWEEN_SECTIONS;
    double thickness = 0;          // axial height of the sheet (m)
    double depth = 0;              // depth across the window (m)
    double axialCoordinate = 0;    // y of the sheet centre (m)
    double coreLoopReluctance = 0; // Zhang R_c1/2 + R_c2 of one half-core loop (A/Wb)
    struct Crossing {
        int side = 0;
        double spanInWindow = 0;
        double windowWidth = 0;
        double depthInWindow = 0;
        // Extent of the sheet in side coordinates (distance from the main column axis, m).
        double sheetStart = 0;
        double sheetEnd = 0;
        // Side gaps (shunt branch), used when the sheet stops short of a column.
        std::optional<double> innerSideGap;
        std::optional<double> outerSideGap;
        // Column overlaps (core loops), used when the sheet enters a column through its gap.
        double innerOverlapWidth = 0;
        double innerOverlapDepth = 0;
        double innerOverlapAirPerSide = 0;
        double outerOverlapWidth = 0;
        double outerOverlapDepth = 0;
        double outerOverlapAirPerSide = 0;
        // Segmented sheet: total piece length inside the window and the gaps between pieces.
        std::vector<double> segmentGaps;
    };
    std::vector<Crossing> crossings;
};

struct MagneticShuntContribution {
    std::string name;
    double relativePermeability = 0;
    // Net ampere-turns enclosed above the sheet per ampere (peak) of source current.
    double enclosedMagnetomotiveForcePerAmpere = 0;
    std::vector<MagneticShuntWindowCrossing> crossings;
    // Inductance stored by the shunt network (sheet, its gaps, the core loops), per the source winding.
    double networkLeakageInductance = 0;
    // Inductance the Energy method attributes to the air the sheet displaces (already counted there).
    double displacedAirLeakageInductance = 0;
    // Peak flux density in the in-window part of the sheet, per ampere of source peak current (T/A).
    double magneticFluxDensityPerAmpere = 0;
    // Peak field strength in the in-window part of the sheet, per ampere of source peak current.
    double magneticFieldStrengthPerAmpere = 0;
    // Sheet losses per squared ampere of source peak current (W/A^2): 0.5 omega mu0 mu'' H0^2 V
    // (Li eq. (27)). Present only when the material publishes a complex permeability whose measured
    // span covers the frequency; the data are not extrapolated.
    std::optional<double> lossesPerAmpereSquared;
};

struct MagneticShuntLeakageResult {
    // Energy-method (air) leakage of the windings, computed without the shunts.
    double windingLeakageInductance = 0;
    double leakageInductance = 0;
    std::vector<MagneticShuntContribution> shunts;
};

class MagneticShuntModel {
  public:
    // Resolves the sheet material (inline record or database name). Throws when the name is not in
    // the core-material database or when the record has no initial permeability.
    static CoreMaterial resolve_material(const MagneticShunt& shunt);

    // Relative permeability of the sheet at the given conditions, from the material's initial
    // permeability (MKF InitialPermeability, frequency and temperature modifiers included).
    static double get_relative_permeability(const CoreMaterial& material, double temperature, double frequency);

    // True when the shunt carries leakage flux (placement inWindow or betweenSections).
    static bool is_leakage_shunt(const MagneticShunt& shunt);
    static bool has_leakage_shunts(const Magnetic& magnetic);

    // Throws for placements that no model handles (outsideWindow).
    static void check_supported_placements(const Magnetic& magnetic);

    // Geometry of shunt shuntIndex against the core, validated. Throws on unsupported cores, on a
    // sheet that collides with a column outside that column's gap, on gapToColumns inconsistent with
    // the drawn position, and on segment lengths that do not add up to the sheet width.
    static MagneticShuntNetworkInputs extract_network_inputs(Magnetic& magnetic, size_t shuntIndex, double temperature, double frequency);

    // Net ampere-turns above the sheet per ampere of source peak current, for the leakage excitation
    // (source +1 A, destination balancing ampere-turns). Throws when a turn intersects the sheet.
    static double calculate_enclosed_magnetomotive_force(Magnetic& magnetic, const MagneticShuntNetworkInputs& inputs, size_t sourceIndex, size_t destinationIndex);

    static MagneticShuntContribution solve_network(const MagneticShuntNetworkInputs& inputs, double relativePermeability, double enclosedMagnetomotiveForcePerAmpere);

    // Adds 0.5 omega mu0 mu'' H^2 V when the material's complex permeability covers the frequency.
    static void add_losses(MagneticShuntContribution& contribution, const MagneticShuntNetworkInputs& inputs, const CoreMaterial& material, double frequency);

    // Core with every gap that holds a shunt sheet reduced to its air: g -> g - t + t / mus (a
    // permeable slab in series with the remaining air). Gaps without a sheet are unchanged.
    // relativePermeabilityPerShunt as in LeakageInductance::calculate_shunt_leakage.
    static Core apply_shunts_to_gapping(Magnetic magnetic, double temperature, double frequency,
                                        std::optional<std::vector<double>> relativePermeabilityPerShunt = std::nullopt);

    // Shunt method total for a given winding (air) leakage: windingLeakageInductance minus the air every
    // in-window sheet displaces plus every sheet network's inductance. See
    // LeakageInductance::calculate_shunt_leakage for relativePermeabilityPerShunt.
    static MagneticShuntLeakageResult assemble_leakage(Magnetic& magnetic, double windingLeakageInductance, double frequency,
                                                       size_t sourceIndex, size_t destinationIndex,
                                                       std::optional<std::vector<double>> relativePermeabilityPerShunt = std::nullopt);

    // Sizing helper (used by no adviser unless Settings::get_coil_adviser_size_magnetic_shunts() is on).
    // Places one betweenSections sheet of `material` in each winding window, in the axial clearance
    // between the source and destination windings' conduction sections that face each other, spanning
    // the window with the given gaps to the columns and the core's depth, and solves its thickness so
    // that the Shunt leakage method returns targetLeakageInductance. Returns the two sheets (+x and -x
    // windows). Throws when the sections are not stacked along the column axis, when the target is not
    // above the shunt-less leakage, when the full clearance cannot reach it, or when the sheets' peak
    // flux density at sourceCurrentPeak is not below 0.7 of the material saturation.
    static std::vector<MagneticShunt> size_shunt_for_leakage(Magnetic magnetic, double targetLeakageInductance, const CoreMaterial& material,
                                                             double innerGapToColumn, double outerGapToColumn, double sourceCurrentPeak,
                                                             double frequency, size_t sourceIndex = 0, size_t destinationIndex = 1);
};

} // namespace OpenMagnetics
