#pragma once
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include <MAS.hpp>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

/**
 * @class ReluctanceNetwork
 * @brief Per-column reluctance network of a core, with one MMF source per wound column.
 *
 * The core is decomposed into one branch per column between the two yokes. Each branch
 * carries the column's share of the (converged) ungapped core reluctance plus the
 * reluctances of the gaps sitting in that column; each winding drives an MMF source in
 * the branch of the column it is wound around (resolved from the coil placement:
 * winding.windingWindow -> winding window column edge -> column).
 *
 * This generalizes the historical lumped single-loop model: with every winding on the
 * main column the network reproduces exactly L = N² / (R_main + parallel(R_returns)) =
 * N² / R_total, and for windings on other columns it yields the per-column flux split,
 * the full magnetizing inductance matrix, and the (im)perfect coupling that the lumped
 * model cannot express.
 *
 * The ungapped-core share per column is split geometrically (column length plus yoke
 * share over column area) and scaled so that the main-column driving-point combination
 * reproduces the lumped IEC effective-parameter reluctance exactly — the shape-constant
 * integrals remain the ground truth for the central-drive case.
 */
class ReluctanceNetwork {
  private:
    Core _core;
    std::vector<double> _columnReluctances;
    // One yoke segment (top plus bottom plate, in series along the loop) per window
    // between x-adjacent legs. FD-validated as non-negligible (ABT #227): collapsing
    // the yokes to perfect conductors overweights the far leg's return share.
    std::vector<double> _yokeSegmentReluctances;
    // Column indexes sorted by x position; the mesh loops run between adjacent legs.
    std::vector<size_t> _orderedColumnIndexes;
    size_t _mainColumnIndex = 0;

    // Mesh solve of the leg chain: per-ORIGINAL-column-index magnetomotive forces in,
    // per-original-column-index fluxes out (positive along the common leg direction).
    std::vector<double> solve_column_fluxes_for_magnetomotive_forces(const std::vector<double>& columnMagnetomotiveForces) const;

  public:
    /**
     * @param core The processed core.
     * @param ungappedCoreReluctance Lumped ungapped core reluctance at the (converged)
     *        permeability, straight out of ReluctanceModel.
     * @param reluctancePerGap Per-gap reluctances aligned with the core's gapping
     *        order, straight out of ReluctanceModel::get_core_reluctance.
     */
    ReluctanceNetwork(Core core, double ungappedCoreReluctance, std::vector<AirGapReluctanceOutput> reluctancePerGap);

    const std::vector<double>& get_column_reluctances() const { return _columnReluctances; }
    const std::vector<double>& get_yoke_segment_reluctances() const { return _yokeSegmentReluctances; }
    const std::vector<size_t>& get_ordered_column_indexes() const { return _orderedColumnIndexes; }
    size_t get_main_column_index() const { return _mainColumnIndex; }

    /**
     * @brief Column each winding is wound around, resolved from the coil placement.
     * Absent placement fields resolve to the schema default (window 0 / main column).
     */
    static std::vector<size_t> resolve_winding_column_indexes(Magnetic magnetic);

    /**
     * @brief True when any winding resolves to a column other than the main one.
     * Cheap and throw-free when no winding carries a placement field.
     */
    static bool has_non_main_placement(Magnetic magnetic);

    /**
     * @brief Sign of the MMF each winding drives into its column per ampere of its reference
     * current: +1 for a winding wound CLOCKWISE around its column, -1 for COUNTER_CLOCKWISE
     * (ABT #1346). Read from the turns' orientation when the coil is wound (every turn of a
     * winding must carry one, and the same one); for a coil not yet wound, the sense the winder
     * gives it in the window it is placed in (Coil::get_turn_orientation_in_winding_window: a
     * window on the negative-x side is mirrored into place and wound the opposite way).
     */
    static std::vector<double> resolve_winding_polarities(Magnetic magnetic);

    /**
     * @brief Magnetizing inductance matrix L_ij from the network (window leakage not
     * included). L_ij = s_i·s_j·N_i·N_j·(δ_{ci,cj}·G_ci − G_ci·G_cj / ΣG), with G_b the branch
     * admittance of column b, c_w the column of winding w and s_w its polarity
     * (resolve_winding_polarities). The sign of a mutual follows how the two windings are
     * wound: two windings on opposite legs wound in opposite senses (what the winder does for
     * a mirrored window) drive the flux the same way round the core and couple positively;
     * wound in the same sense, they couple negatively (ABT #1346; ABT #396 took the topological
     * sign alone, always negative for leg-separated windings, which 3-D FEM contradicts:
     * M/L11 +0.10876 on an E 42/21/20 with both windings on the lateral legs).
     */
    std::vector<std::vector<double>> calculate_magnetizing_inductance_matrix(Magnetic magnetic) const;

    /**
     * @brief Per-column magnetic flux for the given winding currents (positive along
     * the common branch orientation; the fluxes of all columns sum to zero). Each winding
     * drives s_w·N_w·I_w (resolve_winding_polarities).
     */
    std::vector<double> calculate_column_magnetic_fluxes(Magnetic magnetic, const std::vector<double>& windingCurrents) const;

    /**
     * @brief Per-column magnetic flux density for the given winding currents.
     */
    std::vector<double> calculate_column_magnetic_flux_densities(Magnetic magnetic, const std::vector<double>& windingCurrents) const;
};

} // namespace OpenMagnetics
