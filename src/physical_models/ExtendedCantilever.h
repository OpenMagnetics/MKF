#pragma once
#include "constructive_models/Magnetic.h"
#include <vector>

namespace OpenMagnetics {

/**
 * @brief Erickson & Maksimovic extended-cantilever model of an N-winding transformer.
 *
 * Reference: R. W. Erickson and D. Maksimovic, "A Multiple-Winding Magnetics Model Having
 * Directly Measurable Parameters", Colorado Power Electronics Center.
 *
 * The model represents an N-winding transformer with directly-measurable parameters:
 *   - a single shunt magnetizing inductance L11 (self-inductance of the reference winding),
 *   - (N-1) effective turns ratios n_k (open-circuit voltage ratios; n[0] = 1),
 *   - N(N-1)/2 effective leakage inductances l_jk between every winding pair.
 *
 * All quantities are referred to winding 0 (the paper's "winding 1", for which n = 1). The
 * effective leakage inductances may be negative (the paper observes this for side-by-side
 * windings), so signs must be preserved.
 *
 * This is the simulation-ready parameterization the circuit-simulator export uses for
 * 3+ winding components; it avoids the ill-conditioned full inductance matrix by deriving
 * everything from the well-conditioned leakage matrix Λ plus the magnetizing inductance.
 */
struct ExtendedCantileverModel {
    double magnetizingInductance = 0.0;                            ///< L11 (H), referred to winding 0
    std::vector<double> turnsRatios;                              ///< n_k (n[0] = 1), size N
    std::vector<std::vector<double>> effectiveLeakageInductances; ///< l_jk (H), symmetric N×N, diagonal unused

    size_t number_windings() const { return turnsRatios.size(); }
};

class ExtendedCantilever {
public:
    /**
     * @brief Build the extended-cantilever model from a processed magnetic.
     *
     * Computes the leakage inductance matrix Λ (LeakageInductance), the magnetizing
     * inductance L11 and the turns ratios, then derives the cantilever parameters.
     */
    static ExtendedCantileverModel calculate(Magnetic magnetic, double frequency);

    /**
     * @brief Build the model from a precomputed leakage matrix + magnetizing inductance + turns ratios.
     *
     * @param leakageMatrix Symmetric N×N leakage inductance matrix Λ (henries), per-winding current basis.
     * @param magnetizingInductance L11 (H) referred to winding 0.
     * @param turnsRatios n_k (n[0] must be 1), size N.
     */
    static ExtendedCantileverModel build_from_leakage_matrix(
        const std::vector<std::vector<double>>& leakageMatrix,
        double magnetizingInductance,
        const std::vector<double>& turnsRatios);

    /**
     * @brief Diagonal element b_jj of the inverse inductance matrix B = L⁻¹ (paper eq. 5).
     *
     * b_jj = (1/n_j²) [ Σ_{i≠j} 1/l_ji + 1/l_jj ],  with 1/l_jj = 1/L11 for the reference
     * winding (index 0) and 0 otherwise.
     */
    static double inverse_inductance_diagonal(const ExtendedCantileverModel& model, size_t windingIndex);

    /**
     * @brief Short-circuit output inductance at winding k (all other windings shorted) = 1 / b_kk.
     *
     * This is the measured "output impedance" reported in Table 2 of the paper, and equals the
     * N-port series leakage inductance L_ok for the secondary windings.
     */
    static double short_circuit_output_inductance(const ExtendedCantileverModel& model, size_t windingIndex);
};

} // namespace OpenMagnetics
