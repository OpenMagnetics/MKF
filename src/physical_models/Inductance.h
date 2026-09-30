#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/LeakageInductance.h"
#include <MAS.hpp>
#include <vector>
#include <map>
#include <string>

using namespace MAS;

namespace OpenMagnetics {

/**
 * @class Inductance
 * @brief Calculates the complete inductance matrix for multi-winding transformers.
 * 
 * This class computes the inductance matrix [L] that relates terminal voltages and currents
 * according to the transformer equations from Spreen (1990):
 * 
 *   u₁ = (jωL₁₁)i₁ + (jωL₁₂)i₂ + ...
 *   u₂ = (jωL₂₁)i₁ + (jωL₂₂)i₂ + ...
 *   ...
 * 
 * The inductance matrix elements are:
 * - Diagonal elements Lᵢᵢ: Self inductance of winding i (magnetizing + leakage)
 * - Off-diagonal elements Lᵢⱼ: Mutual inductance between windings i and j
 * 
 * All of them come from one energy formulation. The flux of the magnetic splits into the
 * core (magnetizing) flux, with permeance P_m, and the window (leakage) field, whose energy
 * for the winding currents i is W = ½·iᵀΛi (LeakageInductance::calculate_leakage_inductance_matrix).
 * The stored energy ½·iᵀLi is their sum, so
 *
 *   L_ij = N_i·N_j·P_m + Λ_ij          (Lm_i = N_i²·P_m, the magnetizing inductance of winding i)
 *   k_ij = L_ij / √(L_ii·L_jj)
 *
 * and for ampere-turn balanced currents (i_j = −r·i_i, r = N_i/N_j) the magnetizing term
 * cancels, leaving the pairwise leakage referred to winding i:
 *
 *   L_ii + r²·L_jj − 2r·L_ij = Λ_ii + r²·Λ_jj − 2r·Λ_ij
 *
 * (for a 1:1 common-mode choke: L_DM = L11 + L22 − 2·L12). calculate_self_inductance,
 * calculate_mutual_inductance and calculate_coupling_coefficient read the entries of the one
 * matrix calculate_inductance_matrix_values builds, so they cannot disagree with it. For
 * windings on different columns the magnetizing term is the reluctance network's matrix
 * instead of the rank-1 N_i·N_j·P_m (ABT #396).
 */
class Inductance {
private:
    std::string _reluctanceModel;
    
public:
    /**
     * @brief Default constructor using default reluctance model.
     */
    Inductance() {
        _reluctanceModel = to_string(Defaults().reluctanceModelDefault);
    }

    /**
     * @brief Constructor with specified reluctance model.
     * @param model The reluctance model to use for magnetizing inductance calculations.
     */
    Inductance(ReluctanceModels model) {
        _reluctanceModel = to_string(model);
    }

    /**
     * @brief Constructor with string reluctance model.
     * @param model The reluctance model name string.
     */
    Inductance(std::string model) : _reluctanceModel(model) {}

    virtual ~Inductance() = default;

    /**
     * @brief Calculate the complete inductance matrix for a magnetic component.
     * 
     * Computes the inductance matrix at the specified frequency, including:
     * - Self inductances (diagonal elements)
     * - Mutual inductances (off-diagonal elements)
     * 
     * @param magnetic The magnetic component (core + coil).
     * @param frequency Operating frequency in Hz.
     * @param operatingPoint Optional operating point for temperature-dependent calculations.
     * @return ScalarMatrixAtFrequency containing the inductance matrix [L].
     */
    ScalarMatrixAtFrequency calculate_inductance_matrix(
        Magnetic magnetic,
        double frequency,
        OperatingPoint* operatingPoint = nullptr);

    /**
     * @brief Calculate inductance matrices at multiple frequencies.
     * 
     * @param magnetic The magnetic component.
     * @param frequencies Vector of frequencies in Hz.
     * @param operatingPoint Optional operating point.
     * @return Vector of ScalarMatrixAtFrequency, one per frequency.
     */
    std::vector<ScalarMatrixAtFrequency> calculate_inductance_matrix_per_frequency(
        Magnetic magnetic,
        std::vector<double> frequencies,
        OperatingPoint* operatingPoint = nullptr);

    /**
     * @brief Calculate the leakage inductance matrix for a magnetic component.
     * 
     * Builds a matrix [Llk] where each element Llk(i,j) is the leakage inductance
     * between winding i (source) and winding j (destination) referred to winding i.
     *
     * - Diagonal elements are 0 (no leakage from a winding into itself).
     * - Off-diagonal elements are generally NOT symmetric because they are referred
     *   to the source winding.
     *
     * Example for a 3-winding transformer: returns a 3x3 matrix with zeros on the
     * diagonal and Llk(0,1) being the leakage between winding 0 and 1 referred to 0.
     *
     * @param magnetic The magnetic component (core + coil).
     * @param frequency Operating frequency in Hz.
     * @return ScalarMatrixAtFrequency containing the leakage inductance matrix [Llk].
     */
    ScalarMatrixAtFrequency calculate_leakage_inductance_matrix(
        Magnetic magnetic,
        double frequency);


    /**
     * @brief Numeric inductance matrix L = M_mag + Λ (henries), windings in coil order.
     *
     * The single source of truth for every self inductance, mutual inductance and coupling
     * coefficient this class reports; calculate_inductance_matrix is this matrix keyed by
     * winding name.
     */
    std::vector<std::vector<double>> calculate_inductance_matrix_values(
        Magnetic magnetic,
        double frequency,
        OperatingPoint* operatingPoint = nullptr);

    /**
     * @brief Calculate the mutual inductance between two windings.
     *
     * The off-diagonal entry of the inductance matrix: the magnetizing mutual
     * N_i·N_j·P_m = √(Lm_i·Lm_j) plus the mutual leakage Λ_ij of the window field, which is
     * why it needs the frequency.
     *
     * @param magnetic The magnetic component.
     * @param sourceIndex Index of the first winding.
     * @param destinationIndex Index of the second winding.
     * @param frequency Frequency for the leakage (window field) calculation.
     * @param operatingPoint Optional operating point.
     * @return Mutual inductance value in Henries.
     */
    double calculate_mutual_inductance(
        Magnetic magnetic,
        size_t sourceIndex,
        size_t destinationIndex,
        double frequency,
        OperatingPoint* operatingPoint = nullptr);

    /**
     * @brief Calculate the self inductance of a winding.
     * 
     * The diagonal of the inductance matrix: the magnetizing inductance of the winding
     * plus its self-leakage:
     *   Lᵢᵢ = Lₘᵢ + Λᵢᵢ
     * 
     * @param magnetic The magnetic component.
     * @param windingIndex Index of the winding.
     * @param frequency Frequency for leakage inductance calculation.
     * @param operatingPoint Optional operating point.
     * @return Self inductance value in Henries.
     */
    double calculate_self_inductance(
        Magnetic magnetic,
        size_t windingIndex,
        double frequency,
        OperatingPoint* operatingPoint = nullptr);

    /**
     * @brief Calculate the coupling coefficient between two windings.
     * 
     * The coupling coefficient k is defined as:
     *   k = L₁₂ / √(L₁₁ · L₂₂)
     * with all three entries taken from the same inductance matrix, so the mutual carries
     * the mutual leakage Λ₁₂ as well as the magnetizing term. It is signed: leg-separated
     * windings can couple negatively. |k| > 1 throws.
     * 
     * @param magnetic The magnetic component.
     * @param sourceIndex Index of the first winding.
     * @param destinationIndex Index of the second winding.
     * @param frequency Frequency for calculations.
     * @param operatingPoint Optional operating point.
     * @return Coupling coefficient (dimensionless, |k| ≤ 1).
     */
    double calculate_coupling_coefficient(
        Magnetic magnetic,
        size_t sourceIndex,
        size_t destinationIndex,
        double frequency,
        OperatingPoint* operatingPoint = nullptr);

    /**
     * @brief Calculate the magnetizing inductance referred to a specific winding.
     * 
     * The magnetizing inductance scales with the square of turns:
     *   Lₘᵢ = Lₘ_ref · (Nᵢ / N_ref)²
     * 
     * @param magnetic The magnetic component.
     * @param windingIndex Index of the winding to refer to.
     * @param operatingPoint Optional operating point.
     * @return Magnetizing inductance in Henries referred to the specified winding.
     */
    double calculate_magnetizing_inductance_referred_to_winding(
        Magnetic magnetic,
        size_t windingIndex,
        OperatingPoint* operatingPoint = nullptr);

    /**
     * @brief Calculate the leakage inductance between two windings.
     * 
     * Returns the leakage inductance as seen from the source winding
     * when the destination winding carries opposing ampere-turns.
     * 
     * @param magnetic The magnetic component.
     * @param sourceIndex Index of the source winding.
     * @param destinationIndex Index of the destination winding.
     * @param frequency Frequency for calculation.
     * @return Leakage inductance in Henries.
     */
    double calculate_leakage_inductance(
        Magnetic magnetic,
        size_t sourceIndex,
        size_t destinationIndex,
        double frequency);

    /**
     * @brief ABT #396: the ONE place that decides magnetizing coupling for a magnetic
     * whose windings may not all share the main column. Returns the per-column
     * reluctance network's magnetizing inductance matrix when any winding sits off the
     * main column, and nullopt when they all share it (in which case the rank-1
     * sqrt(Lm_i*Lm_j)/turns-ratio closed form used by calculate_inductance_matrix_values
     * and every other consumer -- including
     * ExtendedCantilever::calculate_inductance_matrix, ABT #227.5 -- is already exact and
     * not worth paying the field solve for). Every consumer that assembles a coupling or
     * inductance matrix must call this rather than assuming rank-1 coupling, so a magnetic
     * cannot report two different couplings depending on which entry point was asked.
     *
     * @param magnetic The magnetic component (core + coil).
     * @param magnetizingOutput Output of calculate_inductance_from_number_turns_and_gapping,
     * supplying the ungapped core reluctance and per-gap reluctances the network needs.
     * @return The magnetizing inductance matrix, or nullopt when every winding is on the
     * main column and the caller's closed form already applies.
     */
    static std::optional<std::vector<std::vector<double>>> magnetizing_coupling_matrix(
        Magnetic& magnetic,
        const MagnetizingInductanceOutput& magnetizingOutput);

private:
    /**
     * @brief Calculate magnetizing inductance for the primary winding.
     */
    MagnetizingInductanceOutput calculate_magnetizing_inductance(
        Magnetic magnetic,
        OperatingPoint* operatingPoint);

    /**
     * @brief Throw when a winding index is out of range for the magnetic.
     */
    static void check_winding_index(Magnetic& magnetic, size_t windingIndex);

    /**
     * @brief Get winding name from index for matrix keys.
     */
    std::string get_winding_name(Magnetic& magnetic, size_t windingIndex);
};

} // namespace OpenMagnetics
