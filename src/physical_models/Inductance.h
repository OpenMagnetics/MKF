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
 * The self and mutual inductances are derived from:
 * - Magnetizing inductance (Lₘ): Common flux linking all windings
 * - Leakage inductances (Lₗᵢⱼ): Flux that links only specific windings
 * 
 * For a two-winding transformer with turns N₁ and N₂:
 * - L₁₁ = Lₘ₁ + Lₗ₁₂ (self inductance of primary, referred to primary)
 * - L₂₂ = Lₘ₂ + Lₗ₂₁ (self inductance of secondary, referred to secondary)
 * - L₁₂ = L₂₁ = M = k·√(Lₘ₁·Lₘ₂) (mutual inductance)
 * 
 * Where Lₘ₁ and Lₘ₂ are the magnetizing inductances referred to each winding,
 * related by: Lₘ₂ = Lₘ₁·(N₂/N₁)²
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
     * @brief Calculate the mutual inductance between two windings.
     * 
     * The mutual inductance M is derived from the magnetizing inductance:
     *   M = k · √(Lₘ₁ · Lₘ₂)
     * 
     * For ideal transformers with perfect coupling (k=1):
     *   M = Lₘ₁ · (N₂/N₁) = Lₘ₂ · (N₁/N₂)
     * 
     * @param magnetic The magnetic component.
     * @param sourceIndex Index of the first winding.
     * @param destinationIndex Index of the second winding.
     * @param operatingPoint Optional operating point.
     * @return Mutual inductance value in Henries.
     */
    double calculate_mutual_inductance(
        Magnetic magnetic,
        size_t sourceIndex,
        size_t destinationIndex,
        OperatingPoint* operatingPoint = nullptr);

    /**
     * @brief Calculate the self inductance of a winding.
     * 
     * Self inductance includes the magnetizing inductance contribution
     * and the leakage inductance:
     *   Lᵢᵢ = Lₘᵢ + Lₗᵢ
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
     *   k = M / √(L₁₁ · L₂₂)
     * 
     * Where k = 1 for perfect coupling and k < 1 for real transformers.
     * 
     * @param magnetic The magnetic component.
     * @param sourceIndex Index of the first winding.
     * @param destinationIndex Index of the second winding.
     * @param frequency Frequency for calculations.
     * @param operatingPoint Optional operating point.
     * @return Coupling coefficient (dimensionless, 0 < k ≤ 1).
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

private:
    /**
     * @brief Calculate magnetizing inductance for the primary winding.
     */
    MagnetizingInductanceOutput calculate_magnetizing_inductance(
        Magnetic magnetic,
        OperatingPoint* operatingPoint);

    /**
     * @brief Get winding name from index for matrix keys.
     */
    std::string get_winding_name(Magnetic& magnetic, size_t windingIndex);
};

} // namespace OpenMagnetics
