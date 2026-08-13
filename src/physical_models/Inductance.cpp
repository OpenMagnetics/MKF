#include "physical_models/Inductance.h"
#include "physical_models/ReluctanceNetwork.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/Reluctance.h"
#include "support/Utils.h"
#include "support/Exceptions.h"

#include <cmath>
#include <numbers>

namespace OpenMagnetics {

std::string Inductance::get_winding_name(Magnetic& magnetic, size_t windingIndex) {
    return magnetic.get_coil().get_functional_description()[windingIndex].get_name();
}

MagnetizingInductanceOutput Inductance::calculate_magnetizing_inductance(
    Magnetic magnetic,
    OperatingPoint* operatingPoint) {
    
    ReluctanceModels reluctanceModel;
    from_json(_reluctanceModel, reluctanceModel);
    MagnetizingInductance magnetizingInductanceModel(reluctanceModel);
    
    return magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(
        magnetic, operatingPoint);
}

double Inductance::calculate_magnetizing_inductance_referred_to_winding(
    Magnetic magnetic,
    size_t windingIndex,
    OperatingPoint* operatingPoint) {
    
    auto magnetizingInductanceOutput = calculate_magnetizing_inductance(magnetic, operatingPoint);
    double magnetizingInductancePrimary = magnetizingInductanceOutput.get_magnetizing_inductance().get_nominal().value();
    
    // Magnetizing inductance referred to primary (winding 0)
    // To refer to another winding: Lm_i = Lm_0 * (N_i / N_0)²
    if (windingIndex == 0) {
        return magnetizingInductancePrimary;
    }
    
    double primaryTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();
    double targetTurns = magnetic.get_coil().get_functional_description()[windingIndex].get_number_turns();
    double turnsRatio = targetTurns / primaryTurns;
    
    return magnetizingInductancePrimary * turnsRatio * turnsRatio;
}

double Inductance::calculate_leakage_inductance(
    Magnetic magnetic,
    size_t sourceIndex,
    size_t destinationIndex,
    double frequency) {
    
    if (sourceIndex == destinationIndex) {
        return 0.0;
    }
    
    LeakageInductance leakageInductanceModel;
    auto leakageOutput = leakageInductanceModel.calculate_leakage_inductance(
        magnetic, frequency, sourceIndex, destinationIndex);
    
    return leakageOutput.get_leakage_inductance_per_winding()[0].get_nominal().value();
}

// ABT #396: see the declaration in Inductance.h for the full rationale. On a 3-column
// E 42 with the primary on the centre leg and the secondary on an outer leg,
// calculate_coupling_coefficient said 0.999 where the matrix said 0.624 before this was
// unified; the latter is the real flux divider, since the secondary only links the share
// of primary flux that returns through ITS leg.
std::optional<std::vector<std::vector<double>>> Inductance::magnetizing_coupling_matrix(
    Magnetic& magnetic,
    const MagnetizingInductanceOutput& magnetizingOutput) {

    if (!ReluctanceNetwork::has_non_main_placement(magnetic)) {
        return std::nullopt;
    }
    ReluctanceNetwork magneticCircuit(
        magnetic.get_core(),
        magnetizingOutput.get_ungapped_core_reluctance().value(),
        magnetizingOutput.get_reluctance_per_gap().value_or(std::vector<AirGapReluctanceOutput>{}));
    return magneticCircuit.calculate_magnetizing_inductance_matrix(magnetic);
}

double Inductance::calculate_mutual_inductance(
    Magnetic magnetic,
    size_t sourceIndex,
    size_t destinationIndex,
    OperatingPoint* operatingPoint) {
    
    if (sourceIndex == destinationIndex) {
        // Self inductance, not mutual
        throw std::invalid_argument("Cannot calculate mutual inductance between a winding and itself");
    }

    // Leg-separated windings do not couple ideally: ask the reluctance network, the same
    // source calculate_inductance_matrix uses (ABT #396). The sign is meaningful and is
    // kept — the network returns a negative mutual in the common branch orientation, where
    // flux up one leg comes down the other.
    auto magnetizingOutput = calculate_magnetizing_inductance(magnetic, operatingPoint);
    if (auto networkMatrix = magnetizing_coupling_matrix(magnetic, magnetizingOutput)) {
        return (*networkMatrix)[sourceIndex][destinationIndex];
    }

    // All windings share the main column, so the closed form is exact: with a single
    // magnetizing flux linking every winding, M = sqrt(Lm_source * Lm_dest), equivalently
    // Lm_primary * (N_dest / N_source) for a two-winding transformer.
    double Lm_source = calculate_magnetizing_inductance_referred_to_winding(magnetic, sourceIndex, operatingPoint);
    double Lm_dest = calculate_magnetizing_inductance_referred_to_winding(magnetic, destinationIndex, operatingPoint);
    return std::sqrt(Lm_source * Lm_dest);
}

double Inductance::calculate_self_inductance(
    Magnetic magnetic,
    size_t windingIndex,
    double frequency,
    OperatingPoint* operatingPoint) {
    
    // Self inductance L_ii = Lm_i + Ll_i
    // Where Ll_i is the total leakage inductance as seen from winding i
    
    // ABT #396: the magnetizing part of L_ii has to come from the same place the mutual
    // does, or the coupling coefficient built from them is a ratio of two different
    // models. For leg-separated windings the driving-point magnetizing inductance of
    // winding i is the network's diagonal — its own column in series with the parallel
    // combination of the others — not the rank-1 value referred from the primary, which
    // assumes every winding sees the main column's flux path.
    auto magnetizingOutput = calculate_magnetizing_inductance(magnetic, operatingPoint);
    double Lm_i;
    if (auto networkMatrix = magnetizing_coupling_matrix(magnetic, magnetizingOutput)) {
        Lm_i = (*networkMatrix)[windingIndex][windingIndex];
    }
    else {
        Lm_i = calculate_magnetizing_inductance_referred_to_winding(magnetic, windingIndex, operatingPoint);
    }

    // Self-leakage of winding i is the diagonal of the energy-method leakage
    // matrix (Λ_ii = 4·W(e_i)), the single source of leakage shared with
    // calculate_inductance_matrix. The previous code summed "the maximum
    // leakage across the other windings" but an unconditional break made the
    // loop dead (only winding 0/1 ever used), and it added the full
    // short-circuit pair leakage on the diagonal — which, combined with the
    // ideal-k mutual term, double-counted leakage (ABT #104).
    size_t numWindings = magnetic.get_coil().get_functional_description().size();
    double selfLeakage = 0.0;
    if (numWindings > 1) {
        LeakageInductance leakageModel;
        auto leakageMatrix = leakageModel.calculate_leakage_inductance_matrix(magnetic, frequency);
        selfLeakage = leakageMatrix[windingIndex][windingIndex];
    }

    return Lm_i + selfLeakage;
}

double Inductance::calculate_coupling_coefficient(
    Magnetic magnetic,
    size_t sourceIndex,
    size_t destinationIndex,
    double frequency,
    OperatingPoint* operatingPoint) {
    
    if (sourceIndex == destinationIndex) {
        return 1.0; // Perfect coupling with itself
    }
    
    double L11 = calculate_self_inductance(magnetic, sourceIndex, frequency, operatingPoint);
    double L22 = calculate_self_inductance(magnetic, destinationIndex, frequency, operatingPoint);
    double M = calculate_mutual_inductance(magnetic, sourceIndex, destinationIndex, operatingPoint);
    
    // k = M / sqrt(L11 * L22)
    double denominator = std::sqrt(L11 * L22);
    if (denominator < 1e-15) {
        return 0.0;
    }
    
    double k = M / denominator;

    // ABT #396: the clamp that used to sit here, min(1, max(0, k)), silently destroyed two
    // different things. It turned the legitimately NEGATIVE coupling of leg-separated
    // windings — flux up one leg comes down the other, which the reluctance network reports
    // with a sign — into a flat 0, i.e. "these windings do not couple at all". And it
    // capped |k| > 1, which is not a value to be tidied away but a contradiction: a mutual
    // inductance exceeding sqrt(L11*L22) violates the energy bound, so it means one of the
    // three inputs is wrong. Report the sign, and refuse the impossible.
    constexpr double couplingBoundTolerance = 1e-6;
    if (std::abs(k) > 1 + couplingBoundTolerance) {
        throw InvalidInputException(ErrorCode::CALCULATION_INVALID_RESULT,
            "Coupling coefficient between windings " + std::to_string(sourceIndex) + " and " +
            std::to_string(destinationIndex) + " is " + std::to_string(k) +
            ", which exceeds the energy bound |k| <= 1: the mutual inductance " +
            std::to_string(M) + " H is larger than sqrt(L11*L22) = " + std::to_string(denominator) +
            " H (L11 = " + std::to_string(L11) + " H, L22 = " + std::to_string(L22) + " H)");
    }
    return k;
}

ScalarMatrixAtFrequency Inductance::calculate_leakage_inductance_matrix(
    Magnetic magnetic,
    double frequency) {
    
    auto& functionalDescription = magnetic.get_coil().get_functional_description();
    size_t numWindings = functionalDescription.size();
    
    if (numWindings == 0) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS,
            "Cannot calculate leakage inductance matrix: no windings defined");
    }
    
    ScalarMatrixAtFrequency result;
    result.set_frequency(frequency);
    
    std::map<std::string, std::map<std::string, DimensionWithTolerance>> magnitude;
    
    LeakageInductance leakageInductanceModel;
    
    // Row i: leakage inductances referred to winding i
    for (size_t i = 0; i < numWindings; ++i) {
        std::string windingName_i = get_winding_name(magnetic, i);
        
        auto leakageOutput = leakageInductanceModel.calculate_leakage_inductance_all_windings(
            magnetic, frequency, i);
        auto leakagePerWinding = leakageOutput.get_leakage_inductance_per_winding();
        
        for (size_t j = 0; j < numWindings; ++j) {
            std::string windingName_j = get_winding_name(magnetic, j);
            
            DimensionWithTolerance value;
            if (i == j) {
                // By definition: no leakage from a winding into itself
                value.set_nominal(0.0);
            } else if (j < leakagePerWinding.size() && leakagePerWinding[j].get_nominal()) {
                value.set_nominal(leakagePerWinding[j].get_nominal().value());
            } else {
                value.set_nominal(0.0);
            }
            
            magnitude[windingName_i][windingName_j] = value;
        }
    }
    
    result.set_magnitude(magnitude);
    return result;
}

ScalarMatrixAtFrequency Inductance::calculate_inductance_matrix(
    Magnetic magnetic,
    double frequency,
    OperatingPoint* operatingPoint) {
    
    auto& functionalDescription = magnetic.get_coil().get_functional_description();
    size_t numWindings = functionalDescription.size();
    
    if (numWindings == 0) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS,
            "Cannot calculate inductance matrix: no windings defined");
    }
    // The magnetizing-inductance/reluctance path below reads the core's processed
    // description and resolved shape. On an unprocessed magnetic (shape/material
    // still name strings, no effective parameters) this otherwise fails deep with
    // a cryptic "std::get: wrong index for variant". Surface it clearly instead.
    if (!magnetic.get_core().get_processed_description()) {
        throw CoreNotProcessedException(
            "Cannot calculate inductance matrix: the core has no processed description "
            "(effective parameters/shape unresolved). Run magnetic autocomplete / process the core first.");
    }

    ScalarMatrixAtFrequency result;
    result.set_frequency(frequency);
    
    std::map<std::string, std::map<std::string, DimensionWithTolerance>> magnitude;
    
    // Calculate magnetizing inductance (referred to primary)
    auto magnetizingOutput = calculate_magnetizing_inductance(magnetic, operatingPoint);
    double Lm_primary = magnetizingOutput.get_magnetizing_inductance().get_nominal().value();
    double N_primary = functionalDescription[0].get_number_turns();

    // Full inductance matrix = rank-1 ideal-coupling magnetizing term + energy-method
    // leakage matrix, i.e. L_ij = Lm·(Ni/Np)(Nj/Np) + Λ_ij.
    //
    // ABT #104: the previous hand-rolled matrix put the FULL short-circuit pair
    // leakage on both diagonals (L_ii = Lm_i + Ll_i) while the off-diagonal
    // assumed perfect coupling (M_ij = Lm·(Ni/Np)(Nj/Np), k=1). Those two
    // choices double-count the leakage: the short-circuit inductance
    // L11 - M²/L22 came out ≈ 2·Ll (twice the leakage the model itself
    // computes). The energy-method Λ from LeakageInductance already carries the
    // correct per-winding self-leakage (diagonal) and mutual leakage
    // (off-diagonal), so L = M_mag + Λ reproduces the standalone short-circuit
    // leakage (validated to <0.02% on the ETD 39 40:20 case). Delegating to the
    // model class also removes the parallel hand-rolled magnetics math.
    // A single winding has no leakage (nothing to leak against): keep it free of
    // the energy-method field simulation so single-winding inductors — the hot
    // path through WindingLosses — cost exactly what they did before.
    std::vector<std::vector<double>> leakageMatrix;
    if (numWindings > 1) {
        LeakageInductance leakageModel;
        leakageMatrix = leakageModel.calculate_leakage_inductance_matrix(magnetic, frequency);
    }

    // Multi-column winding placement: replace the rank-1 ideal-coupling magnetizing
    // term with the per-column reluctance network matrix. The network reproduces the
    // rank-1 values exactly when every winding shares the main column, and the
    // imperfect coupling of leg-separated windings when they don't. The energy-method
    // window leakage is only meaningful between windings sharing a column;
    // cross-column pairs take their coupling from the network alone until the window
    // field solvers become window-aware.
    auto networkMatrixIfPlaced = magnetizing_coupling_matrix(magnetic, magnetizingOutput);
    bool multiColumnPlacement = networkMatrixIfPlaced.has_value();
    std::vector<std::vector<double>> networkMatrix;
    std::vector<size_t> columnIndexPerWinding;
    if (multiColumnPlacement) {
        networkMatrix = networkMatrixIfPlaced.value();
        columnIndexPerWinding = ReluctanceNetwork::resolve_winding_column_indexes(magnetic);
    }

    // Build the inductance matrix (symmetric, so compute upper triangular + diagonal)
    for (size_t i = 0; i < numWindings; ++i) {
        std::string windingName_i = get_winding_name(magnetic, i);
        double turns_i = functionalDescription[i].get_number_turns();

        for (size_t j = i; j < numWindings; ++j) {
            std::string windingName_j = get_winding_name(magnetic, j);
            double turns_j = functionalDescription[j].get_number_turns();

            double M_mag_ij = multiColumnPlacement ? networkMatrix[i][j]
                                                   : Lm_primary * (turns_i / N_primary) * (turns_j / N_primary);
            double leakage_ij = (i < leakageMatrix.size() && j < leakageMatrix[i].size())
                                    ? leakageMatrix[i][j]
                                    : 0.0;
            if (multiColumnPlacement && columnIndexPerWinding[i] != columnIndexPerWinding[j]) {
                leakage_ij = 0.0;
            }

            DimensionWithTolerance inductanceValue;
            inductanceValue.set_nominal(M_mag_ij + leakage_ij);

            magnitude[windingName_i][windingName_j] = inductanceValue;
            if (i != j) {
                // Symmetric: L_ji = L_ij
                magnitude[windingName_j][windingName_i] = inductanceValue;
            }
        }
    }

    result.set_magnitude(magnitude);
    return result;
}

std::vector<ScalarMatrixAtFrequency> Inductance::calculate_inductance_matrix_per_frequency(
    Magnetic magnetic,
    std::vector<double> frequencies,
    OperatingPoint* operatingPoint) {
    
    std::vector<ScalarMatrixAtFrequency> results;
    results.reserve(frequencies.size());
    
    for (double freq : frequencies) {
        results.push_back(calculate_inductance_matrix(magnetic, freq, operatingPoint));
    }
    
    return results;
}

} // namespace OpenMagnetics
