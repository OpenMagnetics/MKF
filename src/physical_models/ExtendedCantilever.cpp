#include "physical_models/ExtendedCantilever.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/Exceptions.h"

#include <cmath>
#include <vector>

namespace OpenMagnetics {

namespace {

// Invert a small dense N×N matrix by Gauss-Jordan elimination with partial pivoting.
// N is the number of windings (<= ~8 in practice), so an explicit O(N³) inverse is fine.
std::vector<std::vector<double>> invert_matrix(const std::vector<std::vector<double>>& matrix) {
    size_t n = matrix.size();
    std::vector<std::vector<double>> a(matrix);
    std::vector<std::vector<double>> inverse(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        inverse[i][i] = 1.0;
    }

    for (size_t column = 0; column < n; ++column) {
        // Partial pivoting: find the row with the largest magnitude in this column.
        size_t pivotRow = column;
        double pivotMagnitude = std::abs(a[column][column]);
        for (size_t row = column + 1; row < n; ++row) {
            if (std::abs(a[row][column]) > pivotMagnitude) {
                pivotMagnitude = std::abs(a[row][column]);
                pivotRow = row;
            }
        }
        if (pivotMagnitude == 0.0) {
            throw InvalidInputException(ErrorCode::CALCULATION_ERROR,
                "Cannot build extended-cantilever model: inductance matrix is singular");
        }
        if (pivotRow != column) {
            std::swap(a[pivotRow], a[column]);
            std::swap(inverse[pivotRow], inverse[column]);
        }

        double pivot = a[column][column];
        for (size_t k = 0; k < n; ++k) {
            a[column][k] /= pivot;
            inverse[column][k] /= pivot;
        }
        for (size_t row = 0; row < n; ++row) {
            if (row == column) {
                continue;
            }
            double factor = a[row][column];
            for (size_t k = 0; k < n; ++k) {
                a[row][k] -= factor * a[column][k];
                inverse[row][k] -= factor * inverse[column][k];
            }
        }
    }

    return inverse;
}

} // anonymous namespace

ExtendedCantileverModel ExtendedCantilever::build_from_leakage_matrix(
    const std::vector<std::vector<double>>& leakageMatrix,
    double magnetizingInductance,
    const std::vector<double>& turnsRatios) {

    size_t numberWindings = turnsRatios.size();
    if (numberWindings == 0) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS,
            "Cannot build extended-cantilever model: no windings");
    }
    if (leakageMatrix.size() != numberWindings) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS,
            "Cannot build extended-cantilever model: leakage matrix size does not match number of windings");
    }
    if (magnetizingInductance <= 0.0) {
        throw InvalidInputException(ErrorCode::CALCULATION_ERROR,
            "Cannot build extended-cantilever model: magnetizing inductance must be positive");
    }

    // Full inductance matrix L = M + Λ, referred to winding 0:
    //   M_jk = L11 · n_j · n_k   (rank-1 ideal magnetizing coupling)
    //   Λ                         (leakage matrix from the energy method)
    // With the real leakage matrix added, L is only mildly conditioned (~L11/leakage) and
    // inverts cleanly — unlike the leakage-free ideal matrix, which is singular.
    std::vector<std::vector<double>> inductanceMatrix(numberWindings, std::vector<double>(numberWindings, 0.0));
    for (size_t j = 0; j < numberWindings; ++j) {
        for (size_t k = 0; k < numberWindings; ++k) {
            inductanceMatrix[j][k] = magnetizingInductance * turnsRatios[j] * turnsRatios[k] + leakageMatrix[j][k];
        }
    }

    auto inverseInductanceMatrix = invert_matrix(inductanceMatrix);

    // Effective leakage inductances from the inverse inductance matrix (paper eq. 4):
    //   l_jk = −1 / (n_j · n_k · b_jk),  j ≠ k
    ExtendedCantileverModel model;
    model.magnetizingInductance = magnetizingInductance;
    model.turnsRatios = turnsRatios;
    model.effectiveLeakageInductances.assign(numberWindings, std::vector<double>(numberWindings, 0.0));
    for (size_t j = 0; j < numberWindings; ++j) {
        for (size_t k = 0; k < numberWindings; ++k) {
            if (j == k) {
                continue;
            }
            double denominator = turnsRatios[j] * turnsRatios[k] * inverseInductanceMatrix[j][k];
            if (denominator == 0.0) {
                throw InvalidInputException(ErrorCode::CALCULATION_ERROR,
                    "Cannot build extended-cantilever model: degenerate inverse-inductance coupling between windings");
            }
            model.effectiveLeakageInductances[j][k] = -1.0 / denominator;
        }
    }

    return model;
}

double ExtendedCantilever::inverse_inductance_diagonal(const ExtendedCantileverModel& model, size_t windingIndex) {
    size_t numberWindings = model.number_windings();
    if (windingIndex >= numberWindings) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS,
            "Cannot evaluate inverse-inductance diagonal: winding index out of range");
    }

    // b_jj = (1/n_j²) [ Σ_{i≠j} 1/l_ji + 1/l_jj ],  with 1/l_jj = 1/L11 for the reference
    // winding (index 0) and 0 otherwise (paper eq. 5).
    double sum = 0.0;
    for (size_t i = 0; i < numberWindings; ++i) {
        if (i == windingIndex) {
            continue;
        }
        sum += 1.0 / model.effectiveLeakageInductances[windingIndex][i];
    }
    if (windingIndex == 0) {
        sum += 1.0 / model.magnetizingInductance;
    }

    double turnsRatio = model.turnsRatios[windingIndex];
    return sum / (turnsRatio * turnsRatio);
}

double ExtendedCantilever::short_circuit_output_inductance(const ExtendedCantileverModel& model, size_t windingIndex) {
    return 1.0 / inverse_inductance_diagonal(model, windingIndex);
}

ExtendedCantileverModel ExtendedCantilever::calculate(Magnetic magnetic, double frequency) {
    auto& functionalDescription = magnetic.get_coil().get_functional_description();
    size_t numberWindings = functionalDescription.size();
    if (numberWindings == 0) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS,
            "Cannot build extended-cantilever model: no windings defined");
    }

    auto leakageMatrix = LeakageInductance().calculate_leakage_inductance_matrix(magnetic, frequency);

    double magnetizingInductance = MagnetizingInductance()
        .calculate_inductance_from_number_turns_and_gapping(magnetic)
        .get_magnetizing_inductance().get_nominal().value();

    // Turns ratios referred to winding 0 (n[0] = 1).
    double referenceTurns = functionalDescription[0].get_number_turns();
    std::vector<double> turnsRatios(numberWindings);
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        turnsRatios[windingIndex] = functionalDescription[windingIndex].get_number_turns() / referenceTurns;
    }

    return build_from_leakage_matrix(leakageMatrix, magnetizingInductance, turnsRatios);
}

} // namespace OpenMagnetics
