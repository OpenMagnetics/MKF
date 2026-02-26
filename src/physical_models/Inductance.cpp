#include "physical_models/Inductance.h"
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

double Inductance::calculate_mutual_inductance(
    Magnetic magnetic,
    size_t sourceIndex,
    size_t destinationIndex,
    OperatingPoint* operatingPoint) {
    
    if (sourceIndex == destinationIndex) {
        // Self inductance, not mutual
        throw std::invalid_argument("Cannot calculate mutual inductance between a winding and itself");
    }
    
    // Mutual inductance M = k * sqrt(Lm1 * Lm2)
    // For ideal transformer with k ≈ 1:
    // M = Lm_primary * (N_dest / N_source) = sqrt(Lm_source * Lm_dest)
    
    double Lm_source = calculate_magnetizing_inductance_referred_to_winding(magnetic, sourceIndex, operatingPoint);
    double Lm_dest = calculate_magnetizing_inductance_referred_to_winding(magnetic, destinationIndex, operatingPoint);
    
    // For ideal coupling (k = 1), M = sqrt(Lm1 * Lm2)
    // This is equivalent to Lm_primary * (N2/N1) for a two-winding transformer
    return std::sqrt(Lm_source * Lm_dest);
}

double Inductance::calculate_self_inductance(
    Magnetic magnetic,
    size_t windingIndex,
    double frequency,
    OperatingPoint* operatingPoint) {
    
    // Self inductance L_ii = Lm_i + Ll_i
    // Where Ll_i is the total leakage inductance as seen from winding i
    
    double Lm_i = calculate_magnetizing_inductance_referred_to_winding(magnetic, windingIndex, operatingPoint);
    
    // Calculate leakage inductance contribution
    // For a multi-winding transformer, we sum the leakage contributions from all other windings
    size_t numWindings = magnetic.get_coil().get_functional_description().size();
    double totalLeakage = 0.0;
    
    if (numWindings > 1) {
        // Use leakage inductance to the other windings
        // For simplicity, we use the leakage to the "main" other winding (typically secondary)
        for (size_t otherIndex = 0; otherIndex < numWindings; ++otherIndex) {
            if (otherIndex != windingIndex) {
                double leakage = calculate_leakage_inductance(magnetic, windingIndex, otherIndex, frequency);
                // For multiple secondaries, we'd need a more sophisticated approach
                // Here we take the maximum leakage as the effective leakage
                totalLeakage = std::max(totalLeakage, leakage);
                break; // Use only the first non-self winding for now
            }
        }
    }
    
    return Lm_i + totalLeakage;
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
    
    // Clamp to valid range [0, 1]
    return std::min(1.0, std::max(0.0, k));
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
    
    ScalarMatrixAtFrequency result;
    result.set_frequency(frequency);
    
    std::map<std::string, std::map<std::string, DimensionWithTolerance>> magnitude;
    
    // Calculate magnetizing inductance (referred to primary)
    auto magnetizingOutput = calculate_magnetizing_inductance(magnetic, operatingPoint);
    double Lm_primary = magnetizingOutput.get_magnetizing_inductance().get_nominal().value();
    
    // Build the inductance matrix (symmetric, so compute upper triangular + diagonal)
    for (size_t i = 0; i < numWindings; ++i) {
        std::string windingName_i = get_winding_name(magnetic, i);
        double turns_i = functionalDescription[i].get_number_turns();
        
        for (size_t j = i; j < numWindings; ++j) {
            std::string windingName_j = get_winding_name(magnetic, j);
            double turns_j = functionalDescription[j].get_number_turns();
            
            DimensionWithTolerance inductanceValue;
            
            if (i == j) {
                // Diagonal element: Self inductance L_ii = Lm_i + Ll_i
                double Lm_i = Lm_primary * std::pow(turns_i / functionalDescription[0].get_number_turns(), 2);
                
                double Ll_i = 0.0;
                if (numWindings > 1) {
                    // Get leakage to the first other winding
                    size_t otherWinding = (i == 0) ? 1 : 0;
                    Ll_i = calculate_leakage_inductance(magnetic, i, otherWinding, frequency);
                }
                
                inductanceValue.set_nominal(Lm_i + Ll_i);
                magnitude[windingName_i][windingName_j] = inductanceValue;
            } else {
                // Off-diagonal element: Mutual inductance M_ij
                // M_ij = Lm_primary * (N_i / N_primary) * (N_j / N_primary)
                // This follows from M = k * sqrt(Lm_i * Lm_j) with k = 1
                double N_primary = functionalDescription[0].get_number_turns();
                double M_ij = Lm_primary * (turns_i / N_primary) * (turns_j / N_primary);
                
                inductanceValue.set_nominal(M_ij);
                
                // Set both M_ij and M_ji (symmetric)
                magnitude[windingName_i][windingName_j] = inductanceValue;
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
