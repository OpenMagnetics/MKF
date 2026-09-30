// Coupling coefficient, mutual inductance and self inductance must be the entries of ONE
// inductance matrix. From the energy formulation MKF uses (core flux with permeance P_m plus
// the window field, whose energy is W = ½·iᵀΛi):
//
//   L_ij = N_i·N_j·P_m + Λ_ij,     k_ij = L_ij / √(L_ii·L_jj)
//
// and for ampere-turn balanced currents (i_j = −r·i_i, r = N_i/N_j) the magnetizing term
// cancels exactly, so the pairwise leakage referred to winding i is
//
//   L_pair = L_ii + r²·L_jj − 2r·L_ij = Λ_ii + r²·Λ_jj − 2r·Λ_ij.
//
// calculate_coupling_coefficient used to build M from the magnetizing term alone,
// M = √(Lm_i·Lm_j), dropping the mutual leakage Λ_ij that the matrix (and so the loop leakage)
// contains. On a sector-wound toroidal common-mode choke that gave k = 0.99596 against the
// 0.99745 the matrix implies.
#include <source_location>
#include "physical_models/Inductance.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/ExtendedCantilever.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Wire.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Tolerances. Where the test repeats the model's own arithmetic on the same doubles, the result
// is the same double and the check is exact (==). Where it compares two independent solves
// that agree only in exact arithmetic, the allowed gap is derived from rounding:
//
// - Every inductance MKF computes is a quadrature of the field energy. A chain of n correctly
//   rounded operations has relative error at most γ_n = n·u/(1 − n·u), u = ε/2 (Higham,
//   Accuracy and Stability of Numerical Algorithms, 2nd ed., Lemma 3.1). γ_n <= √ε holds for
//   n <= 2/(√ε·(1 + √ε)) ≈ 1.3e8 roundings, orders more than any of the field solvers performs,
//   so √ε is the error budget of each independently computed inductance.
// - The test's own combination of those inductances adds a few roundings, bounded by γ_m for
//   the m operations counted at the use site.
// - A difference of large, nearly equal terms amplifies their errors by the condition number
//   of the cancellation, Σ|terms| / |result| (Higham §1.7): that factor is computed from the
//   case's own values, not chosen.
const double unitRoundoff = std::numeric_limits<double>::epsilon() / 2;
const double perInductanceRelativeError = std::sqrt(std::numeric_limits<double>::epsilon());
double roundingBound(int roundings) {
    return roundings * unitRoundoff / (1 - roundings * unitRoundoff);
}

auto& couplingSettings = Settings::GetInstance();

// A 1:1 common-mode choke, sector wound: each winding on its own half of a T 20/10/7 ring of a
// 10 000-permeability MnZn. The two windings do not overlap, so the window field of one partly
// threads the other: the mutual leakage Λ01 is not negligible against the self-leakage.
OpenMagnetics::Magnetic sector_wound_choke() {
    couplingSettings.set_use_toroidal_cores(true);
    couplingSettings.set_coil_delimit_and_compact(false);
    couplingSettings.set_coil_try_rewind(false);
    std::vector<int64_t> numberTurns = {19, 19};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string coreShape = "T 20/10/7";
    auto emptyGapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, 1, "3E6");
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, 1,
                                                     WindingOrientation::CONTIGUOUS, WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}

OpenMagnetics::Magnetic gapped_transformer(const std::string& shapeName, std::vector<int64_t> numberTurns) {
    std::vector<int64_t> numberParallels(numberTurns.size(), 1);
    std::vector<OpenMagnetics::Wire> wires;
    for (size_t i = 0; i < numberTurns.size(); ++i) {
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    }
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, "3C97", gapping);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}

// Every pair of windings: the scalar entry points agree with the matrix, and the loop leakage
// the matrix implies is the pairwise leakage LeakageInductance solves for independently (an
// ampere-turn balanced excitation, not the polarization the matrix is assembled from).
void check_consistency(OpenMagnetics::Magnetic magnetic, double frequency) {
    Inductance inductance;
    auto matrix = inductance.calculate_inductance_matrix_values(magnetic, frequency);
    auto named = inductance.calculate_inductance_matrix(magnetic, frequency).get_magnitude();
    auto& functionalDescription = magnetic.get_coil().get_functional_description();
    size_t numberWindings = functionalDescription.size();
    REQUIRE(matrix.size() == numberWindings);

    // The simulator exporters read ExtendedCantilever's matrix: it must be the same one.
    auto cantileverMatrix = ExtendedCantilever::calculate_inductance_matrix(magnetic, frequency);

    for (size_t i = 0; i < numberWindings; ++i) {
        auto name_i = functionalDescription[i].get_name();
        // Same matrix, read back: the same doubles.
        CHECK(inductance.calculate_self_inductance(magnetic, i, frequency) == matrix[i][i]);
        CHECK(named[name_i][name_i].get_nominal().value() == matrix[i][i]);
        for (size_t j = 0; j < numberWindings; ++j) {
            CHECK(cantileverMatrix[i][j] == matrix[i][j]);
            if (i == j) {
                continue;
            }
            auto name_j = functionalDescription[j].get_name();
            double mutual = inductance.calculate_mutual_inductance(magnetic, i, j, frequency);
            double coupling = inductance.calculate_coupling_coefficient(magnetic, i, j, frequency);
            double couplingFromMatrix = matrix[i][j] / std::sqrt(matrix[i][i] * matrix[j][j]);
            INFO("windings " << i << "-" << j << ": M " << mutual << " H, k " << coupling
                 << ", matrix k " << couplingFromMatrix);
            // Same entries and the model's own operation order (M / √(L_ii·L_jj)): the same doubles.
            CHECK(mutual == matrix[i][j]);
            CHECK(named[name_i][name_j].get_nominal().value() == matrix[i][j]);
            CHECK(coupling == couplingFromMatrix);

            double turnsRatio = double(functionalDescription[i].get_number_turns()) / functionalDescription[j].get_number_turns();
            double loopLeakage = matrix[i][i] + turnsRatio * turnsRatio * matrix[j][j] - 2 * turnsRatio * mutual;
            double pairwiseLeakage = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, i, j)
                                         .get_leakage_inductance_per_winding()[0].get_nominal().value();
            INFO("loop leakage from the matrix " << loopLeakage << " H, pairwise leakage " << pairwiseLeakage << " H");
            REQUIRE(pairwiseLeakage > 0);
            // L_ii + r²·L_jj − 2r·M = Λ_ii + r²·Λ_jj − 2r·Λ_ij is exact in exact arithmetic (the
            // magnetizing terms cancel identically and the field model is linear in the currents).
            // In floating point each of the four inductances carries up to √ε relative error, and
            // the cancellation turns that into an absolute error of √ε·Σ|terms|; forming r, r², the
            // three products and the two sums adds at most γ_7 of Σ|terms|.
            double sumOfMagnitudes = matrix[i][i] + turnsRatio * turnsRatio * matrix[j][j] + 2 * turnsRatio * std::abs(mutual);
            double loopTolerance = perInductanceRelativeError * (sumOfMagnitudes + pairwiseLeakage) + roundingBound(7) * sumOfMagnitudes;
            INFO("cancellation condition number " << sumOfMagnitudes / pairwiseLeakage << ", tolerance " << loopTolerance << " H");
            CHECK_THAT(loopLeakage, WithinAbs(pairwiseLeakage, loopTolerance));
        }
    }
}

}  // namespace

// Analytical case: a 1:1 choke. With N1 = N2 the magnetizing terms of L11, L22 and L12 are all
// the same Lm, so the differential-mode loop inductance (series opposing, what a DM measurement
// sees) is L_DM = L11 + L22 − 2·L12 and therefore
//
//   k = L12/√(L11·L22) = (L11 + L22 − L_DM) / (2·√(L11·L22)),
//
// the coupling a bench measurement of the two open-circuit self inductances and the DM
// inductance gives. L_DM here is LeakageInductance's balanced-excitation solve, an independent
// path from the one calculate_coupling_coefficient reads.
TEST_CASE("Test_Coupling_Coefficient_Of_1to1_Choke_Follows_From_Self_And_DM_Inductance", "[physical-model][inductance][coupling][cmc][toroidal][smoke-test]") {
    couplingSettings.reset();
    clear_databases();
    auto magnetic = sector_wound_choke();
    double frequency = 100000;

    Inductance inductance;
    double L11 = inductance.calculate_self_inductance(magnetic, 0, frequency);
    double L22 = inductance.calculate_self_inductance(magnetic, 1, frequency);
    double differentialModeInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1)
                                            .get_leakage_inductance_per_winding()[0].get_nominal().value();
    double expectedCoupling = (L11 + L22 - differentialModeInductance) / (2 * std::sqrt(L11 * L22));
    double coupling = inductance.calculate_coupling_coefficient(magnetic, 0, 1, frequency);

    // What the coupling was with the mutual leakage left out of M.
    double magnetizingOnlyCoupling = std::sqrt(inductance.calculate_magnetizing_inductance_referred_to_winding(magnetic, 0) *
                                               inductance.calculate_magnetizing_inductance_referred_to_winding(magnetic, 1)) /
                                     std::sqrt(L11 * L22);
    INFO("L11 " << L11 << " H, L22 " << L22 << " H, L_DM " << differentialModeInductance << " H, k " << coupling
         << ", expected " << expectedCoupling << ", magnetizing-only " << magnetizingOnlyCoupling);

    // Tolerance on 1 − k. The model's k and expectedCoupling share the denominator √(L11·L22)
    // (the same doubles, the same operations), so they differ only through the numerator:
    // L12 from the matrix against (L11 + L22 − L_DM)/2 from the independent DM solve, equal in
    // exact arithmetic. With √ε relative error on each of L11, L22 and L_DM the numerator is off
    // by at most √ε·(L11 + L22 + L_DM), and the two sums, the product, the root and the division
    // add at most γ_6 of it. Divided by 2·√(L11·L22) that is an error in k, and relative to 1 − k
    // it is amplified by 1/(1 − k): the condition number of reading the leakage out of k.
    double numeratorMagnitude = L11 + L22 + differentialModeInductance;
    double couplingTolerance = (perInductanceRelativeError + roundingBound(6)) * numeratorMagnitude /
                               (2 * std::sqrt(L11 * L22) * (1 - expectedCoupling));
    INFO("relative tolerance on 1 - k " << couplingTolerance);

    // The case must be able to tell the two apart: the old formula falls outside that tolerance.
    REQUIRE(std::abs((1 - magnetizingOnlyCoupling) - (1 - expectedCoupling)) > couplingTolerance * (1 - expectedCoupling));

    // Compare 1 − k, the part that carries the information (k itself is 0.99…).
    CHECK_THAT(1 - coupling, WithinRel(1 - expectedCoupling, couplingTolerance));
    CHECK(coupling < 1.0);

    // The parallel-connected common-mode inductance a choke datasheet states, (L11 + L22 + 2·L12)/4,
    // must follow from the same k: k is M / √(L11·L22) in the model's operation order.
    double mutual = inductance.calculate_mutual_inductance(magnetic, 0, 1, frequency);
    CHECK(coupling == mutual / std::sqrt(L11 * L22));
    couplingSettings.reset();
}

TEST_CASE("Test_Coupling_Mutual_And_Self_Inductance_Agree_With_Inductance_Matrix_Choke", "[physical-model][inductance][coupling][cmc][toroidal][smoke-test]") {
    couplingSettings.reset();
    clear_databases();
    check_consistency(sector_wound_choke(), 100000);
    couplingSettings.reset();
}

TEST_CASE("Test_Coupling_Mutual_And_Self_Inductance_Agree_With_Inductance_Matrix_Transformer", "[physical-model][inductance][coupling][smoke-test]") {
    couplingSettings.reset();
    clear_databases();
    check_consistency(gapped_transformer("ETD 39", {40, 20}), 100000);
    couplingSettings.reset();
}

TEST_CASE("Test_Coupling_Mutual_And_Self_Inductance_Agree_With_Inductance_Matrix_Three_Windings", "[physical-model][inductance][coupling][multi-winding][smoke-test]") {
    couplingSettings.reset();
    clear_databases();
    check_consistency(gapped_transformer("PQ 35/35", {30, 15, 10}), 100000);
    couplingSettings.reset();
}

TEST_CASE("Test_Coupling_Winding_Index_Out_Of_Range_Throws", "[physical-model][inductance][coupling][smoke-test]") {
    couplingSettings.reset();
    clear_databases();
    auto magnetic = gapped_transformer("ETD 39", {40, 20});
    Inductance inductance;
    CHECK_THROWS(inductance.calculate_coupling_coefficient(magnetic, 0, 2, 100000));
    CHECK_THROWS(inductance.calculate_mutual_inductance(magnetic, 2, 0, 100000));
    CHECK_THROWS(inductance.calculate_self_inductance(magnetic, 2, 100000));
    couplingSettings.reset();
}

// Leg-separated windings: the cantilever model the simulator exporters use must not keep its
// own copy of the assembly. It used to add the window leakage between windings on different
// columns, which Inductance::calculate_inductance_matrix leaves out (the window field solvers
// are not window-aware), so the exported coupling differed from calculate_coupling_coefficient.
TEST_CASE("Test_Coupling_Extended_Cantilever_Matrix_Is_Inductance_Matrix_Multi_Column", "[physical-model][inductance][coupling][extended-cantilever][multi-column]") {
    couplingSettings.reset();
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("testData").append("multicolumn_e42_transformer.json");
    std::ifstream masFile(path);
    REQUIRE(masFile.good());
    OpenMagnetics::Magnetic magnetic(json::parse(masFile)["magnetic"]);
    double frequency = 100000;

    auto matrix = Inductance().calculate_inductance_matrix_values(magnetic, frequency);
    auto cantileverMatrix = ExtendedCantilever::calculate_inductance_matrix(magnetic, frequency);
    REQUIRE(cantileverMatrix.size() == matrix.size());
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix.size(); ++j) {
            INFO("entry " << i << "," << j << ": Inductance " << matrix[i][j] << " H, cantilever " << cantileverMatrix[i][j] << " H");
            CHECK(cantileverMatrix[i][j] == matrix[i][j]);
        }
    }
    double coupling = Inductance().calculate_coupling_coefficient(magnetic, 0, 1, frequency);
    CHECK(coupling == cantileverMatrix[0][1] / std::sqrt(cantileverMatrix[0][0] * cantileverMatrix[1][1]));
    couplingSettings.reset();
}
