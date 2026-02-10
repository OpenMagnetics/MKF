#include <source_location>
#include "physical_models/Inductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/LeakageInductance.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "processors/Inputs.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

namespace {

static double maximumError = 0.1;

/**
 * Helper function to create a two-winding transformer magnetic component
 */
OpenMagnetics::Magnetic create_two_winding_magnetic(
    const std::string& shapeName,
    const std::string& coreMaterial,
    std::vector<int64_t> numberTurns,
    std::vector<int64_t> numberParallels) {

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    return magnetic;
}

/**
 * Helper function to create a single-winding inductor magnetic component
 */
OpenMagnetics::Magnetic create_single_winding_magnetic(
    const std::string& shapeName,
    const std::string& coreMaterial,
    int64_t numberTurns,
    double gapLength = 0.001) {

    std::vector<int64_t> turns({numberTurns});
    std::vector<int64_t> parallels({1});

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, turns, parallels, wires);

    auto gapping = OpenMagnetics::Core::create_ground_gapping(gapLength, 1);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    return magnetic;
}


// ============================================================================
// Basic Inductance Matrix Tests
// ============================================================================

TEST_CASE("Test inductance matrix dimensions for single winding", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    auto magnetic = create_single_winding_magnetic("ETD 39", "3C97", 50, 0.001);
    double frequency = 100000;

    Inductance inductance;
    auto matrix = inductance.calculate_inductance_matrix(magnetic, frequency);

    // Check that matrix has correct frequency
    CHECK(matrix.get_frequency() == frequency);

    // Check that matrix has one element (1x1 for single winding)
    auto& magnitude = matrix.get_magnitude();
    CHECK(magnitude.size() == 1);

    // Check that self inductance is positive
    auto windingName = magnetic.get_coil().get_functional_description()[0].get_name();
    double selfInductance = magnitude.at(windingName).at(windingName).get_nominal().value();
    CHECK(selfInductance > 0);

    settings.reset();
}

TEST_CASE("Test inductance matrix dimensions for two windings", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    Inductance inductance;
    auto matrix = inductance.calculate_inductance_matrix(magnetic, frequency);

    // Check that matrix has 2x2 elements
    auto& magnitude = matrix.get_magnitude();
    CHECK(magnitude.size() == 2);

    auto windingName0 = magnetic.get_coil().get_functional_description()[0].get_name();
    auto windingName1 = magnetic.get_coil().get_functional_description()[1].get_name();

    // Check all elements exist
    CHECK(magnitude.count(windingName0) == 1);
    CHECK(magnitude.count(windingName1) == 1);
    CHECK(magnitude.at(windingName0).count(windingName0) == 1);
    CHECK(magnitude.at(windingName0).count(windingName1) == 1);
    CHECK(magnitude.at(windingName1).count(windingName0) == 1);
    CHECK(magnitude.at(windingName1).count(windingName1) == 1);

    settings.reset();
}

TEST_CASE("Test self inductance equals magnetizing plus leakage", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    // Calculate using Inductance class
    Inductance inductance;
    double selfInductance = inductance.calculate_self_inductance(magnetic, 0, frequency);

    // Calculate magnetizing and leakage separately
    MagnetizingInductance magnetizingModel("ZHANG");
    auto magnetizingOutput = magnetizingModel.calculate_inductance_from_number_turns_and_gapping(magnetic);
    double Lm = magnetizingOutput.get_magnetizing_inductance().get_nominal().value();

    LeakageInductance leakageModel;
    auto leakageOutput = leakageModel.calculate_leakage_inductance(magnetic, frequency, 0, 1);
    double Ll = leakageOutput.get_leakage_inductance_per_winding()[0].get_nominal().value();

    // Self inductance should be Lm + Ll
    double expectedSelfInductance = Lm + Ll;

    CHECK_THAT(selfInductance, WithinRel(expectedSelfInductance, maximumError));

    settings.reset();
}

TEST_CASE("Test mutual inductance symmetry", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    Inductance inductance;
    auto matrix = inductance.calculate_inductance_matrix(magnetic, frequency);

    auto windingName0 = magnetic.get_coil().get_functional_description()[0].get_name();
    auto windingName1 = magnetic.get_coil().get_functional_description()[1].get_name();

    auto& magnitude = matrix.get_magnitude();
    double M12 = magnitude.at(windingName0).at(windingName1).get_nominal().value();
    double M21 = magnitude.at(windingName1).at(windingName0).get_nominal().value();

    // Mutual inductance should be symmetric: M12 = M21
    CHECK_THAT(M12, WithinRel(M21, 0.001));

    settings.reset();
}

TEST_CASE("Test mutual inductance from turns ratio", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);

    Inductance inductance;
    double M = inductance.calculate_mutual_inductance(magnetic, 0, 1);

    // Get magnetizing inductance for primary
    MagnetizingInductance magnetizingModel("ZHANG");
    auto magnetizingOutput = magnetizingModel.calculate_inductance_from_number_turns_and_gapping(magnetic);
    double Lm_primary = magnetizingOutput.get_magnetizing_inductance().get_nominal().value();

    // For ideal coupling, M = Lm_primary * (N2/N1)
    double turnsRatio = double(numberTurns[1]) / double(numberTurns[0]);
    double expectedM = Lm_primary * turnsRatio;

    CHECK_THAT(M, WithinRel(expectedM, maximumError));

    settings.reset();
}

TEST_CASE("Test coupling coefficient less than or equal to 1", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    Inductance inductance;
    double k = inductance.calculate_coupling_coefficient(magnetic, 0, 1, frequency);

    // Coupling coefficient should be between 0 and 1
    CHECK(k >= 0.0);
    CHECK(k <= 1.0);

    // For a well-coupled transformer, k should be close to 1
    CHECK(k > 0.9);

    settings.reset();
}

TEST_CASE("Test magnetizing inductance referred to winding scales with turns squared", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);

    Inductance inductance;
    double Lm_primary = inductance.calculate_magnetizing_inductance_referred_to_winding(magnetic, 0);
    double Lm_secondary = inductance.calculate_magnetizing_inductance_referred_to_winding(magnetic, 1);

    // Lm_secondary = Lm_primary * (N2/N1)^2
    double turnsRatio = double(numberTurns[1]) / double(numberTurns[0]);
    double expectedLm_secondary = Lm_primary * turnsRatio * turnsRatio;

    CHECK_THAT(Lm_secondary, WithinRel(expectedLm_secondary, 0.001));

    settings.reset();
}

// ============================================================================
// Inductance Matrix for Different Core Shapes
// ============================================================================

TEST_CASE("Test inductance matrix for E core", "[physical-model][inductance]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({69, 69});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("E 42/33/20", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    Inductance inductance;
    auto matrix = inductance.calculate_inductance_matrix(magnetic, frequency);

    auto windingName0 = magnetic.get_coil().get_functional_description()[0].get_name();
    auto windingName1 = magnetic.get_coil().get_functional_description()[1].get_name();

    auto& magnitude = matrix.get_magnitude();

    // Get matrix elements
    double L11 = magnitude.at(windingName0).at(windingName0).get_nominal().value();
    double L22 = magnitude.at(windingName1).at(windingName1).get_nominal().value();
    double M12 = magnitude.at(windingName0).at(windingName1).get_nominal().value();

    // For equal turns, L11 should approximately equal L22
    CHECK_THAT(L11, WithinRel(L22, 0.1));

    // Self inductances should be greater than mutual inductance (due to leakage)
    CHECK(L11 >= M12);
    CHECK(L22 >= M12);

    // All inductances should be positive
    CHECK(L11 > 0);
    CHECK(L22 > 0);
    CHECK(M12 > 0);

    settings.reset();
}

TEST_CASE("Test inductance matrix for PQ core", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({24, 6});
    std::vector<int64_t> numberParallels({1, 1});

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 75));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 225));

    std::string shapeName = "PQ 32/25";
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;

    Inductance inductance;
    auto matrix = inductance.calculate_inductance_matrix(magnetic, frequency);

    auto windingName0 = magnetic.get_coil().get_functional_description()[0].get_name();
    auto windingName1 = magnetic.get_coil().get_functional_description()[1].get_name();

    auto& magnitude = matrix.get_magnitude();

    double L11 = magnitude.at(windingName0).at(windingName0).get_nominal().value();
    double L22 = magnitude.at(windingName1).at(windingName1).get_nominal().value();
    double M12 = magnitude.at(windingName0).at(windingName1).get_nominal().value();

    // L11 should be larger than L22 due to more turns (scales with N^2)
    CHECK(L11 > L22);

    // Mutual inductance should be between the two self inductances
    CHECK(M12 < L11);
    CHECK(M12 > 0);

    settings.reset();
}

// ============================================================================
// Multi-frequency Tests
// ============================================================================

TEST_CASE("Test inductance matrix at multiple frequencies", "[physical-model][inductance]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);

    std::vector<double> frequencies = {10000, 50000, 100000, 500000};

    Inductance inductance;
    auto matrices = inductance.calculate_inductance_matrix_per_frequency(magnetic, frequencies);

    // Check that we get correct number of matrices
    CHECK(matrices.size() == frequencies.size());

    // Check that each matrix has the correct frequency
    for (size_t i = 0; i < frequencies.size(); ++i) {
        CHECK(matrices[i].get_frequency() == frequencies[i]);
    }

    // Magnetizing inductance (and thus mutual inductance) should be relatively constant with frequency
    auto windingName0 = magnetic.get_coil().get_functional_description()[0].get_name();
    auto windingName1 = magnetic.get_coil().get_functional_description()[1].get_name();

    double M_at_10kHz = matrices[0].get_magnitude().at(windingName0).at(windingName1).get_nominal().value();
    double M_at_100kHz = matrices[2].get_magnitude().at(windingName0).at(windingName1).get_nominal().value();

    // Mutual inductance should not change significantly with frequency
    CHECK_THAT(M_at_10kHz, WithinRel(M_at_100kHz, 0.1));

    settings.reset();
}

// ============================================================================
// Consistency Tests
// ============================================================================

TEST_CASE("Test leakage inductance consistency with standalone calculation", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    // Calculate using Inductance class
    Inductance inductance;
    double Ll_from_class = inductance.calculate_leakage_inductance(magnetic, 0, 1, frequency);

    // Calculate using LeakageInductance directly
    LeakageInductance leakageModel;
    double Ll_direct = leakageModel.calculate_leakage_inductance(magnetic, frequency, 0, 1)
                           .get_leakage_inductance_per_winding()[0].get_nominal().value();

    CHECK_THAT(Ll_from_class, WithinRel(Ll_direct, 0.001));

    settings.reset();
}

TEST_CASE("Test leakage inductance matrix for two windings", "[physical-model][inductance][leakage-inductance-matrix][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    Inductance inductance;
    auto llkMatrix = inductance.calculate_leakage_inductance_matrix(magnetic, frequency);

    CHECK(llkMatrix.get_frequency() == frequency);

    auto& magnitude = llkMatrix.get_magnitude();
    CHECK(magnitude.size() == 2);

    auto windingName0 = magnetic.get_coil().get_functional_description()[0].get_name();
    auto windingName1 = magnetic.get_coil().get_functional_description()[1].get_name();

    // Diagonal must be 0 by definition
    double L00 = magnitude.at(windingName0).at(windingName0).get_nominal().value();
    double L11 = magnitude.at(windingName1).at(windingName1).get_nominal().value();
    CHECK_THAT(L00, WithinAbs(0.0, 1e-18));
    CHECK_THAT(L11, WithinAbs(0.0, 1e-18));

    // Off-diagonal should match direct LeakageInductance calculation
    LeakageInductance leakageModel;
    double L01_direct = leakageModel.calculate_leakage_inductance(magnetic, frequency, 0, 1)
                           .get_leakage_inductance_per_winding()[0].get_nominal().value();
    double L10_direct = leakageModel.calculate_leakage_inductance(magnetic, frequency, 1, 0)
                           .get_leakage_inductance_per_winding()[0].get_nominal().value();

    double L01 = magnitude.at(windingName0).at(windingName1).get_nominal().value();
    double L10 = magnitude.at(windingName1).at(windingName0).get_nominal().value();

    CHECK_THAT(L01, WithinRel(L01_direct, 0.001));
    CHECK_THAT(L10, WithinRel(L10_direct, 0.001));

    // Because leakage is referred to the source winding, matrix is generally not symmetric
    CHECK_FALSE(WithinRel(L01, 1e-12).match(L10));

    settings.reset();
}

TEST_CASE("Test leakage inductance matrix for three windings", "[physical-model][inductance][leakage-inductance-matrix][multi-winding]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({30, 15, 10});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "PQ 35/35";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;

    Inductance inductance;
    auto llkMatrix = inductance.calculate_leakage_inductance_matrix(magnetic, frequency);

    auto& magnitude = llkMatrix.get_magnitude();
    CHECK(magnitude.size() == 3);

    std::vector<std::string> windingNames;
    for (size_t i = 0; i < 3; ++i) {
        windingNames.push_back(magnetic.get_coil().get_functional_description()[i].get_name());
    }

    LeakageInductance leakageModel;

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double lij = magnitude.at(windingNames[i]).at(windingNames[j]).get_nominal().value();
            if (i == j) {
                CHECK_THAT(lij, WithinAbs(0.0, 1e-18));
            } else {
                double lij_direct = leakageModel.calculate_leakage_inductance(magnetic, frequency, i, j)
                                        .get_leakage_inductance_per_winding()[0].get_nominal().value();
                CHECK_THAT(lij, WithinRel(lij_direct, 0.001));
                CHECK(lij > 0.0);
            }
        }
    }

    settings.reset();
}

TEST_CASE("Test magnetizing inductance consistency with standalone calculation", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);

    // Calculate using Inductance class
    Inductance inductance;
    double Lm_from_class = inductance.calculate_magnetizing_inductance_referred_to_winding(magnetic, 0);

    // Calculate using MagnetizingInductance directly
    MagnetizingInductance magnetizingModel("ZHANG");
    double Lm_direct = magnetizingModel.calculate_inductance_from_number_turns_and_gapping(magnetic)
                           .get_magnetizing_inductance().get_nominal().value();

    CHECK_THAT(Lm_from_class, WithinRel(Lm_direct, 0.001));

    settings.reset();
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST_CASE("Test inductance with same winding for leakage returns zero", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    Inductance inductance;
    double Ll_same_winding = inductance.calculate_leakage_inductance(magnetic, 0, 0, frequency);

    CHECK(Ll_same_winding == 0.0);

    settings.reset();
}

TEST_CASE("Test coupling coefficient with self returns 1", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);
    double frequency = 100000;

    Inductance inductance;
    double k_self = inductance.calculate_coupling_coefficient(magnetic, 0, 0, frequency);

    CHECK(k_self == 1.0);

    settings.reset();
}

TEST_CASE("Test mutual inductance throws for self-reference", "[physical-model][inductance][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({40, 20});
    std::vector<int64_t> numberParallels({1, 1});
    auto magnetic = create_two_winding_magnetic("ETD 39", "3C97", numberTurns, numberParallels);

    Inductance inductance;

    // Mutual inductance between winding 0 and itself should throw
    CHECK_THROWS_AS(inductance.calculate_mutual_inductance(magnetic, 0, 0), std::invalid_argument);

    settings.reset();
}

// ============================================================================
// Multi-Winding Tests (3+ windings)
// ============================================================================

TEST_CASE("Test inductance matrix for three windings", "[physical-model][inductance][multi-winding]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({30, 15, 10});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "PQ 35/35";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;

    Inductance inductance;
    auto matrix = inductance.calculate_inductance_matrix(magnetic, frequency);

    // Check matrix dimensions: should be 3x3
    auto& magnitude = matrix.get_magnitude();
    CHECK(magnitude.size() == 3);

    auto windingName0 = magnetic.get_coil().get_functional_description()[0].get_name();
    auto windingName1 = magnetic.get_coil().get_functional_description()[1].get_name();
    auto windingName2 = magnetic.get_coil().get_functional_description()[2].get_name();

    // Check all 9 elements exist
    CHECK(magnitude.at(windingName0).size() == 3);
    CHECK(magnitude.at(windingName1).size() == 3);
    CHECK(magnitude.at(windingName2).size() == 3);

    // Get diagonal elements (self inductances)
    double L11 = magnitude.at(windingName0).at(windingName0).get_nominal().value();
    double L22 = magnitude.at(windingName1).at(windingName1).get_nominal().value();
    double L33 = magnitude.at(windingName2).at(windingName2).get_nominal().value();

    // Get off-diagonal elements (mutual inductances)
    double M12 = magnitude.at(windingName0).at(windingName1).get_nominal().value();
    double M13 = magnitude.at(windingName0).at(windingName2).get_nominal().value();
    double M21 = magnitude.at(windingName1).at(windingName0).get_nominal().value();
    double M23 = magnitude.at(windingName1).at(windingName2).get_nominal().value();
    double M31 = magnitude.at(windingName2).at(windingName0).get_nominal().value();
    double M32 = magnitude.at(windingName2).at(windingName1).get_nominal().value();

    // All inductances should be positive
    CHECK(L11 > 0);
    CHECK(L22 > 0);
    CHECK(L33 > 0);
    CHECK(M12 > 0);
    CHECK(M13 > 0);
    CHECK(M23 > 0);

    // Mutual inductances should be symmetric
    CHECK_THAT(M12, WithinRel(M21, 0.001));
    CHECK_THAT(M13, WithinRel(M31, 0.001));
    CHECK_THAT(M23, WithinRel(M32, 0.001));

    // Self inductances should scale with N^2
    // L11/L22 ≈ (N1/N2)^2 = (30/15)^2 = 4
    // L11/L33 ≈ (N1/N3)^2 = (30/10)^2 = 9
    double ratio_L11_L22 = L11 / L22;
    double ratio_L11_L33 = L11 / L33;
    double expected_ratio_12 = std::pow(double(numberTurns[0]) / double(numberTurns[1]), 2);
    double expected_ratio_13 = std::pow(double(numberTurns[0]) / double(numberTurns[2]), 2);

    CHECK_THAT(ratio_L11_L22, WithinRel(expected_ratio_12, 0.2));
    CHECK_THAT(ratio_L11_L33, WithinRel(expected_ratio_13, 0.2));

    settings.reset();
}

TEST_CASE("Test inductance matrix for four windings", "[physical-model][inductance][multi-winding]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({30, 10, 5, 3});
    std::vector<int64_t> numberParallels({1, 1, 1, 1});
    std::string shapeName = "PQ 35/35";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 50));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 200));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 300));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;

    Inductance inductance;
    auto matrix = inductance.calculate_inductance_matrix(magnetic, frequency);

    // Check matrix dimensions: should be 4x4
    auto& magnitude = matrix.get_magnitude();
    CHECK(magnitude.size() == 4);

    // Verify each row has 4 elements
    for (auto& row : magnitude) {
        CHECK(row.second.size() == 4);
    }

    // Get all winding names
    std::vector<std::string> windingNames;
    for (size_t i = 0; i < 4; ++i) {
        windingNames.push_back(magnetic.get_coil().get_functional_description()[i].get_name());
    }

    // Check that all diagonal elements are positive
    for (size_t i = 0; i < 4; ++i) {
        double Lii = magnitude.at(windingNames[i]).at(windingNames[i]).get_nominal().value();
        CHECK(Lii > 0);
    }

    // Check mutual inductance symmetry for all pairs
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i + 1; j < 4; ++j) {
            double Mij = magnitude.at(windingNames[i]).at(windingNames[j]).get_nominal().value();
            double Mji = magnitude.at(windingNames[j]).at(windingNames[i]).get_nominal().value();
            CHECK_THAT(Mij, WithinRel(Mji, 0.001));
            CHECK(Mij > 0);
        }
    }

    settings.reset();
}

TEST_CASE("Test coupling coefficient for three windings", "[physical-model][inductance][multi-winding]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({30, 15, 10});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "PQ 35/35";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;

    Inductance inductance;

    // Calculate coupling coefficients for all pairs
    double k01 = inductance.calculate_coupling_coefficient(magnetic, 0, 1, frequency);
    double k02 = inductance.calculate_coupling_coefficient(magnetic, 0, 2, frequency);
    double k12 = inductance.calculate_coupling_coefficient(magnetic, 1, 2, frequency);

    // All coupling coefficients should be between 0 and 1
    CHECK(k01 >= 0.0);
    CHECK(k01 <= 1.0);
    CHECK(k02 >= 0.0);
    CHECK(k02 <= 1.0);
    CHECK(k12 >= 0.0);
    CHECK(k12 <= 1.0);

    // For well-coupled transformer, all k should be close to 1
    CHECK(k01 > 0.9);
    CHECK(k02 > 0.9);
    CHECK(k12 > 0.9);

    settings.reset();
}

TEST_CASE("Test mutual inductance relationships for three windings", "[physical-model][inductance][multi-winding][smoke-test]") {
    settings.reset();
    clear_databases();

    std::vector<int64_t> numberTurns({30, 15, 10});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "PQ 35/35";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    Inductance inductance;

    // Get mutual inductances
    double M01 = inductance.calculate_mutual_inductance(magnetic, 0, 1);
    double M02 = inductance.calculate_mutual_inductance(magnetic, 0, 2);
    double M12 = inductance.calculate_mutual_inductance(magnetic, 1, 2);

    // Get magnetizing inductances referred to each winding
    double Lm0 = inductance.calculate_magnetizing_inductance_referred_to_winding(magnetic, 0);
    double Lm1 = inductance.calculate_magnetizing_inductance_referred_to_winding(magnetic, 1);
    double Lm2 = inductance.calculate_magnetizing_inductance_referred_to_winding(magnetic, 2);

    // For ideal coupling: M_ij = sqrt(Lm_i * Lm_j)
    double expected_M01 = std::sqrt(Lm0 * Lm1);
    double expected_M02 = std::sqrt(Lm0 * Lm2);
    double expected_M12 = std::sqrt(Lm1 * Lm2);

    CHECK_THAT(M01, WithinRel(expected_M01, 0.01));
    CHECK_THAT(M02, WithinRel(expected_M02, 0.01));
    CHECK_THAT(M12, WithinRel(expected_M12, 0.01));

    // Verify turns ratio relationship for mutual inductances
    // M01/M02 = N1/N2 (since M_0i = Lm0 * Ni/N0)
    double ratio_M01_M02 = M01 / M02;
    double expected_ratio = double(numberTurns[1]) / double(numberTurns[2]);
    CHECK_THAT(ratio_M01_M02, WithinRel(expected_ratio, 0.01));

    settings.reset();
}

// ============================================================================
// Benchmark Tests
// ============================================================================

TEST_CASE("Benchmark inductance matrix calculation for two windings", "[physical-model][inductance][!benchmark]") {
    BENCHMARK_ADVANCED("two winding inductance matrix")(Catch::Benchmark::Chronometer meter) {
        settings.reset();
        clear_databases();

        std::vector<int64_t> numberTurns({40, 20});
        std::vector<int64_t> numberParallels({1, 1});
        std::string shapeName = "ETD 39";

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

        auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

        std::string coreMaterial = "3C97";
        auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
        auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;

        meter.measure([&magnetic, &frequency] {
            Inductance inductance;
            return inductance.calculate_inductance_matrix(magnetic, frequency);
        });
    };
}

TEST_CASE("Benchmark inductance matrix calculation for three windings", "[physical-model][inductance][!benchmark]") {
    BENCHMARK_ADVANCED("three winding inductance matrix")(Catch::Benchmark::Chronometer meter) {
        settings.reset();
        clear_databases();

        std::vector<int64_t> numberTurns({30, 15, 10});
        std::vector<int64_t> numberParallels({1, 1, 1});
        std::string shapeName = "PQ 35/35";

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

        auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

        std::string coreMaterial = "3C97";
        auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
        auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;

        meter.measure([&magnetic, &frequency] {
            Inductance inductance;
            return inductance.calculate_inductance_matrix(magnetic, frequency);
        });
    };
}

TEST_CASE("Benchmark inductance matrix calculation for four windings", "[physical-model][inductance][!benchmark]") {
    BENCHMARK_ADVANCED("four winding inductance matrix")(Catch::Benchmark::Chronometer meter) {
        settings.reset();
        clear_databases();

        std::vector<int64_t> numberTurns({30, 10, 5, 3});
        std::vector<int64_t> numberParallels({1, 1, 1, 1});
        std::string shapeName = "PQ 35/35";

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 50));
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 200));
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 300));

        auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

        std::string coreMaterial = "3C97";
        auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
        auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;

        meter.measure([&magnetic, &frequency] {
            Inductance inductance;
            return inductance.calculate_inductance_matrix(magnetic, frequency);
        });
    };
}

TEST_CASE("Benchmark self inductance calculation", "[physical-model][inductance][!benchmark]") {
    BENCHMARK_ADVANCED("self inductance")(Catch::Benchmark::Chronometer meter) {
        settings.reset();
        clear_databases();

        std::vector<int64_t> numberTurns({40, 20});
        std::vector<int64_t> numberParallels({1, 1});
        std::string shapeName = "ETD 39";

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

        auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

        std::string coreMaterial = "3C97";
        auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
        auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;

        meter.measure([&magnetic, &frequency] {
            Inductance inductance;
            return inductance.calculate_self_inductance(magnetic, 0, frequency);
        });
    };
}

TEST_CASE("Benchmark coupling coefficient calculation", "[physical-model][inductance][!benchmark]") {
    BENCHMARK_ADVANCED("coupling coefficient")(Catch::Benchmark::Chronometer meter) {
        settings.reset();
        clear_databases();

        std::vector<int64_t> numberTurns({40, 20});
        std::vector<int64_t> numberParallels({1, 1});
        std::string shapeName = "ETD 39";

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
        wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

        auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

        std::string coreMaterial = "3C97";
        auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
        auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;

        meter.measure([&magnetic, &frequency] {
            Inductance inductance;
            return inductance.calculate_coupling_coefficient(magnetic, 0, 1, frequency);
        });
    };
}

} // anonymous namespace
