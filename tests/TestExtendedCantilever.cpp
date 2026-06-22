#include <source_location>
#include "physical_models/ExtendedCantilever.h"
#include "physical_models/LeakageInductance.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Mas.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;

static auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

// Validation of the formula layer against real hardware: the four-winding EC41-3C80 flyback
// transformer measured in Erickson & Maksimovic, "A Multiple-Winding Magnetics Model Having
// Directly Measurable Parameters". Given the measured extended-cantilever parameters (Table 1),
// the model predicts the short-circuit output inductance of each winding (Table 2). This checks
// our eq.-5 implementation against numbers obtained from a physical transformer.
TEST_CASE("Extended cantilever reproduces Erickson-Maksimovic flyback Table 2", "[physical-model][extended-cantilever]") {
    ExtendedCantileverModel model;
    model.magnetizingInductance = 220e-6;                       // L11
    model.turnsRatios = {1.0, 0.4165, 0.4154, 0.1402};         // n1..n4 (winding 1 is the reference)

    // Table 1 measured effective leakage inductances (symmetric); note the negative l34.
    double l12 = 4.5e-6, l13 = 14e-6, l14 = 130e-6;
    double l23 = 34e-6, l24 = 13e-6;
    double l34 = -35e-6;
    model.effectiveLeakageInductances = {
        {0.0, l12, l13, l14},
        {l12, 0.0, l23, l24},
        {l13, l23, 0.0, l34},
        {l14, l24, l34, 0.0},
    };

    // Table 2 predicted short-circuit output inductances, in henries.
    std::vector<double> expected = {3.2e-6, 0.52e-6, 2.3e-6, 0.36e-6};

    for (size_t windingIndex = 0; windingIndex < 4; ++windingIndex) {
        double shortCircuitInductance = ExtendedCantilever::short_circuit_output_inductance(model, windingIndex);
        // Within the paper's two-significant-figure measurement rounding (~4%).
        CHECK_THAT(shortCircuitInductance, WithinRel(expected[windingIndex], 0.06));
    }
}

// The model built from an MKF magnetic must be internally consistent: symmetric effective
// leakage inductances and positive, physically-sized short-circuit output inductances.
TEST_CASE("Extended cantilever from MKF magnetic is consistent", "[physical-model][extended-cantilever][multi-winding]") {
    settings.reset();
    std::vector<int64_t> numberTurns({50, 100, 25});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "E 42/21/15";

    std::vector<OpenMagnetics::Wire> wires;
    for (int i = 0; i < 3; ++i) wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 200));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, "3C97", gapping);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    double frequency = 100000;

    auto model = ExtendedCantilever::calculate(magnetic, frequency);

    REQUIRE(model.number_windings() == 3);
    CHECK(model.magnetizingInductance > 0);
    CHECK_THAT(model.turnsRatios[0], WithinRel(1.0, 1e-12));
    CHECK_THAT(model.turnsRatios[1], WithinRel(2.0, 1e-9));   // 100/50
    CHECK_THAT(model.turnsRatios[2], WithinRel(0.5, 1e-9));   // 25/50

    // Symmetric effective leakage inductances.
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            CHECK_THAT(model.effectiveLeakageInductances[i][j], WithinRel(model.effectiveLeakageInductances[j][i], 1e-6));
        }
    }

    // Short-circuit output inductances must be positive and far below the magnetizing inductance.
    for (size_t k = 0; k < 3; ++k) {
        double shortCircuitInductance = ExtendedCantilever::short_circuit_output_inductance(model, k);
        CHECK(shortCircuitInductance > 0);
        CHECK(shortCircuitInductance < model.magnetizingInductance);
    }

    settings.reset();
}

// For a two-winding transformer the primary short-circuit inductance (secondary shorted) is, by
// definition, the leakage inductance referred to the primary. This cross-checks the cantilever
// derivation against the standalone LeakageInductance solver.
TEST_CASE("Extended cantilever short-circuit inductance matches leakage for two windings", "[physical-model][extended-cantilever]") {
    settings.reset();
    std::vector<int64_t> numberTurns({64, 20});
    std::vector<int64_t> numberParallels({1, 1});
    std::string shapeName = "E 42/33/20";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 25));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 225));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, "3C97", gapping);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    double frequency = 100000;

    auto model = ExtendedCantilever::calculate(magnetic, frequency);
    double shortCircuitPrimary = ExtendedCantilever::short_circuit_output_inductance(model, 0);
    double leakageReferredToPrimary = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1)
        .get_leakage_inductance_per_winding()[0].get_nominal().value();

    CHECK_THAT(shortCircuitPrimary, WithinRel(leakageReferredToPrimary, 0.05));

    settings.reset();
}

// End-to-end validation on a real, fully-processed MAS file. concentric_transformer.json is a
// two-winding component with full geometry; loading it and running the whole geometry→Λ→
// cantilever chain must reproduce the standalone leakage solver's primary short-circuit
// inductance. This exercises the actual MAS pipeline, not a programmatic quick-coil.
TEST_CASE("Extended cantilever from a real MAS file matches leakage solver", "[physical-model][extended-cantilever]") {
    settings.reset();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_transformer.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path.string(), mas);
    auto magnetic = mas.get_magnetic();
    double frequency = 100000;

    auto model = ExtendedCantilever::calculate(magnetic, frequency);
    REQUIRE(model.number_windings() == 2);

    double shortCircuitPrimary = ExtendedCantilever::short_circuit_output_inductance(model, 0);
    double leakageReferredToPrimary = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1)
        .get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(shortCircuitPrimary, WithinRel(leakageReferredToPrimary, 0.02));

    settings.reset();
}

// Multi-winding components round-trip through a MAS file: serialize a three-winding magnetic,
// reload it, and confirm the extended-cantilever model is identical to the in-memory result.
// This proves the cantilever derivation works from a MAS-file-described multi-winding part
// (the repo has no pre-processed 3+ winding fixture, and the catalog example files reference
// wires absent from the default database).
TEST_CASE("Extended cantilever round-trips through a MAS file for three windings", "[physical-model][extended-cantilever][multi-winding]") {
    settings.reset();
    std::vector<int64_t> numberTurns({40, 20, 10});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "E 42/21/15";

    std::vector<OpenMagnetics::Wire> wires;
    for (int i = 0; i < 3; ++i) wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 200));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, "3C97", gapping);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    double frequency = 100000;

    auto inMemoryModel = ExtendedCantilever::calculate(magnetic, frequency);

    // Serialize to a MAS file and read it back.
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    auto outFile = outputFilePath;
    outFile.append("Test_Extended_Cantilever_Three_Windings.mas.json");
    std::filesystem::remove(outFile);
    OpenMagnetics::to_file(outFile, mas);

    OpenMagnetics::Mas reloadedMas;
    OpenMagnetics::from_file(outFile.string(), reloadedMas);
    auto reloadedModel = ExtendedCantilever::calculate(reloadedMas.get_magnetic(), frequency);

    REQUIRE(reloadedModel.number_windings() == 3);
    CHECK_THAT(reloadedModel.magnetizingInductance, WithinRel(inMemoryModel.magnetizingInductance, 1e-6));
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            if (i == j) continue;
            CHECK_THAT(reloadedModel.effectiveLeakageInductances[i][j],
                       WithinRel(inMemoryModel.effectiveLeakageInductances[i][j], 1e-4));
        }
    }
    for (size_t k = 0; k < 3; ++k) {
        double shortCircuitInductance = ExtendedCantilever::short_circuit_output_inductance(reloadedModel, k);
        CHECK(shortCircuitInductance > 0);
        CHECK(shortCircuitInductance < reloadedModel.magnetizingInductance);
    }

    settings.reset();
}
