#include <source_location>
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "physical_models/StrayCapacitance.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include "support/Painter.h"
#include "processors/Sweeper.h"
#include "processors/CircuitSimulatorInterface.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using namespace OpenMagneticsTesting;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

static double maximumError = 0.4;
static auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
static bool plot = true;

TEST_CASE("Calculate capacitance among two windings each with 1 turn and 1 parallel", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();
    auto coilJsonStr = R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })";
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

    auto [core, coil] = prepare_core_and_coil_from_json(coreJsonStr, coilJsonStr);

    StrayCapacitance strayCapacitance;

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 2.8708e-12}, {"Secondary", -2.4259e-12}}},
        {"Secondary", {{"Primary", -2.4259e-12}, {"Secondary", 2.887e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Primary"].size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 2);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }
}

TEST_CASE("Tripole is the positive 3-cap pi-model; 6C is the canonical Biela/Kolar network", "[physical-model][stray-capacitance][models]") {
    settings.reset();
    auto coilJsonStr = R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 10, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 10, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })";
    auto coreJsonStr = R"({"name": "core", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";
    auto [core, coil] = prepare_core_and_coil_from_json(coreJsonStr, coilJsonStr);

    auto output = StrayCapacitance().calculate_capacitance(coil);
    auto amongWindings = output.get_capacitance_among_windings().value();
    auto tripole = output.get_tripole_capacitance_per_winding().value().at("Primary").at("Secondary");
    auto sixCap = output.get_six_capacitor_network_per_winding().value().at("Primary").at("Secondary");

    double cselfP = amongWindings.at("Primary").at("Primary");
    double cselfS = amongWindings.at("Secondary").at("Secondary");
    double cinter = amongWindings.at("Primary").at("Secondary");

    // Tripole = the measurable positive 3-capacitor pi-model: {primary self, secondary self, inter}.
    CHECK(tripole.get_c1() == Catch::Approx(cselfP));
    CHECK(tripole.get_c2() == Catch::Approx(cselfS));
    CHECK(tripole.get_c3() == Catch::Approx(cinter));
    // The three tripole caps must be strictly positive (the whole point of the pi-model).
    CHECK(tripole.get_c1() > 0.0);
    CHECK(tripole.get_c2() > 0.0);
    CHECK(tripole.get_c3() > 0.0);

    // Six-capacitor network = Biela/Kolar eq. (30) from C0 = the inter-winding static capacitance:
    // C1=C2=-C0/6, C3=C4=C0/3, C5=C6=C0/6.
    CHECK(sixCap.get_c3() == Catch::Approx(cinter / 3.0));
    CHECK(sixCap.get_c4() == Catch::Approx(sixCap.get_c3()));
    CHECK(sixCap.get_c1() == Catch::Approx(-sixCap.get_c3() / 2.0));   // -C0/6 = -(C0/3)/2
    CHECK(sixCap.get_c2() == Catch::Approx(sixCap.get_c1()));
    CHECK(sixCap.get_c5() == Catch::Approx(sixCap.get_c3() / 2.0));    // C0/6 = (C0/3)/2
    CHECK(sixCap.get_c6() == Catch::Approx(sixCap.get_c5()));
}

TEST_CASE("Calculate capacitance of a winding with 8 turns and 1 parallel", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();
    auto coilJsonStr = R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })";
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

    auto [core, coil] = prepare_core_and_coil_from_json(coreJsonStr, coilJsonStr);

    StrayCapacitance strayCapacitance;

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 10.166e-12 - 9.8337e-12}}},
    };


    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

    auto magnetic = prepare_magnetic_from_json(coreJsonStr, coilJsonStr);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Get_Capacitance_Among_Windings_1_Windings_8_Turns_1_Parallels.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Calculate capacitance among two windings each with 8 turns and 1 parallel", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();
    auto coilJsonStr = R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })";
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

    auto [core, coil] = prepare_core_and_coil_from_json(coreJsonStr, coilJsonStr);

    StrayCapacitance strayCapacitance;

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 10.166e-12}, {"Secondary", -9.8337e-12}}},
        {"Secondary", {{"Primary", -9.8337e-12}, {"Secondary", 10.199e-12}}},
    };


    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Primary"].size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 2);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

    auto magnetic = prepare_magnetic_from_json(coreJsonStr, coilJsonStr);

    auto outFile = outputFilePath;
    outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels.svg");
    std::filesystem::remove(outFile);
    Painter painter(outFile);
    painter.paint_core(magnetic);
    painter.paint_bobbin(magnetic);
    painter.paint_coil_turns(magnetic);
    painter.export_svg();

}
TEST_CASE("Calculate capacitance among two windings one with 16 and another with 8 turns and both 1 parallel", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();
    auto coilJsonStr = R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 16, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })";
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

    auto [core, coil] = prepare_core_and_coil_from_json(coreJsonStr, coilJsonStr);

    StrayCapacitance strayCapacitance;

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 14.777e-12}, {"Secondary", -6.6715e-12}}},
        {"Secondary", {{"Primary", -6.6715e-12}, {"Secondary", 12.791e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Primary"].size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 2);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

}

TEST_CASE("Calculate capacitance among three windings each with 8 turns and 1 parallel", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();
    auto coilJsonStr = R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" }, {"name": "Tertiary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "tertiary", "wire": "Round 1.00 - Grade 1" } ] })";
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

    auto [core, coil] = prepare_core_and_coil_from_json(coreJsonStr, coilJsonStr);

    StrayCapacitance strayCapacitance;

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 15.437e-12}, {"Secondary", -15.115e-12}, {"Tertiary", -0.12966e-12}}},
        {"Secondary", {{"Primary", -15.115e-12}, {"Secondary", 27.457e-12}, {"Tertiary", -12.21e-12}}},
        {"Tertiary", {{"Primary", -0.12966e-12}, {"Secondary", -12.21e-12}, {"Tertiary", 12.791e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 3);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Primary"].size() == 3);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 3);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Tertiary"].size() == 3);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            if (capacitance != 0) {
                CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
            }
            else {
                CHECK_THAT(capacitance, WithinAbs(expectedValues[firstKey][secondKey], 1e-12));
            }
        }
    }

}

// TEST_CASE("Calculate and compare capacitance matrix of one winding versus two windings", "[physical-model][stray-capacitance][smoke-test]") {
//     {
//         std::vector<int64_t> numberTurns = {110};
//         std::vector<int64_t> numberParallels = {1};
//         std::string coreShapeName = "T 12.5/7.5/5";
//         std::vector<OpenMagnetics::Wire> wires;
//         auto wire = find_wire_by_name("Round 0.15 - Grade 1");
//         wires.push_back(wire);
        
//         auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

//         int64_t numberStacks = 1;
//         std::string coreMaterialName = "A07";
//         std::vector<CoreGap> gapping = {};
//         auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
//         OpenMagnetics::Magnetic magnetic;
//         magnetic.set_core(core);
//         magnetic.set_coil(coil);

//         StrayCapacitance strayCapacitance;

//         std::map<std::string, std::map<std::string, double>> expectedValues = {
//             {"winding 0", {{"winding 0", 3.15617e-12 -3.14887e-12}}},
//         };

//         auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
//         CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);
//         for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
//             for (auto [secondKey, capacitanceWithTolerance] : aux) {
//                 auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
//                 CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
//             }
//         }

//         if (plot) {
//             auto outFile = outputFilePath;
//             outFile.append("Test_One_Versus_Two_Windings_One.svg");
//             std::filesystem::remove(outFile);
//             Painter painter(outFile);
//             painter.paint_core(magnetic);
//             painter.paint_bobbin(magnetic);
//             painter.paint_coil_turns(magnetic);
//             painter.export_svg();
//         }
//     }
//     {
//         std::vector<int64_t> numberTurns = {110, 110};
//         std::vector<int64_t> numberParallels = {1, 1};
//         std::string coreShapeName = "T 12.5/7.5/5";
//         std::vector<OpenMagnetics::Wire> wires;
//         auto wire = find_wire_by_name("Round 0.15 - Grade 1");
//         wires.push_back(wire);
//         wires.push_back(wire);

//         auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

//         int64_t numberStacks = 1;
//         std::string coreMaterialName = "A07";
//         std::vector<CoreGap> gapping = {};
//         auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
//         OpenMagnetics::Magnetic magnetic;
//         magnetic.set_core(core);
//         magnetic.set_coil(coil);

//         StrayCapacitance strayCapacitance;

//         std::map<std::string, std::map<std::string, double>> expectedValues = {
//             {"winding 0", {{"winding 0", 3.15617e-12}, {"winding 1", -3.14887e-12}}},
//             {"winding 1", {{"winding 0", -3.14887e-12}, {"winding 1", 3.15674e-12}}},
//         };

//         auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
//         CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
//         CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["winding 0"].size() == 2);
//         CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["winding 1"].size() == 2);
//         for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
//             for (auto [secondKey, capacitanceWithTolerance] : aux) {
//                 auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
//                 CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
//             }
//         }

//         if (plot) {
//             auto outFile = outputFilePath;
//             outFile.append("Test_One_Versus_Two_Windings_Two.svg");
//             std::filesystem::remove(outFile);
//             Painter painter(outFile);
//             painter.paint_core(magnetic);
//             painter.paint_bobbin(magnetic);
//             painter.paint_coil_turns(magnetic);
//             painter.export_svg();
//         }
//     }
// }

TEST_CASE("Benchmakrs stray capacitance calculation", "[physical-model][stray-capacitance][!benchmark]") {
    BENCHMARK_ADVANCED("measures computation time")(Catch::Benchmark::Chronometer meter) {
        std::vector<int64_t> numberTurns = {110, 110};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string coreShapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 0.15 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

        StrayCapacitance strayCapacitance;

        meter.measure([&strayCapacitance, &coil] { return strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value(); });
    };
}

TEST_CASE("Calculate capacitance of an automatic buck produced in OM with three layers", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "buck_inductor three layers.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();

    StrayCapacitance strayCapacitance;

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 2.6849e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Get_Capacitance_Automatic_Buck_OM.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Calculate capacitance of an automatic buck produced in OM with two layers", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "buck_inductor two layers.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();

    StrayCapacitance strayCapacitance;

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 2.8e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Get_Capacitance_Automatic_Buck_OM.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Calculate area between two round turns", "[physical-model][stray-capacitance][smoke-test]") {
    {
        QuickMagneticConfig config;
        config.numberTurns = {1, 1};
        config.numberParallels = {1, 1};
        auto magnetic = create_quick_test_magnetic(config);
        auto coil = magnetic.get_coil();

        StrayCapacitance strayCapacitance;
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[1];
        auto area = strayCapacitance.calculate_area_between_two_turns(firstTurn, secondTurn);
        CHECK(area < 1e-6);
        CHECK(area > 0.86e-6);
    }
}

TEST_CASE("Calculate area between two rectangular turns", "[physical-model][stray-capacitance][smoke-test]") {
    {
        QuickMagneticConfig config;
        config.numberTurns = {1, 1};
        config.numberParallels = {1, 1};
        config.wireNames = {"Rectangular 2x0.80 - Grade 1", "Rectangular 2x0.80 - Grade 1"};
        auto magnetic = create_quick_test_magnetic(config);
        auto coil = magnetic.get_coil();

        StrayCapacitance strayCapacitance;
        auto layerThickness = coil.get_layers_description().value()[1].get_dimensions()[0];
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[1];
        auto area = strayCapacitance.calculate_area_between_two_turns(firstTurn, secondTurn);

        CHECK_THAT(area, WithinRel(layerThickness * 0.0008, maximumError));
    }
}

TEST_CASE("Calculate area between round and rectangular turns", "[physical-model][stray-capacitance][smoke-test]") {
    {
        QuickMagneticConfig config;
        config.numberTurns = {1, 1};
        config.numberParallels = {1, 1};
        config.wireNames = {"Round 2.00 - Grade 1", "Rectangular 2x0.80 - Grade 1"};
        auto magnetic = create_quick_test_magnetic(config);
        auto coil = magnetic.get_coil();

        StrayCapacitance strayCapacitance;
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[1];
        auto area = strayCapacitance.calculate_area_between_two_turns(firstTurn, secondTurn);
        CHECK(area < 0.6e-6);
        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Energy_Between_Two_turns.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
}

TEST_CASE("Calculate area between round and rectangular turns diagonally", "[physical-model][stray-capacitance][smoke-test]") {
    {
        QuickMagneticConfig config;
        config.numberTurns = {2, 1};
        config.numberParallels = {1, 1};
        config.wireNames = {"Round 2.00 - Grade 1", "Rectangular 2x0.80 - Grade 1"};
        auto magnetic = create_quick_test_magnetic(config);
        auto coil = magnetic.get_coil();

        StrayCapacitance strayCapacitance;
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[2];
        auto area = strayCapacitance.calculate_area_between_two_turns(firstTurn, secondTurn);
        CHECK(area < 2e-6);
        CHECK(area > 1.5e-6);
        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Energy_Between_Two_turns.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
}

TEST_CASE("Calculate energy density between two turns", "[physical-model][stray-capacitance][smoke-test]") {
    {
        QuickMagneticConfig config;
        config.numberTurns = {1, 1};
        config.numberParallels = {1, 1};
        auto magnetic = create_quick_test_magnetic(config);
        auto coil = magnetic.get_coil();
        auto wire = find_wire_by_name("Round 2.00 - Grade 1");

        StrayCapacitance strayCapacitance;
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[1];
        [[maybe_unused]] auto energy = strayCapacitance.calculate_energy_density_between_two_turns(firstTurn, wire, secondTurn, wire, 10);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Energy_Between_Two_turns.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
}

TEST_CASE("Calculate capacitance of a tranformers with low filling factor", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "low_filling_transformer.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 8.8291e-12}, {"Secondary", -4.2151e-12}}},
        {"Secondary", {{"Primary", -4.2151e-12}, {"Secondary", 6.6748e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Primary"].size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 2);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }


    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_low_filling_transformer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Calculate capacitance of a one layer inductor", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "low_filling_inductor.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Secondary", {{"Secondary", 0.7e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 1);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 1);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_low_filling_transformer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Calculate capacitance of a simple planar tranformer", "[physical-model][stray-capacitance][planar][smoke-test]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simple_planar.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 38.026e-12}, {"Secondary", -36.948e-12}}},
        {"Secondary", {{"Primary", -36.948e-12}, {"Secondary", 37.991e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Primary"].size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 2);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_low_filling_transformer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Calculate capacitance of a simple planar tranformer with a turns ratio of 5", "[physical-model][stray-capacitance][planar][smoke-test]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simple_planar_5_to_1.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 18.026e-12}, {"Secondary", -16.948e-12}}},
        {"Secondary", {{"Primary", -16.948e-12}, {"Secondary", 17.991e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Primary"].size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 2);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_low_filling_transformer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Calculate capacitance of a simple planar tranformer with imperfect overlapping", "[physical-model][stray-capacitance][planar][smoke-test]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simple_planar_imperfect_overlapping.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 37.003e-12}, {"Secondary", -27.231e-12}}},
        {"Secondary", {{"Primary", -27.231e-12}, {"Secondary", 25.543e-12}}},
    };

    auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Primary"].size() == 2);
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["Secondary"].size() == 2);
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_low_filling_transformer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}


// ============================================================================
// LITZ WIRE STRAY CAPACITANCE TESTS
// ============================================================================

TEST_CASE("Calculate capacitance of litz winding with 8 turns", "[physical-model][stray-capacitance][litz][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {8, 8};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string coreShapeName = "RM 10/I";
    std::vector<OpenMagnetics::Wire> wires;
    
    // Create a litz wire: 50 strands x 0.1mm, Grade 1
    auto litzWire = OpenMagnetics::Wire::create_quick_litz_wire(0.0001, 50);
    wires.push_back(litzWire);
    wires.push_back(litzWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            // Diagonal positive, off-diagonal negative
            if (firstKey == secondKey) {
                CHECK(capacitance > 0);
            } else {
                CHECK(capacitance < 0);
            }
            CHECK(std::abs(capacitance) < 1e-9); // Sanity: sub-nF range
        }
    }
}

TEST_CASE("Calculate capacitance of single litz winding with 8 turns", "[physical-model][stray-capacitance][litz][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {8};
    std::vector<int64_t> numberParallels = {1};
    std::string coreShapeName = "RM 10/I";
    std::vector<OpenMagnetics::Wire> wires;
    
    // Create litz wire: 100 strands x 0.08mm, Grade 1
    auto litzWire = OpenMagnetics::Wire::create_quick_litz_wire(0.00008, 100);
    wires.push_back(litzWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK(capacitance > 0);
            CHECK(capacitance < 100e-12); // Reasonable upper bound for this geometry
        }
    }
}

TEST_CASE("Litz capacitance with Koch model", "[physical-model][stray-capacitance][litz][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {8};
    std::vector<int64_t> numberParallels = {1};
    std::string coreShapeName = "RM 10/I";
    std::vector<OpenMagnetics::Wire> wires;
    
    auto litzWire = OpenMagnetics::Wire::create_quick_litz_wire(0.0001, 50);
    wires.push_back(litzWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    // Test that Koch model also works with litz
    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::KOCH);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK(capacitance > 0);
        }
    }
}

TEST_CASE("Litz capacitance with Massarini model", "[physical-model][stray-capacitance][litz][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {8};
    std::vector<int64_t> numberParallels = {1};
    std::string coreShapeName = "RM 10/I";
    std::vector<OpenMagnetics::Wire> wires;
    
    auto litzWire = OpenMagnetics::Wire::create_quick_litz_wire(0.0001, 50);
    wires.push_back(litzWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::MASSARINI);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK(capacitance > 0);
        }
    }
}

TEST_CASE("Litz capacitance with Duerdoth model", "[physical-model][stray-capacitance][litz][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {8};
    std::vector<int64_t> numberParallels = {1};
    std::string coreShapeName = "RM 10/I";
    std::vector<OpenMagnetics::Wire> wires;
    
    auto litzWire = OpenMagnetics::Wire::create_quick_litz_wire(0.0001, 50);
    wires.push_back(litzWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::DUERDOTH);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK(capacitance > 0);
        }
    }
}

// ============================================================================
// FOIL WIRE STRAY CAPACITANCE TESTS
// ============================================================================

TEST_CASE("Calculate capacitance of foil windings with 4 turns", "[physical-model][stray-capacitance][foil][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {4, 4};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string coreShapeName = "PQ 32/20";
    std::vector<OpenMagnetics::Wire> wires;
    
    // Create a foil wire using the rectangular helper, then override type to FOIL
    auto foilWire = OpenMagnetics::Wire::create_quick_rectangular_wire(0.010, 0.0002);
    foilWire.set_type(WireType::FOIL);
    wires.push_back(foilWire);
    wires.push_back(foilWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            if (firstKey == secondKey) {
                CHECK(capacitance > 0);
            } else {
                CHECK(capacitance < 0);
            }
        }
    }
}

TEST_CASE("Calculate capacitance of single foil winding with 4 turns", "[physical-model][stray-capacitance][foil][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {4};
    std::vector<int64_t> numberParallels = {1};
    std::string coreShapeName = "PQ 32/20";
    std::vector<OpenMagnetics::Wire> wires;
    
    auto foilWire = OpenMagnetics::Wire::create_quick_rectangular_wire(0.008, 0.0003);
    foilWire.set_type(WireType::FOIL);
    wires.push_back(foilWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            CHECK(capacitance > 0);
            CHECK(capacitance < 1e-9); // Sanity upper bound
        }
    }
}

TEST_CASE("Calculate capacitance of rectangular windings with 4 turns", "[physical-model][stray-capacitance][rectangular][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {4, 4};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string coreShapeName = "PQ 32/20";
    std::vector<OpenMagnetics::Wire> wires;
    
    auto rectWire = OpenMagnetics::Wire::create_quick_rectangular_wire(0.005, 0.001);
    wires.push_back(rectWire);
    wires.push_back(rectWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            if (firstKey == secondKey) {
                CHECK(capacitance > 0);
            } else {
                CHECK(capacitance < 0);
            }
        }
    }
}

// ============================================================================
// MIXED WIRE TYPE TESTS (ROUND + LITZ)
// ============================================================================

TEST_CASE("Calculate capacitance between round and litz windings", "[physical-model][stray-capacitance][mixed][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {8, 8};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string coreShapeName = "RM 10/I";
    std::vector<OpenMagnetics::Wire> wires;
    
    auto roundWire = find_wire_by_name("Round 1.00 - Grade 1");
    auto litzWire = OpenMagnetics::Wire::create_quick_litz_wire(0.0001, 50);
    wires.push_back(roundWire);
    wires.push_back(litzWire);

    auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    // Should not throw - mixed round+litz both go through the round-like path
    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix().value();
    CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);

    for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            if (firstKey == secondKey) {
                CHECK(capacitance > 0);
            } else {
                CHECK(capacitance < 0);
            }
        }
    }
}


// =================== Added Litz/Foil Tests ===================

TEST_CASE("Litz winding stray capacitance smoke test", "[stray-capacitance][litz]") {
    settings.reset();

    // Create two identical litz wires and a simple coil
    auto litzWire = OpenMagnetics::Wire::create_quick_litz_wire(0.0001, 50);
    std::vector<int64_t> numberTurns = {8, 8};
    std::vector<int64_t> numberParallels = {1, 1};
    auto coil = OpenMagnetics::Coil::create_quick_coil("RM 10/I", numberTurns, numberParallels, {litzWire, litzWire});

    StrayCapacitance strayCap(OpenMagnetics::StrayCapacitanceModels::ALBACH);
    auto out = strayCap.calculate_capacitance(coil);
    REQUIRE(out.get_maxwell_capacitance_matrix());
}

TEST_CASE("Foil winding stray capacitance smoke test", "[stray-capacitance][foil]") {
    settings.reset();

    auto foilWire = OpenMagnetics::Wire::create_quick_rectangular_wire(0.008, 0.0003);
    foilWire.set_type(WireType::FOIL);
    std::vector<int64_t> numberTurns = {4,4};
    std::vector<int64_t> numberParallels = {1,1};
    auto coil = OpenMagnetics::Coil::create_quick_coil("PQ 32/20", numberTurns, numberParallels, {foilWire, foilWire});

    StrayCapacitance strayCap(OpenMagnetics::StrayCapacitanceModels::ALBACH);
    auto out = strayCap.calculate_capacitance(coil);
    REQUIRE(out.get_maxwell_capacitance_matrix());
}

TEST_CASE("Mixed round + litz capacitance smoke test", "[stray-capacitance][mixed]") {
    settings.reset();

    auto roundWire = OpenMagnetics::find_wire_by_name("Round 1.00 - Grade 1");
    auto litzWire = OpenMagnetics::Wire::create_quick_litz_wire(0.0001, 50);
    std::vector<int64_t> numberTurns = {8,8};
    std::vector<int64_t> numberParallels = {1,1};
    auto coil = OpenMagnetics::Coil::create_quick_coil("RM 10/I", numberTurns, numberParallels, {roundWire, litzWire});

    StrayCapacitance strayCap(OpenMagnetics::StrayCapacitanceModels::ALBACH);
    auto out = strayCap.calculate_capacitance(coil);
    REQUIRE(out.get_maxwell_capacitance_matrix());
}

TEST_CASE("Bug: beta is NAN in capacitance calculation with processed coil", "[stray-capacitance][bug][beta-nan]") {
    settings.reset();
    
    // Load the test data from JSON file
    auto testDataPath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("testData").append("bug_capacitance_error.json");
    
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    
    json masJson;
    file >> masJson;
    file.close();
    
    // Create Magnetic object from the JSON
    auto magneticJson = masJson["magnetic"];
    OpenMagnetics::Magnetic magnetic(magneticJson);
    
    // Get the coil and operating point
    auto coil = magnetic.get_coil();
    auto inputsJson = masJson["inputs"];
    auto operatingPoints = inputsJson["operatingPoints"];
    REQUIRE(operatingPoints.size() > 0);
    auto operatingPointJson = operatingPoints[0];
    OperatingPoint operatingPoint(operatingPointJson);
    
    
    // Check the wire properties
    const auto& windings = coil.get_functional_description();
    for (size_t i = 0; i < windings.size(); ++i) {
        const auto& winding = windings[i];
        auto wireOpt = winding.get_wire();
        if (std::holds_alternative<OpenMagnetics::Wire>(wireOpt)) {
            const auto& wire = std::get<OpenMagnetics::Wire>(wireOpt);
            auto nameOpt = wire.get_name();
            auto outerDiamOpt = wire.get_outer_diameter();
            if (outerDiamOpt) {
            }
            auto coatingOpt = wire.get_coating();
            if (coatingOpt) {
            }
            auto condDiamOpt = wire.get_conducting_diameter();
            if (condDiamOpt) {
            }
            auto condAreaOpt = wire.get_conducting_area();
            if (condAreaOpt) {
            }
            // For LITZ wires, check strand
            if (wire.get_type() == WireType::LITZ) {
                try {
                    // Create non-const copy to call resolve_strand
                    auto wireCopy = wire;
                    auto strand = wireCopy.resolve_strand();
                    auto strandCoating = OpenMagnetics::Wire::resolve_coating(strand);
                    if (strandCoating && strandCoating->get_grade()) {
                    }
                } catch (const std::exception& e) {
                    FAIL_CHECK("resolve_strand/resolve_coating threw for winding " << i << ": " << e.what());
                }
            }
        }
    }

    // Repro point (beta was NAN): every model must produce a Maxwell capacitance
    // matrix whose entries are all finite.
    for (auto model : {OpenMagnetics::StrayCapacitanceModels::ALBACH,
                       OpenMagnetics::StrayCapacitanceModels::KOCH,
                       OpenMagnetics::StrayCapacitanceModels::MASSARINI,
                       OpenMagnetics::StrayCapacitanceModels::DUERDOTH}) {
        INFO("Model: " << magic_enum::enum_name(model));
        try {
            StrayCapacitance strayCapacitance(model);
            auto output = strayCapacitance.calculate_capacitance(coil);

            REQUIRE(output.get_maxwell_capacitance_matrix());
            auto matrix = output.get_maxwell_capacitance_matrix().value();
            CHECK(matrix.size() > 0);
            for (const auto& entry : matrix) {
                for (const auto& [winding1, innerMap] : entry.get_magnitude()) {
                    for (const auto& [winding2, capacitance] : innerMap) {
                        double value = OpenMagnetics::resolve_dimensional_values(capacitance);
                        INFO("Capacitance " << winding1 << " / " << winding2);
                        CHECK(std::isfinite(value));
                    }
                }
            }
        } catch (const std::exception& e) {
            FAIL_CHECK("calculate_capacitance threw for model "
                       << magic_enum::enum_name(model) << ": " << e.what());
        }
    }
}

TEST_CASE("Investigate capacitance calculation from bug_capacitance_error_2.json", "[physical-model][stray-capacitance][bug-investigation]") {
    settings.reset();
    
    // Load the MAS from the test data file
    auto testDataPath = get_test_data_path(std::source_location::current(), "bug_capacitance_error_2.json");
    
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    
    std::string jsonStr((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    
    auto masJson = nlohmann::json::parse(jsonStr);
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_json(masJson, mas);
    auto magnetic = mas.get_mutable_magnetic();
    auto coil = magnetic.get_mutable_coil();
    
    
    if (coil.get_turns_description()) {
    } else {
    }
    
    if (coil.get_layers_description()) {
    } else {
    }
    
    // Try to calculate capacitance
    
    StrayCapacitance strayCapacitance;
    
    try {
        auto output = strayCapacitance.calculate_capacitance(coil);
        
        
        if (output.get_maxwell_capacitance_matrix()) {
            auto matrix = output.get_maxwell_capacitance_matrix().value();
            for (const auto& entry : matrix) {
                for (const auto& [winding1, innerMap] : entry.get_magnitude()) {
                    for (const auto& [winding2, capacitance] : innerMap) {
                        (void)OpenMagnetics::resolve_dimensional_values(capacitance);
                    }
                }
            }
        }
        
        if (output.get_capacitance_matrix()) {
            auto matrix = output.get_capacitance_matrix().value();
            for (const auto& [winding1, innerMap] : matrix) {
                for (const auto& [winding2, scalarMatrix] : innerMap) {
                    // ScalarMatrixAtFrequency contains magnitude which is map<string, map<string, DimensionWithTolerance>>
                    for (const auto& [sub1, subMap] : scalarMatrix.get_magnitude()) {
                        for (const auto& [sub2, dimWithTol] : subMap) {
                            (void)OpenMagnetics::resolve_dimensional_values(dimWithTol);
                        }
                    }
                }
            }
        }
        
    } catch (const std::exception& e) {
        FAIL("Capacitance calculation failed: " << e.what());
    }
    
    REQUIRE(true);
}


TEST_CASE("Investigate negative C3 in bug_capacitance_error_3", "[physical-model][stray-capacitance][bug-investigation]") {
    settings.reset();
    
    auto testDataPath = get_test_data_path(std::source_location::current(), "bug_capacitance_error_3.json");
    
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    
    std::string jsonStr((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    
    auto masJson = nlohmann::json::parse(jsonStr);
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_json(masJson, mas);
    auto magnetic = mas.get_mutable_magnetic();
    auto coil = magnetic.get_mutable_coil();
    
    
    StrayCapacitance strayCapacitance;
    auto output = strayCapacitance.calculate_capacitance(coil);
    
    
    if (output.get_capacitance_matrix()) {
        auto matrix = output.get_capacitance_matrix().value();
        
        double C11 = 0, C22 = 0, C12 = 0;
        std::string w1, w2;
        
        for (const auto& [winding1, innerMap] : matrix) {
            for (const auto& [winding2, scalarMatrix] : innerMap) {
                auto magnitude = scalarMatrix.get_magnitude();
                double total = 0;
                for (const auto& [k1, subMap] : magnitude) {
                    for (const auto& [k2, dimWithTol] : subMap) {
                        total += OpenMagnetics::resolve_dimensional_values(dimWithTol);
                    }
                }
                
                if (w1.empty()) w1 = winding1;
                if (winding1 == w1 && winding2 == w1) C11 = total;
                if (winding1 != w1 && winding2 != w1) {
                    C22 = total;
                    w2 = winding1;
                }
                if (winding1 == w1 && winding2 != w1) C12 = total;
                
            }
        }
        
        
        [[maybe_unused]] double C1 = C11 + C12;
        [[maybe_unused]] double C2 = C22 + C12;
        double C3 = -C12;
        
        
        if (C3 < 0) {
        }
    }
    
    REQUIRE(true);
}
TEST_CASE("Calculate capacitance for toroidal core", "[physical-model][stray-capacitance][toroidal][smoke-test]") {
    settings.reset();

    std::vector<int64_t> numberTurns = {10};
    std::vector<int64_t> numberParallels = {1};
    std::vector<double> turnsRatios = {};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97";

    // Use SPREAD alignment for toroidal cores (like other toroidal tests)
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape,
                                                     interleavingLevel, sectionOrientation, layersOrientation,
                                                     turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, json::array(), numberStacks, coreMaterial);

    // Ensure coil is wound to generate turn coordinates and additionalCoordinates for toroidal
    coil.wind();

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);

    auto result = strayCapacitance.calculate_capacitance(coil);
    auto maxwellCapacitanceMatrix = result.get_maxwell_capacitance_matrix();

    // Check that we have a valid capacitance matrix
    REQUIRE(maxwellCapacitanceMatrix.has_value());
    REQUIRE(maxwellCapacitanceMatrix.value().size() > 0);

    // Print the capacitance values for debugging
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix.value()[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            (void)OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
        }
    }

    // Check that capacitance among turns is not empty
    auto capacitanceAmongTurns = result.get_capacitance_among_turns();
    
    // Print some capacitance values between turns
    int count = 0;
    for (auto [key, innerMap] : capacitanceAmongTurns.value()) {
        for (auto [key2, capacitance] : innerMap) {
            if (count++ < 20) {
            }
        }
    }

    // For a toroidal core with 10 turns, we should have multiple turn pairs
    // If the capacitance matrix is 0, it means get_surrounding_turns is not finding neighbors
    CHECK(capacitanceAmongTurns->size() > 0);

    // Verify that at least some capacitance values are non-zero
    bool hasNonZeroCapacitance = false;
    for (auto [firstKey, aux] : maxwellCapacitanceMatrix.value()[0].get_magnitude()) {
        for (auto [secondKey, capacitanceWithTolerance] : aux) {
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
            if (std::abs(capacitance) > 1e-15) {
                hasNonZeroCapacitance = true;
                break;
            }
        }
        if (hasNonZeroCapacitance) break;
    }

    CHECK(hasNonZeroCapacitance);
}

TEST_CASE("Investigate missing energy between Tertiary and Secondary turns", "[physical-model][stray-capacitance][bug][debug]") {
    auto testDataPath = get_test_data_path(std::source_location::current(), "bug_capacitance_error_4.json");
    
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    
    std::string jsonStr((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    
    auto masJson = nlohmann::json::parse(jsonStr);
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_json(masJson, mas);
    auto magnetic = mas.get_mutable_magnetic();
    auto coil = magnetic.get_mutable_coil();
    
    
    // Get turns description
    auto turns = coil.get_turns_description().value();
    
    // Find specific turns
    int tertiaryTurn0Idx = -1;
    int secondaryTurn2Idx = -1;
    
    for (size_t i = 0; i < turns.size(); ++i) {
        const auto& turn = turns[i];
        (void)turn;
        
        if (turn.get_winding() == "Tertiary" && turn.get_name().find("turn 0") != std::string::npos) {
            tertiaryTurn0Idx = i;
        }
        if (turn.get_winding() == "secondary" && turn.get_name().find("turn 2") != std::string::npos) {
            secondaryTurn2Idx = i;
        }
    }
    
    
    if (tertiaryTurn0Idx >= 0 && secondaryTurn2Idx >= 0) {
        const auto& turn1 = turns[tertiaryTurn0Idx];
        const auto& turn2 = turns[secondaryTurn2Idx];
        
        double x1 = turn1.get_coordinates()[0];
        double y1 = turn1.get_coordinates()[1];
        double x2 = turn2.get_coordinates()[0];
        double y2 = turn2.get_coordinates()[1];
        
        double dx1 = turn1.get_dimensions().value()[0];
        double dy1 = turn1.get_dimensions().value()[1];
        double dx2 = turn2.get_dimensions().value()[0];
        double dy2 = turn2.get_dimensions().value()[1];
        
        double maxDim1 = std::max(dx1, dy1);
        double maxDim2 = std::max(dx2, dy2);
        [[maybe_unused]] double minDimension = std::min(maxDim1, maxDim2);

        double centerDist = hypot(x2 - x1, y2 - y1);
        [[maybe_unused]] double surfaceGap = centerDist - maxDim1 / 2 - maxDim2 / 2;

    }

    // Now calculate capacitance and check energies
    StrayCapacitance strayCapacitance;
    auto output = strayCapacitance.calculate_capacitance(coil);
    
    // Find turn names first (needed for multiple checks)
    std::string tertiaryTurn0Name;
    std::string secondaryTurn2Name;
    for (const auto& turn : turns) {
        if (turn.get_winding() == "Tertiary" && turn.get_name().find("turn 0") != std::string::npos) {
            tertiaryTurn0Name = turn.get_name();
        }
        if (turn.get_winding() == "secondary" && turn.get_name().find("turn 2") != std::string::npos) {
            secondaryTurn2Name = turn.get_name();
        }
    }
    
    auto energyMap = output.get_electric_energy_among_turns();
    if (energyMap) {
        
        
        auto it1 = energyMap->find(tertiaryTurn0Name);
        if (it1 != energyMap->end()) {
            auto it2 = it1->second.find(secondaryTurn2Name);
            if (it2 != it1->second.end()) {
            } else {
            }
        } else {
        }
        
        // Print all energies for Tertiary turn 0
        if (it1 != energyMap->end()) {
            for ([[maybe_unused]] const auto& [otherTurn, energy] : it1->second) {
            }
        }
        
        // Also check capacitance
        auto capacitanceMap = output.get_capacitance_among_turns();
        if (capacitanceMap) {
            auto itCap = capacitanceMap->find(tertiaryTurn0Name);
            if (itCap != capacitanceMap->end()) {
                for ([[maybe_unused]] const auto& [otherTurn, cap] : itCap->second) {
                }
            }
        }
    }
    
    // Check voltage per turn
    auto voltagePerTurn = output.get_voltage_per_turn();
    if (voltagePerTurn) {
        for (size_t i = 0; i < turns.size() && i < voltagePerTurn->size(); ++i) {
        }
    }
    
    // Check voltage drops
    auto voltageDropMap = output.get_voltage_drop_among_turns();
    if (voltageDropMap) {
        auto itVolt = voltageDropMap->find(tertiaryTurn0Name);
        if (itVolt != voltageDropMap->end()) {
            bool foundNonZero = false;
            for (const auto& [otherTurn, voltage] : itVolt->second) {
                if (voltage != 0) foundNonZero = true;
            }
            if (!foundNonZero) {
            }
        } else {
        }
    }
    
    // Check surrounding turns
    if (tertiaryTurn0Idx >= 0) {
        auto surrounding = StrayCapacitance::get_surrounding_turns(turns[tertiaryTurn0Idx], turns, -1.0);
        for ([[maybe_unused]] const auto& [surroundingTurn, idx] : surrounding) {
        }
    }
    
    // Check if turns are in the same or different windings
    if (tertiaryTurn0Idx >= 0 && secondaryTurn2Idx >= 0) {
    }
    
    REQUIRE(true);
}

TEST_CASE("Debug primary-primary intrawinding capacitance", "[physical-model][stray-capacitance][debug][intrawinding]") {
    auto testDataPath = get_test_data_path(std::source_location::current(), "bug_capacitance_error_4.json");
    
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    
    std::string jsonStr((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    
    auto masJson = nlohmann::json::parse(jsonStr);
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_json(masJson, mas);
    auto magnetic = mas.get_mutable_magnetic();
    auto coil = magnetic.get_mutable_coil();
    
    auto turns = coil.get_turns_description().value();
    
    // Find primary turns 0 and 1 (should be adjacent in layer 0)
    const auto& turn0 = turns[0];  // primary parallel 0 turn 0
    const auto& turn1 = turns[1];  // primary parallel 0 turn 1
    
    
    double x0 = turn0.get_coordinates()[0];
    double y0 = turn0.get_coordinates()[1];
    double x1 = turn1.get_coordinates()[0];
    double y1 = turn1.get_coordinates()[1];
    
    double dx0 = turn0.get_dimensions().value()[0];
    double dy0 = turn0.get_dimensions().value()[1];
    double dx1 = turn1.get_dimensions().value()[0];
    double dy1 = turn1.get_dimensions().value()[1];
    
    double maxDim0 = std::max(dx0, dy0);
    double maxDim1 = std::max(dx1, dy1);
    [[maybe_unused]] double minDimension = std::min(maxDim0, maxDim1);

    double centerDist = hypot(x1 - x0, y1 - y0);
    [[maybe_unused]] double surfaceGap = centerDist - maxDim0 / 2 - maxDim1 / 2;
    
    
    // Check surrounding turns for turn 0
    auto surroundingForTurn0 = StrayCapacitance::get_surrounding_turns(turn0, turns, -1.0);
    for ([[maybe_unused]] const auto& [surroundingTurn, idx] : surroundingForTurn0) {
    }
    
    // Check if turn 1 is in the surrounding turns
    [[maybe_unused]] bool foundTurn1 = false;
    for (const auto& [surroundingTurn, idx] : surroundingForTurn0) {
        if (idx == 1) {
            foundTurn1 = true;
            break;
        }
    }
    
    // Calculate capacitance and check
    StrayCapacitance strayCapacitance;
    auto output = strayCapacitance.calculate_capacitance(coil);
    
    auto capacitanceAmongWindings = output.get_capacitance_among_windings();
    if (capacitanceAmongWindings) {
        for ([[maybe_unused]] const auto& [winding1, inner] : *capacitanceAmongWindings) {
            for ([[maybe_unused]] const auto& [winding2, cap] : inner) {
            }
        }
    }
    
    // Check capacitance between primary turn 0 and turn 1
    auto capacitanceMap = output.get_capacitance_among_turns();
    if (capacitanceMap) {
        auto it = capacitanceMap->find(turn0.get_name());
        if (it != capacitanceMap->end()) {
            auto it2 = it->second.find(turn1.get_name());
            if (it2 != it->second.end()) {
            } else {
            }
        }
    }
    
    // Output full JSON for debugging
    nlohmann::json result;
    to_json(result, output);
    if (result.contains("capacitanceAmongWindings")) {
    }
    
    // Check Maxwell capacitance matrix
    if (result.contains("maxwellCapacitanceMatrix")) {
    }
    
    REQUIRE(true);
}

TEST_CASE("Turn_To_Core_Capacitance_Bounded_And_Monotone", "[physical-model][stray-capacitance][coating]") {
    // §7 building block: capacitance of a turn to the equipotential ferrite surface,
    // through the wire-enamel | air | core-coating dielectric stack. It must be FINITE
    // at zero air gap (the coating is the floor — no divergence like the old
    // wire-above-plane helper) and decrease monotonically as the air gap grows.
    const double r = 0.4e-3;        // conducting radius [m]
    const double L = 30e-3;         // turn length [m]
    const double tEnamel = 25e-6;   // wire enamel thickness [m]
    const double epsEnamel = 3.5;
    const double tCoating = 12.7e-6; // parylene core coating [m]
    const double epsCoating = 3.1;

    std::vector<double> gaps = {0.0, 10e-6, 50e-6, 100e-6, 300e-6, 1e-3};
    std::vector<double> caps;
    for (double g : gaps) {
        double c = StrayCapacitance::calculate_turn_to_core_capacitance(r, L, tEnamel, epsEnamel, g, tCoating, epsCoating);
        caps.push_back(c);
        REQUIRE(std::isfinite(c));
        REQUIRE(c > 0.0);
        REQUIRE(c < 1e-9);  // bounded far below the old ~1.3 nF divergence
    }

    // Finite (bounded by the coating) at zero air gap — the divergence is gone.
    REQUIRE(std::isfinite(caps.front()));
    REQUIRE(caps.front() < 1e-9);

    // Monotonic decrease with increasing air gap.
    for (size_t i = 1; i < caps.size(); ++i) {
        REQUIRE(caps[i] < caps[i - 1]);
    }

    // A thicker core coating lowers the contact (zero-gap) capacitance.
    double cThinCoat = StrayCapacitance::calculate_turn_to_core_capacitance(r, L, tEnamel, epsEnamel, 0.0, 12.7e-6, epsCoating);
    double cThickCoat = StrayCapacitance::calculate_turn_to_core_capacitance(r, L, tEnamel, epsEnamel, 0.0, 100e-6, epsCoating);
    REQUIRE(cThickCoat < cThinCoat);

    // Degenerate inputs return 0 rather than NaN/Inf.
    REQUIRE(StrayCapacitance::calculate_turn_to_core_capacitance(0.0, L, tEnamel, epsEnamel, 0.0, tCoating, epsCoating) == 0.0);
}

// ABT #31 repro: sweep_impedance_over_frequency -> "capacitance cannot be NaN".
// Reproduced with the planar flyback MAS the user reported (the same failure also
// hits the WE CMC catalogue via the El Choker tool).
TEST_CASE("ABT31 sweep_impedance_over_frequency capacitance NaN repro", "[physical-model][stray-capacitance][abt31]") {
    settings.reset();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "custom_magnetic_flyback.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();

    for (auto title : {std::string("Impedance over frequency"), std::string("Common-mode impedance")}) {
        try {
            auto sweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1e3, 1e9, 100, "log", title);
            bool finite = true;
            for (auto y : sweep.get_y_points()) if (!std::isfinite(y)) finite = false;
            UNSCOPED_INFO(title << " => sweep finite=" << finite);
            CHECK(finite);
        }
        catch (const std::exception& e) {
            UNSCOPED_INFO(title << " THREW: " << e.what());
            FAIL_CHECK(title << " threw");
        }
    }
}

// Regression for a NaN report on the stray capacitance of a 3-winding PLANAR
// flyback (custom_magnetic_flyback.json). The coil-only and operating-point
// paths must return finite, real capacitances for planar wires (parallel-plate
// model). See the companion "(through-core)" case for the with-core path.
TEST_CASE("Calculate stray capacitance of custom planar flyback", "[physical-model][stray-capacitance]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "custom_magnetic_flyback.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();
    auto inputs = mas.get_inputs();
    auto operatingPoint = inputs.get_operating_point(0);

    StrayCapacitance strayCapacitance;

    // Raw turn-to-turn capacitances (planar parallel-plate model) must be finite.
    auto capacitanceAmongTurns = strayCapacitance.calculate_capacitance_among_turns(coil);
    REQUIRE_FALSE(capacitanceAmongTurns.empty());
    for (auto& [turnsKey, capacitance] : capacitanceAmongTurns) {
        UNSCOPED_INFO("C[turn " << turnsKey.first << "][turn " << turnsKey.second << "] = " << capacitance);
        REQUIRE(std::isfinite(capacitance));
        REQUIRE(capacitance >= 0);
    }

    auto checkMaxwell = [](const std::string& label, const StrayCapacitanceOutput& output) {
        auto maxwellCapacitanceMatrix = output.get_maxwell_capacitance_matrix().value();
        REQUIRE_FALSE(maxwellCapacitanceMatrix.empty());
        for (auto& matrixAtFreq : maxwellCapacitanceMatrix) {
            for (auto& [firstKey, row] : matrixAtFreq.get_magnitude()) {
                for (auto& [secondKey, capacitanceWithTolerance] : row) {
                    auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
                    UNSCOPED_INFO(label << " C[" << firstKey << "][" << secondKey << "] = " << capacitance);
                    REQUIRE(std::isfinite(capacitance));
                }
            }
        }
    };

    checkMaxwell("coil", strayCapacitance.calculate_capacitance(coil));
    checkMaxwell("coil+operatingPoint", strayCapacitance.calculate_capacitance(coil, operatingPoint));
}

// The with-core (through-core) inter-winding path now supports planar/foil/
// rectangular wires (flat conductors modelled as an area-equivalent cylinder in
// the turn-to-core field-spreading integral), so separated planar windings + core
// must return a finite capacitance matrix instead of throwing.
TEST_CASE("Calculate stray capacitance of custom planar flyback (through-core)", "[physical-model][stray-capacitance]") {
    settings.reset();

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "custom_magnetic_flyback.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();
    auto core = magnetic.get_core();

    StrayCapacitance strayCapacitance;
    StrayCapacitanceOutput output;
    REQUIRE_NOTHROW(output = strayCapacitance.calculate_capacitance(coil, core));

    auto maxwell = output.get_maxwell_capacitance_matrix().value();
    REQUIRE_FALSE(maxwell.empty());
    for (auto& matrixAtFreq : maxwell) {
        for (auto& [firstKey, row] : matrixAtFreq.get_magnitude()) {
            for (auto& [secondKey, capacitanceWithTolerance] : row) {
                auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
                UNSCOPED_INFO("C[" << firstKey << "][" << secondKey << "] = " << capacitance);
                REQUIRE(std::isfinite(capacitance));
            }
        }
    }

    // The through-core element for a flat conductor must be a finite, positive, and
    // physically small (sub-nF) per-winding self term.
    auto amongWindings = output.get_capacitance_among_windings().value();
    for (auto& [firstName, row] : amongWindings) {
        for (auto& [secondName, capacitance] : row) {
            UNSCOPED_INFO("Cww[" << firstName << "][" << secondName << "] = " << capacitance);
            REQUIRE(std::isfinite(capacitance));
        }
    }
}

// ABT #173: the boundary turns of a separated-winding toroidal CMC (where the two
// windings face each other across the separation gap) got a NEGATIVE static
// capacitance (-0.57 pF on WE 744834622): the insulation layers credited between
// those turns (4.7 mm, collected along the winding path) exceeded the actual
// 2.8 mm gap, driving distanceThroughAir negative, and the Albach second-order
// correction — wrongly fed the layer stack instead of the wire coating thickness —
// dominated with Z < 0. A negative pair capacitance propagated into a negative
// among-windings off-diagonal, which silently killed the DM resonance and turned
// the wideband leakage tank into a bare series inductor.
TEST_CASE("Static turn capacitances are non-negative on a separated-winding toroidal CMC", "[physical-model][stray-capacitance]") {
    settings.reset();
    auto testDataPath = get_test_data_path(std::source_location::current(), "cmc_redexpert_744834622.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    auto magneticJson = nlohmann::json::parse(file);
    OpenMagnetics::Magnetic magnetic(magneticJson);
    magnetic = magnetic_autocomplete(magnetic);
    auto coil = magnetic.get_coil();
    auto core = magnetic.get_core();
    if (!coil.get_turns_description()) {
        coil.wind();
    }

    auto amongTurns = OpenMagnetics::StrayCapacitance().calculate_capacitance_among_turns(coil);
    for (auto& [key, capacitance] : amongTurns) {
        UNSCOPED_INFO("C(turn " << key.first << ", turn " << key.second << ") = " << capacitance);
        REQUIRE(capacitance >= 0.0);
    }

    auto amongWindings = OpenMagnetics::StrayCapacitance().calculate_capacitance(coil, core).get_capacitance_among_windings().value();
    auto w0 = coil.get_functional_description()[0].get_name();
    auto w1 = coil.get_functional_description()[1].get_name();
    // The adjacent-winding off-diagonal must be a physical, positive lumped capacitance
    // (it damps/places the DM resonance and the wideband leakage tank).
    CHECK(amongWindings[w0][w1] > 0.0);
    CHECK(amongWindings[w0][w1] < 1e-9);
}

// ABT #366: shielded drum (drumRing). The stray-capacitance pipeline must handle a coil wound
// in the drum groove of the new family end-to-end: real wire, autocompleted quick bobbin,
// finite positive Maxwell self-capacitance in a physically sane range for a millimetre part.
TEST_CASE("Test_Stray_Capacitance_Drum_Ring_Smoke", "[physical-model][stray-capacitance][drum-ring]") {
    settings.reset();
    clear_databases();
    auto core = OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "3C90");
    json coilJson;
    coilJson["bobbin"] = "Dummy";
    coilJson["functionalDescription"] = json::array();
    json winding;
    winding["name"] = "winding 0";
    winding["numberTurns"] = 8;
    winding["numberParallels"] = 1;
    winding["isolationSide"] = "primary";
    winding["wire"] = "Round 0.1 - Grade 1";
    coilJson["functionalDescription"].push_back(winding);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
    auto completed = OpenMagnetics::magnetic_autocomplete(magnetic);
    auto coil = completed.get_coil();
    REQUIRE(coil.get_turns_description().has_value());

    auto output = StrayCapacitance().calculate_capacitance(coil);
    REQUIRE(output.get_capacitance_among_windings().has_value());
    auto amongWindings = output.get_capacitance_among_windings().value();
    auto windingName = coil.get_functional_description()[0].get_name();
    double selfCapacitance = amongWindings.at(windingName).at(windingName);
    CHECK(std::isfinite(selfCapacitance));
    CHECK(selfCapacitance > 0);
    CHECK(selfCapacitance < 1e-10);  // a 2.3 mm part cannot carry >100 pF of self-capacitance
    settings.reset();
}

// ABT #366/#362/#357: stray capacitance across ALL FOUR new families. The pipeline has to wind
// a real coil in each family's window (drum groove or molded cavity), then build the Maxwell
// matrix from those turn positions — so this exercises geometry, winding and capacitance
// together. Physical bound rather than a pinned value: a millimetre-scale part cannot carry
// more than ~100 pF of self-capacitance, and it must be strictly positive.
TEST_CASE("Test_Stray_Capacitance_New_Core_Families",
          "[physical-model][stray-capacitance][drum][drum-ring][drum-semishielded][molded]") {
    settings.reset();
    clear_databases();

    auto buildCustomCore = [](json shapeJson, const std::string& coreType, json coating,
                              const std::string& materialName) {
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", coreType}, {"material", materialName}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        if (!coating.is_null()) {
            coreJson["functionalDescription"]["coating"] = coating;
        }
        OpenMagnetics::Core core(coreJson);
        core.process_data();
        core.process_gap();
        return core;
    };
    json drumDimensions = {
        {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
        {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}}};
    json semishieldedDimensions = drumDimensions;
    semishieldedDimensions["J"] = {{"nominal", 0.0040}};
    semishieldedDimensions["K"] = {{"nominal", 0.0040}};
    semishieldedDimensions["L"] = {{"nominal", 0.0018}};

    std::vector<std::pair<std::string, OpenMagnetics::Core>> cores;
    cores.emplace_back("drum", OpenMagneticsTesting::get_quick_core("DRH-14X20-4C", json::array(), 1, "3C90"));
    cores.emplace_back("drumRing", OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "3C90"));
    cores.emplace_back("drumSemishielded", buildCustomCore(
        {{"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "drumSemishielded"},
         {"aliases", json::array()}, {"name", "LQS-like 4018"}, {"dimensions", semishieldedDimensions}},
        "pieceAndPlate", {{"type", "magneticEpoxy"}, {"thickness", 0.0001}, {"material", "Kool M\u00b5 26"}}, "3C90"));
    cores.emplace_back("molded", buildCustomCore(
        {{"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
         {"aliases", json::array()}, {"name", "MAPI-like 4020"},
         {"dimensions", {
             {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
             {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}},
        "closedShape", json(), "Kool M\u00b5 26"));

    for (auto& [label, core] : cores) {
        json coilJson;
        coilJson["bobbin"] = "Dummy";
        coilJson["functionalDescription"] = json::array({{
            {"name", "winding 0"}, {"numberTurns", 8}, {"numberParallels", 1},
            {"isolationSide", "primary"}, {"wire", "Round 0.1 - Grade 1"}}});
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
        auto completed = OpenMagnetics::magnetic_autocomplete(magnetic);
        auto coil = completed.get_coil();
        REQUIRE(coil.get_turns_description().has_value());

        auto output = StrayCapacitance().calculate_capacitance(coil);
        REQUIRE(output.get_capacitance_among_windings().has_value());
        auto amongWindings = output.get_capacitance_among_windings().value();
        auto windingName = coil.get_functional_description()[0].get_name();
        double selfCapacitance = amongWindings.at(windingName).at(windingName);
        UNSCOPED_INFO(label << ": self-capacitance " << selfCapacitance * 1e12 << " pF");
        CHECK(std::isfinite(selfCapacitance));
        CHECK(selfCapacitance > 0);
        CHECK(selfCapacitance < 1e-10);
    }
    settings.reset();
}

// ABT #395: stray capacitance went NaN for ANY winding with numberParallels > 1, which threw out
// of export_magnetic_as_subcircuit and blocked the SPICE path for every parallel-wound part (119 of
// 504 magnetics in the reporter's corpus; the whole LHMI family). The reported symptom was a
// perfect separation on numberParallels with turn count irrelevant, and this pins the mechanism
// behind that separation.
//
// The energy of a turn pair is 0.5*C*dV^2. dV is EXACTLY zero between corresponding turns of two
// parallels, because the per-turn voltage divider indexes within each parallel and therefore hands
// them identical potentials -- they are the same electrical node. So an infinite C, which a single
// parallel would merely propagate as an infinity, became 0.5*inf*0 = NaN as soon as a second
// parallel existed. Hence parallels, not turns, was the discriminator.
TEST_CASE("Corresponding turns of different parallels sit at identical potential", "[physical-model][stray-capacitance][parallels]") {
    settings.reset();
    auto coreJsonStr = R"({"name": "abt395", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

    std::map<int64_t, size_t> pairsAtIdenticalPotentialPerParallels;
    for (int64_t numberParallels : {1, 2, 4}) {
        auto coilJsonStr = std::string(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": )")
            + std::to_string(numberParallels) + R"(, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })";
        auto [core, coil] = prepare_core_and_coil_from_json(coreJsonStr, coilJsonStr);
        auto voltagesPerTurn = StrayCapacitance::calculate_voltages_per_turn(coil, {{"Primary", 100.0}})
                                   .get_voltage_per_turn().value();
        REQUIRE(voltagesPerTurn.size() == static_cast<size_t>(4 * numberParallels));
        size_t pairsAtIdenticalPotential = 0;
        for (size_t firstTurn = 0; firstTurn < voltagesPerTurn.size(); ++firstTurn) {
            for (size_t secondTurn = firstTurn + 1; secondTurn < voltagesPerTurn.size(); ++secondTurn) {
                if (voltagesPerTurn[firstTurn] == voltagesPerTurn[secondTurn]) {
                    pairsAtIdenticalPotential++;
                }
            }
        }
        pairsAtIdenticalPotentialPerParallels[numberParallels] = pairsAtIdenticalPotential;
    }
    UNSCOPED_INFO("pairs at identical potential: 1 parallel -> " << pairsAtIdenticalPotentialPerParallels[1]
                  << ", 2 -> " << pairsAtIdenticalPotentialPerParallels[2]
                  << ", 4 -> " << pairsAtIdenticalPotentialPerParallels[4]);
    // A single parallel never puts two turns at the same potential, which is exactly why the bug
    // could not reach it. Each extra parallel replicates the whole divider.
    CHECK(pairsAtIdenticalPotentialPerParallels[1] == 0);
    CHECK(pairsAtIdenticalPotentialPerParallels[2] == 4);
    CHECK(pairsAtIdenticalPotentialPerParallels[4] == 24);
}

// The other half: an infinite turn-to-turn capacitance must be reported where the geometry is still
// known, naming the turns and the missing insulation, instead of escaping to become an anonymous
// "Energy cannot be nan" for the whole magnetic. Bare wire (outer diameter == conducting diameter,
// so zero coating) wound tight puts two conductors in surface contact: no dielectric, no gap, and
// every model in the file diverges. That is inconsistent input, not something to approximate.
TEST_CASE("Touching bare conductors report the missing insulation, not a NaN energy", "[physical-model][stray-capacitance][parallels]") {
    settings.reset();
    auto coreJsonStr = R"({"name": "abt395", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";
    std::string bareWireJsonStr = R"({"type": "round", "material": "copper", "conductingDiameter": {"nominal": 0.001}, "outerDiameter": {"nominal": 0.001}, "numberConductors": 1})";

    for (int64_t numberParallels : {1, 2}) {
        auto coilJsonStr = std::string(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": )")
            + std::to_string(numberParallels) + R"(, "isolationSide": "primary", "wire": )" + bareWireJsonStr + R"( } ] })";
        auto magnetic = prepare_magnetic_from_json(coreJsonStr, coilJsonStr);
        OpenMagnetics::CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::NGSPICE);
        std::string message;
        try {
            exporter.export_magnetic_as_subcircuit(magnetic, 100000, 25);
            message = "no exception";
        }
        catch (const std::exception& exception) {
            message = exception.what();
        }
        UNSCOPED_INFO("parallels " << numberParallels << ": " << message);
        // Both parallel counts must report the same root cause. Note that 1 parallel used to
        // "succeed" here: the infinity was carried into the exported netlist rather than raising
        // anything, so this also pins that a non-finite capacitance can no longer ship silently.
        CHECK(message.find("Turn-to-turn capacitance") != std::string::npos);
        CHECK(message.find("zero coating thickness") != std::string::npos);
    }
}

// Insulated wire is the normal case and must keep working with parallels: same core, same turns,
// only numberParallels varying, all exporting a finite self-capacitance.
TEST_CASE("Calculate capacitance of a winding with 4 turns and several parallels", "[physical-model][stray-capacitance][parallels]") {
    settings.reset();
    auto coreJsonStr = R"({"name": "abt395", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

    for (int64_t numberParallels : {1, 2, 4}) {
        auto coilJsonStr = std::string(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": )")
            + std::to_string(numberParallels) + R"(, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })";
        auto [core, coil] = prepare_core_and_coil_from_json(coreJsonStr, coilJsonStr);
        auto maxwellCapacitanceMatrix = StrayCapacitance().calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
        REQUIRE(maxwellCapacitanceMatrix[0].get_magnitude().size() == 1);
        for (auto [firstWindingName, capacitancesToOtherWindings] : maxwellCapacitanceMatrix[0].get_magnitude()) {
            for (auto [secondWindingName, capacitanceWithTolerance] : capacitancesToOtherWindings) {
                auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
                UNSCOPED_INFO(numberParallels << " parallels: " << firstWindingName << "/" << secondWindingName
                              << " capacitance " << capacitance);
                CHECK(std::isfinite(capacitance));
            }
        }

        auto magnetic = prepare_magnetic_from_json(coreJsonStr, coilJsonStr);
        OpenMagnetics::CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::NGSPICE);
        std::string subcircuit;
        REQUIRE_NOTHROW(subcircuit = exporter.export_magnetic_as_subcircuit(magnetic, 100000, 25));
        CHECK(subcircuit.size() > 0);
    }
}
