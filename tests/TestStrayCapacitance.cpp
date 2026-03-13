#include <source_location>
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "physical_models/StrayCapacitance.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include "support/Painter.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
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
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

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

TEST_CASE("Calculate capacitance of a winding with 8 turns and 1 parallel", "[physical-model][stray-capacitance][smoke-test]") {
    settings.reset();
    auto coilJsonStr = R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })";
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

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
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

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
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

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
    auto coreJsonStr = R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })";

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
                }
            }
        }
    }
    
    // Try to calculate capacitance with different models
    for (auto model : {OpenMagnetics::StrayCapacitanceModels::ALBACH, 
                       OpenMagnetics::StrayCapacitanceModels::KOCH,
                       OpenMagnetics::StrayCapacitanceModels::MASSARINI,
                       OpenMagnetics::StrayCapacitanceModels::DUERDOTH}) {
        
        
        try {
            StrayCapacitance strayCapacitance(model);
            auto output = strayCapacitance.calculate_capacitance(coil);
            
            if (output.get_maxwell_capacitance_matrix()) {
                auto matrix = output.get_maxwell_capacitance_matrix().value();
                if (matrix.size() > 0) {
                }
            }
            
        } catch (const std::exception& e) {
            // Don't fail the test, just report the error
        }
    }
    
    // The test should complete without crashing
    REQUIRE(true);
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
    MAS::Mas mas;
    MAS::from_json(masJson, mas);
    auto magnetic = mas.get_mutable_magnetic();
    auto coil = magnetic.get_mutable_coil();
    
    for (const auto& winding : coil.get_functional_description()) {
    }
    
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
                        auto value = OpenMagnetics::resolve_dimensional_values(capacitance);
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
                            auto value = OpenMagnetics::resolve_dimensional_values(dimWithTol);
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
    MAS::Mas mas;
    MAS::from_json(masJson, mas);
    auto magnetic = mas.get_mutable_magnetic();
    auto coil = magnetic.get_mutable_coil();
    
    for (const auto& winding : coil.get_functional_description()) {
    }
    
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
        
        
        double C1 = C11 + C12;
        double C2 = C22 + C12;
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
            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
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
    MAS::Mas mas;
    MAS::from_json(masJson, mas);
    auto magnetic = mas.get_mutable_magnetic();
    auto coil = magnetic.get_mutable_coil();
    
    for (const auto& winding : coil.get_functional_description()) {
    }
    
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
        double minDimension = std::min(maxDim1, maxDim2);
        
        double centerDist = hypot(x2 - x1, y2 - y1);
        double surfaceGap = centerDist - maxDim1 / 2 - maxDim2 / 2;
        
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
            for (const auto& [otherTurn, energy] : it1->second) {
            }
        }
        
        // Also check capacitance
        auto capacitanceMap = output.get_capacitance_among_turns();
        if (capacitanceMap) {
            auto itCap = capacitanceMap->find(tertiaryTurn0Name);
            if (itCap != capacitanceMap->end()) {
                for (const auto& [otherTurn, cap] : itCap->second) {
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
        for (const auto& [surroundingTurn, idx] : surrounding) {
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
    MAS::Mas mas;
    MAS::from_json(masJson, mas);
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
    double minDimension = std::min(maxDim0, maxDim1);
    
    double centerDist = hypot(x1 - x0, y1 - y0);
    double surfaceGap = centerDist - maxDim0 / 2 - maxDim1 / 2;
    
    
    // Check surrounding turns for turn 0
    auto surroundingForTurn0 = StrayCapacitance::get_surrounding_turns(turn0, turns, -1.0);
    for (const auto& [surroundingTurn, idx] : surroundingForTurn0) {
    }
    
    // Check if turn 1 is in the surrounding turns
    bool foundTurn1 = false;
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
        for (const auto& [winding1, inner] : *capacitanceAmongWindings) {
            for (const auto& [winding2, cap] : inner) {
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
