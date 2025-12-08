#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "physical_models/StrayCapacitance.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include "support/Painter.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

static double maximumError = 0.4;
static auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
static bool plot = true;

// TEST_CASE("Calculate capacitance between two bare round turns", "[physical-model][stray-capacitance]") {

//     settings->reset();
//     auto firstWire = OpenMagnetics::find_wire_by_name("Round 2.00 - Grade 1");
//     auto secondWire = OpenMagnetics::find_wire_by_name("Round 2.00 - Grade 1");
//     firstWire.set_bare_coating();
//     secondWire.set_bare_coating();

//     MAS::Turn firstTurn;
//     MAS::Turn secondTurn;
//     firstTurn.set_coordinates({0, 0.0015, 0});
//     firstTurn.set_length(1);
//     secondTurn.set_coordinates({0, -0.001, 0});
//     secondTurn.set_length(1);

//     double expectedValue = 35.4e-12;

//     SECTION( "tests Albach method" ) {
//         StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);
//         auto staticCapacitanceBetweenTwoTurns = strayCapacitance.calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire);
//         std::cout << "ALBACH staticCapacitanceBetweenTwoTurns: " << staticCapacitanceBetweenTwoTurns << std::endl;
//         CHECK_THAT(staticCapacitanceBetweenTwoTurns, WithinRel(expectedValue, maximumError));
//     }
//     SECTION( "tests Koch method" ) {
//         StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::KOCH);
//         auto staticCapacitanceBetweenTwoTurns = strayCapacitance.calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire);
//         std::cout << "KOCH staticCapacitanceBetweenTwoTurns: " << staticCapacitanceBetweenTwoTurns << std::endl;
//         CHECK_THAT(staticCapacitanceBetweenTwoTurns, WithinRel(expectedValue, maximumError));
//     }
//     SECTION( "tests Duerdoth method" ) {
//         StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::DUERDOTH);
//         auto staticCapacitanceBetweenTwoTurns = strayCapacitance.calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire);
//         std::cout << "DUERDOTH staticCapacitanceBetweenTwoTurns: " << staticCapacitanceBetweenTwoTurns << std::endl;
//         CHECK_THAT(staticCapacitanceBetweenTwoTurns, WithinRel(expectedValue, maximumError));
//     }
//     SECTION( "tests Massarini method" ) {
//         StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::MASSARINI);
//         auto staticCapacitanceBetweenTwoTurns = strayCapacitance.calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire);
//         std::cout << "MASSARINI staticCapacitanceBetweenTwoTurns: " << staticCapacitanceBetweenTwoTurns << std::endl;
//         CHECK_THAT(staticCapacitanceBetweenTwoTurns, WithinRel(expectedValue, maximumError));
//     }
// }

// TEST_CASE("Calculate capacitance between two enamelled round turns", "[physical-model][stray-capacitance]") {
//     settings->reset();
//     auto firstWire = OpenMagnetics::find_wire_by_name("Round 2.00 - Grade 1");
//     auto secondWire = OpenMagnetics::find_wire_by_name("Round 2.00 - Grade 1");

//     MAS::Turn firstTurn;
//     MAS::Turn secondTurn;
//     firstTurn.set_coordinates({0, 0.0015, 0});
//     firstTurn.set_length(1);
//     secondTurn.set_coordinates({0, -0.001, 0});
//     secondTurn.set_length(1);

//     double expectedValue = 35.4e-12;

//     SECTION( "tests Albach method" ) {
//         StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::ALBACH);
//         auto staticCapacitanceBetweenTwoTurns = strayCapacitance.calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire);
//         std::cout << "ALBACH staticCapacitanceBetweenTwoTurns: " << staticCapacitanceBetweenTwoTurns << std::endl;
//         CHECK_THAT(staticCapacitanceBetweenTwoTurns, WithinRel(expectedValue, maximumError));
//     }
//     SECTION( "tests Koch method" ) {
//         StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::KOCH);
//         auto staticCapacitanceBetweenTwoTurns = strayCapacitance.calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire);
//         std::cout << "KOCH staticCapacitanceBetweenTwoTurns: " << staticCapacitanceBetweenTwoTurns << std::endl;
//         CHECK_THAT(staticCapacitanceBetweenTwoTurns, WithinRel(expectedValue, maximumError));
//     }
//     SECTION( "tests Duerdoth method" ) {
//         StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::DUERDOTH);
//         auto staticCapacitanceBetweenTwoTurns = strayCapacitance.calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire);
//         std::cout << "DUERDOTH staticCapacitanceBetweenTwoTurns: " << staticCapacitanceBetweenTwoTurns << std::endl;
//         CHECK_THAT(staticCapacitanceBetweenTwoTurns, WithinRel(expectedValue, maximumError));
//     }
//     SECTION( "tests Massarini method" ) {
//         StrayCapacitance strayCapacitance(OpenMagnetics::StrayCapacitanceModels::MASSARINI);
//         auto staticCapacitanceBetweenTwoTurns = strayCapacitance.calculate_static_capacitance_between_two_turns(firstTurn, firstWire, secondTurn, secondWire);
//         std::cout << "MASSARINI staticCapacitanceBetweenTwoTurns: " << staticCapacitanceBetweenTwoTurns << std::endl;
//         CHECK_THAT(staticCapacitanceBetweenTwoTurns, WithinRel(expectedValue, maximumError));
//     }
// }

TEST_CASE("Calculate capacitance among two windings each with 1 turn and 1 parallel", "[physical-model][stray-capacitance]") {
    settings->reset();
    json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
    json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    auto core = OpenMagnetics::Core(coreJson);
    auto coil = OpenMagnetics::Coil(coilJson);

    core.process_data();
    core.process_gap();
    coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
    coil.wind();

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
            std::cout << "-------------------: " << std::endl;
            std::cout << "firstKey: " << firstKey << std::endl;
            std::cout << "secondKey: " << secondKey << std::endl;
            std::cout << "capacitance: " << capacitance << std::endl;
            std::cout << "expectedValues[firstKey][secondKey]: " << expectedValues[firstKey][secondKey] << std::endl;
            CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
        }
    }
}

TEST_CASE("Calculate capacitance of a winding with 8 turns and 1 parallel", "[physical-model][stray-capacitance]") {
    settings->reset();
    json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })");
    json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    auto core = OpenMagnetics::Core(coreJson);
    auto coil = OpenMagnetics::Coil(coilJson);

    core.process_data();
    core.process_gap();
    coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
    coil.wind();

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

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);


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

TEST_CASE("Calculate capacitance among two windings each with 8 turns and 1 parallel", "[physical-model][stray-capacitance]") {
    settings->reset();
    json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
    json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    auto core = OpenMagnetics::Core(coreJson);
    auto coil = OpenMagnetics::Coil(coilJson);

    core.process_data();
    core.process_gap();
    coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
    coil.wind();

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

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);


    auto outFile = outputFilePath;
    outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels.svg");
    std::filesystem::remove(outFile);
    Painter painter(outFile);
    painter.paint_core(magnetic);
    painter.paint_bobbin(magnetic);
    painter.paint_coil_turns(magnetic);
    painter.export_svg();

}
TEST_CASE("Calculate capacitance among two windings one with 16 and another with 8 turns, and both 1 parallel", "[physical-model][stray-capacitance]") {
    settings->reset();
    json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 16, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
    json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    auto core = OpenMagnetics::Core(coreJson);
    auto coil = OpenMagnetics::Coil(coilJson);

    core.process_data();
    core.process_gap();
    coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
    coil.wind();

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

TEST_CASE("Calculate capacitance among three windings each with 8 turns and 1 parallel", "[physical-model][stray-capacitance]") {
    settings->reset();
    json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" }, {"name": "Tertiary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "tertiary", "wire": "Round 1.00 - Grade 1" } ] })");
    json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    auto core = OpenMagnetics::Core(coreJson);
    auto coil = OpenMagnetics::Coil(coilJson);

    core.process_data();
    core.process_gap();
    coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
    coil.wind();

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

TEST_CASE("Calculate and compare capacitance matrix of one winding versus two windings", "[physical-model][stray-capacitance]") {
    {
        std::vector<int64_t> numberTurns = {110};
        std::vector<int64_t> numberParallels = {1};
        std::string coreShapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 0.15 - Grade 1");
        wires.push_back(wire);
        
        auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

        int64_t numberStacks = 1;
        std::string coreMaterialName = "A07";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        StrayCapacitance strayCapacitance;

        std::map<std::string, std::map<std::string, double>> expectedValues = {
            {"winding 0", {{"winding 0", 3.15617e-12 -3.14887e-12}}},
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
            outFile.append("Test_One_Versus_Two_Windings_One.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
    {
        std::vector<int64_t> numberTurns = {110, 110};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string coreShapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 0.15 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

        int64_t numberStacks = 1;
        std::string coreMaterialName = "A07";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        StrayCapacitance strayCapacitance;

        std::map<std::string, std::map<std::string, double>> expectedValues = {
            {"winding 0", {{"winding 0", 3.15617e-12}, {"winding 1", -3.14887e-12}}},
            {"winding 1", {{"winding 0", -3.14887e-12}, {"winding 1", 3.15674e-12}}},
        };

        auto maxwellCapacitanceMatrix = strayCapacitance.calculate_capacitance(coil).get_maxwell_capacitance_matrix().value();
        CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude().size() == 2);
        CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["winding 0"].size() == 2);
        CHECK(maxwellCapacitanceMatrix[0].get_mutable_magnitude()["winding 1"].size() == 2);
        for (auto [firstKey, aux] : maxwellCapacitanceMatrix[0].get_magnitude()) {
            for (auto [secondKey, capacitanceWithTolerance] : aux) {
                auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
                CHECK_THAT(capacitance, WithinRel(expectedValues[firstKey][secondKey], maximumError));
            }
        }

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_One_Versus_Two_Windings_Two.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
}

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

    // TEST(Test_Get_Capacitance_Large) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"_interleavingLevel": 1, "_windingOrientation": "overlapping", "_layersOrientation": "overlapping", "_turnsAlignment": "spread", "_sectionAlignment": "inner or top", "bobbin": {"processedDescription": {"columnDepth": 0.0062, "columnShape": "rectangular", "columnThickness": 0.0011999999999999997, "columnWidth": 0.0062, "coordinates": [0, 0, 0 ], "wallThickness": 0.0011999999999999997, "windingWindows": [{"coordinates": [0.00935, 0, 0 ], "height": 0.022600000000000002, "width": 0.006300000000000001 } ] } }, "functionalDescription": [{"connections": [{"pinName": "4", "type": "Pin" }, {"pinName": "5", "type": "Pin" }, {"pinName": "6", "type": "Pin" }, {"pinName": "7", "type": "Pin" }, {"pinName": "8", "type": "Pin" }, {"pinName": "9", "type": "Pin" }, {"pinName": "10", "type": "Pin" }, {"pinName": "11", "type": "Pin" } ], "isolationSide": "primary", "name": "primary", "numberParallels": 16, "numberTurns": 4, "wire": "Round 25.0 - Single Build" }, {"connections": [{"pinName": "13", "type": "Pin" }, {"pinName": "3", "type": "Pin" } ], "isolationSide": "secondary", "name": "secondary", "numberParallels": 2, "numberTurns": 40, "wire": {"conductingDiameter": {"maximum": 0.000457, "minimum": 0.00045000000000000004, "nominal": 0.000455 }, "material": "copper", "outerDiameter": {"maximum": 0.000516, "minimum": 0.000495, "nominal": 0.000505 }, "coating": {"breakdownVoltage": 2370, "grade": 2, "type": "enamelled" }, "manufacturerInfo": {"name": "Elektrisola" }, "name": "Round 25.0 - Heavy Build", "numberConductors": 1, "standard": "NEMA MW 1000 C", "standardName": "25.0", "type": "round" } } ] })");
    //     json coreJson = json::parse(R"({"functionalDescription": {"gapping": [{"area": 0.0001, "coordinates": [0, 0.00055, 0 ], "distanceClosestNormalSurface": 0.01195, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.0011, "sectionDimensions": [0.01, 0.01 ], "shape": "rectangular", "type": "subtractive" }, {"area": 0.000051, "coordinates": [0.015001, 0, 0 ], "distanceClosestNormalSurface": 0.0125, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.000005, "sectionDimensions": [0.005001, 0.01 ], "shape": "rectangular", "type": "residual" }, {"area": 0.000051, "coordinates": [-0.015001, 0, 0 ], "distanceClosestNormalSurface": 0.0125, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.000005, "sectionDimensions": [0.005001, 0.01 ], "shape": "rectangular", "type": "residual" } ], "material": "TP4A", "numberStacks": 1, "shape": {"aliases": ["E 35" ], "dimensions": {"A": {"maximum": 0.0355, "minimum": 0.0345 }, "B": {"maximum": 0.01775, "minimum": 0.01725 }, "C": {"maximum": 0.0103, "minimum": 0.0097 }, "D": {"maximum": 0.01275, "minimum": 0.01225 }, "E": {"maximum": 0.0255, "minimum": 0.0245 }, "F": {"maximum": 0.0103, "minimum": 0.0097 } }, "family": "e", "magneticCircuit": "open", "name": "E 35/18/10", "type": "standard" }, "type": "two-piece set" }, "processedDescription": {"columns": [{"area": 0.0001, "coordinates": [0, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "central", "width": 0.01 }, {"area": 0.000051, "coordinates": [0.015001, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "lateral", "width": 0.005001 }, {"area": 0.000051, "coordinates": [-0.015001, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "lateral", "width": 0.005001 } ], "depth": 0.01, "effectiveParameters": {"effectiveArea": 0.00010000000000000002, "effectiveLength": 0.08070796326794898, "effectiveVolume": 0.000008070796326794898, "minimumArea": 0.0001 }, "height": 0.035, "width": 0.035, "windingWindows": [{"area": 0.00018750000000000003, "coordinates": [0.005, 0 ], "height": 0.025, "width": 0.007500000000000001 } ] } })");

    //     auto core = OpenMagnetics::Core(coreJson);
    //     auto coil = OpenMagnetics::Coil(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
    //     coil.wind();

    //     StrayCapacitance strayCapacitance;

    //     std::map<std::string, std::map<std::string, double>> expectedValues = {
    //         {"Primary", {{"Primary", 3.12e-12}}},
    //         {"Primary", {"Secondary", 2.2e-12}},
    //         {"Secondary", {"Secondary", 3.3e-12}},
    //     };

    //     // auto capacitanceAmongTurns = strayCapacitance.calculate_capacitance_among_turns(coil);
    //     // CHECK(capacitanceAmongTurns.size() == 4);
    //     // for (auto [keys, capacitance] : capacitanceAmongTurns) {
    //     //     std::cout << "Capacitance between turn " << keys.first << " and turn " << keys.second << ": " << capacitance << std::endl;
    //     // }

    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    //     OpenMagnetics::Magnetic magnetic;
    //     magnetic.set_core(core);
    //     magnetic.set_coil(coil);


    //     auto outFile = outputFilePath;
    //     outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels.svg");
    //     std::filesystem::remove(outFile);
    //     Painter painter(outFile);
    //     painter.paint_core(magnetic);
    //     painter.paint_bobbin(magnetic);
    //     painter.paint_coil_turns(magnetic);
    //     painter.export_svg();

    // }


    // TEST(Test_Get_Surrounding_Turns_All) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"_interleavingLevel": 1, "_windingOrientation": "overlapping", "_layersOrientation": "overlapping", "_turnsAlignment": "spread", "_sectionAlignment": "inner or top", "bobbin": {"processedDescription": {"columnDepth": 0.0062, "columnShape": "rectangular", "columnThickness": 0.0011999999999999997, "columnWidth": 0.0062, "coordinates": [0, 0, 0 ], "wallThickness": 0.0011999999999999997, "windingWindows": [{"coordinates": [0.00935, 0, 0 ], "height": 0.022600000000000002, "width": 0.006300000000000001 } ] } }, "functionalDescription": [{"connections": [{"pinName": "4", "type": "Pin" }, {"pinName": "5", "type": "Pin" }, {"pinName": "6", "type": "Pin" }, {"pinName": "7", "type": "Pin" }, {"pinName": "8", "type": "Pin" }, {"pinName": "9", "type": "Pin" }, {"pinName": "10", "type": "Pin" }, {"pinName": "11", "type": "Pin" } ], "isolationSide": "primary", "name": "primary", "numberParallels": 2, "numberTurns": 16, "wire": "Round 19.0 - Heavy Build" }, {"connections": [{"pinName": "13", "type": "Pin" }, {"pinName": "3", "type": "Pin" } ], "isolationSide": "secondary", "name": "secondary", "numberParallels": 2, "numberTurns": 40, "wire": {"conductingDiameter": {"maximum": 0.000457, "minimum": 0.00045000000000000004, "nominal": 0.000455 }, "material": "copper", "outerDiameter": {"maximum": 0.000516, "minimum": 0.000495, "nominal": 0.000505 }, "coating": {"breakdownVoltage": 2370, "grade": 2, "type": "enamelled" }, "manufacturerInfo": {"name": "Elektrisola" }, "name": "Round 25.0 - Heavy Build", "numberConductors": 1, "standard": "NEMA MW 1000 C", "standardName": "25.0", "type": "round" } } ] })");
    //     json coreJson = json::parse(R"({"functionalDescription": {"gapping": [{"area": 0.0001, "coordinates": [0, 0.00055, 0 ], "distanceClosestNormalSurface": 0.01195, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.0011, "sectionDimensions": [0.01, 0.01 ], "shape": "rectangular", "type": "subtractive" }, {"area": 0.000051, "coordinates": [0.015001, 0, 0 ], "distanceClosestNormalSurface": 0.0125, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.000005, "sectionDimensions": [0.005001, 0.01 ], "shape": "rectangular", "type": "residual" }, {"area": 0.000051, "coordinates": [-0.015001, 0, 0 ], "distanceClosestNormalSurface": 0.0125, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.000005, "sectionDimensions": [0.005001, 0.01 ], "shape": "rectangular", "type": "residual" } ], "material": "TP4A", "numberStacks": 1, "shape": {"aliases": ["E 35" ], "dimensions": {"A": {"maximum": 0.0355, "minimum": 0.0345 }, "B": {"maximum": 0.01775, "minimum": 0.01725 }, "C": {"maximum": 0.0103, "minimum": 0.0097 }, "D": {"maximum": 0.01275, "minimum": 0.01225 }, "E": {"maximum": 0.0255, "minimum": 0.0245 }, "F": {"maximum": 0.0103, "minimum": 0.0097 } }, "family": "e", "magneticCircuit": "open", "name": "E 35/18/10", "type": "standard" }, "type": "two-piece set" }, "processedDescription": {"columns": [{"area": 0.0001, "coordinates": [0, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "central", "width": 0.01 }, {"area": 0.000051, "coordinates": [0.015001, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "lateral", "width": 0.005001 }, {"area": 0.000051, "coordinates": [-0.015001, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "lateral", "width": 0.005001 } ], "depth": 0.01, "effectiveParameters": {"effectiveArea": 0.00010000000000000002, "effectiveLength": 0.08070796326794898, "effectiveVolume": 0.000008070796326794898, "minimumArea": 0.0001 }, "height": 0.035, "width": 0.035, "windingWindows": [{"area": 0.00018750000000000003, "coordinates": [0.005, 0 ], "height": 0.025, "width": 0.007500000000000001 } ] } })");

    //     auto core = OpenMagnetics::Core(coreJson);
    //     auto coil = OpenMagnetics::Coil(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
    //     coil.wind();

    //     auto turns = coil.get_turns_description().value();

    //     OpenMagnetics::Magnetic magnetic;
    //     magnetic.set_core(core);
    //     magnetic.set_coil(coil);

    //     size_t index = 0;
    //     for (auto turn : turns) {

    //         auto surroundingTurns = StrayCapacitance::get_surrounding_turns(turn, turns);
    //         std::vector<Turn> newTurns;
    //         for (auto [surroundingTurn, surroundingTurnIndex] : surroundingTurns) {
    //             newTurns.push_back(surroundingTurn);
    //         }
    //         newTurns.push_back(turn);

    //         // magnetic.get_mutable_coil().set_turns_description(newTurns);

    //         auto outFile = outputFilePath;
    //         outFile.append("Test_Get_Surrounding_Turns_" + std::to_string(index) + ".svg");
    //         std::filesystem::remove(outFile);
    //         Painter painter(outFile);
    //         painter.paint_core(magnetic);
    //         painter.paint_bobbin(magnetic);
    //         painter.paint_coil_turns(magnetic);
    //         painter.export_svg();
    //         index++;
    //         break;
    //     }

    // }



TEST_CASE("Calculate capacitance of an automatic buck produced in OM with three layers", "[physical-model][stray-capacitance]") {
    settings->reset();

    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/buck_inductor three layers.json");
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

TEST_CASE("Calculate capacitance of an automatic buck produced in OM with two layers", "[physical-model][stray-capacitance]") {
    settings->reset();

    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/buck_inductor two layers.json");
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

TEST_CASE("Calculate capacitance of an automatic buck produced in OM with one layer", "[physical-model][stray-capacitance]") {
    settings->reset();

    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/buck_inductor one layer.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();

    StrayCapacitance strayCapacitance;

    std::map<std::string, std::map<std::string, double>> expectedValues = {
        {"Primary", {{"Primary", 0.3e-12}}},
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

TEST_CASE("Calculate area between two round turns", "[physical-model][stray-capacitance]") {
    {
        std::vector<int64_t> numberTurns = {1, 1};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string coreShapeName = "E 35";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 2.00 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);
        
        auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

        int64_t numberStacks = 1;
        std::string coreMaterialName = "A07";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        StrayCapacitance strayCapacitance;
        std::cout << "coil.get_turns_description()->size(): " << coil.get_turns_description()->size() << std::endl;
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[1];
        auto area = strayCapacitance.calculate_area_between_two_turns(firstTurn, secondTurn);
        std::cout << "area: " << area << std::endl;
        CHECK(area < 1e-6);
        CHECK(area > 0.86e-6);
    }
}

TEST_CASE("Calculate area between two rectangular turns", "[physical-model][stray-capacitance]") {
    {
        std::vector<int64_t> numberTurns = {1, 1};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string coreShapeName = "E 35";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Rectangular 2x0.80 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);
        
        auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

        int64_t numberStacks = 1;
        std::string coreMaterialName = "A07";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        StrayCapacitance strayCapacitance;
        std::cout << "coil.get_turns_description()->size(): " << coil.get_turns_description()->size() << std::endl;
        auto layerThickness = coil.get_layers_description().value()[1].get_dimensions()[0];
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[1];
        auto area = strayCapacitance.calculate_area_between_two_turns(firstTurn, secondTurn);

        CHECK_THAT(area, WithinRel(layerThickness * 0.0008, maximumError));
    }
}

TEST_CASE("Calculate area between round and rectangular turns", "[physical-model][stray-capacitance]") {
    {
        std::vector<int64_t> numberTurns = {1, 1};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string coreShapeName = "E 35";
        std::vector<OpenMagnetics::Wire> wires;
        auto firstWire = find_wire_by_name("Round 2.00 - Grade 1");
        wires.push_back(firstWire);
        auto secondWire = find_wire_by_name("Rectangular 2x0.80 - Grade 1");
        wires.push_back(secondWire);
        
        auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

        int64_t numberStacks = 1;
        std::string coreMaterialName = "A07";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        StrayCapacitance strayCapacitance;
        std::cout << "coil.get_turns_description()->size(): " << coil.get_turns_description()->size() << std::endl;
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[1];
        auto area = strayCapacitance.calculate_area_between_two_turns(firstTurn, secondTurn);
        std::cout << "area: " << area << std::endl;
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


TEST_CASE("Calculate area between round and rectangular turns diagonally", "[physical-model][stray-capacitance]") {
    {
        std::vector<int64_t> numberTurns = {2, 1};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string coreShapeName = "E 35";
        std::vector<OpenMagnetics::Wire> wires;
        auto firstWire = find_wire_by_name("Round 2.00 - Grade 1");
        wires.push_back(firstWire);
        auto secondWire = find_wire_by_name("Rectangular 2x0.80 - Grade 1");
        wires.push_back(secondWire);
        
        auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

        int64_t numberStacks = 1;
        std::string coreMaterialName = "A07";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        StrayCapacitance strayCapacitance;
        std::cout << "coil.get_turns_description()->size(): " << coil.get_turns_description()->size() << std::endl;
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[2];
        auto area = strayCapacitance.calculate_area_between_two_turns(firstTurn, secondTurn);
        std::cout << "area: " << area << std::endl;
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

TEST_CASE("Calculate energy density between two turns", "[physical-model][stray-capacitance]") {
    {
        std::vector<int64_t> numberTurns = {1, 1};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string coreShapeName = "E 35";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 2.00 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);
        
        auto coil = OpenMagnetics::Coil::create_quick_coil(coreShapeName, numberTurns, numberParallels, wires);

        int64_t numberStacks = 1;
        std::string coreMaterialName = "A07";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagnetics::Core::create_quick_core(coreShapeName, coreMaterialName, gapping, numberStacks);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        StrayCapacitance strayCapacitance;
        std::cout << "coil.get_turns_description()->size(): " << coil.get_turns_description()->size() << std::endl;
        auto firstTurn = coil.get_turns_description().value()[0];
        auto secondTurn = coil.get_turns_description().value()[1];
        auto energy = strayCapacitance.calculate_energy_density_between_two_turns(firstTurn, wire, secondTurn, wire, 10);
        std::cout << "energy: " << energy << std::endl;

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