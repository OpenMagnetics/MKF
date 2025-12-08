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
TEST_CASE("Calculate capacitance among two windings one with 16 and another with 8 turns and both 1 parallel", "[physical-model][stray-capacitance]") {
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

// TEST_CASE("Calculate and compare capacitance matrix of one winding versus two windings", "[physical-model][stray-capacitance]") {
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

TEST_CASE("Calculate capacitance of a tranformers with low filling factor", "[physical-model][stray-capacitance]") {
    settings->reset();

    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/low_filling_transformer.json");
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

TEST_CASE("Calculate capacitance of a one layer inductor", "[physical-model][stray-capacitance]") {
    settings->reset();

    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/low_filling_inductor.json");
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

TEST_CASE("Calculate capacitance of a simple planar tranformer", "[physical-model][stray-capacitance][planar]") {
    settings->reset();

    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/simple_planar.json");
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

TEST_CASE("Calculate capacitance of a simple planar tranformer with a turns ratio of 5", "[physical-model][stray-capacitance][planar]") {
    settings->reset();

    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/simple_planar_5_to_1.json");
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

TEST_CASE("Calculate capacitance of a simple planar tranformer with imperfect overlapping", "[physical-model][stray-capacitance][planar]") {
    settings->reset();

    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/simple_planar_imperfect_overlapping.json");
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
