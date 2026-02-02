#include <source_location>
#include <iomanip>
#include <map>
#include "physical_models/LeakageInductance.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "processors/MagneticSimulator.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

using namespace MAS; 
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

static double maximumError = 0.3;
static auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
static bool plot = false;

TEST_CASE("Calculate leakage inductance for a E core with same number of turns", "[physical-model][leakage-inductance]") {
    settings.reset();
    std::vector<int64_t> numberTurns({69, 69});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "E 42/33/20";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 25));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 25));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 6.7e-6;
    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a E core with different number of turns", "[physical-model][leakage-inductance]") {
    std::vector<int64_t> numberTurns({64, 20});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "E 42/33/20";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 25));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 225));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 13e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a E core with different number of turns and larger Litz wire", "[physical-model][leakage-inductance][smoke-test]") {

    std::vector<int64_t> numberTurns({16, 6});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "E 42/33/20";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 370));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 666));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 4e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_E_2.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a larger E core with different number of turns and larger Litz wire and interleaving", "[physical-model][leakage-inductance][smoke-test]") {

    std::vector<int64_t> numberTurns({36, 26});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "E 65/32/27";
    uint8_t interleavingLevel = 2;
    auto windingOrientation = WindingOrientation::OVERLAPPING;
    auto layersOrientation = WindingOrientation::OVERLAPPING;
    auto turnsAlignment = CoilAlignment::CENTERED;
    auto sectionsAlignment = CoilAlignment::CENTERED;

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 650));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 700));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment, interleavingLevel);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 9e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with with several parallels", "[physical-model][leakage-inductance][smoke-test]") {

    std::vector<int64_t> numberTurns({24, 6});
    std::vector<int64_t> numberParallels({2, 4});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "PQ 32/25";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 75));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 225));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 9e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_Parallels_No_Interleaving.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with with several parallels with interleaving", "[physical-model][leakage-inductance][smoke-test]") {

    std::vector<int64_t> numberTurns({24, 6});
    std::vector<int64_t> numberParallels({2, 4});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "PQ 32/25";
    uint8_t interleavingLevel = 2;
    auto windingOrientation = WindingOrientation::OVERLAPPING;
    auto layersOrientation = WindingOrientation::OVERLAPPING;
    auto turnsAlignment = CoilAlignment::CENTERED;
    auto sectionsAlignment = CoilAlignment::CENTERED;

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 75));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 225));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment, interleavingLevel);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 5e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
    auto leakageMagneticField = LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency);
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_Parallels_Interleaving.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(OperatingPoint(), magnetic, 1, leakageMagneticField);
        painter.paint_core(magnetic);
        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a ETD core", "[physical-model][leakage-inductance]") {

    std::vector<int64_t> numberTurns({60, 59});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "ETD 39";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.000071, 135));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.000071, 69));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 40e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with contiguous winding orientation", "[physical-model][leakage-inductance][smoke-test]") {
    settings.set_coil_try_rewind(false);
    std::vector<int64_t> numberTurns({27, 3});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "PQ 26/25";
    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    auto layersOrientation = WindingOrientation::OVERLAPPING;
    auto turnsAlignment = CoilAlignment::INNER_OR_TOP;
    auto sectionsAlignment = CoilAlignment::SPREAD;

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.000073, 60));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.000055, 370));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;

    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 86e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Leakage_Inductance_Failing_Test.json");
        OpenMagnetics::to_file(outFile, magnetic);
    }
    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_PQ_26_0.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with overlapping winding orientation and large Litz wires", "[physical-model][leakage-inductance][smoke-test]") {
    settings.set_leakage_inductance_grid_precision_level_wound(2);
    std::vector<int64_t> numberTurns({20, 2});
    std::vector<int64_t> numberParallels({1, 3});
    std::vector<double> turnsRatios({10});
    std::string shapeName = "PQ 40/40";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 800));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 1000));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 9.9e-6;
    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with contiguous winding orientation and large Litz wires", "[physical-model][leakage-inductance][smoke-test]") {
    double localMaximumError = 0.4;
    settings.set_coil_try_rewind(false);
    std::vector<int64_t> numberTurns({20, 2});
    std::vector<int64_t> numberParallels({1, 3});
    std::vector<double> turnsRatios({10});
    std::string shapeName = "PQ 40/40";
    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    auto layersOrientation = WindingOrientation::OVERLAPPING;
    auto turnsAlignment = CoilAlignment::INNER_OR_TOP;
    auto sectionsAlignment = CoilAlignment::SPREAD;

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 800));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 1000));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 52e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();

    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, localMaximumError));
    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Leakage_Inductance_Failing_Test.json");
        OpenMagnetics::to_file(outFile, magnetic);
    }
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for an E core with three windings", "[physical-model][leakage-inductance]") {
    // OpenMagnetics::Logger::getInstance().setLevel(OpenMagnetics::LogLevel::DEBUG);
    settings.reset();
    // settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns({50, 100, 25});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1], double(numberTurns[0]) / numberTurns[2]});
    std::string shapeName = "E 42/21/15";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 200));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 200));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 200));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    double frequency = 100000;

    auto leakageInductance_01 = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    auto leakageInductance_10 = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 1, 0).get_leakage_inductance_per_winding()[0].get_nominal().value();

    CHECK_THAT(leakageInductance_01, WithinRel(leakageInductance_10 * pow(double(numberTurns[0]) / numberTurns[1], 2), 0.01));

    auto leakageInductance_02 = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 2).get_leakage_inductance_per_winding()[0].get_nominal().value();
    auto leakageInductance_20 = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 2, 0).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance_02, WithinRel(leakageInductance_20 * pow(double(numberTurns[0]) / numberTurns[2], 2), 0.01));

    auto leakageInductance_12 = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 1, 2).get_leakage_inductance_per_winding()[0].get_nominal().value();
    auto leakageInductance_21 = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 2, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance_12, WithinRel(leakageInductance_21 * pow(double(numberTurns[1]) / numberTurns[2], 2), 0.01));
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_T_1.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Calculate leakage inductance for toroidal cores with contiguous sections", "[physical-model][leakage-inductance][toroidal-core][smoke-test]") {
    settings.reset();
    clear_databases();
    settings.set_coil_try_rewind(false);
    std::vector<int64_t> numberTurns({10, 200});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "T 48/28/16";
    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    auto layersOrientation = WindingOrientation::OVERLAPPING;
    auto turnsAlignment = CoilAlignment::CENTERED;
    auto sectionsAlignment = CoilAlignment::SPREAD;

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::find_wire_by_name("Round 1.00 - Grade 1"));
    wires.push_back(OpenMagnetics::find_wire_by_name("Round 0.80 - Grade 1"));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);

    std::string coreMaterial = "3C97";
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial);

    std::vector<double> proportionPerWinding = {16.185 / (222.42 + 16.185), 222.42 / (222.42 + 16.185)};
    std::vector<size_t> pattern = {0, 1};
    coil.wind(proportionPerWinding, pattern);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 0.00143;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 1, 0).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    if (plot) {
        settings.set_painter_mode(PainterModes::QUIVER);
        auto leakageMagneticField = LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);

        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_T_0.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(OperatingPoint(), magnetic, 1, leakageMagneticField);
        painter.paint_core(magnetic);
        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }


    settings.reset();
}

TEST_CASE("Calculate leakage inductance for toroidal cores with contiguous sections with few turns", "[physical-model][leakage-inductance][toroidal-core][smoke-test]") {
    settings.reset();
    settings.set_coil_try_rewind(false);
    clear_databases();
    std::vector<int64_t> numberTurns({10, 5});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "T 48/28/16";
    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    auto layersOrientation = WindingOrientation::OVERLAPPING;
    auto turnsAlignment = CoilAlignment::CENTERED;
    auto sectionsAlignment = CoilAlignment::SPREAD;

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_rectangular_wire(0.0038, 0.002));
    wires.push_back(OpenMagnetics::Wire::create_quick_rectangular_wire(0.0038, 0.001));
    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);

    std::string coreMaterial = "3C97";
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    double expectedLeakageInductance = 7e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_T_1.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }


    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a complex planar magnetic", "[physical-model][leakage-inductance][planar][smoke-test]") {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "leakage_inductance_planar.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path.string(), mas);
    auto magnetic = mas.get_magnetic();

    double frequency = 100000;
    double expectedLeakageInductance = 1.4e-6;
    // settings.set_magnetic_field_number_points_x(100);
    // settings.set_magnetic_field_number_points_y(100);
    std::vector<double> turnsRatios = magnetic.get_turns_ratios();

    auto operatingPoint = OpenMagnetics::Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, 0.001, 25, turnsRatios, {sqrt(2), sqrt(2), 0});
    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    // auto leakageMagneticField = LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        settings.set_painter_include_fringing(false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        outFile.append("Test_Leakage_Inductance_Planar.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }


    settings.reset();
}

TEST_CASE("Checks that increasing insulation between layers keeps leakage inductance consistent", "[physical-model][leakage-inductance][planar][bug][smoke-test]") {
    settings.set_leakage_inductance_grid_auto_scaling(true);
    settings.set_coil_maximum_layers_planar(60);
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "bug_unstable_leakage.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path.string(), mas);
    auto magnetic = mas.get_magnetic();

    double frequency = 100000;

    std::vector<double> turnsRatios = magnetic.get_turns_ratios();

    double previousLeakageInductance = 0;
    for (size_t index = 0; index < 10; ++index) {

        std::map<std::pair<size_t, size_t>, double> insulationThickness;
        insulationThickness[{0, 0}] = 20e-6 + index * 20e-6;
        insulationThickness[{0, 1}] = 20e-6 + index * 20e-6;
        insulationThickness[{1, 0}] = 20e-6 + index * 20e-6;
        insulationThickness[{1, 1}] = 20e-6 + index * 20e-6;
        magnetic.get_mutable_coil().wind_planar({1,1,0,0,0,1,1,0,0,1,1,0,0,0,1,1}, 100e-6, {}, insulationThickness, 200e-6);

        auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();

        CHECK(previousLeakageInductance < leakageInductance);
        previousLeakageInductance = leakageInductance;
    }

    if (plot) {
        auto operatingPoint = OpenMagnetics::Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, 0.001, 25, turnsRatios, {sqrt(2), sqrt(2), 0});
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        settings.set_painter_include_fringing(false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        outFile.append("Test_Leakage_Inductance_Planar_Bug_Insulation.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }


    settings.reset();
}

TEST_CASE("Benchmarks leakage inductance calculation in planar", "[physical-model][leakage-inductance][!benchmark]") {
    BENCHMARK_ADVANCED("measures computation time")(Catch::Benchmark::Chronometer meter) {

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "leakage_inductance_planar.json");
        OpenMagnetics::Mas mas;
        OpenMagnetics::from_file(path.string(), mas);
        auto magnetic = mas.get_magnetic();
        double frequency = 100000;

        meter.measure([&magnetic, &frequency] { return LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding(); });
    };
}

TEST_CASE("Calculate leakage inductance for a planar magnetic from the web", "[physical-model][leakage-inductance][planar][bug][smoke-test]") {
    auto json_path_631 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "calculate_leakage_inductance_for_a_planar_magnetic_from_the_web_631.json");
    std::ifstream json_file_631(json_path_631);
    OpenMagnetics::Magnetic magnetic(json::parse(json_file_631));

    double frequency = 150000;
    double expectedLeakageInductance = 1.4e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));

    settings.reset();
}

TEST_CASE("Calculate leakage inductance for a planar magnetic from the web 2", "[physical-model][leakage-inductance][planar][bug][smoke-test]") {
    auto json_path_643 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "calculate_leakage_inductance_for_a_planar_magnetic_from_the_web_2_643.json");
    std::ifstream json_file_643(json_path_643);
    OpenMagnetics::Magnetic magnetic(json::parse(json_file_643));
    auto json_path_644 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "calculate_leakage_inductance_for_a_planar_magnetic_from_the_web_2_644.json");
    std::ifstream json_file_644(json_path_644);
    OpenMagnetics::Inputs inputs(json::parse(json_file_644));

    double expectedLeakageInductance = 1.4e-6;

    OpenMagnetics::MagneticSimulator magneticSimulator;
    auto mas = magneticSimulator.simulate(inputs, magnetic);
    auto leakageInductance = resolve_dimensional_values(mas.get_outputs()[0].get_inductance()->get_leakage_inductance()->get_leakage_inductance_per_winding()[0]);
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Leakage_Inductance_Failing_Test_2.json");
        OpenMagnetics::to_file(outFile, magnetic);
    }

    settings.reset();
}


TEST_CASE("Calculate leakage inductance for a simple planar magnetic", "[physical-model][leakage-inductance][planar][smoke-test]") {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simplest_leakage_inductance_planar.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path.string(), mas);
    auto magnetic = mas.get_magnetic();

    double frequency = 100000;
    double expectedLeakageInductance = 2e-9;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));



    settings.reset();
}


TEST_CASE("Calculate leakage inductance for a planar magnetic from the web 3", "[physical-model][leakage-inductance][planar][smoke-test]") {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "OM Oyang paper example.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path.string(), mas);
    auto magnetic = mas.get_magnetic();

    double frequency = 100000;
    double expectedLeakageInductance = 1.5e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));


    settings.reset();
}

// ============================================================================
// MODEL COMPARISON STUDY - Compare all H-field models for leakage inductance
// ============================================================================

struct LeakageTestCase {
    std::string name;
    std::string shapeName;
    std::vector<int64_t> numberTurns;
    std::vector<int64_t> numberParallels;
    int strandDiameter;  // in microns
    std::vector<int> numberStrands;
    double frequency;
    double expectedLeakageInductance;
};

TEST_CASE("Leakage inductance H-field model comparison study", "[physical-model][leakage-inductance][model-comparison]") {
    settings.reset();
    
    std::vector<LeakageTestCase> testCases = {
        {"E42 16:6 Litz", "E 42/33/20", {16, 6}, {1, 1}, 50, {370, 666}, 100000, 4e-6},
        {"E42 69:69 Litz", "E 42/33/20", {69, 69}, {1, 1}, 50, {25, 25}, 100000, 6.7e-6},
        {"E42 64:20 Litz", "E 42/33/20", {64, 20}, {1, 1}, 50, {25, 225}, 100000, 13e-6},
        {"E65 12:6 Litz", "E 65/32/27", {12, 6}, {1, 1}, 50, {450, 450}, 100000, 9e-6},
        {"PQ40 10:10 Litz x2par", "PQ 40/40", {10, 10}, {2, 2}, 100, {150, 150}, 100000, 5e-6},
    };

    std::vector<std::pair<MagneticFieldStrengthModels, std::string>> models = {
        {MagneticFieldStrengthModels::BINNS_LAWRENSON, "BINNS_LAWRENSON"},
        {MagneticFieldStrengthModels::LAMMERANER, "LAMMERANER"},
        {MagneticFieldStrengthModels::ALBACH, "ALBACH"},
    };

    std::cout << "\n====================================================================================" << std::endl;
    std::cout << "                   LEAKAGE INDUCTANCE H-FIELD MODEL COMPARISON" << std::endl;
    std::cout << "====================================================================================" << std::endl;
    std::cout << std::setw(25) << "Test Case" << " | " << std::setw(10) << "Expected" << " | ";
    for (auto& [model, modelName] : models) {
        std::cout << std::setw(18) << modelName << " | ";
    }
    std::cout << std::endl;
    std::cout << std::string(25 + 3 + 10 + 3 + models.size() * 21, '-') << std::endl;

    // Track errors per model
    std::map<std::string, std::vector<double>> errorsPerModel;

    for (auto& tc : testCases) {
        std::vector<OpenMagnetics::Wire> wires;
        for (size_t i = 0; i < tc.numberTurns.size(); ++i) {
            wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(tc.strandDiameter * 1e-6, tc.numberStrands[i]));
        }
        auto coil = OpenMagnetics::Coil::create_quick_coil(tc.shapeName, tc.numberTurns, tc.numberParallels, wires);

        std::string coreMaterial = "3C97";
        auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
        auto core = OpenMagnetics::Core::create_quick_core(tc.shapeName, coreMaterial, gapping);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        std::cout << std::setw(25) << tc.name << " | " 
                  << std::setw(8) << std::fixed << std::setprecision(1) << (tc.expectedLeakageInductance * 1e6) << " µH | ";

        for (auto& [model, modelName] : models) {
            settings.set_magnetic_field_strength_model(model);
            
            try {
                auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, tc.frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
                double error = (leakageInductance - tc.expectedLeakageInductance) / tc.expectedLeakageInductance * 100;
                errorsPerModel[modelName].push_back(fabs(error));
                
                std::string errorStr = (error > 0 ? "+" : "") + std::to_string(static_cast<int>(error)) + "%";
                std::cout << std::setw(8) << std::fixed << std::setprecision(1) << (leakageInductance * 1e6) 
                          << " µH (" << std::setw(5) << errorStr << ") | ";
            } catch (const std::exception& e) {
                std::cout << std::setw(18) << "ERROR" << " | ";
            }
        }
        std::cout << std::endl;
    }

    std::cout << std::string(25 + 3 + 10 + 3 + models.size() * 21, '-') << std::endl;
    std::cout << std::setw(25) << "Average |Error|" << " | " << std::setw(10) << "" << " | ";
    for (auto& [model, modelName] : models) {
        if (errorsPerModel[modelName].size() > 0) {
            double avg = 0;
            for (auto e : errorsPerModel[modelName]) avg += e;
            avg /= errorsPerModel[modelName].size();
            std::cout << std::setw(18) << std::fixed << std::setprecision(1) << avg << "% | ";
        }
    }
    std::cout << std::endl;
    std::cout << "====================================================================================\n" << std::endl;
    
    settings.reset();
}