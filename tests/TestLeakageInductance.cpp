#include "physical_models/LeakageInductance.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "processors/MagneticSimulator.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

using namespace MAS; 
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

static double maximumError = 0.2;
static auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
static bool plot = false;

TEST_CASE("Calculate leakage inductance for a E core with same number of turns", "[physical-model][leakage-inductance]") {
    settings->reset();
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
    settings->reset();
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
    settings->reset();
}

TEST_CASE("Calculate leakage inductance for a E core with different number of turns and larger Litz wire", "[physical-model][leakage-inductance]") {

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
        auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_E_2.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    settings->reset();
}

TEST_CASE("Calculate leakage inductance for a larger E core with different number of turns and larger Litz wire and interleaving", "[physical-model][leakage-inductance]") {

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
    settings->reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with with several parallels", "[physical-model][leakage-inductance]") {

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

    settings->reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with with several parallels with interleaving", "[physical-model][leakage-inductance]") {

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
        auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

    settings->reset();
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
    settings->reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with contiguous winding orientation", "[physical-model][leakage-inductance]") {
    settings->set_coil_try_rewind(false);
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
        auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Leakage_Inductance_Failing_Test.json");
        OpenMagnetics::to_file(outFile, magnetic);
    }
    if (plot) {
        // settings->set_painter_mode(PainterModes::QUIVER);

        // auto leakageMagneticField = LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);

        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_PQ_26_0.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false);

        // painter.paint_magnetic_field(OperatingPoint(), magnetic, 1, leakageMagneticField);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    settings->reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with overlapping winding orientation and large Litz wires", "[physical-model][leakage-inductance]") {
    settings->set_leakage_inductance_grid_precision_level_wound(2);
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
    std::cout << "leakageInductance: " << leakageInductance << std::endl;
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    // if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_Large_Litzs.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    // }
    settings->reset();
}

TEST_CASE("Calculate leakage inductance for a PQ core with contiguous winding orientation and large Litz wires", "[physical-model][leakage-inductance]") {
    double localMaximumError = 0.4;
    settings->set_coil_try_rewind(false);
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
        auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Leakage_Inductance_Failing_Test.json");
        OpenMagnetics::to_file(outFile, magnetic);
    }
    settings->reset();
}

TEST_CASE("Calculate leakage inductance for an E core with three windings", "[physical-model][leakage-inductance]") {
    settings->reset();
    std::vector<int64_t> numberTurns({50, 100, 25});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1], double(numberTurns[0]) / numberTurns[2]});
    std::string shapeName = "E 42/21/15";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 350));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 350));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 350));
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
    settings->reset();
}

TEST_CASE("Calculate leakage inductance for toroidal cores with contiguous sections", "[physical-model][leakage-inductance][toroidal-core]") {
    settings->reset();
    clear_databases();
    settings->set_coil_try_rewind(false);
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
        settings->set_painter_mode(PainterModes::QUIVER);
        auto leakageMagneticField = LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);

        auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_T_0.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(OperatingPoint(), magnetic, 1, leakageMagneticField);
        painter.paint_core(magnetic);
        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }


    settings->reset();
}

TEST_CASE("Calculate leakage inductance for toroidal cores with contiguous sections with few turns", "[physical-model][leakage-inductance][toroidal-core]") {
    settings->reset();
    settings->set_coil_try_rewind(false);
    clear_databases();
    std::vector<int64_t> numberTurns({10, 5});
    std::vector<int64_t> numberParallels({1, 1});
    std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
    std::string shapeName = "T 48/28/16";
    uint8_t interleavingLevel = 1;
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
    double expectedLeakageInductance = 5e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    if (true) {
        // auto operatingPoint = OpenMagnetics::Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, 0.001, 25, turnsRatios, {sqrt(2), sqrt(2)});
        auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Leakage_Inductance_T_1.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        // painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        // painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }


    settings->reset();
}

TEST_CASE("Calculate leakage inductance for a complex planar magnetic", "[physical-model][leakage-inductance][planar]") {
    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/leakage_inductance_planar.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
    auto magnetic = mas.get_magnetic();

    double frequency = 100000;
    double expectedLeakageInductance = 1.4e-6;
    // settings->set_magnetic_field_number_points_x(100);
    // settings->set_magnetic_field_number_points_y(100);
    std::vector<double> turnsRatios = magnetic.get_turns_ratios();

    auto operatingPoint = OpenMagnetics::Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, 0.001, 25, turnsRatios, {sqrt(2), sqrt(2), 0});
    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    // auto leakageMagneticField = LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);
    std::cout << "leakageInductance: " << leakageInductance << std::endl;
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));
    // if (true) {
        auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        settings->set_painter_include_fringing(false);
        settings->set_painter_mode(PainterModes::CONTOUR);
        outFile.append("Test_Leakage_Inductance_Planar.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    // }


    settings->reset();
}

TEST_CASE("Checks that increasing insulation between layers keeps leakage inductance consistent", "[physical-model][leakage-inductance][planar][bug]") {
    settings->set_leakage_inductance_grid_auto_scaling(true);
    settings->set_coil_maximum_layers_planar(60);
    std::string file_path = __FILE__;
    auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/bug_unstable_leakage.json");
    OpenMagnetics::Mas mas;
    OpenMagnetics::from_file(path, mas);
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

        std::cout << "previousLeakageInductance: " << previousLeakageInductance << std::endl;
        std::cout << "leakageInductance: " << leakageInductance << std::endl;
        CHECK(previousLeakageInductance < leakageInductance);
        previousLeakageInductance = leakageInductance;
    }

    if (plot) {
        auto operatingPoint = OpenMagnetics::Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, 0.001, 25, turnsRatios, {sqrt(2), sqrt(2), 0});
        auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        settings->set_painter_include_fringing(false);
        settings->set_painter_mode(PainterModes::CONTOUR);
        outFile.append("Test_Leakage_Inductance_Planar_Bug_Insulation.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }


    settings->reset();
}

TEST_CASE("Benchmarks leakage inductance calculation in planar", "[physical-model][leakage-inductance][!benchmark]") {
    BENCHMARK_ADVANCED("measures computation time")(Catch::Benchmark::Chronometer meter) {

        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/leakage_inductance_planar.json");
        OpenMagnetics::Mas mas;
        OpenMagnetics::from_file(path, mas);
        auto magnetic = mas.get_magnetic();
        double frequency = 100000;

        meter.measure([&magnetic, &frequency] { return LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding(); });
    };
}

TEST_CASE("Calculate leakage inductance for a planar magnetic from the web", "[physical-model][leakage-inductance][planar][bug]") {
    OpenMagnetics::Magnetic magnetic(json::parse(R"({"coil":{"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.0026,"columnShape":"round","columnThickness":0,"columnWidth":0.0026,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":null,"area":0.000011684999999999999,"coordinates":[0.0041375,0,0],"height":0.0038,"radialHeight":null,"sectionsAlignment":null,"sectionsOrientation":"contiguous","shape":"rectangular","width":0.0030749999999999996}]}},"functionalDescription":[{"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":21,"wire":{"coating":null,"conductingArea":null,"conductingDiameter":null,"conductingHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.00013919},"conductingWidth":{"nominal":0.0007999999999999999},"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Planar 139.19 µm","numberConductors":1,"outerDiameter":null,"outerHeight":{"nominal":0.00013919},"outerWidth":{"nominal":0.0007999999999999999},"standard":"IPC-6012","standardName":"4 oz.","strand":null,"type":"planar"}},{"name":"Secondary","numberTurns":15,"numberParallels":1,"isolationSide":"secondary","wire":{"conductingHeight":{"nominal":0.00024359},"conductingWidth":{"nominal":0.00045},"material":"copper","name":"Planar 243.59 µm","numberConductors":1,"outerHeight":{"nominal":0.00024359},"outerWidth":{"nominal":0.00045},"standard":"IPC-6012","standardName":"7 oz.","type":"planar"}},{"name":"Tertiary","numberTurns":7,"numberParallels":1,"isolationSide":"tertiary","wire":{"conductingHeight":{"nominal":0.0000696},"conductingWidth":{"nominal":0.0007999999999999999},"material":"copper","name":"Planar 69.60 µm","numberConductors":1,"outerHeight":{"nominal":0.0000696},"outerWidth":{"nominal":0.0007999999999999999},"standard":"IPC-6012","standardName":"2 oz.","type":"planar"}}],"layersDescription":[{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0014873549999999995],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"insulationMaterial":null,"name":"Primary layer 0","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0013677599999999995],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index0 and 1","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index0 and 1","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0011959649999999996],"dimensions":[0.0028749999999999995,0.00024359],"fillingFactor":0.9531782608695645,"insulationMaterial":null,"name":"Secondary layer 0","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.3333333333333333],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0010241699999999996],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index1 and 2","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index1 and 2","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0009045749999999996],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"insulationMaterial":null,"name":"Primary layer 1","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"section":"Primary section 1","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0007849799999999996],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index2 and 3","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index2 and 3","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0006131849999999996],"dimensions":[0.0028749999999999995,0.00024359],"fillingFactor":0.9531782608695645,"insulationMaterial":null,"name":"Secondary layer 1","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.3333333333333333],"winding":"Secondary"}],"section":"Secondary section 1","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0004413899999999996],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index3 and 4","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index3 and 4","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0003217949999999996],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"insulationMaterial":null,"name":"Primary layer 2","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"section":"Primary section 2","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0002021999999999996],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index4 and 5","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index4 and 5","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,0.000030404999999999625],"dimensions":[0.0028749999999999995,0.00024359],"fillingFactor":0.9531782608695645,"insulationMaterial":null,"name":"Secondary layer 2","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.3333333333333333],"winding":"Secondary"}],"section":"Secondary section 2","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00014139000000000035],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index5 and 6","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index5 and 6","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00026098500000000036],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"insulationMaterial":null,"name":"Primary layer 3","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"section":"Primary section 3","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00038058000000000037],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index6 and 7","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index6 and 7","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00046538000000000037],"dimensions":[0.0028749999999999995,0.0000696],"fillingFactor":0.2905043478260868,"insulationMaterial":null,"name":"Tertiary layer 0","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.42857142857142855],"winding":"Tertiary"}],"section":"Tertiary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0005501800000000004],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index7 and 8","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index7 and 8","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0006697750000000004],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"insulationMaterial":null,"name":"Primary layer 4","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"section":"Primary section 4","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0007893700000000004],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index8 and 9","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index8 and 9","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0008741700000000004],"dimensions":[0.0028749999999999995,0.0000696],"fillingFactor":0.19366956521739118,"insulationMaterial":null,"name":"Tertiary layer 1","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.2857142857142857],"winding":"Tertiary"}],"section":"Tertiary section 1","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0009589700000000004],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index9 and 10","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index9 and 10","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0010785650000000003],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"insulationMaterial":null,"name":"Primary layer 5","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"section":"Primary section 5","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0011981600000000002],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index10 and 11","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index10 and 11","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00128296],"dimensions":[0.0028749999999999995,0.0000696],"fillingFactor":0.19366956521739118,"insulationMaterial":null,"name":"Tertiary layer 2","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.2857142857142857],"winding":"Tertiary"}],"section":"Tertiary section 2","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00136776],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"insulationMaterial":"FR4","name":"Insulation layer between stack index11 and 12","orientation":"contiguous","partialWindings":[],"section":"Insulation section between stack index11 and 12","turnsAlignment":null,"type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0014873549999999999],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"insulationMaterial":null,"name":"Primary layer 6","orientation":"contiguous","partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"section":"Primary section 6","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"}],"sectionsDescription":[{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0014873549999999995],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Primary section 0","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0013677599999999995],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index0 and 1","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0011959649999999996],"dimensions":[0.0028749999999999995,0.00024359],"fillingFactor":0.9531782608695645,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Secondary section 0","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.3333333333333333],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0010241699999999996],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index1 and 2","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0009045749999999996],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Primary section 1","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0007849799999999996],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index2 and 3","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0006131849999999996],"dimensions":[0.0028749999999999995,0.00024359],"fillingFactor":0.9531782608695645,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Secondary section 1","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.3333333333333333],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0004413899999999996],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index3 and 4","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0003217949999999996],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Primary section 2","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.0002021999999999996],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index4 and 5","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,0.000030404999999999625],"dimensions":[0.0028749999999999995,0.00024359],"fillingFactor":0.9531782608695645,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Secondary section 2","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.3333333333333333],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00014139000000000035],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index5 and 6","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00026098500000000036],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Primary section 3","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00038058000000000037],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index6 and 7","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00046538000000000037],"dimensions":[0.0028749999999999995,0.0000696],"fillingFactor":0.2905043478260868,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Tertiary section 0","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.42857142857142855],"winding":"Tertiary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0005501800000000004],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index7 and 8","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0006697750000000004],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Primary section 4","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0007893700000000004],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index8 and 9","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0008741700000000004],"dimensions":[0.0028749999999999995,0.0000696],"fillingFactor":0.19366956521739118,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Tertiary section 1","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.2857142857142857],"winding":"Tertiary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0009589700000000004],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index9 and 10","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0010785650000000003],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Primary section 5","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0011981600000000002],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index10 and 11","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00128296],"dimensions":[0.0028749999999999995,0.0000696],"fillingFactor":0.19366956521739118,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Tertiary section 2","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.2857142857142857],"winding":"Tertiary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.00136776],"dimensions":[0.0028749999999999995,0.0001],"fillingFactor":1,"group":null,"layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Insulation section between stack index11 and 12","numberLayers":null,"partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.0041375,-0.0014873549999999999],"dimensions":[0.0028749999999999995,0.00013919],"fillingFactor":0.5809669565217387,"group":"Default group","layersAlignment":null,"layersOrientation":"contiguous","margin":[0,0],"name":"Primary section 6","numberLayers":null,"partialWindings":[{"connections":null,"parallelsProportion":[0.14285714285714285],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"}],"turnsDescription":[{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,0.0014873549999999995],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 0","length":0.020106192982974672,"name":"Primary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,0.0014873549999999995],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 0","length":0.02569822790636449,"name":"Primary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004979999999999995,0.0014873549999999995],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 0","length":0.031290262829754306,"name":"Primary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0030249999999999995,0.0011959649999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 0","length":0.019006635554218245,"name":"Secondary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0035649999999999983,0.0011959649999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 0","length":0.022399555620095213,"name":"Secondary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0041049999999999975,0.0011959649999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 0","length":0.025792475685972188,"name":"Secondary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004644999999999996,0.0011959649999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 0","length":0.029185395751849155,"name":"Secondary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005184999999999995,0.0011959649999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 0","length":0.032578315817726126,"name":"Secondary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,0.0009045749999999996],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 1","length":0.020106192982974672,"name":"Primary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 1","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,0.0009045749999999996],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 1","length":0.02569822790636449,"name":"Primary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 1","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004979999999999995,0.0009045749999999996],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 1","length":0.031290262829754306,"name":"Primary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 1","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0030249999999999995,0.0006131849999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 1","length":0.019006635554218245,"name":"Secondary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 1","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0035649999999999983,0.0006131849999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 1","length":0.022399555620095213,"name":"Secondary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 1","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0041049999999999975,0.0006131849999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 1","length":0.025792475685972188,"name":"Secondary parallel 0 turn 7","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 1","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004644999999999996,0.0006131849999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 1","length":0.029185395751849155,"name":"Secondary parallel 0 turn 8","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 1","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005184999999999995,0.0006131849999999996],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 1","length":0.032578315817726126,"name":"Secondary parallel 0 turn 9","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 1","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,0.0003217949999999996],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 2","length":0.020106192982974672,"name":"Primary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 2","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,0.0003217949999999996],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 2","length":0.02569822790636449,"name":"Primary parallel 0 turn 7","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 2","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004979999999999995,0.0003217949999999996],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 2","length":0.031290262829754306,"name":"Primary parallel 0 turn 8","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 2","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0030249999999999995,0.000030404999999999625],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 2","length":0.019006635554218245,"name":"Secondary parallel 0 turn 10","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 2","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0035649999999999983,0.000030404999999999625],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 2","length":0.022399555620095213,"name":"Secondary parallel 0 turn 11","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 2","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0041049999999999975,0.000030404999999999625],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 2","length":0.025792475685972188,"name":"Secondary parallel 0 turn 12","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 2","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004644999999999996,0.000030404999999999625],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 2","length":0.029185395751849155,"name":"Secondary parallel 0 turn 13","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 2","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005184999999999995,0.000030404999999999625],"crossSectionalShape":"rectangular","dimensions":[0.00045,0.00024359],"layer":"Secondary layer 2","length":0.032578315817726126,"name":"Secondary parallel 0 turn 14","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 2","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,-0.00026098500000000036],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 3","length":0.020106192982974672,"name":"Primary parallel 0 turn 9","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 3","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,-0.00026098500000000036],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 3","length":0.02569822790636449,"name":"Primary parallel 0 turn 10","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 3","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004979999999999995,-0.00026098500000000036],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 3","length":0.031290262829754306,"name":"Primary parallel 0 turn 11","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 3","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,-0.00046538000000000037],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.0000696],"layer":"Tertiary layer 0","length":0.020106192982974672,"name":"Tertiary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":0,"section":"Tertiary section 0","winding":"Tertiary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,-0.00046538000000000037],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.0000696],"layer":"Tertiary layer 0","length":0.02569822790636449,"name":"Tertiary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":0,"section":"Tertiary section 0","winding":"Tertiary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004979999999999995,-0.00046538000000000037],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.0000696],"layer":"Tertiary layer 0","length":0.031290262829754306,"name":"Tertiary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":0,"section":"Tertiary section 0","winding":"Tertiary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,-0.0006697750000000004],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 4","length":0.020106192982974672,"name":"Primary parallel 0 turn 12","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 4","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,-0.0006697750000000004],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 4","length":0.02569822790636449,"name":"Primary parallel 0 turn 13","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 4","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004979999999999995,-0.0006697750000000004],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 4","length":0.031290262829754306,"name":"Primary parallel 0 turn 14","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 4","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,-0.0008741700000000004],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.0000696],"layer":"Tertiary layer 1","length":0.020106192982974672,"name":"Tertiary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":0,"section":"Tertiary section 1","winding":"Tertiary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,-0.0008741700000000004],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.0000696],"layer":"Tertiary layer 1","length":0.02569822790636449,"name":"Tertiary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":0,"section":"Tertiary section 1","winding":"Tertiary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,-0.0010785650000000003],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 5","length":0.020106192982974672,"name":"Primary parallel 0 turn 15","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 5","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,-0.0010785650000000003],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 5","length":0.02569822790636449,"name":"Primary parallel 0 turn 16","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 5","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004979999999999995,-0.0010785650000000003],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 5","length":0.031290262829754306,"name":"Primary parallel 0 turn 17","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 5","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,-0.00128296],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.0000696],"layer":"Tertiary layer 2","length":0.020106192982974672,"name":"Tertiary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":0,"section":"Tertiary section 2","winding":"Tertiary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,-0.00128296],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.0000696],"layer":"Tertiary layer 2","length":0.02569822790636449,"name":"Tertiary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":0,"section":"Tertiary section 2","winding":"Tertiary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.0031999999999999993,-0.0014873549999999999],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 6","length":0.020106192982974672,"name":"Primary parallel 0 turn 18","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 6","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004089999999999997,-0.0014873549999999999],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 6","length":0.02569822790636449,"name":"Primary parallel 0 turn 19","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 6","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004979999999999995,-0.0014873549999999999],"crossSectionalShape":"rectangular","dimensions":[0.0007999999999999999,0.00013919],"layer":"Primary layer 6","length":0.031290262829754306,"name":"Primary parallel 0 turn 20","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 6","winding":"Primary"}],"groupsDescription":[{"coordinateSystem":"cartesian","coordinates":[0.0041375,-2.168404344971009e-19],"dimensions":[0.0028749999999999995,0.0031138999999999993],"name":"Default group","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"},{"connections":null,"parallelsProportion":[1],"winding":"Secondary"},{"connections":null,"parallelsProportion":[1],"winding":"Tertiary"}],"sectionsOrientation":"contiguous","type":"Printed"}]},"core":{"functionalDescription":{"gapping":[{"area":0.000021237,"coordinates":[0,0.00025,0],"distanceClosestNormalSurface":0.00165,"distanceClosestParallelSurface":0.0030749999999999996,"length":0.0005,"sectionDimensions":[0.0052,0.0052],"shape":"round","type":"subtractive"},{"area":0.00001125,"coordinates":[0.0063,0,0],"distanceClosestNormalSurface":0.0019,"distanceClosestParallelSurface":0.0030749999999999996,"length":0.000005,"sectionDimensions":[0.00125,0.009],"shape":"irregular","type":"residual"},{"area":0.00001125,"coordinates":[-0.0063,0,0],"distanceClosestNormalSurface":0.0019,"distanceClosestParallelSurface":0.0030749999999999996,"length":0.000005,"sectionDimensions":[0.00125,0.009],"shape":"irregular","type":"residual"}],"material":{"coerciveForce":[{"magneticField":10,"magneticFluxDensity":0,"temperature":100},{"magneticField":15,"magneticFluxDensity":0,"temperature":25}],"curieTemperature":230,"density":4850,"family":"3C","heatCapacity":{"maximum":800,"minimum":700},"heatConductivity":{"maximum":5,"minimum":3.5},"manufacturerInfo":{"datasheetUrl":"https://www.ferroxcube.com/upload/media/product/file/MDS/3c98.pdf","name":"Ferroxcube"},"material":"ferrite","materialComposition":"MnZn","name":"3C98","permeability":{"complex":{"imaginary":[{"frequency":100000,"value":16.305},{"frequency":120000,"value":17.615},{"frequency":150000,"value":19.165},{"frequency":200000,"value":22.715},{"frequency":250000,"value":27.415},{"frequency":300000,"value":33.015},{"frequency":450000,"value":57.26},{"frequency":400000,"value":47.565},{"frequency":500000,"value":69.09},{"frequency":550000,"value":84.59},{"frequency":600000,"value":105.35},{"frequency":650000,"value":133.35000000000002},{"frequency":700000,"value":171.45},{"frequency":750000,"value":223.3},{"frequency":800000,"value":293.1},{"frequency":850000,"value":384.1},{"frequency":900000,"value":497.85},{"frequency":950000,"value":633.45},{"frequency":1000000,"value":786.7},{"frequency":1100000,"value":1123.5},{"frequency":1200000,"value":1473},{"frequency":1300000,"value":1823.5},{"frequency":1400000,"value":2160.5},{"frequency":1500000,"value":2466},{"frequency":1600000,"value":2709},{"frequency":1700000,"value":2862.5},{"frequency":1800000,"value":2920},{"frequency":1900000,"value":2891.5},{"frequency":2000000,"value":2798.5},{"frequency":2200000,"value":2516.5},{"frequency":2400000,"value":2209},{"frequency":2600000,"value":1928.5},{"frequency":2800000,"value":1689},{"frequency":3000000,"value":1489.5},{"frequency":3500000,"value":1124},{"frequency":4000000,"value":886.65},{"frequency":4500000,"value":724.5999999999999},{"frequency":5000000,"value":608.95},{"frequency":6000000,"value":457.25},{"frequency":7000000,"value":365.1},{"frequency":8000000,"value":304.65},{"frequency":9000000,"value":262.65},{"frequency":10000000,"value":231.7},{"frequency":15000000,"value":145.85000000000002},{"frequency":20000000,"value":102.85},{"frequency":25000000,"value":77.84},{"frequency":30000000,"value":61.745000000000005}],"real":[{"frequency":100000,"value":2284},{"frequency":120000,"value":2285},{"frequency":150000,"value":2289.5},{"frequency":200000,"value":2300},{"frequency":250000,"value":2315},{"frequency":300000,"value":2335},{"frequency":450000,"value":2426},{"frequency":400000,"value":2389},{"frequency":500000,"value":2466.5},{"frequency":550000,"value":2514.5},{"frequency":600000,"value":2571},{"frequency":650000,"value":2636},{"frequency":700000,"value":2710},{"frequency":750000,"value":2791.5},{"frequency":800000,"value":2878},{"frequency":850000,"value":2965.5},{"frequency":900000,"value":3048.5},{"frequency":950000,"value":3119},{"frequency":1000000,"value":3172},{"frequency":1100000,"value":3215},{"frequency":1200000,"value":3181},{"frequency":1300000,"value":3074},{"frequency":1400000,"value":2892.5},{"frequency":1500000,"value":2632.5},{"frequency":1600000,"value":2301.5},{"frequency":1700000,"value":1923.5},{"frequency":1800000,"value":1538.5},{"frequency":1900000,"value":1178},{"frequency":2000000,"value":864.3},{"frequency":2200000,"value":397.54999999999995},{"frequency":2400000,"value":108.80000000000001},{"frequency":2600000,"value":-61.405},{"frequency":2800000,"value":-158.9},{"frequency":3000000,"value":-212.45},{"frequency":3500000,"value":-251.9},{"frequency":4000000,"value":-238.89999999999998},{"frequency":4500000,"value":-212.45},{"frequency":5000000,"value":-184.5},{"frequency":6000000,"value":-135.55},{"frequency":7000000,"value":-98.46000000000001},{"frequency":8000000,"value":-71.25999999999999},{"frequency":9000000,"value":-51.66},{"frequency":10000000,"value":-37.7},{"frequency":15000000,"value":-7.8025},{"frequency":20000000,"value":4.0805},{"frequency":25000000,"value":11.735},{"frequency":30000000,"value":17.295}]},"initial":[{"frequency":10000,"temperature":-40,"value":1482},{"frequency":10000,"temperature":-30,"value":1568},{"frequency":10000,"temperature":-20,"value":1696},{"frequency":10000,"temperature":-10,"value":1821},{"frequency":10000,"temperature":0,"value":1958},{"frequency":10000,"temperature":10,"value":2119},{"frequency":10000,"temperature":20,"value":2286},{"frequency":10000,"temperature":30,"value":2507},{"frequency":10000,"temperature":40,"value":2712},{"frequency":10000,"temperature":50,"value":2933},{"frequency":10000,"temperature":60,"value":3170},{"frequency":10000,"temperature":70,"value":3399},{"frequency":10000,"temperature":80,"value":3599},{"frequency":10000,"temperature":90,"value":3796},{"frequency":10000,"temperature":100,"value":3958},{"frequency":10000,"temperature":110,"value":4042},{"frequency":10000,"temperature":120,"value":3975},{"frequency":10000,"temperature":130,"value":3883},{"frequency":10000,"temperature":140,"value":3807},{"frequency":10000,"temperature":150,"value":3694},{"frequency":10000,"temperature":160,"value":3681},{"frequency":10000,"temperature":170,"value":3702},{"frequency":10000,"temperature":180,"value":3752}]},"remanence":[{"magneticField":0,"magneticFluxDensity":0.075,"temperature":100},{"magneticField":0,"magneticFluxDensity":0.14,"temperature":25}],"resistivity":[{"temperature":25,"value":8}],"saturation":[{"magneticField":1200,"magneticFluxDensity":0.44,"temperature":100},{"magneticField":1200,"magneticFluxDensity":0.53,"temperature":25}],"type":"commercial","volumetricLosses":{"default":[{"coefficients":{"excessLossesCoefficient":1.7496785399999968e-21,"resistivityFrequencyCoefficient":1.3082062100000002e-33,"resistivityMagneticFluxDensityCoefficient":0.009345135289999001,"resistivityOffset":0.17955289800000002,"resistivityTemperatureCoefficient":0.000874429148},"method":"roshen"},{"method":"steinmetz","ranges":[{"alpha":1.4399617740036297,"beta":2.8799235480079868,"ct0":3.9697573532165897,"ct1":0.05947093997074883,"ct2":0.00030207891326041445,"k":1.7235233786759798,"maximumFrequency":150000,"minimumFrequency":1},{"alpha":2.1653857592296553,"beta":2.400214766515894,"ct0":1.8680824033596795,"ct1":0.021082654440597735,"ct2":0.00013308219288068598,"k":0.00011523322392865425,"maximumFrequency":1000000,"minimumFrequency":150000},{"alpha":1.96724637893554,"beta":1.6610598090559008,"ct0":1.4273752725474464,"ct1":0.011386084146774614,"ct2":0.00007723962136151744,"k":0.00016592859039323793,"maximumFrequency":1000000000,"minimumFrequency":1000000}]}]}},"numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":{"maximum":0.0141,"minimum":0.0136},"B":{"maximum":0.0032,"minimum":0.0031},"C":{"maximum":0.0092,"minimum":0.0088},"D":{"maximum":0.002,"minimum":0.0018},"E":{"maximum":0.0115,"minimum":0.0112},"F":{"maximum":0.0053,"minimum":0.0051},"G":{"maximum":0.0115,"minimum":0.0112}},"family":"planar er","magneticCircuit":"open","name":"ER 14/4.5/9","type":"standard"},"type":"two-piece set","magneticCircuit":"open"},"geometricalDescription":[{"coordinates":[0,0,0],"machining":[{"coordinates":[0,0.00025,0],"length":0.0005}],"material":"3C98","rotation":[3.141592653589793,3.141592653589793,0],"shape":{"aliases":[],"dimensions":{"A":{"maximum":0.0141,"minimum":0.0136},"B":{"maximum":0.0032,"minimum":0.0031},"C":{"maximum":0.0092,"minimum":0.0088},"D":{"maximum":0.002,"minimum":0.0018},"E":{"maximum":0.0115,"minimum":0.0112},"F":{"maximum":0.0053,"minimum":0.0051},"G":{"maximum":0.0115,"minimum":0.0112}},"family":"planar er","magneticCircuit":"open","name":"ER 14/4.5/9","type":"standard"},"type":"half set"},{"coordinates":[0,0,0],"material":"3C98","rotation":[0,0,0],"shape":{"aliases":[],"dimensions":{"A":{"maximum":0.0141,"minimum":0.0136},"B":{"maximum":0.0032,"minimum":0.0031},"C":{"maximum":0.0092,"minimum":0.0088},"D":{"maximum":0.002,"minimum":0.0018},"E":{"maximum":0.0115,"minimum":0.0112},"F":{"maximum":0.0053,"minimum":0.0051},"G":{"maximum":0.0115,"minimum":0.0112}},"family":"planar er","magneticCircuit":"open","name":"ER 14/4.5/9","type":"standard"},"type":"half set"}],"name":"Custom","processedDescription":{"columns":[{"area":0.000021237,"coordinates":[0,0,0],"depth":0.0052,"height":0.0038,"shape":"round","type":"central","width":0.0052},{"area":0.00001125,"coordinates":[0.0063,0,0],"depth":0.009,"height":0.0038,"minimumWidth":0.00125,"shape":"irregular","type":"lateral","width":0.00125},{"area":0.00001125,"coordinates":[-0.0063,0,0],"depth":0.009,"height":0.0038,"minimumWidth":0.00125,"shape":"irregular","type":"lateral","width":0.00125}],"depth":0.009000000000000001,"effectiveParameters":{"effectiveArea":0.000022112549567984874,"effectiveLength":0.01911852657299317,"effectiveVolume":4.2275936651214745e-7,"minimumArea":0.000021237166338267},"height":0.0063,"width":0.01385,"windingWindows":[{"area":0.000011684999999999999,"coordinates":[0.0026,0],"height":0.0038,"width":0.0030749999999999996}]}},"manufacturerInfo":{"name":"OpenMagnetics","reference":"My custom magnetic"}})"));

    double frequency = 150000;
    double expectedLeakageInductance = 1.4e-6;

    auto leakageInductance = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    auto leakageMagneticField = LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);
    CHECK_THAT(leakageInductance, WithinRel(expectedLeakageInductance, maximumError));

    settings->reset();
}
