#include <source_location>
#include <algorithm>
#include <fstream>
#include "support/Painter.h"
#include "processors/Sweeper.h"
#include "physical_models/Impedance.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/StrayCapacitance.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

double maximumError = 0.25;
TEST_CASE("Test_Impedance_0", "[physical-model][impedance][smoke-test]") {

    std::vector<int64_t> numberTurns = {54, 54};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "T 17/10.7/6.8";
    std::vector<OpenMagnetics::Wire> wires;
    auto wire = find_wire_by_name("Round 0.15 - Grade 1");
    wires.push_back(wire);
    wires.push_back(wire);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment,
                                                     wires,
                                                     false);

    int64_t numberStacks = 1;
    std::string coreMaterial = "80";
    std::vector<CoreGap> gapping = {};
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    double expectedSelfResonantFrequency = 1400000;
    settings._debug = true;
    auto selfResonantFrequency = OpenMagnetics::Impedance().calculate_self_resonant_frequency(magnetic);
    REQUIRE_THAT(expectedSelfResonantFrequency, Catch::Matchers::WithinAbs(selfResonantFrequency, expectedSelfResonantFrequency * maximumError));
    settings._debug = false;

    {
        auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 1000);

        auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;

        outFile.append("Test_Impedance_0.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        #ifdef ENABLE_MATPLOTPP
        painter.paint_curve(impedanceSweep, true);
        #else
        #endif

        #ifdef ENABLE_MATPLOTPP
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        #else
        #endif


    }
    {
        auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_Impedance_0_magnetic.svg";
        outFile.append(filename);
        settings.set_painter_include_fringing(false);
        Painter painter(outFile);

        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        // painter.paint_coil_turns(magnetic);
        #ifdef ENABLE_MATPLOTPP
        painter.export_svg();
        #else
        #endif

    }
}

TEST_CASE("Test_Impedance_Many_Turns", "[physical-model][impedance][smoke-test]") {

    std::vector<int64_t> numberTurns = {110, 110};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "T 12.5/7.5/5";
    std::vector<OpenMagnetics::Wire> wires;
    auto wire = find_wire_by_name("Round 0.15 - Grade 1");
    wires.push_back(wire);
    wires.push_back(wire);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment,
                                                     wires,
                                                     false);

    int64_t numberStacks = 1;
    std::string coreMaterial = "A07";
    std::vector<CoreGap> gapping = {};
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    std::map<double, double> expectedImpedances = {
        {2000, 558},
        {5000, 1350},
        {10000, 2690},
        {25000, 6900},
        {50000, 15900}
    };
    settings._debug = true;

    for (auto [frequency, expectedImpedance] : expectedImpedances) {
        auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
        REQUIRE_THAT(expectedImpedance, Catch::Matchers::WithinAbs(abs(impedance), expectedImpedance * maximumError));
    }

    // {
    //     auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 400000, 1000);

    //     auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
    //     auto outFile = outputFilePath;

    //     outFile.append("Test_Impedance_Many_Turns.svg");
    //     std::filesystem::remove(outFile);
    //     Painter painter(outFile);
    //     painter.paint_curve(impedanceSweep, true);
    //     painter.export_svg();
    //     REQUIRE(std::filesystem::exists(outFile));

    // }
}

TEST_CASE("Test_Self_Resonant_Frequency_Many_Turns", "[physical-model][impedance][smoke-test]") {

    std::vector<int64_t> numberTurns = {110, 110};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "T 12.5/7.5/5";
    std::vector<OpenMagnetics::Wire> wires;
    auto wire = find_wire_by_name("Round 0.15 - Grade 1");
    wires.push_back(wire);
    wires.push_back(wire);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment,
                                                     wires,
                                                     false);

    int64_t numberStacks = 1;
    std::string coreMaterial = "A07";
    std::vector<CoreGap> gapping = {};
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    double expectedSelfResonantFrequency = 180000;
    auto selfResonantFrequency = OpenMagnetics::Impedance().calculate_self_resonant_frequency(magnetic);
    REQUIRE_THAT(expectedSelfResonantFrequency, Catch::Matchers::WithinAbs(selfResonantFrequency, expectedSelfResonantFrequency * maximumError));

}
TEST_CASE("Test_Impedance_Few_Turns", "[physical-model][impedance][smoke-test]") {

    std::vector<int64_t> numberTurns = {18, 18};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "T 12.5/7.5/5";
    std::vector<OpenMagnetics::Wire> wires;
    auto wire = find_wire_by_name("Round 0.425 - Grade 1");
    wires.push_back(wire);
    wires.push_back(wire);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment,
                                                     wires,
                                                     false);

    int64_t numberStacks = 1;
    std::string coreMaterial = "A07";
    std::vector<CoreGap> gapping = {};
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    std::map<double, double> expectedImpedances = {
        {2000, 12.7},
        {5000, 31.8},
        {10000, 62.6},
        {25000, 153},
        {50000, 305},
    };

    for (auto [frequency, expectedImpedance] : expectedImpedances) {
        auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
        REQUIRE_THAT(expectedImpedance, Catch::Matchers::WithinAbs(abs(impedance), expectedImpedance * maximumError));
    }

}

TEST_CASE("Test_Impedance_Many_Turns_Larger_Core", "[physical-model][impedance][smoke-test]") {

    std::vector<int64_t> numberTurns = {9, 9};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "T 36/23/15";
    std::vector<OpenMagnetics::Wire> wires;
    auto wire = find_wire_by_name("Round 2.50 - Grade 1");
    wires.push_back(wire);
    wires.push_back(wire);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment,
                                                     wires,
                                                     false);

    int64_t numberStacks = 1;
    std::string coreMaterial = "A05";
    std::vector<CoreGap> gapping = {};
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    std::map<double, double> expectedImpedances = {
        {2000, 7.49},
        {5000, 19},
        {10000, 37.9},
        {25000, 93.9},
        {50000, 188},
    };

    for (auto [frequency, expectedImpedance] : expectedImpedances) {
        auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
        REQUIRE_THAT(expectedImpedance, Catch::Matchers::WithinAbs(abs(impedance), expectedImpedance * maximumError));
    }
}

TEST_CASE("Test_Impedance_Few_Turns_Larger_Core", "[physical-model][impedance][smoke-test]") {

    std::vector<int64_t> numberTurns = {17, 17};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "T 36/23/15";
    std::vector<OpenMagnetics::Wire> wires;
    auto wire = find_wire_by_name("Round 1.40 - Grade 1");
    wires.push_back(wire);
    wires.push_back(wire);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment,
                                                     wires,
                                                     false);

    int64_t numberStacks = 1;
    std::string coreMaterial = "A05";
    std::vector<CoreGap> gapping = {};
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    std::map<double, double> expectedImpedances = {
        {2000, 21.6},
        {5000, 54.1},
        {10000, 108},
        {25000, 300},
        {50000, 600},
    };

    for (auto [frequency, expectedImpedance] : expectedImpedances) {
        auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
        REQUIRE_THAT(expectedImpedance, Catch::Matchers::WithinAbs(abs(impedance), expectedImpedance * maximumError));
    }

}

TEST_CASE("Test_Differential_Mode_Impedance", "[physical-model][impedance][cmc]") {
    // Two-winding toroid = a common-mode choke. In differential mode the core
    // flux cancels, so the impedance is set by the (much smaller) leakage
    // inductance resonating with the inter-winding capacitance. We assert the
    // qualitative physics: DM impedance is finite, far below the common-mode
    // impedance at low frequency, and inductive (rising) there.
    std::vector<int64_t> numberTurns = {54, 54};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "T 17/10.7/6.8";
    std::vector<OpenMagnetics::Wire> wires;
    auto wire = find_wire_by_name("Round 0.15 - Grade 1");
    wires.push_back(wire);
    wires.push_back(wire);

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                                                     WindingOrientation::CONTIGUOUS, WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires, false);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, std::vector<CoreGap>{}, 1, "80");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double frequency = 100000;
    auto commonMode = abs(OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency));
    auto differentialMode = abs(OpenMagnetics::Impedance().calculate_differential_mode_impedance(magnetic, frequency));

    // DM impedance must be a finite, positive number...
    REQUIRE(std::isfinite(differentialMode));
    REQUIRE(differentialMode > 0);
    // ...and well below the common-mode impedance (leakage ≪ magnetizing L).
    REQUIRE(differentialMode < commonMode);

    // Below its own resonance the DM branch is inductive: |Z| rises with f.
    auto dmLow = abs(OpenMagnetics::Impedance().calculate_differential_mode_impedance(magnetic, 1e4));
    auto dmHigh = abs(OpenMagnetics::Impedance().calculate_differential_mode_impedance(magnetic, 1e5));
    REQUIRE(dmHigh > dmLow);

    // The sweep helper returns a curve of the same length as requested.
    auto sweep = Sweeper().sweep_differential_mode_impedance_over_frequency(magnetic, 1000, 1e8, 50);
    REQUIRE(sweep.get_x_points().size() == 50);
    REQUIRE(sweep.get_y_points().size() == 50);
}

TEST_CASE("Test_Through_Core_Inter_Winding_Capacitance", "[physical-model][impedance][cmc][stray-capacitance]") {
    // The inter-winding capacitance of a separated-winding common-mode choke is the
    // through-core path (turn -> core -> turn): the only capacitive coupling when the
    // two windings have no adjacent turns. Assert the energy/core-potential model gives
    // a finite, positive, physically bounded value, and crucially LESS than the naive
    // parallel-sum series of the two winding-to-core capacitances — the per-turn
    // potential weighting (CPSS core-potential method) must reduce the overestimate.
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "T 20/10/7";
    auto wire = find_wire_by_name("Round 0.15 - Grade 1");
    std::vector<OpenMagnetics::Wire> wires = {wire, wire};

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                                                     WindingOrientation::CONTIGUOUS, WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires, false);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, std::vector<CoreGap>{}, 1, "3C97");
    coil.wind();

    auto primaryName = coil.get_functional_description()[0].get_name();
    auto secondaryName = coil.get_functional_description()[1].get_name();
    std::map<std::string, double> voltageRmsPerWinding = {{primaryName, 10.0}, {secondaryName, 10.0}};
    auto voltagesPerTurn = StrayCapacitance::calculate_voltages_per_turn(coil, voltageRmsPerWinding).get_voltage_per_turn().value();

    double throughCore = StrayCapacitance::calculate_through_core_capacitance(coil, core, primaryName, secondaryName, voltagesPerTurn);
    double primaryToCore = StrayCapacitance::calculate_winding_to_core_capacitance(coil, core, primaryName);
    double secondaryToCore = StrayCapacitance::calculate_winding_to_core_capacitance(coil, core, secondaryName);
    double naiveSeries = (primaryToCore * secondaryToCore) / (primaryToCore + secondaryToCore);

    REQUIRE(std::isfinite(throughCore));
    REQUIRE(throughCore > 0.0);
    REQUIRE(throughCore < 1e-9);             // physically bounded (sub-nF), no divergence
    REQUIRE(throughCore < naiveSeries);      // potential weighting reduces vs the naive sum
}

// ABT #167: the datasheet/REDEXPERT "typical impedance" of a common-mode choke is
// the CM measurement (windings driven in parallel), which excites no leakage
// resonance. sweep_impedance_over_frequency (one winding driven, other floating)
// adds a leakage tank damped only by the mΩ winding resistance, producing a
// near-undamped MΩ spike that measured CM curves do not show. The CM sweep must
// place its single (magnetizing) resonance near the measured one.
// Reference: WE 744834622 measured CM peak 56 kOhm at 0.299 MHz (REDEXPERT .s4p).
TEST_CASE("Test_Impedance_Common_Mode_No_Leakage_Spike", "[physical-model][impedance]") {
    settings.reset();
    auto testDataPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_redexpert_744834405.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    auto magneticJson = nlohmann::json::parse(file);
    OpenMagnetics::Magnetic magnetic(magneticJson);
    magnetic = magnetic_autocomplete(magnetic);

    auto curve = Sweeper::sweep_common_mode_impedance_over_frequency(magnetic, 10e3, 1e9, 300);
    auto frequencies = curve.get_x_points();
    auto impedances = curve.get_y_points();

    size_t peakIndex = 0;
    for (size_t i = 0; i < impedances.size(); ++i) {
        if (impedances[i] > impedances[peakIndex]) {
            peakIndex = i;
        }
    }

    // Single magnetizing resonance near the measured 0.299 MHz (material data is
    // reverse-engineered, so the window is generous but excludes the 1.68 MHz
    // leakage-tank artifact the terminal sweep used to serve as the peak).
    CHECK(frequencies[peakIndex] > 0.15e6);
    CHECK(frequencies[peakIndex] < 0.7e6);
    // No near-undamped leakage spike (the terminal sweep peaked at 1.17 MOhm).
    CHECK(impedances[peakIndex] < 500e3);
    // Capacitive rolloff after the resonance: the top decade must be far below the peak.
    CHECK(impedances.back() < 0.01 * impedances[peakIndex]);
}

// ABT #167: complex-permeability splines are interpolators — past the last measured
// point they diverged polynomially (A10's data ends at 1.3 MHz; µ'' extrapolated to
// -2.5e6 at 1 GHz, i.e. an ACTIVE element), making CMC impedance sweeps rise
// monotonically to 1 GHz instead of rolling off. Out-of-span queries must clamp to
// the nearest measured endpoint: µ'' stays non-negative and the CM curve rolls off.
TEST_CASE("Test_Impedance_Complex_Permeability_No_Extrapolation", "[physical-model][impedance]") {
    settings.reset();
    OpenMagnetics::ComplexPermeability complexPermeabilityModel;
    for (double frequency : {2e6, 1e7, 1e8, 1e9}) {
        auto [real, imaginary] = complexPermeabilityModel.get_complex_permeability(std::string("A10"), frequency);
        CHECK(real >= 1.0);
        CHECK(imaginary >= 0.0);
    }

    auto testDataPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_redexpert_744834622.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    auto magneticJson = nlohmann::json::parse(file);
    OpenMagnetics::Magnetic magnetic(magneticJson);
    magnetic = magnetic_autocomplete(magnetic);

    auto curve = Sweeper::sweep_common_mode_impedance_over_frequency(magnetic, 10e3, 1e9, 300);
    auto impedances = curve.get_y_points();
    double peak = *std::max_element(impedances.begin(), impedances.end());
    // Used to rise monotonically to the 1 GHz end of the sweep; now it must roll off.
    CHECK(impedances.back() < 0.2 * peak);
}

}  // namespace
