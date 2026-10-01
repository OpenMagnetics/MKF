#include <source_location>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <fstream>
#include <sstream>
#include "support/Painter.h"
#include "processors/Sweeper.h"
#include "constructive_models/Bobbin.h"
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
    // CHARACTERISATION, not a physics anchor (re-pinned 2026-10-02). The 1.4 MHz this test used to
    // check was measured on a part whose core coating is unknown, while the coating thickness sets
    // the turn-to-core gap that dominates this self-capacitance; so the measurement cannot
    // discriminate the capacitance physics. This core declares no coating, so it is wound on bare
    // ferrite (an undeclared coating is 0), which puts the turns closest to the core and gives the
    // lowest self-resonance. Value: MKF output with the screened turn-to-core elements and the
    // 2*eps0*lt*Y1 Albach pair element on that bare core.
    // The measured anchor for a toroidal CMC with a documented coating is
    // "Toroidal CMC common-mode resonance against its s4p measurement (WE 744824220)".
    double expectedSelfResonantFrequency = 0.72966e6;
    settings._debug = true;
    auto selfResonantFrequency = OpenMagnetics::Impedance().calculate_self_resonant_frequency(magnetic);
    REQUIRE_THAT(selfResonantFrequency, Catch::Matchers::WithinRel(expectedSelfResonantFrequency, 0.02));
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

    // The real part: WE-CMB 744821039, 110 + 110 turns of Round 0.15 - Grade 1 on an ACME A07
    // T 12.7/7.92/4.9 ring with a 0.6 mm epoxy coating, read from its stored MAS (the same file the
    // coated-anchor test below winds).
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_we_coated_anchor/744821039.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    OpenMagnetics::Magnetic magnetic(nlohmann::json::parse(file));
    magnetic = magnetic_autocomplete(magnetic);
    REQUIRE(magnetic.get_coil().get_turns_description());

    // Anchor: 137.246 kHz, the |Z_CM| peak of the part's WE S-parameter file (744821039.s4p, rev20a,
    // 26/02/2020, mixed-mode Scc21; the value the coated-anchor table below uses). The single-winding
    // SRF with the other winding open is the common-mode resonance: unity coupling makes the open
    // winding mirror the driven one (see Impedance::calculate_self_resonant_frequency).
    //
    // Tolerance, from the measurements alone. The part has a second record: the Heimdall measurement
    // DB (measurements.measurement, origin CMChokesPower, common_mode_impedance, two bit-identical
    // rows) peaks at 125.893 kHz. Nothing records whether the two are the same sample, so their gap,
    // |ln(137246/125893)| = 0.0863 (9.0%), is the measured spread of this part's SRF. Both are read
    // on the same logarithmic grid of ratio 1.01742 (1.74% per point), each peak +-half a step, so
    // the two readings add one full step, ln(1.01742) = 0.0173. The band is ln(f/137246) within
    // +-0.1036, i.e. 123.7 to 152.2 kHz; it contains the measDB peak by construction.
    const double s4pPeak = 137246;
    const double measurementDatabasePeak = 125893;
    const double gridStep = std::log(1.01742);
    const double tolerance = std::abs(std::log(s4pPeak / measurementDatabasePeak)) + gridStep;
    auto selfResonantFrequency = OpenMagnetics::Impedance().calculate_self_resonant_frequency(magnetic);
    INFO("model SRF " << selfResonantFrequency << " Hz, ln(model/anchor) " << std::log(selfResonantFrequency / s4pPeak) << ", tolerance " << tolerance);
    REQUIRE(std::abs(std::log(selfResonantFrequency / s4pPeak)) <= tolerance);

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
    // Expected |Z| are MKF values, re-pinned when ring cores moved to the IEC 60205 clause 5.1 effective
    // parameters (ABT #1502): le/Ae of T 36/23/15 went 0.9506 -> 0.9350 mm^-1 (TDK publishes 0.94 for R36),
    // so L and |Z| rose 1.67% at every point (before: 26.94 / 67.13 / 134.3 / 335.1 / 675.5 Ohm).
    // The measured references were 21.6 / 54.1 / 108 / 300 / 600 Ohm at 2 / 5 / 10 / 25 / 50 kHz;
    // MKF is +27% above the 2 kHz measurement and +14% above the 50 kHz one.
    std::map<double, double> expectedImpedances = {
        {2000, 27.39},
        {5000, 68.25},
        {10000, 136.5},
        {25000, 340.7},
        {50000, 686.8},
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


TEST_CASE("Test_Complex_Permeability_Keeps_Falling_Above_The_Anchor", "[physical-model][impedance]") {
    // ABT #843. get_complex_permeability CLAMPS its interpolation to the tabulated
    // frequency range, so a table that stopped at 100x the anchor handed back the same
    // mu' and mu'' for every frequency above it. Nanoperm 80000 anchors at 18.8 kHz, so
    // everything above 1.88 MHz read as mu' = 19210.6, mu'' = 13647.3 — identical at
    // 4.33 MHz, 10 MHz, 50 MHz and 100 MHz, while the material's own initial-permeability
    // table falls to 156 by 25.6 MHz. Downstream that moved a common-mode choke's
    // modelled self-resonance from ~50 MHz to 4.33 MHz and overstated its peak impedance
    // 12x.
    settings.reset();
    OpenMagnetics::ComplexPermeability complexPermeabilityModel;
    auto coreMaterial = OpenMagnetics::find_core_material_by_name("Nanoperm 80000");

    std::vector<double> frequencies = {1e6, 4.33e6, 1e7, 2.56e7, 5e7, 1e8};
    std::vector<double> magnitudes;
    for (auto frequency : frequencies) {
        auto [real, imaginary] = complexPermeabilityModel.get_complex_permeability(coreMaterial, frequency);
        magnitudes.push_back(sqrt(real * real + imaginary * imaginary));
    }

    // It must keep falling, not freeze: every step strictly below the one before it.
    for (size_t index = 1; index < magnitudes.size(); ++index) {
        CHECK(magnitudes[index] < magnitudes[index - 1]);
    }

    // And it must track the material's own curve, not merely differ. Nanoperm 80000's
    // table falls 2,378 -> 156 between 1 MHz and 25.6 MHz; |mu| now reads 2,822 -> 184
    // across the same span (the excess over mu' is mu''), against 30,116 -> 23,565 before.
    CHECK(magnitudes.back() < 0.1 * magnitudes.front());
    CHECK(magnitudes[3] < 400.0);   // 25.6 MHz, the material's last tabulated point (156)
}


TEST_CASE("PROBE_CMC622_CM_Curve_Dump", "[probe622]") {
    settings.reset();
    auto testDataPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_redexpert_744834622.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    auto magneticJson = nlohmann::json::parse(file);
    OpenMagnetics::Magnetic magnetic(magneticJson);
    magnetic = magnetic_autocomplete(magnetic);
    auto curve = Sweeper::sweep_common_mode_impedance_over_frequency(magnetic, 1e4, 1e8, 120);
    std::ofstream out("/tmp/mkf_zcm_622.txt");
    auto fs = curve.get_x_points();
    auto zs = curve.get_y_points();
    for (size_t i = 0; i < fs.size(); ++i) out << fs[i] << " " << zs[i] << "\n";
    out.close();
    CHECK(fs.size() == 120);
}

// ABT #383: a core carrying only its functionalDescription is a normal thing to hand around —
// it is how MAS files are written when the constructive description is the source of truth and
// the processed values are meant to be derived. Impedance reads the PROCESSED effective area and
// length, and without them it returned garbage that looked like an answer: NaN for drumRing,
// exactly 0 for drumSemishielded, ~1e-10 ohm for molded. The NaN then surfaced far away as
// "Waveform data contains NaN" out of the waveform processor, which points at the wrong
// component entirely — the reporter filed a model bug against impedance because of it.
TEST_CASE("Test_Impedance_Refuses_Unprocessed_Core", "[physical-model][impedance]") {
    settings.reset();
    clear_databases();

    auto processedCore = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::array(), 1, "3C97");
    auto coil = OpenMagneticsTesting::get_quick_coil({10}, {1}, "E 42/21/20");
    OpenMagnetics::Magnetic processedMagnetic;
    processedMagnetic.set_core(processedCore);
    processedMagnetic.set_coil(coil);

    // With the processed description present it answers, and the answer is finite and positive.
    OpenMagnetics::Impedance impedanceModel;
    auto processedImpedance = impedanceModel.calculate_impedance(processedMagnetic, 100000);
    UNSCOPED_INFO("processed core impedance " << std::abs(processedImpedance));
    CHECK(std::isfinite(std::abs(processedImpedance)));
    CHECK(std::abs(processedImpedance) > 0);

    // Strip the processed description and it must refuse rather than answer with NaN/0/1e-10.
    auto unprocessedCore = processedCore;
    unprocessedCore.set_processed_description(std::nullopt);
    OpenMagnetics::Magnetic unprocessedMagnetic;
    unprocessedMagnetic.set_core(unprocessedCore);
    unprocessedMagnetic.set_coil(coil);
    std::string message;
    try {
        impedanceModel.calculate_impedance(unprocessedMagnetic, 100000);
        message = "no exception";
    }
    catch (const std::exception& exception) {
        message = exception.what();
    }
    UNSCOPED_INFO("unprocessed core -> " << message);
    CHECK(message.find("processed description") != std::string::npos);
    settings.reset();
}

// Data-driven competitor cable-core impedance export. Reads /tmp/fr_parts.json
// [{mpn, material, od, id, h}] (toroid dims in mm), computes |Z|(f) 1 MHz..1 GHz
// for a 1-turn core in MKF, writes /tmp/fr_curves.json [{mpn, f[], z[]}]. Lets the
// ingest compute vendor-agnostic curves from material + geometry (the OM way).
TEST_CASE("Competitor_Cable_Core_Impedance_Export", "[cable-core-export]") {
    std::ifstream in("/tmp/fr_parts.json");
    if (!in.good()) { WARN("no /tmp/fr_parts.json — skipping"); return; }
    json parts = json::parse(in);
    json out = json::array();
    for (auto& p : parts) {
        json c;
        c["mpn"] = p.value("mpn", std::string(""));
        try {
            double od = p.at("od").get<double>() / 1000.0;  // mm -> m
            double id = p.at("id").get<double>() / 1000.0;
            double h  = p.at("h").get<double>() / 1000.0;
            std::string material = p.at("material").get<std::string>();
            int turns = p.value("turns", 1);
            json shapeJson = {
                {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "t"},
                {"aliases", json::array()}, {"name", "cable core"},
                {"dimensions", {{"A", {{"nominal", od}}}, {"B", {{"nominal", id}}}, {"C", {{"nominal", h}}}}}
            };
            json coreJson;
            coreJson["functionalDescription"] = {
                {"type", "toroidal"}, {"material", material}, {"shape", shapeJson},
                {"gapping", json::array()}, {"numberStacks", 1}
            };
            Core core(coreJson);
            auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);
            json coilJson;
            to_json(coilJson["bobbin"], bobbin);
            coilJson["functionalDescription"] = json::array();
            json wj;
            wj["name"] = "cable"; wj["numberTurns"] = turns; wj["numberParallels"] = 1;
            wj["isolationSide"] = "primary"; wj["wire"] = "Round 0.475 - Grade 1";
            coilJson["functionalDescription"].push_back(wj);
            OpenMagnetics::Coil coil(coilJson);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            auto curve = Sweeper().sweep_impedance_over_frequency(magnetic, 1e6, 1e9, 150);
            c["f"] = curve.get_x_points();
            c["z"] = curve.get_y_points();
        } catch (const std::exception& e) {
            c["error"] = e.what();
        }
        out.push_back(c);
    }
    std::ofstream of("/tmp/fr_curves.json");
    of << out.dump();
    of.close();
    CHECK(out.size() == parts.size());
}

}  // namespace

TEST_CASE("Test_Core_Dimensional_Attenuation_From_Real_Permittivity", "[physical-model][impedance][smoke-test]") {
    // ABT #848: a MnZn ferrite core is a lossy dielectric with eps' ~ 1e5 (Ferroxcube handbook
    // Table 5) and a conductor (resistivity), so across a 15 mm cross-section the wave is
    // attenuated by 2-4x in the 0.3-1 MHz band. The factor must be ~1 at audio/low-RF, fall to
    // the 0.2-0.5 range around 1 MHz on 15 mm, and be EXACTLY 1 for a material without
    // permittivity data (no invented correction).
    auto a07 = find_core_material_by_name("A07");
    REQUIRE(a07.get_permittivity());
    std::vector<double> dims = {0.005, 0.015};  // T25x15x15 cross-section: 5 mm radial, 15 mm high
    auto mu10k = OpenMagnetics::ComplexPermeability().get_complex_permeability(a07, 1e4);
    auto f10k = OpenMagnetics::Impedance::core_dimensional_attenuation(a07, 1e4, std::complex<double>(mu10k.first, mu10k.second), dims);
    CHECK(std::abs(f10k) > 0.97);
    auto mu1m = OpenMagnetics::ComplexPermeability().get_complex_permeability(a07, 1e6);
    auto f1m = OpenMagnetics::Impedance::core_dimensional_attenuation(a07, 1e6, std::complex<double>(mu1m.first, mu1m.second), dims);
    // Physical check: on a 5 x 15 mm section the field enters mainly through the faces 5 mm
    // apart, and the 15 mm faces only ADD penetration — so the 2D result must be at least the
    // slab factor of the 5 mm thickness (which ignores the extra faces), well above the slab
    // factor of the 15 mm thickness, and below 1.
    auto slab5 = OpenMagnetics::Impedance::core_dimensional_attenuation(a07, 1e6, std::complex<double>(mu1m.first, mu1m.second), {0.005});
    auto slab15 = OpenMagnetics::Impedance::core_dimensional_attenuation(a07, 1e6, std::complex<double>(mu1m.first, mu1m.second), {0.015});
    CHECK(std::abs(f1m) >= std::abs(slab5) * 0.98);
    CHECK(std::abs(f1m) > std::abs(slab15));
    CHECK(std::abs(f1m) < 1.0);
    // monotone: stronger attenuation at higher frequency
    auto mu300k = OpenMagnetics::ComplexPermeability().get_complex_permeability(a07, 3e5);
    auto f300k = OpenMagnetics::Impedance::core_dimensional_attenuation(a07, 3e5, std::complex<double>(mu300k.first, mu300k.second), dims);
    CHECK(std::abs(f300k) > std::abs(f1m));
    // a material WITHOUT permittivity data gets no correction at all
    auto n87 = find_core_material_by_name("N87");
    REQUIRE(!n87.get_permittivity());
    auto fNone = OpenMagnetics::Impedance::core_dimensional_attenuation(n87, 1e6, std::complex<double>(1000, 100), dims);
    CHECK(std::abs(fNone - 1.0) < 1e-12);
}

// Toroidal CMC inter-winding (DM) capacitance with neighbour-screened turn-to-core elements.
// WE 744822222 (T 14/8/9 A07 MnZn in a 0.6 mm case, 2 x 18 turns of 0.5 mm, sectored). Measured
// (REDEXPERT .s4p): DM resonance 40.27 MHz with L_DM 9.93 uH at 1 MHz -> 1.57 pF; a complex
// (R+jwL)||C||Rp fit around the resonance gives 1.97 pF. The isolated (Smythe) turn-to-core
// element gave 4.54 pF (2.3-2.9x high): the close-wound bore row (pitch = wire OD) screens each
// turn's flanks, 23.6 vs 50.9 pF/m. Screened: 2.55 pF (1.3-1.6x). The CM tank capacitance moves
// with it (8.43 -> 6.79 pF), CM peak 0.73 -> 0.98 MHz against 0.93 MHz measured.
TEST_CASE("Toroidal CMC inter-winding capacitance with screened turn-to-core elements (WE 744822222)", "[physical-model][impedance][cmc][stray-capacitance]") {
    auto testDataPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_we_744822222_enriched.json");
    std::ifstream file(testDataPath);
    OpenMagnetics::Magnetic magnetic(nlohmann::json::parse(file));
    magnetic = magnetic_autocomplete(magnetic);
    auto parameters = OpenMagnetics::Impedance().calculate_differential_mode_parameters(magnetic, 1e6);
    CHECK_THAT(parameters.interWindingCapacitance, Catch::Matchers::WithinRel(2.551e-12, 0.02));
    // Within 2x of both measured values, and well below the unscreened 4.54 pF.
    CHECK(parameters.interWindingCapacitance < 2 * 1.57e-12);
    CHECK(parameters.interWindingCapacitance > 0.5 * 1.97e-12);

    auto commonModeModel = OpenMagnetics::Impedance().build_common_mode_impedance_model(magnetic);
    // Characterisation (re-pinned 2026-10-01): 6.795 pF before the Albach pair element was corrected
    // from (2/3)*eps0*lt*Y1 to 2*eps0*lt*Y1 (Albach 2017 eq. 3.14) and toroid turn pairs were split
    // between the bore and the outer-crossing gap. The s4p-derived C_cm of this part is 13.4 pF,
    // but its CM peak is flagged weak/low-Q, so it is not used as an anchor here.
    CHECK_THAT(commonModeModel.tanks[0].capacitance, Catch::Matchers::WithinRel(14.033e-12, 0.02));

    // The DM path honours the magnetic's core electrical reference: a GROUNDED core diverts the
    // through-core path to the reference, and these sectored windings have no adjacent turns,
    // so nothing couples them directly.
    CoreElectricalReference groundedReference;
    groundedReference.set_type(CoreElectricalReferenceType::GROUNDED);
    auto grounded = magnetic;
    grounded.set_core_electrical_reference(groundedReference);
    auto groundedParameters = OpenMagnetics::Impedance().calculate_differential_mode_parameters(grounded, 1e6);
    CHECK(groundedParameters.interWindingCapacitance == 0.0);
    // Floating (absent reference) through the Core/Coil overload equals the Magnetic overload.
    auto coreCoilParameters = OpenMagnetics::Impedance().calculate_differential_mode_parameters(magnetic.get_core(), magnetic.get_coil(), 1e6);
    CHECK_THAT(coreCoilParameters.interWindingCapacitance, Catch::Matchers::WithinRel(parameters.interWindingCapacitance, 1e-12));
}

// MEASURED ANCHOR for the common-mode self-capacitance of epoxy-coated toroidal CMCs.
//
// Parts: every WE-CMB catalogue part that (a) the WE requirements sheet gives a 0.6 mm epoxy core
// coating, (b) has a WE s4p measurement whose CM peak is a clean LC resonance (measDB
// cm_quality_flag 'good': no weak/low-Q/multi-peak/mu-rolloff flag), and (c) MKF winds. 17 parts
// meet (a)+(b); 744821240 is left out by (c) (magnetic_autocomplete leaves it unwound). The core
// is ACME A07 in all 16. No part was selected by its result.
// Measured values (cmc_impedance_whitepaper/work/measdb/ciw_measured.csv): CM impedance peak
// frequency and CM inductance at 10 kHz, both from the WE s4p file. The s4p grid is logarithmic
// with ratio 1.01742 and every peak sits on a grid point: +-0.87% (half a step) reading error.
//
// What is compared: f_res ~ 1/sqrt(L*C), so the modelled peak is rescaled by sqrt(L_model/L_meas)
// at 10 kHz. That removes the A07 initial-permeability tolerance (7000 +-25%, ACME catalogue,
// i.e. up to +-15% on f) using each sample's own measured inductance and leaves the frequency
// error to the capacitance (and to the shape of mu(f), assumed common to the lot).
//
// Coating: the requirements sheet states 0.6 mm; ACME's catalogue states "0.6 mm max" for T9 and
// above. No nominal is documented. A thinner coating raises the capacitance and LOWERS f, so the
// documented bounds are 0 (no coating) and 0.6 mm, and the model must bracket the measurement:
//   median ln(f(0.6 mm)/f_meas) >= -tol   and   median ln(f(0 mm)/f_meas) <= +tol.
// Tolerance: two standard errors of the median of the per-part log ratios (1.2533 s/sqrt(n), s the
// sample standard deviation; ~95% for a normal spread), plus the s4p half-step. The per-part
// scatter (s ~0.13) is not covered by any documented input uncertainty; it is reported, not tested.
// The bracket alone does not reject the pre-fix element (its 0 mm median lands on the
// measurement); the third check, at the requirements-sheet value, does: with (2/3)*eps0*lt*Y1 and
// unsplit toroid pairs the 0.6 mm median sits 28% high (19% before the turn-to-core screening),
// against 3% low now, with a tolerance of ~9%.
TEST_CASE("Toroidal CMC common-mode resonance against its s4p measurement (WE-CMB, 0.6 mm epoxy coating)", "[physical-model][impedance][cmc][stray-capacitance][measured-anchor]") {
    struct MeasuredPart { std::string partNumber; double peakFrequency; double inductanceAt10kHz; };
    const std::vector<MeasuredPart> parts = {
        {"744821039", 137246, 0.044589},   {"744821110", 331131, 0.00976392}, {"744821120", 215030, 0.0211613},
        {"744821150", 444120, 0.00500236}, {"744822110", 331131, 0.00985155}, {"744822120", 197242, 0.0214343},
        {"744823210", 226464, 0.0105858},  {"744823220", 174783, 0.023393},   {"744823305", 336899, 0.00523203},
        {"744823333", 467735, 0.00360972}, {"744824220", 128086, 0.0214685},  {"744824310", 278612, 0.00816787},
        {"744824407", 273842, 0.00692025}, {"744824433", 501187, 0.00311928}, {"744825320", 152230, 0.0202871},
        {"744825433", 113501, 0.0295728},
    };
    const double s4pHalfStep = std::log(std::sqrt(1.01742));

    auto peakOf = [](const Curve2D& curve) {
        auto x = curve.get_x_points();
        auto y = curve.get_y_points();
        return x[std::max_element(y.begin(), y.end()) - y.begin()];
    };
    // ln(f_model * sqrt(L_model / L_meas) / f_meas) for one part at one coating thickness.
    auto logRatio = [&](const MeasuredPart& part, double coatingThickness) {
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_we_coated_anchor/" + part.partNumber + ".json");
        std::ifstream file(path);
        REQUIRE(file.good());
        auto json = nlohmann::json::parse(file);
        json["core"]["functionalDescription"]["coating"]["thickness"] = {{"nominal", coatingThickness}};
        OpenMagnetics::Magnetic magnetic(json);
        magnetic = magnetic_autocomplete(magnetic);
        REQUIRE(magnetic.get_coil().get_turns_description());
        double coarsePeak = peakOf(Sweeper::sweep_common_mode_impedance_over_frequency(magnetic, 1e3, 1e9, 400, "log"));
        double peak = peakOf(Sweeper::sweep_common_mode_impedance_over_frequency(magnetic, coarsePeak / 1.1, coarsePeak * 1.1, 401, "linear"));
        auto lowFrequency = Sweeper::sweep_common_mode_impedance_over_frequency(magnetic, 1e4, 2e4, 2, "linear");
        double modelInductance = lowFrequency.get_y_points()[0] / (2 * std::numbers::pi * lowFrequency.get_x_points()[0]);
        return std::log(peak * std::sqrt(modelInductance / part.inductanceAt10kHz) / part.peakFrequency);
    };
    auto median = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        size_t n = v.size();
        return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    };
    auto tolerance = [&](const std::vector<double>& v) {
        double mean = 0;
        for (double x : v) mean += x;
        mean /= static_cast<double>(v.size());
        double variance = 0;
        for (double x : v) variance += (x - mean) * (x - mean);
        double sampleDeviation = std::sqrt(variance / static_cast<double>(v.size() - 1));
        return 2.0 * 1.2533 * sampleDeviation / std::sqrt(static_cast<double>(v.size())) + s4pHalfStep;
    };

    std::vector<double> atSheetValue, uncoated;
    for (const auto& part : parts) {
        atSheetValue.push_back(logRatio(part, 0.6e-3));
        uncoated.push_back(logRatio(part, 0.0));
        UNSCOPED_INFO(part.partNumber << ": model/measured " << std::exp(atSheetValue.back()) << " at 0.6 mm, " << std::exp(uncoated.back()) << " uncoated");
    }
    double medianAtSheetValue = median(atSheetValue);
    double medianUncoated = median(uncoated);
    double toleranceAtSheetValue = tolerance(atSheetValue);
    double toleranceUncoated = tolerance(uncoated);
    UNSCOPED_INFO("median model/measured: " << std::exp(medianAtSheetValue) << " at 0.6 mm (tol " << toleranceAtSheetValue
                  << "), " << std::exp(medianUncoated) << " uncoated (tol " << toleranceUncoated << ")");
    // The documented coating bounds bracket the measurement.
    CHECK(medianAtSheetValue >= -toleranceAtSheetValue);
    CHECK(medianUncoated <= toleranceUncoated);
    // At the requirements-sheet value itself.
    CHECK(std::abs(medianAtSheetValue) <= toleranceAtSheetValue);
}

// The DM port of a two-winding choke (windings in series opposition, far ends joined) charges
// more than the inter-winding path: the turns of each winding sit at the DM ramp, so the
// turn-to-turn capacitances INSIDE each winding store energy too, and so does every turn's element
// to the core. calculate_differential_mode_capacitance must be exactly the energy method over all
// of them, C_DM = 2 W / V_port^2, with the first winding's turns at +v_i, the second's at -v_j:
//   floating core: W = sum_pairs 1/2 C_ij dV^2 + calculate_winding_pair_to_core_energy (offset 0);
//   grounded core: W = sum_pairs 1/2 C_ij dV^2 + 1/2 sum_i C_i V_i^2 (each winding's self energy
//   against a core held at 0, the sign of -v_j squaring away).
// And the DM impedance must resonate the leakage against C_DM, not against C_iw.
TEST_CASE("Test_Impedance_Differential_Mode_Capacitance_Is_The_DM_Port_Energy", "[physical-model][impedance][cmc][stray-capacitance][captot]") {
    auto testDataPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_we_744822222_enriched.json");
    std::ifstream file(testDataPath);
    OpenMagnetics::Magnetic magnetic(nlohmann::json::parse(file));
    magnetic = magnetic_autocomplete(magnetic);
    auto coil = magnetic.get_coil();
    auto core = magnetic.get_core();
    if (!coil.get_turns_description()) {
        coil.wind();
    }
    auto windings = coil.get_functional_description();
    REQUIRE(windings.size() == 2);
    auto firstName = windings[0].get_name();
    auto secondName = windings[1].get_name();
    StrayCapacitance strayCapacitance(Settings::GetInstance().get_stray_capacitance_model());

    std::map<std::string, double> voltageRmsPerWinding;
    for (auto& winding : windings) {
        voltageRmsPerWinding[winding.get_name()] = 10.0 * winding.get_number_turns() / windings[0].get_number_turns();
    }
    auto voltages = StrayCapacitance::calculate_voltages_per_turn(coil, voltageRmsPerWinding).get_voltage_per_turn().value();
    auto turns = coil.get_turns_description().value();
    std::vector<double> portPotential(turns.size());
    double maximumFirst = std::numeric_limits<double>::lowest();
    double maximumSecond = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < turns.size(); ++i) {
        bool first = turns[i].get_winding() == firstName;
        portPotential[i] = first ? voltages[i] : -voltages[i];
        (first ? maximumFirst : maximumSecond) = std::max(first ? maximumFirst : maximumSecond, voltages[i]);
    }
    const double portVoltage = maximumFirst + maximumSecond;
    double turnToTurnEnergy = 0;
    size_t intraWindingPairs = 0;
    for (const auto& [key, capacitance] : strayCapacitance.calculate_capacitance_among_turns(coil)) {
        if (key.first < key.second) {
            turnToTurnEnergy += 0.5 * capacitance * std::pow(portPotential[key.first] - portPotential[key.second], 2);
            intraWindingPairs += turns[key.first].get_winding() == turns[key.second].get_winding();
        }
    }
    // The DM ramp puts real energy inside each winding: this is what C_iw leaves out.
    REQUIRE(intraWindingPairs > 0);
    REQUIRE(turnToTurnEnergy > 0);

    double floatingCoreEnergy = StrayCapacitance::calculate_winding_pair_to_core_energy(coil, core, firstName, secondName, voltages, 0.0);
    double floating = strayCapacitance.calculate_differential_mode_capacitance(coil, core);
    CHECK_THAT(floating, Catch::Matchers::WithinRel(2 * (turnToTurnEnergy + floatingCoreEnergy) / (portVoltage * portVoltage), 1e-9));

    CoreElectricalReference groundedReference;
    groundedReference.set_type(CoreElectricalReferenceType::GROUNDED);
    double groundedCoreEnergy = StrayCapacitance::calculate_winding_to_core_self_energy(coil, core, firstName, voltages, std::nullopt, 0.0) +
                                StrayCapacitance::calculate_winding_to_core_self_energy(coil, core, secondName, voltages, std::nullopt, 0.0);
    double grounded = strayCapacitance.calculate_differential_mode_capacitance(coil, core, std::nullopt, groundedReference);
    CHECK_THAT(grounded, Catch::Matchers::WithinRel(2 * (turnToTurnEnergy + groundedCoreEnergy) / (portVoltage * portVoltage), 1e-9));

    // The DM parameters carry C_DM (second pass at its own resonance), keep C_iw beside it, and
    // the DM impedance is the leakage branch in parallel with C_DM.
    auto parameters = OpenMagnetics::Impedance().calculate_differential_mode_parameters(magnetic, 1e6);
    double resonance = 1.0 / (2 * std::numbers::pi * std::sqrt(parameters.leakageInductance *
                                                              strayCapacitance.calculate_differential_mode_capacitance(coil, core)));
    CHECK_THAT(parameters.differentialModeCapacitance,
               Catch::Matchers::WithinRel(strayCapacitance.calculate_differential_mode_capacitance(coil, core, resonance), 1e-9));
    CHECK(parameters.differentialModeCapacitance > parameters.interWindingCapacitance);
    for (double frequency : {1e6, 3e7, 1e8}) {
        double omega = 2 * std::numbers::pi * frequency;
        std::complex<double> inductive(parameters.windingResistance, omega * parameters.leakageInductance);
        std::complex<double> capacitive(0, -1.0 / (omega * parameters.differentialModeCapacitance));
        auto expected = 1.0 / (1.0 / inductive + 1.0 / capacitive);
        CHECK_THAT(std::abs(OpenMagnetics::Impedance().differential_mode_impedance_from_parameters(parameters, frequency)),
                   Catch::Matchers::WithinRel(std::abs(expected), 1e-12));
    }
    std::cout << "744822222: C_iw " << parameters.interWindingCapacitance * 1e12 << " pF, C_DM " << parameters.differentialModeCapacitance * 1e12
              << " pF (measured DM resonance 40.27 MHz -> 1.57 pF; fit 1.97 pF)" << std::endl;
}
