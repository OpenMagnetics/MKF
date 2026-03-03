#include <source_location>
#include "support/Painter.h"
#include "converter_models/CommonModeChoke.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "processors/NgspiceRunner.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>
#include <numeric>
#include <cmath>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}
    .parent_path().append("..").append("output");
[[maybe_unused]] double maximumError = 0.1;

// ─── helper: build a minimal valid CMC json ───────────────────────────
static json makeCmcJson(
    double voltageNominal        = 230.0,
    double current               = 5.0,
    double lineFreq              = 50.0,
    int    numWindings           = 2,
    double impedanceAtFreq_Z     = 1000.0,
    double impedanceAtFreq_f     = 150e3)
{
    json j;
    j["operatingVoltage"]  = {{"nominal", voltageNominal}};
    j["operatingCurrent"]  = current;
    j["lineFrequency"]     = lineFreq;
    j["ambientTemperature"]= 25.0;
    j["lineImpedance"]     = 50.0;
    j["numberOfWindings"]  = numWindings;
    j["minimumImpedance"]  = json::array();
    j["minimumImpedance"].push_back({{"frequency", impedanceAtFreq_f},
                                     {"impedance",  impedanceAtFreq_Z}});
    return j;
}

// ═══════════════════════════════════════════════════════════════════════
// Static helper unit tests
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_StaticHelpers", "[converter-model][cmc-topology][unit-test]") {

    SECTION("impedanceToInductance: Z=1000Ω, f=150kHz → L≈1.061 mH") {
        double L = Cmc::impedanceToInductance(1000.0, 150e3);
        double expected = 1000.0 / (2.0 * M_PI * 150e3);
        CHECK_THAT(L, Catch::Matchers::WithinRel(expected, 1e-6));
    }

    SECTION("inductanceToImpedance round-trip") {
        double L = Cmc::impedanceToInductance(750.0, 200e3);
        double Z = Cmc::inductanceToImpedance(L, 200e3);
        CHECK_THAT(Z, Catch::Matchers::WithinRel(750.0, 1e-6));
    }

    SECTION("insertionLossToImpedance: 20 dB, 50 Ω → 50*(10-1)=450 Ω") {
        double Z = Cmc::insertionLossToImpedance(20.0, 50.0);
        CHECK_THAT(Z, Catch::Matchers::WithinRel(450.0, 1e-4));
    }

    SECTION("insertionLossToImpedance: 40 dB → 50*(100-1)=4950 Ω") {
        double Z = Cmc::insertionLossToImpedance(40.0, 50.0);
        CHECK_THAT(Z, Catch::Matchers::WithinRel(4950.0, 1e-4));
    }
}

// ═══════════════════════════════════════════════════════════════════════
// run_checks
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_RunChecks", "[converter-model][cmc-topology]") {

    SECTION("Valid input passes") {
        Cmc cmc(makeCmcJson());
        CHECK(cmc.run_checks(false) == true);
    }

    SECTION("Empty minimumImpedance fails") {
        json j = makeCmcJson();
        j["minimumImpedance"] = json::array();
        Cmc cmc(j);
        CHECK(cmc.run_checks(false) == false);
    }

    SECTION("Zero operating current fails") {
        json j = makeCmcJson();
        j["operatingCurrent"] = 0;
        Cmc cmc(j);
        CHECK(cmc.run_checks(false) == false);
    }

    SECTION("Negative operating current fails") {
        json j = makeCmcJson();
        j["operatingCurrent"] = -1.0;
        Cmc cmc(j);
        CHECK(cmc.run_checks(false) == false);
    }

    SECTION("Zero line frequency fails") {
        json j = makeCmcJson();
        j["lineFrequency"] = 0;
        Cmc cmc(j);
        CHECK(cmc.run_checks(false) == false);
    }

    SECTION("numberOfWindings=5 fails") {
        json j = makeCmcJson();
        j["numberOfWindings"] = 5;
        Cmc cmc(j);
        CHECK(cmc.run_checks(false) == false);
    }

    SECTION("numberOfWindings=1 fails") {
        json j = makeCmcJson();
        j["numberOfWindings"] = 1;
        Cmc cmc(j);
        CHECK(cmc.run_checks(false) == false);
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Single-phase 2-winding CMC
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_SinglePhase", "[converter-model][cmc-topology][smoke-test]") {
    Cmc cmc(makeCmcJson(230, 5, 50, 2, 1000, 150e3));

    SECTION("run_checks passes") {
        CHECK(cmc.run_checks(false) == true);
    }

    SECTION("Design requirements: 1 turns ratio = 1.0") {
        auto req = cmc.process_design_requirements();
        REQUIRE(req.get_turns_ratios().size() == 1);
        double n = resolve_dimensional_values(req.get_turns_ratios()[0]);
        CHECK_THAT(n, Catch::Matchers::WithinRel(1.0, maximumError));
    }

    SECTION("Computed inductance ≥ L(1000Ω, 150kHz)") {
        double expectedL = Cmc::impedanceToInductance(1000.0, 150e3);
        CHECK(cmc.get_computed_inductance() >= expectedL * 0.99);
    }

    SECTION("Dominant frequency = 150 kHz") {
        CHECK_THAT(cmc.get_dominant_frequency(), Catch::Matchers::WithinRel(150e3, 1e-6));
    }

    SECTION("process() produces 1 operating point with 2 excitations") {
        Cmc cmc2(makeCmcJson());
        auto inputs = cmc2.process();
        CHECK(inputs.get_operating_points().size() == 1);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 2);
    }

    SECTION("Winding names are Line and Neutral") {
        auto names = Cmc::windingNames(2);
        CHECK(names[0] == "Line");
        CHECK(names[1] == "Neutral");
    }

    SECTION("Plot Line current waveform") {
        Cmc cmc2(makeCmcJson());
        auto inputs = cmc2.process();
        auto outFile = outputFilePath;
        outFile.append("Test_Cmc_SinglePhase_Line_Current.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(inputs.get_operating_points()[0]
            .get_excitations_per_winding()[0].get_current()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Plot Line voltage waveform") {
        Cmc cmc2(makeCmcJson());
        auto inputs = cmc2.process();
        auto outFile = outputFilePath;
        outFile.append("Test_Cmc_SinglePhase_Line_Voltage.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(inputs.get_operating_points()[0]
            .get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
        painter.export_svg();
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Multiple impedance points — hardest constraint wins
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_MultipleImpedancePoints", "[converter-model][cmc-topology]") {
    json j = makeCmcJson();
    // Add a second, less demanding point
    j["minimumImpedance"].push_back({{"frequency", 1e6}, {"impedance", 500}});
    j["minimumImpedance"].push_back({{"frequency", 10e6}, {"impedance", 100}});

    Cmc cmc(j);

    double L_150k = Cmc::impedanceToInductance(1000, 150e3);   // most demanding
    double L_1M   = Cmc::impedanceToInductance(500,  1e6);
    double L_10M  = Cmc::impedanceToInductance(100,  10e6);
    double expectedMax = std::max({L_150k, L_1M, L_10M});

    CHECK_THAT(cmc.get_computed_inductance(), Catch::Matchers::WithinRel(expectedMax, 1e-6));
    // 150 kHz must dominate
    CHECK_THAT(cmc.get_dominant_frequency(), Catch::Matchers::WithinRel(150e3, 1e-6));
}

// ═══════════════════════════════════════════════════════════════════════
// Insertion-loss specification
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_InsertionLossSpec", "[converter-model][cmc-topology]") {
    json j = makeCmcJson();
    j["minimumImpedance"] = json::array();   // no direct impedance
    j["targetInsertionLoss"] = json::array();
    j["targetInsertionLoss"].push_back({{"frequency", 150e3}, {"insertionLoss", 30}});

    // 30 dB, 50 Ω → Z = 50*(10^1.5 - 1) ≈ 1531 Ω → L ≈ 1.626 mH
    double expectedZ = Cmc::insertionLossToImpedance(30.0, 50.0);
    double expectedL = Cmc::impedanceToInductance(expectedZ, 150e3);

    Cmc cmc(j);
    CHECK_THAT(cmc.get_computed_inductance(), Catch::Matchers::WithinRel(expectedL, 1e-4));

    auto inputs = cmc.process();
    CHECK(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 2);
}

// ═══════════════════════════════════════════════════════════════════════
// Three-phase delta (3 windings)
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_ThreePhaseDelta", "[converter-model][cmc-topology]") {
    Cmc cmc(makeCmcJson(400, 10, 50, 3, 600, 150e3));

    SECTION("3 excitations per operating point") {
        auto inputs = cmc.process();
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 3);
    }

    SECTION("2 turns ratios all = 1.0") {
        auto req = cmc.process_design_requirements();
        REQUIRE(req.get_turns_ratios().size() == 2);
        for (const auto& tr : req.get_turns_ratios())
            CHECK_THAT(resolve_dimensional_values(tr), Catch::Matchers::WithinRel(1.0, maximumError));
    }

    SECTION("Winding names: Phase A, B, C") {
        auto names = Cmc::windingNames(3);
        CHECK(names[0] == "Phase A");
        CHECK(names[1] == "Phase B");
        CHECK(names[2] == "Phase C");
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Three-phase + neutral wye (4 windings)
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_ThreePhaseNeutral", "[converter-model][cmc-topology]") {
    Cmc cmc(makeCmcJson(400, 16, 50, 4, 800, 150e3));

    SECTION("4 excitations per operating point") {
        auto inputs = cmc.process();
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 4);
    }

    SECTION("3 turns ratios all = 1.0") {
        auto req = cmc.process_design_requirements();
        REQUIRE(req.get_turns_ratios().size() == 3);
        for (const auto& tr : req.get_turns_ratios())
            CHECK_THAT(resolve_dimensional_values(tr), Catch::Matchers::WithinRel(1.0, maximumError));
    }

    SECTION("Last winding name is Neutral") {
        auto names = Cmc::windingNames(4);
        CHECK(names[3] == "Neutral");
    }
}

// ═══════════════════════════════════════════════════════════════════════
// AdvancedCmc — user-specified inductance
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_AdvancedCmc", "[converter-model][cmc-topology]") {
    json j = makeCmcJson();
    j["desiredInductance"] = 2e-3;   // 2 mH — override computed value

    AdvancedCmc adv(j);
    adv._assertErrors = false;
    auto inputs = adv.process();

    SECTION("Magnetizing inductance = user-specified 2 mH") {
        double Lm = resolve_dimensional_values(
            inputs.get_design_requirements().get_magnetizing_inductance());
        CHECK_THAT(Lm, Catch::Matchers::WithinRel(2e-3, maximumError));
    }

    SECTION("2 windings in operating point") {
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 2);
    }

    SECTION("Operating point name is Nominal") {
        CHECK(inputs.get_operating_points()[0].get_name().value_or("") == "Nominal");
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Realistic CMC Simulation with Line + Switching Noise
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_RealisticSimulation", "[converter-model][cmc-topology][simulation]") {
    // Skip if ngspice is not available
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available");
    }

    json j = makeCmcJson(230.0, 5.0, 50.0, 2, 1000.0, 150e3);
    CommonModeChoke cmc(j);
    
    double inductance = 1e-3; // 1 mH
    double parasiticCap_pF = 100.0; // 100 pF
    double dvdt_V_ns = 1.0; // 1 V/ns
    
    SECTION("Realistic simulation completes without timeout") {
        // This should complete in under 10 seconds now
        auto operatingPoints = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        
        CHECK(!operatingPoints.empty());
        CHECK(operatingPoints[0].get_excitations_per_winding().size() == 2);
    }
    
    SECTION("Excitations have both voltage and current") {
        auto operatingPoints = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!operatingPoints.empty());
        
        auto excitations = operatingPoints[0].get_excitations_per_winding();
        for (const auto& exc : excitations) {
            CHECK(exc.get_voltage().has_value());
            CHECK(exc.get_current().has_value());
            
            // Voltage should be around line voltage (230V RMS ~ 325V peak)
            if (exc.get_voltage().has_value() && exc.get_voltage().value().get_waveform()) {
                auto vWf = exc.get_voltage().value().get_waveform().value();
                auto data = vWf.get_data();
                if (!data.empty()) {
                    double vPeak = *std::max_element(data.begin(), data.end());
                    CHECK(vPeak > 100.0); // Should be at least 100V
                }
            }
        }
    }
    
    SECTION("Waveforms have proper time data") {
        auto operatingPoints = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!operatingPoints.empty());
        
        auto excitations = operatingPoints[0].get_excitations_per_winding();
        for (const auto& exc : excitations) {
            if (exc.get_current().has_value() && exc.get_current().value().get_waveform()) {
                auto iWf = exc.get_current().value().get_waveform().value();
                auto data = iWf.get_data();
                auto time = iWf.get_time();
                
                CHECK(!data.empty());
                CHECK(time.has_value());
                CHECK(time.value().size() == data.size());
            }
        }
    }
}

} // anonymous namespace
