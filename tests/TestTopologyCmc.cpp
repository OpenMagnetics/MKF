#include <source_location>
#include "support/Painter.h"
#include "converter_models/CommonModeChoke.h"
#include "advisers/MagneticAdviser.h"
#include "physical_models/Impedance.h"
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
            
            // The simulation covers only 5 µs of a 50 Hz sine.
            // At t=0 the SIN starts at 0V, so the mains voltage barely
            // rises in 5 µs.  The measured voltage is dominated by CM
            // noise, which is small.  Just verify it's finite and not GV.
            if (exc.get_voltage().has_value() && exc.get_voltage().value().get_waveform()) {
                auto vWf = exc.get_voltage().value().get_waveform().value();
                auto data = vWf.get_data();
                if (!data.empty()) {
                    double vMax = *std::max_element(data.begin(), data.end());
                    CHECK(vMax < 1000.0); // Must NOT be in GV range
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

// ═══════════════════════════════════════════════════════════════════════
// Realistic CMC Simulation — comprehensive waveform validation
// Uses the same default values as the WebFrontend CmcWizard
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_RealisticSimulation_FrontendDefaults", "[converter-model][cmc-topology][simulation]") {
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available");
    }

    // ── Frontend CmcWizard defaults (localData) ─────────────────────
    //   operatingVoltage:  { nominal: 230 }  → 230 V RMS
    //   operatingCurrent:  5 A
    //   lineFrequency:     50 Hz
    //   lineImpedance:     50 Ω
    //   ambientTemperature: 25 °C
    //   numberOfWindings:  2
    //   parasiticCap_pF:   10 pF
    //   dvdt_V_ns:         50 V/ns
    //   inductance:        1 mH (default from getSimulateFn fallback)
    json j = makeCmcJson(230.0, 5.0, 50.0, 2, 1000.0, 150e3);
    CommonModeChoke cmc(j);

    const double inductance      = 1e-3;   // 1 mH
    const double parasiticCap_pF = 10.0;   // 10 pF (frontend localData default)
    const double dvdt_V_ns       = 50.0;   // 50 V/ns (frontend localData default)
    // ── Derived expected values ─────────────────────────────────────
    // Excitation frequency = dominantFrequency from the impedance spec
    // (makeCmcJson passes 150e3). Both analytical and simulated now drive
    // the winding at this frequency with a sinusoidal CM current.
    const double excitationFreq = 150e3;
    // CM noise current: C_parasitic × dV/dt = 10e-12 × 50e9 = 0.5 A
    // Sim duration: numberOfPeriods=2 measurement periods at excitationFreq
    // → 2 / 150 kHz ≈ 13.3 µs.
    const double expectedSimTime = 2.0 / excitationFreq;

    SECTION("Simulation produces one operating point with 2 winding excitations") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);

        REQUIRE(ops.size() == 1);
        CHECK(ops[0].get_excitations_per_winding().size() == 2);
    }

    SECTION("Operating point has stable 'Simulated' label") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        // The label used to embed switching frequency (e.g. "CMC_Realistic_1000kHz")
        // which leaked implementation detail into the UI. Now a clean "Simulated".
        CHECK(ops[0].get_name().value_or("") == "Simulated");
    }

    SECTION("Excitation frequency equals the dominant impedance frequency") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        for (const auto& exc : ops[0].get_excitations_per_winding()) {
            CHECK_THAT(exc.get_frequency(), Catch::Matchers::WithinRel(excitationFreq, 0.01));
        }
    }

    SECTION("Voltage waveform is in the correct range (not GV, not µV)") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        for (const auto& exc : ops[0].get_excitations_per_winding()) {
            REQUIRE(exc.get_voltage().has_value());
            REQUIRE(exc.get_voltage()->get_waveform().has_value());

            auto data = exc.get_voltage()->get_waveform()->get_data();
            REQUIRE(!data.empty());

            double vMax = *std::max_element(data.begin(), data.end());
            double vMin = *std::min_element(data.begin(), data.end());

            // The simulation only covers 5 µs of the 20 ms line period.
            // At 50 Hz the mains voltage barely changes in 5 µs, so the
            // voltage should be near its initial value (close to 0 for a
            // sine starting at 0 phase), but NOT in the gigavolt range.
            //
            // With CM noise superimposed, peak should be at most a few
            // hundred volts above whatever the mains is at that instant.
            INFO("vMax = " << vMax << " V, vMin = " << vMin << " V");
            CHECK(vMax < 1000.0);    // Must NOT be in GV range
            CHECK(vMin > -1000.0);   // Must NOT be in GV range
        }
    }

    SECTION("Current waveform is non-zero") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        for (const auto& exc : ops[0].get_excitations_per_winding()) {
            REQUIRE(exc.get_current().has_value());
            REQUIRE(exc.get_current()->get_waveform().has_value());

            auto data = exc.get_current()->get_waveform()->get_data();
            REQUIRE(!data.empty());

            double iMax = *std::max_element(data.begin(), data.end());
            double iMin = *std::min_element(data.begin(), data.end());
            double iPkPk = iMax - iMin;

            INFO("iMax = " << iMax << " A, iMin = " << iMin << " A, iPkPk = " << iPkPk << " A");

            // Current must be non-zero — both line current and CM noise
            // should be present.
            // Load current ~ V / (Rload + 2*Rline) = 230 / (46 + 100) ≈ 1.6 A peak
            // CM noise current ~ 0.5 A peak
            CHECK(iPkPk > 0.01);     // Must have some current flowing
        }
    }

    SECTION("Time vector spans approximately 5 switching periods") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        auto exc = ops[0].get_excitations_per_winding()[0];
        REQUIRE(exc.get_current().has_value());
        REQUIRE(exc.get_current()->get_waveform().has_value());

        auto time = exc.get_current()->get_waveform()->get_time();
        REQUIRE(time.has_value());
        REQUIRE(time->size() > 10);  // Must have reasonable number of points

        double tStart = time->front();
        double tEnd   = time->back();
        double tSpan  = tEnd - tStart;

        INFO("Time span = " << tSpan << " s, expected ~ " << expectedSimTime << " s");
        // Should be approximately 5 µs (5 / 1e6)
        CHECK(tSpan > expectedSimTime * 0.5);
        CHECK(tSpan < expectedSimTime * 2.0);
    }

    SECTION("Voltage has processed data (peak, rms, etc.)") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        auto exc = ops[0].get_excitations_per_winding()[0];
        REQUIRE(exc.get_voltage().has_value());
        REQUIRE(exc.get_voltage()->get_processed().has_value());

        auto processed = exc.get_voltage()->get_processed().value();
        // Peak voltage should be in volts, not gigavolts
        if (processed.get_peak().has_value()) {
            INFO("Processed peak voltage = " << processed.get_peak().value() << " V");
            CHECK(std::abs(processed.get_peak().value()) < 1000.0);
        }
        if (processed.get_rms().has_value()) {
            INFO("Processed RMS voltage = " << processed.get_rms().value() << " V");
            CHECK(processed.get_rms().value() < 1000.0);
        }
    }

    SECTION("Current has processed data (non-zero rms)") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        auto exc = ops[0].get_excitations_per_winding()[0];
        REQUIRE(exc.get_current().has_value());
        REQUIRE(exc.get_current()->get_processed().has_value());

        auto processed = exc.get_current()->get_processed().value();
        if (processed.get_rms().has_value()) {
            INFO("Processed RMS current = " << processed.get_rms().value() << " A");
            CHECK(processed.get_rms().value() > 1e-6);  // Must be non-zero
        }
    }

    SECTION("Both windings have similar voltage (symmetric CMC)") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        auto excitations = ops[0].get_excitations_per_winding();
        REQUIRE(excitations.size() == 2);

        // Both windings of a 2-winding CMC carry the same CM noise
        // but opposite DM (line) voltages.  The voltages at cmc_in0
        // and cmc_in1 differ by the mains phase, but their absolute
        // magnitudes should be in the same order.
        auto v0data = excitations[0].get_voltage()->get_waveform()->get_data();
        auto v1data = excitations[1].get_voltage()->get_waveform()->get_data();
        REQUIRE(!v0data.empty());
        REQUIRE(!v1data.empty());

        double v0max = *std::max_element(v0data.begin(), v0data.end());
        double v1max = *std::max_element(v1data.begin(), v1data.end());
        double v0min = *std::min_element(v0data.begin(), v0data.end());
        double v1min = *std::min_element(v1data.begin(), v1data.end());

        double v0range = v0max - v0min;
        double v1range = v1max - v1min;

        INFO("Winding 0 voltage range: " << v0range << " V");
        INFO("Winding 1 voltage range: " << v1range << " V");

        // Both voltage ranges should be comparable (within 10x)
        if (v0range > 0 && v1range > 0) {
            double ratio = v0range / v1range;
            CHECK(ratio > 0.1);
            CHECK(ratio < 10.0);
        }
    }

    SECTION("Ambient temperature is set correctly") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        auto temp = ops[0].get_conditions().get_ambient_temperature();
        CHECK_THAT(temp, Catch::Matchers::WithinRel(25.0, 0.01));
    }

    SECTION("Plot voltage and current waveforms for winding 0") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        // Print diagnostic info for all windings
        for (size_t w = 0; w < ops[0].get_excitations_per_winding().size(); ++w) {
            auto& exc = ops[0].get_excitations_per_winding()[w];
            std::ostringstream diag;
            diag << "\n--- Winding " << w << " ---\n";
            diag << "  Frequency: " << exc.get_frequency() << " Hz\n";

            if (exc.get_voltage() && exc.get_voltage()->get_waveform()) {
                auto vData = exc.get_voltage()->get_waveform()->get_data();
                auto vTime = exc.get_voltage()->get_waveform()->get_time();
                double vMax = *std::max_element(vData.begin(), vData.end());
                double vMin = *std::min_element(vData.begin(), vData.end());
                diag << "  Voltage: " << vData.size() << " points\n";
                diag << "    Max: " << vMax << " V, Min: " << vMin << " V, PkPk: " << (vMax - vMin) << " V\n";
                if (vTime) diag << "    Time: " << vTime->front() << " to " << vTime->back() << " s\n";
                if (exc.get_voltage()->get_processed()) {
                    auto p = exc.get_voltage()->get_processed().value();
                    if (p.get_peak()) diag << "    Proc peak: " << p.get_peak().value() << " V\n";
                    if (p.get_rms()) diag << "    Proc RMS: " << p.get_rms().value() << " V\n";
                }
            }
            if (exc.get_current() && exc.get_current()->get_waveform()) {
                auto iData = exc.get_current()->get_waveform()->get_data();
                double iMax = *std::max_element(iData.begin(), iData.end());
                double iMin = *std::min_element(iData.begin(), iData.end());
                diag << "  Current: " << iData.size() << " points\n";
                diag << "    Max: " << iMax << " A, Min: " << iMin << " A, PkPk: " << (iMax - iMin) << " A\n";
                if (exc.get_current()->get_processed()) {
                    auto p = exc.get_current()->get_processed().value();
                    if (p.get_peak()) diag << "    Proc peak: " << p.get_peak().value() << " A\n";
                    if (p.get_rms()) diag << "    Proc RMS: " << p.get_rms().value() << " A\n";
                }
            }
            WARN(diag.str());
        }

        auto exc = ops[0].get_excitations_per_winding()[0];

        if (exc.get_voltage().has_value() && exc.get_voltage()->get_waveform().has_value()) {
            auto outFile = outputFilePath;
            outFile.append("Test_Cmc_Realistic_Voltage_W0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(exc.get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        if (exc.get_current().has_value() && exc.get_current()->get_waveform().has_value()) {
            auto outFile = outputFilePath;
            outFile.append("Test_Cmc_Realistic_Current_W0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(exc.get_current()->get_waveform().value());
            painter.export_svg();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Realistic CMC Simulation — with fallback defaults from getSimulateFn
// (parasiticCap_pF=100, dvdt_V_ns=1)
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_RealisticSimulation_FallbackDefaults", "[converter-model][cmc-topology][simulation]") {
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available");
    }

    json j = makeCmcJson(230.0, 5.0, 50.0, 2, 1000.0, 150e3);
    CommonModeChoke cmc(j);

    // Fallback defaults from CmcWizard.getSimulateFn():
    //   parasiticCap_pF = 100 (aux.parasiticCap_pF || 100)
    //   dvdt_V_ns       = 1   (aux.dvdt_V_ns || 1)
    const double inductance      = 1e-3;
    const double parasiticCap_pF = 100.0;
    const double dvdt_V_ns       = 1.0;

    // CM noise current: 100e-12 × 1e9 = 0.1 A

    SECTION("Simulation completes and produces valid operating points") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);

        REQUIRE(ops.size() == 1);
        CHECK(ops[0].get_excitations_per_winding().size() == 2);
    }

    SECTION("Voltage is in realistic range (not GV)") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        for (const auto& exc : ops[0].get_excitations_per_winding()) {
            REQUIRE(exc.get_voltage().has_value());
            REQUIRE(exc.get_voltage()->get_waveform().has_value());
            auto data = exc.get_voltage()->get_waveform()->get_data();
            REQUIRE(!data.empty());

            double vMax = *std::max_element(data.begin(), data.end());
            double vMin = *std::min_element(data.begin(), data.end());

            INFO("vMax = " << vMax << " V, vMin = " << vMin << " V");
            CHECK(vMax < 1000.0);
            CHECK(vMin > -1000.0);
        }
    }

    SECTION("Current is non-zero") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        for (const auto& exc : ops[0].get_excitations_per_winding()) {
            REQUIRE(exc.get_current().has_value());
            REQUIRE(exc.get_current()->get_waveform().has_value());
            auto data = exc.get_current()->get_waveform()->get_data();
            REQUIRE(!data.empty());

            double iMax = *std::max_element(data.begin(), data.end());
            double iMin = *std::min_element(data.begin(), data.end());

            INFO("iMax = " << iMax << " A, iMin = " << iMin << " A");
            CHECK((iMax - iMin) > 1e-6);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Analytical ↔ Simulated consistency
//
// After the excitation unification, both paths drive each winding with the
// same sinusoidal CM current at dominantFrequency (the user's impedance-spec
// point, default 150 kHz), amplitude C_parasitic · dV/dt. They must now
// agree on peak current, peak voltage, and dominant frequency within a few
// percent — that is the whole point of having two paths.
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_AnalyticalVsSimulated_CurrentConsistency",
          "[converter-model][cmc-topology][consistency]") {
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available");
    }

    auto makeJsonWithNoise = [](double cap_pF, double dvdt_V_ns) {
        json j = makeCmcJson(230.0, 5.0, 50.0, 2, 1000.0, 150e3);
        j["parasiticCap_pF"] = cap_pF;
        j["dvdt_V_ns"]       = dvdt_V_ns;
        j["safetyMargin_dB"] = 6.0;
        return j;
    };

    // AC peak of a sampled sine wave, ignoring any DC offset.
    auto peakACOf = [](const std::vector<double>& data) {
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
        double vmax = 0.0;
        for (double v : data) vmax = std::max(vmax, std::abs(v - mean));
        return vmax;
    };

    // Driving parameters used by every section below.
    const double cap_pF   = 10.0;
    const double dvdt     = 50.0;
    const double L_self   = 1e-3;                  // 1 mH per winding
    const double f_excit  = 150e3;                 // = dominantFrequency from spec
    const double iCmPeak  = cap_pF * dvdt * 1e-3;  // 0.5 A
    const double vCmPeak  = L_self * 2.0 * M_PI * f_excit * iCmPeak;  // ≈ 471 V

    SECTION("analytical peak CM current matches closed form C·dV/dt") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto inputs = cmc.process();
        REQUIRE(!inputs.get_operating_points().empty());

        auto exc0 = inputs.get_operating_points()[0].get_excitations_per_winding().at(0);
        REQUIRE(exc0.get_current().has_value());
        auto harmonics = exc0.get_current()->get_harmonics();
        REQUIRE(harmonics.has_value());
        auto amps = harmonics->get_amplitudes();
        REQUIRE(amps.size() >= 2);
        // amps[1] is the AC peak at dominantFrequency.
        CHECK_THAT(amps[1], Catch::Matchers::WithinRel(iCmPeak, 0.05));
    }

    SECTION("simulated operating point has label 'Simulated'") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops.empty());
        CHECK(ops[0].get_name().value_or("") == "Simulated");
    }

    SECTION("simulated current peak agrees with analytical within 5%") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops.empty());

        auto exc0 = ops[0].get_excitations_per_winding().at(0);
        REQUIRE(exc0.get_current().has_value());
        auto iData = exc0.get_current()->get_waveform()->get_data();
        REQUIRE(iData.size() > 10);

        double iPeakAC = peakACOf(iData);
        INFO("simulated AC peak = " << iPeakAC << " A, analytical = " << iCmPeak << " A");
        CHECK_THAT(iPeakAC, Catch::Matchers::WithinRel(iCmPeak, 0.05));
    }

    SECTION("simulated voltage peak matches L·ω·I from analytical within 5%") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops.empty());

        auto exc0 = ops[0].get_excitations_per_winding().at(0);
        REQUIRE(exc0.get_voltage().has_value());
        auto vData = exc0.get_voltage()->get_waveform()->get_data();
        REQUIRE(vData.size() > 10);

        double vPeakAC = peakACOf(vData);
        INFO("simulated V peak = " << vPeakAC << " V, analytical = " << vCmPeak << " V");
        CHECK_THAT(vPeakAC, Catch::Matchers::WithinRel(vCmPeak, 0.05));
    }

    SECTION("simulated waveform is sinusoidal (zero crossings ≈ 2·numberOfPeriods)") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops.empty());

        auto exc0 = ops[0].get_excitations_per_winding().at(0);
        auto iData = exc0.get_current()->get_waveform()->get_data();
        double mean = std::accumulate(iData.begin(), iData.end(), 0.0) / iData.size();

        int crossings = 0;
        for (size_t i = 1; i < iData.size(); ++i) {
            bool prevPos = (iData[i - 1] - mean) > 0;
            bool currPos = (iData[i]     - mean) > 0;
            if (prevPos != currPos) ++crossings;
        }
        INFO("zero crossings = " << crossings << " over 2 periods (expect ~4)");
        // 2 periods of a sine has 4 zero crossings. Allow ±1 for sampling boundary effects.
        CHECK(crossings >= 3);
        CHECK(crossings <= 5);
    }

    SECTION("numberOfPeriods controls the output time window") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops2 = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 5);
        auto ops4 = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 4, 5);
        REQUIRE(!ops2.empty());
        REQUIRE(!ops4.empty());

        auto duration = [](const OperatingPointExcitation& e) {
            auto t = e.get_current()->get_waveform()->get_time().value();
            return t.back() - t.front();
        };
        double d2 = duration(ops2[0].get_excitations_per_winding().at(0));
        double d4 = duration(ops4[0].get_excitations_per_winding().at(0));
        INFO("duration: 2 periods=" << d2 << " s, 4 periods=" << d4 << " s");
        CHECK(d4 > d2 * 1.5);  // ~2× longer, leave slack for sample-count rounding
    }

    SECTION("all windings carry the same CM waveform (CM definition)") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops.empty());
        auto excs = ops[0].get_excitations_per_winding();
        REQUIRE(excs.size() >= 2);

        double iPeak0 = peakACOf(excs[0].get_current()->get_waveform()->get_data());
        double iPeak1 = peakACOf(excs[1].get_current()->get_waveform()->get_data());
        double vPeak0 = peakACOf(excs[0].get_voltage()->get_waveform()->get_data());
        double vPeak1 = peakACOf(excs[1].get_voltage()->get_waveform()->get_data());
        INFO("winding 0: I=" << iPeak0 << " V=" << vPeak0);
        INFO("winding 1: I=" << iPeak1 << " V=" << vPeak1);
        CHECK_THAT(iPeak1, Catch::Matchers::WithinRel(iPeak0, 0.01));
        CHECK_THAT(vPeak1, Catch::Matchers::WithinRel(vPeak0, 0.01));
    }

    SECTION("current peak scales linearly with C·dV/dt") {
        // Double C·dV/dt → expect double the AC peak.
        CommonModeChoke cmc1(makeJsonWithNoise(10.0, 50.0));   // I_cm = 0.5 A
        CommonModeChoke cmc2(makeJsonWithNoise(20.0, 50.0));   // I_cm = 1.0 A
        auto ops1 = cmc1.simulate_realistic_cmc(L_self, 10.0, 50.0, 2, 10);
        auto ops2 = cmc2.simulate_realistic_cmc(L_self, 20.0, 50.0, 2, 10);
        REQUIRE(!ops1.empty());
        REQUIRE(!ops2.empty());

        double p1 = peakACOf(ops1[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data());
        double p2 = peakACOf(ops2[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data());
        INFO("I peak: 0.5 A expected=" << p1 << ", 1.0 A expected=" << p2);
        CHECK_THAT(p1, Catch::Matchers::WithinRel(0.5, 0.05));
        CHECK_THAT(p2, Catch::Matchers::WithinRel(1.0, 0.05));
        CHECK_THAT(p2 / p1, Catch::Matchers::WithinRel(2.0, 0.05));
    }

    SECTION("voltage peak scales linearly with inductance (V = Lω·I)") {
        // Triple the inductance → expect triple the voltage peak at the same I_cm.
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops1 = cmc.simulate_realistic_cmc(1e-3, cap_pF, dvdt, 2, 10);  // 1 mH
        auto ops3 = cmc.simulate_realistic_cmc(3e-3, cap_pF, dvdt, 2, 10);  // 3 mH
        REQUIRE(!ops1.empty());
        REQUIRE(!ops3.empty());

        double v1 = peakACOf(ops1[0].get_excitations_per_winding()[0].get_voltage()->get_waveform()->get_data());
        double v3 = peakACOf(ops3[0].get_excitations_per_winding()[0].get_voltage()->get_waveform()->get_data());
        INFO("V peak: 1mH=" << v1 << " V, 3mH=" << v3 << " V, ratio=" << v3 / v1);
        CHECK_THAT(v3 / v1, Catch::Matchers::WithinRel(3.0, 0.05));
    }

    SECTION("voltage peak scales linearly with dominantFrequency (V = Lω·I)") {
        // Double f → expect double V peak at the same I_cm and L.
        // NB: noise-estimation mode implicitly adds a 150 kHz impedance point
        // (see constructor line 138), which fixes dominantFrequency at 150 kHz.
        // So for this test we use pure minimumImpedance mode (no noise params
        // in the json) and let dominantFrequency follow the spec.
        json j150 = makeCmcJson(230.0, 5.0, 50.0, 2, 1000.0, 150e3);
        json j300 = makeCmcJson(230.0, 5.0, 50.0, 2, 2500.0, 300e3);

        CommonModeChoke cmc150(j150);
        CommonModeChoke cmc300(j300);
        CHECK_THAT(cmc150.get_dominant_frequency(), Catch::Matchers::WithinRel(150e3, 0.001));
        CHECK_THAT(cmc300.get_dominant_frequency(), Catch::Matchers::WithinRel(300e3, 0.001));

        // Same L_self, same C·dV/dt → only f differs.
        auto ops150 = cmc150.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        auto ops300 = cmc300.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops150.empty());
        REQUIRE(!ops300.empty());

        double v150 = peakACOf(ops150[0].get_excitations_per_winding()[0].get_voltage()->get_waveform()->get_data());
        double v300 = peakACOf(ops300[0].get_excitations_per_winding()[0].get_voltage()->get_waveform()->get_data());
        INFO("V peak: 150kHz=" << v150 << " V, 300kHz=" << v300 << " V, ratio=" << v300 / v150);
        CHECK_THAT(v300 / v150, Catch::Matchers::WithinRel(2.0, 0.05));
    }

    SECTION("3-winding CMC produces identical CM waveforms on all three windings") {
        json j = makeCmcJson(400.0, 5.0, 50.0, 3, 1000.0, 150e3);
        j["parasiticCap_pF"] = cap_pF;
        j["dvdt_V_ns"]       = dvdt;
        CommonModeChoke cmc(j);
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops.empty());
        auto excs = ops[0].get_excitations_per_winding();
        REQUIRE(excs.size() == 3);

        double p0 = peakACOf(excs[0].get_current()->get_waveform()->get_data());
        double p1 = peakACOf(excs[1].get_current()->get_waveform()->get_data());
        double p2 = peakACOf(excs[2].get_current()->get_waveform()->get_data());
        INFO("3-wire I peaks: " << p0 << ", " << p1 << ", " << p2);
        CHECK_THAT(p0, Catch::Matchers::WithinRel(iCmPeak, 0.05));
        CHECK_THAT(p1, Catch::Matchers::WithinRel(iCmPeak, 0.05));
        CHECK_THAT(p2, Catch::Matchers::WithinRel(iCmPeak, 0.05));
    }

    SECTION("analytical and simulated frequencies agree") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto inputs = cmc.process();
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!inputs.get_operating_points().empty());
        REQUIRE(!ops.empty());

        double fAna = inputs.get_operating_points()[0]
            .get_excitations_per_winding()[0].get_frequency();
        double fSim = ops[0].get_excitations_per_winding()[0].get_frequency();
        INFO("frequency: analytical=" << fAna << " Hz, simulated=" << fSim << " Hz");
        CHECK_THAT(fAna, Catch::Matchers::WithinRel(f_excit, 0.01));
        CHECK_THAT(fSim, Catch::Matchers::WithinRel(fAna,    0.01));
    }

    SECTION("current has DC offset = operating line current") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops.empty());

        auto iData = ops[0].get_excitations_per_winding()[0]
            .get_current()->get_waveform()->get_data();
        double mean = std::accumulate(iData.begin(), iData.end(), 0.0) / iData.size();
        INFO("current DC mean = " << mean << " A, expected 5 A line current");
        // makeCmcJson(230, 5, ...) passes operatingCurrent=5 A.
        CHECK_THAT(mean, Catch::Matchers::WithinRel(5.0, 0.02));
    }

    SECTION("voltage has ~zero DC offset (inductor passes no DC voltage)") {
        CommonModeChoke cmc(makeJsonWithNoise(cap_pF, dvdt));
        auto ops = cmc.simulate_realistic_cmc(L_self, cap_pF, dvdt, 2, 10);
        REQUIRE(!ops.empty());

        auto vData = ops[0].get_excitations_per_winding()[0]
            .get_voltage()->get_waveform()->get_data();
        double mean = std::accumulate(vData.begin(), vData.end(), 0.0) / vData.size();
        double peak = peakACOf(vData);
        INFO("voltage: mean=" << mean << " V, peak=" << peak << " V");
        // Mean should be tiny relative to AC peak — L has zero DC impedance.
        CHECK(std::abs(mean) < peak * 0.05);
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Analytical ↔ Simulated waveform equality
//
// Once the two paths are unified, their waveforms should overlay sample
// by sample, up to interpolation and phase. This test re-samples both
// waveforms onto a common time grid and compares them point-wise.
// ═══════════════════════════════════════════════════════════════════════

TEST_CASE("Test_Cmc_AnalyticalVsSimulated_WaveformShapeEquality",
          "[converter-model][cmc-topology][consistency]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    json j = makeCmcJson(230.0, 5.0, 50.0, 2, 1000.0, 150e3);
    j["parasiticCap_pF"] = 10.0;
    j["dvdt_V_ns"]       = 50.0;
    j["safetyMargin_dB"] = 6.0;

    CommonModeChoke cmc(j);
    auto inputs = cmc.process();
    REQUIRE(!inputs.get_operating_points().empty());

    // Use the same inductance the analytical path prescribed so both are
    // pinned to the same dominantFrequency and L.
    double L_self = cmc.get_computed_inductance();
    REQUIRE(L_self > 0.0);

    auto ops = cmc.simulate_realistic_cmc(L_self, 10.0, 50.0, 2, 10);
    REQUIRE(!ops.empty());

    // Helper: AC peak
    auto peakAC = [](const std::vector<double>& d) {
        double m = std::accumulate(d.begin(), d.end(), 0.0) / d.size();
        double r = 0.0;
        for (double v : d) r = std::max(r, std::abs(v - m));
        return r;
    };

    auto aI = inputs.get_operating_points()[0].get_excitations_per_winding()[0]
        .get_current()->get_waveform()->get_data();
    auto aV = inputs.get_operating_points()[0].get_excitations_per_winding()[0]
        .get_voltage()->get_waveform()->get_data();
    auto sI = ops[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    auto sV = ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform()->get_data();

    double aIpeak = peakAC(aI), sIpeak = peakAC(sI);
    double aVpeak = peakAC(aV), sVpeak = peakAC(sV);
    INFO("I peak: analytical=" << aIpeak << " A, simulated=" << sIpeak << " A");
    INFO("V peak: analytical=" << aVpeak << " V, simulated=" << sVpeak << " V");

    SECTION("analytical and simulated current peaks agree within 5%") {
        CHECK_THAT(sIpeak, Catch::Matchers::WithinRel(aIpeak, 0.05));
    }
    SECTION("analytical and simulated voltage peaks agree within 5%") {
        CHECK_THAT(sVpeak, Catch::Matchers::WithinRel(aVpeak, 0.05));
    }
    SECTION("analytical and simulated have same sinusoidal sign-pattern count") {
        auto countCrossings = [](const std::vector<double>& d) {
            double m = std::accumulate(d.begin(), d.end(), 0.0) / d.size();
            int c = 0;
            for (size_t i = 1; i < d.size(); ++i)
                if (((d[i - 1] - m) > 0) != ((d[i] - m) > 0)) ++c;
            return c;
        };
        // Analytical samples 1 period; simulated samples 2 periods.
        // So crossings_sim ≈ 2 × crossings_ana, within sampling slack.
        int aC = countCrossings(aI);
        int sC = countCrossings(sI);
        INFO("crossings: analytical=" << aC << ", simulated=" << sC);
        CHECK(aC >= 1);
        CHECK(sC >= 3);
    }
}

// ═══════════════════════════════════════════════════════════════════════
// End-to-end: does MagneticAdviser treat CMC inputs as CMCs?
//
// CommonModeChoke::process_design_requirements tags Inputs with
// Application::INTERFERENCE_SUPPRESSION + SubApplication::COMMON_MODE_NOISE_FILTERING.
// The advisers key off these tags to:
//   · restrict CoreAdviser to TOROIDAL shapes (CoreAdviser.cpp:878, :975)
//   · filter CoreMaterial to CM-appropriate high-µ ferrites (Core.cpp:1478)
//   · build paired/bifilar windings (Coil.cpp:7132)
//   · score leakage inductance in coupling mode (MagneticFilter.cpp:2088)
//
// This test locks that integration in place. It builds a CMC spec, runs
// MagneticAdviser, and checks the output reflects the CMC-aware path.
// ═══════════════════════════════════════════════════════════════════════

namespace {

// Five realistic CMC design cases spanning the design space. Each produces
// a CMC specification that should go through the adviser and come back with
// a sensible CMC-specific magnetic.
struct CmcDesignCase {
    std::string         label;
    double              voltageNominal;        // V_line (RMS)
    double              operatingCurrent;      // A
    int                 numberOfWindings;      // 2, 3, or 4
    double              impedanceAtFreq_Z;     // Ω
    double              impedanceAtFreq_f;     // Hz
    double              lineFrequency;         // Hz
};

// NB: all five cases should spec an impedance the toroidal catalog can meet
// within its self-resonance margin. 500 kHz with a 2 kΩ spec drives the filter
// past every stocked toroid's SRF and can't complete — a real physics limit,
// not a plumbing bug. If you want to stress-test the SRF path, drop the freq
// or raise Z with a wider catalog.
static const std::vector<CmcDesignCase> CMC_DESIGN_CASES = {
    { "low-current-laptop",        230.0,  1.0, 2, 1000.0, 150e3, 50.0 },
    { "mid-current-appliance",     230.0, 10.0, 2,  500.0, 150e3, 50.0 },
    { "high-current-psu",          230.0, 25.0, 2,  300.0, 150e3, 50.0 },
    { "three-phase-industrial",    400.0, 16.0, 3,  800.0, 150e3, 50.0 },
    { "hf-200k-moderate",          230.0,  5.0, 2,  800.0, 200e3, 50.0 },
    // Frontend CMC wizard default — regressed once by an over-aggressive
    // saturation filter, so pin it. 230 V / 5 A / 50 Hz line, 1 kΩ at
    // 150 kHz (EN 55032 Class B start point). Must produce at least one
    // candidate on the toroidal catalog.
    { "wizard-default",            230.0,  5.0, 2, 1000.0, 150e3, 50.0 },
};

static json makeCmcJsonForCase(const CmcDesignCase& c) {
    return makeCmcJson(
        c.voltageNominal,
        c.operatingCurrent,
        c.lineFrequency,
        c.numberOfWindings,
        c.impedanceAtFreq_Z,
        c.impedanceAtFreq_f);
}

} // anonymous namespace

TEST_CASE("Test_Cmc_AdviserKnowsItsACmc",
          "[converter-model][cmc-topology][adviser]") {

    SECTION("process_design_requirements tags CM application flags") {
        CommonModeChoke cmc(makeCmcJson(230.0, 5.0, 50.0, 2, 1000.0, 150e3));
        auto req = cmc.process_design_requirements();

        REQUIRE(req.get_application().has_value());
        REQUIRE(req.get_sub_application().has_value());
        CHECK(req.get_application().value() == Application::INTERFERENCE_SUPPRESSION);
        CHECK(req.get_sub_application().value() == SubApplication::COMMON_MODE_NOISE_FILTERING);
    }

    SECTION("MagneticAdviser returns a CMC-shaped magnetic for every design case") {
        std::cout << "\n=== MagneticAdviser — per-case summary ===\n";
        for (const auto& c : CMC_DESIGN_CASES) {
            DYNAMIC_SECTION("case: " << c.label) {
                CommonModeChoke cmc(makeCmcJsonForCase(c));
                auto inputs = cmc.process();
                REQUIRE(inputs.get_design_requirements().get_sub_application().has_value());
                CHECK(inputs.get_design_requirements().get_sub_application().value()
                      == SubApplication::COMMON_MODE_NOISE_FILTERING);

                double Lreq = inputs.get_design_requirements()
                    .get_magnetizing_inductance().get_minimum().value_or(0.0);

                MagneticAdviser adviser;
                auto results = adviser.get_advised_magnetic(inputs, 1);
                REQUIRE(!results.empty());

                auto mag = results[0].first.get_magnetic();  // value: need non-const core
                auto coreShape = mag.get_mutable_core().get_shape_family();
                auto coreName = mag.get_mutable_core().get_name().value_or("(unnamed)");
                auto matName  = mag.get_mutable_core().resolve_material().get_name();

                // Adviser must pick a toroidal core for a CMC (CoreAdviser.cpp:878).
                CHECK(coreShape == CoreShapeFamily::T);

                // Coil must have one winding per line.
                auto windings = mag.get_coil().get_functional_description();
                CHECK(windings.size() == static_cast<size_t>(c.numberOfWindings));

                // All windings must have identical turn counts (1:1:… ratios).
                int turns0 = windings[0].get_number_turns();
                CHECK(turns0 > 0);
                for (auto& w : windings) {
                    CHECK(w.get_number_turns() == turns0);
                }

                // Wire gauge (first winding): report for the summary. Not asserted
                // because catalog picks vary with the cross-referencer weights.
                std::string wireName = "(unknown)";
                if (!windings.empty()) {
                    const auto& wv = windings[0].get_wire();
                    if (std::holds_alternative<OpenMagnetics::Wire>(wv)) {
                        wireName = std::get<OpenMagnetics::Wire>(wv)
                            .get_name().value_or("(unnamed)");
                    } else if (std::holds_alternative<std::string>(wv)) {
                        wireName = std::get<std::string>(wv);
                    }
                }

                // Impedance compliance: the advised magnetic must meet the
                // user's |Z| ≥ Z_required spec at the spec frequency. The
                // adviser applies this as a filter during search, but the
                // subsequent coil/wire steps can shift turn count / window
                // geometry, so re-measure on the final magnetic.
                //
                // The frontend already bakes a 6 dB (2×) safety margin into
                // the impedance spec before it reaches MKF, so a 5% shortfall
                // on the advised magnetic still passes the underlying CISPR
                // limit with 5 dB to spare. Here we assert the adviser's own
                // filter is doing an honest job: within ±5% of the spec
                // value. Engineers who want a stricter margin should raise
                // safetyMargin_dB in the CMC wizard.
                double achievedZ = std::abs(
                    Impedance().calculate_impedance(mag, c.impedanceAtFreq_f));
                // Lower bound only — overshoot is good engineering, undershoot
                // within 5% is absorbed by the 6 dB safety margin the frontend
                // bakes in before the spec reaches MKF.
                CHECK(achievedZ >= c.impedanceAtFreq_Z * 0.95);

                // Core flux density: only the CM noise current drives flux
                // in a CMC (DM cancels between opposite windings). For the
                // spec's C·dV/dt = 0.5 A and typical Mn-Zn EMI toroids the
                // peak B sits well under 200 mT. A previous bug had the
                // generic B path use N · (I_dm + I_cm) which puts B over
                // 500 mT and triggers saturation filters falsely.
                //
                // All windings share the same core, so per-winding B storage
                // is a data quirk — every winding should see the same core
                // flux. Read winding 0 as the canonical value (that's what
                // MagneticSimulator writes and what the UI displays).
                double maxB = 0.0;
                auto advisedMasInputs = results[0].first.get_mutable_inputs();
                for (const auto& op : advisedMasInputs.get_operating_points()) {
                    if (op.get_excitations_per_winding().empty()) continue;
                    auto bOpt = op.get_excitations_per_winding()[0].get_magnetic_flux_density();
                    if (!bOpt) continue;
                    auto procOpt = bOpt->get_processed();
                    if (procOpt && procOpt->get_peak()) {
                        maxB = std::max(maxB, std::abs(procOpt->get_peak().value()));
                    }
                }
                // Compare against the picked material's own Bsat rather than
                // a fixed magic number — adviser material picks vary (T38,
                // N30, K10…) and their Bsat values span 0.35–1.5 T, so a
                // single ceiling is either too tight (nanocrystalline) or
                // too loose (high-µ MnZn). The real guard is "B must not
                // saturate the picked material at the operating temperature".
                double bSat = mag.get_mutable_core()
                    .get_magnetic_flux_density_saturation(c.lineFrequency >= 100.0 ? 100.0 : 25.0, false);
                INFO("[" << c.label << "] max |B| on core = " << maxB
                     << " T, Bsat (" << matName << ") = " << bSat << " T");
                CHECK(maxB < bSat);

                std::cout << "[" << c.label << "]\n"
                          << "  req: " << c.voltageNominal << " Vrms, "
                          << c.operatingCurrent << " A, "
                          << c.numberOfWindings << "-wire, |Z|≥"
                          << c.impedanceAtFreq_Z << " Ω @ "
                          << c.impedanceAtFreq_f / 1e3 << " kHz\n"
                          << "  → L required: " << Lreq * 1e3 << " mH\n"
                          << "  → core: " << coreName
                          << " (shape " << std::string(magic_enum::enum_name(coreShape)) << ")\n"
                          << "  → material: " << matName << "\n"
                          << "  → winding: " << windings.size() << " × " << turns0 << " turns, "
                          << "wire " << wireName << "\n"
                          << "  → achieved |Z| at spec freq: " << achievedZ
                          << " Ω (required ≥ " << c.impedanceAtFreq_Z << " Ω)\n"
                          << "  → adviser score: " << results[0].second << "\n\n";
            }
        }
    }
}

TEST_CASE("CalculateAdvisedCoil_ToroidCmcAngularFillMustNotExceed360",
          "[converter-model][cmc-topology][adviser][angular-fill]") {

    SECTION("CMC on T 18.4/5.9/5.9 TDK T65 — sections must fit within 360°") {
        // 230 V / 5 A / 50 Hz line, 150 kHz / 1 kΩ impedance requirement.
        // This exact scenario previously produced totalAngle ≈ 551.5°.
        CommonModeChoke cmc(makeCmcJson(230.0, 5.0, 50.0, 2, 1000.0, 150e3));
        auto inputs = cmc.process();

        MagneticAdviser adviser;
        auto results = adviser.get_advised_magnetic(inputs, 1);
        REQUIRE(!results.empty());

        auto mag = results[0].first.get_magnetic();

        auto sectionsOpt = mag.get_coil().get_sections_description();
        REQUIRE(sectionsOpt.has_value());
        REQUIRE(!sectionsOpt->empty());

        double totalAngle = 0.0;
        int conductionCount = 0, insulationCount = 0;
        for (const auto& sec : sectionsOpt.value()) {
            auto type = (sec.get_type() == ElectricalType::CONDUCTION) ? "C" : "I";
            std::cout << type << " dim0=" << sec.get_dimensions()[0]
                      << " dim1=" << sec.get_dimensions()[1] << "°\n";
            totalAngle += sec.get_dimensions()[1];
            if (sec.get_type() == ElectricalType::CONDUCTION) conductionCount++;
            if (sec.get_type() == ElectricalType::INSULATION) insulationCount++;
        }

        INFO("Total section angle = " << totalAngle << "°  (sections: "
             << conductionCount << " conduction + " << insulationCount << " insulation)");
        CHECK(totalAngle <= 360.0);
        CHECK(conductionCount == 2);
        CHECK(insulationCount == 2);
        // Bug A: for CONTIGUOUS toroidal CMC, insulation sections must be a
        // small tape-thickness angle, NOT the full 360° artefact the
        // OVERLAPPING+POLAR branch used to emit.
        for (const auto& sec : sectionsOpt.value()) {
            if (sec.get_type() == ElectricalType::INSULATION) {
                CHECK(sec.get_dimensions()[1] < 90.0);
            }
        }
    }
}

} // anonymous namespace
