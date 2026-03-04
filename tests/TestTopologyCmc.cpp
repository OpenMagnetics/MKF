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
    // V_peak = 230 × √2 ≈ 325.3 V
    [[maybe_unused]] const double vPeak = 230.0 * std::sqrt(2.0);
    // Switching frequency: dV/dt [V/s] / V_operating = 50e9 / 230 ~ 217 MHz
    //   → capped at 1 MHz
    const double switchingFrequency = 1e6;
    // CM noise current: C_parasitic × dV/dt = 10e-12 × 50e9 = 0.5 A
    // Load: R_load = V_operating / I_operating = 230 / 5 = 46 Ω
    // Sim time: 5 / f_sw = 5 µs
    const double expectedSimTime    = 5.0 / switchingFrequency;

    SECTION("Simulation produces one operating point with 2 winding excitations") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);

        REQUIRE(ops.size() == 1);
        CHECK(ops[0].get_excitations_per_winding().size() == 2);
    }

    SECTION("Operating point name contains switching frequency") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        auto name = ops[0].get_name().value_or("");
        // Should contain "1000kHz" (1 MHz)
        CHECK(name.find("1000kHz") != std::string::npos);
    }

    SECTION("Excitation frequency equals the switching frequency (not line frequency)") {
        auto ops = cmc.simulate_realistic_cmc(inductance, parasiticCap_pF, dvdt_V_ns);
        REQUIRE(!ops.empty());

        for (const auto& exc : ops[0].get_excitations_per_winding()) {
            CHECK_THAT(exc.get_frequency(), Catch::Matchers::WithinRel(switchingFrequency, 0.01));
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

} // anonymous namespace
