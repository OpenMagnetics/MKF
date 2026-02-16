#include <source_location>
#include "support/Painter.h"
#include "converter_models/Llc.h"
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
#include <chrono>
#include <cmath>

using namespace OpenMagnetics;

namespace {
    auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    // ═══════════════════════════════════════════════════════════════════════
    // Basic Design Tests
    // ═══════════════════════════════════════════════════════════════════════

    TEST_CASE("Test_Llc_HalfBridge_Design", "[converter-model][llc-topology][smoke-test]") {
        json llcJson;
        json inputVoltage;

        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 370.0;
        inputVoltage["maximum"] = 410.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = false;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);

        SECTION("Input validation - valid configuration") {
            CHECK(llc.run_checks(false) == true);
        }

        SECTION("Turns ratio calculation") {
            auto req = llc.process_design_requirements();
            // Formula: n = (Vin_min * k_bridge) / Vout
            // For half-bridge with Vin_min=370V, k_bridge=0.5, Vout=12V:
            // n = (370 * 0.5) / 12 = 15.42
            double expectedN = 370.0 * 0.5 / 12.0;
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
        }

        SECTION("Magnetizing inductance is positive") {
            auto req = llc.process_design_requirements();
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            CHECK(Lm > 0);
            CHECK(Lm > 10e-6);   // At least 10 µH
            CHECK(Lm < 100e-3);  // Less than 100 mH
        }

        SECTION("Resonant tank design") {
            llc.process_design_requirements();
            double Lr = llc.get_computed_resonant_inductance();
            double Cr = llc.get_computed_resonant_capacitance();

            CHECK(Lr > 0);
            CHECK(Cr > 0);

            // Verify resonant frequency: fr = 1 / (2*pi*sqrt(Lr*Cr))
            double fr = 1.0 / (2.0 * M_PI * std::sqrt(Lr * Cr));
            double expectedFr = llc.get_effective_resonant_frequency();
            REQUIRE_THAT(fr, Catch::Matchers::WithinAbs(expectedFr, expectedFr * 0.01));
        }

        SECTION("Inductance ratio") {
            llc.set_computed_inductance_ratio(7.0);
            auto req = llc.process_design_requirements();
            double Lr = llc.get_computed_resonant_inductance();
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            REQUIRE_THAT(Lm / Lr, Catch::Matchers::WithinAbs(7.0, 0.01));
        }

        SECTION("Integrated resonant inductor") {
            llc.set_integrated_resonant_inductor(true);
            auto req = llc.process_design_requirements();
            
            REQUIRE(req.get_leakage_inductance().has_value());
            double leakageTarget = req.get_leakage_inductance().value()[0].get_nominal().value();
            double Lr = llc.get_computed_resonant_inductance();
            REQUIRE_THAT(leakageTarget, Catch::Matchers::WithinAbs(Lr, Lr * 0.01));
        }
    }

    TEST_CASE("Test_Llc_Analytical_Waveforms_NonZero", "[converter-model][llc-topology][analytical-waveforms]") {
        json llcJson;
        json inputVoltage;

        // Use default values from WebFrontend
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 370.0;
        inputVoltage["maximum"] = 410.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = false;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);
        auto req = llc.process_design_requirements();

        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        auto ops = llc.process_operating_points(turnsRatios, Lm);
        REQUIRE(!ops.empty());
        REQUIRE(ops[0].get_excitations_per_winding().size() >= 2);

        // Check primary waveforms
        auto& primaryExc = ops[0].get_excitations_per_winding()[0];
        REQUIRE(primaryExc.get_current().has_value());
        REQUIRE(primaryExc.get_voltage().has_value());

        auto priCurrent = primaryExc.get_current()->get_waveform().value().get_data();
        auto priVoltage = primaryExc.get_voltage()->get_waveform().value().get_data();

        CHECK(!priCurrent.empty());
        CHECK(!priVoltage.empty());

        // Check secondary waveforms
        auto& secondaryExc = ops[0].get_excitations_per_winding()[1];
        REQUIRE(secondaryExc.get_current().has_value());
        REQUIRE(secondaryExc.get_voltage().has_value());

        auto secCurrent = secondaryExc.get_current()->get_waveform().value().get_data();
        auto secVoltage = secondaryExc.get_voltage()->get_waveform().value().get_data();

        CHECK(!secCurrent.empty());
        CHECK(!secVoltage.empty());

        // CRITICAL: Secondary current should NOT be all zeros
        double secCurrentMax = *std::max_element(secCurrent.begin(), secCurrent.end());
        double secCurrentMin = *std::min_element(secCurrent.begin(), secCurrent.end());

        std::cout << "Primary current max: " << *std::max_element(priCurrent.begin(), priCurrent.end()) << std::endl;
        std::cout << "Secondary current max: " << secCurrentMax << std::endl;
        std::cout << "Secondary current min: " << secCurrentMin << std::endl;
        std::cout << "Turns ratio: " << turnsRatios[0] << std::endl;

        // Debug: Print first few values of primary and secondary currents
        std::cout << "DEBUG: Primary current first 5 values: ";
        for (size_t i = 0; i < std::min(size_t(5), priCurrent.size()); ++i) {
            std::cout << priCurrent[i] << " ";
        }
        std::cout << std::endl;
        
        std::cout << "DEBUG: Secondary current first 5 values: ";
        for (size_t i = 0; i < std::min(size_t(5), secCurrent.size()); ++i) {
            std::cout << secCurrent[i] << " ";
        }
        std::cout << std::endl;

        // Secondary current should have some non-zero values
        // The max should be significantly greater than 0 (at least 0.1A for a 10A output)
        CHECK(secCurrentMax > 0.5);  // Should be around 10A/8.33 = 1.2A at least

        // Secondary voltage should also not be all zeros
        double secVoltageMax = *std::max_element(secVoltage.begin(), secVoltage.end());
        CHECK(secVoltageMax > 1.0);  // Should be around 12V
    }

    TEST_CASE("Test_Llc_FullBridge_Design", "[converter-model][llc-topology][smoke-test]") {
        json llcJson;
        json inputVoltage;

        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 370.0;
        inputVoltage["maximum"] = 410.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Full Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.35;
        llcJson["integratedResonantInductor"] = false;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {48.0};
            LlcOperatingPointJson["outputCurrents"] = {5.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);

        SECTION("Turns ratio for full-bridge") {
            auto req = llc.process_design_requirements();
            // n = Vin * 1.0 / (2 * Vout) for full-bridge
            double expectedN = 400.0 * 1.0 / (2.0 * 48.0);
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.02));
        }

        SECTION("Bridge voltage factor") {
            double factor = llc.get_bridge_voltage_factor();
            REQUIRE_THAT(factor, Catch::Matchers::WithinAbs(1.0, 0.001));
        }

        SECTION("Full-bridge primary voltage amplitude") {
            auto req = llc.process_design_requirements();
            
            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());
            
            // Copy waveform to avoid dangling reference
            auto priVwfm = ops[0].get_excitations_per_winding()[0]
                                .get_voltage()->get_waveform().value();
            double Vp = priVwfm.get_data()[1];  // positive peak
            REQUIRE_THAT(std::abs(Vp), Catch::Matchers::WithinAbs(400.0, 1.0));  // Vin_nom = 400 V
        }

        SECTION("Waveform plotting - Full Bridge") {
            auto req = llc.process_design_requirements();

            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto ops = llc.process_operating_points(turnsRatios, Lm);

            // Plot primary current waveform
            auto outFile = outputFilePath;
            outFile.append("Test_Llc_FullBridge_Primary_Current_Waveform.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();

            // Plot primary voltage waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_FullBridge_Primary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
                painter.export_svg();
            }
        }
    }

    TEST_CASE("Test_Llc_OperatingPoints_Generation", "[converter-model][llc-topology][smoke-test]") {
        json llcJson;
        json inputVoltage;

        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 370.0;
        inputVoltage["maximum"] = 410.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);
        auto req = llc.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        SECTION("Multiple input voltages") {
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            // Should have 3 OPs: nominal, minimum, maximum input voltage
            CHECK(ops.size() == 3);
        }

        SECTION("Waveform structure") {
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());

            auto& op = ops[0];
            // 1 primary + 1 secondary = 2 excitations per winding
            CHECK(op.get_excitations_per_winding().size() == 2);

            auto& priExc = op.get_excitations_per_winding()[0];
            REQUIRE(priExc.get_current().has_value());
            REQUIRE(priExc.get_voltage().has_value());
            CHECK(priExc.get_frequency() == 100e3);

            // Copy waveforms to avoid dangling references
            auto currentWfm = priExc.get_current()->get_waveform().value();
            CHECK(currentWfm.get_data().size() > 10);

            auto voltageWfm = priExc.get_voltage()->get_waveform().value();
            CHECK(voltageWfm.get_data().size() == 6);
        }

        SECTION("Primary voltage symmetry") {
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            auto& nomOp = ops[0];
            
            // Copy waveform to avoid dangling reference
            auto priVwfm = nomOp.get_excitations_per_winding()[0]
                                .get_voltage()->get_waveform().value();
            auto vData = priVwfm.get_data();

            // Bipolar rectangular should be symmetric around zero
            // vData[1] = +Vp, vData[3] = -Vp
            REQUIRE_THAT(vData[1], Catch::Matchers::WithinAbs(-vData[3], 1e-6));
        }

        SECTION("Secondary current non-negative") {
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            auto& op = ops[0];
            auto& secExc = op.get_excitations_per_winding()[1];
            
            REQUIRE(secExc.get_current().has_value());
            
            // Copy waveform and data to avoid dangling references
            auto secIWfm = secExc.get_current()->get_waveform().value();
            auto secIData = secIWfm.get_data();
            
            for (double val : secIData) {
                CHECK(val >= -1e-10);  // Allow tiny numerical noise
            }
        }

        SECTION("Waveform plotting - OP Generation") {
            auto ops = llc.process_operating_points(turnsRatios, Lm);

            // Plot primary current waveform
            auto outFile = outputFilePath;
            outFile.append("Test_Llc_OP_Generation_Primary_Current_Waveform.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();

            // Plot primary voltage waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_OP_Generation_Primary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
                painter.export_svg();
            }
        }
    }

    TEST_CASE("Test_Llc_Operating_Modes", "[converter-model][llc-topology][smoke-test]") {
        SECTION("Below resonance (boost mode)") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["resonantFrequency"] = 120000;
            llcJson["operatingPoints"] = json::array();

            {
                json LlcOperatingPointJson;
                LlcOperatingPointJson["ambientTemperature"] = 25.0;
                LlcOperatingPointJson["outputVoltages"] = {12.0};
                LlcOperatingPointJson["outputCurrents"] = {10.0};
                LlcOperatingPointJson["switchingFrequency"] = 90000;
                llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
            }

            OpenMagnetics::Llc llc(llcJson);
            auto req = llc.process_design_requirements();
            
            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());
            
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_current().has_value());
            REQUIRE(priExc.get_voltage().has_value());
        }

        SECTION("Waveform plotting - Below Resonance") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["resonantFrequency"] = 120000;
            llcJson["operatingPoints"] = json::array();

            {
                json LlcOperatingPointJson;
                LlcOperatingPointJson["ambientTemperature"] = 25.0;
                LlcOperatingPointJson["outputVoltages"] = {12.0};
                LlcOperatingPointJson["outputCurrents"] = {10.0};
                LlcOperatingPointJson["switchingFrequency"] = 90000;
                llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
            }

            OpenMagnetics::Llc llc(llcJson);
            auto req = llc.process_design_requirements();

            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto ops = llc.process_operating_points(turnsRatios, Lm);

            // Plot primary current waveform
            auto outFile = outputFilePath;
            outFile.append("Test_Llc_Below_Resonance_Primary_Current_Waveform.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();

            // Plot primary voltage waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_Below_Resonance_Primary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
                painter.export_svg();
            }
        }

        SECTION("Above resonance (buck mode)") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["resonantFrequency"] = 80000;
            llcJson["operatingPoints"] = json::array();

            {
                json LlcOperatingPointJson;
                LlcOperatingPointJson["ambientTemperature"] = 25.0;
                LlcOperatingPointJson["outputVoltages"] = {12.0};
                LlcOperatingPointJson["outputCurrents"] = {10.0};
                LlcOperatingPointJson["switchingFrequency"] = 150000;
                llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
            }

            OpenMagnetics::Llc llc(llcJson);
            auto req = llc.process_design_requirements();
            
            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());
            
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_current().has_value());
        }

        SECTION("Waveform plotting - Above Resonance") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["resonantFrequency"] = 80000;
            llcJson["operatingPoints"] = json::array();

            {
                json LlcOperatingPointJson;
                LlcOperatingPointJson["ambientTemperature"] = 25.0;
                LlcOperatingPointJson["outputVoltages"] = {12.0};
                LlcOperatingPointJson["outputCurrents"] = {10.0};
                LlcOperatingPointJson["switchingFrequency"] = 150000;
                llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
            }

            OpenMagnetics::Llc llc(llcJson);
            auto req = llc.process_design_requirements();

            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto ops = llc.process_operating_points(turnsRatios, Lm);

            // Plot primary current waveform
            auto outFile = outputFilePath;
            outFile.append("Test_Llc_Above_Resonance_Primary_Current_Waveform.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();

            // Plot primary voltage waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_Above_Resonance_Primary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
                painter.export_svg();
            }
        }
    }

    TEST_CASE("Test_Llc_Multiple_Outputs", "[converter-model][llc-topology][smoke-test]") {
        json llcJson;
        json inputVoltage;
        inputVoltage["nominal"] = 400.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0, 5.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0, 3.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);
        auto req = llc.process_design_requirements();

        SECTION("Turns ratios for each output") {
            CHECK(req.get_turns_ratios().size() == 2);
        }

        SECTION("Operating points with multiple secondaries") {
            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            CHECK(ops.size() == 1);  // 1 vin (nominal only) * 1 opPoint
            
            // 1 primary + 2 secondaries = 3 excitations
            CHECK(ops[0].get_excitations_per_winding().size() == 3);
        }

        SECTION("Waveform plotting - Multiple Outputs") {
            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto ops = llc.process_operating_points(turnsRatios, Lm);

            // Plot primary current waveform
            auto outFile = outputFilePath;
            outFile.append("Test_Llc_Multiple_Outputs_Primary_Current_Waveform.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();

            // Plot primary voltage waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_Multiple_Outputs_Primary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
                painter.export_svg();
            }
        }
    }

    TEST_CASE("Test_Llc_Multiple_OperatingPoints", "[converter-model][llc-topology][smoke-test]") {
        json llcJson;
        json inputVoltage;
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 370.0;
        inputVoltage["maximum"] = 410.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["operatingPoints"] = json::array();

        {
            json op1;
            op1["ambientTemperature"] = 25.0;
            op1["outputVoltages"] = {12.0};
            op1["outputCurrents"] = {10.0};
            op1["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(op1);

            json op2;
            op2["ambientTemperature"] = 50.0;
            op2["outputVoltages"] = {12.0};
            op2["outputCurrents"] = {5.0};
            op2["switchingFrequency"] = 130000;
            llcJson["operatingPoints"].push_back(op2);
        }

        OpenMagnetics::Llc llc(llcJson);
        auto req = llc.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
        
        auto ops = llc.process_operating_points(turnsRatios, Lm);
        // 3 vin (nom, min, max) * 2 opPoints = 6
        CHECK(ops.size() == 6);

        SECTION("Waveform plotting - Multiple OP") {
            // Plot primary current waveform for the first operating point
            auto outFile = outputFilePath;
            outFile.append("Test_Llc_Multiple_OP_Primary_Current_Waveform.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();

            // Plot primary voltage waveform for the first operating point
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_Multiple_OP_Primary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
                painter.export_svg();
            }
        }
    }

    TEST_CASE("Test_AdvancedLlc_User_Defined", "[converter-model][llc-topology][advanced][smoke-test]") {
        json llcJson;
        json inputVoltage;
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 370.0;
        inputVoltage["maximum"] = 410.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        llcJson["desiredTurnsRatios"] = {8.33};
        llcJson["desiredMagnetizingInductance"] = 500e-6;
        llcJson["desiredResonantInductance"] = 100e-6;

        OpenMagnetics::AdvancedLlc llc(llcJson);
        
        // Use explicit namespace to avoid ambiguity with MAS::Inputs
        OpenMagnetics::Inputs inputs = llc.process();

        SECTION("Design requirements match user values") {
            CHECK(inputs.get_design_requirements().get_turns_ratios().size() == 1);
            double Lm = inputs.get_design_requirements().get_magnetizing_inductance()
                            .get_nominal().value();
            REQUIRE_THAT(Lm, Catch::Matchers::WithinAbs(500e-6, 1e-9));
        }

        SECTION("Leakage inductance request") {
            REQUIRE(inputs.get_design_requirements().get_leakage_inductance().has_value());
            double Lr = inputs.get_design_requirements().get_leakage_inductance()
                            .value()[0].get_nominal().value();
            REQUIRE_THAT(Lr, Catch::Matchers::WithinAbs(100e-6, 1e-9));
        }

        SECTION("Operating points generated") {
            // 3 Vin * 1 OP = 3 operating points
            CHECK(inputs.get_operating_points().size() == 3);
        }

        SECTION("Waveform plotting - Advanced LLC") {
            // Plot primary current waveform for the first operating point
            {
                auto outFile = outputFilePath;
                outFile.append("Test_AdvancedLlc_Primary_Current_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
                painter.export_svg();
            }

            // Plot primary voltage waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_AdvancedLlc_Primary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
                painter.export_svg();
            }
        }
    }

    TEST_CASE("Test_Llc_EndToEnd_Process", "[converter-model][llc-topology][integration][smoke-test]") {
        json llcJson;
        json inputVoltage;
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 370.0;
        inputVoltage["maximum"] = 410.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);
        llc._assertErrors = true;

        // Use explicit namespace
        OpenMagnetics::Inputs inputs = llc.process();

        SECTION("Design requirements populated") {
            CHECK(inputs.get_operating_points().size() >= 1);
            CHECK(inputs.get_design_requirements().get_turns_ratios().size() == 1);
            
            double Lm = inputs.get_design_requirements().get_magnetizing_inductance()
                            .get_nominal().value();
            CHECK(Lm > 0);
        }

        SECTION("Waveform plotting - End to End") {
            // Plot primary current waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_EndToEnd_Primary_Current_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
                painter.export_svg();
            }

            // Plot secondary current waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_EndToEnd_Secondary_Current_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
                painter.export_svg();
            }

            // Plot primary voltage waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_EndToEnd_Primary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
                painter.export_svg();
            }

            // Plot secondary voltage waveform
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Llc_EndToEnd_Secondary_Voltage_Waveform.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
                painter.export_svg();
            }
        }
    }

    TEST_CASE("Test_Llc_Ngspice_Circuit", "[converter-model][llc-topology][ngspice-simulation]") {
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        json llcJson;
        json inputVoltage;
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 370.0;
        inputVoltage["maximum"] = 410.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);
        auto req = llc.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        SECTION("Half-bridge netlist contains required elements") {
            std::string netlist = llc.generate_ngspice_circuit(turnsRatios, Lm, 0, 0);

            // Input and switches
            CHECK(netlist.find("Vdc") != std::string::npos);
            CHECK(netlist.find("S_hi") != std::string::npos);
            CHECK(netlist.find("S_lo") != std::string::npos);
            // Resonant tank
            CHECK(netlist.find("Lr") != std::string::npos);  // Series resonant inductor
            CHECK(netlist.find("Cr") != std::string::npos);  // Series resonant capacitor
            // Transformer
            CHECK(netlist.find("Lpri") != std::string::npos);  // Primary magnetizing inductance
            CHECK(netlist.find("Lsec1") != std::string::npos); // Secondary winding 1
            CHECK(netlist.find("Lsec2") != std::string::npos); // Secondary winding 2 (center-tapped)
            // Coupling statements
            CHECK(netlist.find("K1") != std::string::npos);
            CHECK(netlist.find("K2") != std::string::npos);
            CHECK(netlist.find("K3") != std::string::npos);
            // Output rectifier
            CHECK(netlist.find("D1") != std::string::npos);
            CHECK(netlist.find("D2") != std::string::npos);
            // Output filter
            CHECK(netlist.find("Cout") != std::string::npos);
            CHECK(netlist.find("Rload") != std::string::npos);
            // Simulation commands
            CHECK(netlist.find(".tran") != std::string::npos);
            CHECK(netlist.find(".end") != std::string::npos);
            CHECK(netlist.find("Half") != std::string::npos);
        }

        SECTION("Full-bridge netlist") {
            json fbJson;
            json fbInputVoltage;
            fbInputVoltage["nominal"] = 400.0;
            fbJson["inputVoltage"] = fbInputVoltage;
            fbJson["bridgeType"] = "Full Bridge";
            fbJson["minSwitchingFrequency"] = 80000;
            fbJson["maxSwitchingFrequency"] = 200000;
            fbJson["qualityFactor"] = 0.35;
            fbJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                op["outputVoltages"] = {48.0};
                op["outputCurrents"] = {5.0};
                op["switchingFrequency"] = 100000;
                fbJson["operatingPoints"].push_back(op);
            }

            OpenMagnetics::Llc fbLlc(fbJson);
            auto fbReq = fbLlc.process_design_requirements();
            
            std::vector<double> fbTurnsRatios;
            for (auto& tr : fbReq.get_turns_ratios()) {
                fbTurnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double fbLm = resolve_dimensional_values(fbReq.get_magnetizing_inductance());
            
            std::string netlist = fbLlc.generate_ngspice_circuit(fbTurnsRatios, fbLm, 0, 0);

            CHECK(netlist.find("Full") != std::string::npos);
            CHECK(netlist.find("S1") != std::string::npos);
            CHECK(netlist.find("S3") != std::string::npos);
            CHECK(netlist.find("S4") != std::string::npos);
            CHECK(netlist.find("bridge_a") != std::string::npos);
            CHECK(netlist.find("bridge_b") != std::string::npos);
        }
    }

    TEST_CASE("Test_Llc_Input_Validation", "[converter-model][llc-topology][validation][smoke-test]") {
        SECTION("Missing operating points") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["operatingPoints"] = json::array();

            OpenMagnetics::Llc llc(llcJson);
            CHECK(llc.run_checks(false) == false);
        }

        SECTION("Missing input voltage") {
            json llcJson;
            llcJson["inputVoltage"] = json::object();
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["operatingPoints"] = json::array();

            {
                json LlcOperatingPointJson;
                LlcOperatingPointJson["ambientTemperature"] = 25.0;
                LlcOperatingPointJson["outputVoltages"] = {12.0};
                LlcOperatingPointJson["outputCurrents"] = {10.0};
                LlcOperatingPointJson["switchingFrequency"] = 100000;
                llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
            }

            OpenMagnetics::Llc llc(llcJson);
            CHECK(llc.run_checks(false) == false);
        }
    }

    TEST_CASE("Test_Llc_Edge_Cases", "[converter-model][llc-topology][edge-case][smoke-test]") {
        SECTION("High-power full-bridge 1 kW") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 390.0;
            inputVoltage["minimum"] = 340.0;
            inputVoltage["maximum"] = 420.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Full Bridge";
            llcJson["minSwitchingFrequency"] = 60000;
            llcJson["maxSwitchingFrequency"] = 150000;
            llcJson["qualityFactor"] = 0.3;
            llcJson["operatingPoints"] = json::array();

            {
                json LlcOperatingPointJson;
                LlcOperatingPointJson["ambientTemperature"] = 40.0;
                LlcOperatingPointJson["outputVoltages"] = {48.0};
                LlcOperatingPointJson["outputCurrents"] = {20.0};
                LlcOperatingPointJson["switchingFrequency"] = 100000;
                llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
            }

            OpenMagnetics::Llc llc(llcJson);
            auto req = llc.process_design_requirements();
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            CHECK(Lm > 0);

            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }

            auto ops = llc.process_operating_points(turnsRatios, Lm);
            CHECK(ops.size() >= 1);
        }

        SECTION("Low-voltage output 5V USB") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 250000;
            llcJson["operatingPoints"] = json::array();

            {
                json LlcOperatingPointJson;
                LlcOperatingPointJson["ambientTemperature"] = 25.0;
                LlcOperatingPointJson["outputVoltages"] = {5.0};
                LlcOperatingPointJson["outputCurrents"] = {6.0};
                LlcOperatingPointJson["switchingFrequency"] = 100000;
                llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
            }

            OpenMagnetics::Llc llc(llcJson);
            auto req = llc.process_design_requirements();

            // n = 400 * 0.5 / (2 * 5) = 20
            double expectedN = 400.0 * 0.5 / (2.0 * 5.0);
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.03));
        }
    }

    TEST_CASE("Test_Llc_Helper_Functions", "[converter-model][llc-topology][unit][smoke-test]") {
        SECTION("Bridge voltage factor - half-bridge") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                op["outputVoltages"] = {12.0};
                op["outputCurrents"] = {10.0};
                op["switchingFrequency"] = 100000;
                llcJson["operatingPoints"].push_back(op);
            }

            OpenMagnetics::Llc llc(llcJson);
            REQUIRE_THAT(llc.get_bridge_voltage_factor(), Catch::Matchers::WithinAbs(0.5, 0.001));
        }

        SECTION("Effective resonant frequency - user provided") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["resonantFrequency"] = 120000;
            llcJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                op["outputVoltages"] = {12.0};
                op["outputCurrents"] = {10.0};
                op["switchingFrequency"] = 100000;
                llcJson["operatingPoints"].push_back(op);
            }

            OpenMagnetics::Llc llc(llcJson);
            REQUIRE_THAT(llc.get_effective_resonant_frequency(), Catch::Matchers::WithinAbs(120000.0, 0.001));
        }

        SECTION("Effective resonant frequency - default geometric mean") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                op["outputVoltages"] = {12.0};
                op["outputCurrents"] = {10.0};
                op["switchingFrequency"] = 100000;
                llcJson["operatingPoints"].push_back(op);
            }

            OpenMagnetics::Llc llc(llcJson);
            double expected = std::sqrt(80000.0 * 200000.0);
            REQUIRE_THAT(llc.get_effective_resonant_frequency(), Catch::Matchers::WithinAbs(expected, 0.01));
        }

        SECTION("Magnetizing current peak sanity") {
            json llcJson;
            json inputVoltage;
            inputVoltage["nominal"] = 400.0;
            llcJson["inputVoltage"] = inputVoltage;
            llcJson["bridgeType"] = "Half Bridge";
            llcJson["minSwitchingFrequency"] = 80000;
            llcJson["maxSwitchingFrequency"] = 200000;
            llcJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                op["outputVoltages"] = {12.0};
                op["outputCurrents"] = {10.0};
                op["switchingFrequency"] = 100000;
                llcJson["operatingPoints"].push_back(op);
            }

            OpenMagnetics::Llc llc(llcJson);
            auto req = llc.process_design_requirements();

            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            double Vin_nom = 400.0;
            double Vp = Vin_nom * 0.5;  // Half-bridge
            double fsw = 100e3;

            // Expected: Im_pk = Vp / (4 * Lm * fsw)
            double expectedImPk = Vp / (4.0 * Lm * fsw);

            CHECK(expectedImPk > 0);
            CHECK(expectedImPk < 10.0);  // Reasonable for this spec
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Simulation Performance Tests
    // ═══════════════════════════════════════════════════════════════════════

    TEST_CASE("Test_Llc_Simulation_Performance_Default_Inputs", "[converter-model][llc-topology][ngspice-simulation][performance]") {
        // Default inputs from LlcWizard.vue
        json llcJson;
        json inputVoltage;
        
        inputVoltage["nominal"] = 400.0;
        inputVoltage["tolerance"] = 0.1;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Full Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 120000;
        llcJson["resonantFrequency"] = 100000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = true;
        llcJson["magnetizingInductance"] = 200e-6;
        llcJson["turnsRatio"] = 4.0;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {48.0};  // outputVoltage
            LlcOperatingPointJson["outputCurrents"] = {10.4166666667};  // outputPower / outputVoltage = 500 / 48
            LlcOperatingPointJson["switchingFrequency"] = 100000;  // resonantFrequency
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        // Set simulation parameters like the wizard does
        size_t numberOfPeriods = 2;
        size_t numberOfSteadyStatePeriods = 1;

        OpenMagnetics::Llc llc(llcJson);
        auto designRequirements = llc.process_design_requirements();
        
        double magnetizingInductance = 200e-6;
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) {
                turnsRatios.push_back(tr.get_nominal().value());
            }
        }

        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        SECTION("simulate_and_extract_topology_waveforms completes in under 15 seconds") {
            auto start = std::chrono::high_resolution_clock::now();
            
            auto topologyWaveforms = llc.simulate_and_extract_topology_waveforms(
                turnsRatios, magnetizingInductance, numberOfPeriods);
            
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            
            INFO("Topology waveforms simulation took " << duration << " seconds");
            CHECK(duration < 15.0);
            
            REQUIRE(!topologyWaveforms.empty());
        }

        SECTION("simulate_and_extract_operating_points completes in under 15 seconds") {
            auto start = std::chrono::high_resolution_clock::now();
            
            // Set number of periods using the class method
            llc.set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
            auto operatingPoints = llc.simulate_and_extract_operating_points(
                turnsRatios, magnetizingInductance);
            
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            
            INFO("Operating points simulation took " << duration << " seconds");
            CHECK(duration < 15.0);
            
            REQUIRE(!operatingPoints.empty());
        }

        SECTION("Combined simulation (topology + operating points) completes in under 15 seconds") {
            auto start = std::chrono::high_resolution_clock::now();
            
            // Run both simulations like libMKF does
            llc.set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
            auto topologyWaveforms = llc.simulate_and_extract_topology_waveforms(
                turnsRatios, magnetizingInductance, numberOfPeriods);
            auto operatingPoints = llc.simulate_and_extract_operating_points(
                turnsRatios, magnetizingInductance);
            
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            
            INFO("Combined simulation took " << duration << " seconds");
            CHECK(duration < 15.0);
            
            REQUIRE(!topologyWaveforms.empty());
            REQUIRE(!operatingPoints.empty());
            
            // Verify we got magnetic waveforms data
            REQUIRE(!operatingPoints[0].get_excitations_per_winding().empty());
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // NaN Debug Tests
    // ═══════════════════════════════════════════════════════════════════════

    TEST_CASE("Test_Llc_Wizard_NaN_Debug", "[converter-model][llc-topology][debug]") {
        std::cout << "\n=== LLC WIZARD NaN DEBUG TEST ===" << std::endl;
        
        json llcJson;
        json inputVoltage;
        
        // Exact values from the wizard
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 360.0;  
        inputVoltage["maximum"] = 440.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Full Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 120000;
        llcJson["resonantFrequency"] = 100000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = true;
        llcJson["operatingPoints"] = json::array();
        
        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {48.0};
            LlcOperatingPointJson["outputCurrents"] = {10.4167};  // 500W / 48V
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }
        
        std::cout << "Creating LLC with wizard parameters..." << std::endl;
        OpenMagnetics::Llc llc(llcJson);
        
        std::cout << "Processing design requirements..." << std::endl;
        auto designRequirements = llc.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal().has_value()) {
                turnsRatios.push_back(tr.get_nominal().value());
                std::cout << "Turns Ratio: " << tr.get_nominal().value() << std::endl;
            }
        }
        
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(
            designRequirements.get_magnetizing_inductance());
        std::cout << "Magnetizing Inductance: " << magnetizingInductance << std::endl;
        
        std::cout << "Processing operating points..." << std::endl;
        auto operatingPoints = llc.process_operating_points(turnsRatios, magnetizingInductance);
        
        std::cout << "Got " << operatingPoints.size() << " operating points" << std::endl;
        
        for (size_t opIdx = 0; opIdx < operatingPoints.size(); ++opIdx) {
            const auto& op = operatingPoints[opIdx];
            std::cout << "\n--- Operating Point " << opIdx << " ---" << std::endl;
            
            auto nameOpt = op.get_name();
            if (nameOpt.has_value()) {
                std::cout << "Name: " << nameOpt.value() << std::endl;
            }
            
            std::cout << "Excitations: " << op.get_excitations_per_winding().size() << std::endl;
            
            for (size_t widx = 0; widx < op.get_excitations_per_winding().size(); ++widx) {
                auto exc = op.get_excitations_per_winding()[widx];
                std::cout << "  Winding " << widx << ": ";
                
                if (exc.get_current() && exc.get_current()->get_waveform()) {
                    auto wf = exc.get_current()->get_waveform().value();
                    auto data = wf.get_data();
                    
                    // Count NaN values
                    int nanCount = 0;
                    int firstNan = -1;
                    for (size_t i = 0; i < data.size(); ++i) {
                        if (std::isnan(data[i])) {
                            if (firstNan == -1) firstNan = i;
                            nanCount++;
                        }
                    }
                    
                    std::cout << "Data size=" << data.size();
                    if (nanCount > 0) {
                        std::cout << " NaN count=" << nanCount 
                                  << " first at index=" << firstNan;
                        if (!data.empty() && !std::isnan(data[0])) {
                            std::cout << " first_valid=" << data[0];
                        }
                        if (firstNan > 0) {
                            std::cout << " last_valid_before_nan=" << data[firstNan-1];
                        }
                    } else {
                        auto maxIt = std::max_element(data.begin(), data.end(),
                            [](double a, double b) { return std::abs(a) < std::abs(b); });
                        if (maxIt != data.end()) {
                            std::cout << " max_abs=" << *maxIt;
                        }
                    }
                    std::cout << std::endl;
                    
                    // Print first 10 values
                    std::cout << "    First 10: ";
                    for (size_t i = 0; i < std::min(size_t(10), data.size()); ++i) {
                        if (std::isnan(data[i])) std::cout << "NaN ";
                        else std::cout << data[i] << " ";
                    }
                    std::cout << std::endl;
                    
                } else {
                    std::cout << "No current waveform" << std::endl;
                }
            }
        }
        
        // Check that no NaN values exist in any waveform
        bool hasNan = false;
        for (auto op : operatingPoints) {
            for (auto exc : op.get_excitations_per_winding()) {
                if (exc.get_current() && exc.get_current()->get_waveform()) {
                    auto data = exc.get_current()->get_waveform()->get_data();
                    for (double val : data) {
                        if (std::isnan(val)) {
                            hasNan = true;
                            break;
                        }
                    }
                }
            }
        }
        
        std::cout << "\n=== RESULT ===" << std::endl;
        if (hasNan) {
            std::cout << "FAIL: Found NaN values in waveforms!" << std::endl;
        } else {
            std::cout << "PASS: No NaN values found!" << std::endl;
        }
        
        CHECK_FALSE(hasNan);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // LLC Classes Tests
    // ═══════════════════════════════════════════════════════════════════════

    TEST_CASE("Test_Llc_AutoDesign_HalfBridge", "[converter-model][llc-topology][auto-design]") {
        json llcJson;
        json inputVoltage;

        // Use proper input voltage range for 12V output
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 350.0;
        inputVoltage["maximum"] = 450.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = false;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);

        SECTION("Input validation - valid configuration") {
            CHECK(llc.run_checks(false) == true);
        }

        SECTION("Design requirements calculated correctly") {
            auto req = llc.process_design_requirements();
            
            // Check turns ratio: n = (Vin * k_bridge) / Vout for LLC
            // Half-bridge: Vo = Vin * 0.5 = 200V, so n = 200V / 12V = 16.67
            double Vin = 400.0;
            double Vout = 12.0;
            double k_bridge = 0.5;  // Half-bridge
            double expectedN = (Vin * k_bridge) / Vout;
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
            
            // Check magnetizing inductance is positive and reasonable
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            CHECK(Lm > 0);
            CHECK(Lm > 10e-6);   // At least 10 µH
            CHECK(Lm < 10e-3);   // Less than 10 mH
            
            // Check resonant tank values are computed
            double Lr = llc.get_computed_resonant_inductance();
            double Cr = llc.get_computed_resonant_capacitance();
            CHECK(Lr > 0);
            CHECK(Cr > 0);
        }

        SECTION("Operating points structure is valid") {
            auto req = llc.process_design_requirements();
            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            
            auto ops = llc.process_operating_points(turnsRatios, Lm);
            
            REQUIRE(ops.size() > 0);
            REQUIRE(ops[0].get_excitations_per_winding().size() >= 2);
        }
    }

    TEST_CASE("Test_Llc_AutoDesign_FullBridge", "[converter-model][llc-topology][auto-design]") {
        json llcJson;
        json inputVoltage;

        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 350.0;
        inputVoltage["maximum"] = 450.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Full Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.35;
        llcJson["integratedResonantInductor"] = false;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {48.0};
            LlcOperatingPointJson["outputCurrents"] = {5.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        OpenMagnetics::Llc llc(llcJson);

        SECTION("Full-bridge turns ratio calculation") {
            auto req = llc.process_design_requirements();
            // Full-bridge: Vo = Vin * 1.0 = 400V, so n = 400V / 48V = 8.33
            double expectedN = 400.0 * 1.0 / 48.0;
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
        }

        SECTION("Bridge voltage factor for full-bridge") {
            double factor = llc.get_bridge_voltage_factor();
            REQUIRE_THAT(factor, Catch::Matchers::WithinAbs(1.0, 0.001));
        }
    }

    TEST_CASE("Test_AdvancedLlc_UserSpecified", "[converter-model][llc-topology][advanced][user-specified]") {
        json llcJson;
        json inputVoltage;

        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 350.0;
        inputVoltage["maximum"] = 450.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = false;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        // User-specified values
        llcJson["desiredTurnsRatios"] = {15.0};  // Custom turns ratio
        llcJson["desiredMagnetizingInductance"] = 150e-6;  // 150 µH

        OpenMagnetics::AdvancedLlc advancedLlc(llcJson);

        SECTION("AdvancedLlc uses user-specified turns ratio") {
            auto inputs = advancedLlc.process();
            auto req = inputs.get_design_requirements();
            
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(15.0, 0.01));
        }

        SECTION("AdvancedLlc uses user-specified magnetizing inductance") {
            auto inputs = advancedLlc.process();
            auto req = inputs.get_design_requirements();
            
            double computedLm = resolve_dimensional_values(req.get_magnetizing_inductance());
            REQUIRE_THAT(computedLm, Catch::Matchers::WithinAbs(150e-6, 150e-6 * 0.01));
        }

        SECTION("AdvancedLlc generates valid operating points") {
            auto inputs = advancedLlc.process();
            auto ops = inputs.get_operating_points();
            
            REQUIRE(ops.size() > 0);
            
            // Check all operating points have valid excitations
            for (const auto& op : ops) {
                REQUIRE(op.get_excitations_per_winding().size() >= 2);
            }
        }
    }

    TEST_CASE("Test_AdvancedLlc_WithResonantComponents", "[converter-model][llc-topology][advanced][resonant-components]") {
        json llcJson;
        json inputVoltage;

        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 350.0;
        inputVoltage["maximum"] = 450.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 200000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = false;
        llcJson["operatingPoints"] = json::array();

        {
            json LlcOperatingPointJson;
            LlcOperatingPointJson["ambientTemperature"] = 25.0;
            LlcOperatingPointJson["outputVoltages"] = {12.0};
            LlcOperatingPointJson["outputCurrents"] = {10.0};
            LlcOperatingPointJson["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(LlcOperatingPointJson);
        }

        // User-specified resonant tank components
        llcJson["desiredTurnsRatios"] = {15.0};
        llcJson["desiredMagnetizingInductance"] = 150e-6;
        llcJson["desiredResonantInductance"] = 25e-6;      // 25 µH
        llcJson["desiredResonantCapacitance"] = 100e-9;    // 100 nF

        OpenMagnetics::AdvancedLlc advancedLlc(llcJson);

        SECTION("AdvancedLlc with user resonant components generates waveforms") {
            auto inputs = advancedLlc.process();
            auto ops = inputs.get_operating_points();
            
            REQUIRE(ops.size() > 0);
            
            // Check secondary current is calculated
            const auto& secondaryExcitation = ops[0].get_excitations_per_winding()[1];
            auto secondaryCurrent = secondaryExcitation.get_current().value();
            
            double maxSecondary = 0;
            if (secondaryCurrent.get_waveform()) {
                auto secondaryWf = secondaryCurrent.get_waveform().value();
                for (auto& v : secondaryWf.get_data()) {
                    if (std::abs(v) > maxSecondary) maxSecondary = std::abs(v);
                }
            }
            
            CHECK(maxSecondary > 0.1);
        }
    }

    TEST_CASE("Test_Llc_InvalidOperatingPoint_Warning", "[converter-model][llc-topology][validation]") {
        json llcJson;
        json inputVoltage;

        // This configuration has Vi < Vo at minimum input (the wizard issue)
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 360.0;  // Too low for 48V output with n=4.17
        inputVoltage["maximum"] = 440.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 120000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = true;
        llcJson["operatingPoints"] = json::array();

        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {48.0};  // 48V output
            op["outputCurrents"] = {10.4167};  // 500W / 48V
            op["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Llc llc(llcJson);
        auto req = llc.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : req.get_turns_ratios()) {
            if (tr.get_nominal().has_value()) {
                turnsRatios.push_back(tr.get_nominal().value());
            }
        }
        
        double magnetizingInductance = resolve_dimensional_values(req.get_magnetizing_inductance());
        auto operatingPoints = llc.process_operating_points(turnsRatios, magnetizingInductance);
        
        // At minimum input voltage, secondary current should be 0 (warning should be printed)
        REQUIRE(operatingPoints.size() == 3);  // min, nom, max
        
        // Check minimum input operating point has 0 secondary current
        const auto& minOp = operatingPoints[0];  // Should be 360V
        const auto& secondaryExcitation = minOp.get_excitations_per_winding()[1];
        auto secondaryCurrent = secondaryExcitation.get_current().value();
        
        double maxSecondary = 0;
        if (secondaryCurrent.get_waveform()) {
            auto secondaryWf = secondaryCurrent.get_waveform().value();
            for (auto& v : secondaryWf.get_data()) {
                if (std::abs(v) > maxSecondary) maxSecondary = std::abs(v);
            }
        }
        
        // At minimum input where Vi < Vo, secondary current is 0
        CHECK(maxSecondary == 0.0);
        
        // But at maximum input, secondary current should be non-zero
        // Nominal may also have Vi ≈ Vo, so check max input instead
        const auto& maxOp = operatingPoints[2];  // Should be 440V
        const auto& secondaryExcitationMax = maxOp.get_excitations_per_winding()[1];
        auto secondaryCurrentMax = secondaryExcitationMax.get_current().value();
        
        double maxSecondaryMax = 0;
        if (secondaryCurrentMax.get_waveform()) {
            auto secondaryWfMax = secondaryCurrentMax.get_waveform().value();
            for (auto& v : secondaryWfMax.get_data()) {
                if (std::abs(v) > maxSecondaryMax) maxSecondaryMax = std::abs(v);
            }
        }
        
        CHECK(maxSecondaryMax > 0.1);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Wizard Tests
    // ═══════════════════════════════════════════════════════════════════════

    TEST_CASE("Test_Llc_Wizard_Default_48V", "[converter-model][llc-topology][wizard]") {
        std::cout << "\n=== Testing LLC with Wizard Defaults (48V output) ===" << std::endl;
        
        json llcJson;
        json inputVoltage;
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 360.0;
        inputVoltage["maximum"] = 440.0;
        llcJson["inputVoltage"] = inputVoltage;
        llcJson["bridgeType"] = "Half Bridge";
        llcJson["minSwitchingFrequency"] = 80000;
        llcJson["maxSwitchingFrequency"] = 120000;
        llcJson["qualityFactor"] = 0.4;
        llcJson["integratedResonantInductor"] = true;
        llcJson["magnetizingInductance"] = 200e-6;
        llcJson["operatingPoints"] = json::array();
        
        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {48.0};  // 48V output like wizard
            op["outputCurrents"] = {10.4167};  // 500W / 48V
            op["switchingFrequency"] = 100000;
            llcJson["operatingPoints"].push_back(op);
        }
        
        OpenMagnetics::Llc llc(llcJson);
        auto designRequirements = llc.process_design_requirements();
        
        std::cout << "  Input Voltage: 400V" << std::endl;
        std::cout << "  Output Voltage: 48V" << std::endl;
        std::cout << "  Output Current: 10.42A" << std::endl;
        std::cout << "  Output Power: 500W" << std::endl;
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal().has_value()) {
                turnsRatios.push_back(tr.get_nominal().value());
                std::cout << "  Computed Turns Ratio: " << tr.get_nominal().value() << std::endl;
            }
        }
        
        double magnetizingInductance = 200e-6;
        auto operatingPoints = llc.process_operating_points(turnsRatios, magnetizingInductance);
        
        REQUIRE(operatingPoints.size() > 0);
        
        // Check maximum input voltage case (should have non-zero secondary current)
        // operatingPoints[0] = minimum input (360V) - may have 0 current if design is invalid
        // operatingPoints[1] = nominal input (400V) - should have current
        // operatingPoints[2] = maximum input (440V) - should have current
        size_t opIndex = operatingPoints.size() > 2 ? 2 : 0;
        REQUIRE(operatingPoints[opIndex].get_excitations_per_winding().size() >= 2);
        
        const auto& secondaryExcitation = operatingPoints[opIndex].get_excitations_per_winding()[1];
        auto secondaryCurrent = secondaryExcitation.get_current().value();
        
        double maxSecondary = 0;
        if (secondaryCurrent.get_waveform()) {
            auto secondaryWf = secondaryCurrent.get_waveform().value();
            for (auto& v : secondaryWf.get_data()) {
                if (std::abs(v) > maxSecondary) maxSecondary = std::abs(v);
            }
        }
        
        std::cout << "  Secondary Current Max: " << maxSecondary << " A" << std::endl;
        
        // This should fail if the bug exists
        REQUIRE(maxSecondary > 0.1);
        
        std::cout << "  ✓ PASS: Secondary current is " << maxSecondary << " A" << std::endl;
    }
}
