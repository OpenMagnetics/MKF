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

using namespace OpenMagnetics;

namespace {
    auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
    double maximumError = 0.1;

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
            double expectedN = 400.0 * 0.5 / (2.0 * 12.0);  // n = Vin * k_bridge / (2 * Vout)
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.02));
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

            CHECK(netlist.find("Vin") != std::string::npos);
            CHECK(netlist.find("Lr") != std::string::npos);
            CHECK(netlist.find("Cr") != std::string::npos);
            CHECK(netlist.find("Lpri") != std::string::npos);
            CHECK(netlist.find("Lsec0") != std::string::npos);
            CHECK(netlist.find("Kpri_sec0") != std::string::npos);
            CHECK(netlist.find(".tran") != std::string::npos);
            CHECK(netlist.find(".end") != std::string::npos);
            CHECK(netlist.find("Half") != std::string::npos);
            CHECK(netlist.find("S1") != std::string::npos);
            CHECK(netlist.find("S2") != std::string::npos);
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
}
