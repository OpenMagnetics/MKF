#include <source_location>
#include "support/Painter.h"
#include "converter_models/IsolatedBuckBoost.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    TEST_CASE("Test_IsolatedBuckBoost", "[converter-model][isolated-buck-boost-topology][smoke-test]") {
        json isolatedBuckBoostInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 10;
        inputVoltage["maximum"] = 30;
        isolatedBuckBoostInputsJson["inputVoltage"] = inputVoltage;
        isolatedBuckBoostInputsJson["diodeVoltageDrop"] = 0;
        isolatedBuckBoostInputsJson["maximumSwitchCurrent"] = 2.5;
        isolatedBuckBoostInputsJson["efficiency"] = 1;
        isolatedBuckBoostInputsJson["operatingPoints"] = json::array();

        {
            json isolatedBuckBoostOperatingPointJson;
            isolatedBuckBoostOperatingPointJson["outputVoltages"] = {6, 5, 5};
            isolatedBuckBoostOperatingPointJson["outputCurrents"] = {0.01, 1, 0.3};
            isolatedBuckBoostOperatingPointJson["switchingFrequency"] = 400000;
            isolatedBuckBoostOperatingPointJson["ambientTemperature"] = 42;
            isolatedBuckBoostInputsJson["operatingPoints"].push_back(isolatedBuckBoostOperatingPointJson);
        }
        OpenMagnetics::IsolatedBuckBoost isolatedBuckBoostInputs(isolatedBuckBoostInputsJson);
        isolatedBuckBoostInputs._assertErrors = true;

        auto inputs = isolatedBuckBoostInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding().size() == isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"].size());
        REQUIRE_THAT(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() != 0);

        REQUIRE_THAT(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset(), 0.01));

        REQUIRE_THAT(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset(), 0.01));

        REQUIRE_THAT(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset(), 0.01));
    }

    TEST_CASE("Test_IsolatedBuckBoost_Ngspice_Simulation", "[converter-model][isolated-buck-boost-topology][ngspice]") {
        json isolatedBuckBoostInputsJson;
        json inputVoltage;

        inputVoltage["minimum"] = 10;
        inputVoltage["maximum"] = 30;
        isolatedBuckBoostInputsJson["inputVoltage"] = inputVoltage;
        isolatedBuckBoostInputsJson["diodeVoltageDrop"] = 0.7;
        isolatedBuckBoostInputsJson["maximumSwitchCurrent"] = 5;
        isolatedBuckBoostInputsJson["efficiency"] = 0.9;
        isolatedBuckBoostInputsJson["operatingPoints"] = json::array();

        {
            json isolatedBuckBoostOperatingPointJson;
            isolatedBuckBoostOperatingPointJson["outputVoltages"] = {12};
            isolatedBuckBoostOperatingPointJson["outputCurrents"] = {1};
            isolatedBuckBoostOperatingPointJson["switchingFrequency"] = 100000;
            isolatedBuckBoostOperatingPointJson["ambientTemperature"] = 25;
            isolatedBuckBoostInputsJson["operatingPoints"].push_back(isolatedBuckBoostOperatingPointJson);
        }

        OpenMagnetics::IsolatedBuckBoost isolatedBuckBoost(isolatedBuckBoostInputsJson);
        isolatedBuckBoost._assertErrors = true;
        
        // First process to get design requirements with inductance
        auto inputs = isolatedBuckBoost.process();
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance());
        
        // Get turns ratios
        std::vector<double> turnsRatios;
        for (const auto& tr : inputs.get_design_requirements().get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        
        // Now run the ngspice simulation
        auto topologyWaveforms = isolatedBuckBoost.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(topologyWaveforms.size() >= 1);
        
        for (size_t opIndex = 0; opIndex < topologyWaveforms.size(); opIndex++) {
            auto& wf = topologyWaveforms[opIndex];
            
            // Check that time vector has reasonable values
            REQUIRE(wf.time.size() > 0);
            REQUIRE(wf.time[0] >= 0);
            
            // Check primary current waveform
            REQUIRE(wf.primaryCurrent.size() == wf.time.size());
            
            // Check primary voltage waveform
            REQUIRE(wf.primaryVoltage.size() == wf.time.size());
            
            // Check output voltages if available
            if (wf.outputVoltages.size() >= 1 && wf.outputVoltages[0].size() == wf.time.size()) {
                // Verify output voltage is close to expected
                double avgOutputVoltage = 0;
                for (double v : wf.outputVoltages[0]) {
                    avgOutputVoltage += v;
                }
                avgOutputVoltage /= wf.outputVoltages[0].size();
                REQUIRE_THAT(std::abs(avgOutputVoltage), Catch::Matchers::WithinAbs(12.0, 5.0));  // Within 5V of expected 12V output
            }
            
            // Paint waveforms for visual inspection
            {
                auto outFile = outputFilePath;
                outFile.append("Test_IsolatedBuckBoost_Ngspice_PrimaryCurrent_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform currentWaveform;
                currentWaveform.set_time(wf.time);
                currentWaveform.set_data(wf.primaryCurrent);
                painter.paint_waveform(currentWaveform);
                painter.export_svg();
            }
            {
                auto outFile = outputFilePath;
                outFile.append("Test_IsolatedBuckBoost_Ngspice_PrimaryVoltage_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform voltageWaveform;
                voltageWaveform.set_time(wf.time);
                voltageWaveform.set_data(wf.primaryVoltage);
                painter.paint_waveform(voltageWaveform);
                painter.export_svg();
            }
            if (wf.outputVoltages.size() > 0 && wf.outputVoltages[0].size() > 0) {
                auto outFile = outputFilePath;
                outFile.append("Test_IsolatedBuckBoost_Ngspice_OutputVoltage_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform outputWaveform;
                outputWaveform.set_time(wf.time);
                outputWaveform.set_data(wf.outputVoltages[0]);
                painter.paint_waveform(outputWaveform);
                painter.export_svg();
            }
            if (wf.secondaryCurrents.size() > 0 && wf.secondaryCurrents[0].size() > 0) {
                auto outFile = outputFilePath;
                outFile.append("Test_IsolatedBuckBoost_Ngspice_SecondaryCurrent_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform currentWaveform;
                currentWaveform.set_time(wf.time);
                currentWaveform.set_data(wf.secondaryCurrents[0]);
                painter.paint_waveform(currentWaveform);
                painter.export_svg();
            }
        }
    }

}  // namespace
