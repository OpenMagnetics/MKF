#include <source_location>
#include "support/Painter.h"
#include "converter_models/Boost.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
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

    TEST_CASE("Test_Boost", "[converter-model][boost-topology][smoke-test]") {
        json boostInputsJson;
        json inputVoltage;

        inputVoltage["minimum"] = 12;
        inputVoltage["maximum"] = 24;
        boostInputsJson["inputVoltage"] = inputVoltage;
        boostInputsJson["diodeVoltageDrop"] = 0.7;
        boostInputsJson["efficiency"] = 1;
        boostInputsJson["maximumSwitchCurrent"] = 8;
        boostInputsJson["operatingPoints"] = json::array();
        {
            json boostOperatingPointJson;
            boostOperatingPointJson["outputVoltage"] = 50;
            boostOperatingPointJson["outputCurrent"] = 1;
            boostOperatingPointJson["switchingFrequency"] = 100000;
            boostOperatingPointJson["ambientTemperature"] = 42;
            boostInputsJson["operatingPoints"].push_back(boostOperatingPointJson);
        }

        OpenMagnetics::Boost boostInputs(boostInputsJson);

        auto inputs = boostInputs.process();
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Boost_Primary.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Boost_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Boost_Primary_Maximum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Boost_Primary_Voltage_Maximum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        REQUIRE_THAT(double(boostInputsJson["operatingPoints"][0]["outputVoltage"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(boostInputsJson["operatingPoints"][0]["outputVoltage"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(boostInputsJson["operatingPoints"][0]["outputVoltage"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(boostInputsJson["operatingPoints"][0]["outputVoltage"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);
    }

    TEST_CASE("Test_Boost_Ngspice_Simulation", "[converter-model][boost-topology][ngspice]") {
        json boostInputsJson;
        json inputVoltage;

        inputVoltage["minimum"] = 12;
        inputVoltage["maximum"] = 24;
        boostInputsJson["inputVoltage"] = inputVoltage;
        boostInputsJson["diodeVoltageDrop"] = 0.7;
        boostInputsJson["efficiency"] = 0.9;
        boostInputsJson["maximumSwitchCurrent"] = 8;
        boostInputsJson["operatingPoints"] = json::array();
        {
            json boostOperatingPointJson;
            boostOperatingPointJson["outputVoltage"] = 48;
            boostOperatingPointJson["outputCurrent"] = 1;
            boostOperatingPointJson["switchingFrequency"] = 100000;
            boostOperatingPointJson["ambientTemperature"] = 25;
            boostInputsJson["operatingPoints"].push_back(boostOperatingPointJson);
        }

        OpenMagnetics::Boost boost(boostInputsJson);
        boost._assertErrors = true;
        
        // First process to get design requirements with inductance
        auto inputs = boost.process();
        double inductance = OpenMagnetics::resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance());
        
        // Now run the ngspice simulation
        auto topologyWaveforms = boost.simulate_and_extract_topology_waveforms(inductance);
        
        REQUIRE(topologyWaveforms.size() >= 1);
        
        for (size_t opIndex = 0; opIndex < topologyWaveforms.size(); opIndex++) {
            auto& wf = topologyWaveforms[opIndex];
            
            // Check that time vector has reasonable values
            REQUIRE(wf.time.size() > 0);
            REQUIRE(wf.time[0] >= 0);
            
            // Check inductor voltage waveform
            REQUIRE(wf.inductorVoltage.size() == wf.time.size());
            
            // Check inductor current waveform
            REQUIRE(wf.inductorCurrent.size() == wf.time.size());
            
            // Check output voltage waveform
            REQUIRE(wf.outputVoltage.size() == wf.time.size());
            
            // Verify output voltage is close to expected
            double avgOutputVoltage = 0;
            for (double v : wf.outputVoltage) {
                avgOutputVoltage += v;
            }
            avgOutputVoltage /= wf.outputVoltage.size();
            REQUIRE_THAT(avgOutputVoltage, Catch::Matchers::WithinAbs(48.0, 10.0));  // Within 10V of expected 48V output (includes diode drop and ripple)
            
            // Paint waveforms for visual inspection
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Boost_Ngspice_InductorCurrent_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform currentWaveform;
                currentWaveform.set_time(wf.time);
                currentWaveform.set_data(wf.inductorCurrent);
                painter.paint_waveform(currentWaveform);
                painter.export_svg();
            }
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Boost_Ngspice_InductorVoltage_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform voltageWaveform;
                voltageWaveform.set_time(wf.time);
                voltageWaveform.set_data(wf.inductorVoltage);
                painter.paint_waveform(voltageWaveform);
                painter.export_svg();
            }
            {
                auto outFile = outputFilePath;
                outFile.append("Test_Boost_Ngspice_OutputVoltage_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform outputWaveform;
                outputWaveform.set_time(wf.time);
                outputWaveform.set_data(wf.outputVoltage);
                painter.paint_waveform(outputWaveform);
                painter.export_svg();
            }
        }
    }

}  // namespace
