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
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    TEST_CASE("Test_Boost", "[converter-model][boost-topology]") {
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

}  // namespace
