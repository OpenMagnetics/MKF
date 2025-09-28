#include "support/Painter.h"
#include "converter_models/Boost.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(TopologyBoost) {
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    TEST(Test_Boost) {
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

        CHECK_CLOSE(double(boostInputsJson["operatingPoints"][0]["outputVoltage"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(boostInputsJson["operatingPoints"][0]["outputVoltage"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(boostInputsJson["operatingPoints"][0]["outputVoltage"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(boostInputsJson["operatingPoints"][0]["outputVoltage"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);
    }
}