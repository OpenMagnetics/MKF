#include "support/Painter.h"
#include "converter_models/IsolatedBuckBoost.h"
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

SUITE(TopologyIsolatedBuckBoost) {
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    TEST(Test_IsolatedBuckBoost) {
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

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding().size() == isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"].size());
        CHECK_CLOSE(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() != 0);

        CHECK_CLOSE(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK_CLOSE(0, inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset(), 0.01);

        CHECK_CLOSE(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        CHECK_CLOSE(0, inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset(), 0.01);

        CHECK_CLOSE(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK_CLOSE(0, inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset(), 0.01);
    }

}