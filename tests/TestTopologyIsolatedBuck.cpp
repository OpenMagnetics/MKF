#include "support/Painter.h"
#include "converter_models/IsolatedBuck.h"
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

SUITE(TopologyIsolatedBuck) {
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    TEST(Test_IsolatedBuck) {
        json isolatedbuckInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 72;
        isolatedbuckInputsJson["inputVoltage"] = inputVoltage;
        isolatedbuckInputsJson["diodeVoltageDrop"] = 0.7;
        isolatedbuckInputsJson["maximumSwitchCurrent"] = 0.7;
        isolatedbuckInputsJson["efficiency"] = 1;
        isolatedbuckInputsJson["operatingPoints"] = json::array();

        {
            json isolatedbuckOperatingPointJson;
            isolatedbuckOperatingPointJson["outputVoltages"] = {10, 10};
            isolatedbuckOperatingPointJson["outputCurrents"] = {0.02, 0.1};
            isolatedbuckOperatingPointJson["switchingFrequency"] = 750000;
            isolatedbuckOperatingPointJson["ambientTemperature"] = 42;
            isolatedbuckInputsJson["operatingPoints"].push_back(isolatedbuckOperatingPointJson);
        }
        OpenMagnetics::IsolatedBuck isolatedbuckInputs(isolatedbuckInputsJson);
        isolatedbuckInputs._assertErrors = true;

        auto inputs = isolatedbuckInputs.process();

        // {
        //     auto outFile = outputFilePath;
        //     outFile.append("Test_IsolatedBuck_Primary.svg");
        //     std::filesystem::remove(outFile);   
        //     Painter painter(outFile, false, true);
        //     painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
        //     painter.export_svg();
        // }
        // {
        //     auto outFile = outputFilePath;
        //     outFile.append("Test_IsolatedBuck_Secondary.svg");
        //     std::filesystem::remove(outFile);   
        //     Painter painter(outFile, false, true);
        //     painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
        //     painter.export_svg();
        // }

        // {
        //     auto outFile = outputFilePath;
        //     outFile.append("Test_IsolatedBuck_Primary_Voltage.svg");
        //     std::filesystem::remove(outFile);   
        //     Painter painter(outFile, false, true);
        //     painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
        //     painter.export_svg();
        // }
        // {
        //     auto outFile = outputFilePath;
        //     outFile.append("Test_IsolatedBuck_Secondary_Voltage.svg");
        //     std::filesystem::remove(outFile);   
        //     Painter painter(outFile, false, true);
        //     painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
        //     painter.export_svg();
        // }

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding().size() == isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"].size());
        CHECK_CLOSE(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(isolatedbuckInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["minimum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(isolatedbuckInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["minimum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(isolatedbuckInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(isolatedbuckInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

}