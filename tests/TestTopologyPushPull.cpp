#include "support/Painter.h"
#include "converter_models/PushPull.h"
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

SUITE(TopologyPushPull) {
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    TEST(Test_PushPull_CCM) {
        json pushPullInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 20;
        inputVoltage["maximum"] = 30;
        pushPullInputsJson["inputVoltage"] = inputVoltage;
        pushPullInputsJson["diodeVoltageDrop"] = 0.6;
        pushPullInputsJson["efficiency"] = 0.9;
        pushPullInputsJson["maximumSwitchCurrent"] = 1;
        pushPullInputsJson["currentRippleRatio"] = 0.3;
        pushPullInputsJson["dutyCycle"] = 0.45;
        pushPullInputsJson["operatingPoints"] = json::array();

        {
            json PushPullOperatingPointJson;
            PushPullOperatingPointJson["outputVoltages"] = {48, 5, 9};
            PushPullOperatingPointJson["outputCurrents"] = {0.7, 0.01, 0.01};
            PushPullOperatingPointJson["switchingFrequency"] = 500000;
            PushPullOperatingPointJson["ambientTemperature"] = 42;
            pushPullInputsJson["operatingPoints"].push_back(PushPullOperatingPointJson);
        }
        OpenMagnetics::PushPull pushPullInputs(pushPullInputsJson);
        pushPullInputs._assertErrors = true;

        auto inputs = pushPullInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Auxiliary_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[4].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Auxiliary_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[5].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Auxiliary_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[4].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Auxiliary_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[5].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        CHECK_CLOSE(double(pushPullInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(pushPullInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(pushPullInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(pushPullInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(pushPullInputsJson["operatingPoints"][0]["outputCurrents"][0]) / 2, inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(pushPullInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(pushPullInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_rms().value(), double(pushPullInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(pushPullInputsJson["operatingPoints"][0]["outputCurrents"][0]) / 2, inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_average().value(), double(pushPullInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(pushPullInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_processed()->get_rms().value(), double(pushPullInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_offset() > 0);

    }

    TEST(Test_PushPull_DCM) {
        json pushPullInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 500;
        inputVoltage["maximum"] = 501;
        pushPullInputsJson["inputVoltage"] = inputVoltage;
        pushPullInputsJson["diodeVoltageDrop"] = 0.6;
        pushPullInputsJson["efficiency"] = 0.9;
        pushPullInputsJson["maximumSwitchCurrent"] = 1;
        pushPullInputsJson["currentRippleRatio"] = 0.3;
        pushPullInputsJson["desiredInductance"] = 1e-3;
        pushPullInputsJson["desiredOutputInductance"] = 5e-6;
        pushPullInputsJson["desiredTurnsRatios"] = {15};
        pushPullInputsJson["dutyCycle"] = 0.45;
        pushPullInputsJson["operatingPoints"] = json::array();

        {
            json PushPullOperatingPointJson;
            PushPullOperatingPointJson["outputVoltages"] = {15};
            PushPullOperatingPointJson["outputCurrents"] = {0.8};
            PushPullOperatingPointJson["switchingFrequency"] = 500000;
            PushPullOperatingPointJson["ambientTemperature"] = 42;
            pushPullInputsJson["operatingPoints"].push_back(PushPullOperatingPointJson);
        }
        OpenMagnetics::AdvancedPushPull pushPullInputs(pushPullInputsJson);
        pushPullInputs._assertErrors = true;

        auto inputs = pushPullInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_First_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_Second_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_First_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_Second_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_First_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_Second_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_First_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_Second_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        CHECK_CLOSE(double(pushPullInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(pushPullInputsJson["inputVoltage"]["maximum"]) * maximumError);

        CHECK_CLOSE(double(pushPullInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(pushPullInputsJson["inputVoltage"]["maximum"]) * maximumError);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_offset() > 0);

    }

}