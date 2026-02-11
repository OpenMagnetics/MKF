#include <source_location>
#include "support/Painter.h"
#include "converter_models/SingleSwitchForward.h"
#include "converter_models/TwoSwitchForward.h"
#include "converter_models/ActiveClampForward.h"
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

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    TEST_CASE("Test_SingleSwitchForward_CCM", "[converter-model][single-switch-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 0.3;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {5};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::SingleSwitchForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Demagnetization_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Demagnetization_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(forwardInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(forwardInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(forwardInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(forwardInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        double outputCurrent = (inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_peak().value() + inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset()) / 2;
        REQUIRE_THAT(double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(outputCurrent, double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_SingleSwitchForward_DCM", "[converter-model][single-switch-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 2;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["desiredInductance"] = 1e-3;
        forwardInputsJson["desiredOutputInductances"] = {5e-6};
        forwardInputsJson["desiredTurnsRatios"] = {1, 2};
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {1};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::AdvancedSingleSwitchForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Demagnetization_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Demagnetization_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
    }

    TEST_CASE("Test_ActiveClampForward_CCM", "[converter-model][active-clamp-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 0.3;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {5};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::ActiveClampForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_CCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_CCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_CCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_CCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(forwardInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(forwardInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        double outputCurrent = (inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_peak().value() + inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset()) / 2;
        REQUIRE_THAT(double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(outputCurrent, double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_ActiveClampForward_DCM", "[converter-model][active-clamp-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 2;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["desiredInductance"] = 1e-3;
        forwardInputsJson["desiredOutputInductances"] = {5e-6};
        forwardInputsJson["desiredTurnsRatios"] = {0.5};
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {1};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::AdvancedActiveClampForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_DCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_DCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_DCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_DCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);
    }

    TEST_CASE("Test_TwoSwitchForward_CCM", "[converter-model][two-switch-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 0.3;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {5};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::TwoSwitchForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_CCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_CCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_CCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_CCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(forwardInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(forwardInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        double outputCurrent = (inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_peak().value() + inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset()) / 2;
        REQUIRE_THAT(double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(outputCurrent, double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_TwoSwitchForward_DCM", "[converter-model][two-switch-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 2;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["desiredInductance"] = 1e-3;
        forwardInputsJson["desiredOutputInductances"] = {5e-6};
        forwardInputsJson["desiredTurnsRatios"] = {0.5};
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {1};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::AdvancedTwoSwitchForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_DCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_DCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_DCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_DCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
    }

    TEST_CASE("Test_SingleSwitchForward_Ngspice_Simulation", "[converter-model][single-switch-forward-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create a Single Switch Forward converter specification
        OpenMagnetics::SingleSwitchForward forward;
        
        // Input voltage
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(36.0);
        inputVoltage.set_maximum(72.0);
        forward.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        forward.set_diode_voltage_drop(0.5);
        
        // Efficiency
        forward.set_efficiency(0.9);
        
        // Current ripple ratio
        forward.set_current_ripple_ratio(0.3);
        
        // Duty cycle
        forward.set_duty_cycle(0.4);
        
        // Operating point: 5V @ 5A output, 200kHz
        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({5.0});
        opPoint.set_output_currents({5.0});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = forward.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Single Switch Forward - Turns ratio: " << turnsRatios[0]);
        INFO("Single Switch Forward - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto converterWaveforms = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!converterWaveforms.empty());
        
        // Verify we have excitations
        // ConverterWaveforms doesn't have excitations_per_winding - check input/output instead
    REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());
        
        // Get primary excitation
        // Primary excitation data is now directly in ConverterWaveforms
        // Voltage data is directly available in ConverterWaveforms
        // Current data is directly available in ConverterWaveforms
        
        // Extract waveform data
        auto priVoltageData = converterWaveforms[0].get_input_voltage().get_data();
        auto priCurrentData = converterWaveforms[0].get_input_current().get_data();
        
        // Calculate statistics
        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priI_max = *std::max_element(priCurrentData.begin(), priCurrentData.end());
        double priI_min = *std::min_element(priCurrentData.begin(), priCurrentData.end());
        
        INFO("Primary voltage max: " << priV_max << " V");
        INFO("Primary current max: " << priI_max << " A");
        INFO("Primary current min: " << priI_min << " A");
        
        // Validate primary voltage: should be close to input voltage during ON time
        CHECK(priV_max > 30.0);  // Should be around 48V
        CHECK(priV_max < 80.0);
        
        // In CCM, current should stay positive (allow for numerical noise)
        CHECK(priI_min > -0.001);
        
        // Validate we have reasonable peak current
        CHECK(priI_max > 0.5);
        CHECK(priI_max < 10.0);
        
        INFO("Single Switch Forward ngspice simulation test passed");
    }

    TEST_CASE("Test_TwoSwitchForward_Ngspice_Simulation", "[converter-model][two-switch-forward-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create a Two Switch Forward converter specification
        OpenMagnetics::TwoSwitchForward forward;
        
        // Input voltage
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(36.0);
        inputVoltage.set_maximum(72.0);
        forward.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        forward.set_diode_voltage_drop(0.5);
        
        // Efficiency
        forward.set_efficiency(0.9);
        
        // Current ripple ratio
        forward.set_current_ripple_ratio(0.3);
        
        // Duty cycle
        forward.set_duty_cycle(0.4);
        
        // Operating point: 5V @ 5A output, 200kHz
        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({5.0});
        opPoint.set_output_currents({5.0});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = forward.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Two Switch Forward - Turns ratio: " << turnsRatios[0]);
        INFO("Two Switch Forward - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto converterWaveforms = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!converterWaveforms.empty());
        
        // Verify we have excitations
        // ConverterWaveforms doesn't have excitations_per_winding - check input/output instead
    REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());
        
        // Get primary excitation
        // Primary excitation data is now directly in ConverterWaveforms
        // Voltage data is directly available in ConverterWaveforms
        // Current data is directly available in ConverterWaveforms
        
        // Extract waveform data
        auto priVoltageData = converterWaveforms[0].get_input_voltage().get_data();
        
        // Calculate statistics
        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        
        INFO("Primary voltage max: " << priV_max << " V");
        
        // Validate primary voltage: should be close to input voltage during ON time
        CHECK(priV_max > 30.0);  // Should be around 48V
        CHECK(priV_max < 80.0);
        
        INFO("Two Switch Forward ngspice simulation test passed");
    }

    TEST_CASE("Test_ActiveClampForward_Ngspice_Simulation", "[converter-model][active-clamp-forward-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create an Active Clamp Forward converter specification
        OpenMagnetics::ActiveClampForward forward;
        
        // Input voltage
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(36.0);
        inputVoltage.set_maximum(72.0);
        forward.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        forward.set_diode_voltage_drop(0.5);
        
        // Efficiency
        forward.set_efficiency(0.9);
        
        // Current ripple ratio
        forward.set_current_ripple_ratio(0.3);
        
        // Duty cycle
        forward.set_duty_cycle(0.45);
        
        // Operating point: 5V @ 5A output, 200kHz
        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({5.0});
        opPoint.set_output_currents({5.0});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = forward.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Active Clamp Forward - Turns ratio: " << turnsRatios[0]);
        INFO("Active Clamp Forward - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto converterWaveforms = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!converterWaveforms.empty());
        
        // Verify we have excitations
        // ConverterWaveforms doesn't have excitations_per_winding - check input/output instead
    REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());
        
        // Get primary excitation
        // Primary excitation data is now directly in ConverterWaveforms
        // Voltage data is directly available in ConverterWaveforms
        // Current data is directly available in ConverterWaveforms
        
        // Extract waveform data
        auto priVoltageData = converterWaveforms[0].get_input_voltage().get_data();
        auto priCurrentData = converterWaveforms[0].get_input_current().get_data();
        
        // Calculate statistics
        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priV_min = *std::min_element(priVoltageData.begin(), priVoltageData.end());
        double priI_avg = std::accumulate(priCurrentData.begin(), priCurrentData.end(), 0.0) / priCurrentData.size();
        
        INFO("Primary voltage max: " << priV_max << " V");
        INFO("Primary voltage min: " << priV_min << " V");
        INFO("Primary current avg: " << priI_avg << " A");
        
        // Validate primary voltage: should be close to input voltage during ON time
        CHECK(priV_max > 30.0);  // Should be around 48V
        CHECK(priV_max < 80.0);
        
        // Active clamp should have negative voltage during reset
        CHECK(priV_min < 0.0);
        
        INFO("Active Clamp Forward ngspice simulation test passed");
    }

}  // namespace
