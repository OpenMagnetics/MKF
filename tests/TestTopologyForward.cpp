#include <source_location>
#include "support/Painter.h"
#include "converter_models/SingleSwitchForward.h"
#include "converter_models/TwoSwitchForward.h"
#include "converter_models/ActiveClampForward.h"
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

    TEST_CASE("Test_Forward_Ngspice_Simulation", "[converter-model][forward-topology][ngspice][!mayfail]") {
        json forwardInputsJson;
        json inputVoltage;

        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 3;
        forwardInputsJson["currentRippleRatio"] = 0.3;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {12};
            ForwardOperatingPointJson["outputCurrents"] = {2};
            ForwardOperatingPointJson["switchingFrequency"] = 100000;
            ForwardOperatingPointJson["ambientTemperature"] = 25;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }

        OpenMagnetics::SingleSwitchForward forward(forwardInputsJson);
        forward._assertErrors = true;
        
        // First process to get design requirements with inductance
        auto inputs = forward.process();
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance());
        
        // Get turns ratios
        std::vector<double> turnsRatios;
        for (const auto& tr : inputs.get_design_requirements().get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        
        // Now run the ngspice simulation
        auto topologyWaveforms = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(topologyWaveforms.size() >= 1);
        
        for (size_t opIndex = 0; opIndex < topologyWaveforms.size(); opIndex++) {
            auto& wf = topologyWaveforms[opIndex];
            
            // Check that time vector has reasonable values
            REQUIRE(wf.time.size() > 0);
            REQUIRE(wf.time[0] >= 0);
            
            // Check primary current waveform
            REQUIRE(wf.primaryCurrent.size() == wf.time.size());
            
            // Check demagnetization current waveform
            REQUIRE(wf.demagCurrent.size() == wf.time.size());
            
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
                outFile.append("Test_Forward_Ngspice_PrimaryCurrent_OP" + std::to_string(opIndex) + ".svg");
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
                outFile.append("Test_Forward_Ngspice_DemagCurrent_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform currentWaveform;
                currentWaveform.set_time(wf.time);
                currentWaveform.set_data(wf.demagCurrent);
                painter.paint_waveform(currentWaveform);
                painter.export_svg();
            }
            if (wf.outputVoltages.size() > 0 && wf.outputVoltages[0].size() > 0) {
                auto outFile = outputFilePath;
                outFile.append("Test_Forward_Ngspice_OutputVoltage_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform outputWaveform;
                outputWaveform.set_time(wf.time);
                outputWaveform.set_data(wf.outputVoltages[0]);
                painter.paint_waveform(outputWaveform);
                painter.export_svg();
            }
            if (wf.outputInductorCurrents.size() > 0 && wf.outputInductorCurrents[0].size() > 0) {
                auto outFile = outputFilePath;
                outFile.append("Test_Forward_Ngspice_OutputInductorCurrent_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                Waveform currentWaveform;
                currentWaveform.set_time(wf.time);
                currentWaveform.set_data(wf.outputInductorCurrents[0]);
                painter.paint_waveform(currentWaveform);
                painter.export_svg();
            }
        }
    }

}  // namespace
