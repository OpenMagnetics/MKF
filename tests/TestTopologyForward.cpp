#include <source_location>
#include "support/Painter.h"
#include "converter_models/SingleSwitchForward.h"
#include "converter_models/TwoSwitchForward.h"
#include "converter_models/ActiveClampForward.h"
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

}  // namespace
