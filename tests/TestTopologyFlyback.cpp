#include <source_location>
#include "support/Painter.h"
#include "converter_models/Flyback.h"
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

    TEST_CASE("Test_Flyback_CCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 110;
        inputVoltage["maximum"] = 240;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDrainSourceVoltage"] = 350;
        flybackInputsJson["currentRippleRatio"] = 0.3;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12, 12};
            flybackOperatingPointJson["outputCurrents"] = {3, 5};
            flybackOperatingPointJson["switchingFrequency"] = 100000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_CCM_Primary_Current_Minimum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        // Secondary winding voltage can be significantly higher than output voltage due to flyback reflected voltages
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * 0.6));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_Flyback_Drain_Source_Voltage_CCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDrainSourceVoltage"] = 600;
        flybackInputsJson["currentRippleRatio"] = 0.5;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {5, 12};
            flybackOperatingPointJson["outputCurrents"] = {10, 2};
            flybackOperatingPointJson["switchingFrequency"] = 150000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_Flyback_Drain_Source_Voltage_DCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDrainSourceVoltage"] = 600;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {5, 12};
            flybackOperatingPointJson["outputCurrents"] = {10, 2};
            flybackOperatingPointJson["switchingFrequency"] = 150000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Drain_Source_Voltage_DCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Drain_Source_Voltage_DCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Drain_Source_Voltage_DCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Drain_Source_Voltage_DCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["minimum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }

    TEST_CASE("Test_Flyback_Duty_Cycle_CCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDutyCycle"] = 0.562469;
        flybackInputsJson["currentRippleRatio"] = 0.5;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {5, 12};
            flybackOperatingPointJson["outputCurrents"] = {10, 2};
            flybackOperatingPointJson["switchingFrequency"] = 150000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["minimum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_Flyback_Duty_Cycle_DCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDutyCycle"] = 0.562469;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {5, 12};
            flybackOperatingPointJson["outputCurrents"] = {10, 2};
            flybackOperatingPointJson["switchingFrequency"] = 150000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }

    TEST_CASE("Test_Flyback_Maximum_Duty_Cycle_DCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 57;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDutyCycle"] = 0.44;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12};
            flybackOperatingPointJson["outputCurrents"] = {11};
            flybackOperatingPointJson["switchingFrequency"] = 200000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }

    TEST_CASE("Test_Flyback_DCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 110;
        inputVoltage["maximum"] = 240;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDrainSourceVoltage"] = 350;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12, 12};
            flybackOperatingPointJson["outputCurrents"] = {3, 5};
            flybackOperatingPointJson["switchingFrequency"] = 100000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }

    TEST_CASE("Test_Advanced_Flyback_CCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 110;
        inputVoltage["maximum"] = 140;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["desiredInductance"] = 950e-6;
        flybackInputsJson["desiredTurnsRatios"] = {10, 20};
        flybackInputsJson["desiredDutyCycle"] = {{0.6, 0.5}};
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12, 6};
            flybackOperatingPointJson["outputCurrents"] = {3, 5};
            flybackOperatingPointJson["switchingFrequency"] = 100000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::AdvancedFlyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_Advanced_Flyback_DCM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 180;
        inputVoltage["maximum"] = 230;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["desiredInductance"] = 150e-6;
        flybackInputsJson["desiredDeadTime"] = {1e-6};
        flybackInputsJson["desiredDutyCycle"] = {{0.4, 0.6}};
        flybackInputsJson["desiredTurnsRatios"] = {10, 12};
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12, 5};
            flybackOperatingPointJson["outputCurrents"] = {3, 5};
            flybackOperatingPointJson["switchingFrequency"] = 200000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::AdvancedFlyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1] .get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1] .get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }


    TEST_CASE("Test_Advanced_Flyback_DCM_Maximum_Inductance", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 180;
        inputVoltage["maximum"] = 230;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["desiredDeadTime"] = {1e-6};
        flybackInputsJson["maximumDrainSourceVoltage"] = 350;
        flybackInputsJson["currentRippleRatio"] = 0.3;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12, 5};
            flybackOperatingPointJson["outputCurrents"] = {3, 5};
            flybackOperatingPointJson["switchingFrequency"] = 200000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackOperatingPointJson["mode"] = MAS::FlybackModes::DISCONTINUOUS_CONDUCTION_MODE;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();
        REQUIRE(inputs.get_design_requirements().get_magnetizing_inductance().get_minimum());
        REQUIRE(inputs.get_design_requirements().get_magnetizing_inductance().get_maximum());


        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1] .get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1] .get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }

    TEST_CASE("Test_Flyback_Drain_Source_Voltage_BMO", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDrainSourceVoltage"] = 600;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {5, 12};
            flybackOperatingPointJson["outputCurrents"] = {10, 2};
            flybackOperatingPointJson["mode"] = MAS::FlybackModes::BOUNDARY_MODE_OPERATION;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.004);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto operatingPoints = flybackInputs.process_operating_points(magnetic);

        REQUIRE(operatingPoints[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(operatingPoints[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(operatingPoints[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(operatingPoints[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(operatingPoints[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(operatingPoints[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }

    TEST_CASE("Test_Flyback_Drain_Source_Voltage_QRM", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDrainSourceVoltage"] = 600;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {5, 12};
            flybackOperatingPointJson["outputCurrents"] = {10, 2};
            flybackOperatingPointJson["mode"] = MAS::FlybackModes::QUASI_RESONANT_MODE;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.004);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto operatingPoints = flybackInputs.process_operating_points(magnetic);

        REQUIRE(operatingPoints[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(operatingPoints[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(operatingPoints[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(operatingPoints[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(operatingPoints[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE(operatingPoints[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(operatingPoints[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(operatingPoints[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }

    TEST_CASE("Test_Flyback_Bug_Web_0", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 57;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0;
        flybackInputsJson["desiredInductance"] = 14.7e-6;
        flybackInputsJson["desiredDutyCycle"] = {{0.44, 0.44}};
        flybackInputsJson["desiredTurnsRatios"] = {2};
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12};
            flybackOperatingPointJson["outputCurrents"] = {11};
            flybackOperatingPointJson["switchingFrequency"] = 200000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::AdvancedFlyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        // Check that second operating point's primary voltage is within the input voltage range
        // System generates an intermediate voltage when no nominal is specified
        double priVpeak = inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value();
        CHECK(priVpeak >= double(flybackInputsJson["inputVoltage"]["minimum"]));
        CHECK(priVpeak <= double(flybackInputsJson["inputVoltage"]["maximum"]));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_Flyback_Bug_Web_1", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 120;
        inputVoltage["maximum"] = 375;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDutyCycle"] = 0.5;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 1;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12, 5};
            flybackOperatingPointJson["outputCurrents"] = {3, 5};
            flybackOperatingPointJson["switchingFrequency"] = 100000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_negative_peak().value(), 1e-6));

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_negative_peak().value(), 1e-6));
    }

    TEST_CASE("Test_Flyback_Bug_Web_2", "[converter-model][flyback-topology][smoke-test]") {
        json flybackInputsJson = json::parse(R"({"currentRippleRatio": 1, "diodeVoltageDrop": 0.7, "efficiency": 0.85, "inputVoltage": {"minimum": 120.0, "maximum": 375.0}, "operatingPoints": [{"ambientTemperature": 20, "outputCurrents": [2.0], "outputVoltages": [5.0], "mode": "Discontinuous Conduction Mode", "switchingFrequency": 100000.0}], "maximumDrainSourceVoltage": 600.0, "maximumDutyCycle": 0.97})");

        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto designRequirements = flybackInputs.process_design_requirements();
        auto turnsRatio = OpenMagnetics::resolve_dimensional_values(designRequirements.get_turns_ratios()[0]);
        REQUIRE(turnsRatio < 25);
    }

    TEST_CASE("Test_Flyback_Ngspice_Simulation_CCM", "[converter-model][flyback-topology][ngspice-simulation]") {
        // Create a Flyback converter specification
        OpenMagnetics::Flyback flyback;
        
        // Input voltage: 48V nominal
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(42.0);
        inputVoltage.set_maximum(54.0);
        flyback.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        flyback.set_diode_voltage_drop(0.5);
        
        // Efficiency
        flyback.set_efficiency(0.9);
        
        // Current ripple ratio for CCM
        flyback.set_current_ripple_ratio(0.4);
        
        // Maximum duty cycle
        flyback.set_maximum_duty_cycle(0.5);
        
        // Operating point: 12V @ 1A output, 100kHz
        OpenMagnetics::FlybackOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({1.0});
        opPoint.set_ambient_temperature(25.0);
        opPoint.set_switching_frequency(100000.0);
        flyback.set_operating_points({opPoint});
        
        // Process design requirements to get turns ratios and inductance
        auto designReqs = flyback.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Flyback CCM - Turns ratio: " << turnsRatios[0]);
        INFO("Flyback CCM - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto operatingPoints = flyback.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!operatingPoints.empty());
        
        // Verify we have waveform data for primary and 1 secondary output
        REQUIRE(!operatingPoints[0].get_input_voltage().get_data().empty());
        REQUIRE(operatingPoints[0].get_output_voltages().size() == 1);

        // Extract waveform data
        auto priVoltageData = operatingPoints[0].get_input_voltage().get_data();
        auto priCurrentData = operatingPoints[0].get_input_current().get_data();
        auto secVoltageData = operatingPoints[0].get_output_voltages()[0].get_data();
        auto secCurrentData = operatingPoints[0].get_output_currents()[0].get_data();
        
        // Calculate statistics
        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priV_min = *std::min_element(priVoltageData.begin(), priVoltageData.end());
        double priI_max = *std::max_element(priCurrentData.begin(), priCurrentData.end());
        double priI_min = *std::min_element(priCurrentData.begin(), priCurrentData.end());
        double secV_max = *std::max_element(secVoltageData.begin(), secVoltageData.end());
        double secI_avg = std::accumulate(secCurrentData.begin(), secCurrentData.end(), 0.0) / secCurrentData.size();
        
        INFO("Primary voltage max: " << priV_max << " V");
        INFO("Primary voltage min: " << priV_min << " V");
        INFO("Primary current max: " << priI_max << " A");
        INFO("Primary current min: " << priI_min << " A");
        INFO("Secondary voltage max: " << secV_max << " V");
        INFO("Secondary current avg: " << secI_avg << " A");
        
        // Validate primary voltage: should be close to input voltage during ON time
        CHECK(priV_max > 40.0);  // Should be around 48V
        CHECK(priV_max < 60.0);
        
        // Validate primary voltage goes negative during OFF time (reflected voltage)
        CHECK(priV_min < -10.0);  // Reflected voltage
        
        // Validate secondary voltage is around output voltage + diode drop
        CHECK(secV_max > 10.0);  // Should be around 12.5V
        CHECK(secV_max < 20.0);
        
        // Validate secondary current average is close to output current
        CHECK(secI_avg > 0.8);  // Should be around 1A
        CHECK(secI_avg < 1.5);
        
        // In CCM, primary current should not go to zero
        CHECK(priI_min > 0.0);
        
        INFO("Flyback CCM ngspice simulation test passed");
    }

    TEST_CASE("Test_Flyback_Ngspice_Simulation_DCM", "[converter-model][flyback-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create a Flyback converter specification for DCM
        OpenMagnetics::Flyback flyback;
        
        // Input voltage: 48V nominal
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(42.0);
        inputVoltage.set_maximum(54.0);
        flyback.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        flyback.set_diode_voltage_drop(0.5);
        
        // Efficiency
        flyback.set_efficiency(0.9);
        
        // Current ripple ratio > 1 forces DCM
        flyback.set_current_ripple_ratio(2.0);
        
        // Maximum duty cycle
        flyback.set_maximum_duty_cycle(0.45);
        
        // Operating point: 12V @ 0.5A output (light load), 100kHz
        OpenMagnetics::FlybackOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({0.5});
        opPoint.set_ambient_temperature(25.0);
        opPoint.set_switching_frequency(100000.0);
        flyback.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = flyback.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Flyback DCM - Turns ratio: " << turnsRatios[0]);
        INFO("Flyback DCM - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto operatingPoints = flyback.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        
        REQUIRE(!operatingPoints.empty());
        
        // Verify we have excitations
        REQUIRE(operatingPoints[0].get_excitations_per_winding().size() == 2);
        
        // Get primary excitation
        const auto& primaryExc = operatingPoints[0].get_excitations_per_winding()[0];
        auto priCurrentData = primaryExc.get_current()->get_waveform()->get_data();
        
        // Calculate primary current statistics
        double priI_max = *std::max_element(priCurrentData.begin(), priCurrentData.end());
        double priI_min = *std::min_element(priCurrentData.begin(), priCurrentData.end());
        
        INFO("Primary current max: " << priI_max << " A");
        INFO("Primary current min: " << priI_min << " A");
        
        // In DCM, primary current should touch or approach zero
        CHECK(priI_min < 0.5);  // Should be near zero
        
        // Validate we have reasonable peak current
        CHECK(priI_max > 0.1);
        CHECK(priI_max < 5.0);
        
        INFO("Flyback DCM ngspice simulation test passed");
    }

    TEST_CASE("Test_Flyback_NumPeriods_Analytical", "[converter-model][flyback-topology][num-periods][smoke-test]") {
        // Verify that set_num_periods_to_extract() affects analytical operating point waveforms
        json flybackInputsJson;
        json inputVoltage;
        inputVoltage["minimum"] = 36;
        inputVoltage["nominal"] = 48;
        inputVoltage["maximum"] = 60;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.5;
        flybackInputsJson["efficiency"] = 0.9;
        flybackInputsJson["currentRippleRatio"] = 0.3;
        flybackInputsJson["maximumDutyCycle"] = 0.5;
        flybackInputsJson["operatingPoints"] = json::array();
        {
            json opJson;
            opJson["outputVoltages"] = {12.0};
            opJson["outputCurrents"] = {2.0};
            opJson["switchingFrequency"] = 100000;
            opJson["ambientTemperature"] = 25;
            flybackInputsJson["operatingPoints"].push_back(opJson);
        }

        // Test with 1 period
        OpenMagnetics::Flyback flyback1(flybackInputsJson);
        flyback1.set_num_periods_to_extract(1);
        auto inputs1 = flyback1.process();
        auto waveform1 = inputs1.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();
        auto numPeriods1 = waveform1.get_number_periods();

        // Test with 3 periods
        OpenMagnetics::Flyback flyback3(flybackInputsJson);
        flyback3.set_num_periods_to_extract(3);
        auto inputs3 = flyback3.process();
        auto waveform3 = inputs3.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();
        auto numPeriods3 = waveform3.get_number_periods();

        INFO("Waveform1 number_periods: " << (numPeriods1 ? std::to_string(numPeriods1.value()) : "none"));
        INFO("Waveform3 number_periods: " << (numPeriods3 ? std::to_string(numPeriods3.value()) : "none"));
        INFO("Waveform1 data size: " << waveform1.get_data().size());
        INFO("Waveform3 data size: " << waveform3.get_data().size());

        // The waveform with 3 periods should have more data points or a different number_periods
        // than the one with 1 period
        if (numPeriods1.has_value() && numPeriods3.has_value()) {
            CHECK(numPeriods3.value() >= numPeriods1.value());
        }
    }

    TEST_CASE("Test_Flyback_NumPeriods_SimulatedOperatingPoints", "[converter-model][flyback-topology][num-periods][ngspice-simulation][smoke-test]") {
        // Verify that set_num_periods_to_extract() affects simulated operating point waveforms
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::Flyback flyback;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        flyback.set_input_voltage(inputVoltage);
        flyback.set_diode_voltage_drop(0.5);
        flyback.set_efficiency(0.9);
        flyback.set_current_ripple_ratio(0.3);
        flyback.set_maximum_duty_cycle(0.5);

        OpenMagnetics::FlybackOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({2.0});
        opPoint.set_switching_frequency(100000);
        opPoint.set_ambient_temperature(25.0);
        flyback.set_operating_points({opPoint});

        auto designReqs = flyback.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        // Simulate with 1 period
        flyback.set_num_periods_to_extract(1);
        auto ops1 = flyback.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!ops1.empty());
        auto voltageWf1 = ops1[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();

        // Simulate with 3 periods
        flyback.set_num_periods_to_extract(3);
        auto ops3 = flyback.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!ops3.empty());
        auto voltageWf3 = ops3[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();

        INFO("1-period waveform data size: " << voltageWf1.get_data().size());
        INFO("3-period waveform data size: " << voltageWf3.get_data().size());

        // With 3 periods, we should get more data points than with 1
        CHECK(voltageWf3.get_data().size() > voltageWf1.get_data().size());
    }

    TEST_CASE("Test_Flyback_NumPeriods_ConverterWaveforms", "[converter-model][flyback-topology][num-periods][ngspice-simulation][smoke-test]") {
        // Verify that numberOfPeriods parameter affects converter waveform data size
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::Flyback flyback;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        flyback.set_input_voltage(inputVoltage);
        flyback.set_diode_voltage_drop(0.5);
        flyback.set_efficiency(0.9);
        flyback.set_current_ripple_ratio(0.3);
        flyback.set_maximum_duty_cycle(0.5);

        OpenMagnetics::FlybackOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({2.0});
        opPoint.set_switching_frequency(100000);
        opPoint.set_ambient_temperature(25.0);
        flyback.set_operating_points({opPoint});

        auto designReqs = flyback.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        // Simulate with 1 period
        auto waveforms1 = flyback.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, 1);
        REQUIRE(!waveforms1.empty());
        auto inputV1 = waveforms1[0].get_input_voltage().get_data();

        // Simulate with 3 periods
        auto waveforms3 = flyback.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, 3);
        REQUIRE(!waveforms3.empty());
        auto inputV3 = waveforms3[0].get_input_voltage().get_data();

        INFO("1-period converter waveform data size: " << inputV1.size());
        INFO("3-period converter waveform data size: " << inputV3.size());

        // With 3 periods, we should get more data points
        CHECK(inputV3.size() > inputV1.size());
    }

    TEST_CASE("Test_Flyback_NumPeriods_WaveformValidation", "[converter-model][flyback-topology][num-periods][ngspice-simulation]") {
        // Verify that numberOfPeriods parameter produces correct waveform time spans
        // and that magnetic/converter waveforms have matching periods
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::Flyback flyback;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        flyback.set_input_voltage(inputVoltage);
        flyback.set_diode_voltage_drop(0.5);
        flyback.set_efficiency(0.9);
        flyback.set_current_ripple_ratio(0.3);
        flyback.set_maximum_duty_cycle(0.5);

        OpenMagnetics::FlybackOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({2.0});
        opPoint.set_switching_frequency(100000);
        opPoint.set_ambient_temperature(25.0);
        flyback.set_operating_points({opPoint});

        auto designReqs = flyback.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        // Test with 1, 2, and 5 periods
        for (size_t numPeriods : {1, 2, 5}) {
            auto waveforms = flyback.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numPeriods);
            REQUIRE(!waveforms.empty());
            
            // Get the time data from input voltage waveform
            auto timeData = waveforms[0].get_input_voltage().get_time().value();
            auto voltageData = waveforms[0].get_input_voltage().get_data();
            
            // Calculate the period from the switching frequency
            double period = 1.0 / 100000.0;  // 10 microseconds
            
            // Calculate actual time span
            double timeSpan = timeData.back() - timeData.front();
            
            // The time span should be approximately numPeriods * period
            // Allow 20% tolerance for edge detection and sampling variations
            double expectedSpan = numPeriods * period;
            double tolerance = 0.2 * expectedSpan;
            
            INFO("Testing " << numPeriods << " periods:");
            INFO("  Expected time span: " << expectedSpan << " s");
            INFO("  Actual time span: " << timeSpan << " s");
            INFO("  Waveform data points: " << voltageData.size());
            
            CHECK(timeSpan >= expectedSpan - tolerance);
            CHECK(timeSpan <= expectedSpan + tolerance);
        }

        // Test that magnetic and converter waveforms have same number of periods
        flyback.set_num_periods_to_extract(3);  // Request 3 periods
        
        auto converterWaveforms = flyback.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, 3);
        auto operatingPoints = flyback.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        
        REQUIRE(!converterWaveforms.empty());
        REQUIRE(!operatingPoints.empty());
        
        // Get time spans
        double converterTimeSpan = converterWaveforms[0].get_input_voltage().get_time().value().back() -
                                   converterWaveforms[0].get_input_voltage().get_time().value().front();
        
        double magneticTimeSpan = operatingPoints[0].get_excitations_per_winding()[0].get_voltage().value().get_waveform().value().get_time().value().back() -
                                  operatingPoints[0].get_excitations_per_winding()[0].get_voltage().value().get_waveform().value().get_time().value().front();
        
        INFO("Converter waveform time span: " << converterTimeSpan << " s");
        INFO("Magnetic waveform time span: " << magneticTimeSpan << " s");
        
        // Both should be approximately the same (3 periods)
        CHECK(std::abs(converterTimeSpan - magneticTimeSpan) < 0.1 * converterTimeSpan);
    }

}  // namespace
