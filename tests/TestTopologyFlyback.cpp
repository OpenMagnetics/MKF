#include "support/Painter.h"
#include "converter_models/Flyback.h"
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

    TEST_CASE("Test_Flyback_CCM", "[converter-model][flyback-topology]") {
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
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_Flyback_Drain_Source_Voltage_CCM", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_Drain_Source_Voltage_DCM", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_Duty_Cycle_CCM", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_Duty_Cycle_DCM", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_Maximum_Duty_Cycle_DCM", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_DCM", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Advanced_Flyback_CCM", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Advanced_Flyback_DCM", "[converter-model][flyback-topology]") {
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


    TEST_CASE("Test_Advanced_Flyback_DCM_Maximum_Inductance", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_Drain_Source_Voltage_BMO", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_Drain_Source_Voltage_QRM", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_Bug_Web_0", "[converter-model][flyback-topology]") {
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

        REQUIRE_THAT(double(flybackInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_Flyback_Bug_Web_1", "[converter-model][flyback-topology]") {
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

    TEST_CASE("Test_Flyback_Bug_Web_2", "[converter-model][flyback-topology]") {
        json flybackInputsJson = json::parse(R"({"currentRippleRatio": 1, "diodeVoltageDrop": 0.7, "efficiency": 0.85, "inputVoltage": {"minimum": 120.0, "maximum": 375.0}, "operatingPoints": [{"ambientTemperature": 20, "outputCurrents": [2.0], "outputVoltages": [5.0], "mode": "Discontinuous Conduction Mode", "switchingFrequency": 100000.0}], "maximumDrainSourceVoltage": 600.0, "maximumDutyCycle": 0.97})");

        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto designRequirements = flybackInputs.process_design_requirements();
        auto turnsRatio = OpenMagnetics::resolve_dimensional_values(designRequirements.get_turns_ratios()[0]);
        REQUIRE(turnsRatio < 25);
    }

}  // namespace
