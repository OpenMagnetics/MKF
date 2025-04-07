#include "Painter.h"
#include "Topology.h"
#include "Utils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

SUITE(Topology) {
    double maximumError = 0.1;

    TEST(Test_Flyback_CCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 110;
        inputVoltage["maximum"] = 240;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDrainSourceVoltage"] = 350;
        flybackInputsJson["currentRippleRatio"] = 0.3;
        flybackInputsJson["efficiency"] = 0.8;
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

        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);
    }

    TEST(Test_Flyback_Maximum_Duty_Cycle_CCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 57;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDutyCycle"] = 0.44;
        flybackInputsJson["currentRippleRatio"] = 0.2;
        flybackInputsJson["efficiency"] = 0.85;
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

        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);
    }

    TEST(Test_Flyback_Maximum_Duty_Cycle_DCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 57;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDutyCycle"] = 0.44;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 0.85;
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

        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
    }

    TEST(Test_Flyback_DCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 110;
        inputVoltage["maximum"] = 240;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDrainSourceVoltage"] = 350;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 0.85;
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

        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
    }

    TEST(Test_Advanced_Flyback_CCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 110;
        inputVoltage["maximum"] = 140;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["desiredInductance"] = 950e-6;
        flybackInputsJson["desiredTurnsRatios"] = {10, 20};
        flybackInputsJson["desiredDutyCycle"] = {0.6};
        flybackInputsJson["efficiency"] = 0.8;
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

        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);
    }

    TEST(Test_Advanced_Flyback_DCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 210;
        inputVoltage["maximum"] = 240;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["desiredInductance"] = 150e-6;
        flybackInputsJson["desiredDeadTime"] = {1e-6};
        flybackInputsJson["desiredDutyCycle"] = {0.4};
        flybackInputsJson["desiredTurnsRatios"] = {10, 12};
        flybackInputsJson["efficiency"] = 0.8;
        flybackInputsJson["operatingPoints"] = json::array();

        {
            json flybackOperatingPointJson;
            flybackOperatingPointJson["outputVoltages"] = {12, 5};
            flybackOperatingPointJson["outputCurrents"] = {3, 5};
            flybackOperatingPointJson["switchingFrequency"] = 100000;
            flybackOperatingPointJson["ambientTemperature"] = 42;
            flybackInputsJson["operatingPoints"].push_back(flybackOperatingPointJson);
        }
        OpenMagnetics::AdvancedFlyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto inputs = flybackInputs.process();

        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1] .get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
    }

    TEST(Test_Flyback_Bug_Web0) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 57;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0;
        flybackInputsJson["desiredInductance"] = 14.7e-6;
        flybackInputsJson["desiredDutyCycle"] = {0.44};
        flybackInputsJson["desiredTurnsRatios"] = {2};
        flybackInputsJson["efficiency"] = 0.85;
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

        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_PRIMARY);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY);
    }
}