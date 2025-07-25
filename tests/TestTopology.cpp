#include "support/Painter.h"
#include "converter_models/Topology.h"
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST(Test_Flyback_Drain_Source_Voltage_CCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDrainSourceVoltage"] = 600;
        flybackInputsJson["currentRippleRatio"] = 0.5;
        flybackInputsJson["efficiency"] = 0.7;
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST(Test_Flyback_Drain_Source_Voltage_DCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDrainSourceVoltage"] = 600;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 0.7;
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);
    }

    TEST(Test_Flyback_Duty_Cycle_CCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDutyCycle"] = 0.562469;
        flybackInputsJson["currentRippleRatio"] = 0.5;
        flybackInputsJson["efficiency"] = 0.7;
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST(Test_Flyback_Duty_Cycle_DCM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDutyCycle"] = 0.562469;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 0.7;
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);
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
        flybackInputsJson["desiredDutyCycle"] = {{0.6, 0.5}};
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

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST(Test_Advanced_Flyback_DCM) {
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
        flybackInputsJson["efficiency"] = 0.8;
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

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1] .get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1] .get_current()->get_processed()->get_offset() == 0);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);
    }

    TEST(Test_Flyback_Drain_Source_Voltage_BMO) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDrainSourceVoltage"] = 600;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 0.7;
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), operatingPoints[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(operatingPoints[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), operatingPoints[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(operatingPoints[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), operatingPoints[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(operatingPoints[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), operatingPoints[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(operatingPoints[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), operatingPoints[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(operatingPoints[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), operatingPoints[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(operatingPoints[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);
    }

    TEST(Test_Flyback_Drain_Source_Voltage_QRM) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 127;
        inputVoltage["maximum"] = 382;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.6;
        flybackInputsJson["maximumDrainSourceVoltage"] = 600;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 0.7;
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), operatingPoints[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(operatingPoints[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), operatingPoints[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(operatingPoints[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), operatingPoints[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(operatingPoints[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(operatingPoints[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), operatingPoints[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(operatingPoints[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), operatingPoints[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(operatingPoints[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), operatingPoints[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(operatingPoints[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(operatingPoints[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);
    }

    TEST(Test_Flyback_Bug_Web_0) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 57;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0;
        flybackInputsJson["desiredInductance"] = 14.7e-6;
        flybackInputsJson["desiredDutyCycle"] = {{0.44, 0.44}};
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

    TEST(Test_Flyback_Bug_Web_1) {
        json flybackInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 120;
        inputVoltage["maximum"] = 375;
        flybackInputsJson["inputVoltage"] = inputVoltage;
        flybackInputsJson["diodeVoltageDrop"] = 0.7;
        flybackInputsJson["maximumDutyCycle"] = 0.5;
        flybackInputsJson["currentRippleRatio"] = 1;
        flybackInputsJson["efficiency"] = 0.85;
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

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(flybackInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][0]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() == 0);

        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(flybackInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError);
        CHECK_CLOSE(double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"]), inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value(), (double(flybackInputsJson["operatingPoints"][0]["outputVoltages"][1]) + double(flybackInputsJson["diodeVoltageDrop"])) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() == 0);
    }

    TEST(Test_Flyback_Bug_Web_2) {
        json flybackInputsJson = json::parse(R"({"currentRippleRatio": 1, "diodeVoltageDrop": 0.7, "efficiency": 0.85, "inputVoltage": {"minimum": 120.0, "maximum": 375.0}, "operatingPoints": [{"ambientTemperature": 20, "outputCurrents": [2.0], "outputVoltages": [5.0], "mode": "Discontinuous Conduction Mode", "switchingFrequency": 100000.0}], "maximumDrainSourceVoltage": 600.0, "maximumDutyCycle": 0.97})");

        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto designRequirements = flybackInputs.process_design_requirements();
        auto turnsRatio = OpenMagnetics::resolve_dimensional_values(designRequirements.get_turns_ratios()[0]);
        CHECK(turnsRatio < 25);
    }

    TEST(Test_Flyback_Bug_Web_3) {
        json flybackInputsJson = json::parse(R"({"currentRippleRatio": 1, "diodeVoltageDrop": 0.7, "efficiency": 0.85, "inputVoltage": {"minimum": 120.0, "maximum": 375.0}, "operatingPoints": [{"ambientTemperature": 20, "outputCurrents": [2.0], "outputVoltages": [5.0], "mode": "Quasi Resonant Mode"}], "maximumDrainSourceVoltage": 600.0, "maximumDutyCycle": 0.97})");

        OpenMagnetics::Flyback flybackInputs(flybackInputsJson);
        flybackInputs._assertErrors = true;

        auto designRequirements = flybackInputs.process_design_requirements();
        auto turnsRatio = OpenMagnetics::resolve_dimensional_values(designRequirements.get_turns_ratios()[0]);
        CHECK(turnsRatio < 25);
    }
}