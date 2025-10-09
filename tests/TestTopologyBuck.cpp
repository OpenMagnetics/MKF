#include "support/Painter.h"
#include "converter_models/Buck.h"
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

SUITE(TopologyBuck) {
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    TEST(Test_Buck) {
        json buckInputsJson;
        json inputVoltage;

        inputVoltage["minimum"] = 20;
        inputVoltage["maximum"] = 240;
        buckInputsJson["inputVoltage"] = inputVoltage;
        buckInputsJson["diodeVoltageDrop"] = 0.7;
        buckInputsJson["efficiency"] = 0.9;
        buckInputsJson["maximumSwitchCurrent"] = 8;
        buckInputsJson["operatingPoints"] = json::array();
        {
            json buckOperatingPointJson;
            buckOperatingPointJson["outputVoltage"] = 12;
            buckOperatingPointJson["outputCurrent"] = 3;
            buckOperatingPointJson["switchingFrequency"] = 100000;
            buckOperatingPointJson["ambientTemperature"] = 42;
            buckInputsJson["operatingPoints"].push_back(buckOperatingPointJson);
        }

        OpenMagnetics::Buck buckInputs(buckInputsJson);

        auto inputs = buckInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Buck_Primary_Minimum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Buck_Primary_Voltage_Minimum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Buck_Primary_Maximum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Buck_Primary_Voltage_Maximum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        CHECK_CLOSE(double(buckInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(buckInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(buckInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(buckInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);
    }
    TEST(Test_Buck_Drain_Source_Voltage_BMO) {
        json buckInputsJson;
        json inputVoltage;

        inputVoltage["minimum"] = 20;
        inputVoltage["maximum"] = 240;
        buckInputsJson["inputVoltage"] = inputVoltage;
        buckInputsJson["diodeVoltageDrop"] = 0.7;
        buckInputsJson["efficiency"] = 0.9;
        buckInputsJson["maximumSwitchCurrent"] = 8;
        buckInputsJson["operatingPoints"] = json::array();
        {
            json buckOperatingPointJson;
            buckOperatingPointJson["outputVoltage"] = 12;
            buckOperatingPointJson["outputCurrent"] = 3;
            buckOperatingPointJson["switchingFrequency"] = 100000;
            buckOperatingPointJson["ambientTemperature"] = 42;
            buckInputsJson["operatingPoints"].push_back(buckOperatingPointJson);
        }
        OpenMagnetics::Buck buckInputs(buckInputsJson);
        buckInputs._assertErrors = true;

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

        auto operatingPoints = buckInputs.process_operating_points(magnetic);

    }
    TEST(Test_Buck_Web_0) {
        json buckInputsJson = json::parse(R"({"inputVoltage":{"minimum":10,"maximum":12},"diodeVoltageDrop":0.7,"efficiency":0.85,"currentRippleRatio":0.4,"operatingPoints":[{"outputVoltage":5,"outputCurrent":2,"switchingFrequency":100000,"ambientTemperature":25}]})");
        // json buckInputsJson;
        // json inputVoltage;

        // inputVoltage["minimum"] = 10;
        // inputVoltage["maximum"] = 12;
        // buckInputsJson["inputVoltage"] = inputVoltage;
        // buckInputsJson["diodeVoltageDrop"] = 0.7;
        // buckInputsJson["efficiency"] = 0.85;
        // buckInputsJson["maximumSwitchCurrent"] = 8;
        // buckInputsJson["operatingPoints"] = json::array();
        // {
        //     json buckOperatingPointJson;
        //     buckOperatingPointJson["outputVoltage"] = 5;
        //     buckOperatingPointJson["outputCurrent"] = 2;
        //     buckOperatingPointJson["switchingFrequency"] = 100000;
        //     buckOperatingPointJson["ambientTemperature"] = 42;
        //     buckInputsJson["operatingPoints"].push_back(buckOperatingPointJson);
        // }

        OpenMagnetics::Buck buckInputs(buckInputsJson);

        auto inputs = buckInputs.process();
        auto currentProcessed = inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed().value();

        std::cout << "currentProcessed.get_duty_cycle().value()" << std::endl;
        std::cout << currentProcessed.get_duty_cycle().value() << std::endl;
        std::cout << currentProcessed.get_duty_cycle().value() << std::endl;
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Buck_Primary_Minimum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Buck_Primary_Voltage_Minimum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Buck_Primary_Maximum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Buck_Primary_Voltage_Maximum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        CHECK_CLOSE(double(buckInputsJson["inputVoltage"]["minimum"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(buckInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        CHECK_CLOSE(double(buckInputsJson["inputVoltage"]["maximum"]), inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(buckInputsJson["inputVoltage"]["maximum"]) * maximumError);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR_WITH_DEADTIME);
        CHECK(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);
    }
}