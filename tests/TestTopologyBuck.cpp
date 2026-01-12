#include <source_location>
#include "support/Painter.h"
#include "converter_models/Buck.h"
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

    TEST_CASE("Test_Buck", "[converter-model][buck-topology][smoke-test]") {
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

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);
    }

    TEST_CASE("Test_Buck_Drain_Source_Voltage_BMO", "[converter-model][buck-topology][smoke-test]") {
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

    TEST_CASE("Test_Buck_Web_0", "[converter-model][buck-topology][smoke-test]") {
        json buckInputsJson = json::parse(R"({"inputVoltage":{"minimum":10,"maximum":12},"diodeVoltageDrop":0.7,"efficiency":0.85,"currentRippleRatio":0.4,"operatingPoints":[{"outputVoltage":5,"outputCurrent":2,"switchingFrequency":100000,"ambientTemperature":25}]})");
        OpenMagnetics::Buck buckInputs(buckInputsJson);

        auto inputs = buckInputs.process();
        auto currentProcessed = inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed().value();

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

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);
    }
}  // namespace
