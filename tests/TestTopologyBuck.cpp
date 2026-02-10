#include <source_location>
#include "support/Painter.h"
#include "converter_models/Buck.h"
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
    [[maybe_unused]] double maximumError = 0.1;

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

        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_Buck_Ngspice_Simulation", "[converter-model][buck-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create a Buck converter specification
        OpenMagnetics::Buck buck;
        
        // Input voltage: 24V nominal (18-32V range)
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(24.0);
        inputVoltage.set_minimum(18.0);
        inputVoltage.set_maximum(32.0);
        buck.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        buck.set_diode_voltage_drop(0.5);
        
        // Efficiency
        buck.set_efficiency(0.95);
        
        // Operating point: 5V @ 2A output, 100kHz
        BuckOperatingPoint opPoint;
        opPoint.set_output_voltage(5.0);
        opPoint.set_output_current(2.0);
        opPoint.set_switching_frequency(100e3);
        opPoint.set_ambient_temperature(25.0);
        
        buck.set_operating_points({opPoint});
        buck.set_current_ripple_ratio(0.4);
        
        // Process design requirements to get inductance
        auto designReqs = buck.process_design_requirements();
        double inductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Buck - Inductance: " << (inductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto operatingPoints = buck.simulate_and_extract_topology_waveforms(inductance);
        
        REQUIRE(!operatingPoints.empty());
        
        // Verify we have excitations
        REQUIRE(!operatingPoints[0].get_excitations_per_winding().empty());
        
        // Get primary (inductor) excitation
        const auto& primaryExc = operatingPoints[0].get_excitations_per_winding()[0];
        REQUIRE(primaryExc.get_voltage().has_value());
        REQUIRE(primaryExc.get_current().has_value());
        
        // Extract waveform data
        auto voltageData = primaryExc.get_voltage()->get_waveform()->get_data();
        auto currentData = primaryExc.get_current()->get_waveform()->get_data();
        
        // Calculate statistics
        double v_max = *std::max_element(voltageData.begin(), voltageData.end());
        double v_min = *std::min_element(voltageData.begin(), voltageData.end());
        double i_max = *std::max_element(currentData.begin(), currentData.end());
        double i_min = *std::min_element(currentData.begin(), currentData.end());
        double i_avg = std::accumulate(currentData.begin(), currentData.end(), 0.0) / currentData.size();
        
        INFO("Inductor voltage max: " << v_max << " V");
        INFO("Inductor voltage min: " << v_min << " V");
        INFO("Inductor current max: " << i_max << " A");
        INFO("Inductor current min: " << i_min << " A");
        INFO("Inductor current avg: " << i_avg << " A");
        
        // For Buck, inductor voltage swings between (Vin - Vout) and -Vout
        // Vin = 24V, Vout = 5V, so voltage should be around 19V and -5V
        CHECK(v_max > 15.0);  // Should be around 19V during switch ON
        CHECK(v_max < 25.0);
        
        CHECK(v_min < 0.0);  // Should be negative during switch OFF
        CHECK(v_min > -10.0);
        
        // Average inductor current should be close to output current
        CHECK(i_avg > 1.5);  // Should be around 2A
        CHECK(i_avg < 2.5);
        
        // In CCM, current should not go to zero
        CHECK(i_min > 0.0);
        
        INFO("Buck ngspice simulation test passed");
    }
}  // namespace
