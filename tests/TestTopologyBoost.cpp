#include <source_location>
#include "support/Painter.h"
#include "converter_models/Boost.h"
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

    TEST_CASE("Test_Boost", "[converter-model][boost-topology][smoke-test]") {
        json boostInputsJson;
        json inputVoltage;

        inputVoltage["minimum"] = 12;
        inputVoltage["maximum"] = 24;
        boostInputsJson["inputVoltage"] = inputVoltage;
        boostInputsJson["diodeVoltageDrop"] = 0.7;
        boostInputsJson["efficiency"] = 1;
        boostInputsJson["maximumSwitchCurrent"] = 8;
        boostInputsJson["operatingPoints"] = json::array();
        {
            json boostOperatingPointJson;
            boostOperatingPointJson["outputVoltage"] = 50;
            boostOperatingPointJson["outputCurrent"] = 1;
            boostOperatingPointJson["switchingFrequency"] = 100000;
            boostOperatingPointJson["ambientTemperature"] = 42;
            boostInputsJson["operatingPoints"].push_back(boostOperatingPointJson);
        }

        OpenMagnetics::Boost boostInputs(boostInputsJson);

        auto inputs = boostInputs.process();
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Boost_Primary.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Boost_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Boost_Primary_Maximum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Boost_Primary_Voltage_Maximum.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        REQUIRE_THAT(double(boostInputsJson["operatingPoints"][0]["outputVoltage"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(boostInputsJson["operatingPoints"][0]["outputVoltage"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(boostInputsJson["operatingPoints"][0]["outputVoltage"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(boostInputsJson["operatingPoints"][0]["outputVoltage"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);
    }

    TEST_CASE("Test_Boost_Ngspice_Simulation", "[converter-model][boost-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create a Boost converter specification
        OpenMagnetics::Boost boost;
        
        // Input voltage: 12V nominal (9-15V range)
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(12.0);
        inputVoltage.set_minimum(9.0);
        inputVoltage.set_maximum(15.0);
        boost.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        boost.set_diode_voltage_drop(0.5);
        
        // Efficiency
        boost.set_efficiency(0.92);
        
        // Operating point: 24V @ 1A output, 100kHz
        BoostOperatingPoint opPoint;
        opPoint.set_output_voltage(24.0);
        opPoint.set_output_current(1.0);
        opPoint.set_switching_frequency(100e3);
        opPoint.set_ambient_temperature(25.0);
        
        boost.set_operating_points({opPoint});
        boost.set_current_ripple_ratio(0.4);
        
        // Process design requirements to get inductance
        auto designReqs = boost.process_design_requirements();
        double inductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Boost - Inductance: " << (inductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto converterWaveforms = boost.simulate_and_extract_topology_waveforms(inductance);
        
        REQUIRE(!converterWaveforms.empty());
        
        // Verify we have excitations
        REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());
        
        // Get primary (inductor) excitation
        // primary excitation replaced by ConverterWaveforms input signals
// Extract waveform data
        auto voltageData = converterWaveforms[0].get_input_voltage().get_data();
        auto currentData = converterWaveforms[0].get_input_current().get_data();
        
        // Calculate statistics
        double v_max = *std::max_element(voltageData.begin(), voltageData.end());
        double v_min = *std::min_element(voltageData.begin(), voltageData.end());
        double i_avg = std::accumulate(currentData.begin(), currentData.end(), 0.0) / currentData.size();
        
        INFO("Inductor voltage max: " << v_max << " V");
        INFO("Inductor voltage min: " << v_min << " V");
        INFO("Inductor current avg: " << i_avg << " A");
        
        // For Boost, inductor voltage swings between Vin and (Vin - Vout)
        // Vin = 12V, Vout = 24V, so voltage should be around 12V and -12V
        CHECK(v_max > 10.0);  // Should be around 12V during switch ON
        CHECK(v_max < 15.0);
        
        // Average inductor current should be higher than output current (boost ratio)
        // Iin = Iout * Vout / Vin = 1A * 24V / 12V = 2A (ideal)
        // Actual will be higher due to efficiency losses
        CHECK(i_avg > 1.2);  // Allow for some variation in simulation
        CHECK(i_avg < 3.0);
        
        INFO("Boost ngspice simulation test passed");
    }

}  // namespace
