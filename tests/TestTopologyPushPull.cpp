#include <source_location>
#include "support/Painter.h"
#include "converter_models/PushPull.h"
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

    TEST_CASE("Test_PushPull_CCM", "[converter-model][push-pull-topology][smoke-test]") {
        json pushPullInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 20;
        inputVoltage["maximum"] = 30;
        pushPullInputsJson["inputVoltage"] = inputVoltage;
        pushPullInputsJson["diodeVoltageDrop"] = 0.6;
        pushPullInputsJson["efficiency"] = 0.9;
        pushPullInputsJson["maximumSwitchCurrent"] = 1;
        pushPullInputsJson["currentRippleRatio"] = 0.3;
        pushPullInputsJson["dutyCycle"] = 0.45;
        pushPullInputsJson["operatingPoints"] = json::array();

        {
            json PushPullOperatingPointJson;
            PushPullOperatingPointJson["outputVoltages"] = {48, 5, 9};
            PushPullOperatingPointJson["outputCurrents"] = {0.7, 0.01, 0.01};
            PushPullOperatingPointJson["switchingFrequency"] = 500000;
            PushPullOperatingPointJson["ambientTemperature"] = 42;
            pushPullInputsJson["operatingPoints"].push_back(PushPullOperatingPointJson);
        }
        OpenMagnetics::PushPull pushPullInputs(pushPullInputsJson);
        pushPullInputs._assertErrors = true;

        auto inputs = pushPullInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Auxiliary_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[4].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Auxiliary_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[5].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_First_Auxiliary_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[4].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_CCM_Second_Auxiliary_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[5].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(pushPullInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(pushPullInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(pushPullInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(pushPullInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(pushPullInputsJson["operatingPoints"][0]["outputCurrents"][0]) / 2, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_average().value(), double(pushPullInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(pushPullInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_rms().value(), double(pushPullInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(pushPullInputsJson["operatingPoints"][0]["outputCurrents"][0]) / 2, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_average().value(), double(pushPullInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(pushPullInputsJson["operatingPoints"][0]["outputVoltages"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_processed()->get_rms().value(), double(pushPullInputsJson["operatingPoints"][0]["outputVoltages"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_offset() > 0);

    }

    TEST_CASE("Test_PushPull_DCM", "[converter-model][push-pull-topology][smoke-test]") {
        json pushPullInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 500;
        inputVoltage["maximum"] = 501;
        pushPullInputsJson["inputVoltage"] = inputVoltage;
        pushPullInputsJson["diodeVoltageDrop"] = 0.6;
        pushPullInputsJson["efficiency"] = 0.9;
        pushPullInputsJson["maximumSwitchCurrent"] = 1;
        pushPullInputsJson["currentRippleRatio"] = 0.3;
        pushPullInputsJson["desiredInductance"] = 1e-3;
        pushPullInputsJson["desiredOutputInductance"] = 5e-6;
        pushPullInputsJson["desiredTurnsRatios"] = {15};
        pushPullInputsJson["dutyCycle"] = 0.45;
        pushPullInputsJson["operatingPoints"] = json::array();

        {
            json PushPullOperatingPointJson;
            PushPullOperatingPointJson["outputVoltages"] = {15};
            PushPullOperatingPointJson["outputCurrents"] = {0.8};
            PushPullOperatingPointJson["switchingFrequency"] = 500000;
            PushPullOperatingPointJson["ambientTemperature"] = 42;
            pushPullInputsJson["operatingPoints"].push_back(PushPullOperatingPointJson);
        }
        OpenMagnetics::AdvancedPushPull pushPullInputs(pushPullInputsJson);
        pushPullInputs._assertErrors = true;

        auto inputs = pushPullInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_First_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_Second_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_First_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_Second_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_First_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_Second_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_First_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_PushPull_DCM_Second_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(pushPullInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(pushPullInputsJson["inputVoltage"]["maximum"]) * maximumError));

        REQUIRE_THAT(double(pushPullInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(pushPullInputsJson["inputVoltage"]["maximum"]) * maximumError));

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[3].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_PushPull_Ngspice_Simulation", "[converter-model][push-pull-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create a Push-Pull converter specification
        OpenMagnetics::PushPull pushPull;
        
        // Input voltage
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(24.0);
        inputVoltage.set_minimum(18.0);
        inputVoltage.set_maximum(32.0);
        pushPull.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        pushPull.set_diode_voltage_drop(0.5);
        
        // Efficiency
        pushPull.set_efficiency(0.9);
        
        // Current ripple ratio
        pushPull.set_current_ripple_ratio(0.3);
        
        // Operating point: 12V @ 5A output, 100kHz
        PushPullOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({5.0});
        opPoint.set_switching_frequency(100000.0);
        opPoint.set_ambient_temperature(25.0);
        pushPull.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = pushPull.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Push-Pull - Turns ratio: " << turnsRatios[0]);
        INFO("Push-Pull - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto operatingPoints = pushPull.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!operatingPoints.empty());
        
        // Verify we have excitations (primary center-tapped + secondary)
        REQUIRE(operatingPoints[0].get_excitations_per_winding().size() >= 3);
        
        // Get primary excitation
        const auto& primaryExc = operatingPoints[0].get_excitations_per_winding()[0];
        REQUIRE(primaryExc.get_voltage().has_value());
        REQUIRE(primaryExc.get_current().has_value());
        
        // Extract waveform data
        auto priVoltageData = primaryExc.get_voltage()->get_waveform()->get_data();
        
        // Calculate statistics
        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priV_min = *std::min_element(priVoltageData.begin(), priVoltageData.end());
        
        INFO("Primary voltage max: " << priV_max << " V");
        INFO("Primary voltage min: " << priV_min << " V");
        
        // For push-pull, primary voltage swings around input voltage
        // Note: Simulation may show voltage spikes due to leakage inductance
        CHECK(priV_max > 15.0);  // Should be around input voltage (24V)
        CHECK(priV_max < 100.0);
        CHECK(priV_min < -15.0);  // Should have negative swing
        CHECK(priV_min > -700.0);  // Allow for voltage spikes
        
        INFO("Push-Pull ngspice simulation test passed");
    }

// End of SUITE

}  // namespace
