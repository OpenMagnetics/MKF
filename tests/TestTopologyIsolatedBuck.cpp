#include "support/Painter.h"
#include "converter_models/IsolatedBuck.h"
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

    TEST_CASE("Test_IsolatedBuck", "[converter-model][isolated-buck-topology]") {
        json isolatedbuckInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 72;
        isolatedbuckInputsJson["inputVoltage"] = inputVoltage;
        isolatedbuckInputsJson["diodeVoltageDrop"] = 0.7;
        isolatedbuckInputsJson["maximumSwitchCurrent"] = 0.7;
        isolatedbuckInputsJson["efficiency"] = 1;
        isolatedbuckInputsJson["operatingPoints"] = json::array();

        {
            json isolatedbuckOperatingPointJson;
            isolatedbuckOperatingPointJson["outputVoltages"] = {10, 10};
            isolatedbuckOperatingPointJson["outputCurrents"] = {0.02, 0.1};
            isolatedbuckOperatingPointJson["switchingFrequency"] = 750000;
            isolatedbuckOperatingPointJson["ambientTemperature"] = 42;
            isolatedbuckInputsJson["operatingPoints"].push_back(isolatedbuckOperatingPointJson);
        }
        OpenMagnetics::IsolatedBuck isolatedbuckInputs(isolatedbuckInputsJson);
        isolatedbuckInputs._assertErrors = true;

        auto inputs = isolatedbuckInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuck_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuck_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuck_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuck_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        REQUIRE_THAT(10e-6, Catch::Matchers::WithinAbs(OpenMagnetics::resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance()), 10e-6 * 0.1));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding().size() == isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"].size());
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["minimum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(-inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_negative_peak().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["minimum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset(), 0.01));

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(-inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_negative_peak().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_IsolatedBuck_CurrentRippleRatio", "[converter-model][isolated-buck-topology]") {
        json isolatedbuckInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 72;
        isolatedbuckInputsJson["inputVoltage"] = inputVoltage;
        isolatedbuckInputsJson["diodeVoltageDrop"] = 0.7;
        isolatedbuckInputsJson["currentRippleRatio"] = 0.8;
        isolatedbuckInputsJson["efficiency"] = 1;
        isolatedbuckInputsJson["operatingPoints"] = json::array();

        {
            json isolatedbuckOperatingPointJson;
            isolatedbuckOperatingPointJson["outputVoltages"] = {10, 10};
            isolatedbuckOperatingPointJson["outputCurrents"] = {0.02, 0.1};
            isolatedbuckOperatingPointJson["switchingFrequency"] = 750000;
            isolatedbuckOperatingPointJson["ambientTemperature"] = 42;
            isolatedbuckInputsJson["operatingPoints"].push_back(isolatedbuckOperatingPointJson);
        }
        OpenMagnetics::IsolatedBuck isolatedbuckInputs(isolatedbuckInputsJson);
        isolatedbuckInputs._assertErrors = true;

        auto inputs = isolatedbuckInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuck_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuck_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuck_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuck_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        REQUIRE_THAT(110e-6, Catch::Matchers::WithinAbs(OpenMagnetics::resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance()), 110e-6 * 0.1));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding().size() == isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"].size());
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["minimum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(-inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_negative_peak().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["minimum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(-inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_negative_peak().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

}  // namespace
