#include <source_location>
#include "support/Painter.h"
#include "converter_models/IsolatedBuck.h"
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
#include <cmath>
#include <limits>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    // ── PtP waveform comparison helpers ────────────────────────────────────
    std::vector<double> ptp_interp(const std::vector<double>& t, const std::vector<double>& d, int N) {
        std::vector<double> out(N);
        double T = t.back();
        for (int i = 0; i < N; ++i) {
            double ti = T * i / (N - 1);
            int lo = 0;
            for (int k = 0; k + 1 < (int)t.size(); ++k) if (t[k] <= ti) lo = k;
            int hi = std::min(lo + 1, (int)t.size() - 1);
            double dt = t[hi] - t[lo];
            out[i] = (dt < 1e-20) ? d[hi] : d[lo] + (ti - t[lo]) / dt * (d[hi] - d[lo]);
        }
        return out;
    }

    double ptp_nrmse(const std::vector<double>& ref, const std::vector<double>& sim) {
        // Mean-subtracted, scale-normalized NRMSE with circular phase alignment.
        // Measures shape similarity only — invariant to DC offset and amplitude scale.
        int N = (int)ref.size();
        double ref_mean = 0.0, sim_mean = 0.0;
        for (int i = 0; i < N; ++i) { ref_mean += ref[i]; sim_mean += sim[i]; }
        ref_mean /= N; sim_mean /= N;
        std::vector<double> r(N), s(N);
        double rAC = 0.0, sAC = 0.0;
        for (int i = 0; i < N; ++i) {
            r[i] = ref[i] - ref_mean; s[i] = sim[i] - sim_mean;
            rAC += r[i] * r[i]; sAC += s[i] * s[i];
        }
        rAC = std::sqrt(rAC / N); sAC = std::sqrt(sAC / N);
        if (rAC < 1e-10 || sAC < 1e-10) return 1.0;
        for (int i = 0; i < N; ++i) { r[i] /= rAC; s[i] /= sAC; }
        int ns = std::min(N, 64);
        double best = std::numeric_limits<double>::max();
        for (int ss = 0; ss < ns; ++ss) {
            int sh = ss * N / ns;
            double ssd = 0.0;
            for (int k = 0; k < N; ++k) { double e = r[k] - s[(k + sh) % N]; ssd += e * e; }
            if (ssd < best) best = ssd;
        }
        return std::sqrt(best / N);
    }

    std::vector<double> ptp_current(const OperatingPoint& op, size_t wi = 0) {
        if (wi >= op.get_excitations_per_winding().size()) return {};
        auto& e = op.get_excitations_per_winding()[wi];
        if (!e.get_current() || !e.get_current()->get_waveform()) return {};
        return e.get_current()->get_waveform()->get_data();
    }

    std::vector<double> ptp_current_time(const OperatingPoint& op, size_t wi = 0) {
        if (wi >= op.get_excitations_per_winding().size()) return {};
        auto& e = op.get_excitations_per_winding()[wi];
        if (!e.get_current() || !e.get_current()->get_waveform()) return {};
        auto tv = e.get_current()->get_waveform()->get_time();
        return tv ? tv.value() : std::vector<double>{};
    }

    TEST_CASE("Test_IsolatedBuck", "[converter-model][isolated-buck-topology][smoke-test]") {
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

    TEST_CASE("Test_IsolatedBuck_CurrentRippleRatio", "[converter-model][isolated-buck-topology][smoke-test]") {
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

    TEST_CASE("Test_IsolatedBuck_Ngspice_Simulation", "[converter-model][isolated-buck-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create an Isolated Buck (Forward with synchronous rectification) converter
        OpenMagnetics::IsolatedBuck isolatedBuck;
        
        // Input voltage
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(36.0);
        inputVoltage.set_maximum(72.0);
        isolatedBuck.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        isolatedBuck.set_diode_voltage_drop(0.5);
        
        // Efficiency
        isolatedBuck.set_efficiency(0.9);
        
        // Current ripple ratio
        isolatedBuck.set_current_ripple_ratio(0.3);
        
        // Operating point: 5V @ 5A output on secondary, 200kHz
        // For Isolated Buck: output_voltages[0] = primary voltage, output_voltages[1] = secondary voltage
        // output_currents[0] = primary current, output_currents[1] = secondary current
        IsolatedBuckOperatingPoint opPoint;
        opPoint.set_output_voltages({5.0, 5.0});  // primary voltage ~5V, secondary voltage 5V
        opPoint.set_output_currents({0.6, 5.0});  // primary current ~0.6A (reflected from secondary), secondary current 5A
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        isolatedBuck.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = isolatedBuck.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Isolated Buck - Turns ratio: " << turnsRatios[0]);
        INFO("Isolated Buck - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto operatingPoints = isolatedBuck.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!operatingPoints.empty());
        
        // Verify we have waveform data
        REQUIRE(!operatingPoints[0].get_input_voltage().get_data().empty());
        REQUIRE(!operatingPoints[0].get_input_current().get_data().empty());

        // Extract waveform data
        auto priVoltageData = operatingPoints[0].get_input_voltage().get_data();
        auto priCurrentData = operatingPoints[0].get_input_current().get_data();
        
        // Calculate statistics
        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priI_avg = std::accumulate(priCurrentData.begin(), priCurrentData.end(), 0.0) / priCurrentData.size();
        
        INFO("Primary voltage max: " << priV_max << " V");
        INFO("Primary current avg: " << priI_avg << " A");
        
        // Validate primary voltage: should be close to input voltage during ON time
        CHECK(priV_max > 30.0);  // Should be around 48V
        CHECK(priV_max < 80.0);

        INFO("Isolated Buck ngspice simulation test passed");

        SECTION("Waveform shape: primary voltage is a pulse (ON ≈ Vin, OFF ≈ 0)") {
            double Vin = 48.0;
            int count_on = 0;
            for (double v : priVoltageData) {
                if (v > 0.5 * Vin) count_on++;
            }
            double frac_on = double(count_on) / priVoltageData.size();
            INFO("Fraction at Vin: " << frac_on << " (expect ≈ D ≈ 0.2-0.5)");
            CHECK(frac_on > 0.10);
            CHECK(frac_on < 0.70);
        }

        SECTION("Waveform shape: primary current ramps up during ON time") {
            auto peak_it = std::max_element(priCurrentData.begin(), priCurrentData.end());
            int peak_idx = (int)std::distance(priCurrentData.begin(), peak_it);
            int N = (int)priCurrentData.size();

            double rising_sum = 0.0;
            for (int k = 1; k <= peak_idx; ++k)
                rising_sum += priCurrentData[k] - priCurrentData[k-1];

            double falling_sum = 0.0;
            for (int k = peak_idx + 1; k < N; ++k)
                falling_sum += priCurrentData[k] - priCurrentData[k-1];

            INFO("Rising sum: " << rising_sum << " A, falling sum: " << falling_sum << " A");
            CHECK(rising_sum > 0.0);
            CHECK(falling_sum < 0.0);
        }
    }

    TEST_CASE("Test_IsolatedBuck_Ngspice_TwoSecondaries_Validation", "[converter-model][isolated-buck-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create an Isolated Buck (Flybuck) converter with TWO secondaries
        // First secondary: non-isolated buck output
        // Second secondary: isolated flyback output
        OpenMagnetics::IsolatedBuck isolatedBuck;
        
        // Input voltage: 48V nominal
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        isolatedBuck.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        isolatedBuck.set_diode_voltage_drop(0.5);
        
        // Efficiency
        isolatedBuck.set_efficiency(0.9);
        
        // Current ripple ratio
        isolatedBuck.set_current_ripple_ratio(0.3);
        
        // Operating point with TWO outputs:
        // output_voltages[0] = primary (non-isolated buck) output voltage
        // output_voltages[1] = secondary 1 (non-isolated, tied to primary ground) voltage  
        // output_voltages[2] = secondary 2 (isolated) voltage
        // output_currents[0] = primary output current
        // output_currents[1] = secondary 1 current
        // output_currents[2] = secondary 2 current
        IsolatedBuckOperatingPoint opPoint;
        double expectedPriVoltage = 5.0;    // Non-isolated buck output: 5V
        double expectedSec1Voltage = 5.0;   // Secondary 1: 5V (tied to same ground)
        double expectedSec2Voltage = 12.0;  // Secondary 2 (isolated): 12V
        double expectedPriCurrent = 1.0;    // Primary output current: 1A
        double expectedSec1Current = 2.0;   // Secondary 1 current: 2A
        double expectedSec2Current = 0.5;   // Secondary 2 current: 0.5A
        
        opPoint.set_output_voltages({expectedPriVoltage, expectedSec1Voltage, expectedSec2Voltage});
        opPoint.set_output_currents({expectedPriCurrent, expectedSec1Current, expectedSec2Current});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        isolatedBuck.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = isolatedBuck.process_design_requirements();
        
        std::vector<double> turnsRatios;
        // Debug: Check what's in designReqs
        REQUIRE(designReqs.get_turns_ratios().size() == 2);
        for (size_t i = 0; i < designReqs.get_turns_ratios().size(); ++i) {
            const auto& tr = designReqs.get_turns_ratios()[i];
            // Try to get the value - should have nominal set
            double val = 1.0;
            if (tr.get_nominal()) {
                val = tr.get_nominal().value();
            } else if (tr.get_minimum()) {
                val = tr.get_minimum().value();
            }
            turnsRatios.push_back(val);
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        // Verify turns ratios were extracted correctly
        REQUIRE(turnsRatios.size() == 2);
        
        INFO("Flybuck Two Secondaries Test");
        INFO("Turns ratios: TR[0]=" << turnsRatios[0] << ", TR[1]=" << turnsRatios[1]);
        INFO("Input voltage: 48V");
        INFO("Expected outputs:");
        INFO("  Primary (non-isolated): " << expectedPriVoltage << "V @ " << expectedPriCurrent << "A");
        INFO("  Secondary 1 (non-isolated): " << expectedSec1Voltage << "V @ " << expectedSec1Current << "A");
        INFO("  Secondary 2 (isolated): " << expectedSec2Voltage << "V @ " << expectedSec2Current << "A");
        INFO("Number of turns ratios from designReqs: " << designReqs.get_turns_ratios().size());
        INFO("Turns ratios:");
        for (size_t i = 0; i < turnsRatios.size(); ++i) {
            INFO("  Secondary " << (i+1) << ": " << turnsRatios[i]);
        }
        INFO("Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Ensure we have turns ratios for the secondaries
        REQUIRE(turnsRatios.size() == 2);  // We expect 2 secondaries
        
        // Run ngspice simulation to extract operating points
        auto operatingPoints = isolatedBuck.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        
        REQUIRE(!operatingPoints.empty());
        REQUIRE(operatingPoints[0].get_excitations_per_winding().size() == 3);  // Primary + 2 secondaries
        
        // Helper lambda to get average and RMS from waveform
        auto getWaveformStats = [](const Waveform& wf) -> std::pair<double, double> {
            const auto& data = wf.get_data();
            if (data.empty()) return {0.0, 0.0};
            
            double avg = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
            double sumSq = 0.0;
            for (double v : data) {
                sumSq += v * v;
            }
            double rms = std::sqrt(sumSq / data.size());
            return {avg, rms};
        };
        
        // Validate primary winding (non-isolated buck output)
        // Winding 0: Primary - voltage should be switching waveform at pri_in node
        const auto& primaryExcitation = operatingPoints[0].get_excitations_per_winding()[0];
        REQUIRE(primaryExcitation.get_voltage().has_value());
        REQUIRE(primaryExcitation.get_current().has_value());
        
        auto [priVoltageAvg, priVoltageRms] = getWaveformStats(primaryExcitation.get_voltage()->get_waveform().value());
        auto [priCurrentAvg, priCurrentRms] = getWaveformStats(primaryExcitation.get_current()->get_waveform().value());
        
        INFO("Primary Winding (Winding 0):");
        INFO("  Voltage - Avg: " << priVoltageAvg << " V, RMS: " << priVoltageRms << " V");
        INFO("  Current - Avg: " << priCurrentAvg << " A, RMS: " << priCurrentRms << " A");
        
        // Primary current should not be thousands of amps (the bug we're fixing)
        CHECK(std::abs(priCurrentAvg) < 100.0);
        CHECK(std::abs(priCurrentRms) < 1000.0);
        
        // Primary current should be in reasonable range (couple of amps for this converter)
        CHECK(std::abs(priCurrentAvg) < 50.0);
        CHECK(std::abs(priCurrentRms) < 100.0);
        
        // Validate secondary 1 (non-isolated)
        const auto& sec1Excitation = operatingPoints[0].get_excitations_per_winding()[1];
        REQUIRE(sec1Excitation.get_voltage().has_value());
        REQUIRE(sec1Excitation.get_current().has_value());
        
        auto [sec1VoltageAvg, sec1VoltageRms] = getWaveformStats(sec1Excitation.get_voltage()->get_waveform().value());
        auto [sec1CurrentAvg, sec1CurrentRms] = getWaveformStats(sec1Excitation.get_current()->get_waveform().value());
        
        INFO("Secondary 1 (Winding 1 - non-isolated):");
        INFO("  Voltage - Avg: " << sec1VoltageAvg << " V, RMS: " << sec1VoltageRms << " V");
        INFO("  Current - Avg: " << sec1CurrentAvg << " A, RMS: " << sec1CurrentRms << " A");
        
        // Secondary 1 current should not be thousands of amps
        CHECK(std::abs(sec1CurrentAvg) < 100.0);
        CHECK(std::abs(sec1CurrentRms) < 1000.0);
        
        // Validate secondary 2 (isolated)
        const auto& sec2Excitation = operatingPoints[0].get_excitations_per_winding()[2];
        REQUIRE(sec2Excitation.get_voltage().has_value());
        REQUIRE(sec2Excitation.get_current().has_value());
        
        auto [sec2VoltageAvg, sec2VoltageRms] = getWaveformStats(sec2Excitation.get_voltage()->get_waveform().value());
        auto [sec2CurrentAvg, sec2CurrentRms] = getWaveformStats(sec2Excitation.get_current()->get_waveform().value());
        
        INFO("Secondary 2 (Winding 2 - isolated):");
        INFO("  Voltage - Avg: " << sec2VoltageAvg << " V, RMS: " << sec2VoltageRms << " V");
        INFO("  Current - Avg: " << sec2CurrentAvg << " A, RMS: " << sec2CurrentRms << " A");
        
        // Secondary 2 current should not be thousands of amps
        CHECK(std::abs(sec2CurrentAvg) < 100.0);
        CHECK(std::abs(sec2CurrentRms) < 1000.0);
        
        // Validate secondary voltages are approximately correct (within 20% of expected)
        // Primary output (non-isolated buck): ~5V
        CHECK(std::abs(priVoltageAvg - expectedPriVoltage) < expectedPriVoltage * 0.2);
        // Secondary 1: ~5V
        CHECK(std::abs(sec1VoltageAvg - expectedSec1Voltage) < expectedSec1Voltage * 0.2);
        // Secondary 2 (isolated): ~12V
        CHECK(std::abs(sec2VoltageAvg - expectedSec2Voltage) < expectedSec2Voltage * 0.2);
        
        // Now run topology waveforms simulation for additional validation
        auto topologyWaveforms = isolatedBuck.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!topologyWaveforms.empty());
        
        // Validate topology waveforms were extracted (input/output signals)
        REQUIRE(!topologyWaveforms[0].get_input_voltage().get_data().empty());
        REQUIRE(!topologyWaveforms[0].get_input_current().get_data().empty());
        
        // The key validation: currents should be in reasonable range (not thousands of amps)
        // This was the main bug we fixed
        auto inputCurrentData = topologyWaveforms[0].get_input_current().get_data();
        double inputI_avg = std::accumulate(inputCurrentData.begin(), inputCurrentData.end(), 0.0) / inputCurrentData.size();
        INFO("Input current average: " << inputI_avg << " A");
        CHECK(std::abs(inputI_avg) < 50.0);  // Should be in single-digit or low tens of amps
        CHECK(std::abs(inputI_avg) < 1000.0);  // Definitely not thousands of amps
        
        INFO("Flybuck Two Secondaries ngspice simulation test completed successfully");
    }

    TEST_CASE("Test_IsolatedBuck_PtP_AnalyticalVsNgspice", "[converter-model][isolated-buck-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        OpenMagnetics::IsolatedBuck isolatedBuck;
        DimensionWithTolerance iv; iv.set_nominal(48.0); iv.set_minimum(36.0); iv.set_maximum(72.0);
        isolatedBuck.set_input_voltage(iv);
        isolatedBuck.set_diode_voltage_drop(0.5);
        isolatedBuck.set_efficiency(0.9);
        isolatedBuck.set_current_ripple_ratio(0.3);

        IsolatedBuckOperatingPoint op;
        op.set_output_voltages({5.0, 5.0}); op.set_output_currents({0.6, 5.0});
        op.set_switching_frequency(200e3); op.set_ambient_temperature(25.0);
        isolatedBuck.set_operating_points({op});

        auto designReqs = isolatedBuck.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
        double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto analyticalOps = isolatedBuck.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        REQUIRE(!aData.empty());
        REQUIRE(!aTime.empty());
        auto aResampled = ptp_interp(aTime, aData, 256);

        isolatedBuck.set_num_periods_to_extract(1);
        auto simOps = isolatedBuck.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        REQUIRE(!sData.empty());
        REQUIRE(!sTime.empty());
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        INFO("Isolated Buck primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.35);
    }

}  // namespace
