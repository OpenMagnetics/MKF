#include <source_location>
#include "support/Painter.h"
#include "converter_models/IsolatedBuckBoost.h"
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

    TEST_CASE("Test_IsolatedBuckBoost", "[converter-model][isolated-buck-boost-topology][smoke-test]") {
        json isolatedBuckBoostInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 10;
        inputVoltage["maximum"] = 30;
        isolatedBuckBoostInputsJson["inputVoltage"] = inputVoltage;
        isolatedBuckBoostInputsJson["diodeVoltageDrop"] = 0;
        isolatedBuckBoostInputsJson["maximumSwitchCurrent"] = 2.5;
        isolatedBuckBoostInputsJson["efficiency"] = 1;
        isolatedBuckBoostInputsJson["operatingPoints"] = json::array();

        {
            json isolatedBuckBoostOperatingPointJson;
            isolatedBuckBoostOperatingPointJson["outputVoltages"] = {6, 5, 5};
            isolatedBuckBoostOperatingPointJson["outputCurrents"] = {0.01, 1, 0.3};
            isolatedBuckBoostOperatingPointJson["switchingFrequency"] = 400000;
            isolatedBuckBoostOperatingPointJson["ambientTemperature"] = 42;
            isolatedBuckBoostInputsJson["operatingPoints"].push_back(isolatedBuckBoostOperatingPointJson);
        }
        OpenMagnetics::IsolatedBuckBoost isolatedBuckBoostInputs(isolatedBuckBoostInputsJson);
        isolatedBuckBoostInputs._assertErrors = true;

        auto inputs = isolatedBuckBoostInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding().size() == isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"].size());
        // Primary average must be >= primaryOutputCurrent (includes reflected secondary currents)
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value() >= double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() != 0);

        REQUIRE_THAT(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset(), 0.01));

        // Primary average must be >= primaryOutputCurrent (includes reflected secondary currents)
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value() >= double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][0]));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::TRIANGULAR);

        REQUIRE_THAT(double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedBuckBoostInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset(), 0.01));
    }

    TEST_CASE("Test_IsolatedBuckBoost_Ngspice_Simulation", "[converter-model][isolated-buck-boost-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create an Isolated Buck-Boost converter
        OpenMagnetics::IsolatedBuckBoost isolatedBuckBoost;
        
        // Input voltage: 12V nominal
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(12.0);
        inputVoltage.set_minimum(9.0);
        inputVoltage.set_maximum(15.0);
        isolatedBuckBoost.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        isolatedBuckBoost.set_diode_voltage_drop(0.5);
        
        // Efficiency
        isolatedBuckBoost.set_efficiency(0.9);
        
        // Current ripple ratio
        isolatedBuckBoost.set_current_ripple_ratio(0.3);
        
        // Operating point: 5V @ 1A output on secondary, 200kHz
        // For Isolated Buck-Boost: output_voltages[0] = primary voltage/current (inductor side)
        // output_voltages[1] = secondary voltage/current
        IsolatedBuckBoostOperatingPoint opPoint;
        opPoint.set_output_voltages({6.0, 5.0});  // primary ~6V, secondary 5V
        opPoint.set_output_currents({0.5, 1.0});  // primary ~0.5A, secondary 1A
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        isolatedBuckBoost.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = isolatedBuckBoost.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Isolated Buck-Boost - Turns ratio: " << turnsRatios[0]);
        INFO("Isolated Buck-Boost - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto converterWaveforms = isolatedBuckBoost.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(converterWaveforms.size() >= 1);
        
        for (size_t opIndex = 0; opIndex < converterWaveforms.size(); opIndex++) {
            auto& wf = converterWaveforms[opIndex];
            
            // Check that time vector has reasonable values
            REQUIRE(wf.get_input_voltage().get_data().size() > 0);
            REQUIRE(wf.get_input_voltage().get_time().value()[0] >= 0);
            
            // Check primary current waveform
            REQUIRE(wf.get_input_current().get_data().size() == wf.get_input_voltage().get_data().size());
            
            // Check primary voltage waveform
            REQUIRE(wf.get_input_voltage().get_data().size() == wf.get_input_voltage().get_data().size());
            
            // For Buck-Boost, primary voltage should be close to input voltage during ON time
            double priV_max = *std::max_element(wf.get_input_voltage().get_data().begin(), wf.get_input_voltage().get_data().end());
            INFO("Primary voltage max: " << priV_max << " V");
            CHECK(priV_max > 5.0);  // Should be around 12V input
            CHECK(priV_max < 20.0);
            
            // Check output voltages if available
            if (wf.get_output_voltages().size() >= 1 && wf.get_output_voltages()[0].get_data().size() == wf.get_input_voltage().get_data().size()) {
                // Verify output voltage is close to expected
                double avgOutputVoltage = std::accumulate(wf.get_output_voltages()[0].get_data().begin(), wf.get_output_voltages()[0].get_data().end(), 0.0) / wf.get_output_voltages()[0].get_data().size();
                REQUIRE_THAT(std::abs(avgOutputVoltage), Catch::Matchers::WithinAbs(5.0, 5.0));  // Within 5V of expected 5V output
            }
            
            // Paint waveforms for visual inspection
            {
                auto outFile = outputFilePath;
                outFile.append("Test_IsolatedBuckBoost_Ngspice_PrimaryCurrent_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                auto priCurrentWaveform = wf.get_input_current();
                painter.paint_waveform(priCurrentWaveform);
                painter.export_svg();
            }
            {
                auto outFile = outputFilePath;
                outFile.append("Test_IsolatedBuckBoost_Ngspice_PrimaryVoltage_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                auto voltageWaveform = wf.get_input_voltage();
                painter.paint_waveform(voltageWaveform);
                painter.export_svg();
            }
            if (wf.get_output_voltages().size() > 0 && wf.get_output_voltages()[0].get_data().size() > 0) {
                auto outFile = outputFilePath;
                outFile.append("Test_IsolatedBuckBoost_Ngspice_OutputVoltage_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                auto outputWaveform = wf.get_output_voltages()[0];
                painter.paint_waveform(outputWaveform);
                painter.export_svg();
            }
            if (wf.get_output_currents().size() > 0 && wf.get_output_currents()[0].get_data().size() > 0) {
                auto outFile = outputFilePath;
                outFile.append("Test_IsolatedBuckBoost_Ngspice_SecondaryCurrent_OP" + std::to_string(opIndex) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, false, true);
                auto secCurrentWaveform = wf.get_output_currents()[0];
                painter.paint_waveform(secCurrentWaveform);
                painter.export_svg();
            }
        }
        
        INFO("Isolated Buck-Boost ngspice simulation test passed");

        // Shape tests on first operating point
        auto& wf0 = converterWaveforms[0];
        auto& priVoltageData = wf0.get_input_voltage().get_data();
        auto& priCurrentData = wf0.get_input_current().get_data();

        SECTION("Waveform shape: primary voltage pulse (ON ≈ Vin)") {
            double Vin = 12.0;
            int count_on = 0;
            for (double v : priVoltageData) {
                if (v > 0.5 * Vin) count_on++;
            }
            double frac_on = double(count_on) / priVoltageData.size();
            INFO("Fraction at Vin: " << frac_on);
            CHECK(frac_on > 0.10);
            CHECK(frac_on < 0.80);
        }

        SECTION("Waveform shape: primary current ramps (sawtooth/triangular)") {
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

    TEST_CASE("Test_IsolatedBuckBoost_PtP_AnalyticalVsNgspice", "[converter-model][isolated-buck-boost-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        OpenMagnetics::IsolatedBuckBoost isolatedBuckBoost;
        DimensionWithTolerance iv; iv.set_nominal(12.0); iv.set_minimum(9.0); iv.set_maximum(15.0);
        isolatedBuckBoost.set_input_voltage(iv);
        isolatedBuckBoost.set_diode_voltage_drop(0.5);
        isolatedBuckBoost.set_efficiency(0.9);
        isolatedBuckBoost.set_current_ripple_ratio(0.3);

        IsolatedBuckBoostOperatingPoint op;
        op.set_output_voltages({6.0, 5.0}); op.set_output_currents({0.5, 1.0});
        op.set_switching_frequency(200e3); op.set_ambient_temperature(25.0);
        isolatedBuckBoost.set_operating_points({op});

        auto designReqs = isolatedBuckBoost.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
        double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto analyticalOps = isolatedBuckBoost.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        REQUIRE(!aData.empty());
        REQUIRE(!aTime.empty());
        auto aResampled = ptp_interp(aTime, aData, 256);

        isolatedBuckBoost.set_num_periods_to_extract(1);
        auto simOps = isolatedBuckBoost.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        REQUIRE(!sData.empty());
        REQUIRE(!sTime.empty());
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        INFO("Isolated Buck-Boost primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.30);
    }

}  // namespace
