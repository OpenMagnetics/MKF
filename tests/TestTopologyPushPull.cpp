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
    // ───────────────────────────────────────────────────────────────────────

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
        auto converterWaveforms = pushPull.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!converterWaveforms.empty());
        
        // Verify we have excitations (primary center-tapped + secondary)
        // ConverterWaveforms doesn't have excitations_per_winding - check input/output instead
    REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());
        
        // Get primary excitation
        // Primary excitation data is now directly in ConverterWaveforms
        // Voltage data is directly available in ConverterWaveforms
        // Current data is directly available in ConverterWaveforms
        
        // Extract waveform data
        auto priVoltageData = converterWaveforms[0].get_input_voltage().get_data();
        
        // Calculate statistics
        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priV_min = *std::min_element(priVoltageData.begin(), priVoltageData.end());
        
        INFO("Primary voltage max: " << priV_max << " V");
        INFO("Primary voltage min: " << priV_min << " V");
        
        // For push-pull, pri_top is a single-ended node (top half of center-tapped primary)
        // When S1 is ON: pri_top pulled low (~0V)
        // When S2 is ON: reflected voltage drives pri_top high (~2*Vin = 48V)
        // Node is never negative — it swings 0V to ~2*Vin
        CHECK(priV_max > 30.0);  // Should reach ~2*Vin = 48V during opposite switch ON
        CHECK(priV_max < 100.0);
        CHECK(priV_min > -2.0);   // Single-ended node: minimum is ~0V

        INFO("Push-Pull ngspice simulation test passed");

        SECTION("Waveform shape: bimodal primary voltage (alternates between 0 and ~2*Vin)") {
            // pri_top alternates: ~0V when S1 conducting, ~2*Vin=48V when S2 conducting
            double Vin = 24.0;
            int count_low = 0, count_high = 0;
            for (double v : priVoltageData) {
                if (v < Vin) count_low++;
                else count_high++;
            }
            double frac_low = double(count_low) / priVoltageData.size();
            double frac_high = double(count_high) / priVoltageData.size();
            INFO("Fraction near 0V: " << frac_low << ", fraction near 2*Vin: " << frac_high);
            // Both states must occupy significant time for a balanced push-pull (D ≈ 0.5)
            CHECK(frac_low > 0.25);
            CHECK(frac_high > 0.25);
        }

        SECTION("Waveform shape: volt-second balance (avg primary voltage ≈ Vin)") {
            // avg(v_pri_top) = 0*D + 2*Vin*(1-D) ≈ Vin for D ≈ 0.5 push-pull
            double priV_avg = std::accumulate(priVoltageData.begin(), priVoltageData.end(), 0.0) / priVoltageData.size();
            INFO("Primary voltage avg: " << priV_avg << " V (expect ≈ 24V = Vin for D=0.5)");
            CHECK(priV_avg > 10.0);
            CHECK(priV_avg < 35.0);
        }

        SECTION("Waveform shape: primary current is non-zero and positive (CCM operation)") {
            auto& priCurrentData = converterWaveforms[0].get_input_current().get_data();
            if (!priCurrentData.empty()) {
                double priI_avg = std::accumulate(priCurrentData.begin(), priCurrentData.end(), 0.0) / priCurrentData.size();
                double priI_min = *std::min_element(priCurrentData.begin(), priCurrentData.end());
                double priI_max = *std::max_element(priCurrentData.begin(), priCurrentData.end());
                INFO("Primary current: avg=" << priI_avg << " A, min=" << priI_min << " A, max=" << priI_max << " A");
                CHECK(priI_avg > 0.5);   // Non-trivial average current
                CHECK(priI_max > 1.0);   // Some peak current
            }
        }
    }

    TEST_CASE("Test_PushPull_Waveform_Polarity", "[converter-model][push-pull-topology][ngspice-simulation][smoke-test]") {
        // Verify Push-Pull converter has correct waveform polarity:
        // - Primary voltage (v(pri_top)) should alternate between positive (during ON) and negative (during opposite switch ON)
        // - Output voltage should be stable around target value
        // In Push-Pull, alternating switches drive the center-tapped primary
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::PushPull pushpull;

        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(24.0);
        inputVoltage.set_minimum(18.0);
        inputVoltage.set_maximum(32.0);
        pushpull.set_input_voltage(inputVoltage);

        pushpull.set_diode_voltage_drop(0.5);
        pushpull.set_efficiency(0.85);
        pushpull.set_current_ripple_ratio(0.3);

        PushPullOperatingPoint opPoint;
        opPoint.set_output_voltages({48.0});
        opPoint.set_output_currents({2.0});
        opPoint.set_switching_frequency(100e3);
        opPoint.set_ambient_temperature(25.0);
        pushpull.set_operating_points({opPoint});

        auto designReqs = pushpull.process_design_requirements();

        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        INFO("Turns ratios count: " << turnsRatios.size());
        INFO("Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");

        // Run simulation and extract operating points (winding-level waveforms)
        auto operatingPoints = pushpull.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!operatingPoints.empty());
        REQUIRE(operatingPoints[0].get_excitations_per_winding().size() >= 2);

        auto& primaryExcitation = operatingPoints[0].get_excitations_per_winding()[0];
        REQUIRE(primaryExcitation.get_voltage().has_value());
        REQUIRE(primaryExcitation.get_voltage()->get_waveform().has_value());

        auto primaryVoltageData = primaryExcitation.get_voltage()->get_waveform()->get_data();

        double priV_max = *std::max_element(primaryVoltageData.begin(), primaryVoltageData.end());
        double priV_min = *std::min_element(primaryVoltageData.begin(), primaryVoltageData.end());

        INFO("Primary voltage max: " << priV_max << " V, min: " << priV_min << " V");

        // Push-Pull primary: voltage should be positive during one half and negative during the other
        // v(pri_top) sees +Vin during S1 ON, and goes negative when S2 ON (via magnetic coupling)
        CHECK(priV_max > 15.0);   // Should be around Vin during S1 ON
        CHECK(priV_min < -5.0);   // Should go negative during S2 ON

        // Also verify converter-level waveforms
        auto converterWaveforms = pushpull.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        REQUIRE(!converterWaveforms.empty());

        auto& cwf = converterWaveforms[0];
        auto cwfInputVoltage = cwf.get_input_voltage().get_data();
        REQUIRE(!cwfInputVoltage.empty());

        double cwfV_max = *std::max_element(cwfInputVoltage.begin(), cwfInputVoltage.end());
        double cwfV_min = *std::min_element(cwfInputVoltage.begin(), cwfInputVoltage.end());

        INFO("Converter input voltage (pri_top) max: " << cwfV_max << " V, min: " << cwfV_min << " V");
        // Note: pri_top node voltage swings between ~0V (S1 ON) and ~Vin (S2 ON via transformer coupling)
        // It doesn't go significantly negative due to the snubber resistor to ground
        CHECK(cwfV_max > 15.0);
        CHECK(cwfV_min < 5.0);  // Allow near-zero minimum (was < -5.0, but pri_top doesn't go negative)

        // Output voltage should be around 48V (stable)
        REQUIRE(!cwf.get_output_voltages().empty());
        auto outputVoltageData = cwf.get_output_voltages()[0].get_data();
        if (!outputVoltageData.empty()) {
            double outV_avg = std::accumulate(outputVoltageData.begin(), outputVoltageData.end(), 0.0) / outputVoltageData.size();
            INFO("Output voltage average: " << outV_avg << " V");
            CHECK(outV_avg > 38.0);
            CHECK(outV_avg < 55.0);
        }
    }

    TEST_CASE("Test_PushPull_NumPeriods_SimulatedOperatingPoints", "[converter-model][push-pull-topology][num-periods][ngspice-simulation][smoke-test]") {
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::PushPull pushpull;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(24.0);
        pushpull.set_input_voltage(inputVoltage);
        pushpull.set_diode_voltage_drop(0.5);
        pushpull.set_efficiency(0.85);
        pushpull.set_current_ripple_ratio(0.3);

        PushPullOperatingPoint opPoint;
        opPoint.set_output_voltages({48.0});
        opPoint.set_output_currents({2.0});
        opPoint.set_switching_frequency(100e3);
        opPoint.set_ambient_temperature(25.0);
        pushpull.set_operating_points({opPoint});

        auto designReqs = pushpull.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        // Simulate with 1 period
        pushpull.set_num_periods_to_extract(1);
        auto ops1 = pushpull.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!ops1.empty());
        auto voltageWf1 = ops1[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();

        // Simulate with 3 periods
        pushpull.set_num_periods_to_extract(3);
        auto ops3 = pushpull.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!ops3.empty());
        auto voltageWf3 = ops3[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();

        INFO("1-period waveform data size: " << voltageWf1.get_data().size());
        INFO("3-period waveform data size: " << voltageWf3.get_data().size());

        CHECK(voltageWf3.get_data().size() > voltageWf1.get_data().size());
    }

    TEST_CASE("Test_PushPull_NumPeriods_ConverterWaveforms", "[converter-model][push-pull-topology][num-periods][ngspice-simulation][smoke-test]") {
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::PushPull pushpull;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(24.0);
        pushpull.set_input_voltage(inputVoltage);
        pushpull.set_diode_voltage_drop(0.5);
        pushpull.set_efficiency(0.85);
        pushpull.set_current_ripple_ratio(0.3);

        PushPullOperatingPoint opPoint;
        opPoint.set_output_voltages({48.0});
        opPoint.set_output_currents({2.0});
        opPoint.set_switching_frequency(100e3);
        opPoint.set_ambient_temperature(25.0);
        pushpull.set_operating_points({opPoint});

        auto designReqs = pushpull.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        // Simulate with 1 period
        auto waveforms1 = pushpull.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, 1);
        REQUIRE(!waveforms1.empty());
        auto inputV1 = waveforms1[0].get_input_voltage().get_data();

        // Simulate with 3 periods
        auto waveforms3 = pushpull.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, 3);
        REQUIRE(!waveforms3.empty());
        auto inputV3 = waveforms3[0].get_input_voltage().get_data();

        INFO("1-period converter waveform data size: " << inputV1.size());
        INFO("3-period converter waveform data size: " << inputV3.size());

        CHECK(inputV3.size() > inputV1.size());
    }

    TEST_CASE("Test_PushPull_Debug_Circuit", "[converter-model][push-pull-topology][ngspice-simulation][debug]") {
        // Debug test: use frontend defaults and print circuit + waveform analysis
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        // Test with min/max voltage range to see dead-time freewheeling at max Vin
        OpenMagnetics::PushPull pushpull;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_minimum(20.0);
        inputVoltage.set_maximum(30.0);
        pushpull.set_input_voltage(inputVoltage);

        pushpull.set_diode_voltage_drop(0.7);
        pushpull.set_efficiency(0.9);
        pushpull.set_current_ripple_ratio(0.3);
        pushpull.set_maximum_switch_current(1.0);

        PushPullOperatingPoint opPoint;
        opPoint.set_output_voltages({48.0});
        opPoint.set_output_currents({0.7});
        opPoint.set_switching_frequency(100000.0);
        opPoint.set_ambient_temperature(25.0);
        pushpull.set_operating_points({opPoint});

        auto designReqs = pushpull.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        std::cout << "\n=== PUSH-PULL DEBUG ===" << std::endl;
        std::cout << "Turns ratios: ";
        for (auto tr : turnsRatios) std::cout << tr << " ";
        std::cout << std::endl;
        std::cout << "Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH" << std::endl;
        std::cout << "Duty cycle: " << pushpull.get_maximum_duty_cycle() << std::endl;

        // Generate and print the circuit
        std::string circuit = pushpull.generate_ngspice_circuit(turnsRatios, magnetizingInductance, 0, 0);
        std::cout << "\n=== GENERATED CIRCUIT ===" << std::endl;
        std::cout << circuit << std::endl;

        // Run simulation - analyze ALL operating points (different input voltages)
        pushpull.set_num_periods_to_extract(1);
        auto operatingPoints = pushpull.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!operatingPoints.empty());

        std::cout << "\nNumber of operating points: " << operatingPoints.size() << std::endl;

        // Print waveform data for EACH operating point
        for (size_t opIdx = 0; opIdx < operatingPoints.size(); ++opIdx) {
            auto& ops = operatingPoints[opIdx];
            std::cout << "\n=== OPERATING POINT " << opIdx << ": " << ops.get_name().value_or("unnamed") << " ===" << std::endl;
            for (size_t w = 0; w < ops.get_excitations_per_winding().size(); ++w) {
                auto& exc = ops.get_excitations_per_winding()[w];
                std::string windingName = exc.get_name().value_or("Winding " + std::to_string(w));
                if (exc.get_current().has_value() && exc.get_current()->get_waveform().has_value()) {
                    auto& iData = exc.get_current()->get_waveform()->get_data();
                    double iMax = *std::max_element(iData.begin(), iData.end());
                    double iMin = *std::min_element(iData.begin(), iData.end());
                    double iAvg = std::accumulate(iData.begin(), iData.end(), 0.0) / iData.size();
                    double iRms = 0;
                    for (auto v : iData) iRms += v * v;
                    iRms = std::sqrt(iRms / iData.size());
                    std::cout << "  " << windingName << " I: min=" << iMin << " max=" << iMax << " avg=" << iAvg << " rms=" << iRms << " pts=" << iData.size() << std::endl;
                    // Print 20 evenly-spaced current samples for ONE period
                    size_t iStep = iData.size() / 20;
                    if (iStep < 1) iStep = 1;
                    std::cout << "    I samples: ";
                    for (size_t i = 0; i < iData.size(); i += iStep) {
                        std::cout << std::fixed << std::setprecision(4) << iData[i] << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }

        // Also get converter waveforms for first operating point
        auto cwfs = pushpull.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        REQUIRE(!cwfs.empty());
        auto& cwf = cwfs[0];

        std::cout << "\n=== CONVERTER WAVEFORMS ===" << std::endl;
        auto outV = cwf.get_output_voltages()[0].get_data();
        double outVAvg = std::accumulate(outV.begin(), outV.end(), 0.0) / outV.size();
        double outVMax = *std::max_element(outV.begin(), outV.end());
        double outVMin = *std::min_element(outV.begin(), outV.end());
        std::cout << "Output V: min=" << outVMin << " max=" << outVMax << " avg=" << outVAvg << std::endl;
        std::cout << "Expected: " << 48.0 << "V" << std::endl;

        // Print expected analytical values
        std::cout << "\n=== EXPECTED ANALYTICAL VALUES ===" << std::endl;
        std::cout << "Turns ratio N (Ns/Np): " << turnsRatios[1] << std::endl;
        std::cout << "Expected primary I avg (reflected load): " << (0.7 / turnsRatios[1]) << " A" << std::endl;
        std::cout << "Peak magnetizing current at Vin_min: " << (20.0 * 4.9e-6 / magnetizingInductance) << " A" << std::endl;
        std::cout << "Expected secondary I avg: 0.7 A" << std::endl;

        CHECK(true);  // Always pass - this is a debug test
    }

    TEST_CASE("Test_PushPull_PtP_AnalyticalVsNgspice", "[converter-model][push-pull-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        OpenMagnetics::PushPull pp;
        DimensionWithTolerance iv; iv.set_nominal(24.0); iv.set_minimum(18.0); iv.set_maximum(32.0);
        pp.set_input_voltage(iv);
        pp.set_diode_voltage_drop(0.5);
        pp.set_efficiency(0.9);
        pp.set_current_ripple_ratio(0.3);

        PushPullOperatingPoint op;
        op.set_output_voltages({12.0}); op.set_output_currents({4.0});
        op.set_switching_frequency(100e3); op.set_ambient_temperature(25.0);
        pp.set_operating_points({op});

        auto designReqs = pp.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
        double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto analyticalOps = pp.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        REQUIRE(!aData.empty());
        REQUIRE(!aTime.empty());
        auto aResampled = ptp_interp(aTime, aData, 256);

        pp.set_num_periods_to_extract(1);
        auto simOps = pp.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        REQUIRE(!sData.empty());
        REQUIRE(!sTime.empty());
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        INFO("PushPull primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.35);
    }

// End of SUITE

}  // namespace
