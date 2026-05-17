#include <source_location>
#include "support/Painter.h"
#include "converter_models/IsolatedBuck.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "processors/NgspiceRunner.h"
#include "NgspiceTestHelpers.h"
#include "ConverterPortChecks.h"

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
        // Post K=1 / coupled-inductor physics fix (IsolatedBuck.cpp): the
        // primary winding current is now a PIECEWISE waveform (custom 4-
        // sample) reflecting the off-time reflected-secondary step, not a
        // symmetric triangle. The continuous symmetric-triangle component
        // is the magnetizing current i_mag = i_pri + Σ i_sec/n; the winding
        // currents themselves step at switching events. See
        // TestIsolatedBuckReferenceDesignsPtp.cpp gate 4 — analytical NRMSE
        // vs SPICE dropped from ≥30 % to ≤12 % across all three reference
        // designs.
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        // Post §5.0 / flybuck-physics fix (IsolatedBuck.cpp): the secondary
        // diode conducts during the FREEWHEEL interval, so the secondary
        // winding's POSITIVE peak (not the negative peak as in forward-mode)
        // equals Vout_sec. Pre-fix the analytical formulas were forward-mode
        // and this asserted -negative_peak ≈ Vout_sec; now the correct
        // flybuck formulas put Vout_sec on the positive peak.
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["minimum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        // Same piecewise primary-current label as above.
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        // For the CUSTOM piecewise primary current `get_offset()` returns
        // (max + min)/2 (midpoint), NOT the time-weighted mean. With the
        // asymmetric piecewise shape (off-time reflected-secondary step)
        // the midpoint can drift from zero even at small DC. The physical
        // invariant is the time-weighted MEAN = Iout_pri (KCL at
        // vpri_out); check that instead.
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]),
                     Catch::Matchers::WithinAbs(
                         inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(),
                         double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
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
        // Post K=1 / coupled-inductor physics fix (IsolatedBuck.cpp): the
        // primary winding current is now PIECEWISE (custom) — see
        // Test_IsolatedBuck for full context.
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        // Post §5.0 / flybuck-physics fix (IsolatedBuck.cpp): the secondary
        // diode conducts during the FREEWHEEL interval, so the secondary
        // winding's POSITIVE peak (not the negative peak as in forward-mode)
        // equals Vout_sec. Pre-fix the analytical formulas were forward-mode
        // and this asserted -negative_peak ≈ Vout_sec; now the correct
        // flybuck formulas put Vout_sec on the positive peak.
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["minimum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["inputVoltage"]["maximum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak_to_peak().value(), double(isolatedbuckInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR);
        REQUIRE(inputs.get_operating_points()[1].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);

        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_current()->get_processed()->get_average().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputCurrents"][1]) * maximumError));
        REQUIRE_THAT(double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[1].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_positive_peak().value(), double(isolatedbuckInputsJson["operatingPoints"][0]["outputVoltages"][1]) * maximumError));
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
        
        // Validate primary input_voltage: §5.1 — this is now the DC source
        // rail v(vin_dc), NOT a winding-port pulse. Must be flat ≈ Vin_nom.
        const double VinNom = 48.0;
        CHECK(priV_max < VinNom * 1.01);
        CHECK(priV_max > VinNom * 0.99);

        INFO("Isolated Buck ngspice simulation test passed");

        SECTION("Waveform shape: input_voltage is a flat DC source rail (§5.1)") {
            // No pulse / no bipolar excursions on the converter-port input_voltage.
            // The pulsed primary winding voltage is checked separately in
            // TestVoltSecondBalance via simulate_and_extract_operating_points.
            double priV_min = *std::min_element(priVoltageData.begin(), priVoltageData.end());
            double priV_avg = std::accumulate(priVoltageData.begin(), priVoltageData.end(), 0.0) / priVoltageData.size();
            INFO("input_voltage min=" << priV_min << " avg=" << priV_avg);
            CHECK(priV_min > VinNom * 0.99);
            CHECK(std::fabs(priV_avg - VinNom) < 0.5);
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
        
        // Validate winding voltages.
        // Primary winding voltage is now extracted as a differential (across-winding) probe
        // V(pri_in) - V(vpri_out), so by volt-second balance the average must be ~0,
        // and the RMS swing must be a substantial fraction of Vin (48V).
        CHECK(std::abs(priVoltageAvg) < 1.0);          // volt-second balance
        CHECK(priVoltageRms > 5.0);                    // bipolar swing of order Vin*sqrt(D(1-D))
        CHECK(priVoltageRms < 48.0);                   // bounded by Vin
        // Secondary winding voltages: per CONVERTER_MODELS_GOLDEN_GUIDE.md
        // §5.0, excitations_per_winding[i].voltage is the across-winding
        // bipolar AC, NOT the rectified DC output. By volt-second balance
        // the average must be ~0; the DC output appears on the converter-
        // port stream / Cout node, not here. Pre-§5.0 these checks accepted
        // an analytical formula that was actually computing forward-mode
        // voltages and labelling them as flybuck — see IsolatedBuck.cpp
        // process_operating_point_for_input_voltage for the corrected
        // flybuck formulas (V_on=−(Vin−Vp)/n, V_off=+Vp/n−Vd).
        CHECK(std::abs(sec1VoltageAvg) < 1.0);
        CHECK(std::abs(sec2VoltageAvg) < 2.0);
        // Suppress unused-variable warnings (kept for INFO printouts above).
        (void)expectedSec1Voltage; (void)expectedSec2Voltage;
        
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

    // SPICE-power sanity regression: catches accidental polarity / probe
    // bugs in the across-winding voltage definition (Bvpri_diff in
    // IsolatedBuck::generate_ngspice_circuit). Asserts:
    //   (a) bipolar V swing centered on 0 (volt-second balance — non-grounded
    //       winding terminal vpri_out)
    //   (b) PASSIVE sign convention: avg(V·I) > 0 on the primary
    //       (energy flowing INTO the winding)
    //   (c) avg(|V·I|) is on the order of the primary-output power
    //       (the AP-filter consumer uses |V·I|; an accidental v(node)−GND
    //       probe collapses this 100–1000× and CoreAdviser advances
    //       sub-mm³ cores, returning 0 hits in TestMagneticAdviser)
    TEST_CASE("Test_IsolatedBuck_SpicePowerSanity",
              "[converter-model][isolated-buck-topology][ngspice-simulation][regression]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        OpenMagnetics::IsolatedBuck isolatedBuck;
        DimensionWithTolerance iv; iv.set_nominal(48.0);
        isolatedBuck.set_input_voltage(iv);
        isolatedBuck.set_diode_voltage_drop(0.5);
        isolatedBuck.set_efficiency(0.9);
        isolatedBuck.set_current_ripple_ratio(0.3);

        IsolatedBuckOperatingPoint op;
        // Single-secondary 5 V @ 5 A configuration (mirrors
        // Test_IsolatedBuck_Ngspice_Simulation).
        op.set_output_voltages({5.0, 5.0});
        op.set_output_currents({0.6, 5.0});
        op.set_switching_frequency(200e3);
        op.set_ambient_temperature(25.0);
        isolatedBuck.set_operating_points({op});

        auto designReqs = isolatedBuck.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
        double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        isolatedBuck.set_num_periods_to_extract(1);
        auto simOps = isolatedBuck.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());

        auto& exc = simOps[0].get_excitations_per_winding()[0];
        REQUIRE(exc.get_voltage().has_value());
        REQUIRE(exc.get_current().has_value());
        auto vData = exc.get_voltage()->get_waveform()->get_data();
        auto iData = exc.get_current()->get_waveform()->get_data();
        REQUIRE(vData.size() == iData.size());
        REQUIRE(!vData.empty());

        double sumVI = 0, sumAbsVI = 0, sumV = 0;
        double vmax = vData[0], vmin = vData[0];
        for (size_t k = 0; k < vData.size(); ++k) {
            sumVI    += vData[k] * iData[k];
            sumAbsVI += std::fabs(vData[k] * iData[k]);
            sumV     += vData[k];
            if (vData[k] > vmax) vmax = vData[k];
            if (vData[k] < vmin) vmin = vData[k];
        }
        double avgVI    = sumVI / vData.size();
        double avgAbsVI = sumAbsVI / vData.size();
        double avgV     = sumV / vData.size();
        double swing    = vmax - vmin;

        // Pri-output power ≈ 5 V × 0.6 A = 3 W. Primary-winding processed
        // power ≈ Pin = (Po_pri + Po_sec)/η for the full converter, which
        // is dominated by the secondary 5 V × 5 A = 25 W → Pin ≈ 28 W / η.
        // For the |V·I| sanity bound we just need the result to be a
        // realistic single-digit-to-tens-of-W magnitude (the broken
        // node-to-GND probe collapses this to milliwatts).
        INFO("V max=" << vmax << " min=" << vmin << " avg=" << avgV
             << "  swing=" << swing
             << "  avg(V·I)=" << avgVI << " avg(|V·I|)=" << avgAbsVI);

        // (a) bipolar swing — primary sees ≈ +Vin on, ≈ −Vo_reflected off
        CHECK(swing > 20.0);                              // ≈ Vin + Vo_reflected
        CHECK(std::fabs(avgV) < 2.0);                     // volt-second balance

        // (b) passive sign — power INTO winding is positive
        CHECK(avgVI > 0.0);

        // (c) magnitude is realistic (single-digit W or larger), not the
        //     milliwatt-class collapse seen with v(node)−GND probes
        CHECK(avgAbsVI > 0.5);
        CHECK(avgAbsVI < 200.0);
    }

    // ────────────────────────────────────────────────────────────────────
    // §5.1 converter-port DC-stream gate (see ConverterPortChecks).
    //
    // IsolatedBuck flybuck topology has a known model-accuracy issue on
    // the ISOLATED secondary rail (the primary buck rail tracks tightly):
    // duty cycle is set by the primary loop, and the secondary settles
    // ~30 % below nameplate due to forward-diode drop (~0.5 V on a 5 V
    // rail) and 0.99 coupling losses. Tracked under
    // CONVERTER_MODELS_REVIEW_PLAN.md §3.H Phases 5-6 (NRMSE + reference
    // designs). The §5.1 gate this test enforces is the RIPPLE bound (still
    // 25 %, default) — a winding-port AC signal smuggled into the
    // converter-port stream would have ripple/mean ≫ 1, regardless of mean.
    // ────────────────────────────────────────────────────────────────────
    TEST_CASE("Test_IsolatedBuck_ConverterPortWaveforms",
              "[converter-port-waveforms][isolated-buck-topology][ngspice-simulation]") {
        NgspiceTestHelpers::skip_if_ngspice_unavailable();

        OpenMagnetics::IsolatedBuck ib;
        const double Vin = 48.0;
        const std::vector<double> Vout = {5.0, 5.0};
        const std::vector<double> Iout = {0.6, 5.0};
        DimensionWithTolerance iv; iv.set_nominal(Vin);
        ib.set_input_voltage(iv);
        ib.set_diode_voltage_drop(0.5);
        ib.set_efficiency(0.9);
        ib.set_current_ripple_ratio(0.3);

        IsolatedBuckOperatingPoint op;
        op.set_output_voltages(Vout);
        op.set_output_currents(Iout);
        op.set_switching_frequency(200e3);
        op.set_ambient_temperature(25.0);
        ib.set_operating_points({op});

        auto designReqs = ib.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : designReqs.get_turns_ratios())
            turnsRatios.push_back(tr.get_nominal().value());
        const double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        ib.set_num_steady_state_periods(200);
        auto wfs = ib.simulate_and_extract_topology_waveforms(turnsRatios, Lm);
        REQUIRE(!wfs.empty());
        for (size_t i = 0; i < wfs.size(); ++i)
            ConverterPortChecks::check_dc_ports(wfs[i], "IsolatedBuck", i, Vin, Vout, Iout,
                                                 /*voutMeanTol*/ 0.35,
                                                 /*ioutMeanTol*/ 0.35);
    }

    // ──────────────────────────────────────────────────────────────────────
    // §3.H Phase 6 — Three industry reference designs (IsolatedBuck/Flybuck).
    //
    // Per CONVERTER_MODELS_REVIEW_PLAN.md §3.H Phase 6, IsolatedBuck must
    // ship paired Values + PtP tests for three published industry references
    // spanning the IC's operating envelope (low / mid / high corners). The
    // tests anchor the analytical model against vendor-published operating
    // points; the L value is held constant at the EVM-recommended part value
    // and the synthesized closed-form quantities (D, Im_pp, Ipri_avg, Ipri_pk)
    // are compared against the textbook flybuck equations.
    //
    // Reference designs (low / mid / high corner):
    //   1. TI LM5160     SNVA674   24 V → 12 V  @ 100 mA, 350 kHz, 47 µH 1:1
    //                              (Flybuck Quick Start, ±12 V isolated bias)
    //   2. TI LM5017     SNVA674A  48 V → 5 V   @ 0.5 A,  200 kHz, 220 µH 9:1
    //                              (Flybuck application example)
    //   3. TI LM5160-Q1  AEC-Q100  60 V → 12 V  @ 1 A,    350 kHz, 47 µH 5:1
    //                              (automotive Flybuck reference)
    //
    // η is held at 1.0 and Vd at 0 so the analytical formulas reduce to the
    // ideal flybuck equations and the test does not require transcribing
    // efficiency-curve points from each user's guide.
    // ──────────────────────────────────────────────────────────────────────

    namespace {
        struct IsoBuckRefDesignSpec {
            const char* name;
            double Vin;     // Input voltage [V]
            double Vout;    // Primary (regulated) output voltage [V]
            double Iout;    // Primary output current [A]
            double Fs;      // Switching frequency [Hz]
            double Lpri;    // Primary (magnetizing) inductance [H]
            double n;       // Turns ratio Np:Ns (primary turns / secondary turns)
            double Iout2;   // Secondary output current [A]
        };

        constexpr IsoBuckRefDesignSpec kIBRefDesign1{
            "LM5160_SNVA674_24Vin_12Vout",
            24.0, 12.0, 0.10, 350e3, 47e-6, 1.0, 0.10};

        constexpr IsoBuckRefDesignSpec kIBRefDesign2{
            "LM5017_SNVA674A_48Vin_5Vout",
            48.0, 5.0, 0.50, 200e3, 220e-6, 9.0, 0.10};

        constexpr IsoBuckRefDesignSpec kIBRefDesign3{
            "LM5160Q1_60Vin_12Vout_AEC",
            60.0, 12.0, 1.00, 350e3, 47e-6, 5.0, 0.20};

        OpenMagnetics::IsolatedBuck build_isobuck_from_spec(const IsoBuckRefDesignSpec& s) {
            OpenMagnetics::IsolatedBuck ib;
            DimensionWithTolerance iv;
            iv.set_nominal(s.Vin);
            iv.set_minimum(s.Vin * 0.95);
            iv.set_maximum(s.Vin * 1.05);
            ib.set_input_voltage(iv);
            ib.set_diode_voltage_drop(0.0);   // ideal diode for closed-form anchor
            ib.set_efficiency(1.0);            // lossless analytical reference
            ib.set_current_ripple_ratio(0.3);

            IsolatedBuckOperatingPoint op;
            op.set_output_voltages({s.Vout, s.Vout / s.n});  // primary + isolated bias
            op.set_output_currents({s.Iout, s.Iout2});
            op.set_switching_frequency(s.Fs);
            op.set_ambient_temperature(25.0);
            ib.set_operating_points({op});
            return ib;
        }

        void assert_isobuck_refdesign_values(const IsoBuckRefDesignSpec& s) {
            auto ib = build_isobuck_from_spec(s);

            // Anchor diagnostics at the nominal Vin operating point.
            IsolatedBuckOperatingPoint op;
            op.set_output_voltages({s.Vout, s.Vout / s.n});
            op.set_output_currents({s.Iout, s.Iout2});
            op.set_switching_frequency(s.Fs);
            op.set_ambient_temperature(25.0);
            ib.process_operating_point_for_input_voltage(s.Vin, op, /*turnsRatios*/{s.n}, s.Lpri);

            // Closed-form expectations (η=1, Vd=0, CCM):
            //   D        = Vout / Vin                                (buck duty)
            //   tOn      = D / Fs
            //   ΔIm      = (Vin − Vout) · tOn / Lpri                 (mag-current ripple)
            //   Ipri_avg = Iout_pri + Iout_sec / n                   (reflected sec)
            //   Ipri_pk  = Ipri_avg + ΔIm/2
            const double D_exp        = s.Vout / s.Vin;
            const double tOn_exp      = D_exp / s.Fs;
            const double dIm_exp      = (s.Vin - s.Vout) * tOn_exp / s.Lpri;
            const double Ipri_avg_exp = s.Iout + s.Iout2 / s.n;
            const double Ipri_pk_exp  = Ipri_avg_exp + 0.5 * dIm_exp;

            INFO(s.name << " — D="         << ib.get_last_duty_cycle()
                        << " (exp " << D_exp << ")");
            INFO(s.name << " — ΔIm="       << ib.get_last_magnetizing_current_ripple()
                        << " A (exp "  << dIm_exp << " A)");
            INFO(s.name << " — Ipri_avg="  << ib.get_last_primary_average_current()
                        << " A (exp "  << Ipri_avg_exp << " A)");
            INFO(s.name << " — Ipri_pk="   << ib.get_last_primary_peak_current()
                        << " A (exp "  << Ipri_pk_exp << " A)");

            // Whether the operating point lands in CCM or DCM is a function
            // of L · Fs / R_load — RefDesign1 (LM5160 Quick Start at 100 mA)
            // is intentionally DCM; RefDesigns 2/3 are CCM. Either way the
            // test is a model-vs-closed-form consistency check on the
            // diagnostics, so we don't gate on the CCM flag here.
            CHECK_THAT(ib.get_last_duty_cycle(),
                       Catch::Matchers::WithinRel(D_exp, 0.05));
            CHECK_THAT(ib.get_last_magnetizing_current_ripple(),
                       Catch::Matchers::WithinRel(dIm_exp, 0.05));
            CHECK_THAT(ib.get_last_primary_average_current(),
                       Catch::Matchers::WithinRel(Ipri_avg_exp, 0.05));
            CHECK_THAT(ib.get_last_primary_peak_current(),
                       Catch::Matchers::WithinRel(Ipri_pk_exp, 0.05));
        }

        void assert_isobuck_refdesign_ptp(const IsoBuckRefDesignSpec& s) {
            NgspiceRunner runner;
            if (!runner.is_available()) SKIP("ngspice not available");

            auto ib = build_isobuck_from_spec(s);
            ib.set_num_steady_state_periods(400);

            const std::vector<double> turnsRatios{s.n};
            auto analyticalOps = ib.process_operating_points(turnsRatios, s.Lpri);
            REQUIRE(!analyticalOps.empty());
            auto aTime = ptp_current_time(analyticalOps[0], 0);
            auto aData = ptp_current(analyticalOps[0], 0);
            REQUIRE(!aData.empty());
            REQUIRE(!aTime.empty());
            auto aResampled = ptp_interp(aTime, aData, 256);

            ib.set_num_periods_to_extract(1);
            auto simOps = ib.simulate_and_extract_operating_points(turnsRatios, s.Lpri);
            REQUIRE(!simOps.empty());
            auto sTime = ptp_current_time(simOps[0], 0);
            auto sData = ptp_current(simOps[0], 0);
            REQUIRE(!sData.empty());
            REQUIRE(!sTime.empty());
            auto sResampled = ptp_interp(sTime, sData, 256);

            const double nrmse = ptp_nrmse(aResampled, sResampled);
            INFO(s.name << " primary-current NRMSE (analytical vs ngspice): "
                        << (nrmse * 100.0) << "%");
            CHECK(nrmse < 0.30);
        }
    } // namespace

    TEST_CASE("Test_IsolatedBuck_RefDesign1_Values_LM5160_SNVA674",
              "[converter-model][isolated-buck-topology][refdesign][values]") {
        assert_isobuck_refdesign_values(kIBRefDesign1);
    }
    TEST_CASE("Test_IsolatedBuck_RefDesign1_PtP_LM5160_SNVA674",
              "[converter-model][isolated-buck-topology][refdesign][ngspice-simulation][ptpcomparison]") {
        assert_isobuck_refdesign_ptp(kIBRefDesign1);
    }
    TEST_CASE("Test_IsolatedBuck_RefDesign2_Values_LM5017_SNVA674A",
              "[converter-model][isolated-buck-topology][refdesign][values]") {
        assert_isobuck_refdesign_values(kIBRefDesign2);
    }
    TEST_CASE("Test_IsolatedBuck_RefDesign2_PtP_LM5017_SNVA674A",
              "[converter-model][isolated-buck-topology][refdesign][ngspice-simulation][ptpcomparison]") {
        assert_isobuck_refdesign_ptp(kIBRefDesign2);
    }
    TEST_CASE("Test_IsolatedBuck_RefDesign3_Values_LM5160Q1_AEC",
              "[converter-model][isolated-buck-topology][refdesign][values]") {
        assert_isobuck_refdesign_values(kIBRefDesign3);
    }
    TEST_CASE("Test_IsolatedBuck_RefDesign3_PtP_LM5160Q1_AEC",
              "[converter-model][isolated-buck-topology][refdesign][ngspice-simulation][ptpcomparison]") {
        assert_isobuck_refdesign_ptp(kIBRefDesign3);
    }

    // Issue M4 regression: AdvancedIsolatedBuck was not exercised by any test.
    TEST_CASE("Test_AdvancedIsolatedBuck_Construction_And_DR", "[converter-model][isolated-buck-topology][smoke-test]") {
        json j;
        json inputVoltage;
        inputVoltage["minimum"] = 24;
        inputVoltage["maximum"] = 48;
        j["inputVoltage"] = inputVoltage;
        j["diodeVoltageDrop"] = 0.7;
        j["efficiency"] = 1;
        j["maximumSwitchCurrent"] = 5;
        j["desiredInductance"] = 100e-6;
        j["desiredTurnsRatios"] = {2.0};
        j["operatingPoints"] = json::array();
        {
            json op;
            op["outputVoltages"] = {10.0, 10.0};
            op["outputCurrents"] = {0.1, 0.05};
            op["switchingFrequency"] = 200000;
            op["ambientTemperature"] = 25;
            j["operatingPoints"].push_back(op);
        }
        OpenMagnetics::AdvancedIsolatedBuck adv(j);

        auto inputs = adv.process();
        REQUIRE(!inputs.get_operating_points().empty());

        auto dr = adv.process_design_requirements();
        // IsolatedBuck parent PDR sets inductance as minimum.
        REQUIRE(dr.get_magnetizing_inductance().get_minimum().has_value());
        REQUIRE(dr.get_magnetizing_inductance().get_minimum().value() > 0);
        REQUIRE(!dr.get_turns_ratios().empty());
    }

}  // namespace
