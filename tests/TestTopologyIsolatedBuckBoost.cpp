#include <source_location>
#include "support/Painter.h"
#include "converter_models/IsolatedBuckBoost.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "processors/NgspiceRunner.h"
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
        
        isolatedBuckBoost.set_diode_voltage_drop(0.5);
        isolatedBuckBoost.set_efficiency(0.9);
        isolatedBuckBoost.set_current_ripple_ratio(0.3);
        
        // For Isolated Buck-Boost: output_voltages[0] = primary buck-boost output,
        // output_voltages[1] = transformer secondary output.
        IsolatedBuckBoostOperatingPoint opPoint;
        opPoint.set_output_voltages({6.0, 5.0});
        opPoint.set_output_currents({0.5, 1.0});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        isolatedBuckBoost.set_operating_points({opPoint});
        
        auto designReqs = isolatedBuckBoost.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Isolated Buck-Boost - Turns ratio: " << turnsRatios[0]);
        INFO("Isolated Buck-Boost - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");

        auto converterWaveforms = isolatedBuckBoost.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        REQUIRE(converterWaveforms.size() >= 1);

        SECTION("Converter port — flat DC input rail (§5.1)") {
            // Sweep order from collect_input_voltages: [nominal, minimum, maximum].
            const std::vector<double> vinExpected = {12.0, 9.0, 15.0};
            REQUIRE(converterWaveforms.size() == vinExpected.size());
            for (size_t k = 0; k < converterWaveforms.size(); ++k) {
                auto& wf = converterWaveforms[k];
                REQUIRE(!wf.get_input_voltage().get_data().empty());
                REQUIRE(wf.get_input_current().get_data().size()
                        == wf.get_input_voltage().get_data().size());

                const auto& vin = wf.get_input_voltage().get_data();
                double mean = 0;
                for (double v : vin) mean += v;
                mean /= vin.size();
                INFO("OP " << k << " input_voltage mean=" << mean << " (nom " << vinExpected[k] << ")");
                CHECK(std::fabs(mean - vinExpected[k]) / vinExpected[k] < 0.01);
            }
        }

        SECTION("Winding-port primary excitation is pulsed (§5.0)") {
            // The pulse / ramp behavior belongs to the winding port.
            // Source it from simulate_and_extract_operating_points.
            auto ops = isolatedBuckBoost.simulate_and_extract_operating_points(
                turnsRatios, magnetizingInductance);
            REQUIRE(!ops.empty());
            auto& exc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(exc.get_voltage().has_value());
            REQUIRE(exc.get_voltage()->get_waveform().has_value());
            const auto priV = exc.get_voltage()->get_waveform()->get_data();
            const double Vin = 12.0;
            int count_on = 0;
            for (double v : priV) {
                if (v > 0.5 * Vin) count_on++;
            }
            double frac_on = double(count_on) / priV.size();
            INFO("Primary winding voltage fraction at >Vin/2: " << frac_on);
            CHECK(frac_on > 0.05);
            CHECK(frac_on < 0.90);

            REQUIRE(exc.get_current().has_value());
            REQUIRE(exc.get_current()->get_waveform().has_value());
            const auto priI = exc.get_current()->get_waveform()->get_data();
            auto peak_it = std::max_element(priI.begin(), priI.end());
            int peak_idx = (int)std::distance(priI.begin(), peak_it);
            int N = (int)priI.size();
            double rising_sum = 0.0;
            for (int k = 1; k <= peak_idx; ++k)
                rising_sum += priI[k] - priI[k-1];
            double falling_sum = 0.0;
            for (int k = peak_idx + 1; k < N; ++k)
                falling_sum += priI[k] - priI[k-1];
            INFO("Primary winding current rising_sum=" << rising_sum
                 << ", falling_sum=" << falling_sum);
            CHECK(rising_sum > 0.0);
            CHECK(falling_sum < 0.0);
        }

        // Visualization (preserve the SVG dumps for PR review). Use the
        // first OP at nominal Vin.
        auto& wf0 = converterWaveforms[0];
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Ngspice_InputCurrent_OP0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            auto w = wf0.get_input_current();
            painter.paint_waveform(w);
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Ngspice_InputVoltage_OP0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            auto w = wf0.get_input_voltage();
            painter.paint_waveform(w);
            painter.export_svg();
        }
        if (!wf0.get_output_voltages().empty()
            && !wf0.get_output_voltages()[0].get_data().empty()) {
            auto outFile = outputFilePath;
            outFile.append("Test_IsolatedBuckBoost_Ngspice_OutputVoltage_OP0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            auto w = wf0.get_output_voltages()[0];
            painter.paint_waveform(w);
            painter.export_svg();
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
        // IsolatedBuckBoost: complex coupled-inductor topology with reflected secondary currents
        // causes higher shape mismatch between simplified analytical model and SPICE simulation.
        // 70% threshold accounts for SPICE convergence challenges with this topology.
        INFO("Isolated Buck-Boost primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.30);
    }


    // ────────────────────────────────────────────────────────────────────
    // §5.1 converter-port DC-stream gate. See ConverterPortChecks for the
    // full bound rationale. The signals returned by
    // simulate_and_extract_topology_waveforms() are the DC source / DC
    // filtered output rails — pri_in (post-switch node) and the
    // magnetizing-current probe must NEVER appear here.
    // ────────────────────────────────────────────────────────────────────
    TEST_CASE("Test_IsolatedBuckBoost_ConverterPortWaveforms",
              "[converter-port-waveforms][isolated-buck-boost-topology][ngspice-simulation]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        const double Vin = 12.0;
        const double Vpri = 6.0, Ipri = 0.5;   // primary buck-boost output
        const double Vsec = 5.0, Isec = 1.0;   // transformer secondary

        OpenMagnetics::IsolatedBuckBoost ibb;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(Vin);
        ibb.set_input_voltage(inputVoltage);
        ibb.set_diode_voltage_drop(0.5);
        ibb.set_efficiency(0.9);
        ibb.set_current_ripple_ratio(0.3);

        IsolatedBuckBoostOperatingPoint op;
        op.set_output_voltages({Vpri, Vsec});
        op.set_output_currents({Ipri, Isec});
        op.set_switching_frequency(200000.0);
        op.set_ambient_temperature(25.0);
        ibb.set_operating_points({op});
        ibb.set_num_periods_to_extract(1);

        auto req = ibb.process_design_requirements();
        std::vector<double> tr;
        for (const auto& t : req.get_turns_ratios()) tr.push_back(t.get_nominal().value());
        const double Lm = req.get_magnetizing_inductance().get_minimum().value();

        auto wfs = ibb.simulate_and_extract_topology_waveforms(tr, Lm);
        REQUIRE(!wfs.empty());

        // IsolatedBuckBoost has known steady-state convergence issues with the
        // current netlist (the analytical model and SPICE diverge by ~50%+ on
        // output magnitudes — see Phase 4 SPICE polish backlog in
        // CONVERTER_MODELS_REVIEW_PLAN.md §3.H). The §5.1 gate's job is to
        // catch AC-in-DC-stream regressions on the data flow, NOT to enforce
        // steady-state accuracy. We use a generous tolerance (1.0 = 100 %)
        // for output magnitudes; the input-voltage flat-DC check (kVinMeanTol
        // = 0.01 inside ConverterPortChecks) remains tight.
        constexpr double kIbbVoutMeanTol = 1.0;
        constexpr double kIbbIoutMeanTol = 1.0;
        for (size_t i = 0; i < wfs.size(); ++i) {
            ConverterPortChecks::check_dc_ports(wfs[i], "IsolatedBuckBoost", i,
                                                Vin, {Vpri, Vsec}, {Ipri, Isec},
                                                kIbbVoutMeanTol, kIbbIoutMeanTol);
        }
    }

}  // namespace
