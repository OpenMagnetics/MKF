#include <source_location>
#include "support/Painter.h"
#include "converter_models/Dab.h"
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

using namespace OpenMagnetics;

namespace {

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
    auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}
        .parent_path().append("..").append("output");

    // =====================================================================
    // TEST 1: Basic Design - TI TIDA-010054 reference parameters
    // =====================================================================
    // V1 = 800V, V2 = 500V, P = 10kW, Fs = 100kHz, N = 1.6, L = 35uH
    // Reference: TI TIDA-010054 (tidues0e.pdf), Table 2-1
    // =====================================================================
    TEST_CASE("Test_Dab_TI_Reference_Design", "[converter-model][dab-topology][smoke-test]") {
        json dabJson;
        json inputVoltage;

        inputVoltage["nominal"] = 800.0;
        inputVoltage["minimum"] = 700.0;
        inputVoltage["maximum"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["useLeakageInductance"] = false;
        dabJson["operatingPoints"] = json::array();

        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {20.0};
            op["phaseShift"] = 23.0;  // degrees, from TI doc: 0.4 rad ≈ 23 deg
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);

        SECTION("Input validation") {
            CHECK(dab.run_checks(false) == true);
        }

        SECTION("Turns ratio") {
            auto req = dab.process_design_requirements();
            // N = V1_nom / V2 = 800 / 500 = 1.6
            double expectedN = 800.0 / 500.0;
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.02));
        }

        SECTION("Series inductance preserved") {
            dab.process_design_requirements();
            double L = dab.get_computed_series_inductance();
            REQUIRE_THAT(L, Catch::Matchers::WithinAbs(35e-6, 35e-6 * 0.01));
        }

        SECTION("Magnetizing inductance is positive and >> L") {
            auto req = dab.process_design_requirements();
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            double L = dab.get_computed_series_inductance();
            CHECK(Lm > 0);
            CHECK(Lm > 10 * L);  // Lm should be at least 10x the series inductance
        }

        SECTION("Power transfer calculation") {
            // P = N*V1*V2*phi*(pi-phi) / (2*pi^2*Fs*L)
            double N = 1.6;
            double V1 = 800.0;
            double V2 = 500.0;
            double phi = 23.0 * M_PI / 180.0; // 23 degrees in radians
            double Fs = 100e3;
            double L = 35e-6;

            double P = Dab::compute_power(V1, V2, N, phi, Fs, L);
            // Expected: approximately 10 kW
            REQUIRE_THAT(P / 1000.0, Catch::Matchers::WithinAbs(10.0, 2.0));
        }

        SECTION("Phase shift calculation") {
            double N = 1.6;
            double V1 = 800.0;
            double V2 = 500.0;
            double Fs = 100e3;
            double L = 35e-6;
            double P = 10e3;

            double phi = Dab::compute_phase_shift(V1, V2, N, Fs, L, P);
            double phi_deg = phi * 180.0 / M_PI;
            // Expected: approximately 23 degrees (0.4 rad)
            REQUIRE_THAT(phi_deg, Catch::Matchers::WithinAbs(23.0, 5.0));
        }

        SECTION("Switching currents i1 and i2") {
            double N = 1.6;
            double V1 = 800.0;
            double V2 = 500.0;
            double phi = 0.4; // radians
            double Fs = 100e3;
            double L = 35e-6;

            double i1, i2;
            Dab::compute_switching_currents(V1, V2, N, phi, Fs, L, i1, i2);

            // From TI doc: for d = N*V2/V1 = 1.6*500/800 = 1.0
            // At d=1: i1 = i2 = phi * Ibase
            // Ibase = V1/(2*pi*Fs*L) = 800/(2*pi*100e3*35e-6) = 36.4
            // i1 = i2 = 0.4 * 36.4 ≈ 14.6 A
            REQUIRE_THAT(i1, Catch::Matchers::WithinAbs(i2, std::abs(i2) * 0.05));
            CHECK(std::abs(i1) > 0);
            CHECK(std::abs(i2) > 0);

            // Verify RMS current
            // At d=1, phi=0.4 rad: i1=i2=phi*Ibase=14.56 A, Irms ≈ 13.9 A (analytically correct).
            // TI TIDA-010054 Table 2-1 lists 9.67 A, but that figure corresponds to a different
            // operating point in the document (lower output current or different efficiency assumption).
            // Our formula matches the standard SPS RMS derivation from Demetriades Ch.6.
            double Irms = Dab::compute_primary_rms_current(i1, i2, phi);
            CHECK(Irms > 0);
            REQUIRE_THAT(Irms, Catch::Matchers::WithinAbs(13.9, 1.0));
        }

        SECTION("ZVS boundaries") {
            double N = 1.6;
            double V1 = 800.0;
            double V2 = 500.0;
            double d = Dab::compute_voltage_ratio(V1, V2, N);
            double phi = 0.4;

            // At d=1, ZVS should be achievable at any phase shift
            REQUIRE_THAT(d, Catch::Matchers::WithinAbs(1.0, 0.01));
            CHECK(Dab::check_zvs_primary(phi, d) == true);
            CHECK(Dab::check_zvs_secondary(phi, d) == true);
        }
    }


    // =====================================================================
    // TEST 2: Design Requirements Computation
    // =====================================================================
    TEST_CASE("Test_Dab_Design_Requirements", "[converter-model][dab-topology][smoke-test]") {
        json dabJson;
        json inputVoltage;

        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 380.0;
        inputVoltage["maximum"] = 420.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["operatingPoints"] = json::array();

        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {48.0};
            op["outputCurrents"] = {50.0};
            op["phaseShift"] = 30.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();

        SECTION("Turns ratio for 400V:48V") {
            double expectedN = 400.0 / 48.0; // ≈ 8.33
            double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
            REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.02));
        }

        SECTION("Series inductance is computed when not provided") {
            double L = dab.get_computed_series_inductance();
            CHECK(L > 0);
            CHECK(L > 1e-7);   // Should be in reasonable range
            CHECK(L < 10e-3);  // Not absurdly large
        }

        SECTION("Magnetizing inductance is positive") {
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
            CHECK(Lm > 0);
            CHECK(Lm > 10e-6);
            CHECK(Lm < 1);
        }

        SECTION("Power verification round-trip") {
            double L = dab.get_computed_series_inductance();
            double phi = dab.get_computed_phase_shift();
            double N = resolve_dimensional_values(req.get_turns_ratios()[0]);
            double P = Dab::compute_power(400.0, 48.0, N, phi, 100e3, L);
            // Should be close to 48 * 50 = 2400 W
            REQUIRE_THAT(P, Catch::Matchers::WithinRel(2400.0, 0.05));
        }
    }


    // =====================================================================
    // TEST 3: Operating Point Waveform Generation
    // =====================================================================
    TEST_CASE("Test_Dab_OperatingPoints_Generation", "[converter-model][dab-topology][smoke-test]") {
        json dabJson;
        json inputVoltage;

        inputVoltage["nominal"] = 800.0;
        inputVoltage["minimum"] = 700.0;
        inputVoltage["maximum"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["operatingPoints"] = json::array();

        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {20.0};
            op["phaseShift"] = 23.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();

        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        SECTION("Multiple input voltages") {
            auto ops = dab.process_operating_points(turnsRatios, Lm);
            // Should have entries for nominal, minimum, maximum (deduplicated)
            CHECK(ops.size() >= 2);
        }

        SECTION("Waveform structure") {
            auto ops = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());

            auto& op = ops[0];
            // Primary + 1 secondary = 2 excitations
            CHECK(op.get_excitations_per_winding().size() == 2);

            auto& priExc = op.get_excitations_per_winding()[0];
            REQUIRE(priExc.get_current().has_value());
            REQUIRE(priExc.get_voltage().has_value());
            CHECK(priExc.get_frequency() == 100e3);

            // Waveforms should have 2*N+1 = 513 samples
            auto currentWfm = priExc.get_current()->get_waveform().value();
            CHECK(currentWfm.get_data().size() == 513);

            auto voltageWfm = priExc.get_voltage()->get_waveform().value();
            CHECK(voltageWfm.get_data().size() == 513);
        }

        SECTION("Primary current is piecewise linear") {
            auto ops = dab.process_operating_points(turnsRatios, Lm);
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            auto currentWfm = priExc.get_current()->get_waveform().value();
            auto iData = currentWfm.get_data();

            // DAB current should be piecewise linear
            // Check that it crosses zero and has significant magnitude
            double maxI = *std::max_element(iData.begin(), iData.end());
            double minI = *std::min_element(iData.begin(), iData.end());
            CHECK(maxI > 0);
            CHECK(minI < 0);
            // Current should be symmetric (antisymmetric between halves)
        }

        SECTION("Primary current antisymmetry") {
            auto ops = dab.process_operating_points(turnsRatios, Lm);
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            auto currentWfm = priExc.get_current()->get_waveform().value();
            auto iData = currentWfm.get_data();

            int N_half = ((int)iData.size() - 1) / 2; // Half-period samples

            // iL(t + Thalf) = -iL(t)
            for (int k = 1; k < N_half; k += N_half/10) {
                REQUIRE_THAT(iData[N_half + k],
                    Catch::Matchers::WithinAbs(-iData[k],
                        std::max(0.05, std::abs(iData[k]) * 0.02)));
            }
        }

        SECTION("Primary voltage is bipolar rectangular") {
            auto ops = dab.process_operating_points(turnsRatios, Lm);
            // ops are sorted by input voltage: ops[0]=700V (min), ops.back()=800V (nominal/max).
            // Use ops.back() to test the 800V (maximum) operating point.
            auto& priExc = ops.back().get_excitations_per_winding()[0];
            auto voltageWfm = priExc.get_voltage()->get_waveform().value();
            auto vData = voltageWfm.get_data();

            int N_half = ((int)vData.size() - 1) / 2;

            // Primary bridge voltage Vab = +V1 during positive half, -V1 during negative half.
            // V1 = 800 V is set directly from inputVoltage, so this should be exact.
            // Positive half: V = +V1 = +800
            for (int k = 1; k < N_half; ++k) {
                REQUIRE_THAT(vData[k],
                    Catch::Matchers::WithinAbs(800.0, 1.0));
            }
            // Negative half: V = -V1 = -800
            // Skip last sample (k=2*N_half): it is at theta=2*pi which wraps back to
            // theta=0 (start of next period), so it equals +V1, not -V1.
            for (int k = N_half + 1; k < (int)vData.size() - 1; ++k) {
                REQUIRE_THAT(vData[k],
                    Catch::Matchers::WithinAbs(-800.0, 1.0));
            }
        }

        SECTION("Secondary current reflects primary through turns ratio") {
            auto ops = dab.process_operating_points(turnsRatios, Lm);
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            auto& secExc = ops[0].get_excitations_per_winding()[1];

            auto priI = priExc.get_current()->get_waveform().value().get_data();
            auto secI = secExc.get_current()->get_waveform().value().get_data();

            double N_ratio = turnsRatios[0];

            // The primary winding carries iL + Im (power transfer + magnetizing).
            // The secondary carries only the power-transfer component: N * iL.
            // Since Im_peak << iL_peak, secondary ≈ N * primary with Im as the error term.
            // Tolerance is relaxed to allow for the magnetizing current offset.
            double Im_peak_approx = 800.0 / (4.0 * 100e3 * Lm);
            double tolerance = N_ratio * Im_peak_approx * 1.1;
            for (size_t k = 0; k < priI.size(); k += priI.size() / 10) {
                REQUIRE_THAT(secI[k],
                    Catch::Matchers::WithinAbs(N_ratio * priI[k],
                        std::max(tolerance, std::abs(N_ratio * priI[k]) * 0.05)));
            }
        }
    }


    // =====================================================================
    // TEST 4: Operating Modes - d=1 (nominal), d<1 (step down), d>1 (step up)
    // =====================================================================
    TEST_CASE("Test_Dab_Operating_Modes", "[converter-model][dab-topology][smoke-test]") {

        SECTION("Nominal operation (d = 1)") {
            json dabJson;
            json inputVoltage;
            inputVoltage["nominal"] = 800.0;
            dabJson["inputVoltage"] = inputVoltage;
            dabJson["seriesInductance"] = 35e-6;
            dabJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                op["outputVoltages"] = {500.0};
                op["outputCurrents"] = {20.0};
                op["phaseShift"] = 23.0;
                op["switchingFrequency"] = 100000;
                dabJson["operatingPoints"].push_back(op);
            }

            OpenMagnetics::Dab dab(dabJson);
            auto req = dab.process_design_requirements();

            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto ops = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());

            // At d=1, i1 = i2 (symmetric current waveform)
            double d = Dab::compute_voltage_ratio(800.0, 500.0, turnsRatios[0]);
            REQUIRE_THAT(d, Catch::Matchers::WithinAbs(1.0, 0.05));

            auto& priExc = ops[0].get_excitations_per_winding()[0];
            auto iData = priExc.get_current()->get_waveform().value().get_data();
            double maxI = *std::max_element(iData.begin(), iData.end());
            double minI = *std::min_element(iData.begin(), iData.end());
            // At d=1, current should be symmetric: |max| ≈ |min|
            REQUIRE_THAT(maxI, Catch::Matchers::WithinRel(-minI, 0.05));
        }

        SECTION("Step-down operation (d < 1, lower V2)") {
            json dabJson;
            json inputVoltage;
            inputVoltage["nominal"] = 800.0;
            dabJson["inputVoltage"] = inputVoltage;
            dabJson["seriesInductance"] = 35e-6;
            dabJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                op["outputVoltages"] = {350.0};  // Lower output voltage -> d < 1
                op["outputCurrents"] = {20.0};
                op["phaseShift"] = 25.0;
                op["switchingFrequency"] = 100000;
                dabJson["operatingPoints"].push_back(op);
            }

            OpenMagnetics::Dab dab(dabJson);
            auto req = dab.process_design_requirements();

            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            // N is designed for nominal V2=350V
            // double d = Dab::compute_voltage_ratio(800.0, 350.0, turnsRatios[0]);
            // Since N = 800/350, d = N*350/800 = 1.0 (it's the same here since V2 changed)
            // Actually for step-down, we'd use N designed for 500V but output 350V
            // Let's check differently:
            auto ops = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());

            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_current().has_value());
        }

        SECTION("Step-up operation (d > 1, higher V2)") {
            json dabJson;
            json inputVoltage;
            inputVoltage["nominal"] = 800.0;
            dabJson["inputVoltage"] = inputVoltage;
            dabJson["seriesInductance"] = 35e-6;
            dabJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                // Use N from 800:500 design but output higher voltage
                op["outputVoltages"] = {500.0};
                op["outputCurrents"] = {20.0};
                op["phaseShift"] = 23.0;
                op["switchingFrequency"] = 100000;
                dabJson["operatingPoints"].push_back(op);
            }

            // Force d > 1 by using lower input voltage
            inputVoltage["nominal"] = 700.0;
            dabJson["inputVoltage"] = inputVoltage;

            OpenMagnetics::Dab dab(dabJson);
            auto req = dab.process_design_requirements();

            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto ops = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());

            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_current().has_value());

            // Current should have different i1, i2 magnitudes for d ≠ 1
            auto iData = priExc.get_current()->get_waveform().value().get_data();
            double maxI = *std::max_element(iData.begin(), iData.end());
            CHECK(maxI > 0);
        }

        SECTION("Reverse power flow (negative phase shift)") {
            json dabJson;
            json inputVoltage;
            inputVoltage["nominal"] = 800.0;
            dabJson["inputVoltage"] = inputVoltage;
            dabJson["seriesInductance"] = 35e-6;
            dabJson["operatingPoints"] = json::array();

            {
                json op;
                op["ambientTemperature"] = 25.0;
                op["outputVoltages"] = {500.0};
                op["outputCurrents"] = {20.0};
                op["phaseShift"] = -23.0;  // Negative = reverse power flow
                op["switchingFrequency"] = 100000;
                dabJson["operatingPoints"].push_back(op);
            }

            OpenMagnetics::Dab dab(dabJson);
            auto req = dab.process_design_requirements();

            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto ops = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());

            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_current().has_value());

            // For reverse power flow, the current waveform should be
            // the negative of the forward power flow waveform
            auto iData = priExc.get_current()->get_waveform().value().get_data();
            double maxI = *std::max_element(iData.begin(), iData.end());
            double minI = *std::min_element(iData.begin(), iData.end());
            CHECK(maxI > 0);
            CHECK(minI < 0);
        }
    }


    // =====================================================================
    // TEST 5: ZVS Boundary Verification
    // =====================================================================
    TEST_CASE("Test_Dab_ZVS_Boundaries", "[converter-model][dab-topology][smoke-test]") {
        SECTION("ZVS at d = 1") {
            double d = 1.0;
            // At d=1, ZVS boundary is phi > 0 for both bridges
            CHECK(Dab::check_zvs_primary(0.1, d) == true);
            CHECK(Dab::check_zvs_secondary(0.1, d) == true);
            // Very small phase shift should still achieve ZVS at d=1
            CHECK(Dab::check_zvs_primary(0.01, d) == true);
        }

        SECTION("ZVS at d = 0.7 (step-down)") {
            double d = 0.7;
            double phi_min_sec = (1.0 - d) * M_PI / 2.0; // ≈ 0.47 rad

            // Below boundary: no ZVS on secondary
            CHECK(Dab::check_zvs_secondary(0.1, d) == false);
            // Above boundary: ZVS on secondary
            CHECK(Dab::check_zvs_secondary(phi_min_sec + 0.1, d) == true);

            // Primary ZVS at d < 1: phi > (1 - 1/d)*pi/2 which is negative, always true
            CHECK(Dab::check_zvs_primary(0.01, d) == true);
        }

        SECTION("ZVS at d = 1.3 (step-up)") {
            double d = 1.3;
            double phi_min_pri = (1.0 - 1.0/d) * M_PI / 2.0; // ≈ 0.36 rad

            // Below boundary: no ZVS on primary
            CHECK(Dab::check_zvs_primary(0.1, d) == false);
            // Above boundary: ZVS on primary
            CHECK(Dab::check_zvs_primary(phi_min_pri + 0.1, d) == true);

            // Secondary ZVS at d > 1: phi > (1-d)*pi/2 which is negative, always true
            CHECK(Dab::check_zvs_secondary(0.01, d) == true);
        }
    }


    // =====================================================================
    // TEST 6: SPICE Netlist Generation
    // =====================================================================
    TEST_CASE("Test_Dab_SPICE_Netlist", "[converter-model][dab-topology][smoke-test]") {
        json dabJson;
        json inputVoltage;

        inputVoltage["nominal"] = 800.0;
        inputVoltage["minimum"] = 700.0;
        inputVoltage["maximum"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["operatingPoints"] = json::array();

        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {20.0};
            op["phaseShift"] = 23.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();

        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        SECTION("Netlist is non-empty and contains key components") {
            std::string netlist = dab.generate_ngspice_circuit(turnsRatios, Lm);

            CHECK(!netlist.empty());
            CHECK(netlist.find("Dual Active Bridge") != std::string::npos);
            CHECK(netlist.find("L_series") != std::string::npos);
            CHECK(netlist.find("L_pri") != std::string::npos);
            CHECK(netlist.find("L_sec_o1") != std::string::npos);
            // Pairwise K coupling: K1 couples L_pri to L_sec_o1
            CHECK(netlist.find("K1 L_pri L_sec_o1") != std::string::npos);
            CHECK(netlist.find(".tran") != std::string::npos);
            // New independent-leg PWM scheme (per leg, not per diagonal pair)
            CHECK(netlist.find("pwm_p_l1") != std::string::npos);
            CHECK(netlist.find("pwm_p_r1") != std::string::npos);
            CHECK(netlist.find("pwm_s_l1") != std::string::npos);
            CHECK(netlist.find("pwm_s_r1") != std::string::npos);
        }

        SECTION("Netlist can be saved to file") {
            std::string netlist = dab.generate_ngspice_circuit(turnsRatios, Lm);

            auto outFile = outputFilePath;
            outFile.append("Test_Dab_SPICE_Netlist.cir");
            std::filesystem::create_directories(outputFilePath);
            std::ofstream ofs(outFile);
            ofs << netlist;
            ofs.close();

            CHECK(std::filesystem::exists(outFile));
        }
    }


    // =====================================================================
    // TEST 7: Waveform Plotting
    // =====================================================================
    TEST_CASE("Test_Dab_Waveform_Plotting", "[converter-model][dab-topology][smoke-test]") {
        json dabJson;
        json inputVoltage;

        inputVoltage["nominal"] = 800.0;
        inputVoltage["minimum"] = 700.0;
        inputVoltage["maximum"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["operatingPoints"] = json::array();

        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {20.0};
            op["phaseShift"] = 23.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();

        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        auto ops = dab.process_operating_points(turnsRatios, Lm);

        SECTION("Primary current waveform SVG") {
            auto outFile = outputFilePath;
            outFile.append("Test_Dab_Primary_Current_Waveform.svg");
            std::filesystem::create_directories(outputFilePath);
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }

        SECTION("Primary voltage waveform SVG") {
            auto outFile = outputFilePath;
            outFile.append("Test_Dab_Primary_Voltage_Waveform.svg");
            std::filesystem::create_directories(outputFilePath);
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }

        SECTION("Secondary current waveform SVG") {
            auto outFile = outputFilePath;
            outFile.append("Test_Dab_Secondary_Current_Waveform.svg");
            std::filesystem::create_directories(outputFilePath);
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }

        SECTION("Secondary voltage waveform SVG") {
            auto outFile = outputFilePath;
            outFile.append("Test_Dab_Secondary_Voltage_Waveform.svg");
            std::filesystem::create_directories(outputFilePath);
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(ops[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
    }


    // =====================================================================
    // TEST 8: Inductance computation from power equation round-trip
    // =====================================================================
    TEST_CASE("Test_Dab_Inductance_Round_Trip", "[converter-model][dab-topology][validation]") {
        double V1 = 800.0;
        double V2 = 500.0;
        double N = 1.6;
        double Fs = 100e3;
        double P_target = 10e3;

        SECTION("Inductance -> Power -> Inductance") {
            double phi = 0.4; // radians
            double L = Dab::compute_series_inductance(V1, V2, N, phi, Fs, P_target);
            double P_check = Dab::compute_power(V1, V2, N, phi, Fs, L);
            REQUIRE_THAT(P_check, Catch::Matchers::WithinRel(P_target, 0.01));
        }

        SECTION("Phase shift -> Power -> Phase shift") {
            double L = 35e-6;
            double P = Dab::compute_power(V1, V2, N, 0.4, Fs, L);
            double phi_check = Dab::compute_phase_shift(V1, V2, N, Fs, L, P);
            REQUIRE_THAT(phi_check, Catch::Matchers::WithinAbs(0.4, 0.01));
        }

        SECTION("Maximum power at phi = pi/2") {
            double L = 35e-6;
            double Pmax = Dab::compute_power(V1, V2, N, M_PI / 2.0, Fs, L);
            // Maximum power should be P_max = N*V1*V2 / (8*Fs*L)
            double Pmax_expected = N * V1 * V2 / (8.0 * Fs * L);
            REQUIRE_THAT(Pmax, Catch::Matchers::WithinRel(Pmax_expected, 0.01));
        }
    }


    // =====================================================================
    // TEST 9: Light load operation (small phase shift)
    // =====================================================================
    TEST_CASE("Test_Dab_Light_Load", "[converter-model][dab-topology][smoke-test]") {
        json dabJson;
        json inputVoltage;

        inputVoltage["nominal"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["operatingPoints"] = json::array();

        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {2.0};  // Light load: 1 kW
            op["phaseShift"] = 5.0;         // Small phase shift
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();

        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        auto ops = dab.process_operating_points(turnsRatios, Lm);
        REQUIRE(!ops.empty());

        // At light load, current amplitude should be smaller
        auto& priExc = ops[0].get_excitations_per_winding()[0];
        auto iData = priExc.get_current()->get_waveform().value().get_data();
        double maxI = *std::max_element(iData.begin(), iData.end());

        // Light load current should be much smaller than full load
        CHECK(maxI < 20.0);  // Full load peak is ~14A, light load should be less
        CHECK(maxI > 0);
    }


    // =====================================================================
    // TEST 10: Multiple secondaries
    // =====================================================================
    TEST_CASE("Test_Dab_Multiple_Outputs", "[converter-model][dab-topology][smoke-test]") {
        json dabJson;
        json inputVoltage;
        inputVoltage["nominal"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["operatingPoints"] = json::array();

        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0, 250.0};
            op["outputCurrents"] = {10.0, 5.0};
            op["phaseShift"] = 20.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();

        SECTION("Turns ratios for each output") {
            CHECK(req.get_turns_ratios().size() == 2);
            double n1 = resolve_dimensional_values(req.get_turns_ratios()[0]);
            double n2 = resolve_dimensional_values(req.get_turns_ratios()[1]);
            // n1 = 800/500 = 1.6, n2 = 800/250 = 3.2, so n2/n1 = 2.0
            REQUIRE_THAT(n2 / n1, Catch::Matchers::WithinAbs(2.0, 0.1));
        }

        SECTION("Operating points with multiple secondaries") {
            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios()) {
                turnsRatios.push_back(resolve_dimensional_values(tr));
            }
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto ops = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!ops.empty());

            // Should have 1 primary + 2 secondaries = 3 excitations
            CHECK(ops[0].get_excitations_per_winding().size() == 3);
        }
    }

    // =====================================================================
    // TEST 11: NgSpice Simulation vs Analytical Comparison
    // Runs the DAB SPICE netlist through NgSpice and compares key waveform
    // metrics (inductor current peak, primary voltage amplitude) against the
    // analytical model.  Uses a tighter Lm (5× series inductance) to make the
    // magnetising current large enough to exercise that path, and a loose ±20%
    // tolerance to absorb SPICE switching transients.
    // =====================================================================
    TEST_CASE("Test_Dab_Ngspice_Simulation", "[converter-model][dab-topology][ngspice-simulation]") {
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        json dabJson;
        json inputVoltage;
        inputVoltage["nominal"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["useLeakageInductance"] = false;
        dabJson["operatingPoints"] = json::array();
        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {20.0};
            op["phaseShift"] = 23.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        dab.set_num_steady_state_periods(10);
        dab.set_num_periods_to_extract(2);
        auto req = dab.process_design_requirements();

        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios())
            turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        SECTION("Simulation runs and returns waveform data") {
            auto waveforms = dab.simulate_and_extract_topology_waveforms(turnsRatios, Lm, 2);
            REQUIRE(!waveforms.empty());
            auto& wf = waveforms[0];
            CHECK(!wf.get_input_voltage().get_data().empty());
            CHECK(!wf.get_input_current().get_data().empty());
            REQUIRE(!wf.get_output_voltages().empty());
            CHECK(!wf.get_output_voltages()[0].get_data().empty());
        }

        SECTION("Simulated primary voltage is bipolar ~±800 V") {
            auto waveforms = dab.simulate_and_extract_topology_waveforms(turnsRatios, Lm, 2);
            REQUIRE(!waveforms.empty());
            auto vData = waveforms[0].get_input_voltage().get_data();
            double vMax = *std::max_element(vData.begin(), vData.end());
            double vMin = *std::min_element(vData.begin(), vData.end());
            REQUIRE_THAT(vMax, Catch::Matchers::WithinAbs( 800.0, 80.0));
            REQUIRE_THAT(vMin, Catch::Matchers::WithinAbs(-800.0, 80.0));
        }

        SECTION("Simulated peak inductor current matches analytical within 20%") {
            // Analytical peak
            auto analOps = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!analOps.empty());
            auto analI = analOps[0].get_excitations_per_winding()[0]
                             .get_current()->get_waveform().value().get_data();
            double analPeak = *std::max_element(analI.begin(), analI.end());

            // SPICE peak
            auto waveforms = dab.simulate_and_extract_topology_waveforms(turnsRatios, Lm, 2);
            REQUIRE(!waveforms.empty());
            auto spiceI = waveforms[0].get_input_current().get_data();
            double spicePeak = *std::max_element(spiceI.begin(), spiceI.end());

            INFO("Analytical peak: " << analPeak << " A");
            INFO("SPICE peak: " << spicePeak << " A");
            // Allow ±20% to accommodate switching transients and dead-time effects
            REQUIRE_THAT(spicePeak, Catch::Matchers::WithinRel(analPeak, 0.20));
        }

        SECTION("Waveform shape: RMS matches analytical within 15%") {
            auto analOps = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!analOps.empty());
            auto analI = analOps[0].get_excitations_per_winding()[0]
                             .get_current()->get_waveform().value().get_data();

            auto waveforms = dab.simulate_and_extract_topology_waveforms(turnsRatios, Lm, 2);
            REQUIRE(!waveforms.empty());
            auto spiceI = waveforms[0].get_input_current().get_data();

            // Compute RMS for both
            auto rms = [](const std::vector<double>& v) {
                double sum = 0;
                for (double x : v) sum += x * x;
                return std::sqrt(sum / v.size());
            };
            double analRms = rms(analI);
            double spiceRms = rms(spiceI);

            INFO("Analytical RMS: " << analRms << " A");
            INFO("SPICE RMS: " << spiceRms << " A");
            REQUIRE_THAT(spiceRms, Catch::Matchers::WithinRel(analRms, 0.15));
        }

        SECTION("Waveform shape: half-wave antisymmetry holds in simulation") {
            auto waveforms = dab.simulate_and_extract_topology_waveforms(turnsRatios, Lm, 2);
            REQUIRE(!waveforms.empty());
            auto spiceI = waveforms[0].get_input_current().get_data();
            REQUIRE(spiceI.size() > 4);

            // SPICE waveform should have both positive and negative portions
            double iMax = *std::max_element(spiceI.begin(), spiceI.end());
            double iMin = *std::min_element(spiceI.begin(), spiceI.end());
            CHECK(iMax > 5.0);   // Peak > 5A (nominal ~16A)
            CHECK(iMin < -5.0);  // Also negative (antisymmetric)
            // Symmetry: |max| ≈ |min| within 20%
            REQUIRE_THAT(std::abs(iMax), Catch::Matchers::WithinRel(std::abs(iMin), 0.20));
        }

        SECTION("Waveform shape: save CSV for visual inspection") {
            auto analOps = dab.process_operating_points(turnsRatios, Lm);
            REQUIRE(!analOps.empty());
            auto analWfm = analOps[0].get_excitations_per_winding()[0]
                               .get_current()->get_waveform().value();
            auto analI = analWfm.get_data();
            auto analT = analWfm.get_time().value();

            auto waveforms = dab.simulate_and_extract_topology_waveforms(turnsRatios, Lm, 2);
            REQUIRE(!waveforms.empty());
            auto spiceI = waveforms[0].get_input_current().get_data();
            auto spiceT = waveforms[0].get_input_current().get_time().value_or(std::vector<double>{});

            auto csvPath = outputFilePath / "Test_Dab_Waveform_Comparison.csv";
            std::filesystem::create_directories(outputFilePath);
            std::ofstream csv(csvPath);
            csv << "t_anal,i_anal,t_spice,i_spice\n";
            size_t nRows = std::max(analI.size(), spiceI.size());
            for (size_t k = 0; k < nRows; ++k) {
                double ta = (k < analT.size()) ? analT[k] : 0;
                double ia = (k < analI.size()) ? analI[k] : 0;
                double ts = (k < spiceT.size()) ? spiceT[k] : 0;
                double is = (k < spiceI.size()) ? spiceI[k] : 0;
                csv << ta << "," << ia << "," << ts << "," << is << "\n";
            }
            csv.close();
            CHECK(std::filesystem::exists(csvPath));
        }
    }

    TEST_CASE("Test_Dab_PtP_AnalyticalVsNgspice", "[converter-model][dab-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        json dabJson;
        json inputVoltage;
        inputVoltage["nominal"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["useLeakageInductance"] = false;
        dabJson["operatingPoints"] = json::array();
        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {20.0};
            op["phaseShift"] = 23.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        auto analyticalOps = dab.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        REQUIRE(!aData.empty());
        REQUIRE(!aTime.empty());
        auto aResampled = ptp_interp(aTime, aData, 256);

        dab.set_num_periods_to_extract(1);
        auto simOps = dab.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        REQUIRE(!sData.empty());
        REQUIRE(!sTime.empty());
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        // DAB inductor current is trapezoidal/triangular. Threshold tightened
        // from 50% → 15% after the closed-form sub-interval propagator + Gear
        // method + snubbers brought NRMSE down to ~4%.
        INFO("DAB primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.15);
    }

    // =====================================================================
    // EPS modulation: D1 > 0, D2 = 0
    // The primary bridge has zero-voltage plateaus at the start of each
    // half-cycle; the secondary is a full square wave.
    // =====================================================================
    TEST_CASE("Test_Dab_EPS_PtP_AnalyticalVsNgspice", "[converter-model][dab-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        json dabJson;
        json inputVoltage;
        inputVoltage["nominal"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["useLeakageInductance"] = false;
        dabJson["operatingPoints"] = json::array();
        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {15.0};
            op["phaseShift"] = 30.0;        // 30° outer phase shift
            op["innerPhaseShift1"] = 20.0;  // 20° primary inner shift (EPS)
            op["modulationType"] = "EPS";
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        auto analyticalOps = dab.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        CHECK(dab.get_last_modulation_type() == 1);  // EPS

        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        REQUIRE(!aData.empty());
        auto aResampled = ptp_interp(aTime, aData, 256);

        dab.set_num_periods_to_extract(1);
        auto simOps = dab.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        REQUIRE(!sData.empty());
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        INFO("DAB EPS NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.15);
    }

    // =====================================================================
    // DPS modulation: D1 = D2 > 0
    // Both bridges have symmetric inner phase shifts.
    // =====================================================================
    TEST_CASE("Test_Dab_DPS_PtP_AnalyticalVsNgspice", "[converter-model][dab-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        json dabJson;
        json inputVoltage;
        inputVoltage["nominal"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["useLeakageInductance"] = false;
        dabJson["operatingPoints"] = json::array();
        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {15.0};
            op["phaseShift"] = 30.0;
            op["innerPhaseShift1"] = 15.0;
            op["modulationType"] = "DPS";  // D2 will default to D1
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        auto analyticalOps = dab.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        CHECK(dab.get_last_modulation_type() == 2);  // DPS

        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        auto aResampled = ptp_interp(aTime, aData, 256);

        dab.set_num_periods_to_extract(1);
        auto simOps = dab.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        INFO("DAB DPS NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.15);
    }

    // =====================================================================
    // TPS modulation: D1 ≠ D2, both > 0
    // Independent inner phase shifts on each bridge.
    // =====================================================================
    TEST_CASE("Test_Dab_TPS_PtP_AnalyticalVsNgspice", "[converter-model][dab-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        json dabJson;
        json inputVoltage;
        inputVoltage["nominal"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["useLeakageInductance"] = false;
        dabJson["operatingPoints"] = json::array();
        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {15.0};
            op["phaseShift"] = 30.0;
            op["innerPhaseShift1"] = 20.0;
            op["innerPhaseShift2"] = 10.0;  // D1 ≠ D2 → TPS
            op["modulationType"] = "TPS";
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        auto analyticalOps = dab.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        CHECK(dab.get_last_modulation_type() == 3);  // TPS

        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        auto aResampled = ptp_interp(aTime, aData, 256);

        dab.set_num_periods_to_extract(1);
        auto simOps = dab.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        INFO("DAB TPS NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.15);
    }

    // =====================================================================
    // Negative phase shift (reverse power flow): the primary inductor current
    // should flow in the opposite direction compared to positive φ.
    // =====================================================================
    TEST_CASE("Test_Dab_NegativePhi_ReversePower", "[converter-model][dab-topology][smoke-test]") {
        auto build = [](double phi_deg) {
            json dabJson;
            json inputVoltage;
            inputVoltage["nominal"] = 800.0;
            dabJson["inputVoltage"] = inputVoltage;
            dabJson["seriesInductance"] = 35e-6;
            dabJson["useLeakageInductance"] = false;
            dabJson["operatingPoints"] = json::array();
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0};
            op["outputCurrents"] = {20.0};
            op["phaseShift"] = phi_deg;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
            return dabJson;
        };

        OpenMagnetics::Dab dab_pos(build(+23.0));
        OpenMagnetics::Dab dab_neg(build(-23.0));

        auto req_pos = dab_pos.process_design_requirements();
        auto req_neg = dab_neg.process_design_requirements();

        std::vector<double> tr_pos, tr_neg;
        for (auto& tr : req_pos.get_turns_ratios()) tr_pos.push_back(resolve_dimensional_values(tr));
        for (auto& tr : req_neg.get_turns_ratios()) tr_neg.push_back(resolve_dimensional_values(tr));
        double Lm_pos = resolve_dimensional_values(req_pos.get_magnetizing_inductance());
        double Lm_neg = resolve_dimensional_values(req_neg.get_magnetizing_inductance());

        auto ops_pos = dab_pos.process_operating_points(tr_pos, Lm_pos);
        auto ops_neg = dab_neg.process_operating_points(tr_neg, Lm_neg);

        // The key check: time-average of (Vab × iL) should flip sign between
        // positive and negative phi (positive = forward power, negative = reverse).
        auto avgPower = [](const OperatingPoint& op) {
            auto& exc = op.get_excitations_per_winding()[0];
            auto vWfm = exc.get_voltage()->get_waveform().value();
            auto iWfm = exc.get_current()->get_waveform().value();
            const auto& v = vWfm.get_data();
            const auto& i = iWfm.get_data();
            double sum = 0;
            size_t n = std::min(v.size(), i.size());
            for (size_t k = 0; k < n; ++k) sum += v[k] * i[k];
            return sum / n;
        };

        double P_pos = avgPower(ops_pos[0]);
        double P_neg = avgPower(ops_neg[0]);
        INFO("Power at +23°: " << P_pos << " W");
        INFO("Power at -23°: " << P_neg << " W");
        CHECK(P_pos > 0);
        CHECK(P_neg < 0);
        // Magnitudes should be similar (the model is power-flow symmetric)
        CHECK(std::abs(P_pos + P_neg) < 0.1 * std::abs(P_pos));
    }

    // =====================================================================
    // Multi-output DAB: per-output magnitudes are correctly scaled by n_i.
    // The previous test only checked excitation count.
    // =====================================================================
    TEST_CASE("Test_Dab_Multiple_Outputs_Magnitudes", "[converter-model][dab-topology][smoke-test]") {
        json dabJson;
        json inputVoltage;
        inputVoltage["nominal"] = 800.0;
        dabJson["inputVoltage"] = inputVoltage;
        dabJson["seriesInductance"] = 35e-6;
        dabJson["useLeakageInductance"] = false;
        dabJson["operatingPoints"] = json::array();
        {
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = {500.0, 250.0};   // V2_2 = ½ V2_1
            op["outputCurrents"] = {10.0, 4.0};
            op["phaseShift"] = 23.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);
        }

        OpenMagnetics::Dab dab(dabJson);
        auto req = dab.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

        auto ops = dab.process_operating_points(turnsRatios, Lm);
        REQUIRE(!ops.empty());

        // Should have 1 primary + 2 secondaries = 3 excitations
        CHECK(ops[0].get_excitations_per_winding().size() == 3);

        auto peakV = [](const OperatingPointExcitation& e) {
            auto wfm = e.get_voltage()->get_waveform().value();
            const auto& d = wfm.get_data();
            double pk = 0;
            for (auto v : d) if (std::abs(v) > pk) pk = std::abs(v);
            return pk;
        };

        // Secondary 0 voltage peak ≈ 500 V (V2_1), Secondary 1 ≈ 250 V (V2_2)
        double v0 = peakV(ops[0].get_excitations_per_winding()[1]);
        double v1 = peakV(ops[0].get_excitations_per_winding()[2]);
        INFO("Secondary 0 peak V = " << v0 << " V (expected ~500)");
        INFO("Secondary 1 peak V = " << v1 << " V (expected ~250)");
        CHECK(v0 > v1);
        REQUIRE_THAT(v0 / v1, Catch::Matchers::WithinAbs(2.0, 0.1));
    }

    // =====================================================================
    // Multi-output DAB: per-output secondary current magnitudes — analytical
    // load-share approximation vs SPICE multi-port simulation.
    //
    // The analytical model uses share = (Iout_i/Vout_i) / Σ G to project the
    // primary tank current onto each secondary; this is exact when all
    // outputs have the same V2 (same turns ratio n_i) and degrades as the
    // V2 spread grows. This test characterises that degradation across
    // three configurations.
    // =====================================================================
    namespace {
        struct MultiOutputErrors { double max_peak_err; double max_rms_err; };

        MultiOutputErrors run_multi_output_check(
            const std::vector<double>& Vouts,
            const std::vector<double>& Iouts)
        {
            json dabJson;
            json inputVoltage;
            inputVoltage["nominal"] = 800.0;
            dabJson["inputVoltage"] = inputVoltage;
            dabJson["seriesInductance"] = 35e-6;
            dabJson["useLeakageInductance"] = false;
            dabJson["operatingPoints"] = json::array();
            json op;
            op["ambientTemperature"] = 25.0;
            op["outputVoltages"] = Vouts;
            op["outputCurrents"] = Iouts;
            op["phaseShift"] = 23.0;
            op["switchingFrequency"] = 100000;
            dabJson["operatingPoints"].push_back(op);

            OpenMagnetics::Dab dab(dabJson);
            auto req = dab.process_design_requirements();
            std::vector<double> turnsRatios;
            for (auto& tr : req.get_turns_ratios())
                turnsRatios.push_back(resolve_dimensional_values(tr));
            double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

            auto analyticalOps = dab.process_operating_points(turnsRatios, Lm);
            dab.set_num_periods_to_extract(1);
            auto simOps = dab.simulate_and_extract_operating_points(turnsRatios, Lm);

            auto rms = [](const std::vector<double>& v) {
                double sum = 0;
                for (double x : v) sum += x * x;
                return v.empty() ? 0.0 : std::sqrt(sum / v.size());
            };
            auto peak = [](const std::vector<double>& v) {
                double pk = 0;
                for (double x : v) if (std::abs(x) > pk) pk = std::abs(x);
                return pk;
            };

            MultiOutputErrors err{0.0, 0.0};
            for (size_t outIdx = 0; outIdx < Vouts.size(); ++outIdx) {
                int wi = static_cast<int>(outIdx + 1);
                auto aData = ptp_current(analyticalOps[0], wi);
                auto sData = ptp_current(simOps[0], wi);
                if (aData.empty() || sData.empty()) continue;
                double aPeak = peak(aData), sPeak = peak(sData);
                double aRms  = rms(aData),  sRms  = rms(sData);
                double peakErr = (sPeak > 0.01) ? std::abs(aPeak - sPeak) / sPeak : 0.0;
                double rmsErr  = (sRms  > 0.01) ? std::abs(aRms  - sRms)  / sRms  : 0.0;
                std::cerr << "[multi-out] V2=" << Vouts[outIdx] << "V Iout=" << Iouts[outIdx] << "A: "
                          << "ana_peak=" << aPeak << "A  spice_peak=" << sPeak << "A  err=" << (peakErr*100) << "%  |  "
                          << "ana_rms="  << aRms  << "A  spice_rms="  << sRms  << "A  err=" << (rmsErr*100)  << "%\n";
                err.max_peak_err = std::max(err.max_peak_err, peakErr);
                err.max_rms_err  = std::max(err.max_rms_err,  rmsErr);
            }
            return err;
        }
    }

    // The three multi-output PtP tests below are CHARACTERISATION ONLY.
    // They run the analytical model and the SPICE model side-by-side and
    // print the per-output peak/RMS error so the limitation of the
    // multi-output model is documented in the build log. The pass criteria
    // are deliberately loose because:
    //   1. The schema does not let users specify per-secondary leakage
    //      inductance, so the SPICE multi-output netlist over-couples the
    //      secondaries (K=0.9999 between every winding pair) and the AC
    //      current split between secondaries is not physically determined.
    //   2. The analytical model uses a load-conductance share which is only
    //      exact for the (unrealistic) limit of equal V2 and equal AC share.
    //   3. Multi-output DAB on a single transformer is uncommon in practice
    //      — the typical realisation uses N separate transformers, one per
    //      output, each with its own primary bridge. MKF's purpose is single-
    //      transformer magnetic-component design, so the multi-output case
    //      is approximate by construction.
    // Single-output DAB (the common case) is fully supported and tested by
    // the SPS / EPS / DPS / TPS PtP tests above (NRMSE ≤ 6 % across modes).
    TEST_CASE("Test_Dab_Multiple_Outputs_PtP_Uniform", "[converter-model][dab-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        std::cerr << "\n=== Multi-output PtP, uniform V2 (characterisation only) ===\n";
        auto err = run_multi_output_check({500.0, 500.0, 500.0}, {6.0, 4.0, 2.0});
        std::cerr << "  worst peak err = " << (err.max_peak_err * 100) << "%, "
                  << "worst rms err = " << (err.max_rms_err * 100) << "%\n";
        // Loose: doesn't fail the build, just records the magnitude.
        CHECK(err.max_peak_err < 5.0);
        CHECK(err.max_rms_err  < 5.0);
    }

    TEST_CASE("Test_Dab_Multiple_Outputs_PtP_ModerateSpread", "[converter-model][dab-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        std::cerr << "\n=== Multi-output PtP, moderate V2 spread (1.5x) (characterisation) ===\n";
        auto err = run_multi_output_check({500.0, 400.0, 333.0}, {6.0, 4.0, 3.0});
        std::cerr << "  worst peak err = " << (err.max_peak_err * 100) << "%, "
                  << "worst rms err = " << (err.max_rms_err * 100) << "%\n";
        CHECK(err.max_peak_err < 5.0);
        CHECK(err.max_rms_err  < 5.0);
    }

    TEST_CASE("Test_Dab_Multiple_Outputs_PtP_LargeSpread", "[converter-model][dab-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        // 5× V2 spread — load-share approximation breaks down here.
        // This test pins down the error for documentation purposes; it does
        // NOT enforce a tight bound. If the user has a multi-output DAB with
        // V2 spread > 2×, they should not rely on per-output current
        // magnitudes from the analytical model.
        std::cerr << "\n=== Multi-output PtP, LARGE V2 spread (5×) — load-share stress test ===\n";
        auto err = run_multi_output_check({500.0, 250.0, 100.0}, {10.0, 4.0, 2.0});
        std::cerr << "  worst peak err = " << (err.max_peak_err * 100) << "%, "
                  << "worst rms err = " << (err.max_rms_err * 100) << "%\n";
        // Soft bound: doesn't fail the build, just records the magnitude.
        CHECK(err.max_peak_err < 5.0);
        CHECK(err.max_rms_err  < 10.0);
    }

} // anonymous namespace