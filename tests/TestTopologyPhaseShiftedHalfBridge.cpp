#include <source_location>
#include "support/Painter.h"
#include "converter_models/PhaseShiftedHalfBridge.h"
#include "converter_models/PhaseShiftedFullBridge.h"
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

using namespace OpenMagnetics;

namespace {

inline std::vector<double> ptp_interp(const std::vector<double>& t,
                                      const std::vector<double>& d, int N) {
    std::vector<double> out(N);
    if (t.empty()) return out;
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

inline double ptp_nrmse(const std::vector<double>& ref,
                        const std::vector<double>& sim) {
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

inline std::vector<double> ptp_current(const OperatingPoint& op, size_t wi = 0) {
    if (wi >= op.get_excitations_per_winding().size()) return {};
    auto& e = op.get_excitations_per_winding()[wi];
    if (!e.get_current() || !e.get_current()->get_waveform()) return {};
    return e.get_current()->get_waveform()->get_data();
}

inline std::vector<double> ptp_current_time(const OperatingPoint& op, size_t wi = 0) {
    if (wi >= op.get_excitations_per_winding().size()) return {};
    auto& e = op.get_excitations_per_winding()[wi];
    if (!e.get_current() || !e.get_current()->get_waveform()) return {};
    auto tv = e.get_current()->get_waveform()->get_time();
    return tv ? tv.value() : std::vector<double>{};
}

auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}
    .parent_path().append("..").append("output");


// =========================================================================
// Helper: create a typical PSHB JSON (400V -> 12V, 25A, 100 kHz, center-tapped)
// Note: PSHB is typically lower power than PSFB due to halved primary voltage
// =========================================================================
json make_pshb_json(double Vin_nom = 400.0, double Vin_min = 370.0, double Vin_max = 410.0,
                    double Vo = 12.0, double Io = 25.0, double Fs = 100000.0,
                    double phaseShift = 135.0,
                    std::string rectType = "centerTapped") {
    json j;
    json inputVoltage;
    inputVoltage["nominal"] = Vin_nom;
    inputVoltage["minimum"] = Vin_min;
    inputVoltage["maximum"] = Vin_max;
    j["inputVoltage"] = inputVoltage;
    j["rectifierType"] = rectType;

    j["operatingPoints"] = json::array();
    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {Vo};
        op["outputCurrents"] = {Io};
        op["switchingFrequency"] = Fs;
        op["phaseShift"] = phaseShift;
        j["operatingPoints"].push_back(op);
    }
    return j;
}


// =========================================================================
// TEST 1: PSHB Design - Center Tapped, 400V -> 12V, 300W
// =========================================================================
TEST_CASE("Test_Pshb_CenterTapped_Design", "[converter-model][pshb-topology][smoke-test]") {
    auto pshbJson = make_pshb_json();
    OpenMagnetics::Pshb pshb(pshbJson);

    SECTION("Input validation") {
        CHECK(pshb.run_checks(false) == true);
    }

    SECTION("Bridge voltage factor is 0.5") {
        CHECK(pshb.get_bridge_voltage_factor() == 0.5);
    }

    SECTION("Turns ratio - half of PSFB for same conditions") {
        auto req = pshb.process_design_requirements();
        double n = resolve_dimensional_values(req.get_turns_ratios()[0]);
        // n = (Vin/2) * Deff / (Vo + Vd)
        // Deff = 135/180 = 0.75
        // n = (400/2) * 0.75 / (12 + 0.6) = 200*0.75 / 12.6 = 150/12.6 ≈ 11.905
        double expectedN = (400.0 / 2.0) * (135.0 / 180.0) / (12.0 + 0.6);
        REQUIRE_THAT(n, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
        CHECK(n > 0);
    }

    SECTION("Turns ratio is roughly half of equivalent PSFB") {
        // PSFB with same conditions: n_fb = Vin * Deff / (Vo + Vd)
        // PSHB: n_hb = (Vin/2) * Deff / (Vo + Vd) = n_fb / 2
        auto req = pshb.process_design_requirements();
        double n_hb = resolve_dimensional_values(req.get_turns_ratios()[0]);
        double n_fb_expected = 400.0 * (135.0 / 180.0) / (12.0 + 0.6);
        REQUIRE_THAT(n_hb, Catch::Matchers::WithinAbs(n_fb_expected / 2.0, n_fb_expected * 0.05));
    }

    SECTION("Magnetizing inductance is positive and reasonable") {
        auto req = pshb.process_design_requirements();
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
        CHECK(Lm > 0);
        CHECK(Lm > 10e-6);
        CHECK(Lm < 500e-3);
    }

    SECTION("Computed series inductance is positive") {
        pshb.process_design_requirements();
        double Lr = pshb.get_computed_series_inductance();
        CHECK(Lr > 0);
        CHECK(Lr < 1e-3);
    }

    SECTION("Computed output inductance is positive") {
        pshb.process_design_requirements();
        double Lo = pshb.get_computed_output_inductance();
        CHECK(Lo > 0);
    }

    SECTION("Effective duty cycle (post duty-cycle-loss)") {
        pshb.process_design_requirements();
        double Deff = pshb.get_computed_effective_duty_cycle();
        CHECK(Deff > 0);
        CHECK(Deff < 1.0);
        // D_cmd = 135°/180° = 0.75; PSHB uses Vbus = Vin/2, so the same
        // turns ratio iteration produces a slightly larger ΔD than PSFB.
        REQUIRE_THAT(Deff, Catch::Matchers::WithinAbs(0.69, 0.07));
    }
}


// =========================================================================
// TEST 2: PSHB Operating Points Generation
// =========================================================================
TEST_CASE("Test_Pshb_OperatingPoints_Generation", "[converter-model][pshb-topology][smoke-test]") {
    auto pshbJson = make_pshb_json();
    OpenMagnetics::Pshb pshb(pshbJson);
    auto req = pshb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    SECTION("Multiple input voltages") {
        auto ops = pshb.process_operating_points(turnsRatios, Lm);
        CHECK(ops.size() == 3);
    }

    SECTION("Waveform structure") {
        auto ops = pshb.process_operating_points(turnsRatios, Lm);
        REQUIRE(!ops.empty());

        auto& op = ops[0];
        // Default rectifier is centerTapped → Pri + Sec0a + Sec0b = 3 windings.
        CHECK(op.get_excitations_per_winding().size() == 3);

        auto& priExc = op.get_excitations_per_winding()[0];
        REQUIRE(priExc.get_current().has_value());
        REQUIRE(priExc.get_voltage().has_value());
        CHECK(priExc.get_frequency() == 100e3);

        auto currentWfm = priExc.get_current()->get_waveform().value();
        CHECK(currentWfm.get_data().size() == 1024);

        auto voltageWfm = priExc.get_voltage()->get_waveform().value();
        CHECK(voltageWfm.get_data().size() == 1024);
    }

    SECTION("Primary voltage is 3-level at HALF amplitude") {
        auto ops = pshb.process_operating_points(turnsRatios, Lm);
        auto& priExc = ops[0].get_excitations_per_winding()[0];
        auto vData = priExc.get_voltage()->get_waveform().value().get_data();

        double vMax = *std::max_element(vData.begin(), vData.end());
        double vMin = *std::min_element(vData.begin(), vData.end());

        // Primary voltage should swing to approximately +(Vin/2) and -(Vin/2)
        double Vin_min = 370.0;
        double Vhb_expected = Vin_min / 2.0;  // 185V
        REQUIRE_THAT(vMax, Catch::Matchers::WithinAbs(Vhb_expected, Vhb_expected * 0.05));
        REQUIRE_THAT(vMin, Catch::Matchers::WithinAbs(-Vhb_expected, Vhb_expected * 0.05));

        // Should have zero-voltage intervals (freewheeling)
        int zeroCount = 0;
        for (double v : vData) {
            if (std::abs(v) < 1.0) zeroCount++;
        }
        CHECK(zeroCount > 0);
    }

    SECTION("PSHB voltage is half of equivalent PSFB voltage") {
        auto ops = pshb.process_operating_points(turnsRatios, Lm);
        auto& priExc = ops[0].get_excitations_per_winding()[0];
        auto vData = priExc.get_voltage()->get_waveform().value().get_data();
        double vMax = *std::max_element(vData.begin(), vData.end());

        // For the min input voltage (first OP = Vin_min = 370V)
        // Half-bridge peak voltage should be ~370/2 = 185V
        double Vin_min = 370.0;
        CHECK(vMax < Vin_min);  // Must be less than full Vin
        REQUIRE_THAT(vMax, Catch::Matchers::WithinAbs(Vin_min / 2.0, Vin_min * 0.05));
    }

    SECTION("Primary current antisymmetry") {
        auto ops = pshb.process_operating_points(turnsRatios, Lm);
        auto& priExc = ops[0].get_excitations_per_winding()[0];
        auto iData = priExc.get_current()->get_waveform().value().get_data();

        size_t N = iData.size();
        size_t half = N / 2;
        double asymmetrySum = 0;
        for (size_t k = 1; k < half; k++) {
            asymmetrySum += std::abs(iData[k] + iData[half + k]);
        }
        double avgAsymmetry = asymmetrySum / half;
        double iPeak = *std::max_element(iData.begin(), iData.end());
        CHECK(avgAsymmetry / iPeak < 0.05);
    }

    SECTION("Secondary winding excitation exists") {
        auto ops = pshb.process_operating_points(turnsRatios, Lm);
        auto& secExc = ops[0].get_excitations_per_winding()[1];
        REQUIRE(secExc.get_current().has_value());
        REQUIRE(secExc.get_voltage().has_value());
    }
}


// =========================================================================
// TEST 3: PSHB Waveform Plotting
// =========================================================================
TEST_CASE("Test_Pshb_Waveform_Plotting", "[converter-model][pshb-topology][visualization]") {
    auto pshbJson = make_pshb_json();
    OpenMagnetics::Pshb pshb(pshbJson);
    auto req = pshb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto ops = pshb.process_operating_points(turnsRatios, Lm);
    REQUIRE(!ops.empty());

    SECTION("Primary current waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Pshb_Primary_Current_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Primary voltage waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Pshb_Primary_Voltage_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Secondary current waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Pshb_Secondary_Current_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
        painter.export_svg();
    }
}


// =========================================================================
// TEST 4: PSHB SPICE Netlist Generation
// =========================================================================
TEST_CASE("Test_Pshb_Spice_Netlist", "[converter-model][pshb-topology][spice]") {
    auto pshbJson = make_pshb_json();
    OpenMagnetics::Pshb pshb(pshbJson);
    // Netlist string assertions in this test target the SW1 + NPC clamp
    // form (S1..S4, DC1, DC2). Pin SW1 mode so the test isn't sensitive
    // to the BEHAVIORAL_PULSE default introduced for ngspice WASM perf.
    pshb.set_bridge_simulation_mode(OpenMagnetics::BridgeSimulationMode::VOLTAGE_CONTROLLED_SWITCH);
    auto req = pshb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    SECTION("Netlist is non-empty and contains key elements") {
        std::string netlist = pshb.generate_ngspice_circuit(turnsRatios, Lm);
        CHECK(!netlist.empty());
        CHECK(netlist.find("Phase-Shifted Half Bridge") != std::string::npos);
        CHECK(netlist.find("Pinheiro-Barbi") != std::string::npos);
        CHECK(netlist.find("L_pri") != std::string::npos);
        CHECK(netlist.find("L_sec") != std::string::npos);
        CHECK(netlist.find("L_out") != std::string::npos);
        CHECK(netlist.find("R_load") != std::string::npos);
        CHECK(netlist.find(".tran") != std::string::npos);
        // Single 3-level NPC leg: 4 stacked S-switches (S1..S4) with body diodes
        CHECK(netlist.find(".model SW1 SW") != std::string::npos);
        CHECK(netlist.find(".model DIDEAL D") != std::string::npos);
        // §8a.5: S1 now drains from qa_drain (downstream of the Vq1_sense
        // 0-V ammeter) rather than directly from vin_dc.
        CHECK(netlist.find("Vq1_sense vin_dc qa_drain") != std::string::npos);
        CHECK(netlist.find("S1 qa_drain nH pwm_S1") != std::string::npos);
        CHECK(netlist.find("S2 nH bridge_a pwm_S2") != std::string::npos);
        CHECK(netlist.find("S3 bridge_a nL pwm_S3") != std::string::npos);
        CHECK(netlist.find("S4 nL 0 pwm_S4") != std::string::npos);
        CHECK(netlist.find("D1 nH vin_dc DIDEAL") != std::string::npos);
        CHECK(netlist.find("D4 0 nL DIDEAL") != std::string::npos);
        // Two NPC clamp diodes anchoring inner nodes to mid_cap
        CHECK(netlist.find("DC1 mid_cap nH DIDEAL") != std::string::npos);
        CHECK(netlist.find("DC2 nL mid_cap DIDEAL") != std::string::npos);
        // Per-switch RC snubbers (DAB/PSFB recipe)
        CHECK(netlist.find("Rsnub_S1") != std::string::npos);
        CHECK(netlist.find("Csnub_S4") != std::string::npos);
        // 4 PWM gate sources (one per switch)
        CHECK(netlist.find("Vpwm_S1 pwm_S1 0 PULSE") != std::string::npos);
        CHECK(netlist.find("Vpwm_S4 pwm_S4 0 PULSE") != std::string::npos);
        CHECK(netlist.find("Vpri_sense bridge_a pri_lr") != std::string::npos);
        // Split capacitor providing Vin/2 reference (point 'b' for vab)
        CHECK(netlist.find("C_split_hi") != std::string::npos);
        CHECK(netlist.find("C_split_lo") != std::string::npos);
        // Robustness features
        CHECK(netlist.find("METHOD=GEAR") != std::string::npos);
        // Differential probe vab = v(bridge_a) - v(mid_cap)
        CHECK(netlist.find("Evab vab 0 bridge_a mid_cap 1") != std::string::npos);
    }

    SECTION("Netlist saved to file") {
        std::string netlist = pshb.generate_ngspice_circuit(turnsRatios, Lm);
        auto outFile = outputFilePath;
        outFile.append("Test_Pshb_Netlist.cir");
        std::ofstream ofs(outFile);
        ofs << netlist;
        ofs.close();
        CHECK(std::filesystem::exists(outFile));
    }
}

// =========================================================================
// REGRESSION: SPICE_PROBE_HANDOFF.md Bug 2 — Vq1_sense must sit between
// `vin_dc` and a switch drain in BEHAVIORAL_PULSE mode, not between the
// behavioral leg source and the bridge node.
//
// Status at time of writing (commit 38e3947b):
//   BEHAVIORAL_PULSE branch emits `Vq1_sense leg_a_src bridge_a 0`,
//   which measures the behavioral source's leg current and misses the
//   DC-link split-cap mid-point currents. The half-bridge offset makes
//   avg(input_current) drift further from output_power / Vin than in
//   the PSFB case.
//
// `[!shouldfail]` until PhaseShiftedHalfBridge.cpp:~729-741 is fixed.
// =========================================================================
TEST_CASE("Test_Pshb_Spice_BehavioralPulse_Vq1Sense_On_VinDc",
          "[converter-model][pshb-topology][spice]") {
    auto pshbJson = make_pshb_json();
    OpenMagnetics::Pshb pshb(pshbJson);
    pshb.set_bridge_simulation_mode(OpenMagnetics::BridgeSimulationMode::BEHAVIORAL_PULSE);
    auto req = pshb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    std::string netlist = pshb.generate_ngspice_circuit(turnsRatios, Lm);

    // Per the SW1 reference at PhaseShiftedHalfBridge.cpp:~772, the correct
    // form is `Vq1_sense vin_dc qa_drain 0` — the ammeter must sit on the
    // DC rail, not between the behavioral source and the bridge node.
    INFO("netlist excerpt around Vq1_sense:");
    auto pos = netlist.find("Vq1_sense");
    if (pos != std::string::npos) {
        INFO(netlist.substr(pos, 80));
    }
    REQUIRE(netlist.find("Vq1_sense vin_dc") != std::string::npos);
    REQUIRE(netlist.find("Vq1_sense leg_a_src") == std::string::npos);
}


// =========================================================================
// TEST 5: PSHB Multiple Outputs
// =========================================================================
TEST_CASE("Test_Pshb_Multiple_Outputs", "[converter-model][pshb-topology][smoke-test]") {
    json pshbJson;
    json inputVoltage;
    inputVoltage["nominal"] = 400.0;
    pshbJson["inputVoltage"] = inputVoltage;
    pshbJson["rectifierType"] = "centerTapped";
    pshbJson["operatingPoints"] = json::array();

    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0, 5.0};
        op["outputCurrents"] = {25.0, 5.0};
        op["switchingFrequency"] = 100000;
        op["phaseShift"] = 135.0;
        pshbJson["operatingPoints"].push_back(op);
    }

    OpenMagnetics::Pshb pshb(pshbJson);
    auto req = pshb.process_design_requirements();

    CHECK(req.get_turns_ratios().size() == 2);

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto ops = pshb.process_operating_points(turnsRatios, Lm);
    REQUIRE(!ops.empty());

    // Default CT rectifier: Pri + 2 outputs * 2 CT half-windings = 5
    CHECK(ops[0].get_excitations_per_winding().size() == 5);
}


// =========================================================================
// TEST 6: PSHB Static calculations
// =========================================================================
TEST_CASE("Test_Pshb_Static_Calculations", "[converter-model][pshb-topology][unit]") {

    SECTION("Effective duty cycle") {
        CHECK(Pshb::compute_effective_duty_cycle(0.0) == 0.0);
        CHECK(Pshb::compute_effective_duty_cycle(90.0) == 0.5);
        CHECK(Pshb::compute_effective_duty_cycle(180.0) == 1.0);
        REQUIRE_THAT(Pshb::compute_effective_duty_cycle(135.0),
                     Catch::Matchers::WithinAbs(0.75, 1e-6));
    }

    SECTION("Output voltage - center tapped (includes Vin/2 factor)") {
        // Vo = (Vin/2) * Deff / n - Vd
        double Vo = Pshb::compute_output_voltage(400.0, 0.75, 11.0, 0.6,
                                                  BRectifierType::CENTER_TAPPED);
        double expected = (400.0 / 2.0) * 0.75 / 11.0 - 0.6;
        REQUIRE_THAT(Vo, Catch::Matchers::WithinAbs(expected, 0.01));
    }

    SECTION("Output voltage - full bridge rectifier") {
        double Vo = Pshb::compute_output_voltage(400.0, 0.75, 11.0, 0.6,
                                                  BRectifierType::FULL_BRIDGE);
        double expected = (400.0 / 2.0) * 0.75 / 11.0 - 2.0 * 0.6;
        REQUIRE_THAT(Vo, Catch::Matchers::WithinAbs(expected, 0.01));
    }

    SECTION("Turns ratio round-trip") {
        double n = Pshb::compute_turns_ratio(400.0, 12.0, 0.75, 0.6,
                                              BRectifierType::CENTER_TAPPED);
        double Vo_check = Pshb::compute_output_voltage(400.0, 0.75, n, 0.6,
                                                        BRectifierType::CENTER_TAPPED);
        REQUIRE_THAT(Vo_check, Catch::Matchers::WithinAbs(12.0, 0.01));
    }

    SECTION("PSHB turns ratio is half of PSFB turns ratio") {
        // For same Vin, Deff, Vo, Vd, rectType:
        // n_fb = Vin * Deff / (Vo+Vd)
        // n_hb = (Vin/2) * Deff / (Vo+Vd) = n_fb / 2
        double n_hb = Pshb::compute_turns_ratio(400.0, 12.0, 0.7, 0.6,
                                                  BRectifierType::CENTER_TAPPED);
        double n_fb_expected = 400.0 * 0.7 / (12.0 + 0.6);
        REQUIRE_THAT(n_hb, Catch::Matchers::WithinAbs(n_fb_expected / 2.0, 0.01));
    }
}


// =========================================================================
// TEST 7: PSHB Current Doubler rectifier
// =========================================================================
TEST_CASE("Test_Pshb_CurrentDoubler_Design", "[converter-model][pshb-topology][smoke-test]") {
    auto pshbJson = make_pshb_json(400.0, 370.0, 410.0, 12.0, 25.0, 100000.0,
                                   135.0, "currentDoubler");
    OpenMagnetics::Pshb pshb(pshbJson);
    auto req = pshb.process_design_requirements();

    double n = resolve_dimensional_values(req.get_turns_ratios()[0]);
    // Current doubler: n = (Vin/2) * Deff / (2*(Vo+Vd))
    double expectedN = (400.0 / 2.0) * 0.75 / (2.0 * (12.0 + 0.6));
    REQUIRE_THAT(n, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
}


// =========================================================================
// TEST 8: AdvancedPshb JSON round-trip
// =========================================================================
TEST_CASE("Test_AdvancedPshb_Process", "[converter-model][pshb-topology][advanced]") {
    json advJson;
    json inputVoltage;
    inputVoltage["nominal"] = 400.0;
    inputVoltage["minimum"] = 370.0;
    inputVoltage["maximum"] = 410.0;
    advJson["inputVoltage"] = inputVoltage;
    advJson["rectifierType"] = "centerTapped";
    advJson["operatingPoints"] = json::array();

    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0};
        op["outputCurrents"] = {25.0};
        op["switchingFrequency"] = 100000;
        op["phaseShift"] = 135.0;
        advJson["operatingPoints"].push_back(op);
    }

    advJson["desiredTurnsRatios"] = {11.0};
    advJson["desiredMagnetizingInductance"] = 500e-6;

    OpenMagnetics::AdvancedPshb advPshb(advJson);
    auto inputs = advPshb.process();

    CHECK(inputs.get_design_requirements().get_turns_ratios().size() == 1);
    double n = resolve_dimensional_values(inputs.get_design_requirements().get_turns_ratios()[0]);
    REQUIRE_THAT(n, Catch::Matchers::WithinAbs(11.0, 0.01));

    double Lm = resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance());
    REQUIRE_THAT(Lm, Catch::Matchers::WithinAbs(500e-6, 1e-7));

    CHECK(!inputs.get_operating_points().empty());
}


// =========================================================================
// TEST 9: Comparison PSHB vs PSFB for same operating conditions
// =========================================================================
TEST_CASE("Test_Pshb_vs_Psfb_Comparison", "[converter-model][pshb-topology][psfb-topology][comparison]") {
    // Both converters: 400V -> 12V, 25A, 100kHz, 135 deg phase shift, center-tapped
    auto pshbJson = make_pshb_json(400.0, 370.0, 410.0, 12.0, 25.0, 100000.0,
                                   135.0, "centerTapped");

    // PSFB with identical operating conditions
    json psfbJson = pshbJson;  // Same JSON (both use PhaseShiftFullBridge schema)

    OpenMagnetics::Pshb pshb(pshbJson);
    OpenMagnetics::Psfb psfb(psfbJson);

    auto reqHb = pshb.process_design_requirements();
    auto reqFb = psfb.process_design_requirements();

    double n_hb = resolve_dimensional_values(reqHb.get_turns_ratios()[0]);
    double n_fb = resolve_dimensional_values(reqFb.get_turns_ratios()[0]);

    SECTION("PSHB turns ratio is half of PSFB") {
        REQUIRE_THAT(n_hb, Catch::Matchers::WithinAbs(n_fb / 2.0, n_fb * 0.02));
    }

    SECTION("PSHB primary voltage amplitude is half of PSFB") {
        std::vector<double> tr_hb, tr_fb;
        for (auto& tr : reqHb.get_turns_ratios())
            tr_hb.push_back(resolve_dimensional_values(tr));
        for (auto& tr : reqFb.get_turns_ratios())
            tr_fb.push_back(resolve_dimensional_values(tr));

        double Lm_hb = resolve_dimensional_values(reqHb.get_magnetizing_inductance());
        double Lm_fb = resolve_dimensional_values(reqFb.get_magnetizing_inductance());

        auto opsHb = pshb.process_operating_points(tr_hb, Lm_hb);
        auto opsFb = psfb.process_operating_points(tr_fb, Lm_fb);

        // Compare at nominal voltage (middle OP if 3 exist, or index 0)
        size_t idxHb = (opsHb.size() > 1) ? 1 : 0;
        size_t idxFb = (opsFb.size() > 1) ? 1 : 0;

        auto vHb = opsHb[idxHb].get_excitations_per_winding()[0]
                       .get_voltage()->get_waveform().value().get_data();
        auto vFb = opsFb[idxFb].get_excitations_per_winding()[0]
                       .get_voltage()->get_waveform().value().get_data();

        double vMaxHb = *std::max_element(vHb.begin(), vHb.end());
        double vMaxFb = *std::max_element(vFb.begin(), vFb.end());

        // HB peak should be ~half of FB peak
        REQUIRE_THAT(vMaxHb, Catch::Matchers::WithinAbs(vMaxFb / 2.0, vMaxFb * 0.05));
    }
}


// =========================================================================
// TEST: PSHB duty-cycle loss diagnostic (Sabate's formula with Vbus = Vin/2)
// =========================================================================
TEST_CASE("Test_Pshb_DutyCycleLoss", "[converter-model][pshb-topology][unit]") {
    json pshbJson;
    json inputVoltage; inputVoltage["nominal"] = 400.0;
    pshbJson["inputVoltage"] = inputVoltage;
    pshbJson["rectifierType"] = "fullBridge";
    pshbJson["seriesInductance"] = 5e-6;
    pshbJson["operatingPoints"] = json::array();
    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0};
        op["outputCurrents"] = {50.0};
        op["switchingFrequency"] = 100000.0;
        op["phaseShift"] = 126.0;
        pshbJson["operatingPoints"].push_back(op);
    }
    OpenMagnetics::Pshb pshb(pshbJson);

    auto req = pshb.process_design_requirements();
    std::vector<double> tr;
    for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto ops = pshb.process_operating_points(tr, Lm);
    REQUIRE(!ops.empty());

    double dcl = pshb.get_last_duty_cycle_loss();
    double n = tr[0];
    // Vbus = Vin/2, so ΔD scales accordingly.
    double Vbus = 400.0 / 2.0;
    double dcl_expected = 4.0 * 5e-6 * 50.0 * 100e3 / (n * Vbus);
    REQUIRE_THAT(dcl, Catch::Matchers::WithinRel(dcl_expected, 0.05));
    CHECK(pshb.get_last_effective_duty_cycle() < 0.7);
    CHECK(pshb.get_last_effective_duty_cycle() > 0.0);
}


// =========================================================================
// TEST: PSHB ZVS diagnostic — non-zero margin and resonant time
// =========================================================================
TEST_CASE("Test_Pshb_ZvsDiagnostics", "[converter-model][pshb-topology][unit]") {
    auto pshbJson = make_pshb_json();
    pshbJson["seriesInductance"] = 50e-6;
    OpenMagnetics::Pshb pshb(pshbJson);
    auto req = pshb.process_design_requirements();
    std::vector<double> tr;
    for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
    pshb.process_operating_points(tr,
        resolve_dimensional_values(req.get_magnetizing_inductance()));

    // With 50 µH, ZVS margin should be comfortably positive.
    CHECK(pshb.get_last_zvs_margin_lagging() > 0.0);
    CHECK(pshb.get_last_resonant_transition_time() > 0.0);
    CHECK(pshb.get_last_resonant_transition_time() < 1e-6);  // sub-µs
    CHECK(pshb.get_last_primary_peak_current() > 0.0);
}


// =========================================================================
// TEST: PSHB analytical vs NgSpice point-to-point
// =========================================================================
TEST_CASE("Test_Pshb_PtP_AnalyticalVsNgspice",
          "[converter-model][pshb-topology][ngspice-simulation][ptpcomparison]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    auto j = make_pshb_json(400.0, 380.0, 420.0, 12.0, 50.0, 100000.0,
                            126.0, "fullBridge");
    j["seriesInductance"] = 5e-6;
    OpenMagnetics::Pshb pshb(j);
    pshb.set_num_periods_to_extract(1);
    // Analytical model is calibrated against the SW1 + NPC clamp bridge.
    // Pin SW1 mode so the PtP NRMSE gate stays meaningful regardless of
    // the BEHAVIORAL_PULSE default.
    pshb.set_bridge_simulation_mode(OpenMagnetics::BridgeSimulationMode::VOLTAGE_CONTROLLED_SWITCH);

    auto req = pshb.process_design_requirements();
    std::vector<double> tr;
    for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto analytical = pshb.process_operating_points(tr, Lm);
    REQUIRE(!analytical.empty());
    auto aTime = ptp_current_time(analytical[0], 0);
    auto aData = ptp_current(analytical[0], 0);
    REQUIRE(!aData.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    auto sim = pshb.simulate_and_extract_operating_points(tr, Lm);
    REQUIRE(!sim.empty());
    auto sTime = ptp_current_time(sim[0], 0);
    auto sData = ptp_current(sim[0], 0);
    if (!sData.empty()) {
        auto sResampled = ptp_interp(sTime, sData, 256);
        double nrmse = ptp_nrmse(aResampled, sResampled);

        // Debug: print key stats
        double aMean = 0, aMax = -1e30, aMin = 1e30, aRms = 0;
        double sMean = 0, sMax = -1e30, sMin = 1e30, sRms = 0;
        for (double v : aResampled) {
            aMean += v; aMax = std::max(aMax, v); aMin = std::min(aMin, v); aRms += v*v;
        }
        for (double v : sResampled) {
            sMean += v; sMax = std::max(sMax, v); sMin = std::min(sMin, v); sRms += v*v;
        }
        aMean /= aResampled.size(); sMean /= sResampled.size();
        aRms = std::sqrt(aRms / aResampled.size()); sRms = std::sqrt(sRms / sResampled.size());
        INFO("Analytical: mean=" << aMean << "A, rms=" << aRms << "A, min=" << aMin << "A, max=" << aMax << "A");
        INFO("SPICE:      mean=" << sMean << "A, rms=" << sRms << "A, min=" << sMin << "A, max=" << sMax << "A");
        INFO("Design: Dcmd=" << pshb.get_computed_effective_duty_cycle()
            << ", Lr=" << (pshb.get_computed_series_inductance()*1e6) << "uH"
            << ", Lm=" << (resolve_dimensional_values(req.get_magnetizing_inductance())*1e6) << "uH"
            << ", n=" << tr[0]);

        INFO("PSHB primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        // DAB-quality threshold (Golden Guide §15: ≤ 0.15).  Currently
        // measured at ~0.112 with the physics-derived freewheel τ =
        // Lr / (2·RON_mosfet) replacing the old empirical 1.0·dur_fw.
        CHECK(nrmse < 0.15);
    }
}


// =========================================================================
// TEST: PSHB at light load (10 % of rated)
// =========================================================================
TEST_CASE("Test_Pshb_Light_Load", "[converter-model][pshb-topology][unit]") {
    auto j = make_pshb_json(400.0, 380.0, 420.0, 12.0, 2.5, 100000.0,
                            135.0, "fullBridge");
    j["seriesInductance"] = 5e-6;
    OpenMagnetics::Pshb pshb(j);
    auto req = pshb.process_design_requirements();
    std::vector<double> tr;
    for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
    pshb.process_operating_points(tr,
        resolve_dimensional_values(req.get_magnetizing_inductance()));

    double threshold = pshb.get_last_zvs_load_threshold();
    INFO("Lagging-leg ZVS minimum-load threshold (A): " << threshold);
    CHECK(threshold > 0.0);
    CHECK(2.5 < threshold);
    CHECK(pshb.get_last_resonant_transition_time() > 0.0);
}


// =========================================================================
// TEST: PSHB Lr / Lo round-trip
// =========================================================================
TEST_CASE("Test_Pshb_Inductance_Round_Trip", "[converter-model][pshb-topology][unit]") {
    auto j = make_pshb_json();
    const double user_Lr = 7.5e-6;
    j["seriesInductance"] = user_Lr;
    OpenMagnetics::Pshb pshb(j);
    pshb.process_design_requirements();
    REQUIRE_THAT(pshb.get_computed_series_inductance(),
                 Catch::Matchers::WithinRel(user_Lr, 1e-9));

    json j2 = make_pshb_json();
    const double user_Lo = 50e-6;
    j2["outputInductance"] = user_Lo;
    OpenMagnetics::Pshb pshb2(j2);
    pshb2.process_design_requirements();
    REQUIRE_THAT(pshb2.get_computed_output_inductance(),
                 Catch::Matchers::WithinRel(user_Lo, 1e-9));
}


// ────────────────────────────────────────────────────────────────────────────
// §5.1 converter-port DC-stream gate. See ConverterPortChecks for the full
// bound rationale. The signals returned by
// simulate_and_extract_topology_waveforms() are the DC source / DC filtered
// output rails — vab (bridge midpoint) and i(L_pri) must NEVER appear here.
// ────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pshb_ConverterPortWaveforms",
          "[converter-port-waveforms][pshb-topology][ngspice-simulation]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    const double Vin = 400.0, Vout = 12.0, Iout = 50.0;
    json j;
    json inputVoltage; inputVoltage["nominal"] = Vin;
    j["inputVoltage"] = inputVoltage;
    j["rectifierType"] = "fullBridge";
    j["seriesInductance"] = 5e-6;
    j["operatingPoints"] = json::array();
    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {Vout};
        op["outputCurrents"] = {Iout};
        op["switchingFrequency"] = 100000.0;
        op["phaseShift"] = 126.0;
        j["operatingPoints"].push_back(op);
    }

    OpenMagnetics::Pshb pshb(j);
    pshb.set_num_periods_to_extract(1);
    auto req = pshb.process_design_requirements();
    std::vector<double> tr;
    for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
    const double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto wfs = pshb.simulate_and_extract_topology_waveforms(tr, Lm);
    REQUIRE(!wfs.empty());

    // Same envelope as PSFB: phase-shift modulation reaches the primary
    // as a *reduced* effective duty (Lr + transformer K=0.999 freewheel
    // notches eat ~10-30 % of the volt-second product). The §5.1 gate's
    // job is to catch AC-in-DC-stream regressions, NOT to enforce
    // steady-state accuracy.
    constexpr double kPshbVoutMeanTol = 0.35;
    constexpr double kPshbIoutMeanTol = 0.35;
    for (size_t i = 0; i < wfs.size(); ++i) {
        ConverterPortChecks::check_dc_ports(wfs[i], "PhaseShiftedHalfBridge", i,
                                            Vin, {Vout}, {Iout},
                                            kPshbVoutMeanTol, kPshbIoutMeanTol);
    }
}


// =========================================================================
// Front-end-default reproduction — mirrors the PSFB version. Replays
// PshbWizard.vue defaults verbatim so we can see whether the perceived
// "Simulated" hang on the frontend's PSHB wizard is a missing SPICE
// model, a 2x duplicated-netlist issue, or a real convergence problem.
// =========================================================================
TEST_CASE("Test_Pshb_FrontendDefaults_TimingDiagnostic",
          "[converter-model][pshb-topology][advanced][timing]") {
    using clk = std::chrono::steady_clock;
    auto ms = [](clk::time_point a, clk::time_point b) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
    };

    // ---- Replicate PshbWizard.vue defaults exactly ------------------
    //   inputVoltage:           { nominal: 400, tolerance: 0.1 }
    //   outputVoltage:          12 V
    //   outputPower:            200 W   (→ outputCurrent = 200 / 12 ≈ 16.67 A)
    //   switchingFreq:          100000 Hz
    //   phaseShift:             72 deg
    //   maxPhaseShift:          144 deg
    //   efficiency:             0.95
    //   seriesInductance:       0   (left out → useLeakageInductance branch)
    //   useLeakageInductance:   true
    //   rectifierType:          "fullBridge"
    //   magnetizingInductance:  1e-3 H
    //   ambientTemperature:     25 C
    //   designMode:             "Help me with the design"  → desiredTurnsRatios = [Vin/Vout] = [400/12]
    //   numberOfPeriods:        2
    //   numberOfSteadyStatePeriods: 5
    json advJson;
    {
        json inputVoltage;
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 360.0;
        inputVoltage["maximum"] = 440.0;
        advJson["inputVoltage"] = inputVoltage;
    }
    advJson["rectifierType"] = "fullBridge";
    advJson["useLeakageInductance"] = true;
    advJson["maximumPhaseShift"] = 144.0;
    advJson["desiredTurnsRatios"] = { 400.0 / 12.0 };
    advJson["desiredMagnetizingInductance"] = 1e-3;
    advJson["operatingPoints"] = json::array();
    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0};
        op["outputCurrents"] = {200.0 / 12.0};
        op["switchingFrequency"] = 100000.0;
        op["phaseShift"] = 72.0;
        advJson["operatingPoints"].push_back(op);
    }

    auto t0 = clk::now();
    OpenMagnetics::AdvancedPshb advPshb(advJson);
    auto inputs = advPshb.process();
    auto t1 = clk::now();
    std::cout << "[PSHB-DEFAULT-TIMING] process()                                       took "
              << ms(t0, t1) << " ms" << std::endl;

    auto& designReq = inputs.get_design_requirements();
    std::vector<double> turnsRatios;
    for (const auto& tr : designReq.get_turns_ratios()) {
        if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
    }
    REQUIRE(!turnsRatios.empty());

    double Lm = 0.0;
    if (designReq.get_magnetizing_inductance().get_minimum()) {
        Lm = designReq.get_magnetizing_inductance().get_minimum().value();
    } else if (designReq.get_magnetizing_inductance().get_nominal()) {
        Lm = designReq.get_magnetizing_inductance().get_nominal().value();
    }
    REQUIRE(Lm > 0.0);

    advPshb.set_num_periods_to_extract(2);
    advPshb.set_num_steady_state_periods(5);

    OpenMagnetics::NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping the simulation timing portion of this diagnostic");
        return;
    }

    auto t2 = clk::now();
    auto topologyWaveforms = advPshb.simulate_and_extract_topology_waveforms(turnsRatios, Lm, /*numberOfPeriods*/ 2);
    auto t3 = clk::now();
    std::cout << "[PSHB-DEFAULT-TIMING] simulate_and_extract_topology_waveforms()      took "
              << ms(t2, t3) << " ms"
              << "  (" << topologyWaveforms.size() << " waveform sets)" << std::endl;

    auto t4 = clk::now();
    auto operatingPoints = advPshb.simulate_and_extract_operating_points(turnsRatios, Lm);
    auto t5 = clk::now();
    std::cout << "[PSHB-DEFAULT-TIMING] simulate_and_extract_operating_points()        took "
              << ms(t4, t5) << " ms"
              << "  (" << operatingPoints.size() << " operating points)" << std::endl;

    std::cout << "[PSHB-DEFAULT-TIMING] total simulation wall time                     "
              << ms(t2, t5) << " ms" << std::endl;

    REQUIRE(!topologyWaveforms.empty());
    REQUIRE(!operatingPoints.empty());
}


// =========================================================================
// Dumps the actual secondary waveform shape produced by the PSHB wizard
// defaults so we can see why the front-end plot "looks like crap".
//
// Investigation strategy:
//   1. Replay wizard defaults verbatim.
//   2. Walk the simulated operating point, isolate the secondary
//      voltage + current waveforms.
//   3. Print samples / stats and check the shape against expected:
//        Secondary voltage: bipolar ±Vin/n square (≈±33.3 V),
//        zero during freewheel (1 - Deff) interval.
//        Secondary current: trapezoidal ramps matching primary/n,
//        approximately Io = 16.67 A peak.
//
//   If V or I show NaN / Inf / mostly-zero / wrong-magnitude /
//   wrong-shape, we know which probe is wrong.
// =========================================================================
TEST_CASE("Test_Pshb_FrontendDefaults_SecondaryWaveformShape",
          "[converter-model][pshb-topology][advanced][secondary-shape]") {
    using namespace OpenMagnetics;

    json advJson;
    {
        json inputVoltage;
        inputVoltage["nominal"] = 400.0;
        inputVoltage["minimum"] = 360.0;
        inputVoltage["maximum"] = 440.0;
        advJson["inputVoltage"] = inputVoltage;
    }
    advJson["rectifierType"] = "fullBridge";
    advJson["useLeakageInductance"] = true;
    advJson["maximumPhaseShift"] = 144.0;
    // Half-bridge primary sees ±Vin/2 (DC-blocking cap mid-point sits at
    // Vin/2). Correct turns ratio for Vout = 12 V is therefore
    // (Vin/2)/Vout = 400/(2·12) = 16.67. PshbWizard.vue currently sends
    // Vin/Vout = 33.33 (the FULL-bridge formula); with n=33.33 the
    // reflected secondary peak is Vin/(2·n) = 6 V, *below* Vout = 12 V,
    // so the secondary diodes never forward-bias and the simulated
    // secondary current is ~0. This test uses the CORRECT half-bridge
    // ratio to confirm the netlist itself is fine — the bug is purely
    // in the wizard's buildParams formula.
    // Replicate PshbWizard.vue.buildParams() in "Help me with the design"
    // mode VERBATIM so this test will catch any future regression of the
    // same shape. The fixed formula (post-bug) is `Vin / (4·Vout)` —
    // half-bridge with ~50 % duty headroom. The buggy historical formula
    // was `Vin / Vout`, which uses the full-bridge formula on a half-
    // bridge primary and is the bug this test was created to catch.
    const double n_from_wizard = 400.0 / (4.0 * 12.0);
    advJson["desiredTurnsRatios"] = { n_from_wizard };
    advJson["desiredMagnetizingInductance"] = 1e-3;
    advJson["operatingPoints"] = json::array();
    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0};
        op["outputCurrents"] = {200.0 / 12.0};
        op["switchingFrequency"] = 100000.0;
        op["phaseShift"] = 72.0;
        advJson["operatingPoints"].push_back(op);
    }

    AdvancedPshb advPshb(advJson);
    auto inputs = advPshb.process();

    std::vector<double> turnsRatios;
    for (const auto& tr : inputs.get_design_requirements().get_turns_ratios()) {
        if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
    }
    double Lm = inputs.get_design_requirements().get_magnetizing_inductance()
                    .get_nominal().value_or(1e-3);

    advPshb.set_num_periods_to_extract(2);
    advPshb.set_num_steady_state_periods(5);

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping shape check");
        return;
    }

    auto operatingPoints = advPshb.simulate_and_extract_operating_points(turnsRatios, Lm);
    REQUIRE(!operatingPoints.empty());

    auto dumpStats = [](const std::string& label,
                        const std::vector<double>& xs) {
        double mn = std::numeric_limits<double>::infinity();
        double mx = -std::numeric_limits<double>::infinity();
        double sum = 0.0;
        size_t nan = 0, inf = 0, nonzero = 0;
        for (double v : xs) {
            if (std::isnan(v)) { nan++; continue; }
            if (std::isinf(v)) { inf++; continue; }
            mn = std::min(mn, v);
            mx = std::max(mx, v);
            sum += v;
            if (std::abs(v) > 1e-9) nonzero++;
        }
        double mean = xs.empty() ? 0.0 : sum / xs.size();
        std::cout << "  " << std::left << std::setw(34) << label
                  << "  N="     << std::setw(6) << xs.size()
                  << "  min="   << std::setw(14) << mn
                  << "  max="   << std::setw(14) << mx
                  << "  mean="  << std::setw(14) << mean
                  << "  nonzero=" << nonzero
                  << "  nan=" << nan << "  inf=" << inf << std::endl;
    };

    auto dumpFirstAndLast = [](const std::string& label,
                                const std::vector<double>& xs, size_t k = 8) {
        std::cout << "  " << label << " first " << k << ":";
        for (size_t i = 0; i < std::min(k, xs.size()); ++i) std::cout << " " << xs[i];
        std::cout << std::endl;
        std::cout << "  " << label << " last  " << k << ":";
        size_t start = (xs.size() > k) ? xs.size() - k : 0;
        for (size_t i = start; i < xs.size(); ++i) std::cout << " " << xs[i];
        std::cout << std::endl;
    };

    std::cout << std::scientific << std::setprecision(4);
    for (size_t opi = 0; opi < operatingPoints.size(); ++opi) {
        const auto& op = operatingPoints[opi];
        std::cout << "\n========== OP " << opi << " : " << op.get_name().value_or("?") << " ==========" << std::endl;

        const auto& exc = op.get_excitations_per_winding();
        for (size_t w = 0; w < exc.size(); ++w) {
            std::cout << "  --- winding " << w << " ---" << std::endl;

            if (exc[w].get_voltage()) {
                auto V = exc[w].get_voltage().value();
                if (V.get_waveform()) {
                    auto wf = V.get_waveform().value();
                    dumpStats("Voltage data", wf.get_data());
                    if (wf.get_time()) dumpStats("Voltage time", wf.get_time().value());
                    if (w == 1) dumpFirstAndLast("Voltage data", wf.get_data(), 12);
                } else {
                    std::cout << "  Voltage waveform: <missing>" << std::endl;
                }
            } else {
                std::cout << "  Voltage: <none>" << std::endl;
            }

            if (exc[w].get_current()) {
                auto I = exc[w].get_current().value();
                if (I.get_waveform()) {
                    auto wf = I.get_waveform().value();
                    dumpStats("Current data", wf.get_data());
                    if (wf.get_time()) dumpStats("Current time", wf.get_time().value());
                    if (w == 1) dumpFirstAndLast("Current data", wf.get_data(), 12);
                } else {
                    std::cout << "  Current waveform: <missing>" << std::endl;
                }
            } else {
                std::cout << "  Current: <none>" << std::endl;
            }
        }
    }
    std::cout.unsetf(std::ios::scientific);

    // Sanity gates on the secondary waveform at the wizard's default
    // scenario. These are the cheapest assertions that would have caught
    // the original "looks like crap" bug at HEAD time:
    //
    //   (1) Secondary voltage |peak| must EXCEED the rectifier-output
    //       voltage by enough headroom for diode forward-bias. We chose
    //       1.3·Vout as the lower gate (200 V / 8.3 ≈ 24 V vs Vout = 12 V,
    //       so plenty of slack — but 200 V / 33.3 ≈ 6 V (the bug) would
    //       fall below 1.3·12 = 15.6 V and trip this gate immediately).
    //   (2) Secondary current |peak| must be a non-trivial fraction of
    //       Iout. With the wizard's bug, current was 0.37 A on a 16.7 A
    //       nominal — 2 % of design. A 20 % floor (3.3 A here) is well
    //       above the 2 % bug and well below the ~30 % steady-state value
    //       we see post-fix, leaving healthy margin in both directions.
    const double Vout_nominal = 12.0;
    const double Iout_nominal = 200.0 / 12.0;
    for (size_t opi = 0; opi < operatingPoints.size(); ++opi) {
        const auto& exc = operatingPoints[opi].get_excitations_per_winding();
        REQUIRE(exc.size() >= 2);

        if (exc[1].get_voltage() && exc[1].get_voltage()->get_waveform()) {
            auto wf = exc[1].get_voltage()->get_waveform().value();
            double absMax = 0.0;
            for (double v : wf.get_data()) absMax = std::max(absMax, std::abs(v));
            INFO("OP " << opi << " secondary |V|max = " << absMax
                       << " vs gate 1.3*Vout = " << (1.3 * Vout_nominal));
            CHECK(absMax > 1.3 * Vout_nominal);
            CHECK(absMax < 100.0);
        }

        if (exc[1].get_current() && exc[1].get_current()->get_waveform()) {
            auto wf = exc[1].get_current()->get_waveform().value();
            double absMax = 0.0;
            for (double v : wf.get_data()) absMax = std::max(absMax, std::abs(v));
            INFO("OP " << opi << " secondary |I|max = " << absMax
                       << " vs gate 0.2*Iout = " << (0.2 * Iout_nominal));
            CHECK(absMax > 0.2 * Iout_nominal);
        }
    }
}


} // anonymous namespace
