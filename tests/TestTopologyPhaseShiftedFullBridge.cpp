#include <source_location>
#include "support/Painter.h"
#include "converter_models/PhaseShiftedFullBridge.h"
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

using namespace OpenMagnetics;

namespace {

// ── PtP waveform comparison helpers (mirror of TestTopologyDab helpers) ──
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
// Helper: create a typical PSFB JSON (400V -> 12V, 50A, 100 kHz, center-tapped)
// =========================================================================
json make_psfb_json(double Vin_nom = 400.0, double Vin_min = 370.0, double Vin_max = 410.0,
                    double Vo = 12.0, double Io = 50.0, double Fs = 100000.0,
                    double phaseShift = 126.0,
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
// TEST 1: PSFB Design - Center Tapped, 400V -> 12V, 600W
// =========================================================================
TEST_CASE("Test_Psfb_CenterTapped_Design", "[converter-model][psfb-topology][smoke-test]") {
    auto psfbJson = make_psfb_json();
    OpenMagnetics::Psfb psfb(psfbJson);

    SECTION("Input validation") {
        CHECK(psfb.run_checks(false) == true);
    }

    SECTION("Turns ratio") {
        auto req = psfb.process_design_requirements();
        double n = resolve_dimensional_values(req.get_turns_ratios()[0]);
        // n = Vin * Deff / (Vo + Vd)
        // Deff = 126/180 = 0.7
        // n = 400 * 0.7 / (12 + 0.6) = 280 / 12.6 ≈ 22.22
        double expectedN = 400.0 * (126.0 / 180.0) / (12.0 + 0.6);
        REQUIRE_THAT(n, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
        CHECK(n > 0);
    }

    SECTION("Magnetizing inductance is positive and reasonable") {
        auto req = psfb.process_design_requirements();
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
        CHECK(Lm > 0);
        CHECK(Lm > 10e-6);
        CHECK(Lm < 500e-3);
    }

    SECTION("Computed series inductance is positive") {
        psfb.process_design_requirements();
        double Lr = psfb.get_computed_series_inductance();
        CHECK(Lr > 0);
        CHECK(Lr < 1e-3);
    }

    SECTION("Computed output inductance is positive") {
        psfb.process_design_requirements();
        double Lo = psfb.get_computed_output_inductance();
        CHECK(Lo > 0);
    }

    SECTION("Effective duty cycle (post duty-cycle-loss)") {
        psfb.process_design_requirements();
        double Deff = psfb.get_computed_effective_duty_cycle();
        CHECK(Deff > 0);
        CHECK(Deff < 1.0);
        // D_cmd = 126°/180° = 0.70. With the default Lr seed and the
        // duty-cycle-loss correction (Sabate 1990) the secondary sees
        // D_eff ≈ D_cmd − ΔD with ΔD typically 3–8 % at full load.
        REQUIRE_THAT(Deff, Catch::Matchers::WithinAbs(0.66, 0.05));
    }
}


// =========================================================================
// TEST 2: PSFB Operating Points Generation
// =========================================================================
TEST_CASE("Test_Psfb_OperatingPoints_Generation", "[converter-model][psfb-topology][smoke-test]") {
    auto psfbJson = make_psfb_json();
    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    SECTION("Multiple input voltages") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        // Vin has nominal, min, max → 3 operating points
        CHECK(ops.size() == 3);
    }

    SECTION("Waveform structure") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        REQUIRE(!ops.empty());

        auto& op = ops[0];
        // Primary + 2 CT secondary half-windings = 3 windings
        // (rectifier-type-aware: CT splits into two half-windings).
        CHECK(op.get_excitations_per_winding().size() == 3);

        auto& priExc = op.get_excitations_per_winding()[0];
        REQUIRE(priExc.get_current().has_value());
        REQUIRE(priExc.get_voltage().has_value());
        CHECK(priExc.get_frequency() == 100e3);

        // Waveforms should have 2*256+1 = 513 samples
        auto currentWfm = priExc.get_current()->get_waveform().value();
        CHECK(currentWfm.get_data().size() == 513);

        auto voltageWfm = priExc.get_voltage()->get_waveform().value();
        CHECK(voltageWfm.get_data().size() == 513);
    }

    SECTION("Primary voltage is 3-level") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        auto& priExc = ops[0].get_excitations_per_winding()[0];
        auto vData = priExc.get_voltage()->get_waveform().value().get_data();

        double vMax = *std::max_element(vData.begin(), vData.end());
        double vMin = *std::min_element(vData.begin(), vData.end());

        // Primary voltage should swing to approximately +Vin and -Vin
        double Vin_min = 370.0;
        REQUIRE_THAT(vMax, Catch::Matchers::WithinAbs(Vin_min, Vin_min * 0.05));
        REQUIRE_THAT(vMin, Catch::Matchers::WithinAbs(-Vin_min, Vin_min * 0.05));

        // Should have zero-voltage intervals (freewheeling)
        int zeroCount = 0;
        for (double v : vData) {
            if (std::abs(v) < 1.0) zeroCount++;
        }
        CHECK(zeroCount > 0);
    }

    SECTION("Primary current antisymmetry") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        auto& priExc = ops[0].get_excitations_per_winding()[0];
        auto iData = priExc.get_current()->get_waveform().value().get_data();

        // PSFB primary current should be antisymmetric: i(t + T/2) = -i(t)
        size_t N = iData.size();
        size_t half = N / 2;
        double asymmetrySum = 0;
        for (size_t k = 1; k < half; k++) {
            asymmetrySum += std::abs(iData[k] + iData[half + k]);
        }
        double avgAsymmetry = asymmetrySum / half;
        double iPeak = *std::max_element(iData.begin(), iData.end());
        // Antisymmetry error should be small relative to peak
        CHECK(avgAsymmetry / iPeak < 0.05);
    }

    SECTION("Secondary winding excitation exists") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        auto& secExc = ops[0].get_excitations_per_winding()[1];
        REQUIRE(secExc.get_current().has_value());
        REQUIRE(secExc.get_voltage().has_value());
    }
}


// =========================================================================
// TEST 3: PSFB Waveform Plotting
// =========================================================================
TEST_CASE("Test_Psfb_Waveform_Plotting", "[converter-model][psfb-topology][visualization]") {
    auto psfbJson = make_psfb_json();
    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto ops = psfb.process_operating_points(turnsRatios, Lm);
    REQUIRE(!ops.empty());

    SECTION("Primary current waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_Primary_Current_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Primary voltage waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_Primary_Voltage_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Secondary current waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_Secondary_Current_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
        painter.export_svg();
    }
}


// =========================================================================
// TEST 4: PSFB SPICE Netlist Generation
// =========================================================================
TEST_CASE("Test_Psfb_Spice_Netlist", "[converter-model][psfb-topology][spice]") {
    auto psfbJson = make_psfb_json();
    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    SECTION("Netlist is non-empty and contains key elements") {
        std::string netlist = psfb.generate_ngspice_circuit(turnsRatios, Lm);
        CHECK(!netlist.empty());
        CHECK(netlist.find("Phase-Shifted Full Bridge") != std::string::npos);
        CHECK(netlist.find("L_pri") != std::string::npos);
        CHECK(netlist.find("L_sec") != std::string::npos);
        // Multi-secondary K-matrix uses K1, K2, ... (DAB-style), not K_trafo.
        CHECK(netlist.find("K1 L_pri L_sec_o1") != std::string::npos);
        CHECK(netlist.find("L_out") != std::string::npos);
        CHECK(netlist.find("R_load") != std::string::npos);
        CHECK(netlist.find(".tran") != std::string::npos);
        // Full bridge: 4 switches (SA, SB, SC, SD)
        CHECK(netlist.find("SA ") != std::string::npos);
        CHECK(netlist.find("SB ") != std::string::npos);
        CHECK(netlist.find("SC ") != std::string::npos);
        CHECK(netlist.find("SD ") != std::string::npos);
    }

    SECTION("Netlist saved to file") {
        std::string netlist = psfb.generate_ngspice_circuit(turnsRatios, Lm);
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_Netlist.cir");
        std::ofstream ofs(outFile);
        ofs << netlist;
        ofs.close();
        CHECK(std::filesystem::exists(outFile));
    }
}


// =========================================================================
// TEST 5: PSFB Multiple Outputs
// =========================================================================
TEST_CASE("Test_Psfb_Multiple_Outputs", "[converter-model][psfb-topology][smoke-test]") {
    json psfbJson;
    json inputVoltage;
    inputVoltage["nominal"] = 400.0;
    psfbJson["inputVoltage"] = inputVoltage;
    psfbJson["rectifierType"] = "centerTapped";
    psfbJson["operatingPoints"] = json::array();

    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0, 5.0};
        op["outputCurrents"] = {50.0, 10.0};
        op["switchingFrequency"] = 100000;
        op["phaseShift"] = 126.0;
        psfbJson["operatingPoints"].push_back(op);
    }

    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    CHECK(req.get_turns_ratios().size() == 2);

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto ops = psfb.process_operating_points(turnsRatios, Lm);
    REQUIRE(!ops.empty());

    // Default CT rectifier: Primary + 2 outputs * 2 CT half-windings = 5
    CHECK(ops[0].get_excitations_per_winding().size() == 5);
}


// =========================================================================
// TEST 6: PSFB Static calculations
// =========================================================================
TEST_CASE("Test_Psfb_Static_Calculations", "[converter-model][psfb-topology][unit]") {

    SECTION("Effective duty cycle") {
        CHECK(Psfb::compute_effective_duty_cycle(0.0) == 0.0);
        CHECK(Psfb::compute_effective_duty_cycle(90.0) == 0.5);
        CHECK(Psfb::compute_effective_duty_cycle(180.0) == 1.0);
        REQUIRE_THAT(Psfb::compute_effective_duty_cycle(126.0),
                     Catch::Matchers::WithinAbs(0.7, 1e-6));
    }

    SECTION("Output voltage - center tapped") {
        // Vo = Vin * Deff / n - Vd
        double Vo = Psfb::compute_output_voltage(400.0, 0.7, 22.0, 0.6,
                                                  PsfbRectifierType::CENTER_TAPPED);
        double expected = 400.0 * 0.7 / 22.0 - 0.6;
        REQUIRE_THAT(Vo, Catch::Matchers::WithinAbs(expected, 0.01));
    }

    SECTION("Turns ratio round-trip") {
        double n = Psfb::compute_turns_ratio(400.0, 12.0, 0.7, 0.6,
                                              PsfbRectifierType::CENTER_TAPPED);
        double Vo_check = Psfb::compute_output_voltage(400.0, 0.7, n, 0.6,
                                                        PsfbRectifierType::CENTER_TAPPED);
        REQUIRE_THAT(Vo_check, Catch::Matchers::WithinAbs(12.0, 0.01));
    }
}


// =========================================================================
// TEST 7: PSFB Current Doubler rectifier
// =========================================================================
TEST_CASE("Test_Psfb_CurrentDoubler_Design", "[converter-model][psfb-topology][smoke-test]") {
    auto psfbJson = make_psfb_json(400.0, 370.0, 410.0, 12.0, 50.0, 100000.0,
                                   126.0, "currentDoubler");
    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    double n = resolve_dimensional_values(req.get_turns_ratios()[0]);
    // Current doubler: n = Vin * Deff / (2*(Vo+Vd))
    double expectedN = 400.0 * 0.7 / (2.0 * (12.0 + 0.6));
    REQUIRE_THAT(n, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
}


// =========================================================================
// TEST 8: AdvancedPsfb JSON round-trip
// =========================================================================
TEST_CASE("Test_AdvancedPsfb_Process", "[converter-model][psfb-topology][advanced]") {
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
        op["outputCurrents"] = {50.0};
        op["switchingFrequency"] = 100000;
        op["phaseShift"] = 126.0;
        advJson["operatingPoints"].push_back(op);
    }

    advJson["desiredTurnsRatios"] = {22.0};
    advJson["desiredMagnetizingInductance"] = 1e-3;

    OpenMagnetics::AdvancedPsfb advPsfb(advJson);
    auto inputs = advPsfb.process();

    CHECK(inputs.get_design_requirements().get_turns_ratios().size() == 1);
    double n = resolve_dimensional_values(inputs.get_design_requirements().get_turns_ratios()[0]);
    REQUIRE_THAT(n, Catch::Matchers::WithinAbs(22.0, 0.01));

    double Lm = resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance());
    REQUIRE_THAT(Lm, Catch::Matchers::WithinAbs(1e-3, 1e-6));

    CHECK(!inputs.get_operating_points().empty());
}


// =========================================================================
// TEST: Duty-cycle loss diagnostic populated and matches Sabate's formula
// =========================================================================
TEST_CASE("Test_Psfb_DutyCycleLoss", "[converter-model][psfb-topology][unit]") {
    // Single input-voltage so the lastDutyCycleLoss diagnostic unambiguously
    // refers to the same Vin as the test's expected formula.
    json psfbJson;
    json inputVoltage; inputVoltage["nominal"] = 400.0;
    psfbJson["inputVoltage"] = inputVoltage;
    psfbJson["rectifierType"] = "fullBridge";
    psfbJson["seriesInductance"] = 5e-6;
    psfbJson["operatingPoints"] = json::array();
    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0};
        op["outputCurrents"] = {50.0};
        op["switchingFrequency"] = 100000.0;
        op["phaseShift"] = 126.0;
        psfbJson["operatingPoints"].push_back(op);
    }

    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios())
        turnsRatios.push_back(resolve_dimensional_values(tr));
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto ops = psfb.process_operating_points(turnsRatios, Lm);
    REQUIRE(!ops.empty());

    double dcl = psfb.get_last_duty_cycle_loss();
    double Deff_actual = psfb.get_last_effective_duty_cycle();
    double n = turnsRatios[0];

    // Sabate ΔD = 4·Lk·Io·Fs / (n·Vin)
    double dcl_expected = 4.0 * 5e-6 * 50.0 * 100e3 / (n * 400.0);
    REQUIRE_THAT(dcl, Catch::Matchers::WithinRel(dcl_expected, 0.05));
    CHECK(Deff_actual < 0.7);
    CHECK(Deff_actual > 0.0);
}


// =========================================================================
// TEST: ZVS-margin diagnostic — sign flips when Lk crosses the threshold
// =========================================================================
TEST_CASE("Test_Psfb_ZvsBoundary", "[converter-model][psfb-topology][unit]") {
    // ZVS at light load is hard; we sweep Lr at fixed load and look for the
    // sign flip. Vin=400V, Coss=200pF, Ip_pk ≈ Io/n + Im_peak.
    auto psfbJson = make_psfb_json(400.0, 370.0, 410.0, 12.0, 50.0, 100000.0,
                                   126.0, "fullBridge");
    psfbJson["seriesInductance"] = 1e-6;  // small Lk → likely no ZVS
    OpenMagnetics::Psfb psfb_low(psfbJson);
    auto req_low = psfb_low.process_design_requirements();
    std::vector<double> tr_low;
    for (auto& t : req_low.get_turns_ratios()) tr_low.push_back(resolve_dimensional_values(t));
    psfb_low.process_operating_points(tr_low, resolve_dimensional_values(req_low.get_magnetizing_inductance()));
    double margin_low = psfb_low.get_last_zvs_margin_lagging();

    psfbJson["seriesInductance"] = 50e-6;  // large Lk → ZVS easily
    OpenMagnetics::Psfb psfb_high(psfbJson);
    auto req_high = psfb_high.process_design_requirements();
    std::vector<double> tr_high;
    for (auto& t : req_high.get_turns_ratios()) tr_high.push_back(resolve_dimensional_values(t));
    psfb_high.process_operating_points(tr_high, resolve_dimensional_values(req_high.get_magnetizing_inductance()));
    double margin_high = psfb_high.get_last_zvs_margin_lagging();

    // Margin grows monotonically with Lk (E_inductor = ½·Lk·Ip² scales linearly
    // with Lk; Ip stays roughly the same).
    CHECK(margin_high > margin_low);
    CHECK(psfb_low.get_last_resonant_transition_time() > 0.0);
    CHECK(psfb_high.get_last_resonant_transition_time() > psfb_low.get_last_resonant_transition_time());
}


// =========================================================================
// TEST: Rectifier-type secondary winding count
// =========================================================================
TEST_CASE("Test_Psfb_RectifierType_WindingCount", "[converter-model][psfb-topology][unit]") {
    auto run = [](const std::string& rect) {
        auto j = make_psfb_json(400.0, 370.0, 410.0, 12.0, 50.0, 100000.0,
                                126.0, rect);
        OpenMagnetics::Psfb psfb(j);
        auto req = psfb.process_design_requirements();
        std::vector<double> tr;
        for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
        auto ops = psfb.process_operating_points(tr,
            resolve_dimensional_values(req.get_magnetizing_inductance()));
        return ops[0].get_excitations_per_winding().size();
    };
    // CT  : Pri + Sec0a + Sec0b   = 3
    // CD  : Pri + Sec0             = 2 (single secondary, two output inductors)
    // FB  : Pri + Sec0             = 2
    CHECK(run("centerTapped")    == 3);
    CHECK(run("currentDoubler")  == 2);
    CHECK(run("fullBridge")      == 2);
}


// =========================================================================
// TEST: Current-doubler emits a second output inductor in extra components
// =========================================================================
TEST_CASE("Test_Psfb_CurrentDoubler_TwoOutputInductors",
          "[converter-model][psfb-topology][unit]") {
    auto j = make_psfb_json(400.0, 370.0, 410.0, 12.0, 50.0, 100000.0,
                            126.0, "currentDoubler");
    OpenMagnetics::Psfb psfb(j);
    auto req = psfb.process_design_requirements();
    std::vector<double> tr;
    for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
    psfb.process_operating_points(tr,
        resolve_dimensional_values(req.get_magnetizing_inductance()));

    auto extras = psfb.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    // Lo, Lo2, (optional Lr)
    CHECK(extras.size() >= 2);

    // Confirm the two output inductors are named differently
    int loCount = 0;
    for (auto& v : extras) {
        if (std::holds_alternative<OpenMagnetics::Inputs>(v)) {
            auto name = std::get<OpenMagnetics::Inputs>(v).get_design_requirements().get_name().value_or("");
            if (name == "outputInductor")  loCount++;
            if (name == "outputInductor2") loCount++;
        }
    }
    CHECK(loCount == 2);
}


// =========================================================================
// TEST: Analytical vs NgSpice point-to-point (FB rectifier, single output)
// =========================================================================
TEST_CASE("Test_Psfb_PtP_AnalyticalVsNgspice", "[converter-model][psfb-topology][ngspice-simulation][ptpcomparison]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    // Single input voltage so analytical and simulate_and_extract orderings
    // match (the latter uses collect_input_voltages which is nominal-first;
    // the former sorts ascending — they only line up when there's one Vin).
    json j;
    json inputVoltage; inputVoltage["nominal"] = 400.0;
    j["inputVoltage"] = inputVoltage;
    j["rectifierType"] = "fullBridge";
    j["seriesInductance"] = 5e-6;
    j["operatingPoints"] = json::array();
    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0};
        op["outputCurrents"] = {50.0};
        op["switchingFrequency"] = 100000.0;
        op["phaseShift"] = 126.0;
        j["operatingPoints"].push_back(op);
    }

    OpenMagnetics::Psfb psfb(j);
    psfb.set_num_periods_to_extract(1);

    auto req = psfb.process_design_requirements();
    std::vector<double> tr;
    for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto analytical = psfb.process_operating_points(tr, Lm);
    REQUIRE(!analytical.empty());
    auto aTime = ptp_current_time(analytical[0], 0);
    auto aData = ptp_current(analytical[0], 0);
    REQUIRE(!aData.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // Dump the SPICE netlist for offline inspection.
    {
        auto netFile = outputFilePath;
        netFile.append("Test_Psfb_PtP_netlist.cir");
        std::filesystem::remove(netFile);
        std::ofstream nf(netFile);
        nf << psfb.generate_ngspice_circuit(tr, Lm, 0, 0);
    }

    auto sim = psfb.simulate_and_extract_operating_points(tr, Lm);
    REQUIRE(!sim.empty());
    auto sTime = ptp_current_time(sim[0], 0);
    auto sData = ptp_current(sim[0], 0);
    if (!sData.empty()) {
        auto sResampled = ptp_interp(sTime, sData, 256);

        // Dump CSV for offline inspection.
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_PtP_compare.csv");
        std::filesystem::remove(outFile);
        std::ofstream f(outFile);
        f << "k,t_a,a_raw,t_s,s_raw\n";
        size_t nA = std::min(aData.size(), aTime.size());
        size_t nS = std::min(sData.size(), sTime.size());
        size_t nMax = std::max(nA, nS);
        for (size_t k = 0; k < nMax; ++k) {
            f << k << ",";
            if (k < nA) f << aTime[k] << "," << aData[k];
            else f << ",";
            f << ",";
            if (k < nS) f << sTime[k] << "," << sData[k];
            else f << ",";
            f << "\n";
        }

        double nrmse = ptp_nrmse(aResampled, sResampled);
        INFO("PSFB primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        // DAB-quality threshold (≤ 0.15 + small margin for PSFB's harder
        // freewheel modelling). Achieved: 0.15 with the corrected
        // Im_peak (D_cmd, not Deff), correct ILo trough/peak geometry,
        // and the empirical 30 % freewheel-decay model.
        CHECK(nrmse < 0.18);
    }
}


// =========================================================================
// TEST: Generated SPICE netlist contains required robustness features
// =========================================================================
TEST_CASE("Test_Psfb_SpiceNetlistRobustness", "[converter-model][psfb-topology][spice]") {
    auto j = make_psfb_json();
    OpenMagnetics::Psfb psfb(j);
    auto req = psfb.process_design_requirements();
    std::vector<double> tr;
    for (auto& t : req.get_turns_ratios()) tr.push_back(resolve_dimensional_values(t));
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    auto net = psfb.generate_ngspice_circuit(tr, Lm, 0, 0);

    // Snubbers
    CHECK(net.find("Rsnub_QA") != std::string::npos);
    CHECK(net.find("Csnub_QA") != std::string::npos);
    // Solver options
    CHECK(net.find("METHOD=GEAR") != std::string::npos);
    CHECK(net.find("ITL1=500") != std::string::npos);
    // Bridge differential probe
    CHECK(net.find("Evab vab 0 mid_A mid_C 1") != std::string::npos);
    // IC pre-charge
    CHECK(net.find(".ic v(out_node_o1)") != std::string::npos);
    // Rectifier-type-aware (default centerTapped)
    CHECK(net.find("Center-tapped rectifier") != std::string::npos);
}


} // anonymous namespace
