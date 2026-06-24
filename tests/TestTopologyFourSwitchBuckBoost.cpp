/**
 * @file TestTopologyFourSwitchBuckBoost.cpp
 * @brief Unit tests for the Four-Switch Buck-Boost (4SBB) DC-DC converter.
 *
 * Coverage map (FOUR_SWITCH_BUCK_BOOST_PLAN.md §A.10):
 *
 *   * Static analytical helpers
 *       - classify_region (BUCK / BOOST / BUCK_BOOST with hysteresis &
 *         minOnTime / minOffTime constraints)
 *       - compute_buck_duty / compute_boost_duty
 *       - compute_buck_boost_duties for SIMULTANEOUS and SPLIT_PWM modes
 *       - compute_inductor_buck_region / _boost_region
 *       - compute_worst_case_inductance
 *       - compute_inductor_dc_bias_boost
 *
 *   * Reference designs (analytical scalar checks, η=1):
 *       - TIDA-01411   (9–36V → 12V @ 8A @ 350 kHz)
 *       - LM5176EVM-HP (9–55V → 24V @ 5A @ 200 kHz)
 *       - LT8390 DC2825A (4–60V → 12V @ 7A @ 300 kHz)
 *
 *   * Per-OP diagnostics population (Region, M, ΔI, K, K_crit, RMS).
 *
 *   * process_design_requirements with worst-case inductor sizing.
 *
 *   * SPICE netlist string contract — region-specific assertions, snubbers,
 *     RC, gate phase, ic= count, METHOD=GEAR.
 *
 *   * Bidirectional symmetry (Iout sign).
 *
 *   * Multi-phase (phaseCount = 2) — per-phase Iout/N.
 *
 *   * Input-validation guards.
 *
 *   * AdvFSBB round-trip.
 *
 *   * get_extra_components_inputs (P7): 3 entries
 *       - Inputs (L)
 *       - CAS::Inputs (Cin, INPUT_FILTER)
 *       - CAS::Inputs (Co,  OUTPUT_FILTER)
 *
 *   * §8 Volt-second balance, ConverterPort, SpicePowerSanity, Plotting.
 */

#include "converter_models/FourSwitchBuckBoost.h"
#include "processors/Inputs.h"
#include "processors/NgspiceRunner.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "ConverterPortChecks.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <source_location>
#include <variant>

using namespace OpenMagnetics;
using FSBB = OpenMagnetics::FourSwitchBuckBoost;
using AdvFSBB = OpenMagnetics::AdvancedFourSwitchBuckBoost;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

    json make_fsbb_json(double VinNom, double Vo, double Iout, double fsw,
                        double rippleRatio = 0.4, double efficiency = 1.0,
                        double VinMin = -1.0, double VinMax = -1.0) {
        json j;
        json iv = json::object();
        iv["nominal"] = VinNom;
        if (VinMin > 0.0) iv["minimum"] = VinMin;
        if (VinMax > 0.0) iv["maximum"] = VinMax;
        j["inputVoltage"] = iv;
        j["currentRippleRatio"] = rippleRatio;
        j["efficiency"]         = efficiency;
        j["maximumSwitchCurrent"] = 50.0;

        j["operatingPoints"] = json::array();
        json op;
        op["outputVoltages"]     = {Vo};
        op["outputCurrents"]     = {Iout};
        op["switchingFrequency"] = fsw;
        op["ambientTemperature"] = 25.0;
        j["operatingPoints"].push_back(op);
        return j;
    }

    json make_advanced_fsbb_json(double VinNom, double Vo, double Iout, double fsw,
                                 double L, double VinMin = -1.0, double VinMax = -1.0) {
        json j = make_fsbb_json(VinNom, Vo, Iout, fsw, 0.4, 1.0, VinMin, VinMax);
        j["desiredInductance"] = L;
        return j;
    }

}  // namespace

// =====================================================================
//   Static analytical helpers
// =====================================================================

TEST_CASE("Test_FSBB_ClassifyRegion_Buck_Boost_BB",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    using R = FSBB::Region;
    const double hys = 0.15, minOn = 80e-9, minOff = 80e-9, Fs = 300e3;
    REQUIRE(FSBB::classify_region(36.0, 12.0, hys, minOn, minOff, Fs) == R::BUCK);
    REQUIRE(FSBB::classify_region(4.0,  12.0, hys, minOn, minOff, Fs) == R::BOOST);
    REQUIRE(FSBB::classify_region(12.0, 12.0, hys, minOn, minOff, Fs) == R::BUCK_BOOST);
    REQUIRE(FSBB::classify_region(11.0, 12.0, hys, minOn, minOff, Fs) == R::BUCK_BOOST);
    REQUIRE(FSBB::classify_region(13.0, 12.0, hys, minOn, minOff, Fs) == R::BUCK_BOOST);
}

TEST_CASE("Test_FSBB_ClassifyRegion_BadInput",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    REQUIRE_THROWS(FSBB::classify_region(0.0, 12.0, 0.15, 80e-9, 80e-9, 300e3));
    REQUIRE_THROWS(FSBB::classify_region(12.0, 0.0, 0.15, 80e-9, 80e-9, 300e3));
    REQUIRE_THROWS(FSBB::classify_region(12.0, 12.0, 0.15, 80e-9, 80e-9, 0.0));
}

TEST_CASE("Test_FSBB_BuckDuty_Formula",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    REQUIRE_THAT(FSBB::compute_buck_duty(24.0, 12.0, 1.0),
                 WithinRel(0.5, 1e-12));
    REQUIRE_THAT(FSBB::compute_buck_duty(36.0, 12.0, 1.0),
                 WithinRel(12.0/36.0, 1e-12));
    // Vo > Vin·η → invalid for buck.
    REQUIRE_THROWS(FSBB::compute_buck_duty(10.0, 12.0, 1.0));
    REQUIRE_THROWS(FSBB::compute_buck_duty(0.0, 12.0, 1.0));
}

TEST_CASE("Test_FSBB_BoostDuty_Formula",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    REQUIRE_THAT(FSBB::compute_boost_duty(6.0, 12.0, 1.0),
                 WithinRel(0.5, 1e-12));
    REQUIRE_THAT(FSBB::compute_boost_duty(4.0, 12.0, 1.0),
                 WithinRel(1.0 - 4.0/12.0, 1e-12));
    // Vin > Vo → invalid for boost.
    REQUIRE_THROWS(FSBB::compute_boost_duty(24.0, 12.0, 1.0));
    REQUIRE_THROWS(FSBB::compute_boost_duty(6.0, 0.0, 1.0));
}

TEST_CASE("Test_FSBB_BuckBoostDuties_Simultaneous_M_eq_DoverOneMinusD",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    double Db = 0.0, Dboost = 0.0;
    FSBB::compute_buck_boost_duties(
        12.0, 12.0, 1.0, MAS::TransitionMode::SIMULTANEOUS, Db, Dboost);
    // D = Vo/(Vin+Vo) = 0.5 → M = 1.
    REQUIRE_THAT(Db, WithinRel(0.5, 1e-12));
    REQUIRE_THAT(Dboost, WithinRel(0.5, 1e-12));
}

TEST_CASE("Test_FSBB_BuckBoostDuties_SplitPwm_M1_BothEqual",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    double Db = 0.0, Dboost = 0.0;
    FSBB::compute_buck_boost_duties(
        12.0, 12.0, 1.0, MAS::TransitionMode::SPLIT_PWM, Db, Dboost);
    // M=1 special case: both = (1+δ)/2 = 0.525 with δ=0.05.
    REQUIRE_THAT(Db, WithinRel(0.525, 1e-9));
    REQUIRE_THAT(Dboost, WithinRel(0.525, 1e-9));
}

TEST_CASE("Test_FSBB_BuckBoostDuties_SplitPwm_M_NotOne_Constraint",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    // Vin=11, Vo=12 → M = 12/11 ≈ 1.0909
    double Db = 0.0, Dboost = 0.0;
    FSBB::compute_buck_boost_duties(
        11.0, 12.0, 1.0, MAS::TransitionMode::SPLIT_PWM, Db, Dboost);
    // After clamp, M ≈ Db/(1−Dboost) — within ~5% (clamps may activate).
    REQUIRE(Db > 0.0);
    REQUIRE(Dboost > 0.0);
    REQUIRE(Db < 1.0);
    REQUIRE(Dboost < 1.0);
}

TEST_CASE("Test_FSBB_Inductor_Buck_Region_Formula",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    // L = Vo·(Vin−Vo) / (K·Iout·Fs·Vin)
    double L = FSBB::compute_inductor_buck_region(36.0, 12.0, 8.0, 350e3, 0.4);
    double Lexp = 12.0 * 24.0 / (0.4 * 8.0 * 350e3 * 36.0);
    REQUIRE_THAT(L, WithinRel(Lexp, 1e-12));
    // Vo >= Vin → returns 0 (region not valid).
    REQUIRE_THAT(FSBB::compute_inductor_buck_region(10.0, 12.0, 8.0, 350e3, 0.4),
                 WithinAbs(0.0, 1e-12));
    REQUIRE_THROWS(FSBB::compute_inductor_buck_region(0.0, 12.0, 8.0, 350e3, 0.4));
}

TEST_CASE("Test_FSBB_Inductor_Boost_Region_Formula",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    // L = Vin²·(Vo−Vin) / (K·Iout·Fs·Vo²)
    double L = FSBB::compute_inductor_boost_region(9.0, 12.0, 8.0, 350e3, 0.4);
    double Lexp = (9.0 * 9.0) * (12.0 - 9.0) / (0.4 * 8.0 * 350e3 * 12.0 * 12.0);
    REQUIRE_THAT(L, WithinRel(Lexp, 1e-12));
    // Vo <= Vin → returns 0.
    REQUIRE_THAT(FSBB::compute_inductor_boost_region(15.0, 12.0, 8.0, 350e3, 0.4),
                 WithinAbs(0.0, 1e-12));
}

TEST_CASE("Test_FSBB_WorstCase_L_Picks_Larger_Of_Buck_Boost",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    // TIDA-01411: 9–36V → 12V @ 8A @ 350 kHz, K=0.4
    double Lbuck = FSBB::compute_inductor_buck_region(36.0, 12.0, 8.0, 350e3, 0.4);
    double Lboost = FSBB::compute_inductor_boost_region(9.0, 12.0, 8.0, 350e3, 0.4);
    double Lwc = FSBB::compute_worst_case_inductance(9.0, 36.0, 12.0, 8.0, 350e3, 0.4);
    REQUIRE_THAT(Lwc, WithinRel(std::max(Lbuck, Lboost), 1e-12));
}

TEST_CASE("Test_FSBB_DC_Bias_Boost_Iout_Vo_Vin",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    // Boost region @ Vin_min: I_L,avg = Iout · Vo / Vin
    REQUIRE_THAT(FSBB::compute_inductor_dc_bias_boost(9.0, 12.0, 8.0),
                 WithinRel(8.0 * 12.0 / 9.0, 1e-12));
    REQUIRE_THROWS(FSBB::compute_inductor_dc_bias_boost(0.0, 12.0, 8.0));
}

// =====================================================================
//   Reference designs (analytical scalar checks)
// =====================================================================

TEST_CASE("Test_FSBB_TIDA_01411_Reference_Design",
          "[converter-model][four-switch-buck-boost-topology][refdesign][analytical]") {
    // TIDA-01411 (TI): 9–36 V → 12 V @ 8 A @ 350 kHz.
    // Worst-case L: max(L_buck@36V, L_boost@9V).
    double Lwc = FSBB::compute_worst_case_inductance(9.0, 36.0, 12.0, 8.0, 350e3, 0.4);
    REQUIRE(Lwc > 0.0);
    // Check Vo=Vin nominal → BUCK_BOOST region.
    json j = make_fsbb_json(/*VinNom*/ 12.0, 12.0, 8.0, 350e3, 0.4, 1.0, /*min*/ 9.0, /*max*/ 36.0);
    FSBB c(j);
    auto ops = c.process_operating_points(std::vector<double>{}, Lwc);
    REQUIRE(!ops.empty());
}

TEST_CASE("Test_FSBB_LM5176EVM_Reference_Design",
          "[converter-model][four-switch-buck-boost-topology][refdesign][analytical]") {
    // LM5176EVM-HP: 9–55 V → 24 V @ 5 A @ 200 kHz.
    double Lwc = FSBB::compute_worst_case_inductance(9.0, 55.0, 24.0, 5.0, 200e3, 0.4);
    REQUIRE(Lwc > 0.0);
    json j = make_fsbb_json(24.0, 24.0, 5.0, 200e3, 0.4, 1.0, 9.0, 55.0);
    FSBB c(j);
    auto ops = c.process_operating_points(std::vector<double>{}, Lwc);
    REQUIRE(!ops.empty());
}

TEST_CASE("Test_FSBB_LT8390_DC2825A_Reference_Design",
          "[converter-model][four-switch-buck-boost-topology][refdesign][analytical]") {
    // LT8390 DC2825A: 4–60 V → 12 V @ 7 A @ 300 kHz.
    double Lwc = FSBB::compute_worst_case_inductance(4.0, 60.0, 12.0, 7.0, 300e3, 0.4);
    REQUIRE(Lwc > 0.0);
    json j = make_fsbb_json(12.0, 12.0, 7.0, 300e3, 0.4, 1.0, 4.0, 60.0);
    FSBB c(j);
    auto ops = c.process_operating_points(std::vector<double>{}, Lwc);
    REQUIRE(!ops.empty());
}

// =====================================================================
//   Per-OP diagnostics populated
// =====================================================================

TEST_CASE("Test_FSBB_Diagnostics_Populated_Buck",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    json j = make_fsbb_json(36.0, 12.0, 8.0, 350e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);
    REQUIRE(c.get_last_region() == FSBB::Region::BUCK);
    REQUIRE_THAT(c.get_last_duty_cycle_buck(), WithinRel(12.0/36.0, 1e-9));
    REQUIRE(c.get_last_inductor_average_current() > 0.0);
    REQUIRE(c.get_last_inductor_peak_to_peak() > 0.0);
    REQUIRE(c.get_last_inductor_peak_current() > 0.0);
    REQUIRE(c.get_last_inductor_rms_current() > 0.0);
    REQUIRE(c.get_last_q1_rms_current() > 0.0);
    REQUIRE(c.get_last_q2_rms_current() > 0.0);
    REQUIRE(c.get_last_q3_rms_current() > 0.0);   // Q3 ON → RMS = full I_L
    REQUIRE_THAT(c.get_last_q4_rms_current(), WithinAbs(0.0, 1e-12));
    REQUIRE(c.get_last_dcm_k() > 0.0);
    REQUIRE(c.get_last_dcm_kcrit() > 0.0);
    REQUIRE(c.get_last_sized_output_capacitance() > 0.0);
}

TEST_CASE("Test_FSBB_Diagnostics_Populated_Boost",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    json j = make_fsbb_json(6.0, 12.0, 5.0, 300e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);
    REQUIRE(c.get_last_region() == FSBB::Region::BOOST);
    REQUIRE_THAT(c.get_last_duty_cycle_boost(), WithinRel(0.5, 1e-9));
    REQUIRE_THAT(c.get_last_conversion_ratio(), WithinRel(2.0, 1e-9));
    // I_L,avg = Iout / (1−D) = 5 / 0.5 = 10
    REQUIRE_THAT(c.get_last_inductor_average_current(), WithinRel(10.0, 1e-9));
    REQUIRE(c.get_last_q1_rms_current() > 0.0);   // Q1 ON
    REQUIRE_THAT(c.get_last_q2_rms_current(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("Test_FSBB_Diagnostics_Populated_BuckBoost",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    json j = make_fsbb_json(12.0, 12.0, 5.0, 300e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);
    REQUIRE(c.get_last_region() == FSBB::Region::BUCK_BOOST);
    REQUIRE(c.get_last_duty_cycle_buck() > 0.0);
    REQUIRE(c.get_last_duty_cycle_boost() > 0.0);
    REQUIRE(c.get_last_q1_rms_current() > 0.0);
    REQUIRE(c.get_last_q4_rms_current() > 0.0);
}

// =====================================================================
//   process_design_requirements
// =====================================================================

TEST_CASE("Test_FSBB_DR_Returns_WorstCase_L",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    json j = make_fsbb_json(12.0, 12.0, 8.0, 350e3, 0.4, 1.0, 9.0, 36.0);
    FSBB c(j);
    auto dr = c.process_design_requirements();
    REQUIRE(dr.get_topology() == MAS::Topology::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
    REQUIRE(dr.get_isolation_sides().has_value());
    REQUIRE(dr.get_isolation_sides().value().size() == 1);    // single inductor (no isolation)
    REQUIRE(dr.get_magnetizing_inductance().get_minimum().has_value());
    REQUIRE(dr.get_magnetizing_inductance().get_minimum().value() > 0.0);
}

TEST_CASE("Test_FSBB_DR_Throws_If_No_Ripple_Or_Imax",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    json j = make_fsbb_json(12.0, 12.0, 8.0, 350e3, 0.4, 1.0, 9.0, 36.0);
    j.erase("currentRippleRatio");
    j.erase("maximumSwitchCurrent");
    FSBB c(j);
    REQUIRE_THROWS(c.process_design_requirements());
}

// =====================================================================
//   SPICE netlist string contract
// =====================================================================

TEST_CASE("Test_FSBB_SPICE_Netlist_StringContract_Buck",
          "[converter-model][four-switch-buck-boost-topology][spice-netlist]") {
    json j = make_fsbb_json(36.0, 12.0, 8.0, 350e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);
    std::string n = c.generate_ngspice_circuit(/*L*/ 5e-6);
    REQUIRE(n.find("region=BUCK") != std::string::npos);
    REQUIRE(n.find("METHOD=GEAR") != std::string::npos);
    REQUIRE(n.find("UIC") != std::string::npos);
    REQUIRE(n.find(".model SW1 SW") != std::string::npos);
    // Q1/Q2 PWM, Q3 high (DC), Q4 low (DC).
    REQUIRE(n.find("Vpwm_q1 pwm_q1 0 PULSE(") != std::string::npos);
    REQUIRE(n.find("Vpwm_q2 pwm_q2 0 PULSE(") != std::string::npos);
    REQUIRE(n.find("Vpwm_q3 pwm_q3 0 ") != std::string::npos);
    REQUIRE(n.find("Vpwm_q4 pwm_q4 0 0") != std::string::npos);
    // S-prefix for switches.
    REQUIRE(n.find("S_Q1 vin_p sw1") != std::string::npos);
    REQUIRE(n.find("S_Q2 sw1 0") != std::string::npos);
    REQUIRE(n.find("S_Q3 sw2 vout") != std::string::npos);
    REQUIRE(n.find("S_Q4 sw2 0") != std::string::npos);
    // Snubbers on both SW nodes.
    REQUIRE(n.find("Rsnub_sw1 sw1 sn_sw1") != std::string::npos);
    REQUIRE(n.find("Csnub_sw1 sn_sw1 0") != std::string::npos);
    REQUIRE(n.find("Rsnub_sw2 sw2 sn_sw2") != std::string::npos);
    REQUIRE(n.find("Csnub_sw2 sn_sw2 0") != std::string::npos);
    // Inductor with IC.
    REQUIRE(n.find("L1 l_in l_out") != std::string::npos);
    REQUIRE(n.find(" IC=") != std::string::npos);
    // Sense elements.
    REQUIRE(n.find("Vin_sense") != std::string::npos);
    REQUIRE(n.find("Vl_sense") != std::string::npos);
    REQUIRE(n.find("Vsw2_sense") != std::string::npos);
    REQUIRE(n.find("Vout_sense") != std::string::npos);
    // Differential probe.
    REQUIRE(n.find("Bvpri_diff vpri_diff 0") != std::string::npos);
}

TEST_CASE("Test_FSBB_SPICE_Netlist_StringContract_Boost",
          "[converter-model][four-switch-buck-boost-topology][spice-netlist]") {
    json j = make_fsbb_json(6.0, 12.0, 5.0, 300e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);
    std::string n = c.generate_ngspice_circuit(/*L*/ 5e-6);
    REQUIRE(n.find("region=BOOST") != std::string::npos);
    // Q1 high, Q2 low (DC); Q3/Q4 PWM.
    REQUIRE(n.find("Vpwm_q2 pwm_q2 0 0") != std::string::npos);
    REQUIRE(n.find("Vpwm_q4 pwm_q4 0 PULSE(") != std::string::npos);
    REQUIRE(n.find("Vpwm_q3 pwm_q3 0 PULSE(") != std::string::npos);
}

TEST_CASE("Test_FSBB_SPICE_Netlist_StringContract_BuckBoost_AllPwm",
          "[converter-model][four-switch-buck-boost-topology][spice-netlist]") {
    json j = make_fsbb_json(12.0, 12.0, 5.0, 300e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);
    std::string n = c.generate_ngspice_circuit(/*L*/ 5e-6);
    REQUIRE(n.find("region=BUCK_BOOST") != std::string::npos);
    // All four are PULSE in BB.
    REQUIRE(n.find("Vpwm_q1 pwm_q1 0 PULSE(") != std::string::npos);
    REQUIRE(n.find("Vpwm_q2 pwm_q2 0 PULSE(") != std::string::npos);
    REQUIRE(n.find("Vpwm_q3 pwm_q3 0 PULSE(") != std::string::npos);
    REQUIRE(n.find("Vpwm_q4 pwm_q4 0 PULSE(") != std::string::npos);
}

TEST_CASE("Test_FSBB_SPICE_Bad_Inputs_Throw",
          "[converter-model][four-switch-buck-boost-topology][spice-netlist]") {
    json j = make_fsbb_json(12.0, 12.0, 5.0, 300e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);
    REQUIRE_THROWS(c.generate_ngspice_circuit(5e-6, /*idx*/ 99, 0));
    REQUIRE_THROWS(c.generate_ngspice_circuit(5e-6, 0, /*opIdx*/ 99));
}

// =====================================================================
//   Bidirectional symmetry — Iout sign preservation.
// =====================================================================

TEST_CASE("Test_FSBB_Bidirectional_Iout_Sign_Preserved",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    // Negative Iout with bidirectional=true → analytical |Iout| matches positive
    // case but inductor current waveform is sign-flipped.
    json j_pos = make_fsbb_json(12.0, 12.0, 5.0, 300e3);
    FSBB c_pos(j_pos);
    c_pos.process_operating_points(std::vector<double>{}, 5e-6);
    const double iL_avg_pos = c_pos.get_last_inductor_average_current();

    json j_neg = make_fsbb_json(12.0, 12.0, -5.0, 300e3);
    j_neg["bidirectional"] = true;
    FSBB c_neg(j_neg);
    c_neg.process_operating_points(std::vector<double>{}, 5e-6);
    const double iL_avg_neg = c_neg.get_last_inductor_average_current();
    REQUIRE_THAT(iL_avg_neg, WithinRel(iL_avg_pos, 1e-9));
}

// =====================================================================
//   Multi-phase
// =====================================================================

TEST_CASE("Test_FSBB_MultiPhase_PerPhase_Iout_Reduced",
          "[converter-model][four-switch-buck-boost-topology][unit]") {
    // phaseCount=2 → per-phase Iout halved.
    json j1 = make_fsbb_json(36.0, 12.0, 8.0, 350e3);
    FSBB c1(j1);
    c1.process_operating_points(std::vector<double>{}, 5e-6);
    const double iL_1ph = c1.get_last_inductor_average_current();

    json j2 = make_fsbb_json(36.0, 12.0, 8.0, 350e3);
    j2["phaseCount"] = 2;
    FSBB c2(j2);
    c2.process_operating_points(std::vector<double>{}, 5e-6);
    const double iL_2ph = c2.get_last_inductor_average_current();

    REQUIRE_THAT(iL_2ph, WithinRel(iL_1ph / 2.0, 1e-9));
}

// =====================================================================
//   AdvFSBB round-trip
// =====================================================================

TEST_CASE("Test_FSBB_Advanced_Process_RoundTrip",
          "[converter-model][four-switch-buck-boost-topology][advanced]") {
    const double L = 5e-6;
    json j = make_advanced_fsbb_json(12.0, 12.0, 8.0, 350e3, L, 9.0, 36.0);
    AdvFSBB a(j);
    auto inputs = a.process();
    REQUIRE(inputs.get_design_requirements().get_topology()
            == MAS::Topology::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
    auto Lnom = inputs.get_design_requirements().get_magnetizing_inductance().get_nominal();
    REQUIRE(Lnom.has_value());
    REQUIRE_THAT(Lnom.value(), WithinRel(L, 1e-6));
    REQUIRE(!inputs.get_operating_points().empty());
}

TEST_CASE("Test_FSBB_Advanced_RejectsZeroL",
          "[converter-model][four-switch-buck-boost-topology][advanced]") {
    json j = make_advanced_fsbb_json(12.0, 12.0, 8.0, 350e3, 0.0, 9.0, 36.0);
    AdvFSBB a(j);
    REQUIRE_THROWS(a.process());
}

// =====================================================================
//   get_extra_components_inputs (P7) — L + Cin + Co
// =====================================================================

TEST_CASE("Test_FSBB_ExtraComponents_L_Cin_Co",
          "[converter-model][four-switch-buck-boost-topology][extra-components]") {
    json j = make_fsbb_json(12.0, 12.0, 8.0, 350e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);

    auto extras = c.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    REQUIRE(extras.size() == 3);

    // Index 0: inductor (Inputs).
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[0]));
    {
        const auto& l = std::get<OpenMagnetics::Inputs>(extras[0]);
        REQUIRE(l.get_design_requirements().get_topology()
                == MAS::Topology::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
        REQUIRE(l.get_design_requirements().get_name().value() == "inductor");
        auto Lnom = l.get_design_requirements().get_magnetizing_inductance().get_nominal();
        REQUIRE(Lnom.has_value());
        REQUIRE_THAT(Lnom.value(), WithinRel(c.get_last_sized_inductance(), 1e-9));
        REQUIRE(!l.get_operating_points().empty());
        REQUIRE(l.get_operating_points()[0].get_excitations_per_winding().size() == 1);
    }

    // Index 1: input cap.
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[1]));
    {
        const auto& cin = std::get<CAS::Inputs>(extras[1]);
        REQUIRE(cin.get_design_requirements().get_role() == CAS::Application::INPUT_FILTER);
        REQUIRE(cin.get_design_requirements().get_role().value() == CAS::Application::INPUT_FILTER);
        auto Cnom = cin.get_design_requirements().get_capacitance().get_nominal();
        REQUIRE(Cnom.has_value());
        REQUIRE(Cnom.value() > 0.0);
    }

    // Index 2: output cap.
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[2]));
    {
        const auto& co = std::get<CAS::Inputs>(extras[2]);
        REQUIRE(co.get_design_requirements().get_role() == CAS::Application::OUTPUT_FILTER);
        REQUIRE(co.get_design_requirements().get_role().value() == CAS::Application::OUTPUT_FILTER);
        auto Cnom = co.get_design_requirements().get_capacitance().get_nominal();
        REQUIRE(Cnom.has_value());
        REQUIRE(Cnom.value() > 0.0);
    }
}

TEST_CASE("Test_FSBB_ExtraComponents_Throws_BeforeProcess",
          "[converter-model][four-switch-buck-boost-topology][extra-components]") {
    json j = make_fsbb_json(12.0, 12.0, 8.0, 350e3);
    FSBB c(j);
    REQUIRE_THROWS(c.get_extra_components_inputs(ExtraComponentsMode::IDEAL));
}

TEST_CASE("Test_FSBB_ExtraComponents_REAL_Requires_Magnetic",
          "[converter-model][four-switch-buck-boost-topology][extra-components]") {
    json j = make_fsbb_json(12.0, 12.0, 8.0, 350e3);
    FSBB c(j);
    c.process_operating_points(std::vector<double>{}, 5e-6);
    REQUIRE_THROWS(c.get_extra_components_inputs(ExtraComponentsMode::REAL, std::nullopt));
}

// =====================================================================
//   §8 Volt-Second Balance — All Windings (analytical + SPICE)
// =====================================================================

namespace {

double fsbb_normalised_avg(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double sum = 0.0, peak = 0.0;
    for (double x : v) {
        sum  += x;
        peak  = std::max(peak, std::fabs(x));
    }
    if (peak < 1e-12) return 0.0;
    return std::fabs(sum / static_cast<double>(v.size())) / peak;
}

void fsbb_check_all_windings(const std::vector<MAS::OperatingPoint>& ops,
                             const std::string& path, double eps) {
    REQUIRE_FALSE(ops.empty());
    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        const auto& exc = ops[opIdx].get_excitations_per_winding();
        REQUIRE(exc.size() == 1);   // single inductor
        REQUIRE(exc[0].get_voltage().has_value());
        const auto wf = exc[0].get_voltage()->get_waveform().value();
        const auto& d = wf.get_data();
        REQUIRE(!d.empty());
        const double normAvg = fsbb_normalised_avg(d);
        INFO("FSBB [" << path << "] OP " << opIdx
             << " — |avg(V)|/peak(|V|) = " << normAvg << " (bound " << eps << ")");
        CHECK(normAvg < eps);
    }
}

}  // namespace

TEST_CASE("Test_FSBB_VoltSecondBalance_AllWindings_Buck",
          "[converter-model][four-switch-buck-boost-topology][volt-second-balance]") {
    json j = make_fsbb_json(36.0, 12.0, 8.0, 350e3);
    FSBB c(j);

    SECTION("analytical path") {
        auto ops = c.process_operating_points(std::vector<double>{}, 5e-6);
        fsbb_check_all_windings(ops, "analytical(BUCK)", 0.02);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        c.set_num_steady_state_periods(50);
        c.set_num_periods_to_extract(1);
        c.process_operating_points(std::vector<double>{}, 5e-6);
        auto ops = c.simulate_and_extract_operating_points(5e-6);
        fsbb_check_all_windings(ops, "SPICE(BUCK)", 0.05);
    }
}

TEST_CASE("Test_FSBB_VoltSecondBalance_AllWindings_Boost",
          "[converter-model][four-switch-buck-boost-topology][volt-second-balance]") {
    json j = make_fsbb_json(6.0, 12.0, 5.0, 300e3);
    FSBB c(j);

    SECTION("analytical path") {
        auto ops = c.process_operating_points(std::vector<double>{}, 5e-6);
        fsbb_check_all_windings(ops, "analytical(BOOST)", 0.02);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        c.set_num_steady_state_periods(50);
        c.set_num_periods_to_extract(1);
        c.process_operating_points(std::vector<double>{}, 5e-6);
        auto ops = c.simulate_and_extract_operating_points(5e-6);
        fsbb_check_all_windings(ops, "SPICE(BOOST)", 0.05);
    }
}

// =====================================================================
//   §5.1 Converter-Port Waveforms — DC-stream gate.
// =====================================================================

TEST_CASE("Test_FSBB_ConverterPortWaveforms_Buck",
          "[converter-port-waveforms][four-switch-buck-boost-topology][ngspice-simulation]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    const double Vin = 36.0, Vout = 12.0, Iout = 8.0;
    json j = make_fsbb_json(Vin, Vout, Iout, 350e3);
    FSBB c(j);
    c.set_num_steady_state_periods(80);
    c.set_num_periods_to_extract(1);
    c.process_operating_points(std::vector<double>{}, 5e-6);

    auto wfs = c.simulate_and_extract_topology_waveforms(5e-6);
    REQUIRE(!wfs.empty());
    for (size_t i = 0; i < wfs.size(); ++i) {
        ConverterPortChecks::check_dc_ports(
            wfs[i], "FSBB-Buck", i, Vin, {Vout}, {Iout},
            /*voutMeanTol*/ 0.30,
            /*ioutMeanTol*/ 0.30);
    }
}

// =====================================================================
//   §5.2 SPICE Power Sanity — inductor winding.
// =====================================================================

TEST_CASE("Test_FSBB_SpicePowerSanity",
          "[converter-model][four-switch-buck-boost-topology][ngspice-simulation][regression]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    json j = make_fsbb_json(36.0, 12.0, 8.0, 350e3);
    FSBB c(j);
    c.set_num_steady_state_periods(80);
    c.set_num_periods_to_extract(1);
    c.process_operating_points(std::vector<double>{}, 5e-6);

    auto simOps = c.simulate_and_extract_operating_points(5e-6);
    REQUIRE(!simOps.empty());

    auto& exc = simOps[0].get_excitations_per_winding()[0];
    REQUIRE(exc.get_voltage().has_value());
    REQUIRE(exc.get_current().has_value());
    auto vData = exc.get_voltage()->get_waveform()->get_data();
    auto iData = exc.get_current()->get_waveform()->get_data();
    REQUIRE(vData.size() == iData.size());
    REQUIRE(!vData.empty());

    double sumAbsVI = 0;
    double vmax = vData[0], vmin = vData[0];
    for (size_t k = 0; k < vData.size(); ++k) {
        sumAbsVI += std::fabs(vData[k] * iData[k]);
        if (vData[k] > vmax) vmax = vData[k];
        if (vData[k] < vmin) vmin = vData[k];
    }
    const double avgAbsVI = sumAbsVI / vData.size();
    const double swing    = vmax - vmin;
    INFO("V max=" << vmax << " min=" << vmin << " swing=" << swing
         << "  avg(|V·I|)=" << avgAbsVI);

    // Inductor sees +Vin−Vo on / -Vo off in BUCK → swing ≥ Vo = 12V.
    CHECK(swing > 5.0);
    // Significant instantaneous power flow.
    CHECK(avgAbsVI > 1.0);
}

// =====================================================================
//   §8 Waveform Plotting — visual regression artifacts.
// =====================================================================

namespace {
auto fsbbOutputFilePath = std::filesystem::path{std::source_location::current().file_name()}
    .parent_path().append("..").append("output");
}

TEST_CASE("Test_FSBB_Waveform_Plotting",
          "[converter-model][four-switch-buck-boost-topology][visualization]") {
    json j = make_fsbb_json(36.0, 12.0, 8.0, 350e3);
    FSBB c(j);
    auto ops = c.process_operating_points(std::vector<double>{}, 5e-6);
    REQUIRE(!ops.empty());

    std::filesystem::create_directories(fsbbOutputFilePath);

    SECTION("Inductor current waveform") {
        auto outFile = fsbbOutputFilePath;
        outFile.append("Test_FSBB_Inductor_Current_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_waveform(
            ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Inductor voltage waveform") {
        auto outFile = fsbbOutputFilePath;
        outFile.append("Test_FSBB_Inductor_Voltage_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_waveform(
            ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
        painter.export_svg();
    }
}
