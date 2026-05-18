/**
 * @file TestTopologyWeinberg.cpp
 * @brief Unit tests for the Weinberg DC-DC converter — V1 classic
 *        push-pull primary, V2 H-bridge primary variant, and optional
 *        synchronous-rectifier secondary.
 *
 * Coverage map (WEINBERG_PLAN.md §9, mirrors Zeta §9):
 *
 *   * Static analytical helpers
 *       - calculate_conversion_ratio_boost / _buck
 *       - calculate_overlap_fraction = max(0, 2D-1)
 *       - calculate_switch_peak_voltage_classic / _bridge
 *       - calculate_l1_min, calculate_dcm_K_boost, calculate_rhp_zero_freq
 *       - detect_operating_regime (0=buck, 1=boundary, 2=boost)
 *
 *   * Reference designs (analytical scalar checks)
 *       - Weinberg-Schreuders 1985 (50V→150V@1.5kW@50kHz, n=1/3)
 *       - Yadav-Gowrishankara 2016 (42V→100V@500W@100kHz, n=1/1.4)
 *       - IJRTE BDR (24V→50V@250W@100kHz, n=1/2)
 *
 *   * Per-OP diagnostics population
 *
 *   * Volt-second balance — analytical + SPICE, V1 + V2
 *
 *   * Input-validation guards: Vin≤0, Vo≤0, n≤0, D≥0.95, Iout=0
 *
 *   * SPICE netlist string contract (WEINBERG_PLAN.md §6.2)
 *       - METHOD=GEAR, UIC, ≥3 .ic= occurrences
 *       - K_in 0.999 (L1 coupling), K_xfmr 0.99 (main xfmr coupling, NEVER 1.0)
 *       - Llk_pa/Llk_pb/Llk_sa/Llk_sb (V1) or Llk_pri/Llk_sa/Llk_sb (V2)
 *       - V1: Q1/Q2 push-pull + D3a/D3b energy-recovery diodes
 *       - V2: QA1/QA2/QB1/QB2 H-bridge, NO D3
 *       - Sync rectifier: SW2 model + S_pos/S_neg, no D_pos/D_neg DIDEAL
 *       - Bvab probe (V_ab differential)
 *
 *   * V1 vs V2 switch peak voltage: 2·Vin/(1-D) vs Vin
 *
 *   * V1 default DR: 2 windings + turnsRatio[0] > 0
 *
 *   * AdvancedWeinberg round-trip
 *
 *   * get_extra_components_inputs (P7): 2 entries (L1 Inputs, Co CAS::Inputs OUTPUT_FILTER)
 */

#include "converter_models/Weinberg.h"
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
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

    // Build a minimal Weinberg JSON for one (Vin, Vo, Iout, fsw) operating point.
    json make_weinberg_json(double Vin, double Vo, double Iout, double fsw,
                            double diodeVdrop = 0.7, double efficiency = 0.85,
                            double rippleRatio = 0.30) {
        json j;
        j["inputVoltage"]         = {{"nominal", Vin}};
        j["diodeVoltageDrop"]     = diodeVdrop;
        j["currentRippleRatio"]   = rippleRatio;
        j["efficiency"]           = efficiency;
        j["maximumSwitchCurrent"] = 200.0;

        j["operatingPoints"] = json::array();
        json op;
        op["outputVoltages"]      = {Vo};
        op["outputCurrents"]      = {Iout};
        op["switchingFrequency"]  = fsw;
        op["ambientTemperature"]  = 25.0;
        j["operatingPoints"].push_back(op);
        return j;
    }

    json make_advanced_weinberg_json(double Vin, double Vo, double Iout, double fsw,
                                     double desiredL1, double desiredN) {
        json j = make_weinberg_json(Vin, Vo, Iout, fsw);
        j["desiredInductance"]  = desiredL1;
        j["desiredTurnsRatio"]  = desiredN;
        return j;
    }

}  // namespace

// =====================================================================
//   Static analytical helpers
// =====================================================================

TEST_CASE("Test_Weinberg_ConversionRatio_Boost_Formula",
          "[converter-model][weinberg-topology][unit]") {
    // Boost regime: M = 1/(2·n·(1-D))
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_conversion_ratio_boost(0.6, 1.0/3.0),
                 WithinRel(1.0 / (2.0 * (1.0/3.0) * 0.4), 1e-12));
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_conversion_ratio_boost(0.75, 0.5),
                 WithinRel(1.0 / (2.0 * 0.5 * 0.25), 1e-12));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_conversion_ratio_boost(0.6, /*n*/ 0.0));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_conversion_ratio_boost(/*D*/ 1.0, 0.5));
}

TEST_CASE("Test_Weinberg_ConversionRatio_Buck_Formula",
          "[converter-model][weinberg-topology][unit]") {
    // Buck regime: M = 2D/n  (full-wave CT-PP rectification, continuous with
    // boost branch at D=0.5 where M = 1/n).
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_conversion_ratio_buck(0.4, 0.5),
                 WithinRel(2.0 * 0.4 / 0.5, 1e-12));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_conversion_ratio_buck(0.4, 0.0));
}

TEST_CASE("Test_Weinberg_OverlapFraction",
          "[converter-model][weinberg-topology][unit]") {
    // overlap = max(0, 2D-1)
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_overlap_fraction(0.3), WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_overlap_fraction(0.5), WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_overlap_fraction(0.6), WithinRel(0.2, 1e-12));
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_overlap_fraction(0.9), WithinRel(0.8, 1e-12));
}

TEST_CASE("Test_Weinberg_SwitchPeakVoltage_V1_vs_V2",
          "[converter-model][weinberg-topology][unit]") {
    // V1 classic: 2·Vin/(1-D); V2 bridge: Vin (halved)
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_switch_peak_voltage_classic(50.0, 0.6),
                 WithinRel(2.0 * 50.0 / 0.4, 1e-12));
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_switch_peak_voltage_bridge(50.0),
                 WithinRel(50.0, 1e-12));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_switch_peak_voltage_classic(50.0, /*D*/ 1.0));
}

TEST_CASE("Test_Weinberg_DetectOperatingRegime",
          "[converter-model][weinberg-topology][unit]") {
    REQUIRE(OpenMagnetics::Weinberg::detect_operating_regime(0.3) == 0);   // buck
    REQUIRE(OpenMagnetics::Weinberg::detect_operating_regime(0.5) == 1);   // boundary
    REQUIRE(OpenMagnetics::Weinberg::detect_operating_regime(0.7) == 2);   // boost
}

TEST_CASE("Test_Weinberg_L1_Min_Sizing",
          "[converter-model][weinberg-topology][unit]") {
    // Boost regime (overlap dominates): L1 = Vin · (2D-1) / (4·ΔIL1·fsw)
    double L1 = OpenMagnetics::Weinberg::calculate_l1_min(/*Vin*/ 50.0, /*D*/ 0.6,
                                            /*ΔIL1*/ 1.0, /*fsw*/ 50e3);
    REQUIRE_THAT(L1, WithinRel(50.0 * 0.2 / (4.0 * 1.0 * 50e3), 1e-12));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_l1_min(50.0, 0.6, /*ΔIL1*/ 0.0, 50e3));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_l1_min(50.0, 0.6, 1.0, /*fsw*/ 0.0));
}

TEST_CASE("Test_Weinberg_DCM_K_Boost_Formula",
          "[converter-model][weinberg-topology][unit]") {
    // K = 2·L1·fsw·(1-D)² / R
    double L1 = 200e-6, fsw = 50e3, D = 0.6, R = 15.0;
    double Kexp = 2.0 * L1 * fsw * (1.0 - D) * (1.0 - D) / R;
    REQUIRE_THAT(OpenMagnetics::Weinberg::calculate_dcm_K_boost(L1, fsw, D, R),
                 WithinRel(Kexp, 1e-12));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_dcm_K_boost(/*L1*/ 0.0, fsw, D, R));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_dcm_K_boost(L1, fsw, D, /*R*/ 0.0));
}

// =====================================================================
//   Input-validation guards
// =====================================================================

TEST_CASE("Test_Weinberg_DutyCycle_BadInput_Throws",
          "[converter-model][weinberg-topology][unit]") {
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_duty_cycle(/*Vin*/ 0.0, 150.0, 1.0/3.0, 0.7, 0.85));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_duty_cycle(/*Vin*/ -1.0, 150.0, 1.0/3.0, 0.7, 0.85));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_duty_cycle(50.0, /*Vo*/ 0.0, 1.0/3.0, 0.7, 0.85));
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_duty_cycle(50.0, 150.0, /*n*/ 0.0, 0.7, 0.85));
}

TEST_CASE("Test_Weinberg_DutyCycle_Clamp_Throws_Above_095",
          "[converter-model][weinberg-topology][unit]") {
    // Boost branch: D = 1 - 1/(2·n·M).  With Vin small, Vo huge, n small,
    // we force D very close to 1. Vin=1, Vo=1000, n=1, η=1, Vd=0:
    //   M = 1000  → D = 1 - 1/2000 = 0.9995  → throws.
    REQUIRE_THROWS(OpenMagnetics::Weinberg::calculate_duty_cycle(/*Vin*/ 1.0, /*Vo*/ 1000.0,
                                                  /*n*/ 1.0, 0.0, 1.0));
}

TEST_CASE("Test_Weinberg_DivideByZero_Iout_Guard",
          "[converter-model][weinberg-topology][unit]") {
    json j = make_weinberg_json(50.0, 150.0, /*Iout*/ 10.0, 50e3);
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{1.0/3.0}, /*L1*/ 200e-6);

    auto ops = w.get_mutable_operating_points();
    ops[0].set_output_currents({0.0});
    w.set_operating_points(ops);
    REQUIRE_THROWS(w.generate_ngspice_circuit(/*n*/ 1.0/3.0, /*L1*/ 200e-6));
}

// =====================================================================
//   Reference designs (analytical scalar checks — no SPICE)
// =====================================================================

TEST_CASE("Test_Weinberg_Schreuders_1985_Reference_Design",
          "[converter-model][weinberg-topology][refdesign][analytical]") {
    // Weinberg-Schreuders PESC ~1985: 50 V → 150 V @ 1.5 kW @ 50 kHz, n = 1/3.
    // Lossless (η=1, Vd=0) → M = 3, D = 1 − 1/(2·(1/3)·3) = 1 − 1/2 = 0.5.
    // (This is the boost-buck boundary.)
    json j = make_weinberg_json(/*Vin*/ 50.0, /*Vo*/ 150.0, /*Iout*/ 10.0,
                                /*fsw*/ 50e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{1.0/3.0}, /*L1*/ 200e-6);

    double D = w.get_last_duty_cycle();
    REQUIRE_THAT(D, WithinRel(0.5, 1e-9));
    // M at boundary = 1/n = 3 (boundary uses buck branch since regime==1).
    REQUIRE_THAT(w.get_last_conversion_ratio(), WithinRel(3.0, 1e-9));
    // V1 switch peak voltage at D=0.5: 2·50/0.5 = 200 V
    REQUIRE_THAT(w.get_last_switch_peak_voltage(), WithinRel(200.0, 1e-9));
    // Diode peak reverse voltage = 2·Vout = 300 V (CT-FW).
    REQUIRE_THAT(w.get_last_diode_peak_reverse_voltage(), WithinRel(300.0, 1e-9));
}

TEST_CASE("Test_Weinberg_Yadav_Gowrishankara_2016_Reference_Design",
          "[converter-model][weinberg-topology][refdesign][analytical]") {
    // Yadav-Gowrishankara 2016: 42 V → 100 V @ 500 W @ 100 kHz, n = 1/1.4.
    // Lossless: M = 100/42 ≈ 2.381; boost branch D = 1 − 1/(2·(1/1.4)·M).
    double n = 1.0/1.4, Vin = 42.0, Vo = 100.0;
    double Mexp = Vo / Vin;
    double Dexp = 1.0 - 1.0 / (2.0 * n * Mexp);
    REQUIRE(Dexp > 0.5);

    json j = make_weinberg_json(Vin, Vo, /*Iout*/ 5.0,
                                /*fsw*/ 100e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{n}, /*L1*/ 100e-6);

    REQUIRE_THAT(w.get_last_duty_cycle(),       WithinRel(Dexp, 1e-9));
    REQUIRE_THAT(w.get_last_conversion_ratio(), WithinRel(Mexp, 1e-9));
    REQUIRE(w.get_last_operating_regime() == 2);   // boost
}

TEST_CASE("Test_Weinberg_IJRTE_BDR_Reference_Design",
          "[converter-model][weinberg-topology][refdesign][analytical]") {
    // IJRTE 2020-style Weinberg BDR: 24 V → 50 V @ 250 W @ 100 kHz, n = 1/2.
    double n = 0.5, Vin = 24.0, Vo = 50.0;
    double Mexp = Vo / Vin;
    double Dexp = 1.0 - 1.0 / (2.0 * n * Mexp);

    json j = make_weinberg_json(Vin, Vo, /*Iout*/ 5.0,
                                /*fsw*/ 100e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{n}, /*L1*/ 80e-6);

    REQUIRE_THAT(w.get_last_duty_cycle(),       WithinRel(Dexp, 1e-9));
    REQUIRE_THAT(w.get_last_conversion_ratio(), WithinRel(Mexp, 1e-9));
}

// =====================================================================
//   Per-OP diagnostics populated
// =====================================================================

TEST_CASE("Test_Weinberg_Diagnostics_Populated",
          "[converter-model][weinberg-topology][unit]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{1.0/3.0}, /*L1*/ 200e-6);

    REQUIRE(w.get_last_duty_cycle() > 0.0);
    REQUIRE(w.get_last_duty_cycle() < 1.0);
    REQUIRE(w.get_last_conversion_ratio() > 0.0);
    REQUIRE(w.get_last_overlap_fraction() >= 0.0);
    REQUIRE(w.get_last_switch_peak_voltage() > 0.0);
    REQUIRE(w.get_last_switch_peak_current() > 0.0);
    REQUIRE(w.get_last_diode_peak_reverse_voltage() > 0.0);
    REQUIRE(w.get_last_diode_peak_current() > 0.0);
    REQUIRE(w.get_last_input_inductor_average() > 0.0);
    REQUIRE(w.get_last_input_inductor_ripple() > 0.0);
    REQUIRE(w.get_last_dcm_k() > 0.0);
    REQUIRE(w.get_last_dcm_kcrit() > 0.0);
    REQUIRE(w.get_last_sized_l1() > 0.0);
    REQUIRE(w.get_last_sized_co() > 0.0);
    REQUIRE(w.get_last_output_voltage_ripple() > 0.0);
}

// =====================================================================
//   process_design_requirements
// =====================================================================

TEST_CASE("Test_Weinberg_V1_Default_DR_Has_Two_Windings",
          "[converter-model][weinberg-topology][v1-default][unit]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    OpenMagnetics::Weinberg w(j);
    auto dr = w.process_design_requirements();

    REQUIRE(dr.get_topology() == Topologies::WEINBERG_CONVERTER);
    REQUIRE(dr.get_isolation_sides().has_value());
    REQUIRE(dr.get_isolation_sides().value().size() == 2);
    REQUIRE(dr.get_turns_ratios().size() == 1);
    auto Nratio = dr.get_turns_ratios()[0].get_nominal();
    REQUIRE(Nratio.has_value());
    REQUIRE(Nratio.value() > 0.0);
}

TEST_CASE("Test_Weinberg_DR_Throws_If_No_Ripple_Or_Imax",
          "[converter-model][weinberg-topology][unit]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    j.erase("currentRippleRatio");
    j.erase("maximumSwitchCurrent");
    OpenMagnetics::Weinberg w(j);
    REQUIRE_THROWS(w.process_design_requirements());
}

// =====================================================================
//   SPICE netlist string contract — V1 classic (default)
// =====================================================================

TEST_CASE("Test_Weinberg_V1_SPICE_Netlist_StringContract",
          "[converter-model][weinberg-topology][spice-netlist][v1]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{1.0/3.0}, /*L1*/ 200e-6);

    std::string netlist = w.generate_ngspice_circuit(/*n*/ 1.0/3.0, /*L1*/ 200e-6);

    // GEAR + UIC for .ic
    REQUIRE(netlist.find("METHOD=GEAR") != std::string::npos);
    REQUIRE(netlist.find("UIC") != std::string::npos);
    REQUIRE(netlist.find("V1 classic push-pull") != std::string::npos);
    REQUIRE(netlist.find("/CT-FW diode") != std::string::npos);

    // Couplings
    REQUIRE(netlist.find("K_in L1a L1b 0.999") != std::string::npos);
    // V1 main xfmr: 6 pairwise K statements at k=0.99 (NEVER 1.0)
    REQUIRE(netlist.find("K_pa_pb Lpri_a Lpri_b 0.99") != std::string::npos);
    REQUIRE(netlist.find("K_pa_sa Lpri_a Lsec_a 0.99") != std::string::npos);
    REQUIRE(netlist.find("K_pb_sb Lpri_b Lsec_b 0.99") != std::string::npos);
    REQUIRE(netlist.find("K_sa_sb Lsec_a Lsec_b 0.99") != std::string::npos);
    REQUIRE(netlist.find("K_pa_pb Lpri_a Lpri_b 1") == std::string::npos);

    // Two L1 windings
    REQUIRE(netlist.find("L1a vin_p") != std::string::npos);
    REQUIRE(netlist.find("L1b vin_p") != std::string::npos);

    // V1 push-pull primary windings + leakage
    REQUIRE(netlist.find("Lpri_a") != std::string::npos);
    REQUIRE(netlist.find("Lpri_b") != std::string::npos);
    REQUIRE(netlist.find("Llk_pa") != std::string::npos);
    REQUIRE(netlist.find("Llk_pb") != std::string::npos);
    // CT-FW secondary windings + leakage
    REQUIRE(netlist.find("Lsec_a") != std::string::npos);
    REQUIRE(netlist.find("Lsec_b") != std::string::npos);
    REQUIRE(netlist.find("Llk_sa") != std::string::npos);
    REQUIRE(netlist.find("Llk_sb") != std::string::npos);

    // V1 push-pull switches (S1 phase=0, S2 phase=Tsw/2) + snubbers.
    // Use "S" prefix (ngspice voltage-controlled SW), not "Q" (BJT).
    REQUIRE(netlist.find("S1 drainQ1 0 pwm1 0 SW1") != std::string::npos);
    REQUIRE(netlist.find("S2 drainQ2 0 pwm2 0 SW1") != std::string::npos);
    REQUIRE(netlist.find("Vpwm1") != std::string::npos);
    REQUIRE(netlist.find("Vpwm2") != std::string::npos);
    REQUIRE(netlist.find("Rsnub_s1") != std::string::npos);
    REQUIRE(netlist.find("Csnub_s2") != std::string::npos);

    // D3a/D3b (Schreuders' literature drain→Vin clamp) are intentionally
    // OMITTED from V1 SPICE — the literal anode-at-drain topology clamps the
    // legitimate 2·Vin/(1−D) push-pull swing to Vin+0.7 and breaks power
    // transfer. The literature uses an aux-winding clamp reference that we
    // do not model. Snubbers handle the leakage spike adequately.
    REQUIRE(netlist.find("D3a drainQ1 vin_p DIDEAL") == std::string::npos);
    REQUIRE(netlist.find("D3b drainQ2 vin_p DIDEAL") == std::string::npos);

    // Bvab probe
    REQUIRE(netlist.find("Bvab vab 0") != std::string::npos);

    // CT-FW diodes (default rectifier)
    REQUIRE(netlist.find("D_pos sec_pos_in out_node DIDEAL") != std::string::npos);
    REQUIRE(netlist.find("D_neg sec_neg_in out_node DIDEAL") != std::string::npos);
    // No SR (default)
    REQUIRE(netlist.find(".model SW2") == std::string::npos);

    // ≥3 .ic= occurrences (L1a, L1b, Co)
    size_t icCount = 0; size_t pos = 0;
    while ((pos = netlist.find(" ic=", pos)) != std::string::npos) { ++icCount; ++pos; }
    REQUIRE(icCount >= 3);
}

// =====================================================================
//   SPICE netlist string contract — V2 H-bridge primary
// =====================================================================

TEST_CASE("Test_Weinberg_V2_Bridge_SPICE_Netlist_StringContract",
          "[converter-model][weinberg-topology][spice-netlist][v2]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    j["variant"] = "bridge";
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{1.0/3.0}, /*L1*/ 200e-6);

    std::string netlist = w.generate_ngspice_circuit(/*n*/ 1.0/3.0, /*L1*/ 200e-6);

    REQUIRE(netlist.find("V2 bridge primary") != std::string::npos);

    // Couplings (3 pairwise K statements for V2; single primary Lpri)
    REQUIRE(netlist.find("K_in L1a L1b 0.999") != std::string::npos);
    REQUIRE(netlist.find("K_pri_sa Lpri Lsec_a 0.99") != std::string::npos);
    REQUIRE(netlist.find("K_pri_sb Lpri Lsec_b 0.99") != std::string::npos);
    REQUIRE(netlist.find("K_sa_sb Lsec_a Lsec_b 0.99") != std::string::npos);

    // 4 H-bridge switches (S-prefix for ngspice SW elements).
    REQUIRE(netlist.find("SA1 priCT_a priLeft  pwm1 0 SW1") != std::string::npos);
    REQUIRE(netlist.find("SA2 priLeft 0       pwm2 0 SW1") != std::string::npos);
    REQUIRE(netlist.find("SB1 priCT_a priRight pwm2 0 SW1") != std::string::npos);
    REQUIRE(netlist.find("SB2 priRight 0       pwm1 0 SW1") != std::string::npos);

    // V2: NO push-pull S1/S2, NO D3 energy-recovery diodes
    REQUIRE(netlist.find("D3a drainQ1") == std::string::npos);
    REQUIRE(netlist.find("D3b drainQ2") == std::string::npos);
    REQUIRE(netlist.find("\nS1 drainQ1") == std::string::npos);

    // V2: single primary winding Lpri (no Lpri_a/Lpri_b split)
    REQUIRE(netlist.find("Lpri pri_top priRight") != std::string::npos);
    REQUIRE(netlist.find("Llk_pri") != std::string::npos);
    // CT-FW secondary windings still present
    REQUIRE(netlist.find("Lsec_a") != std::string::npos);
    REQUIRE(netlist.find("Lsec_b") != std::string::npos);
}

// =====================================================================
//   Synchronous rectifier
// =====================================================================

TEST_CASE("Test_Weinberg_Synchronous_Netlist_Replaces_Diodes",
          "[converter-model][weinberg-topology][synchronous][unit]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    j["synchronousRectifier"] = true;
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{1.0/3.0}, /*L1*/ 200e-6);

    std::string netlist = w.generate_ngspice_circuit(/*n*/ 1.0/3.0, /*L1*/ 200e-6);

    // No CT-FW DIDEAL D_pos/D_neg
    REQUIRE(netlist.find("D_pos sec_pos_in out_node DIDEAL") == std::string::npos);
    REQUIRE(netlist.find("D_neg sec_neg_in out_node DIDEAL") == std::string::npos);

    // SR: SW2 model + S_pos / S_neg
    REQUIRE(netlist.find(".model SW2") != std::string::npos);
    REQUIRE(netlist.find("S_pos sec_pos_in out_node pwm1 0 SW2") != std::string::npos);
    REQUIRE(netlist.find("S_neg sec_neg_in out_node pwm2 0 SW2") != std::string::npos);
    REQUIRE(netlist.find("/synchronous SR") != std::string::npos);
}

// =====================================================================
//   V1 vs V2 switch peak voltage diagnostic comparison
// =====================================================================

TEST_CASE("Test_Weinberg_V1_vs_V2_SwitchPeakVoltage",
          "[converter-model][weinberg-topology][unit]") {
    // V1 switch sees 2·Vin/(1-D); V2 switch sees Vin only.
    {
        json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3,
                                    /*Vd*/ 0.0, /*η*/ 1.0);
        OpenMagnetics::Weinberg v1(j);
        v1.process_operating_points(std::vector<double>{1.0/3.0}, 200e-6);
        REQUIRE_THAT(v1.get_last_switch_peak_voltage(),
                     WithinRel(2.0 * 50.0 / (1.0 - v1.get_last_duty_cycle()), 1e-9));
    }
    {
        json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3,
                                    /*Vd*/ 0.0, /*η*/ 1.0);
        j["variant"] = "bridge";
        OpenMagnetics::Weinberg v2(j);
        v2.process_operating_points(std::vector<double>{1.0/3.0}, 200e-6);
        REQUIRE_THAT(v2.get_last_switch_peak_voltage(), WithinRel(50.0, 1e-9));
    }
}

// =====================================================================
//   AdvancedWeinberg round-trip
// =====================================================================

TEST_CASE("Test_Weinberg_AdvancedWeinberg_Process_RoundTrip",
          "[converter-model][weinberg-topology][advanced]") {
    const double L1 = 200e-6;
    const double n  = 1.0/3.0;
    json j = make_advanced_weinberg_json(50.0, 150.0, 10.0, 50e3, L1, n);
    OpenMagnetics::AdvancedWeinberg aw(j);
    auto inputs = aw.process();
    REQUIRE(inputs.get_design_requirements().get_topology() == Topologies::WEINBERG_CONVERTER);
    auto Lnom = inputs.get_design_requirements().get_magnetizing_inductance().get_nominal();
    REQUIRE(Lnom.has_value());
    REQUIRE_THAT(Lnom.value(), WithinRel(L1, 1e-6));
    REQUIRE(inputs.get_design_requirements().get_turns_ratios().size() == 1);
    auto Nnom = inputs.get_design_requirements().get_turns_ratios()[0].get_nominal();
    REQUIRE(Nnom.has_value());
    REQUIRE_THAT(Nnom.value(), WithinRel(n, 1e-3));
    REQUIRE(inputs.get_operating_points().size() == 1);
    // process_design_requirements declares 2 isolation sides (combined Pri +
    // combined Sec), so process_operating_points emits 2 excitations.
    REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 2);
}

// =====================================================================
//   get_extra_components_inputs (P7) — L1 + Co
// =====================================================================

TEST_CASE("Test_Weinberg_ExtraComponents_L1_Co",
          "[converter-model][weinberg-topology][extra-components]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{1.0/3.0}, /*L1*/ 200e-6);

    auto extras = w.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    REQUIRE(extras.size() == 2);

    // Index 0: L1 input coupled inductor (2 windings 1:1).
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[0]));
    {
        const auto& l1 = std::get<OpenMagnetics::Inputs>(extras[0]);
        REQUIRE(l1.get_design_requirements().get_topology() == Topologies::WEINBERG_CONVERTER);
        REQUIRE(l1.get_design_requirements().get_name().value() == "inputCoupledInductor");
        auto Lnom = l1.get_design_requirements().get_magnetizing_inductance().get_nominal();
        REQUIRE(Lnom.has_value());
        REQUIRE_THAT(Lnom.value(), WithinRel(w.get_last_sized_l1(), 1e-9));
        REQUIRE(l1.get_design_requirements().get_turns_ratios().size() == 1);
        auto Nratio = l1.get_design_requirements().get_turns_ratios()[0].get_nominal();
        REQUIRE(Nratio.has_value());
        REQUIRE_THAT(Nratio.value(), WithinRel(1.0, 1e-12));
        REQUIRE(l1.get_design_requirements().get_isolation_sides().has_value());
        REQUIRE(l1.get_design_requirements().get_isolation_sides().value().size() == 2);
        REQUIRE(l1.get_operating_points().size() == 1);
        REQUIRE(l1.get_operating_points()[0].get_excitations_per_winding().size() == 2);
    }

    // Index 1: Co output cap (OUTPUT_FILTER).
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[1]));
    {
        const auto& co = std::get<CAS::Inputs>(extras[1]);
        REQUIRE(co.get_design_requirements().get_name().value() == "outputCapacitor");
        auto Cnom = co.get_design_requirements().get_capacitance().get_nominal();
        REQUIRE(Cnom.has_value());
        REQUIRE_THAT(Cnom.value(), WithinRel(w.get_last_sized_co(), 1e-9));
        REQUIRE(co.get_design_requirements().get_role() == CAS::Application::OUTPUT_FILTER);
        REQUIRE(co.get_operating_points().size() == 1);
    }
}

TEST_CASE("Test_Weinberg_ExtraComponents_Throws_BeforeProcess",
          "[converter-model][weinberg-topology][extra-components]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    OpenMagnetics::Weinberg w(j);
    REQUIRE_THROWS(w.get_extra_components_inputs(ExtraComponentsMode::IDEAL));
}

TEST_CASE("Test_Weinberg_ExtraComponents_REAL_Requires_Magnetic",
          "[converter-model][weinberg-topology][extra-components]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3);
    OpenMagnetics::Weinberg w(j);
    w.process_operating_points(std::vector<double>{1.0/3.0}, /*L1*/ 200e-6);
    REQUIRE_THROWS(w.get_extra_components_inputs(ExtraComponentsMode::REAL,
                                                  /*magnetic*/ std::nullopt));
}

// =====================================================================
//   §8 Volt-Second Balance — All Windings (analytical + SPICE)
// =====================================================================

namespace {

double weinberg_normalised_avg(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double sum = 0.0, peak = 0.0;
    for (double x : v) {
        sum  += x;
        peak  = std::max(peak, std::fabs(x));
    }
    if (peak < 1e-12) return 0.0;
    return std::fabs(sum / static_cast<double>(v.size())) / peak;
}

void weinberg_check_all_windings(const std::vector<MAS::OperatingPoint>& ops,
                                 const std::string& path,
                                 double eps,
                                 size_t expectedWindings) {
    REQUIRE_FALSE(ops.empty());
    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        const auto& exc = ops[opIdx].get_excitations_per_winding();
        REQUIRE(exc.size() == expectedWindings);
        for (size_t wi = 0; wi < exc.size(); ++wi) {
            REQUIRE(exc[wi].get_voltage().has_value());
            const auto wf = exc[wi].get_voltage()->get_waveform().value();
            const auto& d = wf.get_data();
            REQUIRE(!d.empty());
            const double normAvg = weinberg_normalised_avg(d);
            INFO("Weinberg [" << path << "] OP " << opIdx
                 << " winding " << wi
                 << " — |avg(V)|/peak(|V|) = " << normAvg
                 << " (bound " << eps << ")");
            CHECK(normAvg < eps);
        }
    }
}

}  // namespace

TEST_CASE("Test_Weinberg_VoltSecondBalance_AllWindings_V1",
          "[converter-model][weinberg-topology][volt-second-balance]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Weinberg w(j);

    SECTION("analytical path") {
        auto ops = w.process_operating_points(std::vector<double>{1.0/3.0}, 200e-6);
        weinberg_check_all_windings(ops, "analytical(V1)", 0.02, /*windings*/ 2);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        w.set_num_steady_state_periods(400);
        w.set_num_periods_to_extract(1);
        auto ops = w.simulate_and_extract_operating_points(1.0/3.0, 200e-6);
        weinberg_check_all_windings(ops, "SPICE(V1)", 0.05, /*windings*/ 2);
    }
}

TEST_CASE("Test_Weinberg_VoltSecondBalance_AllWindings_V2_Bridge",
          "[converter-model][weinberg-topology][volt-second-balance]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3, /*Vd*/ 0.0, /*η*/ 1.0);
    j["variant"] = "bridge";
    OpenMagnetics::Weinberg w(j);

    SECTION("analytical path") {
        auto ops = w.process_operating_points(std::vector<double>{1.0/3.0}, 200e-6);
        weinberg_check_all_windings(ops, "analytical(V2)", 0.02, /*windings*/ 2);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        w.set_num_steady_state_periods(400);
        w.set_num_periods_to_extract(1);
        auto ops = w.simulate_and_extract_operating_points(1.0/3.0, 200e-6);
        weinberg_check_all_windings(ops, "SPICE(V2)", 0.05, /*windings*/ 2);
    }
}

// =====================================================================
//   §5.1 Converter-Port Waveforms — DC-stream gate.
// =====================================================================

TEST_CASE("Test_Weinberg_ConverterPortWaveforms",
          "[converter-port-waveforms][weinberg-topology][ngspice-simulation]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    const double Vin = 50.0, Vout = 150.0, Iout = 10.0;
    json j = make_weinberg_json(Vin, Vout, Iout, 50e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Weinberg w(j);
    w.set_num_steady_state_periods(400);
    w.set_num_periods_to_extract(1);

    auto wfs = w.simulate_and_extract_topology_waveforms(1.0/3.0, 200e-6);
    REQUIRE(!wfs.empty());
    for (size_t i = 0; i < wfs.size(); ++i) {
        ConverterPortChecks::check_dc_ports(
            wfs[i], "Weinberg", i, Vin, {Vout}, {Iout},
            /*voutMeanTol*/ 0.35,   // open-loop fixed-D, Weinberg settles slowly
            /*ioutMeanTol*/ 0.35);
    }
}

// =====================================================================
//   §5.2 SPICE Power Sanity — primary winding.
// =====================================================================

TEST_CASE("Test_Weinberg_SpicePowerSanity",
          "[converter-model][weinberg-topology][ngspice-simulation][regression]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Weinberg w(j);
    w.set_num_steady_state_periods(400);
    w.set_num_periods_to_extract(1);

    auto simOps = w.simulate_and_extract_operating_points(1.0/3.0, 200e-6);
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
    const double avgVI    = sumVI    / vData.size();
    const double avgAbsVI = sumAbsVI / vData.size();
    const double swing    = vmax - vmin;

    INFO("V max=" << vmax << " min=" << vmin << " swing=" << swing
         << "  avg(V·I)=" << avgVI << " avg(|V·I|)=" << avgAbsVI);

    // (a) bipolar primary winding swing — V_ab = ±Vin (push-pull) → ≥ 50 V.
    CHECK(swing > 50.0);

    // (b) passive sign — power into the modelled primary winding.
    CHECK(std::fabs(avgVI) < avgAbsVI);    // not collapsed to constant sign

    // (c) realistic magnitude for 1.5 kW fixture.
    CHECK(avgAbsVI > 5.0);
}

// =====================================================================
//   §8 Waveform Plotting — visual regression artifacts.
// =====================================================================

namespace {
auto weinbergOutputFilePath = std::filesystem::path{std::source_location::current().file_name()}
    .parent_path().append("..").append("output");
}

TEST_CASE("Test_Weinberg_Waveform_Plotting",
          "[converter-model][weinberg-topology][visualization]") {
    json j = make_weinberg_json(50.0, 150.0, 10.0, 50e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Weinberg w(j);
    auto ops = w.process_operating_points(std::vector<double>{1.0/3.0}, 200e-6);
    REQUIRE(!ops.empty());

    std::filesystem::create_directories(weinbergOutputFilePath);

    SECTION("Primary current waveform plot") {
        auto outFile = weinbergOutputFilePath;
        outFile.append("Test_Weinberg_Primary_Current_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Primary voltage waveform plot") {
        auto outFile = weinbergOutputFilePath;
        outFile.append("Test_Weinberg_Primary_Voltage_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
        painter.export_svg();
    }
}
