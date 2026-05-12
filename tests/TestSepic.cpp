/**
 * @file TestSepic.cpp
 * @brief Unit tests for the SEPIC (Single-Ended Primary-Inductor Converter)
 *        DC-DC converter model — V1 uncoupled, V2 coupled-inductor and
 *        synchronous-rectifier variants. CCM analytical + DCM-flag coverage,
 *        plus structural SPICE netlist contract.
 *
 * Coverage map (SEPIC_PLAN.md §9):
 *
 *   * Static analytical helpers
 *       - Conversion ratio M(D) = +D/(1-D)
 *       - Coupling-cap voltage VCs = Vin (steady state)
 *       - L1 / L2 / Cs sizing inverses
 *       - K(D) (DCM boundary)
 *
 *   * Reference designs (analytical scalar checks)
 *       - TI SNVA168E §6 worked example  (5 V → 12 V @ 0.5 A, 600 kHz)
 *       - LTC1871 SEPIC (3.3 V → 5 V @ 1 A, 250 kHz)
 *
 *   * Per-OP diagnostics population
 *
 *   * Volt-second balance on L1
 *
 *   * Input-validation guards
 *       - Vin ≤ 0           → throws
 *       - |Vo| ≤ 0          → throws
 *       - D ≥ 0.95          → throws
 *       - Iout = 0 in SPICE → throws
 *
 *   * SPICE netlist string contract (SEPIC_PLAN.md §6.2)
 *       - METHOD=GEAR, .ic on every reactive element, snubbers, probes.
 *       - L2 bottom terminal at GND (the SEPIC distinguishing feature).
 *       - D1 anode at node_B side, cathode at vout (non-inverting output).
 *
 *   * OpenMagnetics::AdvancedSepic round-trip
 *
 *   * get_extra_components_inputs (P7)
 *       - 3 entries: L2 (Inputs), Cs (CAS::Inputs, DC_LINK), Co (CAS::Inputs, OUTPUT_FILTER)
 *
 *   * V2 coupled-inductor: DR shape + K_L1L2 in netlist
 *
 *   * Synchronous rectifier: D1 replaced by S2 + complementary PWM
 */

#include "converter_models/Sepic.h"
#include "processors/Inputs.h"
#include "support/Utils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <variant>

using namespace OpenMagnetics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

    // Build a minimal SEPIC JSON for one (Vin, Vo, Iout, fsw) operating point.
    json make_sepic_json(double Vin, double Vo, double Iout, double fsw,
                         double diodeVdrop = 0.7, double efficiency = 0.85,
                         double rippleRatio = 0.30) {
        json j;
        j["inputVoltage"]         = {{"nominal", Vin}};
        j["diodeVoltageDrop"]     = diodeVdrop;
        j["currentRippleRatio"]   = rippleRatio;
        j["efficiency"]           = efficiency;
        j["maximumSwitchCurrent"] = 100.0;

        j["operatingPoints"] = json::array();
        json op;
        op["outputVoltages"]      = {Vo};
        op["outputCurrents"]      = {Iout};
        op["switchingFrequency"]  = fsw;
        op["ambientTemperature"]  = 25.0;
        j["operatingPoints"].push_back(op);
        return j;
    }

    json make_advanced_sepic_json(double Vin, double Vo, double Iout, double fsw,
                                  double desiredL1) {
        json j = make_sepic_json(Vin, Vo, Iout, fsw);
        j["desiredInductance"] = desiredL1;
        return j;
    }

}  // namespace

// =====================================================================
//   Static analytical helpers
// =====================================================================

TEST_CASE("Test_Sepic_ConversionRatio_Formula",
          "[converter-model][sepic-topology][unit]") {
    // M(D) = +D/(1-D)  (NON-inverting; cf. Cuk which is -D/(1-D)).
    REQUIRE_THAT(OpenMagnetics::Sepic::calculate_conversion_ratio(0.2), WithinRel(0.25, 1e-12));
    REQUIRE_THAT(OpenMagnetics::Sepic::calculate_conversion_ratio(0.5), WithinRel(1.0,  1e-12));
    REQUIRE_THAT(OpenMagnetics::Sepic::calculate_conversion_ratio(0.7), WithinRel(7.0/3.0, 1e-12));
}

TEST_CASE("Test_Sepic_L1_L2_Cs_Inverse_Sizing",
          "[converter-model][sepic-topology][unit]") {
    // L1 = Vin·D / (ΔIL1·fsw)
    double L1 = OpenMagnetics::Sepic::calculate_l1_min(/*Vin*/ 5.0, /*D*/ 0.6, /*ΔIL1*/ 0.5, /*fsw*/ 600e3);
    REQUIRE_THAT(L1, WithinRel(5.0 * 0.6 / (0.5 * 600e3), 1e-12));

    // L2 ALSO uses Vin (not Vo) — distinguishing SEPIC feature versus Cuk.
    double L2 = OpenMagnetics::Sepic::calculate_l2_min(/*Vin*/ 5.0, /*D*/ 0.6, /*ΔIL2*/ 0.3, /*fsw*/ 600e3);
    REQUIRE_THAT(L2, WithinRel(5.0 * 0.6 / (0.3 * 600e3), 1e-12));

    double Cs = OpenMagnetics::Sepic::calculate_cs_min(/*Iout*/ 0.5, /*D*/ 0.6, /*ΔVCs*/ 0.25, /*fsw*/ 600e3);
    REQUIRE_THAT(Cs, WithinRel(0.5 * 0.6 / (0.25 * 600e3), 1e-12));

    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_l1_min(5.0, 0.6, /*ΔIL1*/ 0.0, 600e3));
    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_l2_min(5.0, 0.6, 0.3, /*fsw*/ 0.0));
    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_cs_min(0.5, 0.6, /*ΔVCs*/ -1.0, 600e3));
}

TEST_CASE("Test_Sepic_DCM_K_Formula",
          "[converter-model][sepic-topology][unit]") {
    // K(D) = 2·Le·fsw/R, Le = L1·L2/(L1+L2)
    double L1 = 47e-6, L2 = 47e-6, fsw = 600e3, R = 24.0;
    double Le = L1 * L2 / (L1 + L2);
    double Kexpected = 2.0 * Le * fsw / R;
    REQUIRE_THAT(OpenMagnetics::Sepic::calculate_dcm_K(L1, L2, fsw, R), WithinRel(Kexpected, 1e-12));

    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_dcm_K(0.0, 0.0, fsw, R));
    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_dcm_K(L1,  L2,  fsw, /*R*/ 0.0));
}

// =====================================================================
//   Input-validation guards
// =====================================================================

TEST_CASE("Test_Sepic_DutyCycle_BadInput_Throws",
          "[converter-model][sepic-topology][unit]") {
    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_duty_cycle(/*Vin*/ 0.0, 12.0, 0.7, 0.85));
    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_duty_cycle(/*Vin*/ -1.0, 12.0, 0.7, 0.85));
    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_duty_cycle(5.0, /*Vo*/ 0.0, 0.7, 0.85));
}

TEST_CASE("Test_Sepic_DutyCycle_Clamp_Throws_Above_095",
          "[converter-model][sepic-topology][unit]") {
    // D = (Vo+Vd)/(Vin·η + Vo + Vd). For Vin=1, η=1, Vd=0: D=Vo/(Vo+1).
    // D ≥ 0.95 ⇒ Vo ≥ 19.
    REQUIRE_THROWS(OpenMagnetics::Sepic::calculate_duty_cycle(/*Vin*/ 1.0, /*Vo*/ 100.0, 0.0, 1.0));
}

TEST_CASE("Test_Sepic_DivideByZero_Iout_Guard",
          "[converter-model][sepic-topology][unit]") {
    json j = make_sepic_json(5.0, 12.0, /*Iout*/ 0.5, 600e3);
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);

    auto ops = sepic.get_mutable_operating_points();
    ops[0].set_output_currents({0.0});
    sepic.set_operating_points(ops);
    REQUIRE_THROWS(sepic.generate_ngspice_circuit(/*L1*/ 22e-6));
}

// =====================================================================
//   Reference designs (analytical scalar checks — no SPICE)
// =====================================================================

TEST_CASE("Test_Sepic_TI_SNVA168E_Reference_Design",
          "[converter-model][sepic-topology][refdesign][analytical]") {
    // TI SNVA168E §6 (Designing A SEPIC Converter):
    //   Vin = 5 V, Vo = 12 V, Iout = 0.5 A, fsw = 600 kHz.
    // Use ideal η=1, Vd=0 so we hit the textbook D = 12/17 ≈ 0.7059.
    json j = make_sepic_json(/*Vin*/ 5.0, /*Vo*/ 12.0, /*Iout*/ 0.5,
                             /*fsw*/ 600e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);

    double D = sepic.get_last_duty_cycle();
    REQUIRE_THAT(D, WithinRel(12.0 / 17.0, 1e-9));
    REQUIRE_THAT(sepic.get_last_conversion_ratio(),       WithinRel(D / (1.0 - D), 1e-9));
    REQUIRE_THAT(sepic.get_last_coupling_cap_voltage(),   WithinRel(5.0, 1e-9));   // VCs = Vin
    REQUIRE_THAT(sepic.get_last_output_inductor_average(),WithinRel(0.5, 1e-9));   // = Iout
    REQUIRE_THAT(sepic.get_last_input_inductor_average(), WithinRel(0.5 * D / (1.0 - D), 1e-9));
    REQUIRE_THAT(sepic.get_last_switch_peak_voltage(),    WithinRel(17.0, 1e-9));  // Vin + Vo
    REQUIRE_THAT(sepic.get_last_diode_peak_reverse_voltage(), WithinRel(17.0, 1e-9));
}

TEST_CASE("Test_Sepic_LTC1871_Reference_Design",
          "[converter-model][sepic-topology][refdesign][analytical]") {
    // ADI LTC1871-style SEPIC: Vin=3.3V, Vo=5V, Iout=1A, fsw=250kHz, η≈0.85.
    json j = make_sepic_json(/*Vin*/ 3.3, /*Vo*/ 5.0, /*Iout*/ 1.0,
                             /*fsw*/ 250e3, /*Vd*/ 0.0, /*η*/ 0.85);
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 10e-6);

    double D = sepic.get_last_duty_cycle();
    // D = 5/(3.3·0.85+5) = 5/7.805 ≈ 0.6406
    REQUIRE(D > 0.55);
    REQUIRE(D < 0.75);
    REQUIRE_THAT(sepic.get_last_coupling_cap_voltage(), WithinRel(3.3, 1e-9));
    REQUIRE_THAT(sepic.get_last_switch_peak_voltage(),  WithinRel(8.3, 1e-9));
    REQUIRE(sepic.get_last_conversion_ratio() > 0.0);   // non-inverting
}

// =====================================================================
//   Per-OP diagnostics populated
// =====================================================================

TEST_CASE("Test_Sepic_Diagnostics_Populated",
          "[converter-model][sepic-topology][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);

    REQUIRE(sepic.get_last_duty_cycle() > 0.0);
    REQUIRE(sepic.get_last_duty_cycle() < 1.0);
    REQUIRE(sepic.get_last_conversion_ratio() > 0.0);   // non-inverting
    REQUIRE(sepic.get_last_coupling_cap_voltage() > 0.0);
    REQUIRE(sepic.get_last_input_inductor_average() > 0.0);
    REQUIRE(sepic.get_last_output_inductor_average() > 0.0);
    REQUIRE(sepic.get_last_input_inductor_ripple() > 0.0);
    REQUIRE(sepic.get_last_output_inductor_ripple() > 0.0);
    REQUIRE(sepic.get_last_switch_peak_voltage() > 0.0);
    REQUIRE(sepic.get_last_switch_peak_current() > 0.0);
    REQUIRE(sepic.get_last_diode_peak_reverse_voltage() > 0.0);
    REQUIRE(sepic.get_last_diode_peak_current() > 0.0);
    REQUIRE(sepic.get_last_coupling_cap_rms_current() > 0.0);
    REQUIRE(sepic.get_last_dcm_k() > 0.0);
    REQUIRE(sepic.get_last_dcm_kcrit() > 0.0);
    REQUIRE(sepic.get_last_sized_l2() > 0.0);
    REQUIRE(sepic.get_last_sized_cs() > 0.0);
    REQUIRE(sepic.get_last_sized_co() > 0.0);
}

// =====================================================================
//   Volt-second balance on L1
// =====================================================================

TEST_CASE("Test_Sepic_OperatingPoints_VoltSecondBalance",
          "[converter-model][sepic-topology][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Sepic sepic(j);
    auto ops = sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0].get_excitations_per_winding().size() == 1);

    auto& exc = ops[0].get_excitations_per_winding()[0];
    REQUIRE(exc.get_voltage().has_value());
    auto vWfm = exc.get_voltage().value().get_waveform().value();
    auto& vData = vWfm.get_data();
    auto vTime  = vWfm.get_time().value();

    double area = 0.0;
    for (size_t k = 1; k < vData.size(); ++k) {
        double dt = vTime[k] - vTime[k - 1];
        area += 0.5 * (vData[k] + vData[k - 1]) * dt;
    }
    REQUIRE(std::abs(area) < 1e-5);
}

// =====================================================================
//   SPICE netlist string contract
// =====================================================================

TEST_CASE("Test_Sepic_SPICE_Netlist_StringContract",
          "[converter-model][sepic-topology][spice-netlist]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);

    std::string netlist = sepic.generate_ngspice_circuit(/*L1*/ 22e-6);

    // GEAR integration + UIC for .ic
    REQUIRE(netlist.find("METHOD=GEAR") != std::string::npos);
    REQUIRE(netlist.find("UIC") != std::string::npos);

    // All four reactive elements present
    REQUIRE(netlist.find("L1 ")   != std::string::npos);
    REQUIRE(netlist.find("L2 ")   != std::string::npos);
    REQUIRE(netlist.find("Cs ")   != std::string::npos);
    REQUIRE(netlist.find("Cout ") != std::string::npos);

    // ≥4 .ic= occurrences (one per reactive element)
    size_t icCount = 0; size_t pos = 0;
    while ((pos = netlist.find(" ic=", pos)) != std::string::npos) { ++icCount; ++pos; }
    REQUIRE(icCount >= 4);

    // L2 BOTTOM TERMINAL AT GND — the SEPIC distinguishing feature.
    REQUIRE(netlist.find("L2 0 l2_top") != std::string::npos);

    // Switch (low-side) and snubbers
    REQUIRE(netlist.find("S1 node_A_int 0 pwm_ctrl 0 SW1") != std::string::npos);
    REQUIRE(netlist.find("Rsnub_s1") != std::string::npos);
    REQUIRE(netlist.find("Csnub_s1") != std::string::npos);
    REQUIRE(netlist.find("Rsnub_d1") != std::string::npos);

    // Diode oriented toward positive vout (non-inverting).
    REQUIRE(netlist.find("D1 rect_in vout DIDEAL") != std::string::npos);

    // Probes
    REQUIRE(netlist.find("Bvpri_diff") != std::string::npos);
    REQUIRE(netlist.find("Bsw") != std::string::npos);
    REQUIRE(netlist.find("Vin_sense") != std::string::npos);
    REQUIRE(netlist.find("Vl1_sense") != std::string::npos);
    REQUIRE(netlist.find("Vl2_sense") != std::string::npos);
    REQUIRE(netlist.find("Vrect_sense") != std::string::npos);
}

// =====================================================================
//   OpenMagnetics::AdvancedSepic round-trip
// =====================================================================

TEST_CASE("Test_Sepic_AdvancedSepic_Process_RoundTrip",
          "[converter-model][sepic-topology][advanced]") {
    const double L1 = 22e-6;
    json j = make_advanced_sepic_json(5.0, 12.0, 0.5, 600e3, L1);
    OpenMagnetics::AdvancedSepic asepic(j);
    auto inputs = asepic.process();
    REQUIRE(inputs.get_design_requirements().get_topology() == Topologies::SEPIC_CONVERTER);
    auto Lnom = inputs.get_design_requirements().get_magnetizing_inductance().get_nominal();
    REQUIRE(Lnom.has_value());
    REQUIRE_THAT(Lnom.value(), WithinRel(L1, 1e-6));
    REQUIRE(inputs.get_operating_points().size() == 1);
    REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 1);
}

// =====================================================================
//   get_extra_components_inputs (P7) — L2 / Cs / Co
// =====================================================================

TEST_CASE("Test_Sepic_ExtraComponents_L2_Cs_Co",
          "[converter-model][sepic-topology][extra-components]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);

    auto extras = sepic.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    REQUIRE(extras.size() == 3);

    // Index 0: L2 inductor
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[0]));
    {
        const auto& l2 = std::get<OpenMagnetics::Inputs>(extras[0]);
        REQUIRE(l2.get_design_requirements().get_topology() == Topologies::SEPIC_CONVERTER);
        REQUIRE(l2.get_design_requirements().get_name().value() == "outputInductor");
        auto Lnom = l2.get_design_requirements().get_magnetizing_inductance().get_nominal();
        REQUIRE(Lnom.has_value());
        REQUIRE_THAT(Lnom.value(), WithinRel(sepic.get_last_sized_l2(), 1e-9));
        REQUIRE(l2.get_operating_points().size() == 1);
        REQUIRE(l2.get_operating_points()[0].get_excitations_per_winding().size() == 1);
    }

    // Index 1: Cs coupling cap (DC_LINK)
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[1]));
    {
        const auto& cs = std::get<CAS::Inputs>(extras[1]);
        REQUIRE(cs.get_design_requirements().get_name().value() == "couplingCapacitor");
        auto Cnom = cs.get_design_requirements().get_capacitance().get_nominal();
        REQUIRE(Cnom.has_value());
        REQUIRE_THAT(Cnom.value(), WithinRel(sepic.get_last_sized_cs(), 1e-9));
        REQUIRE(cs.get_design_requirements().get_role() == CAS::Application::DC_LINK);
        // Rated voltage > peak Cs voltage (with ripple, +20 % margin).
        double Vrated = cs.get_design_requirements().get_rated_voltage();
        REQUIRE(Vrated > sepic.get_last_coupling_cap_voltage());
        REQUIRE(cs.get_operating_points().size() == 1);
    }

    // Index 2: Co output cap (OUTPUT_FILTER)
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[2]));
    {
        const auto& co = std::get<CAS::Inputs>(extras[2]);
        REQUIRE(co.get_design_requirements().get_name().value() == "outputCapacitor");
        auto Cnom = co.get_design_requirements().get_capacitance().get_nominal();
        REQUIRE(Cnom.has_value());
        REQUIRE_THAT(Cnom.value(), WithinRel(sepic.get_last_sized_co(), 1e-9));
        REQUIRE(co.get_design_requirements().get_role() == CAS::Application::OUTPUT_FILTER);
        REQUIRE(co.get_operating_points().size() == 1);
    }
}

TEST_CASE("Test_Sepic_ExtraComponents_Throws_BeforeProcess",
          "[converter-model][sepic-topology][extra-components]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    OpenMagnetics::Sepic sepic(j);
    REQUIRE_THROWS(sepic.get_extra_components_inputs(ExtraComponentsMode::IDEAL));
}

TEST_CASE("Test_Sepic_ExtraComponents_REAL_Requires_Magnetic",
          "[converter-model][sepic-topology][extra-components]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);
    REQUIRE_THROWS(sepic.get_extra_components_inputs(ExtraComponentsMode::REAL,
                                                     /*magnetic*/ std::nullopt));
}

// =====================================================================
//   CCM/DCM boundary detection
// =====================================================================

TEST_CASE("Test_Sepic_CCM_DCM_Boundary_Flips_With_Load",
          "[converter-model][sepic-topology][unit]") {
    auto run_at = [](double Iout) -> bool {
        json j = make_sepic_json(/*Vin*/ 5.0, /*Vo*/ 12.0, Iout,
                                 /*fsw*/ 600e3, /*Vd*/ 0.0, /*η*/ 1.0);
        OpenMagnetics::Sepic sepic(j);
        sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);
        return sepic.get_last_is_ccm();
    };
    // Heavy vs light load — at least one must flip CCM/DCM regime.
    bool isCcmHeavy = run_at(/*Iout*/ 5.0);
    bool isCcmLight = run_at(/*Iout*/ 0.001);
    REQUIRE(isCcmHeavy != isCcmLight);
}

// =====================================================================
//   V2 coupled-inductor: DR + netlist contract
// =====================================================================

TEST_CASE("Test_Sepic_V2_CoupledInductor_DR_Has_Two_Windings",
          "[converter-model][sepic-topology][v2-coupled][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    j["coupledInductor"] = true;
    OpenMagnetics::Sepic sepic(j);

    auto dr = sepic.process_design_requirements();
    REQUIRE(dr.get_topology() == Topologies::SEPIC_CONVERTER);
    REQUIRE(dr.get_isolation_sides().has_value());
    REQUIRE(dr.get_isolation_sides().value().size() == 2);
    REQUIRE(dr.get_turns_ratios().size() == 1);
    auto Nratio = dr.get_turns_ratios()[0].get_nominal();
    REQUIRE(Nratio.has_value());
    REQUIRE_THAT(Nratio.value(), WithinRel(1.0, 1e-12));
}

TEST_CASE("Test_Sepic_V2_CoupledInductor_Netlist_Has_K_L1L2",
          "[converter-model][sepic-topology][v2-coupled][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    j["coupledInductor"]     = true;
    j["couplingCoefficient"] = 0.98;
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);

    std::string netlist = sepic.generate_ngspice_circuit(/*L1*/ 22e-6);

    REQUIRE(netlist.find("K_L1L2 L1 L2") != std::string::npos);
    REQUIRE(netlist.find("0.98") != std::string::npos);
    REQUIRE(netlist.find("V2 coupled-inductor") != std::string::npos);
}

TEST_CASE("Test_Sepic_V2_CoupledInductor_Default_K_Is_Tight",
          "[converter-model][sepic-topology][v2-coupled][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    j["coupledInductor"] = true;
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);
    std::string netlist = sepic.generate_ngspice_circuit(/*L1*/ 22e-6);
    REQUIRE(netlist.find("K_L1L2 L1 L2 0.999") != std::string::npos);
}

TEST_CASE("Test_Sepic_V2_CoupledInductor_Throws_On_Bad_K",
          "[converter-model][sepic-topology][v2-coupled][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    j["coupledInductor"]     = true;
    j["couplingCoefficient"] = 1.5;
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);
    REQUIRE_THROWS(sepic.generate_ngspice_circuit(/*L1*/ 22e-6));
}

TEST_CASE("Test_Sepic_V1_Default_DR_Has_Single_Winding",
          "[converter-model][sepic-topology][v1-default][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    OpenMagnetics::Sepic sepic(j);

    auto dr = sepic.process_design_requirements();
    REQUIRE(dr.get_isolation_sides().has_value());
    REQUIRE(dr.get_isolation_sides().value().size() == 1);
    REQUIRE(dr.get_turns_ratios().empty());
}

TEST_CASE("Test_Sepic_V1_Default_Netlist_Has_No_Coupling",
          "[converter-model][sepic-topology][v1-default][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);
    std::string netlist = sepic.generate_ngspice_circuit(/*L1*/ 22e-6);

    REQUIRE(netlist.find("K_L1L2") == std::string::npos);
    REQUIRE(netlist.find("V1 uncoupled") != std::string::npos);
    REQUIRE(netlist.find(".model DIDEAL") != std::string::npos);
    REQUIRE(netlist.find("D1 rect_in vout DIDEAL") != std::string::npos);
    REQUIRE(netlist.find(".model SW2") == std::string::npos);
}

// =====================================================================
//   Synchronous rectifier: D1 → S2 with complementary PWM
// =====================================================================

TEST_CASE("Test_Sepic_Synchronous_Netlist_Replaces_Diode",
          "[converter-model][sepic-topology][synchronous][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    j["synchronousRectifier"] = true;
    OpenMagnetics::Sepic sepic(j);
    sepic.process_operating_points(std::vector<double>{}, /*L1*/ 22e-6);

    std::string netlist = sepic.generate_ngspice_circuit(/*L1*/ 22e-6);

    // No D1 model / element
    REQUIRE(netlist.find("D1 rect_in vout DIDEAL") == std::string::npos);
    REQUIRE(netlist.find(".model DIDEAL") == std::string::npos);

    // S2 + complementary PWM source
    REQUIRE(netlist.find(".model SW2") != std::string::npos);
    REQUIRE(netlist.find("S2 rect_in vout pwm_ctrl_inv 0 SW2") != std::string::npos);
    REQUIRE(netlist.find("Vpwm_inv pwm_ctrl_inv") != std::string::npos);
    REQUIRE(netlist.find("/synchronous") != std::string::npos);

    // Primary PWM unchanged
    REQUIRE(netlist.find("S1 node_A_int 0 pwm_ctrl 0 SW1") != std::string::npos);
}

// =====================================================================
//   process_design_requirements: missing both ripple/Imax → throws
// =====================================================================

TEST_CASE("Test_Sepic_DR_Throws_If_No_Ripple_Or_Imax",
          "[converter-model][sepic-topology][unit]") {
    json j = make_sepic_json(5.0, 12.0, 0.5, 600e3);
    j.erase("currentRippleRatio");
    j.erase("maximumSwitchCurrent");
    OpenMagnetics::Sepic sepic(j);
    REQUIRE_THROWS(sepic.process_design_requirements());
}
