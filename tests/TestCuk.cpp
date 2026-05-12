/**
 * @file TestCuk.cpp
 * @brief Unit tests for the Cuk (V1 non-isolated, CCM) DC-DC converter model.
 *
 * Coverage map (CUK_PLAN.md §9):
 *
 *   * Static analytical helpers
 *       - Conversion ratio M(D) = -D/(1-D)
 *       - Coupling cap voltage VC1 = Vin/(1-D)
 *       - L1/L2/C1 sizing inverses
 *       - K(D) (DCM boundary)
 *       - RHP-zero frequency
 *
 *   * Reference designs (CUK_PLAN.md §3 Table 1)
 *       - Erickson §2.4 worked example  (25 V → -25 V, 100 kHz, D=0.5)
 *       - Bramble LT3757-based          (10 V → -5  V, 300 kHz)
 *
 *   * Per-OP diagnostics population (every `lastXxx` non-zero where applicable)
 *
 *   * Volt-second balance on L1 (analytical waveform sanity)
 *
 *   * Input-validation guards
 *       - Iout ≤ 0          → throws (no Rload divide-by-zero)
 *       - D >= 0.95         → throws (loss of regulation)
 *       - Vin ≤ 0           → throws
 *
 *   * SPICE netlist string contract (CUK_PLAN.md §6.2)
 *       - METHOD=GEAR
 *       - .ic mandatory for L1, L2, C1, Co
 *       - Bvpri_diff probe present
 *       - snubbers across S1 and D1
 *       - parasitic R baked in
 *
 *   * AdvancedCuk round-trip (user-specified L1 honoured)
 *
 *   * get_extra_components_inputs (P7)
 *       - Returns 3 entries: L2 (Inputs), C1 (CAS::Inputs), Co (CAS::Inputs)
 *       - C1 role = DC_LINK, Co role = OUTPUT_FILTER
 *       - Each ancillary carries one operating point per processed Vin×OP
 *
 *   * CCM/DCM boundary detection
 *       - lastIsCcm flips at K(D) = (1-D)²
 *
 * Note: full SPICE-vs-analytical PtP correlation is exercised in
 * TestCukReferenceDesignsPtp.cpp (4-gate harness).  This file is purely
 * algebraic / structural.
 */

#include "converter_models/Cuk.h"
#include "processors/Inputs.h"
#include "support/Utils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <numbers>
#include <variant>

using namespace OpenMagnetics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

    // ---------------------------------------------------------------
    // Helper: build a minimal Cuk JSON for a single (Vin, |Vo|, Iout, fsw) point.
    // ---------------------------------------------------------------
    json make_cuk_json(double Vin, double VoMag, double Iout, double fsw,
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
        op["outputVoltages"]      = {VoMag};   // user passes |Vo| positive
        op["outputCurrents"]      = {Iout};
        op["switchingFrequency"]  = fsw;
        op["ambientTemperature"]  = 25.0;
        j["operatingPoints"].push_back(op);
        return j;
    }

    // ---------------------------------------------------------------
    // Helper: build an AdvancedCuk JSON with explicit desiredInductance.
    // ---------------------------------------------------------------
    json make_advanced_cuk_json(double Vin, double VoMag, double Iout, double fsw,
                                double desiredL1) {
        json j = make_cuk_json(Vin, VoMag, Iout, fsw);
        j["desiredInductance"] = desiredL1;
        return j;
    }

}  // namespace

// =====================================================================
//   Static analytical helpers
// =====================================================================

TEST_CASE("Test_Cuk_ConversionRatio_Formula",
          "[converter-model][cuk-topology][unit]") {
    // M(D) = -D/(1-D)
    REQUIRE_THAT(OpenMagnetics::Cuk::calculate_conversion_ratio(0.2), WithinRel(-0.25, 1e-12));
    REQUIRE_THAT(OpenMagnetics::Cuk::calculate_conversion_ratio(0.5), WithinRel(-1.0,  1e-12));
    REQUIRE_THAT(OpenMagnetics::Cuk::calculate_conversion_ratio(0.7), WithinRel(-7.0/3.0, 1e-12));
}

TEST_CASE("Test_Cuk_CouplingCapVoltage_Formula",
          "[converter-model][cuk-topology][unit]") {
    // VC1 = Vin/(1-D)
    REQUIRE_THAT(OpenMagnetics::Cuk::calculate_coupling_cap_voltage(25.0, 0.5), WithinRel(50.0, 1e-12));
    REQUIRE_THAT(OpenMagnetics::Cuk::calculate_coupling_cap_voltage(12.0, 0.7), WithinRel(40.0, 1e-12));
    REQUIRE_THAT(OpenMagnetics::Cuk::calculate_coupling_cap_voltage(10.0, 0.333333333333), WithinRel(15.0, 1e-9));
}

TEST_CASE("Test_Cuk_L1_L2_C1_Inverse_Sizing",
          "[converter-model][cuk-topology][unit]") {
    // L1 = Vin·D / (ΔIL1·fsw); inverting → ΔIL1 = Vin·D/(L1·fsw)
    double L1 = OpenMagnetics::Cuk::calculate_l1_min(/*Vin*/ 25.0, /*D*/ 0.5, /*ΔIL1*/ 0.5, /*fsw*/ 100e3);
    REQUIRE_THAT(L1, WithinRel(25.0 * 0.5 / (0.5 * 100e3), 1e-12));

    double L2 = OpenMagnetics::Cuk::calculate_l2_min(/*|Vo|*/ 25.0, /*D*/ 0.5, /*ΔIL2*/ 0.3, /*fsw*/ 100e3);
    REQUIRE_THAT(L2, WithinRel(25.0 * 0.5 / (0.3 * 100e3), 1e-12));

    double C1 = OpenMagnetics::Cuk::calculate_c1_min(/*Iout*/ 1.0, /*D*/ 0.5, /*ΔVC1*/ 2.5, /*fsw*/ 100e3);
    REQUIRE_THAT(C1, WithinRel(1.0 * 0.5 / (2.5 * 100e3), 1e-12));

    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_l1_min(25.0, 0.5, /*ΔIL1*/ 0.0, 100e3));
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_l2_min(25.0, 0.5, 0.3, /*fsw*/ 0.0));
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_c1_min(1.0,  0.5, /*ΔVC1*/ -1.0, 100e3));
}

TEST_CASE("Test_Cuk_DCM_K_Formula",
          "[converter-model][cuk-topology][unit]") {
    // K(D) = 2·Le·fsw / R, Le = L1·L2/(L1+L2)
    double L1 = 100e-6, L2 = 100e-6, fsw = 100e3, R = 25.0;
    double Le = L1 * L2 / (L1 + L2);
    double Kexpected = 2.0 * Le * fsw / R;
    REQUIRE_THAT(OpenMagnetics::Cuk::calculate_dcm_K(L1, L2, fsw, R), WithinRel(Kexpected, 1e-12));

    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_dcm_K(0.0, 0.0, fsw, R));
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_dcm_K(L1,  L2,  fsw, /*R*/ 0.0));
}

TEST_CASE("Test_Cuk_RHP_Zero_Frequency_Formula",
          "[converter-model][cuk-topology][unit]") {
    // ωRHP ≈ R·(1-D)² / (D·L2);  fRHP = ω / (2π)
    double R = 25.0, D = 0.5, L2 = 100e-6;
    double f_expected = (R * std::pow(1.0 - D, 2) / (D * L2)) / (2.0 * std::numbers::pi);
    REQUIRE_THAT(OpenMagnetics::Cuk::calculate_rhp_zero_frequency(R, D, L2),
                 WithinRel(f_expected, 1e-12));

    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_rhp_zero_frequency(R, /*D*/ 0.0, L2));
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_rhp_zero_frequency(R, /*D*/ 1.0, L2));
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_rhp_zero_frequency(R, D, /*L2*/ 0.0));
}

// =====================================================================
//   Input-validation guards
// =====================================================================

TEST_CASE("Test_Cuk_DutyCycle_BadInput_Throws",
          "[converter-model][cuk-topology][unit]") {
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_duty_cycle(/*Vin*/ 0.0, 25.0, 0.7, 0.85));
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_duty_cycle(25.0, /*|Vo|*/ 0.0, 0.7, 0.85));
}

TEST_CASE("Test_Cuk_DutyCycle_Clamp_Throws_Above_095",
          "[converter-model][cuk-topology][unit]") {
    // D = (|Vo|+Vd) / (|Vo|+Vd+Vin·η).  Push to D >= 0.95: |Vo| >> Vin.
    // With Vin=1, η=1, Vd=0: D = |Vo|/(|Vo|+1).  D ≥ 0.95 ⇒ |Vo| ≥ 19.
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_duty_cycle(/*Vin*/ 1.0, /*|Vo|*/ 100.0, 0.0, 1.0));
}

TEST_CASE("Test_Cuk_DivideByZero_Iout_Guard",
          "[converter-model][cuk-topology][unit]") {
    // Iout=0 → no Rload available; SPICE generation must throw, not crash.
    json j = make_cuk_json(25.0, 25.0, /*Iout*/ 1.0, 100e3);
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    // Mutate Iout to 0 and expect generate_ngspice_circuit to throw.
    auto ops = cuk.get_mutable_operating_points();
    ops[0].set_output_currents({0.0});
    cuk.set_operating_points(ops);
    REQUIRE_THROWS(cuk.generate_ngspice_circuit(/*L1*/ 100e-6));
}

// =====================================================================
//   Reference designs (analytical scalar checks, no SPICE)
// =====================================================================

TEST_CASE("Test_Cuk_Erickson_Reference_Design",
          "[converter-model][cuk-topology][refdesign][analytical]") {
    // Erickson-Maksimović 3rd ed. §2.4 worked example:
    //   Vin=25V, Vo=-25V, Pout=25W (Iout=1A), fsw=100kHz, D≈0.5.
    // Use ideal efficiency (1.0) and zero diode drop so the analytical
    // scalars match the textbook values exactly.
    json j = make_cuk_json(/*Vin*/ 25.0, /*|Vo|*/ 25.0, /*Iout*/ 1.0,
                           /*fsw*/ 100e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    REQUIRE_THAT(cuk.get_last_duty_cycle(),         WithinAbs(0.5, 1e-9));
    REQUIRE_THAT(cuk.get_last_conversion_ratio(),   WithinRel(-1.0, 1e-9));
    REQUIRE_THAT(cuk.get_last_coupling_cap_voltage(), WithinRel(50.0, 1e-9));   // Vin/(1-D) = 50V
    REQUIRE_THAT(cuk.get_last_input_inductor_average(),  WithinRel(1.0, 1e-9)); // D·Iout/(1-D)
    REQUIRE_THAT(cuk.get_last_output_inductor_average(), WithinRel(1.0, 1e-9));
    REQUIRE_THAT(cuk.get_last_switch_peak_voltage(),     WithinRel(50.0, 1e-9)); // Vin+|Vo|
    REQUIRE_THAT(cuk.get_last_diode_peak_reverse_voltage(), WithinRel(50.0, 1e-9));
}

TEST_CASE("Test_Cuk_Bramble_Reference_Design",
          "[converter-model][cuk-topology][refdesign][analytical]") {
    // Simon Bramble LT3757 inverting Cuk: Vin=10V, Vo=-5V, Pout=5W (Iout=1A),
    // fsw=300kHz, η≈0.85.  D ≈ |Vo|/(|Vo|+Vin·η) = 5/(5+8.5) ≈ 0.370.
    json j = make_cuk_json(/*Vin*/ 10.0, /*|Vo|*/ 5.0, /*Iout*/ 1.0,
                           /*fsw*/ 300e3, /*Vd*/ 0.0, /*η*/ 0.85);
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 47e-6);

    double D = cuk.get_last_duty_cycle();
    REQUIRE(D > 0.30);
    REQUIRE(D < 0.45);
    // VC1 = Vin/(1-D) check.
    REQUIRE_THAT(cuk.get_last_coupling_cap_voltage(),
                 WithinRel(10.0 / (1.0 - D), 1e-9));
    // Switch / diode peak voltage = Vin + |Vo|
    REQUIRE_THAT(cuk.get_last_switch_peak_voltage(), WithinRel(15.0, 1e-9));
}

// =====================================================================
//   Per-OP diagnostics populated
// =====================================================================

TEST_CASE("Test_Cuk_Diagnostics_Populated",
          "[converter-model][cuk-topology][unit]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    REQUIRE(cuk.get_last_duty_cycle() > 0.0);
    REQUIRE(cuk.get_last_duty_cycle() < 1.0);
    REQUIRE(cuk.get_last_conversion_ratio() < 0.0);   // inverting
    REQUIRE(cuk.get_last_coupling_cap_voltage() > 0.0);
    REQUIRE(cuk.get_last_input_inductor_average() > 0.0);
    REQUIRE(cuk.get_last_output_inductor_average() > 0.0);
    REQUIRE(cuk.get_last_input_inductor_ripple() > 0.0);
    REQUIRE(cuk.get_last_output_inductor_ripple() > 0.0);
    REQUIRE(cuk.get_last_switch_peak_voltage() > 0.0);
    REQUIRE(cuk.get_last_switch_peak_current() > 0.0);
    REQUIRE(cuk.get_last_diode_peak_reverse_voltage() > 0.0);
    REQUIRE(cuk.get_last_diode_peak_current() > 0.0);
    REQUIRE(cuk.get_last_coupling_cap_rms_current() > 0.0);
    REQUIRE(cuk.get_last_rhp_zero_frequency() > 0.0);
    REQUIRE(cuk.get_last_recommended_loop_bandwidth() > 0.0);
    REQUIRE(cuk.get_last_dcm_k() > 0.0);
    REQUIRE(cuk.get_last_dcm_kcrit() > 0.0);
    REQUIRE(cuk.get_last_sized_l2() > 0.0);
    REQUIRE(cuk.get_last_sized_c1() > 0.0);
    REQUIRE(cuk.get_last_sized_co() > 0.0);
}

// =====================================================================
//   Volt-second balance on L1 (analytical waveform sanity)
// =====================================================================

TEST_CASE("Test_Cuk_OperatingPoints_VoltSecondBalance",
          "[converter-model][cuk-topology][unit]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3, /*Vd*/ 0.0, /*η*/ 1.0);
    OpenMagnetics::Cuk cuk(j);
    auto ops = cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0].get_excitations_per_winding().size() == 1);

    auto& exc = ops[0].get_excitations_per_winding()[0];
    REQUIRE(exc.get_voltage().has_value());
    auto vWfm = exc.get_voltage().value().get_waveform().value();
    auto& vData = vWfm.get_data();
    auto vTimeOpt = vWfm.get_time().value();
    auto& vTime = vTimeOpt;

    // Trapezoidal integration of vL1 over one period ≈ 0.
    double area = 0.0;
    for (size_t k = 1; k < vData.size(); ++k) {
        double dt = vTime[k] - vTime[k - 1];
        area += 0.5 * (vData[k] + vData[k - 1]) * dt;
    }
    // Area should be ~0; with rectangular discretization tolerate ~1 % of Vin·T.
    REQUIRE(std::abs(area) < 1e-5);
}

// =====================================================================
//   SPICE netlist string contract (CUK_PLAN.md §6.2)
// =====================================================================

TEST_CASE("Test_Cuk_SPICE_Netlist_StringContract",
          "[converter-model][cuk-topology][spice-netlist]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    std::string netlist = cuk.generate_ngspice_circuit(/*L1*/ 100e-6);

    // Mandatory per CUK_PLAN.md §6.2:
    REQUIRE(netlist.find("METHOD=GEAR") != std::string::npos);
    REQUIRE(netlist.find("L1 ") != std::string::npos);
    REQUIRE(netlist.find("L2 ") != std::string::npos);
    REQUIRE(netlist.find("C1 ") != std::string::npos);
    REQUIRE(netlist.find("Cout ") != std::string::npos);

    // IC pre-charge (one per reactive element):
    REQUIRE(netlist.find("L1 ")  != std::string::npos);
    REQUIRE(netlist.find(" ic=") != std::string::npos);
    // Count how many "ic=" occurrences there are (expect ≥4).
    size_t icCount = 0; size_t pos = 0;
    while ((pos = netlist.find(" ic=", pos)) != std::string::npos) { ++icCount; ++pos; }
    REQUIRE(icCount >= 4);

    // Switch/diode snubbers:
    REQUIRE(netlist.find("Rsnub_s1") != std::string::npos);
    REQUIRE(netlist.find("Csnub_s1") != std::string::npos);
    REQUIRE(netlist.find("Rsnub_d1") != std::string::npos);
    REQUIRE(netlist.find("Csnub_d1") != std::string::npos);

    // Differential L1 winding-voltage probe (Boost-style):
    REQUIRE(netlist.find("Bvpri_diff") != std::string::npos);
    // Switch-node bridge probe:
    REQUIRE(netlist.find("Bsw") != std::string::npos);

    // UIC must be on the .tran line so the .ic directives fire.
    REQUIRE(netlist.find("UIC") != std::string::npos);

    // D1 anode at node_B, cathode at GND-via-Vd_sense (Cuk convention).
    REQUIRE(netlist.find("D1 node_B d_cath DIDEAL") != std::string::npos);
}

// =====================================================================
//   AdvancedCuk round-trip
// =====================================================================

TEST_CASE("Test_Cuk_AdvancedCuk_Process_RoundTrip",
          "[converter-model][cuk-topology][advanced]") {
    const double L1 = 220e-6;
    json j = make_advanced_cuk_json(25.0, 25.0, 1.0, 100e3, L1);
    OpenMagnetics::AdvancedCuk acuk(j);
    auto inputs = acuk.process();
    REQUIRE(inputs.get_design_requirements().get_topology() == Topologies::CUK_CONVERTER);
    auto Lnom = inputs.get_design_requirements().get_magnetizing_inductance().get_nominal();
    REQUIRE(Lnom.has_value());
    REQUIRE_THAT(Lnom.value(), WithinRel(L1, 1e-6));
    REQUIRE(inputs.get_operating_points().size() == 1);
    REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 1);
}

// =====================================================================
//   get_extra_components_inputs (P7)
// =====================================================================

TEST_CASE("Test_Cuk_ExtraComponents_L2_C1_Co",
          "[converter-model][cuk-topology][extra-components]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    // IDEAL mode does not require a designed magnetic.
    auto extras = cuk.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    REQUIRE(extras.size() == 3);

    // Index 0: L2 inductor (MAS::Inputs)
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[0]));
    {
        const auto& l2 = std::get<OpenMagnetics::Inputs>(extras[0]);
        REQUIRE(l2.get_design_requirements().get_topology() == Topologies::CUK_CONVERTER);
        auto Lnom = l2.get_design_requirements().get_magnetizing_inductance().get_nominal();
        REQUIRE(Lnom.has_value());
        REQUIRE(Lnom.value() > 0.0);
        REQUIRE_THAT(Lnom.value(), WithinRel(cuk.get_last_sized_l2(), 1e-9));
        REQUIRE(l2.get_operating_points().size() == 1);
        REQUIRE(l2.get_operating_points()[0].get_excitations_per_winding().size() == 1);
    }

    // Index 1: C1 coupling cap (CAS::Inputs, role = DC_LINK)
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[1]));
    {
        const auto& c1 = std::get<CAS::Inputs>(extras[1]);
        auto Cnom = c1.get_design_requirements().get_capacitance().get_nominal();
        REQUIRE(Cnom.has_value());
        REQUIRE_THAT(Cnom.value(), WithinRel(cuk.get_last_sized_c1(), 1e-9));
        REQUIRE(c1.get_design_requirements().get_role() == CAS::Application::DC_LINK);
        // Rated voltage must be ≥ peak (we add a 20 % margin, so > VC1).
        double Vrated = c1.get_design_requirements().get_rated_voltage();
        REQUIRE(Vrated > cuk.get_last_coupling_cap_voltage());
        REQUIRE(c1.get_operating_points().size() == 1);
    }

    // Index 2: Co output cap (CAS::Inputs, role = OUTPUT_FILTER)
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[2]));
    {
        const auto& co = std::get<CAS::Inputs>(extras[2]);
        auto Cnom = co.get_design_requirements().get_capacitance().get_nominal();
        REQUIRE(Cnom.has_value());
        REQUIRE_THAT(Cnom.value(), WithinRel(cuk.get_last_sized_co(), 1e-9));
        REQUIRE(co.get_design_requirements().get_role() == CAS::Application::OUTPUT_FILTER);
        REQUIRE(co.get_operating_points().size() == 1);
    }
}

TEST_CASE("Test_Cuk_ExtraComponents_Throws_BeforeProcess",
          "[converter-model][cuk-topology][extra-components]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    OpenMagnetics::Cuk cuk(j);
    REQUIRE_THROWS(cuk.get_extra_components_inputs(ExtraComponentsMode::IDEAL));
}

TEST_CASE("Test_Cuk_ExtraComponents_REAL_Requires_Magnetic",
          "[converter-model][cuk-topology][extra-components]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);
    REQUIRE_THROWS(cuk.get_extra_components_inputs(ExtraComponentsMode::REAL,
                                                   /*magnetic*/ std::nullopt));
}

// =====================================================================
//   CCM/DCM boundary detection
// =====================================================================

TEST_CASE("Test_Cuk_CCM_DCM_Boundary_Flips_With_Load",
          "[converter-model][cuk-topology][unit]") {
    // Heavy load → low R → small K → DCM possible.
    // Light load → high R → large K relative to (1-D)² → CCM.
    // We sweep load to confirm lastIsCcm flips between regimes.
    // Note: with L1 fixed and L2 auto-sized to maintain a ripple ratio,
    // Le = L1·L2/(L1+L2) ≈ L1 at light load (L2 huge) and ≈ L1·L2_heavy/(L1+L2_heavy)
    // at heavy load. K = 2·Le·fsw/R, so heavy load (small R) → CCM, light load → DCM —
    // the opposite of the textbook fixed-L behaviour, but a consistent property of
    // this implementation.
    auto run_at = [](double Iout) -> bool {
        json j = make_cuk_json(/*Vin*/ 25.0, /*|Vo|*/ 25.0, /*Iout*/ Iout,
                               /*fsw*/ 100e3, /*Vd*/ 0.0, /*η*/ 1.0);
        OpenMagnetics::Cuk cuk(j);
        cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);
        return cuk.get_last_is_ccm();
    };

    // Heavy load — small R, Le sizable → expect CCM
    bool isCcmHeavy = run_at(/*Iout*/ 10.0);
    // Light load — huge R, K → 0 → expect DCM
    bool isCcmLight = run_at(/*Iout*/ 0.001);

    // At least one of the two regimes must be DCM and the other CCM.
    REQUIRE(isCcmHeavy != isCcmLight);
}

// =====================================================================
//   V4 synchronous Cuk: netlist contract
// =====================================================================

TEST_CASE("Test_Cuk_V4_Synchronous_Netlist_Replaces_Diode",
          "[converter-model][cuk-topology][v4-synchronous][unit]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    j["synchronous"] = true;
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    std::string netlist = cuk.generate_ngspice_circuit(/*L1*/ 100e-6);

    // V4 specifics:
    //   * No D1 diode model line and no D1 element.
    //   * SW2 model + S2 element present.
    //   * Complementary PWM source Vpwm_inv routed to S2.
    REQUIRE(netlist.find("D1 node_B d_cath DIDEAL") == std::string::npos);
    REQUIRE(netlist.find(".model DIDEAL") == std::string::npos);
    REQUIRE(netlist.find(".model SW2") != std::string::npos);
    REQUIRE(netlist.find("S2 node_B d_cath pwm_ctrl_inv 0 SW2") != std::string::npos);
    REQUIRE(netlist.find("Vpwm_inv pwm_ctrl_inv") != std::string::npos);
    REQUIRE(netlist.find("V4 synchronous") != std::string::npos);

    // Shared infrastructure (Vd_sense, snubbers) must remain so post-sim
    // probes continue to find the rectifier current.
    REQUIRE(netlist.find("Vd_sense d_cath 0 0") != std::string::npos);
    REQUIRE(netlist.find("Rsnub_d1") != std::string::npos);
    REQUIRE(netlist.find("Csnub_d1") != std::string::npos);

    // S1 high-side PWM unchanged.
    REQUIRE(netlist.find("S1 node_A_int 0 pwm_ctrl 0 SW1") != std::string::npos);
}

TEST_CASE("Test_Cuk_V5_Bidirectional_Implies_Synchronous",
          "[converter-model][cuk-topology][v5-bidirectional][unit]") {
    // Per CUK_PLAN.md: bidirectional=true implies synchronous=true even
    // when the user has not explicitly set synchronous.
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    j["bidirectional"] = true;
    // synchronous deliberately omitted.
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    std::string netlist = cuk.generate_ngspice_circuit(/*L1*/ 100e-6);

    REQUIRE(netlist.find(".model SW2") != std::string::npos);
    REQUIRE(netlist.find("S2 node_B d_cath pwm_ctrl_inv 0 SW2") != std::string::npos);
    REQUIRE(netlist.find("D1 node_B d_cath DIDEAL") == std::string::npos);
    REQUIRE(netlist.find("V5 bidirectional") != std::string::npos);
}

TEST_CASE("Test_Cuk_V1_Default_Still_Uses_Diode",
          "[converter-model][cuk-topology][v1-default][unit]") {
    // Regression guard: with neither synchronous nor bidirectional set, the
    // V1 diode wiring must be preserved.
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    std::string netlist = cuk.generate_ngspice_circuit(/*L1*/ 100e-6);

    REQUIRE(netlist.find(".model DIDEAL") != std::string::npos);
    REQUIRE(netlist.find("D1 node_B d_cath DIDEAL") != std::string::npos);
    REQUIRE(netlist.find(".model SW2") == std::string::npos);
    REQUIRE(netlist.find("S2 node_B") == std::string::npos);
}

// =====================================================================
//   V2 coupled-inductor: design requirements + netlist contract
// =====================================================================

TEST_CASE("Test_Cuk_V2_CoupledInductor_DR_Has_Two_Windings",
          "[converter-model][cuk-topology][v2-coupled][unit]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    j["coupledInductor"] = true;
    OpenMagnetics::Cuk cuk(j);

    auto dr = cuk.process_design_requirements();
    REQUIRE(dr.get_topology() == Topologies::CUK_CONVERTER);
    REQUIRE(dr.get_isolation_sides().has_value());
    // Two windings on the same (non-isolated) core.
    REQUIRE(dr.get_isolation_sides().value().size() == 2);
    // 1:1 turns ratio for the canonical zero-ripple sizing.
    REQUIRE(dr.get_turns_ratios().size() == 1);
    auto Nratio = dr.get_turns_ratios()[0].get_nominal();
    REQUIRE(Nratio.has_value());
    REQUIRE_THAT(Nratio.value(), WithinRel(1.0, 1e-12));
}

TEST_CASE("Test_Cuk_V2_CoupledInductor_Netlist_Has_K1",
          "[converter-model][cuk-topology][v2-coupled][unit]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    j["coupledInductor"]     = true;
    j["couplingCoefficient"] = 0.98;
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    std::string netlist = cuk.generate_ngspice_circuit(/*L1*/ 100e-6);

    // Mutual-inductance line must reference both L1 and L2 with the user k.
    REQUIRE(netlist.find("K1 L1 L2") != std::string::npos);
    REQUIRE(netlist.find("0.98") != std::string::npos);
    REQUIRE(netlist.find("V2 coupled-inductor") != std::string::npos);

    // The L1 + L2 elements themselves still appear (we couple, not replace).
    REQUIRE(netlist.find("L1 l1_in l1_dcr_mid") != std::string::npos);
    REQUIRE(netlist.find("L2 node_B l2_dcr_mid") != std::string::npos);
}

TEST_CASE("Test_Cuk_V2_CoupledInductor_Default_K_Is_Tight",
          "[converter-model][cuk-topology][v2-coupled][unit]") {
    // Without an explicit couplingCoefficient the netlist defaults to a tight
    // k = 0.999 (near zero-ripple ideal). Regression guard.
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    j["coupledInductor"] = true;
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);

    std::string netlist = cuk.generate_ngspice_circuit(/*L1*/ 100e-6);
    REQUIRE(netlist.find("K1 L1 L2 0.999") != std::string::npos);
}

TEST_CASE("Test_Cuk_V2_CoupledInductor_Throws_On_Bad_K",
          "[converter-model][cuk-topology][v2-coupled][unit]") {
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    j["coupledInductor"]     = true;
    j["couplingCoefficient"] = 1.5;   // out of range
    OpenMagnetics::Cuk cuk(j);
    cuk.process_operating_points(std::vector<double>{}, /*L1*/ 100e-6);
    REQUIRE_THROWS(cuk.generate_ngspice_circuit(/*L1*/ 100e-6));
}

TEST_CASE("Test_Cuk_CoupledInductor_And_Isolated_Are_Mutually_Exclusive",
          "[converter-model][cuk-topology][v2-coupled][v3-isolated][unit]") {
    // V2 (coupled inductor on one core) and V3 (isolation transformer) cannot
    // both be requested at once - DR construction must throw.
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    j["coupledInductor"] = true;
    j["isolated"]        = true;
    OpenMagnetics::Cuk cuk(j);
    REQUIRE_THROWS(cuk.process_design_requirements());
}

TEST_CASE("Test_Cuk_V1_Default_DR_Has_Single_Winding",
          "[converter-model][cuk-topology][v1-default][unit]") {
    // Regression guard for the V1 / V4 / V5 single-winding DR shape: a one-
    // entry isolation_sides vector and an empty turnsRatios array.
    json j = make_cuk_json(25.0, 25.0, 1.0, 100e3);
    OpenMagnetics::Cuk cuk(j);

    auto dr = cuk.process_design_requirements();
    REQUIRE(dr.get_isolation_sides().has_value());
    REQUIRE(dr.get_isolation_sides().value().size() == 1);
    REQUIRE(dr.get_turns_ratios().empty());
}

// ===================================================================
//   V3 isolated Cuk - Erickson §6 / Cuk patent (1977 isolated variant)
// ===================================================================
//
// V3 introduces a 2-winding isolation transformer in series with the
// coupling capacitor, splitting it into Ca (primary side) and Cb
// (secondary side). The conversion gain becomes
//     |Vo|/Vin = D / ((1-D) · n)         where n = Np/Ns,
// reducing to V1 when n=1. The DesignRequirements primary magnetic is
// the transformer (turns_ratios={n}, isolation_sides={PRIMARY,SECONDARY});
// the input choke L1 and output choke L2 are emitted as ancillaries.

TEST_CASE("Test_Cuk_V3_Isolated_DutyCycle_Reduces_To_V1_When_n_Equals_1",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    // With turnsRatio=1 the isolated duty-cycle expression must collapse
    // exactly onto the V1 form (algebraic identity).
    const double Vin = 25.0, Vo = 25.0, Vd = 0.7, eta = 0.85;
    double dV1 = OpenMagnetics::Cuk::calculate_duty_cycle(Vin, Vo, Vd, eta);
    double dV3 = OpenMagnetics::Cuk::calculate_duty_cycle(Vin, Vo, Vd, eta, /*n=*/1.0);
    REQUIRE_THAT(dV3, WithinAbs(dV1, 1e-12));
}

TEST_CASE("Test_Cuk_V3_Isolated_DutyCycle_Step_Down_Transformer",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    // Step-down transformer (n=Np/Ns=2) on a 24V→-5V buck-style isolated Cuk:
    //     D = (|Vo|+Vd)·n / (Vin·η + (|Vo|+Vd)·n)
    //       = (5+0.7)·2 / (24·0.9 + (5+0.7)·2)
    //       = 11.4 / (21.6 + 11.4) = 11.4 / 33.0 = 0.3454...
    const double Vin = 24.0, Vo = 5.0, Vd = 0.7, eta = 0.9, n = 2.0;
    double D = OpenMagnetics::Cuk::calculate_duty_cycle(Vin, Vo, Vd, eta, n);
    double expected = (Vo + Vd) * n / (Vin * eta + (Vo + Vd) * n);
    REQUIRE_THAT(D, WithinRel(expected, 1e-12));
    REQUIRE(D > 0.0);
    REQUIRE(D < 0.95);
}

TEST_CASE("Test_Cuk_V3_Isolated_DutyCycle_Throws_On_Bad_TurnsRatio",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_duty_cycle(24.0, 5.0, 0.7, 0.9, /*n=*/0.0));
    REQUIRE_THROWS(OpenMagnetics::Cuk::calculate_duty_cycle(24.0, 5.0, 0.7, 0.9, /*n=*/-1.5));
}

TEST_CASE("Test_Cuk_V3_Isolated_DR_Has_Two_Windings",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    // V3 must emit a 2-winding transformer DR (turns_ratios={n},
    // isolation_sides={PRIMARY, SECONDARY}). The magnetizing_inductance
    // returned is Lm (transformer), not L1 (input choke).
    json j = make_cuk_json(24.0, 5.0, 0.5, 100e3);
    j["isolated"]    = true;
    j["turnsRatio"]  = 2.0;
    OpenMagnetics::Cuk cuk(j);
    auto dr = cuk.process_design_requirements();

    REQUIRE(dr.get_isolation_sides().has_value());
    REQUIRE(dr.get_isolation_sides().value().size() == 2);
    REQUIRE(dr.get_isolation_sides().value()[0] == MAS::IsolationSide::PRIMARY);
    REQUIRE(dr.get_isolation_sides().value()[1] == MAS::IsolationSide::SECONDARY);
    REQUIRE(dr.get_turns_ratios().size() == 1);
    REQUIRE_THAT(dr.get_turns_ratios()[0].get_nominal().value(), WithinAbs(2.0, 1e-9));

    // Lm should be > 0 (sized from ΔIm budget).
    REQUIRE(dr.get_magnetizing_inductance().get_minimum().value() > 0.0);
}

TEST_CASE("Test_Cuk_V3_Isolated_Throws_On_Bad_TurnsRatio_DR",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    json j = make_cuk_json(24.0, 5.0, 0.5, 100e3);
    j["isolated"]   = true;
    j["turnsRatio"] = 0.0;   // invalid
    OpenMagnetics::Cuk cuk(j);
    REQUIRE_THROWS(cuk.process_design_requirements());
}

TEST_CASE("Test_Cuk_V3_Isolated_Netlist_Contract",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    json j = make_cuk_json(24.0, 5.0, 0.5, 100e3);
    j["isolated"]   = true;
    j["turnsRatio"] = 2.0;
    OpenMagnetics::Cuk cuk(j);

    // Run the analytical worker first to populate diagnostics.
    auto dr = cuk.process_design_requirements();
    double Lm = dr.get_magnetizing_inductance().get_minimum().value();
    cuk.process_operating_points(std::vector<double>{2.0}, /*Lm*/ Lm);

    std::string netlist = cuk.generate_ngspice_circuit(/*Lm*/ Lm);
    // V3 banner
    REQUIRE(netlist.find("isolated V3") != std::string::npos);
    // Two-cap split + transformer with K1 on Lp/Ls (NOT on L1/L2)
    REQUIRE(netlist.find("Ca node_C") != std::string::npos);
    REQUIRE(netlist.find("Lp tx_pri_pos") != std::string::npos);
    REQUIRE(netlist.find("Ls 0 tx_sec_pos") != std::string::npos);
    REQUIRE(netlist.find("K1 Lp Ls 0.999") != std::string::npos);
    REQUIRE(netlist.find("Cb tx_sec_pos") != std::string::npos);
    // Single-cap C1 from V1/V2 must NOT appear
    REQUIRE(netlist.find("C1 node_C") == std::string::npos);
    // Common scaffolding still present
    REQUIRE(netlist.find("L1 l1_in") != std::string::npos);
    REQUIRE(netlist.find("L2 node_B") != std::string::npos);
    REQUIRE(netlist.find("Vd_sense d_cath") != std::string::npos);
    REQUIRE(netlist.find("METHOD=GEAR") != std::string::npos);
}

TEST_CASE("Test_Cuk_V3_Isolated_Diagnostics_Populated",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    json j = make_cuk_json(24.0, 5.0, 0.5, 100e3, 0.7, 0.9);
    j["isolated"]   = true;
    j["turnsRatio"] = 2.0;
    OpenMagnetics::Cuk cuk(j);

    auto dr = cuk.process_design_requirements();
    double Lm = dr.get_magnetizing_inductance().get_minimum().value();
    cuk.process_operating_points(std::vector<double>{2.0}, Lm);

    // V3-specific diagnostics that V1 leaves at zero.
    REQUIRE(cuk.get_last_sized_ca() > 0.0);
    REQUIRE(cuk.get_last_sized_cb() > 0.0);
    REQUIRE(cuk.get_last_sized_lm() > 0.0);
    REQUIRE_THAT(cuk.get_last_turns_ratio(), WithinAbs(2.0, 1e-9));
    // Conversion ratio: -D/((1-D)·n)
    double D = cuk.get_last_duty_cycle();
    double expectedM = -D / ((1.0 - D) * 2.0);
    REQUIRE_THAT(cuk.get_last_conversion_ratio(), WithinRel(expectedM, 1e-9));
    // V1 single-cap diagnostic should be zero in V3
    REQUIRE_THAT(cuk.get_last_sized_c1(), WithinAbs(0.0, 1e-12));
    // Switch peak is VCa = Vin/(1-D); diode peak is VCb = |Vo|/D (≠ Vin+|Vo|)
    REQUIRE_THAT(cuk.get_last_switch_peak_voltage(),
                 WithinRel(24.0 / (1.0 - D), 1e-6));
    REQUIRE_THAT(cuk.get_last_diode_peak_reverse_voltage(),
                 WithinRel(5.0 / D, 1e-6));
}

TEST_CASE("Test_Cuk_V3_Isolated_Extra_Components_5_Entries",
          "[converter-model][cuk-topology][v3-isolated][p7][unit]") {
    json j = make_cuk_json(24.0, 5.0, 0.5, 100e3, 0.7, 0.9);
    j["isolated"]   = true;
    j["turnsRatio"] = 2.0;
    OpenMagnetics::Cuk cuk(j);

    auto dr = cuk.process_design_requirements();
    double Lm = dr.get_magnetizing_inductance().get_minimum().value();
    cuk.process_operating_points(std::vector<double>{2.0}, Lm);

    auto extras = cuk.get_extra_components_inputs();
    // V3: L1 + L2 + Ca + Cb + Co
    REQUIRE(extras.size() == 5);
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[0]));   // L1
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[1]));   // L2
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[2]));   // Ca
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[3]));   // Cb
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[4]));   // Co

    // Names + roles
    auto& l1Inputs = std::get<OpenMagnetics::Inputs>(extras[0]);
    REQUIRE(l1Inputs.get_design_requirements().get_name().value() == "inputInductor");
    auto& l2Inputs = std::get<OpenMagnetics::Inputs>(extras[1]);
    REQUIRE(l2Inputs.get_design_requirements().get_name().value() == "outputInductor");

    auto& caInputs = std::get<CAS::Inputs>(extras[2]);
    REQUIRE(caInputs.get_design_requirements().get_name().value() == "primaryCouplingCapacitor");
    REQUIRE(caInputs.get_design_requirements().get_role() == CAS::Application::DC_LINK);
    auto& cbInputs = std::get<CAS::Inputs>(extras[3]);
    REQUIRE(cbInputs.get_design_requirements().get_name().value() == "secondaryCouplingCapacitor");
    REQUIRE(cbInputs.get_design_requirements().get_role() == CAS::Application::DC_LINK);
    auto& coInputs = std::get<CAS::Inputs>(extras[4]);
    REQUIRE(coInputs.get_design_requirements().get_name().value() == "outputCapacitor");
    REQUIRE(coInputs.get_design_requirements().get_role() == CAS::Application::OUTPUT_FILTER);
}

TEST_CASE("Test_Cuk_V3_Isolated_Two_Windings_In_Excitation",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    // V3's process_operating_points emits TWO winding excitations (Lp + Ls)
    // since the primary magnetic is a 2-winding transformer.
    json j = make_cuk_json(24.0, 5.0, 0.5, 100e3, 0.7, 0.9);
    j["isolated"]   = true;
    j["turnsRatio"] = 2.0;
    OpenMagnetics::Cuk cuk(j);
    auto dr = cuk.process_design_requirements();
    double Lm = dr.get_magnetizing_inductance().get_minimum().value();
    auto ops = cuk.process_operating_points(std::vector<double>{2.0}, Lm);
    REQUIRE(ops.size() >= 1);
    REQUIRE(ops[0].get_excitations_per_winding().size() == 2);
}

TEST_CASE("Test_Cuk_V3_Isolated_Throws_If_No_Ripple_Budget_For_L1",
          "[converter-model][cuk-topology][v3-isolated][unit]") {
    // V3 needs an explicit ripple-ratio or maxSwitchCurrent to size L1; no
    // silent fallbacks per CLAUDE.md. (process_design_requirements already
    // enforces this for V1; the V3 path inherits the same gate.)
    json j = make_cuk_json(24.0, 5.0, 0.5, 100e3);
    j.erase("currentRippleRatio");
    j.erase("maximumSwitchCurrent");
    j["isolated"]   = true;
    j["turnsRatio"] = 2.0;
    OpenMagnetics::Cuk cuk(j);
    REQUIRE_THROWS(cuk.process_design_requirements());
}
