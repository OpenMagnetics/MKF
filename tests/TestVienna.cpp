/**
 * @file TestVienna.cpp
 * @brief Phase 1+2 unit tests for the Vienna rectifier (3-phase 3-level PFC).
 *
 * Coverage (VIENNA_PLAN.md §A.10, pragmatic Phase 1+2 scope):
 *   - JSON parsing and run_checks for valid / malformed inputs
 *   - Static analytical helpers (Kolar 1994 closed forms)
 *   - process_design_requirements: M, I_pk, L derivation
 *   - process_operating_points: returns ONE OP with THREE windings
 *   - Switch voltage stress = Vdc/2 (the 3-level signature)
 *   - Throws on over-modulation, Vdc <= sqrt(2)*V_LL, unsupported flags
 *   - AdvancedVienna user override of L
 */
#include "converter_models/Vienna.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace OpenMagnetics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Build a minimal 3-phase Vienna JSON.  Defaults: 400 V_LL EU grid, 750 V dc bus,
// 70 kHz Fsw (typical SiC), 10 kW.
json make_vienna_json(double V_LL = 400.0,
                      double Vdc  = 750.0,
                      double Fsw  = 70e3,
                      double P    = 10e3,
                      double rippleRatio = 0.25) {
    json j;
    j["lineToLineVoltage"]   = {{"nominal", V_LL}};
    j["lineFrequency"]       = 50.0;
    j["outputDcVoltage"]     = Vdc;
    j["switchingFrequency"]  = Fsw;
    j["currentRippleRatio"]  = rippleRatio;
    j["efficiency"]          = 0.97;
    j["powerFactor"]         = 0.99;

    // OP: outputVoltages[0] = Vdc, outputCurrents[0] = total Idc (so V*I = P).
    json op;
    op["outputVoltages"]     = {Vdc};
    op["outputCurrents"]     = {P / Vdc};
    op["switchingFrequency"] = Fsw;
    op["ambientTemperature"] = 25.0;
    j["operatingPoints"]     = json::array({op});
    return j;
}

} // namespace


TEST_CASE("Vienna: JSON parses and run_checks accepts a valid 10 kW design", "[vienna-topology]") {
    auto j = make_vienna_json();
    Vienna v(j);
    REQUIRE(v.run_checks(false));
    REQUIRE_FALSE(v.is_bridge_topology());
    REQUIRE(v.topology_kind() == MAS::Topologies::VIENNA_RECTIFIER_CONVERTER);
}


TEST_CASE("Vienna: phase-peak / modulation-index helpers (Kolar 1994 conventions)", "[vienna-topology]") {
    // V_phase_peak = sqrt(2)*400/sqrt(3) = 326.6 V
    double Vpk = Vienna::compute_phase_peak_voltage(400.0);
    REQUIRE_THAT(Vpk, WithinRel(std::sqrt(2.0) * 400.0 / std::sqrt(3.0), 1e-12));

    // M = 326.6 / (750/2) = 0.871
    double M = Vienna::compute_modulation_index(Vpk, 750.0);
    REQUIRE_THAT(M, WithinRel(Vpk / 375.0, 1e-12));
    REQUIRE(M < 1.0);

    // Hartmann ETH 19755 (2011) switch RMS: I_pk * sqrt(1/4 - 2*M/(3*pi))
    double Irms = Vienna::compute_switch_rms(20.0, 0.871);
    double expected = 20.0 * std::sqrt(0.25 - 2.0 * 0.871 / (3.0 * M_PI));
    REQUIRE_THAT(Irms, WithinRel(expected, 1e-9));

    // Diode avg = I_pk / pi
    REQUIRE_THAT(Vienna::compute_diode_avg(31.4159265), WithinRel(10.0, 1e-6));
}


TEST_CASE("Vienna: process_design_requirements derives L, I_pk and switch stress", "[vienna-topology]") {
    auto j = make_vienna_json(/*V_LL*/400, /*Vdc*/750, /*Fsw*/70e3,
                              /*P*/10e3, /*ripple*/0.25);
    Vienna v(j);
    auto dr = v.process_design_requirements();

    REQUIRE(dr.get_topology() == Topologies::VIENNA_RECTIFIER_CONVERTER);
    REQUIRE(dr.get_turns_ratios().size() == 1);
    REQUIRE_THAT(dr.get_turns_ratios()[0].get_nominal().value(), WithinRel(1.0, 1e-9));

    // I_pk = sqrt(2)*P / (3 * V_phase_rms * eta * PF)
    //      = sqrt(2)*10000 / (3 * 230.94 * 0.97 * 0.99) = ~21.3 A
    double V_pk_phase  = std::sqrt(2.0) * 400.0 / std::sqrt(3.0);
    double V_rms_phase = V_pk_phase / std::sqrt(2.0);
    double I_pk_expected = std::sqrt(2.0) * 10e3 / (3.0 * V_rms_phase * 0.97 * 0.99);
    REQUIRE_THAT(v.get_computed_line_peak_current(),
                 WithinRel(I_pk_expected, 1e-6));

    // Switch stress = Vdc/2 = 375 V (the 3-level signature).
    REQUIRE_THAT(v.get_computed_switch_voltage_stress(), WithinRel(375.0, 1e-9));

    // M < 1 (no over-modulation).
    REQUIRE(v.get_computed_modulation_index() < 1.0);
    REQUIRE(v.get_computed_modulation_index() > 0.5);

    // L from peak-to-peak ripple target.
    //   L = V_phase_peak * (1 - M) / (rippleRatio * I_pk * Fsw)
    double M_expected   = V_pk_phase / (750.0 / 2.0);
    double L_expected   = V_pk_phase * (1.0 - M_expected) /
                          (0.25 * I_pk_expected * 70e3);
    REQUIRE_THAT(v.get_computed_boost_inductance(),
                 WithinRel(L_expected, 1e-6));
    REQUIRE(v.get_computed_boost_inductance() > 0);
}


TEST_CASE("Vienna: process_operating_points returns one OP with three phase windings", "[vienna-topology]") {
    auto j = make_vienna_json();
    Vienna v(j);
    auto ops = v.process_operating_points(std::vector<double>{}, 0.0);
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0].get_excitations_per_winding().size() == 3);

    // Diagnostics populated.
    REQUIRE(v.get_last_inductor_peak_current()      > 0);
    REQUIRE(v.get_last_inductor_ripple_peak_to_peak() > 0);
    REQUIRE(v.get_last_switch_voltage_stress()      > 0);
    REQUIRE(v.get_last_switch_rms_current()         > 0);
    REQUIRE(v.get_last_diode_avg_current()          > 0);

    // Switch V-stress = Vdc/2 (3-level signature).
    REQUIRE_THAT(v.get_last_switch_voltage_stress(), WithinRel(375.0, 1e-9));

    // Diode avg = I_pk / pi.
    double I_pk = v.get_computed_line_peak_current();
    REQUIRE_THAT(v.get_last_diode_avg_current(),
                 WithinRel(I_pk / M_PI, 1e-6));
}


TEST_CASE("Vienna: throws on outputDcVoltage below sqrt(2)*V_LL", "[vienna-topology]") {
    // 400 V_LL → sqrt(2)*400 = 565 V; setting Vdc = 500 V is below boost-PFC requirement.
    auto j = make_vienna_json(/*V_LL*/400, /*Vdc*/500);
    Vienna v(j);
    REQUIRE_THROWS_AS(v.process_design_requirements(), std::runtime_error);
}


TEST_CASE("Vienna: rejects unsupported viennaII / non-tType / synchronousRectifier in Phase 1+2",
          "[vienna-topology]") {
    {
        auto j = make_vienna_json();
        j["viennaVariant"] = "viennaII";
        Vienna v(j);
        REQUIRE_THROWS_AS(v.run_checks(true), std::runtime_error);
    }
    {
        auto j = make_vienna_json();
        j["switchType"] = "backToBackMosfet";
        Vienna v(j);
        REQUIRE_THROWS_AS(v.run_checks(true), std::runtime_error);
    }
    {
        auto j = make_vienna_json();
        j["synchronousRectifier"] = true;
        Vienna v(j);
        REQUIRE_THROWS_AS(v.run_checks(true), std::runtime_error);
    }
    {
        auto j = make_vienna_json();
        j["phaseCount"] = 2;
        Vienna v(j);
        REQUIRE_THROWS_AS(v.run_checks(true), std::runtime_error);
    }
    {
        auto j = make_vienna_json();
        j["samplingStrategy"] = "fullLineCycle";
        Vienna v(j);
        REQUIRE_THROWS_AS(v.run_checks(true), std::runtime_error);
    }
}


TEST_CASE("Vienna: missing operating points fails run_checks", "[vienna-topology]") {
    auto j = make_vienna_json();
    j["operatingPoints"] = json::array();
    Vienna v(j);
    REQUIRE_FALSE(v.run_checks(false));
    REQUIRE_THROWS_AS(v.run_checks(true), std::runtime_error);
}


TEST_CASE("Vienna: AdvancedVienna user override of L takes precedence", "[vienna-topology]") {
    auto j = make_vienna_json();
    j["desiredBoostInductance"] = 250e-6;
    AdvancedVienna adv(j);
    auto dr = adv.process_design_requirements();
    REQUIRE_THAT(adv.get_computed_boost_inductance(), WithinRel(250e-6, 1e-12));
    REQUIRE_THAT(dr.get_magnetizing_inductance().get_nominal().value(),
                 WithinRel(250e-6, 1e-9));
}


TEST_CASE("Vienna: switch RMS, diode avg and peak ripple match closed-form predictions",
          "[vienna-topology]") {
    auto j = make_vienna_json(/*V_LL*/400, /*Vdc*/750, /*Fsw*/70e3,
                              /*P*/10e3, /*ripple*/0.25);
    Vienna v(j);
    v.process_operating_points(std::vector<double>{}, 0.0);

    double V_pk_phase = std::sqrt(2.0) * 400.0 / std::sqrt(3.0);
    double M          = V_pk_phase / 375.0;
    double I_pk       = v.get_computed_line_peak_current();

    // Switch RMS via Hartmann ETH 19755 (2011) closed form.
    double Irms_expected = I_pk * std::sqrt(0.25 - 2.0 * M / (3.0 * M_PI));
    REQUIRE_THAT(v.get_last_switch_rms_current(),
                 WithinRel(Irms_expected, 1e-6));

    // Peak inductor current = I_pk + DeltaI/2.
    double L          = v.get_computed_boost_inductance();
    double DeltaI_pp  = V_pk_phase * (1.0 - M) / (L * 70e3);
    REQUIRE_THAT(v.get_last_inductor_ripple_peak_to_peak(),
                 WithinRel(DeltaI_pp, 1e-6));
    REQUIRE_THAT(v.get_last_inductor_peak_current(),
                 WithinRel(I_pk + DeltaI_pp / 2.0, 1e-6));
}
