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
#include <catch2/catch_approx.hpp>
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
        // Phase 3 (this file) lifted fullLineCycle out of the deferred list.
        // peakOfLinePlusSectors is still gated.
        auto j = make_vienna_json();
        j["samplingStrategy"] = "peakOfLinePlusSectors";
        Vienna v(j);
        REQUIRE_THROWS_AS(v.run_checks(true), std::runtime_error);
    }
}


// =========================================================================
// Phase 3 — fullLineCycle sampling.
//
// peakOfLineOnly (default) emits three identical switching-period snapshots
// because the analytical model represents each phase at its own peak-of-
// line, which is the same waveform by symmetry. The wizard plot of those
// three windings looked indistinguishable. fullLineCycle emits the full
// line cycle per phase, with the three phases shifted by 120° in time, so
// the plot now shows a recognisable 3-phase pattern.
//
// This test pins:
//   (1) run_checks accepts fullLineCycle.
//   (2) Each winding carries a line-cycle-long waveform (period = 1/F_line,
//       not 1/Fsw).
//   (3) Phase A current at angle θ ≈ Phase B current at angle (θ − 120°),
//       and Phase C at (θ + 120°). This is the heart of the "phases must
//       not be identical" complaint.
//   (4) Peak inductor current is within tolerance of I_pk computed from
//       the analytical helper.
// =========================================================================
TEST_CASE("Vienna: fullLineCycle emits 120-degree-shifted phase waveforms",
          "[vienna-topology][phase3-line-cycle]") {
    auto j = make_vienna_json();
    j["samplingStrategy"] = "fullLineCycle";
    Vienna v(j);

    // (1) run_checks must now accept fullLineCycle.
    REQUIRE(v.run_checks(false));

    auto inputs = v.process();
    const auto& ops = inputs.get_operating_points();
    REQUIRE(!ops.empty());
    const auto& exc = ops[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 3);

    const double F_line = 50.0;  // make_vienna_json() default
    const double T_line = 1.0 / F_line;

    // (2) Each winding's waveform spans one line period.
    auto get_current = [&](size_t w) {
        REQUIRE(exc[w].get_current().has_value());
        auto sig = exc[w].get_current().value();
        REQUIRE(sig.get_waveform().has_value());
        return sig.get_waveform().value();
    };
    Waveform wfA = get_current(0);
    Waveform wfB = get_current(1);
    Waveform wfC = get_current(2);

    REQUIRE(wfA.get_time().has_value());
    REQUIRE(wfA.get_data().size() == wfA.get_time().value().size());
    REQUIRE(wfA.get_time().value().back() ==
            Catch::Approx(T_line).epsilon(0.001));

    // (3) Time-domain phase relationship. Sample A at t0 ≈ T/4 (= π/2, peak
    // of phase A), then locate the same value on phase B at t = t0 + T/3
    // (120° later, since B trails A by 120° in positive-sequence). Use
    // index lookup so the test isn't fragile to the exact sample count.
    //
    // A leads B by 120° means: A(t) = B(t + T/3) = C(t − T/3).
    const size_t N = wfA.get_data().size();
    const size_t iAtQuarter      = N / 4;                          // ≈ T/4
    const size_t iAtQuarter_plus = (iAtQuarter + N / 3) % N;       // ≈ T/4 + T/3
    const size_t iAtQuarter_minus= (iAtQuarter + N - N / 3) % N;   // ≈ T/4 − T/3

    double Ipeak = v.get_last_inductor_peak_current();
    REQUIRE(Ipeak > 0.0);

    // The line-cycle envelope is exact but the switching ripple at the
    // 4096-sample / 100 kHz Fsw / 50 Hz F_line ratio is heavily aliased
    // (≈2 samples per switching cycle). That makes the exact value at any
    // individual sample drift by up to ~ΔI_pp/2 (≈10-15 % of I_pk for the
    // default ripple ratio). Use a generous 25 %·I_pk margin — still far
    // tighter than the bug we're protecting against, where A = B = C
    // exactly, and where the line envelope at iAtQuarter_plus would be
    // shifted by ≈sqrt(3)/2 · I_pk ≈ 0.87·I_pk from A's reading.
    const double margin = 0.25 * Ipeak;

    // A at its own peak — envelope ≈ Ipeak * sin(π/2) = Ipeak.
    CHECK(wfA.get_data()[iAtQuarter] == Catch::Approx(Ipeak).margin(margin));
    // B trails A by T/3 — at the same wall-clock time A is +T/3 ahead of
    // B's peak, so B(t0+T/3) ≈ A(t0).
    CHECK(wfB.get_data()[iAtQuarter_plus] ==
          Catch::Approx(wfA.get_data()[iAtQuarter]).margin(margin));
    CHECK(wfC.get_data()[iAtQuarter_minus] ==
          Catch::Approx(wfA.get_data()[iAtQuarter]).margin(margin));

    // (4) The three waveforms must NOT be element-wise identical — the bug
    // we are fixing. At iAtQuarter (A's peak), B is at angle (π/2 − 2π/3)
    // = −π/6, sin = −0.5, so B's value is ≈ −0.5·Ipeak — that's >0.5·Ipeak
    // away from A's +Ipeak reading. Same for C. Use a generous 0.5·Ipeak
    // floor that still trips on A=B=C.
    CHECK(std::abs(wfA.get_data()[iAtQuarter] - wfB.get_data()[iAtQuarter]) >
          0.5 * Ipeak);
    CHECK(std::abs(wfA.get_data()[iAtQuarter] - wfC.get_data()[iAtQuarter]) >
          0.5 * Ipeak);
}


TEST_CASE("Vienna: peakOfLineOnly default keeps identical-phases behaviour",
          "[vienna-topology][phase3-line-cycle]") {
    // Regression guard — make sure the new fullLineCycle path didn't change
    // the default behaviour. Three windings, all identical, sampled at the
    // switching period.
    auto j = make_vienna_json();
    // Explicit peakOfLineOnly (== absent default).
    j["samplingStrategy"] = "peakOfLineOnly";
    Vienna v(j);
    auto inputs = v.process();
    const auto& exc = inputs.get_operating_points()[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 3);

    auto get_current_data = [&](size_t w) {
        return exc[w].get_current().value().get_waveform().value().get_data();
    };
    auto da = get_current_data(0);
    auto db = get_current_data(1);
    auto dc = get_current_data(2);
    REQUIRE(da.size() == db.size());
    REQUIRE(da.size() == dc.size());
    for (size_t i = 0; i < da.size(); ++i) {
        CHECK(da[i] == Catch::Approx(db[i]).epsilon(1e-9));
        CHECK(da[i] == Catch::Approx(dc[i]).epsilon(1e-9));
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
