// =============================================================================
// TestTopologyClllc.cpp — CLLLC bidirectional symmetric resonant converter.
//
// CURRENT SCOPE (matches the v1 skeleton in Clllc.cpp):
//   - Construction / from_json
//   - run_checks (accept/reject paths)
//   - Static analytical helpers (compute_*)
//   - process_design_requirements (turns ratio, Lr1/Cr1/Lm sizing,
//     symmetric Lr2/Cr2 derivation, resonant-frequency identities)
//   - Stubbed methods throw with the documented "not yet implemented" tag
//   - AdvancedClllc JSON round-trip
//
// OUT OF SCOPE (pending CLLLC_PLAN.md §A.5–A.10 — v2):
//   - Full 5-state Nielsen time-domain solver
//   - Direction-aware 8-switch SPICE netlist
//   - DAB-quality reference-design PtP gates
//   - VoltSecondBalance / ConverterPort / SpicePowerSanity tests
//
// This file exists primarily as a regression net for the parts of Clllc
// that DO work today (so the static helpers and design solver do not
// silently regress while the rest of the topology is built out). When
// the v2 solver lands, append additional TEST_CASEs in this same file
// rather than splitting into a sibling.
// =============================================================================

#include "converter_models/Clllc.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <nlohmann/json.hpp>

#include <cmath>

using nlohmann::json;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

using CLLLC    = OpenMagnetics::Clllc;
using AdvCLLLC = OpenMagnetics::AdvancedClllc;
using OpenMagnetics::DimensionWithTolerance;

namespace {

// Minimal valid CLLLC JSON modelled after the TI TIDM-02002 6.6 kW OBC:
// HV bus 400 V, LV bus 400 V (1:1 nominal), 250–500 kHz fsw range, single
// forward-direction OP delivering 6.6 kW at 400 V / 16.5 A.
json minimal_valid_clllc() {
    return json{
        {"highVoltageBusVoltage", {{"nominal", 400.0}}},
        {"lowVoltageBusVoltage",  {{"nominal", 400.0}}},
        {"minSwitchingFrequency", 250e3},
        {"maxSwitchingFrequency", 500e3},
        {"primaryResonantFrequency", 350e3},
        {"qualityFactor", 0.4},
        {"inductanceRatioK", 6.0},
        {"operatingPoints", json::array({
            json{
                {"ambientTemperature", 25.0},
                {"switchingFrequency", 350e3},
                {"outputVoltages", {400.0}},
                {"outputCurrents", {16.5}},
                {"powerFlowDirection", "forward"}
            }
        })}
    };
}

}  // namespace


// ─────────────────────────────────────────────────────────────────────────
// Construction
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("CLLLC: minimal valid config constructs without throwing",
          "[clllc-topology][construction]") {
    REQUIRE_NOTHROW(CLLLC(minimal_valid_clllc()));
    CLLLC c(minimal_valid_clllc());
    CHECK(c.get_high_voltage_bus_voltage().get_nominal().value() == 400.0);
    CHECK(c.get_low_voltage_bus_voltage().get_nominal().value()  == 400.0);
    CHECK(c.get_min_switching_frequency() == 250e3);
    CHECK(c.get_max_switching_frequency() == 500e3);
    CHECK(c.get_operating_points().size() == 1);
}


// ─────────────────────────────────────────────────────────────────────────
// run_checks
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("CLLLC: run_checks accepts a well-formed config",
          "[clllc-topology][run-checks]") {
    CLLLC c(minimal_valid_clllc());
    CHECK(c.run_checks(false) == true);
    CHECK_NOTHROW(c.run_checks(true));
}

TEST_CASE("CLLLC: run_checks rejects empty operating-point list",
          "[clllc-topology][run-checks]") {
    CLLLC c;
    DimensionWithTolerance hv; hv.set_nominal(400.0);
    DimensionWithTolerance lv; lv.set_nominal(400.0);
    c.set_high_voltage_bus_voltage(hv);
    c.set_low_voltage_bus_voltage(lv);
    c.set_min_switching_frequency(250e3);
    c.set_max_switching_frequency(500e3);
    c.set_operating_points({});
    CHECK(c.run_checks(false) == false);
    CHECK_THROWS(c.run_checks(true));
}

TEST_CASE("CLLLC: run_checks rejects fsw outside [min, max]",
          "[clllc-topology][run-checks]") {
    json j = minimal_valid_clllc();
    j["operatingPoints"][0]["switchingFrequency"] = 1e6;  // > maxSwitchingFrequency
    CLLLC c(j);
    CHECK(c.run_checks(false) == false);
    CHECK_THROWS(c.run_checks(true));
}

TEST_CASE("CLLLC: run_checks rejects voltage/current count mismatch",
          "[clllc-topology][run-checks]") {
    json j = minimal_valid_clllc();
    j["operatingPoints"][0]["outputVoltages"] = {400.0, 400.0};
    j["operatingPoints"][0]["outputCurrents"] = {16.5};
    CLLLC c(j);
    CHECK(c.run_checks(false) == false);
    CHECK_THROWS(c.run_checks(true));
}


// ─────────────────────────────────────────────────────────────────────────
// Static analytical helpers — values
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("CLLLC: compute_primary_resonant_frequency matches 1/(2π√(LC))",
          "[clllc-topology][static-helpers]") {
    // Lr1=20µH, Cr1=10nF → fr = 1/(2π√(2e-13)) ≈ 355.88 kHz.
    const double Lr1 = 20e-6, Cr1 = 10e-9;
    const double expected = 1.0 / (2.0 * M_PI * std::sqrt(Lr1 * Cr1));
    CHECK_THAT(CLLLC::compute_primary_resonant_frequency(Lr1, Cr1),
               WithinRel(expected, 1e-12));
}

TEST_CASE("CLLLC: compute_inductance_ratio_K = Lm / Lr1",
          "[clllc-topology][static-helpers]") {
    CHECK_THAT(CLLLC::compute_inductance_ratio_K(120e-6, 20e-6),
               WithinRel(6.0, 1e-12));
}

TEST_CASE("CLLLC: compute_quality_factor = √(L/C) / R",
          "[clllc-topology][static-helpers]") {
    // Lr1=20µH, Cr1=10nF, R=100 Ω → Q = √(20e-6/10e-9)/100 ≈ 0.4472.
    const double Lr1 = 20e-6, Cr1 = 10e-9, R = 100.0;
    const double expected = std::sqrt(Lr1 / Cr1) / R;
    CHECK_THAT(CLLLC::compute_quality_factor(Lr1, Cr1, R),
               WithinRel(expected, 1e-12));
}

TEST_CASE("CLLLC: compute_zvs_lm_max = T_dead / (16·Coss·Fs)",
          "[clllc-topology][static-helpers]") {
    // Huang SLUP263 §6 form. T_dead=200ns, Coss=200pF, Fs=350kHz.
    const double T_dead = 200e-9, Coss = 200e-12, Fs = 350e3;
    const double expected = T_dead / (16.0 * Coss * Fs);
    CHECK_THAT(CLLLC::compute_zvs_lm_max(T_dead, Coss, Fs),
               WithinRel(expected, 1e-12));
}

TEST_CASE("CLLLC: compute_turns_ratio = V_HV / V_LV",
          "[clllc-topology][static-helpers]") {
    CHECK_THAT(CLLLC::compute_turns_ratio(400.0, 200.0),
               WithinRel(2.0, 1e-12));
    CHECK_THAT(CLLLC::compute_turns_ratio(400.0, 400.0),
               WithinRel(1.0, 1e-12));
}

TEST_CASE("CLLLC: symmetric Lr2 = Lr1 / n² and Cr2 = Cr1 · n²",
          "[clllc-topology][static-helpers]") {
    const double Lr1 = 20e-6, Cr1 = 10e-9;
    for (double n : {1.0, 1.5, 2.0, 3.5}) {
        CHECK_THAT(CLLLC::compute_symmetric_lr2(Lr1, n),
                   WithinRel(Lr1 / (n * n), 1e-12));
        CHECK_THAT(CLLLC::compute_symmetric_cr2(Cr1, n),
                   WithinRel(Cr1 * n * n, 1e-12));
    }
}

TEST_CASE("CLLLC: compute_fha_gain is unity at f_n=1 (symmetric)",
          "[clllc-topology][static-helpers]") {
    // At fn=1 (operating at fr), the (1 - 1/fn²) and (fn - 1/fn) terms
    // vanish → denom = 1 → gain = 1. Holds for any K, Q.
    for (double K : {3.0, 5.0, 10.0}) {
        for (double Q : {0.2, 0.4, 0.8}) {
            CHECK_THAT(CLLLC::compute_fha_gain(1.0, K, Q, 1.0, 1.0),
                       WithinRel(1.0, 1e-9));
        }
    }
}

TEST_CASE("CLLLC: compute_fha_gain is monotonically decreasing for f_n > 1",
          "[clllc-topology][static-helpers]") {
    const double K = 6.0, Q = 0.4;
    double prev = CLLLC::compute_fha_gain(1.0, K, Q, 1.0, 1.0);
    for (double fn : {1.1, 1.2, 1.5, 2.0, 3.0}) {
        double g = CLLLC::compute_fha_gain(fn, K, Q, 1.0, 1.0);
        INFO("fn=" << fn << " gain=" << g << " prev=" << prev);
        CHECK(g < prev);
        prev = g;
    }
}


// ─────────────────────────────────────────────────────────────────────────
// Static analytical helpers — argument validation (no silent fallbacks)
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("CLLLC: static helpers throw on non-positive inputs",
          "[clllc-topology][static-helpers]") {
    CHECK_THROWS_AS(CLLLC::compute_primary_resonant_frequency(0,    1e-9), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_primary_resonant_frequency(1e-6, 0   ), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_inductance_ratio_K(0,    1e-6), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_inductance_ratio_K(1e-6, 0   ), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_quality_factor(0,    1e-9, 100), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_quality_factor(1e-6, 0,    100), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_quality_factor(1e-6, 1e-9, 0  ), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_zvs_lm_max(0, 1e-12, 1e5), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_turns_ratio(0,   1.0), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_turns_ratio(1.0, 0  ), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_symmetric_lr2(0,    1.0), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_symmetric_lr2(1e-6, 0  ), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_symmetric_cr2(0,    1.0), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_symmetric_cr2(1e-9, 0  ), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_fha_gain(0,   5.0, 0.4, 1.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(CLLLC::compute_fha_gain(1.0, 0,   0.4, 1.0, 1.0), std::invalid_argument);
}


// ─────────────────────────────────────────────────────────────────────────
// process_design_requirements — sizing identities
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("CLLLC: process_design_requirements derives expected design values",
          "[clllc-topology][design]") {
    json j = minimal_valid_clllc();
    CLLLC c(j);
    auto dr = c.process_design_requirements();

    // Turns ratio at nominal: Vhv·k_pri / (Vlv·k_sec) = 400/400 = 1.0
    // (both bridges default to fullBridge ⇒ k_pri = k_sec = 1.0).
    const double n = c.get_computed_turns_ratio();
    CHECK_THAT(n, WithinRel(1.0, 1e-12));

    // K = 6.0 → Lm = 6 · Lr1.
    const double Lr1 = c.get_computed_primary_series_inductance();
    const double Lm  = c.get_computed_magnetizing_inductance();
    REQUIRE(Lr1 > 0);
    REQUIRE(Lm  > 0);
    CHECK_THAT(Lm, WithinRel(6.0 * Lr1, 1e-9));

    // Symmetric tank at n=1: Lr2 == Lr1, Cr2 == Cr1.
    const double Lr2 = c.get_computed_secondary_series_inductance();
    const double Cr1 = c.get_computed_primary_resonant_capacitance();
    const double Cr2 = c.get_computed_secondary_resonant_capacitance();
    REQUIRE(Lr2 > 0);
    REQUIRE(Cr1 > 0);
    REQUIRE(Cr2 > 0);
    CHECK_THAT(Lr2, WithinRel(Lr1, 1e-12));
    CHECK_THAT(Cr2, WithinRel(Cr1, 1e-12));

    // f_r1 = 1/(2π√(Lr1·Cr1)) — should match the configured fr (350 kHz).
    const double fr_computed = c.get_computed_primary_resonant_frequency();
    CHECK_THAT(fr_computed, WithinRel(350e3, 1e-3));

    // f_m1 = 1/(2π√((Lr1+Lm)·Cr1)) < f_r1 (always, since Lm > 0).
    const double fm_fwd = c.get_computed_magnetizing_resonant_freq_fwd();
    CHECK(fm_fwd > 0);
    CHECK(fm_fwd < fr_computed);

    // f_m2 (reverse) — at n=1 with symmetric tank, f_m2 == f_m1 to tight tol.
    const double fm_rev = c.get_computed_magnetizing_resonant_freq_rev();
    CHECK_THAT(fm_rev, WithinRel(fm_fwd, 1e-9));

    // DesignRequirements.topology must carry the CLLLC enum.
    CHECK(dr.get_topology() == MAS::Topologies::CLLLC_RESONANT_CONVERTER);
    REQUIRE(dr.get_magnetizing_inductance().get_nominal().has_value());
    CHECK_THAT(dr.get_magnetizing_inductance().get_nominal().value(),
               WithinRel(Lm, 1e-12));
}


TEST_CASE("CLLLC: process_design_requirements honours user-supplied Lr1/Cr1",
          "[clllc-topology][design]") {
    json j = minimal_valid_clllc();
    CLLLC c(j);
    c.set_user_primary_series_inductance(15e-6);
    c.set_user_primary_resonant_capacitance(12e-9);
    c.process_design_requirements();
    CHECK_THAT(c.get_computed_primary_series_inductance(), WithinRel(15e-6, 1e-12));
    CHECK_THAT(c.get_computed_primary_resonant_capacitance(), WithinRel(12e-9, 1e-12));
}


TEST_CASE("CLLLC: process_design_requirements rejects empty OP list",
          "[clllc-topology][design]") {
    CLLLC c;
    DimensionWithTolerance hv; hv.set_nominal(400.0);
    DimensionWithTolerance lv; lv.set_nominal(400.0);
    c.set_high_voltage_bus_voltage(hv);
    c.set_low_voltage_bus_voltage(lv);
    c.set_min_switching_frequency(250e3);
    c.set_max_switching_frequency(500e3);
    c.set_operating_points({});
    CHECK_THROWS(c.process_design_requirements());
}


// ─────────────────────────────────────────────────────────────────────────
// Phase A solver — process_operating_points works for forward + reverse
// (4-state linear-ODE one-shot affine-propagator solver, see Clllc.cpp
// section "4-state linear-ODE solver — Phase A core" for theory).
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("CLLLC Phase A: process_operating_points returns one OP per HV input voltage",
          "[clllc-topology][phase-a]") {
    CLLLC c(minimal_valid_clllc());
    auto req = c.process_design_requirements();
    auto ops = c.process_operating_points(
        {c.get_computed_turns_ratio()}, c.get_computed_magnetizing_inductance());
    REQUIRE(ops.size() >= 1);
    CHECK(ops[0].get_excitations_per_winding().size() == 2);  // primary + secondary
}

TEST_CASE("CLLLC Phase A: half-wave antisymmetry residual is essentially zero",
          "[clllc-topology][phase-a][solver]") {
    CLLLC c(minimal_valid_clllc());
    c.process_design_requirements();
    c.process_operating_points({c.get_computed_turns_ratio()},
                                c.get_computed_magnetizing_inductance());
    // The one-shot affine-propagator solver enforces antisymmetry by
    // construction (modulo RK4 truncation). Residual should be << peak.
    double residual = c.get_last_steady_state_residual();
    double iPk      = c.get_last_primary_peak_current();
    REQUIRE(iPk > 0);
    CHECK(residual < 1e-3 * iPk);  // residual is a 4-state norm, peak is scalar
}

TEST_CASE("CLLLC Phase A: peak / RMS currents and Cr peak voltages are sane",
          "[clllc-topology][phase-a][diagnostics]") {
    CLLLC c(minimal_valid_clllc());
    c.process_design_requirements();
    c.process_operating_points({c.get_computed_turns_ratio()},
                                c.get_computed_magnetizing_inductance());

    double iPriPk  = c.get_last_primary_peak_current();
    double iSecPk  = c.get_last_secondary_peak_current();
    double iLmPk   = c.get_last_magnetizing_peak_current();
    double iPriRms = c.get_last_primary_rms_current();
    double iSecRms = c.get_last_secondary_rms_current();
    double vCr1Pk  = c.get_last_cr1_peak_voltage();
    double vCr2Pk  = c.get_last_cr2_peak_voltage();

    CHECK(iPriPk > 0);
    CHECK(iSecPk > 0);
    CHECK(iLmPk  > 0);
    CHECK(iPriRms > 0);
    CHECK(iSecRms > 0);
    CHECK(iPriRms <= iPriPk);   // RMS ≤ peak by definition
    CHECK(iSecRms <= iSecPk);
    CHECK(vCr1Pk > 0);
    CHECK(vCr2Pk > 0);

    // Magnetizing-current peak: physically reasonable bound is roughly
    //   I_Lm_pk ≈ V_pri · T/2 / (4·Lm)  (triangular ramp clamped to ±V_pri)
    // For 400 V / 350 kHz / Lm = K·Lr1 = 6 · Lr1, this should be < primary
    // peak (i.e. magnetizing is a smaller fraction of total tank current).
    CHECK(iLmPk < iPriPk * 2.0);  // generous bound; tighter gates are gated
                                   // behind reference-design tests in next phase.
}

TEST_CASE("CLLLC Phase A: symmetric tank → current-sharing ratio close to 1",
          "[clllc-topology][phase-a][symmetry]") {
    CLLLC c(minimal_valid_clllc());
    c.process_design_requirements();
    c.process_operating_points({c.get_computed_turns_ratio()},
                                c.get_computed_magnetizing_inductance());
    double ratio = c.get_last_current_sharing_ratio();
    // 1.0 ± 25 % — the headline diagnostic that justifies a separate model
    // from LLC. Tight ±5 % gate gated behind a reference-design test.
    CHECK_THAT(ratio, WithinAbs(1.0, 0.25));
}

TEST_CASE("CLLLC Phase A: ZVS-margin diagnostics are populated for BOTH bridges",
          "[clllc-topology][phase-a][zvs]") {
    CLLLC c(minimal_valid_clllc());
    c.process_design_requirements();
    c.process_operating_points({c.get_computed_turns_ratio()},
                                c.get_computed_magnetizing_inductance());
    // Both populated (any finite value, including negative which means
    // "ZVS lost" — that's a valid diagnostic).
    double zPri = c.get_last_zvs_margin_primary_lagging();
    double zSec = c.get_last_zvs_margin_secondary_lagging();
    CHECK(std::isfinite(zPri));
    CHECK(std::isfinite(zSec));
    CHECK(c.get_last_resonant_transition_time() > 0);
    // For 1:1 turns ratio at fsw = fr, the ZVS margins should both be
    // close to the magnetizing peak current minus the load contribution.
    CHECK(std::abs(zPri) < c.get_last_primary_peak_current() * 2);
    CHECK(std::abs(zSec) < c.get_last_secondary_peak_current() * 2);
}

TEST_CASE("CLLLC Phase A: operating region classified by fsw / fr",
          "[clllc-topology][phase-a][region]") {
    using Region = CLLLC::OperatingRegion;
    auto run_at = [](double fsw_test) {
        json j = minimal_valid_clllc();
        j["operatingPoints"][0]["switchingFrequency"] = fsw_test;
        // Make sure fsw is inside [min, max]
        j["minSwitchingFrequency"] = std::min(fsw_test, 100e3);
        j["maxSwitchingFrequency"] = std::max(fsw_test, 1e6);
        CLLLC c(j);
        c.process_design_requirements();
        c.process_operating_points({c.get_computed_turns_ratio()},
                                    c.get_computed_magnetizing_inductance());
        return std::pair<Region, double>{c.get_last_operating_region(),
                                          c.get_computed_primary_resonant_frequency()};
    };

    auto [r_at, fr1]   = run_at(350e3);          // at fr (design)
    auto [r_above, _2] = run_at(550e3);          // above fr
    auto [r_below, _3] = run_at(150e3);          // below fr (Lm dominates)
    (void)_2; (void)_3;

    CHECK(r_at    == Region::AT_RESONANCE);
    CHECK(r_above == Region::ABOVE_RESONANCE);
    CHECK(r_below == Region::BELOW_RESONANCE);
    CHECK(fr1 > 0);
}

TEST_CASE("CLLLC Phase A: forward direction is reported correctly",
          "[clllc-topology][phase-a][direction]") {
    CLLLC c(minimal_valid_clllc());
    c.process_design_requirements();
    c.process_operating_points({c.get_computed_turns_ratio()},
                                c.get_computed_magnetizing_inductance());
    CHECK(c.get_last_power_flow_direction() == CLLLC::PowerFlowDirection::FORWARD);
    CHECK(c.get_last_mode_forward() == 1);
    CHECK(c.get_last_mode_reverse() == 0);
}

TEST_CASE("CLLLC Phase A: reverse direction reuses solver via symmetric reflection",
          "[clllc-topology][phase-a][direction][reverse]") {
    json j = minimal_valid_clllc();
    j["operatingPoints"][0]["powerFlowDirection"] = "reverse";
    CLLLC c(j);
    c.process_design_requirements();
    c.process_operating_points({c.get_computed_turns_ratio()},
                                c.get_computed_magnetizing_inductance());
    CHECK(c.get_last_power_flow_direction() == CLLLC::PowerFlowDirection::REVERSE);
    CHECK(c.get_last_mode_forward() == 0);
    CHECK(c.get_last_mode_reverse() == 1);
    CHECK(c.get_last_primary_peak_current() > 0);  // physically meaningful
}

TEST_CASE("CLLLC Phase A: forward/reverse symmetry — peak currents within 5 % at 1:1",
          "[clllc-topology][phase-a][symmetry][headline]") {
    // Headline test: at a 1:1 turns ratio with a symmetric tank, the same
    // operating point run forward and reverse must produce essentially the
    // same peak primary tank current (i_Lr1 fwd ≡ i_Lr2_reflected rev). This
    // is the property that justifies CLLLC as a separate model from CLLC.
    json jf = minimal_valid_clllc();
    json jr = minimal_valid_clllc();
    jr["operatingPoints"][0]["powerFlowDirection"] = "reverse";

    CLLLC cf(jf);
    cf.process_design_requirements();
    cf.process_operating_points({cf.get_computed_turns_ratio()},
                                 cf.get_computed_magnetizing_inductance());
    double fwdPri = cf.get_last_primary_peak_current();
    double fwdSec = cf.get_last_secondary_peak_current();

    CLLLC cr(jr);
    cr.process_design_requirements();
    cr.process_operating_points({cr.get_computed_turns_ratio()},
                                 cr.get_computed_magnetizing_inductance());
    double revPri = cr.get_last_primary_peak_current();   // reverse-active side
    double revSec = cr.get_last_secondary_peak_current(); // reverse-passive side

    // At 1:1 the symmetric tank means swapping direction is a relabel: the
    // active-side peak in reverse mirrors the active-side peak in forward.
    CHECK_THAT(revPri, WithinRel(fwdSec, 0.05));   // Lr2 reverse ↔ Lr2 forward
    CHECK_THAT(revSec, WithinRel(fwdPri, 0.05));   // Lr1 reverse ↔ Lr1 forward
}


// ─────────────────────────────────────────────────────────────────────────
// Stubbed v2 surface still throws (Phase B+ work)
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("CLLLC: Phase B stubs (SPICE + extras) still throw std::logic_error",
          "[clllc-topology][v2-pending]") {
    CLLLC c(minimal_valid_clllc());
    CHECK_THROWS_AS(c.generate_ngspice_circuit({1.0}, 120e-6),
                    std::logic_error);
    CHECK_THROWS_AS(c.simulate_and_extract_operating_points({1.0}, 120e-6),
                    std::logic_error);
    CHECK_THROWS_AS(c.simulate_and_extract_topology_waveforms({1.0}, 120e-6),
                    std::logic_error);
    CHECK_THROWS_AS(c.get_extra_components_inputs(
                        OpenMagnetics::ExtraComponentsMode::IDEAL),
                    std::logic_error);
}


// ─────────────────────────────────────────────────────────────────────────
// AdvancedClllc — JSON round-trip
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("CLLLC: AdvancedClllc JSON round-trip preserves desired fields",
          "[clllc-topology][advanced]") {
    json j = minimal_valid_clllc();
    j["desiredTurnsRatios"]              = {1.0};
    j["desiredMagnetizingInductance"]    = 120e-6;
    j["desiredPrimarySeriesInductance"]  = 20e-6;
    j["desiredPrimaryResonantCapacitance"] = 10e-9;

    AdvCLLLC adv(j);
    CHECK(adv.get_desired_turns_ratios().size() == 1);
    CHECK_THAT(adv.get_desired_turns_ratios()[0], WithinRel(1.0, 1e-12));
    CHECK_THAT(adv.get_desired_magnetizing_inductance(),
               WithinRel(120e-6, 1e-12));
    REQUIRE(adv.get_desired_primary_series_inductance().has_value());
    CHECK_THAT(adv.get_desired_primary_series_inductance().value(),
               WithinRel(20e-6, 1e-12));
    REQUIRE(adv.get_desired_primary_resonant_capacitance().has_value());
    CHECK_THAT(adv.get_desired_primary_resonant_capacitance().value(),
               WithinRel(10e-9, 1e-12));

    // process() is a v2 stub — must throw rather than silently return junk.
    CHECK_THROWS_AS(adv.process(), std::logic_error);
}
