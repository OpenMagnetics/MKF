/**
 * @file TestDiagnosticsNonZero.cpp
 * @brief Sanity-check that every wizard-exposed `*Diagnostics` field for the
 *        five topologies whose getters were added in commit
 *        10c595e4 (PFC, ActiveClampForward, SingleSwitchForward,
 *        TwoSwitchForward, DifferentialModeChoke) is populated with a
 *        non-zero value on a healthy reference design.
 *
 * A zero value here means the C++ `process_*` path silently skipped its
 * assignment — that's a bug, not a normal state. Users were burned by
 * exactly this on the SRC wizard (`lastGainM = 0.0` initializer surviving
 * to the UI), so each new diagnostic gets a regression gate here.
 *
 * Reference operating points are minimal designs that exercise the
 * primary `process_*` path; specific numeric values are deliberately not
 * asserted (those live in the dedicated test files like
 * TestActiveClampForwardReferenceDesignsPtp.cpp). What we assert here is
 * the weaker "computed something" invariant.
 */
#include "converter_models/PowerFactorCorrection.h"
#include "converter_models/ActiveClampForward.h"
#include "converter_models/SingleSwitchForward.h"
#include "converter_models/TwoSwitchForward.h"
#include "converter_models/DifferentialModeChoke.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

// --------------------------------------------------------------------------
// PowerFactorCorrection
// --------------------------------------------------------------------------
TEST_CASE("PFC: diagnostics are non-zero on a healthy CCM design", "[diagnostics][pfc]") {
    nlohmann::json j = {
        {"inputVoltage", {{"minimum", 85.0}, {"maximum", 265.0}, {"nominal", 230.0}}},
        {"outputVoltage", 400.0},
        {"outputPower", 300.0},
        {"switchingFrequency", 100000.0},
        {"lineFrequency", 50.0},
        {"currentRippleRatio", 0.3},
        {"efficiency", 0.95},
        {"mode", "continuousConductionMode"},
        {"diodeVoltageDrop", 0.7},
        {"ambientTemperature", 25.0},
    };
    OpenMagnetics::PowerFactorCorrection pfc(j);
    auto dr = pfc.process_design_requirements();

    // process_design_requirements populates these without needing OPs.
    REQUIRE(pfc.get_computed_inductance() > 0.0);
    REQUIRE_FALSE(pfc.get_computed_actual_mode().empty());

    auto ops = pfc.process_operating_points({}, pfc.get_computed_inductance());
    REQUIRE_FALSE(ops.empty());

    REQUIRE(pfc.get_last_duty_cycle_peak()       > 0.0);
    REQUIRE(pfc.get_last_duty_cycle_peak()       < 1.0);
    REQUIRE(pfc.get_last_peak_inductor_current() > 0.0);
    REQUIRE(pfc.get_last_inductor_ripple()       > 0.0);
    REQUIRE(pfc.get_last_line_rms_current()      > 0.0);
    REQUIRE(pfc.get_last_input_power()           > 0.0);
}

// --------------------------------------------------------------------------
// SingleSwitchForward
// --------------------------------------------------------------------------
TEST_CASE("SingleSwitchForward: diagnostics are non-zero on a healthy design", "[diagnostics][ssf]") {
    nlohmann::json j = {
        {"inputVoltage", {{"minimum", 36.0}, {"maximum", 75.0}, {"nominal", 48.0}}},
        {"diodeVoltageDrop", 0.7},
        {"maximumDutyCycle", 0.45},
        {"currentRippleRatio", 0.3},
        {"efficiency", 0.92},
        {"operatingPoints", nlohmann::json::array({
            {
                {"outputVoltages", nlohmann::json::array({12.0})},
                {"outputCurrents", nlohmann::json::array({5.0})},
                {"switchingFrequency", 100000.0},
                {"ambientTemperature", 25.0},
            }
        })},
    };
    OpenMagnetics::SingleSwitchForward ssf(j);
    auto dr = ssf.process_design_requirements();

    REQUIRE(ssf.get_last_maximum_duty_cycle()              > 0.0);
    REQUIRE(ssf.get_last_computed_magnetizing_inductance() > 0.0);
    REQUIRE(ssf.get_last_computed_secondary_turns_ratio()  > 0.0);
    REQUIRE(ssf.get_last_reset_voltage()                   > 0.0);

    std::vector<double> turnsRatios;
    for (const auto& tr : dr.get_turns_ratios()) {
        if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
    }
    double Lm = dr.get_magnetizing_inductance().get_minimum().value();
    auto ops = ssf.process_operating_points(turnsRatios, Lm);
    REQUIRE_FALSE(ops.empty());

    REQUIRE(ssf.get_last_primary_peak_current()     > 0.0);
    REQUIRE(ssf.get_last_secondary_peak_current()   > 0.0);
    REQUIRE(ssf.get_last_magnetizing_peak_current() > 0.0);
}

// --------------------------------------------------------------------------
// TwoSwitchForward
// --------------------------------------------------------------------------
TEST_CASE("TwoSwitchForward: diagnostics are non-zero on a healthy design", "[diagnostics][tsf]") {
    nlohmann::json j = {
        {"inputVoltage", {{"minimum", 300.0}, {"maximum", 400.0}, {"nominal", 350.0}}},
        {"diodeVoltageDrop", 0.7},
        {"maximumDutyCycle", 0.45},
        {"currentRippleRatio", 0.3},
        {"efficiency", 0.92},
        {"operatingPoints", nlohmann::json::array({
            {
                {"outputVoltages", nlohmann::json::array({24.0})},
                {"outputCurrents", nlohmann::json::array({4.0})},
                {"switchingFrequency", 100000.0},
                {"ambientTemperature", 25.0},
            }
        })},
    };
    OpenMagnetics::TwoSwitchForward tsf(j);
    auto dr = tsf.process_design_requirements();

    REQUIRE(tsf.get_last_maximum_duty_cycle()              > 0.0);
    REQUIRE(tsf.get_last_computed_magnetizing_inductance() > 0.0);
    REQUIRE(tsf.get_last_computed_secondary_turns_ratio()  > 0.0);

    std::vector<double> turnsRatios;
    for (const auto& tr : dr.get_turns_ratios()) {
        if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
    }
    double Lm = dr.get_magnetizing_inductance().get_minimum().value();
    auto ops = tsf.process_operating_points(turnsRatios, Lm);
    REQUIRE_FALSE(ops.empty());

    REQUIRE(tsf.get_last_primary_peak_current()     > 0.0);
    REQUIRE(tsf.get_last_secondary_peak_current()   > 0.0);
    REQUIRE(tsf.get_last_magnetizing_peak_current() > 0.0);
}

// --------------------------------------------------------------------------
// ActiveClampForward
// --------------------------------------------------------------------------
TEST_CASE("ActiveClampForward: diagnostics are non-zero on a healthy design", "[diagnostics][acf]") {
    nlohmann::json j = {
        {"inputVoltage", {{"minimum", 36.0}, {"maximum", 75.0}, {"nominal", 48.0}}},
        {"diodeVoltageDrop", 0.7},
        {"maximumDutyCycle", 0.55},
        {"currentRippleRatio", 0.3},
        {"efficiency", 0.93},
        {"operatingPoints", nlohmann::json::array({
            {
                {"outputVoltages", nlohmann::json::array({3.3})},
                {"outputCurrents", nlohmann::json::array({20.0})},
                {"switchingFrequency", 200000.0},
                {"ambientTemperature", 25.0},
            }
        })},
    };
    OpenMagnetics::ActiveClampForward acf(j);
    auto dr = acf.process_design_requirements();

    REQUIRE(acf.get_last_maximum_duty_cycle()              > 0.0);
    REQUIRE(acf.get_last_computed_magnetizing_inductance() > 0.0);
    REQUIRE(acf.get_last_computed_secondary_turns_ratio()  > 0.0);
    REQUIRE(acf.get_last_clamp_cap_voltage()               > 0.0);  // V_clamp = D/(1-D)·Vin

    std::vector<double> turnsRatios;
    for (const auto& tr : dr.get_turns_ratios()) {
        if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
    }
    double Lm = dr.get_magnetizing_inductance().get_minimum().value();
    auto ops = acf.process_operating_points(turnsRatios, Lm);
    REQUIRE_FALSE(ops.empty());

    REQUIRE(acf.get_last_primary_peak_current()     > 0.0);
    REQUIRE(acf.get_last_secondary_peak_current()   > 0.0);
    REQUIRE(acf.get_last_magnetizing_peak_current() > 0.0);
}

// --------------------------------------------------------------------------
// DifferentialModeChoke
// --------------------------------------------------------------------------
TEST_CASE("DMC: diagnostics are non-zero on an impedance-spec design", "[diagnostics][dmc]") {
    nlohmann::json j = {
        {"configuration", "singlePhase"},
        {"lineFrequency", 50.0},
        {"ratedCurrent", 16.0},
        {"operatingCurrent", 10.0},
        {"ambientTemperature", 25.0},
        {"inputVoltage", {{"nominal", 230.0}, {"minimum", 207.0}, {"maximum", 253.0}}},
        {"minimumImpedance", nlohmann::json::array({
            {{"frequency", 150000.0},  {"impedance", {{"magnitude", 500.0}, {"phase", 90.0}}}},
            {{"frequency", 1000000.0}, {"impedance", {{"magnitude", 2000.0}, {"phase", 90.0}}}},
        })},
        {"insulationType", "no"},
        {"operatingPoints", nlohmann::json::array()},
    };
    OpenMagnetics::DifferentialModeChoke dmc(j);
    auto dr = dmc.process_design_requirements();

    REQUIRE(dmc.get_computed_inductance()           > 0.0);
    REQUIRE(dmc.get_computed_min_frequency()        > 0.0);
    REQUIRE(dmc.get_computed_max_frequency()        > dmc.get_computed_min_frequency());
    REQUIRE(dmc.get_computed_impedance_at_min_freq() > 0.0);
    REQUIRE(dmc.get_computed_number_windings()      >= 1);
}
