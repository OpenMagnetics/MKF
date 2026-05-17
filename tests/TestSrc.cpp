/**
 * @file TestSrc.cpp
 * @brief Phase 1+2 unit tests for the Series Resonant Converter (SRC).
 *
 * Coverage (SRC_PLAN.md §11 Phase 1+2):
 *   - JSON parsing and run_checks for valid / malformed inputs
 *   - process_design_requirements: turns ratio, Lr/Cr derivation
 *     (Steigerwald 1988 / FHA conventions)
 *   - process_operating_points: steady-state at resonance (Λ=1) gives
 *     M = 1 (HB) and Vout matches the expected n·Vbridge_dc
 *   - Above-resonance CCM: gain decreases monotonically with Λ (FHA)
 *   - Below-resonance: throws (deferred to Phase 3)
 *   - User overrides for Lr / Cr take precedence over Q-derived values
 *   - AdvancedSrc round-trip sets desired turns ratio and resonant components
 */
#include "converter_models/Src.h"
#include "converter_models/Topology.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <variant>

using namespace OpenMagnetics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Build a minimal half-bridge SRC JSON. fr is implicit via fmin/fmax bounds
// (geometric mean) unless resonantFrequency is set.
json make_src_json(double Vin, double Vout, double Iout,
                   double fsw, double fr,
                   double Q = 2.0,
                   const std::string& bridge = "halfBridge") {
    json j;
    j["inputVoltage"]          = {{"nominal", Vin}};
    j["minSwitchingFrequency"] = fr * 0.5;
    j["maxSwitchingFrequency"] = fr * 2.0;
    j["resonantFrequency"]     = fr;
    j["qualityFactor"]         = Q;
    j["bridgeType"]            = bridge;
    j["isolated"]              = true;

    json op;
    op["outputVoltages"]      = {Vout};
    op["outputCurrents"]      = {Iout};
    op["switchingFrequency"]  = fsw;
    op["ambientTemperature"]  = 25.0;
    j["operatingPoints"]      = json::array({op});
    return j;
}

} // namespace


TEST_CASE("SRC: JSON parses and run_checks accepts a valid HB design", "[src-topology]") {
    auto j = make_src_json(/*Vin*/400, /*Vout*/48, /*Iout*/10,
                           /*fsw*/100e3, /*fr*/100e3);
    Src src(j);
    REQUIRE(src.run_checks(/*assert*/false));
    REQUIRE(src.is_bridge_topology());
    REQUIRE(src.topology_kind() == MAS::Topologies::SERIES_RESONANT_CONVERTER);
    REQUIRE(src.get_bridge_voltage_factor() == 0.5);  // HB
}


TEST_CASE("SRC: full-bridge factor is 1.0", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0, "fullBridge");
    Src src(j);
    REQUIRE(src.get_bridge_voltage_factor() == 1.0);
}


TEST_CASE("SRC: process_design_requirements derives turns ratio and tank", "[src-topology]") {
    // HB, fr=100kHz, Q=2, Vin=400, Vout=48, Iout=10A → Rload=4.8Ω
    // n = k_bridge·Vin/Vout = 0.5·400/48 = 4.1667
    // Rac = (8·n²/π²)·Rload = 8·17.36·4.8/π² ≈ 67.5
    // Zr = Q·Rac = 135;  Lr = Zr/(2π·fr) = 215 µH;  Cr = 1/(2π·fr·Zr) = 11.8 nF
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, /*Q*/2.0);
    Src src(j);
    auto dr = src.process_design_requirements();
    REQUIRE(dr.get_topology() == Topologies::SERIES_RESONANT_CONVERTER);

    REQUIRE(dr.get_turns_ratios().size() == 1);
    auto n = dr.get_turns_ratios()[0].get_nominal().value();
    REQUIRE_THAT(n, WithinRel(0.5 * 400.0 / 48.0, 0.02));

    double Lr = src.get_computed_resonant_inductance();
    double Cr = src.get_computed_resonant_capacitance();
    REQUIRE(Lr > 0);
    REQUIRE(Cr > 0);
    double fr_computed = 1.0 / (2.0 * M_PI * std::sqrt(Lr * Cr));
    REQUIRE_THAT(fr_computed, WithinRel(100e3, 0.01));
}


TEST_CASE("SRC: at resonance (Λ=1) FHA gain M = 1 (HB)", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, /*fsw*/100e3, /*fr*/100e3, 2.0);
    Src src(j);
    auto ops = src.process_operating_points(std::vector<double>{}, /*Lm*/0);
    REQUIRE(ops.size() >= 1);

    // FHA gain at Λ=1: M = 1 / sqrt(1 + Q²·(Λ-1/Λ)²) = 1
    REQUIRE_THAT(src.get_last_gain(), WithinRel(1.0, 1e-6));
    REQUIRE_THAT(src.get_last_normalized_fsw(), WithinRel(1.0, 1e-6));
    REQUIRE(src.get_last_is_above_resonance());
    REQUIRE(src.get_last_ir_peak()  > 0);
    REQUIRE(src.get_last_vcr_peak() > 0);
}


TEST_CASE("SRC: gain decreases monotonically above resonance", "[src-topology]") {
    // Same design, sweep fsw above resonance and verify M(Λ) is monotone ↓.
    double fr = 100e3;
    double prev_M = 2.0;
    for (double Lambda : {1.0, 1.1, 1.3, 1.6, 2.0}) {
        auto j = make_src_json(400, 48, 10, /*fsw*/Lambda * fr, fr, 2.0);
        Src src(j);
        src.process_operating_points(std::vector<double>{}, 0);
        double M = src.get_last_gain();
        REQUIRE(M <= prev_M + 1e-9);
        prev_M = M;
    }
    REQUIRE(prev_M < 1.0);  // far above resonance must give M < 1
}


TEST_CASE("SRC: below-resonance operation throws (Phase 3 deferred)", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, /*fsw*/80e3, /*fr*/100e3, 2.0);
    Src src(j);
    REQUIRE_THROWS_AS(
        src.process_operating_points(std::vector<double>{}, 0),
        std::runtime_error);
}


TEST_CASE("SRC: user override of Lr/Cr replaces Q-derived values", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src src(j);
    src.set_user_resonant_inductance(150e-6);
    src.set_user_resonant_capacitance(10e-9);
    src.process_design_requirements();
    REQUIRE_THAT(src.get_computed_resonant_inductance(),  WithinRel(150e-6, 1e-9));
    REQUIRE_THAT(src.get_computed_resonant_capacitance(), WithinRel(10e-9,  1e-9));
}


TEST_CASE("SRC: missing operating points fails run_checks", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3);
    j["operatingPoints"] = json::array();
    Src src(j);
    REQUIRE(!src.run_checks(false));
    REQUIRE_THROWS_AS(src.run_checks(true), std::runtime_error);
}


TEST_CASE("SRC: AdvancedSrc round-trips desired turns ratio and Lr/Cr", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3);
    j["desiredTurnsRatios"]         = {3.5};
    j["desiredResonantInductance"]  = 200e-6;
    j["desiredResonantCapacitance"] = 12e-9;
    AdvancedSrc adv(j);
    auto dr = adv.process_design_requirements();
    REQUIRE(dr.get_turns_ratios().size() == 1);
    REQUIRE_THAT(dr.get_turns_ratios()[0].get_nominal().value(), WithinRel(3.5, 1e-6));
    REQUIRE_THAT(adv.get_computed_resonant_inductance(),  WithinRel(200e-6, 1e-9));
    REQUIRE_THAT(adv.get_computed_resonant_capacitance(), WithinRel(12e-9,  1e-9));
}


TEST_CASE("SRC: rejects unsupported rectifier types in Phase 2", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3);
    j["rectifierType"] = "centerTappedDiode";
    Src src(j);
    REQUIRE_THROWS_AS(
        src.process_operating_points(std::vector<double>{}, 0),
        std::runtime_error);
}


// =====================================================================
// get_extra_components_inputs — Cr (CAS) + Lr (MAS Inputs)
// =====================================================================

TEST_CASE("SRC: get_extra_components_inputs throws before processing", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src src(j);
    REQUIRE_THROWS_AS(
        src.get_extra_components_inputs(ExtraComponentsMode::IDEAL),
        std::runtime_error);
}


TEST_CASE("SRC: get_extra_components_inputs REAL requires magnetic", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src src(j);
    src.process_operating_points(std::vector<double>{}, 0);
    REQUIRE_THROWS_AS(
        src.get_extra_components_inputs(ExtraComponentsMode::REAL, std::nullopt),
        std::invalid_argument);
}


TEST_CASE("SRC: get_extra_components_inputs IDEAL returns Cr + Lr", "[src-topology]") {
    // HB, fr=100kHz, Q=2, Vin=400, Vout=48, Iout=10A
    auto j = make_src_json(400, 48, 10, /*fsw*/100e3, /*fr*/100e3, /*Q*/2.0);
    Src src(j);
    src.process_operating_points(std::vector<double>{}, 0);

    auto extras = src.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    REQUIRE(extras.size() == 2);

    // [0] = Cr as CAS::Inputs (RESONANT application)
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[0]));
    auto& crInputs = std::get<CAS::Inputs>(extras[0]);
    auto& crDr = crInputs.get_design_requirements();
    REQUIRE(crDr.get_capacitance().get_nominal().has_value());
    double Cr_emit = crDr.get_capacitance().get_nominal().value();
    REQUIRE_THAT(Cr_emit, WithinRel(src.get_computed_resonant_capacitance(), 1e-9));
    REQUIRE(crDr.get_rated_voltage() > 0);
    REQUIRE(crDr.get_role().has_value());
    REQUIRE(crDr.get_role().value() == CAS::Application::RESONANT);
    REQUIRE(!crInputs.get_operating_points().empty());
    {
        auto& exc = crInputs.get_operating_points()[0].get_excitation();
        REQUIRE(exc.get_voltage().has_value());
        REQUIRE(exc.get_current().has_value());
        REQUIRE(!exc.get_voltage()->get_waveform()->get_data().empty());
        REQUIRE(!exc.get_current()->get_waveform()->get_data().empty());
    }

    // [1] = Lr as MAS Inputs (always external in IDEAL)
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[1]));
    auto& lrInputs = std::get<OpenMagnetics::Inputs>(extras[1]);
    auto& lrDr = lrInputs.get_design_requirements();
    REQUIRE(lrDr.get_magnetizing_inductance().get_nominal().has_value());
    double Lr_emit = lrDr.get_magnetizing_inductance().get_nominal().value();
    REQUIRE_THAT(Lr_emit, WithinRel(src.get_computed_resonant_inductance(), 1e-9));
    REQUIRE(lrDr.get_topology().has_value());
    REQUIRE(lrDr.get_topology().value() == Topologies::SERIES_RESONANT_CONVERTER);
    REQUIRE(!lrInputs.get_operating_points().empty());
    {
        auto& exc = lrInputs.get_operating_points()[0].get_excitations_per_winding()[0];
        REQUIRE(exc.get_current().has_value());
        REQUIRE(exc.get_voltage().has_value());
    }
}


TEST_CASE("SRC: get_extra_components_inputs round-trips fr = 1/(2π√(LrCr))",
          "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src src(j);
    src.process_operating_points(std::vector<double>{}, 0);
    auto extras = src.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    REQUIRE(extras.size() == 2);
    double Cr = std::get<CAS::Inputs>(extras[0])
                    .get_design_requirements().get_capacitance()
                    .get_nominal().value();
    double Lr = std::get<OpenMagnetics::Inputs>(extras[1])
                    .get_design_requirements().get_magnetizing_inductance()
                    .get_nominal().value();
    double fr = 1.0 / (2.0 * M_PI * std::sqrt(Lr * Cr));
    REQUIRE_THAT(fr, WithinRel(100e3, 1e-3));
}
