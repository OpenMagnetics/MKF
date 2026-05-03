// =============================================================================
// TestTopologyAsymmetricHalfBridge.cpp — Asymmetric Half-Bridge (AHB)
//
// P1 (stub) coverage only: constructor, run_checks, "not implemented" guards.
// Phases P2-P12 will add static-helper, operating-point, NRMSE, and SPICE
// tests as they land. See ASYMMETRIC_HALF_BRIDGE_PLAN.md.
// =============================================================================

#include "converter_models/AsymmetricHalfBridge.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <nlohmann/json.hpp>

using nlohmann::json;
// Do NOT `using namespace OpenMagnetics;` — both OpenMagnetics::AsymmetricHalfBridge
// and MAS::AsymmetricHalfBridge are visible (the former inherits from the latter
// and the header pulls in `using namespace MAS;` inside OpenMagnetics). Qualify
// each call site explicitly via type aliases, mirroring TestTopologyDab.cpp.
using AHB    = OpenMagnetics::AsymmetricHalfBridge;
using AdvAHB = OpenMagnetics::AdvancedAsymmetricHalfBridge;
using OpenMagnetics::AhbRectifierType;
using OpenMagnetics::DimensionWithTolerance;

namespace {

// Minimal valid AHB JSON: 100 V -> 5 V / 20 A, 200 kHz, D=0.45.
// Mirrors the TI SLUP223 worked example identified in the AHB plan §3.
// baseOperatingPoint fields are inlined flat into AhbOperatingPoint by
// quicktype's allOf handling — no nested `conditions` object.
json minimal_valid_ahb() {
    return json{
        {"inputVoltage", {{"nominal", 100.0}}},
        {"operatingPoints", json::array({
            json{
                {"ambientTemperature", 25.0},
                {"switchingFrequency", 200000.0},
                {"outputVoltages", {5.0}},
                {"outputCurrents", {20.0}},
                {"dutyCycle", 0.45}
            }
        })}
    };
}

// Helper: capture the message of a thrown std::exception.
template <typename Fn>
std::string capture_throw_message(Fn&& fn) {
    try {
        fn();
    } catch (const std::exception& e) {
        return e.what();
    }
    return "";
}

} // namespace


TEST_CASE("AHB P1: minimal valid config constructs without throwing", "[ahb-topology][P1]") {
    json j = minimal_valid_ahb();
    REQUIRE_NOTHROW(AHB(j));

    AHB ahb(j);
    CHECK(ahb.get_operating_points().size() == 1);
    CHECK(ahb.get_input_voltage().get_nominal().has_value());
    CHECK(ahb.get_input_voltage().get_nominal().value() == 100.0);
    CHECK(ahb.get_operating_points()[0].get_duty_cycle() == 0.45);
}


TEST_CASE("AHB P1: run_checks accepts a well-formed config", "[ahb-topology][P1]") {
    AHB ahb(minimal_valid_ahb());
    CHECK(ahb.run_checks(false) == true);
    CHECK_NOTHROW(ahb.run_checks(true));
}


TEST_CASE("AHB P1: run_checks rejects empty operating-point list", "[ahb-topology][P1]") {
    AHB ahb;
    DimensionWithTolerance vIn;
    vIn.set_nominal(100.0);
    ahb.set_input_voltage(vIn);
    ahb.set_operating_points({});
    CHECK(ahb.run_checks(false) == false);
    CHECK_THROWS_AS(ahb.run_checks(true), std::runtime_error);
}


TEST_CASE("AHB P1: run_checks rejects missing input-voltage spec", "[ahb-topology][P1]") {
    AHB ahb(minimal_valid_ahb());
    DimensionWithTolerance vIn;  // no nominal/min/max
    ahb.set_input_voltage(vIn);
    CHECK(ahb.run_checks(false) == false);
    CHECK_THROWS_AS(ahb.run_checks(true), std::runtime_error);
}


TEST_CASE("AHB P1: schema + setters reject out-of-range duty cycle", "[ahb-topology][P1]") {
    // The schema bounds dutyCycle to [0,1] and quicktype enforces this in
    // BOTH from_json AND set_duty_cycle. Verify both layers reject 1.5.
    json j = minimal_valid_ahb();
    j["operatingPoints"][0]["dutyCycle"] = 1.5;
    CHECK_THROWS(AHB(j));

    AHB ahb(minimal_valid_ahb());
    auto ops = ahb.get_operating_points();
    CHECK_THROWS(ops[0].set_duty_cycle(1.5));
}


TEST_CASE("AHB P1: run_checks rejects non-positive switching frequency", "[ahb-topology][P1]") {
    json j = minimal_valid_ahb();
    j["operatingPoints"][0]["switchingFrequency"] = 0.0;
    AHB ahb(j);
    CHECK(ahb.run_checks(false) == false);
    CHECK_THROWS_AS(ahb.run_checks(true), std::runtime_error);
}


TEST_CASE("AHB P1: run_checks rejects out-of-range maximumDutyCycle", "[ahb-topology][P1]") {
    json j = minimal_valid_ahb();
    j["maximumDutyCycle"] = 0.0;
    AHB ahb(j);
    CHECK(ahb.run_checks(false) == false);
    CHECK_THROWS_AS(ahb.run_checks(true), std::runtime_error);
}


TEST_CASE("AHB P1: rectifierType round-trips through JSON", "[ahb-topology][P1]") {
    for (auto rt : {AhbRectifierType::CENTER_TAPPED,
                    AhbRectifierType::CURRENT_DOUBLER,
                    AhbRectifierType::FULL_BRIDGE,
                    AhbRectifierType::AHB_FLYBACK}) {
        json j = minimal_valid_ahb();
        json rtJson; to_json(rtJson, rt);
        j["rectifierType"] = rtJson;
        AHB ahb(j);
        REQUIRE(ahb.get_rectifier_type().has_value());
        CHECK(ahb.get_rectifier_type().value() == rt);
        CHECK(ahb.run_checks(false) == true);
    }
}


TEST_CASE("AHB P1: stubbed methods throw with plan reference", "[ahb-topology][P1]") {
    AHB ahb(minimal_valid_ahb());

    CHECK_THROWS_AS(ahb.process_design_requirements(), std::runtime_error);
    CHECK(capture_throw_message([&]{ ahb.process_design_requirements(); })
              .find("P6") != std::string::npos);

    CHECK_THROWS_AS(ahb.process_operating_points({1.0, 10.0}, 50e-6),
                    std::runtime_error);
    CHECK(capture_throw_message([&]{
        ahb.process_operating_points({1.0, 10.0}, 50e-6);
    }).find("P3") != std::string::npos);

    CHECK_THROWS_AS(ahb.generate_ngspice_circuit({1.0, 10.0}, 50e-6),
                    std::runtime_error);
    CHECK(capture_throw_message([&]{
        ahb.generate_ngspice_circuit({1.0, 10.0}, 50e-6);
    }).find("P7") != std::string::npos);
}


// =============================================================================
// P2: Static analytical helpers
//
// Coverage targets (per ASYMMETRIC_HALF_BRIDGE_PLAN.md §10):
//   - V_Cb=(1-D)Vin (E1)
//   - non-monotonic gain peaking at D=0.5 with M_max=1/(2n) for CT/FB
//     and 1/(4n) for CD; AHB-Flyback is monotonic
//   - duty cycle ↔ conversion ratio inverse round-trip on lower branch
//   - turns-ratio sizing equal to gain-formula inverse when Vd=0, eta=1
//   - Lo, Lm, Cb sizing identities
//   - ZVS energy balance sign change
//   - Dead-time = pi*sqrt(2*Llk*Coss)
//   - Steady-state flux excursion peaks at D=0.5
//   - Transient flux excursion >= steady-state
//   - All helpers reject malformed inputs with std::invalid_argument
// =============================================================================

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("AHB P2: V_Cb = (1-D)*Vin (E1)", "[ahb-topology][P2]") {
    const double Vin = 100.0;
    for (double D : {0.2, 0.4, 0.5, 0.6, 0.8}) {
        CHECK_THAT(AHB::compute_dc_blocking_cap_voltage(Vin, D),
                   WithinRel((1.0 - D) * Vin, 1e-12));
    }
}


TEST_CASE("AHB P2: gain curve is non-monotonic with peak at D=0.5",
          "[ahb-topology][P2]") {
    const double n = 5.0;

    SECTION("CENTER_TAPPED peaks at 1/(2n)") {
        const double M_peak = AHB::compute_conversion_ratio(
            0.5, n, AhbRectifierType::CENTER_TAPPED);
        CHECK_THAT(M_peak, WithinRel(1.0 / (2.0 * n), 1e-12));
        for (double D : {0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9}) {
            const double M = AHB::compute_conversion_ratio(
                D, n, AhbRectifierType::CENTER_TAPPED);
            CHECK(M < M_peak);
            // Symmetry: M(D) == M(1-D)
            const double M_mirror = AHB::compute_conversion_ratio(
                1.0 - D, n, AhbRectifierType::CENTER_TAPPED);
            CHECK_THAT(M, WithinRel(M_mirror, 1e-12));
        }
    }

    SECTION("FULL_BRIDGE matches CENTER_TAPPED") {
        for (double D : {0.2, 0.4, 0.5, 0.6}) {
            const double M_ct = AHB::compute_conversion_ratio(
                D, n, AhbRectifierType::CENTER_TAPPED);
            const double M_fb = AHB::compute_conversion_ratio(
                D, n, AhbRectifierType::FULL_BRIDGE);
            CHECK_THAT(M_ct, WithinRel(M_fb, 1e-12));
        }
    }

    SECTION("CURRENT_DOUBLER peaks at 1/(4n)") {
        const double M_peak = AHB::compute_conversion_ratio(
            0.5, n, AhbRectifierType::CURRENT_DOUBLER);
        CHECK_THAT(M_peak, WithinRel(1.0 / (4.0 * n), 1e-12));
    }

    SECTION("AHB_FLYBACK is monotonic and singular at D->1") {
        double prev = 0.0;
        for (double D : {0.1, 0.3, 0.5, 0.7, 0.9}) {
            const double M = AHB::compute_conversion_ratio(
                D, n, AhbRectifierType::AHB_FLYBACK);
            CHECK(M > prev);   // strictly monotonic
            prev = M;
        }
    }
}


TEST_CASE("AHB P2: duty-cycle ↔ conversion-ratio round-trips on lower branch",
          "[ahb-topology][P2]") {
    const double Vin = 100.0;
    const double n   = 5.0;

    for (auto rect : {AhbRectifierType::CENTER_TAPPED,
                      AhbRectifierType::FULL_BRIDGE,
                      AhbRectifierType::CURRENT_DOUBLER,
                      AhbRectifierType::AHB_FLYBACK}) {
        for (double D : {0.10, 0.20, 0.30, 0.40, 0.45}) {
            const double M  = AHB::compute_conversion_ratio(D, n, rect);
            const double Vo = M * Vin;
            const double D_back = AHB::compute_duty_cycle(
                Vin, Vo, n, /*Vd=*/0.0, /*eta=*/1.0, rect);
            CHECK_THAT(D_back, WithinAbs(D, 1e-9));
        }
    }
}


TEST_CASE("AHB P2: compute_duty_cycle rejects targets above the AHB peak",
          "[ahb-topology][P2]") {
    // CT peak gain = 1/(2n) = 0.1 for n=5. Asking for Vo/Vin = 0.2 is impossible.
    CHECK_THROWS_AS(
        AHB::compute_duty_cycle(100.0, 20.0, 5.0, 0.0, 1.0,
                                AhbRectifierType::CENTER_TAPPED),
        std::invalid_argument);
}


TEST_CASE("AHB P2: turns-ratio sizing matches conversion ratio (Vd=0)",
          "[ahb-topology][P2]") {
    const double Vin_min = 80.0;
    const double Vo      = 5.0;

    for (auto rect : {AhbRectifierType::CENTER_TAPPED,
                      AhbRectifierType::FULL_BRIDGE,
                      AhbRectifierType::CURRENT_DOUBLER,
                      AhbRectifierType::AHB_FLYBACK}) {
        for (double D_max : {0.20, 0.35, 0.45}) {
            const double n = AHB::compute_turns_ratio(
                Vin_min, Vo, D_max, /*Vd=*/0.0, rect);
            // Reverse: M(D_max, n) * Vin_min should equal Vo.
            const double M = AHB::compute_conversion_ratio(D_max, n, rect);
            CHECK_THAT(M * Vin_min, WithinRel(Vo, 1e-10));
        }
    }
}


TEST_CASE("AHB P2: Lo_min sizing identity (S2)", "[ahb-topology][P2]") {
    const double Vo  = 5.0;
    const double dI  = 1.0;     // 1 A pp
    const double fsw = 200e3;
    for (double D : {0.2, 0.4, 0.5, 0.6}) {
        const double Lo = AHB::compute_lo_min(Vo, D, dI, fsw);
        CHECK_THAT(Lo * dI * fsw,
                   WithinRel(Vo * (1.0 - 2.0 * D * (1.0 - D)), 1e-12));
    }
}


TEST_CASE("AHB P2: Lm_min for ZVS (S4)", "[ahb-topology][P2]") {
    const double Vin = 100.0, D = 0.45, Tsw = 5e-6, Im = 0.5;
    const double Lm = AHB::compute_lm_min_for_zvs(Vin, D, Tsw, Im);
    CHECK_THAT(Lm, WithinRel((1.0 - D) * Vin * D * Tsw / (2.0 * Im), 1e-12));
}


TEST_CASE("AHB P2: Cb_min sizing (S3, AN-4153 eq. 16)",
          "[ahb-topology][P2]") {
    const double Ipri = 5.0, D = 0.45, dV = 5.0, fsw = 200e3;
    const double Cb = AHB::compute_cb_min(Ipri, D, dV, fsw);
    CHECK_THAT(Cb, WithinRel((Ipri * D) / (fsw * dV), 1e-12));
}


TEST_CASE("AHB P2: ZVS energy balance changes sign with Llk",
          "[ahb-topology][P2]") {
    const double Lm_refl = 0, Ipri = 1.0, Coss = 200e-12, Vin = 400.0;
    // Threshold: 0.5*Llk*Ipri^2 = 2*Coss*Vin^2 => Llk_t = 4*Coss*Vin^2/Ipri^2
    const double Llk_threshold = 4.0 * Coss * Vin * Vin / (Ipri * Ipri);
    CHECK(AHB::compute_zvs_energy_balance(
              Llk_threshold * 0.5, Lm_refl, Ipri, Coss, Vin) < 0.0);
    CHECK(AHB::compute_zvs_energy_balance(
              Llk_threshold * 2.0, Lm_refl, Ipri, Coss, Vin) > 0.0);
    CHECK_THAT(AHB::compute_zvs_energy_balance(
                   Llk_threshold, Lm_refl, Ipri, Coss, Vin),
               WithinAbs(0.0, 1e-15));
}


TEST_CASE("AHB P2: dead-time = pi*sqrt(2*Llk*Coss) (Z2)",
          "[ahb-topology][P2]") {
    const double Llk = 1e-6, Coss = 200e-12;
    CHECK_THAT(AHB::compute_dead_time(Llk, Coss),
               WithinRel(M_PI * std::sqrt(2.0 * Llk * Coss), 1e-12));
}


TEST_CASE("AHB P2: steady-state flux excursion peaks at D=0.5 (F4)",
          "[ahb-topology][P2]") {
    const double Vin = 100.0, Tsw = 5e-6, Np = 20, Ae = 50e-6;

    const double B_peak = AHB::compute_steady_state_flux_excursion(
        Vin, 0.5, Tsw, Np, Ae);
    CHECK_THAT(B_peak,
               WithinRel(0.25 * Vin * Tsw / (2.0 * Np * Ae), 1e-12));

    for (double D : {0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9}) {
        const double B = AHB::compute_steady_state_flux_excursion(
            Vin, D, Tsw, Np, Ae);
        CHECK(B < B_peak);
        // Symmetric in D.
        const double B_mirror = AHB::compute_steady_state_flux_excursion(
            Vin, 1.0 - D, Tsw, Np, Ae);
        CHECK_THAT(B, WithinRel(B_mirror, 1e-12));
    }
}


TEST_CASE("AHB P2: transient flux excursion >= steady-state when dVin>0",
          "[ahb-topology][P2]") {
    const double Vin = 100.0, D = 0.4, Tsw = 5e-6, Np = 20, Ae = 50e-6;
    const double Lm = 50e-6, Rdcr = 0.1;

    const double B_ss = AHB::compute_steady_state_flux_excursion(
        Vin, D, Tsw, Np, Ae);
    const double B_tr_zero = AHB::compute_transient_flux_excursion(
        Vin, /*dVin=*/0.0, D, Tsw, Np, Ae, Lm, Rdcr);
    CHECK_THAT(B_tr_zero, WithinRel(B_ss, 1e-12));   // dVin=0 -> equality

    for (double dVin : {1.0, 5.0, 20.0}) {
        const double B_tr = AHB::compute_transient_flux_excursion(
            Vin, dVin, D, Tsw, Np, Ae, Lm, Rdcr);
        CHECK(B_tr > B_ss);
    }
}


TEST_CASE("AHB P2: helpers reject malformed inputs",
          "[ahb-topology][P2]") {
    using R = AhbRectifierType;

    CHECK_THROWS_AS(AHB::compute_dc_blocking_cap_voltage(-1.0, 0.5),
                    std::invalid_argument);
    CHECK_THROWS_AS(AHB::compute_dc_blocking_cap_voltage(100.0, 1.5),
                    std::invalid_argument);
    CHECK_THROWS_AS(AHB::compute_dc_blocking_cap_voltage(
                        std::nan(""), 0.5), std::invalid_argument);

    CHECK_THROWS_AS(AHB::compute_conversion_ratio(0.5, -1.0, R::CENTER_TAPPED),
                    std::invalid_argument);
    CHECK_THROWS_AS(AHB::compute_conversion_ratio(1.0, 5.0, R::AHB_FLYBACK),
                    std::invalid_argument);   // singular

    CHECK_THROWS_AS(AHB::compute_turns_ratio(0.0, 5.0, 0.4, 0.0,
                                             R::CENTER_TAPPED),
                    std::invalid_argument);
    CHECK_THROWS_AS(AHB::compute_turns_ratio(100.0, 5.0, 0.0, 0.0,
                                             R::CENTER_TAPPED),
                    std::invalid_argument);   // DDc=0

    CHECK_THROWS_AS(AHB::compute_duty_cycle(100.0, 5.0, 5.0, -1.0, 1.0,
                                            R::CENTER_TAPPED),
                    std::invalid_argument);   // Vd<0
    CHECK_THROWS_AS(AHB::compute_duty_cycle(100.0, 5.0, 5.0, 0.0, 1.5,
                                            R::CENTER_TAPPED),
                    std::invalid_argument);   // eta>1

    CHECK_THROWS_AS(AHB::compute_lo_min(5.0, 0.4, 0.0, 200e3),
                    std::invalid_argument);   // dILo=0
    CHECK_THROWS_AS(AHB::compute_lm_min_for_zvs(100.0, 0.4, 5e-6, 0.0),
                    std::invalid_argument);   // Im=0
    CHECK_THROWS_AS(AHB::compute_cb_min(5.0, 0.4, 0.0, 200e3),
                    std::invalid_argument);   // dVCb=0
    CHECK_THROWS_AS(AHB::compute_zvs_energy_balance(1e-6, -1.0, 1.0,
                                                    200e-12, 400.0),
                    std::invalid_argument);   // Lm_refl<0
    CHECK_THROWS_AS(AHB::compute_zvs_energy_balance(1e-6, 0.0, -1.0,
                                                    200e-12, 400.0),
                    std::invalid_argument);   // Ipri<0
    CHECK_THROWS_AS(AHB::compute_dead_time(0.0, 200e-12),
                    std::invalid_argument);
    CHECK_THROWS_AS(AHB::compute_steady_state_flux_excursion(
                        100.0, 0.4, 5e-6, 0.0, 50e-6),
                    std::invalid_argument);   // Np=0
    CHECK_THROWS_AS(AHB::compute_transient_flux_excursion(
                        100.0, -1.0, 0.4, 5e-6, 20, 50e-6, 50e-6, 0.1),
                    std::invalid_argument);   // dVin<0
}


TEST_CASE("AHB P1: AdvancedAsymmetricHalfBridge constructs and round-trips", "[ahb-topology][P1]") {
    json j = minimal_valid_ahb();
    j["desiredTurnsRatios"] = json::array({10.0});
    j["desiredMagnetizingInductance"] = 50e-6;
    j["desiredLeakageInductance"] = 1e-6;
    j["desiredOutputInductance"] = 1e-6;
    j["desiredDcBlockingCapacitance"] = 0.47e-6;
    j["desiredOutputCapacitance"] = 220e-6;

    REQUIRE_NOTHROW(AdvAHB(j));
    AdvAHB ahb(j);

    CHECK(ahb.get_desired_turns_ratios().size() == 1);
    CHECK(ahb.get_desired_turns_ratios()[0] == 10.0);
    CHECK(ahb.get_desired_magnetizing_inductance() == 50e-6);
    REQUIRE(ahb.get_desired_leakage_inductance().has_value());
    CHECK(ahb.get_desired_leakage_inductance().value() == 1e-6);

    json out;
    to_json(out, ahb);
    CHECK(out["desiredMagnetizingInductance"] == 50e-6);
    CHECK(out["desiredTurnsRatios"][0] == 10.0);

    CHECK_THROWS_AS(ahb.process(), std::runtime_error);
    CHECK(capture_throw_message([&]{ ahb.process(); })
              .find("P12") != std::string::npos);
}
