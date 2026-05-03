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

    CHECK_THROWS_AS(ahb.generate_ngspice_circuit({10.0}, 50e-6),
                    std::runtime_error);
    CHECK(capture_throw_message([&]{
        ahb.generate_ngspice_circuit({10.0}, 50e-6);
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


// =============================================================================
// P3: Operating-point processing (CCM analytical model, V1 CT secondary)
//
// Coverage targets per ASYMMETRIC_HALF_BRIDGE_PLAN.md §10:
//   - OperatingPoints_Generation: process_operating_points returns the right
//     number of OPs with non-empty waveforms.
//   - Asymmetric_Switch_RMS: I_Q1,rms / I_Q2,rms ≃ √(D / (1-D)) at D=0.3
//     within 5 % (Imbertson-Mohan §IV.A).
//   - Steady_State_Flux_Excursion diagnostic populated.
//   - Diagnostics_Populated: every per-OP diagnostic non-zero where the
//     physics dictates non-zero.
//   - Volt-second balance on v_pri and v_Lo.
//   - SLUP223 worked example (Vin=100 V, Vo=5 V/20 A, n=10, fsw=200 kHz,
//     D=0.45) lands in the analytical ballpark.
// =============================================================================

#include "converter_models/AsymmetricHalfBridge.h"

namespace {

// Trapezoidal integral of y(t) over the entire waveform.
double trapz(const std::vector<double>& t, const std::vector<double>& y) {
    REQUIRE(t.size() == y.size());
    double sum = 0.0;
    for (size_t k = 0; k + 1 < t.size(); ++k)
        sum += 0.5 * (y[k] + y[k + 1]) * (t[k + 1] - t[k]);
    return sum;
}

double mean_of(const std::vector<double>& t, const std::vector<double>& y) {
    const double duration = t.back() - t.front();
    REQUIRE(duration > 0);
    return trapz(t, y) / duration;
}

// Build a minimal valid AHB with a chosen duty and pre-sized turns ratio for
// a centre-tapped secondary (Vo = 2D(1-D)Vin/n at this duty).
json ahb_with_duty(double Vin, double Vo, double Io, double D, double n,
                   double fsw = 200000.0) {
    return json{
        {"inputVoltage", {{"nominal", Vin}}},
        {"rectifierType", "centerTapped"},
        {"operatingPoints", json::array({
            json{
                {"ambientTemperature", 25.0},
                {"switchingFrequency", fsw},
                {"outputVoltages", {Vo}},
                {"outputCurrents", {Io}},
                {"dutyCycle", D}
            }
        })}
    };
}

} // namespace


TEST_CASE("AHB P3: process_operating_points returns one OP per (Vin × OP)",
          "[ahb-topology][P3]") {
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;   // = 4.95 V
    AHB ahb(ahb_with_duty(Vin, Vo, 20.0, D, n));

    auto ops = ahb.process_operating_points({n}, 50e-6);
    REQUIRE(ops.size() == 1);                            // 1 Vin × 1 OP
    auto& exc = ops[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 3);                            // pri + sec_a + sec_b
    CHECK(exc[0].get_current()->get_waveform()->get_data().size() > 100);
    CHECK(exc[0].get_voltage()->get_waveform()->get_data().size() > 100);
}


TEST_CASE("AHB P3: primary voltage volt-second balance is exact",
          "[ahb-topology][P3]") {
    const double Vin = 100.0, n = 10.0;
    for (double D : {0.20, 0.30, 0.45, 0.50}) {
        const double Vo = 2.0 * D * (1.0 - D) * Vin / n;
        AHB ahb(ahb_with_duty(Vin, Vo, 10.0, D, n));
        auto ops = ahb.process_operating_points({n}, 50e-6);
        const auto& vpri = ops[0].get_excitations_per_winding()[0]
                              .get_voltage()->get_waveform();
        const double vs = trapz(vpri->get_time().value(),
                                vpri->get_data());
        CHECK_THAT(vs, WithinAbs(0.0, 1e-9));
    }
}


TEST_CASE("AHB P3: output-inductor voltage volt-second balance is exact",
          "[ahb-topology][P3]") {
    const double Vin = 100.0, n = 10.0, D = 0.4;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_duty(Vin, Vo, 10.0, D, n));
    auto ops = ahb.process_operating_points({n}, 50e-6);
    // Lo is the first extra component.
    // We don't have a public getter for extras, but the op's last excitations
    // are the 2 secondaries. Lo waveform balance is implicit in the sec
    // voltages — verify v_pri/n - Vo integrates to ~0 at the rectified node:
    const auto& vpri = ops[0].get_excitations_per_winding()[0]
                          .get_voltage()->get_waveform();
    const auto time = vpri->get_time().value();
    const auto& data = vpri->get_data();
    std::vector<double> vLo;
    vLo.reserve(data.size());
    for (double v : data) vLo.push_back(std::abs(v) / n - Vo);
    CHECK_THAT(trapz(time, vLo), WithinAbs(0.0, 1e-9));
}


TEST_CASE("AHB P3: asymmetric switch RMS — I_Q1/I_Q2 = √(D/(1-D))",
          "[ahb-topology][P3]") {
    const double D = 0.3, n = 10.0, Vin = 100.0;
    const double Vo = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_duty(Vin, Vo, 10.0, D, n));
    // Use a generous Lm so the magnetizing ripple is small relative to
    // the reflected load current — this is the regime in which the
    // textbook √(D/(1-D)) ratio holds (Imbertson-Mohan §IV.A).
    ahb.process_operating_points({n}, 5e-3);

    const double Iq1 = ahb.get_last_switch_rms_current_q1();
    const double Iq2 = ahb.get_last_switch_rms_current_q2();
    REQUIRE(Iq1 > 0.0);
    REQUIRE(Iq2 > 0.0);
    const double ratio_observed = Iq1 / Iq2;
    const double ratio_textbook = std::sqrt(D / (1.0 - D));
    CHECK_THAT(ratio_observed, WithinRel(ratio_textbook, 0.05));
}


TEST_CASE("AHB P3: per-OP diagnostics are populated",
          "[ahb-topology][P3]") {
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_duty(Vin, Vo, 20.0, D, n));
    ahb.process_operating_points({n}, 50e-6);

    CHECK_THAT(ahb.get_last_duty_cycle(),                 WithinRel(D, 1e-9));
    CHECK_THAT(ahb.get_last_conversion_ratio(),           WithinRel(Vo / Vin, 1e-9));
    CHECK_THAT(ahb.get_last_dc_blocking_cap_voltage(),    WithinRel((1.0 - D) * Vin, 1e-9));
    CHECK_THAT(ahb.get_last_primary_peak_voltage_positive(),
               WithinRel((1.0 - D) * Vin, 1e-9));
    CHECK_THAT(ahb.get_last_primary_peak_voltage_negative(),
               WithinRel(-D * Vin, 1e-9));
    CHECK_THAT(ahb.get_last_switch_peak_voltage_q1(),     WithinRel(Vin, 1e-12));
    CHECK_THAT(ahb.get_last_switch_peak_voltage_q2(),     WithinRel(Vin, 1e-12));
    CHECK(ahb.get_last_switch_rms_current_q1() > 0.0);
    CHECK(ahb.get_last_switch_rms_current_q2() > 0.0);
    CHECK(ahb.get_last_resonant_transition_time() > 0.0);
    CHECK(ahb.get_last_steady_state_flux_excursion() > 0.0);
    CHECK(ahb.get_last_magnetizing_current_ripple() > 0.0);
    CHECK(ahb.get_last_output_inductor_ripple() > 0.0);
    CHECK(ahb.get_last_operating_mode() == 0);   // CCM at this design
    CHECK(ahb.get_last_rectifier_type() == 0);   // CT
}


TEST_CASE("AHB P3: steady-state flux excursion matches analytical formula",
          "[ahb-topology][P3]") {
    const double Vin = 100.0, n = 10.0, fsw = 200e3;
    for (double D : {0.20, 0.30, 0.45, 0.50}) {
        const double Vo = 2.0 * D * (1.0 - D) * Vin / n;
        AHB ahb(ahb_with_duty(Vin, Vo, 10.0, D, n, fsw));
        ahb.process_operating_points({n}, 50e-6);
        const double Tsw    = 1.0 / fsw;
        const double expected = D * (1.0 - D) * Vin * Tsw;   // V·s per turn
        CHECK_THAT(ahb.get_last_steady_state_flux_excursion(),
                   WithinRel(expected, 1e-12));
    }
}


TEST_CASE("AHB P3: nominal+min+max input voltages each produce one OP",
          "[ahb-topology][P3]") {
    const double n = 10.0, D = 0.45;
    const double Vin_nom = 100.0;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin_nom / n;
    json j = ahb_with_duty(Vin_nom, Vo, 20.0, D, n);
    j["inputVoltage"]["minimum"] = 80.0;
    j["inputVoltage"]["maximum"] = 120.0;
    AHB ahb(j);

    auto ops = ahb.process_operating_points({n}, 50e-6);
    CHECK(ops.size() == 3);   // min, nom, max
}


TEST_CASE("AHB P3: rejects malformed processing inputs",
          "[ahb-topology][P3]") {
    AHB ahb(minimal_valid_ahb());
    // turnsRatios empty
    CHECK_THROWS_AS(ahb.process_operating_points({}, 50e-6),
                    std::invalid_argument);
    // turnsRatios non-positive
    CHECK_THROWS_AS(ahb.process_operating_points({-1.0}, 50e-6),
                    std::invalid_argument);
    // magnetizingInductance non-positive
    CHECK_THROWS_AS(ahb.process_operating_points({10.0}, 0.0),
                    std::invalid_argument);
    // multi-secondary not yet supported (P4/P11)
    CHECK_THROWS_AS(ahb.process_operating_points({10.0, 5.0}, 50e-6),
                    std::invalid_argument);
}


TEST_CASE("AHB P3: SLUP223-like 100V→5V/20A worked example sanity",
          "[ahb-topology][P3]") {
    // Vin=100 V, Vo=5 V, Io=20 A, n=10, fsw=200 kHz, D=0.45.
    // Expected steady-state values (analytical):
    //   V_Cb        = 0.55 * 100   = 55 V
    //   M_target    = 5 / 100      = 0.05
    //   M_actual    = 2*0.45*0.55/10 = 0.0495 → Vo_actual = 4.95 V
    //   Im_pp       = 0.55*100*0.45*5e-6 / 50e-6 = 2.475 A
    //   td_resonant = π·sqrt(1e-6 * 2 * 200e-12) ≈ 1.987 ns
    AHB ahb(ahb_with_duty(100.0, 4.95, 20.0, 0.45, 10.0));
    ahb.process_operating_points({10.0}, 50e-6);

    CHECK_THAT(ahb.get_last_dc_blocking_cap_voltage(),
               WithinRel(55.0, 1e-9));
    CHECK_THAT(ahb.get_last_magnetizing_current_ripple(),
               WithinRel(0.55 * 100.0 * 0.45 * 5e-6 / 50e-6, 1e-9));
    // Resonant transition time uses default Coss=200 pF, default Llk=1 µH.
    // Per the AHB physics (P5 refinement) the resonant tank during dead-
    // time is Llk + Lm because the secondary diodes block (Mappus 2014,
    // Imbertson-Mohan §V); the dead-time helper is fed Llk+Lm = 51 µH
    // here.
    const double td_expected =
        M_PI * std::sqrt((1e-6 + 50e-6) * 2.0 * 200e-12);
    CHECK_THAT(ahb.get_last_resonant_transition_time(),
               WithinRel(td_expected, 1e-9));
}


TEST_CASE("AHB P3: Lo current waveform mean equals Io (extras populated)",
          "[ahb-topology][P3]") {
    // We don't have a public extras getter, but we can verify indirectly
    // via secondary-winding current means: in CT, sec_a+sec_b carry I_Lo
    // alternately, so their sum equals I_Lo at every instant. The mean
    // of (i_sec_a + i_sec_b) over one period must equal Io.
    const double Vin = 100.0, n = 10.0, D = 0.45, Io = 20.0;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_duty(Vin, Vo, Io, D, n));
    auto ops = ahb.process_operating_points({n}, 50e-6);

    const auto& iSecA = ops[0].get_excitations_per_winding()[1]
                          .get_current()->get_waveform();
    const auto& iSecB = ops[0].get_excitations_per_winding()[2]
                          .get_current()->get_waveform();
    const auto tA = iSecA->get_time().value();
    const auto& dA = iSecA->get_data();
    const auto& dB = iSecB->get_data();
    REQUIRE(dA.size() == dB.size());
    std::vector<double> sum(dA.size());
    for (size_t k = 0; k < dA.size(); ++k) sum[k] = dA[k] + dB[k];
    CHECK_THAT(mean_of(tA, sum), WithinRel(Io, 0.02));
}


// =============================================================================
// P4: Rectifier-type-aware processing — Full-Bridge, Current-Doubler, Flyback.
//
// Coverage targets per ASYMMETRIC_HALF_BRIDGE_PLAN.md §10:
//   - FB rectifier emits exactly 1 secondary winding (vs CT's 2 half-windings),
//     volt-second balance still exact, lastRectifierType == 2.
//   - CD rectifier emits 1 secondary, lastRectifierType == 1, vSec balance
//     exact, asymmetric Lo1/Lo2 inferred via Io_per_inductor split.
//   - At D=0.5 the CD design is symmetric: |dILo1| ≃ |dILo2|.
//   - AHB_FLYBACK throws with a "P10" tag (storage topology, separate scope).
// =============================================================================

namespace {

json ahb_with_duty_rect(double Vin, double Vo, double Io, double D, double n,
                        const std::string& rectifierType,
                        double fsw = 200000.0) {
    return json{
        {"inputVoltage", {{"nominal", Vin}}},
        {"rectifierType", rectifierType},
        {"operatingPoints", json::array({
            json{
                {"ambientTemperature", 25.0},
                {"switchingFrequency", fsw},
                {"outputVoltages", {Vo}},
                {"outputCurrents", {Io}},
                {"dutyCycle", D}
            }
        })}
    };
}

} // namespace


TEST_CASE("AHB P4: full-bridge rectifier emits 1 secondary winding",
          "[ahb-topology][P4]") {
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_duty_rect(Vin, Vo, 20.0, D, n, "fullBridge"));
    auto ops = ahb.process_operating_points({n}, 50e-6);
    REQUIRE(ops.size() == 1);
    const auto& exc = ops[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 2);                         // pri + single sec
    CHECK(exc[1].get_name() == "Secondary 0");
    CHECK(ahb.get_last_rectifier_type() == 2);        // FB
}


TEST_CASE("AHB P4: full-bridge primary volt-second balance is exact",
          "[ahb-topology][P4]") {
    const double Vin = 100.0, n = 10.0;
    for (double D : {0.20, 0.30, 0.45, 0.50}) {
        const double Vo = 2.0 * D * (1.0 - D) * Vin / n;
        AHB ahb(ahb_with_duty_rect(Vin, Vo, 10.0, D, n, "fullBridge"));
        auto ops = ahb.process_operating_points({n}, 50e-6);
        const auto& vpri = ops[0].get_excitations_per_winding()[0]
                              .get_voltage()->get_waveform();
        CHECK_THAT(trapz(vpri->get_time().value(), vpri->get_data()),
                   WithinAbs(0.0, 1e-9));
        const auto& vsec = ops[0].get_excitations_per_winding()[1]
                              .get_voltage()->get_waveform();
        CHECK_THAT(trapz(vsec->get_time().value(), vsec->get_data()),
                   WithinAbs(0.0, 1e-9));
    }
}


TEST_CASE("AHB P4: full-bridge secondary current is bipolar with mean 0",
          "[ahb-topology][P4]") {
    // FB single secondary winding sees +i_Lo in interval A and -i_Lo in
    // interval C. With i_Lo ≃ Io for low ripple, the period-mean is
    //   D*Io + (1-D)*(-Io) = (2D-1)*Io
    // — i.e. zero only at D=0.5; non-zero otherwise. Verify D=0.5.
    const double Vin = 100.0, n = 10.0, D = 0.5;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    const double Io  = 10.0;
    AHB ahb(ahb_with_duty_rect(Vin, Vo, Io, D, n, "fullBridge"));
    auto ops = ahb.process_operating_points({n}, 5e-3);
    const auto& iSec = ops[0].get_excitations_per_winding()[1]
                          .get_current()->get_waveform();
    const auto t = iSec->get_time().value();
    const auto& d = iSec->get_data();
    CHECK_THAT(mean_of(t, d), WithinAbs(0.0, 0.05));
    // Peak amplitude ≃ Io (large Lm so Lo ripple dominates within ±5 %).
    double peak = 0.0;
    for (double v : d) peak = std::max(peak, std::abs(v));
    CHECK_THAT(peak, WithinRel(Io, 0.10));
}


TEST_CASE("AHB P4: current-doubler rectifier emits 1 secondary winding",
          "[ahb-topology][P4]") {
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_duty_rect(Vin, Vo, 20.0, D, n, "currentDoubler"));
    auto ops = ahb.process_operating_points({n}, 50e-6);
    REQUIRE(ops.size() == 1);
    const auto& exc = ops[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 2);
    CHECK(exc[1].get_name() == "Secondary 0");
    CHECK(ahb.get_last_rectifier_type() == 1);        // CD
}


TEST_CASE("AHB P4: current-doubler primary volt-second balance is exact",
          "[ahb-topology][P4]") {
    const double Vin = 100.0, n = 10.0;
    for (double D : {0.30, 0.45, 0.50}) {
        const double Vo = 2.0 * D * (1.0 - D) * Vin / n;
        AHB ahb(ahb_with_duty_rect(Vin, Vo, 10.0, D, n, "currentDoubler"));
        auto ops = ahb.process_operating_points({n}, 50e-6);
        const auto& vpri = ops[0].get_excitations_per_winding()[0]
                              .get_voltage()->get_waveform();
        CHECK_THAT(trapz(vpri->get_time().value(), vpri->get_data()),
                   WithinAbs(0.0, 1e-9));
        const auto& vsec = ops[0].get_excitations_per_winding()[1]
                              .get_voltage()->get_waveform();
        CHECK_THAT(trapz(vsec->get_time().value(), vsec->get_data()),
                   WithinAbs(0.0, 1e-9));
    }
}


TEST_CASE("AHB P4: current-doubler primary current asymmetry encodes Io/2 split",
          "[ahb-topology][P4]") {
    // In CD each output inductor carries Io/2. The reflected primary
    // current peak is therefore (Io/2)/n (plus i_Lm), not Io/n. Check
    // by comparing CT vs CD primary peaks at the same operating point
    // and confirming the CD primary peak is roughly half (within Im
    // contribution).
    const double Vin = 100.0, n = 10.0, D = 0.45, Io = 20.0;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    const double Lm  = 5e-3;   // large -> i_Lm small vs reflected load

    AHB ahbCT(ahb_with_duty_rect(Vin, Vo, Io, D, n, "centerTapped"));
    auto opsCT = ahbCT.process_operating_points({n}, Lm);
    AHB ahbCD(ahb_with_duty_rect(Vin, Vo, Io, D, n, "currentDoubler"));
    auto opsCD = ahbCD.process_operating_points({n}, Lm);

    auto peak_abs = [](const std::optional<OpenMagnetics::Waveform>& w) {
        double p = 0.0;
        for (double v : w->get_data()) p = std::max(p, std::abs(v));
        return p;
    };
    const double iPriPkCT = peak_abs(opsCT[0].get_excitations_per_winding()[0]
                                        .get_current()->get_waveform());
    const double iPriPkCD = peak_abs(opsCD[0].get_excitations_per_winding()[0]
                                        .get_current()->get_waveform());
    CHECK_THAT(iPriPkCD, WithinRel(iPriPkCT / 2.0, 0.10));
}


TEST_CASE("AHB P4: AHB_FLYBACK rectifier is deferred to P10 (throws)",
          "[ahb-topology][P4]") {
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_duty_rect(Vin, Vo, 20.0, D, n, "ahbFlyback"));
    CHECK_THROWS_AS(ahb.process_operating_points({n}, 50e-6),
                    std::runtime_error);
    CHECK(capture_throw_message([&]{
              ahb.process_operating_points({n}, 50e-6);
          }).find("P10") != std::string::npos);
}


// =============================================================================
// P5: Transient flux excursion + ZVS boundary + DCM detection.
//
// Coverage targets per ASYMMETRIC_HALF_BRIDGE_PLAN.md §10:
//   - Transient_Flux_Excursion: lastTransientFluxExcursionEstimate = 0
//     when inputVoltageStepRange unset; equals lambda_ss + D·dVin·Tsw
//     when set; rejects negative dVin; exceeds steady-state for any
//     positive dVin.
//   - ZVS_Boundaries: lastZvsMargin sweeps from negative (low Lm/Llk) to
//     positive as the resonant-tank energy crosses 2·Coss·Vin² and
//     lastResonantTransitionTime grows with sqrt(Llk + Lm).
//   - DCM_Detection: forcing a small Lo with a non-zero pinned value
//     drives ILo_min < 0 and toggles lastOperatingMode to 1.
// =============================================================================

namespace {

json ahb_with_step_range(double Vin, double Vo, double Io, double D, double n,
                         double inputVoltageStepRange,
                         double fsw = 200000.0) {
    json j = ahb_with_duty(Vin, Vo, Io, D, n, fsw);
    j["inputVoltageStepRange"] = inputVoltageStepRange;
    return j;
}

} // namespace


TEST_CASE("AHB P5: lastTransientFluxExcursionEstimate is zero when step range unset",
          "[ahb-topology][P5]") {
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_duty(Vin, Vo, 10.0, D, n));
    ahb.process_operating_points({n}, 50e-6);
    CHECK(ahb.get_last_transient_flux_excursion_estimate() == 0.0);
}


TEST_CASE("AHB P5: transient flux excursion equals lambda_ss + D·dVin·Tsw",
          "[ahb-topology][P5]") {
    const double Vin = 100.0, n = 10.0, D = 0.40, fsw = 200e3;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    const double dVin = 30.0;
    AHB ahb(ahb_with_step_range(Vin, Vo, 10.0, D, n, dVin, fsw));
    ahb.process_operating_points({n}, 50e-6);
    const double Tsw      = 1.0 / fsw;
    const double expected = D * (1.0 - D) * Vin * Tsw + D * dVin * Tsw;
    CHECK_THAT(ahb.get_last_transient_flux_excursion_estimate(),
               WithinRel(expected, 1e-12));
}


TEST_CASE("AHB P5: transient flux excursion strictly exceeds steady-state for dVin>0",
          "[ahb-topology][P5]") {
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_step_range(Vin, Vo, 10.0, D, n, /*dVin=*/20.0));
    ahb.process_operating_points({n}, 50e-6);
    CHECK(ahb.get_last_transient_flux_excursion_estimate() >
          ahb.get_last_steady_state_flux_excursion());
}


TEST_CASE("AHB P5: negative inputVoltageStepRange is rejected",
          "[ahb-topology][P5]") {
    // Schema enforces type=number with no minimum, so the JSON parses;
    // the rejection happens inside the operating-point processor.
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    AHB ahb(ahb_with_step_range(Vin, Vo, 10.0, D, n, /*dVin=*/-5.0));
    CHECK_THROWS_AS(ahb.process_operating_points({n}, 50e-6),
                    std::invalid_argument);
}


TEST_CASE("AHB P5: ZVS margin grows monotonically with leakage inductance",
          "[ahb-topology][P5]") {
    // With Lm and load fixed, the primary current at Q1-off is fixed,
    // so E_stored = 0.5*(Llk+Lm)*Ipri² is monotonically increasing in
    // Llk. Sweep three Llk values and verify.
    const double Vin = 100.0, n = 10.0, D = 0.45, Io = 5.0;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;

    auto margin_for_Llk = [&](double Llk) {
        AHB ahb(ahb_with_duty(Vin, Vo, Io, D, n));
        ahb.set_computed_leakage_inductance(Llk);
        ahb.process_operating_points({n}, 200e-6);
        return ahb.get_last_zvs_margin();
    };
    const double m1 = margin_for_Llk(1e-7);
    const double m2 = margin_for_Llk(1e-6);
    const double m3 = margin_for_Llk(1e-5);
    CHECK(m2 > m1);
    CHECK(m3 > m2);
}


TEST_CASE("AHB P5: resonant transition time grows with sqrt(Llk+Lm)",
          "[ahb-topology][P5]") {
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    auto td_for = [&](double Llk, double Lm) {
        AHB ahb(ahb_with_duty(Vin, Vo, 10.0, D, n));
        ahb.set_computed_leakage_inductance(Llk);
        ahb.process_operating_points({n}, Lm);
        return ahb.get_last_resonant_transition_time();
    };
    const double Llk = 1e-6, Coss = 200e-12;   // Coss is class default
    const double td1 = td_for(Llk, 50e-6);
    const double td2 = td_for(Llk, 200e-6);
    CHECK_THAT(td1,
               WithinRel(M_PI * std::sqrt((Llk + 50e-6)  * 2.0 * Coss), 1e-9));
    CHECK_THAT(td2,
               WithinRel(M_PI * std::sqrt((Llk + 200e-6) * 2.0 * Coss), 1e-9));
    CHECK(td2 > td1);
}


TEST_CASE("AHB P5: ZVS margin sign flips between heavy load and light load",
          "[ahb-topology][P5]") {
    // ZVS is a load-dependent threshold: at heavy load the reflected Io
    // dominates the primary peak current and 0.5*(Llk+Lm)*I_pri² >> 2·Coss·Vin²
    // (margin > 0). At very light load with a large Lm (so the magnetizing
    // ripple is also small), I_pri collapses to ~Io/n and the stored
    // energy can fall below the Coss-charge requirement (margin < 0).
    // This is the classic "AHB loses ZVS at light load" condition.
    const double Vin = 100.0, n = 10.0, D = 0.45;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;
    const double Lm  = 50e-3;   // 50 mH -> i_lm,pk ~ µA at this OP

    auto margin_for_Io = [&](double Io) {
        AHB ahb(ahb_with_duty(Vin, Vo, Io, D, n));
        ahb.set_computed_leakage_inductance(1e-6);
        ahb.process_operating_points({n}, Lm);
        return ahb.get_last_zvs_margin();
    };
    CHECK(margin_for_Io(0.01) < 0.0);   // light  -> ZVS lost
    CHECK(margin_for_Io(20.0) > 0.0);   // heavy  -> ZVS achieved
}


TEST_CASE("AHB P5: small pinned Lo drives DCM and toggles operating-mode flag",
          "[ahb-topology][P5]") {
    // Pin a deliberately undersized Lo so the analytical I_Lo,min < 0,
    // which triggers the DCM flag. We don't model true DCM waveforms in
    // P5 (that's P9 territory); we just surface the flag so a user can
    // see the design has crossed the CCM/DCM boundary.
    const double Vin = 100.0, n = 10.0, D = 0.45, Io = 1.0;
    const double Vo  = 2.0 * D * (1.0 - D) * Vin / n;

    AHB ahb_ccm(ahb_with_duty(Vin, Vo, Io, D, n));
    ahb_ccm.process_operating_points({n}, 50e-6);
    REQUIRE(ahb_ccm.get_last_operating_mode() == 0);   // CCM at default Lo

    AHB ahb_dcm(ahb_with_duty(Vin, Vo, Io, D, n));
    ahb_dcm.set_computed_output_inductance(1e-7);      // 0.1 µH -> giant ripple
    ahb_dcm.process_operating_points({n}, 50e-6);
    CHECK(ahb_dcm.get_last_operating_mode() == 1);     // DCM flagged
    CHECK(ahb_dcm.get_last_output_inductor_ripple() > 2.0 * Io);
}
