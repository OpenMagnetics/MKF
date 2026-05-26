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


// (No "rejects unsupported" cases remain — Phase 3 lifted every gate that
// Phase 1+2 had. phaseCount < 1 is still rejected as a validation error,
// not a "deferred to Phase X" gate.)


// =========================================================================
// Phase 3 — Item 6: phaseCount > 1 interleaved channels.
//
// N parallel boost stages per phase. Each channel sees I_pk/N average
// current and L_single/N inductance, so the per-channel ripple is N×
// the single-channel ripple but the SUM-current ripple at the input
// cancels to ~1/N (ideal). The adviser sizes ONE channel inductor;
// the user builds 3·N physical units.
//
// Gates:
//   (1) run_checks accepts phaseCount=N for N >= 1.
//   (2) computedBoostInductance with phaseCount=2 is half the
//       single-channel computed inductance (at the same operating
//       point and ripple-ratio target).
//   (3) process_operating_points emits 3·N windings per OP.
//   (4) Winding names tag the channel ("Phase A ch 1", ...).
// =========================================================================
TEST_CASE("Vienna: phaseCount=2 halves computed L and emits 6 windings",
          "[vienna-topology][phase3]") {
    auto jN1 = make_vienna_json();
    Vienna vN1(jN1);
    (void) vN1.process();
    const double L_single = vN1.get_computed_boost_inductance();

    auto jN2 = make_vienna_json();
    jN2["phaseCount"] = 2;
    Vienna vN2(jN2);
    REQUIRE(vN2.run_checks(false));  // (1)
    auto inputs = vN2.process();

    // (2)
    CHECK(vN2.get_computed_boost_inductance() ==
          Catch::Approx(L_single / 2.0).epsilon(1e-6));

    // (3) 1 input × 1 OP × 3 phases × 2 channels = 6 windings
    REQUIRE(inputs.get_operating_points().size() == 1);
    const auto& exc = inputs.get_operating_points()[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 6);

    // (4) names — first six should be Phase A ch 1, Phase A ch 2, Phase B
    //     ch 1, Phase B ch 2, Phase C ch 1, Phase C ch 2.
    REQUIRE(exc[0].get_name().has_value());
    CHECK(exc[0].get_name().value() == "Phase A ch 1");
    CHECK(exc[1].get_name().value() == "Phase A ch 2");
    CHECK(exc[2].get_name().value() == "Phase B ch 1");
    CHECK(exc[5].get_name().value() == "Phase C ch 2");
}


TEST_CASE("Vienna: phaseCount=3 per-channel current is 1/3 of single-channel",
          "[vienna-topology][phase3]") {
    auto jN1 = make_vienna_json();
    Vienna vN1(jN1);
    (void) vN1.process();
    double I_pk_phase = vN1.get_computed_line_peak_current();

    auto jN3 = make_vienna_json();
    jN3["phaseCount"] = 3;
    Vienna vN3(jN3);
    auto inputs = vN3.process();
    const auto& exc = inputs.get_operating_points()[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 9);  // 3 phases × 3 channels

    // Per-channel triangular current — at peakOfLineOnly (default), Phase A
    // ch 1 should have I_avg ≈ +I_pk_phase / 3 (since sin(π/2) = 1).
    auto offsetA1 = exc[0].get_current().value().get_processed().value().get_offset();
    CHECK(offsetA1 ==
          Catch::Approx(I_pk_phase / 3.0).epsilon(0.05));
}


// =========================================================================
// Phase 3 — Item 3: viennaII variant.
//
// Vienna II splits each leg into two back-to-back MOSFETs (bidirectional
// clamp), each conducting one polarity half-cycle. Per-switch RMS / avg
// are HALF their Vienna I counterparts; the fast-rectifier diodes are
// replaced by body diodes of the OFF-polarity switch.
//
// Gates:
//   (1) run_checks accepts viennaII.
//   (2) Closed-form check: I_sw_rms² for Vienna II = ½ · I_sw_rms² for
//       Vienna I at the same (I_pk, M).
//   (3) Inductor waveforms are UNCHANGED — the variant difference is on
//       the converter switch arrangement, not the boost inductor itself.
// =========================================================================
TEST_CASE("Vienna: viennaII switch RMS is half Vienna I in energy terms",
          "[vienna-topology][phase3]") {
    // (2) closed-form sanity check.
    const double I_pk = 10.0;
    const double M    = 0.85;
    double rmsI    = Vienna::compute_switch_rms(I_pk, M);
    double rmsII   = Vienna::compute_switch_rms_vienna_ii(I_pk, M);
    // I²_II / I²_I = 0.5 (exact, by integral construction).
    CHECK((rmsII * rmsII) ==
          Catch::Approx(0.5 * rmsI * rmsI).epsilon(1e-9));

    double avgII = Vienna::compute_switch_avg_vienna_ii(I_pk, M);
    // I_avg_ii = I_pk · (1/π − M/4) — quick spot check.
    CHECK(avgII == Catch::Approx(I_pk * (1.0 / M_PI - M / 4.0)).epsilon(1e-9));
}


// =========================================================================
// Phase 3 — Item 5: synchronousRectifier.
//
// The boost diode and its sync-rec MOSFET replacement carry identical
// current — the inductor doesn't care which rectifier device is in
// circuit. What changes is the caller's conduction-loss model:
//   diode  → P = Vf · I_avg
//   syncRec → P = Rds · I_rms²
//
// So sync rec primarily needs an I_rms diagnostic to multiply by Rds.
// The closed-form (Hartmann ETH 19755 §3.2.3, half-cycle integration):
//   I_avg = I_pk · M / 4
//   I_rms = I_pk · sqrt(2M / (3π))
// =========================================================================
TEST_CASE("Vienna: boost-diode avg/rms closed forms",
          "[vienna-topology][phase3]") {
    const double I_pk = 25.0;
    const double M    = 0.85;
    CHECK(Vienna::compute_boost_diode_avg(I_pk, M) ==
          Catch::Approx(I_pk * M / 4.0).epsilon(1e-9));
    CHECK(Vienna::compute_boost_diode_rms(I_pk, M) ==
          Catch::Approx(I_pk * std::sqrt(2.0 * M / (3.0 * M_PI))).epsilon(1e-9));
    // I_rms > I_avg always (RMS ≥ avg by Cauchy-Schwarz on positive signal).
    CHECK(Vienna::compute_boost_diode_rms(I_pk, M) >
          Vienna::compute_boost_diode_avg(I_pk, M));
}


TEST_CASE("Vienna: synchronousRectifier=true accepted + boost-diode RMS populated",
          "[vienna-topology][phase3]") {
    auto j = make_vienna_json();
    j["synchronousRectifier"] = true;
    Vienna v(j);
    REQUIRE(v.run_checks(false));
    auto inputs = v.process();
    REQUIRE(!inputs.get_operating_points().empty());

    double I_pk = v.get_computed_line_peak_current();
    double M    = v.get_last_modulation_index();
    REQUIRE(I_pk > 0.0);

    CHECK(v.get_last_boost_diode_avg_current() ==
          Catch::Approx(Vienna::compute_boost_diode_avg(I_pk, M)).epsilon(1e-6));
    CHECK(v.get_last_boost_diode_rms_current() ==
          Catch::Approx(Vienna::compute_boost_diode_rms(I_pk, M)).epsilon(1e-6));

    // Inductor waveform is unchanged from the diode-rectifier variant —
    // sync rec only changes the rectifier device, not the boost current.
    auto j2 = make_vienna_json();
    j2["synchronousRectifier"] = false;
    Vienna v2(j2);
    (void) v2.process();
    CHECK(v.get_last_inductor_peak_current() ==
          Catch::Approx(v2.get_last_inductor_peak_current()).epsilon(1e-6));
}


// =========================================================================
// Phase 3 — Item 4: alternative switch types.
//
// switchType=backToBackMosfet shares Vienna II's two-switches-per-leg
// arrangement, so the per-switch RMS / avg numbers must equal the
// viennaVariant=viennaII case at the same operating point.
//
// switchType=singleMosfetIn4DiodeBridge is single-switch-per-leg (like
// Vienna I tType) but with 4 forward diode drops in the conduction path
// — the per-switch RMS / avg are unchanged from Vienna I; the diode-loss
// model is the caller's responsibility.
// =========================================================================
TEST_CASE("Vienna: backToBackMosfet switchType matches viennaII diagnostics",
          "[vienna-topology][phase3]") {
    auto jBb = make_vienna_json();
    jBb["switchType"] = "backToBackMosfet";
    Vienna vBb(jBb);
    REQUIRE(vBb.run_checks(false));
    (void) vBb.process();

    auto jII = make_vienna_json();
    jII["viennaVariant"] = "viennaII";
    Vienna vII(jII);
    (void) vII.process();

    CHECK(vBb.get_last_switch_rms_current() ==
          Catch::Approx(vII.get_last_switch_rms_current()).epsilon(1e-9));
    CHECK(vBb.get_last_diode_avg_current() ==
          Catch::Approx(vII.get_last_diode_avg_current()).epsilon(1e-9));
}


TEST_CASE("Vienna: singleMosfetIn4DiodeBridge switchType keeps Vienna I diagnostics",
          "[vienna-topology][phase3]") {
    auto jSm = make_vienna_json();
    jSm["switchType"] = "singleMosfetIn4DiodeBridge";
    Vienna vSm(jSm);
    REQUIRE(vSm.run_checks(false));
    (void) vSm.process();

    auto jI = make_vienna_json();  // tType + viennaI (defaults)
    Vienna vI(jI);
    (void) vI.process();

    CHECK(vSm.get_last_switch_rms_current() ==
          Catch::Approx(vI.get_last_switch_rms_current()).epsilon(1e-9));
    CHECK(vSm.get_last_diode_avg_current() ==
          Catch::Approx(vI.get_last_diode_avg_current()).epsilon(1e-9));
}


TEST_CASE("Vienna: viennaII run_checks + diagnostics dispatch",
          "[vienna-topology][phase3]") {
    auto j = make_vienna_json();
    j["viennaVariant"] = "viennaII";
    Vienna v(j);

    // (1)
    REQUIRE(v.run_checks(false));

    auto inputs = v.process();
    REQUIRE(!inputs.get_operating_points().empty());

    // Diagnostics should reflect Vienna II's per-switch numbers, NOT
    // Vienna I's. Compute both and confirm we got the II flavour.
    double I_pk = v.get_computed_line_peak_current();
    double M    = v.get_last_modulation_index();
    REQUIRE(I_pk > 0.0);
    REQUIRE(M > 0.0);

    double expectedII = Vienna::compute_switch_rms_vienna_ii(I_pk, M);
    CHECK(v.get_last_switch_rms_current() ==
          Catch::Approx(expectedII).epsilon(1e-6));

    // (3) Inductor waveform is unchanged: Vienna I and Vienna II produce
    // the same inductor current envelope. Compare peak-current diagnostic
    // across the two variants at the same OP.
    auto j2 = make_vienna_json();
    j2["viennaVariant"] = "viennaI";
    Vienna vI(j2);
    (void) vI.process();
    CHECK(v.get_last_inductor_peak_current() ==
          Catch::Approx(vI.get_last_inductor_peak_current()).epsilon(1e-6));
}


// =========================================================================
// Phase 3 — peakOfLinePlusSectors sampling.
//
// Emits 6 OPs per input op (DPWM sector mid-points at 30°/90°/150°/210°/
// 270°/330°). Each OP carries 3 windings — at the same wall-clock instant
// each phase is at its own local sin·θ. So within a single sector OP the
// three windings are NOT identical (different |V_phase|, different duty,
// different ripple), and across the 6 sectors the per-phase stress varies
// from zero (sin near 0) to peak (sin = ±1).
//
// Gates:
//   (1) run_checks accepts peakOfLinePlusSectors.
//   (2) process_operating_points returns 6 OPs per input op.
//   (3) Each OP carries 3 windings (one per phase).
//   (4) The sector at θ=90° is identical to peakOfLineOnly's snapshot for
//       Phase A (sin(π/2)=1 — worst-case), validating the helper is the
//       same building block both paths share.
//   (5) At sector θ=30° (sin=0.5), Phase A's local |V|·d ripple is
//       roughly half of the peak-of-line value (sanity check on the
//       per-angle duty/ripple computation).
// =========================================================================
TEST_CASE("Vienna: peakOfLinePlusSectors emits 6 OPs with per-sector physics",
          "[vienna-topology][phase3]") {
    auto j = make_vienna_json();
    j["samplingStrategy"] = "peakOfLinePlusSectors";
    Vienna v(j);

    // (1)
    REQUIRE(v.run_checks(false));

    auto inputs = v.process();
    const auto& ops = inputs.get_operating_points();

    // (2) 1 input × 6 sectors = 6 OPs
    REQUIRE(ops.size() == 6);

    // (3)
    for (const auto& op : ops) {
        REQUIRE(op.get_excitations_per_winding().size() == 3);
    }

    // (4) Peak-of-line cross-check: sector @ 90° (index 1) — Phase A there
    //     equals what peakOfLineOnly emits, modulo the names/labels.
    auto j2 = make_vienna_json();
    j2["samplingStrategy"] = "peakOfLineOnly";
    Vienna v2(j2);
    auto inputs2 = v2.process();
    REQUIRE(inputs2.get_operating_points().size() == 1);
    const auto& peakOp     = inputs2.get_operating_points()[0];
    const auto& sectorOp90 = ops[1];

    // Average from the Processed.offset field — Inputs::calculate_processed_data
    // computes the DC offset robustly from the FFT (bin 0). Summing the raw
    // sampled-waveform data array directly is unreliable: power-of-2 resampling
    // of the 3-point triangular can land sample points unevenly across the
    // triangle's apex, biasing a naive arithmetic mean. The FFT DC term is
    // immune to that.
    auto get_iavg = [](const OperatingPoint& op, size_t w) {
        auto sig = op.get_excitations_per_winding()[w].get_current().value();
        REQUIRE(sig.get_processed().has_value());
        return sig.get_processed().value().get_offset();
    };
    double iPeakA   = get_iavg(peakOp,     0);
    double iSec90A  = get_iavg(sectorOp90, 0);
    CHECK(iSec90A == Catch::Approx(iPeakA).epsilon(0.05));

    // (5) Sector @ 30° (sin=0.5) — Phase A's |I_avg| ≈ 0.5·I_pk;
    //     Phase B at 30°-120° = -90° (sin=-1) — Phase B's |I_avg| ≈ I_pk
    //     with negative sign (the inductor is at its -peak).
    double Ipeak = v.get_last_inductor_peak_current();
    double Ipeak_envelope = v.get_computed_line_peak_current();
    REQUIRE(Ipeak_envelope > 0.0);

    double iSec30A = get_iavg(ops[0], 0);
    double iSec30B = get_iavg(ops[0], 1);
    double iSec30C = get_iavg(ops[0], 2);

    CHECK(std::abs(iSec30A) ==
          Catch::Approx(0.5 * Ipeak_envelope).margin(0.20 * Ipeak_envelope));
    CHECK(iSec30B < -0.5 * Ipeak_envelope);  // at its negative peak
    CHECK(std::abs(iSec30C) ==
          Catch::Approx(0.5 * Ipeak_envelope).margin(0.20 * Ipeak_envelope));

    (void) Ipeak;  // suppress unused warning if any
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
          "[vienna-topology][phase3]") {
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
          "[vienna-topology][phase3]") {
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


TEST_CASE("Vienna: missing operating points fails at construction", "[vienna-topology]") {
    // Vienna validates spec shape eagerly in the ctor (validate_vienna_spec_shape
    // at Vienna.cpp:86). Empty operatingPoints arrays are rejected before
    // run_checks ever runs, which is the desired fail-fast contract.
    auto j = make_vienna_json();
    j["operatingPoints"] = json::array();
    REQUIRE_THROWS_AS(Vienna(j), std::runtime_error);
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
