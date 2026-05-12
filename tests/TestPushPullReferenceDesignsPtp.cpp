// Phase 8j — Golden PtP regression tests for three Push-Pull reference
// designs (TI SN6501 1.2 W, SN6505B 3.3 W, SN6507 5 W).
//
// Mirrors the canonical Buck/Boost/ACF/Flyback PtP harness pattern (4 gates):
//
//   wall_time            < tol_walltime
//   load consistency     |Vout_spice/Iout_spice − Rload_nom|/Rload_nom
//                        < tol_rload_pct
//   power balance        0 ≤ (Pin − Pout_spice)/Pin ≤ tol_loss_max
//   primary-current NRMSE analytical vs SPICE < tol_nrmse
//
// Specs reuse the three reference designs already characterised in
// TestTopologyPushPull.cpp (kPushPullRefDesign1..3); duplicated here as
// struct literals so each PtP test is self-describing.
//
// Tags: [converter-model][push-pull-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "converter_models/PushPull.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

using namespace OpenMagnetics::Testing;

struct RefDesignSpec {
    const char* name;
    double Vin;             // input voltage [V]
    double Vout;            // output voltage [V]   (nominal — load wiring)
    double Iout;            // output current [A]   (nominal — load wiring)
    double Fs;              // switching frequency [Hz]
    double N;               // turns ratio Np:Ns
    double ripple;          // currentRippleRatio (CCM)
    double tol_walltime;    // wall-time gate [s]
    double tol_rload_pct;   // |Vout_spice/Iout_spice − Rload_nom|/Rload_nom [%]
    double tol_loss_max;    // (Pin − Pout_spice)/Pin upper bound [fraction]
    double tol_nrmse;       // primary-current NRMSE [fraction, 0..1]
};

OpenMagnetics::PushPull build(const RefDesignSpec& s) {
    OpenMagnetics::PushPull pp;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    pp.set_input_voltage(iv);
    pp.set_diode_voltage_drop(0.0);   // lossless analytical reference
    pp.set_efficiency(1.0);
    pp.set_current_ripple_ratio(s.ripple);
    PushPullOperatingPoint op;
    op.set_output_voltages({s.Vout});
    op.set_output_currents({s.Iout});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    pp.set_operating_points(std::vector<PushPullOperatingPoint>{op});
    return pp;
}

// Lm chosen so the magnetizing peak-to-peak equals ripple·Iout/N — keeps
// SPICE ripple consistent with the analytical currentRippleRatio:
//   ΔImag = Vin·t1/Lm = ripple·Iout/N,
//   t1 = (Vout·N) / (2·Vin·Fs)
//   ⇒ Lm = Vout·N² / (2·Fs·ripple·Iout).
double consistent_lm(const RefDesignSpec& s) {
    return s.Vout * s.N * s.N / (2.0 * s.Fs * s.ripple * s.Iout);
}

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== PushPull PtP — " << s.name << " ==========\n";

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    auto pp = build(s);
    const double Lm = consistent_lm(s);
    // turnsRatios layout: [secondPrimary=1, mainSec=N, ...]
    const std::vector<double> tr{1.0, s.N, s.N};

    // ── Analytical pass — anchors primary-current shape ───────────────
    auto analyticalOps = pp.process_operating_points(tr, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // ── SPICE switching pass ──────────────────────────────────────────
    pp.set_num_steady_state_periods(500);
    pp.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = pp.simulate_and_extract_topology_waveforms(tr, Lm, 1);
    auto simOps = pp.simulate_and_extract_operating_points(tr, Lm);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency.
    const auto& voutW = wfs[0].get_output_voltages()[0];
    const auto& ioutW = wfs[0].get_output_currents()[0];
    const auto& voutD = voutW.get_data();
    const auto& ioutD = ioutW.get_data();
    const auto voutT_opt = voutW.get_time();
    const auto ioutT_opt = ioutW.get_time();
    REQUIRE(!voutD.empty());
    REQUIRE(!ioutD.empty());
    REQUIRE(voutT_opt.has_value());
    REQUIRE(ioutT_opt.has_value());
    const auto& voutT = voutT_opt.value();
    const auto& ioutT = ioutT_opt.value();
    const double voutMean = ConverterPortChecks::time_weighted_mean(voutT, voutD);
    const double ioutMean = ConverterPortChecks::time_weighted_mean(ioutT, ioutD);
    const double rloadNom = s.Vout / s.Iout;
    const double rloadSpc = voutMean / ioutMean;
    const double rloadErr = (rloadSpc - rloadNom) / rloadNom;
    std::cout << "  Vout_spice = " << voutMean << " V, "
              << "Iout_spice = " << ioutMean << " A, "
              << "Rload_spice = " << rloadSpc << " Ω "
              << "(Rload_nom = " << rloadNom << " Ω, err "
              << 100.0 * rloadErr << " %, gate "
              << s.tol_rload_pct << " %)\n";
    CHECK(std::fabs(rloadErr) < s.tol_rload_pct / 100.0);

    // Gate 3 — power balance.
    const auto& vinW = wfs[0].get_input_voltage();
    const auto& iinW = wfs[0].get_input_current();
    const auto& vinD = vinW.get_data();
    const auto& iinD = iinW.get_data();
    const auto vinT_opt = vinW.get_time();
    REQUIRE(!vinD.empty());
    REQUIRE(!iinD.empty());
    REQUIRE(vinT_opt.has_value());
    const auto& vinT = vinT_opt.value();
    const double pin  = ConverterPortChecks::time_weighted_mean_product(
        vinT, vinD, iinD);
    const double pout = ConverterPortChecks::time_weighted_mean_product(
        voutT, voutD, ioutD);
    const double loss = (pin - pout) / pin;
    std::cout << "  Pin = " << pin << " W, Pout = " << pout << " W, "
              << "loss = " << 100.0 * loss << " %   (gate 0..."
              << 100.0 * s.tol_loss_max << " %)\n";
    CHECK(loss >= 0.0);
    CHECK(loss <= s.tol_loss_max);

    // Gate 4 — primary-current shape NRMSE.
    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty());
    REQUIRE(!sTime.empty());
    auto sResampled = ptp_interp(sTime, sData, 256);
    const double nrmse = ptp_nrmse(aResampled, sResampled);
    std::cout << "  iPri NRMSE (analytical vs SPICE) = "
              << 100.0 * nrmse << " %   (gate "
              << 100.0 * s.tol_nrmse << " %)\n";
    CHECK(nrmse < s.tol_nrmse);
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────
// Per-design TEST_CASEs.  Specs match kPushPullRefDesign<1,2,3> in
// TestTopologyPushPull.cpp (canonical source).
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("PushPull reference design PtP — TI SN6501 (1.2 W, 5 V→3.3 V)",
          "[converter-model][push-pull-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "SN6501",
        /*Vin*/          5.0,
        /*Vout*/         3.3,
        /*Iout*/         0.35,
        /*Fs*/           410e3,
        /*N*/            1.0,
        /*ripple*/       0.4,
        /*tol_walltime*/ 6.0,
        /*tol_rload_pct*/ 0.5,
        /*tol_loss_max*/ 0.40,
        /*tol_nrmse*/    0.30
    };
    run_ptp_gates(s);
}

TEST_CASE("PushPull reference design PtP — TI SN6505B (3.3 W, 5 V→3.3 V)",
          "[converter-model][push-pull-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "SN6505B",
        /*Vin*/          5.0,
        /*Vout*/         3.3,
        /*Iout*/         1.0,
        /*Fs*/           420e3,
        /*N*/            1.0,
        /*ripple*/       0.4,
        /*tol_walltime*/ 6.0,
        /*tol_rload_pct*/ 0.5,
        /*tol_loss_max*/ 0.40,
        /*tol_nrmse*/    0.30
    };
    run_ptp_gates(s);
}

TEST_CASE("PushPull reference design PtP — TI SN6507 (5 W, 12 V→5 V)",
          "[converter-model][push-pull-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "SN6507",
        /*Vin*/          12.0,
        /*Vout*/         5.0,
        /*Iout*/         1.0,
        /*Fs*/           200e3,
        /*N*/            1.0,
        /*ripple*/       0.4,
        /*tol_walltime*/ 6.0,
        /*tol_rload_pct*/ 0.5,
        /*tol_loss_max*/ 0.40,
        /*tol_nrmse*/    0.30
    };
    run_ptp_gates(s);
}
