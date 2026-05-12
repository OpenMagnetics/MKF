// SEPIC reference-design PtP regression tests — three published designs
// spanning the topology's realistic envelope.
//
// Mirrors TestCukReferenceDesignsPtp.cpp / TestBuckReferenceDesignsPtp.cpp —
// the canonical 4-gate PtP pattern documented in tests/PtpHelpers.h.
//
// Each test:
//   1. Builds a SEPIC model with the ref-design specs (η=1, Vd=0 — lossless
//      analytical anchor).
//   2. Runs the analytical pass to anchor the expected L1-current shape.
//   3. Runs the ngspice switching netlist for one steady-state period
//      (after ≥400 settling periods to drain the 4th-order LC tank tail).
//   4. Asserts the four PtP gates:
//        Gate 1  walltime          < tol_walltime
//        Gate 2  load consistency  |Vout_spice/Iout_spice − Rload_nom| / Rload_nom
//                                  < tol_rload_pct/100
//                                  (SEPIC: Vout positive, Iout positive ⇒
//                                   Rload positive)
//        Gate 3  power balance     0 ≤ (Pin − Pout)/Pin ≤ tol_loss_max
//        Gate 4  L1-current NRMSE  ptp_nrmse(analytical, SPICE) < tol_nrmse
//
// SEPIC-specific notes:
//   * Vo is positive (non-inverting) — opposite sign convention to Cuk.
//     Iout drawn into the load is also positive in our SPICE wiring
//     (Vout_sense vout vout_load 0; Rload vout_load 0 R), so
//     Rload_spice = +Vout_mean / +Iout_mean is positive.
//   * Both inductors see ±Vin during the ON sub-interval (NOT ±Vo as in
//     Cuk's L2). The L1-current shape is therefore Vin-driven on both
//     half-periods.
//   * 4th-order LC tank — same settling characteristic as Cuk (≥400
//     periods needed at all corners).
//
// Sources (SEPIC_PLAN.md §3 reference-design table):
//   #1 TI SNVA168E §6 worked example (5 V → 12 V @ 0.5 A, 600 kHz) —
//      low-corner step-up.
//   #2 ADI LTC1871-style SEPIC (3.3 V → 5 V @ 1 A, 250 kHz) — mid-corner.
//   #3 TI TIDA-00781 12 W reference design (12 V → 12 V @ 1 A, 250 kHz) —
//      high-corner buck-boost (D ≈ 0.5).
//
// Tags: [converter-model][sepic-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "converter_models/Sepic.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;

namespace {

struct SepicRefDesignSpec {
    const char* name;
    double Vin;            // input voltage [V]
    double Vout;           // output voltage [V]   (non-inverting, positive)
    double Iout;           // load current [A]
    double Fs;             // switching frequency [Hz]
    double L1;             // input inductor [H]
    double tol_walltime;   // wall-time gate [s]
    double tol_rload_pct;  // |Rload_spice − Rload_nom|/Rload_nom [%]
    double tol_loss_max;   // (Pin − Pout)/Pin upper bound [fraction]
    double tol_nrmse;      // L1-current shape NRMSE [fraction, 0..1]
};

OpenMagnetics::Sepic build(const SepicRefDesignSpec& s) {
    OpenMagnetics::Sepic c;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    c.set_input_voltage(iv);
    c.set_diode_voltage_drop(0.0);   // ideal-diode analytical reference
    c.set_efficiency(1.0);            // lossless anchor
    c.set_current_ripple_ratio(0.4);
    MAS::TopologyExcitation op;       // SEPIC uses unified TopologyExcitation
    op.set_output_voltages({s.Vout});
    op.set_output_currents({s.Iout});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    c.set_operating_points({op});
    return c;
}

using namespace OpenMagnetics::Testing;

void run_ptp_gates(const SepicRefDesignSpec& s) {
    std::cout << "\n========== SEPIC PtP — " << s.name << " ==========\n";

    OpenMagnetics::NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    // ── Analytical pass ────────────────────────────────────────────────
    auto sepic = build(s);
    auto analyticalOps = sepic.process_operating_points(std::vector<double>{}, s.L1);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // ── SPICE switching pass ───────────────────────────────────────────
    sepic.set_num_steady_state_periods(400);
    sepic.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = sepic.simulate_and_extract_topology_waveforms(s.L1);
    auto simOps = sepic.simulate_and_extract_operating_points(s.L1);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency (Vout, Iout both positive in SEPIC).
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
    REQUIRE(std::fabs(ioutMean) > 1e-9);
    const double rloadSpc = voutMean / ioutMean;
    const double rloadErr = (rloadSpc - rloadNom) / rloadNom;
    std::cout << "  Vout_spice = " << voutMean << " V, "
              << "Iout_spice = " << ioutMean << " A, "
              << "Rload_spice = " << rloadSpc << " Ω "
              << "(Rload_nom = " << rloadNom << " Ω, err "
              << 100.0 * rloadErr << " %, gate "
              << s.tol_rload_pct << " %)\n";
    CHECK(std::fabs(rloadErr) < s.tol_rload_pct / 100.0);

    // Gate 3 — power balance: Pin must dominate Pout. Both positive in SEPIC.
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
    CHECK(loss >= -0.05);   // small negative tolerance for numerical noise
    CHECK(loss <= s.tol_loss_max);

    // Gate 4 — L1-current waveform NRMSE.
    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty());
    REQUIRE(!sTime.empty());
    auto sResampled = ptp_interp(sTime, sData, 256);
    const double nrmse = ptp_nrmse(aResampled, sResampled);
    std::cout << "  iL1 NRMSE (analytical vs SPICE) = "
              << 100.0 * nrmse << " %   (gate "
              << 100.0 * s.tol_nrmse << " %)\n";
    CHECK(nrmse < s.tol_nrmse);
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────
// Per-design TEST_CASEs.
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("SEPIC reference design PtP — TI SNVA168E (6 W, 600 kHz)",
          "[converter-model][sepic-topology][refdesign][ptp][slow]") {
    // TI SNVA168E "Designing A SEPIC Converter" §6 worked example:
    //   Vin = 5 V, Vo = 12 V, Iout = 0.5 A, fsw = 600 kHz.
    // D ≈ 12/17 ≈ 0.706 (lossless). Step-up SEPIC, low-corner.
    SepicRefDesignSpec s{
        /*name*/          "TI SNVA168E",
        /*Vin*/           5.0,
        /*Vout*/          12.0,
        /*Iout*/          0.5,
        /*Fs*/            600e3,
        /*L1*/            22e-6,
        /*tol_walltime*/  6.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_loss_max*/  0.65,
        /*tol_nrmse*/     0.10
    };
    run_ptp_gates(s);
}

TEST_CASE("SEPIC reference design PtP — ADI LTC1871 (5 W, 250 kHz)",
          "[converter-model][sepic-topology][refdesign][ptp][slow]") {
    // ADI LTC1871-style SEPIC controller reference:
    //   Vin = 3.3 V, Vo = 5 V, Iout = 1 A, fsw = 250 kHz.
    // D ≈ 5/8.3 ≈ 0.602. Mid-corner step-up.
    SepicRefDesignSpec s{
        /*name*/          "ADI LTC1871",
        /*Vin*/           3.3,
        /*Vout*/          5.0,
        /*Iout*/          1.0,
        /*Fs*/            250e3,
        /*L1*/            10e-6,
        /*tol_walltime*/  4.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_loss_max*/  0.65,
        /*tol_nrmse*/     0.10
    };
    run_ptp_gates(s);
}

TEST_CASE("SEPIC reference design PtP — TI TIDA-00781 (12 W, 250 kHz)",
          "[converter-model][sepic-topology][refdesign][ptp][slow]") {
    // TI TIDA-00781 12 W SEPIC reference design:
    //   Vin = 12 V (nominal), Vo = 12 V, Iout = 1 A, fsw = 250 kHz.
    // D = 0.5 — buck-boost regulation point. High-corner test
    // since the converter must hold regulation across Vin droop / surge.
    SepicRefDesignSpec s{
        /*name*/          "TI TIDA-00781",
        /*Vin*/           12.0,
        /*Vout*/          12.0,
        /*Iout*/          1.0,
        /*Fs*/            250e3,
        /*L1*/            22e-6,
        /*tol_walltime*/  4.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_loss_max*/  0.65,
        /*tol_nrmse*/     0.10
    };
    run_ptp_gates(s);
}
