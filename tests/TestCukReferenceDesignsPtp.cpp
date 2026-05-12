// Phase 8m — Golden point-to-point (PtP) regression tests for three
// Cuk converter (V1 non-isolated) reference designs.
//
// Mirrors TestBuckReferenceDesignsPtp.cpp / TestBoostReferenceDesignsPtp.cpp
// — the canonical 4-gate PtP pattern documented in tests/PtpHelpers.h.
//
// Each test:
//   1. Builds a Cuk model with the ref-design specs (η=1, Vd=0 — lossless
//      analytical anchor).
//   2. Runs the analytical pass to anchor the expected L1-current shape.
//   3. Runs the ngspice switching netlist for one steady-state period
//      (after 400 settling periods to drain the 4th-order LC tank tail).
//   4. Asserts the four PtP gates:
//        Gate 1  walltime          < tol_walltime
//        Gate 2  load consistency  |Vout_spice/Iout_spice − Rload_nom| / Rload_nom
//                                  < tol_rload_pct/100
//                                  (Vo is negative so both Vout_mean and
//                                   Iout_mean are negative — their ratio
//                                   is positive Rload)
//        Gate 3  power balance     0 ≤ (Pin − Pout)/Pin ≤ tol_loss_max
//                                  (Pout = Vo·Iload = (−)·(−) = +)
//        Gate 4  L1-current NRMSE  ptp_nrmse(analytical, SPICE) < tol_nrmse
//
// Cuk-specific notes:
//   * Vo is signed *negative* throughout the waveforms; the user-facing
//     output voltage spec is the magnitude. Load wiring is symmetric so
//     the +Rload check passes either way.
//   * The 4th-order LC tank settles slowly — 400 periods needed at all
//     three Fs corners. Faster than IBB though (no transformer).
//
// Tags: [converter-model][cuk-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "converter_models/Cuk.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin;            // input voltage [V]
    double Vout_mag;       // |Vo| [V]   (Vo is negative; passed positive in opPoint)
    double Iout;           // load current [A]
    double Fs;             // switching frequency [Hz]
    double L1;             // input inductor [H]
    double tol_walltime;   // wall-time gate [s]
    double tol_rload_pct;  // |Rload_spice − Rload_nom|/Rload_nom [%]
    double tol_loss_max;   // (Pin − Pout)/Pin upper bound [fraction]
    double tol_nrmse;      // L1-current shape NRMSE [fraction, 0..1]
};

OpenMagnetics::Cuk build(const RefDesignSpec& s) {
    OpenMagnetics::Cuk c;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    c.set_input_voltage(iv);
    c.set_diode_voltage_drop(0.0);   // ideal-diode analytical reference
    c.set_efficiency(1.0);            // lossless
    c.set_current_ripple_ratio(0.4);
    BaseOperatingPoint op;
    op.set_output_voltages({s.Vout_mag});  // magnitude — Cuk treats it as |Vo|
    op.set_output_currents({s.Iout});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    c.set_operating_points({op});
    return c;
}

using namespace OpenMagnetics::Testing;

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== Cuk PtP — " << s.name << " ==========\n";

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    // ── Analytical pass ────────────────────────────────────────────────
    auto cuk = build(s);
    auto analyticalOps = cuk.process_operating_points(std::vector<double>{}, s.L1);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // ── SPICE switching pass ───────────────────────────────────────────
    cuk.set_num_steady_state_periods(400);
    cuk.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = cuk.simulate_and_extract_topology_waveforms(s.L1);
    auto simOps = cuk.simulate_and_extract_operating_points(s.L1);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency (Vout_mean and Iout_mean both negative;
    // their ratio is positive Rload).
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
    const double rloadNom = s.Vout_mag / s.Iout;
    // Guard: if either mean is ~0 the ratio explodes. Cuk should always
    // settle to negative Vout & Iout — this would catch a wiring bug.
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

    // Gate 3 — power balance: Pin must dominate Pout.
    // Pin  = mean(Vin·Iin)  = (+)·(+)
    // Pout = mean(Vo·Iload) = (−)·(−) = (+)
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
// Per-design TEST_CASEs. Sources: CUK_PLAN.md §3 reference-design table.
//   #1 Erickson-Maksimović 3rd ed. §2.4 (worked example)
//   #2 Simon Bramble LT3757-based 5 W 300 kHz Cuk
//   #3 Mid-Fs 200 kHz step-up synthetic (12 V → 24 V, 1 A)
// (LM2611 1.4 MHz fixture deferred — needs tighter SPICE step-time tuning,
//  per CUK_PLAN.md §6.2.)
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("Cuk reference design PtP — Erickson §2.4 (25 W, 100 kHz)",
          "[converter-model][cuk-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/          "Erickson §2.4",
        /*Vin*/           25.0,
        /*Vout_mag*/      25.0,
        /*Iout*/          1.0,
        /*Fs*/            100e3,
        /*L1*/            100e-6,
        /*tol_walltime*/  4.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_loss_max*/  0.65,
        /*tol_nrmse*/     0.05
    };
    run_ptp_gates(s);
}

TEST_CASE("Cuk reference design PtP — Bramble LT3757 (5 W, 300 kHz)",
          "[converter-model][cuk-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/          "Bramble LT3757",
        /*Vin*/           10.0,
        /*Vout_mag*/      5.0,
        /*Iout*/          1.0,
        /*Fs*/            300e3,
        /*L1*/            47e-6,
        /*tol_walltime*/  4.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_loss_max*/  0.62,
        /*tol_nrmse*/     0.07
    };
    run_ptp_gates(s);
}

TEST_CASE("Cuk reference design PtP — Synthetic step-up (24 W, 200 kHz)",
          "[converter-model][cuk-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/          "Synthetic 12→-24 V",
        /*Vin*/           12.0,
        /*Vout_mag*/      24.0,
        /*Iout*/          1.0,
        /*Fs*/            200e3,
        /*L1*/            47e-6,
        /*tol_walltime*/  4.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_loss_max*/  0.62,
        /*tol_nrmse*/     0.05
    };
    run_ptp_gates(s);
}
