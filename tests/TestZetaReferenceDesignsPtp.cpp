// Zeta reference-design PtP regression tests — three published designs
// spanning the topology's realistic envelope.
//
// Mirrors TestSepicReferenceDesignsPtp.cpp / TestCukReferenceDesignsPtp.cpp —
// the canonical 4-gate PtP pattern documented in tests/PtpHelpers.h.
//
// Each test:
//   1. Builds a Zeta model with the ref-design specs (η=1, Vd=0 — lossless
//      analytical anchor).
//   2. Runs the analytical pass to anchor the expected L1-current shape.
//   3. Runs the ngspice switching netlist for one steady-state period
//      (after ≥400 settling periods to drain the 4th-order LC tank tail).
//   4. Asserts the four PtP gates:
//        Gate 1  walltime          < tol_walltime
//        Gate 2  load consistency  |Vout_spice/Iout_spice − Rload_nom| / Rload_nom
//                                  < tol_rload_pct/100
//                                  (Zeta: Vout positive, Iout positive ⇒
//                                   Rload positive)
//        Gate 3  power balance     0 ≤ (Pin − Pout)/Pin ≤ tol_loss_max
//        Gate 4  L1-current NRMSE  ptp_nrmse(analytical, SPICE) < tol_nrmse
//
// Zeta-specific notes:
//   * Vo is positive (non-inverting) — same sign convention as SEPIC,
//     opposite to Cuk. Iout drawn into the load is also positive in our
//     SPICE wiring (Vout_sense vout vout_load 0; Rload vout_load 0 R), so
//     Rload_spice = +Vout_mean / +Iout_mean is positive.
//   * Both inductors see +Vin during the ON sub-interval and −Vo during
//     OFF — same shape as SEPIC's L1 voltage. The L1-current shape is
//     therefore Vin/Vo-driven on each half-period.
//   * 4th-order LC tank — same settling characteristic as SEPIC and Cuk
//     (≥400 periods needed at all corners). The L2-Co output filter
//     yields LOW output ripple (vs SEPIC's pulsed-diode output), which
//     can speed convergence at the load port.
//
// Sources (ZETA_PLAN.md §3 reference-design table):
//   #1 TI PMP9581-style step-down (12 V → 5 V @ 1 A, 600 kHz) —
//      low-corner buck-mode regulation.
//   #2 TI LM5085 SNVA608-style buck-boost (12 V → 12 V @ 1 A, 300 kHz) —
//      mid-corner automotive 12 V rail regulation (D=0.5).
//   #3 TI PMP9581-style step-up (5 V → 12 V @ 0.5 A, 600 kHz) —
//      high-corner step-up (D ≈ 0.706, mirrors SEPIC SNVA168E spec).
//
// Tags: [converter-model][zeta-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "converter_models/Zeta.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;

namespace {

struct ZetaRefDesignSpec {
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

OpenMagnetics::Zeta build(const ZetaRefDesignSpec& s) {
    OpenMagnetics::Zeta c;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    c.set_input_voltage(iv);
    c.set_diode_voltage_drop(0.0);   // ideal-diode analytical reference
    c.set_efficiency(1.0);            // lossless anchor
    c.set_current_ripple_ratio(0.4);
    MAS::TopologyExcitation op;       // Zeta uses unified TopologyExcitation
    op.set_output_voltages({s.Vout});
    op.set_output_currents({s.Iout});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    c.set_operating_points({op});
    return c;
}

using namespace OpenMagnetics::Testing;

void run_ptp_gates(const ZetaRefDesignSpec& s) {
    std::cout << "\n========== Zeta PtP — " << s.name << " ==========\n";

    OpenMagnetics::NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    // ── Analytical pass ────────────────────────────────────────────────
    auto zeta = build(s);
    auto analyticalOps = zeta.process_operating_points(std::vector<double>{}, s.L1);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // ── SPICE switching pass ───────────────────────────────────────────
    zeta.set_num_steady_state_periods(400);
    zeta.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = zeta.simulate_and_extract_topology_waveforms(s.L1);
    auto simOps = zeta.simulate_and_extract_operating_points(s.L1);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency (Vout, Iout both positive in Zeta).
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

    // Gate 3 — power balance: Pin must dominate Pout. Both positive in Zeta.
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

TEST_CASE("Zeta reference design PtP — TI PMP9581 step-down (5 W, 600 kHz)",
          "[converter-model][zeta-topology][refdesign][ptp][slow]") {
    // TI PMP9581-style Zeta step-down:
    //   Vin = 12 V, Vo = 5 V, Iout = 1 A, fsw = 600 kHz.
    // D ≈ 5/17 ≈ 0.294 (lossless). Buck-mode Zeta, low-corner.
    ZetaRefDesignSpec s{
        /*name*/          "TI PMP9581 step-down",
        /*Vin*/           12.0,
        /*Vout*/          5.0,
        /*Iout*/          1.0,
        /*Fs*/            600e3,
        /*L1*/            22e-6,
        /*tol_walltime*/  6.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_loss_max*/  0.65,
        /*tol_nrmse*/     0.10
    };
    run_ptp_gates(s);
}

TEST_CASE("Zeta reference design PtP — TI LM5085 SNVA608 (12 W, 300 kHz)",
          "[converter-model][zeta-topology][refdesign][ptp][slow]") {
    // TI LM5085 SNVA608-style Zeta buck-boost (automotive 12 V rail):
    //   Vin = 12 V, Vo = 12 V, Iout = 1 A, fsw = 300 kHz.
    // D = 0.5 — buck-boost regulation point. Mid-corner.
    ZetaRefDesignSpec s{
        /*name*/          "TI LM5085 SNVA608",
        /*Vin*/           12.0,
        /*Vout*/          12.0,
        /*Iout*/          1.0,
        /*Fs*/            300e3,
        /*L1*/            33e-6,
        /*tol_walltime*/  4.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_loss_max*/  0.65,
        /*tol_nrmse*/     0.10
    };
    run_ptp_gates(s);
}

TEST_CASE("Zeta reference design PtP — Step-up 5→12V (6 W, 600 kHz)",
          "[converter-model][zeta-topology][refdesign][ptp][slow]") {
    // Zeta high-corner step-up (mirror of SEPIC SNVA168E spec to exercise
    // the same D ≈ 0.706 regime through the Zeta sub-interval physics):
    //   Vin = 5 V, Vo = 12 V, Iout = 0.5 A, fsw = 600 kHz.
    // D ≈ 12/17 ≈ 0.706. High-corner — large IL1 bias = Iout·D/(1-D) ≈ 1.2 A.
    ZetaRefDesignSpec s{
        /*name*/          "Zeta step-up 5V→12V",
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
