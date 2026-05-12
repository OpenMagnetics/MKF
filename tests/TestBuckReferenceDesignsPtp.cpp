// Phase 7 — Golden point-to-point (PtP) regression tests for the three
// Buck reference designs (TPS54202EVM-716 10 W, LMR33630ADDAEVM 15 W,
// LM5146-Q1-EVM12V 96 W).
//
// This file mirrors `TestPfcReferenceDesignsPtp.cpp` for the boost-PFC
// trio — same "single-purpose golden PtP" structure, but applied to the
// DC-DC Buck converter.  The per-design Values + raw analytical-vs-SPICE
// NRMSE checks already live in `TestTopologyBuck.cpp` (the
// `Test_Buck_RefDesign<N>_*` cases).  This file consolidates the four
// acceptance gates that file is missing — wall-time, DC-port regulation,
// power balance, and waveform NRMSE — under the canonical
// [refdesign][ptp][slow] tag set so a single suite filter exercises the
// full per-design contract for every topology.
//
// Each test:
//   1. Builds a Buck model with the ref-design specs (synchronous, η=1,
//      Vd=0 — the actual EVM topology, lossless analytical reference).
//   2. Runs the analytical pass to anchor the inductor-current shape.
//   3. Runs the ngspice switching netlist for one steady-state period
//      (after 400 settling periods to drain the Cout RC tail).
//   4. Asserts the four PtP gates:
//        wall-time         < 30 s          (Buck designs typically <1 s)
//        load consistency  |Vout_spice/Iout_spice − Rload_nom|/Rload_nom
//                          within tol_rload_pct  (verifies resistive load
//                          is wired to Vout_nom/Iout_nom; the open-loop
//                          power stage settles wherever conduction loss
//                          puts it, but the *ratio* must hold)
//        power balance     0 ≤ (Pin − Pout_spice)/Pin ≤ tol_loss_max
//                          (Pin must dominate Pout — η ∈ [1 − tol_loss_max, 1])
//        inductor NRMSE    analytical vs SPICE shape NRMSE < 15 %
//
//   Gates 2 and 3 are intentionally SPICE-internal (self-consistent),
//   not vs nominal Vout/Iout/Pout: the Buck SPICE netlist is open-loop
//   (no Vout regulator), so Vout settles below nominal in proportion to
//   the conduction/switching loss.  Gating against nominal would make
//   the test a duty-cycle pre-compensation exercise rather than a model-
//   faithfulness check.  The PFC PtP analog gates against nominal
//   because PFC SPICE includes a closed-loop Vbus controller.
//
// Tags: [converter-model][buck-topology][refdesign][ptp][slow]
//   The [slow] tag matches the PFC PtP convention.  Buck PtP is actually
//   sub-second per design — well under the 30 s wall-time gate — but the
//   tag keeps the tri-topology PtP suite uniformly filterable.

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "converter_models/Buck.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin;            // input voltage [V]
    double Vout;           // output voltage [V]   (nominal — load wiring)
    double Iout;           // output current [A]   (nominal — load wiring)
    double Fs;             // switching frequency [Hz]
    double L;              // inductor value [H]
    double tol_walltime;   // wall-time gate [s]
    double tol_rload_pct;  // |Vout_spice/Iout_spice − Rload_nom|/Rload_nom [%]
    double tol_loss_max;   // (Pin − Pout_spice)/Pin upper bound [fraction]
    double tol_nrmse;      // inductor-current shape NRMSE [fraction, 0..1]
};

OpenMagnetics::Buck build(const RefDesignSpec& s) {
    OpenMagnetics::Buck b;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    b.set_input_voltage(iv);
    b.set_diode_voltage_drop(0.0);   // synchronous buck (matches EVMs)
    b.set_efficiency(1.0);            // lossless analytical reference
    b.set_current_ripple_ratio(0.4);
    BaseOperatingPoint op;
    op.set_output_voltages({s.Vout});
    op.set_output_currents({s.Iout});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    b.set_operating_points({op});
    return b;
}

using namespace OpenMagnetics::Testing;

// ── PtP harness ────────────────────────────────────────────────────────
void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== Buck PtP — " << s.name << " ==========\n";

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    // ── Analytical pass (anchors the expected inductor-current shape) ─
    auto buck = build(s);
    auto analyticalOps = buck.process_operating_points(std::vector<double>{}, s.L);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // ── SPICE switching pass — the wall-clock-gated workhorse ─────────
    // The Buck class hardcodes a 100 µF Cout, forming an RC tail of
    // τ = (Vout/Iout) · 100 µF.  At the high-Fs / low-Iout corners the
    // default 50 settling periods is < 8·τ; bump to 400 so the
    // extraction window opens after the load step has fully decayed.
    buck.set_num_steady_state_periods(400);
    buck.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = buck.simulate_and_extract_topology_waveforms(s.L);
    auto simOps = buck.simulate_and_extract_operating_points(s.L);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency.  The Buck SPICE netlist is open-loop
    // (no Vout regulator); Vout settles below nominal in proportion to
    // the conduction/switching loss.  But the resistive load is wired
    // for Rload_nom = Vout_nom / Iout_nom, so Vout_spice / Iout_spice
    // must reproduce that ratio.  This catches load-wiring bugs without
    // making the test sensitive to the absolute loss budget.
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

    // Gate 3 — power balance: Pin must dominate Pout (lossy but causal).
    // Compare SPICE-internal Pin vs Pout (NOT vs nominal Pout, which is
    // unreachable open-loop).  η = Pout / Pin must lie in
    // [1 − tol_loss_max, 1].  Negative loss (Pin < Pout) indicates a
    // sign error or extraction-window bug.
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

    // Gate 4 — inductor-current waveform NRMSE (analytical vs SPICE).
    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty());
    REQUIRE(!sTime.empty());
    auto sResampled = ptp_interp(sTime, sData, 256);
    const double nrmse = ptp_nrmse(aResampled, sResampled);
    std::cout << "  iL NRMSE (analytical vs SPICE) = "
              << 100.0 * nrmse << " %   (gate "
              << 100.0 * s.tol_nrmse << " %)\n";
    CHECK(nrmse < s.tol_nrmse);
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────
// Per-design TEST_CASEs.  Specs match `kBuckRefDesign<1,2,3>` in
// TestTopologyBuck.cpp (canonical source); duplicated here as struct
// literals so each PtP test is self-describing.
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("Buck reference design PtP — TPS54202EVM-716 (10 W)",
          "[converter-model][buck-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "TPS54202EVM-716",
        /*Vin*/          12.0,
        /*Vout*/         5.0,
        /*Iout*/         2.0,
        /*Fs*/           500e3,
        /*L*/            22e-6,
        /*tol_walltime*/ 30.0,
        /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.50,
        /*tol_nrmse*/    0.06
    };
    run_ptp_gates(s);
}

TEST_CASE("Buck reference design PtP — LMR33630ADDAEVM (15 W)",
          "[converter-model][buck-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "LMR33630ADDAEVM",
        /*Vin*/          12.0,
        /*Vout*/         5.0,
        /*Iout*/         3.0,
        /*Fs*/           400e3,
        /*L*/            10e-6,
        /*tol_walltime*/ 30.0,
        /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.50,
        /*tol_nrmse*/    0.06
    };
    run_ptp_gates(s);
}

TEST_CASE("Buck reference design PtP — LM5146-Q1-EVM12V (96 W)",
          "[converter-model][buck-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "LM5146-Q1-EVM12V",
        /*Vin*/          24.0,
        /*Vout*/         12.0,
        /*Iout*/         8.0,
        /*Fs*/           400e3,
        /*L*/            22e-6,
        /*tol_walltime*/ 30.0,
        /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.50,
        /*tol_nrmse*/    0.06
    };
    run_ptp_gates(s);
}
