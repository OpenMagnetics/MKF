// Phase 8i — Golden PtP regression tests for three Flyback reference
// designs (TI PMP30817 1.2 W, TI LM5180EVM-DUAL 3 W, TI TIDA-00709 33 W).
//
// Mirrors the canonical Buck/Boost/ACF PtP harness pattern (4 gates):
//
//   wall_time            < tol_walltime
//   load consistency     |Vout_spice/Iout_spice − Rload_nom|/Rload_nom
//                        < tol_rload_pct
//   power balance        0 ≤ (Pin − Pout_spice)/Pin ≤ tol_loss_max
//   primary-current NRMSE analytical vs SPICE < tol_nrmse
//
// Specs reuse the three reference designs already characterised in
// TestTopologyFlyback.cpp (kFlybackRefDesign1..3); duplicated here as
// struct literals so each PtP test is self-describing.
//
// Tags: [converter-model][flyback-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "converter_models/Flyback.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"
#include "WaveformDumpHelpers.h"

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

OpenMagnetics::Flyback build(const RefDesignSpec& s) {
    OpenMagnetics::Flyback f;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    f.set_input_voltage(iv);
    f.set_diode_voltage_drop(0.0);   // lossless analytical reference
    f.set_efficiency(1.0);
    f.set_current_ripple_ratio(s.ripple);
    f.set_maximum_duty_cycle(0.5);
    f.set_maximum_drain_source_voltage(s.Vin + s.N * s.Vout + 50.0);
    OpenMagnetics::FlybackOperatingPoint op;
    op.set_output_voltages({s.Vout});
    op.set_output_currents({s.Iout});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    f.set_operating_points(std::vector<OpenMagnetics::FlybackOperatingPoint>{op});
    return f;
}

// Build a Flyback identical to `build(s)` but with output current
// overridden to the SPICE-settled value (`Iout_spice`) and the current
// ripple ratio explicitly set to its physics value
//   ΔI_pp / (2·Iavg)  with  ΔI_pp = Vin·D/(Lm·Fs)
// rather than the nominal-design 0.4. Used for the NRMSE shape
// comparison so the analytical and SPICE primary-current waveforms
// reference the same operating point: open-loop SPICE settles Vout
// below nominal due to diode/switch losses, which lowers the actual
// primary-current offset (Iavg = Iout/N/(1-D)) and raises the actual
// ripple-ratio (ΔI/Iavg) — both effects need to be reflected in the
// analytical model before the shape comparison is meaningful.
OpenMagnetics::Flyback build_for_shape_match(const RefDesignSpec& s,
                                             double Iout_spice,
                                             double Lm) {
    // D is fixed by the analytical-derived PWM and unchanged by Vout
    // settling because the SPICE netlist hard-codes tOn = D_nominal·T.
    double Vor = s.N * s.Vout;
    double D = Vor / (s.Vin + Vor);
    double Iavg = Iout_spice / s.N / (1.0 - D);
    double dIpp = s.Vin * D / (Lm * s.Fs);
    double ripple = dIpp / (2.0 * Iavg);   // (peak-trough)/(2·Iavg) — matches
                                           // how Flyback.cpp consumes the
                                           // currentRippleRatio field:
                                           //   pp = center · ripple · 2

    OpenMagnetics::Flyback f;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    f.set_input_voltage(iv);
    f.set_diode_voltage_drop(0.0);
    f.set_efficiency(1.0);
    f.set_current_ripple_ratio(ripple);
    f.set_maximum_duty_cycle(0.5);
    f.set_maximum_drain_source_voltage(s.Vin + s.N * s.Vout + 50.0);
    OpenMagnetics::FlybackOperatingPoint op;
    op.set_output_voltages({s.Vout});         // nominal Vout keeps D invariant
    op.set_output_currents({Iout_spice});     // settled Iout from SPICE
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    f.set_operating_points(std::vector<OpenMagnetics::FlybackOperatingPoint>{op});
    return f;
}

// Lm chosen so the analytical primary-current ripple matches what the
// SPICE primary current produces — keeps PtP shape comparison fair.
//   Vor = N·Vout, D = Vor/(Vin+Vor), I_center = (Iout/N)/(1-D),
//   ΔIL_pp = 2·ripple·I_center,  Lm = Vin·D / (Fs·ΔIL_pp).
double consistent_lm(const RefDesignSpec& s) {
    double Vor    = s.N * s.Vout;
    double D      = Vor / (s.Vin + Vor);
    double center = (s.Iout / s.N) / (1.0 - D);
    double dIL    = 2.0 * s.ripple * center;
    return s.Vin * D / (s.Fs * dIL);
}

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== Flyback PtP — " << s.name << " ==========\n";

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    auto fly = build(s);
    const double Lm = consistent_lm(s);
    const std::vector<double> tr{s.N};

    // ── Analytical pass (nominal OP) — used for diagnostic only ────────
    // The NRMSE comparison further down uses a second analytical pass
    // anchored on the SPICE-settled operating point (see
    // `build_for_shape_match`) because open-loop SPICE settles Vout
    // below nominal and the analytical primary-current
    // offset/ripple-ratio shape is highly sensitive to that settling.
    auto analyticalOpsNominal = fly.process_operating_points(tr, Lm);
    REQUIRE(!analyticalOpsNominal.empty());

    // ── SPICE switching pass ──────────────────────────────────────────
    // Cout RC tail at high-Fs / low-Iout corners is long; mirror Buck
    // pattern at 400 settling periods.
    fly.set_num_steady_state_periods(400);
    fly.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = fly.simulate_and_extract_topology_waveforms(tr, Lm, 1);
    auto simOps = fly.simulate_and_extract_operating_points(tr, Lm);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency.  Open-loop SPICE settles Vout below
    // nominal proportional to loss, but the resistive load (wired for
    // Rload_nom = Vout_nom/Iout_nom) preserves the ratio.
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

    // Gate 3 — power balance (Pin must dominate Pout, lossy but causal).
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

    // Gate 4 — primary-current shape NRMSE (analytical vs SPICE).
    // Re-run analytical with SPICE-settled Iout and auto-physics ripple
    // ratio so both waveforms reference the same operating point.
    auto flyShape = build_for_shape_match(s, ioutMean, Lm);
    auto analyticalOps = flyShape.process_operating_points(tr, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty());
    REQUIRE(!sTime.empty());
    auto sResampled = ptp_interp(sTime, sData, 256);
    // 256-shift NRMSE — the sub-period commutation of the secondary
    // diode (with the SPICE Csnub charging through the primary-current
    // sense as the drain voltage swings) introduces a ~40-100 ns lag
    // between SPICE's primary-current step and the analytical ideal
    // step; the default 64-shift phase grid (period/64 ≈ 60-220 ns at
    // these Fs) is too coarse to align that lag and leaves a ~5 pp
    // spurious NRMSE residual. 256 shifts resolves it.
    const double nrmse = ptp_nrmse(aResampled, sResampled, 256);
    std::cout << "  iPri NRMSE (analytical vs SPICE) = "
              << 100.0 * nrmse << " %   (gate "
              << 100.0 * s.tol_nrmse << " %)\n";
    CHECK(nrmse < s.tol_nrmse);

    OpenMagnetics::Testing::dump_waveforms_csv(
        std::string("flyback_") + s.name, analyticalOps[0], simOps[0]);
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────
// Per-design TEST_CASEs.  Specs match kFlybackRefDesign<1,2,3> in
// TestTopologyFlyback.cpp (canonical source).
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("Flyback reference design PtP — TI PMP30817 (1.2 W primary-side-regulated)",
          "[converter-model][flyback-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "PMP30817",
        /*Vin*/          24.0,
        /*Vout*/         6.0,
        /*Iout*/         0.2,
        /*Fs*/           250e3,
        /*N*/            3.0,
        /*ripple*/       0.4,
        /*tol_walltime*/ 5.0,
        /*tol_rload_pct*/ 2.0,
        /*tol_loss_max*/ 0.60,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("Flyback reference design PtP — TI LM5180EVM-DUAL (3 W dual-output)",
          "[converter-model][flyback-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "LM5180EVM-DUAL",
        /*Vin*/          24.0,
        /*Vout*/         15.0,
        /*Iout*/         0.2,
        /*Fs*/           200e3,
        /*N*/            1.5,
        /*ripple*/       0.4,
        /*tol_walltime*/ 5.0,
        /*tol_rload_pct*/ 2.0,
        /*tol_loss_max*/ 0.60,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("Flyback reference design PtP — TI TIDA-00709 (33 W QR-flyback)",
          "[converter-model][flyback-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "TIDA-00709",
        /*Vin*/          120.0,
        /*Vout*/         12.0,
        /*Iout*/         2.75,
        /*Fs*/           70e3,
        /*N*/            8.0,
        /*ripple*/       0.4,
        /*tol_walltime*/ 5.0,
        /*tol_rload_pct*/ 2.0,
        /*tol_loss_max*/ 0.60,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}
