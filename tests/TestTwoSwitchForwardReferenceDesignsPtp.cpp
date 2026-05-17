// Phase 8 — Golden PtP regression tests for three TwoSwitchForward
// (TSF) reference operating points.  TSF is the classic two-switch
// forward (TI SLUP089, Pressman ch.6 §6.3): high-side + low-side
// switches in series with the primary, body diodes provide automatic
// flux reset (no tertiary winding, no RCD clamp).  Switch stress is
// clamped to Vin (vs single-switch forward where the switch sees
// 2·Vin during reset).  D ≤ 0.5.
//
//   1. SLUP089-class 100 W telecom    — 48 V → 12 V @ 8 A,  250 kHz
//   2. UCC2897-class 50 W (2SF mode)  — 48 V → 5 V  @ 10 A, 200 kHz
//   3. Erickson §6.3-class 25 W       — 24 V → 5 V  @ 5 A,  150 kHz
//        Erickson-Maksimović "Fundamentals of Power Electronics" 3rd ed.
//
// NRMSE gate: < 70 %.  The existing
// Test_TwoSwitchForward_PtP_AnalyticalVsNgspice gates at < 35 % at
// a single lightly-loaded operating point; these three ref designs
// run at higher load and the open-loop SPICE Vout droops 25–40 %
// below nominal, shifting the actual duty cycle off the analytical
// model.  70 % accommodates this droop regime.
//
// Tags: [converter-model][two-switch-forward-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "converter_models/TwoSwitchForward.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin, Vin_min, Vin_max;
    double Vout, Iout, Fs, D;
    double tol_walltime, tol_vin_pct;
    double tol_loss_neg, tol_loss_max;
    double tol_nrmse;
};

OpenMagnetics::TwoSwitchForward build(const RefDesignSpec& s) {
    OpenMagnetics::TwoSwitchForward fwd;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin_min);
    iv.set_maximum(s.Vin_max);
    fwd.set_input_voltage(iv);
    fwd.set_diode_voltage_drop(0.5);
    fwd.set_efficiency(0.9);
    fwd.set_current_ripple_ratio(0.3);
    fwd.set_duty_cycle(s.D);
    ForwardOperatingPoint op;
    op.set_output_voltages({s.Vout});
    op.set_output_currents({s.Iout});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    fwd.set_operating_points({op});
    return fwd;
}

using namespace OpenMagnetics::Testing;

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== TSF PtP — " << s.name << " ==========\n";
    NgspiceRunner runner;
    if (!runner.is_available()) { WARN("ngspice not available"); return; }

    auto fwd = build(s);
    auto designReqs = fwd.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : designReqs.get_turns_ratios())
        turnsRatios.push_back(tr.get_nominal().value());
    REQUIRE(!turnsRatios.empty());
    const double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();
    REQUIRE(Lm > 0.0);
    std::cout << "  N=" << turnsRatios[0] << "  Lm=" << 1e6*Lm << " uH\n";

    auto analyticalOps = fwd.process_operating_points(turnsRatios, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty()); REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    fwd.set_num_steady_state_periods(400);
    fwd.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = fwd.simulate_and_extract_topology_waveforms(turnsRatios, Lm);
    auto simOps = fwd.simulate_and_extract_operating_points(turnsRatios, Lm);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    REQUIRE(!wfs.empty()); REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — Vin sanity
    const auto& vinW = wfs[0].get_input_voltage();
    const auto& iinW = wfs[0].get_input_current();
    const auto& vinD = vinW.get_data();
    const auto& iinD = iinW.get_data();
    const auto vinT_opt = vinW.get_time();
    REQUIRE(!vinD.empty()); REQUIRE(!iinD.empty()); REQUIRE(vinT_opt.has_value());
    const auto& vinT = vinT_opt.value();
    const double vinMean = ConverterPortChecks::time_weighted_mean(vinT, vinD);
    const double vinErr = (vinMean - s.Vin) / s.Vin;
    std::cout << "  Vin_spice = " << vinMean << " V (err "
              << 100.0*vinErr << " %, gate " << s.tol_vin_pct << " %)\n";
    CHECK(std::fabs(vinErr) < s.tol_vin_pct/100.0);

    // Gate 3 — Power balance, Pout reconstructed from Vout²/Rload_nom
    const auto& voutW = wfs[0].get_output_voltages()[0];
    const auto& voutD = voutW.get_data();
    const auto voutT_opt = voutW.get_time();
    REQUIRE(!voutD.empty()); REQUIRE(voutT_opt.has_value());
    const auto& voutT = voutT_opt.value();
    const double pin = ConverterPortChecks::time_weighted_mean_product(vinT, vinD, iinD);
    const double voutMs = ConverterPortChecks::time_weighted_mean_product(voutT, voutD, voutD);
    const double rloadNom = s.Vout / s.Iout;
    const double poutRec = voutMs / rloadNom;
    const double loss = (pin - poutRec) / pin;
    std::cout << "  Pin = " << pin << " W, Pout_recon = " << poutRec
              << " W, loss = " << 100.0*loss << " %   (gate "
              << -100.0*s.tol_loss_neg << "..." << 100.0*s.tol_loss_max << " %)\n";
    CHECK(pin > 0.0);
    CHECK(loss >= -s.tol_loss_neg);
    CHECK(loss <=  s.tol_loss_max);

    // Gate 4 — primary current NRMSE
    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty()); REQUIRE(!sTime.empty());
    auto sResampled = ptp_interp(sTime, sData, 256);
    const double nrmse = ptp_nrmse(aResampled, sResampled);
    std::cout << "  iPri NRMSE = " << 100.0*nrmse << " %   (gate "
              << 100.0*s.tol_nrmse << " %)\n";
    CHECK(nrmse < s.tol_nrmse);
}

}  // namespace

TEST_CASE("TSF reference design PtP — SLUP089-class 100 W telecom "
          "(48V→12V/8A 250 kHz)",
          "[converter-model][two-switch-forward-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"SLUP089-100W", 48.0, 36.0, 72.0, 12.0, 8.0, 250e3, 0.40,
                    30.0, 2.0, 0.05, 0.60, 0.15};
    run_ptp_gates(s);
}

TEST_CASE("TSF reference design PtP — UCC2897-class 50 W (2SF mode) "
          "(48V→5V/10A 200 kHz)",
          "[converter-model][two-switch-forward-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"UCC2897-50W-2SF", 48.0, 36.0, 72.0, 5.0, 10.0, 200e3, 0.40,
                    30.0, 2.0, 0.05, 0.60, 0.15};
    run_ptp_gates(s);
}

TEST_CASE("TSF reference design PtP — Erickson §6.3-class 25 W "
          "(24V→5V/5A 150 kHz)",
          "[converter-model][two-switch-forward-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"Erickson-25W", 24.0, 20.0, 30.0, 5.0, 5.0, 150e3, 0.40,
                    30.0, 2.0, 0.05, 0.60, 0.15};
    run_ptp_gates(s);
}
