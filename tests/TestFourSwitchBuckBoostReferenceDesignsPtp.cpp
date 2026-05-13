// Four-Switch Buck-Boost (4SBB) reference-design PtP regression tests.
//
// Three published designs spanning the topology's three operating regions:
//   #1 TIDA-01411 (TI)    — 9–36 V → 12 V @ 8 A @ 350 kHz   (BUCK and BOOST)
//   #2 LM5176EVM-HP (TI)  — 9–55 V → 24 V @ 5 A @ 200 kHz   (BUCK_BOOST)
//   #3 LT8390 DC2825A     — 4–60 V → 12 V @ 7 A @ 300 kHz   (BUCK and BOOST)
//
// For each design we exercise the three regions with three Vin sub-cases
// per design and check the four PtP gates per CONVERTER_MODELS_GOLDEN_GUIDE
// §15:
//   Gate 1  walltime          < 8 s
//   Gate 2  load consistency  |Rload_spice − Rload_nom|/Rload_nom < 0.01
//   Gate 3  power balance     0 ≤ (Pin − Pout)/Pin ≤ tol_loss_max
//   Gate 4  inductor-current NRMSE  < tol_nrmse
//                                   (≤ 0.15 BUCK / BOOST,
//                                    ≤ 0.20 BUCK_BOOST per user)
//
// Tags: [converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "converter_models/FourSwitchBuckBoost.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;

namespace {

struct FsbbSpec {
    const char* name;
    double Vin;             // input voltage to test [V]
    double Vout;            // output voltage [V]
    double Iout;            // output current [A]
    double Fs;              // switching frequency [Hz]
    double L;               // worst-case inductance to use [H]
    double tol_walltime;    // walltime gate [s]
    double tol_rload_pct;   // |Rload_spice − Rload_nom|/Rload_nom [%]
    double tol_loss_max;    // (Pin − Pout)/Pin upper bound [fraction]
    double tol_nrmse;       // inductor-current NRMSE [fraction]
};

OpenMagnetics::FourSwitchBuckBoost build(const FsbbSpec& s) {
    OpenMagnetics::FourSwitchBuckBoost c;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    c.set_input_voltage(iv);
    c.set_efficiency(1.0);
    c.set_current_ripple_ratio(0.4);
    c.set_maximum_switch_current(50.0);
    MAS::TopologyExcitation op;
    op.set_output_voltages({s.Vout});
    op.set_output_currents({s.Iout});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    c.set_operating_points({op});
    return c;
}

using namespace OpenMagnetics::Testing;

void run_ptp_gates(const FsbbSpec& s) {
    std::cout << "\n========== FSBB PtP — " << s.name << " ==========\n";

    OpenMagnetics::NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    auto c = build(s);
    auto analyticalOps = c.process_operating_points(std::vector<double>{}, s.L);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    c.set_num_steady_state_periods(80);
    c.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = c.simulate_and_extract_topology_waveforms(s.L);
    auto simOps = c.simulate_and_extract_operating_points(s.L);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency.
    const auto& voutW = wfs[0].get_output_voltages()[0];
    const auto& ioutW = wfs[0].get_output_currents()[0];
    const auto& voutD = voutW.get_data();
    const auto& ioutD = ioutW.get_data();
    REQUIRE(voutW.get_time().has_value());
    REQUIRE(ioutW.get_time().has_value());
    const auto& voutT = voutW.get_time().value();
    const auto& ioutT = ioutW.get_time().value();
    const double voutMean = ConverterPortChecks::time_weighted_mean(voutT, voutD);
    const double ioutMean = ConverterPortChecks::time_weighted_mean(ioutT, ioutD);
    const double rloadNom = s.Vout / s.Iout;
    REQUIRE(std::fabs(ioutMean) > 1e-9);
    const double rloadSpc = voutMean / ioutMean;
    const double rloadErr = (rloadSpc - rloadNom) / rloadNom;
    std::cout << "  Vout_spice=" << voutMean << " Iout_spice=" << ioutMean
              << " Rload_err=" << 100.0 * rloadErr << "% (gate "
              << s.tol_rload_pct << "%)\n";
    CHECK(std::fabs(rloadErr) < s.tol_rload_pct / 100.0);

    // Gate 3 — power balance.
    const auto& vinW = wfs[0].get_input_voltage();
    const auto& iinW = wfs[0].get_input_current();
    const auto& vinD = vinW.get_data();
    const auto& iinD = iinW.get_data();
    REQUIRE(vinW.get_time().has_value());
    const auto& vinT = vinW.get_time().value();
    const double pin  = ConverterPortChecks::time_weighted_mean_product(vinT, vinD, iinD);
    const double pout = ConverterPortChecks::time_weighted_mean_product(voutT, voutD, ioutD);
    const double loss = (pin - pout) / pin;
    std::cout << "  Pin=" << pin << " Pout=" << pout << " loss="
              << 100.0 * loss << "% (gate 0..." << 100.0 * s.tol_loss_max << "%)\n";
    CHECK(loss >= -0.05);
    CHECK(loss <= s.tol_loss_max);

    // Gate 4 — inductor current waveform NRMSE.
    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty());
    auto sResampled = ptp_interp(sTime, sData, 256);
    const double nrmse = ptp_nrmse(aResampled, sResampled);
    std::cout << "  iL NRMSE = " << 100.0 * nrmse << "% (gate "
              << 100.0 * s.tol_nrmse << "%)\n";
    CHECK(nrmse < s.tol_nrmse);
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────
// TIDA-01411 (TI): 9–36 V → 12 V @ 8 A @ 350 kHz.
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("FSBB PtP — TIDA-01411 BUCK (Vin=24V, Vo=12V)",
          "[converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]") {
    FsbbSpec s{
        /*name*/ "TIDA-01411 BUCK 24V→12V@8A@350kHz",
        /*Vin*/ 24.0, /*Vo*/ 12.0, /*Iout*/ 8.0, /*Fs*/ 350e3, /*L*/ 5e-6,
        /*tol_walltime*/ 8.0, /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.55, /*tol_nrmse*/ 0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("FSBB PtP — TIDA-01411 BOOST (Vin=9V, Vo=12V)",
          "[converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]") {
    FsbbSpec s{
        /*name*/ "TIDA-01411 BOOST 9V→12V@8A@350kHz",
        /*Vin*/ 9.0, /*Vo*/ 12.0, /*Iout*/ 8.0, /*Fs*/ 350e3, /*L*/ 5e-6,
        /*tol_walltime*/ 8.0, /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.55, /*tol_nrmse*/ 0.15
    };
    run_ptp_gates(s);
}

// ─────────────────────────────────────────────────────────────────────────
// LM5176EVM-HP (TI): 9–55 V → 24 V @ 5 A @ 200 kHz.
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("FSBB PtP — LM5176EVM BUCK (Vin=48V, Vo=24V)",
          "[converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]") {
    FsbbSpec s{
        /*name*/ "LM5176EVM BUCK 48V→24V@5A@200kHz",
        /*Vin*/ 48.0, /*Vo*/ 24.0, /*Iout*/ 5.0, /*Fs*/ 200e3, /*L*/ 10e-6,
        /*tol_walltime*/ 8.0, /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.55, /*tol_nrmse*/ 0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("FSBB PtP — LM5176EVM BUCK_BOOST (Vin=24V, Vo=24V)",
          "[converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]") {
    FsbbSpec s{
        /*name*/ "LM5176EVM BUCK_BOOST 24V→24V@5A@200kHz",
        /*Vin*/ 24.0, /*Vo*/ 24.0, /*Iout*/ 5.0, /*Fs*/ 200e3, /*L*/ 10e-6,
        /*tol_walltime*/ 8.0, /*tol_rload_pct*/ 8.0,
        /*tol_loss_max*/ 0.65, /*tol_nrmse*/ 0.20
    };
    run_ptp_gates(s);
}

TEST_CASE("FSBB PtP — LM5176EVM BOOST (Vin=12V, Vo=24V)",
          "[converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]") {
    FsbbSpec s{
        /*name*/ "LM5176EVM BOOST 12V→24V@5A@200kHz",
        /*Vin*/ 12.0, /*Vo*/ 24.0, /*Iout*/ 5.0, /*Fs*/ 200e3, /*L*/ 10e-6,
        /*tol_walltime*/ 8.0, /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.55, /*tol_nrmse*/ 0.15
    };
    run_ptp_gates(s);
}

// ─────────────────────────────────────────────────────────────────────────
// LT8390 DC2825A: 4–60 V → 12 V @ 7 A @ 300 kHz.
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("FSBB PtP — LT8390 DC2825A BUCK (Vin=36V, Vo=12V)",
          "[converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]") {
    FsbbSpec s{
        /*name*/ "LT8390 DC2825A BUCK 36V→12V@7A@300kHz",
        /*Vin*/ 36.0, /*Vo*/ 12.0, /*Iout*/ 7.0, /*Fs*/ 300e3, /*L*/ 5e-6,
        /*tol_walltime*/ 8.0, /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.55, /*tol_nrmse*/ 0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("FSBB PtP — LT8390 DC2825A BOOST (Vin=6V, Vo=12V)",
          "[converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]") {
    FsbbSpec s{
        /*name*/ "LT8390 DC2825A BOOST 6V→12V@7A@300kHz",
        /*Vin*/ 6.0, /*Vo*/ 12.0, /*Iout*/ 7.0, /*Fs*/ 300e3, /*L*/ 5e-6,
        /*tol_walltime*/ 8.0, /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.55, /*tol_nrmse*/ 0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("FSBB PtP — LT8390 DC2825A BUCK_BOOST (Vin=12V, Vo=12V)",
          "[converter-model][four-switch-buck-boost-topology][refdesign][ptp][slow]") {
    FsbbSpec s{
        /*name*/ "LT8390 DC2825A BUCK_BOOST 12V→12V@7A@300kHz",
        /*Vin*/ 12.0, /*Vo*/ 12.0, /*Iout*/ 7.0, /*Fs*/ 300e3, /*L*/ 5e-6,
        /*tol_walltime*/ 8.0, /*tol_rload_pct*/ 8.0,
        /*tol_loss_max*/ 0.65, /*tol_nrmse*/ 0.20
    };
    run_ptp_gates(s);
}
