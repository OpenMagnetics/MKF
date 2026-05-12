// Phase 8f — Golden PtP regression tests for three Llc reference
// operating points spanning common LLC application classes:
//
//   1. Telecom 120 W server brick     — 400 V → 12 V @ 10 A, 100 kHz
//        Half-bridge LLC at fr.  Q=0.4, k=Lm/Ls=5.0.  STMicro AN2644.
//   2. ATX 240 W AC-DC rail           — 400 V → 24 V @ 10 A, 100 kHz
//        Half-bridge LLC at fr.  Q=0.4, k=5.0.  Infineon AN2008-03.
//   3. EV onboard 1 kW charger        — 400 V → 48 V @ 20 A, 100 kHz
//        Half-bridge LLC at fr.  Q=0.4, k=5.0.  ON Semi NCP1399 class.
//
// LLC primary (resonant tank) current is approximately sinusoidal at fr
// with a triangular magnetizing component.
//
// NRMSE gate: < 35 %.  The existing Test_Llc_PtP_AnalyticalVsNgspice
// (TestTopologyLlc.cpp:1326) gates at < 30 % at a single 12 V/10 A
// operating point; we allow 35 % to absorb the wider Vout range.
//
// Tags: [converter-model][llc-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include "converter_models/Llc.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "Definitions.h"
#include "PtpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin, Vin_min, Vin_max;
    double Vout, Iout, Fs, Fmin, Fmax;
    double tol_walltime, tol_vin_pct;
    double tol_vout_max;
    double tol_nrmse;
};

nlohmann::json build_fixture(const RefDesignSpec& s) {
    return nlohmann::json{
        {"inputVoltage", {{"nominal", s.Vin}, {"minimum", s.Vin_min}, {"maximum", s.Vin_max}}},
        {"bridgeType", "halfBridge"},
        {"minSwitchingFrequency", s.Fmin},
        {"maxSwitchingFrequency", s.Fmax},
        {"qualityFactor", 0.4},
        {"inductanceRatio", 5.0},
        {"operatingPoints", nlohmann::json::array({
            nlohmann::json{
                {"ambientTemperature", 25.0},
                {"outputVoltages", {s.Vout}},
                {"outputCurrents", {s.Iout}},
                {"switchingFrequency", s.Fs}
            }
        })}
    };
}

using namespace OpenMagnetics::Testing;

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== LLC PtP — " << s.name << " ==========\n";
    NgspiceRunner runner;
    if (!runner.is_available()) { WARN("ngspice not available"); return; }

    OpenMagnetics::Llc llc(build_fixture(s));
    auto req = llc.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios())
        turnsRatios.push_back(resolve_dimensional_values(tr));
    REQUIRE(!turnsRatios.empty());
    const double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    REQUIRE(Lm > 0.0);
    std::cout << "  n=" << turnsRatios[0] << "  Lm=" << 1e6*Lm << " uH\n";

    auto analyticalOps = llc.process_operating_points(turnsRatios, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty()); REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    llc.set_num_steady_state_periods(50);
    llc.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = llc.simulate_and_extract_topology_waveforms(turnsRatios, Lm, 1);
    auto simOps = llc.simulate_and_extract_operating_points(turnsRatios, Lm, 1);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    REQUIRE(!wfs.empty()); REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — Vin sanity
    const auto& vinW = wfs[0].get_input_voltage();
    const auto& vinD = vinW.get_data();
    const auto vinT_opt = vinW.get_time();
    REQUIRE(!vinD.empty()); REQUIRE(vinT_opt.has_value());
    const auto& vinT = vinT_opt.value();
    const double vinMean = ConverterPortChecks::time_weighted_mean(vinT, vinD);
    const double vinErr = (vinMean - s.Vin) / s.Vin;
    std::cout << "  Vin_spice = " << vinMean << " V (err "
              << 100.0*vinErr << " %, gate " << s.tol_vin_pct << " %)\n";
    CHECK(std::fabs(vinErr) < s.tol_vin_pct/100.0);

    // Gate 3 — Vout regulation (LLC at fr is open-loop unity gain;
    // input_current is bursty bridge-leg current per ConverterPortChecks
    // §5.1 guidance which omits an Iin bound).
    const auto& voutW = wfs[0].get_output_voltages()[0];
    const auto& voutD = voutW.get_data();
    const auto voutT_opt = voutW.get_time();
    REQUIRE(!voutD.empty()); REQUIRE(voutT_opt.has_value());
    const auto& voutT = voutT_opt.value();
    const double voutMean = ConverterPortChecks::time_weighted_mean(voutT, voutD);
    const double voutErr = (voutMean - s.Vout) / s.Vout;
    std::cout << "  Vout_spice = " << voutMean << " V (err "
              << 100.0*voutErr << " %, gate ±" << 100.0*s.tol_vout_max << " %)\n";
    CHECK(std::fabs(voutErr) < s.tol_vout_max);

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

TEST_CASE("LLC reference design PtP — Telecom 120 W server brick "
          "(400V→12V/10A 100 kHz half-bridge)",
          "[converter-model][llc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"Telecom-120W", 400.0, 370.0, 410.0, 12.0, 10.0,
                    100e3, 80e3, 200e3,
                    30.0, 2.0, 0.30, 0.30};
    run_ptp_gates(s);
}

TEST_CASE("LLC reference design PtP — ATX 240 W AC-DC rail "
          "(400V→24V/10A 100 kHz half-bridge)",
          "[converter-model][llc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"ATX-240W", 400.0, 370.0, 410.0, 24.0, 10.0,
                    100e3, 80e3, 200e3,
                    30.0, 2.0, 0.30, 0.30};
    run_ptp_gates(s);
}

TEST_CASE("LLC reference design PtP — EV onboard 1 kW charger "
          "(400V→48V/20A 100 kHz half-bridge)",
          "[converter-model][llc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"EV-1kW", 400.0, 370.0, 410.0, 48.0, 20.0,
                    100e3, 80e3, 200e3,
                    30.0, 2.0, 0.30, 0.30};
    run_ptp_gates(s);
}
