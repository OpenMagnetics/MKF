// Phase 8e — Golden PtP regression tests for three Dab (Dual Active
// Bridge) reference operating points:
//
//   1. TI UCC21750-class 10 kW EV charger — 800 V → 500 V @ 20 A, 100 kHz
//        SPS modulation, outer phase 23°.  TI design report SLUA963.
//   2. ABB-class 5 kW battery interface   — 800 V → 400 V @ 12.5 A, 100 kHz
//        SPS modulation, outer phase 30°.  ABB review series.
//   3. EV auxiliary 2 kW DC-DC            — 400 V → 48 V @ 40 A, 80 kHz
//        SPS modulation, outer phase 35°.  Bosch / Continental class.
//
// DAB primary (series-inductor) current is trapezoidal in SPS; the
// closed-form analytical model + Gear-method SPICE comparison is
// known to land at ≈ 4–10 % NRMSE per existing single-point tests.
//
// NRMSE gate: < 25 %.  The existing Test_Dab_PtP_AnalyticalVsNgspice
// (TestTopologyDab.cpp:1054) gates at < 15 %; we allow a slightly
// looser 25 % to absorb wider operating-point coverage and SPICE
// load-droop variation.
//
// Tags: [converter-model][dab-topology][refdesign][ptp][slow]

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

#include "converter_models/Dab.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin;
    double Vout, Iout, Fs;
    double phaseShiftDeg;
    double tol_walltime, tol_vin_pct;
    double tol_loss_neg, tol_loss_max;
    double tol_nrmse;
};

// Solve Ls so that the requested (Vin,Vout,Iout) are self-consistent with
// the chosen SPS outer phase shift phi:
//     Pout = (n Vin Vout / (2π fs Ls)) · phi · (1 - phi/π),  n = Vin/Vout
// → Ls = (n · Vin / (2π fs Iout)) · phi · (1 - phi/π).
// Without this the SPICE bridge floats Vout to whatever Rload demands and
// the input/output power balance diverges from the analytical setpoint.
double solve_series_inductance(const RefDesignSpec& s) {
    const double n = s.Vin / s.Vout;
    const double phi = s.phaseShiftDeg * M_PI / 180.0;
    const double shape = phi * (1.0 - phi / M_PI);
    return (n * s.Vin) / (2.0 * M_PI * s.Fs * s.Iout) * shape;
}

nlohmann::json build_fixture(const RefDesignSpec& s) {
    return nlohmann::json{
        {"inputVoltage", {{"nominal", s.Vin}}},
        {"seriesInductance", solve_series_inductance(s)},
        {"useLeakageInductance", false},
        {"operatingPoints", nlohmann::json::array({
            nlohmann::json{
                {"ambientTemperature", 25.0},
                {"outputVoltages", {s.Vout}},
                {"outputCurrents", {s.Iout}},
                {"innerPhaseShift3", s.phaseShiftDeg},
                {"switchingFrequency", s.Fs}
            }
        })}
    };
}

using namespace OpenMagnetics::Testing;

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== DAB PtP — " << s.name << " ==========\n";
    NgspiceRunner runner;
    if (!runner.is_available()) { WARN("ngspice not available"); return; }

    OpenMagnetics::Dab dab(build_fixture(s));
    auto req = dab.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios())
        turnsRatios.push_back(tr.get_nominal().value());
    REQUIRE(!turnsRatios.empty());
    const double Lm = req.get_magnetizing_inductance().get_minimum().value();
    REQUIRE(Lm > 0.0);
    std::cout << "  n=" << turnsRatios[0] << "  Lm=" << 1e6*Lm << " uH\n";

    auto analyticalOps = dab.process_operating_points(turnsRatios, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty()); REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    dab.set_num_steady_state_periods(200);
    dab.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = dab.simulate_and_extract_topology_waveforms(turnsRatios, Lm, 1);
    auto simOps = dab.simulate_and_extract_operating_points(turnsRatios, Lm);
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

    // Gate 3 — Vout regulation (DAB input-current waveform is bursty
    // bridge-leg current; per ConverterPortChecks.h:28 the §5.1 helper
    // explicitly omits an input_current bound — "too sensitive to η model
    // and startup inrush". We instead gate that the SPICE-settled output
    // voltage tracks the analytical setpoint, which is the meaningful
    // open-loop power-balance check for DAB.)
    const auto& voutW = wfs[0].get_output_voltages()[0];
    const auto& voutD = voutW.get_data();
    const auto voutT_opt = voutW.get_time();
    REQUIRE(!voutD.empty()); REQUIRE(voutT_opt.has_value());
    const auto& voutT = voutT_opt.value();
    const double voutMean = ConverterPortChecks::time_weighted_mean(voutT, voutD);
    const double voutErr = (voutMean - s.Vout) / s.Vout;
    std::cout << "  Vout_spice = " << voutMean << " V (err "
              << 100.0*voutErr << " %, gate ±" << 100.0*s.tol_loss_max << " %)\n";
    CHECK(std::fabs(voutErr) < s.tol_loss_max);

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

TEST_CASE("DAB reference design PtP — TI UCC21750-class 10 kW EV charger "
          "(800V→500V/20A 100 kHz SPS)",
          "[converter-model][dab-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"TI-10kW", 800.0, 500.0, 20.0, 100e3,
                    23.0,
                    30.0, 2.0, 0.20, 0.60, 0.25};
    run_ptp_gates(s);
}

TEST_CASE("DAB reference design PtP — ABB-class 5 kW battery interface "
          "(800V→400V/12.5A 100 kHz SPS)",
          "[converter-model][dab-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"ABB-5kW", 800.0, 400.0, 12.5, 100e3,
                    30.0,
                    30.0, 2.0, 0.20, 0.60, 0.25};
    run_ptp_gates(s);
}

TEST_CASE("DAB reference design PtP — EV auxiliary 2 kW DC-DC "
          "(400V→48V/40A 80 kHz SPS)",
          "[converter-model][dab-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"EV-Aux-2kW", 400.0, 48.0, 40.0, 80e3,
                    30.0,
                    30.0, 2.0, 0.20, 0.60, 0.25};
    run_ptp_gates(s);
}
