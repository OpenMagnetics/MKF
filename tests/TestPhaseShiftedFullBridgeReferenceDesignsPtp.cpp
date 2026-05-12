// Phase 8g — Golden PtP regression tests for three Psfb (Phase-Shifted
// Full Bridge) reference operating points spanning common ZVS-PSFB
// applications:
//
//   1. Telecom 600 W 12 V brick      — 400 V → 12 V @ 50 A, 100 kHz
//        Full-bridge rectifier.  TI UCC28950 reference, Pressman ch.8.
//   2. Server 1.2 kW 24 V rail       — 400 V → 24 V @ 50 A, 100 kHz
//        Full-bridge rectifier.  TI SLUA560 worked example.
//   3. EV onboard 1 kW 48 V auxiliary — 400 V → 48 V @ 21 A, 100 kHz
//        Full-bridge rectifier.  ON Semi NCP1395 class.
//
// PSFB primary current is a trapezoidal lagging-leg waveform with a
// freewheeling interval where the primary recirculates through the
// rectifier diodes.
//
// NRMSE gate: < 50 %.  The existing Test_Psfb_PtP_AnalyticalVsNgspice
// (TestTopologyPhaseShiftedFullBridge.cpp:615) gates at < 18 % at a
// single 12 V/50 A operating point; we allow 50 % to absorb wider Vout
// coverage and droop-induced phase-mismatch.
//
// Tags: [converter-model][psfb-topology][refdesign][ptp][slow]

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

#include "converter_models/PhaseShiftedFullBridge.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "Definitions.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin;
    double Vout, Iout, Fs;
    double seriesInductance;
    double phaseShiftDeg;
    double tol_walltime, tol_vin_pct;
    double tol_vout_max;
    double tol_nrmse;
};

nlohmann::json build_fixture(const RefDesignSpec& s) {
    return nlohmann::json{
        {"inputVoltage", {{"nominal", s.Vin}}},
        {"rectifierType", "fullBridge"},
        {"seriesInductance", s.seriesInductance},
        {"operatingPoints", nlohmann::json::array({
            nlohmann::json{
                {"ambientTemperature", 25.0},
                {"outputVoltages", {s.Vout}},
                {"outputCurrents", {s.Iout}},
                {"switchingFrequency", s.Fs},
                {"phaseShift", s.phaseShiftDeg}
            }
        })}
    };
}

// Local PtP helper functions (self-contained file, no header)
std::vector<double> ptp_interp(const std::vector<double>& t,
                                const std::vector<double>& d, int N) {
    std::vector<double> out(N);
    double T = t.back();
    for (int i = 0; i < N; ++i) {
        double ti = T * i / (N - 1);
        int lo = 0;
        for (int k = 0; k + 1 < (int)t.size(); ++k) if (t[k] <= ti) lo = k;
        int hi = std::min(lo + 1, (int)t.size() - 1);
        double dt = t[hi] - t[lo];
        out[i] = (dt < 1e-20) ? d[hi]
                              : d[lo] + (ti - t[lo]) / dt * (d[hi] - d[lo]);
    }
    return out;
}
double ptp_nrmse(const std::vector<double>& ref, const std::vector<double>& sim) {
    int N = (int)ref.size();
    double rm = 0, sm = 0;
    for (int i = 0; i < N; ++i) { rm += ref[i]; sm += sim[i]; }
    rm /= N; sm /= N;
    std::vector<double> r(N), s(N);
    double rA = 0, sA = 0;
    for (int i = 0; i < N; ++i) {
        r[i] = ref[i] - rm; s[i] = sim[i] - sm;
        rA += r[i]*r[i]; sA += s[i]*s[i];
    }
    rA = std::sqrt(rA/N); sA = std::sqrt(sA/N);
    if (rA < 1e-10 || sA < 1e-10) return 1.0;
    for (int i = 0; i < N; ++i) { r[i] /= rA; s[i] /= sA; }
    int ns = std::min(N, 64);
    double best = std::numeric_limits<double>::max();
    for (int ss = 0; ss < ns; ++ss) {
        int sh = ss * N / ns;
        double ssd = 0;
        for (int k = 0; k < N; ++k) {
            double e = r[k] - s[(k + sh) % N];
            ssd += e * e;
        }
        if (ssd < best) best = ssd;
    }
    return std::sqrt(best/N);
}
std::vector<double> ptp_current(const OperatingPoint& op, size_t wi = 0) {
    if (wi >= op.get_excitations_per_winding().size()) return {};
    auto& e = op.get_excitations_per_winding()[wi];
    if (!e.get_current() || !e.get_current()->get_waveform()) return {};
    return e.get_current()->get_waveform()->get_data();
}
std::vector<double> ptp_current_time(const OperatingPoint& op, size_t wi = 0) {
    if (wi >= op.get_excitations_per_winding().size()) return {};
    auto& e = op.get_excitations_per_winding()[wi];
    if (!e.get_current() || !e.get_current()->get_waveform()) return {};
    auto tv = e.get_current()->get_waveform()->get_time();
    return tv ? tv.value() : std::vector<double>{};
}

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== PSFB PtP — " << s.name << " ==========\n";
    NgspiceRunner runner;
    if (!runner.is_available()) { WARN("ngspice not available"); return; }

    OpenMagnetics::Psfb psfb(build_fixture(s));
    auto req = psfb.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios())
        turnsRatios.push_back(resolve_dimensional_values(tr));
    REQUIRE(!turnsRatios.empty());
    const double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    REQUIRE(Lm > 0.0);
    std::cout << "  n=" << turnsRatios[0] << "  Lm=" << 1e6*Lm << " uH\n";

    auto analyticalOps = psfb.process_operating_points(turnsRatios, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty()); REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    psfb.set_num_steady_state_periods(200);
    psfb.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = psfb.simulate_and_extract_topology_waveforms(turnsRatios, Lm);
    auto simOps = psfb.simulate_and_extract_operating_points(turnsRatios, Lm);
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

    // Gate 3 — Vout regulation (PSFB input current is bursty bridge-leg
    // current per ConverterPortChecks §5.1 guidance which omits an Iin
    // bound; we instead gate that SPICE-settled Vout tracks setpoint).
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

TEST_CASE("PSFB reference design PtP — Telecom 600 W 12 V brick "
          "(400V→12V/50A 100 kHz)",
          "[converter-model][psfb-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"Telecom-600W", 400.0, 12.0, 50.0, 100e3,
                    5e-6, 126.0,
                    30.0, 2.0, 0.30, 0.50};
    run_ptp_gates(s);
}

TEST_CASE("PSFB reference design PtP — Server 1.2 kW 24 V rail "
          "(400V→24V/50A 100 kHz)",
          "[converter-model][psfb-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"Server-1.2kW", 400.0, 24.0, 50.0, 100e3,
                    5e-6, 126.0,
                    30.0, 2.0, 0.30, 0.50};
    run_ptp_gates(s);
}

TEST_CASE("PSFB reference design PtP — EV onboard 1 kW 48 V auxiliary "
          "(400V→48V/21A 100 kHz)",
          "[converter-model][psfb-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"EV-Aux-1kW", 400.0, 48.0, 21.0, 100e3,
                    5e-6, 126.0,
                    30.0, 2.0, 0.30, 0.50};
    run_ptp_gates(s);
}
