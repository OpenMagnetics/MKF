// =====================================================================
// Golden PtP regression tests for the Vienna 3-phase 3-level PFC.
// Mirrors tests/TestSrcReferenceDesignsPtp.cpp; gates per
// CONVERTER_MODELS_GOLDEN_GUIDE.md §15 (NRMSE ≤ 0.15, peak ratio band,
// power balance, Vin sanity, wall time).
//
// Vienna current capability (Phase 1+2 + Phase-1 SPICE, May 2026):
//   - viennaI variant, T-type switch, diode rectifier only
//   - peakOfLineOnly sampling strategy
//   - Single-phase boost emulation in SPICE (one phase modelled as a
//     stand-alone DC-DC boost at its peak-of-line operating point; the
//     resulting waveform is replicated to all three windings by
//     analytical symmetry).
//
// FIXME-vienna-1 (single-phase emulation):
//   Phase-1 SPICE models ONE phase as a stand-alone boost converter at
//   frozen line peak (V_phase_peak DC source, switch-to-neutral, diode
//   to upper half-bus). The three Vienna phase inductors are identical
//   by analytical symmetry, so we duplicate the Phase-A waveform across
//   "Phase A/B/C" windings in the extractor. The full 3-phase netlist
//   with all three phases at frozen-DC line-angle voltages, T-type
//   switches, split-bus caps + neutral-point modelling (per
//   VIENNA_PLAN.md §A.6) is deferred to Phase 3+.
//
// FIXME-vienna-2 (Vin sanity gate uses V_phase_peak):
//   Because the single-phase emulation has no real 3-phase bus, the
//   Vin-sanity gate checks input_voltage against V_phase_peak (the
//   frozen-DC source rail), not against V_LL.
//
// FIXME-vienna-3 (PtP shape vs amplitude):
//   ptp_nrmse is mean-subtracted AND RMS-normalised → shape-only.
//   We pair it with a peak-amplitude gate (current band [0.5, 1.7])
//   matching the SRC / CLLC convention.
//
// Tags: [converter-model][vienna-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <nlohmann/json.hpp>

#include "converter_models/Vienna.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"
#include "WaveformDumpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;
using namespace OpenMagnetics::Testing;

namespace {

struct ViennaRefSpec {
    const char* name;
    double V_LL;                // line-to-line RMS
    double Vdc;                 // bus
    double Fsw;                 // switching frequency
    double P_total;             // total 3-phase output power
    double rippleRatio;         // current ripple as fraction of I_pk
    double tol_walltime;
    double tol_vin_pct;
    double tol_loss_neg;
    double tol_loss_max;
    double tol_nrmse;
    double tol_peak_ratio_min = 0.5;
    double tol_peak_ratio_max = 1.7;
};

nlohmann::json build_fixture(const ViennaRefSpec& s) {
    nlohmann::json j;
    j["lineToLineVoltage"]   = {{"nominal", s.V_LL}};
    j["lineFrequency"]       = 50.0;
    j["outputDcVoltage"]     = s.Vdc;
    j["switchingFrequency"]  = s.Fsw;
    j["currentRippleRatio"]  = s.rippleRatio;
    j["efficiency"]          = 0.97;
    j["powerFactor"]         = 0.99;

    nlohmann::json op;
    op["outputVoltages"]     = nlohmann::json::array({s.Vdc});
    op["outputCurrents"]     = nlohmann::json::array({s.P_total / s.Vdc});
    op["switchingFrequency"] = s.Fsw;
    op["ambientTemperature"] = 25.0;
    j["operatingPoints"]     = nlohmann::json::array({op});
    return j;
}

void run_ptp_gates(const ViennaRefSpec& s) {
    std::cout << "\n========== Vienna PtP — " << s.name << " ==========\n";
    NgspiceRunner runner;
    if (!runner.is_available()) { WARN("ngspice not available"); return; }

    Vienna v(build_fixture(s));
    auto dr = v.process_design_requirements();
    REQUIRE(!dr.get_turns_ratios().empty());
    const double n  = dr.get_turns_ratios()[0].get_nominal().value();
    const double L  = resolve_dimensional_values(dr.get_magnetizing_inductance());
    REQUIRE(n  > 0.0);
    REQUIRE(L  > 0.0);
    const double V_phase_peak = std::sqrt(2.0) * s.V_LL / std::sqrt(3.0);
    std::cout << "  M=" << v.get_computed_modulation_index()
              << "  I_pk=" << v.get_computed_line_peak_current() << " A"
              << "  L=" << 1e6*v.get_computed_boost_inductance() << " uH"
              << "  V_phase_peak=" << V_phase_peak << " V\n";

    auto analyticalOps = v.process_operating_points({n}, L);
    REQUIRE(!analyticalOps.empty());
    REQUIRE(analyticalOps[0].get_excitations_per_winding().size() == 3);
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty()); REQUIRE(!aTime.empty());

    v.set_num_steady_state_periods(50);
    v.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = v.simulate_and_extract_topology_waveforms({n}, L, 1);
    auto simOps = v.simulate_and_extract_operating_points({n}, L);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    REQUIRE(!wfs.empty()); REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";
    CHECK(wallTime < s.tol_walltime);

    // Gate — Vin sanity: input_voltage stream is v(vphase) = V_phase_peak
    // (single-phase emulation; see FIXME-vienna-2).
    const auto& vinW = wfs[0].get_input_voltage();
    const auto& iinW = wfs[0].get_input_current();
    const auto& vinD = vinW.get_data();
    const auto& iinD = iinW.get_data();
    const auto vinT_opt = vinW.get_time();
    REQUIRE(!vinD.empty()); REQUIRE(!iinD.empty()); REQUIRE(vinT_opt.has_value());
    const auto& vinT = vinT_opt.value();
    const double vinMean = ConverterPortChecks::time_weighted_mean(vinT, vinD);
    const double vinErr  = (vinMean - V_phase_peak) / V_phase_peak;
    std::cout << "  Vin_spice = " << vinMean << " V (err "
              << 100.0*vinErr << " %, gate " << s.tol_vin_pct << " %)\n";
    CHECK(std::fabs(vinErr) < s.tol_vin_pct/100.0);

    // Gate — per-phase power balance. iinD is the reconstructed per-phase
    // DC equivalent current at V_phase_peak; Pout_per_phase = P_total/3.
    const auto& voutW = wfs[0].get_output_voltages()[0];
    const auto& voutD = voutW.get_data();
    const auto voutT_opt = voutW.get_time();
    REQUIRE(!voutD.empty()); REQUIRE(voutT_opt.has_value());
    const auto& voutT = voutT_opt.value();
    const double pin_per_phase = ConverterPortChecks::time_weighted_mean_product(vinT, vinD, iinD);
    const double voutMs        = ConverterPortChecks::time_weighted_mean_product(voutT, voutD, voutD);
    const double rloadNom      = s.Vdc / (s.P_total / s.Vdc);
    const double pout_total    = voutMs / rloadNom;
    const double pout_per_phase = pout_total / 3.0;
    const double loss   = (pin_per_phase > 0.0) ? (pin_per_phase - pout_per_phase) / pin_per_phase : 0.0;
    std::cout << "  Pin/phase = " << pin_per_phase << " W, Pout/phase = " << pout_per_phase
              << " W, loss = " << 100.0*loss << " %   (gate "
              << -100.0*s.tol_loss_neg << "..." << 100.0*s.tol_loss_max << " %)\n";
    CHECK(pin_per_phase > 0.0);
    CHECK(loss >= -s.tol_loss_neg);
    CHECK(loss <=  s.tol_loss_max);

    // Gate — phase-A inductor current NRMSE.
    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty()); REQUIRE(!sTime.empty());
    const double T_common = std::min(aTime.back(), sTime.back());
    auto resample_to = [](const std::vector<double>& t,
                          const std::vector<double>& d,
                          double Tmax, int N) {
        std::vector<double> out(N);
        for (int i = 0; i < N; ++i) {
            double ti = Tmax * i / (N - 1);
            int lo = 0;
            for (int k = 0; k + 1 < (int)t.size(); ++k) if (t[k] <= ti) lo = k;
            int hi = std::min(lo + 1, (int)t.size() - 1);
            double dt = t[hi] - t[lo];
            out[i] = (dt < 1e-20) ? d[hi]
                                  : d[lo] + (ti - t[lo]) / dt * (d[hi] - d[lo]);
        }
        return out;
    };
    auto aCommon = resample_to(aTime, aData, T_common, 256);
    auto sCommon = resample_to(sTime, sData, T_common, 256);
    const double nrmse = ptp_nrmse(aCommon, sCommon);
    auto stats = [](const std::vector<double>& v) {
        double mn = v[0], mx = v[0], mean = 0.0;
        for (double x : v) { mn = std::min(mn, x); mx = std::max(mx, x); mean += x; }
        return std::tuple<double,double,double>(mn, mx, mean / v.size());
    };
    auto [a_min, a_max, a_mean] = stats(aData);
    auto [s_min, s_max, s_mean] = stats(sData);
    std::cout << "  analytical: t_span=" << (aTime.back()-aTime.front())*1e6
              << " us, N=" << aData.size()
              << ", i_range=[" << a_min << ", " << a_max << "], mean=" << a_mean << "\n";
    std::cout << "  spice:      t_span=" << (sTime.back()-sTime.front())*1e6
              << " us, N=" << sData.size()
              << ", i_range=[" << s_min << ", " << s_max << "], mean=" << s_mean << "\n";
    std::cout << "  T_common = " << T_common*1e6 << " us\n";
    std::cout << "  iL NRMSE = " << 100.0*nrmse << " %   (gate "
              << 100.0*s.tol_nrmse << " %)\n";
    CHECK(nrmse < s.tol_nrmse);

    // Gate — peak-amplitude ratio (regression detector, see FIXME-vienna-3).
    const double a_pkpk = a_max - a_min;
    const double s_pkpk = s_max - s_min;
    const double peak_ratio = (a_pkpk > 1e-12) ? s_pkpk / a_pkpk : 0.0;
    std::cout << "  iL peak ratio (spice/anal) = " << peak_ratio
              << "   (gate " << s.tol_peak_ratio_min << ".."
              << s.tol_peak_ratio_max << ")\n";
    CHECK(peak_ratio > s.tol_peak_ratio_min);
    CHECK(peak_ratio < s.tol_peak_ratio_max);

    dump_waveforms_csv(std::string("vienna_") + s.name,
                       analyticalOps[0], simOps[0]);
}

}  // namespace


// First design: TI TIDM-1000-style 10 kW Vienna, 400 V_LL → 800 V/50 kHz.
// Mid-range power and Fsw; M ≈ 0.816, d ≈ 0.184.
TEST_CASE("Vienna reference design PtP — TIDM-1000-like 10 kW "
          "(400V_LL→800V, Fsw=50kHz)",
          "[converter-model][vienna-topology][refdesign][ptp][slow]") {
    ViennaRefSpec s{"TIDM-1000-10kW",
                    400.0,                // V_LL
                    800.0,                // Vdc
                    50e3,                 // Fsw
                    10e3,                 // P
                    0.25,                 // ripple ratio
                    60.0,                 // wall time
                    2.0,                  // Vin sanity %
                    0.30, 0.60,           // loss neg, loss max
                    0.15};                // NRMSE
    run_ptp_gates(s);
}


// Second design: STDES-VIENNARECT 15 kW, 400 V_LL → 800 V/70 kHz.
// Higher power, higher Fsw — same M but larger I_pk, smaller L.
TEST_CASE("Vienna reference design PtP — STDES-VIENNARECT 15 kW "
          "(400V_LL→800V, Fsw=70kHz)",
          "[converter-model][vienna-topology][refdesign][ptp][slow]") {
    ViennaRefSpec s{"STDES-VIENNARECT-15kW",
                    400.0, 800.0, 70e3, 15e3, 0.25,
                    60.0, 2.0, 0.30, 0.60, 0.15};
    run_ptp_gates(s);
}


// Third design: Infineon REF_11KW_PFC_SIC_QD-like 11 kW, 65 kHz. Adds
// coverage at an intermediate power point on a slightly different Fsw.
TEST_CASE("Vienna reference design PtP — Infineon 11 kW SiC "
          "(400V_LL→800V, Fsw=65kHz)",
          "[converter-model][vienna-topology][refdesign][ptp][slow]") {
    ViennaRefSpec s{"Infineon-11kW-SiC",
                    400.0, 800.0, 65e3, 11e3, 0.25,
                    60.0, 2.0, 0.30, 0.60, 0.15};
    run_ptp_gates(s);
}
