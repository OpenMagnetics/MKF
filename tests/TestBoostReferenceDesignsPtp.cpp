// Phase 7 — Golden point-to-point (PtP) regression tests for the three
// Boost reference designs (TPS61089EVM-742 18 W, TPS61178EVM-792 32 W,
// LM5122EVM-1PH 108 W).
//
// Mirror image of `TestBuckReferenceDesignsPtp.cpp`, applied to the
// Boost converter.  See that file's header for the rationale; the
// per-design Values + raw analytical-vs-SPICE NRMSE checks already live
// in `TestTopologyBoost.cpp` (the `Test_Boost_RefDesign<N>_*` cases),
// and this file consolidates the four PFC-style acceptance gates under
// the canonical `[refdesign][ptp][slow]` tag set so a single suite
// filter exercises the full per-design contract for every topology.
//
// Each test:
//   1. Builds a Boost model with the ref-design specs (synchronous,
//      η=1, Vd=0 — matches the EVM topology, lossless analytical ref).
//   2. Runs the analytical pass to anchor the inductor-current shape.
//   3. Runs the ngspice switching netlist for one steady-state period
//      (after 400 settling periods to drain the Cout RC tail).
//   4. Asserts the four PtP gates:
//        wall-time         < 30 s          (Boost designs typically <1 s)
//        load consistency  |Vout_spice/Iout_spice − Rload_nom|/Rload_nom
//                          within tol_rload_pct (open-loop SPICE; load
//                          ratio is fixed even if absolute Vout drifts)
//        power balance     0 ≤ (Pin − Pout_spice)/Pin ≤ tol_loss_max
//        inductor NRMSE    analytical vs SPICE shape NRMSE < 15 %
//
//   Gates 2 and 3 are SPICE-internal (self-consistent), not vs nominal:
//   the Boost SPICE netlist is open-loop, so Vout settles below nominal
//   in proportion to conduction/switching loss.  See TestBuckRefDesignsPtp
//   for the full rationale.
//
// Tags: [converter-model][boost-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "converter_models/Boost.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin;            // input voltage [V]
    double Vout;           // output voltage [V]
    double Iout;           // output current [A]
    double Fs;             // switching frequency [Hz]
    double L;              // inductor value [H]
    double tol_walltime;   // wall-time gate [s]
    double tol_rload_pct;  // |Vout_spice/Iout_spice − Rload_nom|/Rload_nom [%]
    double tol_loss_max;   // (Pin − Pout_spice)/Pin upper bound [fraction]
    double tol_nrmse;      // inductor-current shape NRMSE [fraction, 0..1]
};

OpenMagnetics::Boost build(const RefDesignSpec& s) {
    OpenMagnetics::Boost b;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    b.set_input_voltage(iv);
    b.set_diode_voltage_drop(0.0);   // synchronous boost (matches EVMs)
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

// ── Local copies of the ptp_* helpers used by TestTopologyBoost.cpp ───
// Inlined so this PtP file is self-contained.
std::vector<double> ptp_interp(const std::vector<double>& t,
                                const std::vector<double>& d,
                                int N) {
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

double ptp_nrmse(const std::vector<double>& ref,
                 const std::vector<double>& sim) {
    int N = (int)ref.size();
    double ref_mean = 0.0, sim_mean = 0.0;
    for (int i = 0; i < N; ++i) { ref_mean += ref[i]; sim_mean += sim[i]; }
    ref_mean /= N; sim_mean /= N;
    std::vector<double> r(N), s(N);
    double rAC = 0.0, sAC = 0.0;
    for (int i = 0; i < N; ++i) {
        r[i] = ref[i] - ref_mean; s[i] = sim[i] - sim_mean;
        rAC += r[i] * r[i]; sAC += s[i] * s[i];
    }
    rAC = std::sqrt(rAC / N); sAC = std::sqrt(sAC / N);
    if (rAC < 1e-10 || sAC < 1e-10) return 1.0;
    for (int i = 0; i < N; ++i) { r[i] /= rAC; s[i] /= sAC; }
    int ns = std::min(N, 64);
    double best = std::numeric_limits<double>::max();
    for (int ss = 0; ss < ns; ++ss) {
        int sh = ss * N / ns;
        double ssd = 0.0;
        for (int k = 0; k < N; ++k) {
            double e = r[k] - s[(k + sh) % N];
            ssd += e * e;
        }
        if (ssd < best) best = ssd;
    }
    return std::sqrt(best / N);
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

// ── PtP harness ────────────────────────────────────────────────────────
void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== Boost PtP — " << s.name << " ==========\n";

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    auto boost = build(s);
    auto analyticalOps = boost.process_operating_points(std::vector<double>{}, s.L);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // 400 settling periods to drain the Cout RC tail (same rationale as
    // assert_boost_refdesign_ptp in TestTopologyBoost.cpp).
    boost.set_num_steady_state_periods(400);
    boost.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = boost.simulate_and_extract_topology_waveforms(s.L);
    auto simOps = boost.simulate_and_extract_operating_points(s.L);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency.  Boost SPICE is open-loop (no Vout
    // regulator); Vout settles below nominal in proportion to the
    // conduction/switching loss.  But Rload = Vout_nom/Iout_nom is
    // fixed, so Vout_spice / Iout_spice must reproduce that ratio.
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

    // Gate 4 — inductor-current waveform NRMSE
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


TEST_CASE("Boost reference design PtP — TPS61089EVM-742 (18 W)",
          "[converter-model][boost-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "TPS61089EVM-742",
        /*Vin*/          5.0,
        /*Vout*/         9.0,
        /*Iout*/         2.0,
        /*Fs*/           400e3,
        /*L*/            4.7e-6,
        /*tol_walltime*/ 30.0,
        /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.50,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("Boost reference design PtP — TPS61178EVM-792 (32 W)",
          "[converter-model][boost-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "TPS61178EVM-792",
        /*Vin*/          7.2,
        /*Vout*/         16.0,
        /*Iout*/         2.0,
        /*Fs*/           300e3,
        /*L*/            6.8e-6,
        /*tol_walltime*/ 30.0,
        /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.50,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("Boost reference design PtP — LM5122EVM-1PH (108 W)",
          "[converter-model][boost-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "LM5122EVM-1PH",
        /*Vin*/          12.0,
        /*Vout*/         24.0,
        /*Iout*/         4.5,
        /*Fs*/           250e3,
        /*L*/            10e-6,
        /*tol_walltime*/ 30.0,
        /*tol_rload_pct*/ 5.0,
        /*tol_loss_max*/ 0.50,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}
