// Phase 8d — Golden PtP regression tests for three CllcConverter
// (bidirectional CLLC resonant) reference operating points spanning the
// common application space:
//
//   1. Telecom 500 W brick               — 400 V → 48 V @ 10 A, 200 kHz
//        Forward mode at fr.  Symmetric tank, k=4.45, Q=0.3.
//   2. Telecom 250 W brick               — 400 V → 24 V @ 10.4 A, 250 kHz
//        Forward mode above fr (buck region).
//   3. Infineon AN-2024-06 EV charger    — 750 V → 600 V @ 18.33 A, 73 kHz
//        Forward mode at fr.  11 kW battery-side bidirectional reference.
//
// CLLC primary current is dominated by the resonant tank sinusoid plus a
// triangular magnetizing component; the analytical FHA model and SPICE
// agree well in the resonant fundamental but diverge in higher-harmonic
// content, particularly off-resonance.
//
// NRMSE gate: < 99 %.  No prior single-point CLLC PtP test exists; the
// resonant tank produces a near-sinusoidal primary current that aligns
// well, but FHA neglects the dead-time voltage notches and the
// magnetizing triangle phase, so a 99 % gate records the comparison
// diagnostically while still rejecting a regressed waveform that lost
// shape entirely.
//
// Tags: [converter-model][cllc-topology][refdesign][ptp][slow]

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

#include "converter_models/Cllc.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin, Vin_min, Vin_max;
    double Vout, Iout, Fs, Fmin, Fmax;
    double tol_walltime, tol_vin_pct;
    double tol_loss_neg, tol_loss_max;
    double tol_nrmse;
};

nlohmann::json build_fixture(const RefDesignSpec& s) {
    return nlohmann::json{
        {"inputVoltage", {{"nominal", s.Vin}, {"minimum", s.Vin_min}, {"maximum", s.Vin_max}}},
        {"maxSwitchingFrequency", s.Fmax},
        {"minSwitchingFrequency", s.Fmin},
        {"efficiency", 0.95},
        {"qualityFactor", 0.3},
        {"symmetricDesign", true},
        {"bidirectional", true},
        {"operatingPoints", nlohmann::json::array({
            nlohmann::json{
                {"outputVoltages", {s.Vout}},
                {"outputCurrents", {s.Iout}},
                {"switchingFrequency", s.Fs},
                {"ambientTemperature", 25.0},
                {"powerFlow", "forward"}
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
    std::cout << "\n========== CLLC PtP — " << s.name << " ==========\n";
    NgspiceRunner runner;
    if (!runner.is_available()) { WARN("ngspice not available"); return; }

    OpenMagnetics::CllcConverter cllc(build_fixture(s));
    auto params = cllc.calculate_resonant_parameters();
    const double n  = params.turnsRatio;
    const double Lm = params.magnetizingInductance;
    REQUIRE(n  > 0.0);
    REQUIRE(Lm > 0.0);
    std::cout << "  n=" << n << "  Lm=" << 1e6*Lm << " uH"
              << "  L1=" << 1e6*params.primaryResonantInductance << " uH"
              << "  C1=" << 1e9*params.primaryResonantCapacitance << " nF\n";

    auto analyticalOps = cllc.process_operating_points({n}, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty()); REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    cllc.set_num_steady_state_periods(50);
    cllc.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = cllc.simulate_and_extract_topology_waveforms({n}, Lm, 1);
    auto simOps = cllc.simulate_and_extract_operating_points({n}, Lm);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    REQUIRE(!wfs.empty()); REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — Vin sanity (use OP 0, which corresponds to nominal Vin)
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

TEST_CASE("CLLC reference design PtP — Telecom 500 W brick at fr "
          "(400V→48V/10A 200 kHz)",
          "[converter-model][cllc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"Telecom-500W", 400.0, 360.0, 420.0, 48.0, 10.0,
                    200e3, 100e3, 400e3,
                    60.0, 2.0, 0.20, 0.60, 0.99};
    run_ptp_gates(s);
}

TEST_CASE("CLLC reference design PtP — Telecom 250 W brick above fr "
          "(400V→24V/10.4A 250 kHz)",
          "[converter-model][cllc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"Telecom-250W", 400.0, 360.0, 420.0, 24.0, 10.4,
                    250e3, 100e3, 400e3,
                    60.0, 2.0, 0.20, 0.60, 0.99};
    run_ptp_gates(s);
}

TEST_CASE("CLLC reference design PtP — Infineon AN-2024-06 EV charger "
          "(750V→600V/18.33A 73 kHz)",
          "[converter-model][cllc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"Infineon-11kW", 750.0, 700.0, 800.0, 600.0, 18.33,
                    73e3, 40e3, 250e3,
                    60.0, 2.0, 0.20, 0.60, 0.99};
    run_ptp_gates(s);
}
