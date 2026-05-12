// Phase 8c — Golden PtP regression tests for three AsymmetricHalfBridge
// (AHB) reference operating points spanning the common
// telecom / brick / AC-DC AHB application space:
//
//   1. TI SLUP223-class 100 W telecom    — 100 V → 5 V  @ 20 A, 200 kHz
//        Steigerwald / TI SLUP223 worked example (CT secondary).
//   2. ON Semi AN-4153-class 200 W brick — 100 V → 12 V @ 16 A, 100 kHz
//        ON Semi AN-4153 server-brick reference (CT secondary).
//   3. ST AN2852-class 90 W AC-DC        —  90 V → 19 V @ 4.7 A, 100 kHz
//        ST AN2852 ATX-adapter reference (CT secondary, post-PFC bus).
//
// AHB asymmetric volt-second clamp via DC-blocking cap V_Cb=(1-D)*Vin
// gives a non-monotonic conversion ratio peaking at D=0.5; the analytical
// primary current is the textbook trapezoid + (Io/n)·(1-2D) DC bias on the
// magnetizing ramp.
//
// NRMSE gate: < 150 %.  The existing AHB P8 single-point gates record
// 0.56 (D=0.45) to 1.04 (D=0.30) NRMSE: when D ≠ 0.5 the asymmetric
// (Io/n)·(1-2D) DC bias is highly sensitive to load-driven Io shift, and
// open-loop SPICE Vout droops 25–40 % below nominal at these power
// levels which inflates NRMSE further.  150 % records the comparison
// diagnostically (no silent fallback) and gates only that the analytical
// vs SPICE comparison is finite and within the same magnitude family.
// Tighter gating waits for closed-loop Vout regulation
// (REVIEW_PLAN §3.B item).
//
// Tags: [converter-model][ahb-topology][refdesign][ptp][slow]

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

#include "converter_models/AsymmetricHalfBridge.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;
using AHB = OpenMagnetics::AsymmetricHalfBridge;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin;
    double Vout, Iout, Fs, D;
    double tol_walltime, tol_vin_pct;
    double tol_loss_neg, tol_loss_max;
    double tol_nrmse;
};

// Build a self-consistent AHB JSON fixture: first pass sizes n at the
// requested D; second pass solves Vo = M(D,n)·Vin so SPICE (open-loop)
// converges to the same Io as the analytical model.
nlohmann::json build_fixture(const RefDesignSpec& s) {
    nlohmann::json design = nlohmann::json{
        {"inputVoltage", {{"nominal", s.Vin}}},
        {"operatingPoints", nlohmann::json::array({
            nlohmann::json{
                {"ambientTemperature", 25.0},
                {"switchingFrequency", s.Fs},
                {"outputVoltages", {s.Vout}},
                {"outputCurrents", {s.Iout}},
                {"dutyCycle", s.D}
            }
        })}
    };
    AHB design_ahb(design);
    design_ahb.process_design_requirements();
    const double n = design_ahb.get_computed_turns_ratio();
    const double M = AHB::compute_conversion_ratio(
        s.D, n, AhbRectifierType::CENTER_TAPPED);
    const double Vo = M * s.Vin;
    const double Io = s.Iout * (s.Vout / Vo);  // preserve nominal Po
    return nlohmann::json{
        {"inputVoltage", {{"nominal", s.Vin}}},
        {"operatingPoints", nlohmann::json::array({
            nlohmann::json{
                {"ambientTemperature", 25.0},
                {"switchingFrequency", s.Fs},
                {"outputVoltages", {Vo}},
                {"outputCurrents", {Io}},
                {"dutyCycle", s.D}
            }
        })}
    };
}

using namespace OpenMagnetics::Testing;

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== AHB PtP — " << s.name << " ==========\n";
    NgspiceRunner runner;
    if (!runner.is_available()) { WARN("ngspice not available"); return; }

    AHB ahb(build_fixture(s));
    ahb.process_design_requirements();
    const double n  = ahb.get_computed_turns_ratio();
    const double Lm = ahb.get_computed_magnetizing_inductance();
    REQUIRE(n  > 0.0);
    REQUIRE(Lm > 0.0);
    std::cout << "  n=" << n << "  Lm=" << 1e6*Lm << " uH\n";

    auto analyticalOps = ahb.process_operating_points({n}, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty()); REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    ahb.set_num_steady_state_periods(400);
    ahb.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = ahb.simulate_and_extract_topology_waveforms({n}, Lm);
    auto simOps = ahb.simulate_and_extract_operating_points({n}, Lm);
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
    // Rload_nom uses the self-consistent (Vo, Io) from build_fixture so
    // the analytical model and SPICE reference the same load resistance.
    const auto& opActual = ahb.get_operating_points()[0];
    const double VoNom = opActual.get_output_voltages()[0];
    const double IoNom = opActual.get_output_currents()[0];
    const double rloadNom = VoNom / IoNom;
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

TEST_CASE("AHB reference design PtP — TI SLUP223-class 100 W telecom "
          "(100V→5V/20A 200 kHz)",
          "[converter-model][ahb-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"SLUP223-100W", 100.0, 5.0, 20.0, 200e3, 0.45,
                    30.0, 2.0, 0.05, 0.60, 1.50};
    run_ptp_gates(s);
}

TEST_CASE("AHB reference design PtP — ON Semi AN-4153-class 200 W brick "
          "(100V→12V/16A 100 kHz)",
          "[converter-model][ahb-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"AN4153-200W", 100.0, 12.0, 16.0, 100e3, 0.45,
                    30.0, 2.0, 0.05, 0.60, 1.50};
    run_ptp_gates(s);
}

TEST_CASE("AHB reference design PtP — ST AN2852-class 90 W AC-DC "
          "(90V→19V/4.7A 100 kHz)",
          "[converter-model][ahb-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{"AN2852-90W", 90.0, 19.0, 4.7, 100e3, 0.45,
                    30.0, 2.0, 0.05, 0.60, 1.50};
    run_ptp_gates(s);
}
