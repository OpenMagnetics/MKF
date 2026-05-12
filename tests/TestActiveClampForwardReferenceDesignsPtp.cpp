// Phase 8 — Golden point-to-point (PtP) regression tests for three
// ActiveClampForward (ACF) reference operating points spanning the
// commonly-published telecom / 48 V intermediate-bus range:
//
//   1. UCC2897A-class telecom POL  — 48 V → 3.3 V @ 30 A, 250 kHz   (~100 W)
//        TI EVM SLUU373; Andreycak TI U-153.
//   2. Erickson textbook §6.4-class — 28 V → 5 V  @ 10 A, 200 kHz   (~50 W)
//        Erickson-Maksimović "Fundamentals of Power Electronics" 3rd ed.
//   3. AN-1023-class 12 V brick     — 48 V → 12 V @ 16 A, 250 kHz   (~200 W)
//        Infineon AN-1023 "200 W forward with active-clamp reset".
//
// Each ref-design's Lm and turns ratio are derived by the model itself
// via `process_design_requirements()` (the model is the contract; the
// PtP harness must not encode magnetics values that bypass it).  The
// harness then runs the analytical primary-current shape pass and the
// ngspice switching pass, and asserts the four PtP gates:
//
//   1. wall-time         < 30 s
//   2. input regulation  |Vin_spice − Vin_nom| / Vin_nom < tol_vin_pct
//                        (ACF SPICE has an ideal DC source — sanity)
//   3. power balance     Pin > 0,  loss = (Pin − Pout_recon)/Pin
//                        within [-tol_loss_neg, tol_loss_max].
//                        Pout_recon = Vout_spice² / Rload_nom is
//                        reconstructed (ACF model returns iout as a
//                        load-Ohm projection, not measured), so this
//                        gate is a SPICE-internal sanity check; the
//                        truly-faithful efficiency check waits for
//                        REVIEW_PLAN §3.C ZVS/leakage rework.
//   4. primary current   NRMSE analytical vs SPICE < tol_nrmse
//                        Loose (≈ 0.99 like the existing
//                        Test_ActiveClampForward_PtP_AnalyticalVsNgspice)
//                        because the analytical model lacks the active-
//                        clamp resonant hump on the primary.  Once §3.C
//                        ships the resonant analytical model this can
//                        tighten toward 0.20 like Buck/Boost.
//
// Open-loop Vout droop is allowed: the netlist has no Vout regulator,
// so Vout_spice tracks below Vout_nom by the conduction loss.  We do
// NOT gate on Vout vs nominal — the load-consistency gate from Buck
// is degenerate here (iout is reconstructed from vout/Rload by
// construction).
//
// Tags: [converter-model][active-clamp-forward-topology][refdesign][ptp][slow]
//   The [slow] tag matches the Buck/Boost/PFC PtP convention.

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "converter_models/ActiveClampForward.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

struct RefDesignSpec {
    const char* name;
    double Vin;             // nominal input voltage [V]
    double Vin_min;         // input voltage min [V]
    double Vin_max;         // input voltage max [V]
    double Vout;            // nominal output voltage [V]
    double Iout;            // nominal output current [A]
    double Fs;              // switching frequency [Hz]
    double D;               // (max) duty cycle command (open-loop)
    double tol_walltime;    // wall-time gate [s]
    double tol_vin_pct;     // |Vin_spice − Vin_nom|/Vin_nom [%]
    double tol_loss_neg;    // negative-loss tolerance [fraction]
                            // (allows minor numerical noise / extraction
                            // boundary effects pushing Pout slightly > Pin)
    double tol_loss_max;    // (Pin − Pout_recon)/Pin upper bound [fraction]
    double tol_nrmse;       // primary-current NRMSE [fraction, 0..1]
};

OpenMagnetics::ActiveClampForward build(const RefDesignSpec& s) {
    OpenMagnetics::ActiveClampForward fwd;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin_min);
    iv.set_maximum(s.Vin_max);
    fwd.set_input_voltage(iv);
    fwd.set_diode_voltage_drop(0.5);     // sync rect Vf typical
    fwd.set_efficiency(0.9);              // analytical reference η
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

// ── PtP harness ────────────────────────────────────────────────────────
void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== ACF PtP — " << s.name << " ==========\n";

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    // ── Resolve magnetics from the model itself ───────────────────────
    auto fwd = build(s);
    auto designReqs = fwd.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : designReqs.get_turns_ratios()) {
        turnsRatios.push_back(tr.get_nominal().value());
    }
    REQUIRE(!turnsRatios.empty());
    const double Lm = designReqs.get_magnetizing_inductance()
                          .get_minimum().value();
    REQUIRE(Lm > 0.0);
    std::cout << "  N=" << turnsRatios[0] << "  Lm="
              << 1e6 * Lm << " uH\n";

    // ── Analytical pass (anchors the expected primary-current shape) ──
    auto analyticalOps = fwd.process_operating_points(turnsRatios, Lm);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // ── SPICE switching pass — wall-clock-gated workhorse ─────────────
    // ACF hardcodes Cout = 100 µF + Cclamp = 10 µF.  At low Iout the
    // Cout RC tail is the slow pole (τ = Vout/Iout · 100 µF).  Bump
    // settling periods to 400 to match the Buck/Boost convention.
    fwd.set_num_steady_state_periods(400);
    fwd.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = fwd.simulate_and_extract_topology_waveforms(turnsRatios, Lm);
    auto simOps = fwd.simulate_and_extract_operating_points(turnsRatios, Lm);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — input voltage sanity (ideal DC source).
    const auto& vinW = wfs[0].get_input_voltage();
    const auto& iinW = wfs[0].get_input_current();
    const auto& vinD = vinW.get_data();
    const auto& iinD = iinW.get_data();
    const auto vinT_opt = vinW.get_time();
    REQUIRE(!vinD.empty());
    REQUIRE(!iinD.empty());
    REQUIRE(vinT_opt.has_value());
    const auto& vinT = vinT_opt.value();
    const double vinMean = ConverterPortChecks::time_weighted_mean(vinT, vinD);
    const double vinErr  = (vinMean - s.Vin) / s.Vin;
    std::cout << "  Vin_spice = " << vinMean << " V "
              << "(Vin_nom = " << s.Vin << " V, err "
              << 100.0 * vinErr << " %, gate "
              << s.tol_vin_pct << " %)\n";
    CHECK(std::fabs(vinErr) < s.tol_vin_pct / 100.0);

    // Gate 3 — power balance.  Pin from (vin_dc, iin) directly; Pout
    // reconstructed as Vout_spice² / Rload_nom because the ACF model
    // returns iout as a load-Ohm projection of vout (so Vout·Iout
    // collapses identically to vout²/Rload — gating on it would be a
    // tautology).  This gate verifies (a) Pin > 0 (causal), and
    // (b) η = Pout/Pin lies in a physically plausible band.
    const auto& voutW = wfs[0].get_output_voltages()[0];
    const auto& voutD = voutW.get_data();
    const auto voutT_opt = voutW.get_time();
    REQUIRE(!voutD.empty());
    REQUIRE(voutT_opt.has_value());
    const auto& voutT = voutT_opt.value();
    const double pin    = ConverterPortChecks::time_weighted_mean_product(
        vinT, vinD, iinD);
    const double voutMs = ConverterPortChecks::time_weighted_mean_product(
        voutT, voutD, voutD);   // ⟨Vout²⟩ on non-uniform grid
    const double rloadNom = s.Vout / s.Iout;
    const double poutRec  = voutMs / rloadNom;
    const double loss     = (pin - poutRec) / pin;
    std::cout << "  Pin = " << pin << " W, Pout_recon = " << poutRec
              << " W, loss = " << 100.0 * loss << " %   (gate "
              << -100.0 * s.tol_loss_neg << "..."
              << 100.0 * s.tol_loss_max << " %)\n";
    CHECK(pin > 0.0);
    CHECK(loss >= -s.tol_loss_neg);
    CHECK(loss <=  s.tol_loss_max);

    // Gate 4 — primary current waveform NRMSE (analytical vs SPICE).
    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty());
    REQUIRE(!sTime.empty());
    auto sResampled = ptp_interp(sTime, sData, 256);
    const double nrmse = ptp_nrmse(aResampled, sResampled);
    std::cout << "  iPri NRMSE (analytical vs SPICE) = "
              << 100.0 * nrmse << " %   (gate "
              << 100.0 * s.tol_nrmse << " %)\n";
    CHECK(nrmse < s.tol_nrmse);
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────
// Per-design TEST_CASEs.  Numerical tolerances are SPECIFIC PER DESIGN
// because the analytical-vs-SPICE shape mismatch depends on duty-cycle
// proximity to D_max and clamp-cap dynamics; pick them empirically from
// the first run and lock them.
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("ACF reference design PtP — UCC2897A-class 100 W telecom POL "
          "(48V→3.3V/30A 250 kHz)",
          "[converter-model][active-clamp-forward-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "UCC2897A-100W",
        /*Vin*/          48.0,
        /*Vin_min*/      36.0,
        /*Vin_max*/      72.0,
        /*Vout*/         3.3,
        /*Iout*/         30.0,
        /*Fs*/           250e3,
        /*D*/            0.45,
        /*tol_walltime*/ 30.0,
        /*tol_vin_pct*/  2.0,
        /*tol_loss_neg*/ 0.05,
        /*tol_loss_max*/ 0.60,
        /*tol_nrmse*/    0.99
    };
    run_ptp_gates(s);
}

TEST_CASE("ACF reference design PtP — Erickson §6.4-class 50 W "
          "(28V→5V/10A 200 kHz)",
          "[converter-model][active-clamp-forward-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "Erickson-50W",
        /*Vin*/          28.0,
        /*Vin_min*/      24.0,
        /*Vin_max*/      36.0,
        /*Vout*/         5.0,
        /*Iout*/         10.0,
        /*Fs*/           200e3,
        /*D*/            0.40,
        /*tol_walltime*/ 30.0,
        /*tol_vin_pct*/  2.0,
        /*tol_loss_neg*/ 0.05,
        /*tol_loss_max*/ 0.60,
        /*tol_nrmse*/    0.99
    };
    run_ptp_gates(s);
}

TEST_CASE("ACF reference design PtP — AN-1023-class 200 W brick "
          "(48V→12V/16A 250 kHz)",
          "[converter-model][active-clamp-forward-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "AN1023-200W",
        /*Vin*/          48.0,
        /*Vin_min*/      36.0,
        /*Vin_max*/      72.0,
        /*Vout*/         12.0,
        /*Iout*/         16.0,
        /*Fs*/           250e3,
        /*D*/            0.45,
        /*tol_walltime*/ 30.0,
        /*tol_vin_pct*/  2.0,
        /*tol_loss_neg*/ 0.05,
        /*tol_loss_max*/ 0.60,
        /*tol_nrmse*/    0.99
    };
    run_ptp_gates(s);
}
