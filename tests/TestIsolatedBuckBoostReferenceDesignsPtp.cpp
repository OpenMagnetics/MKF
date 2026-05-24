// Phase 8l — Golden PtP regression tests for three IsolatedBuckBoost
// (isolated flyback-class buck-boost) reference designs:
//   MAX17498A 12V→5V, LM5180 SNVA827 24V→12V, Erickson Sec.6.1.3 24V→36V.
//
// Mirrors the canonical Buck/Boost/ACF/Flyback/PushPull/IsolatedBuck PtP
// harness pattern (4 gates):
//
//   wall_time              < tol_walltime
//   load consistency       |Vout_spice/Iout_spice − Rload_nom|/Rload_nom
//                          < tol_rload_pct       (primary regulated output)
//   Pout balance           |Pout_total_spice − Pout_total_nom|/Pout_total_nom
//                          < tol_pout_pct        (Pout sums primary + bias)
//                          Mirrors the IsolatedBuck PtP pattern: Pin-vs-Pout
//                          comparison is unreliable for flyback-class
//                          topologies because the single-period extraction
//                          window only sees Vin·Iin during MOSFET-on time
//                          while output ports integrate V·I across the full
//                          period.  Comparing Pout_spice against the
//                          nominal V·I design point cross-checks the open-
//                          loop operating point.  For low-Vout designs
//                          (e.g. 12→5 V), real ~0.6 V diode drops eat
//                          ~12 % of volt-seconds in the buck-boost transfer
//                          function and shift the settled Vout below
//                          nominal; the tol_pout gate is sized per-design.
//   primary-current NRMSE  analytical vs SPICE < tol_nrmse
//
// Specs reuse the three reference designs already characterised in
// TestTopologyIsolatedBuckBoost.cpp (kIBBRefDesign1..3); duplicated here
// as struct literals so each PtP test is self-describing.
//
// Tags: [converter-model][isolated-buck-boost-topology][refdesign][ptp][slow]

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "converter_models/IsolatedBuckBoost.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"
#include "PtpHelpers.h"
#include "WaveformDumpHelpers.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

using namespace OpenMagnetics::Testing;

struct RefDesignSpec {
    const char* name;
    double Vin;             // input voltage [V]
    double Vout;            // primary (regulated) output voltage [V]
    double Iout;            // primary output current [A]
    double Fs;              // switching frequency [Hz]
    double Lpri;            // primary (magnetizing) inductance [H]
    double n;               // turns ratio Np:Ns
    double Iout2;           // isolated bias output current [A]
    double tol_walltime;    // wall-time gate [s]
    double tol_rload_pct;   // |Vout_spice/Iout_spice − Rload_nom|/Rload_nom [%]
    double tol_pout_pct;    // |Pout_total_spice − Pout_total_nom|/Pout_total_nom [%]
    double tol_nrmse;       // primary-current NRMSE [fraction, 0..1]
};

OpenMagnetics::IsolatedBuckBoost build(const RefDesignSpec& s) {
    OpenMagnetics::IsolatedBuckBoost ibb;
    DimensionWithTolerance iv;
    iv.set_nominal(s.Vin);
    iv.set_minimum(s.Vin * 0.95);
    iv.set_maximum(s.Vin * 1.05);
    ibb.set_input_voltage(iv);
    ibb.set_diode_voltage_drop(0.0);   // ideal diode for closed-form anchor
    ibb.set_efficiency(1.0);            // lossless analytical reference
    ibb.set_current_ripple_ratio(0.3);
    IsolatedBuckBoostOperatingPoint op;
    op.set_output_voltages({s.Vout, s.Vout / s.n});  // primary + isolated bias
    op.set_output_currents({s.Iout, s.Iout2});
    op.set_switching_frequency(s.Fs);
    op.set_ambient_temperature(25.0);
    ibb.set_operating_points(std::vector<IsolatedBuckBoostOperatingPoint>{op});
    return ibb;
}

void run_ptp_gates(const RefDesignSpec& s) {
    std::cout << "\n========== IsolatedBuckBoost PtP — " << s.name << " ==========\n";

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not available — skipping " << s.name);
        return;
    }

    auto ibb = build(s);
    const std::vector<double> tr{s.n};

    // ── Analytical pass — anchors primary-current shape ───────────────
    auto analyticalOps = ibb.process_operating_points(tr, s.Lpri);
    REQUIRE(!analyticalOps.empty());
    auto aTime = ptp_current_time(analyticalOps[0], 0);
    auto aData = ptp_current(analyticalOps[0], 0);
    REQUIRE(!aData.empty());
    REQUIRE(!aTime.empty());
    auto aResampled = ptp_interp(aTime, aData, 256);

    // ── SPICE switching pass ──────────────────────────────────────────
    // Flyback-class converters with multi-output rectification need long
    // settling. Output caps are pre-charged via `IC=Vout` (see
    // IsolatedBuckBoost.cpp:~531, ~544) and `.nodeset v(vout*)=Vout`
    // (~582), so cap voltage doesn't have to charge from 0. Despite
    // that, Gate 3 (Pout vs nominal) still requires 1200+ periods:
    //
    //   600  periods  → Pout err 33% / 9% / 21%   (MAX / LM5180 / Erickson)
    //   900  periods  → Pout err 33% / 9% / 21%   (identical — not settled)
    //   1200 periods  → all three within ±10-30% gates
    //
    // The plateau between 600 and 900 indicates this isn't an RC tail
    // (which would decay monotonically). It's a slow LC ring between the
    // primary inductor and the multi-output cap network that doesn't
    // damp until ~1000+ periods. IC=/.nodeset can't shortcut a resonant
    // mode the netlist itself excites.
    //
    // If this needs to be faster, the right fix is in the SPICE netlist:
    // add a small series-R damping element or shrink Cout (reducing the
    // RC time constant), not raise settling. See
    // src/converter_models/IsolatedBuckBoost.cpp:541-545.
    ibb.set_num_steady_state_periods(1200);
    ibb.set_num_periods_to_extract(1);

    const auto t0 = std::chrono::steady_clock::now();
    auto wfs    = ibb.simulate_and_extract_topology_waveforms(tr, s.Lpri, 1);
    auto simOps = ibb.simulate_and_extract_operating_points(tr, s.Lpri);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    REQUIRE(!wfs.empty());
    REQUIRE(!simOps.empty());

    std::cout << "  wall_time = " << wallTime << " s   (gate "
              << s.tol_walltime << " s)\n";

    // Gate 1 — wall time
    CHECK(wallTime < s.tol_walltime);

    // Gate 2 — load consistency on the primary (regulated) output.
    REQUIRE(wfs[0].get_output_voltages().size() >= 1);
    REQUIRE(wfs[0].get_output_currents().size() >= 1);
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

    // Gate 3 — Pout(total) cross-check against nominal V·I design point.
    // See header comment for why Pin-vs-Pout is unreliable here.  Vout is
    // signed (Fly-Buck-Boost inverts) so Pout = V·I is positive whether
    // measured against +5 V/+0.5 A nominal or −5 V/−0.5 A settled.
    double pout = ConverterPortChecks::time_weighted_mean_product(
        voutT, voutD, ioutD);
    if (wfs[0].get_output_voltages().size() >= 2 &&
        wfs[0].get_output_currents().size() >= 2) {
        const auto& v2W = wfs[0].get_output_voltages()[1];
        const auto& i2W = wfs[0].get_output_currents()[1];
        const auto& v2D = v2W.get_data();
        const auto& i2D = i2W.get_data();
        const auto v2T_opt = v2W.get_time();
        if (!v2D.empty() && !i2D.empty() && v2T_opt.has_value()) {
            pout += ConverterPortChecks::time_weighted_mean_product(
                v2T_opt.value(), v2D, i2D);
        }
    }
    const double poutNomTotal = s.Vout * s.Iout + (s.Vout / s.n) * s.Iout2;
    const double poutErr = (pout - poutNomTotal) / poutNomTotal;
    std::cout << "  Pout(total)_spice = " << pout << " W, "
              << "Pout(total)_nom = " << poutNomTotal << " W, "
              << "err = " << 100.0 * poutErr << " %   (gate ±"
              << s.tol_pout_pct << " %)\n";
    CHECK(std::fabs(poutErr) <= s.tol_pout_pct / 100.0);

    // Gate 4 — primary-current shape NRMSE.  η=1, Vd=0 ideal cases expose
    // resonant ringing on primary current at switching transitions that
    // the analytical model approximates as ideal triangles; gate is sized
    // accordingly (see existing PtP rationale in
    // TestTopologyIsolatedBuckBoost.cpp:543-556).
    auto sTime = ptp_current_time(simOps[0], 0);
    auto sData = ptp_current(simOps[0], 0);
    REQUIRE(!sData.empty());
    REQUIRE(!sTime.empty());
    auto sResampled = ptp_interp(sTime, sData, 256);
    const double nrmse = ptp_nrmse(aResampled, sResampled, 256);
    std::cout << "  iPri NRMSE (analytical vs SPICE) = "
              << 100.0 * nrmse << " %   (gate "
              << 100.0 * s.tol_nrmse << " %)\n";
    CHECK(nrmse < s.tol_nrmse);

    OpenMagnetics::Testing::dump_waveforms_csv(
        std::string("isolated_buck_boost_") + s.name, analyticalOps[0], simOps[0]);
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────
// Per-design TEST_CASEs.  Specs match kIBBRefDesign<1,2,3> in
// TestTopologyIsolatedBuckBoost.cpp (canonical source).
// ─────────────────────────────────────────────────────────────────────────

TEST_CASE("IsolatedBuckBoost reference design PtP — MAX17498A (12V→5V)",
          "[converter-model][isolated-buck-boost-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "MAX17498A_12Vin_5Vout",
        /*Vin*/          12.0,
        /*Vout*/         5.0,
        /*Iout*/         0.50,
        /*Fs*/           250e3,
        /*Lpri*/         100e-6,
        /*n*/            1.5,
        /*Iout2*/        0.10,
        /*tol_walltime*/ 14.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_pout_pct*/  30.0,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("IsolatedBuckBoost reference design PtP — TI LM5180 SNVA827 (24V→12V)",
          "[converter-model][isolated-buck-boost-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "LM5180_SNVA827_24Vin_12Vout",
        /*Vin*/          24.0,
        /*Vout*/         12.0,
        /*Iout*/         1.00,
        /*Fs*/           100e3,
        /*Lpri*/         47e-6,
        /*n*/            1.0,
        /*Iout2*/        0.10,
        /*tol_walltime*/ 24.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_pout_pct*/  15.0,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}

TEST_CASE("IsolatedBuckBoost reference design PtP — Erickson §6.1.3 (24V→36V step-up)",
          "[converter-model][isolated-buck-boost-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "Erickson_Sec6p1p3_24Vin_36Vout",
        /*Vin*/          24.0,
        /*Vout*/         36.0,
        /*Iout*/         0.50,
        /*Fs*/           100e3,
        /*Lpri*/         220e-6,
        /*n*/            2.0,
        /*Iout2*/        0.10,
        /*tol_walltime*/ 34.0,
        /*tol_rload_pct*/ 0.01,
        /*tol_pout_pct*/  10.0,
        /*tol_nrmse*/    0.15
    };
    run_ptp_gates(s);
}
