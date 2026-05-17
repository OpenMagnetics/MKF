// Phase 5 — Golden point-to-point (PtP) regression tests for the three
// boost-PFC reference designs (NCP1654 100 W, UCC28180 360 W, L4981 1000 W).
//
// Each test:
//   1. Builds a PowerFactorCorrection model with the ref-design specs.
//   2. Drives the switching netlist through ngspice.
//   3. Asserts the analytical-vs-SPICE acceptance gates from §3.D Phase 5:
//        wall-time      < 60 s
//        bus regulation within ±6 % of nominal Vbus
//        Pin balance    within ±5 % of nominal Pout
//        inductor envelope NRMSE within per-design tolerance
//          (NCP1654 < 10 %, UCC28180 < 30 %, L4981 < 60 %)
//      The envelope tolerances reflect measured Phase-4 reality:
//      mid/high-power designs suffer the well-known PFC zero-crossing
//      tracking distortion (post-ZC inductor overshoot + pre-ZC saturation
//      of the duty command at D = 1).  A future per-design feed-forward
//      pass can tighten these.
//
// Tags: [converter-model][pfc-topology][refdesign][ptp][slow]
//   The [slow] tag keeps these out of the snappy default suite — they
//   each take several seconds of ngspice wall-time.

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "converter_models/PowerFactorCorrection.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"

using namespace OpenMagnetics;
using namespace MAS;

namespace {

struct RefDesignSpec {
    const char* name;
    double vrms;          // RMS line voltage [V]
    double vbus;          // nominal bus voltage [V]
    double pout;          // nominal output power [W]
    double fsw;           // switching frequency [Hz]
    double cbus;          // bulk capacitance [F]
    double L;             // boost inductance [H]
    int    cycles;        // # of full line cycles to simulate
    double tol_vbus_pct;  // bus regulation tolerance [%]
    double tol_pin_pct;   // power balance tolerance [%]
    double tol_envelope;  // inductor envelope NRMSE tolerance [fraction, 0..1]
    double tol_walltime;  // wall-time gate [s]
};

OpenMagnetics::PowerFactorCorrection build(const RefDesignSpec& s) {
    OpenMagnetics::PowerFactorCorrection pfc;
    DimensionWithTolerance iv;
    iv.set_nominal(s.vrms);
    iv.set_minimum(s.vrms);
    iv.set_maximum(s.vrms);
    pfc.set_input_voltage(iv);
    pfc.set_output_voltage(s.vbus);
    pfc.set_output_power(s.pout);
    pfc.set_switching_frequency(s.fsw);
    pfc.set_line_frequency(50.0);
    pfc.set_efficiency(1.0);
    pfc.set_diode_voltage_drop(0.0);
    pfc.set_current_ripple_ratio(0.3);
    pfc.set_bulk_capacitance(s.cbus);
    pfc.set_mode(PfcModes::CONTINUOUS_CONDUCTION_MODE);
    pfc.set_ambient_temperature(25.0);
    return pfc;
}

double interp(const std::vector<double>& t,
              const std::vector<double>& v,
              double t0) {
    if (t0 <= t.front()) return v.front();
    if (t0 >= t.back())  return v.back();
    auto it = std::lower_bound(t.begin(), t.end(), t0);
    size_t i = std::distance(t.begin(), it);
    if (i == 0) return v.front();
    double a = (t0 - t[i-1]) / (t[i] - t[i-1]);
    return v[i-1] + a * (v[i] - v[i-1]);
}

// Trapezoidal integration of v(t) over [t0,t1] on the simulator's own
// non-uniform time grid (subset to the window).
double mean_over(const std::vector<double>& t,
                 const std::vector<double>& v,
                 double t0, double t1) {
    double area = 0.0;
    double tprev = t0, vprev = interp(t, v, t0);
    for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] <= t0) continue;
        if (t[i] >= t1) break;
        area += 0.5 * (v[i] + vprev) * (t[i] - tprev);
        tprev = t[i]; vprev = v[i];
    }
    double vend = interp(t, v, t1);
    area += 0.5 * (vend + vprev) * (t1 - tprev);
    return area / (t1 - t0);
}

double mean_product(const std::vector<double>& t,
                    const std::vector<double>& a,
                    const std::vector<double>& b,
                    double t0, double t1) {
    double area = 0.0;
    double tprev = t0;
    double aprev = interp(t, a, t0);
    double bprev = interp(t, b, t0);
    for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] <= t0) continue;
        if (t[i] >= t1) break;
        double pcur  = a[i] * b[i];
        double pprev = aprev * bprev;
        area += 0.5 * (pcur + pprev) * (t[i] - tprev);
        tprev = t[i]; aprev = a[i]; bprev = b[i];
    }
    double aend = interp(t, a, t1);
    double bend = interp(t, b, t1);
    area += 0.5 * (aend*bend + aprev*bprev) * (t1 - tprev);
    return area / (t1 - t0);
}

// Run ngspice on the ref-design netlist and assert the §3.D Phase 5 PtP
// acceptance gates.
void run_ptp_gates(const RefDesignSpec& s) {
    INFO("Reference design: " << s.name);

    const auto tBuild0 = std::chrono::steady_clock::now();
    auto pfc = build(s);
    const std::string netlist = pfc.generate_ngspice_switching_circuit(
        s.L, s.cycles);
    const double tBuild = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - tBuild0).count();

    NgspiceRunner runner;
    if (!runner.is_available()) {
        WARN("ngspice not installed on this host; skipping PtP test");
        return;
    }

    SimulationConfig cfg;
    cfg.frequency        = 50.0;
    cfg.extractOnePeriod = false;
    cfg.numberOfPeriods  = static_cast<size_t>(s.cycles);
    cfg.keepTempFiles    = false;
    cfg.timeout          = 180.0;

    const auto t0 = std::chrono::steady_clock::now();
    const auto sim = runner.run_simulation(netlist, cfg);
    const double wallTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    std::cout << "[PFC PtP " << s.name << "] tBuild=" << tBuild
              << "s tSim=" << wallTime << "s\n";

    const auto tPostSim0 = std::chrono::steady_clock::now();

    INFO("wall_time=" << wallTime << " s   success=" << sim.success
         << "   waveforms=" << sim.waveformNames.size());
    REQUIRE(sim.success);
    REQUIRE(wallTime < s.tol_walltime);

    // Resolve waveforms by name (case-insensitive substring match).
    auto find_by = [&](const std::string& key) -> const Waveform* {
        for (size_t k = 0; k < sim.waveformNames.size(); ++k) {
            std::string n = sim.waveformNames[k];
            std::transform(n.begin(), n.end(), n.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            if (n.find(key) != std::string::npos) return &sim.waveforms[k];
        }
        return nullptr;
    };
    const Waveform* wfTime = find_by("time");
    const Waveform* wfVbus = find_by("vbus");
    const Waveform* wfVin  = find_by("vin_rect");
    const Waveform* wfIL   = find_by("vl_sense");
    REQUIRE(wfTime != nullptr);
    REQUIRE(wfVbus != nullptr);
    REQUIRE(wfVin  != nullptr);
    REQUIRE(wfIL   != nullptr);

    const auto tvec_opt = wfTime->get_time();
    REQUIRE(tvec_opt.has_value());
    const auto tvec = tvec_opt.value();
    const auto vbus = wfVbus->get_data();
    const auto vin  = wfVin->get_data();
    const auto iL   = wfIL->get_data();

    const double Tline = 1.0 / 50.0;
    const double t_a   = tvec.back() - Tline;
    const double t_b   = tvec.back();

    // ---- Gate 1: bus regulation -----------------------------------------
    const double vbus_mean = mean_over(tvec, vbus, t_a, t_b);
    const double vbus_err  = 100.0 * (vbus_mean - s.vbus) / s.vbus;
    INFO("vbus_mean=" << vbus_mean << " V (target " << s.vbus
         << ", err " << vbus_err << " %, tol ±" << s.tol_vbus_pct << " %)");
    REQUIRE(std::fabs(vbus_err) < s.tol_vbus_pct);

    // ---- Gate 2: power balance ------------------------------------------
    const double pin     = mean_product(tvec, vin, iL, t_a, t_b);
    const double pin_err = 100.0 * (pin - s.pout) / s.pout;
    INFO("Pin=" << pin << " W (target " << s.pout
         << ", err " << pin_err << " %, tol ±" << s.tol_pin_pct << " %)");
    REQUIRE(std::fabs(pin_err) < s.tol_pin_pct);

    // ---- Gate 3: inductor-current envelope NRMSE ------------------------
    // Compare a one-Tsw sliding mean of iL against the ideal Iin_pk·|sin| line
    // envelope.  The sliding mean strips the switching ripple so this metric
    // exposes line-frequency tracking error only.
    //
    // Implementation: precompute the trapezoidal cumulative integral
    // F[i] = ∫_{tvec[0]}^{tvec[i]} iL(τ) dτ once in O(N); each sliding-mean
    // query then reduces to (F_at(b) - F_at(a)) / (b - a) where F_at(t) for
    // t inside segment [tvec[k-1], tvec[k]] is F[k-1] + 0.5·(iL[k-1] +
    // iL_lin(t))·(t − tvec[k-1]).  Window-edge lookup is O(log N) via
    // binary search.  Replaces the previous O(N²) inner mean_over loop.
    const double Iin_rms = s.pout / s.vrms;
    const double Iin_pk  = std::sqrt(2.0) * Iin_rms;
    const double Tsw     = 1.0 / s.fsw;
    std::vector<double> Fcum(tvec.size(), 0.0);
    for (size_t i = 1; i < tvec.size(); ++i) {
        Fcum[i] = Fcum[i-1] + 0.5 * (iL[i-1] + iL[i]) * (tvec[i] - tvec[i-1]);
    }
    auto F_at = [&](double t) -> double {
        if (t <= tvec.front()) return 0.0;
        if (t >= tvec.back())  return Fcum.back();
        auto it = std::upper_bound(tvec.begin(), tvec.end(), t);
        size_t k = std::distance(tvec.begin(), it);   // tvec[k-1] < t <= tvec[k]
        const double dt   = tvec[k] - tvec[k-1];
        const double frac = (t - tvec[k-1]) / dt;
        const double iLt  = iL[k-1] + frac * (iL[k] - iL[k-1]);
        return Fcum[k-1] + 0.5 * (iL[k-1] + iLt) * (t - tvec[k-1]);
    };
    auto sliding_mean = [&](double tc) {
        double a = std::max(tc - 0.5*Tsw, tvec.front());
        double b = std::min(tc + 0.5*Tsw, tvec.back());
        return (F_at(b) - F_at(a)) / (b - a);
    };
    double sse = 0.0, sm = 0.0;
    size_t n = 0;
    for (size_t i = 0; i < tvec.size(); ++i) {
        if (tvec[i] < t_a || tvec[i] > t_b) continue;
        const double phase = 2.0 * M_PI * 50.0 * tvec[i];
        const double iref  = Iin_pk * std::fabs(std::sin(phase));
        const double e     = sliding_mean(tvec[i]) - iref;
        sse += e * e;
        sm  += iref * iref;
        ++n;
    }
    REQUIRE(n > 0);
    const double nrmse_env = std::sqrt(sse / n) / std::sqrt(sm / n);
    INFO("envelope NRMSE=" << 100.0*nrmse_env << " % (tol "
         << 100.0*s.tol_envelope << " %)");

    // Optional waveform dump (env var MKF_DUMP_WAVEFORMS=1) — emits
    // /tmp/mkf_wf_pfc_<name>.csv with time, ideal sinusoidal reference,
    // sliding-mean SPICE iL, raw SPICE iL.  Used for debugging analytical
    // models against the SPICE truth.
    if (const char* dump = std::getenv("MKF_DUMP_WAVEFORMS"); dump && dump[0] == '1') {
        std::string fname = "/tmp/mkf_wf_pfc_";
        fname += s.name;
        fname += ".csv";
        std::ofstream f(fname);
        f << "time_s,vin_V,iL_ref_A,iL_avg_A,iL_raw_A\n";
        for (size_t i = 0; i < tvec.size(); ++i) {
            if (tvec[i] < t_a || tvec[i] > t_b) continue;
            const double phase = 2.0 * M_PI * 50.0 * tvec[i];
            const double iref  = Iin_pk * std::fabs(std::sin(phase));
            f << tvec[i] << ',' << vin[i] << ',' << iref << ','
              << sliding_mean(tvec[i]) << ',' << iL[i] << '\n';
        }
        std::cout << "[PFC PtP " << s.name << "] dumped to " << fname << "\n";
    }

    REQUIRE(nrmse_env < s.tol_envelope);

    // ---- Gate 4: Phase-6 ConverterPortChecks switching-leg helper -------
    // Reuse the SAME simulation already executed above — build a minimal
    // OperatingPoint from sim.waveforms and feed it to the canonical
    // check_pfc_switching_ports helper.  This avoids a SECOND full ngspice
    // transient run (each one takes 1.5x to 2x the wall-time of a Buck/
    // Boost test, dominating the PtP suite when duplicated).  Web-frontend
    // users can still call simulate_with_ngspice_switching directly with
    // any cycle count they want — the cost is theirs to bear there.
    //
    // check_pfc_switching_ports expects the OperatingPoint to be trimmed
    // to a whole number of full line cycles (it computes the line-cycle
    // mean / RMS of vin_rect to compare against analytical 0.9·Vrms / Vrms).
    // simulate_with_ngspice_switching trims to the last full cycle by
    // default; we replicate that here by slicing [t_a, t_b].
    auto sliceToLastCycle = [&](const std::vector<double>& src)
        -> std::vector<double> {
        std::vector<double> out;
        out.reserve(src.size());
        for (size_t i = 0; i < tvec.size(); ++i) {
            if (tvec[i] >= t_a && tvec[i] <= t_b) out.push_back(src[i]);
        }
        return out;
    };
    const std::vector<double> tvecLC = sliceToLastCycle(tvec);
    const std::vector<double> vinLC  = sliceToLastCycle(vin);
    const std::vector<double> vbusLC = sliceToLastCycle(vbus);
    const std::vector<double> iLLC   = sliceToLastCycle(iL);

    auto makeWfWithTime = [&](const std::vector<double>& data,
                               const std::vector<double>& time) -> Waveform {
        Waveform wf;
        wf.set_data(data);
        wf.set_time(time);
        return wf;
    };
    auto makeSig = [&](const Waveform& wf) -> SignalDescriptor {
        SignalDescriptor sig;
        sig.set_waveform(wf);
        return sig;
    };
    OperatingPoint op;
    OperatingPointExcitation ipx;   // InputPort:  V=vin_rect, I=iL
    ipx.set_name("InputPort");
    ipx.set_voltage(makeSig(makeWfWithTime(vinLC, tvecLC)));
    ipx.set_current(makeSig(makeWfWithTime(iLLC,  tvecLC)));
    OperatingPointExcitation psx;   // PowerStage: V=vbus,    I=iL
    psx.set_name("PowerStage");
    psx.set_voltage(makeSig(makeWfWithTime(vbusLC, tvecLC)));
    psx.set_current(makeSig(makeWfWithTime(iLLC,   tvecLC)));
    op.set_excitations_per_winding({ipx, psx});

    ConverterPortChecks::check_pfc_switching_ports(
        op,
        s.name,
        s.vrms,
        s.vbus,
        s.pout,
        /*vinTol*/      0.05,
        /*voutMeanTol*/ s.tol_vbus_pct / 100.0,
        /*iinMeanTol*/  s.tol_pin_pct  / 100.0);
    const double tPostSim = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - tPostSim0).count();
    std::cout << "[PFC PtP " << s.name << "] tPostSim=" << tPostSim << "s\n";
}

} // namespace


TEST_CASE("PFC reference design PtP — NCP1654 100 W",
          "[converter-model][pfc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "NCP1654-100W",
        /*vrms*/         230.0,
        /*vbus*/         400.0,
        /*pout*/         100.0,
        /*fsw*/          100e3,
        /*cbus*/         100e-6,
        /*L*/            3.30e-3,
        /*cycles*/       3,        // NCP1654 needs 3 cycles for Pin
                                   // balance to settle within ±5 %; the
                                   // 100 W design has the smallest Cbus
                                   // ripple-vs-Pout ratio of the trio.
        /*tol_vbus_pct*/ 6.0,
        /*tol_pin_pct*/  5.0,
        /*tol_envelope*/ 0.10,
        /*tol_walltime*/ 30.0
    };
    run_ptp_gates(s);
}

TEST_CASE("PFC reference design PtP — UCC28180 360 W",
          "[converter-model][pfc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "UCC28180-360W",
        /*vrms*/         230.0,
        /*vbus*/         390.0,
        /*pout*/         360.0,
        /*fsw*/          65e3,
        /*cbus*/         470e-6,
        /*L*/            1.25e-3,
        /*cycles*/       2,
        /*tol_vbus_pct*/ 6.0,
        /*tol_pin_pct*/  5.0,
        /*tol_envelope*/ 0.30,
        /*tol_walltime*/ 10.0
    };
    run_ptp_gates(s);
}

TEST_CASE("PFC reference design PtP — L4981 1000 W",
          "[converter-model][pfc-topology][refdesign][ptp][slow]") {
    RefDesignSpec s{
        /*name*/         "L4981-1000W",
        /*vrms*/         230.0,
        /*vbus*/         400.0,
        /*pout*/         1000.0,
        /*fsw*/          50e3,
        /*cbus*/         1500e-6,
        /*L*/            659e-6,
        /*cycles*/       2,
        /*tol_vbus_pct*/ 6.0,
        /*tol_pin_pct*/  5.0,
        /*tol_envelope*/ 0.60,
        /*tol_walltime*/ 30.0
    };
    run_ptp_gates(s);
}
