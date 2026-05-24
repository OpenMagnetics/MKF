// Phase 4 exploratory driver — runs the switching netlist for each ref
// design through ngspice and prints the convergence + acceptance-gate
// metrics (wall-time, bus settle, line-cycle current envelope, power
// balance). NOT a regression test (no asserts beyond success). Tag:
// [phase4-explore][.] so it never runs in the default Catch invocation.

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <source_location>
#include <string>
#include <vector>

#include "converter_models/PowerFactorCorrection.h"
#include "processors/NgspiceRunner.h"
#include "support/Painter.h"

using namespace OpenMagnetics;

namespace {

struct Spec {
    const char* name;
    double vrms, vbus, pout, fsw, cbus, L;
    int    cycles;
};

OpenMagnetics::PowerFactorCorrection build(const Spec& s) {
    OpenMagnetics::PowerFactorCorrection pfc;
    DimensionWithTolerance iv;
    iv.set_nominal(s.vrms); iv.set_minimum(s.vrms); iv.set_maximum(s.vrms);
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

// Linear-interpolated lookup of vector v(t) at time t0 (t monotonic).
double interp(const std::vector<double>& t,
              const std::vector<double>& v,
              double t0) {
    if (t0 <= t.front()) return v.front();
    if (t0 >= t.back())  return v.back();
    auto it = std::lower_bound(t.begin(), t.end(), t0);
    size_t i = std::distance(t.begin(), it);
    if (i == 0) return v.front();
    double t1 = t[i-1], t2 = t[i];
    double v1 = v[i-1], v2 = v[i];
    double a = (t0 - t1) / (t2 - t1);
    return v1 + a * (v2 - v1);
}

// Mean of v(t) over [t0, t1] using trapezoidal rule on the simulation's
// own non-uniform time grid (subset to the window).
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
    // a and b share the same time axis (both come from same sim).
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
    double pend = aend * bend;
    double pprev = aprev * bprev;
    area += 0.5 * (pend + pprev) * (t1 - tprev);
    return area / (t1 - t0);
}

void run_one(const Spec& s) {
    std::cout << "\n========== " << s.name << " ==========\n";
    auto pfc = build(s);
    const std::string netlist = pfc.generate_ngspice_switching_circuit(s.L, s.cycles);
    std::ofstream("/tmp/pfc_" + std::string(s.name) + ".cir") << netlist;

    NgspiceRunner runner;
    REQUIRE(runner.is_available());
    SimulationConfig cfg;
    cfg.frequency        = 50.0;
    cfg.extractOnePeriod = false;
    cfg.numberOfPeriods  = static_cast<size_t>(s.cycles);
    cfg.keepTempFiles    = true;
    cfg.timeout          = 180.0;

    auto t0 = std::chrono::steady_clock::now();
    auto sim = runner.run_simulation(netlist, cfg);
    auto dt = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::cout << "wall_time=" << dt << " s   success=" << sim.success
              << "   waveforms=" << sim.waveformNames.size() << "\n";
    if (!sim.success) {
        std::cout << "  ERROR: " << sim.errorMessage << "\n";
        return;
    }

    // Resolve waveforms by name (case-insensitive prefix/suffix match).
    auto find_by = [&](const std::string& key) -> const Waveform* {
        for (size_t k = 0; k < sim.waveformNames.size(); ++k) {
            std::string n = sim.waveformNames[k];
            std::transform(n.begin(), n.end(), n.begin(),
                           [](unsigned char c){return std::tolower(c);});
            if (n.find(key) != std::string::npos) return &sim.waveforms[k];
        }
        return nullptr;
    };
    auto* wfTime  = find_by("time");
    auto* wfVbus  = find_by("vbus");
    auto* wfVin   = find_by("vin_rect");
    auto* wfIL    = find_by("vl_sense");
    if (!wfTime || !wfVbus || !wfVin || !wfIL) {
        std::cout << "  Missing required waveform.\n";
        return;
    }
    auto tvec_opt = wfTime->get_time();
    REQUIRE(tvec_opt.has_value());
    const auto tvec = tvec_opt.value();
    const auto vbus = wfVbus->get_data();
    const auto vin  = wfVin->get_data();
    const auto iL   = wfIL->get_data();
    std::cout << "  time samples=" << tvec.size()
              << "  t_end=" << tvec.back()*1e3 << " ms\n";

    const double Tline = 1.0 / 50.0;
    // Last full line cycle for steady-state metrics.
    const double t_a = tvec.back() - Tline;
    const double t_b = tvec.back();

    double vbus_mean   = mean_over(tvec, vbus, t_a, t_b);
    double vbus_min    = *std::min_element(vbus.begin(), vbus.end());
    double vbus_max    = *std::max_element(vbus.begin(), vbus.end());
    double iL_mean     = mean_over(tvec, iL,   t_a, t_b);
    double pin         = mean_product(tvec, vin, iL, t_a, t_b);

    std::cout << "  vbus mean(last cycle) = " << vbus_mean
              << "  V   target=" << s.vbus
              << "   err=" << 100.0*(vbus_mean - s.vbus)/s.vbus << " %\n";
    std::cout << "  vbus min/max overall   = " << vbus_min << " / " << vbus_max << "\n";
    std::cout << "  iL  mean(last cycle)   = " << iL_mean << " A\n";
    std::cout << "  Pin = mean(vin*iL)     = " << pin
              << " W   target=" << s.pout
              << "   err=" << 100.0*(pin - s.pout)/s.pout << " %\n";

    // Current-envelope NRMSE vs ideal sin envelope at line freq.
    //
    // The raw inductor current carries 100 kHz switching ripple riding on the
    // 50 Hz line-frequency envelope.  Comparing instantaneous iL to the ideal
    // half-sine envelope gives an NRMSE floor set by the ripple itself
    // (ΔI = Vin·D/(L·fsw) at the line peak), which for L4981's 680 µH choke
    // is ≈12 % — leaving no room for tracking error.  Compute *both*: the raw
    // value (for visibility) and an envelope-vs-envelope value where iL is
    // first averaged over a sliding window of one switching period (Tsw =
    // 1/fsw), which collapses the ripple and exposes line-frequency tracking
    // error only.
    const double Iin_rms = s.pout / s.vrms;
    const double Iin_pk  = std::sqrt(2.0) * Iin_rms;
    const double Tsw     = 1.0 / s.fsw;
    auto sliding_mean = [&](double t0)->double{
        // Centred moving average of iL over [t0-Tsw/2, t0+Tsw/2].
        double a = t0 - 0.5*Tsw;
        double b = t0 + 0.5*Tsw;
        if (a < tvec.front()) a = tvec.front();
        if (b > tvec.back())  b = tvec.back();
        return mean_over(tvec, iL, a, b);
    };
    double sse_raw=0.0, sse_env=0.0, sm=0.0; size_t n=0;
    for (size_t i=0;i<tvec.size();++i) {
        if (tvec[i] < t_a || tvec[i] > t_b) continue;
        double phase = 2.0*M_PI*50.0*tvec[i];
        double iref  = Iin_pk * std::fabs(std::sin(phase));
        double e_raw = iL[i] - iref;
        double e_env = sliding_mean(tvec[i]) - iref;
        sse_raw += e_raw*e_raw;
        sse_env += e_env*e_env;
        sm      += iref*iref;
        ++n;
    }
    double nrmse_raw = std::sqrt(sse_raw/std::max<size_t>(n,1)) /
                       std::sqrt(sm/std::max<size_t>(n,1));
    double nrmse_env = std::sqrt(sse_env/std::max<size_t>(n,1)) /
                       std::sqrt(sm/std::max<size_t>(n,1));
    std::cout << "  iL NRMSE raw / envelope vs ideal sin = "
              << 100.0*nrmse_raw << " % / "
              << 100.0*nrmse_env << " %\n";

    // Stack-plot the 3 diagnostic windings (PowerStage / VoltageLoop /
    // CurrentLoop) packed by simulate_with_ngspice_switching.  The default
    // trim-to-last-line-cycle keeps SVGs to a single 20 ms window of the
    // steady-state output.
    //
    // SVG emission is opt-in via the `MKF_PFC_EMIT_SVG` environment
    // variable.  The stacked plots run 8–16 MB each (NCP1654 ≈ 16 MB,
    // UCC28180 ≈ 12 MB, L4981 ≈ 8 MB) and re-running this exploratory
    // harness should not silently bloat the workspace.  Set
    //   `MKF_PFC_EMIT_SVG=1` (any non-empty value) to write the files;
    // leave unset to skip the second simulate-and-paint pass entirely.
    const char* emitSvg = std::getenv("MKF_PFC_EMIT_SVG");
    if (emitSvg && emitSvg[0] != '\0') {
        auto outDir = std::filesystem::path{
            std::source_location::current().file_name()
        }.parent_path() / ".." / "output" / "pfc_phase4";
        std::filesystem::create_directories(outDir);
        auto outFile = outDir / (std::string(s.name) + "_waveforms.svg");

        auto op = pfc.simulate_with_ngspice_switching(s.L, s.cycles);
        Painter painter(outFile);
        painter.paint_operating_point_waveforms(
            op,
            std::string(s.name) + " — last line cycle (steady state)",
            1400.0, 1000.0);
        painter.export_svg();
        std::cout << "  SVG: " << outFile.string() << "\n";
    } else {
        std::cout << "  SVG: skipped (set MKF_PFC_EMIT_SVG=1 to emit)\n";
    }
}

} // namespace

TEST_CASE("Phase4_Explore_All", "[phase4-explore][.]") {
    Spec s1{"NCP1654-100W",  230.0, 400.0,  100.0, 100e3,  100e-6, 3.30e-3, 3};
    Spec s2{"UCC28180-360W", 230.0, 390.0,  360.0,  65e3,  470e-6, 1.25e-3, 3};
    Spec s3{"L4981-1000W",   230.0, 400.0, 1000.0,  50e3, 1500e-6, 659e-6,  3};
    run_one(s1);
    run_one(s2);
    run_one(s3);
}
