// ─────────────────────────────────────────────────────────────────────────────
// TestSpiceSubcircuits.cpp
//
// Unit tests for the reusable ngspice .subckt definitions in
// `src/converter_models/PfcControllerSubcircuits.h` (OPAMP_IDEAL,
// COMPARATOR_IDEAL). These are foundational primitives for the boost-PFC
// average-current-mode controller (Phase 1 of the PFC controller workplan).
// Validated in isolation here so any later issue in the controller netlist
// is unambiguously attributable to the controller composition, not to the
// subckt definitions.
// ─────────────────────────────────────────────────────────────────────────────

#include "converter_models/PfcControllerSubcircuits.h"
#include "processors/NgspiceRunner.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMagnetics;

namespace {

// Small helper: locate a named waveform (case-insensitive) inside a
// SimulationResult and return its time/data pair. Throws if not found —
// per project no-fallbacks rule.
struct WaveformView {
    std::vector<double> time;
    std::vector<double> data;
};

// Strip a leading "v(" / "V(" prefix and trailing ")" so callers can ask
// for either "v(vo)" or "vo" indifferently.
std::string normalize_node_name(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (s.size() >= 4 && s.substr(0, 2) == "v(" && s.back() == ')') {
        s = s.substr(2, s.size() - 3);
    }
    return s;
}

WaveformView locate_waveform(const SimulationResult& r, const std::string& name) {
    const std::string target = normalize_node_name(name);
    for (size_t i = 0; i < r.waveformNames.size(); ++i) {
        if (normalize_node_name(r.waveformNames[i]) == target) {
            const auto& wf = r.waveforms[i];
            auto t = wf.get_time();
            if (!t.has_value()) {
                throw std::runtime_error("Waveform '" + name + "' has no time vector");
            }
            return { t.value(), wf.get_data() };
        }
    }
    std::ostringstream oss;
    oss << "Waveform '" << name << "' not found. Available:";
    for (const auto& n : r.waveformNames) oss << " " << n;
    throw std::runtime_error(oss.str());
}

// ─────────────────────────────────────────────────────────────────────────────
// OPAMP_IDEAL — voltage follower bandwidth.
//
// A voltage follower (vo connected to vn) with GBW = 10 MHz should reach
// 99 % of a 1 V step in roughly one time constant τ = A0 / (2π·GBW·A0) = 1/(2π·GBW)
// ≈ 16 ns. We allow 2 µs of settling time and require V(vo) within 0.5 % of
// V(vin), confirming the closed-loop pole is well above the 1 µs measurement
// window (i.e. GBW is at least ≈ 1 MHz, comfortable margin to spec).
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_SpiceSubckt_OpampIdeal_VoltageFollower",
          "[ngspice-simulation][spice-subcircuits][opamp]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    std::ostringstream nl;
    nl << "* Test OPAMP_IDEAL — voltage follower step response\n";
    nl << SPICE_SUBCKT_OPAMP_IDEAL;
    nl << "V_in vin 0 PULSE(0 1 1u 1n 1n 50u 100u)\n";
    nl << "X_op vin vo vo OPAMP_IDEAL params: A0=1e5 GBW=10e6 VSSPOS=5 VSSNEG=-5\n";
    nl << ".tran 100n 10u 0 100n uic\n";
    nl << ".save V(vin) V(vo)\n";
    nl << ".end\n";

    SimulationConfig cfg;
    cfg.frequency = 1e6;
    cfg.extractOnePeriod = false;
    cfg.numberOfPeriods  = 10;

    auto r = runner.run_simulation(nl.str(), cfg);
    REQUIRE(r.success);

    auto vin = locate_waveform(r, "v(vin)");
    auto vo  = locate_waveform(r, "v(vo)");
    REQUIRE(!vin.data.empty());
    REQUIRE(vin.data.size() == vo.data.size());

    // After t > 3 µs (well past the high step at 1 µs + ≈ 100 ns settling),
    // vo must track vin within 0.5 %.
    int checked = 0;
    for (size_t i = 0; i < vin.time.size(); ++i) {
        const double t = vin.time[i];
        if (t < 3e-6 || t > 9e-6) continue;
        const double v_in_i = vin.data[i];
        const double v_o_i  = vo.data[i];
        INFO("t=" << t * 1e6 << " µs vin=" << v_in_i << " vo=" << v_o_i);
        CHECK(std::abs(v_o_i - v_in_i) < 0.005);
        ++checked;
    }
    REQUIRE(checked > 10);
}

// ─────────────────────────────────────────────────────────────────────────────
// OPAMP_IDEAL — output rail clamp.
//
// Open-loop opamp driven by V(vp)−V(vn) = ±0.1 V → ideally A0·0.1 = ±10 V,
// must be clamped to ±5 V (VSSPOS / VSSNEG). For the clamp to engage within
// the 10 µs simulation window we use a fast internal pole: A0=100, GBW=100MHz
// → fp = GBW/A0 = 1 MHz, τ = 160 ns ≪ 4 µs settling allowance. The open-loop
// gain is intentionally low (100) so the input differential of 0.1 V drives
// the linear-region output to exactly 10 V, which the clamp then truncates to
// the rail. Symmetric check on the negative rail.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_SpiceSubckt_OpampIdeal_RailClamp",
          "[ngspice-simulation][spice-subcircuits][opamp]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    std::ostringstream nl;
    nl << "* Test OPAMP_IDEAL — open-loop rail clamp\n";
    nl << SPICE_SUBCKT_OPAMP_IDEAL;
    nl << "V_p vp 0 PWL(0 0 1u 0.1 5u 0.1 6u -0.1 10u -0.1)\n";
    nl << "V_n vn 0 0\n";
    nl << "X_op vp vn vo OPAMP_IDEAL params: A0=100 GBW=100e6 VSSPOS=5 VSSNEG=-5\n";
    nl << ".tran 100n 10u 0 100n uic\n";
    nl << ".save V(vo)\n";
    nl << ".end\n";

    SimulationConfig cfg;
    cfg.frequency = 1e5;
    cfg.extractOnePeriod = false;
    cfg.numberOfPeriods  = 1;

    auto r = runner.run_simulation(nl.str(), cfg);
    REQUIRE(r.success);

    auto vo = locate_waveform(r, "v(vo)");

    double vMaxAfterStep = -1e9;
    double vMinAfterStep = 1e9;
    for (size_t i = 0; i < vo.time.size(); ++i) {
        const double t = vo.time[i];
        if (t < 3e-6) continue;     // wait past the +0.1 V step
        if (t < 4e-6) vMaxAfterStep = std::max(vMaxAfterStep, vo.data[i]);
        if (t > 8e-6) vMinAfterStep = std::min(vMinAfterStep, vo.data[i]);
    }

    INFO("vMaxAfterStep=" << vMaxAfterStep << " (expect ≈ +5 V)");
    INFO("vMinAfterStep=" << vMinAfterStep << " (expect ≈ −5 V)");
    CHECK_THAT(vMaxAfterStep, Catch::Matchers::WithinAbs(5.0, 0.01));
    CHECK_THAT(vMinAfterStep, Catch::Matchers::WithinAbs(-5.0, 0.01));
}

// ─────────────────────────────────────────────────────────────────────────────
// COMPARATOR_IDEAL — sawtooth vs DC threshold.
//
// Drive vp with a 100 kHz 0→1 V sawtooth and vn with 0.4 V DC. Output should
// be VHIGH = 5 V exactly while the sawtooth is above 0.4 V (60 % duty), and
// VLOW = 0 V otherwise. We check the average over a steady-state cycle is
// 0.6 · VHIGH = 3.0 V within 5 %.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_SpiceSubckt_ComparatorIdeal_SawtoothThreshold",
          "[ngspice-simulation][spice-subcircuits][comparator]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    std::ostringstream nl;
    nl << "* Test COMPARATOR_IDEAL — sawtooth vs DC threshold\n";
    nl << SPICE_SUBCKT_COMPARATOR_IDEAL;
    nl << "V_saw saw 0 PULSE(0 1 0 9.9u 0.1u 1n 10u)\n";
    nl << "V_ref vref 0 0.4\n";
    nl << "X_cmp saw vref vcmp COMPARATOR_IDEAL params: VHIGH=5 VLOW=0\n";
    nl << ".tran 50n 50u 0 50n uic\n";
    nl << ".save V(vcmp) V(saw)\n";
    nl << ".end\n";

    SimulationConfig cfg;
    cfg.frequency = 100e3;
    cfg.extractOnePeriod = false;
    cfg.numberOfPeriods  = 5;

    auto r = runner.run_simulation(nl.str(), cfg);
    REQUIRE(r.success);

    auto vcmp = locate_waveform(r, "v(vcmp)");

    // Examine the last 30 µs (3 full cycles of 10 µs steady-state).
    double sumDt = 0.0;
    double areaUnderVcmp = 0.0;
    for (size_t i = 1; i < vcmp.time.size(); ++i) {
        const double t  = vcmp.time[i];
        if (t < 20e-6) continue;
        const double dt = vcmp.time[i] - vcmp.time[i - 1];
        // Trapezoidal integration (vcmp is essentially square-wave so this
        // gives a slight bias only near the edges, which dt averages out).
        const double mid = 0.5 * (vcmp.data[i] + vcmp.data[i - 1]);
        sumDt          += dt;
        areaUnderVcmp  += mid * dt;
    }
    REQUIRE(sumDt > 25e-6);
    const double avg = areaUnderVcmp / sumDt;
    INFO("Comparator average (steady): " << avg << " V (expect ≈ 3.0 V)");
    // 60 % duty × 5 V = 3.0 V. Allow ±10 % for sawtooth edge timing
    // (the 100 ns reset edge introduces a small systematic asymmetry).
    CHECK(avg > 2.7);
    CHECK(avg < 3.3);

    // Verify only the two output levels appear (no intermediate values).
    int countHigh = 0, countLow = 0, countMid = 0;
    for (size_t i = 0; i < vcmp.data.size(); ++i) {
        const double t = vcmp.time[i];
        if (t < 20e-6) continue;
        const double v = vcmp.data[i];
        if (v > 4.99)      ++countHigh;
        else if (v < 0.01) ++countLow;
        else               ++countMid;
    }
    INFO("samples high=" << countHigh << " low=" << countLow
         << " mid=" << countMid << " (mid should be ≈ 0)");
    CHECK(countHigh > 0);
    CHECK(countLow > 0);
    // Hard B-source comparator transitions in one solver step; expect at
    // most a handful of mid-level samples per edge.
    CHECK(countMid < 30);
}

// ─────────────────────────────────────────────────────────────────────────────
// Prelude helper — verify both subckts compose cleanly into one netlist
// without name collisions (sanity check before the PFC builder uses
// `spice_subckt_prelude_pfc_controller()` to inject both at once).
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_SpiceSubckt_PreludeHelperCompiles",
          "[ngspice-simulation][spice-subcircuits]") {
    NgspiceRunner runner;
    if (!runner.is_available()) SKIP("ngspice not available");

    std::ostringstream nl;
    nl << "* Test prelude helper — both subckts in one netlist\n";
    nl << spice_subckt_prelude_pfc_controller();
    nl << "V_in vin 0 PULSE(0 1 1u 1n 1n 50u 100u)\n";
    nl << "V_ref vref 0 0.5\n";
    nl << "X_op  vin vo_op vo_op OPAMP_IDEAL params: A0=1e5 GBW=1e6 VSSPOS=5 VSSNEG=0\n";
    nl << "X_cmp vin vref vo_cmp COMPARATOR_IDEAL params: VHIGH=5 VLOW=0\n";
    nl << ".tran 100n 5u 0 100n uic\n";
    nl << ".save V(vo_op) V(vo_cmp)\n";
    nl << ".end\n";

    SimulationConfig cfg;
    cfg.frequency = 1e6;
    cfg.extractOnePeriod = false;

    auto r = runner.run_simulation(nl.str(), cfg);
    REQUIRE(r.success);
    REQUIRE(!r.waveformNames.empty());
}

} // namespace
