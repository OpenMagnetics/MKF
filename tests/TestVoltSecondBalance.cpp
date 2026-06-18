// ─────────────────────────────────────────────────────────────────────────────
// TestVoltSecondBalance.cpp
//
// Cross-topology gate enforcing the §5.0 rule of
// CONVERTER_MODELS_GOLDEN_GUIDE.md:
//
//   "Every signal that ends up in excitations_per_winding[i].voltage MUST be
//    measured at the magnetic component's winding ports, never at a rectified
//    output cap, DC link, or any node downstream of the rectifier."
//
// Volt-second balance on a true winding voltage in CCM steady state requires
//
//   | avg( v_winding(t) ) |  /  peak( | v_winding(t) | )   <  ε
//
// (we use ε = 2 % for the analytical path, 5 % for the SPICE path to absorb
// startup transients and discretization). A signal that violates this bound is
// almost certainly the rectified output node smuggled into the winding stream,
// or a v(node)−GND probe on a winding terminal that is not GND.
//
// Each topology gets its own TEST_CASE so the failure message names the
// topology + winding index + path that violates the rule. The test
// covers every winding (primary + every secondary), in both the analytical
// and SPICE paths.
//
// Output rails (Vout / ILout) are NOT checked here — they are intentionally
// DC by design and live in the separate converter-port stream returned by
// simulate_and_extract_topology_waveforms (see Test_<Topo>_ConverterPortWaveforms
// — to be added per topology).
//
// ─────────────────────────────────────────────────────────────────────────────
// BASELINE AUDIT (recorded at the commit that introduced this file).
// These failures define the §3.H+ work order — they are NOT to be relaxed.
//
//   PushPull       analytical ✅   SPICE ✅   (clean — fixed in 4ffd3c28)
//   DAB            analytical ✅   SPICE ❌   secondary probe wrong (48 %)
//   Buck           analytical ✅   SPICE ❌   primary probe wrong (26 %)
//   Boost          analytical ✅   SPICE ❌   primary probe wrong (100 % DC)
//   Flyback        analytical ✅   SPICE ❌   both probes wrong (15 %)
//   IsolatedBuck   analytical ❌   SPICE ❌   secondary probe wrong (89/98 %)
//                                              — bug is in BOTH paths, not
//                                                only SPICE
//
// Each fix follows the §5.0 template proven on PushPull (4ffd3c28) and the
// IsolatedBuck primary (22bfee9e):
//   1. Add Bv<winding>_diff <name>_diff 0 V=V(<dot>)−V(<other>)
//   2. Add v(<name>_diff) to .save
//   3. Update waveformMapping for that winding
//   4. Verify avg(V·I) > 0 (passive sign) — flip the B source if not
// And, when the analytical path also leaks output-rail V into the winding
// (IsolatedBuck), fix process_operating_point_for_input_voltage to
// synthesise the differential winding voltage symbolically rather than
// reuse the rectified output cap voltage.
// ─────────────────────────────────────────────────────────────────────────────
// ─────────────────────────────────────────────────────────────────────────────

#include "converter_models/Buck.h"
#include "converter_models/Boost.h"
#include "converter_models/Flyback.h"
#include "converter_models/Dab.h"
#include "converter_models/PushPull.h"
#include "converter_models/IsolatedBuck.h"
#include "processors/NgspiceRunner.h"

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

constexpr double kAnalyticalEps = 0.02;   // 2 % of peak — strict
constexpr double kSpiceEps      = 0.05;   // 5 % of peak — allow startup / discretization

// Compute |avg(V)| / peak(|V|). Returns 0 if the waveform is uniformly zero.
double normalised_avg(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double sum = 0.0, peak = 0.0;
    for (double x : v) {
        sum  += x;
        peak  = std::max(peak, std::fabs(x));
    }
    if (peak < 1e-12) return 0.0;
    return std::fabs(sum / static_cast<double>(v.size())) / peak;
}

// Per-winding check: for every excitation in the OP, fetch the voltage
// waveform and assert volt-second balance.
void check_all_windings(const std::vector<OperatingPoint>& ops,
                        const std::string&                 topoName,
                        const std::string&                 path,
                        double                             eps,
                        size_t                             expectedWindings) {
    REQUIRE_FALSE(ops.empty());
    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        const auto& exc = ops[opIdx].get_excitations_per_winding();
        REQUIRE(exc.size() == expectedWindings);
        for (size_t w = 0; w < exc.size(); ++w) {
            REQUIRE(exc[w].get_voltage().has_value());
            const auto wf = exc[w].get_voltage()->get_waveform().value();
            const auto& d  = wf.get_data();
            REQUIRE(!d.empty());
            double normAvg = normalised_avg(d);
            INFO(topoName << " [" << path << "] OP " << opIdx
                 << " winding " << w
                 << " — |avg(V)|/peak(|V|) = " << normAvg
                 << " (bound " << eps << ")");
            CHECK(normAvg < eps);
        }
    }
}

}  // namespace

// ────────────────────────────── DAB ─────────────────────────────────────────
TEST_CASE("Test_VoltSecondBalance_Dab", "[volt-second-balance][dab-topology]") {
    OpenMagnetics::Dab dab;
    DimensionWithTolerance iv; iv.set_nominal(400.0);
    dab.set_input_voltage(iv);
    dab.set_efficiency(0.95);

    MAS::DabOperatingPoint op;
    op.set_output_voltages({48.0});
    op.set_output_currents({10.0});
    op.set_switching_frequency(100e3);
    op.set_ambient_temperature(25.0);
    dab.set_operating_points({op});

    auto designReqs = dab.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
    double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

    SECTION("analytical path") {
        auto ops = dab.process_operating_points(turnsRatios, Lm);
        check_all_windings(ops, "DAB", "analytical", kAnalyticalEps, /*expectedWindings*/ 2);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        dab.set_num_periods_to_extract(1);
        auto ops = dab.simulate_and_extract_operating_points(turnsRatios, Lm);
        check_all_windings(ops, "DAB", "SPICE", kSpiceEps, /*expectedWindings*/ 2);
    }
}

// ────────────────────────────── Buck ────────────────────────────────────────
TEST_CASE("Test_VoltSecondBalance_Buck", "[volt-second-balance][buck-topology]") {
    OpenMagnetics::Buck buck;
    DimensionWithTolerance iv; iv.set_nominal(48.0);
    buck.set_input_voltage(iv);
    buck.set_efficiency(0.92);
    buck.set_current_ripple_ratio(0.3);

    BaseOperatingPoint op;
    op.set_output_voltages({12.0});
    op.set_output_currents({2.0});
    op.set_switching_frequency(250e3);
    op.set_ambient_temperature(25.0);
    buck.set_operating_points({op});

    double Lm = buck.process_design_requirements().get_magnetizing_inductance().get_minimum().value();

    SECTION("analytical path") {
        auto ops = buck.process_operating_points(std::vector<double>{}, Lm);
        check_all_windings(ops, "Buck", "analytical", kAnalyticalEps, /*expectedWindings*/ 1);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        auto ops = buck.simulate_and_extract_operating_points(Lm);
        check_all_windings(ops, "Buck", "SPICE", kSpiceEps, /*expectedWindings*/ 1);
    }
}

// ────────────────────────────── Boost ───────────────────────────────────────
TEST_CASE("Test_VoltSecondBalance_Boost", "[volt-second-balance][boost-topology]") {
    OpenMagnetics::Boost boost;
    DimensionWithTolerance iv; iv.set_nominal(12.0);
    boost.set_input_voltage(iv);
    boost.set_efficiency(0.92);
    boost.set_current_ripple_ratio(0.3);

    BaseOperatingPoint op;
    op.set_output_voltages({24.0});
    op.set_output_currents({1.0});
    op.set_switching_frequency(250e3);
    op.set_ambient_temperature(25.0);
    boost.set_operating_points({op});

    double Lm = boost.process_design_requirements().get_magnetizing_inductance().get_minimum().value();

    SECTION("analytical path") {
        auto ops = boost.process_operating_points(std::vector<double>{}, Lm);
        check_all_windings(ops, "Boost", "analytical", kAnalyticalEps, /*expectedWindings*/ 1);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        auto ops = boost.simulate_and_extract_operating_points(Lm);
        check_all_windings(ops, "Boost", "SPICE", kSpiceEps, /*expectedWindings*/ 1);
    }
}

// ────────────────────────────── Flyback ─────────────────────────────────────
TEST_CASE("Test_VoltSecondBalance_Flyback", "[volt-second-balance][flyback-topology]") {
    OpenMagnetics::Flyback fb;
    DimensionWithTolerance iv; iv.set_nominal(48.0);
    fb.set_input_voltage(iv);
    fb.set_diode_voltage_drop(0.5);
    fb.set_efficiency(0.9);
    fb.set_current_ripple_ratio(0.3);
    fb.set_maximum_duty_cycle(0.5);

    OpenMagnetics::FlybackOperatingPoint op;
    op.set_output_voltages({12.0});
    op.set_output_currents({1.0});
    op.set_switching_frequency(100e3);
    op.set_ambient_temperature(25.0);
    fb.set_operating_points({op});

    auto designReqs = fb.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
    double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

    SECTION("analytical path") {
        auto ops = fb.process_operating_points(turnsRatios, Lm);
        check_all_windings(ops, "Flyback", "analytical", kAnalyticalEps, /*expectedWindings*/ 2);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        fb.set_num_steady_state_periods(200);
        fb.set_num_periods_to_extract(1);
        auto ops = fb.simulate_and_extract_operating_points(turnsRatios, Lm);
        check_all_windings(ops, "Flyback", "SPICE", kSpiceEps, /*expectedWindings*/ 2);
    }
}

// ────────────────────────────── PushPull ────────────────────────────────────
TEST_CASE("Test_VoltSecondBalance_PushPull", "[volt-second-balance][push-pull-topology]") {
    OpenMagnetics::PushPull pp;
    // Vin must be high enough that for the chosen N (set by process_design_requirements
    // from Vin/Vout) we still have T1 ≤ T/2 in process_operating_points. Vin=24, Vout=5
    // → effective D ≈ Vo/(η·Vin) ≈ 0.23, comfortably under 0.5.
    DimensionWithTolerance iv; iv.set_nominal(24.0);
    pp.set_input_voltage(iv);
    pp.set_diode_voltage_drop(0.5);
    pp.set_efficiency(0.9);
    pp.set_current_ripple_ratio(0.4);

    OpenMagnetics::PushPullOperatingPoint op;
    op.set_output_voltages({5.0});
    op.set_output_currents({1.0});
    op.set_switching_frequency(420e3);
    op.set_ambient_temperature(25.0);
    pp.set_operating_points({op});

    auto designReqs = pp.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
    double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

    // PushPull has 2 primary halves (center-tapped) + 1 secondary = 3 windings.
    // Confirm by inspecting the analytical OP's excitation count.
    auto analytical = pp.process_operating_points(turnsRatios, Lm);
    REQUIRE(!analytical.empty());
    size_t expectedWindings = analytical[0].get_excitations_per_winding().size();

    SECTION("analytical path") {
        check_all_windings(analytical, "PushPull", "analytical", kAnalyticalEps, expectedWindings);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        pp.set_num_periods_to_extract(1);
        auto ops = pp.simulate_and_extract_operating_points(turnsRatios, Lm);
        check_all_windings(ops, "PushPull", "SPICE", kSpiceEps, expectedWindings);
    }
}

// ────────────────────────────── IsolatedBuck ────────────────────────────────
//
// NOTE — KNOWN §5.0 VIOLATION (tracked in REVIEW_PLAN §3.H Phase 4 remainders):
//   The secondary voltage probe in IsolatedBuck::generate_ngspice_circuit
//   currently saves v(sec<N>_rect) — the rectified output cap node — instead
//   of a Bvsec<N>_diff differential probe across the secondary winding. The
//   primary differential probe was added in commit 22bfee9e but the secondary
//   was left as a cosmetic remainder. This test will fail on the SPICE-path
//   secondary winding(s) until that is fixed; failure is the intended signal
//   driving the fix, NOT a reason to relax the bound.
TEST_CASE("Test_VoltSecondBalance_IsolatedBuck", "[volt-second-balance][isolated-buck-topology]") {
    OpenMagnetics::IsolatedBuck ib;
    DimensionWithTolerance iv; iv.set_nominal(48.0);
    ib.set_input_voltage(iv);
    ib.set_diode_voltage_drop(0.5);
    ib.set_efficiency(0.9);
    ib.set_current_ripple_ratio(0.3);

    IsolatedBuckOperatingPoint op;
    op.set_output_voltages({5.0, 5.0});
    op.set_output_currents({0.6, 5.0});
    op.set_switching_frequency(200e3);
    op.set_ambient_temperature(25.0);
    ib.set_operating_points({op});

    auto designReqs = ib.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
    double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

    auto analytical = ib.process_operating_points(turnsRatios, Lm);
    REQUIRE(!analytical.empty());
    size_t expectedWindings = analytical[0].get_excitations_per_winding().size();

    SECTION("analytical path") {
        check_all_windings(analytical, "IsolatedBuck", "analytical", kAnalyticalEps, expectedWindings);
    }
    SECTION("SPICE path") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");
        ib.set_num_periods_to_extract(1);
        auto ops = ib.simulate_and_extract_operating_points(turnsRatios, Lm);
        check_all_windings(ops, "IsolatedBuck", "SPICE", kSpiceEps, expectedWindings);
    }
}
