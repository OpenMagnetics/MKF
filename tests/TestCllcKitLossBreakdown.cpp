// Diagnostic — KIT-20kW-Asymmetric CLLC: identify where SPICE's 33 % "loss"
// (Pin=29.4 kW vs Pout=19.7 kW) is being dissipated. Compares against the
// analytical lossless FHA peak primary current (~39 A) which differs from
// SPICE's ~54 A precisely because SPICE delivers ~10 kW more from the
// source — almost certainly into switch RON, snubbers, or SR cross-
// conduction during the 10 ns gate transitions.
//
// Hidden behind [.diag][cllc-loss-breakdown] so it does NOT run in the
// normal CI matrix; invoke explicitly:
//   ./MKF_tests "[cllc-loss-breakdown]"

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include "converter_models/Cllc.h"
#include "processors/NgspiceRunner.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

// Augment the netlist's .save line with branch currents for every switch
// and snubber resistor, then add an extra full-period transient window so
// we can integrate steady-state powers cleanly.
std::string augment_saves(const std::string& netlist) {
    // Branches to add (ngspice voltage-controlled-switch current = i(sN);
    // resistor current = i(rN)).
    const char* extras =
        "+ i(s1) i(s2) i(s3) i(s4) i(sa) i(sb) i(sc) i(sd)\n"
        "+ i(rsn_s1) i(rsn_s2) i(rsn_s3) i(rsn_s4)\n"
        "+ i(rsn_sa) i(rsn_sb) i(rsn_sc) i(rsn_sd)\n"
        "+ v(node_a) v(node_b) v(node_c) v(node_d)\n"
        "+ v(ns1) v(ns2) v(ns3) v(ns4) v(nsa) v(nsb) v(nsc) v(nsd)\n";
    // Find the existing .save block (starts with ".save" and continues with
    // '+' continuation lines). Inject extras immediately before the blank
    // line that terminates the block.
    std::string out = netlist;
    auto pos = out.find(".save");
    if (pos == std::string::npos) return out;
    // Walk forward over continuation lines to find the end of the .save block.
    auto end = pos;
    while (end < out.size()) {
        auto eol = out.find('\n', end);
        if (eol == std::string::npos) break;
        auto next_start = eol + 1;
        if (next_start >= out.size()) { end = next_start; break; }
        char first = out[next_start];
        if (first != '+') { end = next_start; break; }
        end = next_start;
    }
    out.insert(end, extras);
    return out;
}

// Trapezoidal time-weighted mean of f(t)*g(t) over [t_lo, t_hi].
double mean_product_window(const std::vector<double>& t,
                           const std::vector<double>& f,
                           const std::vector<double>& g,
                           double t_lo, double t_hi) {
    if (t.size() < 2 || t.size() != f.size() || t.size() != g.size()) return 0.0;
    double sum = 0.0, span = 0.0;
    for (size_t k = 1; k < t.size(); ++k) {
        double a = std::max(t[k-1], t_lo), b = std::min(t[k], t_hi);
        if (b <= a) continue;
        // Linear interp f,g at a and b.
        double dt = t[k] - t[k-1];
        if (dt <= 0) continue;
        auto interp = [&](double tt, const std::vector<double>& v) {
            return v[k-1] + (tt - t[k-1]) / dt * (v[k] - v[k-1]);
        };
        double fa = interp(a, f), fb = interp(b, f);
        double ga = interp(a, g), gb = interp(b, g);
        double pa = fa * ga, pb = fb * gb;
        sum += 0.5 * (pa + pb) * (b - a);
        span += (b - a);
    }
    return (span > 0) ? sum / span : 0.0;
}

double mean_square_window(const std::vector<double>& t,
                          const std::vector<double>& f,
                          double t_lo, double t_hi) {
    return mean_product_window(t, f, f, t_lo, t_hi);
}

} // namespace

TEST_CASE("CLLC diag — KIT 20 kW asymmetric loss breakdown",
          "[cllc-loss-breakdown][.diag]") {
    nlohmann::json fixture{
        {"inputVoltage", {{"nominal", 800.0}, {"minimum", 700.0}, {"maximum", 900.0}}},
        {"maxSwitchingFrequency", 300e3},
        {"minSwitchingFrequency", 50e3},
        {"efficiency", 0.95},
        {"qualityFactor", 0.3},
        {"symmetricDesign", false},
        {"resonantInductorRatio", 0.95},
        {"resonantCapacitorRatio", 1.052},
        {"bidirectional", true},
        {"operatingPoints", nlohmann::json::array({
            nlohmann::json{
                {"outputVoltages", {800.0}},
                {"outputCurrents", {25.0}},
                {"switchingFrequency", 100e3},
                {"ambientTemperature", 25.0},
                {"powerFlow", "forward"}
            }
        })}
    };

    CllcConverter cllc(fixture);
    auto params = cllc.calculate_resonant_parameters();
    const double n = params.turnsRatio;

    std::string netlist = cllc.generate_ngspice_circuit(n, params, 0, 0);
    netlist = augment_saves(netlist);

    NgspiceRunner runner;
    if (!runner.is_available()) { WARN("ngspice not available"); return; }

    SimulationConfig cfg;
    cfg.frequency = 100e3;
    cfg.extractOnePeriod = false;
    cfg.numberOfPeriods = 4;       // collect several cycles of steady state
    cfg.steadyStateCycles = 30;    // skip transient
    cfg.timeout = 120.0;
    cfg.keepTempFiles = false;

    auto sim = runner.run_simulation(netlist, cfg);
    REQUIRE(sim.success);

    // Resolve waveforms by name (case-insensitive substring match).
    auto find_by = [&](const std::string& key) -> const Waveform* {
        for (size_t k = 0; k < sim.waveformNames.size(); ++k) {
            std::string nm = sim.waveformNames[k];
            std::transform(nm.begin(), nm.end(), nm.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            if (nm.find(key) != std::string::npos) return &sim.waveforms[k];
        }
        return nullptr;
    };
    const Waveform* wfTime = find_by("time");
    REQUIRE(wfTime != nullptr);
    const auto tvec_opt = wfTime->get_time();
    REQUIRE(tvec_opt.has_value());
    const auto t = tvec_opt.value();
    REQUIRE(t.size() > 100);

    const double Tperiod = 1.0 / 100e3;
    const double t_hi = t.back();
    const double t_lo = std::max(t.front(), t_hi - 2.0 * Tperiod);
    std::cout << "[KIT-loss] integration window: " << t_lo*1e6 << " us to "
              << t_hi*1e6 << " us  (" << (t_hi-t_lo)/Tperiod << " cycles)\n";

    // Helper: dual lookup since ngspice saves nodes as bare names (no v())
    // and source branch currents as "<name>#branch".
    auto wf_node = [&](const std::string& node) -> const Waveform* {
        const Waveform* w = find_by(node);
        return w;
    };
    auto wf_branch = [&](const std::string& src) -> const Waveform* {
        return find_by(src + "#branch");
    };

    // Pin = <V(vin_p) · -i(vin)>  (i(vin) flows from + to - INSIDE source,
    // so external current INTO vin_p is -i(vin)).
    {
        const Waveform* wfI = wf_branch("vin");
        const Waveform* wfV = wf_node("vin_p");
        REQUIRE(wfI != nullptr);
        REQUIRE(wfV != nullptr);
        auto iD = wfI->get_data();
        auto vD = wfV->get_data();
        std::vector<double> iext(iD.size());
        for (size_t k = 0; k < iD.size(); ++k) iext[k] = -iD[k];
        double pin = mean_product_window(t, vD, iext, t_lo, t_hi);
        std::cout << "[KIT-loss] Pin = " << pin << " W (lossless analytical = 20 kW)\n";
    }

    // Pout = <V(vout_load) · I(Vout_sense)>  — but easier: use V(vout_p)²/Rload.
    {
        const Waveform* wfVout = wf_node("vout_p");
        const Waveform* wfIout = wf_branch("vout_sense");
        if (wfVout && wfIout) {
            auto vD = wfVout->get_data();
            auto iD = wfIout->get_data();
            // i(Vout_sense) flows from vout_p TO vout_load (the +V side of the
            // load sense source). Positive when current leaves vout_p into the
            // load. Power INTO load = V(vout_p)·i(Vout_sense).
            double pout = mean_product_window(t, vD, iD, t_lo, t_hi);
            std::cout << "[KIT-loss] Pout (load) = " << pout << " W\n";
        }
    }

    // Snubber power: each R-snubber R=1kΩ. Compute V_R = V(node_in) - V(ns_n),
    // then P = mean(V_R²) / R.
    struct Snub { const char* name; const char* a; const char* b; };
    std::vector<Snub> snubs = {
        {"Rsn_S1", "vin_p", "ns1"},
        {"Rsn_S2", "node_a", "ns2"},
        {"Rsn_S3", "vin_p", "ns3"},
        {"Rsn_S4", "node_b", "ns4"},
        {"Rsn_Sa", "node_c", "nsa"},
        {"Rsn_Sb", "vout_n", "nsb"},
        {"Rsn_Sc", "node_d", "nsc"},
        {"Rsn_Sd", "vout_n", "nsd"},
    };
    double p_snub_total = 0.0;
    for (const auto& s : snubs) {
        const Waveform* wfA = wf_node(s.a);
        const Waveform* wfB = wf_node(s.b);
        if (!wfA || !wfB) {
            std::cout << "[KIT-loss]   " << s.name << " : missing nodes\n";
            continue;
        }
        auto vA = wfA->get_data();
        auto vB = wfB->get_data();
        std::vector<double> vR(vA.size());
        for (size_t k = 0; k < vR.size(); ++k) vR[k] = vA[k] - vB[k];
        double v2 = mean_square_window(t, vR, t_lo, t_hi);
        double p = v2 / 1000.0;
        std::cout << "[KIT-loss]   " << s.name << " : V_R_rms=" << std::sqrt(v2)
                  << " V, P=" << p << " W\n";
        p_snub_total += p;
    }
    std::cout << "[KIT-loss] TOTAL snubber loss: " << p_snub_total << " W\n";

    // Switch loss via KCL: for each leg, the bridge current entering the
    // midpoint goes through one of the two leg switches. We can compute the
    // RESISTIVE conduction loss in each switch from V_ds·i(switch), but we
    // don't have i(switch). Approximate aggregate: total switch loss =
    // Pin - Pout - p_snub - (tank loss, ~0 lossless caps/inductors).
    // (The K=0.9999 leakage adds negligible loss; caps and inductors are
    // lossless in the netlist.)
    std::cout << "[KIT-loss] note: switch loss = Pin − Pout − snubbers\n";
}
