// Mutual (cross-coupling) resistance ngspice export — abt #50/#72/#76.
//
// The off-diagonal Hesterman winding-loss term. The historical n==2 realization used a
// coupled auxiliary winding whose emitted R-L ladder topology did NOT reproduce its own
// fitted model (its DC impedance was the last-stage R, not dcMutualResistance — abt #76),
// and the n>=3 auxiliary realization left the coupled-L matrix non-positive-definite so
// ngspice collapsed. Both are now unified onto the PD-safe BEHAVIOURAL realization: each
// winding carries a series voltage drop sum_{j!=i} R_ij(f_export)*I_j built from linear
// current-controlled voltage sources (Hmut CCVS) + a current sense (Vsns), adding NO
// inductor and NO coupling, so the Lmag matrix stays positive-definite. These tests pin
// that realization and verify the emitted Hmut gains reproduce the fitted mutual
// resistance R_ij(f) at the DC and HF limits.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "processors/CircuitSimulatorInterface.h"
#include "support/Settings.h"
#include "TestingUtils.h"

// Kirchhoff's in-process libngspice runner, compiled into the test binary by
// tests/KirchhoffNgspiceRunnerShim.cpp (the libKirchhoffApi.so facade hides it).
#include "NgspiceRunner.hpp"

using namespace OpenMagnetics;

namespace {

// Emit an n-winding transformer subcircuit with the mutual-resistance term on, at a chosen
// export frequency (the frequency at which the behavioural realization freezes R_ij).
std::string export_transformer_subckt(const std::vector<int64_t>& numberTurns,
                                      double frequency = 100000.0) {
    auto gapping = OpenMagneticsTesting::get_residual_gap();
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 39", gapping, numberTurns, 1, "3C97");
    // get_quick_magnetic leaves the coil unwound; the subcircuit exporter needs a
    // processed coil (turns description) to compute per-winding + mutual resistances.
    auto coil = magnetic.get_coil();
    coil.wind();
    magnetic.set_coil(coil);
    auto& settings = Settings::GetInstance();
    bool previous = settings.get_circuit_simulator_include_mutual_resistance();
    settings.set_circuit_simulator_include_mutual_resistance(true);
    CircuitSimulatorExporterNgspiceModel exporter;
    std::string subckt = exporter.export_magnetic_as_subcircuit(magnetic, frequency, 25.0);
    settings.set_circuit_simulator_include_mutual_resistance(previous);
    return subckt;
}

// Parse the subcircuit name from the emitted ".subckt <name> ..." line.
std::string parse_subckt_name(const std::string& subckt) {
    std::istringstream in(subckt);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string tok0, name;
        ls >> tok0 >> name;
        if (tok0 == ".subckt") return name;
    }
    return "";
}

// Run the emitted subcircuit through Kirchhoff's in-process libngspice runner and measure
// the mutual trans-impedance V(P2+,P2-)/I1 with winding 2 OPEN and a DC current I0 driven
// through winding 1. At DC steady state every inductive term is inert (dI/dt = 0, so the
// K-coupled transformer voltage is zero and every ladder inductor is a short) and I2 = 0,
// so the ONLY contribution to V(P2+,P2-) is the Hmut_2_1 CCVS drop = R12(f_export) * I0.
// This exercises the FULL emitted deck end-to-end (parse + solve), and reads back exactly
// the mutual resistance the emission realized.
double measure_mutual_transimpedance_ngspice(const std::string& subckt, double I0) {
    std::string name = parse_subckt_name(subckt);
    REQUIRE_FALSE(name.empty());
    std::string deck;
    deck += "* mutual-resistance fit/emit verification (abt #76)\n";
    deck += subckt + "\n";
    deck += "I1 0 in1 DC " + std::to_string(I0) + "\n";
    deck += "X1 in1 0 out2 0 " + name + "\n";
    deck += ".tran 1u 20u\n";
    deck += ".end\n";
    auto result = Kirchhoff::run_ngspice_in_process(deck);
    INFO("ngspice error: " << result.error);
    REQUIRE(result.success);
    auto v = result.average("v(out2)", 5e-6, 2e-5);
    REQUIRE(v.has_value());
    return *v / I0;
}

// Parse the gain (last token) of an emitted element line whose first token is `elementName`.
// Returns NaN if the element is not present.
double parse_element_gain(const std::string& subckt, const std::string& elementName) {
    std::istringstream in(subckt);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string first;
        ls >> first;
        if (first != elementName) continue;
        std::string tok, last;
        while (ls >> tok) last = tok;
        try {
            return std::stod(last);
        } catch (...) {
            return std::nan("");
        }
    }
    return std::nan("");
}

}  // namespace

TEST_CASE("mutual resistance: 2-winding uses the PD-safe behavioural realization (abt #76)",
          "[circuit][mutual][export][smoke-test]") {
    auto subckt = export_transformer_subckt({20, 20});
    // The n==2 path is now the SAME PD-safe behavioural realization as n>=3 (abt #76): the
    // old coupled auxiliary winding (LA_/KA_) is gone — it emitted a ladder that did not
    // reproduce its own fitted mutual resistance. The mutual drop is a linear CCVS chain.
    REQUIRE(subckt.find("PD-safe behavioural") != std::string::npos);
    CHECK(subckt.find("Hmut_") != std::string::npos);       // CCVS mutual drop
    CHECK(subckt.find("Vsns_") != std::string::npos);       // winding-current sense
    CHECK(subckt.find("Node_Wtop_") != std::string::npos);  // series routing node
    CHECK(subckt.find("Rdc1 Node_Wtop_1") != std::string::npos);

    // No coupled auxiliary winding and no mutual ladder inductor: the behavioural realization
    // adds NO inductor and NO coupling, so the Lmag matrix stays positive-definite.
    CHECK(subckt.find("LA_") == std::string::npos);
    CHECK(subckt.find("KA_") == std::string::npos);
    CHECK(subckt.find("Lmut_") == std::string::npos);
    CHECK(subckt.find("Fmut_") == std::string::npos);
}

TEST_CASE("mutual resistance: 2-winding emitted Hmut reproduces the fitted model at DC and HF (abt #76)",
          "[circuit][mutual][export][smoke-test]") {
    // The behavioural realization emits, per winding pair, Hmut CCVS whose gain is the
    // fitted mutual resistance R_ij evaluated at the EXPORT frequency:
    //     R_ij(f) = R0 + parallel(L1, R1 + parallel(L2, R2 + parallel(L3, R3)))
    // with R0 = dcMutualResistance the SERIES DC term (coefficients = [R1,L1,R2,L2,R3,L3]).
    // We verify the emitted gain against the fit's two closed-form limits, computed
    // independently from the public coefficients:
    //   * f -> 0  : every parallel(Lk, .) -> 0 (inductor short)  => R_ij -> R0 = dcMutualResistance
    //   * f -> inf: every parallel(Lk, X) -> X (inductor open)   => R_ij -> R0 + R1 + R2 + R3
    // The OLD auxiliary ladder failed the DC limit (its DC impedance was the last-stage R,
    // not dcMutualResistance) — this test pins that the emission now matches its own fit.
    auto magnetic = [] {
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto m = OpenMagneticsTesting::get_quick_magnetic("ETD 39", gapping, {20, 20}, 1, "3C97");
        auto coil = m.get_coil();
        coil.wind();
        m.set_coil(coil);
        return m;
    }();

    // Same fit the exporter uses internally (deterministic levmar, fixed seeds).
    auto coeffs = CircuitSimulatorExporter::calculate_mutual_resistance_coefficients(magnetic, 25.0);
    REQUIRE(coeffs.size() == 1);                 // exactly one pair for 2 windings
    REQUIRE(coeffs[0].coefficients.size() >= 6);
    double dc = coeffs[0].dcMutualResistance;
    double hf = dc + coeffs[0].coefficients[0] + coeffs[0].coefficients[2] + coeffs[0].coefficients[4];
    INFO("fitted dcMutualResistance = " << dc << " Ohm; HF asymptote = " << hf << " Ohm");

    // DC limit: export at a very low frequency so R_ij -> R0.
    auto subLow = export_transformer_subckt({20, 20}, 1.0);
    double gainLow = parse_element_gain(subLow, "Hmut_1_2");
    REQUIRE(std::isfinite(gainLow));
    INFO("emitted Hmut_1_2 @ 1 Hz     = " << gainLow << " Ohm (expect ~ dcMutualResistance)");
    CHECK(std::abs(gainLow - dc) <= 1e-3 * std::abs(dc) + 1e-9);

    // The symmetric drop in the other winding carries the same R_ij.
    double gainLow21 = parse_element_gain(subLow, "Hmut_2_1");
    REQUIRE(std::isfinite(gainLow21));
    CHECK(std::abs(gainLow21 - gainLow) <= 1e-9);

    // HF asymptote: export at a very high frequency so R_ij -> R0 + sum(Rk).
    auto subHigh = export_transformer_subckt({20, 20}, 1e12);
    double gainHigh = parse_element_gain(subHigh, "Hmut_1_2");
    REQUIRE(std::isfinite(gainHigh));
    INFO("emitted Hmut_1_2 @ 1e12 Hz  = " << gainHigh << " Ohm (expect ~ R0 + R1 + R2 + R3)");
    CHECK(std::abs(gainHigh - hf) <= 1e-2 * std::abs(hf) + 1e-9);

    // End-to-end: run BOTH emitted subcircuits through the in-process libngspice engine and
    // measure the realized mutual trans-impedance (winding 2 open, DC drive in winding 1 —
    // see measure_mutual_transimpedance_ngspice). The simulated R12 must match the fitted
    // model at the low-frequency and HF-asymptote export points. The old auxiliary ladder
    // realized a DIFFERENT transfer function than the fit — this is the abt #76 regression.
    double simLow = measure_mutual_transimpedance_ngspice(subLow, 1e4);
    INFO("ngspice-simulated R12 (low-frequency export)  = " << simLow << " Ohm");
    CHECK(std::abs(simLow - dc) <= 1e-2 * std::abs(dc) + 1e-9);

    double simHigh = measure_mutual_transimpedance_ngspice(subHigh, 1.0);
    INFO("ngspice-simulated R12 (HF-asymptote export)   = " << simHigh << " Ohm");
    CHECK(std::abs(simHigh - hf) <= 1e-2 * std::abs(hf) + 1e-9);
}

TEST_CASE("mutual resistance: 3-winding uses the PD-safe behavioural realization", "[circuit][mutual][export][smoke-test]") {
    auto subckt = export_transformer_subckt({20, 20, 20});

    // The cross-coupling term is KEPT for n>=3, not skipped ...
    REQUIRE(subckt.find("PD-safe behavioural") != std::string::npos);
    CHECK(subckt.find("MUTUAL RESISTANCE NETWORK: SKIPPED") == std::string::npos);

    // ... via linear current-controlled voltage sources (Hmut) driven by a
    // per-winding current sense (Vsns), routed through the series top node.
    CHECK(subckt.find("Hmut_") != std::string::npos);       // CCVS mutual drop
    CHECK(subckt.find("Vsns_") != std::string::npos);       // winding-current sense
    CHECK(subckt.find("Node_Wtop_") != std::string::npos);  // series routing node

    // PD-safety BY CONSTRUCTION: the mutual network adds NO inductor and NO coupling.
    // The only magnetically-coupled inductors remain the Lmag windings (K12/K13/K23);
    // there must be no auxiliary winding (LA_/KA_) and no mutual ladder inductor (Lmut_)
    // — adding coupling/inductance is exactly what made the n>=3 deck non-PD.
    CHECK(subckt.find("LA_") == std::string::npos);
    CHECK(subckt.find("KA_") == std::string::npos);
    CHECK(subckt.find("Lmut_") == std::string::npos);
    CHECK(subckt.find("Fmut_") == std::string::npos);

    // The winding self-resistance is routed through the behavioural top node.
    CHECK(subckt.find("Rdc1 Node_Wtop_1") != std::string::npos);
}
