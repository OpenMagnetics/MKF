// 2-winding transformer coupling coefficient ngspice export — abt #56/#61.
//
// The ngspice subcircuit exporter emits the magnetic coupling between the two winding
// magnetizing inductors as a `K Lmag_1 Lmag_2 <k>` statement. That coupling used to be
// HARD-CAPPED at 0.98 "for stability", which injected ~2% ARTIFICIAL series leakage: for a
// well-coupled transformer the field-solver coupling is ~0.999, and forcing it to 0.98
// strangled simulated power transfer (AHB vout capped ~3 V instead of 12 V; PSFB/PSHB decks
// transferred ~0 power). ngspice only needs K strictly below 1 (k=1.0 makes the mutual
// matrix singular), so the cap is now a numerical-stability-only clamp just under 1. This
// test pins that the emitted K reflects the COMPUTED coupling, not the old 0.98 cap.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "processors/CircuitSimulatorInterface.h"
#include "support/Settings.h"
#include "TestingUtils.h"

using namespace OpenMagnetics;

namespace {

// Parse `K Lmag_1 Lmag_2 <value>` — the emitted 2-winding coupling statement.
double parse_coupling_K(const std::string& subckt) {
    std::istringstream in(subckt);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string tok0, tok1, tok2, val;
        ls >> tok0 >> tok1 >> tok2 >> val;
        if (tok0 == "K" && tok1 == "Lmag_1" && tok2 == "Lmag_2") {
            try { return std::stod(val); } catch (...) { return std::nan(""); }
        }
    }
    return std::nan("");
}

// Parse the emitted `.param CouplingCoefficient_12_Value=<value>` — the coupling the field
// solver computed from the measured leakage (independent of the K-statement clamp).
double parse_coupling_param(const std::string& subckt) {
    const std::string key = ".param CouplingCoefficient_12_Value=";
    auto pos = subckt.find(key);
    if (pos == std::string::npos) return std::nan("");
    pos += key.size();
    auto end = subckt.find_first_of("\r\n", pos);
    try { return std::stod(subckt.substr(pos, end - pos)); } catch (...) { return std::nan(""); }
}

}  // namespace

TEST_CASE("2-winding coupling: emitted K reflects the computed coupling, not the old 0.98 cap (abt #61)",
          "[circuit][export][smoke-test]") {
    // A tightly-wound 2-winding transformer on an (essentially ungapped) ETD 39 — a
    // well-coupled part whose field-solver coupling is comfortably above the old 0.98 cap.
    auto gapping = OpenMagneticsTesting::get_residual_gap();
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 39", gapping, {20, 20}, 1, "3C97");
    auto coil = magnetic.get_coil();
    coil.wind();
    magnetic.set_coil(coil);

    CircuitSimulatorExporterNgspiceModel exporter;
    std::string subckt = exporter.export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0);

    double emittedK = parse_coupling_K(subckt);
    double computed = parse_coupling_param(subckt);
    REQUIRE(std::isfinite(emittedK));
    REQUIRE(std::isfinite(computed));
    INFO("emitted K = " << emittedK << ", computed CouplingCoefficient_12 = " << computed);

    // The emitted K must equal the field-solver coupling (the only allowed adjustment is the
    // numerical clamp just under 1) — NOT the discarded 0.98 cap. The K statement is printed
    // via std::to_string (6 decimal places), so compare at that print precision.
    CHECK(std::abs(emittedK - std::min(0.999999, computed)) <= 5e-7);

    // For a well-coupled transformer the computed coupling is far above 0.98, so the emitted
    // K must be too. If this ever regressed to the 0.98 cap, both checks below would fail.
    CHECK(computed > 0.98);
    CHECK(emittedK > 0.98);
    CHECK(std::abs(emittedK - 0.98) > 1e-6);

    // And it must stay strictly below 1 (k=1.0 => singular mutual matrix in ngspice).
    CHECK(emittedK < 1.0);
}
