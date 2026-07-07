// Mutual (cross-coupling) resistance ngspice export — abt #50/#72.
//
// The off-diagonal Hesterman winding-loss term is realized for a 2-winding
// transformer with a coupled auxiliary winding (positive-definite, works). For
// n>=3 that auxiliary realization leaves the coupled-L matrix non-positive-
// definite and ngspice collapses, so it used to be SKIPPED. These tests pin the
// PD-safe BEHAVIOURAL realization that keeps the term for n>=3 without adding any
// coupled inductors (grounded uncoupled ladders + sense/behavioural sources).
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "processors/CircuitSimulatorInterface.h"
#include "support/Settings.h"
#include "TestingUtils.h"

using namespace OpenMagnetics;

namespace {

// Emit an n-winding transformer subcircuit with the mutual-resistance term on.
std::string export_transformer_subckt(const std::vector<int64_t>& numberTurns) {
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
    std::string subckt = exporter.export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0);
    settings.set_circuit_simulator_include_mutual_resistance(previous);
    return subckt;
}

}  // namespace

TEST_CASE("mutual resistance: 2-winding keeps the auxiliary-winding realization", "[circuit][mutual][export][smoke-test]") {
    auto subckt = export_transformer_subckt({20, 20});
    // The 2-winding path is positive-definite already and left untouched (abt #72):
    // it still uses the coupled auxiliary winding LA_ and its KA_ couplings, and it
    // does NOT route the winding through the behavioural series node.
    REQUIRE(subckt.find("MUTUAL RESISTANCE NETWORK") != std::string::npos);
    CHECK(subckt.find("LA_12") != std::string::npos);
    CHECK(subckt.find("KA_") != std::string::npos);
    CHECK(subckt.find("PD-safe behavioural") == std::string::npos);
    CHECK(subckt.find("Node_Wtop_") == std::string::npos);
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
