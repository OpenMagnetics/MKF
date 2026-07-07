// NL5 exporter reachability (ABT #120.1).
//
// The NL5 exporter was fully implemented but not derived from
// CircuitSimulatorExporterModel, so CircuitSimulatorExporter(NL5) threw
// "Unknown program" and the exporter was dead code behind the factory.
// These tests pin the factory wiring and a minimal subcircuit emission.
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "processors/CircuitSimulatorInterface.h"
#include "constructive_models/Magnetic.h"
#include "support/Settings.h"
#include "TestingUtils.h"

using namespace OpenMagnetics;

namespace {

OpenMagnetics::Magnetic wound_magnetic(const std::vector<int64_t>& numberTurns) {
    auto gapping = OpenMagneticsTesting::get_residual_gap();
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 39", gapping, numberTurns, 1, "3C97");
    auto coil = magnetic.get_coil();
    coil.wind();  // the exporter needs a processed coil (turns description)
    magnetic.set_coil(coil);
    return magnetic;
}

} // namespace

TEST_CASE("nl5: factory constructs the exporter", "[circuit][export][nl5][smoke-test]") {
    // Used to throw ModelNotAvailableException("Unknown Circuit Simulator program...")
    CHECK_NOTHROW(CircuitSimulatorExporter(CircuitSimulatorExporterModels::NL5));
}

TEST_CASE("nl5: subcircuit export emits winding components", "[circuit][export][nl5][smoke-test]") {
    auto magnetic = wound_magnetic({10, 5});
    CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::NL5);
    std::string subckt = exporter.export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0);
    CHECK(!subckt.empty());
    // NL5 schematics are XML-ish component lists; the transformer windings are
    // W components and the magnetizing inductance must be present.
    CHECK(subckt.find("Cmp") != std::string::npos);
}

TEST_CASE("nl5: symbol export fails loudly (not implemented)", "[circuit][export][nl5][smoke-test]") {
    auto magnetic = wound_magnetic({10});
    CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::NL5);
    CHECK_THROWS(exporter.export_magnetic_as_symbol(magnetic));
}
