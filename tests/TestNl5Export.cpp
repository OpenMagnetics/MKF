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

// ABT #1456 (Alf's decision): the NL5 exporter built the core fractional-pole network inside
// try{}catch(...){}, so a network that could not be built was silently left out. N22 carries only a loss
// factor, no Steinmetz, and the core fracpole network is anchored on Steinmetz: the FRACPOLE export must
// now throw, while the ladder export (loss-factor core resistance) still works.
TEST_CASE("nl5: fracpole export throws when the core network cannot be built", "[circuit][export][nl5][abt-1456]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 39", OpenMagneticsTesting::get_residual_gap(),
                                                             std::vector<int64_t>{10}, 1, "N22");
    auto coil = magnetic.get_coil();
    coil.wind();
    magnetic.set_coil(coil);
    CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::NL5);
    CHECK_NOTHROW(exporter.export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0));
    CHECK_THROWS(exporter.export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0, std::nullopt, std::nullopt,
                                                        CircuitSimulatorExporterCurveFittingModes::FRACPOLE));
}

// ABT #1623: the part reference is free text written into every export. In NL5 XML it must be
// escaped; in SPICE and PLECS a line break would end the comment and run the rest as netlist or
// script on open (PLECS init commands are executed), and a quote would end PLECS's quoted name.
TEST_CASE("exporters: a hostile reference is escaped or refused, never written raw", "[circuit][export][nl5][smoke-test][abt-1623]") {
    auto magnetic = wound_magnetic({10});
    auto setReference = [&](const std::string& reference) {
        MAS::MagneticManufacturerInfo manufacturerInfo;
        manufacturerInfo.set_name("test");
        manufacturerInfo.set_reference(reference);
        magnetic.set_manufacturer_info(manufacturerInfo);
    };
    setReference("Part <A&B> \"x\"");
    std::string nl5 = CircuitSimulatorExporter(CircuitSimulatorExporterModels::NL5).export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0);
    CHECK(nl5.find("Component: Part &lt;A&amp;B&gt; &quot;x&quot;</c>") != std::string::npos);
    CHECK(nl5.find("<A&B>") == std::string::npos);
    CHECK_THROWS(CircuitSimulatorExporter(CircuitSimulatorExporterModels::PLECS).export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0));

    setReference("Part\nV1 in 0 1000");
    for (auto model : {CircuitSimulatorExporterModels::NGSPICE, CircuitSimulatorExporterModels::LTSPICE, CircuitSimulatorExporterModels::PLECS}) {
        CHECK_THROWS(CircuitSimulatorExporter(model).export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0));
    }

    setReference("Plain part 42");
    for (auto model : {CircuitSimulatorExporterModels::NGSPICE, CircuitSimulatorExporterModels::LTSPICE, CircuitSimulatorExporterModels::PLECS}) {
        CHECK_NOTHROW(CircuitSimulatorExporter(model).export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0));
    }
}
