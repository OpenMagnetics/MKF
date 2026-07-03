// Stray/parasitic capacitance ngspice export.
//
// The StrayCapacitance energy method (turn-to-turn -> per-winding) is fully modeled
// in MKF but was never emitted into the SPICE subcircuit. These tests pin the export
// of the per-winding self capacitance (self-resonance) and inter-winding capacitance
// (CM coupling) as lumped positive terminal caps, gated by
// circuit_simulator_include_stray_capacitance (default off).
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

std::string export_with_capacitance(const std::vector<int64_t>& numberTurns, bool includeCaps) {
    auto magnetic = wound_magnetic(numberTurns);
    auto& settings = Settings::GetInstance();
    bool previous = settings.get_circuit_simulator_include_stray_capacitance();
    settings.set_circuit_simulator_include_stray_capacitance(includeCaps);
    CircuitSimulatorExporterNgspiceModel exporter;
    std::string subckt = exporter.export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0);
    settings.set_circuit_simulator_include_stray_capacitance(previous);
    return subckt;
}

std::string export_ltspice_with_capacitance(const std::vector<int64_t>& numberTurns, bool includeCaps) {
    auto magnetic = wound_magnetic(numberTurns);
    auto& settings = Settings::GetInstance();
    bool previous = settings.get_circuit_simulator_include_stray_capacitance();
    settings.set_circuit_simulator_include_stray_capacitance(includeCaps);
    CircuitSimulatorExporterLtspiceModel exporter;
    std::string subckt = exporter.export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0);
    settings.set_circuit_simulator_include_stray_capacitance(previous);
    return subckt;
}

}  // namespace

TEST_CASE("stray capacitance: OFF by default (no Cself_/Cwind_ emitted)", "[circuit][capacitance][export]") {
    auto subckt = export_with_capacitance({20, 20}, false);
    CHECK(subckt.find("STRAY CAPACITANCE NETWORK") == std::string::npos);
    CHECK(subckt.find("Cself_") == std::string::npos);
    CHECK(subckt.find("Cwind_") == std::string::npos);
}

TEST_CASE("stray capacitance: 2-winding emits self + inter-winding caps", "[circuit][capacitance][export]") {
    auto subckt = export_with_capacitance({20, 20}, true);

    REQUIRE(subckt.find("STRAY CAPACITANCE NETWORK") != std::string::npos);
    // Self-capacitance across each winding — the self-resonant element.
    CHECK(subckt.find("Cself_1 P1+ P1-") != std::string::npos);
    CHECK(subckt.find("Cself_2 P2+ P2-") != std::string::npos);
    // Inter-winding capacitance between the two windings — the CM coupling path.
    CHECK(subckt.find("Cwind_1_2 P1+ P2+") != std::string::npos);
}

TEST_CASE("stray capacitance: single inductor emits a self-capacitance", "[circuit][capacitance][export]") {
    auto subckt = export_with_capacitance({30}, true);
    // A one-winding inductor still gets a self-capacitance (its self-resonance) and no
    // inter-winding term.
    CHECK(subckt.find("Cself_1 P1+ P1-") != std::string::npos);
    CHECK(subckt.find("Cwind_") == std::string::npos);
}

// The frontend's SPICE downloads use the LTspice-format exporter, so it must emit the same
// positive 3-cap network as ngspice (shared emit_stray_capacitance_spice).
TEST_CASE("stray capacitance: LTspice OFF by default", "[circuit][capacitance][export]") {
    auto subckt = export_ltspice_with_capacitance({20, 20}, false);
    CHECK(subckt.find("Cself_") == std::string::npos);
    CHECK(subckt.find("Cwind_") == std::string::npos);
}

TEST_CASE("stray capacitance: LTspice 2-winding emits self + inter-winding caps", "[circuit][capacitance][export]") {
    auto subckt = export_ltspice_with_capacitance({20, 20}, true);
    REQUIRE(subckt.find("STRAY CAPACITANCE NETWORK") != std::string::npos);
    CHECK(subckt.find("Cself_1 P1+ P1-") != std::string::npos);
    CHECK(subckt.find("Cself_2 P2+ P2-") != std::string::npos);
    CHECK(subckt.find("Cwind_1_2 P1+ P2+") != std::string::npos);
}

TEST_CASE("stray capacitance: LTspice single inductor emits a self-capacitance", "[circuit][capacitance][export]") {
    auto subckt = export_ltspice_with_capacitance({30}, true);
    CHECK(subckt.find("Cself_1 P1+ P1-") != std::string::npos);
    CHECK(subckt.find("Cwind_") == std::string::npos);
}
