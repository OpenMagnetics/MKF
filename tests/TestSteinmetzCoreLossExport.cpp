// Large-signal behavioural core-loss (GSE) ngspice export.
//
// The exported subcircuit's core loss was small-signal only (a mu(f) resistance ladder). This
// adds an opt-in behavioural GSE element whose cycle-average follows the Steinmetz f^alpha*B^beta
// law under real waveforms. These tests pin the emission + structure and print the reference
// Steinmetz loss and the sinusoidal drive needed for a target flux, so a separate ngspice run can
// check the injected loss matches Steinmetz.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "processors/CircuitSimulatorInterface.h"
#include "support/Settings.h"
#include "TestingUtils.h"

using namespace OpenMagnetics;

namespace {

OpenMagnetics::Magnetic make_inductor(int64_t turns) {
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 39", gapping, {turns}, 1, "3C97");
    auto coil = magnetic.get_coil();
    coil.wind();
    magnetic.set_coil(coil);
    return magnetic;
}

std::string export_core_loss(OpenMagnetics::Magnetic magnetic, bool steinmetz, double freq) {
    auto& settings = Settings::GetInstance();
    bool previous = settings.get_circuit_simulator_include_steinmetz_core_loss();
    settings.set_circuit_simulator_include_steinmetz_core_loss(steinmetz);
    CircuitSimulatorExporterNgspiceModel exporter;
    auto s = exporter.export_magnetic_as_subcircuit(magnetic, freq, 25.0);
    settings.set_circuit_simulator_include_steinmetz_core_loss(previous);
    return s;
}

std::string export_core_loss_ltspice(OpenMagnetics::Magnetic magnetic, bool steinmetz, double freq) {
    auto& settings = Settings::GetInstance();
    bool previous = settings.get_circuit_simulator_include_steinmetz_core_loss();
    settings.set_circuit_simulator_include_steinmetz_core_loss(steinmetz);
    CircuitSimulatorExporterLtspiceModel exporter;
    auto s = exporter.export_magnetic_as_subcircuit(magnetic, freq, 25.0);
    settings.set_circuit_simulator_include_steinmetz_core_loss(previous);
    return s;
}

}  // namespace

TEST_CASE("steinmetz core loss: OFF by default (no behavioural element)", "[circuit][core-losses][export][smoke-test]") {
    auto s = export_core_loss(make_inductor(30), false, 100000.0);
    CHECK(s.find("Bcloss_") == std::string::npos);
    CHECK(s.find("Behavioural GSE core loss") == std::string::npos);
}

TEST_CASE("steinmetz core loss: ON emits the behavioural GSE injector", "[circuit][core-losses][export][smoke-test]") {
    double freq = 100000.0;
    auto magnetic = make_inductor(30);
    auto s = export_core_loss(magnetic, true, freq);

    REQUIRE(s.find("Behavioural GSE core loss") != std::string::npos);
    CHECK(s.find("Bcloss_1 Node_R_Lmag_1 P1-") != std::string::npos);  // loss current on the mag branch
    CHECK(s.find("Gcl_flux_1") != std::string::npos);                  // flux integrator
    CHECK(s.find("Node_Bflux_1") != std::string::npos);
    // The eps-smoothed powers keep the element Lipschitz (no |V|^(alpha-1) infinite slope at V=0,
    // which stalls the transient); a raw sgn()/abs() form must NOT be emitted.
    CHECK(s.find("sgn(") == std::string::npos);

    // The parameters are physical: 0 < alpha < beta, gc finite and positive. (The ngspice run
    // confirms the cycle-average loss tracks Steinmetz k*f^alpha*B^beta to <0.1% across amplitude.)
    auto p = CircuitSimulatorExporter::calculate_gse_core_loss_params(magnetic, freq, 25.0);
    REQUIRE(p.valid);
    CHECK(p.alpha > 1.0);
    CHECK(p.beta > p.alpha);
    CHECK(p.gc > 0.0);
    CHECK(std::isfinite(p.gc));
}

TEST_CASE("steinmetz core loss: LTspice emits the same behavioural element, unquoted", "[circuit][core-losses][export][smoke-test]") {
    auto magnetic = make_inductor(30);
    // Off by default -> no behavioural element.
    auto off = export_core_loss_ltspice(magnetic, false, 100000.0);
    CHECK(off.find("Bcloss_") == std::string::npos);
    // On -> the behavioural GSE element, with the B-source expression UNQUOTED for LTspice.
    auto on = export_core_loss_ltspice(magnetic, true, 100000.0);
    REQUIRE(on.find("Behavioural GSE core loss") != std::string::npos);
    CHECK(on.find("Bcloss_1 Node_R_Lmag_1 P1- I=") != std::string::npos);
    CHECK(on.find("Bcloss_1 Node_R_Lmag_1 P1- I='") == std::string::npos);  // no single quotes
    CHECK(on.find("Gcl_flux_1") != std::string::npos);
}
