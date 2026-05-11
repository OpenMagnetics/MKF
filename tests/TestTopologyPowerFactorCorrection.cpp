// ─────────────────────────────────────────────────────────────────────────────
// TestTopologyPowerFactorCorrection.cpp
//
// Per CONVERTER_MODELS_GOLDEN_GUIDE.md §5.0 / §5.1, every converter topology
// must expose `simulate_and_extract_topology_waveforms()` returning
// `ConverterWaveforms` whose external-port signals satisfy the appropriate
// converter-port gate.
//
// PFC is the one topology in the family whose input port is fundamentally AC
// (rectified line voltage), so the standard `check_dc_ports` helper does NOT
// apply to its input. A PFC-specific helper, `check_pfc_ports`, validates:
//   - input_voltage line-cycle MEAN ≈ 0.9·Vrms and RMS ≈ Vrms (rectified |sin|)
//   - output_voltage MEAN ≈ Vbus_nominal (DC bus)
//   - output_current MEAN ≈ Iload_nominal (DC load)
//
// This file mirrors the structure of TestTopologyBuck.cpp / TestTopologyBoost.cpp
// for the standard golden-tier set: design-requirements, operating-point
// generation, and the §5.1 converter-port gate.
// ─────────────────────────────────────────────────────────────────────────────

#include "converter_models/PowerFactorCorrection.h"
#include "ConverterPortChecks.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

using namespace OpenMagnetics;

namespace {

// Universal PFC nameplate (230 Vrms / 50 Hz / 400 V bus / 500 W).
constexpr double kVinRmsNom = 230.0;
constexpr double kVinRmsMin = 185.0;
constexpr double kVinRmsMax = 265.0;
constexpr double kVbus      = 400.0;
constexpr double kPout      = 500.0;
constexpr double kFsw       = 100e3;
constexpr double kFline     = 50.0;
constexpr double kEta       = 0.95;
constexpr double kRipple    = 0.3;

OpenMagnetics::PowerFactorCorrection make_default_pfc(bool sweepInputVoltage = true) {
    OpenMagnetics::PowerFactorCorrection pfc;
    DimensionWithTolerance iv;
    iv.set_nominal(kVinRmsNom);
    if (sweepInputVoltage) {
        iv.set_minimum(kVinRmsMin);
        iv.set_maximum(kVinRmsMax);
    }
    pfc.set_input_voltage(iv);
    pfc.set_output_voltage(kVbus);
    pfc.set_output_power(kPout);
    pfc.set_switching_frequency(kFsw);
    pfc.set_line_frequency(kFline);
    pfc.set_efficiency(kEta);
    pfc.set_current_ripple_ratio(kRipple);
    pfc.set_mode(PfcModes::CONTINUOUS_CONDUCTION_MODE);
    pfc.set_ambient_temperature(25.0);
    return pfc;
}

// ─────────────────────────────────────────────────────────────────────────────
// Design requirements: process_design_requirements() must produce a finite,
// positive boost-inductor inductance window (and a single primary-side winding).
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_Design_Requirements",
          "[converter-model][pfc-topology][design-requirements]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
    auto dr  = pfc.process_design_requirements();

    REQUIRE((dr.get_magnetizing_inductance().get_nominal().has_value() ||
             dr.get_magnetizing_inductance().get_minimum().has_value()));

    const double L = dr.get_magnetizing_inductance().get_nominal().value_or(
                     dr.get_magnetizing_inductance().get_minimum().value_or(0.0));
    INFO("PFC computed boost inductance: " << L * 1e6 << " µH");
    CHECK(L > 0.0);
    CHECK(std::isfinite(L));
}

// ─────────────────────────────────────────────────────────────────────────────
// Operating point generation: process_operating_points() must return at least
// one OperatingPoint with non-empty current/voltage waveforms on a valid time
// grid spanning at least one mains period.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_OperatingPoints_Generation",
          "[converter-model][pfc-topology][operating-points]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
    pfc.set_num_periods_to_extract(1);
    const double L = pfc.calculate_inductance_ccm();
    REQUIRE(L > 0.0);

    auto ops = pfc.process_operating_points({}, L);
    REQUIRE(!ops.empty());

    const auto& exc = ops[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 1);

    auto iWfOpt = exc[0].get_current().value().get_waveform();
    auto vWfOpt = exc[0].get_voltage().value().get_waveform();
    REQUIRE(iWfOpt.has_value());
    REQUIRE(vWfOpt.has_value());
    const auto& iWf = iWfOpt.value();
    const auto& vWf = vWfOpt.value();
    REQUIRE(!iWf.get_data().empty());
    REQUIRE(!vWf.get_data().empty());
    REQUIRE(iWf.get_time().has_value());
    REQUIRE(vWf.get_time().has_value());

    const auto iTime = iWf.get_time().value();
    const double tSpan = iTime.back() - iTime.front();
    INFO("Operating-point time span: " << tSpan * 1e3 << " ms (≥ "
         << (1.0 / kFline) * 1e3 << " ms expected for 1 line cycle)");
    CHECK(tSpan >= 0.99 / kFline);
}

// ─────────────────────────────────────────────────────────────────────────────
// §5.1 converter-port gate (PFC-flavoured).
//
// Sweeps nominal/min/max input RMS; validates that on each ConverterWaveforms:
//   - rectified input voltage line-cycle MEAN ≈ 0.9·Vrms, RMS ≈ Vrms
//   - DC bus output voltage MEAN ≈ Vbus
//   - DC load current MEAN ≈ Pout/Vbus
// per ConverterPortChecks::check_pfc_ports.
//
// This is the analytical PFC waveform synthesis (no ngspice call needed); the
// underlying SPICE netlist is exercised separately by TestNgspiceRunner.cpp.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_ConverterPortWaveforms",
          "[converter-port-waveforms][pfc-topology]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/true);
    const double L = pfc.calculate_inductance_ccm();
    REQUIRE(L > 0.0);

    auto wfs = pfc.simulate_and_extract_topology_waveforms(L, /*numberOfCycles=*/1);
    REQUIRE(!wfs.empty());
    REQUIRE(wfs.size() == 3); // nominal, minimum, maximum (collect_input_voltages order)

    const std::vector<double> vinSweep = {kVinRmsNom, kVinRmsMin, kVinRmsMax};
    const double iLoadNom = kPout / kVbus; // 1.25 A

    for (size_t i = 0; i < wfs.size(); ++i) {
        ConverterPortChecks::check_pfc_ports(wfs[i], "PFC", i,
                                             vinSweep[i], kVbus, iLoadNom);
    }
}

} // namespace
