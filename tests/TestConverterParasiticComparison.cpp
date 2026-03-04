#include "processors/CircuitSimulatorInterface.h"
#include "TestingUtils.h"
#include "physical_models/WindingLosses.h"
#include "converter_models/Flyback.h"
#include "converter_models/SingleSwitchForward.h"
#include "converter_models/PushPull.h"
#include "converter_models/Llc.h"
#include "converter_models/Buck.h"
#include "converter_models/Boost.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/NgspiceRunner.h"
#include "support/Painter.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

auto outputFilePath = std::filesystem::path{__FILE__}.parent_path().append("..").append("output");

// Helper to compare ideal vs parasitic for any converter
template<typename ConverterType, typename OperatingPointType>
void compare_ideal_vs_parasitic(
    const std::string& converterName,
    ConverterType& converter,
    const OperatingPointType& opPoint,
    OpenMagnetics::Magnetic magnetic,
    double turnsRatio,
    double magnetizingInductance) {
    
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Calculate winding resistances
    double primaryResistance = 0;
    double secondaryResistance = 0;
    auto coil = magnetic.get_coil();
    if (coil.get_functional_description().size() > 0) {
        primaryResistance = WindingLosses::calculate_effective_resistance_of_winding(
            magnetic, 0, opPoint.get_switching_frequency().value_or(100000), 100);
    }
    if (coil.get_functional_description().size() > 1) {
        secondaryResistance = WindingLosses::calculate_effective_resistance_of_winding(
            magnetic, 1, opPoint.get_switching_frequency().value_or(100000), 100);
    }
    
    INFO(converterName << " - Turns ratio: " << turnsRatio);
    INFO(converterName << " - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
    INFO(converterName << " - Primary resistance: " << primaryResistance << " Ohm");
    
    // Run ideal simulation
    auto idealOps = converter.simulate_and_extract_operating_points({turnsRatio}, magnetizingInductance);
    REQUIRE(idealOps.size() == 1);
    
    // Run real magnetic simulation
    auto realOps = converter.simulate_with_magnetic_and_extract_operating_points(magnetic);
    REQUIRE(realOps.size() == 1);
    
    // Both simulations should produce valid waveforms
    REQUIRE(idealOps[0].get_excitations_per_winding().size() >= 1);
    REQUIRE(realOps[0].get_excitations_per_winding().size() >= 1);
    
    // Extract primary current waveforms
    auto idealPriI = idealOps[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    auto realPriI = realOps[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    
    // Calculate RMS currents
    double idealPriIrms = std::sqrt(std::inner_product(idealPriI.begin(), idealPriI.end(), 
                                                       idealPriI.begin(), 0.0) / idealPriI.size());
    double realPriIrms = std::sqrt(std::inner_product(realPriI.begin(), realPriI.end(), 
                                                      realPriI.begin(), 0.0) / realPriI.size());
    
    INFO(converterName << " - Ideal primary RMS current: " << idealPriIrms << " A");
    INFO(converterName << " - Real primary RMS current: " << realPriIrms << " A");
    
    // Real current should be within reasonable bounds of ideal (allowing for parasitics)
    CHECK(realPriIrms > idealPriIrms * 0.7);  // Real should not be too much lower
    CHECK(realPriIrms < idealPriIrms * 2.5);  // Real can be up to 2.5x higher with parasitics
    
    // Calculate copper losses
    double realCopperLoss = realPriIrms * realPriIrms * primaryResistance;
    INFO(converterName << " - Real copper losses: " << realCopperLoss << " W");
    
    // Real copper loss should be positive if we have winding resistance
    if (primaryResistance > 0.01) {
        CHECK(realCopperLoss > 0.001);
    }
}

// ============================================================================
// Flyback Tests
// ============================================================================

TEST_CASE("Converter_Flyback_Ideal_vs_Parasitic_CCM", "[converter][flyback][ideal-vs-parasitic][smoke-test]") {
    OpenMagnetics::Flyback flyback;
    
    DimensionWithTolerance inputVoltage;
    inputVoltage.set_nominal(24.0);
    flyback.set_input_voltage(inputVoltage);
    flyback.set_diode_voltage_drop(0.5);
    
    OpenMagnetics::FlybackOperatingPoint opPoint;
    opPoint.set_output_voltages({5.0});
    opPoint.set_output_currents({2.0});
    opPoint.set_ambient_temperature(25.0);
    opPoint.set_switching_frequency(100e3);
    flyback.set_operating_points({opPoint});
    
    // Create magnetic
    std::vector<int64_t> numberTurns = {20, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "E 25/13/7");
    auto core = OpenMagneticsTesting::get_quick_core("E 25/13/7", 
        OpenMagneticsTesting::get_distributed_gap(0.0004, 1), 1, "3C90");
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    double turnsRatio = 4.0;
    double magnetizingInductance = 75e-6;  // 75 uH
    
    compare_ideal_vs_parasitic("Flyback_CCM", flyback, opPoint, magnetic, turnsRatio, magnetizingInductance);
}

// ============================================================================
// Additional tests for converters that don't support magnetic comparison
// ============================================================================

TEST_CASE("Converter_Basic_Simulation_Tests", "[converter][basic-simulation][smoke-test]") {
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Test Buck
    {
        OpenMagnetics::Buck buck;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(12.0);
        buck.set_input_voltage(inputVoltage);
        buck.set_diode_voltage_drop(0.5);
        
        auto opPoint = BuckOperatingPoint();
        opPoint.set_output_voltage(5.0);
        opPoint.set_output_current(3.0);
        opPoint.set_ambient_temperature(25.0);
        opPoint.set_switching_frequency(500e3);
        buck.set_operating_points({opPoint});
        
        double inductance = 10e-6;
        auto ops = buck.simulate_and_extract_operating_points(inductance);
        REQUIRE(ops.size() == 1);
        REQUIRE(ops[0].get_excitations_per_winding().size() >= 1);
    }
    
    // Test Boost
    {
        OpenMagnetics::Boost boost;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(5.0);
        boost.set_input_voltage(inputVoltage);
        boost.set_diode_voltage_drop(0.5);
        
        auto opPoint = BoostOperatingPoint();
        opPoint.set_output_voltage(12.0);
        opPoint.set_output_current(1.0);
        opPoint.set_ambient_temperature(25.0);
        opPoint.set_switching_frequency(300e3);
        boost.set_operating_points({opPoint});
        
        double inductance = 22e-6;
        auto ops = boost.simulate_and_extract_operating_points(inductance);
        REQUIRE(ops.size() == 1);
        REQUIRE(ops[0].get_excitations_per_winding().size() >= 1);
    }
    
    INFO("Basic converter simulation tests passed");
}

}  // namespace
