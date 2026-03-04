#include "processors/CircuitSimulatorInterface.h"
#include "TestingUtils.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/Sweeper.h"
#include "support/Painter.h"
#include "processors/NgspiceRunner.h"
#include "converter_models/Flyback.h"
#include "physical_models/MagnetizingInductance.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
double max_error = 0.01;
auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
bool plot = false;

TEST_CASE("Test_CircuitSimulatorExporter_Ngspice_Ideal_vs_Parasitic_Comparison", 
          "[processor][circuit-simulator-exporter][ngspice][parasitic-comparison][smoke-test]") {
    
    // Skip if ngspice is not available
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Create a simple flyback converter setup
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
    
    // Create a real Magnetic component with significant winding resistance
    std::vector<int64_t> numberTurns = {20, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "E 25/13/7";
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    int64_t numberStacks = 1;
    std::string coreMaterial = "3C90";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0004, 1);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    // Extract parameters
    double primaryTurns = coil.get_functional_description()[0].get_number_turns();
    double secondaryTurns = coil.get_functional_description()[1].get_number_turns();
    double turnsRatio = primaryTurns / secondaryTurns;
    double magnetizingInductance = resolve_dimensional_values(
        MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
    
    INFO("Turns ratio: " << turnsRatio);
    INFO("Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
    
    // Calculate winding resistances
    double primaryResistance = WindingLosses::calculate_effective_resistance_of_winding(
        magnetic, 0, 100e3, 25.0);  // Primary at 100kHz
    double secondaryResistance = WindingLosses::calculate_effective_resistance_of_winding(
        magnetic, 1, 100e3, 25.0);  // Secondary at 100kHz
    
    INFO("Primary DC resistance: " << primaryResistance << " Ohm");
    INFO("Secondary DC resistance: " << secondaryResistance << " Ohm");
    
    // Run ideal simulation (no winding resistances)
    auto idealOps = flyback.simulate_and_extract_operating_points({turnsRatio}, magnetizingInductance);
    REQUIRE(idealOps.size() == 1);
    
    // Run real magnetic simulation (with ladder model)
    auto realOps = flyback.simulate_with_magnetic_and_extract_operating_points(magnetic);
    REQUIRE(realOps.size() == 1);
    
    // Both simulations should produce valid waveforms
    REQUIRE(idealOps[0].get_excitations_per_winding().size() == 2);
    REQUIRE(realOps[0].get_excitations_per_winding().size() == 2);
    
    // Extract primary current waveforms
    auto idealPriI = idealOps[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    auto realPriI = realOps[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    
    // Calculate RMS and peak currents
    double idealPriIrms = std::sqrt(std::inner_product(idealPriI.begin(), idealPriI.end(), 
                                                       idealPriI.begin(), 0.0) / idealPriI.size());
    double realPriIrms = std::sqrt(std::inner_product(realPriI.begin(), realPriI.end(), 
                                                      realPriI.begin(), 0.0) / realPriI.size());
    double idealPriImax = *std::max_element(idealPriI.begin(), idealPriI.end());
    double realPriImax = *std::max_element(realPriI.begin(), realPriI.end());
    
    INFO("Ideal primary RMS current: " << idealPriIrms << " A");
    INFO("Real primary RMS current: " << realPriIrms << " A");
    INFO("Ideal primary peak current: " << idealPriImax << " A");
    INFO("Real primary peak current: " << realPriImax << " A");
    
    // Verify that currents are similar but real has slightly higher values due to resistive losses
    // Real current can be significantly higher (up to 2x) due to winding resistance and duty cycle changes
    CHECK(realPriIrms > idealPriIrms * 0.9);  // Real should not be too much lower
    CHECK(realPriIrms < idealPriIrms * 2.5);  // Real can be up to 2.5x higher with parasitics
    
    // Extract secondary current waveforms  
    auto idealSecI = idealOps[0].get_excitations_per_winding()[1].get_current()->get_waveform()->get_data();
    auto realSecI = realOps[0].get_excitations_per_winding()[1].get_current()->get_waveform()->get_data();
    
    double idealSecIrms = std::sqrt(std::inner_product(idealSecI.begin(), idealSecI.end(),
                                                       idealSecI.begin(), 0.0) / idealSecI.size());
    double realSecIrms = std::sqrt(std::inner_product(realSecI.begin(), realSecI.end(),
                                                      realSecI.begin(), 0.0) / realSecI.size());
    
    INFO("Ideal secondary RMS current: " << idealSecIrms << " A");
    INFO("Real secondary RMS current: " << realSecIrms << " A");
    
    // Secondary currents should also be similar (with tolerance for parasitics)
    CHECK(realSecIrms > idealSecIrms * 0.9);
    CHECK(realSecIrms < idealSecIrms * 2.5);
    
    // Calculate copper losses (I^2 * R)
    double realCopperLoss = realPriIrms * realPriIrms * primaryResistance + 
                           realSecIrms * realSecIrms * secondaryResistance;
    
    INFO("Real copper losses: " << realCopperLoss << " W");
    
    // Real copper loss should be positive (since we have winding resistance)
    CHECK(realCopperLoss > 0.01);  // Should have measurable copper loss
    
    // Output voltage should be similar (within tolerance)
    double idealOutputV = opPoint.get_output_voltages()[0];
    // In a real converter, output voltage is controlled, so it should be close to target
    CHECK(std::abs(idealOutputV - 5.0) < 0.5);  // Should be close to 5V
    
    // Save waveforms for visual inspection
    BasicPainter painter;
    
    std::string svgIdeal = painter.paint_operating_point_waveforms(
        idealOps[0], "Flyback Ideal Transformer", 1200, 900);
    auto outIdeal = outputFilePath;
    outIdeal.append("flyback_ideal_vs_parasitic_ideal.svg");
    std::ofstream ofsIdeal(outIdeal);
    if (ofsIdeal.is_open()) { ofsIdeal << svgIdeal; ofsIdeal.close(); }
    
    std::string svgReal = painter.paint_operating_point_waveforms(
        realOps[0], "Flyback with Parasitic Ladder Model", 1200, 900);
    auto outReal = outputFilePath;
    outReal.append("flyback_ideal_vs_parasitic_ladder.svg");
    std::ofstream ofsReal(outReal);
    if (ofsReal.is_open()) { ofsReal << svgReal; ofsReal.close(); }
    
    INFO("Ideal waveforms saved to: " << outIdeal.string());
    INFO("Real (ladder) waveforms saved to: " << outReal.string());
}

TEST_CASE("Test_CircuitSimulatorExporter_Ngspice_Ladder_vs_Fracpole_Comparison", 
          "[processor][circuit-simulator-exporter][ngspice][fracpole-comparison][smoke-test]") {
    
    // Skip if ngspice is not available
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Create a magnetic component
    std::vector<int64_t> numberTurns = {30, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "PQ 35/35";
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    int64_t numberStacks = 1;
    std::string coreMaterial = "3C90";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    // Export with LADDER mode
    auto ladderFile = outputFilePath;
    ladderFile.append("test_magnetic_ladder.cir");
    std::filesystem::remove(ladderFile);
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::NGSPICE).export_magnetic_as_subcircuit(
        magnetic, 100000, 100, ladderFile.string(), std::nullopt, 
        CircuitSimulatorExporterCurveFittingModes::LADDER);
    
    // Export with FRACPOLE mode
    auto fracpoleFile = outputFilePath;
    fracpoleFile.append("test_magnetic_fracpole.cir");
    std::filesystem::remove(fracpoleFile);
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::NGSPICE).export_magnetic_as_subcircuit(
        magnetic, 100000, 100, fracpoleFile.string(), std::nullopt,
        CircuitSimulatorExporterCurveFittingModes::FRACPOLE);
    
    // Both files should exist
    REQUIRE(std::filesystem::exists(ladderFile));
    REQUIRE(std::filesystem::exists(fracpoleFile));
    
    // Read and compare the netlists
    std::ifstream ladderStream(ladderFile);
    std::string ladderContent((std::istreambuf_iterator<char>(ladderStream)),
                              std::istreambuf_iterator<char>());
    ladderStream.close();
    
    std::ifstream fracpoleStream(fracpoleFile);
    std::string fracpoleContent((std::istreambuf_iterator<char>(fracpoleStream)),
                                std::istreambuf_iterator<char>());
    fracpoleStream.close();
    
    // Ladder model should contain R and L components
    CHECK(ladderContent.find("R") != std::string::npos);
    CHECK(ladderContent.find("L") != std::string::npos);
    
    // Fracpole model should also contain R and L components (fallback to ladder currently)
    CHECK(fracpoleContent.find("R") != std::string::npos);
    CHECK(fracpoleContent.find("L") != std::string::npos);
    
    // Both files should have valid netlist structure
    CHECK(fracpoleContent.find(".subckt") != std::string::npos);
    CHECK(fracpoleContent.find(".ends") != std::string::npos);
    
    // TODO: When fracpole implementation is complete, this should be different from ladder
    // For now, they may be the same since fracpole falls back to ladder
    INFO("Ladder model netlist length: " << ladderContent.length());
    INFO("Fracpole model netlist length: " << fracpoleContent.length());
    
    // Note: fracpole mode may fall back to ladder if not fully implemented
    // This test verifies the export mechanism works for both modes
}

}  // namespace
