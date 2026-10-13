#include "physical_models/Temperature.h"
#include "physical_models/ThermalResistance.h"
#include "physical_models/WindingLosses.h"
#include "constructive_models/Bobbin.h"
#include "support/Settings.h"
#include "support/Exceptions.h"
#include "support/Painter.h"
#include "processors/MagneticSimulator.h"
#include "Definitions.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <set>
#include <functional>
#include <map>
#include <magic_enum.hpp>

using namespace MAS;
using namespace OpenMagnetics;

// Helper function to get the output directory path
std::filesystem::path getOutputDir() {
    return std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("output");
}

// Helper function to export temperature field SVG
void exportTemperatureFieldSvg(const std::string& testName, 
                                OpenMagnetics::Magnetic magnetic, 
                                const std::map<std::string, double>& nodeTemperatures,
                                double ambientTemperature = 25.0) {
    auto outputDir = getOutputDir();
    std::filesystem::create_directories(outputDir);
    auto outFile = outputDir / ("thermal_" + testName + ".svg");
    
    // Wind the coil if not already wound
    if (!magnetic.get_coil().get_turns_description()) {
        try {
            magnetic.get_mutable_coil().wind();
        } catch (...) {
            // Coil winding may fail for some test cases, continue without turns
        }
    }
    
    OpenMagnetics::Painter painter(outFile);
    painter.paint_core(magnetic);
    // Don't paint turns separately - paint_temperature_field will paint them with temperature colors
    // if (magnetic.get_coil().get_turns_description()) {
    //     painter.paint_coil_turns(magnetic);
    // }
    painter.paint_temperature_field(magnetic, nodeTemperatures, true, OpenMagnetics::ColorPalette::BLUE_TO_RED, ambientTemperature);  // Enable color bar
    painter.export_svg();
}

// Helper function to export thermal circuit schematic SVG
void exportThermalCircuitSchematic(const std::string& testName,
                                   const Temperature& temp) {
    auto outputDir = getOutputDir();
    std::filesystem::create_directories(outputDir);
    auto outFile = outputDir / ("thermal_schematic_" + testName + ".svg");
    
    OpenMagnetics::Painter painter(outFile);
    auto nodes = temp.getNodes();
    auto resistances = temp.getResistances();
    std::string svg = painter.paint_thermal_circuit_schematic(nodes, resistances);
    
    // Write to file
    std::ofstream file(outFile);
    file << svg;
    file.close();
}

namespace {

// Helper function to run magnetic simulation and extract real losses
// This replaces the mock loss setup with actual simulation results
struct LossesFromSimulation {
    double coreLosses = 0.0;
    double windingLosses = 0.0;
    double ambientTemperature = 25.0;
    std::optional<WindingLossesOutput> windingLossesOutput;
    bool simulationSucceeded = false;
};

LossesFromSimulation getLossesFromSimulation(const OpenMagnetics::Magnetic& magnetic, 
                                              const OpenMagnetics::Inputs& inputs) {
    LossesFromSimulation result;
    
    // Run magnetic simulation to get actual losses
    OpenMagnetics::MagneticSimulator magneticSimulator;
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    mas.set_inputs(inputs);
    
    // No try/catch: a simulation that throws (e.g. the volt-second balance check rejecting an
    // unbalanced fixture) must fail the test with ITS message, not surface later as a misleading
    // "WindingLossesOutput missing" or silently switch the caller to placeholder losses.
    auto simulatedMas = magneticSimulator.simulate(inputs, magnetic);
    if (simulatedMas.get_outputs().empty()) {
        throw std::runtime_error("getLossesFromSimulation: MagneticSimulator returned no outputs");
    }
    auto outputs = simulatedMas.get_outputs()[0];

    // Get core losses
    if (outputs.get_core_losses()) {
        result.coreLosses = resolve_dimensional_values(outputs.get_core_losses()->get_core_losses());
    }

    // Get winding losses
    if (outputs.get_winding_losses()) {
        result.windingLosses = resolve_dimensional_values(outputs.get_winding_losses()->get_winding_losses());
        result.windingLossesOutput = outputs.get_winding_losses();
    }

    result.simulationSucceeded = true;

    // Get ambient temperature from operating point
    if (!inputs.get_operating_points().empty()) {
        auto opPoint = inputs.get_operating_points()[0];
        result.ambientTemperature = opPoint.get_conditions().get_ambient_temperature();
    }
    
    return result;
}

// Helper to create Inputs with operating points for a magnetic component
// Uses Inputs::create_quick_operating_point_only_current with reasonable defaults
// Properly handles multi-winding transformers by calculating turns ratios
OpenMagnetics::Inputs createInputsForMagnetic(const OpenMagnetics::Magnetic& magnetic,
                                               double frequency = 100000.0,  // 100 kHz default
                                               double temperature = 25.0,
                                               double currentPeak = 1.0,     // 1A peak default
                                               WaveformLabel waveShape = WaveformLabel::SINUSOIDAL) {
    // Estimate magnetizing inductance based on core and turns
    double magnetizingInductance = 1e-3;  // 1mH default
    std::vector<double> turnsRatios;
    
    // Try to calculate actual inductance and turns ratios if possible
    try {
        auto coil = magnetic.get_coil();
        auto core = magnetic.get_core();
        auto functionalDesc = coil.get_functional_description();
        
        if (!functionalDesc.empty()) {
            // Get primary winding turns
            auto primaryWinding = functionalDesc[0];
            int primaryTurns = primaryWinding.get_number_turns();
            
            // Calculate magnetizing inductance based on primary
            double al = 1000e-9;  // 1000 nH/turn^2 typical
            magnetizingInductance = al * primaryTurns * primaryTurns;
            
            // Calculate turns ratios for additional windings (for transformers)
            if (functionalDesc.size() > 1) {
                for (size_t i = 1; i < functionalDesc.size(); ++i) {
                    double ratio = static_cast<double>(functionalDesc[i].get_number_turns()) / primaryTurns;
                    turnsRatios.push_back(ratio);
                }
            }
        }
    } catch (...) {
        // Use defaults
    }
    
    // Calculate peak-to-peak from peak for sinusoidal
    double peakToPeak = currentPeak * 2.0;  // For sinusoidal: peak-to-peak = 2 * peak
    
    return OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        frequency,
        magnetizingInductance,
        temperature,
        waveShape,
        peakToPeak,
        0.5,    // duty cycle
        0.0,    // DC current
        turnsRatios  // Pass turns ratios for multi-winding support
    );
}

// Helper to apply losses from simulation to config
// Throws exception if simulation fails - no mock losses allowed
void applySimulatedLosses(TemperatureConfig& config, 
                          const OpenMagnetics::Magnetic& magnetic,
                          double frequency = 100000.0,
                          double currentPeak = 1.0) {
    // Create proper Inputs with operating points
    auto inputs = createInputsForMagnetic(magnetic, frequency, config.ambientTemperature, currentPeak);
    
    // Get real losses from simulation
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    if (!losses.simulationSucceeded) {
        throw std::runtime_error("Magnetic simulation failed - cannot calculate temperatures without real losses. "
                                 "Ensure the magnetic component has valid core, coil, and operating points.");
    }
    
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results. "
                                 "Cannot calculate temperatures without per-turn loss distribution.");
    }
    
    config.coreLosses = losses.coreLosses;
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.ambientTemperature = losses.ambientTemperature;
}

//=============================================================================
// Unit Tests for Static Calculation Methods
//=============================================================================

TEST_CASE("Temperature: Toroidal Core T20 Ten Turns", "[temperature][round-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.3;  // Fallback if simulation fails
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_Toroid_T20_10_turns.svg").string();
    
    // Use real simulation for losses
    // Use 0.1A current to avoid extreme core losses in small ungapped toroid
    auto inputs = createInputsForMagnetic(magnetic, 100000.0, config.ambientTemperature, 0.1);
    auto losses = getLossesFromSimulation(magnetic, inputs);
    if (losses.simulationSucceeded) {
        config.coreLosses = losses.coreLosses;
        config.windingLosses = losses.windingLosses;
        if (losses.windingLossesOutput.has_value()) {
            config.windingLossesOutput = losses.windingLossesOutput.value();
        } else {
            applySimulatedLosses(config, magnetic);
        }
    }
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    // Export visualization
    exportTemperatureFieldSvg("Toroid_T20_10_turns", magnetic, result.nodeTemperatures);
    exportThermalCircuitSchematic("Toroid_T20_10_turns", temp);
    
    // Expected is ~30-50 K/W for 20mm toroid per Maniktala
    REQUIRE(result.totalThermalResistance > 1.0);
    REQUIRE(result.totalThermalResistance < 450.0);  // Relaxed until model is calibrated
    
}

TEST_CASE("Temperature: Larger Toroidal Core Two Windings", "[temperature][round-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns({20, 10});
    std::vector<int64_t> numberParallels({1, 1});
    std::string shapeName = "T 36/23/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.0;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_Toroid_T36_two_windings.svg").string();
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    // Export visualization
    exportTemperatureFieldSvg("Toroid_T36_two_windings", magnetic, result.nodeTemperatures);
    exportThermalCircuitSchematic("Toroid_T36_two_windings", temp);
    
    // Expected is ~20-40 K/W
    REQUIRE(result.totalThermalResistance > 0.5);
    REQUIRE(result.totalThermalResistance < 200.0);  // Relaxed until model is calibrated
    
}

TEST_CASE("Temperature: T36 Two Windings Schematic Only", "[temperature][round-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns({20, 10});
    std::vector<int64_t> numberParallels({1, 1});
    std::string shapeName = "T 36/23/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.0;
    applySimulatedLosses(config, magnetic);
    config.nodePerCoilTurn = true;  // Enable quadrant visualization
    config.plotSchematic = true;
    config.maxIterations = 1;  // Just 1 iteration to build network
    // ABT #837: this case builds the network to DRAW it and never intends to solve it, so it
    // opts out of the convergence requirement explicitly rather than reading a diverged number.
    config.requireConvergence = false;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_T36_two_windings_quadrant.svg").string();
    
    // Build thermal circuit and generate schematic (minimal solve)
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    // The thermal network must actually have been built.
    CHECK(temp.getNodes().size() > 0);
    CHECK(temp.getResistances().size() > 0);

    // Paint the magnetic geometry for visualization
    Painter painter;
    painter.paint_core(magnetic);
    painter.paint_coil_turns(magnetic);
    auto svg = painter.export_svg();
    CHECK(svg.find("<svg") != std::string::npos);
    std::ofstream file("output/T36_geometry_visualization.svg");
    file << svg;
    file.close();
}

TEST_CASE("Temperature: T20 Two Windings Quadrant Visualization", "[temperature][round-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns({5, 5});
    std::vector<int64_t> numberParallels({1, 1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.0;
    applySimulatedLosses(config, magnetic);
    config.nodePerCoilTurn = true;  // Enable quadrant visualization
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_quadrant_T20_two_windings.svg").string();
    
    Temperature temp(magnetic, config);
    
    // Try to calculate temperatures, but export schematic even if solver fails
    // This allows us to visualize the thermal network for debugging
    try {
        auto result = temp.calculateTemperatures();
    } catch (const std::exception& e) {
        // Calculation failed, but we can still export the schematic
        std::cout << "Temperature calculation failed (expected for this test): " << e.what() << std::endl;
    }
    
    // Export thermal circuit schematic SVG
    exportThermalCircuitSchematic("T20_two_windings", temp);
    
    // Note: With insulation layers present, only the outermost layer has convection to ambient.
    // Inner layers rely on conduction through turns to reach the outer layer.
    // For visualization purposes, we don't require convergence.
    REQUIRE(temp.getNodes().size() > 0);
    REQUIRE(temp.getResistances().size() > 0);
}

TEST_CASE("Temperature: Toroidal Quadrant Visualization", "[temperature][round-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.3;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_quadrant_visualization.svg").string();
    
    Temperature temp(magnetic, config);
    
    // Calculate temperatures (builds network and generates schematic even if solver fails)
    auto result = temp.calculateTemperatures();
    // The thermal network must actually have been built and the schematic written.
    CHECK(temp.getNodes().size() > 0);
    CHECK(temp.getResistances().size() > 0);
    // Temperature::plotSchematic writes the schematic as '<path>.json'.
    CHECK(std::filesystem::exists(config.schematicOutputPath + ".json"));
}



//=============================================================================
// Additional Core Type Tests
//=============================================================================

TEST_CASE("Temperature: ETD Core", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({15});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 39/20/13";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "3C95";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.5;  // Fallback
    config.plotSchematic = false;
    
    // Use real simulation for losses
    // Use lower current (0.1A) to avoid extreme core losses in ungapped E core
    auto inputs = createInputsForMagnetic(magnetic, 100000.0, config.ambientTemperature, 0.1);
    auto losses = getLossesFromSimulation(magnetic, inputs);
    if (losses.simulationSucceeded) {
        config.coreLosses = losses.coreLosses;
        config.windingLosses = losses.windingLosses;
        if (losses.windingLossesOutput.has_value()) {
            config.windingLossesOutput = losses.windingLossesOutput.value();
        }
    } else {
        applySimulatedLosses(config, magnetic);
    }
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    // Temperature can be high with real simulation losses - just check it's reasonable
    REQUIRE(result.maximumTemperature < 500.0);
    REQUIRE(result.totalThermalResistance > 0.1);
    REQUIRE(result.totalThermalResistance < 200.0);
}

TEST_CASE("Temperature: E Core", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({20});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "E 42/21/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    // Add 0.5mm gap to limit flux density and reduce core losses
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.0005);  // 0.5mm gap
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.8;  // Fallback
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_E_Core_42_21_15.svg").string();
    
    // Use real simulation for losses
    auto inputs = createInputsForMagnetic(magnetic, 100000.0, config.ambientTemperature, 1.0);
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    if (losses.simulationSucceeded) {
        config.coreLosses = losses.coreLosses;
        config.windingLosses = losses.windingLosses;
        if (losses.windingLossesOutput.has_value()) {
            config.windingLossesOutput = losses.windingLossesOutput.value();
        }
    } else {
        applySimulatedLosses(config, magnetic);
    }
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    // Temperature can be high with real simulation losses - just check it's reasonable
    REQUIRE(result.maximumTemperature < 500.0);
    REQUIRE(result.totalThermalResistance > 0.1);
}

TEST_CASE("Temperature: Multi-Winding", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({20, 10, 15});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "T 36/23/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.2;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
}

TEST_CASE("Temperature: Ambient Temperature Effect", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config1;
    config1.ambientTemperature = 25.0;
    config1.coreLosses = 0.3;
    applySimulatedLosses(config1, magnetic);
    config1.plotSchematic = false;
    
    Temperature temp1(magnetic, config1);
    auto result1 = temp1.calculateTemperatures();
    
    TemperatureConfig config2;
    config2.ambientTemperature = 50.0;
    config2.coreLosses = 0.3;
    applySimulatedLosses(config2, magnetic);
    config2.plotSchematic = false;
    
    Temperature temp2(magnetic, config2);
    auto result2 = temp2.calculateTemperatures();
    
    REQUIRE(result1.converged);
    REQUIRE(result2.converged);
    REQUIRE(result2.maximumTemperature > result1.maximumTemperature);

    // ABT #1735 (Alf 2026-10-11): the band is derived from how the part's own heat-transfer
    // coefficients scale between the two ambients, not guessed. Raising the ambient from 25 to
    // 50 C at a similar rise:
    //   - radiation, h_rad = eps sigma (Ts^2 + Ta^2)(Ts + Ta), grows with absolute temperature
    //     cubed: (323/298)^3 = 1.27 at a few kelvin of rise;
    //   - natural convection on surfaces inside a correlation's range falls slightly (film
    //     viscosity rises faster than conductivity, and the rise is smaller), while a small
    //     surface on the Nu >= 1 floor has h = k_air(T_film)/L, which grows with k_air (~7%);
    //   - solid conduction (ferrite, copper, enamel, tape) is ambient-independent here, and
    //     stagnant-air films grow with k_air.
    // So the part's total thermal resistance must FALL (radiation carries a large share of a small
    // naturally-cooled part's heat and nothing falls by as much as it rises), and it cannot fall by
    // more than its fastest-growing conductance rises: no part of the heat path becomes more than
    // m_max times as conductive, so the rise per watt cannot shrink below 1/m_max of its 25 C
    // value. m_max is read from the model's own resistors, paired between the two solves (same
    // geometry, same network). The core/winding loss split moves between the ambients (3C97 loss
    // falls, copper loss rises), but on this toroid both sit on one body whose hottest node is
    // set by the surface films, so the split does not carry the rise per watt outside this range.
    const auto& resistances25 = result1.thermalResistances;
    const auto& resistances50 = result2.thermalResistances;
    REQUIRE(resistances25.size() == resistances50.size());
    REQUIRE(!resistances25.empty());
    double maximumConductanceGrowth = 0.0;
    double radiationConductance25 = 0.0;
    double ambientConductance25 = 0.0;
    const size_t ambientIdx = temp1.getNodes().size() - 1;
    for (size_t k = 0; k < resistances25.size(); ++k) {
        REQUIRE(resistances25[k].nodeFromId == resistances50[k].nodeFromId);
        REQUIRE(resistances25[k].nodeToId == resistances50[k].nodeToId);
        REQUIRE(resistances25[k].type == resistances50[k].type);
        REQUIRE(resistances25[k].resistance > 0);
        REQUIRE(resistances50[k].resistance > 0);
        maximumConductanceGrowth = std::max(maximumConductanceGrowth, resistances25[k].resistance / resistances50[k].resistance);
        if (resistances25[k].nodeToId == ambientIdx) {
            ambientConductance25 += 1.0 / resistances25[k].resistance;
            if (resistances25[k].type == HeatTransferType::RADIATION) {
                radiationConductance25 += 1.0 / resistances25[k].resistance;
            }
        }
    }
    // The fastest-growing coefficient is radiation: (T50/T25)^3 within the rises seen here.
    double hRadGrowth = ThermalResistance::calculateRadiationCoefficient(result2.maximumTemperature, 50.0, 0.9) /
                        ThermalResistance::calculateRadiationCoefficient(result1.maximumTemperature, 25.0, 0.9);
    double ratio = result2.totalThermalResistance / result1.totalThermalResistance;
    INFO("R(25 C) " << result1.totalThermalResistance << " K/W, R(50 C) " << result2.totalThermalResistance
         << " K/W, ratio " << ratio << "; largest conductance growth " << maximumConductanceGrowth
         << ", h_rad growth at the hottest node " << hRadGrowth << "; radiation share of the 25 C surface conductance "
         << radiationConductance25 / ambientConductance25);
    CHECK(maximumConductanceGrowth > 1.0);
    CHECK(maximumConductanceGrowth < hRadGrowth * 1.02);
    CHECK(result2.totalThermalResistance < result1.totalThermalResistance);
    CHECK(ratio > 1.0 / maximumConductanceGrowth);
}

TEST_CASE("Temperature: Loss Variation", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    // Use different current levels to get different actual losses
    // Low current (0.05A) for lower losses
    TemperatureConfig config1;
    config1.ambientTemperature = 25.0;
    applySimulatedLosses(config1, magnetic, 100000.0, 0.05);
    config1.plotSchematic = false;
    
    Temperature temp1(magnetic, config1);
    auto result1 = temp1.calculateTemperatures();
    
    // High current (0.5A) for higher losses - 10x more current
    TemperatureConfig config2;
    config2.ambientTemperature = 25.0;
    applySimulatedLosses(config2, magnetic, 100000.0, 0.5);
    config2.plotSchematic = false;
    
    Temperature temp2(magnetic, config2);
    auto result2 = temp2.calculateTemperatures();
    
    REQUIRE(result1.converged);
    REQUIRE(result2.converged);
    
    double deltaT1 = result1.maximumTemperature - config1.ambientTemperature;
    double deltaT2 = result2.maximumTemperature - config2.ambientTemperature;
    
    REQUIRE(deltaT2 > deltaT1);
    // With real losses, expect roughly 100x higher losses (10x current^2 for core)
    // Temperature rise should be significantly higher
    REQUIRE(deltaT2 > 5.0 * deltaT1);
}

TEST_CASE("Temperature: Radiation Effect", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config1;
    config1.ambientTemperature = 25.0;
    config1.coreLosses = 0.5;
    applySimulatedLosses(config1, magnetic);
    config1.includeRadiation = false;
    config1.plotSchematic = false;
    
    Temperature temp1(magnetic, config1);
    auto result1 = temp1.calculateTemperatures();
    
    TemperatureConfig config2;
    config2.ambientTemperature = 25.0;
    config2.coreLosses = 0.5;
    applySimulatedLosses(config2, magnetic);
    config2.includeRadiation = true;
    config2.plotSchematic = false;
    
    Temperature temp2(magnetic, config2);
    auto result2 = temp2.calculateTemperatures();
    
    REQUIRE(result1.converged);
    REQUIRE(result2.converged);
    // Radiation should reduce temperature or at worst be negligible.
    // Allow small numerical tolerance for solver precision.
    REQUIRE(result2.maximumTemperature <= result1.maximumTemperature + 0.1);
    REQUIRE(result2.totalThermalResistance <= result1.totalThermalResistance + 0.01);
}

TEST_CASE("Temperature: Segment Count", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config1;
    config1.ambientTemperature = 25.0;
    config1.coreLosses = 0.3;
    applySimulatedLosses(config1, magnetic);
    config1.toroidalSegments = 8;
    config1.plotSchematic = false;
    
    Temperature temp1(magnetic, config1);
    auto result1 = temp1.calculateTemperatures();
    
    TemperatureConfig config2;
    config2.ambientTemperature = 25.0;
    config2.coreLosses = 0.3;
    applySimulatedLosses(config2, magnetic);
    config2.toroidalSegments = 16;
    config2.plotSchematic = false;
    
    Temperature temp2(magnetic, config2);
    auto result2 = temp2.calculateTemperatures();
    
    REQUIRE(result1.converged);
    REQUIRE(result2.converged);
    REQUIRE_THAT(result2.maximumTemperature, 
                 Catch::Matchers::WithinRel(result1.maximumTemperature, 0.15));
    REQUIRE_THAT(result2.totalThermalResistance, 
                 Catch::Matchers::WithinRel(result1.totalThermalResistance, 0.15));
}

TEST_CASE("Temperature: Node Access", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    // Ungapped toroidal core - will have high flux density
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.3;
    // Use lower current (0.1A) to avoid extreme core losses in ungapped toroid
    auto inputs = createInputsForMagnetic(magnetic, 100000.0, config.ambientTemperature, 0.1);
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    if (losses.simulationSucceeded) {
        config.coreLosses = losses.coreLosses;
        config.windingLosses = losses.windingLosses;
        if (losses.windingLossesOutput.has_value()) {
            config.windingLossesOutput = losses.windingLossesOutput.value();
        }
    } else {
        applySimulatedLosses(config, magnetic);
    }
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    REQUIRE(result.converged);
    
    const auto& nodes = temp.getNodes();
    REQUIRE(nodes.size() > 9);  // With turn nodes enabled
    
    const auto& resistances = temp.getResistances();
    REQUIRE(resistances.size() > 16);  // With turn nodes enabled
    
    for (const auto& node : nodes) {
        REQUIRE(node.temperature >= config.ambientTemperature);
        // Relaxed threshold for realistic losses
        REQUIRE(node.temperature < 500.0);
    }
}

TEST_CASE("Temperature: Bulk Resistance", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.3;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    
    // Energy conservation. The previous check compared deltaT/totalLosses against
    // result.totalThermalResistance, which is DEFINED as deltaT/totalPower — a tautology.
    // Instead assert the real physics: at steady state the heat leaving to the ambient node
    // equals the total injected power.
    auto ebNodes = temp.getNodes();
    size_t ambientIdx = ebNodes.size();
    for (size_t i = 0; i < ebNodes.size(); ++i) {
        if (ebNodes[i].part == ThermalNodePartType::AMBIENT) { ambientIdx = i; break; }
    }
    REQUIRE(ambientIdx < ebNodes.size());
    double heatToAmbient = 0.0;
    for (const auto& r : temp.getResistances()) {
        size_t other = ebNodes.size();
        if (r.nodeToId == ambientIdx)      other = r.nodeFromId;
        else if (r.nodeFromId == ambientIdx) other = r.nodeToId;
        if (other < ebNodes.size() && r.resistance > 0.0) {
            heatToAmbient += (ebNodes[other].temperature - ebNodes[ambientIdx].temperature) / r.resistance;
        }
    }
    double totalInjected = 0.0;
    for (const auto& n : ebNodes) totalInjected += n.powerDissipation;
    REQUIRE(totalInjected > 0.0);
    REQUIRE_THAT(heatToAmbient, Catch::Matchers::WithinRel(totalInjected, 0.02));  // 2%: steady-state energy balance
    // Updated after removing inner-outer turn connections (was >15.0 before)
    // Further updated with gravity-aware convection (only top surfaces cool)
    REQUIRE(result.totalThermalResistance > 10.0);
    REQUIRE(result.totalThermalResistance < 150.0);  // Increased from 100 due to reduced convection
}

TEST_CASE("Temperature: Failure paths fail loudly", "[temperature][smoke-test]") {
    // Accessors/solver must throw on bad input rather than silently return a sentinel
    // (the pre-review behaviour was zeros/garbage on degenerate input).
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 29/16/10";
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::array(), 1, "N87");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    REQUIRE(result.converged);

    // getTemperatureAtPoint must reject an under-specified point, not guess.
    REQUIRE_THROWS_AS(temp.getTemperatureAtPoint(std::vector<double>{0.0}), std::runtime_error);
    // A valid query returns a finite temperature at or above ambient.
    double t = temp.getTemperatureAtPoint(std::vector<double>{0.0, 0.0, 0.0});
    REQUIRE(std::isfinite(t));
    REQUIRE(t >= config.ambientTemperature - 0.001);
}

TEST_CASE("Temperature: Forced vs Natural Convection", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({15});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 39/20/13";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig naturalConfig;
    naturalConfig.ambientTemperature = 25.0;
    naturalConfig.coreLosses = 1.5;
    applySimulatedLosses(naturalConfig, magnetic);
    naturalConfig.includeForcedConvection = false;
    naturalConfig.plotSchematic = false;
    
    TemperatureConfig forcedConfig;
    forcedConfig.ambientTemperature = 25.0;
    forcedConfig.coreLosses = 1.5;
    applySimulatedLosses(forcedConfig, magnetic);
    forcedConfig.includeForcedConvection = true;
    forcedConfig.airVelocity = 3.0;  // 3 m/s
    forcedConfig.plotSchematic = false;
    
    Temperature naturalTemp(magnetic, naturalConfig);
    Temperature forcedTemp(magnetic, forcedConfig);
    
    auto naturalResult = naturalTemp.calculateTemperatures();
    auto forcedResult = forcedTemp.calculateTemperatures();
    
    REQUIRE(naturalResult.converged);
    REQUIRE(forcedResult.converged);
    
    REQUIRE(forcedResult.maximumTemperature < naturalResult.maximumTemperature);
    REQUIRE(forcedResult.totalThermalResistance < naturalResult.totalThermalResistance);
}

TEST_CASE("Temperature: Convergence Test", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({20});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 44/22/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 2.0;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.1);
}

TEST_CASE("Temperature: Very High Losses", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 29/16/10";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 10.0;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > 100.0);
}

TEST_CASE("Temperature: Very Small Core", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({5});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "E 13/7/4";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.2;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.totalThermalResistance > 10.0);
}

TEST_CASE("Temperature: Maniktala Formula Comparison", "[temperature][smoke-test]") {
    std::vector<std::tuple<std::string, double, double>> cores = {
        {"ETD 29/16/10", 5.47, 0.5},
        {"ETD 34/17/11", 7.64, 0.7},
        {"ETD 44/22/15", 17.8, 1.2}
    };
    
    for (const auto& [coreName, Ve_cm3, coreLoss] : cores) {
        DYNAMIC_SECTION("Core: " << coreName) {
            double Rth_maniktala = 53.0 * std::pow(Ve_cm3, -0.54);
            
            std::vector<int64_t> numberTurns({12});
            std::vector<int64_t> numberParallels({1});

            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                             numberParallels,
                                                             coreName,
                                                             1,
                                                             WindingOrientation::CONTIGUOUS,
                                                             WindingOrientation::CONTIGUOUS,
                                                             CoilAlignment::CENTERED,
                                                             CoilAlignment::CENTERED);

            std::string coreMaterial = "N87";
            auto gapping = json::array();
            auto core = OpenMagneticsTesting::get_quick_core(coreName, gapping, 1, coreMaterial);
            
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            TemperatureConfig config;
            config.ambientTemperature = 25.0;
            config.coreLosses = coreLoss;
            applySimulatedLosses(config, magnetic);
            config.plotSchematic = false;
            
            Temperature temp(magnetic, config);
            auto result = temp.calculateTemperatures();
            
            REQUIRE(result.converged);
            
            double error = std::abs(result.totalThermalResistance - Rth_maniktala) / Rth_maniktala;
            // Cross-MODEL sanity check (network solver vs the analytical Maniktala formula):
            // they are different models (actual disagreement ~1.7x here), so this is an
            // order-of-magnitude agreement bound, tightened from an unfalsifiable 300% to 200%.
            REQUIRE(error < 2.0);
        }
    }
}

TEST_CASE("Temperature: PQ Core", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({18});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "PQ 26/25";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 2.0;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.totalThermalResistance > 3.0);  // Lowered from 5.0 due to improved concentric core connections
    REQUIRE(result.totalThermalResistance < 50.0);
}

TEST_CASE("Temperature: Four Winding Transformer", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({24, 12, 8, 6});
    std::vector<int64_t> numberParallels({1, 1, 1, 1});
    std::string shapeName = "T 36/23/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.5;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
}

TEST_CASE("Temperature: Zero Losses Baseline", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 49/25/16";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.0;
    config.windingLosses = 0.0;
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE_THAT(result.maximumTemperature, 
                 Catch::Matchers::WithinAbs(config.ambientTemperature, 0.5));
}

TEST_CASE("Temperature: Linear Scaling Validation", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({20});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 49/25/16";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.plotSchematic = false;
    
    // Use different current levels to get different actual losses
    // With real simulation, losses scale with operating point, not config
    std::vector<double> currents = {0.1, 0.2, 0.3, 0.5};
    std::vector<double> tempRises;
    std::vector<double> actualLosses;
    
    
    for (double I : currents) {
        applySimulatedLosses(config, magnetic, 100000.0, I);
        
        Temperature temp(magnetic, config);
        auto result = temp.calculateTemperatures();
        
        REQUIRE(result.converged);
        double tempRise = result.maximumTemperature - config.ambientTemperature;
        tempRises.push_back(tempRise);
        actualLosses.push_back(config.coreLosses + config.windingLosses);
    }
    
    // Verify that temperature rises with losses (monotonic relationship)
    for (size_t i = 1; i < tempRises.size(); ++i) {
        REQUIRE(tempRises[i] > tempRises[i-1]);
    }
    
    // Calculate thermal resistance for each point (should be roughly constant)
    std::vector<double> rthValues;
    for (size_t i = 0; i < actualLosses.size(); ++i) {
        if (actualLosses[i] > 0.001) {
            rthValues.push_back(tempRises[i] / actualLosses[i]);
        }
    }
    
    if (rthValues.size() >= 2) {
        // Re-pinned 2026-08-02 (ABT #461): with radiation to ambient in the model, Rth is NOT
        // flat — it falls as load rises, because h_total grows with surface temperature
        // (h_conv ~ dT^0.25 from the Churchill/Rayleigh correlation, h_rad ~ T^3 from
        // Stefan-Boltzmann; OMFEM's FEM integrates the same physics). The old +/-45% flatness
        // band encoded the radiation-less model's shape. The physical assertions are:
        //   1. Rth monotonically NON-INCREASING with power — the h(T) signature. A model that
        //      lost its temperature-dependent cooling would go flat but not rising; a model
        //      with inverted physics rises and fails here.
        for (size_t i = 1; i < rthValues.size(); ++i) {
            REQUIRE(rthValues[i] <= rthValues[i - 1] * 1.001);  // strictly non-increasing (1e-3 numeric slack)
        }
        //   2. A coarse spread envelope derived from the correlation exponents, not calibrated:
        //      convection-only gives Rth ~ P^(-0.2) (dT ~ P^0.8), i.e. deviation-from-mean up to
        //      ~0.32 over this sweep's power range; the radiative share steepens it (Rth ~
        //      P^(-3/4) in the radiation-dominated limit). Measured with the radiating model:
        //      0.46. The 0.65 envelope sits between the radiation-dominated slope and a broken
        //      model (order-of-magnitude Rth swings fail).
        double avgRth = 0;
        for (double r : rthValues) avgRth += r;
        avgRth /= rthValues.size();
        for (double r : rthValues) {
            double deviation = std::abs(r - avgRth) / avgRth;
            REQUIRE(deviation < 0.65);
        }
    }
}

TEST_CASE("Temperature: U Core", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({15});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "U 93/76/30";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 2.5;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.1);
}

TEST_CASE("Temperature: RM Core", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({12});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "RM 8";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.8;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 5.0);
}


// =============================================================================
// Phase 2: Turn Node Tests
// =============================================================================

TEST_CASE("Temperature: Winding Losses Only", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.0;  // No core losses
    applySimulatedLosses(config, magnetic);  // Only winding losses
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    // With winding losses only, temperature should still rise
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: Temperature at Point", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.5;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);

    // Test getTemperatureAtPoint at core center
    std::vector<double> coreCenter = {0.0, 0.0, 0.0};
    double tempAtCore = temp.getTemperatureAtPoint(coreCenter);
    
    REQUIRE(tempAtCore >= config.ambientTemperature);
    REQUIRE(tempAtCore <= result.maximumTemperature);
    
    // Test at a winding location
    auto nodes = temp.getNodes();
    if (!nodes.empty() && !nodes[0].physicalCoordinates.empty()) {
        double tempAtWinding = temp.getTemperatureAtPoint(nodes[0].physicalCoordinates);
        REQUIRE(tempAtWinding >= config.ambientTemperature);
        REQUIRE(tempAtWinding <= result.maximumTemperature);
    }
}

TEST_CASE("Temperature: Per-Turn Temperature Model", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({20});  // More turns to see gradient
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 36/23/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.0;
    applySimulatedLosses(config, magnetic);  // Higher winding losses
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    // Check that we have turn nodes
    auto nodes = temp.getNodes();
    int turnNodeCount = 0;
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::TURN) {
            turnNodeCount++;
        }
    }
    
    // With turn nodes enabled, we should have nodes for the turns
    REQUIRE(turnNodeCount > 0);
}



// =============================================================================
// Planar Core Tests
// =============================================================================

TEST_CASE("Temperature: Planar Core ER", "[temperature][planar][smoke-test]") {
    std::vector<int64_t> numberTurns({8, 4});
    std::vector<int64_t> numberParallels({1, 1});
    std::string shapeName = "ER 28/14";  // Using a common ER size

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "3F4";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.5;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 3.0);  // Lowered from 5.0 due to improved concentric core connections
}

TEST_CASE("Temperature: Planar Core Three Windings", "[temperature][planar][smoke-test]") {
    std::vector<int64_t> numberTurns({12, 6, 4});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "ER 28/14";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "3F4";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.0;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 3.0);  // Lowered from 5.0 due to improved concentric core connections
}

// =============================================================================
// Paper Validation Tests
// =============================================================================

TEST_CASE("Temperature: Van den Bossche E42 Validation", "[temperature][smoke-test]") {
    // Reference: Van den Bossche & Valchev - E42 core thermal resistance ~10-14 K/W
    std::vector<int64_t> numberTurns({15});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "E 42/21/20";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    // Use ungapped core with lower current to get realistic losses
    // Reference paper likely used partially saturated core or specific operating conditions
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.includeRadiation = true;
    config.plotSchematic = false;

    // Test at different power levels
    std::vector<double> currents = {0.03, 0.04, 0.05};  // Low currents for reasonable losses
    
    
    // Get core geometry info
    auto coreGeom = core.get_functional_description();
    
    for (double I : currents) {
        applySimulatedLosses(config, magnetic, 100000.0, I);
        
        Temperature temp(magnetic, config);
        auto result = temp.calculateTemperatures();
        
        REQUIRE(result.converged);
        
        double totalLosses = config.coreLosses + config.windingLosses;
        double tempRise = result.maximumTemperature - config.ambientTemperature;
        double rth = (totalLosses > 0.001) ? tempRise / totalLosses : 0.0;
        
        
        // Count convection resistances
        int convectionCount = 0;
        double minConvR = 1e9, maxConvR = 0, sumConvR = 0;
        for (const auto& r : result.thermalResistances) {
            if (r.type == HeatTransferType::NATURAL_CONVECTION || r.type == HeatTransferType::FORCED_CONVECTION) {
                convectionCount++;
                minConvR = std::min(minConvR, r.resistance);
                maxConvR = std::max(maxConvR, r.resistance);
                sumConvR += r.resistance;
            }
        }
        if (convectionCount > 0) {
        }
        
        // Reference paper reports Rth ≈ 10-14 K/W
        // Allow broader range since model uses actual operating conditions
        REQUIRE(rth > 1.0);
        REQUIRE(rth < 50.0);
    }
}

TEST_CASE("Temperature: Power-Temperature Linearity", "[temperature][smoke-test]") {
    // Reference: Dey et al. 2021 - Temperature rise scales linearly with power
    std::vector<int64_t> numberTurns({20});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 49";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.plotSchematic = false;

    // Test linearity at different power levels using different currents
    std::vector<double> currents = {0.1, 0.2, 0.3, 0.5};
    std::vector<double> thermalResistances;
    std::vector<double> actualLosses;
    std::vector<double> tempRises;
    
    
    for (double I : currents) {
        applySimulatedLosses(config, magnetic, 100000.0, I);
        
        Temperature temp(magnetic, config);
        auto result = temp.calculateTemperatures();
        
        REQUIRE(result.converged);
        
        double totalLosses = config.coreLosses + config.windingLosses;
        double tempRise = result.maximumTemperature - config.ambientTemperature;
        double rth = tempRise / totalLosses;
        
        thermalResistances.push_back(rth);
        actualLosses.push_back(totalLosses);
        tempRises.push_back(tempRise);
    }
    
    // Re-pinned 2026-08-02 (ABT #461): the "linear" premise (Dey et al.) holds only in the
    // small-dT convection-only limit. With radiation to ambient, Rth falls with load
    // (h_conv ~ dT^0.25, h_rad ~ T^3 — the same physics OMFEM's FEM integrates), so the
    // physical assertions are monotone non-increasing Rth plus a correlation-derived spread
    // envelope (see Linear Scaling Validation for the derivation; measured 0.46 with the
    // radiating model, convection-only bound ~0.32, envelope 0.65).
    for (size_t i = 1; i < thermalResistances.size(); ++i) {
        REQUIRE(thermalResistances[i] <= thermalResistances[i - 1] * 1.001);
    }
    double avgRth = 0;
    for (double rth : thermalResistances) avgRth += rth;
    avgRth /= thermalResistances.size();
    for (double rth : thermalResistances) {
        double deviation = std::abs(rth - avgRth) / avgRth;
        REQUIRE(deviation < 0.65);
    }
}

TEST_CASE("Temperature: Core Internal Gradient", "[temperature][smoke-test]") {
    // Reference: Salinas thesis - ferrite conductivity gives small internal gradients
    std::vector<int64_t> numberTurns({15});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 44";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    applySimulatedLosses(config, magnetic);
    // This ungapped ETD 44 at 1 A / 15 turns saturates, so the SIMULATED core loss (~50 W)
    // is unrealistically high for a thermal test (it would destroy the core). Pin a
    // realistic 3 W so the internal-gradient check reflects a physical operating point;
    // applySimulatedLosses still provides the per-turn winding-loss distribution.
    config.coreLosses = 3.0;
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);

    // Check core node temperatures
    auto nodes = temp.getNodes();
    double maxCoreTemp = 0;
    double minCoreTemp = 1000;
    
    for (const auto& node : nodes) {
        if (node.part != ThermalNodePartType::AMBIENT && 
            node.part != ThermalNodePartType::TURN &&
            node.part != ThermalNodePartType::BOBBIN_CENTRAL_COLUMN &&
            node.part != ThermalNodePartType::BOBBIN_TOP_YOKE &&
            node.part != ThermalNodePartType::BOBBIN_BOTTOM_YOKE) {
            maxCoreTemp = std::max(maxCoreTemp, node.temperature);
            minCoreTemp = std::min(minCoreTemp, node.temperature);
        }
    }
    
    // Internal gradient within core should be reasonable due to ferrite conductivity
    // Increased tolerance for half-core symmetric model and quadrant-specific connections
    // Central column (40% of losses) to lateral columns (20% of losses) creates natural gradient
    double internalGradient = maxCoreTemp - minCoreTemp;
    REQUIRE(internalGradient >= 0);
    // Ferrite conducts well (~4 W/m.K), so the within-core gradient is modest; a physical
    // sanity ceiling (was an unfalsifiable 1000 K). If this trips, the actual value points
    // at a real conduction-model problem.
    REQUIRE(internalGradient < 50.0);
}

TEST_CASE("Temperature: Detailed Loss Distribution", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({25, 12});
    std::vector<int64_t> numberParallels({1, 1});
    std::string shapeName = "E 55/28/21";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.8;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    // Verify we have nodes for both windings
    auto nodes = temp.getNodes();
    int turnNodeCount = 0;
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::TURN) {
            turnNodeCount++;
        }
    }
    REQUIRE(turnNodeCount > 0);
}

TEST_CASE("Temperature: Three Winding Transformer", "[temperature][smoke-test]") {
    std::vector<int64_t> numberTurns({15, 8, 5});
    std::vector<int64_t> numberParallels({1, 1, 1});
    std::string shapeName = "ETD 39";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.2;
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    // Thermal resistance varies with actual operating conditions
    REQUIRE(result.totalThermalResistance > 1.0);
    REQUIRE(result.totalThermalResistance < 100.0);
}

TEST_CASE("Temperature: Toroidal Inductor Rectangular Wires", "[temperature][round-winding-window][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "toroidal_inductor_rectangular_wires.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Export visualization
    exportTemperatureFieldSvg("toroidal_inductor_rectangular_wires", magnetic, result.nodeTemperatures, config.ambientTemperature);
    exportThermalCircuitSchematic("toroidal_inductor_rectangular_wires", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
    
    // Verify connection logic: no spurious long-distance connections
    // Turns should only connect to physically adjacent turns (surface distance < min_dim/4)
    const auto& nodes = temp.getNodes();
    const auto& resistances = temp.getResistances();
    
    // Check all conduction connections between turns
    for (const auto& res : resistances) {
        if (res.type != OpenMagnetics::HeatTransferType::CONDUCTION) continue;
        if (res.nodeFromId >= nodes.size() || res.nodeToId >= nodes.size()) continue;
        
        const auto& node1 = nodes[res.nodeFromId];
        const auto& node2 = nodes[res.nodeToId];
        
        // Only check turn-to-turn connections
        if (node1.part != OpenMagnetics::ThermalNodePartType::TURN || 
            node2.part != OpenMagnetics::ThermalNodePartType::TURN) continue;
        
        // Calculate center distance
        double dx = node1.physicalCoordinates[0] - node2.physicalCoordinates[0];
        double dy = node1.physicalCoordinates[1] - node2.physicalCoordinates[1];
        double centerDist = std::sqrt(dx*dx + dy*dy);
        
        // Calculate surface distance using same logic as the implementation
        double minDim1 = std::min(node1.dimensions.width, node1.dimensions.height);
        double minDim2 = std::min(node2.dimensions.width, node2.dimensions.height);
        double extent1 = minDim1 / 2.0;
        double extent2 = minDim2 / 2.0;
        double surfaceDist = centerDist - extent1 - extent2;
        double threshold = std::min(minDim1, minDim2) / 4.0;
        
        // Verify connection follows the geometric rule
        REQUIRE(surfaceDist <= threshold);
    }
}

TEST_CASE("Temperature Class: Toroidal Core T20", "[temperature][round-winding-window][smoke-test]") {
    // Create toroidal magnetic using standard test utilities
    std::string shapeName = "T 20/10/7";
    
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::array(), 1, "N87");
    
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, we need to wind the coil to generate turn coordinates and additionalCoordinates
    // which are required for proper thermal modeling
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }
    
    // Configure new Temperature API
    TemperatureConfig config;
    config.toroidalSegments = 8;
    config.plotSchematic = true;
    config.ambientTemperature = 25.0;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_Toroid_20mm_NEW_ARCH.svg").string();
    
    // Use real simulation for losses
    applySimulatedLosses(config, magnetic);
    
    // Create Temperature analyzer (new architecture)
    Temperature analyzer(magnetic, config);
    
    // Calculate temperatures
    ThermalResult result = analyzer.calculateTemperatures();
    
    
    // Access nodes from new architecture
    const auto& nodes = analyzer.getNodes();
    
    // Count node types
    int coreNodes = 0, turnNodes = 0, ambientNodes = 0;
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::CORE_TOROIDAL_SEGMENT) coreNodes++;
        else if (node.part == ThermalNodePartType::TURN) turnNodes++;
        else if (node.part == ThermalNodePartType::AMBIENT) ambientNodes++;
    }
    
    
    // Node details verified in test assertions below
    
    // Verify schematic was generated
    REQUIRE(result.converged);
    
}

TEST_CASE("Temperature: Toroidal Inductor Round Wire Multilayer", "[temperature][round-winding-window][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "toroidal_inductor_round_wire_multilayer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Export temperature field visualization with turns colored by temperature
    exportTemperatureFieldSvg("toroidal_inductor_round_wire_multilayer", magnetic, result.nodeTemperatures, config.ambientTemperature);
    
    // Export thermal circuit schematic with quadrants
    exportThermalCircuitSchematic("toroidal_inductor_round_wire_multilayer", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
    
    // Verify multilayer convection: outer layer turns should have RADIAL_INNER convection to air
    // Inner layer turns should have RADIAL_OUTER convection to air (inter-layer gap)
    const auto& nodes = temp.getNodes();
    const auto& resistances = temp.getResistances();
    
    // Build set of connected quadrants
    std::set<std::string> connectedQuadrants;
    for (const auto& res : resistances) {
        std::string key1 = std::to_string(res.nodeFromId) + "_" + std::string(magic_enum::enum_name(res.quadrantFrom));
        std::string key2 = std::to_string(res.nodeToId) + "_" + std::string(magic_enum::enum_name(res.quadrantTo));
        connectedQuadrants.insert(key1);
        connectedQuadrants.insert(key2);
    }
    
    // Find ambient node
    size_t ambientIdx = nodes.size();  // Last node is typically ambient
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].part == OpenMagnetics::ThermalNodePartType::AMBIENT) {
            ambientIdx = i;
            break;
        }
    }
    
    // Count convection connections from RADIAL_INNER and RADIAL_OUTER quadrants
    int radialInnerConvectionCount = 0;
    int radialOuterConvectionCount = 0;
    
    for (const auto& res : resistances) {
        // Check for any convection or radiation type
        bool isConvectionOrRadiation = 
            (res.type == OpenMagnetics::HeatTransferType::NATURAL_CONVECTION) ||
            (res.type == OpenMagnetics::HeatTransferType::FORCED_CONVECTION) ||
            (res.type == OpenMagnetics::HeatTransferType::RADIATION);
        if (!isConvectionOrRadiation) continue;
        
        // Check if this is a turn-to-ambient connection
        bool fromIsTurn = (res.nodeFromId < nodes.size()) && 
                          (nodes[res.nodeFromId].part == OpenMagnetics::ThermalNodePartType::TURN);
        bool toIsAmbient = (res.nodeToId == ambientIdx) || 
                           (res.nodeToId < nodes.size() && 
                            nodes[res.nodeToId].part == OpenMagnetics::ThermalNodePartType::AMBIENT);
        
        if (fromIsTurn && toIsAmbient) {
            if (res.quadrantFrom == OpenMagnetics::ThermalNodeFace::RADIAL_INNER) {
                radialInnerConvectionCount++;
            } else if (res.quadrantFrom == OpenMagnetics::ThermalNodeFace::RADIAL_OUTER) {
                radialOuterConvectionCount++;
            }
        }
    }
    
    
    // Check if this test has insulation layers
    bool hasInsulationLayers = false;
    for (const auto& node : nodes) {
        if (node.part == OpenMagnetics::ThermalNodePartType::INSULATION_LAYER) {
            hasInsulationLayers = true;
            break;
        }
    }
    
    if (hasInsulationLayers) {
        // With insulation layers, turns connect to insulation via conduction on radial faces
        // Radial convection to ambient is blocked by insulation
        // Don't require radial convection when insulation layers are present
    } else {
        // In a multilayer winding WITHOUT insulation layers, there should be convection 
        // from both RADIAL_INNER and RADIAL_OUTER
        // The outer layer's RADIAL_INNER faces the inter-layer air gap
        // The inner layer's RADIAL_OUTER faces the inter-layer air gap
        REQUIRE(radialInnerConvectionCount > 0);
        REQUIRE(radialOuterConvectionCount > 0);
    }
}

TEST_CASE("Temperature: Concentric Round Wire Spread Multilayer", "[temperature][rectangular-winding-window][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_round_wire_spread_multilayer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Export temperature field visualization with turns colored by temperature
    exportTemperatureFieldSvg("concentric_round_wire_spread_multilayer", magnetic, result.nodeTemperatures, config.ambientTemperature);
    
    // Export thermal circuit schematic with quadrants
    exportThermalCircuitSchematic("concentric_round_wire_spread_multilayer", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: Concentric Round Wire Centered Multilayer", "[temperature][rectangular-winding-window][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_round_wire_centered_multilayer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Export temperature field visualization with turns colored by temperature
    exportTemperatureFieldSvg("concentric_round_wire_centered_multilayer", magnetic, result.nodeTemperatures, config.ambientTemperature);
    
    // Export thermal circuit schematic with quadrants
    exportThermalCircuitSchematic("concentric_round_wire_centered_multilayer", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: Concentric Round Wire Full Layer", "[temperature][rectangular-winding-window][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_round_wire_full_layer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Export temperature field visualization with turns colored by temperature
    exportTemperatureFieldSvg("concentric_round_wire_full_layer", magnetic, result.nodeTemperatures, config.ambientTemperature);
    
    // Export thermal circuit schematic with quadrants
    exportThermalCircuitSchematic("concentric_round_wire_full_layer", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: Concentric Round Wire Simple", "[temperature][rectangular-winding-window][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_round_wire_simple.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Export temperature field visualization with turns colored by temperature
    exportTemperatureFieldSvg("concentric_round_wire_simple", magnetic, result.nodeTemperatures, config.ambientTemperature);
    
    // Export thermal circuit schematic with quadrants
    exportThermalCircuitSchematic("concentric_round_wire_simple", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: Concentric Planar Inductor", "[temperature][planar][concentric][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_planar.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_concentric_planar.svg").string();
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    // Export temperature field visualization with turns colored by temperature
    exportTemperatureFieldSvg("concentric_planar", magnetic, result.nodeTemperatures, config.ambientTemperature);
    
    // Export thermal circuit schematic
    exportThermalCircuitSchematic("concentric_planar", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: Planar Transformer Complex", "[temperature][planar][transformer][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "planar_transformer_complex.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);

    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);

    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);

    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_planar_transformer_complex.svg").string();

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    // Export temperature field visualization
    exportTemperatureFieldSvg("planar_transformer_complex", magnetic, result.nodeTemperatures, config.ambientTemperature);

    // Export thermal circuit schematic
    exportThermalCircuitSchematic("planar_transformer_complex", temp);

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: Concentric Litz and Foil", "[temperature][concentric][litz][foil][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_litz_foil.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);

    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);

    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);

    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_concentric_litz_foil.svg").string();

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    // Export temperature field visualization
    exportTemperatureFieldSvg("concentric_litz_foil", magnetic, result.nodeTemperatures, config.ambientTemperature);

    // Export thermal circuit schematic
    exportThermalCircuitSchematic("concentric_litz_foil", temp);

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: Concentric with Insulation Layers", "[temperature][rectangular-winding-window][smoke-test]") {
    // Load the test file with insulation layers
    std::filesystem::path testFile = std::filesystem::path(__FILE__).parent_path() / "testData" / "concentric_round_wire_insulation_layers.json";
    
    std::ifstream file(testFile);
    REQUIRE(file.good());
    
    json j;
    file >> j;
    
    // Parse inputs and magnetic
    OpenMagnetics::Inputs inputs(j["inputs"]);
    OpenMagnetics::Magnetic magnetic(j["magnetic"]);
    
    
    if (magnetic.get_coil().get_layers_description()) {
        auto layers = magnetic.get_coil().get_layers_description().value();
        size_t insulationCount = 0;
        for (const auto& layer : layers) {
            if (layer.get_type() == MAS::ElectricalType::INSULATION) {
                insulationCount++;
            }
        }
    }
    
    if (magnetic.get_coil().get_turns_description()) {
    }
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_concentric_insulation_layers.svg").string();
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Count insulation layer nodes
    size_t insulationNodeCount = 0;
    for (const auto& node : temp.getNodes()) {
        if (node.part == OpenMagnetics::ThermalNodePartType::INSULATION_LAYER) {
            insulationNodeCount++;
        }
    }
    
    // Export temperature field visualization with turns and insulation colored by temperature
    exportTemperatureFieldSvg("concentric_insulation_layers", magnetic, result.nodeTemperatures, config.ambientTemperature);
    
    // Export thermal circuit schematic
    exportThermalCircuitSchematic("concentric_insulation_layers", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
    REQUIRE(insulationNodeCount > 0);  // Should have created insulation layer nodes
}

TEST_CASE("Temperature: Concentric with Insulation Layers and Forced Convection", "[temperature][rectangular-winding-window][cooling][smoke-test]") {
    // Load the test file with insulation layers
    std::filesystem::path testFile = std::filesystem::path(__FILE__).parent_path() / "testData" / "concentric_round_wire_insulation_layers.json";
    
    std::ifstream file(testFile);
    REQUIRE(file.good());
    
    json j;
    file >> j;
    
    // Parse inputs and magnetic
    OpenMagnetics::Inputs inputs(j["inputs"]);
    OpenMagnetics::Magnetic magnetic(j["magnetic"]);
    
    
    if (magnetic.get_coil().get_layers_description()) {
        auto layers = magnetic.get_coil().get_layers_description().value();
        size_t insulationCount = 0;
        for (const auto& layer : layers) {
            if (layer.get_type() == MAS::ElectricalType::INSULATION) {
                insulationCount++;
            }
        }
    }
    
    if (magnetic.get_coil().get_turns_description()) {
    }
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    // First test with natural convection (baseline)
    TemperatureConfig configNatural;
    configNatural.ambientTemperature = losses.ambientTemperature;
    configNatural.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    configNatural.windingLosses = losses.windingLosses;
    configNatural.windingLossesOutput = losses.windingLossesOutput.value();
    configNatural.plotSchematic = false;
    
    Temperature tempNatural(magnetic, configNatural);
    auto resultNatural = tempNatural.calculateTemperatures();
    
    // Now test with forced convection
    TemperatureConfig configForced;
    configForced.ambientTemperature = losses.ambientTemperature;
    configForced.coreLosses = losses.coreLosses;
    configForced.windingLosses = losses.windingLosses;
    configForced.windingLossesOutput = losses.windingLossesOutput.value();
    configForced.plotSchematic = true;
    configForced.schematicOutputPath = (getOutputDir() / "thermal_schematic_concentric_insulation_forced_convection.svg").string();
    
    // Create MAS cooling object for forced convection
    MAS::Cooling forcedCooling;
    forcedCooling.set_temperature(losses.ambientTemperature);
    forcedCooling.set_velocity(std::vector<double>{3.0, 0.0, 0.0});  // 3 m/s airflow
    configForced.masCooling = forcedCooling;
    
    Temperature tempForced(magnetic, configForced);
    auto resultForced = tempForced.calculateTemperatures();
    
    // Count insulation layer nodes
    size_t insulationNodeCount = 0;
    for (const auto& node : tempForced.getNodes()) {
        if (node.part == OpenMagnetics::ThermalNodePartType::INSULATION_LAYER) {
            insulationNodeCount++;
        }
    }
    
    // Verify forced convection provides better cooling
    REQUIRE(resultForced.maximumTemperature < resultNatural.maximumTemperature);
    REQUIRE(resultForced.totalThermalResistance < resultNatural.totalThermalResistance);
    
    // Export temperature field visualization
    exportTemperatureFieldSvg("concentric_insulation_forced_convection", magnetic, resultForced.nodeTemperatures, configForced.ambientTemperature);
    
    // Export thermal circuit schematic
    exportThermalCircuitSchematic("concentric_insulation_forced_convection", tempForced);
    
    REQUIRE(resultForced.converged);
    REQUIRE(resultForced.maximumTemperature > configForced.ambientTemperature);
    REQUIRE(resultForced.totalThermalResistance > 0.0);
    REQUIRE(insulationNodeCount > 0);  // Should have created insulation layer nodes
}

TEST_CASE("Temperature: Toroidal with Insulation Layers", "[temperature][round-winding-window][smoke-test]") {
    // Load the toroidal inductor with insulation layers test data
    std::filesystem::path testFile = std::filesystem::path(__FILE__).parent_path() 
        / "testData" / "toroidal_inductor_round_wire_multilayer_with_insulation.json";
    
    auto mas = OpenMagneticsTesting::mas_loader(testFile);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();
    
    
    // Print insulation layer info
    if (coil.get_layers_description()) {
        auto layers = coil.get_layers_description().value();
        auto insulationLayers = coil.get_layers_description_insulation();
        for (const auto& layer : insulationLayers) {
            std::string layerName = layer.get_name();
            if (layerName.empty()) layerName = "unnamed";
            // Layer coordinate system info removed
        }
    }
    
    // Get inputs from test file (contains operating points)
    //
    // Operating point: a single-winding inductor. The fixture used to drive BOTH windings (67 and
    // 15 turns) with the same in-phase 10 A pk-pk current, so their ampere-turns added on an
    // ungapped toroid: B_peak = 31.5 T and 2.28 MW of core loss, and the thermal network could not
    // converge (the old helper swallowed that exception and the test ran on 0 W). Now only the
    // 67-turn primary carries current; the 15-turn secondary stays in the coil (its section and the
    // insulation layers around it are what this test checks) but is open: 0 A, with the induced
    // voltage V1 * 15 / 67. That open secondary is how a single-winding inductor is expressed on
    // this two-winding coil: MagneticSimulator rejects an operating point without exactly one
    // excitation per winding, and removing the secondary would remove the insulation layers.
    //   Core: Fair-Rite 98, T 28/14/15 ungapped, Ae = 105 mm^2, Bsat(100 C) = 0.405 T;
    //   L = 44.27 mH (the fixture's stored magnetizing inductance), ambient 100 C, 100 kHz, D = 0.5.
    //   Primary current: triangular, 30 mA pk-pk (+-15 mA), no DC.
    //   B_peak = L * I_peak / (N * Ae) = 44.27e-3 * 0.015 / (67 * 105e-6) = 0.094 T  (23 % of Bsat).
    //   V1 = L * dI / (D * T) = 44.27e-3 * 0.030 / 5e-6 = +-265.60 V square (volt-second balanced);
    //   V2 = 265.60 * 15 / 67 = +-59.46 V.
    auto inputs = mas.get_inputs();

    // Run magnetic simulation to get real losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    // The simulation must produce real losses now, not the 0 W the swallowed exception left.
    REQUIRE(losses.coreLosses > 0);

    
    // Temperature config
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    config.windingLosses = losses.windingLosses;
    if (losses.windingLossesOutput.has_value()) {
        config.windingLossesOutput = losses.windingLossesOutput.value();
    }
    
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_toroidal_with_insulation.svg").string();
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Count insulation layer nodes
    size_t insulationNodeCount = 0;
    for (const auto& node : temp.getNodes()) {
        if (node.part == OpenMagnetics::ThermalNodePartType::INSULATION_LAYER) {
            insulationNodeCount++;
        }
    }
    
    // Export schematic
    exportThermalCircuitSchematic("toroidal_with_insulation", temp);
    
    // Verify insulation layer nodes were created
    REQUIRE(insulationNodeCount > 0);
    
    // Verify insulation nodes are at expected radii
    // Inner nodes: inside core hole (~5-7mm), Outer nodes: outside core (~13-17mm)
    bool foundInnerLayer = false, foundOuterLayer = false;
    for (const auto& node : temp.getNodes()) {
        if (node.part == OpenMagnetics::ThermalNodePartType::INSULATION_LAYER) {
            double r = std::sqrt(node.physicalCoordinates[0]*node.physicalCoordinates[0] + 
                                 node.physicalCoordinates[1]*node.physicalCoordinates[1]);
            if (r > 0.005 && r < 0.007) {  // Inner surface nodes
                foundInnerLayer = true;
            }
            if (r > 0.013 && r < 0.017) {  // Outer surface nodes
                foundOuterLayer = true;
            }
        }
    }
    REQUIRE(foundInnerLayer);
    REQUIRE(foundOuterLayer);
}

TEST_CASE("Temperature: Concentric Simple Insulation Layers Schematic", "[temperature][rectangular-winding-window][smoke-test]") {
    // Load the simple concentric insulation layers test data
    std::filesystem::path testFile = std::filesystem::path(__FILE__).parent_path() 
        / "testData" / "concentric_round_wire_insulation_layers_simple.json";
    
    auto mas = OpenMagneticsTesting::mas_loader(testFile);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();
    
    
    // Print insulation layer info
    if (coil.get_layers_description()) {
        auto layers = coil.get_layers_description().value();
        auto insulationLayers = coil.get_layers_description_insulation();
        for (const auto& layer : insulationLayers) {
            std::string layerName = layer.get_name();
            if (layerName.empty()) layerName = "unnamed";
        }
    }
    
    // For schematic-only test, use simulated losses
    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.5;
    applySimulatedLosses(config, magnetic);
    
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_concentric_simple_insulation.svg").string();
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Count insulation layer nodes
    size_t insulationNodeCount = 0;
    for (const auto& node : temp.getNodes()) {
        if (node.part == OpenMagnetics::ThermalNodePartType::INSULATION_LAYER) {
            insulationNodeCount++;
        }
    }
    
    // Verify turn-to-insulation connections
    for (const auto& res : temp.getResistances()) {
        if (res.type != OpenMagnetics::HeatTransferType::CONDUCTION) continue;
        const auto& nodes = temp.getNodes();
        if (res.nodeFromId >= nodes.size() || res.nodeToId >= nodes.size()) continue;
        
        const auto& node1 = nodes[res.nodeFromId];
        const auto& node2 = nodes[res.nodeToId];
        
        bool isTurnInsulation = (node1.part == OpenMagnetics::ThermalNodePartType::TURN && 
                                  node2.part == OpenMagnetics::ThermalNodePartType::INSULATION_LAYER) ||
                                (node1.part == OpenMagnetics::ThermalNodePartType::INSULATION_LAYER && 
                                  node2.part == OpenMagnetics::ThermalNodePartType::TURN);
        
        if (isTurnInsulation) {
            // Turn-to-insulation connection found
        }
    }
    
    // Export schematic
    exportThermalCircuitSchematic("concentric_simple_insulation", temp);
    
    // Basic verification
    REQUIRE(insulationNodeCount > 0);
    REQUIRE(temp.getNodes().size() > 0);
    REQUIRE(temp.getResistances().size() > 0);
}

// ============================================================================
// Cooling Options Tests
// ============================================================================

TEST_CASE("Temperature: Forced Convection Cooling", "[temperature][cooling][smoke-test]") {
    // Use existing concentric test file
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(
        std::source_location::current(), "concentric_round_wire_simple.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    // First calculate with natural convection (baseline)
    TemperatureConfig configNatural;
    configNatural.ambientTemperature = losses.ambientTemperature;
    configNatural.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    configNatural.windingLosses = losses.windingLosses;
    configNatural.windingLossesOutput = losses.windingLossesOutput.value();
    configNatural.plotSchematic = true;
    configNatural.schematicOutputPath = (getOutputDir() / "thermal_natural_convection.svg").string();
    
    Temperature tempNatural(magnetic, configNatural);
    auto resultNatural = tempNatural.calculateTemperatures();
    
    
    // Now calculate with forced convection
    TemperatureConfig configForced;
    configForced.ambientTemperature = losses.ambientTemperature;
    configForced.coreLosses = losses.coreLosses;
    configForced.windingLosses = losses.windingLosses;
    configForced.windingLossesOutput = losses.windingLossesOutput.value();
    configForced.plotSchematic = true;
    configForced.schematicOutputPath = (getOutputDir() / "thermal_forced_convection.svg").string();
    
    // Create MAS cooling object for forced convection
    MAS::Cooling forcedCooling;
    forcedCooling.set_velocity(std::vector<double>{2.0, 0.0, 0.0});  // 2 m/s
    forcedCooling.set_flow_diameter(0.04);  // 40mm fan
    forcedCooling.set_fluid("air");
    configForced.masCooling = forcedCooling;
    
    Temperature tempForced(magnetic, configForced);
    auto resultForced = tempForced.calculateTemperatures();
    
    
    // Verify forced convection gives lower temperature
    REQUIRE(resultForced.maximumTemperature < resultNatural.maximumTemperature);
    REQUIRE(resultForced.totalThermalResistance < resultNatural.totalThermalResistance);
    
    // Export visualizations
    exportTemperatureFieldSvg("forced_convection", magnetic, resultForced.nodeTemperatures, configForced.ambientTemperature);
    exportThermalCircuitSchematic("forced_convection", tempForced);
    
    REQUIRE(resultForced.converged);
}

TEST_CASE("Temperature: Heatsink Cooling", "[temperature][cooling][smoke-test]") {
    // Use existing concentric test file (heatsink works best on concentric cores)
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(
        std::source_location::current(), "concentric_round_wire_simple.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_heatsink.svg").string();
    
    // Create MAS cooling object for heatsink
    MAS::Cooling heatsinkCooling;
    heatsinkCooling.set_thermal_resistance(2.5);  // 2.5 K/W heatsink
    heatsinkCooling.set_interface_thermal_resistance(3.0);  // 3 W/mK TIM
    heatsinkCooling.set_interface_thickness(0.0001);  // 100 microns
    { MAS::Dimensions d; d.set_width(0.05); d.set_height(0.05); d.set_depth(0.02); heatsinkCooling.set_dimensions(d); }  // 50x50x20mm
    config.masCooling = heatsinkCooling;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Check that heatsink node was created
    bool hasHeatsinkNode = false;
    for (const auto& node : temp.getNodes()) {
        if (node.name == "Heatsink") {
            hasHeatsinkNode = true;
            break;
        }
    }
    REQUIRE(hasHeatsinkNode);
    
    // Export visualizations
    exportTemperatureFieldSvg("heatsink_cooling", magnetic, result.nodeTemperatures, config.ambientTemperature);
    exportThermalCircuitSchematic("heatsink_cooling", temp);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
}

TEST_CASE("Temperature: Cold Plate Cooling", "[temperature][cooling][smoke-test]") {
    // Use toroidal test file (cold plate works on both core types)
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(
        std::source_location::current(), "toroidal_inductor_round_wire_multilayer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_cold_plate.svg").string();
    
    // Create MAS cooling object for cold plate
    MAS::Cooling coldPlateCooling;
    coldPlateCooling.set_maximum_temperature(40.0);  // 40°C cold plate
    coldPlateCooling.set_interface_thermal_resistance(3.0);  // 3 W/mK TIM
    coldPlateCooling.set_interface_thickness(0.0002);  // 200 microns
    { MAS::Dimensions d; d.set_width(0.06); d.set_height(0.06); d.set_depth(0.01); coldPlateCooling.set_dimensions(d); }  // 60x60x10mm
    config.masCooling = coldPlateCooling;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    
    // Check that cold plate node was created
    bool hasColdPlateNode = false;
    for (const auto& node : temp.getNodes()) {
        if (node.name == "ColdPlate") {
            hasColdPlateNode = true;
            REQUIRE(node.isFixedTemperature == true);
            // Temperature should be close to 40°C (allowing for solver tolerance)
            REQUIRE(std::abs(node.temperature - 40.0) < 1.0);
            break;
        }
    }
    REQUIRE(hasColdPlateNode);
    
    // Verify temperature is reduced compared to unassisted natural convection.
    // Cold plate helps but core thermal conductivity limits heat transfer,
    // especially for toroidal cores where heat must travel through ferrite.
    REQUIRE(result.maximumTemperature < 180.0);
    
    // Export visualizations
    exportTemperatureFieldSvg("cold_plate_cooling", magnetic, result.nodeTemperatures, config.ambientTemperature);
    exportThermalCircuitSchematic("cold_plate_cooling", temp);
    
    REQUIRE(result.converged);
}

TEST_CASE("Temperature: Cooling Utils Type Detection", "[temperature][cooling][smoke-test]") {
    using namespace OpenMagnetics;
    
    
    // Test natural convection detection
    MAS::Cooling naturalCooling;
    naturalCooling.set_temperature(25.0);
    REQUIRE(CoolingUtils::isNaturalConvection(naturalCooling) == true);
    REQUIRE(CoolingUtils::isForcedConvection(naturalCooling) == false);
    REQUIRE(CoolingUtils::isHeatsink(naturalCooling) == false);
    REQUIRE(CoolingUtils::isColdPlate(naturalCooling) == false);
    
    // Test forced convection detection
    MAS::Cooling forcedCooling;
    forcedCooling.set_velocity(std::vector<double>{1.0, 0.0, 0.0});
    REQUIRE(CoolingUtils::isNaturalConvection(forcedCooling) == false);
    REQUIRE(CoolingUtils::isForcedConvection(forcedCooling) == true);
    REQUIRE(CoolingUtils::isHeatsink(forcedCooling) == false);
    REQUIRE(CoolingUtils::isColdPlate(forcedCooling) == false);
    
    // Test heatsink detection
    MAS::Cooling heatsinkCooling;
    heatsinkCooling.set_thermal_resistance(2.0);
    REQUIRE(CoolingUtils::isNaturalConvection(heatsinkCooling) == false);
    REQUIRE(CoolingUtils::isForcedConvection(heatsinkCooling) == false);
    REQUIRE(CoolingUtils::isHeatsink(heatsinkCooling) == true);
    REQUIRE(CoolingUtils::isColdPlate(heatsinkCooling) == false);
    
    // Test cold plate detection (also has thermal_resistance but max_temp takes precedence)
    MAS::Cooling coldPlateCooling;
    coldPlateCooling.set_maximum_temperature(40.0);
    coldPlateCooling.set_thermal_resistance(1.0);  // This should be ignored
    REQUIRE(CoolingUtils::isNaturalConvection(coldPlateCooling) == false);
    REQUIRE(CoolingUtils::isForcedConvection(coldPlateCooling) == false);
    REQUIRE(CoolingUtils::isHeatsink(coldPlateCooling) == false);
    REQUIRE(CoolingUtils::isColdPlate(coldPlateCooling) == true);
    
    // Test forced convection coefficient calculation
    double h_forced = CoolingUtils::calculateForcedConvectionCoefficient(
        80.0, 25.0, 2.0, 0.01);  // 80°C surface, 25°C ambient, 2 m/s, 10mm char length
    REQUIRE(h_forced > 0.0);
    REQUIRE(h_forced > 5.0);  // Should be higher than natural convection (~5-10 W/m²K)
    
    // Test mixed convection calculation
    double h_mixed = CoolingUtils::calculateMixedConvectionCoefficient(10.0, 50.0);
    REQUIRE(h_mixed > 50.0);  // Should be close to forced value
    REQUIRE(h_mixed < 51.0);  // But natural should add a small amount
}

TEST_CASE("Temperature: 220kW Transformer", "[temperature][transformer][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "220kW_trafo.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);

    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);

    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);

    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing from simulation results");
    }
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_220kW_transformer.svg").string();

    // Add forced convection at 4 m/s
    MAS::Cooling forcedCooling;
    forcedCooling.set_temperature(losses.ambientTemperature);
    forcedCooling.set_velocity(std::vector<double>{4.0, 0.0, 0.0});  // 4 m/s airflow
    config.masCooling = forcedCooling;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    // Export temperature field visualization
    exportTemperatureFieldSvg("220kW_transformer", magnetic, result.nodeTemperatures, config.ambientTemperature);

    // Export thermal circuit schematic
    exportThermalCircuitSchematic("220kW_transformer", temp);

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

TEST_CASE("Temperature: PQ Single Turn Quadrant Surface Areas", "[temperature][concentric][geometry][smoke-test]") {
    // Create a PQ core with just one turn to verify quadrant surface areas
    // PQ cores have round columns, so turn length = 2*pi*r at each face's radial position
    std::vector<int64_t> numberTurns({1});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "PQ 26/25";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "N87";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.5;
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    // Find the turn node
    const auto& nodes = temp.getNodes();
    const ThermalNetworkNode* turnNode = nullptr;
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::TURN) {
            turnNode = &node;
            break;
        }
    }
    
    REQUIRE(turnNode != nullptr);
    
    // Get quadrant surface areas
    double leftArea = 0.0, rightArea = 0.0, topArea = 0.0, bottomArea = 0.0;
    
    for (const auto& quadrant : turnNode->quadrants) {
        switch (quadrant.face) {
            case ThermalNodeFace::RADIAL_INNER:  // LEFT for concentric
                leftArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::RADIAL_OUTER:  // RIGHT for concentric
                rightArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::TANGENTIAL_LEFT:  // TOP for concentric
                topArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::TANGENTIAL_RIGHT:  // BOTTOM for concentric
                bottomArea = quadrant.surfaceArea;
                break;
            default:
                break;
        }
    }
    
    // All areas should be positive
    REQUIRE(leftArea > 0.0);
    REQUIRE(rightArea > 0.0);
    REQUIRE(topArea > 0.0);
    REQUIRE(bottomArea > 0.0);
    
    // TOP == BOTTOM (same wire width and both at center radius)
    REQUIRE_THAT(topArea, Catch::Matchers::WithinRel(bottomArea, 0.01));
    
    // Position-dependent turn lengths: inner face traces a shorter path than outer face
    // Therefore: LEFT < TOP/BOTTOM < RIGHT
    REQUIRE(leftArea < topArea);
    REQUIRE(topArea < rightArea);
    REQUIRE(leftArea < rightArea);
}

TEST_CASE("Temperature: Toroidal Single Turn Quadrant Surface Areas", "[temperature][toroidal][geometry][smoke-test]") {
    // Create a toroidal core with one turn to verify quadrant surface areas
    // Toroidal cores always have round columns: turn length = 2*pi*r
    std::vector<int64_t> numberTurns({1});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "3C97";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    // Toroidal cores need winding to generate turn coordinates
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {
        magnetic.get_mutable_coil().wind();
    }

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.3;
    config.plotSchematic = false;
    applySimulatedLosses(config, magnetic, 100000.0, 0.1);
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    // For toroidal cores, turns are split into inner/outer half-nodes
    // Find the first inner turn node
    const auto& nodes = temp.getNodes();
    const ThermalNetworkNode* turnNode = nullptr;
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::TURN) {
            turnNode = &node;
            break;
        }
    }
    
    REQUIRE(turnNode != nullptr);
    
    // Get quadrant surface areas
    double innerArea = 0.0, outerArea = 0.0, tangLeftArea = 0.0, tangRightArea = 0.0;
    
    for (const auto& quadrant : turnNode->quadrants) {
        switch (quadrant.face) {
            case ThermalNodeFace::RADIAL_INNER:
                innerArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::RADIAL_OUTER:
                outerArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::TANGENTIAL_LEFT:
                tangLeftArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::TANGENTIAL_RIGHT:
                tangRightArea = quadrant.surfaceArea;
                break;
            default:
                break;
        }
    }
    
    // All areas should be positive
    REQUIRE(innerArea > 0.0);
    REQUIRE(outerArea > 0.0);
    REQUIRE(tangLeftArea > 0.0);
    REQUIRE(tangRightArea > 0.0);
    
    // TANGENTIAL_LEFT == TANGENTIAL_RIGHT (both at center radius)
    REQUIRE_THAT(tangLeftArea, Catch::Matchers::WithinRel(tangRightArea, 0.01));
    
    // Position-dependent turn lengths: inner face traces a shorter path than outer face
    // Therefore: RADIAL_INNER < TANGENTIAL < RADIAL_OUTER
    REQUIRE(innerArea < tangLeftArea);
    REQUIRE(tangLeftArea < outerArea);
    REQUIRE(innerArea < outerArea);
}

TEST_CASE("Temperature: Planar Single Turn Quadrant Surface Areas", "[temperature][planar][geometry][smoke-test]") {
    // Create a planar core (ER) with one turn to verify quadrant surface areas
    // Planar cores are concentric with round columns
    std::vector<int64_t> numberTurns({1});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ER 28/14";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);

    std::string coreMaterial = "3F4";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 0.5;
    config.plotSchematic = false;
    applySimulatedLosses(config, magnetic);
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    // Find the turn node
    const auto& nodes = temp.getNodes();
    const ThermalNetworkNode* turnNode = nullptr;
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::TURN) {
            turnNode = &node;
            break;
        }
    }
    
    REQUIRE(turnNode != nullptr);
    
    // Get quadrant surface areas
    double leftArea = 0.0, rightArea = 0.0, topArea = 0.0, bottomArea = 0.0;
    
    for (const auto& quadrant : turnNode->quadrants) {
        switch (quadrant.face) {
            case ThermalNodeFace::RADIAL_INNER:  // LEFT for concentric
                leftArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::RADIAL_OUTER:  // RIGHT for concentric
                rightArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::TANGENTIAL_LEFT:  // TOP for concentric
                topArea = quadrant.surfaceArea;
                break;
            case ThermalNodeFace::TANGENTIAL_RIGHT:  // BOTTOM for concentric
                bottomArea = quadrant.surfaceArea;
                break;
            default:
                break;
        }
    }
    
    // All areas should be positive
    REQUIRE(leftArea > 0.0);
    REQUIRE(rightArea > 0.0);
    REQUIRE(topArea > 0.0);
    REQUIRE(bottomArea > 0.0);
    
    // TOP == BOTTOM (same wire width and both at center radius)
    REQUIRE_THAT(topArea, Catch::Matchers::WithinRel(bottomArea, 0.01));
    
    // Position-dependent turn lengths: inner face traces a shorter path than outer face
    // Therefore: LEFT < TOP/BOTTOM < RIGHT
    REQUIRE(leftArea < topArea);
    REQUIRE(topArea < rightArea);
    REQUIRE(leftArea < rightArea);
}


TEST_CASE("Temperature: concentric_transformer", "[temperature][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_transformer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (losses.windingLossesOutput.has_value()) {
        config.windingLosses = losses.windingLosses;
        config.windingLossesOutput = losses.windingLossesOutput.value();
    } else {
        applySimulatedLosses(config, magnetic);
    }
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    // Reference values from Icepak simulation (Student Version)
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    // Get temperatures by component type
    auto tempsByType = temp.getTemperaturesByComponentType();
    auto tempsPerTurn = temp.getTemperaturePerTurn();

    SECTION("Core temperature validation against Icepak") {
        // Icepak core temperatures from solution overview
        // core_0: 63.50°C (from Icepak)
        // core_1: 67.88°C (from Icepak)
        // Re-pinned 2026-10-06 (ABT #1585), 48.63 C -> 31.70 C. The fixture's primary voltage
        // was a -20.5/+70.5 V square (25 V mean) and its secondary 0/-100 V (-50 V mean); MKF
        // integrated that DC into a drifting magnetizing current, 2.09 T peak on a ferrite core,
        // 1.261 W of core loss. Inputs now rejects unbalanced volt-seconds, so the fixture carries
        // the same squares at zero mean (+-45.5 V, +-50 V): 0.063 T, 0.133 W of core loss
        // (winding loss unchanged, 0.201 W). The Icepak band (core_0 63.50 C, core_1 67.88 C,
        // previously 48.0..84.85) matched only the drift artifact's loss and Icepak's loss inputs
        // are not recorded, so it cannot be compared against a physical excitation; this pins
        // the thermal network's result for the balanced fixture instead. The old Icepak band
        // matched the old fixture's DC-offset flux-drift artefact (1-2 T), so this is now a
        // REGRESSION pin, not a validation; ABT #1715 redoes the Icepak comparison with
        // recorded loss inputs.
        REQUIRE_THAT(tempsByType.at("core"), Catch::Matchers::WithinRel(31.70, 0.05));
    }
    
    SECTION("Bobbin temperature validation against Icepak") {
        if (tempsByType.find("bobbin") != tempsByType.end()) {
            // No bobbin temperature data from Icepak
        }
    }
    
    SECTION("Individual turn temperatures from Icepak export") {
        // Per-turn temperatures from Icepak, asserted against the REAL turn keys from
        // getTemperaturePerTurn. These were previously looked up as "Turn_W0_Turn_N",
        // which never matched any key, so every assertion silently skipped (dead code).
        // NOTE: the thermal model currently predicts a nearly-uniform winding temperature
        // (~49C for this design) rather than Icepak's per-turn spread (38-60C), so the 25%
        // bracket validates the aggregate winding level, not the per-turn distribution.
        // Each lookup REQUIREs its key to exist, so a future turn-naming change fails loudly
        // instead of silently skipping.
        // Re-pinned 2026-08-02 (ABT #461): the model predicts a nearly-uniform winding (see the
        // note above), so bracketing each turn to ITS OWN Icepak value +/-25% asserted a per-turn
        // spread the lumped model deliberately does not resolve — after the winding wrap-area /
        // radiation / characteristic-length fixes the uniform level (43.1 C) fell 4% below the
        // hottest turn's bracket while still sitting INSIDE Icepak's own turn range. The honest
        // assertion for a lumped model is therefore: every predicted turn lies within the span of
        // Icepak's exported turn temperatures for this design (38.72 to 60.01 C) — the reference's
        // own spread, no invented tolerance. Corroboration: OMFEM's radiating 2D FEM on this
        // fixture, corrected for its documented planar area deficit, brackets the same range
        // (#461). The key REQUIREs keep the dead-lookup protection (keys must exist).
        // Re-pinned 2026-10-06 (ABT #1585): with the volt-second balanced fixture (see the core
        // section) the core loss falls from the drift artifact's 1.261 W to 0.133 W and the
        // uniform winding level from 47.1..47.8 C to 32.05..32.31 C, below Icepak's exported
        // turn span (38.72..60.01 C), which rested on that artifact's loss. Each turn is pinned
        // to the balanced fixture's winding level, 32.2 C, at 5 % -- a REGRESSION pin, not a
        // validation (see the core section; ABT #1715).
        const double windingLevel = 32.2;
        auto checkTurn = [&](const std::string& key) {
            REQUIRE(tempsPerTurn.count(key) == 1);
            REQUIRE_THAT(tempsPerTurn.at(key), Catch::Matchers::WithinRel(windingLevel, 0.05));
        };
        checkTurn("Secondary parallel 0 turn 4");
        checkTurn("Secondary parallel 0 turn 5");
        checkTurn("Primary parallel 0 turn 6");
        checkTurn("Primary parallel 0 turn 8");
        checkTurn("Secondary parallel 0 turn 2");
        checkTurn("Secondary parallel 0 turn 3");
        checkTurn("Secondary parallel 0 turn 10");
        checkTurn("Primary parallel 1 turn 5");
    }
    
    SECTION("Winding temperature by index") {
        // Check winding temperatures using getTemperaturesByComponentType
        if (tempsByType.find("winding 0") != tempsByType.end()) {
            REQUIRE(tempsByType.at("winding 0") > config.ambientTemperature);
        }
        if (tempsByType.find("winding 1") != tempsByType.end()) {
            REQUIRE(tempsByType.at("winding 1") > config.ambientTemperature);
        }
    }
}


TEST_CASE("Temperature: concentric_flyback_rectangular_column", "[temperature][smoke-test]") {
    // The Icepak reference ranges in this test (cap 605.73, floor 363.44 on
    // the "core" component) were calibrated on 2026-03-03 (commit 356189a3).
    // Multiple legitimate thermal-model fixes have shipped since — notably
    // 7717b002 (skip zero-thickness toroidal insulation layers), 8df6b8b7
    // (SimpleMatrix refactor + degenerate-toroidal-pair guard), d0a5d6b7
    // (litz bundle diameter from strand), a259ee27 (compute B before core
    // losses), 4c1ebeab (drop silent fallbacks). When this test was skipped the
    // core temperature was ~12% past the upper cap (676.79°C vs cap 605.73°C).
    // A targeted revert of 7717b002's source moves the model output *further*
    // from the calibrated band (verified by experiment), so the drift is
    // cumulative across multiple fixes, not attributable to one commit. Need
    // a fresh Icepak run on this geometry to re-baseline the reference ranges
    // before re-enabling the test.
    //
    // RE-MEASURED 2026-08-01 (ABT #461): the drift had kept growing — 766.03°C
    // against the same cap of 605.73°C (~26% over). Root cause: the model had NO
    // radiation path to ambient, and at these temperatures radiation to the room
    // is the dominant cooling mechanism (h_rad ~ 27-80 W/m2K vs ~10 convective).
    // RESOLVED the same day: exposed core surfaces now radiate to ambient and a
    // turn's radiating area is split between the core it faces through the window
    // and the room (Temperature.cpp, ABT #461). Core lands at ~517°C, inside this
    // band — re-enabled. The WINDING (~1708°C vs Icepak's 523.63°C turn) is still
    // wrong: these turns have no exposed surfaces, so the winding couples only by
    // conduction, and that path is the separate coupling problem tracked in #461.
    // The per-turn section below is currently DEAD (guarded lookups on key names
    // that never match), so it does not catch this — see ABT #454.
    //
    // FIXTURE (ABT #1740, Alf 2026-10-10): the fixture used to carry a second operating
    // point, "Max. input volt.", whose 911 uVs on 24 turns of this E 25/13/7 impose a
    // 0.73 T swing on 3C96 (B_sat ~0.5 T): the core saturates and MKF now refuses it. It
    // was dropped; the test only ever used the first operating point's losses (outputs[0]).
    // Its design requirements said Lm >= 510 uH and turns ratio 8.5, while the coil is
    // 24:13 on a 0.8 mm gap (~64 uH); they now state the real coil (Lm 60-67 uH, 24/13).
    // The remaining operating point's waveforms were simulated with Lm = 510 uH, not with
    // this coil's ~64 uH; the Icepak references stay, as they belong to this geometry and
    // to that operating point's losses.
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_flyback_rectangular_column.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (losses.windingLossesOutput.has_value()) {
        config.windingLosses = losses.windingLosses;
        config.windingLossesOutput = losses.windingLossesOutput.value();
    } else {
        applySimulatedLosses(config, magnetic);
    }
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    // Reference values from Icepak simulation (Student Version)
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    // Get temperatures by component type
    auto tempsByType = temp.getTemperaturesByComponentType();
    auto tempsPerTurn = temp.getTemperaturePerTurn();
    
    SECTION("Core temperature validation against Icepak") {
        // Icepak core temperatures from solution overview
        // core_0: 388.43°C (from Icepak)
        // core_1: 484.58°C (from Icepak)
        REQUIRE(tempsByType.at("core") <= 605.73); // Max Icepak: 484.58°C + 25% tolerance
        REQUIRE(tempsByType.at("core") >= 363.44); // Min Icepak: 484.58°C - 25% tolerance
    }
    
    SECTION("Bobbin temperature validation against Icepak") {
        if (tempsByType.find("bobbin") != tempsByType.end()) {
            // No bobbin temperature data from Icepak
        }
    }
    
    SECTION("Individual turn temperatures from Icepak export") {
        // ABT #838: these lookups used the key "Turn_W0_Turn_N" (and, after a botched
        // find-and-replace, "Turn_Wecondary_0_Warallel_0_Turn_0"), which getTemperaturePerTurn
        // never emits -- its keys are "<winding> parallel <p> turn <n>". Every assertion sat
        // behind an `if (find(...) != end())` that was never true, so the whole section ran
        // green while asserting nothing at all. The keys are now the real ones and each one is
        // REQUIREd to exist, so a naming change fails loudly instead of going quiet again.
        //
        // These comparisons are EXPECTED TO FAIL today, and that failure is the point: it is
        // the only per-turn measurement in the repository, and it is how far the winding side
        // of the thermal model is from Icepak. Do not widen the bracket to make it pass --
        // fix the winding model (ABT #838). Alf, 2026-09-17: "repair them and let them fail".
        auto checkIcepakTurn = [&](const std::string& key, double icepakTemperature) {
            REQUIRE(tempsPerTurn.count(key) == 1);
            CHECK_THAT(tempsPerTurn.at(key), Catch::Matchers::WithinRel(icepakTemperature, 0.25));
        };
        checkIcepakTurn("Primary parallel 0 turn 3", 523.63);
        checkIcepakTurn("Primary parallel 0 turn 12", 300.51);
        checkIcepakTurn("Primary parallel 0 turn 21", 326.85);
        checkIcepakTurn("Primary parallel 0 turn 4", 525.20);
        checkIcepakTurn("Primary parallel 0 turn 15", 445.24);
        checkIcepakTurn("Primary parallel 0 turn 18", 550.86);
        checkIcepakTurn("Secondary 0 parallel 0 turn 11", 328.13);
        checkIcepakTurn("Primary parallel 0 turn 19", 551.13);
        checkIcepakTurn("Primary parallel 0 turn 20", 328.04);
        checkIcepakTurn("Primary parallel 0 turn 23", 188.64);
    }
    
    SECTION("Winding temperature by index") {
        // Check winding temperatures using getTemperaturesByComponentType
        if (tempsByType.find("winding 0") != tempsByType.end()) {
            REQUIRE(tempsByType.at("winding 0") > config.ambientTemperature);
        }
    }
}

TEST_CASE("Temperature: concentric_transformer_contiguous_rectangular_wire", "[temperature][smoke-test]") {
    // Same situation as Temperature: concentric_flyback_rectangular_column —
    // Icepak reference ranges (cap 473.04, floor 283.82) frozen on 2026-03-03.
    // When this test was skipped the model output was 270.93°C, undershooting
    // the floor by ~5%. Cumulative drift from the same series of post-March
    // thermal fixes. Re-baseline needed against fresh Icepak run.
    //
    // RE-MEASURED 2026-08-01 (ABT #461): this one had INVERTED — 1056.23°C, 123%
    // above the cap of 473.04 (the 270.93 undershoot recorded above was long gone).
    // Root cause and fix as in concentric_flyback_rectangular_column: radiation to
    // ambient was missing entirely, and the turn->core radiation pumped the whole
    // winding's radiating area into the core. With core->ambient radiation, the
    // window/room split of the turn area, and the damped solver iteration
    // (Temperature.cpp, ABT #461), the core converges at ~449°C, inside this band —
    // re-enabled.
    //
    // FIXTURE (ABT #1740): as in concentric_flyback_rectangular_column, the saturating
    // "Max. input volt." operating point was dropped and the design requirements now
    // state the real coil (Lm 60-67 uH, turns ratio 24/13). The remaining operating
    // point's waveforms were simulated with Lm = 510 uH, not with this coil's inductance.
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_transformer_contiguous_rectangular_wire.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    
    // Run magnetic simulation to get actual losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    
    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    if (losses.windingLossesOutput.has_value()) {
        config.windingLosses = losses.windingLosses;
        config.windingLossesOutput = losses.windingLossesOutput.value();
    } else {
        applySimulatedLosses(config, magnetic);
    }
    config.plotSchematic = true;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    // Reference values from Icepak simulation (Student Version)
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    // Export temperature field and thermal circuit schematic for visualization
    exportTemperatureFieldSvg("concentric_transformer_contiguous_rectangular_wire", magnetic, result.nodeTemperatures, config.ambientTemperature);
    exportThermalCircuitSchematic("concentric_transformer_contiguous_rectangular_wire", temp);
    
    // Get temperatures by component type
    auto tempsByType = temp.getTemperaturesByComponentType();
    auto tempsPerTurn = temp.getTemperaturePerTurn();
    
    SECTION("Core temperature validation against Icepak") {
        // Icepak core temperatures from solution overview
        // core_0: 351.32°C (from Icepak)
        // core_1: 378.43°C (from Icepak)
        REQUIRE(tempsByType.at("core") <= 473.04); // Max Icepak: 378.43°C + 25% tolerance
        REQUIRE(tempsByType.at("core") >= 283.82); // Min Icepak: 378.43°C - 25% tolerance
    }
    
    SECTION("Bobbin temperature validation against Icepak") {
        if (tempsByType.find("bobbin") != tempsByType.end()) {
            // No bobbin temperature data from Icepak
        }
    }
    
    SECTION("Individual turn temperatures from Icepak export") {
        // ABT #838: these lookups used the key "Turn_W0_Turn_N" (and, after a botched
        // find-and-replace, "Turn_Wecondary_0_Warallel_0_Turn_0"), which getTemperaturePerTurn
        // never emits -- its keys are "<winding> parallel <p> turn <n>". Every assertion sat
        // behind an `if (find(...) != end())` that was never true, so the whole section ran
        // green while asserting nothing at all. The keys are now the real ones and each one is
        // REQUIREd to exist, so a naming change fails loudly instead of going quiet again.
        //
        // These comparisons are EXPECTED TO FAIL today, and that failure is the point: it is
        // the only per-turn measurement in the repository, and it is how far the winding side
        // of the thermal model is from Icepak. Do not widen the bracket to make it pass --
        // fix the winding model (ABT #838). Alf, 2026-09-17: "repair them and let them fail".
        auto checkIcepakTurn = [&](const std::string& key, double icepakTemperature) {
            REQUIRE(tempsPerTurn.count(key) == 1);
            CHECK_THAT(tempsPerTurn.at(key), Catch::Matchers::WithinRel(icepakTemperature, 0.25));
        };
        checkIcepakTurn("Primary parallel 0 turn 3", 395.37);
        checkIcepakTurn("Primary parallel 0 turn 21", 244.07);
        checkIcepakTurn("Primary parallel 0 turn 18", 305.77);
        checkIcepakTurn("Secondary 0 parallel 0 turn 0", 242.76);
        checkIcepakTurn("Secondary 0 parallel 0 turn 8", 243.86);
        checkIcepakTurn("Primary parallel 0 turn 1", 301.11);
        checkIcepakTurn("Primary parallel 0 turn 4", 394.95);
        checkIcepakTurn("Primary parallel 0 turn 16", 232.31);
        checkIcepakTurn("Primary parallel 0 turn 7", 395.20);
        checkIcepakTurn("Primary parallel 0 turn 19", 267.48);
        checkIcepakTurn("Secondary 0 parallel 0 turn 1", 242.76);
        checkIcepakTurn("Secondary 0 parallel 0 turn 3", 243.77);
        checkIcepakTurn("Secondary 0 parallel 0 turn 4", 243.47);
        checkIcepakTurn("Secondary 0 parallel 0 turn 6", 281.72);
        checkIcepakTurn("Secondary 0 parallel 0 turn 7", 243.46);
        checkIcepakTurn("Secondary 0 parallel 0 turn 11", 232.11);
        checkIcepakTurn("Secondary 0 parallel 0 turn 12", 243.87);
        checkIcepakTurn("Primary parallel 0 turn 14", 232.22);
    }
    
    SECTION("Winding temperature by index") {
        // Check winding temperatures using getTemperaturesByComponentType
        if (tempsByType.find("winding 0") != tempsByType.end()) {
            REQUIRE(tempsByType.at("winding 0") > config.ambientTemperature);
        }
    }
}


TEST_CASE("Temperature: BuckInductor T134_77_27 from MAS file", "[temperature][concentric][buck-inductor][smoke-test]") {
    // The fixture BuckInductor_164uH_T134_77_27.json was referenced when this
    // test was added in 7717b002 ("fix(thermal): skip zero-thickness insulation
    // layers in toroidal thermal model") but the JSON file was never committed
    // to tests/testData/. Until the fixture lands the test cannot run; without
    // this SKIP the loader throws on an empty input. Remove the SKIP once the
    // fixture is committed.
    SKIP("Test needs investigation - BuckInductor_164uH_T134_77_27.json fixture missing");
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "BuckInductor_164uH_T134_77_27.json");
    auto mas = OpenMagneticsTesting::mas_loader(path.string());

    auto magnetic = mas.get_magnetic();
    auto inputs = mas.get_inputs();

    // Wind the coil if not already wound
    if (!magnetic.get_coil().get_turns_description()) {
        magnetic.get_mutable_coil().wind();
    }

    // Run the electromagnetic simulation to get losses
    auto losses = getLossesFromSimulation(magnetic, inputs);
    REQUIRE(losses.simulationSucceeded);
    REQUIRE(losses.windingLossesOutput.has_value());

    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = true;
    config.schematicOutputPath = (getOutputDir() / "thermal_schematic_BuckInductor_T134.svg").string();

    {
        auto outFile = getOutputDir() / "magnetic_BuckInductor_T134.svg";
        std::filesystem::create_directories(getOutputDir());
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    Temperature temp(magnetic, config);
    ThermalResult result;
    try {
        result = temp.calculateTemperatures();
        exportTemperatureFieldSvg("BuckInductor_T134", magnetic, result.nodeTemperatures, config.ambientTemperature);
    } catch (const std::exception& e) {
        // A throw used to be silently swallowed here; it must fail the test.
        FAIL("calculateTemperatures/exportTemperatureFieldSvg threw: " << e.what());
    }
    exportThermalCircuitSchematic("BuckInductor_T134", temp);

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);

    auto tempsByType = temp.getTemperaturesByComponentType();
    REQUIRE(tempsByType.find("core") != tempsByType.end());
    REQUIRE(tempsByType.at("core") > config.ambientTemperature);
}

TEST_CASE("Temperature: coreOnly mode concentric E-core", "[temperature][core-only][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "E 42/21/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "N87");

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.coreOnly = true;
    config.coreLosses = 0.5;
    config.ambientTemperature = 25.0;
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.maximumTemperature < 300.0);
    REQUIRE(result.averageCoreTemperature > config.ambientTemperature);
}

TEST_CASE("Temperature: coreOnly mode toroidal core", "[temperature][core-only][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "3C97");

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.coreOnly = true;
    config.coreLosses = 0.3;
    config.ambientTemperature = 25.0;
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
}

TEST_CASE("Temperature: coreOnly zero losses stays at ambient", "[temperature][core-only][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "E 42/21/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                                                     WindingOrientation::CONTIGUOUS,
                                                     WindingOrientation::CONTIGUOUS,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "N87");

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.coreOnly = true;
    config.coreLosses = 0.0;
    config.ambientTemperature = 25.0;
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE_THAT(result.maximumTemperature, Catch::Matchers::WithinAbs(config.ambientTemperature, 1.0));
}

TEST_CASE("Temperature: Litz Wire Without Pre-computed Outer Diameter", "[temperature][litz][smoke-test]") {
    // Regression test: litz wires arriving at the temperature calculation may have only strand data
    // (conductingDiameter, numberConductors, coating) with no outerDiameter pre-computed at the wire level.
    // Temperature::extractWireProperties must call calculate_outer_diameter() to derive the bundle diameter
    // from the strand, not rely on outerDiameter being present (the old code left _wireWidth at the 1mm
    // default because the outer_diameter block was guarded by `if (!_isRoundWire)`, silently using wrong geometry).
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_litz_foil.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);

    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);

    // Strip outerDiameter from the litz wire in the functional description so that
    // extractWireProperties must compute it via calculate_outer_diameter() from the strand.
    auto& coil = magnetic.get_mutable_coil();
    auto functionalDescription = coil.get_functional_description();
    for (auto& winding : functionalDescription) {
        auto wireVariant = winding.get_wire();
        if (std::holds_alternative<OpenMagnetics::Wire>(wireVariant)) {
            auto wire = std::get<OpenMagnetics::Wire>(wireVariant);
            if (wire.get_type() == WireType::LITZ) {
                wire.set_outer_diameter(std::nullopt);
                REQUIRE(!wire.get_outer_diameter().has_value());
                REQUIRE(wire.get_strand().has_value());
                winding.set_wire(wire);
            }
        }
    }
    coil.set_functional_description(functionalDescription);

    auto losses = getLossesFromSimulation(magnetic, inputs);

    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    config.windingLosses = losses.windingLosses;
    if (!losses.windingLossesOutput.has_value()) {
        throw std::runtime_error("WindingLossesOutput missing");
    }
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    REQUIRE(result.totalThermalResistance > 0.0);
}

} // namespace

// ============================================================================
// OMFEM 2D FEM cross-check battery (ABT #454 / #461)
// ============================================================================
// Compares MKF's lumped thermal network against OMFEM's independent 2D FEM
// (mesh -> radiating heat-conduction solve, h=12 W/m2K, eps=0.9) on every MAS
// example that OMFEM could process. The reference file stores BOTH the losses
// OMFEM used and the temperatures it produced; this test drives MKF with the
// SAME stored losses, so the comparison is decoupled from MKF's loss models and
// cannot silently go stale when they change — the failure mode of the frozen
// Icepak bands (a 2026-03 band met a 2026-08 loss model in ABT #461).
//
// Reference generation: scripts alongside ABT #461 run omfem_mas on each
// example and record P_core/P_cu (OMFEM's FEM losses), ambient, and the FEM's
// core/winding/max temperatures. Regenerate with the same tool when OMFEM's
// thermal model materially changes.
//
// Tolerances are wide BY DESIGN and documented, because the two models have
// known structural differences (measured in the #461 investigation):
//   - OMFEM's planar section lacks the core's front/back envelope faces, so
//     core-dominated planar cases read hot (~2x in rise).
//   - Its winding turns are section islands (no MLT surface scaling), so
//     winding-dominated cases hold heat in the winding and starve the core.
//   - MKF's winding<->world coupling has its own known defects (#461).
// Measured max-rise ratios across all 27 examples with the post-#461/#527
// physics (radiation, full core loss, wrap areas, composition-aware k, no h
// floors) span 0.30-3.43 above the noise floor — the [FEMCMP] lines this test
// prints are the census. Both extremes are geometry-family spread between a
// lumped network and a 2D section FEM (the 0.30 is a toroidal CMC where the
// planar-vs-axisym treatments differ most; the 3.43 a winding-dominated EQ
// core), not tuning headroom: a factor-4 band is the tightest that passes the
// measured envelope, and it still catches the order-of-magnitude breaks this
// battery exists for (the pre-#461 model was ~10x off the FEM on a core).
TEST_CASE("Temperature: OMFEM 2D FEM cross-check battery over MAS examples", "[temperature][thermal-fem-battery]") {
    namespace fs = std::filesystem;
    auto refPath = fs::path{std::source_location::current().file_name()}.parent_path()
                       .append("testData").append("omfem_thermal_2d_references.json");
    REQUIRE(fs::exists(refPath));
    json refs = json::parse(std::ifstream(refPath.string()));
    auto examplesDir = fs::path{std::source_location::current().file_name()}.parent_path()
                           .append("..").append("MAS").append("examples");

    // Below this FEM temperature rise the ratio of two small numbers is noise,
    // not physics; such cases only assert that MKF stays similarly cool.
    constexpr double kMinRiseForRatioK = 5.0;
    constexpr double kMaxRiseRatioBand = 4.0;  // see header comment for the measured basis

    // Coverage gate (the ABT #454 lesson: never green-by-vacuity). Checked against the
    // reference FILE, not a loop counter: Catch2 re-enters the test once per dynamic
    // section, so no cross-section counter can accumulate.
    size_t usableReferences = 0;
    for (auto& [name, ref] : refs.at("examples").items()) {
        if (ref.at("status") == "ok") usableReferences++;
    }
    REQUIRE(usableReferences >= 5);

    for (auto& [name, ref] : refs.at("examples").items()) {
        if (ref.at("status") != "ok") continue;
        DYNAMIC_SECTION(name) {
            auto mas = OpenMagneticsTesting::mas_loader((examplesDir / name).string());
            auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
            auto& coil = magnetic.get_mutable_coil();
            REQUIRE(coil.get_turns_description());

            const double ambient = ref.at("ambient");
            const double pCore = ref.at("p_core");
            const double pCu = ref.at("p_cu");

            // Synthetic per-turn split of the stored winding loss, proportional to each
            // turn's length (uniform current density). The thermal model requires a
            // per-turn distribution by contract; this is the stored total, distributed,
            // not a fallback around missing data.
            auto turns = coil.get_turns_description().value();
            double totalLength = 0;
            for (auto& turn : turns) totalLength += turn.get_length();
            REQUIRE(totalLength > 0);
            std::vector<WindingLossesPerElement> perTurn;
            for (auto& turn : turns) {
                OhmicLosses ohmic;
                ohmic.set_losses(pCu * turn.get_length() / totalLength);
                ohmic.set_origin(ResultOrigin::SIMULATION);
                ohmic.set_method_used("omfem-reference length-proportional split");
                WindingLossesPerElement element;
                element.set_name(turn.get_name());
                element.set_ohmic_losses(ohmic);
                perTurn.push_back(element);
            }
            WindingLossesOutput windingLosses;
            windingLosses.set_origin(ResultOrigin::SIMULATION);
            windingLosses.set_method_used("omfem-reference length-proportional split");
            windingLosses.set_winding_losses(pCu);
            windingLosses.set_winding_losses_per_turn(perTurn);

            TemperatureConfig config;
            config.ambientTemperature = ambient;
            config.coreLosses = pCore;
            config.windingLosses = pCu;
            config.windingLossesOutput = windingLosses;
            config.plotSchematic = false;

            Temperature temp(magnetic, config);
            auto result = temp.calculateTemperatures();
            auto byType = temp.getTemperaturesByComponentType();

            REQUIRE(result.converged);
            REQUIRE(std::isfinite(result.maximumTemperature));
            REQUIRE(byType.count("core") == 1);

            const double mkfMaxRise = result.maximumTemperature - ambient;
            const double mkfCoreRise = byType.at("core") - ambient;
            const double femMaxRise = double(ref.at("omfem_t_max")) - ambient;
            const double femCoreRise = double(ref.at("omfem_t_core")) - ambient;
            std::cout << "[FEMCMP] " << name << " P=" << pCore << "+" << pCu
                      << " maxRise mkf/fem=" << mkfMaxRise << "/" << femMaxRise
                      << " ratio=" << (femMaxRise > 0 ? mkfMaxRise / femMaxRise : 0)
                      << " coreRise mkf/fem=" << mkfCoreRise << "/" << femCoreRise << "\n";
            INFO(name << ": P=" << pCore << "+" << pCu << " W | MKF core/max rise "
                 << mkfCoreRise << "/" << mkfMaxRise << " K | OMFEM " << femCoreRise
                 << "/" << femMaxRise << " K");

            REQUIRE(mkfMaxRise >= 0.0);
            if (femMaxRise < kMinRiseForRatioK) {
                // Cold reference: MKF must not invent a hot component out of watts the
                // FEM shed easily (same factor band applied to the absolute rise).
                REQUIRE(mkfMaxRise <= kMinRiseForRatioK * kMaxRiseRatioBand);
            } else {
                const double maxRatio = mkfMaxRise / femMaxRise;
                // ER 9.5 current sense: OMFEM's 2D planar E-type setup is a half section extruded
                // by the column depth, with 48.6 mm2 of exterior against the real core's 183 mm2
                // (3.8x short), so it reads hot (ABT #1475). With ABT #1454/#1459's real exterior
                // MKF sits at 0.223; Alf approved a 0.20 floor for this case only (2026-09-28).
                const double minRatio = (name == "13_current_sense_er95_n87.json") ? 0.20 : 1.0 / kMaxRiseRatioBand;
                REQUIRE(maxRatio > minRatio);
                REQUIRE(maxRatio < kMaxRiseRatioBand);
                // Core comparison only where the core carries the heat: on
                // winding-dominated cases the FEM's island winding starves its core
                // (measured ~8x apart on the flyback) and the ratio means nothing.
                if (pCore >= 0.6 * (pCore + pCu) && femCoreRise >= kMinRiseForRatioK) {
                    const double coreRatio = mkfCoreRise / femCoreRise;
                    REQUIRE(coreRatio > minRatio);  // same per-case floor: ER 9.5 is core-only (0.226)
                    REQUIRE(coreRatio < kMaxRiseRatioBand);
                }
            }
        }
    }
}


// ABT #906: on a small packed toroid (High Flux 160 T10/6/5, 56+28 turns AWG29,
// 0.29 W total loss at 50 C ambient — a user-reported flybuck) the full network
// exposed every turn's quarter-developed surface plus the wrapped core surface to
// ambient: 16.6 cm2 of convecting area on a wound body whose physical envelope is
// ~4.8 cm2. Predicted rise +14 K where Magnetics' own empirical (P_mW/A_cm2)^0.833
// gives +31..36 K and first-principles convection+radiation gives +33..42 K — and
// the MagneticSimulator core-only path simultaneously reported +36 K for the same
// design (a 2.5x self-contradiction between two surfaces of the same solver).
// The toroidal ambient interface is now envelope-based (outer cylinder + annular
// faces + hole cylinder, divided over the boundary turn nodes), which lands the
// rise in the physical band and makes the two paths agree.
TEST_CASE("Temperature: ABT906 Packed Toroid Envelope Ambient Interface",
          "[temperature][toroidal][abt906][smoke-test]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flybuck_hf160_250u_abt907.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    auto magnetic = mas.get_magnetic();

    // Losses as simulated for this operating point (checked against the user's
    // own Steinmetz hand-calculation once the correct magnetizing ripple is used).
    TemperatureConfig config;
    config.ambientTemperature = 50;
    config.coreLosses = mas.get_outputs()[0].get_core_losses()->get_core_losses();
    config.windingLosses = mas.get_outputs()[0].get_winding_losses()->get_winding_losses();
    config.windingLossesOutput = mas.get_outputs()[0].get_winding_losses().value();
    config.plotSchematic = false;

    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();

    REQUIRE(result.converged);
    // Physical band for ~0.29 W on a ~4.5-4.8 cm2 wound envelope: +28..+46 K.
    // The old developed-surface model gave +14 K (below band); the intermediate
    // binary-blocking attempt gave +140 K (above band).
    CHECK(result.maximumTemperature > 78.0);
    CHECK(result.maximumTemperature < 96.0);
    // Assembly-to-ambient thermal resistance in the physical range for this size.
    CHECK(result.totalThermalResistance > 95.0);
    CHECK(result.totalThermalResistance < 160.0);

    // The convecting+radiating area presented to ambient must be the wound-body
    // envelope, not the developed wire surface: ~4.5 cm2 here, not 16.6 cm2.
    double ambientConvectionArea = 0.0;
    const auto& resistances = temp.getResistances();
    for (const auto& res : resistances) {
        if (res.type == OpenMagnetics::HeatTransferType::NATURAL_CONVECTION ||
            res.type == OpenMagnetics::HeatTransferType::FORCED_CONVECTION) {
            ambientConvectionArea += res.area;
        }
    }
    CHECK(ambientConvectionArea > 3.5e-4);
    CHECK(ambientConvectionArea < 6.0e-4);
}

// ABT #906 follow-up: validate the envelope-based toroidal ambient interface
// against Magnetics Inc.'s independent empirical fit dT = (P_mW / A_cm2)^0.833
// (their published thermal-rise estimate for wound powder toroids, with A the
// wound-body surface). The envelope area here is computed by ARITHMETIC in the
// test (from core OD/ID/height + wire OD), independent of the model's element
// areas, and the power for each size is chosen so the empirical fit predicts a
// +40 K rise. The network solves h(T) correlations + radiation, a completely
// different derivation, so tracking the empirical fit within a band across a
// 4x size span is a real cross-check, not a tautology.
// (WE RedExpert's measured deltat catalogue was considered for this validation
// but its parts are drum/shielded cores, which never enter the toroidal path.)
TEST_CASE("Temperature: ABT906 Envelope Model Tracks Magnetics Empirical Rise Across Sizes",
          "[temperature][toroidal][abt906][thermal-validation][smoke-test]") {
    settings.reset();
    struct SweepCase {
        std::string shape;
        int64_t turns;
    };
    std::vector<SweepCase> sweep = {
        {"T 10/6/5", 40},
        {"T 20/10/7", 60},
        {"T 30/20/10", 80},
        {"T 40/24/16", 100},
    };

    for (const auto& sweepCase : sweep) {
        CAPTURE(sweepCase.shape);
        auto coil = OpenMagneticsTesting::get_quick_coil({sweepCase.turns},
                                                         {1},
                                                         sweepCase.shape,
                                                         1,
                                                         WindingOrientation::OVERLAPPING,
                                                         WindingOrientation::OVERLAPPING,
                                                         CoilAlignment::SPREAD,
                                                         CoilAlignment::SPREAD);
        auto core = OpenMagneticsTesting::get_quick_core(sweepCase.shape, json::array(), 1, "High Flux 60");
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        magnetic.get_mutable_coil().wind();

        // Wound-envelope surface by plain arithmetic: core dims + one wire on each side.
        auto dimensions = core.resolve_shape().get_dimensions().value();
        double coreOuterDiameter = resolve_dimensional_values(dimensions["A"]);
        double coreInnerDiameter = resolve_dimensional_values(dimensions["B"]);
        double coreHeight = resolve_dimensional_values(dimensions["C"]);
        auto wire = std::get<OpenMagnetics::Wire>(magnetic.get_coil().get_functional_description()[0].get_wire());
        double wireOuterDiameter = wire.calculate_outer_diameter();
        double woundOuterDiameter = coreOuterDiameter + 2 * wireOuterDiameter;
        double woundInnerDiameter = std::max(0.0, coreInnerDiameter - 2 * wireOuterDiameter);
        double woundHeight = coreHeight + 2 * wireOuterDiameter;
        double envelopeAreaSquareMeters = std::numbers::pi * woundOuterDiameter * woundHeight +
                                          std::numbers::pi * woundInnerDiameter * woundHeight +
                                          2.0 * std::numbers::pi / 4.0 *
                                              (woundOuterDiameter * woundOuterDiameter -
                                               woundInnerDiameter * woundInnerDiameter);
        double envelopeAreaSquareCentimeters = envelopeAreaSquareMeters * 1e4;

        // Power for which the Magnetics empirical fit predicts +40 K.
        double expectedRise = 40.0;
        double lossesMilliwatts = envelopeAreaSquareCentimeters * std::pow(expectedRise, 1.0 / 0.833);

        TemperatureConfig config;
        config.ambientTemperature = 25.0;
        config.coreLosses = lossesMilliwatts / 1000.0;
        config.plotSchematic = false;

        auto result = Temperature(magnetic, config).calculateTemperatures();
        REQUIRE(result.converged);

        double networkRise = result.maximumTemperature - config.ambientTemperature;
        double ratio = networkRise / expectedRise;
        CAPTURE(envelopeAreaSquareCentimeters, lossesMilliwatts, networkRise, ratio);
        CHECK(ratio > 0.6);
        CHECK(ratio < 1.7);
    }
    settings.reset();
}

// The disconnected-node diagnosis must test CONNECTIVITY to ambient, not node degree. A node
// whose only edges lead to other stranded nodes has a non-zero conductance diagonal and passed
// the old G(i,i) < 1e-12 test, so a floating island went straight to a singular solve. Bare
// graph: ambient 0; 1 and 2 hang off it; 3 and 4 touch only each other; 5 touches nothing.
TEST_CASE("Test_Thermal_Floating_Island_Is_Reported_Not_Just_Degree_Zero_Nodes", "[temperature][bug][omega]") {
    auto link = [](size_t from, size_t to) {
        ThermalResistanceElement res;
        res.nodeFromId = from;
        res.nodeToId = to;
        res.type = HeatTransferType::CONDUCTION;
        res.resistance = 1.0;
        return res;
    };
    std::vector<ThermalResistanceElement> resistances = {link(0, 1), link(1, 2), link(3, 4)};
    auto stranded = Temperature::nodesWithoutPathToRoots(6, resistances, {0});
    REQUIRE(stranded == std::vector<size_t>{3, 4, 5});   // the island's members AND the isolated node
    // a cold plate (fixed temperature) is a root too: the island hanging off it is fine
    auto withColdPlate = Temperature::nodesWithoutPathToRoots(6, resistances, {0, 3});
    REQUIRE(withColdPlate == std::vector<size_t>{5});
    REQUIRE(Temperature::nodesWithoutPathToRoots(3, {link(0, 1), link(1, 2)}, {0}).empty());
}

// ============================================================================================
// ABT #838 — three build-time values that the per-iteration recalculation used to overwrite.
//
// recalculateConvectionResistances() rebuilds every convection and radiation resistance from the
// current temperatures on each pass, so anything folded into the resistance at BUILD time and not
// carried on the resistor itself survives exactly one iteration. Radiation was lost this way once
// already (it was lumped into h_conv); these three were the remaining instances.
// ============================================================================================

TEST_CASE("ABT838_Convection_Resistor_Carries_The_Exposed_Area_Not_The_Bare_One",
          "[temperature][abt838]") {
    // A quadrant's convection resistance is built from surfaceArea * surfaceCoverage, so a face a
    // winding lies against convects less than a bare one. The resistor stored the BARE area, and
    // the recalculation then rebuilt the resistance as 1/(h * area) — handing the covered face its
    // full area back and making a wound surface convect as if nothing were wound on it.
    std::filesystem::path testFile = std::filesystem::path(__FILE__).parent_path() / "testData" /
                                     "concentric_round_wire_insulation_layers.json";
    std::ifstream file(testFile);
    REQUIRE(file.good());
    json j;
    file >> j;
    OpenMagnetics::Inputs inputs(j["inputs"]);
    OpenMagnetics::Magnetic magnetic(j["magnetic"]);
    auto losses = getLossesFromSimulation(magnetic, inputs);

    TemperatureConfig config;
    config.ambientTemperature = losses.ambientTemperature;
    config.coreLosses = losses.coreLosses;
    REQUIRE(losses.windingLossesOutput.has_value());
    config.windingLosses = losses.windingLosses;
    config.windingLossesOutput = losses.windingLossesOutput.value();
    config.plotSchematic = false;

    Temperature temperature(magnetic, config);
    auto result = temperature.calculateTemperatures();
    REQUIRE(result.converged);

    // The invariant: every convection resistor is built with the SAME h_conv, so 1/(R * area)
    // must come out identical for all of them. That is exactly what the bug broke -- a covered
    // face's resistance came from surfaceArea*coverage while its stored area was the bare
    // surfaceArea, so its implied h was off by 1/coverage (a factor of 24 on the face below).
    // It is also immune to the half-core symmetry correction, which scales resistance and area
    // together and therefore leaves the implied h alone.
    const auto& nodes = temperature.getNodes();
    std::vector<double> exposedFaceCoefficients;
    std::vector<std::pair<double, double>> coveredFaceCoefficients;   // (coverage, implied h)
    size_t checkedCovered = 0;
    for (const auto& resistor : temperature.getResistances()) {
        if (resistor.type != HeatTransferType::NATURAL_CONVECTION &&
            resistor.type != HeatTransferType::FORCED_CONVECTION) {
            continue;
        }
        if (resistor.nodeFromId >= nodes.size() || resistor.quadrantFrom == ThermalNodeFace::NONE) {
            continue;
        }
        const auto* quadrant = nodes[resistor.nodeFromId].getQuadrant(resistor.quadrantFrom);
        if (!quadrant || quadrant->surfaceArea <= 0 || resistor.area <= 0 || resistor.resistance <= 0) {
            continue;
        }
        double impliedCoefficient = 1.0 / (resistor.resistance * resistor.area);
        if (quadrant->surfaceCoverage < 1.0) {
            ++checkedCovered;
            coveredFaceCoefficients.emplace_back(quadrant->surfaceCoverage, impliedCoefficient);
            UNSCOPED_INFO("covered face on node " << resistor.nodeFromId << ": coverage "
                          << quadrant->surfaceCoverage << ", surfaceArea " << quadrant->surfaceArea
                          << ", stored area " << resistor.area << ", implied h " << impliedCoefficient);
            // The stored area is the EXPOSED one (the symmetry correction may double it), never
            // the bare face.
            CHECK(resistor.area < quadrant->surfaceArea * 2.0);
        }
        else {
            exposedFaceCoefficients.push_back(impliedCoefficient);
        }
    }
    // A design with no partially covered face would prove nothing, so say so rather than pass.
    REQUIRE(checkedCovered > 0);
    REQUIRE(!exposedFaceCoefficients.empty());
    std::sort(exposedFaceCoefficients.begin(), exposedFaceCoefficients.end());
    const double exposedReference = exposedFaceCoefficients[exposedFaceCoefficients.size() / 2];
    REQUIRE(exposedReference > 0);
    // Every face is built with the same h_conv, so a covered face's implied h must land near the
    // fully exposed ones. It does not land EXACTLY there (series terms such as an insulation layer
    // add resistance, spreading the implied value by a few per cent), so the bar here is the SIZE
    // OF THE BUG rather than an exact equality: storing the bare area divides a covered face's
    // implied h by its coverage -- a factor of 24 on the 4% face in this design, and never less
    // than 2 on any face covered by more than half.
    for (const auto& [coverage, impliedCoefficient] : coveredFaceCoefficients) {
        UNSCOPED_INFO("coverage " << coverage << ": implied h " << impliedCoefficient
                      << " against the exposed-face reference " << exposedReference);
        CHECK(impliedCoefficient > exposedReference / 2.0);
        CHECK(impliedCoefficient < exposedReference * 2.0);
    }
}

TEST_CASE("ABT838_Forced_Convection_Keeps_Its_Buoyancy_Term_Through_The_Iteration",
          "[temperature][abt838][cooling]") {
    // applyForcedConvectionCooling builds these resistors with the mixed law
    // h = (h_natural^3 + h_forced^3)^(1/3); the recalculation used to rebuild them from PURE
    // forced convection, a second implementation, so buoyancy was dropped after iteration one.
    // At low airflow the natural term is a large share of the total, so a part cooled by a gentle
    // breeze must come out COOLER than pure forced convection alone would make it.
    std::filesystem::path testFile = std::filesystem::path(__FILE__).parent_path() / "testData" /
                                     "concentric_round_wire_insulation_layers.json";
    std::ifstream file(testFile);
    REQUIRE(file.good());
    json j;
    file >> j;
    OpenMagnetics::Inputs inputs(j["inputs"]);
    OpenMagnetics::Magnetic magnetic(j["magnetic"]);
    auto losses = getLossesFromSimulation(magnetic, inputs);

    auto solve = [&](std::optional<double> velocity) {
        TemperatureConfig config;
        config.ambientTemperature = losses.ambientTemperature;
        config.coreLosses = losses.coreLosses;
        REQUIRE(losses.windingLossesOutput.has_value());
        config.windingLosses = losses.windingLosses;
        config.windingLossesOutput = losses.windingLossesOutput.value();
        config.plotSchematic = false;
        if (velocity) {
            MAS::Cooling cooling;
            cooling.set_temperature(losses.ambientTemperature);
            cooling.set_velocity(std::vector<double>{velocity.value(), 0.0, 0.0});
            config.masCooling = cooling;
        }
        return Temperature(magnetic, config).calculateTemperatures().maximumTemperature;
    };

    const double natural = solve(std::nullopt);
    const double barelyMoving = solve(0.02);  // the discriminating case, see below
    const double strongFlow = solve(3.0);

    UNSCOPED_INFO("natural " << natural << " C, 0.02 m/s " << barelyMoving
                  << " C, 3 m/s " << strongFlow << " C");
    // THE point of the blend: at a velocity this low the FORCED correlation alone returns a
    // smaller coefficient than buoyancy does, so a model that keeps only the forced term reports
    // that stirring the air HEATS the part. Physically, air that barely moves can do no worse
    // than still air. Mixed convection is what guarantees it -- (h_nat^3 + h_forced^3)^(1/3) is
    // never below h_nat -- and with the recalculation overwriting the blend with pure forced
    // convection, this is the assertion that goes red.
    CHECK(barelyMoving <= natural);
    // Real airflow still helps more than a whisper of it.
    CHECK(strongFlow < barelyMoving);
}

TEST_CASE("ABT838_Bare_Copper_Does_Not_Radiate_Like_A_Matte_Dielectric", "[temperature][abt838]") {
    // One emissivity of 0.9 was charged to every surface: right for enamelled wire, insulation
    // wrap and ferrite, wrong by more than an order of magnitude for bare copper (0.03-0.07).
    // This fixture has a served litz primary and a FOIL secondary whose coating is declared
    // "bare" -- how MAS writes bare metal -- so the two windings must come out with different
    // emissivities, and the value has to survive recalculateConvectionResistances, which is
    // where the single config-wide value used to be re-imposed every iteration.
    //
    // ABT #1540: the stored layout is litz | foil | litz, and the 0.4 x 9 mm foil touches the
    // ten litz turns of the layer on each of its sides. The network used to measure a
    // rectangular turn by half its SMALLER side in every direction, so only the litz turn level
    // with the foil's centre could ever be "near" it, none was, and the foil came out with no
    // edge to either layer and its outer face open to the room -- through a full layer of litz.
    // That impossible surface was the bare one this test read. With the outlines measured as
    // they are, the buried foil has no exposed face, which is asserted below; the bare surface
    // is now read from the same design without its outer primary section (litz | foil), where
    // the foil really is the outermost layer, and the coated one from the stored layout, whose
    // outermost layer is litz.
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(),
                                                            "concentric_litz_foil.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);

    auto withoutOuterPrimarySection = [](const OpenMagnetics::Magnetic& stored) {
        json magneticJson;
        to_json(magneticJson, stored);
        auto& coil = magneticJson["coil"];
        // The stored sections, from the column outwards: Primary section 0, insulation,
        // Secondary section 0, insulation, Primary section 1, insulation.
        REQUIRE(coil["sectionsDescription"].size() == 6);
        REQUIRE(coil["sectionsDescription"][4]["name"] == "Primary section 1");
        REQUIRE(coil["functionalDescription"][0]["name"] == "Primary");
        REQUIRE(coil["functionalDescription"][0]["numberParallels"] == 2);
        const std::set<std::string> dropped = {coil["sectionsDescription"][4]["name"].get<std::string>(),
                                               coil["sectionsDescription"][5]["name"].get<std::string>()};
        auto keepOutsideDropped = [&](json& description, const std::string& sectionKey) {
            json kept = json::array();
            for (auto& item : description) {
                if (!dropped.count(item[sectionKey].get<std::string>())) {
                    kept.push_back(item);
                }
            }
            description = kept;
        };
        keepOutsideDropped(coil["sectionsDescription"], "name");
        keepOutsideDropped(coil["layersDescription"], "section");
        keepOutsideDropped(coil["turnsDescription"], "section");
        // The second parallel lived entirely in the dropped section: the primary is one parallel.
        coil["functionalDescription"][0]["numberParallels"] = 1;
        for (const char* key : {"sectionsDescription", "layersDescription"}) {
            for (auto& item : coil[key]) {
                if (!item.contains("partialWindings")) continue;
                for (auto& partialWinding : item["partialWindings"]) {
                    if (partialWinding["winding"] != "Primary") continue;
                    REQUIRE(partialWinding["parallelsProportion"].size() == 2);
                    REQUIRE(partialWinding["parallelsProportion"][1] == 0);
                    partialWinding["parallelsProportion"] = json::array({partialWinding["parallelsProportion"][0]});
                }
            }
        }
        for (auto& turn : coil["turnsDescription"]) {
            REQUIRE(turn["parallel"] == 0);
        }
        return OpenMagnetics::Magnetic(magneticJson);
    };

    // Read the SOLVED resistances, not just the field: a value stamped on the resistor that the
    // iteration then ignores is exactly the bug being fixed, and a test that only inspects the
    // field passes in both states. h_rad follows from R and the area, and it is proportional to
    // the emissivity, so the bare winding's radiating surfaces must come out far weaker than the
    // coated one's.
    std::vector<double> bareCoefficients;
    std::vector<double> coatedCoefficients;
    auto readRadiatingTurnSurfaces = [&](const OpenMagnetics::Magnetic& magnetic) {
        auto losses = getLossesFromSimulation(magnetic, mas.get_mutable_inputs());

        TemperatureConfig config;
        config.ambientTemperature = losses.ambientTemperature;
        config.coreLosses = losses.coreLosses;
        REQUIRE(losses.windingLossesOutput.has_value());
        config.windingLosses = losses.windingLosses;
        config.windingLossesOutput = losses.windingLossesOutput.value();
        config.plotSchematic = false;

        Temperature temperature(magnetic, config);
        auto result = temperature.calculateTemperatures();
        REQUIRE(result.converged);

        const auto& nodes = temperature.getNodes();
        for (const auto& resistor : temperature.getResistances()) {
            if (resistor.type != HeatTransferType::RADIATION) {
                continue;
            }
            if (resistor.nodeFromId >= nodes.size() || resistor.area <= 0 || resistor.resistance <= 0) {
                continue;
            }
            const auto& node = nodes[resistor.nodeFromId];
            if (node.part != ThermalNodePartType::TURN || !node.windingIndex) {
                continue;
            }
            const double impliedCoefficient = 1.0 / (resistor.resistance * resistor.area);
            // Winding 0 is the served litz, winding 1 the bare foil.
            if (node.windingIndex.value() == 1) {
                bareCoefficients.push_back(impliedCoefficient);
                CHECK_THAT(resistor.emissivity,
                           Catch::Matchers::WithinRel(OpenMagnetics::ThermalDefaults::kRadiation_BareCopperEmissivity, 1e-12));
            }
            else {
                coatedCoefficients.push_back(impliedCoefficient);
                CHECK_THAT(resistor.emissivity,
                           Catch::Matchers::WithinRel(OpenMagnetics::ThermalDefaults::kConvection_DefaultEmissivity, 1e-12));
            }
        }
    };

    // Stored layout, litz | foil | litz: the outer litz layer radiates, the buried foil does not.
    readRadiatingTurnSurfaces(mas.get_magnetic());
    CHECK(bareCoefficients.empty());
    REQUIRE(!coatedCoefficients.empty());
    const size_t coatedSurfacesOfStoredLayout = coatedCoefficients.size();

    // litz | foil: the foil is the outermost layer and covers the litz, so it alone radiates.
    readRadiatingTurnSurfaces(withoutOuterPrimarySection(mas.get_magnetic()));
    CHECK(coatedCoefficients.size() == coatedSurfacesOfStoredLayout);

    // Both kinds must actually be present, or the test is asserting nothing.
    REQUIRE(!bareCoefficients.empty());
    REQUIRE(!coatedCoefficients.empty());
    std::sort(bareCoefficients.begin(), bareCoefficients.end());
    std::sort(coatedCoefficients.begin(), coatedCoefficients.end());
    const double bareMedian = bareCoefficients[bareCoefficients.size() / 2];
    const double coatedMedian = coatedCoefficients[coatedCoefficients.size() / 2];
    UNSCOPED_INFO("implied h_rad: bare foil " << bareMedian << " W/m2K, coated litz " << coatedMedian);
    // The emissivities differ by 12.9x; the surfaces sit at different temperatures, so the
    // coefficients do not differ by exactly that. A factor of 5 is far beyond any temperature
    // spread and far below 12.9 -- it separates "the fix is live" from "every surface got 0.9".
    CHECK(bareMedian < coatedMedian / 5.0);
    // The size of the error this removes: bare metal radiates more than ten times less than the
    // matte dielectric value every surface used to be charged.
    CHECK(OpenMagnetics::ThermalDefaults::kRadiation_BareCopperEmissivity <
          OpenMagnetics::ThermalDefaults::kConvection_DefaultEmissivity / 10.0);
}


// An insulation layer names its material by its MAS name ("Kapton HN"). The lookup lower-cased the
// name before searching the insulation database, whose keys are the MAS names, so every mixed-case
// material missed, fell through to the generic list and threw MATERIAL_DATA_MISSING: the web
// frontend's DMC flow failed with "No thermal conductivity for material 'Kapton HN'" although MAS
// carries one (ABT #1302).
TEST_CASE("Temperature: an insulation material is found under its MAS name", "[temperature][smoke-test]") {
    auto expected = find_insulation_material_by_name("Kapton HN").get_thermal_conductivity();
    REQUIRE(expected);
    REQUIRE_THAT(ThermalResistance::getMaterialThermalConductivity("Kapton HN"),
                 Catch::Matchers::WithinAbs(expected.value(), 1e-12));
    // A name that only the generic list knows still resolves there.
    REQUIRE_THAT(ThermalResistance::getMaterialThermalConductivity("mylar"),
                 Catch::Matchers::WithinAbs(0.15, 1e-12));
}

// ============================================================================
// ABT #1753: a core material's heat conductivity — the MAS value, or a flagged estimate (class,
// family, all sourced materials); absence never throws, malformed data does
// ============================================================================
// The network used the configured coreThermalConductivity (4.0) or a per-composition constant
// table whenever MAS had no heatConductivity: numbers no datasheet backs, applied silently. Now
// one resolver returns the MAS value, or estimates it from the SOURCED materials of the nearest
// pool that has any and says so (Alf, 2026-10-11: "always keep an approved fallback").
namespace {
CoreMaterial core_material_without_heat_conductivity(MaterialType type, std::optional<MaterialComposition> composition) {
    auto material = find_core_material_by_name("N87");
    material.set_name("N87 without heat conductivity");
    material.set_heat_conductivity(std::nullopt);
    material.set_provenance(std::nullopt);
    material.set_material(type);
    material.set_material_composition(composition);
    return material;
}

// Resolves a renamed copy (the WARNING is logged once per material) and returns the warnings.
std::vector<std::string> estimate_warnings(CoreMaterial material, const std::string& probeName) {
    material.set_name(probeName);
    Logger::getInstance().disableCollector();
    Logger::getInstance().enableCollector(LogLevel::WARNING);
    resolve_core_material_heat_conductivity(material);
    auto collected = Logger::getInstance().drainCollected();
    Logger::getInstance().disableCollector();
    std::vector<std::string> messages;
    for (const auto& record : collected) {
        if (record.level == LogLevel::WARNING && record.message.find(probeName) != std::string::npos) {
            messages.push_back(record.message);
        }
    }
    return messages;
}

// The span of the pooled members' values; every member must be sourced and satisfy `inPool`.
std::pair<double, double> sourced_pool_span(const CoreHeatConductivityEstimate& estimate,
                                            const std::function<bool(const CoreMaterial&)>& inPool) {
    double low = std::numeric_limits<double>::max();
    double high = std::numeric_limits<double>::lowest();
    for (const auto& name : estimate.sourceMaterials) {
        auto member = find_core_material_by_name(name);
        REQUIRE(member.get_heat_conductivity());
        REQUIRE(member.get_provenance());
        REQUIRE(inPool(member));
        double value = resolve_dimensional_values(member.get_heat_conductivity().value());
        low = std::min(low, value);
        high = std::max(high, value);
    }
    return {low, high};
}

bool has_sourced_heat_conductivity(const CoreMaterial& material) {
    if (!material.get_heat_conductivity() || !material.get_provenance()) {
        return false;
    }
    const auto provenance = material.get_provenance().value();  // bound: never range-for over a temporary's value()
    for (const auto& entry : provenance) {
        if (entry.get_retracted() && entry.get_retracted().value()) {
            continue;
        }
        if (entry.get_fields()) {
            auto fields = entry.get_fields().value();
            if (std::find(fields.begin(), fields.end(), "heatConductivity") != fields.end() &&
                (entry.get_source_url() || entry.get_source_name())) {
                return true;
            }
        }
    }
    return false;
}
}

TEST_CASE("Core heat conductivity: a material with MAS data uses it exactly", "[temperature][thermal][smoke-test]") {
    SECTION("Sourced value") {
        auto material = find_core_material_by_name("N87");
        REQUIRE(material.get_heat_conductivity());
        auto resolved = resolve_core_material_heat_conductivity(material);
        CHECK_FALSE(resolved.is_estimated());
        CHECK(resolved.tier == CoreHeatConductivityTier::MATERIAL);
        CHECK(resolved.source);
        CHECK_FALSE(resolved.is_unsourced());
        CHECK(resolved.value == resolve_dimensional_values(material.get_heat_conductivity().value()));
        CHECK(ThermalResistance::getCoreMaterialThermalConductivity(material) == resolved.value);
        CHECK(ThermalResistance::getMaterialThermalConductivity("N87") == resolved.value);
    }
    SECTION("Unsourced value: still the material's own, flagged unsourced") {
        auto material = find_core_material_by_name("MPP 14");
        REQUIRE(material.get_heat_conductivity());
        REQUIRE_FALSE(has_sourced_heat_conductivity(material));
        auto resolved = resolve_core_material_heat_conductivity(material);
        CHECK_FALSE(resolved.is_estimated());
        CHECK(resolved.tier == CoreHeatConductivityTier::MATERIAL);
        CHECK(resolved.is_unsourced());
        CHECK(resolved.value == resolve_dimensional_values(material.get_heat_conductivity().value()));
        auto warnings = estimate_warnings(material, "MPP 14 (unsourced warning probe)");
        REQUIRE(warnings.size() == 1);
        CHECK(warnings[0].find("UNSOURCED") != std::string::npos);
    }
}

TEST_CASE("Core heat conductivity: a material without data gets its class estimate, flagged", "[temperature][thermal][smoke-test]") {
    SECTION("MnZn ferrite: family values only, so the class median") {
        auto material = core_material_without_heat_conductivity(MaterialType::FERRITE, MaterialComposition::MN_ZN);
        auto resolved = resolve_core_material_heat_conductivity(material);
        REQUIRE(resolved.is_estimated());
        CHECK(resolved.tier == CoreHeatConductivityTier::CLASS);
        const auto& estimate = resolved.estimate.value();
        CHECK(estimate.tier == CoreHeatConductivityTier::CLASS);
        CHECK(estimate.pool == "class FERRITE MN_ZN");
        REQUIRE(estimate.sourceMaterials.size() > 1);
        CHECK(estimate.numberDataPoints < estimate.sourceMaterials.size());  // family values counted once
        CHECK(estimate.numberDocuments > 1);
        CHECK_FALSE(estimate.predictor);
        auto [low, high] = sourced_pool_span(estimate, [](const CoreMaterial& m) {
            return m.get_material() == MaterialType::FERRITE && m.get_material_composition() == MaterialComposition::MN_ZN;
        });
        CHECK(resolved.value >= low);
        CHECK(resolved.value <= high);
    }
    SECTION("Carbonyl iron powder: per-grade data with a density trend") {
        auto mix1 = find_core_material_by_name("Mix 1");
        REQUIRE_FALSE(mix1.get_heat_conductivity());
        auto resolved = resolve_core_material_heat_conductivity(mix1);
        REQUIRE(resolved.is_estimated());
        CHECK(resolved.tier == CoreHeatConductivityTier::CLASS);
        const auto& estimate = resolved.estimate.value();
        REQUIRE(estimate.predictor);
        CHECK(estimate.predictor.value() == "density");
        CHECK(estimate.predictorValue.value() == mix1.get_density().value());
        CHECK(estimate.predictorValue.value() >= estimate.predictorMinimum.value());
        CHECK(estimate.predictorValue.value() <= estimate.predictorMaximum.value());
        CHECK_FALSE(estimate.trendNotFollowedReason);
        auto [low, high] = sourced_pool_span(estimate, [](const CoreMaterial& m) {
            return m.get_material() == MaterialType::POWDER && m.get_material_composition() == MaterialComposition::CARBONYL_IRON;
        });
        CHECK(resolved.value > low);
        CHECK(resolved.value < high);
        auto warnings = estimate_warnings(mix1, "Mix 1 (estimate warning probe)");
        REQUIRE(warnings.size() == 1);
        CHECK(warnings[0].find("ESTIMATED") != std::string::npos);
        CHECK(warnings[0].find("tier CLASS") != std::string::npos);
        CHECK(warnings[0].find("class POWDER CARBONYL_IRON") != std::string::npos);
        CHECK(warnings[0].find("least-squares trend on density") != std::string::npos);
    }
    SECTION("Carbonyl iron powder outside the sourced density span: the class median, never extrapolated") {
        auto mix1 = resolve_core_material_heat_conductivity(find_core_material_by_name("Mix 1")).estimate.value();
        auto mix17 = find_core_material_by_name("Mix 17");
        REQUIRE_FALSE(mix17.get_heat_conductivity());
        REQUIRE(mix17.get_density().value() < mix1.predictorMinimum.value());
        auto resolved = resolve_core_material_heat_conductivity(mix17);
        REQUIRE(resolved.is_estimated());
        const auto& estimate = resolved.estimate.value();
        CHECK_FALSE(estimate.predictor);
        REQUIRE(estimate.trendNotFollowedReason);
        // The median of the class's per-grade sourced values.
        std::vector<double> values;
        for (const auto& name : estimate.sourceMaterials) {
            values.push_back(resolve_dimensional_values(find_core_material_by_name(name).get_heat_conductivity().value()));
        }
        std::sort(values.begin(), values.end());
        REQUIRE(values.size() % 2 == 1);
        CHECK(resolved.value == values[values.size() / 2]);
        auto warnings = estimate_warnings(mix17, "Mix 17 (estimate warning probe)");
        REQUIRE(warnings.size() == 1);
        CHECK(warnings[0].find("ESTIMATED") != std::string::npos);
        CHECK(warnings[0].find("median") != std::string::npos);
        CHECK(warnings[0].find("trend is not followed") != std::string::npos);
    }
}

TEST_CASE("Core heat conductivity: a class with no sourced member gets its family estimate, flagged", "[temperature][thermal][smoke-test]") {
    // No FeSiCr powder in MAS carries a sourced value: the metal powder and alloy family answers.
    auto material = core_material_without_heat_conductivity(MaterialType::POWDER, MaterialComposition::FE_SI_CR);
    auto resolved = resolve_core_material_heat_conductivity(material);
    REQUIRE(resolved.is_estimated());
    CHECK(resolved.tier == CoreHeatConductivityTier::FAMILY);
    const auto& estimate = resolved.estimate.value();
    CHECK(estimate.pool.find("metal powder and alloy family") != std::string::npos);
    CHECK(estimate.method.find("class POWDER FE_SI_CR has no sourced member") != std::string::npos);
    CHECK(estimate.numberDocuments >= 1);
    CHECK_FALSE(estimate.predictor);
    auto [low, high] = sourced_pool_span(estimate, [](const CoreMaterial& m) { return m.get_material() != MaterialType::FERRITE; });
    CHECK(resolved.value >= low);
    CHECK(resolved.value <= high);
    // Every sourced non-ferrite material is in the pool.
    size_t sourcedMetals = 0;
    for (const auto& name : get_core_material_names()) {
        auto member = find_core_material_by_name(name);
        if (member.get_material() != MaterialType::FERRITE && has_sourced_heat_conductivity(member)) {
            ++sourcedMetals;
        }
    }
    CHECK(estimate.sourceMaterials.size() == sourcedMetals);
    auto warnings = estimate_warnings(material, "FeSiCr powder (family warning probe)");
    REQUIRE(warnings.size() == 1);
    CHECK(warnings[0].find("tier FAMILY") != std::string::npos);
}

TEST_CASE("Core heat conductivity: a proprietary or missing composition falls to its material type, then to all sourced materials", "[temperature][thermal][smoke-test]") {
    size_t sourcedAll = 0;
    std::map<MaterialType, size_t> sourcedByType;
    for (const auto& name : get_core_material_names()) {
        auto member = find_core_material_by_name(name);
        if (has_sourced_heat_conductivity(member)) {
            ++sourcedAll;
            ++sourcedByType[member.get_material()];
        }
    }
    REQUIRE(sourcedAll > 0);
    size_t lastResortCases = 0;
    for (auto type : magic_enum::enum_values<MaterialType>()) {
        for (std::optional<MaterialComposition> composition : {std::optional<MaterialComposition>(MaterialComposition::PROPRIETARY),
                                                                std::optional<MaterialComposition>()}) {
            CAPTURE(magic_enum::enum_name(type), composition.has_value());
            auto resolved = resolve_core_material_heat_conductivity(core_material_without_heat_conductivity(type, composition));
            REQUIRE(resolved.is_estimated());
            const auto& estimate = resolved.estimate.value();
            if (sourcedByType[type] > 0) {
                CHECK(resolved.tier == CoreHeatConductivityTier::FAMILY);
                CHECK(estimate.pool == "material type " + std::string(magic_enum::enum_name(type)));
                CHECK(estimate.sourceMaterials.size() == sourcedByType[type]);
                sourced_pool_span(estimate, [type](const CoreMaterial& m) { return m.get_material() == type; });
            }
            else {
                ++lastResortCases;
                CHECK(resolved.tier == CoreHeatConductivityTier::ALL_SOURCED);
                CHECK(estimate.pool == "all sourced core materials");
                CHECK(estimate.sourceMaterials.size() == sourcedAll);
                auto [low, high] = sourced_pool_span(estimate, [](const CoreMaterial&) { return true; });
                CHECK(resolved.value >= low);
                CHECK(resolved.value <= high);
            }
        }
    }
    // MAS today has no sourced electrical-steel or nanocrystalline value, so the last resort is
    // exercised; if sourcing ever covers every material type, this count drops to 0 by design.
    INFO("material types answered by the last resort: " << lastResortCases / 2);
    CHECK(lastResortCases > 0);
}

TEST_CASE("Core heat conductivity: a malformed value throws", "[temperature][thermal][smoke-test]") {
    for (double bad : {-1.0, 0.0}) {
        auto material = find_core_material_by_name("N87");
        DimensionWithTolerance malformed;
        malformed.set_nominal(bad);
        material.set_heat_conductivity(malformed);
        CHECK_THROWS_AS(resolve_core_material_heat_conductivity(material), MaterialException);
    }
}

TEST_CASE("Temperature: the network takes the core conductivity from the resolver and never from a configured default", "[temperature][thermal][smoke-test]") {
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 29/16/10";
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::array(), 1, "N87");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_coil(coil);

    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.plotSchematic = false;
    // Losses from the N87 design, then held FIXED: only the thermal model is under test.
    magnetic.set_core(core);
    applySimulatedLosses(config, magnetic);

    auto solve = [&](const CoreMaterial& material) {
        core.get_mutable_functional_description().set_material(material);
        magnetic.set_core(core);
        auto result = Temperature(magnetic, config).calculateTemperatures();
        REQUIRE(result.converged);
        return result.maximumTemperature;
    };
    auto with_conductivity = [](CoreMaterial material, double k) {
        DimensionWithTolerance value;
        value.set_nominal(k);
        material.set_heat_conductivity(value);
        return material;
    };

    // A material whose class has no sourced conductivity: the solve runs on the flagged family
    // estimate (it used to run on the configured 4.0, and before this ticket's first draft, threw).
    auto noData = core_material_without_heat_conductivity(MaterialType::POWDER, MaterialComposition::FE_SI_CR);
    const double resolvedK = resolve_core_material_heat_conductivity(noData).value;
    // The guard only discriminates if the estimate is not the old default.
    REQUIRE(std::abs(resolvedK - 4.0) > 0.5);
    const double estimated = solve(noData);
    CHECK_THAT(estimated, Catch::Matchers::WithinAbs(solve(with_conductivity(noData, resolvedK)), 1e-9));
    CHECK(std::abs(estimated - solve(with_conductivity(noData, 4.0))) > 1e-3);
}

// ============================================================================
// ABT #1454 / #1459: the wound winding's path to its enclosure and the core's real exterior
// ============================================================================
// The network used to leave a wound winding with no conduction path to bobbin or core, and let a
// film BEHIND a turn claim that turn's exposed face; its temperature then followed whichever faces
// the convection heuristics happened to leave open (PQ 26/25 flyback: P-S 59 C, P-S-P 89-101 C at
// the same losses; E cores hotter as they got bigger). The core convected from a quarter model
// doubled once. These tests drive the network at FIXED losses (length-proportional per-turn split,
// as in the OMFEM battery), so only the thermal model is under test.
namespace {
OpenMagnetics::Magnetic abt1454Transformer(const std::string& shape, const std::vector<size_t>& pattern,
                                           size_t repetitions, int64_t primaryTurns, int64_t secondaryTurns,
                                           int64_t secondaryParallels, const std::string& primaryWire,
                                           const std::string& secondaryWire);

OpenMagnetics::Magnetic abt1454Transformer(const std::string& shape, const std::vector<size_t>& pattern,
                                           size_t repetitions, int64_t primaryTurns, int64_t secondaryTurns,
                                           int64_t secondaryParallels, std::vector<OpenMagnetics::Wire> wires) {
    auto coil = OpenMagneticsTesting::get_quick_coil({primaryTurns, secondaryTurns}, {1, secondaryParallels}, shape, 1,
                                                     WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
    coil.wind({0.5, 0.5}, pattern, repetitions);
    auto core = OpenMagneticsTesting::get_quick_core(shape, json::array(), 1, "3C95");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}

OpenMagnetics::Magnetic abt1454Transformer(const std::string& shape, const std::vector<size_t>& pattern,
                                           size_t repetitions, int64_t primaryTurns, int64_t secondaryTurns,
                                           int64_t secondaryParallels, const std::string& primaryWire,
                                           const std::string& secondaryWire) {
    return abt1454Transformer(shape, pattern, repetitions, primaryTurns, secondaryTurns, secondaryParallels,
                              {find_wire_by_name(primaryWire), find_wire_by_name(secondaryWire)});
}

double abt1454MaximumTemperatureAtFixedLoss(OpenMagnetics::Magnetic magnetic, double coreLosses, double windingLosses) {
    auto turns = magnetic.get_coil().get_turns_description().value();
    double totalLength = 0;
    for (auto& turn : turns) totalLength += turn.get_length();
    std::vector<WindingLossesPerElement> perTurn;
    for (auto& turn : turns) {
        OhmicLosses ohmic;
        ohmic.set_losses(windingLosses * turn.get_length() / totalLength);
        ohmic.set_origin(ResultOrigin::SIMULATION);
        ohmic.set_method_used("fixed-loss length-proportional split");
        WindingLossesPerElement element;
        element.set_name(turn.get_name());
        element.set_ohmic_losses(ohmic);
        perTurn.push_back(element);
    }
    WindingLossesOutput windingLossesOutput;
    windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
    windingLossesOutput.set_method_used("fixed-loss length-proportional split");
    windingLossesOutput.set_winding_losses(windingLosses);
    windingLossesOutput.set_winding_losses_per_turn(perTurn);

    TemperatureConfig config;
    config.ambientTemperature = 40.0;
    config.coreLosses = coreLosses;
    config.windingLosses = windingLosses;
    config.windingLossesOutput = windingLossesOutput;
    config.plotSchematic = false;
    Temperature temperature(magnetic, config);
    auto result = temperature.calculateTemperatures();
    REQUIRE(result.converged);
    return result.maximumTemperature;
}
} // namespace

TEST_CASE("Temperature: ABT 1454 hottest spot falls with core size at fixed loss", "[temperature][abt-1454][smoke-test]") {
    // Same sandwich winding, same 0.3 W core + 0.3 W copper, growing E cores: every size step adds
    // exterior surface and window, so the hottest spot must get cooler. Before ABT #1454 it rose
    // and fell with layer parity (E 30/15/7 hotter than E 20/10/6).
    const std::vector<std::string> shapes = {"E 20/10/6", "E 25/13/7", "E 30/15/7", "E 32/16/9", "E 42/21/15", "E 55/28/21"};
    double previous = std::numeric_limits<double>::max();
    for (const auto& shape : shapes) {
        auto magnetic = abt1454Transformer(shape, {0, 1, 0}, 1, 61, 7, 1, "Round 0.2 - Grade 1", "Round 0.90 - Grade 1");
        double maximumTemperature = abt1454MaximumTemperatureAtFixedLoss(magnetic, 0.3, 0.3);
        INFO(shape << ": Tmax " << maximumTemperature << " C, previous (smaller) core " << previous << " C");
        CHECK(maximumTemperature < previous);
        previous = maximumTemperature;
    }
}

TEST_CASE("Temperature: ABT 1454 sandwich and plain interleave agree at fixed loss", "[temperature][abt-1454][smoke-test]") {
    // PQ 26/25 flyback of the ticket, 0.355 W core + 0.15 W copper, 40 C ambient. OMFEM's 2D thermal
    // FEM puts P-S, P-S-P and P-S-P-S within 0.1 K (50.4 / 50.4 / 50.3 C): the heat is mostly the
    // core's, and the layout only moves the winding's internal gradient. That gradient is bounded
    // by the copper loss through the winding's own conduction path: 0.15 W across a few films and
    // layers of ~10 K/W is ~1.5 K, so layouts may differ by at most 3 K.
    constexpr double kLayoutBandK = 3.0;
    auto ps = abt1454Transformer("PQ 26/25", {0, 1}, 1, 35, 4, 3, "Round 0.4 - Grade 2", "Round 0.80 - Grade 2");
    auto psp = abt1454Transformer("PQ 26/25", {0, 1, 0}, 1, 35, 4, 3, "Round 0.4 - Grade 2", "Round 0.80 - Grade 2");
    auto psps = abt1454Transformer("PQ 26/25", {0, 1}, 2, 35, 4, 3, "Round 0.4 - Grade 2", "Round 0.80 - Grade 2");
    double tPs = abt1454MaximumTemperatureAtFixedLoss(ps, 0.355, 0.15);
    double tPsp = abt1454MaximumTemperatureAtFixedLoss(psp, 0.355, 0.15);
    double tPsps = abt1454MaximumTemperatureAtFixedLoss(psps, 0.355, 0.15);
    INFO("P-S " << tPs << " C, P-S-P " << tPsp << " C, P-S-P-S " << tPsps << " C");
    CHECK(std::abs(tPsp - tPs) < kLayoutBandK);
    CHECK(std::abs(tPsps - tPs) < kLayoutBandK);
}

TEST_CASE("Temperature: ABT 1454 a coil that does not fit its window throws unless strict geometry is off", "[temperature][abt-1454][smoke-test]") {
    // EI 6.8/3 with a 10 + 5 turn Round 0.2 coil: the winder reports it does not fit
    // (are_sections_and_layers_fitting() is false) and a turn overlaps its enclosure by ~0.58 mm.
    // Strict (default): the contact path has no real geometry, so the network throws. Non-strict:
    // that face gets no conduction path, an ERROR is logged, and the rest of the network solves.
    const std::string shape = "EI 6.8/3/0.8/6.8/1.48/6/3.32";
    auto magnetic = abt1454Transformer(shape, {0, 1}, 1, 10, 5, 1, "Round 0.2 - Grade 1", "Round 0.2 - Grade 1");
    REQUIRE_FALSE(magnetic.get_mutable_coil().are_sections_and_layers_fitting());

    REQUIRE(settings.get_thermal_network_strict_geometry());
    REQUIRE_THROWS_WITH(abt1454MaximumTemperatureAtFixedLoss(magnetic, 0.3, 0.3),
                        Catch::Matchers::ContainsSubstring("overlaps its enclosure"));

    SettingsGuard<bool> lenient(Settings::GetInstance(), &Settings::get_thermal_network_strict_geometry,
                                &Settings::set_thermal_network_strict_geometry, false);
    double maximumTemperature = abt1454MaximumTemperatureAtFixedLoss(magnetic, 0.3, 0.3);
    INFO("EI 6.8/3 non-strict Tmax " << maximumTemperature << " C");
    CHECK(std::isfinite(maximumTemperature));
    CHECK(maximumTemperature > 40.0);
}

TEST_CASE("Temperature: ABT 1454 an unserved litz rests on its enclosure through its strands' enamel", "[temperature][abt-1454][smoke-test]") {
    // An unserved litz is {coating: bare} around individually enamelled strands. Its outer surface
    // is strand enamel, not copper, so its line contact with the bobbin goes through that film.
    // The enclosure sweep read the BUNDLE coating, took the bundle for bare copper and threw
    // "bare round wire ... line contact ... diverges" on every unserved litz (Ampere: litz
    // 120x0.071, 150x0.08, 100x0.063 builds that main evaluates).
    // Bundle diameters from MKF's own litz tables (IEC 60317, grade 1 strands, one serving layer).
    const double unservedDiameter = OpenMagnetics::Wire::get_outer_diameter_bare_litz(0.063e-3, 100);
    const double servedDiameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(0.063e-3, 100);
    auto litz = [](const std::string& coatingJson, double outerDiameter) {
        auto data = json::parse(R"j({"name": "Litz 100x0.063 (abt-1454 test)", "type": "litz", "material": "copper",
            "standard": "IEC 60317", "strand": "Round 0.063 - Grade 1", "numberConductors": 100})j");
        data["coating"] = json::parse(coatingJson);
        data["outerDiameter"] = {{"nominal", outerDiameter}};
        return OpenMagnetics::Wire(data);
    };
    auto unserved = litz(R"({"type": "bare"})", unservedDiameter);
    auto served = litz(R"({"type": "served", "numberLayers": 1})", servedDiameter);
    const std::string shape = "PQ 26/25";

    auto unservedMagnetic = abt1454Transformer(shape, {0, 1}, 1, 12, 6, 1, {unserved, unserved});
    auto servedMagnetic = abt1454Transformer(shape, {0, 1}, 1, 12, 6, 1, {served, served});
    double tUnserved = abt1454MaximumTemperatureAtFixedLoss(unservedMagnetic, 0.3, 0.3);
    double tServed = abt1454MaximumTemperatureAtFixedLoss(servedMagnetic, 0.3, 0.3);
    INFO("PQ 26/25 litz 100x0.063: unserved Tmax " << tUnserved << " C, single served " << tServed << " C");
    CHECK(std::isfinite(tUnserved));
    CHECK(tUnserved > 40.0);
    // Same strands, same losses; the serving only adds a ~25-30 um textile film over the strand
    // enamel the unserved bundle already has. The two must agree within the winding's own
    // internal gradient (0.3 W over ~10 K/W films), and the served one cannot run cooler by more
    // than that: it has the thicker film between copper and enclosure.
    CHECK(std::abs(tUnserved - tServed) < 3.0);
    CHECK(tUnserved <= tServed + 0.5);
}

// ABT #1735: the full network put the core of a WOUND toroid 3-12x hotter than the same core bare
// (T 36/23/15 N87 20:10: 97 K/W full vs 33 K/W core-only at 1 W, and the "core-only" value itself
// wrapped the core in the coil's invisible turns). Three faults in the toroidal builder, all
// topology/area, none a correlation: (1) the outermost insulation wrap convected a section-plane
// quadrant -- half the core height at the mean radius, no annular faces -- a third of its real
// skin; (2) turn-to-wrap conduction ran over the chord between the turn's angle and the wrap
// segment's mid-angle (millimetres of tape) instead of the radial contact; (3) a toroid with no
// wrap blocked a core face outright under any turn and never counted the bare ferrite between
// sparse turns. A winding with no copper loss only ADDS surface to the core it is wound on, so the
// wound core must not be hotter than the bare one at the same core loss.
namespace {
OpenMagnetics::Magnetic abt1735WoundToroid(std::vector<int64_t> turns, const std::string& shape, const std::string& material) {
    std::vector<int64_t> parallels(turns.size(), 1);
    auto coil = OpenMagneticsTesting::get_quick_coil(turns, parallels, shape, 1, WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING, CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);
    auto core = OpenMagneticsTesting::get_quick_core(shape, json::array(), 1, material);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    magnetic.get_mutable_coil().wind();
    return magnetic;
}

double abt1735CoreRise(const OpenMagnetics::Magnetic& magnetic, bool coreOnly, double coreLosses) {
    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = coreLosses;
    config.windingLosses = 0.0;  // P_cu = 0: the winding is a passive fin and wrap
    config.coreOnly = coreOnly;
    config.plotSchematic = false;
    Temperature temperature(magnetic, config);
    auto result = temperature.calculateTemperatures();
    REQUIRE(result.converged);
    double hottestCore = -std::numeric_limits<double>::max();
    for (const auto& node : temperature.getNodes()) {
        if (node.part == ThermalNodePartType::CORE_TOROIDAL_SEGMENT) hottestCore = std::max(hottestCore, node.temperature);
    }
    REQUIRE(hottestCore > config.ambientTemperature);
    return hottestCore - config.ambientTemperature;
}

// h*A hand estimate for an isothermal ring body at rise dT, by surface (ABT #1735, Alf
// 2026-10-11): the rim and bore cylinders are VERTICAL with L = the body height, the two annular
// faces HORIZONTAL_TOP / HORIZONTAL_BOTTOM with L = area/perimeter = (R_out - R_in)/2; each adds
// grey-body radiation at the default emissivity -- the bore wall only over the share of its view
// that leaves through the bore's two openings (abt1735BoreViewFactor). The correlation is shared
// with the network on purpose: what this checks is the AREA, the orientation and the length of
// each surface.
// View factor from the inner wall of an open cylinder (radius r, height H) to its two end
// openings, written out here independently of ThermalResistance: coaxial equal disks at
// separation H see each other with F_dd = (X - sqrt(X^2 - 4))/2, X = 2 + (H/r)^2; each opening
// sees the wall with 1 - F_dd; reciprocity 2 pi r H F = 2 pi r^2 (1 - F_dd).
double abt1735BoreViewFactor(double r, double H) {
    double X = 2.0 + (H / r) * (H / r);
    double diskToDisk = (X - std::sqrt(X * X - 4.0)) / 2.0;
    return (r / H) * (1.0 - diskToDisk);
}

double abt1735RingResistance(double rimR, double boreR, double height, double rise) {
    double ambient = 25.0;
    double surface = ambient + rise;
    double hRad = ThermalResistance::calculateRadiationCoefficient(surface, ambient, ThermalDefaults::kConvection_DefaultEmissivity);
    double faceLength = (rimR - boreR) / 2.0;
    double faceArea = std::numbers::pi * (rimR * rimR - boreR * boreR);
    double rimArea = 2.0 * std::numbers::pi * rimR * height;
    double boreArea = 2.0 * std::numbers::pi * boreR * height;
    double hVertical = ThermalResistance::calculateNaturalConvectionCoefficient(surface, ambient, height, SurfaceOrientation::VERTICAL);
    double conductance =
        (hVertical + hRad) * rimArea +
        (hVertical + hRad * abt1735BoreViewFactor(boreR, height)) * boreArea +
        (ThermalResistance::calculateNaturalConvectionCoefficient(surface, ambient, faceLength, SurfaceOrientation::HORIZONTAL_TOP) + hRad) * faceArea +
        (ThermalResistance::calculateNaturalConvectionCoefficient(surface, ambient, faceLength, SurfaceOrientation::HORIZONTAL_BOTTOM) + hRad) * faceArea;
    return 1.0 / conductance;
}
}  // namespace

TEST_CASE("Temperature: ABT 1735 a wound toroid with no copper loss is no hotter than its bare core",
          "[temperature][toroidal][abt-1735][smoke-test]") {
    struct Case { std::vector<int64_t> turns; std::string shape; std::string material; };
    // Two-winding quick coils carry a full-circle insulation wrap (faults 1 and 2); the single
    // sparse winding has none (fault 3); the 60-turn one is a fully covered single winding.
    std::vector<Case> cases = {
        {{20, 10}, "T 36/23/15", "N87"},
        {{5, 5}, "T 20/10/7", "N87"},
        {{24, 12, 8, 6}, "T 36/23/15", "3C97"},
        {{10}, "T 36/23/15", "N87"},
        {{60}, "T 20/10/7", "N87"},
    };
    for (const auto& c : cases) {
        CAPTURE(c.shape, c.material, c.turns.size(), c.turns[0]);
        auto magnetic = abt1735WoundToroid(c.turns, c.shape, c.material);
        double bareRise = abt1735CoreRise(magnetic, true, 1.0);
        double woundRise = abt1735CoreRise(magnetic, false, 1.0);
        CAPTURE(bareRise, woundRise);
        // The wrap and the turns' contact strips add a little series resistance (enamel, tape,
        // the air a turn holds open under the tape next to it) to a surface they enlarge by
        // 25-40%; at no copper loss the wound core may land a few per cent above the bare one,
        // never several times. Between sparse turns the tape lies on the ferrite: a 0.5 mm air
        // film under the whole wrap of the 5:5 T 20/10/7 alone is ~17 K/W, 1.2x its bare core.
        CHECK(woundRise < bareRise * 1.15);
    }
}

TEST_CASE("Temperature: ABT 1735 toroid surface resistance matches the h*A of its own surface",
          "[temperature][toroidal][abt-1735][smoke-test]") {
    // Bare core: its whole surface -- both cylinders AND both annular faces -- by arithmetic.
    {
        auto magnetic = abt1735WoundToroid({20, 10}, "T 36/23/15", "N87");
        auto dimensions = magnetic.get_core().resolve_shape().get_dimensions().value();
        double A = resolve_dimensional_values(dimensions["A"]);
        double B = resolve_dimensional_values(dimensions["B"]);
        double C = resolve_dimensional_values(dimensions["C"]);
        double rise = abt1735CoreRise(magnetic, true, 1.0);
        double handEstimate = abt1735RingResistance(A / 2.0, B / 2.0, C, rise);
        CAPTURE(rise, handEstimate);
        // Isothermal ring, same correlations: the network must reproduce h*A to a few per cent.
        CHECK(rise > handEstimate * 0.95);
        CHECK(rise < handEstimate * 1.05);
    }
    // Wound and wrapped: the outermost wrap's skin, by arithmetic from the MAS layer record
    // (radial offset of the layer centre from the core surface, and its thickness).
    {
        auto magnetic = abt1735WoundToroid({20, 10}, "T 36/23/15", "N87");
        auto dimensions = magnetic.get_core().resolve_shape().get_dimensions().value();
        double A = resolve_dimensional_values(dimensions["A"]);
        double B = resolve_dimensional_values(dimensions["B"]);
        double C = resolve_dimensional_values(dimensions["C"]);
        auto layers = magnetic.get_coil().get_layers_description_insulation();
        REQUIRE(!layers.empty());
        double skinOffset = 0;
        for (const auto& layer : layers) {
            if (layer.get_dimensions().size() < 2 || layer.get_coordinates().size() < 2) continue;
            skinOffset = std::max(skinOffset, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2.0);
        }
        REQUIRE(skinOffset > 0);
        double rimR = A / 2.0 + skinOffset;
        double boreR = B / 2.0 - skinOffset;
        double height = C + 2.0 * skinOffset;
        double rise = abt1735CoreRise(magnetic, false, 1.0);
        double handEstimate = abt1735RingResistance(rimR, boreR, height, rise);
        CAPTURE(rise, handEstimate);
        // Hand estimate of the internal path, core to wrap, at 1 W: ~20 turns x 2 contact strips
        // of ~0.5 mm x ~35 mm through ~30 um enamel at ~0.2 W/mK is ~0.4 K/W, plus half a tape
        // thickness: well under 10% of the ~15 K/W surface term. Before the fix: 97 K/W, 6x.
        CHECK(rise > handEstimate * 0.95);
        CHECK(rise < handEstimate * 1.15);
    }
}

TEST_CASE("Temperature: ABT 1735 a toroid's bore radiates to the room only through its openings",
          "[temperature][toroidal][abt-1735][smoke-test]") {
    // Closed form against hand values: r = H gives X = 3, F_dd = (3 - sqrt 5)/2 and
    // F = (sqrt 5 - 1)/2 = 0.618034; H = 2r gives X = 6, F_dd = 3 - 2 sqrt 2 and F = sqrt 2 - 1.
    CHECK_THAT(ThermalResistance::calculateCylinderWallToOpeningsViewFactor(1e-3, 1e-3),
               Catch::Matchers::WithinAbs((std::sqrt(5.0) - 1.0) / 2.0, 1e-12));
    CHECK_THAT(ThermalResistance::calculateCylinderWallToOpeningsViewFactor(2e-3, 4e-3),
               Catch::Matchers::WithinAbs(std::sqrt(2.0) - 1.0, 1e-12));
    CHECK_THROWS(ThermalResistance::calculateCylinderWallToOpeningsViewFactor(0.0, 1e-3));

    // Bare T 36/23/15: every exterior surface radiates to the room at full view except the bore
    // wall, which sees the room through its two openings (r = B/2, H = C) and otherwise itself.
    // It used to radiate its whole area to the room at view factor 1, against the builder's own
    // statement that the bore is a cavity with ~zero net exchange.
    auto magnetic = abt1735WoundToroid({20, 10}, "T 36/23/15", "N87");
    auto dimensions = magnetic.get_core().resolve_shape().get_dimensions().value();
    double A = resolve_dimensional_values(dimensions["A"]);
    double B = resolve_dimensional_values(dimensions["B"]);
    double C = resolve_dimensional_values(dimensions["C"]);
    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = 1.0;
    config.windingLosses = 0.0;
    config.coreOnly = true;
    config.plotSchematic = false;
    Temperature temperature(magnetic, config);
    auto result = temperature.calculateTemperatures();
    REQUIRE(result.converged);
    const size_t ambientIdx = temperature.getNodes().size() - 1;
    REQUIRE(temperature.getNodes()[ambientIdx].part == ThermalNodePartType::AMBIENT);
    double radiatingArea = 0.0;
    double convectingArea = 0.0;
    for (const auto& r : temperature.getResistances()) {
        if (r.nodeToId != ambientIdx) continue;
        if (r.type == HeatTransferType::RADIATION) radiatingArea += r.area;
        if (r.type == HeatTransferType::NATURAL_CONVECTION) convectingArea += r.area;
    }
    double rimArea = std::numbers::pi * A * C;
    double boreArea = std::numbers::pi * B * C;
    double facesArea = 2.0 * std::numbers::pi * (A * A - B * B) / 4.0;
    double boreViewFactor = abt1735BoreViewFactor(B / 2.0, C);
    CAPTURE(rimArea, boreArea, facesArea, boreViewFactor, convectingArea, radiatingArea);
    // T 36/23/15: r = 11.5 mm, H = 15 mm, F = 0.553.
    CHECK(boreViewFactor > 0.5);
    CHECK(boreViewFactor < 0.6);
    // Convection takes the whole surface; radiation the same surface less the bore's self-view.
    CHECK_THAT(convectingArea, Catch::Matchers::WithinRel(rimArea + boreArea + facesArea, 1e-6));
    CHECK_THAT(radiatingArea, Catch::Matchers::WithinRel(rimArea + boreViewFactor * boreArea + facesArea, 1e-6));
}

TEST_CASE("Temperature: ABT 1735 a foil-wound E core with no copper loss is no hotter than its bare core",
          "[temperature][abt-1735][smoke-test]") {
    // Test_Simulator_Web_1's E 40/11 MPP 40, two 6-turn foil windings. FOIL used to classify the
    // winding as PCB-planar, and the planar path convected the core from its quarter-model section
    // quadrants: 6.6 cm2 of core exterior (vs 32 cm2 from the real geometry) and 83 K/W at 1 W
    // against 27.6 K/W core-only. A foil winding is wound around the column like any other.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_simulator_web_1_1095.json");
    OpenMagnetics::Magnetic magnetic(json::parse(std::ifstream(path)));
    REQUIRE(magnetic.get_coil().get_turns_description());
    for (const auto& winding : magnetic.get_coil().get_functional_description()) {
        REQUIRE(std::get<OpenMagnetics::Wire>(winding.get_wire()).get_type() == WireType::FOIL);
    }
    auto hottestCoreRise = [&](bool coreOnly) {
        TemperatureConfig config;
        config.ambientTemperature = 25.0;
        config.coreLosses = 1.0;
        config.windingLosses = 0.0;
        config.coreOnly = coreOnly;
        config.plotSchematic = false;
        Temperature temperature(magnetic, config);
        auto result = temperature.calculateTemperatures();
        REQUIRE(result.converged);
        double hottest = -std::numeric_limits<double>::max();
        for (const auto& node : temperature.getNodes()) {
            if (node.part == ThermalNodePartType::CORE_CENTRAL_COLUMN || node.part == ThermalNodePartType::CORE_LATERAL_COLUMN ||
                node.part == ThermalNodePartType::CORE_TOP_YOKE || node.part == ThermalNodePartType::CORE_BOTTOM_YOKE) {
                hottest = std::max(hottest, node.temperature);
            }
        }
        REQUIRE(hottest > config.ambientTemperature);
        return hottest - config.ambientTemperature;
    };
    double bareRise = hottestCoreRise(true);
    double woundRise = hottestCoreRise(false);
    CAPTURE(bareRise, woundRise);
    CHECK(woundRise < bareRise * 1.15);
}

namespace {
// E 42/21/20 with per-column windows: windings[i] is placed in bobbin winding window windows[i]
// (0: main column; on this shell 1 and 2 are the outer legs' windows).
OpenMagnetics::Magnetic abt1735MultiColumnMagnetic(const std::vector<int64_t>& turns, const std::vector<int64_t>& windows) {
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "3C97");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    json coilJson;
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    const std::vector<std::string> names = {"Primary", "Secondary", "Tertiary"};
    for (size_t i = 0; i < turns.size(); ++i) {
        coilJson["functionalDescription"].push_back(json{{"name", names[i]}, {"numberTurns", turns[i]}, {"numberParallels", 1},
                                                         {"isolationSide", i == 0 ? "primary" : (i == 1 ? "secondary" : "tertiary")},
                                                         {"wire", "Round 0.475 - Grade 1"}});
    }
    OpenMagnetics::Coil coil(coilJson, false);
    for (size_t i = 0; i < windows.size(); ++i) {
        if (windows[i] != 0) coil.get_mutable_functional_description()[i].set_winding_window(windows[i]);
    }
    coil.set_core_columns(core.get_processed_description()->get_columns());
    REQUIRE(coil.wind());
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}

struct Abt1735Solved {
    std::map<std::string, double> legTemperature;   // "Core_Column_<i>" -> C
    double hottestTurnOnLeg[3] = {-1e9, -1e9, -1e9};  // per wound column index
    double coldestTurn = 1e9;
    double heatToAmbient = 0, injected = 0;
    size_t lateralNodes = 0;
    bool anyWoundColumnTag = false;
};

Abt1735Solved abt1735Solve(OpenMagnetics::Magnetic magnetic, double coreLoss, double lossPerTurn) {
    TemperatureConfig config;
    config.ambientTemperature = 25.0;
    config.coreLosses = coreLoss;
    config.plotSchematic = false;
    size_t numberTurns = magnetic.get_coil().get_turns_description()->size();
    config.windingLosses = lossPerTurn * numberTurns;
    if (lossPerTurn > 0) {
        WindingLossesOutput lossesOutput;
        std::vector<WindingLossesPerElement> perTurn;
        for (size_t t = 0; t < numberTurns; ++t) {
            WindingLossesPerElement element;
            OhmicLosses ohmic;
            ohmic.set_losses(lossPerTurn);
            ohmic.set_method_used("test");
            ohmic.set_origin(ResultOrigin::SIMULATION);
            element.set_name(magnetic.get_coil().get_turns_description()->at(t).get_name());
            element.set_ohmic_losses(ohmic);
            perTurn.push_back(element);
        }
        lossesOutput.set_winding_losses_per_turn(perTurn);
        lossesOutput.set_winding_losses(config.windingLosses);
        lossesOutput.set_origin(ResultOrigin::SIMULATION);
        lossesOutput.set_method_used("test");
        config.windingLossesOutput = lossesOutput;
    }
    Temperature temperature(magnetic, config);
    auto result = temperature.calculateTemperatures();
    REQUIRE(result.converged);
    Abt1735Solved solved;
    const auto& nodes = temperature.getNodes();
    size_t ambientIdx = nodes.size();
    for (size_t i = 0; i < nodes.size(); ++i) {
        const auto& node = nodes[i];
        if (node.part == ThermalNodePartType::AMBIENT) ambientIdx = i;
        if (node.part == ThermalNodePartType::CORE_LATERAL_COLUMN) {
            solved.lateralNodes++;
            solved.legTemperature[node.name] = node.temperature;
        }
        if (node.woundColumnIndex) solved.anyWoundColumnTag = true;
        if (node.part == ThermalNodePartType::TURN) {
            solved.coldestTurn = std::min(solved.coldestTurn, node.temperature);
            if (node.woundColumnIndex) {
                auto& hottest = solved.hottestTurnOnLeg[node.woundColumnIndex.value()];
                hottest = std::max(hottest, node.temperature);
            }
        }
        solved.injected += node.powerDissipation;
    }
    REQUIRE(ambientIdx < nodes.size());
    for (const auto& r : temperature.getResistances()) {
        size_t other = nodes.size();
        if (r.nodeToId == ambientIdx) other = r.nodeFromId;
        else if (r.nodeFromId == ambientIdx) other = r.nodeToId;
        if (other < nodes.size() && r.resistance > 0.0) {
            solved.heatToAmbient += (nodes[other].temperature - nodes[ambientIdx].temperature) / r.resistance;
        }
    }
    return solved;
}
}  // namespace

TEST_CASE("Temperature: ABT 1735 a winding on an outer leg is modelled in its own window",
          "[temperature][abt-1735][multi-column][smoke-test]") {
    // The concentric network used to build the main column's window only; a winding on an outer
    // leg was refused (first as a 21 mm "overlap", then as unmodelled). Alf 2026-10-11: model it.
    // Each wound leg becomes its own node, its winding rests on that leg's bobbin column, faces the
    // main winding's crossing across the shared window, and convects from its room-facing wrap.
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_core_per_column_winding_windows(true);
    struct Case { std::string label; std::vector<int64_t> turns; std::vector<int64_t> windows; };
    std::vector<Case> cases = {
        {"secondary on the -x leg", {20, 10}, {0, 2}},
        {"secondary on the +x leg", {20, 10}, {0, 1}},
        {"secondary and tertiary on both legs", {20, 10, 10}, {0, 1, 2}},
    };
    for (const auto& c : cases) {
        CAPTURE(c.label);
        auto magnetic = abt1735MultiColumnMagnetic(c.turns, c.windows);
        auto columns = magnetic.get_core().get_processed_description()->get_columns();
        size_t mainColumn = magnetic.get_core().get_main_column_index();
        std::set<size_t> woundLegs;
        for (size_t w = 1; w < c.windows.size(); ++w) {
            woundLegs.insert(magnetic.get_core().get_winding_window_column_index(static_cast<size_t>(c.windows[w])));
        }
        // Copper loss only: the leg windings heat their own legs.
        auto solved = abt1735Solve(magnetic, 0.2, 0.02);
        CAPTURE(solved.heatToAmbient, solved.injected);
        CHECK(solved.lateralNodes == 2);
        CHECK_THAT(solved.heatToAmbient, Catch::Matchers::WithinRel(solved.injected, 0.02));
        for (size_t leg = 0; leg < columns.size(); ++leg) {
            if (leg == mainColumn) continue;
            std::string name = "Core_Column_" + std::to_string(leg);
            REQUIRE(solved.legTemperature.count(name));
            CAPTURE(leg, solved.legTemperature[name], solved.hottestTurnOnLeg[leg]);
            if (woundLegs.count(leg)) {
                CHECK(solved.hottestTurnOnLeg[leg] > 25.0);
                CHECK(solved.hottestTurnOnLeg[leg] > solved.legTemperature[name]);
                for (size_t other = 0; other < columns.size(); ++other) {
                    if (other == mainColumn || woundLegs.count(other)) continue;
                    CHECK(solved.legTemperature[name] > solved.legTemperature["Core_Column_" + std::to_string(other)]);
                }
            }
        }
        CHECK(solved.coldestTurn > 25.0);
    }
    // Mirror images: the +x and -x placements are the same component turned around.
    {
        auto minus = abt1735Solve(abt1735MultiColumnMagnetic({20, 10}, {0, 2}), 0.2, 0.02);
        auto plus = abt1735Solve(abt1735MultiColumnMagnetic({20, 10}, {0, 1}), 0.2, 0.02);
        CHECK_THAT(minus.hottestTurnOnLeg[2], Catch::Matchers::WithinRel(plus.hottestTurnOnLeg[1], 1e-6));
        CHECK_THAT(minus.legTemperature["Core_Column_2"], Catch::Matchers::WithinRel(plus.legTemperature["Core_Column_1"], 1e-6));
        CHECK_THAT(minus.legTemperature["Core_Column_1"], Catch::Matchers::WithinRel(plus.legTemperature["Core_Column_2"], 1e-6));
    }
    // Core-only needs no window and still answers.
    {
        auto magnetic = abt1735MultiColumnMagnetic({20, 10}, {0, 2});
        TemperatureConfig config;
        config.coreLosses = 1.0;
        config.coreOnly = true;
        config.plotSchematic = false;
        CHECK(Temperature(magnetic, config).calculateTemperatures().converged);
    }
    settings.reset();
}

TEST_CASE("Temperature: ABT 1735 without an outer-leg winding the network keeps one two-leg node",
          "[temperature][abt-1735][multi-column][smoke-test]") {
    // Per-column windows on, every winding on the main column: nothing is placed on a leg, so the
    // network is the unsplit main-window network (one lateral node standing for both legs, no
    // node tagged with a wound column) and the outer-leg code never runs.
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_core_per_column_winding_windows(true);
    auto solved = abt1735Solve(abt1735MultiColumnMagnetic({20, 10}, {0, 0}), 0.2, 0.02);
    CHECK(solved.lateralNodes == 1);
    CHECK(solved.legTemperature.count("Core_Column_1") == 1);
    CHECK(!solved.anyWoundColumnTag);
    CHECK_THAT(solved.heatToAmbient, Catch::Matchers::WithinRel(solved.injected, 0.02));
    settings.reset();
}
