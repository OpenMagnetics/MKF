#include "physical_models/ThermalEquivalentCircuit.h"
#include "physical_models/WindingLosses.h"
#include "support/Painter.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <fstream>

using namespace MAS;
using namespace OpenMagnetics;

// Helper function to export temperature field SVG
void exportTemperatureFieldSvg(const std::string& testName, 
                                OpenMagnetics::Magnetic magnetic, 
                                const OpenMagnetics::ThermalAnalysisOutput& result) {
    std::filesystem::path outFile = "output";
    std::filesystem::create_directories(outFile);
    outFile.append("thermal_" + testName + ".svg");
    
    // Wind the coil if not already wound
    if (!magnetic.get_coil().get_turns_description()) {
        try {
            magnetic.get_mutable_coil().wind();
        } catch (...) {
            // Coil winding may fail for some test cases, continue without turns
        }
    }
    
    OpenMagnetics::BasicPainter painter(outFile);
    painter.paint_core(magnetic);
    if (magnetic.get_coil().get_turns_description()) {
        painter.paint_coil_turns(magnetic);
    }
    painter.paint_temperature_field(magnetic, result.nodeTemperatures);
    painter.export_svg();
}

// Helper function to export thermal circuit schematic SVG
void exportThermalCircuitSchematic(const std::string& testName,
                                   const OpenMagnetics::ThermalEquivalentCircuit& circuit) {
    std::filesystem::path outFile = "output";
    std::filesystem::create_directories(outFile);
    outFile.append("thermal_schematic_" + testName + ".svg");
    
    OpenMagnetics::BasicPainter painter(outFile);
    auto nodes = circuit.getNodes();
    auto resistances = circuit.getResistances();
    std::string svg = painter.paint_thermal_circuit_schematic(nodes, resistances);
    
    // Write to file
    std::ofstream file(outFile);
    file << svg;
    file.close();
}

namespace {

//=============================================================================
// Unit Tests for Static Calculation Methods
//=============================================================================

TEST_CASE("ThermalEquivalentCircuit: Conduction Resistance Calculation", "[thermal][conduction]") {
    // R = L / (k * A)
    // For a 10mm path through copper (k=385 W/m·K) with 1cm² area:
    // R = 0.01 / (385 * 0.0001) = 0.26 K/W
    
    SECTION("Copper conduction") {
        double length = 0.01;  // 10mm
        double k = 385.0;      // Copper
        double area = 0.0001;  // 1 cm²
        
        double R = ThermalEquivalentCircuit::calculateConductionResistance(length, k, area);
        
        REQUIRE_THAT(R, Catch::Matchers::WithinRel(0.2597, 0.01));
    }
    
    SECTION("Ferrite conduction") {
        double length = 0.02;  // 20mm
        double k = 4.0;        // Ferrite
        double area = 0.001;   // 10 cm²
        
        double R = ThermalEquivalentCircuit::calculateConductionResistance(length, k, area);
        
        // R = 0.02 / (4 * 0.001) = 5 K/W
        REQUIRE_THAT(R, Catch::Matchers::WithinRel(5.0, 0.001));
    }
    
    SECTION("Zero length returns zero resistance") {
        double R = ThermalEquivalentCircuit::calculateConductionResistance(0, 385.0, 0.0001);
        REQUIRE(R == 0.0);
    }
    
    SECTION("Invalid parameters throw exception") {
        REQUIRE_THROWS(ThermalEquivalentCircuit::calculateConductionResistance(0.01, 0, 0.0001));
        REQUIRE_THROWS(ThermalEquivalentCircuit::calculateConductionResistance(0.01, 385.0, 0));
        REQUIRE_THROWS(ThermalEquivalentCircuit::calculateConductionResistance(0.01, -1, 0.0001));
    }
}

TEST_CASE("ThermalEquivalentCircuit: Convection Resistance Calculation", "[thermal][convection]") {
    // R = 1 / (h * A)
    
    SECTION("Basic convection resistance") {
        double h = 10.0;       // Typical natural convection
        double area = 0.01;    // 100 cm²
        
        double R = ThermalEquivalentCircuit::calculateConvectionResistance(h, area);
        
        // R = 1 / (10 * 0.01) = 10 K/W
        REQUIRE_THAT(R, Catch::Matchers::WithinRel(10.0, 0.001));
    }
    
    SECTION("Forced convection lower resistance") {
        double h = 100.0;      // Forced convection
        double area = 0.01;    // 100 cm²
        
        double R = ThermalEquivalentCircuit::calculateConvectionResistance(h, area);
        
        // R = 1 / (100 * 0.01) = 1 K/W
        REQUIRE_THAT(R, Catch::Matchers::WithinRel(1.0, 0.001));
    }
    
    SECTION("Invalid parameters throw exception") {
        REQUIRE_THROWS(ThermalEquivalentCircuit::calculateConvectionResistance(0, 0.01));
        REQUIRE_THROWS(ThermalEquivalentCircuit::calculateConvectionResistance(10.0, 0));
    }
}

TEST_CASE("ThermalEquivalentCircuit: Natural Convection Coefficient", "[thermal][convection]") {
    // Natural convection h typically 5-25 W/(m²·K)
    
    SECTION("Vertical surface, moderate temperature difference") {
        double surfaceTemp = 80.0;   // 80°C
        double ambientTemp = 25.0;   // 25°C
        double charLength = 0.05;    // 5cm characteristic length
        
        double h = ThermalEquivalentCircuit::calculateNaturalConvectionCoefficient(
            surfaceTemp, ambientTemp, charLength, SurfaceOrientation::VERTICAL);
        
        // Should be in typical natural convection range
        REQUIRE(h >= 5.0);
        REQUIRE(h <= 30.0);
    }
    
    SECTION("Horizontal top surface has higher h than bottom") {
        double surfaceTemp = 100.0;
        double ambientTemp = 25.0;
        double charLength = 0.05;
        
        double h_top = ThermalEquivalentCircuit::calculateNaturalConvectionCoefficient(
            surfaceTemp, ambientTemp, charLength, SurfaceOrientation::HORIZONTAL_TOP);
        
        double h_bottom = ThermalEquivalentCircuit::calculateNaturalConvectionCoefficient(
            surfaceTemp, ambientTemp, charLength, SurfaceOrientation::HORIZONTAL_BOTTOM);
        
        // Hot surface facing up has better convection than facing down
        REQUIRE(h_top > h_bottom);
    }
    
    SECTION("Higher temperature difference increases h") {
        double ambientTemp = 25.0;
        double charLength = 0.05;
        
        double h_small_dt = ThermalEquivalentCircuit::calculateNaturalConvectionCoefficient(
            40.0, ambientTemp, charLength, SurfaceOrientation::VERTICAL);
        
        double h_large_dt = ThermalEquivalentCircuit::calculateNaturalConvectionCoefficient(
            100.0, ambientTemp, charLength, SurfaceOrientation::VERTICAL);
        
        REQUIRE(h_large_dt > h_small_dt);
    }
    
    SECTION("Small temperature difference still gives valid h") {
        double h = ThermalEquivalentCircuit::calculateNaturalConvectionCoefficient(
            25.5, 25.0, 0.05, SurfaceOrientation::VERTICAL);
        
        // Should return minimum practical value
        REQUIRE(h >= 5.0);
    }
}

TEST_CASE("ThermalEquivalentCircuit: Forced Convection Coefficient", "[thermal][convection]") {
    // Forced convection h typically 25-250+ W/(m²·K)
    
    SECTION("Low velocity air") {
        double h = ThermalEquivalentCircuit::calculateForcedConvectionCoefficient(
            1.0, 0.05, 25.0);  // 1 m/s, 5cm length, 25°C
        
        REQUIRE(h >= 10.0);
        REQUIRE(h <= 100.0);
    }
    
    SECTION("High velocity air") {
        double h = ThermalEquivalentCircuit::calculateForcedConvectionCoefficient(
            10.0, 0.05, 25.0);  // 10 m/s
        
        REQUIRE(h >= 50.0);
        REQUIRE(h <= 500.0);
    }
    
    SECTION("Higher velocity gives higher h") {
        double h_low = ThermalEquivalentCircuit::calculateForcedConvectionCoefficient(
            1.0, 0.05, 25.0);
        
        double h_high = ThermalEquivalentCircuit::calculateForcedConvectionCoefficient(
            5.0, 0.05, 25.0);
        
        REQUIRE(h_high > h_low);
    }
    
    SECTION("Zero velocity falls back to natural convection") {
        double h = ThermalEquivalentCircuit::calculateForcedConvectionCoefficient(
            0.0, 0.05, 25.0);
        
        // Should give natural convection value
        REQUIRE(h >= 5.0);
    }
}

TEST_CASE("ThermalEquivalentCircuit: Radiation Coefficient", "[thermal][radiation]") {
    // h_rad = ε * σ * (Ts² + Ta²) * (Ts + Ta)
    // At 100°C surface, 25°C ambient, ε=0.9: h_rad ≈ 7-8 W/(m²·K)
    
    SECTION("Typical operating temperatures") {
        double surfaceTemp = 100.0;
        double ambientTemp = 25.0;
        double emissivity = 0.9;
        
        double h_rad = ThermalEquivalentCircuit::calculateRadiationCoefficient(
            surfaceTemp, ambientTemp, emissivity);
        
        // Should be in 5-10 W/(m²·K) range for these temperatures
        REQUIRE(h_rad >= 5.0);
        REQUIRE(h_rad <= 12.0);
    }
    
    SECTION("Emissivity affects coefficient proportionally") {
        double surfaceTemp = 100.0;
        double ambientTemp = 25.0;
        
        double h_high_e = ThermalEquivalentCircuit::calculateRadiationCoefficient(
            surfaceTemp, ambientTemp, 0.9);
        
        double h_low_e = ThermalEquivalentCircuit::calculateRadiationCoefficient(
            surfaceTemp, ambientTemp, 0.5);
        
        // h should be roughly proportional to emissivity
        REQUIRE_THAT(h_high_e / h_low_e, Catch::Matchers::WithinRel(0.9 / 0.5, 0.01));
    }
    
    SECTION("Higher temperature increases radiation coefficient") {
        double ambientTemp = 25.0;
        double emissivity = 0.9;
        
        double h_100 = ThermalEquivalentCircuit::calculateRadiationCoefficient(
            100.0, ambientTemp, emissivity);
        
        double h_150 = ThermalEquivalentCircuit::calculateRadiationCoefficient(
            150.0, ambientTemp, emissivity);
        
        REQUIRE(h_150 > h_100);
    }
}

TEST_CASE("ThermalEquivalentCircuit: Material Thermal Conductivity", "[thermal][materials]") {
    SECTION("Known materials return correct values") {
        // Copper: MAS data interpolates to ~399 at 25°C (from 0°C=401, 127°C=392)
        // Traditional textbook value is 385, MAS uses more detailed data
        REQUIRE_THAT(ThermalEquivalentCircuit::getMaterialThermalConductivity("copper"), 
                     Catch::Matchers::WithinRel(399.0, 0.02));
        
        // Aluminium: MAS data is ~237 at 25°C
        REQUIRE_THAT(ThermalEquivalentCircuit::getMaterialThermalConductivity("aluminium"),
                     Catch::Matchers::WithinRel(237.0, 0.02));
        
        REQUIRE_THAT(ThermalEquivalentCircuit::getMaterialThermalConductivity("ferrite"),
                     Catch::Matchers::WithinRel(4.0, 0.01));
    }
    
    SECTION("Case insensitive lookup") {
        REQUIRE_THAT(ThermalEquivalentCircuit::getMaterialThermalConductivity("COPPER"),
                     Catch::Matchers::WithinRel(399.0, 0.02));
        
        REQUIRE_THAT(ThermalEquivalentCircuit::getMaterialThermalConductivity("Ferrite"),
                     Catch::Matchers::WithinRel(4.0, 0.01));
    }
    
    SECTION("Unknown material returns default") {
        double k = ThermalEquivalentCircuit::getMaterialThermalConductivity("unknown_material");
        REQUIRE(k > 0);  // Should return some default value
    }
}

TEST_CASE("ThermalEquivalentCircuit: Fluid Properties", "[thermal][fluids]") {
    SECTION("Air properties at room temperature") {
        FluidProperties air = FluidProperties::getAirProperties(25.0);
        
        // Density around 1.2 kg/m³
        REQUIRE(air.density > 1.0);
        REQUIRE(air.density < 1.4);
        
        // Thermal conductivity around 0.025 W/(m·K)
        REQUIRE(air.thermalConductivity > 0.020);
        REQUIRE(air.thermalConductivity < 0.030);
        
        // Prandtl number around 0.71 for air
        REQUIRE(air.prandtlNumber > 0.65);
        REQUIRE(air.prandtlNumber < 0.75);
    }
    
    SECTION("Air properties change with temperature") {
        FluidProperties air_cold = FluidProperties::getAirProperties(0.0);
        FluidProperties air_hot = FluidProperties::getAirProperties(100.0);
        
        // Density decreases with temperature (ideal gas)
        REQUIRE(air_cold.density > air_hot.density);
        
        // Thermal conductivity increases with temperature
        REQUIRE(air_hot.thermalConductivity > air_cold.thermalConductivity);
        
        // Viscosity increases with temperature
        REQUIRE(air_hot.dynamicViscosity > air_cold.dynamicViscosity);
    }
}

//=============================================================================
// Integration Tests with Magnetic Components
//=============================================================================

TEST_CASE("ThermalEquivalentCircuit: Configuration", "[thermal][config]") {
    ThermalModelConfiguration config;
    
    SECTION("Default configuration values") {
        REQUIRE(config.ambientTemperature == 25.0);
        REQUIRE(config.convergenceTolerance == 0.1);
        REQUIRE(config.maxIterations == 100);
        REQUIRE(config.includeForcedConvection == false);
        REQUIRE(config.includeRadiation == true);
    }
    
    SECTION("Configuration can be modified") {
        config.ambientTemperature = 40.0;
        config.includeForcedConvection = true;
        config.airVelocity = 2.0;
        
        ThermalEquivalentCircuit circuit(config);
        // Should not throw
    }
}

TEST_CASE("ThermalEquivalentCircuit: Factory Creation", "[thermal][factory]") {
    SECTION("Create equivalent circuit model") {
        auto model = ThermalModel::factory(ThermalModel::ModelType::EQUIVALENT_CIRCUIT);
        REQUIRE(model != nullptr);
    }
    
    SECTION("Create simple equivalent circuit model") {
        auto model = ThermalModel::factory(ThermalModel::ModelType::EQUIVALENT_CIRCUIT_SIMPLE);
        REQUIRE(model != nullptr);
    }
    
    SECTION("Default factory creates equivalent circuit") {
        auto model = ThermalModel::factory();
        REQUIRE(model != nullptr);
    }
}

TEST_CASE("ThermalEquivalentCircuit: Simple ETD Core Temperature", "[thermal][integration]") {
    // Create a simple ETD49 magnetic for testing
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 49",                   // shape
        json::array(),              // no gaps
        {10},                       // 10 turns, single winding
        1,                          // 1 stack
        "N87"                       // material
    );
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    SECTION("Zero losses give ambient temperature") {
        auto result = circuit.calculateTemperatures(magnetic, 0.0, 0.0);
        
        REQUIRE(result.converged);
        REQUIRE_THAT(result.maximumTemperature, 
                     Catch::Matchers::WithinAbs(config.ambientTemperature, 0.5));
    }
    
    SECTION("Core losses increase temperature") {
        double coreLosses = 2.0;  // 2W
        double windingLosses = 0.0;
        
        auto result = circuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
        
        REQUIRE(result.converged);
        REQUIRE(result.maximumTemperature > config.ambientTemperature);
        REQUIRE(result.totalThermalResistance > 0);
        
        exportTemperatureFieldSvg("ETD49_core_losses", magnetic, result);
        
        // Temperature rise should be roughly losses * Rth
        double expectedRise = coreLosses * result.totalThermalResistance;
        double actualRise = result.maximumTemperature - config.ambientTemperature;
        REQUIRE_THAT(actualRise, Catch::Matchers::WithinRel(expectedRise, 0.1));
    }
    
    SECTION("Winding losses increase temperature") {
        double coreLosses = 0.0;
        double windingLosses = 1.5;  // 1.5W
        
        auto result = circuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
        
        REQUIRE(result.converged);
        REQUIRE(result.maximumTemperature > config.ambientTemperature);
        
        exportTemperatureFieldSvg("ETD49_winding_losses", magnetic, result);
    }
    
    SECTION("Combined losses") {
        double coreLosses = 1.0;
        double windingLosses = 1.0;
        
        auto result = circuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
        
        REQUIRE(result.converged);
        
        // Compare with separate losses
        auto coreOnlyResult = circuit.calculateTemperatures(magnetic, coreLosses, 0.0);
        auto windingOnlyResult = circuit.calculateTemperatures(magnetic, 0.0, windingLosses);
        
        // Combined should be higher than either alone (but not necessarily sum due to nonlinearity)
        REQUIRE(result.maximumTemperature >= coreOnlyResult.maximumTemperature);
        REQUIRE(result.maximumTemperature >= windingOnlyResult.maximumTemperature);
    }
}

TEST_CASE("ThermalEquivalentCircuit: Temperature at Point", "[thermal][interpolation]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "E 55/28/21",
        json::array(),
        {20, 10},  // Primary and secondary
        1,
        "N87"
    );
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(magnetic, 2.0, 1.5);
    REQUIRE(result.converged);
    
    exportTemperatureFieldSvg("E55_temperature_at_point", magnetic, result);
    
    SECTION("Can query temperature at any point") {
        std::vector<double> point = {0.0, 0.0, 0.0};  // Center
        double temp = circuit.getTemperatureAtPoint(point);
        
        REQUIRE(temp >= config.ambientTemperature);
        REQUIRE(temp <= result.maximumTemperature + 1.0);  // Allow small tolerance
    }
    
    SECTION("Different points may have different temperatures") {
        std::vector<double> center = {0.0, 0.0, 0.0};
        std::vector<double> edge = {0.05, 0.0, 0.0};  // 5cm from center
        
        double tempCenter = circuit.getTemperatureAtPoint(center);
        double tempEdge = circuit.getTemperatureAtPoint(edge);
        
        // Both should be valid temperatures
        REQUIRE(tempCenter >= config.ambientTemperature);
        REQUIRE(tempEdge >= config.ambientTemperature);
    }
}

TEST_CASE("ThermalEquivalentCircuit: Forced vs Natural Convection", "[thermal][convection]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 39",
        json::array(),
        {15},
        1,
        "N87"
    );
    
    double coreLosses = 1.5;
    double windingLosses = 1.0;
    
    ThermalModelConfiguration naturalConfig;
    naturalConfig.ambientTemperature = 25.0;
    naturalConfig.includeForcedConvection = false;
    
    ThermalModelConfiguration forcedConfig;
    forcedConfig.ambientTemperature = 25.0;
    forcedConfig.includeForcedConvection = true;
    forcedConfig.airVelocity = 3.0;  // 3 m/s
    
    ThermalEquivalentCircuit naturalCircuit(naturalConfig);
    ThermalEquivalentCircuit forcedCircuit(forcedConfig);
    
    auto naturalResult = naturalCircuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    auto forcedResult = forcedCircuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    
    REQUIRE(naturalResult.converged);
    REQUIRE(forcedResult.converged);
    
    exportTemperatureFieldSvg("ETD39_natural_convection", magnetic, naturalResult);
    exportTemperatureFieldSvg("ETD39_forced_convection", magnetic, forcedResult);
    
    // Forced convection should result in lower temperature
    REQUIRE(forcedResult.maximumTemperature < naturalResult.maximumTemperature);
    
    // Forced convection should have lower thermal resistance
    REQUIRE(forcedResult.totalThermalResistance < naturalResult.totalThermalResistance);
}

TEST_CASE("ThermalEquivalentCircuit: Radiation Effect", "[thermal][radiation]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 34",
        json::array(),
        {12},
        1,
        "N87"
    );
    
    double coreLosses = 2.0;
    double windingLosses = 0.5;
    
    ThermalModelConfiguration withRadiation;
    withRadiation.ambientTemperature = 25.0;
    withRadiation.includeRadiation = true;
    
    ThermalModelConfiguration withoutRadiation;
    withoutRadiation.ambientTemperature = 25.0;
    withoutRadiation.includeRadiation = false;
    
    ThermalEquivalentCircuit circuitWithRad(withRadiation);
    ThermalEquivalentCircuit circuitWithoutRad(withoutRadiation);
    
    auto resultWithRad = circuitWithRad.calculateTemperatures(magnetic, coreLosses, windingLosses);
    auto resultWithoutRad = circuitWithoutRad.calculateTemperatures(magnetic, coreLosses, windingLosses);
    
    REQUIRE(resultWithRad.converged);
    REQUIRE(resultWithoutRad.converged);
    
    exportTemperatureFieldSvg("ETD34_with_radiation", magnetic, resultWithRad);
    exportTemperatureFieldSvg("ETD34_without_radiation", magnetic, resultWithoutRad);
    
    // Radiation adds another heat transfer path, so temperature should be lower
    REQUIRE(resultWithRad.maximumTemperature < resultWithoutRad.maximumTemperature);
}

TEST_CASE("ThermalEquivalentCircuit: Ambient Temperature Effect", "[thermal][ambient]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "PQ 35/35",
        json::array(),
        {25},
        1,
        "N87"
    );
    
    double coreLosses = 1.5;
    double windingLosses = 1.0;
    
    ThermalModelConfiguration coldAmbient;
    coldAmbient.ambientTemperature = 20.0;
    
    ThermalModelConfiguration hotAmbient;
    hotAmbient.ambientTemperature = 50.0;
    
    ThermalEquivalentCircuit coldCircuit(coldAmbient);
    ThermalEquivalentCircuit hotCircuit(hotAmbient);
    
    auto coldResult = coldCircuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    auto hotResult = hotCircuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    
    REQUIRE(coldResult.converged);
    REQUIRE(hotResult.converged);
    
    exportTemperatureFieldSvg("PQ35_cold_ambient", magnetic, coldResult);
    exportTemperatureFieldSvg("PQ35_hot_ambient", magnetic, hotResult);
    
    // Temperature rise above ambient should be similar
    double coldRise = coldResult.maximumTemperature - coldAmbient.ambientTemperature;
    double hotRise = hotResult.maximumTemperature - hotAmbient.ambientTemperature;
    
    // Within 20% due to temperature-dependent properties
    REQUIRE_THAT(coldRise, Catch::Matchers::WithinRel(hotRise, 0.20));
    
    // Absolute temperatures should differ by roughly the ambient difference
    double ambientDiff = hotAmbient.ambientTemperature - coldAmbient.ambientTemperature;
    double maxTempDiff = hotResult.maximumTemperature - coldResult.maximumTemperature;
    
    REQUIRE_THAT(maxTempDiff, Catch::Matchers::WithinRel(ambientDiff, 0.25));
}

TEST_CASE("ThermalEquivalentCircuit: Node Information", "[thermal][nodes]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "E 42/21/15",
        json::array(),
        {18, 9},
        1,
        "N87"
    );
    
    ThermalModelConfiguration config;
    config.nodePerCoilLayer = true;
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(magnetic, 1.5, 1.0);
    REQUIRE(result.converged);
    
    exportTemperatureFieldSvg("E42_node_info", magnetic, result);
    
    SECTION("Nodes are created") {
        auto nodes = circuit.getNodes();
        REQUIRE(!nodes.empty());
        
        // Should have ambient node
        bool hasAmbient = false;
        for (const auto& node : nodes) {
            if (node.type == ThermalNodeType::AMBIENT) {
                hasAmbient = true;
                break;
            }
        }
        REQUIRE(hasAmbient);
    }
    
    SECTION("Node temperatures are stored") {
        REQUIRE(!result.nodeTemperatures.empty());
        
        // All temperatures should be >= ambient
        for (const auto& [name, temp] : result.nodeTemperatures) {
            REQUIRE(temp >= config.ambientTemperature - 0.1);
        }
    }
    
    SECTION("Thermal resistances are created") {
        auto resistances = circuit.getResistances();
        REQUIRE(!resistances.empty());
        
        // All resistances should be positive
        for (const auto& res : resistances) {
            REQUIRE(res.resistance > 0);
        }
    }
}

TEST_CASE("ThermalEquivalentCircuit: Convergence", "[thermal][convergence]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 44",
        json::array(),
        {20},
        1,
        "N87"
    );
    
    SECTION("Normal case converges quickly") {
        ThermalModelConfiguration config;
        config.maxIterations = 100;
        config.convergenceTolerance = 0.1;
        
        ThermalEquivalentCircuit circuit(config);
        
        auto result = circuit.calculateTemperatures(magnetic, 2.0, 1.0);
        
        REQUIRE(result.converged);
        REQUIRE(result.iterationsToConverge < 50);  // Should converge well before limit
        
        exportTemperatureFieldSvg("ETD44_convergence", magnetic, result);
    }
    
    SECTION("Tight tolerance may need more iterations") {
        ThermalModelConfiguration config;
        config.maxIterations = 200;
        config.convergenceTolerance = 0.001;  // Very tight
        
        ThermalEquivalentCircuit circuit(config);
        
        auto result = circuit.calculateTemperatures(magnetic, 2.0, 1.0);
        
        // May or may not converge with tight tolerance, but should run
        REQUIRE(result.iterationsToConverge > 0);
    }
}

TEST_CASE("ThermalEquivalentCircuit: Bulk Thermal Resistance", "[thermal][bulk]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 49",
        json::array(),
        {15},
        1,
        "N87"
    );
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    double totalLosses = 3.0;
    auto result = circuit.calculateTemperatures(magnetic, 2.0, 1.0);
    
    REQUIRE(result.converged);
    
    exportTemperatureFieldSvg("ETD49_bulk_resistance", magnetic, result);
    
    double bulkRth = circuit.getBulkThermalResistance();
    
    SECTION("Bulk resistance is consistent with results") {
        // Rth = (Tmax - Tambient) / Ptotal
        double expectedRth = (result.maximumTemperature - config.ambientTemperature) / totalLosses;
        
        REQUIRE_THAT(bulkRth, Catch::Matchers::WithinRel(expectedRth, 0.01));
    }
    
    SECTION("Bulk resistance is in reasonable range for this core size") {
        // ETD49 should have Rth around 5-15 K/W depending on model
        REQUIRE(bulkRth > 2.0);
        REQUIRE(bulkRth < 25.0);
    }
}

//=============================================================================
// Edge Cases and Robustness Tests
//=============================================================================

TEST_CASE("ThermalEquivalentCircuit: Detailed Loss Distribution", "[thermal][detailed]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "E 55/28/21",
        json::array(),
        {25, 12},
        1,
        "N87"
    );
    
    ThermalModelConfiguration config;
    
    ThermalEquivalentCircuit circuit(config);
    
    double coreLosses = 0.8;  // Total core losses
    double windingLosses = 1.2;  // Total winding losses (primary + secondary)
    
    auto result = circuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    
    REQUIRE(result.converged);
    REQUIRE(result.methodUsed == "ThermalEquivalentCircuit");
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("E55_detailed_losses", magnetic, result);
}

TEST_CASE("ThermalEquivalentCircuit: Very High Losses", "[thermal][edge]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 29",
        json::array(),
        {10},
        1,
        "N87"
    );
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    config.maxIterations = 150;
    
    ThermalEquivalentCircuit circuit(config);
    
    // Very high losses that would cause extreme heating
    double coreLosses = 10.0;
    double windingLosses = 5.0;
    
    auto result = circuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    
    // Should still converge (may be unrealistic temperature)
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > 100.0);  // Will be quite hot
    
    exportTemperatureFieldSvg("ETD29_high_losses", magnetic, result);
}

TEST_CASE("ThermalEquivalentCircuit: Very Small Core", "[thermal][edge]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "E 13/7/4",
        json::array(),
        {5},
        1,
        "N87"
    );
    
    ThermalModelConfiguration config;
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(magnetic, 0.2, 0.1);
    
    REQUIRE(result.converged);
    
    exportTemperatureFieldSvg("E13_small_core", magnetic, result);
    
    // Small core should have higher thermal resistance
    double rth = result.totalThermalResistance;
    REQUIRE(rth > 10.0);  // Small cores have high Rth
}

//=============================================================================
// Paper-Based Validation Tests
// References:
// - Van den Bossche & Valchev: "Thermal Modeling of E-type Magnetic Components"
// - Dey et al. 2021: "Lumped Parameter Thermal Network Modelling of Power Transformers"
// - Guillermo Salinas PhD Thesis: Thermal Modelling of High-Frequency Magnetic Components
//=============================================================================

TEST_CASE("ThermalEquivalentCircuit: Comparison with Maniktala formula", "[thermal][validation][maniktala]") {
    // Maniktala empirical formula: Rth = 53 * Ve^(-0.54) K/W, where Ve is in cm³
    // This validates our equivalent circuit gives results in the same ballpark
    // (Note: our model is physics-based, so some deviation is expected)
    
    constexpr double MAX_ERROR = 0.40;  // 40% tolerance (methods are fundamentally different)
    
    // Test with several core sizes with known effective volumes
    std::vector<std::tuple<std::string, double>> cores = {
        {"ETD 29", 5.47},    // Effective volume in cm³
        {"ETD 34", 7.64},
        {"ETD 44", 17.8},
        {"ETD 49", 24.0}
    };
    
    for (const auto& [coreName, Ve_cm3] : cores) {
        DYNAMIC_SECTION("Core: " << coreName << " (Ve=" << Ve_cm3 << " cm³)") {
            // Calculate Maniktala reference
            double Rth_maniktala = 53.0 * std::pow(Ve_cm3, -0.54);
            
            auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
                coreName, json::array(), {12}, 1, "N87");
            
            ThermalModelConfiguration config;
            config.ambientTemperature = 25.0;
            
            ThermalEquivalentCircuit circuit(config);
            
            // Use 2W total losses for test
            auto result = circuit.calculateTemperatures(magnetic, 1.4, 0.6);
            
            REQUIRE(result.converged);
            
            // Replace spaces with underscores for filename
            std::string safeCoreName = coreName;
            std::replace(safeCoreName.begin(), safeCoreName.end(), ' ', '_');
            exportTemperatureFieldSvg("Maniktala_" + safeCoreName, magnetic, result);
            
            double error = std::abs(result.totalThermalResistance - Rth_maniktala) / Rth_maniktala;
            
            if (verboseTests) {
                std::cout << coreName << ": Maniktala=" << Rth_maniktala 
                          << " K/W, Circuit=" << result.totalThermalResistance 
                          << " K/W, Error=" << error * 100 << "%" << std::endl;
            }
            
            REQUIRE(error < MAX_ERROR);
        }
    }
}

TEST_CASE("ThermalEquivalentCircuit: Van den Bossche E42 experimental validation", "[thermal][validation][paper]") {
    // Reference: Van den Bossche & Valchev, "Thermal Modeling of E-type Magnetic Components"
    // EE42 core experimental data - thermal resistance approx 10-14 K/W
    
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "E 42/21/20", json::array(), {15}, 1, "N87");
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    config.includeRadiation = true;
    config.includeForcedConvection = false;
    
    ThermalEquivalentCircuit circuit(config);
    
    // Paper reports Rth ≈ 10-14 K/W for E42 cores at natural convection
    // Test at several power levels
    std::vector<std::pair<double, double>> powerVsExpectedTempRise = {
        {1.0, 12.0},   // 1W -> ~12°C rise
        {2.0, 24.0},   // 2W -> ~24°C rise  
        {3.0, 35.0},   // 3W -> ~35°C rise
        {5.0, 55.0}    // 5W -> ~55°C rise
    };
    
    for (const auto& [power, expectedRise] : powerVsExpectedTempRise) {
        DYNAMIC_SECTION("Power = " << power << " W") {
            auto result = circuit.calculateTemperatures(magnetic, power * 0.65, power * 0.35);
            
            REQUIRE(result.converged);
            
            exportTemperatureFieldSvg("VanDenBossche_E42_" + std::to_string(static_cast<int>(power)) + "W", magnetic, result);
            
            double actualRise = result.maximumTemperature - config.ambientTemperature;
            double error = std::abs(actualRise - expectedRise) / expectedRise;
            
            if (verboseTests) {
                std::cout << "E42 @ " << power << "W: Expected rise=" << expectedRise 
                          << "°C, Actual rise=" << actualRise 
                          << "°C, Error=" << error * 100 << "%" << std::endl;
            }
            
            // Allow 30% error for complex thermal phenomena
            REQUIRE(error < 0.30);
        }
    }
}

TEST_CASE("ThermalEquivalentCircuit: Dey2021 LPTN Temperature Scaling", "[thermal][validation][paper]") {
    // Reference: Dey et al. 2021, "Lumped Parameter Thermal Network Modelling of Power Transformers"
    // Key finding: Temperature rise scales approximately linearly with power at steady state
    // We validate this relationship holds in our model
    
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 49", json::array(), {20}, 1, "N87");
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    // Measure temperature rise at different power levels
    std::vector<double> powers = {1.0, 2.0, 3.0, 5.0};
    std::vector<double> tempRises;
    
    for (double P : powers) {
        auto result = circuit.calculateTemperatures(magnetic, P * 0.7, P * 0.3);
        REQUIRE(result.converged);
        tempRises.push_back(result.maximumTemperature - config.ambientTemperature);
    }
    
    // Check approximate linearity: ΔT / P should be roughly constant
    std::vector<double> ratios;
    for (size_t i = 0; i < powers.size(); ++i) {
        ratios.push_back(tempRises[i] / powers[i]);
    }
    
    // Average ratio (effective thermal resistance)
    double avgRatio = 0;
    for (double r : ratios) avgRatio += r;
    avgRatio /= ratios.size();
    
    // All ratios should be within 20% of average (allowing for temp-dependent convection)
    for (double r : ratios) {
        double deviation = std::abs(r - avgRatio) / avgRatio;
        REQUIRE(deviation < 0.20);
    }
    
    if (verboseTests) {
        std::cout << "Dey2021 linearity test - ΔT/P ratios:" << std::endl;
        for (size_t i = 0; i < powers.size(); ++i) {
            std::cout << "  P=" << powers[i] << "W: ΔT=" << tempRises[i] 
                      << "°C, Rth=" << ratios[i] << " K/W" << std::endl;
        }
        std::cout << "  Average Rth: " << avgRatio << " K/W" << std::endl;
    }
}

TEST_CASE("ThermalEquivalentCircuit: Salinas Thesis Ferrite Conductivity", "[thermal][validation][paper]") {
    // Reference: Guillermo Salinas PhD Thesis
    // Validates that ferrite thermal conductivity (~4 W/m·K) gives reasonable internal gradients
    // Internal temperature gradient in core should be small due to low losses density
    
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 44", json::array(), {15}, 1, "N87");
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    config.coreThermalConductivity = 4.0;  // Ferrite k ≈ 4 W/(m·K)
    
    ThermalEquivalentCircuit circuit(config);
    
    // Moderate losses
    auto result = circuit.calculateTemperatures(magnetic, 3.0, 1.5);
    
    REQUIRE(result.converged);
    
    exportTemperatureFieldSvg("Salinas_ETD44_ferrite", magnetic, result);
    
    // Get core node temperatures from the result
    double maxCoreTemp = 0;
    double minCoreTemp = 1000;
    
    for (const auto& [name, temp] : result.nodeTemperatures) {
        if (name.find("Core_") != std::string::npos && name != "Ambient") {
            maxCoreTemp = std::max(maxCoreTemp, temp);
            minCoreTemp = std::min(minCoreTemp, temp);
        }
    }
    
    // Internal gradient within core should be relatively small (< 10°C typically)
    // due to ferrite's thermal conductivity and distributed losses
    double internalGradient = maxCoreTemp - minCoreTemp;
    
    if (verboseTests) {
        std::cout << "Core internal gradient: " << internalGradient << "°C" << std::endl;
        std::cout << "Max core temp: " << maxCoreTemp << "°C, Min: " << minCoreTemp << "°C" << std::endl;
    }
    
    // Gradient should be meaningful but not excessive
    REQUIRE(internalGradient >= 0);
    REQUIRE(internalGradient < 25.0);  // Less than 25°C internal gradient
}

TEST_CASE("ThermalEquivalentCircuit: PQ Core Validation", "[thermal][validation]") {
    // PQ cores are commonly used in power electronics
    // Validate that PQ cores work and have reasonable thermal characteristics
    
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "PQ 26/25", json::array(), {18}, 1, "N87");
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(magnetic, 2.0, 1.0);
    
    REQUIRE(result.converged);
    
    exportTemperatureFieldSvg("PQ26_validation", magnetic, result);
    
    // PQ 26/25 has Ve ≈ 5.1 cm³
    // Maniktala formula: Rth ≈ 53 * 5.1^(-0.54) ≈ 21.7 K/W
    double Rth_maniktala_approx = 53.0 * std::pow(5.1, -0.54);
    
    // Our calculated value should be in reasonable range
    REQUIRE(result.totalThermalResistance > 5.0);   // Not unreasonably low
    REQUIRE(result.totalThermalResistance < 50.0);  // Not unreasonably high
    
    if (verboseTests) {
        std::cout << "PQ 26/25 thermal resistance: " << result.totalThermalResistance 
                  << " K/W (Maniktala ref: " << Rth_maniktala_approx << " K/W)" << std::endl;
    }
}

TEST_CASE("ThermalEquivalentCircuit: Toroidal Core", "[thermal][validation]") {
    // Test that toroidal cores are handled (they have different geometry)
    // Define a toroidal shape with dimensions A (outer diameter), B (inner diameter), C (height)
    // Approximate dimensions for a 20mm OD toroid: OD=20mm, ID=10mm, H=7mm
    
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["functionalDescription"]["name"] = "ToroidThermalTest";
    coreJson["functionalDescription"]["type"] = "toroidal";
    coreJson["functionalDescription"]["material"] = "N87";
    coreJson["functionalDescription"]["shape"] = json();
    coreJson["functionalDescription"]["shape"]["family"] = "t";
    coreJson["functionalDescription"]["shape"]["type"] = "standard";
    coreJson["functionalDescription"]["shape"]["magneticCircuit"] = "closed";
    coreJson["functionalDescription"]["shape"]["name"] = "T 20/10/7";
    coreJson["functionalDescription"]["shape"]["dimensions"] = json();
    coreJson["functionalDescription"]["shape"]["dimensions"]["A"] = {{"nominal", 0.020}};  // 20mm outer diameter
    coreJson["functionalDescription"]["shape"]["dimensions"]["B"] = {{"nominal", 0.010}};  // 10mm inner diameter
    coreJson["functionalDescription"]["shape"]["dimensions"]["C"] = {{"nominal", 0.007}};  // 7mm height
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    
    Core core(coreJson);
    
    // Create coil with 10 turns
    OpenMagnetics::Coil coil;
    std::vector<Winding> windings;
    Winding winding;
    winding.set_number_turns(10);
    winding.set_number_parallels(1);
    winding.set_name("Primary");
    winding.set_isolation_side(IsolationSide::PRIMARY);
    windings.push_back(winding);
    coil.set_functional_description(windings);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    // Should handle toroidal geometry without crashing
    auto result = circuit.calculateTemperatures(magnetic, 0.5, 0.3);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("Toroid_20mm", magnetic, result);
    
    // Check that result temperatures are stored
    REQUIRE(!result.nodeTemperatures.empty());
    
    // Temperature rise should be reasonable for 0.8W in a small toroid
    double tempRise = result.maximumTemperature - config.ambientTemperature;
    REQUIRE(tempRise > 0.0);
    REQUIRE(tempRise < 100.0);  // Should not be extreme
}

//=============================================================================
// Multi-Winding Tests (3-4 windings)
// Based on patterns from TestWindingLosses.cpp
//=============================================================================

TEST_CASE("ThermalEquivalentCircuit: Three Winding Transformer", "[thermal][validation][multi-winding]") {
    // Test thermal model with a three-winding transformer similar to 
    // Test_Resistance_Matrix_Three_Windings in TestWindingLosses.cpp
    
    std::vector<int64_t> numberTurns = {30, 15, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    std::string shapeName = "PQ 35/35";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00005, 100));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(2e-5, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    // Distribute losses across windings: primary has most, secondary and tertiary have less
    double coreLosses = 2.0;
    double windingLosses = 1.5;  // Total winding losses
    
    auto result = circuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("PQ35_three_windings", magnetic, result);
    
    // Check thermal resistance is reasonable for PQ 35/35 (Ve ~ 8.7 cm³)
    // Maniktala: Rth ≈ 53 * 8.7^(-0.54) ≈ 16 K/W
    REQUIRE(result.totalThermalResistance > 5.0);
    REQUIRE(result.totalThermalResistance < 30.0);
    
    if (verboseTests) {
        std::cout << "Three-winding transformer:" << std::endl;
        std::cout << "  Max temperature: " << result.maximumTemperature << "°C" << std::endl;
        std::cout << "  Thermal resistance: " << result.totalThermalResistance << " K/W" << std::endl;
    }
}

TEST_CASE("ThermalEquivalentCircuit: Four Winding Transformer", "[thermal][validation][multi-winding]") {
    // Test thermal model with a four-winding transformer
    // Common in multi-output power supplies
    
    std::vector<int64_t> numberTurns = {24, 12, 8, 6};  // Primary + 3 secondaries
    std::vector<int64_t> numberParallels = {1, 1, 1, 1};
    std::string shapeName = "ETD 49";

    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00004, 80));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00004, 80));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00004, 80));
    wires.push_back(OpenMagnetics::Wire::create_quick_litz_wire(0.00004, 80));

    auto coil = OpenMagnetics::Coil::create_quick_coil(shapeName, numberTurns, numberParallels, wires);

    std::string coreMaterial = "N87";
    auto gapping = OpenMagnetics::Core::create_ground_gapping(1e-4, 3);
    auto core = OpenMagnetics::Core::create_quick_core(shapeName, coreMaterial, gapping);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    // Moderate losses for a multi-output supply
    double coreLosses = 3.0;
    double windingLosses = 2.5;
    
    auto result = circuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("ETD49_four_windings", magnetic, result);
    
    // ETD49 has Ve ≈ 24 cm³, Maniktala Rth ≈ 53 * 24^(-0.54) ≈ 10 K/W
    REQUIRE(result.totalThermalResistance > 3.0);
    REQUIRE(result.totalThermalResistance < 25.0);
    
    // With 5.5W total losses, expect ~50-80°C temperature rise
    double tempRise = result.maximumTemperature - config.ambientTemperature;
    REQUIRE(tempRise > 20.0);
    REQUIRE(tempRise < 150.0);
    
    if (verboseTests) {
        std::cout << "Four-winding transformer:" << std::endl;
        std::cout << "  Max temperature: " << result.maximumTemperature << "°C" << std::endl;
        std::cout << "  Temperature rise: " << tempRise << "°C" << std::endl;
        std::cout << "  Thermal resistance: " << result.totalThermalResistance << " K/W" << std::endl;
    }
}

//=============================================================================
// Toroidal Core Tests (Additional)
// Based on patterns from TestWindingLossesToroidalCores namespace
//=============================================================================

TEST_CASE("ThermalEquivalentCircuit: Toroidal Core Ten Turns", "[thermal][validation][toroidal]") {
    // Test toroidal core with multiple turns, similar to 
    // Test_Winding_Losses_Ten_Turns_Round_Sinusoidal_Toroidal_Core
    
    std::vector<int64_t> numberTurns({10});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "T 20/10/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,  // interleavingLevel
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

    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(magnetic, 0.3, 0.2);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("Toroid_T20_10_turns", magnetic, result);
    
    // Toroidal core thermal resistance
    REQUIRE(result.totalThermalResistance > 1.0);
    REQUIRE(result.totalThermalResistance < 50.0);
    
    if (verboseTests) {
        std::cout << "Toroidal T 20/10/7 ten turns:" << std::endl;
        std::cout << "  Max temperature: " << result.maximumTemperature << "°C" << std::endl;
        std::cout << "  Thermal resistance: " << result.totalThermalResistance << " K/W" << std::endl;
    }
}

TEST_CASE("ThermalEquivalentCircuit: Larger Toroidal Core Two Windings", "[thermal][validation][toroidal]") {
    // Test larger toroidal core with two windings
    
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

    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(magnetic, 1.0, 0.8);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("Toroid_T36_two_windings", magnetic, result);
    
    // Larger toroid thermal resistance
    REQUIRE(result.totalThermalResistance > 0.5);
    REQUIRE(result.totalThermalResistance < 40.0);
    
    if (verboseTests) {
        std::cout << "Toroidal T 36/23/15 two windings:" << std::endl;
        std::cout << "  Max temperature: " << result.maximumTemperature << "°C" << std::endl;
        std::cout << "  Thermal resistance: " << result.totalThermalResistance << " K/W" << std::endl;
    }
}

//=============================================================================
// Planar Core Tests
// Based on patterns from TestWindingLossesPlanar namespace
//=============================================================================

TEST_CASE("ThermalEquivalentCircuit: Planar Core ER 14.5", "[thermal][validation][planar]") {
    // Test planar core thermal model
    // ER cores are commonly used for planar transformers
    
    std::vector<int64_t> numberTurns = {8, 4};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "ER 14.5/3/7";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    std::string coreMaterial = "3F4";  // High frequency material for planar
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    // Planar cores typically have lower thermal resistance due to flat shape
    auto result = circuit.calculateTemperatures(magnetic, 0.5, 0.3);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("ER14_planar", magnetic, result);
    
    // ER 14.5 is a small planar core (may have high Rth due to small size)
    REQUIRE(result.totalThermalResistance > 10.0);
    REQUIRE(result.totalThermalResistance < 200.0);
    
    if (verboseTests) {
        std::cout << "Planar ER 14.5/3/7:" << std::endl;
        std::cout << "  Max temperature: " << result.maximumTemperature << "°C" << std::endl;
        std::cout << "  Thermal resistance: " << result.totalThermalResistance << " K/W" << std::endl;
    }
}

TEST_CASE("ThermalEquivalentCircuit: Planar Core LP 32", "[thermal][validation][planar]") {
    // Test larger LP planar core
    // LP cores are specifically designed for planar magnetics
    
    std::vector<int64_t> numberTurns = {16, 8};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "LP 32/15/22";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    std::string coreMaterial = "3F4";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(magnetic, 1.5, 1.0);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("LP32_planar", magnetic, result);
    
    // LP 32 is a larger planar core with better thermal dissipation
    REQUIRE(result.totalThermalResistance > 5.0);
    REQUIRE(result.totalThermalResistance < 40.0);
    
    if (verboseTests) {
        std::cout << "Planar LP 32/15/22:" << std::endl;
        std::cout << "  Max temperature: " << result.maximumTemperature << "°C" << std::endl;
        std::cout << "  Thermal resistance: " << result.totalThermalResistance << " K/W" << std::endl;
    }
}

TEST_CASE("ThermalEquivalentCircuit: Planar Core Three Windings", "[thermal][validation][planar][multi-winding]") {
    // Test planar core with three windings (typical for isolated DC-DC with auxiliary)
    
    std::vector<int64_t> numberTurns = {12, 6, 4};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    std::string shapeName = "ER 23/3.6/13";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    std::string coreMaterial = "3F4";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(magnetic, 1.0, 0.8);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("ER23_planar_three_windings", magnetic, result);
    
    REQUIRE(result.totalThermalResistance > 5.0);
    REQUIRE(result.totalThermalResistance < 100.0);
    
    if (verboseTests) {
        std::cout << "Planar ER 23 three windings:" << std::endl;
        std::cout << "  Max temperature: " << result.maximumTemperature << "°C" << std::endl;
        std::cout << "  Thermal resistance: " << result.totalThermalResistance << " K/W" << std::endl;
    }
}

TEST_CASE("ThermalEquivalentCircuit: Per-Turn Model with WindingLosses", "[thermal][validation][per-turn]") {
    // Test per-turn thermal model using actual winding losses per turn
    // This should show temperature variation between turns with different losses
    
    std::vector<int64_t> numberTurns = {10, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "ETD 34/17/11";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    std::string coreMaterial = "3C97";
    auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    // Wind the coil to get turns description
    magnetic.get_mutable_coil().wind();
    
    // Verify coil has turns
    REQUIRE(magnetic.get_coil().get_turns_description().has_value());
    size_t numTurns = magnetic.get_coil().get_turns_description()->size();
    REQUIRE(numTurns > 0);
    
    // Create a mock WindingLossesOutput with per-turn losses
    WindingLossesOutput windingLossesOutput;
    windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
    windingLossesOutput.set_winding_losses(1.5);  // Total winding losses
    
    // Create per-turn losses (varying losses to simulate proximity effects)
    std::vector<WindingLossesPerElement> lossesPerTurn;
    for (size_t i = 0; i < numTurns; ++i) {
        WindingLossesPerElement turnLoss;
        
        // Create ohmic losses with large variation - hotspot in the middle
        double baseLoss = 0.01;  // Base low loss
        if (i >= 4 && i <= 6) {
            baseLoss = 0.5;  // 50x higher loss in the middle turns (hotspot)
        }
        OhmicLosses ohmic;
        ohmic.set_losses(baseLoss);
        ohmic.set_origin(ResultOrigin::SIMULATION);
        turnLoss.set_ohmic_losses(ohmic);
        
        lossesPerTurn.push_back(turnLoss);
    }
    windingLossesOutput.set_winding_losses_per_turn(lossesPerTurn);
    
    // Calculate core losses (estimate based on typical ferrite at 100kHz)
    double coreLosses = 0.9;  // Total core losses in watts
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    config.nodePerCoilTurn = true;  // Enable per-turn modeling
    
    ThermalEquivalentCircuit circuit(config);
    
    auto result = circuit.calculateTemperatures(
        magnetic, coreLosses, windingLossesOutput);
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
    exportTemperatureFieldSvg("ETD34_per_turn", magnetic, result);
    
    // Verify we get individual turn temperatures
    // In per-turn mode, there should be more temperature nodes
    size_t turnNodeCount = 0;
    double maxTurnTemp = 0, minTurnTemp = 1000;
    for (const auto& [name, temp] : result.nodeTemperatures) {
        if (name.find("turn") != std::string::npos) {  // Note: lowercase "turn"
            turnNodeCount++;
            maxTurnTemp = std::max(maxTurnTemp, temp);
            minTurnTemp = std::min(minTurnTemp, temp);
            std::cout << "  " << name << ": " << temp << "°C" << std::endl;
        }
    }
    std::cout << "Turn temperature spread: " << (maxTurnTemp - minTurnTemp) << "°C" << std::endl;
    
    // We should have individual turn nodes
    REQUIRE(turnNodeCount > 0);
    
    if (verboseTests) {
        std::cout << "Per-Turn Model ETD 34:" << std::endl;
        std::cout << "  Max temperature: " << result.maximumTemperature << "°C" << std::endl;
        std::cout << "  Turn nodes: " << turnNodeCount << std::endl;
        std::cout << "  Max turn temp: " << maxTurnTemp << "°C" << std::endl;
        std::cout << "  Min turn temp: " << minTurnTemp << "°C" << std::endl;
        std::cout << "  Turn temperature spread: " << (maxTurnTemp - minTurnTemp) << "°C" << std::endl;
        std::cout << "  Thermal resistance: " << result.totalThermalResistance << " K/W" << std::endl;
    }
}

TEST_CASE("ThermalEquivalentCircuit: Circuit Schematic Export", "[thermal][schematic]") {
    // Create a magnetic component
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 34",
        json::array(),
        {15, 8},  // Two windings
        1,
        "N87"
    );
    
    double coreLosses = 1.5;
    double windingLosses = 1.0;
    
    ThermalModelConfiguration config;
    config.ambientTemperature = 25.0;
    
    ThermalEquivalentCircuit circuit(config);
    
    // Run calculation to build the circuit
    auto result = circuit.calculateTemperatures(magnetic, coreLosses, windingLosses);
    
    REQUIRE(result.converged);
    
    // Export schematic
    exportThermalCircuitSchematic("ETD34_two_windings", circuit);
    
    // Verify the schematic was created
    std::filesystem::path outFile = "output/thermal_schematic_ETD34_two_windings.svg";
    REQUIRE(std::filesystem::exists(outFile));
    
    // Read the file and check it contains expected elements
    std::ifstream file(outFile);
    std::string content((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());
    file.close();
    
    // Check for key SVG elements
    REQUIRE(content.find("<svg") != std::string::npos);
    REQUIRE(content.find("Thermal Equivalent Circuit") != std::string::npos);
    REQUIRE(content.find("°C") != std::string::npos);
    REQUIRE(content.find("K/W") != std::string::npos);
    
    std::cout << "Thermal circuit schematic exported to: " << outFile << std::endl;
}

} // anonymous namespace
