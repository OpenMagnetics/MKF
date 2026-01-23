#include "processors/CircuitSimulatorInterface.h"
#include "TestingUtils.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/Sweeper.h"
#include "support/Painter.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
double max_error = 0.01;
auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
bool plot = false;

TEST_CASE("Test_CircuitSimulatorExporter_Simba_Only_Magnetic", "[processor][circuit-simulator-exporter][simba]") {
    std::vector<int64_t> numberTurns = {30, 10, 5, 1};
    std::vector<int64_t> numberParallels = {1, 1, 1, 2};
    std::string shapeName = "PQ 35/35";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);
    coil.get_mutable_functional_description()[3].set_isolation_side(MAS::IsolationSide::PRIMARY);

    int64_t numberStacks = 1;
    std::string coreMaterial = "95";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto flyback_jsimba_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback.jsimba");
    auto jsimbaFile = outputFilePath;
    jsimbaFile.append("./Test_CircuitSimulatorExporter_Simba_Only_Magnetic.jsimba");

    std::filesystem::remove(jsimbaFile);
    CircuitSimulatorExporter().export_magnetic_as_subcircuit(magnetic, 10000, 100, jsimbaFile.string());
    REQUIRE(std::filesystem::exists(jsimbaFile));
}

TEST_CASE("Test_CircuitSimulatorExporter", "[processor][circuit-simulator-exporter][simba]") {
    std::vector<int64_t> numberTurns = {30, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "PQ 35/35";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    int64_t numberStacks = 1;
    std::string coreMaterial = "95";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto flyback_jsimba_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback.jsimba");
    auto jsimbaFile = outputFilePath;
    jsimbaFile.append("./Test_CircuitSimulatorExporter.jsimba");

    std::filesystem::remove(jsimbaFile);
    CircuitSimulatorExporter().export_magnetic_as_subcircuit(magnetic, 10000, 100, jsimbaFile.string(), flyback_jsimba_path);
    REQUIRE(std::filesystem::exists(jsimbaFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Simba_Json_Ur", "[processor][circuit-simulator-exporter][simba]") {
    auto json_path_74 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_circuitsimulatorexporter_simba_json_ur_74.json");
    std::ifstream json_file_74(json_path_74);
    OpenMagnetics::Magnetic magnetic(json::parse(json_file_74));

    auto flyback_jsimba_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback.jsimba");
    auto jsimbaFile = outputFilePath;
    jsimbaFile.append("./Test_CircuitSimulatorExporter_Simba_Json_Ur.jsimba");

    std::filesystem::remove(jsimbaFile);
    CircuitSimulatorExporter().export_magnetic_as_subcircuit(magnetic, 10000, 100, jsimbaFile.string(), flyback_jsimba_path);
    REQUIRE(std::filesystem::exists(jsimbaFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Simba_Json", "[processor][circuit-simulator-exporter][simba]") {
    auto json_path_87 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_circuitsimulatorexporter_simba_json_87.json");
    std::ifstream json_file_87(json_path_87);
    OpenMagnetics::Magnetic magnetic(json::parse(json_file_87));

    auto flyback_jsimba_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback.jsimba");
    auto jsimbaFile = outputFilePath;
    jsimbaFile.append("./Test_CircuitSimulatorExporter_Simba_Json.jsimba");

    std::filesystem::remove(jsimbaFile);
    CircuitSimulatorExporter().export_magnetic_as_subcircuit(magnetic, 10000, 100, jsimbaFile.string(), flyback_jsimba_path);
    REQUIRE(std::filesystem::exists(jsimbaFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Simba_Json_Toroidal_Core", "[processor][circuit-simulator-exporter][simba]") {
    auto json_path_100 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_circuitsimulatorexporter_simba_json_toroidal_core_100.json");
    std::ifstream json_file_100(json_path_100);
    OpenMagnetics::Magnetic magnetic(json::parse(json_file_100));

    auto flyback_jsimba_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback.jsimba");
    auto jsimbaFile = outputFilePath;
    jsimbaFile.append("./Test_CircuitSimulatorExporter_Simba_Json_Toroidal_Core.jsimba");

    std::filesystem::remove(jsimbaFile);
    CircuitSimulatorExporter().export_magnetic_as_subcircuit(magnetic, 10000, 100, jsimbaFile.string(), flyback_jsimba_path);
    REQUIRE(std::filesystem::exists(jsimbaFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Simba_Json_Ep_Core", "[processor][circuit-simulator-exporter][simba]") {
    auto json_path_113 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_circuitsimulatorexporter_simba_json_ep_core_113.json");
    std::ifstream json_file_113(json_path_113);
    OpenMagnetics::Magnetic magnetic(json::parse(json_file_113));

    auto flyback_jsimba_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback.jsimba");
    auto jsimbaFile = outputFilePath;
    jsimbaFile.append("./Test_CircuitSimulatorExporter_Simba_Json_Ep_Core.jsimba");

    std::filesystem::remove(jsimbaFile);
    CircuitSimulatorExporter().export_magnetic_as_subcircuit(magnetic, 10000, 100, jsimbaFile.string(), flyback_jsimba_path);
    REQUIRE(std::filesystem::exists(jsimbaFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Simba_Powder_Core", "[processor][circuit-simulator-exporter][simba]") {
    std::vector<int64_t> numberTurns = {30, 10, 5, 1};
    std::vector<int64_t> numberParallels = {1, 1, 1, 2};
    std::string shapeName = "PQ 35/35";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);
    coil.get_mutable_functional_description()[3].set_isolation_side(MAS::IsolationSide::PRIMARY);

    int64_t numberStacks = 1;
    std::string coreMaterial = "GX 60";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto flyback_jsimba_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback.jsimba");
    auto jsimbaFile = outputFilePath;
    jsimbaFile.append("./Test_CircuitSimulatorExporter_Simba_Only_Magnetic.jsimba");

    std::filesystem::remove(jsimbaFile);
    CircuitSimulatorExporter().export_magnetic_as_subcircuit(magnetic, 10000, 100, jsimbaFile.string());
    REQUIRE(std::filesystem::exists(jsimbaFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Simba_Bug_0_Gap_Length", "[processor][circuit-simulator-exporter][simba]") {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_0_length_gap.json");
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = mas.get_magnetic();

    auto flyback_jsimba_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback.jsimba");
    auto jsimbaFile = outputFilePath;
    jsimbaFile.append("./Test_CircuitSimulatorExporter_Simba_0_Length_Gap.jsimba");

    std::filesystem::remove(jsimbaFile);
    CircuitSimulatorExporter().export_magnetic_as_subcircuit(magnetic, 10000, 100, jsimbaFile.string());
    REQUIRE(std::filesystem::exists(jsimbaFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Ngspice_Only_Magnetic", "[processor][circuit-simulator-exporter][ngspice]") {
    std::vector<int64_t> numberTurns = {30, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "PQ 35/35";
    std::vector<OpenMagnetics::Wire> wires;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    int64_t numberStacks = 1;
    std::string coreMaterial = "95";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto cirFile = outputFilePath;
    cirFile.append("./Test_CircuitSimulatorExporter_Ngspice_Only_Magnetic.cir");

    std::filesystem::remove(cirFile);
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::NGSPICE).export_magnetic_as_subcircuit(magnetic, 10000, 100, cirFile.string());
    REQUIRE(std::filesystem::exists(cirFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Ltspice_Only_Magnetic_Analytical", "[processor][circuit-simulator-exporter][ltspice]") {
    std::vector<int64_t> numberTurns = {30, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "PQ 35/35";
    std::vector<OpenMagnetics::Wire> wires;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    int64_t numberStacks = 1;
    std::string coreMaterial = "95";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto cirFile = outputFilePath;
    cirFile.append("./Custom_component_made_with_OpenMagnetic.cir");
    std::filesystem::remove(cirFile);
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).export_magnetic_as_subcircuit(magnetic, 10000, 100, cirFile.string(), std::nullopt, CircuitSimulatorExporterCurveFittingModes::ANALYTICAL);
    REQUIRE(std::filesystem::exists(cirFile));

    auto asyFile = outputFilePath;
    asyFile.append("./Custom_component_made_with_OpenMagnetic.asy");
    std::filesystem::remove(asyFile);
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).export_magnetic_as_symbol(magnetic, asyFile.string());
    REQUIRE(std::filesystem::exists(asyFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Ltspice_Only_Magnetic_Ladder", "[processor][circuit-simulator-exporter][ltspice]") {
    std::vector<int64_t> numberTurns = {30, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "PQ 35/35";
    std::vector<OpenMagnetics::Wire> wires;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    int64_t numberStacks = 1;
    std::string coreMaterial = "95";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto cirFile = outputFilePath;
    cirFile.append("./Custom_component_made_with_OpenMagnetic.cir");
    std::filesystem::remove(cirFile);
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).export_magnetic_as_subcircuit(magnetic, 10000, 100, cirFile.string());
    REQUIRE(std::filesystem::exists(cirFile));

    auto asyFile = outputFilePath;
    asyFile.append("./Custom_component_made_with_OpenMagnetic.asy");
    std::filesystem::remove(asyFile);
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).export_magnetic_as_symbol(magnetic, asyFile.string());
    REQUIRE(std::filesystem::exists(asyFile));
}

TEST_CASE("Test_CircuitSimulatorExporter_Ac_Resistance_Coefficients_Analytical", "[processor][circuit-simulator-exporter][ltspice]") {
    std::vector<int64_t> numberTurns = {30};
    std::vector<int64_t> numberParallels = {1};
    std::string shapeName = "PQ 35/35";
    std::vector<OpenMagnetics::Wire> wires;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    int64_t numberStacks = 1;
    std::string coreMaterial = "95";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto coefficientsPerWinding = CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).calculate_ac_resistance_coefficients_per_winding(magnetic, 42, CircuitSimulatorExporterCurveFittingModes::ANALYTICAL);

    size_t numberElements = 100;
    size_t windingIndex = 0;
    double startingFrequency = 0.1;
    double endingFrequency = 1000000;

    Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, windingIndex);
    auto frequenciesVector = windingAcResistanceData.get_x_points();
    auto acResistanceVector = windingAcResistanceData.get_y_points();

    double errorAverage = 0;
    for (size_t index = 0; index < acResistanceVector.size(); ++index) {
        std::vector<double> c(coefficientsPerWinding[0].size());
        for (size_t coefficientIndex = 0; coefficientIndex < coefficientsPerWinding[0].size(); ++coefficientIndex) {
            c[coefficientIndex] = coefficientsPerWinding[0][coefficientIndex];
        }
        auto frequency = frequenciesVector[index];
        double modeledAcResistance = CircuitSimulatorExporter::analytical_model(c.data(), frequency);
        double error = fabs(acResistanceVector[index] - modeledAcResistance) / acResistanceVector[index];
        errorAverage += error;
    }

    errorAverage /= acResistanceVector.size();

    // std::cout << "errorAverage: " << errorAverage << std::endl;
    REQUIRE(0.25 > errorAverage);
}

TEST_CASE("Test_CircuitSimulatorExporter_Ac_Resistance_Coefficients_Ladder", "[processor][circuit-simulator-exporter][ltspice]") {
    std::vector<int64_t> numberTurns = {10};
    std::vector<int64_t> numberParallels = {1};
    std::string shapeName = "PQ 35/35";
    std::vector<OpenMagnetics::Wire> wires;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    int64_t numberStacks = 1;
    std::string coreMaterial = "95";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto coefficientsPerWinding = CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).calculate_ac_resistance_coefficients_per_winding(magnetic, 42, CircuitSimulatorExporterCurveFittingModes::LADDER);

    size_t numberElements = 100;
    size_t windingIndex = 0;
    double startingFrequency = 0.1;
    double endingFrequency = 1000000;

    Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, windingIndex);
    auto frequenciesVector = windingAcResistanceData.get_x_points();
    auto acResistanceVector = windingAcResistanceData.get_y_points();

    // for (size_t coefficientIndex = 0; coefficientIndex < coefficientsPerWinding[0].size(); ++coefficientIndex) {
    //     std::cout << coefficientIndex << std::endl;
    //     std::cout << coefficientsPerWinding[0][coefficientIndex] << std::endl;
    // }

    double errorAverage = 0;
    for (size_t index = 0; index < acResistanceVector.size(); ++index) {
        std::vector<double> c(coefficientsPerWinding[0].size());
        for (size_t coefficientIndex = 0; coefficientIndex < coefficientsPerWinding[0].size(); ++coefficientIndex) {
            c[coefficientIndex] = coefficientsPerWinding[0][coefficientIndex];
        }
        auto frequency = frequenciesVector[index];
        double modeledAcResistance = CircuitSimulatorExporter::ladder_model(c.data(), frequency, acResistanceVector[0]);
        double error = fabs(acResistanceVector[index] - modeledAcResistance) / acResistanceVector[index];
        errorAverage += error;
    }

    errorAverage /= acResistanceVector.size();
    // std::cout << "errorAverage: " << errorAverage << std::endl;

    REQUIRE(0.01 > errorAverage);
}

TEST_CASE("Test_CircuitSimulatorExporter_Ac_Resistance_Coefficients_Ladder_Planar", "[processor][circuit-simulator-exporter][ltspice]") {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "ladder_planar.json");
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = mas.get_magnetic();

    auto coefficientsPerWinding = CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).calculate_ac_resistance_coefficients_per_winding(magnetic, 42, CircuitSimulatorExporterCurveFittingModes::LADDER);

    size_t numberElements = 100;
    size_t windingIndex = 0;
    double startingFrequency = 0.1;
    double endingFrequency = 10000000;

    Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, windingIndex);
    auto frequenciesVector = windingAcResistanceData.get_x_points();
    auto acResistanceVector = windingAcResistanceData.get_y_points();

    // for (size_t coefficientIndex = 0; coefficientIndex < coefficientsPerWinding[0].size(); ++coefficientIndex) {
    //     std::cout << coefficientIndex << std::endl;
    //     std::cout << coefficientsPerWinding[0][coefficientIndex] << std::endl;
    // }

    double errorAverage = 0;
    for (size_t index = 0; index < acResistanceVector.size(); ++index) {
        std::vector<double> c(coefficientsPerWinding[0].size());
        for (size_t coefficientIndex = 0; coefficientIndex < coefficientsPerWinding[0].size(); ++coefficientIndex) {
            c[coefficientIndex] = coefficientsPerWinding[0][coefficientIndex];
        }
        auto frequency = frequenciesVector[index];
        double modeledAcResistance = CircuitSimulatorExporter::ladder_model(c.data(), frequency, acResistanceVector[0]);
        double error = fabs(acResistanceVector[index] - modeledAcResistance) / acResistanceVector[index];
        errorAverage += error;
    }

    errorAverage /= acResistanceVector.size();
    // std::cout << "errorAverage: " << errorAverage << std::endl;

    REQUIRE(0.01 > errorAverage);
}

TEST_CASE("Test_CircuitSimulatorExporter_Core_Resistance_Coefficients_Ladder", "[processor][circuit-simulator-exporter][ltspice]") {
    std::vector<int64_t> numberTurns = {10, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "PQ 35/35";
    std::vector<OpenMagnetics::Wire> wires;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName);

    int64_t numberStacks = 1;
    std::string coreMaterial = "3C97";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto coefficients = CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).calculate_core_resistance_coefficients(magnetic);

    size_t numberElements = 20;
    double startingFrequency = 1000;
    double endingFrequency = 300000;

    Curve2D windingCoreResistanceData = Sweeper().sweep_core_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, 25);
    auto frequenciesVector = windingCoreResistanceData.get_x_points();
    auto coreResistanceVector = windingCoreResistanceData.get_y_points();

    {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;

        outFile.append("Test_CircuitSimulatorExporter_Core_Resistance_Coefficients_Ladder_Theory.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(windingCoreResistanceData, true);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));

    }
    // for (size_t coefficientIndex = 0; coefficientIndex < coefficients.size(); ++coefficientIndex) {
    //     std::cout << coefficientIndex << std::endl;
    //     std::cout << coefficients[coefficientIndex] << std::endl;
    // }

    std::vector<double> modeledCoreResistances;
    double errorAverage = 0;
    for (size_t index = 0; index < coreResistanceVector.size(); ++index) {
        std::vector<double> c(coefficients.size());
        for (size_t coefficientIndex = 0; coefficientIndex < coefficients.size(); ++coefficientIndex) {
            c[coefficientIndex] = coefficients[coefficientIndex];
        }
        auto frequency = frequenciesVector[index];
        double modeledCoreResistance = CircuitSimulatorExporter::core_ladder_model(c.data(), frequency, coreResistanceVector[0]);
        modeledCoreResistances.push_back(modeledCoreResistance);
        double error = fabs(coreResistanceVector[index] - modeledCoreResistance) / coreResistanceVector[index];
        errorAverage += error;
    }

    {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_CircuitSimulatorExporter_Core_Resistance_Coefficients_Ladder_Modeled.svg");
        Painter painter(outFile, false, true);
        painter.paint_curve(Curve2D(frequenciesVector, modeledCoreResistances, "meh"));
        painter.export_svg();
    }

    errorAverage /= coreResistanceVector.size();
    // std::cout << "errorAverage: " << errorAverage << std::endl;

    REQUIRE(0.4 > errorAverage);
}

TEST_CASE("Test_CircuitSimulatorExporter_Ltspice_Web_0", "[processor][circuit-simulator-exporter][ltspice]") {

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "bug_dc_resistance_ltspice.json");
    auto mas = OpenMagneticsTesting::mas_loader(path.string());

    auto magnetic = mas.get_magnetic();

    auto cirFile = outputFilePath;
    cirFile.append("./Custom_component_made_with_OpenMagnetic.cir");
    std::filesystem::remove(cirFile);
    auto effectiveResistanceThisWinding_0 = WindingLosses::calculate_effective_resistance_of_winding(magnetic, 0, 0.1, 100);
    auto effectiveResistanceThisWinding_1 = WindingLosses::calculate_effective_resistance_of_winding(magnetic, 1, 0.1, 100);
    auto dcResistancePerWinding = WindingOhmicLosses::calculate_dc_resistance_per_winding(magnetic.get_coil(), 100);
    REQUIRE_THAT(effectiveResistanceThisWinding_0, Catch::Matchers::WithinAbs(dcResistancePerWinding[0], effectiveResistanceThisWinding_0 * max_error));
    REQUIRE_THAT(effectiveResistanceThisWinding_1, Catch::Matchers::WithinAbs(dcResistancePerWinding[1], effectiveResistanceThisWinding_1 * max_error));
}

TEST_CASE("Test_Guess_Periodicity_Simba", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    std::ifstream is(simulation_path);
    std::vector<std::vector<double>> columns;
    std::vector<std::string> column_names;

    if(is.is_open()) {
        std::string line;
        while(getline(is, line)) {
            std::stringstream ss(line);
            std::string token;
            if (column_names.size() == 0) {
                // Getting column names
                while(getline(ss, token, ',')) {
                    column_names.push_back(token);
                    columns.push_back({});
                }
            }
            else {
                size_t currentColumnIndex = 0;
                while(getline(ss, token, ',')) {
                    columns[currentColumnIndex].push_back(stod(token));
                    currentColumnIndex++;
                }

            }
        }
        is.close();
    }
    else {
        throw std::runtime_error("File not found");
    }
    MAS::Waveform waveform;
    waveform.set_time(columns[0]);
    waveform.set_data(columns[1]);
    auto waveformOnePeriod = CircuitSimulationReader().get_one_period(waveform, 100000);
    REQUIRE(128U == waveformOnePeriod.get_data().size());
}

TEST_CASE("Test_Guess_Separator_Commas", "[processor][circuit-simulation-reader]") {
    std::string row = "columns,separated,by,commas";
    REQUIRE(',' == CircuitSimulationReader::guess_separator(row));
}

TEST_CASE("Test_Guess_Separator_Semicolon", "[processor][circuit-simulation-reader]") {
    std::string row = "columns;separated;by;semicolon";
    REQUIRE(';' == CircuitSimulationReader::guess_separator(row));
}

TEST_CASE("Test_Guess_Separator_Tabs", "[processor][circuit-simulation-reader]") {
    std::string row = "columns\tseparated\tby\ttabs";
    REQUIRE('\t' == CircuitSimulationReader::guess_separator(row));
}

TEST_CASE("Test_Guess_Separator_Simba", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    std::ifstream is(simulation_path);
    std::vector<std::vector<double>> columns;
    std::vector<std::string> column_names;

    if(is.is_open()) {
        std::string line;
        while(getline(is, line)) {
            REQUIRE(',' == CircuitSimulationReader::guess_separator(line));
        }
    }
}

TEST_CASE("Test_Guess_Separator_Ltspice", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "ltspice_simulation.txt");

    std::ifstream is(simulation_path);
    std::vector<std::vector<double>> columns;
    std::vector<std::string> column_names;

    if(is.is_open()) {
        std::string line;
        while(getline(is, line)) {
            REQUIRE('\t' == CircuitSimulationReader::guess_separator(line));
        }
    }
}

TEST_CASE("Test_Import_Csv_Rosano_Forward", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "forward_case.csv");

    double frequency = 200000;
    auto reader = CircuitSimulationReader(simulation_path);
    std::vector<std::map<std::string, std::string>> mapColumnNames;
    std::map<std::string, std::string> primaryColumnNames;
    std::map<std::string, std::string> secondaryColumnNames;
    primaryColumnNames["time"] = "Time";
    primaryColumnNames["current"] = "Ipri";
    primaryColumnNames["magnetizingCurrent"] = "Im";
    primaryColumnNames["voltage"] = "Vpri";
    mapColumnNames.push_back(primaryColumnNames);
    secondaryColumnNames["time"] = "Time";
    secondaryColumnNames["current"] = "Isec";
    secondaryColumnNames["voltage"] = "Vsec";
    mapColumnNames.push_back(secondaryColumnNames);
    auto operatingPoint = reader.extract_operating_point(2, frequency, mapColumnNames);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 121e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto primaryMagnetizingCurrent = primaryExcitation.get_magnetizing_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryMagnetizingCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryMagnetizingCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryVoltage.get_waveform().value());
        painter.export_svg();
    }
}

TEST_CASE("Test_Import_Csv_Rosano_Flyback", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback_case.csv");

    double frequency = 200000;
    auto reader = CircuitSimulationReader(simulation_path);
    std::vector<std::map<std::string, std::string>> mapColumnNames;
    std::map<std::string, std::string> primaryColumnNames;
    std::map<std::string, std::string> secondaryColumnNames;
    primaryColumnNames["time"] = "Time";
    primaryColumnNames["current"] = "Ipri";
    primaryColumnNames["magnetizingCurrent"] = "Imag";
    primaryColumnNames["voltage"] = "Vpri";
    mapColumnNames.push_back(primaryColumnNames);
    secondaryColumnNames["time"] = "Time";
    secondaryColumnNames["current"] = "Isec";
    secondaryColumnNames["voltage"] = "Vsec";
    mapColumnNames.push_back(secondaryColumnNames);
    auto operatingPoint = reader.extract_operating_point(2, frequency, mapColumnNames);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 50e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto primaryMagnetizingCurrent = primaryExcitation.get_magnetizing_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryMagnetizingCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryMagnetizingCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryVoltage.get_waveform().value());
        painter.export_svg();
    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());
        painter.export_svg();
    }
}

TEST_CASE("Test_Simba", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    double turnsRatio = 1.0 / 0.3;
    double frequency = 100000;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 220e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto secondaryFrequency = secondaryExcitation.get_frequency();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE(frequency == secondaryFrequency);
    REQUIRE_THAT(2.79694, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 2.79694 * max_error));
    REQUIRE_THAT(primaryCurrent.get_processed().value().get_rms().value() / turnsRatio, Catch::Matchers::WithinAbs(secondaryCurrent.get_processed().value().get_rms().value(), primaryCurrent.get_processed().value().get_rms().value() / turnsRatio * max_error));
    REQUIRE_THAT(13.1204, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 13.1204 * max_error));
    REQUIRE_THAT(primaryVoltage.get_processed().value().get_rms().value() * turnsRatio, Catch::Matchers::WithinAbs(secondaryVoltage.get_processed().value().get_rms().value(), primaryVoltage.get_processed().value().get_rms().value() * turnsRatio * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryVoltage.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());
        painter.export_svg();
    }
}

TEST_CASE("Test_PFC_Only_Current", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "only_pfc_current_waveform.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 110e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryCurrent = primaryExcitation.get_current().value();
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryCurrent.get_waveform().value());
        painter.export_svg();
    }
}

TEST_CASE("Test_Simba_File_Loaded", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    std::string file = "";
    std::string line;
    std::ifstream is(simulation_path);
    if(is.is_open()) {
        while(getline(is, line)) {
            file += line;
            file += "\n";
        }
        is.close();
    }
    else {
        throw std::runtime_error("File not found");
    }

    double turnsRatio = 1.0 / 0.3;
    double frequency = 100000;
    auto reader = CircuitSimulationReader(file);
    auto operatingPoint = reader.extract_operating_point(2, frequency);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 220e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto secondaryFrequency = secondaryExcitation.get_frequency();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE(frequency == secondaryFrequency);
    REQUIRE_THAT(2.79694, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 2.79694 * max_error));
    REQUIRE_THAT(primaryCurrent.get_processed().value().get_rms().value() / turnsRatio, Catch::Matchers::WithinAbs(secondaryCurrent.get_processed().value().get_rms().value(), primaryCurrent.get_processed().value().get_rms().value() / turnsRatio * max_error));
    REQUIRE_THAT(13.1204, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 13.1204 * max_error));
    REQUIRE_THAT(primaryVoltage.get_processed().value().get_rms().value() * turnsRatio, Catch::Matchers::WithinAbs(secondaryVoltage.get_processed().value().get_rms().value(), primaryVoltage.get_processed().value().get_rms().value() * turnsRatio * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryVoltage.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());
        painter.export_svg();
    }
}

TEST_CASE("Test_Ltspice", "[processor][circuit-simulation-reader][ltspice]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "ltspice_simulation.txt");

    double frequency = 372618;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto secondaryFrequency = secondaryExcitation.get_frequency();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE(frequency == secondaryFrequency);
    REQUIRE_THAT(0.0524431, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 0.0524431 * max_error));
    REQUIRE_THAT(0.4, Catch::Matchers::WithinAbs(secondaryCurrent.get_processed().value().get_rms().value(), 0.4 * max_error));
    REQUIRE_THAT(6, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 6 * max_error));
    REQUIRE_THAT(64, Catch::Matchers::WithinAbs(secondaryVoltage.get_processed().value().get_rms().value(), 64 * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryVoltage.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());
        painter.export_svg();
    }
}

TEST_CASE("Test_Plecs", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "plecs_simulation.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(1, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE_THAT(11.3, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 11.3 * max_error));
    REQUIRE_THAT(324, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 324 * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryVoltage.get_waveform().value());
        painter.export_svg();
    }
}

TEST_CASE("Test_Plecs_Missing_Windings", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "wrong_plecs_simulation.csv");

    double frequency = 50;
    {
        auto reader = CircuitSimulationReader(simulation_path);
        auto operatingPoint = reader.extract_operating_point(1, frequency);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

        REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
        REQUIRE(!operatingPoint.get_excitations_per_winding()[0].get_current());
        REQUIRE(!operatingPoint.get_excitations_per_winding()[0].get_voltage());
    }
    {
        auto reader = CircuitSimulationReader(simulation_path);
        std::vector<std::map<std::string, std::string>> mapColumnNames;
        std::map<std::string, std::string> primaryColumnNames;
        primaryColumnNames["time"] = "IHave";
        primaryColumnNames["current"] = "no";
        primaryColumnNames["voltage"] = "idea";
        mapColumnNames.push_back(primaryColumnNames);
        auto operatingPoint = reader.extract_operating_point(1, frequency, mapColumnNames);
        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

        REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
        REQUIRE(operatingPoint.get_excitations_per_winding()[0].get_current());
        REQUIRE(operatingPoint.get_excitations_per_winding()[0].get_voltage());

        auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
        auto primaryFrequency = primaryExcitation.get_frequency();
        auto primaryCurrent = primaryExcitation.get_current().value();
        auto primaryVoltage = primaryExcitation.get_voltage().value();

        REQUIRE(frequency == primaryFrequency);
        REQUIRE_THAT(11.3, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 11.3 * max_error));
        REQUIRE_THAT(324, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 324 * max_error));

        if (plot) {
            auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("primaryCurrent.svg");
            Painter painter(outFile, false, true);
            painter.paint_waveform(primaryCurrent.get_waveform().value());
            painter.export_svg();
        }
        if (plot) {
            auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("primaryVoltage.svg");
            Painter painter(outFile, false, true);
            painter.paint_waveform(primaryVoltage.get_waveform().value());
            painter.export_svg();
        }
    }
}

TEST_CASE("Test_Psim", "[processor][circuit-simulation-reader][psim]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "psim_simulation.csv");

    double frequency = 120000;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 52e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE_THAT(1.25, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 1.25 * max_error));
    REQUIRE_THAT(29.7, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 29.7 * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(primaryVoltage.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());
        painter.export_svg();
    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile, false, true);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());
        painter.export_svg();
    }
}

TEST_CASE("Test_Psim_Harmonics_Size_Error", "[processor][circuit-simulation-reader][psim]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "psim_simulation.csv");

    double frequency = 100000;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 0.0001);

    auto commonHarmonicIndexes = get_main_harmonic_indexes(operatingPoint, 0.05);
    REQUIRE(49U == commonHarmonicIndexes.back());
}

TEST_CASE("Test_Simba_Column_Names", "[processor][circuit-simulation-reader][simba]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    double frequency = 100000;
    auto reader = CircuitSimulationReader(simulation_path);
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time [s]"));
    REQUIRE(!mapColumnNames[0]["current"].compare("TX1 - W2 - Current [A]"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("TX1 - W2 - Voltage [V]"));
    REQUIRE(!mapColumnNames[1]["time"].compare("Time [s]"));
    REQUIRE(!mapColumnNames[1]["current"].compare("TX1 - W5 - Current [A]"));
    REQUIRE(!mapColumnNames[1]["voltage"].compare("TX1 - W5 - Voltage [V]"));
}

TEST_CASE("Test_Ltspice_Column_Names", "[processor][circuit-simulation-reader][ltspice]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "ltspice_simulation.txt");

    double frequency = 372618;
    auto reader = CircuitSimulationReader(simulation_path);
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    REQUIRE(!mapColumnNames[0]["time"].compare("time"));
    REQUIRE(!mapColumnNames[0]["current"].compare("I(L1)"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("V(n001)"));
    REQUIRE(!mapColumnNames[1]["time"].compare("time"));
    REQUIRE(!mapColumnNames[1]["current"].compare("I(L2)"));
    REQUIRE(!mapColumnNames[1]["voltage"].compare("V(n002)"));
}

TEST_CASE("Test_Plecs_Column_Names", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "plecs_simulation.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path);
    auto mapColumnNames = reader.extract_map_column_names(1, frequency);

    REQUIRE(mapColumnNames.size() == 1);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time / s"));
    REQUIRE(!mapColumnNames[0]["current"].compare("L2:Inductor current"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("L2:Inductor voltage"));
}

TEST_CASE("Test_Plecs_Web", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "plecs_simulation.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path); 
    std::string mapColumnNamesString = R"([{"current":"L2:Inductor current","time":"Time / s","voltage":"L2:Inductor voltage"}])";

    std::vector<std::map<std::string, std::string>> mapColumnNames = json::parse(mapColumnNamesString).get<std::vector<std::map<std::string, std::string>>>();


    auto operatingPoint = reader.extract_operating_point(1, frequency, mapColumnNames);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

    REQUIRE(mapColumnNames.size() == 1);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time / s"));
    REQUIRE(!mapColumnNames[0]["current"].compare("L2:Inductor current"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("L2:Inductor voltage"));
}

TEST_CASE("Test_Plecs_Column_Names_Missing_Windings", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "wrong_plecs_simulation.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path);
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    REQUIRE(!mapColumnNames[0]["current"].compare(""));
    REQUIRE(!mapColumnNames[0]["voltage"].compare(""));
    REQUIRE(!mapColumnNames[1]["current"].compare(""));
    REQUIRE(!mapColumnNames[1]["voltage"].compare(""));
}

TEST_CASE("Test_Psim_Column_Names", "[processor][circuit-simulation-reader][psim]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "psim_simulation.csv");

    double frequency = 120000;
    auto reader = CircuitSimulationReader(simulation_path);
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time"));
    REQUIRE(!mapColumnNames[0]["current"].compare("Ipri"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("Vpri"));
    REQUIRE(!mapColumnNames[1]["time"].compare("Time"));
    REQUIRE(!mapColumnNames[1]["current"].compare("Isec"));
    REQUIRE(!mapColumnNames[1]["voltage"].compare("Vsec"));
}

TEST_CASE("Test_Extract_Column_Names_Web_0", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_0.csv");

    std::string file = "";
    std::string line;
    std::ifstream is(simulation_path);
    if(is.is_open()) {
        while(getline(is, line)) {
            file += line;
            file += "\n";
        }
        is.close();
    }
    else {
        throw std::runtime_error("File not found");
    }

    double frequency = 250000;
    auto reader = CircuitSimulationReader(simulation_path);
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time / s"));
    REQUIRE(!mapColumnNames[0]["current"].compare("I_trafo_HV"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("U_trafo_HV"));
    REQUIRE(!mapColumnNames[1]["time"].compare("Time / s"));
    REQUIRE(!mapColumnNames[1]["current"].compare("I_trafo_LV"));
    REQUIRE(!mapColumnNames[1]["voltage"].compare("U_trafo_LV"));
}

TEST_CASE("Test_Import_Csv_Web_0", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_0.csv");

    double frequency = 250000;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 10e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
}

TEST_CASE("Test_Import_Csv_Web_1", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_1.csv");

    double frequency = 919963.201472;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 10e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
}

TEST_CASE("Test_Import_Csv_Web_2", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_2.csv");

    double frequency = 1e6;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 10e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
}

TEST_CASE("Test_Import_Csv_Web_3", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_3.csv");

    double frequency = 50e3;
    auto reader = CircuitSimulationReader(simulation_path);
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 10e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
}

TEST_CASE("Test_Import_Csv_Web_4", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_4.csv");

    double frequency = 239600;
    auto reader = CircuitSimulationReader(simulation_path);
    auto mapColumnNamesJson = json::parse(R"([{"current":"I(L1","time":"time","voltage":"V(Vin_q1_drain"},{"current":"I(L2","time":"time","voltage":"V(q5_drain_q2_drain"}])");
    std::vector<std::map<std::string, std::string>> mapColumnNames = mapColumnNamesJson.get<std::vector<std::map<std::string, std::string>>>();

    auto operatingPoint = reader.extract_operating_point(2, frequency, mapColumnNames);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 10e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
}

TEST_CASE("Test_CircuitSimulatorExporter_Ltspice_LLC_Trafo_First", "[processor][circuit-simulator-exporter][ltspice][llc]") {
    // Load the LLC_trafo_first MAS file that was reported to have LTspice errors
    auto mas_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "LLC_trafo_first.json");
    auto mas = OpenMagneticsTesting::mas_loader(mas_path.string());
    auto magnetic = mas.get_magnetic();

    // Verify the magnetic has 3 windings (Primary, Secondary, Tertiary)
    auto windingsDescription = magnetic.get_coil().get_functional_description();
    REQUIRE(windingsDescription.size() == 3);

    // Export to LTspice subcircuit (.cir file)
    auto cirFile = outputFilePath;
    cirFile.append("./LLC_trafo_first.cir");
    std::filesystem::remove(cirFile);
    
    double frequency = 100000;  // 100 kHz as per the MAS operating point
    double temperature = 100;   // As per the MAS ambient temperature
    
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).export_magnetic_as_subcircuit(
        magnetic, frequency, temperature, cirFile.string(), std::nullopt, CircuitSimulatorExporterCurveFittingModes::LADDER);
    REQUIRE(std::filesystem::exists(cirFile));

    // Export to LTspice symbol (.asy file)
    auto asyFile = outputFilePath;
    asyFile.append("./LLC_trafo_first.asy");
    std::filesystem::remove(asyFile);
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).export_magnetic_as_symbol(magnetic, asyFile.string());
    REQUIRE(std::filesystem::exists(asyFile));

    // Create a simple test circuit that uses the exported magnetic component
    auto testCircuitFile = outputFilePath;
    testCircuitFile.append("./LLC_trafo_first_test.asc");
    std::filesystem::remove(testCircuitFile);

    // Create a minimal LTspice test circuit with the exported component
    // This circuit applies a sine wave to the primary and measures the secondary outputs
    std::ofstream testCircuit(testCircuitFile);
    testCircuit << "Version 4\n";
    testCircuit << "SHEET 1 880 680\n";
    testCircuit << "WIRE 160 160 80 160\n";
    testCircuit << "WIRE 320 160 240 160\n";
    testCircuit << "WIRE 80 240 80 160\n";
    testCircuit << "WIRE 240 240 240 160\n";
    testCircuit << "WIRE 320 240 320 160\n";
    testCircuit << "WIRE 80 320 80 300\n";
    testCircuit << "WIRE 240 320 240 300\n";
    testCircuit << "WIRE 320 320 320 300\n";
    testCircuit << "WIRE 240 320 80 320\n";
    testCircuit << "WIRE 320 320 240 320\n";
    testCircuit << "FLAG 80 320 0\n";
    testCircuit << "SYMBOL voltage 80 140 R0\n";
    testCircuit << "SYMATTR InstName V1\n";
    testCircuit << "SYMATTR Value SINE(0 400 100000)\n";
    testCircuit << "SYMBOL res 224 144 R0\n";
    testCircuit << "SYMATTR InstName R1\n";
    testCircuit << "SYMATTR Value 10\n";
    testCircuit << "SYMBOL res 304 144 R0\n";
    testCircuit << "SYMATTR InstName R2\n";
    testCircuit << "SYMATTR Value 10\n";
    testCircuit << "SYMBOL LLC_trafo_first 160 200 R0\n";
    testCircuit << "SYMATTR InstName X1\n";
    testCircuit << "TEXT 56 344 Left 2 !.tran 100u\n";
    testCircuit << "TEXT 56 376 Left 2 !.lib LLC_trafo_first.cir\n";
    testCircuit.close();
    REQUIRE(std::filesystem::exists(testCircuitFile));

    // Verify the generated .cir file is not empty and has valid content
    std::ifstream cirFileStream(cirFile);
    std::string cirContent((std::istreambuf_iterator<char>(cirFileStream)),
                           std::istreambuf_iterator<char>());
    cirFileStream.close();
    
    // Check that the subcircuit contains expected elements (use lowercase as generated by exporter)
    REQUIRE(cirContent.find(".subckt") != std::string::npos);
    REQUIRE(cirContent.find(".ends") != std::string::npos);
    
    // Check for inductance definitions (Lm for magnetizing inductance)
    REQUIRE(cirContent.find("L") != std::string::npos);
    
    // Verify the symbol file has valid content
    std::ifstream asyFileStream(asyFile);
    std::string asyContent((std::istreambuf_iterator<char>(asyFileStream)),
                           std::istreambuf_iterator<char>());
    asyFileStream.close();
    
    REQUIRE((asyContent.find("SYMDEF") != std::string::npos || asyContent.find("Version") != std::string::npos));
}

TEST_CASE("Test_CircuitSimulatorExporter_Ltspice_LLC_Trafo_First_Analytical", "[processor][circuit-simulator-exporter][ltspice][llc]") {
    // Same test but with ANALYTICAL curve fitting mode
    auto mas_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "LLC_trafo_first.json");
    auto mas = OpenMagneticsTesting::mas_loader(mas_path.string());
    auto magnetic = mas.get_magnetic();

    auto cirFile = outputFilePath;
    cirFile.append("./LLC_trafo_first_analytical.cir");
    std::filesystem::remove(cirFile);
    
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).export_magnetic_as_subcircuit(
        magnetic, 100000, 100, cirFile.string(), std::nullopt, CircuitSimulatorExporterCurveFittingModes::ANALYTICAL);
    REQUIRE(std::filesystem::exists(cirFile));
    
    // Verify the file has content
    std::ifstream cirFileStream(cirFile);
    std::string cirContent((std::istreambuf_iterator<char>(cirFileStream)),
                           std::istreambuf_iterator<char>());
    cirFileStream.close();
    
    REQUIRE(!cirContent.empty());
    REQUIRE(cirContent.find(".subckt") != std::string::npos);
}

TEST_CASE("Test_CircuitSimulatorExporter_Ltspice_LLC_Trafo_Runnable_Netlist", "[processor][circuit-simulator-exporter][ltspice][llc]") {
    // This test creates a complete runnable LTspice netlist that includes both the 
    // exported subcircuit and a test circuit to validate it runs in LTspice
    
    auto mas_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "LLC_trafo_first.json");
    auto mas = OpenMagneticsTesting::mas_loader(mas_path.string());
    auto magnetic = mas.get_magnetic();

    // Export the subcircuit
    auto subcktFile = outputFilePath;
    subcktFile.append("./LLC_trafo_runnable.cir");
    std::filesystem::remove(subcktFile);
    
    CircuitSimulatorExporter(CircuitSimulatorExporterModels::LTSPICE).export_magnetic_as_subcircuit(
        magnetic, 100000, 100, subcktFile.string(), std::nullopt, CircuitSimulatorExporterCurveFittingModes::LADDER);
    REQUIRE(std::filesystem::exists(subcktFile));

    // Read the subcircuit content
    std::ifstream subcktStream(subcktFile);
    std::string subcktContent((std::istreambuf_iterator<char>(subcktStream)),
                              std::istreambuf_iterator<char>());
    subcktStream.close();

    // Create a complete runnable netlist that includes the subcircuit definition
    // and a simple test circuit
    auto runnableFile = outputFilePath;
    runnableFile.append("./LLC_trafo_runnable_test.cir");
    std::filesystem::remove(runnableFile);

    std::ofstream netlist(runnableFile);
    netlist << "* LTspice runnable test netlist for LLC_trafo_first\n";
    netlist << "* Generated by OpenMagnetics MKF test suite\n\n";
    
    // Include the subcircuit definition
    netlist << subcktContent << "\n\n";
    
    // Create a simple test circuit with the transformer
    netlist << "* Test circuit\n";
    netlist << "V1 in 0 SINE(0 400 100000)  ; 400V peak, 100kHz\n";
    netlist << "R_in in P1_pos 0.1  ; Small input resistance to limit current\n";
    netlist << "X1 P1_pos P1_neg P2_pos P2_neg P3_pos P3_neg My_custom_magnetic\n";
    netlist << "R_load1 P2_pos P2_neg 10  ; Load on secondary\n";
    netlist << "R_load2 P3_pos P3_neg 10  ; Load on tertiary\n";
    netlist << "Rgnd P1_neg 0 0.001  ; Connect primary negative to ground\n";
    netlist << "\n";
    netlist << ".tran 0 100u 0 10n  ; Transient analysis for 100us with 10ns max step\n";
    netlist << ".options plotwinsize=0\n";
    netlist << ".end\n";
    netlist.close();

    REQUIRE(std::filesystem::exists(runnableFile));
    
    // Verify the netlist structure is valid
    std::ifstream netlistStream(runnableFile);
    std::string netlistContent((std::istreambuf_iterator<char>(netlistStream)),
                               std::istreambuf_iterator<char>());
    netlistStream.close();
    
    // Check for required components
    REQUIRE(netlistContent.find(".subckt") != std::string::npos);
    REQUIRE(netlistContent.find(".ends") != std::string::npos);
    REQUIRE(netlistContent.find(".tran") != std::string::npos);
    REQUIRE(netlistContent.find(".end") != std::string::npos);
    REQUIRE(netlistContent.find("X1") != std::string::npos);  // Instance of the subcircuit
    
    // Check that mutual coupling K statements have unique names
    size_t k2_pos = netlistContent.find("K2 Lmag_1 Lmag_2");
    size_t k3_pos = netlistContent.find("K3 Lmag_1 Lmag_3");
    REQUIRE(k2_pos != std::string::npos);
    REQUIRE(k3_pos != std::string::npos);
}

}  // namespace
