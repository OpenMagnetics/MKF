#include "processors/CircuitSimulatorInterface.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <filesystem>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
auto outputFilePath = std::filesystem::path{__FILE__}.parent_path().append("..").append("output");

TEST_CASE("Test_CircuitSimulatorExporter_Plecs_Simple_Inductor", "[processor][circuit-simulator-exporter][plecs]") {
    std::vector<int64_t> numberTurns = {20};
    std::vector<int64_t> numberParallels = {1};
    std::string shapeName = "E 32/16/9";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    int64_t numberStacks = 1;
    std::string coreMaterial = "N87";
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto plecsFile = outputFilePath;
    plecsFile.append("Test_Plecs_Simple_Inductor.plecs");
    std::filesystem::remove(plecsFile);

    CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::PLECS);
    exporter.export_magnetic_as_subcircuit(magnetic, 100000, 25, plecsFile.string());

    REQUIRE(std::filesystem::exists(plecsFile));

    std::ifstream f(plecsFile);
    std::string content((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    REQUIRE(content.find("Plecs {") != std::string::npos);
    REQUIRE(content.find("MagneticInterface") != std::string::npos);
    REQUIRE(content.find("P_sat") != std::string::npos);
    REQUIRE(content.find("P_air") != std::string::npos);
    REQUIRE(content.find("InitializationCommands") != std::string::npos);
    REQUIRE(content.find("ACVoltageSource") != std::string::npos);
}

TEST_CASE("Test_CircuitSimulatorExporter_Plecs_Transformer", "[processor][circuit-simulator-exporter][plecs]") {
    std::vector<int64_t> numberTurns = {24, 12};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "E 42/21/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    int64_t numberStacks = 1;
    std::string coreMaterial = "N87";
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto plecsFile = outputFilePath;
    plecsFile.append("Test_Plecs_Transformer.plecs");
    std::filesystem::remove(plecsFile);

    CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::PLECS);
    exporter.export_magnetic_as_subcircuit(magnetic, 100000, 25, plecsFile.string());

    REQUIRE(std::filesystem::exists(plecsFile));

    std::ifstream f(plecsFile);
    std::string content((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    // Transformer should have two MagneticInterface components
    size_t firstMagInt = content.find("MagneticInterface");
    REQUIRE(firstMagInt != std::string::npos);
    size_t secondMagInt = content.find("MagneticInterface", firstMagInt + 1);
    REQUIRE(secondMagInt != std::string::npos);

    // Should have a Resistor (load)
    REQUIRE(content.find("Resistor") != std::string::npos);
}

TEST_CASE("Test_CircuitSimulatorExporter_Plecs_Symbol", "[processor][circuit-simulator-exporter][plecs]") {
    std::vector<int64_t> numberTurns = {24, 12};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "E 42/21/15";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    int64_t numberStacks = 1;
    std::string coreMaterial = "N87";
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto plecsFile = outputFilePath;
    plecsFile.append("Test_Plecs_Symbol.plecs");
    std::filesystem::remove(plecsFile);

    CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::PLECS);
    exporter.export_magnetic_as_symbol(magnetic, plecsFile.string());

    REQUIRE(std::filesystem::exists(plecsFile));

    std::ifstream f(plecsFile);
    std::string content((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    // Symbol should have magnetic components but no test circuit
    REQUIRE(content.find("Plecs {") != std::string::npos);
    REQUIRE(content.find("MagneticInterface") != std::string::npos);
    REQUIRE(content.find("P_sat") != std::string::npos);
    REQUIRE(content.find("ACVoltageSource") == std::string::npos);
    REQUIRE(content.find("Scope") == std::string::npos);
}

} // namespace
