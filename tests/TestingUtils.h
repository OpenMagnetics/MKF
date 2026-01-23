#pragma once

#include "Constants.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "constructive_models/Wire.h"
#include <magic_enum.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <source_location>
#include <vector>

extern bool verboseTests;

using namespace MAS;
using namespace OpenMagnetics;

namespace OpenMagneticsTesting {

// Cross-platform helper to get test data path from source file location
// Usage: get_test_data_path(std::source_location::current(), "myfile.json")
inline std::filesystem::path get_test_data_path(const std::source_location& loc, const std::string& filename) {
    return std::filesystem::path(loc.file_name()).parent_path() / "testData" / filename;
}

// Overload that returns just the testData directory
inline std::filesystem::path get_test_data_dir(const std::source_location& loc) {
    return std::filesystem::path(loc.file_name()).parent_path() / "testData";
}
OpenMagnetics::Coil get_quick_coil(std::vector<int64_t> numberTurns,
                                          std::vector<int64_t> numberParallels,
                                          std::string shapeName,
                                          uint8_t interleavingLevel = 1,
                                          MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING,
                                          MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING,
                                          MAS::CoilAlignment turnsAlignment = MAS::CoilAlignment::CENTERED,
                                          MAS::CoilAlignment sectionsAlignment = MAS::CoilAlignment::CENTERED,
                                          std::vector<OpenMagnetics::Wire> wires = std::vector<OpenMagnetics::Wire>({}),
                                          bool useBobbin = true,
                                          int numberStacks = 1);

OpenMagnetics::Coil get_quick_coil(std::vector<int64_t> numberTurns,
                                          std::vector<int64_t> numberParallels,
                                          double bobbinHeight,
                                          double bobbinWidth,
                                          std::vector<double> bobbinCenterCoodinates,
                                          uint8_t interleavingLevel = 1,
                                          MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING,
                                          MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING,
                                          MAS::CoilAlignment turnsAlignment = MAS::CoilAlignment::CENTERED,
                                          MAS::CoilAlignment sectionsAlignment = MAS::CoilAlignment::CENTERED,
                                          std::vector<OpenMagnetics::Wire> wires = std::vector<OpenMagnetics::Wire>({}));

OpenMagnetics::Coil get_quick_coil_no_compact(std::vector<int64_t> numberTurns,
                                                     std::vector<int64_t> numberParallels,
                                                     double bobbinHeight,
                                                     double bobbinWidth,
                                                     std::vector<double> bobbinCenterCoodinates,
                                                     uint8_t interleavingLevel = 1,
                                                     MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING,
                                                     MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING,
                                                     MAS::CoilAlignment turnsAlignment = MAS::CoilAlignment::CENTERED,
                                                     MAS::CoilAlignment sectionsAlignment = MAS::CoilAlignment::CENTERED,
                                                     std::vector<OpenMagnetics::Wire> wires = std::vector<OpenMagnetics::Wire>({}));

OpenMagnetics::Coil get_quick_toroidal_coil_no_compact(std::vector<int64_t> numberTurns,
                                                              std::vector<int64_t> numberParallels,
                                                              double bobbinRadialHeight,
                                                              double bobbinAngle,
                                                              double columnDepth,
                                                              uint8_t interleavingLevel = 1,
                                                              MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING,
                                                              MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING,
                                                              MAS::CoilAlignment turnsAlignment = MAS::CoilAlignment::CENTERED,
                                                              MAS::CoilAlignment sectionsAlignment = MAS::CoilAlignment::CENTERED,
                                                              std::vector<OpenMagnetics::Wire> wires = std::vector<OpenMagnetics::Wire>({}));


OpenMagnetics::Inputs get_quick_insulation_inputs(MAS::DimensionWithTolerance altitude,
                                                         MAS::Cti cti,
                                                         MAS::InsulationType insulation_type,
                                                         MAS::DimensionWithTolerance main_supply_voltage,
                                                         MAS::OvervoltageCategory overvoltage_category,
                                                         MAS::PollutionDegree pollution_degree,
                                                         std::vector<MAS::InsulationStandards> standards,
                                                         double maximumVoltageRms,
                                                         double maximumVoltagePeak,
                                                         double frequency,
                                                         MAS::WiringTechnology wiringTechnology = MAS::WiringTechnology::WOUND);

MAS::InsulationRequirements get_quick_insulation_requirements(MAS::DimensionWithTolerance altitude,
                                                                        MAS::Cti cti,
                                                                        MAS::InsulationType insulation_type,
                                                                        MAS::DimensionWithTolerance main_supply_voltage,
                                                                        MAS::OvervoltageCategory overvoltage_category,
                                                                        MAS::PollutionDegree pollution_degree,
                                                                        std::vector<MAS::InsulationStandards> standards);

Core get_quick_core(std::string shapeName,
                                          json basicGapping,
                                          int numberStacks = 1,
                                          std::string materialName = "N87");

OpenMagnetics::Magnetic get_quick_magnetic(std::string shapeName,
                                                  json basicGapping,
                                                  std::vector<int64_t> numberTurns,
                                                  int numberStacks = 1,
                                                  std::string materialName = "N87");
json get_ground_gap(double gapLength);
json get_distributed_gap(double gapLength, int numberGaps);
json get_spacer_gap(double gapLength);
json get_residual_gap();

void print(std::vector<double> data);
void print(std::vector<std::vector<double>> data);
void print(std::vector<int64_t> data);
void print(std::vector<uint64_t> data);
void print(std::vector<std::string> data);
void print(double data);
void print(std::string data);
void print_json(json data);


void check_sections_description(OpenMagnetics::Coil coil,
                                std::vector<int64_t> numberTurns,
                                std::vector<int64_t> numberParallels,
                                uint8_t interleavingLevel = 1,
                                MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING);

void check_layers_description(OpenMagnetics::Coil coil,
                                      MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING);


bool check_turns_description(OpenMagnetics::Coil coil);
bool check_wire_standards(OpenMagnetics::Coil coil);
void check_winding_losses(OpenMagnetics::Mas mas);

OpenMagnetics::Mas mas_loader(const std::filesystem::path& path);

// Helper to create Core and Coil from JSON strings, process them, and optionally create magnetic
std::pair<OpenMagnetics::Core, OpenMagnetics::Coil> prepare_core_and_coil_from_json(
    const std::string& coreJsonStr,
    const std::string& coilJsonStr);

OpenMagnetics::Magnetic prepare_magnetic_from_json(
    const std::string& coreJsonStr,
    const std::string& coilJsonStr);

// Struct to hold common painter test configuration
struct PainterTestConfig {
    std::vector<int64_t> numberTurns = {23, 13};
    std::vector<int64_t> numberParallels = {2, 2};
    uint8_t interleavingLevel = 2;
    int numberStacks = 1;
    double voltagePeakToPeak = 2000;
    std::string coreShape = "PQ 26/25";
    std::string coreMaterial = "3C97";
    double gapLength = 0.001;
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    double frequency = 125000;
    double magnetizingInductance = 0.001;
    double temperature = 25;
    WaveformLabel waveformLabel = WaveformLabel::TRIANGULAR;
    double dutyCycle = 0.5;
    double offset = 0;
    std::vector<std::string> wireNames = {};  // If empty, uses default wires
    std::vector<OpenMagnetics::Wire> customWires = {};  // For wires that need modification after lookup
    bool compactCoil = true;  // Whether to call delimit_and_compact
};

// Helper to prepare magnetic and inputs for painter tests
std::pair<OpenMagnetics::Magnetic, OpenMagnetics::Inputs> prepare_painter_test(const PainterTestConfig& config);

// Helper to create and configure a Coil from JSON string with all properties
OpenMagnetics::Coil prepare_coil_from_json(const std::string& coilJsonStr);

// Struct to hold coil winding configuration for complex test patterns
struct CoilWindingConfig {
    std::string coilJsonStr;
    std::vector<size_t> pattern = {};
    std::vector<double> proportionPerWinding = {};
    std::vector<std::vector<double>> marginPairs = {};
    size_t repetitions = 1;
    bool windCoil = true;  // If true, call coil.wind() automatically
};

// Helper to create and wind a coil from JSON with complex settings
// Handles object/array/value forms of _layersOrientation and _turnsAlignment
OpenMagnetics::Coil prepare_and_wind_coil(const CoilWindingConfig& config);

// Configuration for creating quick operating point inputs
struct QuickInputsConfig {
    double frequency = 100000;
    double magnetizingInductance = 100e-6;
    double temperature = 20;
    WaveformLabel label = WaveformLabel::TRIANGULAR;
    double peakToPeak = 2 * 1.73205;
    double dutyCycle = 0.5;
    double offset = 0;
};

// Helper to create quick inputs for testing using common defaults
OpenMagnetics::Inputs create_quick_test_inputs(const QuickInputsConfig& config = QuickInputsConfig{});

// Configuration for quick magnetic setup using create_quick_core and create_quick_coil
struct QuickMagneticConfig {
    std::vector<int64_t> numberTurns = {1, 1};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string coreShapeName = "E 35";
    std::string coreMaterialName = "A07";
    std::vector<std::string> wireNames = {};  // If empty, uses "Round 2.00 - Grade 1" for each winding
    int numberStacks = 1;
};

// Helper to create a quick Magnetic for testing
OpenMagnetics::Magnetic create_quick_test_magnetic(const QuickMagneticConfig& config = QuickMagneticConfig{});

} // namespace OpenMagneticsTesting