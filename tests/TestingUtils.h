#pragma once

#include "Constants.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "constructive_models/Wire.h"
#include "support/Settings.h"
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
                                                         MAS::IsolationClass insulation_type,
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
                                                                        MAS::IsolationClass insulation_type,
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

// Toroids: every turn clears the flat bar on each side of its section. A section's MAS margin is
// the half-thickness of the bar on its side; the bar's centre plane passes through the ring axis at
// the section's edge pushed out by the coil's own bar angle. Each turn's perpendicular distance to
// that plane is at least the margin plus the wire's radius, and neighbouring sections' bars do not
// overlap (a spread window may leave room between them).
void check_turns_clear_toroidal_bars(OpenMagnetics::Coil coil);
void check_winding_losses(OpenMagnetics::Mas mas);

// Asserts that an exported SVG file exists, is non-empty, has an <svg root and
// contains at least one drawn element. Callers should delete any stale file
// before painting so the existence check is meaningful.
void check_svg(const std::filesystem::path& svgPath);

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

// Applies the web engine's settings object (the JSON the WebFrontend hands to
// WebLibMKF set_settings) to OpenMagnetics::Settings, key for key as that binding
// does, so a fixture captured at the web engine proxy replays under the same
// settings in C++. A key the fixture lacks throws (json::at), never a default.
inline void apply_web_engine_settings(const nlohmann::json& s) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.set_magnetizing_inductance_include_air_inductance(s.at("magnetizingInductanceIncludeAirInductance").get<bool>());
    settings.set_coil_allow_margin_tape(s.at("coilAllowMarginTape").get<bool>());
    settings.set_coil_allow_insulated_wire(s.at("coilAllowInsulatedWire").get<bool>());
    settings.set_coil_fill_sections_with_margin_tape(s.at("coilFillSectionsWithMarginTape").get<bool>());
    settings.set_coil_wind_even_if_not_fit(s.at("coilWindEvenIfNotFit").get<bool>());
    settings.set_coil_delimit_and_compact(s.at("coilDelimitAndCompact").get<bool>());
    settings.set_coil_only_one_turn_per_layer_in_contiguous_rectangular(s.at("coilOnlyOneTurnPerLayerInContiguousRectangular").get<bool>());
    settings.set_coil_try_rewind(s.at("coilTryRewind").get<bool>());
    settings.set_coil_maximum_layers_planar(s.at("coilMaximumLayersPlanar").get<int>());
    settings.set_preferred_core_material_ferrite_manufacturer(s.at("preferredCoreMaterialFerriteManufacturer").get<std::string>());
    settings.set_preferred_core_material_powder_manufacturer(s.at("preferredCoreMaterialPowderManufacturer").get<std::string>());
    if (s.at("preferredWireStandard").is_null()) {
        settings.set_preferred_wire_standard(std::nullopt);
    }
    else {
        settings.set_preferred_wire_standard(s.at("preferredWireStandard").get<MAS::WireStandard>());
    }
    settings.set_coil_include_additional_coordinates(s.at("coilIncludeAdditionalCoordinates").get<bool>());
    settings.set_coil_use_real_winding_geometry(s.at("coilUseRealWindingGeometry").get<bool>());
    settings.set_use_only_cores_in_stock(s.at("useOnlyCoresInStock").get<bool>());
    settings.set_magnetic_field_number_points_x(s.at("magneticFieldNumberPointsX").get<int>());
    settings.set_magnetic_field_number_points_y(s.at("magneticFieldNumberPointsY").get<int>());
    settings.set_magnetic_field_mirroring_dimension(s.at("magneticFieldMirroringDimension").get<int>());
    settings.set_magnetic_field_include_fringing(s.at("magneticFieldIncludeFringing").get<bool>());
    settings.set_coil_adviser_maximum_number_wires(s.at("coilAdviserMaximumNumberWires").get<int>());
    settings.set_core_adviser_include_margin(s.at("coreIncludeMargin").get<bool>());
    settings.set_core_adviser_include_stacks(s.at("coreIncludeStacks").get<bool>());
    settings.set_core_adviser_include_distributed_gaps(s.at("coreIncludeDistributedGaps").get<bool>());
    settings.set_use_toroidal_cores(s.at("useToroidalCores").get<bool>());
    settings.set_use_concentric_cores(s.at("useConcentricCores").get<bool>());
    settings.set_magnetic_field_strength_model(static_cast<OpenMagnetics::MagneticFieldStrengthModels>(s.at("magneticFieldStrengthModel").get<int>()));
    settings.set_magnetic_field_strength_fringing_effect_model(static_cast<OpenMagnetics::MagneticFieldStrengthFringingEffectModels>(s.at("magneticFieldStrengthFringingEffectModel").get<int>()));
    settings.set_reluctance_model(static_cast<OpenMagnetics::ReluctanceModels>(s.at("reluctanceModel").get<int>()));
    settings.set_core_temperature_model(static_cast<OpenMagnetics::CoreTemperatureModels>(s.at("coreTemperatureModel").get<int>()));
    settings.set_core_thermal_resistance_model(static_cast<OpenMagnetics::CoreThermalResistanceModels>(s.at("coreThermalResistanceModel").get<int>()));
    settings.set_winding_skin_effect_losses_model(static_cast<OpenMagnetics::WindingSkinEffectLossesModels>(s.at("windingSkinEffectLossesModel").get<int>()));
    settings.set_winding_proximity_effect_losses_model(static_cast<OpenMagnetics::WindingProximityEffectLossesModels>(s.at("windingProximityEffectLossesModel").get<int>()));
    settings.set_stray_capacitance_model(static_cast<OpenMagnetics::StrayCapacitanceModels>(s.at("strayCapacitanceModel").get<int>()));
    settings.set_coil_enable_user_winding_losses_models(s.at("coilEnableUserWindingLossesModels").get<bool>());
    settings.set_core_per_column_winding_windows(s.at("corePerColumnWindingWindows").get<bool>());
    settings.set_coil_adviser_allow_lateral_placement(s.at("coilAdviserAllowLateralPlacement").get<bool>());
    settings.set_thermal_network_strict_geometry(s.at("thermalNetworkStrictGeometry").get<bool>());
    settings.set_allow_material_data_extrapolation(s.at("allowMaterialDataExtrapolation").get<bool>());
}

// Rebuilds a coil the way the web does after an adviser returns it, when the user turns real
// winding geometry on: MVB++ magnetic_autocomplete_safe takes the stored magnetic JSON without
// its layersDescription, builds Coil(json, false) (no wind), and runs magnetic_autocomplete with
// coil_use_real_winding_geometry on and MVB++'s lead bend policy (kRoundCornerBendFactor 1.05,
// no minimum bend radius). Settings are restored on every way out.
inline OpenMagnetics::Magnetic rebuild_with_real_winding_as_web(const OpenMagnetics::Magnetic& stored) {
    nlohmann::json magneticJson;
    to_json(magneticJson, stored);
    magneticJson.at("coil").erase("layersDescription");
    auto& settings = OpenMagnetics::Settings::GetInstance();
    OpenMagnetics::SettingsGuard<bool> realWindingGuard(settings,
        &OpenMagnetics::Settings::get_coil_use_real_winding_geometry,
        &OpenMagnetics::Settings::set_coil_use_real_winding_geometry, true);
    OpenMagnetics::SettingsGuard<std::optional<double>> bendFactorGuard(settings,
        &OpenMagnetics::Settings::get_coil_lead_bend_radius_factor,
        &OpenMagnetics::Settings::set_coil_lead_bend_radius_factor, std::optional<double>(1.05));
    OpenMagnetics::SettingsGuard<std::optional<double>> bendMinimumGuard(settings,
        &OpenMagnetics::Settings::get_coil_lead_minimum_bend_radius,
        &OpenMagnetics::Settings::set_coil_lead_minimum_bend_radius, std::optional<double>());
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagnetics::Core(magneticJson.at("core")));
    magnetic.set_coil(OpenMagnetics::Coil(magneticJson.at("coil"), false));
    return OpenMagnetics::magnetic_autocomplete(magnetic, nlohmann::json{});
}

} // namespace OpenMagneticsTesting