#pragma once
#include "constructive_models/Wire.h"
#include "constructive_models/Core.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Mas.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/InsulationMaterial.h"

#include <MAS.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>
#include <complex>
#include "spline.h"
#include <exception>

#include "Defaults.h"
#include "Constants.h"
#include "support/Settings.h"
#include "Definitions.h"
#include "Cache.h"
#include "support/Logger.h"

using namespace MAS;

namespace OpenMagnetics {

inline const OpenMagnetics::Defaults defaults = OpenMagnetics::Defaults();
inline const OpenMagnetics::Constants constants = OpenMagnetics::Constants();
inline OpenMagnetics::Settings& settings = OpenMagnetics::Settings::GetInstance();

inline std::map<OpenMagnetics::MagneticFilters, std::map<std::string, double>> _scorings;

inline std::vector<OpenMagnetics::Core> coreDatabase;
inline std::map<std::string, MAS::CoreMaterial> coreMaterialDatabase;
inline std::map<std::string, MAS::CoreShape> coreShapeDatabase;
inline std::vector<MAS::CoreShapeFamily> coreShapeFamiliesInDatabase;
inline std::map<std::string, OpenMagnetics::Wire> wireDatabase;
inline std::map<std::string, OpenMagnetics::Bobbin> bobbinDatabase;
inline std::map<std::string, OpenMagnetics::InsulationMaterial> insulationMaterialDatabase;
inline std::map<std::string, MAS::WireMaterial> wireMaterialDatabase;

inline OpenMagnetics::MagneticsCache magneticsCache;

void add_scoring(std::string name, OpenMagnetics::MagneticFilters filter, double scoring);
void clear_scoring();
std::optional<double> get_scoring(std::string name, OpenMagnetics::MagneticFilters filter);

// Legacy logging interface - prefer using Logger directly
// Verbosity levels: 0=ERROR, 1=WARNING, 2=INFO, 3+=DEBUG
void logEntry(std::string entry, std::string module = "", uint8_t entryVerbosity = 1);
std::string read_log();

bool check_requirement(DimensionWithTolerance requirement, double value);
Core find_core_by_name(std::string name);
CoreMaterial find_core_material_by_name(std::string name);
CoreShape find_core_shape_by_name(std::string name);
Wire find_wire_by_name(std::string name);
Wire find_wire_by_dimension(double dimension, std::optional<WireType> wireType=std::nullopt, std::optional<WireStandard> wireStandard=std::nullopt, bool obfuscate=true);
Bobbin find_bobbin_by_name(std::string name);
InsulationMaterial find_insulation_material_by_name(std::string name);
WireMaterial find_wire_material_by_name(std::string name);
CoreShape find_core_shape_by_winding_window_perimeter(double desiredPerimeter, std::optional<CoreShapeFamily> family=std::nullopt);
CoreShape find_core_shape_by_area_product(double desiredAreaProduct, std::optional<CoreShapeFamily> family=std::nullopt);
CoreShape find_core_shape_by_winding_window_area(double desiredWindingWindowArea, std::optional<CoreShapeFamily> family=std::nullopt);
CoreShape find_core_shape_by_winding_window_dimensions(double desiredWidthOrRadius, double desiredHeight, std::optional<CoreShapeFamily> family=std::nullopt);
CoreShape find_core_shape_by_effective_parameters(double desiredEffectiveLength, double desiredEffectiveArea, double desiredEffectiveVolume, std::optional<CoreShapeFamily> family=std::nullopt);
std::pair<MAS::CoreShape, double> find_closest_core_shape(std::vector<std::string> coreShapeCandidates, MagneticCoreSearchElement magneticCoreSearchElement);
double get_error_by_winding_window_perimeter(CoreShape shape, double perimeter);
double get_error_by_area_product(CoreShape shape, double areaProduct);
double get_error_by_winding_window_area(CoreShape shape, double windingWindowArea);
double get_error_by_winding_window_dimensions(CoreShape shape, double widthOrRadius, double height);
double get_error_by_effective_parameters(CoreShape shape, double effectiveLength, double effectiveArea, double effectiveVolume);


void clear_loaded_cores();
void clear_databases();
void load_cores(std::optional<std::string> fileToLoad=std::nullopt);
void load_core_materials(std::optional<std::string> fileToLoad=std::nullopt);
void load_advanced_core_materials(std::string fileToLoad, bool onlyDataFromManufacturer = true);
void load_core_shapes(bool withAliases=true, std::optional<std::string> fileToLoad=std::nullopt);
void load_wires(std::optional<std::string> fileToLoad=std::nullopt);
void load_bobbins();
void load_insulation_materials();
void load_wire_materials();
void load_databases(json data, bool withAliases=true, bool addInternalData=true);

std::vector<std::string> get_core_shape_names(std::string manufacturer);
std::vector<std::string> get_core_material_names(std::optional<std::string> manufacturer=std::nullopt);
std::vector<std::string> get_core_shape_names(CoreShapeFamily family);
std::vector<std::string> get_core_shape_names();
std::vector<std::string> get_shape_family_dimensions(CoreShapeFamily family, std::optional<std::string> familySubtype=std::nullopt);
std::vector<std::string> get_shape_family_subtypes(CoreShapeFamily family);
std::vector<CoreShapeFamily> get_core_shape_families();
std::vector<std::string> get_core_material_families(std::optional<MaterialType> materialType = std::nullopt);
std::vector<std::string> get_wire_names();
std::vector<std::string> get_bobbin_names();
std::vector<std::string> get_insulation_material_names();
std::vector<std::string> get_wire_material_names();

std::vector<CoreMaterial> get_materials(std::optional<std::string> manufacturer=std::nullopt);
std::vector<CoreShape> get_shapes(bool includeToroidal = true);
std::vector<Wire> get_wires(std::optional<WireType> wireType=std::nullopt, std::optional<WireStandard> wireStandard=std::nullopt);
std::vector<Bobbin> get_bobbins();
std::vector<InsulationMaterial> get_insulation_materials();
std::vector<WireMaterial> get_wire_materials();


bool is_size_power_of_2(std::vector<double> data);
size_t round_up_size_to_power_of_2(std::vector<double> data);
size_t round_up_size_to_power_of_2(size_t size);

std::vector<size_t> get_main_harmonic_indexes(OperatingPointExcitation excitation, double windingLossesHarmonicAmplitudeThreshold, std::string signal="current", std::optional<size_t> mainHarmonicIndex=std::nullopt);
std::vector<size_t> get_main_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::string signal="current", std::optional<size_t> mainHarmonicIndex=std::nullopt);
std::vector<size_t> get_main_harmonic_indexes(Harmonics harmonics, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);
std::vector<size_t> get_operating_point_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);
std::vector<size_t> get_excitation_harmonic_indexes(OperatingPointExcitation excitation, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);

double roundFloat(double value, int64_t decimals = 9);
double ceilFloat(double value, size_t decimals);
double floorFloat(double value, size_t decimals);

CoreShape flatten_dimensions(CoreShape shape);
std::map<std::string, double> flatten_dimensions(std::map<std::string, Dimension> dimensions);

double try_get_duty_cycle(Waveform waveform, double frequency);
std::complex<double> modified_bessel_first_kind(double order, std::complex<double> x);
std::complex<double> bessel_first_kind(double order, std::complex<double> z);
double kelvin_function_real(double order, double x);
double kelvin_function_imaginary(double order, double x);
double derivative_kelvin_function_real(double order, double x);
double derivative_kelvin_function_imaginary(double order, double x);
double comp_ellint_1(double x);
double comp_ellint_2(double x);

std::string to_title_case(std::string text);

double wound_distance_to_angle(double distance, double radius);
double angle_to_wound_distance(double angle, double radius);


bool check_collisions(std::map<std::string, std::vector<double>> dimensionsByName, std::map<std::string, std::vector<double>> coordinatesByName, bool roundTurn = false);
std::vector<IsolationSide> get_ordered_isolation_sides();
IsolationSide get_isolation_side_from_index(size_t index);
std::string get_isolation_side_name_from_index(size_t index);

std::vector<std::string> split(std::string s, std::string delimiter);
std::vector<double> linear_spaced_array(double a, double b, size_t N);
std::vector<double> logarithmic_spaced_array(double a, double b, size_t N);

double decibels_to_amplitude(double decibels);
double amplitude_to_decibels(double amplitude);

std::string fix_filename(std::string filename);
Inputs inputs_autocomplete(Inputs inputs, std::optional<Magnetic> magnetic = std::nullopt, json configuration = {});
Magnetic magnetic_autocomplete(Magnetic magnetic, json configuration = {});
Mas mas_autocomplete(Mas mas, bool simulate = true, json configuration = {});

std::map<std::string, double> normalize_scoring(std::map<std::string, double> scoring, double weight, std::map<std::string, bool> filterConfiguration);
std::map<std::string, double> normalize_scoring(std::map<std::string, double> scoring, double weight, bool invert, bool log);
std::vector<double> normalize_scoring(std::vector<double> scoring, double weight, std::map<std::string, bool> filterConfiguration);
std::vector<double> normalize_scoring(std::vector<double> scoring, double weight, bool invert, bool log);
std::map<std::string, double> normalize_scoring(std::map<std::string, double> scoring, OpenMagnetics::MagneticFilterOperation filterConfiguration);
std::vector<double> normalize_scoring(std::vector<double> scoring, OpenMagnetics::MagneticFilterOperation filterConfiguration);
void normalize_scoring(std::vector<std::pair<Mas, double>>* masesWithScoring, std::vector<double> scoring, double weight, std::map<std::string, bool> filterConfiguration);
void normalize_scoring(std::vector<std::pair<Mas, double>>* masesWithScoring, std::vector<double> scoring, MagneticFilterOperation filterConfiguration);

std::string generate_random_string(size_t length = 8);

size_t find_closest_index(std::vector<double> vector, double value);
double get_closest(double val1, double val2, double value);

template <typename Type>
std::string to_string(const Type & model) {
    json modelJson;
    to_json(modelJson, model);
    return modelJson;
}

std::vector<std::string> try_extract_core_shape_names(std::string coreShapeText);
std::optional<MAS::CoreShape> try_extract_core_shape(std::string coreShapeText);
std::vector<std::string> try_extract_core_material_names(std::string coreMaterialText);
std::optional<MAS::CoreMaterial> try_extract_core_material(std::string coreMaterialText);
std::vector<MAS::CoreGap> extract_core_gapping(OpenMagnetics::Core ungappedCore, int64_t numberTurns, double inductance);

Inputs get_defaults_inputs(OpenMagnetics::Magnetic magnetic);


} // namespace OpenMagnetics

