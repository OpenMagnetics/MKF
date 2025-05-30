#pragma once
#include "constructive_models/Wire.h"
#include "constructive_models/Core.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Mas.h"
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


extern OpenMagnetics::Defaults defaults;
extern OpenMagnetics::Constants constants;
extern OpenMagnetics::Settings *settings;

extern std::vector<OpenMagnetics::Core> coreDatabase;
extern std::map<std::string, MAS::CoreMaterial> coreMaterialDatabase;
extern std::map<std::string, MAS::CoreShape> coreShapeDatabase;
extern std::vector<MAS::CoreShapeFamily> coreShapeFamiliesInDatabase;
extern std::map<std::string, OpenMagnetics::Wire> wireDatabase;
extern std::map<std::string, OpenMagnetics::Bobbin> bobbinDatabase;
extern std::map<std::string, OpenMagnetics::InsulationMaterial> insulationMaterialDatabase;
extern std::map<std::string, MAS::WireMaterial> wireMaterialDatabase;

extern tk::spline bobbinFillingFactorInterpWidth;
extern tk::spline bobbinFillingFactorInterpHeight;
extern tk::spline bobbinWindingWindowProportionInterpWidth;
extern tk::spline bobbinWindingWindowProportionInterpHeight;
extern std::map<std::string, tk::spline> wireFillingFactorInterps;
extern std::map<std::string, tk::spline> wirePackingFactorInterps;
extern std::map<std::string, tk::spline> wireCoatingThicknessProportionInterps;
extern std::map<std::string, tk::spline> wireConductingAreaProportionInterps;
extern double minBobbinWidth;
extern double maxBobbinWidth;
extern double minBobbinHeight;
extern double maxBobbinHeight;
extern double minWindingWindowWidth;
extern double maxWindingWindowWidth;
extern double minWindingWindowHeight;
extern double maxWindingWindowHeight;

extern std::map<std::string, std::variant<double, tk::spline>> initialPermeabilityMagneticFieldDcBiasInterps;
extern std::map<std::string, std::variant<double, tk::spline>> initialPermeabilityFrequencyInterps;
extern std::map<std::string, std::variant<double, tk::spline>> initialPermeabilityTemperatureInterps;
extern std::map<std::string, tk::spline> complexPermeabilityRealInterps;
extern std::map<std::string, tk::spline> complexPermeabilityImaginaryInterps;

extern std::map<std::string, tk::spline> lossFactorInterps;

extern bool _addInternalData;

extern std::map<std::string, double> minWireConductingDimensions;
extern std::map<std::string, double> maxWireConductingDimensions;
extern std::map<std::string, int64_t> minLitzWireNumberConductors;
extern std::map<std::string, int64_t> maxLitzWireNumberConductors;

extern OpenMagnetics::Cache magneticsCache;

using namespace MAS;

namespace OpenMagnetics {

double resolve_dimensional_values(Dimension dimensionValue, DimensionalValues preferredValue = DimensionalValues::NOMINAL);
bool check_requirement(DimensionWithTolerance requirement, double value);
Core find_core_by_name(std::string name);
CoreMaterial find_core_material_by_name(std::string name);
CoreShape find_core_shape_by_name(std::string name);
Wire find_wire_by_name(std::string name);
Wire find_wire_by_dimension(double dimension, std::optional<WireType> wireType=std::nullopt, std::optional<WireStandard> wireStandard=std::nullopt, bool obfuscate=true);
Bobbin find_bobbin_by_name(std::string name);
InsulationMaterial find_insulation_material_by_name(std::string name);
WireMaterial find_wire_material_by_name(std::string name);
CoreShape find_core_shape_by_winding_window_perimeter(double desiredPerimeter);


void clear_loaded_cores();
void clear_databases();
void load_cores();
void load_core_materials(std::optional<std::string> fileToLoad=std::nullopt);
void load_core_shapes(bool withAliases=true, std::optional<std::string> fileToLoad=std::nullopt);
void load_wires(std::optional<std::string> fileToLoad=std::nullopt);
void load_bobbins();
void load_insulation_materials();
void load_wire_materials();
void load_databases(json data, bool withAliases=true, bool addInternalData=true);

std::vector<std::string> get_core_shapes_names(std::optional<std::string> manufacturer);
std::vector<std::string> get_material_names(std::optional<std::string> manufacturer);
std::vector<std::string> get_shape_names(CoreShapeFamily family);
std::vector<std::string> get_shape_names();
std::vector<std::string> get_shape_family_dimensions(CoreShapeFamily family, std::optional<std::string> familySubtype=std::nullopt);
std::vector<std::string> get_shape_family_subtypes(CoreShapeFamily family);
std::vector<CoreShapeFamily> get_shape_families();
std::vector<std::string> get_wire_names();
std::vector<std::string> get_bobbin_names();
std::vector<std::string> get_insulation_material_names();
std::vector<std::string> get_wire_material_names();

std::vector<CoreMaterial> get_materials(std::optional<std::string> manufacturer);
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

double roundFloat(double value, int64_t decimals = 6);
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

} // namespace OpenMagnetics

