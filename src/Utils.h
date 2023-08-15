#pragma once
#include "BobbinWrapper.h"
#include "InsulationMaterialWrapper.h"

#include "Constants.h"

#include <MAS.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "spline.h"

extern std::map<std::string, OpenMagnetics::CoreMaterial> coreMaterialDatabase;
extern std::map<std::string, OpenMagnetics::CoreShape> coreShapeDatabase;
extern std::map<std::string, OpenMagnetics::WireS> wireDatabase;
extern std::map<std::string, OpenMagnetics::BobbinWrapper> bobbinDatabase;
extern std::map<std::string, OpenMagnetics::InsulationMaterialWrapper> insulationMaterialDatabase;
extern std::map<std::string, OpenMagnetics::WireMaterial> wireMaterialDatabase;
extern OpenMagnetics::Constants constants;

extern tk::spline bobbinFillingFactorInterpWidth;
extern tk::spline bobbinFillingFactorInterpHeight;
extern tk::spline bobbinWindingWindowInterpWidth;
extern tk::spline bobbinWindingWindowInterpHeight;
extern std::map<std::string, tk::spline> wireFillingFactorInterps;
extern double minBobbinWidth;
extern double maxBobbinWidth;
extern double minBobbinHeight;
extern double maxBobbinHeight;
extern double minWindingWindowWidth;
extern double maxWindingWindowWidth;
extern double minWindingWindowHeight;
extern double maxWindingWindowHeight;

extern std::map<std::string, double> minWireConductingWidths;
extern std::map<std::string, double> maxWireConductingWidths;


namespace OpenMagnetics {

enum class DimensionalValues : int {
    MAXIMUM,
    NOMINAL,
    MINIMUM
};
enum class GappingType : int {
    GRINDED,
    SPACER,
    RESIDUAL,
    DISTRIBUTED
};

double resolve_dimensional_values(OpenMagnetics::Dimension dimensionValue, DimensionalValues preferredValue = DimensionalValues::NOMINAL);
bool check_requirement(DimensionWithTolerance requirement, double value);
OpenMagnetics::CoreMaterial find_core_material_by_name(std::string name);
OpenMagnetics::CoreShape find_core_shape_by_name(std::string name);
OpenMagnetics::WireS find_wire_by_name(std::string name);
OpenMagnetics::Bobbin find_bobbin_by_name(std::string name);
OpenMagnetics::InsulationMaterial find_insulation_material_by_name(std::string name);
OpenMagnetics::WireMaterial find_wire_material_by_name(std::string name);

void load_databases(bool withAliases=true);

std::vector<std::string> get_material_names();
std::vector<std::string> get_shape_names();
std::vector<std::string> get_wire_names();

template<int decimals> double roundFloat(double value);

bool is_size_power_of_2(std::vector<double> data);

double roundFloat(double value, size_t decimals);

CoreShape flatten_dimensions(CoreShape shape);

double try_get_duty_cycle(Waveform waveform, double frequency);

bool check_collisions(std::map<std::string, std::vector<double>> dimensionsByName, std::map<std::string, std::vector<double>> coordinatesByName);
} // namespace OpenMagnetics