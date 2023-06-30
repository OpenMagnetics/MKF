#pragma once
#include "BobbinWrapper.h"
#include "InsulationMaterialWrapper.h"

#include "Constants.h"
#include "json.hpp"

#include <MAS.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json-schema.hpp>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;
#include <libInterpolate/Interpolate.hpp>

extern std::map<std::string, OpenMagnetics::CoreMaterial> coreMaterialDatabase;
extern std::map<std::string, OpenMagnetics::CoreShape> coreShapeDatabase;
extern std::map<std::string, OpenMagnetics::WireS> wireDatabase;
extern std::map<std::string, OpenMagnetics::BobbinWrapper> bobbinDatabase;
extern std::map<std::string, OpenMagnetics::InsulationMaterialWrapper> insulationMaterialDatabase;
extern std::map<std::string, OpenMagnetics::WireMaterial> wireMaterialDatabase;
extern OpenMagnetics::Constants constants;

extern _1D::LinearInterpolator<double> bobbinFillingFactorInterpWidth;
extern _1D::LinearInterpolator<double> bobbinFillingFactorInterpHeight;
extern _1D::LinearInterpolator<double> bobbinWindingWindowInterpWidth;
extern _1D::LinearInterpolator<double> bobbinWindingWindowInterpHeight;
extern std::map<std::string, _1D::LinearInterpolator<double>> wireFillingFactorInterps;
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


template<OpenMagnetics::DimensionalValues preferredValue>
double resolve_dimensional_values(OpenMagnetics::Dimension dimensionValue) {
    double doubleValue = 0;
    if (std::holds_alternative<OpenMagnetics::DimensionWithTolerance>(dimensionValue)) {
        switch (preferredValue) {
            case OpenMagnetics::DimensionalValues::MAXIMUM:
                if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                else if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                break;
            case OpenMagnetics::DimensionalValues::NOMINAL:
                if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_maximum().has_value() &&
                         std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue =
                        (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_maximum().value() +
                         std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_minimum().value()) /
                        2;
                else if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                else if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                break;
            case OpenMagnetics::DimensionalValues::MINIMUM:
                if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                else if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<OpenMagnetics::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                break;
            default:
                throw std::runtime_error("Unknown type of dimension, options are {MAXIMUM, NOMINAL, MINIMUM}");
        }
    }
    else if (std::holds_alternative<double>(dimensionValue)) {
        doubleValue = std::get<double>(dimensionValue);
    }
    else {
        throw std::runtime_error("Unknown variant in dimensionValue, holding variant");
    }
    return doubleValue;
}

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