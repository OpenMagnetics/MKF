#include "Defaults.h"
#include "Constants.h"
#include "Utils.h"
#include "json.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

std::map<std::string, OpenMagnetics::CoreMaterial> coreMaterialDatabase;
std::map<std::string, OpenMagnetics::CoreShape> coreShapeDatabase;
std::map<std::string, OpenMagnetics::WireS> wireDatabase;
std::map<std::string, OpenMagnetics::BobbinWrapper> bobbinDatabase;
std::map<std::string, OpenMagnetics::InsulationMaterialWrapper> insulationMaterialDatabase;
std::map<std::string, OpenMagnetics::WireMaterial> wireMaterialDatabase;
OpenMagnetics::Constants constants = OpenMagnetics::Constants();

namespace OpenMagnetics {

void load_databases(bool withAliases) {
    {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
        auto dataFilePath = masPath + "data/materials.ndjson";
        std::ifstream ndjsonFile(dataFilePath);
        std::string myline;
        while (std::getline(ndjsonFile, myline)) {
            json jf = json::parse(myline);
            CoreMaterial coreMaterial(jf);
            coreMaterialDatabase[jf["name"]] = coreMaterial;
        }
    }

    {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
        auto dataFilePath = masPath + "data/shapes.ndjson";
        std::ifstream ndjsonFile(dataFilePath);
        std::string myline;
        while (std::getline(ndjsonFile, myline)) {
            json jf = json::parse(myline);
            CoreShape coreShape(jf);
            coreShapeDatabase[jf["name"]] = coreShape;
            if (withAliases) {
                for (auto& alias : jf["aliases"]) {
                    coreShapeDatabase[alias] = coreShape;
                }
            }
        }
    }

    {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
        auto dataFilePath = masPath + "data/wires.ndjson";
        std::ifstream ndjsonFile(dataFilePath);
        std::string myline;
        while (std::getline(ndjsonFile, myline)) {
            json jf = json::parse(myline);
            WireS wire(jf);
            wireDatabase[jf["name"]] = wire;
        }
    }

    {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
        auto dataFilePath = masPath + "data/bobbins.ndjson";
        std::ifstream ndjsonFile(dataFilePath);
        std::string myline;
        while (std::getline(ndjsonFile, myline)) {
            json jf = json::parse(myline);
            BobbinWrapper bobbin(jf);
            bobbinDatabase[jf["name"]] = bobbin;
        }
    }

    {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
        auto dataFilePath = masPath + "data/insulation_materials.ndjson";
        std::ifstream ndjsonFile(dataFilePath);
        std::string myline;
        while (std::getline(ndjsonFile, myline)) {
            json jf = json::parse(myline);
            InsulationMaterialWrapper insulationMaterial(jf);
            insulationMaterialDatabase[jf["name"]] = insulationMaterial;
        }
    }

    {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
        auto dataFilePath = masPath + "data/wire_materials.ndjson";
        std::ifstream ndjsonFile(dataFilePath);
        std::string myline;
        while (std::getline(ndjsonFile, myline)) {
            json jf = json::parse(myline);
            WireMaterial wireMaterial(jf);
            wireMaterialDatabase[jf["name"]] = wireMaterial;
        }
    }
}

std::vector<std::string> get_material_names() {
    if (coreMaterialDatabase.empty()) {
        load_databases();
    }

    std::vector<std::string> materialNames;

    for (auto& datum : coreMaterialDatabase) {
        materialNames.push_back(datum.first);
    }

    return materialNames;
}

std::vector<std::string> get_shape_names() {
    if (coreShapeDatabase.empty()) {
        load_databases(true);
    }

    std::vector<std::string> shapeNames;

    for (auto& datum : coreShapeDatabase) {
        shapeNames.push_back(datum.first);
    }

    return shapeNames;
}

OpenMagnetics::CoreMaterial find_core_material_by_name(std::string name) {
    if (coreMaterialDatabase.empty()) {
        load_databases();
    }
    if (coreMaterialDatabase.count(name)) {
        return coreMaterialDatabase[name];
    }
    else {
        throw std::runtime_error("Core material not found: " + name);
    }
}

OpenMagnetics::CoreShape find_core_shape_by_name(std::string name) {
    if (coreShapeDatabase.empty()) {
        load_databases();
    }
    if (coreShapeDatabase.count(name)) {
        return coreShapeDatabase[name];
    }
    else {
        throw std::runtime_error("Core shape not found: " + name);
    }
}

std::vector<std::string> get_wire_names() {
    if (wireDatabase.empty()) {
        load_databases(true);
    }

    std::vector<std::string> wireNames;

    for (auto& datum : wireDatabase) {
        wireNames.push_back(datum.first);
    }

    return wireNames;
}

OpenMagnetics::WireS find_wire_by_name(std::string name) {
    if (wireDatabase.empty()) {
        load_databases();
    }
    if (wireDatabase.count(name)) {
        return wireDatabase[name];
    }
    else {
        return json::parse("{}");
    }
}

OpenMagnetics::Bobbin find_bobbin_by_name(std::string name) {
    if (bobbinDatabase.empty()) {
        load_databases();
    }
    if (bobbinDatabase.count(name)) {
        return bobbinDatabase[name];
    }
    else {
        return json::parse("{}");
    }
}

OpenMagnetics::InsulationMaterial find_insulation_material_by_name(std::string name) {
    if (insulationMaterialDatabase.empty()) {
        load_databases();
    }
    if (insulationMaterialDatabase.count(name)) {
        return insulationMaterialDatabase[name];
    }
    else {
        return json::parse("{}");
    }
}

OpenMagnetics::WireMaterial find_wire_material_by_name(std::string name) {
    if (wireMaterialDatabase.empty()) {
        load_databases();
    }
    if (wireMaterialDatabase.count(name)) {
        return wireMaterialDatabase[name];
    }
    else {
        return json::parse("{}");
    }
}

double resolve_dimensional_values(OpenMagnetics::Dimension dimensionValue, DimensionalValues preferredValue) {
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

bool check_requirement(DimensionWithTolerance requirement, double value){
    if (requirement.get_minimum() && requirement.get_maximum()) {
        return requirement.get_minimum().value() <= value && value <= requirement.get_maximum().value();
    }
    else if (!requirement.get_minimum() && requirement.get_nominal() && requirement.get_maximum()) {
        return requirement.get_nominal().value() <= value && value <= requirement.get_maximum().value();
    }
    else if (requirement.get_minimum() && requirement.get_nominal() && !requirement.get_maximum()) {
        return requirement.get_minimum().value() <= value && value <= requirement.get_nominal().value();
    }
    else if (!requirement.get_minimum() && requirement.get_nominal() && !requirement.get_maximum()) {
        return value == requirement.get_nominal().value();
    }

    return false;
}

bool check_collisions(std::map<std::string, std::vector<double>> dimensionsByName, std::map<std::string, std::vector<double>> coordinatesByName){
    for (auto& leftDimension : dimensionsByName) {
        std::string leftName = leftDimension.first;
        std::vector<double> leftDimensions = leftDimension.second;
        std::vector<double> leftCoordinates = coordinatesByName[leftName];
        for (auto& rightDimension : dimensionsByName) {
            std::string rightName = rightDimension.first;
            if (rightName == leftName) {
                continue;
            }
            std::vector<double> rightDimensions = rightDimension.second;
            std::vector<double> rightCoordinates = coordinatesByName[rightName];

            if (roundFloat(fabs(leftCoordinates[0] - rightCoordinates[0]), 9) < roundFloat(leftDimensions[0] / 2 + rightDimensions[0] / 2, 9) &&
                roundFloat(fabs(leftCoordinates[1] - rightCoordinates[1]), 9) < roundFloat(leftDimensions[1] / 2 + rightDimensions[1] / 2, 9)) {
                return true;
            }
        }
    }
    return false;
}

} // namespace OpenMagnetics
