#include "Defaults.h"
#include "Constants.h"
#include "Utils.h"
#include "json.hpp"

#include <math.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

using json = nlohmann::json;

std::map<std::string, OpenMagnetics::CoreMaterial> coreMaterialDatabase;
std::map<std::string, OpenMagnetics::CoreShape> coreShapeDatabase;
std::map<std::string, OpenMagnetics::Wire> wireDatabase;
std::map<std::string, OpenMagnetics::BobbinWrapper> bobbinDatabase;
std::map<std::string, OpenMagnetics::InsulationMaterialWrapper> insulationMaterialDatabase;
std::map<std::string, OpenMagnetics::WireMaterial> wireMaterialDatabase;
OpenMagnetics::Constants constants = OpenMagnetics::Constants();

namespace OpenMagnetics {

void load_databases(bool withAliases) {
    {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
        auto dataFilePath = masPath + "data/core_materials.ndjson";
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
        auto dataFilePath = masPath + "data/core_shapes.ndjson";
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
            Wire wire(jf);
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

void load_databases(json data, bool withAliases) {
    for (auto& element : data["coreMaterials"].items()) {
        json jf = element.value();
        try {
            CoreMaterial coreMaterial(jf);
            coreMaterialDatabase[jf["name"]] = coreMaterial;
        }
        catch (...) {
            std::cout << "error: " << jf << std::endl;
            continue;
        }
    }

    for (auto& element : data["coreShapes"].items()) {
        json jf = element.value();

        int familySubtype = 1;
        if (jf["familySubtype"] != nullptr) {
            familySubtype = jf["familySubtype"];
        }

        jf["familySubtype"] = std::to_string(familySubtype);
        std::string family = jf["family"];
        std::transform(family.begin(), family.end(), family.begin(), ::toupper);
        std::replace(family.begin(), family.end(), ' ', '_');
        jf["family"] = family;
        CoreShape coreShape;
        from_json(jf, coreShape);

        coreShape.set_family(magic_enum::enum_cast<CoreShapeFamily>(family).value());

        coreShapeDatabase[jf["name"]] = coreShape;

        if (withAliases) {
            for (auto& alias : jf["aliases"]) {
                coreShapeDatabase[alias] = coreShape;
            }
        }
    }

    for (auto& element : data["wires"].items()) {
        json jf = element.value();
        std::string standardName;
        if (jf["standardName"].is_string()){
            standardName = jf["standardName"];
        }
        else {
            standardName = std::to_string(jf["standardName"].get<int>());
        }
        jf["standardName"] = standardName;
        Wire wire(jf);
        wireDatabase[jf["name"]] = wire;
    }

    for (auto& element : data["bobbins"].items()) {
        json jf = element.value();
        BobbinWrapper bobbin(jf);
        bobbinDatabase[jf["name"]] = bobbin;
    }

    for (auto& element : data["insulationMaterials"].items()) {
        json jf = element.value();
        InsulationMaterialWrapper insulationMaterial(jf);
        insulationMaterialDatabase[jf["name"]] = insulationMaterial;
    }

    for (auto& element : data["wireMaterials"].items()) {
        json jf = element.value();
        WireMaterial wireMaterial(jf);
        wireMaterialDatabase[jf["name"]] = wireMaterial;
    }
}

std::vector<std::string> get_material_names(std::optional<std::string> manufacturer) {
    if (coreMaterialDatabase.empty()) {
        load_databases();
    }

    std::vector<std::string> materialNames;

    for (auto& datum : coreMaterialDatabase) {
        if (!manufacturer) {
            materialNames.push_back(datum.first);
        }
        else {
            if (datum.second.get_manufacturer_info().get_name() == manufacturer.value() || manufacturer.value() == "") {
                materialNames.push_back(datum.first);
            }
        }
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

OpenMagnetics::Wire find_wire_by_name(std::string name) {
    if (wireDatabase.empty()) {
        load_databases();
    }
    if (wireDatabase.count(name)) {
        return wireDatabase[name];
    }
    else {
        throw std::runtime_error("wire not found: " + name);
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
        auto defaults = Defaults();
        return requirement.get_nominal().value() * (1 - defaults.magnetizingInductanceThresholdValidity) <= value &&
        value <= requirement.get_nominal().value() * (1 + defaults.magnetizingInductanceThresholdValidity);
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

// double modified_bessel_first_kind(double order, std::complex<double> x) {
//     std::complex<double> integral = 0;
//     std::cout << order << std::endl;
//     std::cout << x << std::endl;
//     double dtetha = 0.000001;
//     for (double tetha = 0; tetha < std::numbers::pi; tetha+=dtetha)
//     {
//         integral += exp(x * cos(tetha)) * cos(order * tetha) * dtetha;
//     }

//     auto bessel = integral / std::numbers::pi;
//     return bessel.real();
// }

std::complex<double> modified_bessel_first_kind(double order, std::complex<double> z) {
    std::complex<double> sum = 0;
    size_t limitK = 1000;
    for (size_t k = 0; k < limitK; ++k)
    {
        double divider = tgammal(k + 1) * tgammal(order + k + 1);
        if (std::isinf(divider)) {
            break;
        }
        sum += pow(0.25 * pow(z, 2), k) / divider;
        // if (std::isnan(sum.real())) {
        //     break;
        // }
    }

    auto bessel = sum * pow(0.5 * z, order);
    return bessel;
}

} // namespace OpenMagnetics
