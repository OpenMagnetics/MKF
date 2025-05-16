#include "support/Settings.h"
#include "Defaults.h"
#include "Constants.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/MagneticSimulator.h"
#include "support/Utils.h"
#include "json.hpp"

#include <math.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cfloat>
#include <limits>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include <cmrc/cmrc.hpp>

CMRC_DECLARE(data);

#ifndef M_PI_2
#    define M_PI_2 (1.57079632679)
#endif

using json = nlohmann::json;

std::vector<OpenMagnetics::CoreWrapper> coreDatabase;
std::map<std::string, OpenMagnetics::CoreMaterial> coreMaterialDatabase;
std::map<std::string, OpenMagnetics::CoreShape> coreShapeDatabase;
std::vector<OpenMagnetics::CoreShapeFamily> coreShapeFamiliesInDatabase;
std::map<std::string, OpenMagnetics::WireWrapper> wireDatabase;
std::map<std::string, OpenMagnetics::BobbinWrapper> bobbinDatabase;
std::map<std::string, OpenMagnetics::InsulationMaterialWrapper> insulationMaterialDatabase;
std::map<std::string, OpenMagnetics::WireMaterial> wireMaterialDatabase;
OpenMagnetics::Constants constants = OpenMagnetics::Constants();

bool _addInternalData = true;

namespace OpenMagnetics {

void load_cores() {
    auto settings = OpenMagnetics::Settings::GetInstance();
    bool includeToroidalCores = settings->get_use_toroidal_cores();
    bool includeConcentricCores = settings->get_use_concentric_cores();
    bool useOnlyCoresInStock = settings->get_use_only_cores_in_stock();

    auto fs = cmrc::data::get_filesystem();
    if (useOnlyCoresInStock && fs.exists("MAS/data/cores_stock.ndjson")) {
        auto data = fs.open("MAS/data/cores_stock.ndjson");
        std::string database = std::string(data.begin(), data.end());
        std::string delimiter = "\n";
        size_t pos = 0;
        std::string token;
        if (database.back() != delimiter.back()) {
            database += delimiter;
        }
        while ((pos = database.find(delimiter)) != std::string::npos) {
            token = database.substr(0, pos);
            json jf = json::parse(token);
            CoreType type;
            from_json(jf["functionalDescription"]["type"], type);
            if ((includeToroidalCores && type == CoreType::TOROIDAL) || (includeConcentricCores && type != CoreType::TOROIDAL)) {
                CoreWrapper core(jf, false, true, false);
                coreDatabase.push_back(core);
            }
            database.erase(0, pos + delimiter.length());
        }
    }
    else {
        auto data = fs.open("MAS/data/cores.ndjson");
        std::string database = std::string(data.begin(), data.end());
        std::string delimiter = "\n";

        if (database.back() == delimiter.back()) {
            database.pop_back();
        }
        database = std::regex_replace(database, std::regex("\n"), ", ");
        database = "[" + database + "]";
        json arr = json::parse(database);
        std::vector<OpenMagnetics::CoreWrapper> tempCoreDatabase;

        for (auto elem : arr) {

            if ((includeToroidalCores && elem["functionalDescription"]["type"] == "toroidal") || (includeConcentricCores && elem["functionalDescription"]["type"] != "toroidal")) {
                tempCoreDatabase.push_back(OpenMagnetics::CoreWrapper(elem));
            }
        }

        if (includeToroidalCores && includeConcentricCores) {
            coreDatabase = tempCoreDatabase;
        }
        else {
            for (auto core : tempCoreDatabase) {
                if ((includeToroidalCores && core.get_type() == CoreType::TOROIDAL) || (includeConcentricCores && core.get_type() != CoreType::TOROIDAL)) {
                    coreDatabase.push_back(core);
                }
            }
        }
    }
}

void clear_loaded_cores() {
    coreDatabase.clear();
}

void clear_databases() {
    coreDatabase.clear();
    coreMaterialDatabase.clear();
    coreShapeDatabase.clear();
    coreShapeFamiliesInDatabase.clear();
    wireDatabase.clear();
    bobbinDatabase.clear();
    insulationMaterialDatabase.clear();
    wireMaterialDatabase.clear();
}

void load_core_materials(std::optional<std::string> fileToLoad) {
    if (!_addInternalData) {
        return;
    }
    auto fs = cmrc::data::get_filesystem();
    {
        std::string database;
        if (fileToLoad) {
            database = fileToLoad.value();
        }
        else {
            auto data = fs.open("MAS/data/core_materials.ndjson");
            database = std::string(data.begin(), data.end());
        }

        std::string delimiter = "\n";
        size_t pos = 0;
        std::string token;
        if (database.back() != delimiter.back()) {
            database += delimiter;
        }
        while ((pos = database.find(delimiter)) != std::string::npos) {
            token = database.substr(0, pos);
            json jf = json::parse(token);
            CoreMaterial coreMaterial(jf);
            coreMaterialDatabase[jf["name"]] = coreMaterial;
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_core_shapes(bool withAliases, std::optional<std::string> fileToLoad) {
    if (!_addInternalData) {
        return;
    }
    auto fs = cmrc::data::get_filesystem();
    {
        auto settings = OpenMagnetics::Settings::GetInstance();
        bool includeToroidalCores = settings->get_use_toroidal_cores();
        bool includeConcentricCores = settings->get_use_concentric_cores();

        std::string database;
        if (fileToLoad) {
            database = fileToLoad.value();
        }
        else {
            auto data = fs.open("MAS/data/core_shapes.ndjson");
            database = std::string(data.begin(), data.end());
        }

        std::string delimiter = "\n";
        size_t pos = 0;
        std::string token;
        if (database.back() != delimiter.back()) {
            database += delimiter;
        }
        while ((pos = database.find(delimiter)) != std::string::npos) {
            token = database.substr(0, pos);
            json jf = json::parse(token);
            CoreShape coreShape(jf);
            if ((includeToroidalCores && coreShape.get_family() == CoreShapeFamily::T) || (includeConcentricCores && coreShape.get_family() != CoreShapeFamily::T)) {
                if (std::find(coreShapeFamiliesInDatabase.begin(), coreShapeFamiliesInDatabase.end(), coreShape.get_family()) == coreShapeFamiliesInDatabase.end()) {
                    coreShapeFamiliesInDatabase.push_back(coreShape.get_family());
                }
                coreShapeDatabase[jf["name"]] = coreShape;
                if (withAliases) {
                    for (auto& alias : jf["aliases"]) {
                        coreShapeDatabase[alias] = coreShape;
                    }
                }
            }
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_wires(std::optional<std::string> fileToLoad) {
    if (!_addInternalData) {
        return;
    }
    auto fs = cmrc::data::get_filesystem();
    {
        std::string database;
        if (fileToLoad) {
            database = fileToLoad.value();
        }
        else {
            auto data = fs.open("MAS/data/wires.ndjson");
            database = std::string(data.begin(), data.end());
        }

        std::string delimiter = "\n";
        size_t pos = 0;
        std::string token;
        if (database.back() != delimiter.back()) {
            database += delimiter;
        }
        while ((pos = database.find(delimiter)) != std::string::npos) {
            token = database.substr(0, pos);
            json jf = json::parse(token);
            WireWrapper wire(jf);
            wireDatabase[jf["name"]] = wire;
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_bobbins() {
    if (!_addInternalData) {
        return;
    }
    auto fs = cmrc::data::get_filesystem();
    {
        auto data = fs.open("MAS/data/bobbins.ndjson");
        std::string database = std::string(data.begin(), data.end());
        std::string delimiter = "\n";
        size_t pos = 0;
        std::string token;
        if (database.back() != delimiter.back()) {
            database += delimiter;
        }
        while ((pos = database.find(delimiter)) != std::string::npos) {
            token = database.substr(0, pos);
            json jf = json::parse(token);
            BobbinWrapper bobbin(jf);
            bobbinDatabase[jf["name"]] = bobbin;
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_insulation_materials() {
    if (!_addInternalData) {
        return;
    }
    auto fs = cmrc::data::get_filesystem();
    {
        auto data = fs.open("MAS/data/insulation_materials.ndjson");
        std::string database = std::string(data.begin(), data.end());
        std::string delimiter = "\n";
        size_t pos = 0;
        std::string token;
        if (database.back() != delimiter.back()) {
            database += delimiter;
        }
        while ((pos = database.find(delimiter)) != std::string::npos) {
            token = database.substr(0, pos);
            json jf = json::parse(token);
            InsulationMaterialWrapper insulationMaterial(jf);
            insulationMaterialDatabase[jf["name"]] = insulationMaterial;
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_wire_materials() {
    if (!_addInternalData) {
        return;
    }
    auto fs = cmrc::data::get_filesystem();
    {
        auto data = fs.open("MAS/data/wire_materials.ndjson");
        std::string database = std::string(data.begin(), data.end());
        std::string delimiter = "\n";
        size_t pos = 0;
        std::string token;
        if (database.back() != delimiter.back()) {
            database += delimiter;
        }
        while ((pos = database.find(delimiter)) != std::string::npos) {
            token = database.substr(0, pos);
            json jf = json::parse(token);
            WireMaterial wireMaterial(jf);
            wireMaterialDatabase[jf["name"]] = wireMaterial;
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_databases(json data, bool withAliases, bool addInternalData) {
    _addInternalData = addInternalData;
    if (addInternalData) {
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        if (coreShapeDatabase.empty()) {
            load_core_shapes();
        }
        if (wireDatabase.empty()) {
            load_wires();
        }
        if (bobbinDatabase.empty()) {
            load_bobbins();
        }
        if (insulationMaterialDatabase.empty()) {
            load_insulation_materials();
        }
        if (wireMaterialDatabase.empty()) {
            load_wire_materials();
        }
    }

    for (auto& element : data["coreMaterials"].items()) {
        json jf = element.value();
        try {
            CoreMaterial coreMaterial(jf);
            coreMaterialDatabase[jf["name"]] = coreMaterial;
        }
        catch (...) {
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
        if (jf.contains("standardName")) {
            if (jf["standardName"].is_string()){
                standardName = jf["standardName"];
            }
            else {
                standardName = std::to_string(jf["standardName"].get<int>());
            }
            jf["standardName"] = standardName;
        }
        WireWrapper wire(jf);
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

OpenMagnetics::CoreWrapper find_core_by_name(std::string name) {
    if (coreDatabase.empty()) {
        load_cores();
    }
    for (auto core : coreDatabase) {
        if (core.get_name() == name) {
            return core;
        }
    }
    throw std::runtime_error("Core not found: " + name);
}

OpenMagnetics::CoreMaterial find_core_material_by_name(std::string name) {
    if (coreMaterialDatabase.empty()) {
        load_core_materials();
    }
    if (coreMaterialDatabase.count(name)) {
        return coreMaterialDatabase[name];
    }
    else {
        for (auto [name, coreMaterial] : coreMaterialDatabase) {
            std::string commercialName;
            if (coreMaterial.get_commercial_name()) {
                commercialName = coreMaterial.get_commercial_name().value();
            }
            else {
                commercialName = coreMaterial.get_manufacturer_info().get_name() + " " + coreMaterial.get_name();
            }

            if (commercialName == name) {
                return coreMaterial;
            }
        }
        throw std::runtime_error("Core material not found: " + name);
    }
}

OpenMagnetics::CoreShape find_core_shape_by_name(std::string name) {
    if (coreShapeDatabase.empty()) {
        load_core_shapes();
    }
    if (coreShapeDatabase.count(name)) {
        return coreShapeDatabase[name];
    }
    else {
        for (const auto& [key, value] : coreShapeDatabase) {
            std::string dbName = key;
            dbName.erase(remove(dbName.begin(), dbName.end(), ' '), dbName.end());
            if (name == dbName) {
                return value;
            }
        }
        throw std::runtime_error("Core shape not found: " + name);
    }
}

std::vector<std::string> get_material_names(std::optional<std::string> manufacturer) {
    if (coreMaterialDatabase.empty()) {
        load_core_materials();
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

std::vector<std::string> get_core_shapes_names(std::optional<std::string> manufacturer) {
    if (coreDatabase.empty()) {
        load_cores();
    }

    std::vector<std::string> coreNames;

    for (auto& core : coreDatabase) {
        std::string coreShapeName = core.get_shape_name();
        if (!manufacturer) {
            if (std::find(coreNames.begin(), coreNames.end(), coreShapeName) == coreNames.end()) {
                coreNames.push_back(coreShapeName);
            }
        }
        else {
            if (!core.get_manufacturer_info()) {
                continue;
            }
            std::string manufacturerName = core.get_manufacturer_info()->get_name();

            if (manufacturerName == manufacturer.value() || manufacturer.value() == "") {
                if (std::find(coreNames.begin(), coreNames.end(), coreShapeName) == coreNames.end()) {
                    coreNames.push_back(coreShapeName);
                }
            }
        }
    }

    return coreNames;
}

std::vector<std::string> get_shape_names(CoreShapeFamily family) {
    if (coreShapeDatabase.empty()) {
        load_core_shapes(true);
    }

    std::vector<std::string> shapeNames;
 
    for (auto& [name, shape] : coreShapeDatabase) {
        if (shape.get_family() == family) {
            shapeNames.push_back(name);
        }
    }

    return shapeNames;
}

std::vector<std::string> get_shape_family_dimensions(CoreShapeFamily family, std::optional<std::string> familySubtype) {
    if (coreShapeDatabase.empty()) {
        load_core_shapes(true);
    }

    std::vector<std::string> distinctDimensions;
 
    for (auto& [name, shape] : coreShapeDatabase) {
        if (shape.get_family() == family) {
            if (familySubtype && shape.get_family_subtype()) {
                if (shape.get_family_subtype().value() != familySubtype.value()) {
                    continue;
                }
            }
            auto dimensions = shape.get_dimensions().value();
            for (auto& [dimensionKey, dimensionValue] : dimensions) {
                if (std::find(distinctDimensions.begin(), distinctDimensions.end(), dimensionKey) == distinctDimensions.end()) {
                    distinctDimensions.push_back(dimensionKey);
                }
            }
        }
    }

    std::sort(distinctDimensions.begin(), distinctDimensions.end(), [](std::string a, std::string b) {return a<b;});

    return distinctDimensions;
}


std::vector<std::string> get_shape_family_subtypes(CoreShapeFamily family) {
    if (coreShapeDatabase.empty()) {
        load_core_shapes(true);
    }

    std::vector<std::string> distinctSubtypes;
 
    for (auto& [name, shape] : coreShapeDatabase) {
        if (shape.get_family() == family) {
            if (shape.get_family_subtype()) {
                std::string familySubtype = shape.get_family_subtype().value();
                if (std::find(distinctSubtypes.begin(), distinctSubtypes.end(), familySubtype) == distinctSubtypes.end()) {
                    distinctSubtypes.push_back(familySubtype);
                }
            }
        }
    }

    std::sort(distinctSubtypes.begin(), distinctSubtypes.end(), [](std::string a, std::string b) {return a<b;});

    return distinctSubtypes;
}

std::vector<std::string> get_shape_names() {
    if (coreShapeDatabase.empty()) {
        load_core_shapes(true);
    }
    auto settings = OpenMagnetics::Settings::GetInstance();
    bool includeToroidalCores = settings->get_use_toroidal_cores();
    bool includeConcentricCores = settings->get_use_concentric_cores();

    std::vector<std::string> shapeNames;
 
    for (auto& [name, shape] : coreShapeDatabase) {
        if ((includeToroidalCores && shape.get_family() == CoreShapeFamily::T) || (includeConcentricCores && shape.get_family() != CoreShapeFamily::T)) {
            shapeNames.push_back(name);
        }
    }

    return shapeNames;
}


std::vector<OpenMagnetics::CoreShapeFamily> get_shape_families() {
    if (coreShapeDatabase.empty()) {
        load_core_shapes(true);
    }

    return coreShapeFamiliesInDatabase;
}

std::vector<std::string> get_wire_names() {
    if (wireDatabase.empty()) {
        load_wires();
    }

    std::vector<std::string> wireNames;

    for (auto& datum : wireDatabase) {
        wireNames.push_back(datum.first);
    }

    return wireNames;
}


std::vector<std::string> get_bobbin_names() {
    if (wireDatabase.empty()) {
        load_bobbins();
    }

    std::vector<std::string> bobbinNames;

    for (auto& datum : bobbinDatabase) {
        bobbinNames.push_back(datum.first);
    }

    return bobbinNames;
}


std::vector<std::string> get_insulation_material_names() {
    if (wireDatabase.empty()) {
        load_insulation_materials();
    }

    std::vector<std::string> insulationMaterialNames;

    for (auto& datum : insulationMaterialDatabase) {
        insulationMaterialNames.push_back(datum.first);
    }

    return insulationMaterialNames;
}


std::vector<std::string> get_wire_material_names() {
    if (wireDatabase.empty()) {
        load_wire_materials();
    }

    std::vector<std::string> wireMaterialNames;

    for (auto& datum : wireMaterialDatabase) {
        wireMaterialNames.push_back(datum.first);
    }

    return wireMaterialNames;
}

std::vector<OpenMagnetics::CoreMaterial> get_materials(std::optional<std::string> manufacturer) {
    if (coreMaterialDatabase.empty()) {
        load_core_materials();
    }

    std::vector<OpenMagnetics::CoreMaterial> materials;

    for (auto& datum : coreMaterialDatabase) {
        if (!manufacturer) {
            materials.push_back(datum.second);
        }
        else {
            if (datum.second.get_manufacturer_info().get_name() == manufacturer.value() || manufacturer.value() == "") {
                materials.push_back(datum.second);
            }
        }
    }

    return materials;
}

std::vector<OpenMagnetics::CoreShape> get_shapes(bool includeToroidal) {
    if (coreShapeDatabase.empty()) {
        load_core_shapes(true);
    }

    std::vector<OpenMagnetics::CoreShape> shapes;

    for (auto& datum : coreShapeDatabase) {
        if (includeToroidal || (datum.second.get_family() != CoreShapeFamily::T)) {
            shapes.push_back(datum.second);
        }
    }

    return shapes;
}


std::vector<OpenMagnetics::WireWrapper> get_wires(std::optional<WireType> wireType, std::optional<WireStandard> wireStandard) {
    if (wireDatabase.empty()) {
        load_wires();
    }

    std::vector<OpenMagnetics::WireWrapper> wires;

    for (auto& datum : wireDatabase) {
        if (wireStandard && !datum.second.get_standard()) {
            continue;
        }

        if (wireStandard && datum.second.get_standard().value() != wireStandard) {
            continue;
        }

        if (wireType && datum.second.get_type() != wireType) {
            continue;
        }
        
        wires.push_back(datum.second);
    }

    return wires;
}


std::vector<OpenMagnetics::BobbinWrapper> get_bobbins() {
    if (wireDatabase.empty()) {
        load_bobbins();
    }

    std::vector<OpenMagnetics::BobbinWrapper> bobbins;

    for (auto& datum : bobbinDatabase) {
        bobbins.push_back(datum.second);
    }

    return bobbins;
}


std::vector<OpenMagnetics::InsulationMaterialWrapper> get_insulation_materials() {
    if (wireDatabase.empty()) {
        load_insulation_materials();
    }

    std::vector<OpenMagnetics::InsulationMaterialWrapper> insulationMaterials;

    for (auto& datum : insulationMaterialDatabase) {
        insulationMaterials.push_back(datum.second);
    }

    return insulationMaterials;
}


std::vector<OpenMagnetics::WireMaterial> get_wire_materials() {
    if (wireDatabase.empty()) {
        load_wire_materials();
    }

    std::vector<OpenMagnetics::WireMaterial> wireMaterials;

    for (auto& datum : wireMaterialDatabase) {
        wireMaterials.push_back(datum.second);
    }

    return wireMaterials;
}

OpenMagnetics::WireWrapper find_wire_by_name(std::string name) {
    if (wireDatabase.empty()) {
        load_wires();
    }
    if (wireDatabase.count(name)) {
        return wireDatabase[name];
    }
    else {
        throw std::runtime_error("wire not found: " + name);
    }
}


OpenMagnetics::WireWrapper find_wire_by_dimension(double dimension, std::optional<WireType> wireType, std::optional<WireStandard> wireStandard, bool obfuscate) {
    if (wireDatabase.empty()) {
        load_wires();
    }

    auto wires = get_wires(wireType, wireStandard);
    double minimumDistance = DBL_MAX;
    OpenMagnetics::WireWrapper chosenWire;
    double minimumDimension = DBL_MAX;
    std::vector<OpenMagnetics::WireWrapper> possibleWires;

    for (auto wire : wires) {
        double distance = 0;

        switch (wire.get_type()) {
            case WireType::LITZ:
                continue;
                // throw std::runtime_error("Not defined for Litz wire");
            case WireType::ROUND:
                {
                    if (!wire.get_conducting_diameter()) {
                        throw std::runtime_error("Missing conducting diameter in round wire");
                    }
                    auto conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
                    distance = fabs(conductingDiameter - dimension);
                    break;
                }
            case WireType::PLANAR:
            case WireType::RECTANGULAR:
                {
                    if (!wire.get_conducting_height()) {
                        throw std::runtime_error("Missing conducting height in rectangular wire");
                    }
                    auto conductingHeight = resolve_dimensional_values(wire.get_conducting_height().value());
                    distance = fabs(conductingHeight - dimension);
                    break;
                }
            case WireType::FOIL:
                {
                    if (!wire.get_conducting_width()) {
                        throw std::runtime_error("Missing conducting width in foil wire");
                    }
                    auto conductingWidth = resolve_dimensional_values(wire.get_conducting_width().value());
                    distance = fabs(conductingWidth - dimension);
                    break;
                }
            default:
                throw std::runtime_error("Unknow type of wire");
        }
        if (distance == 0 || fabs(distance) <= 0.000000001) {
            possibleWires.push_back(wire);
        }

        if (distance < minimumDistance) {
            minimumDistance = distance;
            chosenWire = wire;
        }
        else if (distance == minimumDistance || fabs(distance - minimumDistance) <= 0.000000001) {
            if (wire.get_maximum_outer_dimension() < minimumDimension) {
                minimumDimension = wire.get_maximum_outer_dimension();
                chosenWire = wire;
            }
        }
    }

    double minimumOuterDimension = DBL_MAX;
    for (auto wire : possibleWires) {
        if (minimumOuterDimension > wire.get_maximum_outer_dimension()) {
            minimumOuterDimension = wire.get_maximum_outer_dimension();
            chosenWire = wire;
        }
    }

    if (obfuscate) {
        chosenWire.set_name(std::nullopt);
        chosenWire.set_coating(std::nullopt);
    }

    return chosenWire;
}

OpenMagnetics::BobbinWrapper find_bobbin_by_name(std::string name) {
    if (bobbinDatabase.empty()) {
        load_bobbins();
    }
    if (bobbinDatabase.count(name)) {
        return bobbinDatabase[name];
    }
    else {
        return json::parse("{}");
    }
}

OpenMagnetics::InsulationMaterialWrapper find_insulation_material_by_name(std::string name) {
    if (insulationMaterialDatabase.empty()) {
        load_insulation_materials();
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
        load_wire_materials();
    }
    if (wireMaterialDatabase.count(name)) {
        return wireMaterialDatabase[name];
    }
    else {
        return json::parse("{}");
    }
}

OpenMagnetics::CoreShape find_core_shape_by_winding_window_perimeter(double desiredPerimeter) {
    if (coreShapeDatabase.empty()) {
        load_core_shapes();
    }

    double minimumPerimeterError = DBL_MAX;
    OpenMagnetics::CoreShape closestShape;
    for (auto [name, shape] : coreShapeDatabase) {
        if (shape.get_family() != CoreShapeFamily::PQI && shape.get_family() != CoreShapeFamily::UI && shape.get_family() != CoreShapeFamily::UT) {
            auto corePiece = OpenMagnetics::CorePiece::factory(shape);
            auto mainColumn = corePiece->get_columns()[0];
            double perimeter;
            if (mainColumn.get_shape() == ColumnShape::RECTANGULAR || mainColumn.get_shape() == ColumnShape::IRREGULAR) {
                perimeter = 2 * (mainColumn.get_width() + mainColumn.get_depth());
            }
            else if (mainColumn.get_shape() == ColumnShape::ROUND) {
                perimeter = std::numbers::pi * mainColumn.get_width();
            }
            else if (mainColumn.get_shape() == ColumnShape::OBLONG) {
                perimeter = std::numbers::pi * mainColumn.get_width() + 2 * (mainColumn.get_depth() - mainColumn.get_width());
            }
            else {
                throw std::runtime_error("Unsupported column shape");
            }
 
            double perimeterError = fabs(perimeter - desiredPerimeter) / desiredPerimeter;
            if (perimeterError < minimumPerimeterError) {
                minimumPerimeterError = perimeterError;
                closestShape = shape;
            }
        }
    }
    return closestShape;
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
        if (requirement.get_maximum().value() < requirement.get_minimum().value()) {
            throw std::runtime_error("Minimum requirement cannot be larger than maximum");
        }
        if (requirement.get_nominal()) {
            if (requirement.get_maximum().value() < requirement.get_nominal().value()) {
                throw std::runtime_error("Nominal requirement cannot be larger than maximum");
            }
        }
        return requirement.get_minimum().value() <= value && value <= requirement.get_maximum().value();
    }
    else if (!requirement.get_minimum() && requirement.get_nominal() && requirement.get_maximum()) {
        if (requirement.get_maximum().value() < requirement.get_nominal().value()) {
            throw std::runtime_error("Nominal requirement cannot be larger than maximum");
        }
        return requirement.get_nominal().value() <= value && value <= requirement.get_maximum().value();
    }
    else if (requirement.get_minimum() && requirement.get_nominal() && !requirement.get_maximum()) {
        if (requirement.get_nominal().value() < requirement.get_minimum().value()) {
            throw std::runtime_error("Minimum requirement cannot be larger than nominal");
        }
        return requirement.get_minimum().value() <= value && value <= requirement.get_nominal().value();
    }
    else if (!requirement.get_minimum() && requirement.get_nominal() && !requirement.get_maximum()) {
        auto defaults = Defaults();
        return requirement.get_nominal().value() * (1 - defaults.magnetizingInductanceThresholdValidity) <= value &&
        value <= requirement.get_nominal().value() * (1 + defaults.magnetizingInductanceThresholdValidity);
    }
    else if (requirement.get_minimum() && !requirement.get_nominal() && !requirement.get_maximum()) {
        auto defaults = Defaults();
        return value > requirement.get_minimum().value();
    }
    else if (!requirement.get_minimum() && !requirement.get_nominal() && requirement.get_maximum()) {
        auto defaults = Defaults();
        return value > requirement.get_maximum().value();
    }

    return false;
}

bool check_collisions(std::map<std::string, std::vector<double>> dimensionsByName, std::map<std::string, std::vector<double>> coordinatesByName, bool roundWindingWinw){
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
            double distanceBetweenCenters = roundFloat(sqrt(pow(fabs(leftCoordinates[0] - rightCoordinates[0]), 2) + pow(fabs(leftCoordinates[1] - rightCoordinates[1]), 2)), 9);

            if (roundWindingWinw) {
                if (distanceBetweenCenters - roundFloat(leftDimensions[0] / 2 + rightDimensions[0] / 2, 9) < -1e-8) {
                    std::cout << "leftName: " << leftName << std::endl;
                    std::cout << "rightName: " << rightName << std::endl;
                    std::cout << "distanceBetweenCenters: " << distanceBetweenCenters << std::endl;
                    std::cout << "(distanceBetweenCenters - roundFloat(leftDimensions[0] / 2 + rightDimensions[0] / 2, 9)): " << (distanceBetweenCenters - roundFloat(leftDimensions[0] / 2 + rightDimensions[0] / 2, 9)) << std::endl;
                    std::cout << "(leftDimensions[0] / 2 + rightDimensions[0] / 2): " << (leftDimensions[0] / 2 + rightDimensions[0] / 2) << std::endl;
                    return true;
                }
            }
            else {
                if (roundFloat(fabs(leftCoordinates[0] - rightCoordinates[0]), 9) < roundFloat(leftDimensions[0] / 2 + rightDimensions[0] / 2, 9) &&
                    roundFloat(fabs(leftCoordinates[1] - rightCoordinates[1]), 9) < roundFloat(leftDimensions[1] / 2 + rightDimensions[1] / 2, 9)) {
                    std::cout << "leftName: " << leftName << std::endl;
                    std::cout << "rightName: " << rightName << std::endl;
                    std::cout << "roundFloat(fabs(leftCoordinates[0] - rightCoordinates[0]), 9): " << roundFloat(fabs(leftCoordinates[0] - rightCoordinates[0]), 9) << std::endl;
                    std::cout << "roundFloat(leftDimensions[0] / 2 + rightDimensions[0] / 2, 9): " << roundFloat(leftDimensions[0] / 2 + rightDimensions[0] / 2, 9) << std::endl;
                    std::cout << "roundFloat(fabs(leftCoordinates[1] - rightCoordinates[1]), 9): " << roundFloat(fabs(leftCoordinates[1] - rightCoordinates[1]), 9) << std::endl;
                    std::cout << "roundFloat(leftDimensions[1] / 2 + rightDimensions[1] / 2, 9): " << roundFloat(leftDimensions[1] / 2 + rightDimensions[1] / 2, 9) << std::endl;
                    return true;
                }
            }
        }
    }
    return false;
}

std::complex<double> modified_bessel_first_kind(double order, std::complex<double> z) {
    std::complex<double> sum = 0;
    std::complex<double> inc = 0;
    std::complex<double> aux = 0.25 * pow(z, 2);
    size_t limitK = 1000;
    for (size_t k = 0; k < limitK; ++k)
    {
        double divider = tgammaf(k + 1) * tgammaf(order + k + 1);
        if (std::isinf(divider)) {
            break;
        }
        inc = pow(aux, k) / divider;
        sum += inc;
        if (std::abs(inc) < std::abs(sum) * 0.0001){
            break;
        }
    }

    auto bessel = sum * pow(0.5 * z, order);
    return bessel;
}

std::complex<double> bessel_first_kind(double order, std::complex<double> z) {
    std::complex<double> sum = 0;
    std::complex<double> inc = 0;
    std::complex<double> aux = 0.25 * pow(z, 2);
    size_t limitK = 1000;
    for (size_t k = 0; k < limitK; ++k)
    {
        double divider = tgammaf(k + 1) * tgammaf(order + k + 1);
        if (std::isinf(divider)) {
            break;
        }
        inc = pow(-1, k) * pow(aux, k) / divider;
        sum += inc;
        if (std::abs(inc) < std::abs(sum) * 0.0001){
            break;
        }
    }
    auto bessel = sum * pow(0.5 * z, order);
    return bessel;
}

double kelvin_function_real(double order, double x) {
    return bessel_first_kind(order, x * exp(3.0 / 4 * std::numbers::pi * std::complex<double>{0.0, 1.0})).real();
}

double kelvin_function_imaginary(double order, double x) {
    return bessel_first_kind(order, x * exp(3.0 / 4 * std::numbers::pi * std::complex<double>{0.0, 1.0})).imag();
}

double derivative_kelvin_function_real(double order, double x) {
    return (kelvin_function_real(order + 1.0, x) + kelvin_function_imaginary(order + 1.0, x)) / sqrt(2) + order / x * kelvin_function_real(order, x);
}

double derivative_kelvin_function_imaginary(double order, double x) {
    return (kelvin_function_imaginary(order + 1.0, x) - kelvin_function_real(order + 1.0, x)) / sqrt(2) + order / x * kelvin_function_imaginary(order, x);
}

IsolationSide get_isolation_side_from_index(size_t index) {
    auto orderedIsolationSide = magic_enum::enum_cast<OrderedIsolationSide>(index).value();
    auto orderedIsolationSideString = std::string{magic_enum::enum_name(orderedIsolationSide)};
    return magic_enum::enum_cast<IsolationSide>(orderedIsolationSideString).value();
}

std::string get_isolation_side_name_from_index(size_t index) {
    auto isolationSide = get_isolation_side_from_index(index);
    std::string isolationSideString = std::string{magic_enum::enum_name(isolationSide)};
    return isolationSideString;
}

double comp_ellint_1(double x) {
    double k;      // modulus
    double m;      // the parameter of the elliptic function m = modulus^2
    double a;      // arithmetic mean
    double g;      // geometric mean
    double a_old;  // previous arithmetic mean
    double g_old;  // previous geometric mean
    double two_n;  // power of 2

    if ( x == 0.0 ) { return M_PI_2;    }
    k = fabs(x);
    m = k * k;
    if ( m == 1.0 ) { return std::numeric_limits<double>::quiet_NaN();    }

    a = 1.0;
    g = sqrt(1.0 - m);
    two_n = 1.0;
    for (int i=0;i<100;i++) 
    {
        g_old = g;
        a_old = a;
        a = 0.5 * (g_old + a_old);
        g = g_old * a_old;
        two_n += two_n;
        if ( fabs(a_old - g_old) <= (a_old * DBL_EPSILON) ) break;
        g = sqrt(g);
    } 
    return std::numbers::pi / 2 / a;
}

double comp_ellint_2(double x) {
    double k;      // modulus
    double m;      // the parameter of the elliptic function m = modulus^2
    double a;      // arithmetic mean
    double g;      // geometric mean
    double a_old;  // previous arithmetic mean
    double g_old;  // previous geometric mean
    double two_n;  // power of 2
    double sum;

    if ( x == 0.0 ) { return M_PI_2;    }
    k = fabs(x);
    m = k * k;
    if ( m == 1.0 ) { return 1.0;    }

    a = 1.0;
    g = sqrt(1.0 - m);
    two_n = 1.0;
    sum = 2.0 - m;
    for (int i=0;i<100;i++) 
    {
        g_old = g;
        a_old = a;
        a = 0.5 * (g_old + a_old);
        g = g_old * a_old;
        two_n += two_n;
        sum -= two_n * (a * a - g);
        if ( fabs(a_old - g_old) <= (a_old * DBL_EPSILON) ) break;
        g = sqrt(g);
    } 
    return (std::numbers::pi / 4 / a) * sum;
}


double ceilFloat(double value, size_t decimals) {
    return std::ceil(value * pow(10, decimals)) / pow(10, decimals);
}


double floorFloat(double value, size_t decimals) {
    return std::floor(value * pow(10, decimals)) / pow(10, decimals);
}

std::string to_title_case(std::string text) {
    
    std::string titleText = "";

    size_t pos = 0;
    size_t pre_pos = 0;

    pos = text.find(' ', pre_pos);

    while (pos != std::string::npos)
    {
        std::string sub = "";

        sub = text.substr(pre_pos, (pos - pre_pos));

        if (pre_pos != pos)
        {
            sub = text.substr(pre_pos, (pos - pre_pos));
        }
        else 
        {
            sub = text.substr(pre_pos, 1);
        }
        
        sub[0] = std::toupper(sub[0]);
        titleText += sub + text[pos];

        if (pos < (text.length() - 1))
        {
            pre_pos = (pos + 1);
        }
        else
        {
            pre_pos = pos;
            break;
        }

        pos = text.find(' ', pre_pos);

    }

    std::string sub = text.substr(pre_pos, std::string::npos);
    sub[0] = std::toupper(sub[0]);
    titleText += sub;
    
    return titleText;
}

double wound_distance_to_angle(double distance, double radius) {
    double angle = 2 * asin((distance / 2) / radius) * 180 / std::numbers::pi;
    if (std::isnan(angle)) {
        return 360;
    }
    return angle;
}

double angle_to_wound_distance(double angle, double radius) {
    return 2 * sin(angle / 2 / 180 * std::numbers::pi) * radius;
}

bool is_size_power_of_2(std::vector<double> data) {
    return data.size() > 0 && ((data.size() & (data.size() - 1)) == 0);
}

size_t round_up_size_to_power_of_2(std::vector<double> data) {
    return pow(2, ceil(log(data.size()) / log(2)));
}

size_t round_up_size_to_power_of_2(size_t size) {
    return pow(2, ceil(log(size) / log(2)));
}


std::vector<size_t> get_main_harmonic_indexes(OperatingPointExcitation excitation, double windingLossesHarmonicAmplitudeThreshold, std::string signal, std::optional<size_t> mainHarmonicIndex) {
    std::vector<size_t> mainHarmonicIndexes;
    double mainOrMaximumHarmonicAmplitudeTimesRootFrequency = 0;
    SignalDescriptor signalDescriptor;
    if (signal == "current") {
        if (!excitation.get_current()) {
            throw std::runtime_error("Missing current");
        }
        signalDescriptor = excitation.get_current().value();
    }
    else if (signal == "voltage") {
        if (!excitation.get_voltage()) {
            throw std::runtime_error("Missing voltage");
        }
        signalDescriptor = excitation.get_voltage().value();
    }
    else if (signal == "magnetizingCurrent") {
        if (!excitation.get_magnetizing_current()) {
            return mainHarmonicIndexes;
        }
        signalDescriptor = excitation.get_magnetizing_current().value();
    }
    else {
         throw std::runtime_error("Not supported harmonic common index extraction for " + signal);
    }
    if (!signalDescriptor.get_harmonics()) {
        if (!signalDescriptor.get_waveform()) {
            throw std::runtime_error("Missing harmonics");
        }
    }
    auto harmonics = signalDescriptor.get_harmonics().value();
    size_t maximumCommonIndex = harmonics.get_amplitudes().size();
    if (!mainHarmonicIndex) {
        for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            mainOrMaximumHarmonicAmplitudeTimesRootFrequency = std::max(harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex]), mainOrMaximumHarmonicAmplitudeTimesRootFrequency);
        }
    }
    else {
        mainOrMaximumHarmonicAmplitudeTimesRootFrequency = harmonics.get_amplitudes()[mainHarmonicIndex.value()] * sqrt(harmonics.get_frequencies()[mainHarmonicIndex.value()]);
    }

    if (mainOrMaximumHarmonicAmplitudeTimesRootFrequency == 0) {
        return mainHarmonicIndexes;
    }

    for (size_t harmonicIndex = 1; harmonicIndex < maximumCommonIndex; ++harmonicIndex) {

        if ((harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex])) < mainOrMaximumHarmonicAmplitudeTimesRootFrequency * windingLossesHarmonicAmplitudeThreshold) {
            continue;
        }
        if (std::find(mainHarmonicIndexes.begin(), mainHarmonicIndexes.end(), harmonicIndex) == mainHarmonicIndexes.end()) {

            mainHarmonicIndexes.push_back(harmonicIndex);
        }
    }

    return mainHarmonicIndexes;
}

std::vector<size_t> get_main_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::string signal, std::optional<size_t> mainHarmonicIndex) {
    std::vector<size_t> mainHarmonicIndexesInOperatingPoint;
    for (size_t windingIndex = 0; windingIndex < operatingPoint.get_excitations_per_winding().size(); ++windingIndex) {
        auto mainHarmonicIndexes = get_main_harmonic_indexes(operatingPoint.get_excitations_per_winding()[windingIndex], windingLossesHarmonicAmplitudeThreshold, signal, mainHarmonicIndex);

        for (auto index : mainHarmonicIndexes) {
            if (std::find(mainHarmonicIndexesInOperatingPoint.begin(), mainHarmonicIndexesInOperatingPoint.end(), index) == mainHarmonicIndexesInOperatingPoint.end()) {
                mainHarmonicIndexesInOperatingPoint.push_back(index);
            }
        }
    }

    return mainHarmonicIndexesInOperatingPoint;
}


std::vector<size_t> get_operating_point_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex) {
    std::vector<size_t> commonHarmonicIndexes;
    auto currentCommonHarmonicIndexes = get_main_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold, "current", mainHarmonicIndex);
    auto voltageCommonHarmonicIndexes = get_main_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold, "voltage", mainHarmonicIndex);
    auto magnetizingCurrentCommonHarmonicIndexes = get_main_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold, "magnetizingCurrent", mainHarmonicIndex);

    for (auto index : currentCommonHarmonicIndexes) {
        if (std::find(commonHarmonicIndexes.begin(), commonHarmonicIndexes.end(), index) == commonHarmonicIndexes.end()) {
            commonHarmonicIndexes.push_back(index);
        }
    }
    for (auto index : voltageCommonHarmonicIndexes) {
        if (std::find(commonHarmonicIndexes.begin(), commonHarmonicIndexes.end(), index) == commonHarmonicIndexes.end()) {
            commonHarmonicIndexes.push_back(index);
        }
    }
    for (auto index : magnetizingCurrentCommonHarmonicIndexes) {
        if (std::find(commonHarmonicIndexes.begin(), commonHarmonicIndexes.end(), index) == commonHarmonicIndexes.end()) {
            commonHarmonicIndexes.push_back(index);
        }
    }
    return commonHarmonicIndexes;
}
std::vector<size_t> get_excitation_harmonic_indexes(OperatingPointExcitation excitation, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex) {
    std::vector<size_t> commonHarmonicIndexes;
    auto currentCommonHarmonicIndexes = get_main_harmonic_indexes(excitation, windingLossesHarmonicAmplitudeThreshold, "current", mainHarmonicIndex);
    auto voltageCommonHarmonicIndexes = get_main_harmonic_indexes(excitation, windingLossesHarmonicAmplitudeThreshold, "voltage", mainHarmonicIndex);
    auto magnetizingCurrentCommonHarmonicIndexes = get_main_harmonic_indexes(excitation, windingLossesHarmonicAmplitudeThreshold, "magnetizingCurrent", mainHarmonicIndex);

    for (auto index : currentCommonHarmonicIndexes) {
        if (std::find(commonHarmonicIndexes.begin(), commonHarmonicIndexes.end(), index) == commonHarmonicIndexes.end()) {
            commonHarmonicIndexes.push_back(index);
        }
    }
    for (auto index : voltageCommonHarmonicIndexes) {
        if (std::find(commonHarmonicIndexes.begin(), commonHarmonicIndexes.end(), index) == commonHarmonicIndexes.end()) {
            commonHarmonicIndexes.push_back(index);
        }
    }
    for (auto index : magnetizingCurrentCommonHarmonicIndexes) {
        if (std::find(commonHarmonicIndexes.begin(), commonHarmonicIndexes.end(), index) == commonHarmonicIndexes.end()) {
            commonHarmonicIndexes.push_back(index);
        }
    }
    return commonHarmonicIndexes;
}

std::vector<size_t> get_main_harmonic_indexes(Harmonics harmonics, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex) {
    double mainOrMaximumHarmonicAmplitudeTimesRootFrequency = 0;
    if (!mainHarmonicIndex) {
        for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            mainOrMaximumHarmonicAmplitudeTimesRootFrequency = std::max(harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex]), mainOrMaximumHarmonicAmplitudeTimesRootFrequency);
        }
    }
    else {
        mainOrMaximumHarmonicAmplitudeTimesRootFrequency = harmonics.get_amplitudes()[mainHarmonicIndex.value()] * sqrt(harmonics.get_frequencies()[mainHarmonicIndex.value()]);
    }

    std::vector<size_t> mainHarmonicIndexes;
    for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
        if ((harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex])) < mainOrMaximumHarmonicAmplitudeTimesRootFrequency * windingLossesHarmonicAmplitudeThreshold) {
            continue;
        }

        if (std::find(mainHarmonicIndexes.begin(), mainHarmonicIndexes.end(), harmonicIndex) == mainHarmonicIndexes.end()) {
            mainHarmonicIndexes.push_back(harmonicIndex);
        }
    }
    return mainHarmonicIndexes;
}

std::vector<std::string> split(std::string s, std::string delimiter) {
    std::vector<std::string> tokens;
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        tokens.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    tokens.push_back(s);

    return tokens;
}

std::vector<double> linear_spaced_array(double a, double b, size_t N) {
    double h = (b - a) / static_cast<double>(N-1);
    std::vector<double> xs(N);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}


std::vector<double> logarithmic_spaced_array(double a, double b, size_t N) {
    double logmin = log10(a);
    double logmax = log10(b);

    double h = (logmax - logmin) / static_cast<double>(N-1);
    std::vector<double> xs(N);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = log10(a); x != xs.end(); ++x, val += h) {
        *x = pow(10, val);
    }
    return xs;
}


double decibels_to_amplitude(double decibels) {
    return pow(10, decibels / 20);
}
double amplitude_to_decibels(double amplitude) {
    return 20 * log10(amplitude);
}

std::string fix_filename(std::string filename) {
    filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex(" "), "_")).string();
    filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex("\\,"), "_")).string();
    filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex("\\."), "_")).string();
    filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex("\\:"), "_")).string();
    filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex("\\/"), "_")).string();
    return filename;
}

SignalDescriptor standardize_signal_descriptor(SignalDescriptor signalDescriptor, double frequency) {

    auto standardSignalDescriptor = InputsWrapper::standardize_waveform(signalDescriptor, frequency);
    if (standardSignalDescriptor.get_harmonics()) {
        auto processed = InputsWrapper::calculate_processed_data(standardSignalDescriptor.get_harmonics().value(), standardSignalDescriptor.get_waveform().value(), true);
        standardSignalDescriptor.set_processed(processed);
    }
    else {
        auto processed = InputsWrapper::calculate_processed_data(standardSignalDescriptor.get_waveform().value(), frequency, true);
        standardSignalDescriptor.set_processed(processed);
    }

    return standardSignalDescriptor;
}

OperatingPointExcitation calculate_reflected_secondary(OperatingPointExcitation primaryExcitation, double turnRatio){
    OpenMagnetics::OperatingPointExcitation excitationOfThisWinding(primaryExcitation);
    auto currentSignalDescriptorProcessed = OpenMagnetics::InputsWrapper::calculate_basic_processed_data(primaryExcitation.get_current().value().get_waveform().value());
    auto voltageSignalDescriptorProcessed = OpenMagnetics::InputsWrapper::calculate_basic_processed_data(primaryExcitation.get_voltage().value().get_waveform().value());

    auto voltageSignalDescriptor = OpenMagnetics::InputsWrapper::reflect_waveform(primaryExcitation.get_voltage().value(), 1.0 / turnRatio, voltageSignalDescriptorProcessed.get_label());
    auto currentSignalDescriptor = OpenMagnetics::InputsWrapper::reflect_waveform(primaryExcitation.get_current().value(), turnRatio, currentSignalDescriptorProcessed.get_label());

    auto voltageSampledWaveform = OpenMagnetics::InputsWrapper::calculate_sampled_waveform(voltageSignalDescriptor.get_waveform().value(), excitationOfThisWinding.get_frequency());
    voltageSignalDescriptor.set_harmonics(OpenMagnetics::InputsWrapper::calculate_harmonics_data(voltageSampledWaveform, excitationOfThisWinding.get_frequency()));
    voltageSignalDescriptor.set_processed(OpenMagnetics::InputsWrapper::calculate_processed_data(voltageSignalDescriptor, voltageSampledWaveform, true));

    auto currentSampledWaveform = OpenMagnetics::InputsWrapper::calculate_sampled_waveform(currentSignalDescriptor.get_waveform().value(), excitationOfThisWinding.get_frequency());
    currentSignalDescriptor.set_harmonics(OpenMagnetics::InputsWrapper::calculate_harmonics_data(currentSampledWaveform, excitationOfThisWinding.get_frequency()));
    currentSignalDescriptor.set_processed(OpenMagnetics::InputsWrapper::calculate_processed_data(currentSignalDescriptor, currentSampledWaveform, true));

    excitationOfThisWinding.set_voltage(voltageSignalDescriptor);
    excitationOfThisWinding.set_current(currentSignalDescriptor);

    return excitationOfThisWinding;
}

MasWrapper mas_autocomplete(MasWrapper mas, bool simulate) {

    //Inputs
    size_t numberWindings = mas.get_inputs().get_design_requirements().get_turns_ratios().size() + 1;
    if (!mas.get_inputs().get_design_requirements().get_isolation_sides()) {
        for (size_t i = 0; i < numberWindings; i++) {
            std::vector<IsolationSide> isolationSides;
            isolationSides.push_back(get_isolation_side_from_index(i));
            mas.get_mutable_inputs().get_mutable_design_requirements().set_isolation_sides(isolationSides);
        }
    }
    else {
        std::vector<IsolationSide> isolationSides;
        std::vector<IsolationSide> currentIsolationSides = mas.get_inputs().get_design_requirements().get_isolation_sides().value();
        for (size_t i = 0; i < numberWindings; i++) {
            if (currentIsolationSides.size() <= i) {
                isolationSides.push_back(get_isolation_side_from_index(i));
            }
            else {
                isolationSides.push_back(currentIsolationSides[i]);
            }
        }
        mas.get_mutable_inputs().get_mutable_design_requirements().set_isolation_sides(isolationSides);
    }

    auto [result, resultDescription] = mas.get_mutable_inputs().check_integrity();

    for (size_t operatingPointIndex = 0; operatingPointIndex <  mas.get_mutable_inputs().get_mutable_operating_points().size(); operatingPointIndex++) {
        for (size_t windingIndex = 0; windingIndex < numberWindings; windingIndex++) {
            if (mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_current()) {
                auto current = mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_current().value();
                if (!current.get_waveform() ||
                    !current.get_processed()) {
                    current = standardize_signal_descriptor(mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_current().value(),  mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_frequency());
                }
                if (!current.get_harmonics()) {
                    auto sampledCurrentWaveform = InputsWrapper::calculate_sampled_waveform(current.get_waveform().value(), mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_frequency());
                    auto harmonics = InputsWrapper::calculate_harmonics_data(sampledCurrentWaveform, mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_frequency());
                    current.set_harmonics(harmonics);
                }
                mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_current(current);
            }
            if (mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_voltage()) {
                auto voltage = mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_voltage().value();
                if (!voltage.get_waveform() ||
                    !voltage.get_processed()) {
                    voltage = standardize_signal_descriptor(mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_voltage().value(),  mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_frequency());
                }
                if (!voltage.get_harmonics()) {
                    auto sampledvoltageWaveform = InputsWrapper::calculate_sampled_waveform(voltage.get_waveform().value(), mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_frequency());
                    auto harmonics = InputsWrapper::calculate_harmonics_data(sampledvoltageWaveform, mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].get_frequency());
                    voltage.set_harmonics(harmonics);
                }
                mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_voltage(voltage);
            }
        }
    }

    auto magnetic = magnetic_autocomplete(mas.get_magnetic());
    mas.set_magnetic(magnetic);

    // Magnetizing current

    OpenMagnetics::MagnetizingInductance magnetizingInductanceObj;
    for (size_t operatingPointIndex = 0; operatingPointIndex <  mas.get_mutable_inputs().get_mutable_operating_points().size(); operatingPointIndex++) {
        for (size_t windingIndex = 0; windingIndex < numberWindings; windingIndex++) {
            double magnetizingInductance = magnetizingInductanceObj.calculate_inductance_from_number_turns_and_gapping(mas.get_mutable_magnetic().get_mutable_core(), mas.get_mutable_magnetic().get_mutable_coil(), &mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex]).get_magnetizing_inductance().get_nominal().value();
            auto excitation = mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex];

            auto magnetizingCurrent = OpenMagnetics::InputsWrapper::calculate_magnetizing_current(excitation, magnetizingInductance, true, 0.0);

            if (excitation.get_voltage()) {
                if (excitation.get_voltage().value().get_processed()) {
                    if (excitation.get_voltage().value().get_processed().value().get_duty_cycle()) {
                        auto processed = magnetizingCurrent.get_processed().value();
                        processed.set_duty_cycle(excitation.get_voltage().value().get_processed().value().get_duty_cycle().value());
                        magnetizingCurrent.set_processed(processed);
                    }
                }
            }
            mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_magnetizing_current(magnetizingCurrent);

            if (windingIndex == 0 && numberWindings == 2 && mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding().size() == 1) {
                auto turnRatio = mas.get_mutable_magnetic().get_turns_ratios()[0];
                auto secondaryExcitation = calculate_reflected_secondary(mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0], turnRatio);
                mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding().push_back(secondaryExcitation);
            }
        }
    }

    if (simulate) {
        // Outputs
        mas = MagneticSimulator().simulate(mas.get_inputs(), mas.get_magnetic());
    }

    return mas;
}

MagneticWrapper magnetic_autocomplete(MagneticWrapper magnetic) {
    // Core
    auto shape = magnetic.get_mutable_core().resolve_shape();

    if (magnetic.get_mutable_core().get_shape_family() == CoreShapeFamily::T) {
        magnetic.get_mutable_core().get_mutable_functional_description().set_type(CoreType::TOROIDAL);
        shape.set_magnetic_circuit(MagneticCircuit::CLOSED);
        magnetic.get_mutable_core().get_mutable_functional_description().get_mutable_gapping().clear();
    }
    else {
        magnetic.get_mutable_core().get_mutable_functional_description().set_type(CoreType::TWO_PIECE_SET);
        shape.set_magnetic_circuit(MagneticCircuit::OPEN);
    }
    magnetic.get_mutable_core().get_mutable_functional_description().set_shape(shape);

    auto material = magnetic.get_mutable_core().resolve_material();
    magnetic.get_mutable_core().get_mutable_functional_description().set_material(material);

    if (!magnetic.get_core().get_processed_description()) {
        magnetic.get_mutable_core().process_data();
        magnetic.get_mutable_core().process_gap();
    }

    if (!magnetic.get_core().get_geometrical_description()) {
        auto geometricalDescription = magnetic.get_mutable_core().create_geometrical_description();
        magnetic.get_mutable_core().set_geometrical_description(geometricalDescription);
    }

    // Coil

    for (size_t i = 0; i < magnetic.get_coil().get_functional_description().size(); i++) {
        if (magnetic.get_coil().get_functional_description().size() <= i) {
            CoilFunctionalDescription dummyWinding;
            dummyWinding.set_name(get_isolation_side_name_from_index(i));
            dummyWinding.set_number_turns(1);
            dummyWinding.set_number_parallels(1);
            dummyWinding.set_isolation_side(get_isolation_side_from_index(i));
            dummyWinding.set_wire("Dummy");
            magnetic.get_mutable_coil().get_mutable_functional_description().push_back(dummyWinding);
        }
        else {
            if (magnetic.get_coil().get_functional_description()[i].get_name() == "") {
                magnetic.get_mutable_coil().get_mutable_functional_description()[i].set_name(get_isolation_side_name_from_index(i));
            }
            if (magnetic.get_coil().get_functional_description()[i].get_number_turns() == 0) {
                magnetic.get_mutable_coil().get_mutable_functional_description()[i].set_number_turns(1);
            }
            if (magnetic.get_coil().get_functional_description()[i].get_number_parallels() == 0) {
                magnetic.get_mutable_coil().get_mutable_functional_description()[i].set_number_parallels(1);
            }
            if (std::holds_alternative<std::string>(magnetic.get_coil().get_functional_description()[i].get_wire())) {
                auto wireName = std::get<std::string>(magnetic.get_coil().get_functional_description()[i].get_wire());
                if (wireName == "") {
                    magnetic.get_mutable_coil().get_mutable_functional_description()[i].set_wire("Dummy");
                }
            }
        }
    }

    for (size_t i = 0; i < magnetic.get_coil().get_functional_description().size(); i++) {
        auto wire = magnetic.get_mutable_coil().resolve_wire(i);
        auto coating = wire.resolve_coating();
        InsulationWireCoating insulationWireCoating;
        if (coating) {
            insulationWireCoating = wire.resolve_coating().value();
        }
        else {
            insulationWireCoating.set_type(OpenMagnetics::InsulationWireCoatingType::BARE);
        }

        if (!insulationWireCoating.get_material()) {
            auto coatingType = insulationWireCoating.get_type().value();
            if (coatingType == InsulationWireCoatingType::ENAMELLED) {
                insulationWireCoating.set_material(Defaults().defaultEnamelledInsulationMaterial);
            }
            else {
                insulationWireCoating.set_material(Defaults().defaultInsulationMaterial);
            }
        }

        auto insulationWireCoatingMaterial = wire.resolve_coating_insulation_material();
        insulationWireCoating.set_material(insulationWireCoatingMaterial);
        wire.set_coating(insulationWireCoating);

        if (wire.get_strand()) {
            auto strand = wire.resolve_strand();
            wire.set_strand(strand);
        }

        magnetic.get_mutable_coil().get_mutable_functional_description()[i].set_wire(wire);
    }

    BobbinWrapper bobbin = magnetic.get_mutable_coil().resolve_bobbin();

    if (!bobbin.get_functional_description() && !bobbin.get_processed_description()) {
        if (magnetic.get_mutable_core().get_type() == CoreType::TWO_PIECE_SET && magnetic.get_wire(0).get_type() != WireType::RECTANGULAR && magnetic.get_wire(0).get_type() != WireType::PLANAR) {
            bobbin = BobbinWrapper::create_quick_bobbin(magnetic.get_mutable_core(), false);
        }
        else {
            bobbin = BobbinWrapper::create_quick_bobbin(magnetic.get_mutable_core(), true);
        }

    }

    if (!bobbin.get_processed_description()->get_mutable_winding_windows()[0].get_sections_orientation()) {
        if (magnetic.get_mutable_core().get_type() == CoreType::TWO_PIECE_SET) {
            bobbin.get_processed_description()->get_mutable_winding_windows()[0].set_sections_alignment(CoilAlignment::CENTERED);
            bobbin.get_processed_description()->get_mutable_winding_windows()[0].set_sections_orientation(WindingOrientation::OVERLAPPING);
        }
        else {
            bobbin.get_processed_description()->get_mutable_winding_windows()[0].set_sections_alignment(CoilAlignment::SPREAD);
            bobbin.get_processed_description()->get_mutable_winding_windows()[0].set_sections_orientation(WindingOrientation::CONTIGUOUS);
        }
    }
    magnetic.get_mutable_coil().set_bobbin(bobbin);

    if (!magnetic.get_mutable_coil().get_turns_description()) {
        if (magnetic.get_mutable_core().get_type() == CoreType::TWO_PIECE_SET) {
            magnetic.get_mutable_coil().set_turns_alignment(CoilAlignment::SPREAD);
        }
        else {
            magnetic.get_mutable_coil().set_turns_alignment(CoilAlignment::CENTERED);
        }

        magnetic.get_mutable_coil().wind();
    }

    if (magnetic.get_mutable_coil().get_layers_description()) {
        auto layers = magnetic.get_mutable_coil().get_layers_description().value();
        for (size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
            auto insulationMaterial = CoilWrapper::resolve_insulation_layer_insulation_material(magnetic.get_mutable_coil(), layers[layerIndex].get_name());
            layers[layerIndex].set_insulation_material(insulationMaterial);
        }
        magnetic.get_mutable_coil().set_layers_description(layers);
    }

    return magnetic;
}

} // namespace OpenMagnetics
