#include "Defaults.h"
#include "Constants.h"
#include "Utils.h"
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
std::map<std::string, OpenMagnetics::WireWrapper> wireDatabase;
std::map<std::string, OpenMagnetics::BobbinWrapper> bobbinDatabase;
std::map<std::string, OpenMagnetics::InsulationMaterialWrapper> insulationMaterialDatabase;
std::map<std::string, OpenMagnetics::WireMaterial> wireMaterialDatabase;
OpenMagnetics::Constants constants = OpenMagnetics::Constants();

namespace OpenMagnetics {

void load_cores(bool includeToroids, bool useOnlyCoresInStock) {
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
            CoreWrapper core(jf, false, true, false);
            if (includeToroids || core.get_type() != CoreType::TOROIDAL) {
                coreDatabase.push_back(core);
            }
            database.erase(0, pos + delimiter.length());
        }
    }
    else {
        auto data = fs.open("MAS/data/cores.ndjson");
        std::string database = std::string(data.begin(), data.end());
        std::string delimiter = "\n";
        // size_t pos = 0;
        // std::string token;
        // if (database.back() != delimiter.back()) {
        //     database += delimiter;
        // }
        if (database.back() == delimiter.back()) {
            database.pop_back();
        }
        database = std::regex_replace(database, std::regex("\n"), ", ");
        database = "[" + database + "]";
        json arr = json::parse(database);
        coreDatabase = std::vector<CoreWrapper>(arr);

        // throw std::runtime_error("Ea");
        // while ((pos = database.find(delimiter)) != std::string::npos) {
        //     token = database.substr(0, pos);
        //     json jf = json::parse(token);
        //     CoreWrapper core(jf, false, true, false);
        //     if (includeToroids || core.get_type() != CoreType::TOROIDAL) {
        //         coreDatabase.push_back(core);
        //     }
        //     database.erase(0, pos + delimiter.length());
        // }
    }
}

void clear_loaded_cores() {
    coreDatabase.clear();
}

void load_core_materials() {
    auto fs = cmrc::data::get_filesystem();
    {
        auto data = fs.open("MAS/data/core_materials.ndjson");
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
            CoreMaterial coreMaterial(jf);
            coreMaterialDatabase[jf["name"]] = coreMaterial;
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_core_shapes(bool withAliases) {
    auto fs = cmrc::data::get_filesystem();
    {
        auto data = fs.open("MAS/data/core_shapes.ndjson");
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
            CoreShape coreShape(jf);
            coreShapeDatabase[jf["name"]] = coreShape;
            if (withAliases) {
                for (auto& alias : jf["aliases"]) {
                    coreShapeDatabase[alias] = coreShape;
                }
            }
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_wires() {
    auto fs = cmrc::data::get_filesystem();
    {
        auto data = fs.open("MAS/data/wires.ndjson");
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
            WireWrapper wire(jf);
            wireDatabase[jf["name"]] = wire;
            database.erase(0, pos + delimiter.length());
        }
    }
}

void load_bobbins() {
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

void load_databases(json data, bool withAliases) {
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


OpenMagnetics::CoreMaterial find_core_material_by_name(std::string name) {
    if (coreMaterialDatabase.empty()) {
        load_core_materials();
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

std::vector<std::string> get_shape_names(bool includeToroidal) {
    if (coreShapeDatabase.empty()) {
        load_core_shapes(true);
    }

    std::vector<std::string> shapeNames;

    for (auto& datum : coreShapeDatabase) {
        if (includeToroidal || (datum.second.get_family() != CoreShapeFamily::T)) {
            shapeNames.push_back(datum.first);
        }
    }

    return shapeNames;
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

std::vector<OpenMagnetics::WireWrapper> get_wires() {
    if (wireDatabase.empty()) {
        load_wires();
    }

    std::vector<OpenMagnetics::WireWrapper> wires;

    for (auto& datum : wireDatabase) {
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
        // if (std::isnan(sum.real())) {
        //     break;
        // }
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
        // if (std::isnan(sum.real())) {
        //     break;
        // }
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

} // namespace OpenMagnetics
