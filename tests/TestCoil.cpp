#include "RandomUtils.h"
#include <source_location>
#include "support/Settings.h"
#include "support/Painter.h"
#include "constructive_models/Coil.h"
#include "json.hpp"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <limits>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <time.h>
using json = nlohmann::json;
#include <typeinfo>
#include <chrono>
#include <thread>

using namespace MAS;
using namespace OpenMagnetics;
using namespace OpenMagneticsTesting;

namespace {

auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
bool plot = false;

TEST_CASE("Test_Coil_Json_0", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"bobbin":"Dummy","functionalDescription":[{"isolationSide":"Primary","name":"Primary","numberParallels":1,"numberTurns":23,"wire":"Dummy"}]})";

    auto coilJson = json::parse(coilString);
    auto Coil(coilJson);
}

TEST_CASE("Test_Coil_Json_1", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"_interleavingLevel":3,"_windingOrientation":"contiguous","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":1,"numberTurns":9,"wire":"Round 0.475 - Grade 1"}]})";

    auto coil = prepare_coil_from_json(coilString);
    json coilJson;
    to_json(coilJson, coil);

    REQUIRE(coil.get_functional_description().size() > 0);
    REQUIRE(coilJson["functionalDescription"].size() > 0);

    coil.wind();

    auto section = coil.get_sections_description().value()[0];
    REQUIRE(!std::isnan(section.get_dimensions()[0]));
    REQUIRE(!std::isnan(section.get_dimensions()[1]));
}

TEST_CASE("Test_Coil_Json_2", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"_interleavingLevel":7,"_windingOrientation":"overlapping","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":27,"numberTurns":36,"wire":"Round 0.475 - Grade 1"}]})";
    settings.set_coil_wind_even_if_not_fit(false);

    auto coil = prepare_coil_from_json(coilString);
    coil.wind();

    auto section = coil.get_sections_description().value()[0];
    REQUIRE(!std::isnan(section.get_dimensions()[0]));
    REQUIRE(!std::isnan(section.get_dimensions()[1]));
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Json_2.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Test_Coil_Json_3", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"_interleavingLevel":7,"_windingOrientation":"contiguous","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":88,"numberTurns":1,"wire":"Round 0.475 - Grade 1"}]})";
    settings.set_coil_delimit_and_compact(false);

    auto coil = prepare_coil_from_json(coilString);
    coil.wind();

    auto section = coil.get_sections_description().value()[0];
    REQUIRE(!std::isnan(section.get_dimensions()[0]));
    REQUIRE(!std::isnan(section.get_dimensions()[1]));
    std::vector<int64_t> numberTurns = {1};
    std::vector<int64_t> numberParallels = {88};
    uint8_t interleavingLevel = 7;
    check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Json_3.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Test_Coil_Json_4", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.006,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0032500000000000003,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0002835287369864788,"coordinates":[0.0095,0,0],"height":null,"radialHeight":0.0095,"sectionsAlignment":"outer or bottom","sectionsOrientation":"contiguous","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":27,"wire":{"coating":{"breakdownVoltage":2700,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":4.116868676970209e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000724},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 21.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000757},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"21 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"secondary","name":"Secondary","numberParallels":1,"numberTurns":27,"wire":{"coating":{"breakdownVoltage":5000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":4.620411001469214e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000767},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 20.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000831},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"20.5 AWG","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription": null, "turnsDescription":null,"_turnsAlignment":{"Primary section 0":"spread","Secondary section 0":"spread"},"_layersOrientation":{"Primary section 0":"overlapping","Secondary section 0":"overlapping"}})";

    std::vector<size_t> pattern = {0, 1};
    std::vector<double> proportionPerWinding = {0.5, 0.5};
    size_t repetitions = 2;

    auto coilJson = json::parse(coilString);

    auto coilFunctionalDescription = coilJson["functionalDescription"].get<std::vector<OpenMagnetics::Winding>>();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(coilJson["bobbin"]);
    coil.set_functional_description(coilFunctionalDescription);

    if (coilJson["_layersOrientation"].is_object()) {
        auto layersOrientationPerSection = std::map<std::string, WindingOrientation>(coilJson["_layersOrientation"]);
        for (auto [sectionName, layerOrientation] : layersOrientationPerSection) {
            coil.set_layers_orientation(layerOrientation, sectionName);
        }
    }
    else if (coilJson["_layersOrientation"].is_array()) {
        coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
        if (coil.get_sections_description()) {
            auto sections = coil.get_sections_description_conduction();
            auto layersOrientationPerSection = coilJson["_layersOrientation"].get<std::vector<WindingOrientation>>();
            for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                if (sectionIndex < layersOrientationPerSection.size()) {
                    coil.set_layers_orientation(layersOrientationPerSection[sectionIndex], sections[sectionIndex].get_name());
                }
            }
        }
    }
    else {
        WindingOrientation layerOrientation(coilJson["_layersOrientation"]);
        coil.set_layers_orientation(layerOrientation);

    }
    if (coilJson["_turnsAlignment"].is_object()) {
        auto turnsAlignmentPerSection = std::map<std::string, CoilAlignment>(coilJson["_turnsAlignment"]);
        for (auto [sectionName, turnsAlignment] : turnsAlignmentPerSection) {
            coil.set_turns_alignment(turnsAlignment, sectionName);
        }
    }
    else if (coilJson["_turnsAlignment"].is_array()) {
        coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
        if (coil.get_sections_description()) {
            auto sections = coil.get_sections_description_conduction();
            auto turnsAlignmentPerSection = coilJson["_turnsAlignment"].get<std::vector<CoilAlignment>>();
            for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                if (sectionIndex < turnsAlignmentPerSection.size()) {
                    coil.set_turns_alignment(turnsAlignmentPerSection[sectionIndex], sections[sectionIndex].get_name());
                }
            }
        }
    }
    else {
        CoilAlignment turnsAlignment(coilJson["_turnsAlignment"]);
        coil.set_turns_alignment(turnsAlignment);
    }

    coil.set_bobbin(coilJson["bobbin"]);
    coil.set_functional_description(coilFunctionalDescription);
    coil.wind();
    REQUIRE(bool(coil.get_sections_description()));
    REQUIRE(bool(coil.get_layers_description()));
    REQUIRE(bool(coil.get_turns_description()));
}

TEST_CASE("Test_Coil_Json_5", "[constructive-model][coil][bug][smoke-test]") {
    auto json_path_180 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_json_5_180.json");
    std::ifstream json_file_180(json_path_180);
    std::string coilString((std::istreambuf_iterator<char>(json_file_180)), std::istreambuf_iterator<char>());

    std::vector<size_t> pattern = {0, 1};
    std::vector<double> proportionPerWinding = {0.25, 0.75};
    size_t repetitions = 2;

    auto coilJson = json::parse(coilString);

    auto coilFunctionalDescription = coilJson["functionalDescription"].get<std::vector<OpenMagnetics::Winding>>();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(coilJson["bobbin"]);
    coil.set_functional_description(coilFunctionalDescription);

    if (coilJson["_layersOrientation"].is_object()) {
        auto layersOrientationPerSection = std::map<std::string, WindingOrientation>(coilJson["_layersOrientation"]);
        for (auto [sectionName, layerOrientation] : layersOrientationPerSection) {
            coil.set_layers_orientation(layerOrientation, sectionName);
        }
    }
    else if (coilJson["_layersOrientation"].is_array()) {
        coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
        if (coil.get_sections_description()) {
            auto sections = coil.get_sections_description_conduction();
            auto layersOrientationPerSection = coilJson["_layersOrientation"].get<std::vector<WindingOrientation>>();
            for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                if (sectionIndex < layersOrientationPerSection.size()) {
                    coil.set_layers_orientation(layersOrientationPerSection[sectionIndex], sections[sectionIndex].get_name());
                }
            }
        }
    }
    else {
        WindingOrientation layerOrientation(coilJson["_layersOrientation"]);
        coil.set_layers_orientation(layerOrientation);

    }
    if (coilJson["_turnsAlignment"].is_object()) {
        auto turnsAlignmentPerSection = std::map<std::string, CoilAlignment>(coilJson["_turnsAlignment"]);
        for (auto [sectionName, turnsAlignment] : turnsAlignmentPerSection) {
            coil.set_turns_alignment(turnsAlignment, sectionName);
        }
    }
    else if (coilJson["_turnsAlignment"].is_array()) {
        coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
        if (coil.get_sections_description()) {
            auto sections = coil.get_sections_description_conduction();
            auto turnsAlignmentPerSection = coilJson["_turnsAlignment"].get<std::vector<CoilAlignment>>();
            for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                if (sectionIndex < turnsAlignmentPerSection.size()) {
                    coil.set_turns_alignment(turnsAlignmentPerSection[sectionIndex], sections[sectionIndex].get_name());
                }
            }
        }
    }
    else {
        CoilAlignment turnsAlignment(coilJson["_turnsAlignment"]);
        coil.set_turns_alignment(turnsAlignment);
    }

    coil.set_bobbin(coilJson["bobbin"]);
    coil.set_functional_description(coilFunctionalDescription);
    coil.wind();
    REQUIRE(bool(coil.get_sections_description()));
    REQUIRE(bool(coil.get_layers_description()));
    REQUIRE(bool(coil.get_turns_description()));
    REQUIRE(bool(coil.are_sections_and_layers_fitting()));
}

TEST_CASE("Test_Coil_Json_6", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"_sectionsAlignment":"spread","_turnsAlignment":"centered","bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.0075,"columnShape":"rectangular","columnThickness":0.0,"columnWidth":0.0026249999999999997,"coordinates":[0.0,0.0,0.0],"pins":null,"wallThickness":0.0,"windingWindows":[{"angle":360.0,"area":0.00017203361371057708,"coordinates":[0.0074,0.0,0.0],"height":null,"radialHeight":0.0074,"sectionsAlignment":"spread","sectionsOrientation":"contiguous","shape":"round","width":null}]}},"functionalDescription":[{"isolationSide":"primary","name":"primary","numberParallels":1,"numberTurns":15,"wire":{"coating":{"breakdownVoltage":2700.0,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.00125},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Round 1.25 - Grade 1","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001316},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"1.25 mm","strand":null,"type":"round"}},{"isolationSide":"secondary","name":"secondary","numberParallels":1,"numberTurns":15,"wire":{"coating":{"breakdownVoltage":2700.0,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.00125},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Round 1.25 - Grade 1","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001316},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"1.25 mm","strand":null,"type":"round"}}]})";

    CoilWindingConfig config;
    config.coilJsonStr = coilString;
    config.pattern = {0, 1};
    config.repetitions = 1;

    auto coil = prepare_and_wind_coil(config);

    REQUIRE(coil.get_turns_description().has_value());

    json result;
    to_json(result, coil);
}

TEST_CASE("Test_Coil_Json_7", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.00356,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0022725,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0000637587014444212,"coordinates":[0.004505,0,0],"height":null,"radialHeight":0.004505,"sectionsAlignment":"inner or top","sectionsOrientation":"overlapping","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":3,"numberTurns":55,"wire":{"coating":{"breakdownVoltage":1220,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":8.042477193189871e-8},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000323,"minimum":0.00031800000000000003,"nominal":0.00032},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 28.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000356,"minimum":0.00033800000000000003,"nominal":0.000347},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"28 AWG","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null,"_turnsAlignment":["spread"],"_layersOrientation":["overlapping"]})";

    CoilWindingConfig config;
    config.coilJsonStr = coilString;
    config.pattern = {0};
    config.proportionPerWinding = {1.0};
    config.repetitions = 1;

    auto coil = prepare_and_wind_coil(config);

    REQUIRE(coil.get_turns_description().has_value());

    json result;
    to_json(result, coil);
}

TEST_CASE("Test_Coil_Json_8", "[constructive-model][coil][bug][smoke-test]") {
    auto json_path_282 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_json_8_282.json");
    std::ifstream json_file_282(json_path_282);
    json coilJson = json::parse(json_file_282);
    OpenMagnetics::Coil coil(coilJson, false);
    auto layers = coil.get_layers_description().value();

    for (auto layer : layers) {
        if (layer.get_type() == ElectricalType::INSULATION) {
            auto material = OpenMagnetics::Coil::resolve_insulation_layer_insulation_material(coil, layer.get_name());
            json mierda;
            to_json(mierda, material);
        }

    }
}

TEST_CASE("Test_Coil_Json_9", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.01295, "columnShape": "round", "columnThickness": 0.0, "columnWidth": 0.01295, "coordinates": [0.0, 0.0, 0.0 ], "pins": null, "wallThickness": 0.0, "windingWindows": [{"angle": null, "area": 0.0001596, "coordinates": [0.0196, 0.0, 0.0 ], "height": 0.012, "radialHeight": null, "sectionsAlignment": "centered", "sectionsOrientation": "contiguous", "shape": "rectangular", "width": 0.0133 } ] } }, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "Primary", "numberParallels": 3, "numberTurns": 12, "wire": {"coating": {"breakdownVoltage": null, "grade": null, "material": {"aliases": ["Tefzel ETFE" ], "composition": "Ethylene Tetrafluoroethylene", "dielectricStrength": [{"humidity": null, "temperature": 23.0, "thickness": 2.5e-05, "value": 160000000.0 } ], "manufacturer": "Chemours", "meltingPoint": 220.0, "name": "ETFE", "relativePermittivity": 2.7, "resistivity": [{"temperature": 170.0, "value": 1000000000000000.0 } ], "specificHeat": 1172.0, "temperatureClass": 155.0, "thermalConductivity": 0.24 }, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "bare" }, "conductingArea": null, "conductingDiameter": null, "conductingHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00020999999999999998 }, "conductingWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002 }, "edgeRadius": null, "manufacturerInfo": null, "material": "copper", "name": null, "numberConductors": 1, "outerDiameter": null, "outerHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00021020999999999995 }, "outerWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002002 }, "standard": null, "standardName": null, "strand": null, "type": "rectangular" } }, {"connections": null, "isolationSide": "secondary", "name": "Secondary", "numberParallels": 3, "numberTurns": 15, "wire": {"coating": {"breakdownVoltage": null, "grade": null, "material": {"aliases": [], "composition": "Polyurethane", "dielectricStrength": [{"humidity": null, "temperature": null, "thickness": 0.0001, "value": 25000000.0 } ], "manufacturer": "MWS Wire", "meltingPoint": null, "name": "Polyurethane 155", "relativePermittivity": 3.7, "resistivity": [{"temperature": null, "value": 1e+16 } ], "specificHeat": null, "temperatureClass": 155.0, "thermalConductivity": null }, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled" }, "conductingArea": null, "conductingDiameter": null, "conductingHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00020999999999999998 }, "conductingWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002 }, "edgeRadius": null, "manufacturerInfo": null, "material": "copper", "name": null, "numberConductors": 1, "outerDiameter": null, "outerHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00021020999999999995 }, "outerWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002002 }, "standard": null, "standardName": null, "strand": null, "type": "rectangular" } } ], "layersOrientation": "contiguous", "turnsAlignment": "spread" })";
    auto coilJson = json::parse(coilString);
    size_t repetitions = 1;
    double insulationThickness = 0.10 / 1000;
    std::string proportionPerWindingString = "[]";
    std::string patternString = "[0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]";

    std::vector<double> proportionPerWinding = json::parse(proportionPerWindingString);
    std::vector<size_t> pattern = json::parse(patternString);
    auto coilFunctionalDescription = coilJson["functionalDescription"].get<std::vector<OpenMagnetics::Winding>>();
    OpenMagnetics::Coil coil;

    if (coilJson.contains("_interleavingLevel")) {
        coil.set_interleaving_level(coilJson["_interleavingLevel"]);
    }
    if (coilJson.contains("_windingOrientation")) {
        coil.set_winding_orientation(coilJson["_windingOrientation"]);
    }
    if (coilJson.contains("_layersOrientation")) {
        coil.set_layers_orientation(coilJson["_layersOrientation"]);
    }
    if (coilJson.contains("_turnsAlignment")) {
        coil.set_turns_alignment(coilJson["_turnsAlignment"]);
    }
    if (coilJson.contains("_sectionAlignment")) {
        coil.set_section_alignment(coilJson["_sectionAlignment"]);
    }

    coil.set_bobbin(coilJson["bobbin"]);
    coil.set_functional_description(coilFunctionalDescription);

    if (insulationThickness > 0) {
        coil.calculate_custom_thickness_insulation(insulationThickness);
    }
    if (proportionPerWinding.size() == coilFunctionalDescription.size()) {
        if (pattern.size() > 0 && repetitions > 0) {
            coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
        }
        else if (repetitions > 0) {
            coil.wind_by_sections(repetitions);
        }
        else {
            coil.wind_by_sections();
        }
    }
    else {
        if (pattern.size() > 0 && repetitions > 0) {
            coil.wind_by_sections(pattern, repetitions);
        }
        else if (repetitions > 0) {
            coil.wind_by_sections(repetitions);
        }
        else {
            coil.wind_by_sections();
        }
    }

    json result;
    to_json(result, coil);
}

TEST_CASE("Test_Coil_Json_10", "[constructive-model][coil][bug][smoke-test]") {
    auto json_path_359 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_json_10_359.json");
    std::ifstream json_file_359(json_path_359);
    std::string coilString((std::istreambuf_iterator<char>(json_file_359)), std::istreambuf_iterator<char>());
    auto coilJson = json::parse(coilString);

    auto coilFunctionalDescription = coilJson["functionalDescription"].get<std::vector<OpenMagnetics::Winding>>();
    auto coilSectionsDescription = coilJson["sectionsDescription"].get<std::vector<Section>>();
    auto coilLayersDescription = coilJson["layersDescription"].get<std::vector<Layer>>();
    OpenMagnetics::Coil coil;

    if (coilJson.contains("_interleavingLevel")) {
        coil.set_interleaving_level(coilJson["_interleavingLevel"]);
    }
    if (coilJson.contains("_windingOrientation")) {
        coil.set_winding_orientation(coilJson["_windingOrientation"]);
    }
    if (coilJson.contains("_layersOrientation")) {
        coil.set_layers_orientation(coilJson["_layersOrientation"]);
    }
    if (coilJson.contains("_turnsAlignment")) {
        coil.set_turns_alignment(coilJson["_turnsAlignment"]);
    }
    if (coilJson.contains("_sectionAlignment")) {
        coil.set_section_alignment(coilJson["_sectionAlignment"]);
    }

    coil.set_bobbin(coilJson["bobbin"]);
    coil.set_functional_description(coilFunctionalDescription);
    coil.set_sections_description(coilSectionsDescription);
    coil.set_layers_description(coilLayersDescription);
    coil.wind_by_turns();

    json result;
    to_json(result, coil);
}

TEST_CASE("Test_Coil_Json_11", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.006175,"columnShape":"round","columnThickness":0,"columnWidth":0.006175,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":null,"area":0.000041283000000000004,"coordinates":[0.0098875,0,0],"height":0.00556,"radialHeight":null,"sectionsAlignment":"inner or top","sectionsOrientation":"contiguous","shape":"rectangular","width":0.007425000000000001}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"primary","numberParallels":1,"numberTurns":5,"wire":{"coating":null,"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.2293100000000003e-7},"conductingDiameter":null,"conductingHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000522},"conductingWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0023550000000000008},"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Planar 52.20 µm","numberConductors":1,"outerDiameter":null,"outerHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000522},"outerWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0023550000000000008},"standard":"IPC-6012","standardName":"1.5 oz.","strand":null,"type":"planar"}},{"connections":null,"isolationSide":"secondary","name":"SECONDARY","numberParallels":1,"numberTurns":3,"wire":{"coating":null,"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.2449700000000003e-7},"conductingDiameter":null,"conductingHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000348},"conductingWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.003577500000000001},"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Planar 34.80 µm","numberConductors":1,"outerDiameter":null,"outerHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000348},"outerWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.003577500000000001},"standard":"IPC-6012","standardName":"1 oz.","strand":null,"type":"planar"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null,"_turnsAlignment":["spread","spread","spread","spread"],"_layersOrientation":["contiguous","contiguous","contiguous","contiguous"],"_interlayerInsulationThickness":0,"_intersectionInsulationThickness":0.0001})";

    std::vector<size_t> pattern = {0, 1, 0, 1};
    std::vector<double> proportionPerWinding = {0.5, 0.5};
    size_t repetitions = 1;

    auto coilJson = json::parse(coilString);

    auto coilFunctionalDescription = coilJson["functionalDescription"].get<std::vector<OpenMagnetics::Winding>>();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(coilJson["bobbin"]);
    coil.set_functional_description(coilFunctionalDescription);

    if (coilJson["_layersOrientation"].is_object()) {
        auto layersOrientationPerSection = std::map<std::string, WindingOrientation>(coilJson["_layersOrientation"]);
        for (auto [sectionName, layerOrientation] : layersOrientationPerSection) {
            coil.set_layers_orientation(layerOrientation, sectionName);
        }
    }
    else if (coilJson["_layersOrientation"].is_array()) {
        coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
        if (coil.get_sections_description()) {
            auto sections = coil.get_sections_description_conduction();
            auto layersOrientationPerSection = coilJson["_layersOrientation"].get<std::vector<WindingOrientation>>();
            for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                if (sectionIndex < layersOrientationPerSection.size()) {
                    coil.set_layers_orientation(layersOrientationPerSection[sectionIndex], sections[sectionIndex].get_name());
                }
            }
        }
    }
    else {
        WindingOrientation layerOrientation(coilJson["_layersOrientation"]);
        coil.set_layers_orientation(layerOrientation);

    }
    if (coilJson["_turnsAlignment"].is_object()) {
        auto turnsAlignmentPerSection = std::map<std::string, CoilAlignment>(coilJson["_turnsAlignment"]);
        for (auto [sectionName, turnsAlignment] : turnsAlignmentPerSection) {
            coil.set_turns_alignment(turnsAlignment, sectionName);
        }
    }
    else if (coilJson["_turnsAlignment"].is_array()) {
        coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
        if (coil.get_sections_description()) {
            auto sections = coil.get_sections_description_conduction();
            auto turnsAlignmentPerSection = coilJson["_turnsAlignment"].get<std::vector<CoilAlignment>>();
            for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                if (sectionIndex < turnsAlignmentPerSection.size()) {
                    coil.set_turns_alignment(turnsAlignmentPerSection[sectionIndex], sections[sectionIndex].get_name());
                }
            }
        }
    }
    else {
        CoilAlignment turnsAlignment(coilJson["_turnsAlignment"]);
        coil.set_turns_alignment(turnsAlignment);
    }

    coil.set_bobbin(coilJson["bobbin"]);
    coil.set_functional_description(coilFunctionalDescription);
    coil.wind();
    REQUIRE(bool(coil.get_sections_description()));
    REQUIRE(bool(coil.get_layers_description()));
    REQUIRE(bool(coil.get_turns_description()));
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered", "[constructive-model][coil][margin][smoke-test]") {
    settings.reset();
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[1] == sectionDimensionsAfterMarginNoFill[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth, 0.001));
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 25, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.001;
    
    settings.set_coil_fill_sections_with_margin_tape(false);
    // settings.set_coil_wind_even_if_not_fit(false);
    // settings.set_coil_try_rewind(true);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Centered_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_1[1] > sectionDimensionsAfterMarginNoFill_1[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Top", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE_THAT(sectionDimensionsAfterMarginFill[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill[1], 0.0001));
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth, 0.0001));
    REQUIRE_THAT(marginAfterMarginFill[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill[0], 0.0001));
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Top_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 25, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsBeforeMargin_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginBeforeMargin_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Top_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginAfterMarginNoFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Top_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginAfterMarginFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Top_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_2[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(0 == marginBeforeMargin_2[1]);
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_0[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_0[1], 0.0001));
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_1[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_1[1], 0.0001));
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_2[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_2[1], 0.0001));
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.0001));
    REQUIRE_THAT(marginAfterMarginFill_0[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_0[0], 0.0001));
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE_THAT(marginAfterMarginFill_2[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_2[0], 0.0001));
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_2[1] > marginAfterMarginNoFill_2[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE_THAT(sectionDimensionsBeforeMargin_1[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_1[1], 0.0001));
    REQUIRE_THAT(sectionDimensionsBeforeMargin_2[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_2[1], 0.0001));
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Bottom", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE_THAT(sectionDimensionsAfterMarginFill[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill[1], 0.0001));
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth, 0.0001));
    REQUIRE_THAT(marginAfterMarginFill[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill[1], 0.0001));
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Bottom_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 25, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsBeforeMargin_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginBeforeMargin_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Bottom_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginAfterMarginNoFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Bottom_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginAfterMarginFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Bottom_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_2[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(0 == marginBeforeMargin_2[1]);
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_0[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_0[1], 0.0001));
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_1[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_1[1], 0.0001));
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_2[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_2[1], 0.0001));
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.0001));
    REQUIRE_THAT(marginAfterMarginFill_0[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_0[1], 0.0001));
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    REQUIRE_THAT(marginAfterMarginFill_2[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_2[1], 0.0001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
    REQUIRE(marginAfterMarginFill_2[0] > marginAfterMarginNoFill_2[0]);
    REQUIRE(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE_THAT(sectionDimensionsBeforeMargin_1[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_1[1], 0.0001));
    REQUIRE_THAT(sectionDimensionsBeforeMargin_2[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_2[1], 0.0001));
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Spread", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE_THAT(sectionDimensionsAfterMarginFill[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill[1], 0.0001));
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth, 0.0001));
    REQUIRE_THAT(marginAfterMarginFill[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill[1], 0.0001));
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Spread_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 25, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsBeforeMargin_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginBeforeMargin_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginAfterMarginNoFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginAfterMarginFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_2[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(0 == marginBeforeMargin_2[1]);
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_0[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_0[1], 0.0001));
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_1[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_1[1], 0.0001));
    REQUIRE_THAT(sectionDimensionsAfterMarginFill_2[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_2[1], 0.0001));
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.0001));
    REQUIRE_THAT(marginAfterMarginFill_0[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_0[1], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE_THAT(marginAfterMarginFill_2[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_2[1], 0.0001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
    REQUIRE(marginAfterMarginFill_2[0] > marginAfterMarginNoFill_2[0]);
    REQUIRE(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE_THAT(sectionDimensionsBeforeMargin_1[1], Catch::Matchers::WithinAbs(sectionDimensionsAfterMarginNoFill_1[1], 0.0001));
    REQUIRE(sectionDimensionsBeforeMargin_2[1] > sectionDimensionsAfterMarginNoFill_2[1]);
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Inner_No_Filling_Then_Filling_Horizontal_Centered", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Inner_No_Filling_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Inner_No_Filling_Then_Filling_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[1] == sectionDimensionsAfterMarginNoFill[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth, 0.001));
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Inner_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 25, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.001;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Inner_No_Filling_Horizontal_Centered_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Inner_No_Filling_Horizontal_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Inner_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_1[1] > sectionDimensionsAfterMarginNoFill_1[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Outer_No_Filling_Then_Filling_Horizontal_Centered", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowEndingWidth = windingWindowCoordinates[0] + windingWindowDimensions[0] / 2;
    auto sectionEndingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] + coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Outer_No_Filling_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Outer_No_Filling_Then_Filling_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[1] == sectionDimensionsAfterMarginNoFill[1]);
    REQUIRE_THAT(windingWindowEndingWidth, Catch::Matchers::WithinAbs(sectionEndingWidth, 0.001));
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Outer_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 25, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.001;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Outer_No_Filling_Horizontal_Centered_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Outer_No_Filling_Horizontal_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Outer_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_1[1] > sectionDimensionsAfterMarginNoFill_1[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Horizontal_Centered", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.002;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[1] == sectionDimensionsAfterMarginNoFill[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth, 0.001));
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 25, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.001;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Horizontal_Centered_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Horizontal_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_1[1] > sectionDimensionsAfterMarginNoFill_1[1]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Centered", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Centered_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginAfterMarginNoFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto marginAfterMarginFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
    REQUIRE_THAT(marginAfterMarginFill_2[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_2[1], 0.0001));
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Top", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Top_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Top_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Top_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Top_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Bottom", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Bottom_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Bottom_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Bottom_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Bottom_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Spread", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Spread_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Centered", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Centered_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 3, margin * 0.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Inner", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Inner_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Inner_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Outer", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Outer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Outer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Outer_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Outer_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Spread", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Spread_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Centered", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Centered_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Inner", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Inner_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Inner_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Outer", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Outer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Outer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Outer_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Outer_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Spread", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Spread_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Centered", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Centered_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Inner", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Inner_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Inner_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Outer", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Outer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Outer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Outer_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Outer_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Spread", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {47};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginBeforeMargin = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
    auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin[0]);
    REQUIRE(0 == marginBeforeMargin[1]);
    REQUIRE(sectionDimensionsAfterMarginFill[0] == sectionDimensionsAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
    REQUIRE(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
    REQUIRE(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins", "[constructive-model][coil][margin][smoke-test]") {
    std::vector<int64_t> numberTurns = {34, 12, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    double margin = 0.0005;
    
    settings.set_coil_fill_sections_with_margin_tape(false);

    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);
    auto sectionDimensionsBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginBeforeMargin_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginBeforeMargin_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Spread_Three_Different_Margins_No_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_fill_sections_with_margin_tape(false);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
    auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginNoFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginNoFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
    auto bobbin = coil.resolve_bobbin();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    OpenMagneticsTesting::check_turns_description(coil);

    settings.set_coil_fill_sections_with_margin_tape(true);
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
    auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
    auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
    auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
    auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    REQUIRE(0 == marginBeforeMargin_0[0]);
    REQUIRE(0 == marginBeforeMargin_0[1]);
    REQUIRE(0 == marginBeforeMargin_1[0]);
    REQUIRE(0 == marginBeforeMargin_1[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_0[1] == sectionDimensionsAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsAfterMarginFill_1[1] == sectionDimensionsAfterMarginNoFill_1[1]);
    REQUIRE_THAT(windingWindowStartingWidth, Catch::Matchers::WithinAbs(sectionStartingWidth_0, 0.001));
    REQUIRE(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
    REQUIRE(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
    REQUIRE(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
    REQUIRE(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
    REQUIRE(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Wind_By_Section_Wind_By_Consecutive_Parallels", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {42};
    std::vector<int64_t> numberParallels = {3};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 2;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
}

TEST_CASE("Test_Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {41};
    std::vector<int64_t> numberParallels = {3};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 2;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
}

TEST_CASE("Test_Wind_By_Section_Wind_By_Full_Turns", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {2};
    std::vector<int64_t> numberParallels = {7};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 2;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
}

TEST_CASE("Test_Wind_By_Section_Wind_By_Full_Parallels", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {2};
    std::vector<int64_t> numberParallels = {7};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 7;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
}

TEST_CASE("Test_Wind_By_Section_Wind_By_Full_Parallels_Multiwinding", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {2, 5};
    std::vector<int64_t> numberParallels = {7, 7};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 7;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
}

TEST_CASE("Test_Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced_Vertical", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {41};
    std::vector<int64_t> numberParallels = {3};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 2;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
}

TEST_CASE("Test_Wind_By_Section_Random_0", "[constructive-model][coil][rectangular-winding-window]") {
    std::vector<int64_t> numberTurns = {9};
    std::vector<int64_t> numberParallels = {1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 3;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
}

TEST_CASE("Test_Wind_By_Section_Random_1", "[constructive-model][coil][rectangular-winding-window]") {
    std::vector<int64_t> numberTurns = {6};
    std::vector<int64_t> numberParallels = {2};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 3;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
}

TEST_CASE("Test_Wind_By_Section_Random_2", "[constructive-model][coil][rectangular-winding-window]") {
    std::vector<int64_t> numberTurns = {5};
    std::vector<int64_t> numberParallels = {2};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 3;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
}

TEST_CASE("Test_Wind_By_Section_Random_3", "[constructive-model][coil][rectangular-winding-window]") {
    std::vector<int64_t> numberTurns = {5};
    std::vector<int64_t> numberParallels = {1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 3;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
}

TEST_CASE("Test_Wind_By_Section_Random_4", "[constructive-model][coil][rectangular-winding-window]") {
    std::vector<int64_t> numberTurns = {91};
    std::vector<int64_t> numberParallels = {2};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 3;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
}

TEST_CASE("Test_Wind_By_Section_Random_5", "[constructive-model][coil][rectangular-winding-window]") {
    std::vector<int64_t> numberTurns = {23};
    std::vector<int64_t> numberParallels = {1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 7;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
}

TEST_CASE("Test_Wind_By_Section_Random_6", "[constructive-model][coil][rectangular-winding-window]") {
    std::vector<int64_t> numberTurns = {1};
    std::vector<int64_t> numberParallels = {43};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 5;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
}

TEST_CASE("Test_Wind_By_Section_Random", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_try_rewind(false);
    for (size_t i = 0; i < 1000; ++i)
    {
        std::vector<int64_t> numberTurns = {OpenMagnetics::TestUtils::randomInt64(1, 100 + 1 - 1)};
        std::vector<int64_t> numberParallels = {OpenMagnetics::TestUtils::randomInt64(1, 100 + 1 - 1)};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        int64_t numberPhysicalTurns = numberTurns[0] * numberParallels[0];
        uint8_t interleavingLevel = uint8_t(OpenMagnetics::TestUtils::randomInt(1, 10 + static_cast<int>(1) - 1));
        interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
        auto windingOrientation = OpenMagnetics::TestUtils::randomInt(0, 2 - 1)? WindingOrientation::CONTIGUOUS : WindingOrientation::OVERLAPPING;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, windingOrientation);
    }
    settings.reset();
}

TEST_CASE("Test_Wind_By_Section_Random_Multiwinding", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_try_rewind(false);
    for (size_t i = 0; i < 1000; ++i)
    {
        std::vector<int64_t> numberTurns;
        std::vector<int64_t> numberParallels;
        int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
        for (size_t windingIndex = 0; windingIndex < OpenMagnetics::TestUtils::randomSize(1, 10 + 1 - 1); ++windingIndex)
        {
            numberTurns.push_back(OpenMagnetics::TestUtils::randomInt64(1, 100 + 1 - 1));
            numberParallels.push_back(OpenMagnetics::TestUtils::randomInt64(1, 100 + 1 - 1));
            numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
        }
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        int64_t interleavingLevel = OpenMagnetics::TestUtils::randomInt(1, 10 + static_cast<int>(1) - 1);
        interleavingLevel = std::min(numberPhysicalTurns, interleavingLevel);
        auto windingOrientation = OpenMagnetics::TestUtils::randomInt(0, 2 - 1)? WindingOrientation::CONTIGUOUS : WindingOrientation::OVERLAPPING;
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, windingOrientation);
    }
    settings.reset();
}

TEST_CASE("Test_Wind_By_Section_With_Insulation_Sections", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {23, 42};
    std::vector<int64_t> numberParallels = {2, 1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
    uint8_t interleavingLevel = 2;

    auto wires = std::vector<OpenMagnetics::Wire>({find_wire_by_name("Round 0.014 - Grade 1")});

    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
    double voltagePeakToPeak = 400;
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
    coil.set_inputs(inputs);
    coil.wind();
    auto log = coil.read_log();

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, sectionOrientation);
}

TEST_CASE("Test_Wind_By_Section_Pattern", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {21, 21};
    std::vector<int64_t> numberParallels = {2, 2};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 2;

    std::vector<size_t> pattern = {0, 1};
    size_t repetitions = 2;

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

    coil.wind_by_sections(pattern, repetitions);
    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
}

TEST_CASE("Test_Wind_By_Layers_Wind_One_Section_One_Layer", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 9;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
    auto layersDescription = coil.get_layers_description().value();
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Wind_One_Section_Two_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 6;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
    auto layersDescription = coil.get_layers_description().value();
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Wind_One_Section_One_Layer_Two_Parallels", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {2};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 15;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
    auto layersDescription = coil.get_layers_description().value();
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Wind_One_Section_Two_Layers_Two_Parallels", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {2};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 6;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
    auto layersDescription = coil.get_layers_description().value();
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Two_Sections_Two_Layers_Two_Parallels", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {2};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 6;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 2;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
    auto layersDescription = coil.get_layers_description().value();
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Two_Sections_One_Layer_One_Parallel", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 6;
    int64_t numberMaximumLayers = 1;
    uint8_t interleavingLevel = 2;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
    auto layersDescription = coil.get_layers_description().value();
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Two_Sections_One_Layer_Two_Parallels", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {2};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 6;
    int64_t numberMaximumLayers = 1;
    uint8_t interleavingLevel = 2;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
    auto layersDescription = coil.get_layers_description().value();
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Two_Sections_Two_Layers_One_Parallel", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 2;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 2;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
 
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Vertical_Winding_Horizontal_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 2;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    auto layersOrientation = WindingOrientation::CONTIGUOUS;
    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
    OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Vertical_Winding_Vertical_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 2;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    auto layersOrientation = WindingOrientation::CONTIGUOUS;
    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
 
    OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Horizontal_Winding_Horizontal_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 2;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    auto layersOrientation = WindingOrientation::OVERLAPPING;
    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
 
    OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Horizontal_Winding_Vertical_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 2;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    auto layersOrientation = WindingOrientation::CONTIGUOUS;
    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
 
    OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
}

TEST_CASE("Test_Wind_By_Layers_Wind_Horizontal_Winding", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 2;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

    auto windingOrientation = WindingOrientation::OVERLAPPING;
    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
 
    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Random_0", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {5};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 1;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 2;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Random", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    for (size_t i = 0; i < 1000; ++i)
    {
        std::vector<int64_t> numberTurns = {OpenMagnetics::TestUtils::randomInt64(1, 10 + 1 - 1)};
        std::vector<int64_t> numberParallels = {OpenMagnetics::TestUtils::randomInt64(1, 3 + 1 - 1)};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = OpenMagnetics::TestUtils::randomInt64(1, 4 + 1 - 1);
        int64_t numberMaximumLayers = OpenMagnetics::TestUtils::randomInt64(1, 3 + 1 - 1);
        uint8_t interleavingLevel = OpenMagnetics::TestUtils::randomInt(1, 10 + static_cast<int>(1) - 1);
        int64_t numberPhysicalTurns = numberTurns[0] * numberParallels[0];
        interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        OpenMagneticsTesting::check_layers_description(coil);
    }
}

TEST_CASE("Test_Wind_By_Layers_With_Insulation_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {23, 42};
    std::vector<int64_t> numberParallels = {2, 1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
    uint8_t interleavingLevel = 2;

    auto wires = std::vector<OpenMagnetics::Wire>({find_wire_by_name("Round 0.014 - Grade 1")});

    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
    double voltagePeakToPeak = 400;
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
    coil.set_inputs(inputs);
    coil.wind();
    auto log = coil.read_log();

    OpenMagneticsTesting::check_layers_description(coil);
}

TEST_CASE("Test_External_Insulation_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {

    std::string insulationLayersString = "{\"(0, 1)\":[{\"coordinates\":[0.0035501599999999997,0],\"dimensions\":[2.032e-05,0.0102],\"orientation\":\"overlapping\",\"margin\":[0,0],\"name\":\"section_1_insulation_layer_0\",\"partialWindings\":[],\"type\":\"insulation\"},{\"coordinates\":[0.00709016,0],\"dimensions\":[2.032e-05,0.0102],\"orientation\":\"overlapping\",\"margin\":[0,0],\"name\":\"section_1_insulation_layer_1\",\"partialWindings\":[],\"type\":\"insulation\"}],\"(1, 2)\":[{\"coordinates\":[0.004212799999998001,0],\"dimensions\":[2.032e-05,0.0102],\"orientation\":\"overlapping\",\"margin\":[0,0],\"name\":\"section_3_insulation_layer_0\",\"partialWindings\":[],\"type\":\"insulation\"},{\"coordinates\":[0.008415439999996001,0],\"dimensions\":[2.032e-05,0.0102],\"orientation\":\"overlapping\",\"margin\":[0,0],\"name\":\"section_3_insulation_layer_1\",\"partialWindings\":[],\"type\":\"insulation\"}],\"(2, 0)\":[{\"coordinates\":[0.004423439999998001,0],\"dimensions\":[2.032e-05,0.0102],\"orientation\":\"overlapping\",\"margin\":[0,0],\"name\":\"section_5_insulation_layer_0\",\"partialWindings\":[],\"type\":\"insulation\"},{\"coordinates\":[0.008836719999996,0],\"dimensions\":[2.032e-05,0.0102],\"orientation\":\"overlapping\",\"margin\":[0,0],\"name\":\"section_5_insulation_layer_1\",\"partialWindings\":[],\"type\":\"insulation\"}]}";
    auto insulationLayersJson = json::parse(insulationLayersString);

    std::map<std::pair<size_t, size_t>, std::vector<Layer>> insulationLayers;

    for (auto [key, layersJson] : insulationLayersJson.items()) {
        std::pair<size_t, size_t> windingsMapKey(key[0], key[1]);
        std::vector<Layer> layers;
        for (auto layerJson : layersJson) {
            layers.push_back(Layer(layerJson));
        }
        insulationLayers[windingsMapKey] = layers;
    }

    OpenMagnetics::Coil coil;

    if (insulationLayers.size() > 0) {
        coil.set_insulation_layers(insulationLayers);
    }
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_One_Layer", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    double wireDiameter = 0.000509;
    int64_t numberMaximumTurnsPerLayer = 9;
    int64_t numberMaximumLayers = 2;
    uint8_t interleavingLevel = 1;
    double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
    double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    auto numberReallyTestedWound = std::vector<int>(2, 0);
    for (size_t testIndex = 0; testIndex < 2; ++testIndex) {
        if (testIndex == 0) {
            settings.set_coil_try_rewind(false);
        }
        else {
            settings.set_coil_try_rewind(true);
        }

        for (size_t i = 0; i < 100; ++i)
        {
            std::vector<int64_t> numberTurns;
            std::vector<int64_t> numberParallels;
            int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
            for (size_t windingIndex = 0; windingIndex < OpenMagnetics::TestUtils::randomSize(1, 2 + 1 - 1); ++windingIndex)
            // for (size_t windingIndex = 0; windingIndex < OpenMagnetics::TestUtils::randomSize(1, 10 + 1 - 1); ++windingIndex)
            {
                int64_t numberPhysicalTurnsThisWinding = OpenMagnetics::TestUtils::randomSize(1, 300 + 1 - 1);
                int64_t numberTurnsThisWinding = OpenMagnetics::TestUtils::randomInt64(1, 100 + 1 - 1);
                int64_t numberParallelsThisWinding = std::max(1.0, std::ceil(double(numberPhysicalTurnsThisWinding) / numberTurnsThisWinding));
                numberTurns.push_back(numberTurnsThisWinding);
                numberParallels.push_back(numberParallelsThisWinding);
                numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
            }
            double bobbinHeight = 0.01;
            double bobbinWidth = 0.01;
            std::vector<double> bobbinCenterCoodinates = {0.05, 0, 0};
            uint8_t interleavingLevel = OpenMagnetics::TestUtils::randomInt(1, 10 + static_cast<int>(1) - 1);
            interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
            int windingOrientationIndex = OpenMagnetics::TestUtils::randomInt(0, 2 - 1);
            WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(windingOrientationIndex).value();

            // auto windingOrientation = OpenMagnetics::TestUtils::randomInt(0, 2 - 1)? WindingOrientation::CONTIGUOUS : WindingOrientation::OVERLAPPING;
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                bobbinWidth *= numberTurns.size();
                // bobbinCenterCoodinates[0] += bobbinWidth / 2;
            }
            else {
                bobbinHeight *= numberTurns.size();
            }

            int64_t numberPhysicalTurnsDebug = 0;
            for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                numberPhysicalTurnsDebug += numberTurns[windingIndex] * numberParallels[windingIndex];
            }

            try {
                auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

                if (coil.get_turns_description()) {
                    numberReallyTestedWound[testIndex]++;
                }

                bool result = OpenMagneticsTesting::check_turns_description(coil);

                if (!result) {

                    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                        std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
                    }
                    for (size_t windingIndex = 0; windingIndex < numberParallels.size(); ++windingIndex) {
                        std::cout << "numberParallels: " << numberParallels[windingIndex] << std::endl;
                    }
                    std::cout << "interleavingLevel: " << double(interleavingLevel) << std::endl;
                    std::cout << "windingOrientationIndex: " << windingOrientationIndex << std::endl;
                    return;

                }
            }
            catch (...) {
                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
                }
                for (size_t windingIndex = 0; windingIndex < numberParallels.size(); ++windingIndex) {
                    std::cout << "numberParallels: " << numberParallels[windingIndex] << std::endl;
                }
                std::cout << "interleavingLevel: " << double(interleavingLevel) << std::endl;
                std::cout << "windingOrientationIndex: " << windingOrientationIndex << std::endl;
                return;
            }
        }
    }

    REQUIRE(numberReallyTestedWound[1] > numberReallyTestedWound[0]);

    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding_0", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns;
    std::vector<int64_t> numberParallels;
    int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
    for (size_t windingIndex = 0; windingIndex < 1UL; ++windingIndex)
    {
        numberTurns.push_back(4);
        numberParallels.push_back(12);
        numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
    }
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 10;
    interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
    auto windingOrientation = WindingOrientation::CONTIGUOUS;
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        bobbinWidth *= numberTurns.size();
    }
    else {
        bobbinHeight *= numberTurns.size();
    }

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding_1", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {80};
    std::vector<int64_t> numberParallels = {3};
    int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();

    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
    {
        numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
    }
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 9;
    interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
    auto windingOrientation = WindingOrientation::OVERLAPPING;
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        bobbinWidth *= numberTurns.size();
    }
    else {
        bobbinHeight *= numberTurns.size();
    }

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

    OpenMagneticsTesting::check_layers_description(coil);
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding_2", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {39};
    std::vector<int64_t> numberParallels = {8};
    int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();

    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
    {
        numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
    }
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 7;
    interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        bobbinWidth *= numberTurns.size();
        bobbinCenterCoodinates[0] += bobbinWidth / 2;
    }
    else {
        bobbinHeight *= numberTurns.size();
    }

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

    OpenMagneticsTesting::check_turns_description(coil);
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Wind_By_Turn_Random_Multiwinding_2.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding_3", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {33, 18};
    std::vector<int64_t> numberParallels = {8, 2};
    int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();

    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
    {
        numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
    }
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 3;
    interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        bobbinWidth *= numberTurns.size();
        bobbinCenterCoodinates[0] += bobbinWidth / 2;
    }
    else {
        bobbinHeight *= numberTurns.size();
    }

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

    OpenMagneticsTesting::check_layers_description(coil);
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding_4", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {48, 68};
    std::vector<int64_t> numberParallels = {5, 2};
    int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();

    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
    {
        numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
    }
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 2;
    interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        bobbinWidth *= numberTurns.size();
        bobbinCenterCoodinates[0] += bobbinWidth / 2;
    }
    else {
        bobbinHeight *= numberTurns.size();
    }

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

    OpenMagneticsTesting::check_turns_description(coil);
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Wind_By_Turn_Random_Multiwinding_4.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding_5", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {16};
    std::vector<int64_t> numberParallels = {3};
    int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();

    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
    {
        numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
    }
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 4;
    interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        bobbinWidth *= numberTurns.size();
        bobbinCenterCoodinates[0] += bobbinWidth / 2;
    }
    else {
        bobbinHeight *= numberTurns.size();
    }

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

    OpenMagneticsTesting::check_turns_description(coil);
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Wind_By_Turn_Random_Multiwinding_4.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding_6", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {90, 37};
    std::vector<int64_t> numberParallels = {1, 1};
    int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();

    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
    {
        numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
    }
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 2;
    interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        bobbinWidth *= numberTurns.size();
        bobbinCenterCoodinates[0] += bobbinWidth / 2;
    }
    else {
        bobbinHeight *= numberTurns.size();
    }

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Random_Multiwinding_7", "[constructive-model][coil][rectangular-winding-window]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {1, 8};
    std::vector<int64_t> numberParallels = {7, 30};
    int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();

    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
    {
        numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
    }
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.01;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    uint8_t interleavingLevel = 1;
    interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        bobbinWidth *= numberTurns.size();
        bobbinCenterCoodinates[0] += bobbinWidth / 2;
    }
    else {
        bobbinHeight *= numberTurns.size();
    }

    auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

    OpenMagneticsTesting::check_turns_description(coil);
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Wind_By_Turn_Random_Multiwinding_4.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Rectangular_No_Bobbin", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.0038);
    wire.set_nominal_value_conducting_height(0.00076);
    wire.set_nominal_value_outer_width(0.004);
    wire.set_nominal_value_outer_height(0.0008);
    wire.set_type(WireType::RECTANGULAR);
    wires.push_back(wire);

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment,
                                                     wires,
                                                     false);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {40, 40};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 32/30",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 32/30", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    coil.set_interlayer_insulation(0.0001);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Two_Times", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Two_Times_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_interlayer_insulation(0.0002);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Two_Times_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_interlayer_insulation(0.0001);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Two_Times_After_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_interlayer_insulation(0);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Two_Times_After_After_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Toroidal_Core", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "T 17.3/9.7/12.7",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("T 17.3/9.7/12.7", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Toroidal_Core_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_interlayer_insulation(0.0005);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Toroidal_Core_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Toroidal_Core_Contiguous", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    // settings.set_coil_delimit_and_compact(false);
    // settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "T 17.3/9.7/12.7",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("T 17.3/9.7/12.7", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Toroidal_Core_Contiguous_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }

    coil.set_interlayer_insulation(0.0001);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Toroidal_Core_Contiguous_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterLayers_Only_Primary", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_Only_Primary_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_interlayer_insulation(0.0001, std::nullopt, "winding 0");

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_Only_Primary_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterLayers_Only_Secondary", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_Only_Secondary_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_interlayer_insulation(0.0001, std::nullopt, "winding 1");

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_Only_Secondary_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Contiguous_Layers", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    // settings.set_coil_delimit_and_compact(false);
    // settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Contiguous_Layers_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_interlayer_insulation(0.0001);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Contiguous_Layers_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_intersection_insulation(0.0002, 1);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Interleaved", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 2;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 40/40",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 40/40", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Interleaved_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_intersection_insulation(0.0001, 1);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Interleaved_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Interleaved_After_Sections.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        // painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Contiguous", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 28/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Contiguous_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_intersection_insulation(0.0002, 1);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Contiguous_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Layers_Toroidal_Core", "[constructive-model][coil][round-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 3};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "T 17.3/9.7/12.7",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("T 17.3/9.7/12.7", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Layers_Toroidal_Core_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_intersection_insulation(0.0001, 1);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Layers_Toroidal_Core_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Layers_Toroidal_Core_Contiguous", "[constructive-model][coil][round-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "T 17.3/9.7/12.7",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("T 17.3/9.7/12.7", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Layers_Toroidal_Core_Contiguous_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_intersection_insulation(0.0001, 1);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterSections_All_Layers_Toroidal_Core_Contiguous_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Change_Insulation_InterLayers_And_InterSections_All_Sections", "[constructive-model][coil][rectangular-winding-window][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {50, 50};
    std::vector<int64_t> numberParallels = {3, 2};
    uint8_t interleavingLevel = 2;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 40/40",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 40/40", json::parse("[]"), 1, "Dummy");
    // auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    // coil.set_bobbin(bobbin);
    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_And_InterSections_All_Sections_Before.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }

    coil.set_interlayer_insulation(0.00005);
    coil.set_intersection_insulation(0.0002, 1);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Change_Insulation_InterLayers_And_InterSections_All_Sections_After.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_One_Large_Layer_Toroidal", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {42};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, WindingOrientation::CONTIGUOUS);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Large_Layer_Toroidal.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(1U == coil.get_layers_description().value().size());
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_One_Full_Layer_Toroidal", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {58};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, WindingOrientation::CONTIGUOUS);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Full_Layer_Toroidal.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(1U == coil.get_layers_description().value().size());
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_Two_Layers_Toroidal", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {59};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, WindingOrientation::CONTIGUOUS);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_One_Section_Two_Layers_Toroidal.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {3};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset(); 
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(180, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[1].get_coordinates()[1], 0.001));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Top", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {3};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();

    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 0.5));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Bottom", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {3};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();

    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(357, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[2].get_coordinates()[1], 0.5));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Spread", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {3};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();

    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(60, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 0.5));
    REQUIRE_THAT(180, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[1].get_coordinates()[1], 0.5));
    REQUIRE_THAT(300, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[2].get_coordinates()[1], 0.5));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Contiguous_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {3, 3};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    settings.set_coil_try_rewind(false);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Contiguous_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(90, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[1].get_coordinates()[1], 0.5));
    REQUIRE_THAT(270, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[4].get_coordinates()[1], 0.5));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Overlapping_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {55, 55};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Overlapping_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Turn_Wind_Four_Sections_One_Layer_Toroidal_Overlapping_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_delimit_and_compact(false);
    std::vector<int64_t> numberTurns = {42, 42};
    std::vector<int64_t> numberParallels = {2, 2};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    // settings.set_coil_delimit_and_compact(false);
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turn_Wind_Four_Sections_One_Layer_Toroidal_Overlapping_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    // settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(coil.get_turns_description().value().size() == 135);
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(182, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[59].get_coordinates()[1], 1));
    REQUIRE_THAT(4.25, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(327, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[101].get_coordinates()[1], 1));
    REQUIRE_THAT(5.5, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    REQUIRE_THAT(299, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Bottom", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (true) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(coil.get_turns_description().value().size() == 135);
    REQUIRE_THAT(160, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(357, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[59].get_coordinates()[1], 1));
    REQUIRE_THAT(32, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(356, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[101].get_coordinates()[1], 1));
    REQUIRE_THAT(60, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    REQUIRE_THAT(355, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(81, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(173, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[15].get_coordinates()[1], 1));
    REQUIRE_THAT(180, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[16].get_coordinates()[1], 1));
    REQUIRE_THAT(272, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[31].get_coordinates()[1], 1));
    REQUIRE_THAT(327, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Spread", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    // settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(5, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(353, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[59].get_coordinates()[1], 1));
    REQUIRE_THAT(354, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Top", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[18].get_coordinates()[1], 1));
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[34].get_coordinates()[1], 1));
    REQUIRE_THAT(317, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[119].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Bottom", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(12, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(117, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[17].get_coordinates()[1], 1));
    REQUIRE_THAT(123, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(221, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    OpenMagneticsTesting::check_turns_description(coil);
    // Not clearly what this combination should do, so I check nothing
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Spread", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(117, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[17].get_coordinates()[1], 1));
    REQUIRE_THAT(123, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(243, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Top", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(42, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(147, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[17].get_coordinates()[1], 1));
    REQUIRE_THAT(43, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[34].get_coordinates()[1], 1));
    REQUIRE_THAT(357, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[119].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Bottom", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(42, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(147, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[17].get_coordinates()[1], 1));
    REQUIRE_THAT(44, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[34].get_coordinates()[1], 1));
    REQUIRE_THAT(357, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[119].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Spread", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(117, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[17].get_coordinates()[1], 1));
    REQUIRE_THAT(123, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(243, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Top", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(23, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(177, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[67].get_coordinates()[1], 1));
    REQUIRE_THAT(232, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    REQUIRE_THAT(329, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Bottom", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(23, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(177, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[67].get_coordinates()[1], 1));
    REQUIRE_THAT(232, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    REQUIRE_THAT(336, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(23, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(177, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[67].get_coordinates()[1], 1));
    REQUIRE_THAT(232, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    REQUIRE_THAT(333, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Spread", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(117, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[17].get_coordinates()[1], 1));
    REQUIRE_THAT(123, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(243, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(123, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(243, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Bottom", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(12, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(115, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[59].get_coordinates()[1], 1));
    REQUIRE_THAT(236, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[101].get_coordinates()[1], 1));
    REQUIRE_THAT(356, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Centered", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(7, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(109, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[59].get_coordinates()[1], 1));
    REQUIRE_THAT(223, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[101].get_coordinates()[1], 1));
    REQUIRE_THAT(348, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Spread", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(117, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[17].get_coordinates()[1], 1));
    REQUIRE_THAT(123, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(243, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Different_Wires", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 20, 20};
    std::vector<int64_t> numberParallels = {1, 5, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    // settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    std::vector<OpenMagnetics::Wire> wires;

    wires.push_back({find_wire_by_name("Round 0.335 - Grade 1")});
    wires.push_back({find_wire_by_name("Round 0.1 - Grade 2")});
    wires.push_back({find_wire_by_name("Litz 225x0.04 - Grade 1 - Double Served")});


    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Different_Wires.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(coil.get_turns_description().value().size() == 180);
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Different_Wires", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 20, 20};
    std::vector<int64_t> numberParallels = {1, 5, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    // settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    std::vector<OpenMagnetics::Wire> wires;

    wires.push_back({find_wire_by_name("Round 0.335 - Grade 1")});
    wires.push_back({find_wire_by_name("Round 0.1 - Grade 2")});
    wires.push_back({find_wire_by_name("Litz 225x0.04 - Grade 1 - Double Served")});


    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Different_Wires.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(coil.get_turns_description().value().size() == 180);
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Huge_Wire", "[constructive-model][coil][round-winding-window][smoke-test]") {
    std::vector<int64_t> numberTurns = {3};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    // settings.set_coil_try_rewind(false);
    // settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    std::vector<OpenMagnetics::Wire> wires;

    wires.push_back({find_wire_by_name("Litz 200x0.2 - Grade 2 - Double Served")});

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
    clear_databases();
    settings.set_use_toroidal_cores(true);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Huge_Wire.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(coil.get_turns_description().value().size() == 3);
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Rectangular_Wire", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {11, 90};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    std::vector<OpenMagnetics::Wire> wires;

    wires.push_back({find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
    wires.push_back({find_wire_by_name("Round 0.335 - Grade 1")});

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Rectangular_Wire.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(coil.get_turns_description().value().size() == 101);
    // Check this one manually, checking collision between two rotated rectangles is not worth it
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Rectangular_Wire", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {6, 90};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
    std::vector<OpenMagnetics::Wire> wires;

    wires.push_back({find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
    wires.push_back({find_wire_by_name("Round 0.335 - Grade 1")});

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Rectangular_Wire.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();

    REQUIRE(coil.get_turns_description().value().size() == 96);
    // Check this one manually, checking collision between two rotated rectangles is not worth it
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top_Margin", "[constructive-model][coil][round-winding-window][margin][smoke-test]") {
    settings.set_coil_equalize_margins(false);
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    settings.set_coil_delimit_and_compact(false);
    // settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    double margin = 0.0001;
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE(coil.get_turns_description().value().size() == 135);
    REQUIRE_THAT(3, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(186, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[59].get_coordinates()[1], 1));
    REQUIRE_THAT(4.25, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(175, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[101].get_coordinates()[1], 1));
    REQUIRE_THAT(7, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    REQUIRE_THAT(261, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Top_Margin", "[constructive-model][coil][round-winding-window][margin][smoke-test]") {
    settings.set_coil_equalize_margins(false);
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    double margin = 0.0002;
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Top_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(6, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(161, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(258, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[102].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Top_Margin", "[constructive-model][coil][round-winding-window][margin][smoke-test]") {
    settings.set_coil_equalize_margins(false);
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    double margin = 0.0002;
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Top_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(33, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(188, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(332, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Top_Margin", "[constructive-model][coil][round-winding-window][margin][smoke-test]") {
    settings.set_coil_equalize_margins(false);
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    double margin = 0.0002;
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Top_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(20, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(174, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(318, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Margin", "[constructive-model][coil][round-winding-window][margin][smoke-test]") {
    settings.set_coil_equalize_margins(false);
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    double margin = 0.0002;
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(7, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(131, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(341, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Spread_Margin", "[constructive-model][coil][round-winding-window][margin][smoke-test]") {
    settings.set_coil_equalize_margins(false);
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    double margin = 0.0002;
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
    coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
    coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Spread_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    REQUIRE_THAT(7, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[0].get_coordinates()[1], 1));
    REQUIRE_THAT(131, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[60].get_coordinates()[1], 1));
    REQUIRE_THAT(349, Catch::Matchers::WithinAbs(coil.get_turns_description().value()[134].get_coordinates()[1], 1));
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Additional_Coordinates", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    auto turns = coil.get_turns_description().value();

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Additional_Coordinates.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    for (auto turn : turns) {
        REQUIRE(turn.get_additional_coordinates());
    }
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Spread_Top_Additional_Coordinates", "[constructive-model][coil][round-winding-window][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {60, 42, 33};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);


    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Spread_Top_Additional_Coordinates.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
    }
    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    auto turns = coil.get_turns_description().value();
    for (auto turn : turns) {
        REQUIRE(turn.get_additional_coordinates());
        if (turn.get_additional_coordinates()) {
            auto additionalCoordinates = turn.get_additional_coordinates().value();

            for (auto additionalCoordinate : additionalCoordinates){
                REQUIRE(additionalCoordinate[0] < 0);
            }
        }
    }
    OpenMagneticsTesting::check_turns_description(coil);
}

TEST_CASE("Test_Wind_By_Layers_Planar_One_Layer", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};
    std::vector<size_t> stackUp = {0};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.02;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(bobbinHeight, bobbinWidth);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001);
    wire.set_nominal_value_conducting_height(0.00076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp, {});
    coil.wind_by_planar_layers();
    auto layersDescription = coil.get_layers_description().value();
    REQUIRE(layersDescription.size() == 1);
}

TEST_CASE("Test_Wind_By_Layers_Planar_Two_Layers", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};
    std::vector<size_t> stackUp = {0, 0};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.02;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(bobbinHeight, bobbinWidth);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001);
    wire.set_nominal_value_conducting_height(0.00076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp, {});
    coil.wind_by_planar_layers();
    auto layersDescription = coil.get_layers_description().value();
    REQUIRE(layersDescription.size() == 3);
}

TEST_CASE("Test_Wind_By_Layers_Planar_Two_Windings", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7, 7};
    std::vector<int64_t> numberParallels = {1, 1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
    std::vector<size_t> stackUp = {0, 1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.02;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(bobbinHeight, bobbinWidth);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001);
    wire.set_nominal_value_conducting_height(0.00076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp, {});
    coil.wind_by_planar_layers();
    auto layersDescription = coil.get_layers_description().value();
    REQUIRE(3U == layersDescription.size());
    REQUIRE(1U == layersDescription[0].get_partial_windings().size());
    REQUIRE("PRIMARY" == layersDescription[0].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(1 == layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[2].get_partial_windings().size());
    REQUIRE("SECONDARY" == layersDescription[2].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(1 == layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);
}

TEST_CASE("Test_Wind_By_Layers_Planar_Two_Windings_Two_Layers_No_Interleaved", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {8, 8};
    std::vector<int64_t> numberParallels = {1, 1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
    std::vector<size_t> stackUp = {0, 0, 1, 1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.02;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(bobbinHeight, bobbinWidth);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001);
    wire.set_nominal_value_conducting_height(0.00076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp, {});
    coil.wind_by_planar_layers();
    auto layersDescription = coil.get_layers_description().value();
    REQUIRE(7U == layersDescription.size());
    REQUIRE(1U == layersDescription[0].get_partial_windings().size());
    REQUIRE("PRIMARY" == layersDescription[0].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(0.5 == layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[2].get_partial_windings().size());
    REQUIRE("PRIMARY" == layersDescription[2].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(0.5 == layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[4].get_partial_windings().size());
    REQUIRE("SECONDARY" == layersDescription[4].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[4].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(0.5 == layersDescription[4].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[6].get_partial_windings().size());
    REQUIRE("SECONDARY" == layersDescription[6].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[6].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(0.5 == layersDescription[6].get_partial_windings()[0].get_parallels_proportion()[0]);
}

TEST_CASE("Test_Wind_By_Layers_Planar_Two_Windings_Two_Layers_No_Interleaved_Odd_Turns", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {3, 3};
    std::vector<int64_t> numberParallels = {1, 1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
    std::vector<size_t> stackUp = {0, 0, 1, 1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.02;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(bobbinHeight, bobbinWidth);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001);
    wire.set_nominal_value_conducting_height(0.00076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp, {});
    coil.wind_by_planar_layers();
    auto layersDescription = coil.get_layers_description().value();
    REQUIRE(7U == layersDescription.size());
    REQUIRE(1U == layersDescription[0].get_partial_windings().size());
    REQUIRE("PRIMARY" == layersDescription[0].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(2.0 / 3 == layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[2].get_partial_windings().size());
    REQUIRE("PRIMARY" == layersDescription[2].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(1.0 / 3 == layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[4].get_partial_windings().size());
    REQUIRE("SECONDARY" == layersDescription[4].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[4].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(2.0 / 3 == layersDescription[4].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[6].get_partial_windings().size());
    REQUIRE("SECONDARY" == layersDescription[6].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[6].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(1.0 / 3 == layersDescription[6].get_partial_windings()[0].get_parallels_proportion()[0]);
}

TEST_CASE("Test_Wind_By_Layers_Planar_Two_Windings_Two_Layers_Interleaved_Odd_Turns", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {3, 3};
    std::vector<int64_t> numberParallels = {1, 1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
    std::vector<size_t> stackUp = {0, 1, 0, 1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.02;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(bobbinHeight, bobbinWidth);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001);
    wire.set_nominal_value_conducting_height(0.00076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp, {});
    coil.wind_by_planar_layers();
    auto layersDescription = coil.get_layers_description().value();
    REQUIRE(7U == layersDescription.size());
    REQUIRE(1U == layersDescription[0].get_partial_windings().size());
    REQUIRE("PRIMARY" == layersDescription[0].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(2.0 / 3 == layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[2].get_partial_windings().size());
    REQUIRE("SECONDARY" == layersDescription[2].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(2.0 / 3 == layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[4].get_partial_windings().size());
    REQUIRE("PRIMARY" == layersDescription[4].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[4].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(1.0 / 3 == layersDescription[4].get_partial_windings()[0].get_parallels_proportion()[0]);
    REQUIRE(1U == layersDescription[6].get_partial_windings().size());
    REQUIRE("SECONDARY" == layersDescription[6].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[6].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(1.0 / 3 == layersDescription[6].get_partial_windings()[0].get_parallels_proportion()[0]);
}

TEST_CASE("Test_Wind_By_Layers_Planar_Two_Windings_Two_Layers_Interleaved_Odd_Turns_With_Insulation", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {3, 3};
    std::vector<int64_t> numberParallels = {1, 1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
    std::vector<size_t> stackUp = {0, 1, 0, 1};
    double bobbinHeight = 0.01;
    double bobbinWidth = 0.02;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(bobbinHeight, bobbinWidth);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001);
    wire.set_nominal_value_conducting_height(0.00076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp);
    coil.wind_by_planar_layers();
    auto layersDescription = coil.get_layers_description().value();
    REQUIRE(7U == layersDescription.size());
    REQUIRE(MAS::ElectricalType::CONDUCTION == layersDescription[0].get_type());
    REQUIRE(1U == layersDescription[0].get_partial_windings().size());
    REQUIRE("PRIMARY" == layersDescription[0].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(2.0 / 3 == layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);

    REQUIRE(MAS::ElectricalType::INSULATION == layersDescription[1].get_type());

    REQUIRE(1U == layersDescription[2].get_partial_windings().size());
    REQUIRE(MAS::ElectricalType::CONDUCTION == layersDescription[2].get_type());
    REQUIRE("SECONDARY" == layersDescription[2].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(2.0 / 3 == layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);

    REQUIRE(MAS::ElectricalType::INSULATION == layersDescription[3].get_type());

    REQUIRE(1U == layersDescription[4].get_partial_windings().size());
    REQUIRE(MAS::ElectricalType::CONDUCTION == layersDescription[4].get_type());
    REQUIRE("PRIMARY" == layersDescription[4].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[4].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(1.0 / 3 == layersDescription[4].get_partial_windings()[0].get_parallels_proportion()[0]);

    REQUIRE(MAS::ElectricalType::INSULATION == layersDescription[5].get_type());

    REQUIRE(1U == layersDescription[6].get_partial_windings().size());
    REQUIRE(MAS::ElectricalType::CONDUCTION == layersDescription[6].get_type());
    REQUIRE("SECONDARY" == layersDescription[6].get_partial_windings()[0].get_winding());
    REQUIRE(1U == layersDescription[6].get_partial_windings()[0].get_parallels_proportion().size());
    REQUIRE(1.0 / 3 == layersDescription[6].get_partial_windings()[0].get_parallels_proportion()[0]);
}

TEST_CASE("Test_Wind_By_Turns_Planar_One_Layer", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};
    std::vector<size_t> stackUp = {0};
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto core = OpenMagneticsTesting::get_quick_core("ELP 32/6/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001);
    wire.set_nominal_value_conducting_height(0.000076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp, {});
    coil.wind_by_planar_layers();
    coil.wind_by_planar_turns(0.0002, {{0, 0.0002}});
    coil.delimit_and_compact();
    REQUIRE(coil.get_turns_description());
    auto turnsDescription = coil.get_turns_description().value();
    REQUIRE(turnsDescription.size() == 7);
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turns_Planar_One_Layer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Test_Wind_By_Turns_Planar_Two_Windings_Two_Layers_Interleaved_Odd_Turns_With_Insulation", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {20, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
    std::vector<size_t> stackUp = {0, 1, 0, 1};
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto core = OpenMagneticsTesting::get_quick_core("ELP 38/8/25", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);

    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.0008);
    wire.set_nominal_value_conducting_height(0.000076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);
    wires.push_back(wire);
    wire.set_nominal_value_conducting_width(0.0032);
    wire.set_nominal_value_conducting_height(0.000076);
    wires.push_back(wire);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wires[windingIndex]);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);
    coil.set_strict(false);

    coil.wind_by_planar_sections(stackUp, {{{0, 1}, 0.0005}}, 0.0005);
    coil.wind_by_planar_layers();
    coil.wind_by_planar_turns(0.0002, {{0, 0.0002}, {1, 0.0002}});
    coil.delimit_and_compact();
    REQUIRE(coil.get_turns_description());
    if (coil.get_turns_description()) {
        auto turnsDescription = coil.get_turns_description().value();
        REQUIRE(turnsDescription.size() == 25);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turns_Planar_Two_Windings_Two_Layers_Interleaved_Odd_Turns_With_Insulation.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
}

TEST_CASE("Test_Wind_By_Turns_Planar_Many_Layers", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {20, 5};
    std::vector<int64_t> numberParallels = {4, 4};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
    std::vector<size_t> stackUp = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto core = OpenMagneticsTesting::get_quick_core("ELP 38/8/25", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);

    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.0008);
    wire.set_nominal_value_conducting_height(0.000076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);
    wires.push_back(wire);
    wire.set_nominal_value_conducting_width(0.0032);
    wire.set_nominal_value_conducting_height(0.000076);
    wires.push_back(wire);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wires[windingIndex]);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);
    coil.set_strict(false);

    coil.wind_by_planar_sections(stackUp, {{{0, 1}, 0.0001}}, 0.0001);
    coil.wind_by_planar_layers();
    coil.wind_by_planar_turns(0.0002, {{0, 0.0002}, {1, 0.0002}});
    coil.delimit_and_compact();
    REQUIRE(coil.get_turns_description());
    if (coil.get_turns_description()) {
        auto turnsDescription = coil.get_turns_description().value();
        REQUIRE(turnsDescription.size() == 100);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turns_Planar_Many_Layers.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            // painter.paint_coil_layers(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
}

TEST_CASE("Test_Wind_By_Turns_Planar_One_Layer_Distance_To_Core", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {7};
    std::vector<int64_t> numberParallels = {1};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};
    std::vector<size_t> stackUp = {0};
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto core = OpenMagneticsTesting::get_quick_core("ELP 32/6/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.0005);
    wire.set_nominal_value_conducting_height(0.000076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wire);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);

    coil.wind_by_planar_sections(stackUp, {{{0, 1}, 0.0001}}, 0.001);
    coil.wind_by_planar_layers();
    coil.wind_by_planar_turns(0, {{0, 0.0002}, {1, 0.0002}});
    coil.delimit_and_compact();
    REQUIRE(coil.get_turns_description());
    auto turnsDescription = coil.get_turns_description().value();
    REQUIRE(turnsDescription.size() == 7);
    if (plot) {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turns_Planar_One_Layer_Distance_To_Core.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
}

TEST_CASE("Test_Wind_By_Turns_Planar_Many_Layers_Magnetic_Field", "[constructive-model][coil][planar][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);

    std::vector<int64_t> numberTurns = {20, 5};
    std::vector<int64_t> numberParallels = {4, 4};
    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
    std::vector<size_t> stackUp = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    auto core = OpenMagneticsTesting::get_quick_core("ELP 38/8/25", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);

    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.0008);
    wire.set_nominal_value_conducting_height(0.000076);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);
    wires.push_back(wire);
    wire.set_nominal_value_conducting_width(0.0032);
    wire.set_nominal_value_conducting_height(0.000076);
    wires.push_back(wire);

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
        OpenMagnetics::Winding coilFunctionalDescription; 
        coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
        coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
        coilFunctionalDescription.set_name(std::string{magic_enum::enum_name(isolationSides[windingIndex])});
        coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
        coilFunctionalDescription.set_wire(wires[windingIndex]);
        coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
    }
    coil.set_bobbin(bobbin);
    coil.set_strict(false);

    coil.wind_by_planar_sections(stackUp, {{{0, 1}, 0.0001}}, 0.0001);
    coil.wind_by_planar_layers();
    coil.wind_by_planar_turns(0.0002, {{0, 0.0002}, {1, 0.0002}});
    coil.delimit_and_compact();
    REQUIRE(coil.get_turns_description());
    if (coil.get_turns_description()) {
        auto turnsDescription = coil.get_turns_description().value();
        REQUIRE(turnsDescription.size() == 100);
        if (plot) {
            double voltagePeakToPeak = 2000;
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, {double(numberTurns[0]) / numberTurns[1]});
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turns_Planar_Many_Layers_Magnetic_Field.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            // painter.paint_coil_layers(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
}

TEST_CASE("Test_Get_Round_Wire_From_Dc_Resistance", "[constructive-model][coil][smoke-test]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    std::vector<int64_t> numberTurns = {1, 60};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    std::string coreShape = "EE5";
    std::string coreMaterial = "3C97"; 
    auto emptyGapping = json::array();
    // settings.set_coil_delimit_and_compact(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    std::vector<double> dcResistances = {0.00075, 1.75};
    auto wires = coil.guess_round_wire_from_dc_resistance(dcResistances, 0.01);
    REQUIRE(wires[0].get_name().value() == "Round 0.63 - Grade 1");
    REQUIRE(wires[1].get_name().value() == "Round 0.106 - Grade 1");
}

TEST_CASE("Test_Wind_By_Sections_Two_Windings_Together", "[constructive-model][coil][groups][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {5, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 40/40",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    coil.get_mutable_functional_description()[0].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[1].get_name()});
    coil.get_mutable_functional_description()[0].set_isolation_side(MAS::IsolationSide::PRIMARY);
    coil.get_mutable_functional_description()[1].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[0].get_name()});
    coil.get_mutable_functional_description()[1].set_isolation_side(MAS::IsolationSide::PRIMARY);
    auto core = OpenMagneticsTesting::get_quick_core("PQ 40/40", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    coil.set_bobbin(bobbin);
    coil.wind_by_sections();
    REQUIRE(1 == coil.get_sections_description()->size());
    REQUIRE(2 == coil.get_sections_description().value()[0].get_partial_windings().size());
    REQUIRE("winding 0" == coil.get_sections_description().value()[0].get_partial_windings()[0].get_winding());
    REQUIRE("winding 1" == coil.get_sections_description().value()[0].get_partial_windings()[1].get_winding());
    auto virtualFunctionalDescription = coil.virtualize_functional_description();
    REQUIRE(1 == virtualFunctionalDescription.size());
    REQUIRE(numberTurns[0] + numberTurns[1] == virtualFunctionalDescription[0].get_number_turns());
    REQUIRE(numberParallels[0] == virtualFunctionalDescription[0].get_number_parallels());

    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Sections_Two_Windings_Together.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Sections_Two_Windings_Together_One_Not", "[constructive-model][coil][groups][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {5, 5, 12};
    std::vector<int64_t> numberParallels = {2, 2, 3};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 40/40",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    coil.get_mutable_functional_description()[0].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[1].get_name()});
    coil.get_mutable_functional_description()[0].set_isolation_side(MAS::IsolationSide::PRIMARY);
    coil.get_mutable_functional_description()[1].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[0].get_name()});
    coil.get_mutable_functional_description()[1].set_isolation_side(MAS::IsolationSide::PRIMARY);
    auto core = OpenMagneticsTesting::get_quick_core("PQ 40/40", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    coil.set_bobbin(bobbin);
    coil.wind_by_sections();

    REQUIRE(4 == coil.get_sections_description()->size());
    REQUIRE(2 == coil.get_sections_description().value()[0].get_partial_windings().size());
    REQUIRE(1 == coil.get_sections_description().value()[2].get_partial_windings().size());
    REQUIRE("winding 0" == coil.get_sections_description().value()[0].get_partial_windings()[0].get_winding());
    REQUIRE("winding 1" == coil.get_sections_description().value()[0].get_partial_windings()[1].get_winding());
    REQUIRE("winding 2" == coil.get_sections_description().value()[2].get_partial_windings()[0].get_winding());
    auto virtualFunctionalDescription = coil.virtualize_functional_description();
    REQUIRE(2 == virtualFunctionalDescription.size());
    REQUIRE(numberTurns[0] + numberTurns[1] == virtualFunctionalDescription[0].get_number_turns());
    REQUIRE(numberParallels[0] == virtualFunctionalDescription[0].get_number_parallels());
    REQUIRE(numberTurns[2] == virtualFunctionalDescription[1].get_number_turns());
    REQUIRE(numberParallels[2] == virtualFunctionalDescription[1].get_number_parallels());

    OpenMagneticsTesting::check_turns_description(coil);

    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Sections_Two_Windings_Together_One_Not.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Layers_Two_Windings_Together", "[constructive-model][coil][groups][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {5, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 40/40",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    coil.get_mutable_functional_description()[0].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[1].get_name()});
    coil.get_mutable_functional_description()[0].set_isolation_side(MAS::IsolationSide::PRIMARY);
    coil.get_mutable_functional_description()[1].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[0].get_name()});
    coil.get_mutable_functional_description()[1].set_isolation_side(MAS::IsolationSide::PRIMARY);
    auto core = OpenMagneticsTesting::get_quick_core("PQ 40/40", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    coil.set_bobbin(bobbin);
    coil.wind_by_sections();
    coil.wind_by_layers();
    REQUIRE(1 == coil.get_layers_description()->size());
    REQUIRE(2 == coil.get_layers_description().value()[0].get_partial_windings().size());
    REQUIRE("winding 0" == coil.get_layers_description().value()[0].get_partial_windings()[0].get_winding());
    REQUIRE("winding 1" == coil.get_layers_description().value()[0].get_partial_windings()[1].get_winding());

    OpenMagneticsTesting::check_turns_description(coil);
    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Layers_Two_Windings_Together.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        // painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Layers_Two_Windings_Together_One_Not", "[constructive-model][coil][groups][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {5, 5, 12};
    std::vector<int64_t> numberParallels = {2, 2, 3};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 40/40",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    coil.get_mutable_functional_description()[0].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[1].get_name()});
    coil.get_mutable_functional_description()[0].set_isolation_side(MAS::IsolationSide::PRIMARY);
    coil.get_mutable_functional_description()[1].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[0].get_name()});
    coil.get_mutable_functional_description()[1].set_isolation_side(MAS::IsolationSide::PRIMARY);
    auto core = OpenMagneticsTesting::get_quick_core("PQ 40/40", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    coil.set_bobbin(bobbin);
    coil.wind_by_sections();
    coil.wind_by_layers();
    REQUIRE(4 == coil.get_layers_description()->size());
    REQUIRE(2 == coil.get_layers_description().value()[0].get_partial_windings().size());
    REQUIRE(1 == coil.get_layers_description().value()[2].get_partial_windings().size());
    REQUIRE("winding 0" == coil.get_layers_description().value()[0].get_partial_windings()[0].get_winding());
    REQUIRE("winding 1" == coil.get_layers_description().value()[0].get_partial_windings()[1].get_winding());
    REQUIRE("winding 2" == coil.get_layers_description().value()[2].get_partial_windings()[0].get_winding());

    OpenMagneticsTesting::check_turns_description(coil);
    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Layers_Two_Windings_Together_One_Not.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        // painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turns_Two_Windings_Together", "[constructive-model][coil][groups][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {5, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 40/40",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    coil.get_mutable_functional_description()[0].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[1].get_name()});
    coil.get_mutable_functional_description()[0].set_isolation_side(MAS::IsolationSide::PRIMARY);
    coil.get_mutable_functional_description()[1].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[0].get_name()});
    coil.get_mutable_functional_description()[1].set_isolation_side(MAS::IsolationSide::PRIMARY);
    auto core = OpenMagneticsTesting::get_quick_core("PQ 40/40", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    coil.set_bobbin(bobbin);
    coil.wind_by_sections();
    coil.wind_by_layers();
    coil.wind_by_turns();
    REQUIRE(1 == coil.get_layers_description()->size());
    REQUIRE(2 == coil.get_layers_description().value()[0].get_partial_windings().size());
    REQUIRE("winding 0" == coil.get_layers_description().value()[0].get_partial_windings()[0].get_winding());
    REQUIRE("winding 1" == coil.get_layers_description().value()[0].get_partial_windings()[1].get_winding());

    OpenMagneticsTesting::check_turns_description(coil);
    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turns_Two_Windings_Together.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

TEST_CASE("Test_Wind_By_Turns_Two_Windings_Together_One_Not", "[constructive-model][coil][groups][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(true);
    std::vector<int64_t> numberTurns = {5, 5, 12};
    std::vector<int64_t> numberParallels = {2, 2, 3};
    uint8_t interleavingLevel = 1;
    std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;


    WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     "PQ 40/40",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    coil.get_mutable_functional_description()[0].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[1].get_name()});
    coil.get_mutable_functional_description()[0].set_isolation_side(MAS::IsolationSide::PRIMARY);
    coil.get_mutable_functional_description()[1].set_wound_with(std::vector<std::string>{coil.get_mutable_functional_description()[0].get_name()});
    coil.get_mutable_functional_description()[1].set_isolation_side(MAS::IsolationSide::PRIMARY);
    auto core = OpenMagneticsTesting::get_quick_core("PQ 40/40", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    coil.set_bobbin(bobbin);
    coil.wind();
    REQUIRE(4 == coil.get_layers_description()->size());
    REQUIRE(2 == coil.get_layers_description().value()[0].get_partial_windings().size());
    REQUIRE("winding 0" == coil.get_layers_description().value()[0].get_partial_windings()[0].get_winding());
    REQUIRE("winding 1" == coil.get_layers_description().value()[0].get_partial_windings()[1].get_winding());
    REQUIRE("winding 2" == coil.get_layers_description().value()[2].get_partial_windings()[0].get_winding());

    OpenMagneticsTesting::check_turns_description(coil);
    {
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Wind_By_Turns_Two_Windings_Together_One_Not.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_coil(coil);
        magnetic.set_core(core);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        settings.reset();
    }
}

}  // namespace
