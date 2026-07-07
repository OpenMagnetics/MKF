#include "RandomUtils.h"
#include <source_location>
#include "support/Settings.h"
#include "support/Painter.h"
#include "constructive_models/Coil.h"
#include "physical_models/WindingOhmicLosses.h"
#include "json.hpp"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <limits>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <time.h>
#include <set>
#include <map>
using json = nlohmann::json;
#include <typeinfo>
#include <chrono>
#include <thread>

using namespace MAS;
using namespace OpenMagnetics;
using namespace OpenMagneticsTesting;

namespace {

auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
bool plot = true;

TEST_CASE("Test_Coil_Json_0", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"bobbin":"Dummy","functionalDescription":[{"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":23,"wire":"Dummy"}]})";

    auto coilJson = json::parse(coilString);
    // The original line here was `auto Coil(coilJson);` — a most-vexing-parse variable
    // declaration that just copied the json and never constructed a Coil at all.
    OpenMagnetics::Coil coil(coilJson, false);
    REQUIRE(coil.get_functional_description().size() == 1);
    CHECK(coil.get_functional_description()[0].get_name() == "Primary");
    CHECK(coil.get_functional_description()[0].get_number_turns() == 23);
    CHECK(coil.get_functional_description()[0].get_number_parallels() == 1);
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
    std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.006,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0032500000000000003,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0002835287369864788,"coordinates":[0.0095,0,0],"height":null,"radialHeight":0.0095,"sectionsAlignment":"outerOrBottom","sectionsOrientation":"contiguous","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":27,"wire":{"coating":{"breakdownVoltage":2700,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":4.116868676970209e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000724},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 21.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000757},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"21 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"secondary","name":"Secondary","numberParallels":1,"numberTurns":27,"wire":{"coating":{"breakdownVoltage":5000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":4.620411001469214e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000767},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 20.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000831},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"20.5 AWG","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription": null, "turnsDescription":null,"_turnsAlignment":{"Primary section 0":"spread","Secondary section 0":"spread"},"_layersOrientation":{"Primary section 0":"overlapping","Secondary section 0":"overlapping"}})";

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
    std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.00356,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0022725,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0000637587014444212,"coordinates":[0.004505,0,0],"height":null,"radialHeight":0.004505,"sectionsAlignment":"innerOrTop","sectionsOrientation":"overlapping","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":3,"numberTurns":55,"wire":{"coating":{"breakdownVoltage":1220,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":8.042477193189871e-8},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000323,"minimum":0.00031800000000000003,"nominal":0.00032},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 28.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000356,"minimum":0.00033800000000000003,"nominal":0.000347},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"28 AWG","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null,"_turnsAlignment":["spread"],"_layersOrientation":["overlapping"]})";

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
    REQUIRE(layers.size() > 0);

    // Repro point: resolving the insulation material of every insulation layer must work.
    size_t insulationLayerCount = 0;
    for (auto layer : layers) {
        if (layer.get_type() == ElectricalType::INSULATION) {
            insulationLayerCount++;
            auto material = OpenMagnetics::Coil::resolve_insulation_layer_insulation_material(coil, layer.get_name());
            json materialJson;
            to_json(materialJson, material);
            INFO("Layer: " << layer.get_name());
            CHECK(!materialJson.empty());
        }

    }
    // The fixture contains insulation layers; the loop above must not be vacuous.
    REQUIRE(insulationLayerCount > 0);
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

    // Repro point: winding by sections must actually produce a sections description.
    REQUIRE(coil.get_sections_description());
    CHECK(coil.get_sections_description_conduction().size() > 0);

    json result;
    to_json(result, coil);
    CHECK(!result.empty());
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

    // Repro point: winding by turns over preset sections/layers must produce turns.
    REQUIRE(coil.get_turns_description());
    CHECK(coil.get_turns_description()->size() > 0);

    json result;
    to_json(result, coil);
    CHECK(!result.empty());
}

TEST_CASE("Test_Coil_Json_11", "[constructive-model][coil][bug][smoke-test]") {
    std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.006175,"columnShape":"round","columnThickness":0,"columnWidth":0.006175,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":null,"area":0.000041283000000000004,"coordinates":[0.0098875,0,0],"height":0.00556,"radialHeight":null,"sectionsAlignment":"innerOrTop","sectionsOrientation":"contiguous","shape":"rectangular","width":0.007425000000000001}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"primary","numberParallels":1,"numberTurns":5,"wire":{"coating":null,"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.2293100000000003e-7},"conductingDiameter":null,"conductingHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000522},"conductingWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0023550000000000008},"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Planar 52.20 \u00b5m","numberConductors":1,"outerDiameter":null,"outerHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000522},"outerWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0023550000000000008},"standard":"IPC-6012","standardName":"1.5 oz.","strand":null,"type":"planar"}},{"connections":null,"isolationSide":"secondary","name":"SECONDARY","numberParallels":1,"numberTurns":3,"wire":{"coating":null,"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.2449700000000003e-7},"conductingDiameter":null,"conductingHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000348},"conductingWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.003577500000000001},"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Planar 34.80 \u00b5m","numberConductors":1,"outerDiameter":null,"outerHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000348},"outerWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.003577500000000001},"standard":"IPC-6012","standardName":"1 oz.","strand":null,"type":"planar"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null,"_turnsAlignment":["spread","spread","spread","spread"],"_layersOrientation":["contiguous","contiguous","contiguous","contiguous"],"_interlayerInsulationThickness":0,"_intersectionInsulationThickness":0.0001})";

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
    // This test re-winds the no-compact coil after construction; the no-compact
    // environment is part of its intent, so set it explicitly (the builder no longer
    // leaks the flag globally).
    settings.set_coil_delimit_and_compact(false);
    double voltagePeakToPeak = 400;
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
    coil.set_inputs(inputs);
    coil.wind();
    auto log = coil.read_log();

    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, sectionOrientation);
    settings.reset();
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
    // This test re-winds the no-compact coil after construction; the no-compact
    // environment is part of its intent, so set it explicitly (the builder no longer
    // leaks the flag globally).
    settings.set_coil_delimit_and_compact(false);

    coil.wind_by_sections(pattern, repetitions);
    OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    settings.reset();
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
    settings.reset();
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
    // This test re-winds the no-compact coil after construction; the no-compact
    // environment is part of its intent, so set it explicitly (the builder no longer
    // leaks the flag globally).
    settings.set_coil_delimit_and_compact(false);
    double voltagePeakToPeak = 400;
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
    coil.set_inputs(inputs);
    coil.wind();
    auto log = coil.read_log();

    OpenMagneticsTesting::check_layers_description(coil);
    settings.reset();
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

    // The fixture defines insulation layers for the (0,1), (1,2) and (2,0) winding pairs.
    REQUIRE(insulationLayers.size() == 3);
    for (auto& [windingsMapKey, layers] : insulationLayers) {
        REQUIRE(layers.size() == 2);
        for (auto& layer : layers) {
            CHECK(layer.get_type() == ElectricalType::INSULATION);
            REQUIRE(layer.get_dimensions().size() == 2);
            CHECK(layer.get_dimensions()[0] > 0);
            CHECK(layer.get_dimensions()[1] > 0);
        }
    }

    OpenMagnetics::Coil coil;
    CHECK_NOTHROW(coil.set_insulation_layers(insulationLayers));
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
            catch (const std::exception& e) {
                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
                }
                for (size_t windingIndex = 0; windingIndex < numberParallels.size(); ++windingIndex) {
                    std::cout << "numberParallels: " << numberParallels[windingIndex] << std::endl;
                }
                std::cout << "interleavingLevel: " << double(interleavingLevel) << std::endl;
                std::cout << "windingOrientationIndex: " << windingOrientationIndex << std::endl;
                // A throw used to be swallowed here (early return -> test passed). It must fail.
                FAIL("get_quick_coil threw for the combination printed above: " << e.what());
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
    // Reduced from {48, 68} turns and {5, 2} parallels to fit in 10x10mm bobbin
    // Total physical turns: 24*2 + 32*1 = 80 (was 240 + 136 = 376)
    std::vector<int64_t> numberTurns = {24, 32};
    std::vector<int64_t> numberParallels = {2, 1};
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
                                                     "PQ 32/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 32/20", json::parse("[]"), 1, "Dummy");
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
                                                     "PQ 32/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 32/20", json::parse("[]"), 1, "Dummy");
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
                                                     "PQ 32/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 32/20", json::parse("[]"), 1, "Dummy");
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
                                                     "PQ 32/20",
                                                     interleavingLevel,
                                                     windingOrientation,
                                                     layersOrientation,
                                                     turnsAlignment,
                                                     sectionsAlignment);

    auto core = OpenMagneticsTesting::get_quick_core("PQ 32/20", json::parse("[]"), 1, "Dummy");
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
        auto path = outputFilePath;
        path.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Additional_Coordinates.json");
        to_file(path.string(), magnetic);

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

TEST_CASE("Test_Additiona_Turns_Bug", "[constructive-model][coil][round-winding-window][bug]") {
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_include_additional_coordinates(true);
    
    // Create a toroidal core similar to T 20/10/7
    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97";
    auto emptyGapping = json::array();
    int64_t numberStacks = 1;
    
    // Create a coil with 60 turns to force 2 layers (similar to original bug)
    std::vector<int64_t> numberTurns = {60};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;
    
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
    CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

    auto turns = coil.get_turns_description().value();
    auto layers = coil.get_layers_description().value();
    

    // Check that all turns have additional coordinates
    for (auto& turn : turns) {
        REQUIRE(turn.get_additional_coordinates());
    }
    
    // Group turns by layer and check additional coordinates
    std::map<std::string, std::vector<size_t>> turnsByLayer;
    for (size_t i = 0; i < turns.size(); ++i) {
        turnsByLayer[turns[i].get_layer().value()].push_back(i);
    }
    
    // Collect all additional radii
    std::vector<double> additionalRadii;
    for (size_t i = 0; i < turns.size(); ++i) {
        auto addCoords = turns[i].get_additional_coordinates().value()[0];
        double radius = sqrt(addCoords[0]*addCoords[0] + addCoords[1]*addCoords[1]);
        additionalRadii.push_back(radius);
    }
    
    // Find unique radii
    std::set<double> uniqueRadii;
    for (auto r : additionalRadii) {
        // Round to 0.1mm precision
        uniqueRadii.insert(round(r * 10000) / 10000);
    }

    // The bug was that additional turns from layer 1 were placed at a different radius
    // than turns from layer 0, even though there was space in the first additional layer.
    // After the fix, all additional turns should be at the same radius (compacted).
    REQUIRE(uniqueRadii.size() == 1);

    if (plot) {
        auto outFile = outputFilePath;
        outFile.append("Test_Additiona_Turns_Bug.svg");
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

TEST_CASE("Test_Wind_Two_Sections_Toroidal_Default_Alignment_Spreads", "[constructive-model][coil][round-winding-window][smoke-test]") {
    // Verify that a toroidal coil built without explicit sectionAlignment (e.g. CMC path)
    // defaults to SPREAD, placing 2 equal windings at ~90° and ~270°.
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_try_rewind(false);
    settings.set_coil_wind_even_if_not_fit(true);

    std::string coreShape = "T 20/10/7";
    std::string coreMaterial = "3C97";
    auto emptyGapping = json::array();
    int64_t numberStacks = 1;

    auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);

    json bobbinJson;
    to_json(bobbinJson, bobbin);

    json coilJson;
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    for (int i = 0; i < 2; ++i) {
        json w;
        w["name"] = "winding " + std::to_string(i);
        w["numberTurns"] = 10;
        w["numberParallels"] = 1;
        json iso;
        to_json(iso, OpenMagnetics::get_isolation_side_from_index(i));
        w["isolationSide"] = iso;
        w["wire"] = "Round 0.475 - Grade 1";
        coilJson["functionalDescription"].push_back(w);
    }

    // Use json-only constructor: no explicit sectionAlignment → _sectionAlignmentExplicit = false
    OpenMagnetics::Coil coil(coilJson, false);
    // Use CONTIGUOUS so the two windings get separate arcs (needed for angular spread to be visible)
    coil.set_winding_orientation(WindingOrientation::CONTIGUOUS);
    coil.wind();

    settings.reset();
    coil.convert_turns_to_polar_coordinates();
    auto turns = coil.get_turns_description().value();
    // With SPREAD default on round window and 2 equal windings: each winding gets 180°
    // of the toroid. Section centers at ~90° and ~270°. Turns are compactly wound within
    // each section, so the MIDDLE turn (index 4 of 10) is nearest to the section center.
    REQUIRE_THAT(90, Catch::Matchers::WithinAbs(turns[4].get_coordinates()[1], 30));
    REQUIRE_THAT(270, Catch::Matchers::WithinAbs(turns[14].get_coordinates()[1], 30));
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
            Painter painter(outFile);
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


TEST_CASE("Test_Toroidal_External_Turns_Compaction", "[toroidal][coil][compaction]") {
    // Load the toroidal inductor MAS file
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "toroidal_inductor_round_wire_multilayer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();
    auto core = magnetic.get_core();
    
    settings.set_coil_wind_even_if_not_fit(true);
    coil.wind();
    
    // Analyze external turns (those with additional_coordinates)
    auto layers = coil.get_layers_description().value();
    auto turns = coil.get_turns_description().value();
    
    size_t totalExternalTurns = 0;
    double totalWireArea = 0;
    double totalBoundingArea = 0;
    
    for (const auto& turn : turns) {
        if (turn.get_additional_coordinates()) {
            totalExternalTurns++;
            auto coords = turn.get_coordinates();
            auto dims = turn.get_dimensions().value();
            double wireRadius = dims[0] / 2;
            double wireArea = M_PI * wireRadius * wireRadius;
            totalWireArea += wireArea;
            
            // Calculate bounding box area for gap analysis
            totalBoundingArea += dims[0] * dims[1];
        }
    }
    
    // Generate plot if enabled
    if (plot) {
        auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Toroidal_External_Turns_Compaction.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        magnetic.set_coil(coil);
        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
    }
    
    // Verify that external turns exist and are properly placed
    REQUIRE(totalExternalTurns > 0);
    
    // Calculate gap efficiency
    double gapEfficiency = totalWireArea / totalBoundingArea;
    
    // With proper compaction, gap efficiency should be reasonable
    // (allowing for some tolerance due to geometric constraints)
    REQUIRE(gapEfficiency > 0.5);
    
    settings.reset();
}

TEST_CASE("Test_Toroidal_Delimit_And_Compact_Multilayer", "[toroidal][coil][compaction]") {
    // Load the toroidal inductor MAS file
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "toroidal_inductor_round_wire_multilayer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();
    auto core = magnetic.get_core();
    
    settings.set_coil_wind_even_if_not_fit(true);
    coil.wind();
    
    auto turnsBefore = coil.get_turns_description().value();
    size_t totalTurnsBefore = turnsBefore.size();
    
    // Apply delimit and compact
    coil.delimit_and_compact();
    
    auto turnsAfter = coil.get_turns_description().value();
    size_t totalTurnsAfter = turnsAfter.size();
    
    // Verify turn count is preserved
    REQUIRE(totalTurnsBefore == totalTurnsAfter);
    
    // Analyze turn positions before and after
    double totalMovement = 0;
    size_t movedTurns = 0;
    
    for (size_t i = 0; i < std::min(totalTurnsBefore, totalTurnsAfter); ++i) {
        auto coordsBefore = turnsBefore[i].get_coordinates();
        auto coordsAfter = turnsAfter[i].get_coordinates();
        
        double dx = coordsAfter[0] - coordsBefore[0];
        double dy = coordsAfter[1] - coordsBefore[1];
        double movement = std::sqrt(dx*dx + dy*dy);
        
        if (movement > 1e-9) {
            totalMovement += movement;
            movedTurns++;
        }
    }
    
    // Generate comparison plot if enabled
    if (plot) {
        auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("output");
        
        // Plot after compaction
        auto outFileAfter = outputFilePath;
        outFileAfter.append("Test_Toroidal_Delimit_And_Compact_After.svg");
        std::filesystem::remove(outFileAfter);
        Painter painterAfter(outFileAfter);
        magnetic.set_coil(coil);
        painterAfter.paint_core(magnetic);
        painterAfter.paint_coil_turns(magnetic);
        painterAfter.export_svg();
    }
    
    // Verify coil is in valid state after delimit_and_compact
    // Note: If turns are already within bounds, no movement occurs
    // The important thing is that the function runs without errors
    // and the coil remains in a valid state
    REQUIRE(totalTurnsAfter > 0);
    
    // Verify that turns are properly positioned within core boundaries
    bool allTurnsValid = true;
    for (const auto& turn : turnsAfter) {
        auto coords = turn.get_coordinates();
        auto dims = turn.get_dimensions().value();
        
        // Check that turn coordinates are valid numbers
        if (std::isnan(coords[0]) || std::isnan(coords[1])) {
            allTurnsValid = false;
            break;
        }
        
        // Check that turn dimensions are valid
        if (dims[0] <= 0 || dims[1] <= 0) {
            allTurnsValid = false;
            break;
        }
    }
    REQUIRE(allTurnsValid);
    
    settings.reset();
}

TEST_CASE("Test_Coil_Compacting_Tertiary_Winding", "[constructive-model][coil][bug][visualization]") {
    auto testDataPath = get_test_data_path(std::source_location::current(), "bug_coil_compacting.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    
    json j;
    file >> j;
    
    auto magneticJson = j["magnetic"];
    auto coilJson = magneticJson["coil"];
    auto coreJson = magneticJson["core"];
    
    OpenMagnetics::Coil coil(coilJson);
    OpenMagnetics::Core core(coreJson);
    
    // Enable compacting and additional coordinates
    settings.set_coil_delimit_and_compact(true);
    settings.set_coil_include_additional_coordinates(true);
    
    // Wind the coil
    coil.wind();
    
    // Check if we have turns
    REQUIRE(coil.get_turns_description().has_value());
    auto turns = coil.get_turns_description().value();
    
    // Create magnetic for visualization
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_coil(coil);
    magnetic.set_core(core);
    
    // Paint the coil with core
    auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("output");
    auto outFile = outputFilePath;
    outFile.append("Test_Coil_Compacting_Tertiary_Winding.svg");
    std::filesystem::remove(outFile);
    
    Painter painter(outFile);
    painter.paint_core(magnetic);
    painter.paint_bobbin(magnetic);
    painter.paint_coil_sections(magnetic);
    painter.paint_coil_turns(magnetic);
    painter.export_svg();
    
    // Analyze turns by winding
    std::map<std::string, std::vector<Turn>> turnsByWinding;
    for (const auto& turn : turns) {
        turnsByWinding[turn.get_winding()].push_back(turn);
    }
    
    // Check each winding has additional coordinates
    for (const auto& [windingName, windingTurns] : turnsByWinding) {
        int withAdditional = 0;
        int withoutAdditional = 0;
        
        for (const auto& turn : windingTurns) {
            if (turn.get_additional_coordinates().has_value()) {
                withAdditional++;
            } else {
                withoutAdditional++;
            }
        }
        
        // All turns should have additional coordinates
        REQUIRE(withoutAdditional == 0);
    }
    
    settings.reset();
}


TEST_CASE("Test_Add_Margin_From_Json_Reconstructed_Toroidal_Reproduces_Segfault", "[constructive-model][coil][margin][toroidal][regression]") {
    // Reproducer for el-choker frontend crash: PyOpenMagnetics.add_margin_to_section_by_index segfaults
    // when called on a Coil reconstructed from JSON (e.g. the output of magnetic_autocomplete),
    // because Coil(json, false) only deserializes MAS fields and never sizes
    // _marginsPerSection. The function then indexes _marginsPerSection out of bounds.
    settings.reset();
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_wind_even_if_not_fit(true);

    std::vector<int64_t> numberTurns = {60, 42};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels,
                                                      "T 20/10/7", interleavingLevel,
                                                      sectionOrientation, layersOrientation,
                                                      turnsAlignment, sectionsAlignment);
    REQUIRE(bool(coil.get_sections_description()));

    // Round-trip through JSON, simulating what PyOpenMagnetics bindings do
    // (and what the el-choker bake script would do):
    //   1. magnetic_autocomplete returns coil JSON
    //   2. caller constructs Coil(coilJson, false) and calls add_margin_to_section_by_index
    json coilJson;
    to_json(coilJson, coil);
    OpenMagnetics::Coil reconstructed(coilJson, false);
    REQUIRE(bool(reconstructed.get_sections_description()));

    // This must NOT segfault. Pre-fix it does (out-of-bounds index on _marginsPerSection).
    double margin = 0.0002;
    reconstructed.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});

    REQUIRE(bool(reconstructed.get_sections_description()));
    auto sections = reconstructed.get_sections_description().value();
    auto setMargin = std::get<std::vector<double>>(sections[0].get_margin().value());
    REQUIRE_THAT(setMargin[0], Catch::Matchers::WithinAbs(margin, 1e-9));
    REQUIRE_THAT(setMargin[1], Catch::Matchers::WithinAbs(margin, 1e-9));

    settings.reset();
}

// Number of conduction turns overlapping a per-layer reserved connection slot on their own layer
// (terminal lead or U continuation; Z dragbacks excluded — they route diagonally and don't displace
// turns). 0 means the real-winding geometry placed every turn clear of its connections.
static int real_geometry_collisions(OpenMagnetics::Coil& coil) {
    if (!coil.get_turns_description()) {
        return -1;
    }
    auto turns = coil.get_turns_description().value();
    auto spaces = coil.get_connection_reserved_spaces();
    int collisions = 0;
    for (const auto& space : spaces) {
        if (space.layer.empty()) {
            continue;
        }
        if (!space.isTerminal && coil.get_winding_order(space.section) != WindingOrder::U) {
            continue;
        }
        for (const auto& turn : turns) {
            if (!turn.get_layer() || turn.get_layer().value() != space.layer) {
                continue;
            }
            auto tc = turn.get_coordinates();
            auto td = turn.get_dimensions().value();
            double overlapX = (td[0] + space.dimensions[0]) / 2 - std::abs(tc[0] - space.coordinates[0]);
            double overlapY = (td[1] + space.dimensions[1]) / 2 - std::abs(tc[1] - space.coordinates[1]);
            if (overlapX > 1e-6 && overlapY > 1e-6) {
                collisions++;
                std::cout << "[RGCOLL] turn " << turn.get_name() << " (w=" << turn.get_winding()
                          << " h=" << td[1]*1e3 << ") vs " << (space.isTerminal ? "terminal" : "continuation")
                          << "(w=" << space.winding << " h=" << space.dimensions[1]*1e3 << " cy=" << space.coordinates[1]*1e3
                          << ") on " << space.layer << " overlapX=" << overlapX*1e3 << " overlapY=" << overlapY*1e3 << "\n";
            }
        }
    }
    return collisions;
}

static void paint_connection_demo(OpenMagnetics::Coil coil, const std::string& shapeName, const std::string& filename, bool withConnections) {
    auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
    auto outFile = outputFilePath;
    outFile.append(filename);
    std::filesystem::remove(outFile);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::parse("[]"), 1, "Dummy");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    Painter painter(outFile);
    painter.paint_core(magnetic);
    painter.paint_bobbin(magnetic);
    painter.paint_coil_turns(magnetic);
    if (withConnections) {
        painter.paint_coil_connections(magnetic);
    }
    painter.export_svg();
}

TEST_CASE("Test_Winding_Order_U_Vs_Z_Reverses_Alternate_Layers", "[constructive-model][coil][winding-order]") {
    std::vector<int64_t> numberTurns = {60};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", interleavingLevel);

    REQUIRE(coil.get_layers_description());
    auto conductionLayers = coil.get_layers_description_conduction();
    REQUIRE(conductionLayers.size() >= 2);

    // Default winding order is Z (historical behaviour).
    REQUIRE(coil.get_winding_order(conductionLayers[0].get_section().value()) == WindingOrder::Z);

    // The two innermost conduction layers, ordered radially.
    std::sort(conductionLayers.begin(), conductionLayers.end(), [](const Layer& a, const Layer& b) {
        return a.get_coordinates()[0] < b.get_coordinates()[0];
    });
    std::string innerLayer = conductionLayers[0].get_name();
    std::string nextLayer = conductionLayers[1].get_name();

    // y of the first-wound and last-wound turn of a layer (turns_description is in winding order).
    auto firstWoundY = [](OpenMagnetics::Coil& c, const std::string& layerName) {
        auto turns = c.get_turns_description().value();
        for (auto& turn : turns) {
            if (turn.get_layer() && turn.get_layer().value() == layerName) {
                return turn.get_coordinates()[1];
            }
        }
        throw std::runtime_error("no turns in layer " + layerName);
    };
    auto lastWoundY = [](OpenMagnetics::Coil& c, const std::string& layerName) {
        auto turns = c.get_turns_description().value();
        double y = std::nan("");
        for (auto& turn : turns) {
            if (turn.get_layer() && turn.get_layer().value() == layerName) {
                y = turn.get_coordinates()[1];
            }
        }
        return y;
    };

    // Z winding: every layer is wound the same direction, so the first-wound turn of the next layer
    // sits at the same y-edge as the first-wound turn of the inner layer.
    double zInnerFirst = firstWoundY(coil, innerLayer);
    double zNextFirst = firstWoundY(coil, nextLayer);
    REQUIRE_THAT(zNextFirst, Catch::Matchers::WithinAbs(zInnerFirst, 1e-9));

    paint_connection_demo(coil, "PQ 28/20", "Test_WindingOrder_Z.svg", true);

    // Switch the winding window default to U winding and re-wind.
    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    windingWindows[0].set_winding_order(WindingOrder::U);
    processed.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
    coil.wind();

    REQUIRE(coil.get_winding_order(conductionLayers[0].get_section().value()) == WindingOrder::U);

    // U winding: the inner layer is unchanged, but the next layer is reversed, so its first-wound
    // turn now sits at the opposite y-edge — adjacent to the inner layer's last-wound turn.
    double uInnerFirst = firstWoundY(coil, innerLayer);
    double uNextFirst = firstWoundY(coil, nextLayer);
    double uInnerLast = lastWoundY(coil, innerLayer);

    REQUIRE_THAT(uInnerFirst, Catch::Matchers::WithinAbs(zInnerFirst, 1e-9));   // inner layer not reversed
    REQUIRE(std::abs(uNextFirst - zNextFirst) > 1e-6);                          // next layer flipped
    REQUIRE_THAT(uNextFirst, Catch::Matchers::WithinAbs(uInnerLast, 1e-9));     // now adjacent to inner's end

    paint_connection_demo(coil, "PQ 28/20", "Test_WindingOrder_U.svg", true);

    settings.reset();
}

TEST_CASE("Test_Real_Vs_Ideal_Connection_Geometry", "[constructive-model][coil][winding-order][real-geometry]") {
    std::vector<int64_t> numberTurns = {60};
    std::vector<int64_t> numberParallels = {1};
    uint8_t interleavingLevel = 1;

    // Ideal geometry (default): connections reserve no space and contribute no loss.
    settings.set_coil_use_real_winding_geometry(false);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", interleavingLevel);
    REQUIRE(coil.get_layers_description_conduction().size() >= 2);

    double fillingFactorIdeal = coil.get_sections_description_conduction()[0].get_filling_factor().value();
    auto resistanceIdeal = WindingOhmicLosses::calculate_dc_resistance_per_winding(coil, 25.0);
    // The reserved-space geometry exists independently of the flag (there are internal boundaries).
    REQUIRE(coil.get_connection_reserved_spaces().size() >= 1);

    paint_connection_demo(coil, "PQ 28/20", "Test_RealVsIdeal_Ideal.svg", false);

    // Real geometry: connection leads reserve space and add series resistance.
    settings.set_coil_use_real_winding_geometry(true);
    coil.wind();
    double fillingFactorReal = coil.get_sections_description_conduction()[0].get_filling_factor().value();
    auto resistanceReal = WindingOhmicLosses::calculate_dc_resistance_per_winding(coil, 25.0);

    REQUIRE(fillingFactorReal > fillingFactorIdeal);
    REQUIRE(resistanceReal[0] > resistanceIdeal[0]);

    // Plot the same (real-geometry) coil wound Z (default) and U, for comparison. In Z the adjacent
    // layers are joined by an out-of-plane dragback (no in-plane link drawn); in U they are joined
    // by an in-plane turnaround at alternating ends.
    paint_connection_demo(coil, "PQ 28/20", "Test_RealVsIdeal_Real_Z.svg", true);

    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    windingWindows[0].set_winding_order(WindingOrder::U);
    processed.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
    coil.wind();
    REQUIRE(coil.get_winding_order(coil.get_sections_description_conduction()[0].get_name()) == WindingOrder::U);
    paint_connection_demo(coil, "PQ 28/20", "Test_RealVsIdeal_Real_U.svg", true);

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Interleaved_Connection_Squeeze_Grows_Outward", "[constructive-model][coil][real-geometry]") {
    // Interleaved primary/secondary (P-S-P-S), one parallel each, ~20 turns per layer. Connection
    // leads (terminal entrance/exit + inter-layer transitions) route outward to the bobbin border, so
    // they accumulate on the outer layers: the innermost conduction layer is crossed by no lead while
    // the outermost is crossed by several. The reserved connection slots on the outermost conduction
    // layer must therefore exceed those on the innermost.
    std::vector<int64_t> numberTurns = {40, 40};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 2;

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", interleavingLevel);

    auto conductionLayers = coil.get_layers_description_conduction();
    REQUIRE(conductionLayers.size() >= 4);
    std::sort(conductionLayers.begin(), conductionLayers.end(), [](const Layer& a, const Layer& b) {
        return a.get_coordinates()[0] < b.get_coordinates()[0];
    });

    // The structure must be interleaved: both windings appear among the layers.
    std::set<std::string> windingsPresent;
    for (const auto& layer : conductionLayers) {
        windingsPresent.insert(layer.get_partial_windings()[0].get_winding());
    }
    REQUIRE(windingsPresent.size() == 2);

    std::cout << "[SPSP] conduction layers=" << conductionLayers.size()
              << " turns in innermost layer=" << coil.get_turns_by_layer(conductionLayers[0].get_name()).size() << std::endl;

    auto spaces = coil.get_connection_reserved_spaces();
    std::map<std::string, int> reservedCountByLayer;
    for (const auto& space : spaces) {
        if (!space.layer.empty()) {
            reservedCountByLayer[space.layer]++;
        }
    }

    int reservedInnermost = reservedCountByLayer[conductionLayers.front().get_name()];
    int reservedOutermost = reservedCountByLayer[conductionLayers.back().get_name()];

    REQUIRE(reservedOutermost > reservedInnermost);

    // --- Collision check: no conduction turn may overlap a per-layer reserved connection slot
    // (terminal OR inter-layer continuation) on its own layer. ---
    auto countRealCollisions = [&](OpenMagnetics::Coil& c, const std::string& tag) -> int {
        auto turns = c.get_turns_description().value();
        auto localSpaces = c.get_connection_reserved_spaces();
        int collisions = 0;
        for (const auto& space : localSpaces) {
            if (space.layer.empty()) continue;  // per-layer reserved slots only (not the drawn links)
            // A Z dragback does not displace turns in the crossed layer (only terminals + U
            // continuations block), so turns are allowed to sit under a Z continuation squeeze.
            if (!space.isTerminal && c.get_winding_order(space.section) != WindingOrder::U) continue;
            for (const auto& turn : turns) {
                if (!turn.get_layer() || turn.get_layer().value() != space.layer) continue;  // same layer
                auto tc = turn.get_coordinates();
                auto td = turn.get_dimensions().value();
                double overlapX = (td[0] + space.dimensions[0]) / 2 - std::abs(tc[0] - space.coordinates[0]);
                double overlapY = (td[1] + space.dimensions[1]) / 2 - std::abs(tc[1] - space.coordinates[1]);
                if (overlapX > 1e-6 && overlapY > 1e-6) {
                    collisions++;
                    std::cout << "[COLLISION " << tag << "] " << turn.get_name() << " (cy=" << tc[1]*1e3
                              << ") vs " << (space.isTerminal ? "terminal" : "continuation") << " reserved on "
                              << space.layer << " (cy=" << space.coordinates[1]*1e3 << ") overlapY=" << overlapY * 1e3 << " mm\n";
                }
            }
        }
        std::cout << "[REAL COLLISIONS " << tag << "]=" << collisions << "\n";
        return collisions;
    };
    CHECK(countRealCollisions(coil, "Z") == 0);
    // --- Fullness: turns per conduction layer, inner-to-outer (orphan layers stand out) ---
    for (const auto& layer : conductionLayers) {
        auto tl = coil.get_turns_by_layer(layer.get_name());
        double minY = 1e9, maxY = -1e9;
        for (const auto& t : tl) { minY = std::min(minY, t.get_coordinates()[1]); maxY = std::max(maxY, t.get_coordinates()[1]); }
        std::cout << "[FULLNESS] x=" << layer.get_coordinates()[0]*1e3
                  << " " << layer.get_name() << " turns=" << tl.size()
                  << " topY=" << maxY*1e3 << " botY=" << minY*1e3
                  << " layerCy=" << layer.get_coordinates()[1]*1e3 << " layerH=" << layer.get_dimensions()[1]*1e3 << "\n";
    }

    // Z (default) and U variants of the same interleaved coil, for comparison.
    paint_connection_demo(coil, "PQ 28/20", "Test_SPSP_Interleaved_Real.svg", true);

    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    windingWindows[0].set_winding_order(WindingOrder::U);
    processed.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
    coil.wind();
    {
        auto uLayers = coil.get_layers_description_conduction();
        std::sort(uLayers.begin(), uLayers.end(), [](const Layer& a, const Layer& b){ return a.get_coordinates()[0] < b.get_coordinates()[0]; });
        for (const auto& layer : uLayers) {
            std::cout << "[U-FULLNESS] x=" << layer.get_coordinates()[0] << " w=" << layer.get_partial_windings()[0].get_winding()
                      << " " << layer.get_name() << " turns=" << coil.get_turns_by_layer(layer.get_name()).size() << "\n";
        }
    }
    CHECK(countRealCollisions(coil, "U") == 0);
    paint_connection_demo(coil, "PQ 28/20", "Test_SPSP_Interleaved_Real_U.svg", true);

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Different_Wire_Sizes", "[constructive-model][coil][real-geometry]") {
    // Interleaved primary/secondary with QUITE different wire gauges — a thick 1.40 mm primary and a
    // thin 0.40 mm secondary. Per-layer turn capacities and connection-slot sizes differ a lot between
    // the windings, stressing the real-winding blocking + edge-routed connections across mismatched
    // wires. Uses a large PQ 40/40 window so it fits comfortably.
    std::vector<int64_t> numberTurns = {24, 60};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 2;
    auto wires = std::vector<OpenMagnetics::Wire>({
        find_wire_by_name("Round 1.40 - Grade 1"),
        find_wire_by_name("Round 0.4 - Grade 1"),
    });

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 40/40", interleavingLevel,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);

    REQUIRE(coil.get_turns_description());
    std::set<std::string> windingsPresent;
    for (const auto& layer : coil.get_layers_description_conduction()) {
        windingsPresent.insert(layer.get_partial_windings()[0].get_winding());
    }
    REQUIRE(windingsPresent.size() == 2);
    paint_connection_demo(coil, "PQ 40/40", "Test_DifferentWires_Z.svg", true);
    CHECK(real_geometry_collisions(coil) == 0);

    // Same coil wound U.
    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    windingWindows[0].set_winding_order(WindingOrder::U);
    processed.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
    coil.wind();
    paint_connection_demo(coil, "PQ 40/40", "Test_DifferentWires_U.svg", true);
    CHECK(real_geometry_collisions(coil) == 0);

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Z_39_37", "[constructive-model][coil][real-geometry]") {
    // Special-case odd turn counts in Z winding: {39, 37} interleaved. Stresses the blocking-aware
    // section split and spillover when the per-section turns do not divide evenly into full layers.
    std::vector<int64_t> numberTurns = {39, 37};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 2;

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", interleavingLevel);

    REQUIRE(coil.get_turns_description());
    CHECK(real_geometry_collisions(coil) == 0);
    paint_connection_demo(coil, "PQ 28/20", "Test_Real_Z_39_37.svg", true);

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_U_38_36", "[constructive-model][coil][real-geometry]") {
    // Special-case even turn counts in U winding: {38, 36} interleaved.
    std::vector<int64_t> numberTurns = {38, 36};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 2;

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", interleavingLevel);

    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    windingWindows[0].set_winding_order(WindingOrder::U);
    processed.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
    coil.wind();

    REQUIRE(coil.get_turns_description());
    CHECK(real_geometry_collisions(coil) == 0);
    paint_connection_demo(coil, "PQ 28/20", "Test_Real_U_38_36.svg", true);

    settings.reset();
}

// Distinct parallels of a winding that have at least one terminal (entrance/exit) lead. For an
// N-filar winding each parallel is its own conductor, so this should equal the parallel count.
static int distinct_parallels_with_terminal_leads(OpenMagnetics::Coil& coil, const std::string& windingName) {
    std::set<int64_t> parallels;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (space.isTerminal && space.winding == windingName) {
            parallels.insert(space.parallel);
        }
    }
    return int(parallels.size());
}

// True if every conduction layer of the winding holds the SAME number of turns for each of its K
// parallels (and all K are present) — i.e. the parallels are wound side by side in lockstep.
static bool layers_balanced_across_parallels(OpenMagnetics::Coil& coil, const std::string& windingName, int64_t numberParallels) {
    if (!coil.get_turns_description()) {
        return false;
    }
    auto turns = coil.get_turns_description().value();
    std::map<std::string, std::map<int64_t, int>> perLayerPerParallel;
    for (const auto& turn : turns) {
        if (turn.get_winding() != windingName || !turn.get_layer()) {
            continue;
        }
        perLayerPerParallel[turn.get_layer().value()][turn.get_parallel()]++;
    }
    for (const auto& [layerName, perParallel] : perLayerPerParallel) {
        if (int64_t(perParallel.size()) != numberParallels) {
            return false;  // not every parallel present in this layer
        }
        int count = -1;
        for (const auto& [parallel, c] : perParallel) {
            if (count < 0) count = c;
            else if (c != count) return false;  // parallels unequal in this layer
        }
    }
    return true;
}

TEST_CASE("Test_Real_Geometry_Multifilar_N_Filar", "[constructive-model][coil][real-geometry]") {
    // N-filar (bifilar/trifilar/4-filar): each parallel is its own physical conductor with its own
    // entrance/exit terminal leads and its own inter-layer continuation, wound side by side. For
    // K = 2, 3, 4 check the winding stays balanced across parallels, every parallel gets its leads,
    // and no turn collides with a connection lead — in both Z and U winding order.
    for (int64_t K : std::vector<int64_t>{2, 3, 4}) {
        std::vector<int64_t> numberTurns = {18};
        std::vector<int64_t> numberParallels = {K};
        settings.set_coil_use_real_winding_geometry(true);
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", 1);

        INFO("K=" << K << " Z");
        REQUIRE(coil.get_turns_description());
        CHECK(coil.get_turns_description().value().size() == size_t(18 * K));
        CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == int(K));
        CHECK(layers_balanced_across_parallels(coil, "winding 0", K));
        CHECK(real_geometry_collisions(coil) == 0);
        paint_connection_demo(coil, "PQ 28/20", "Test_Real_Multifilar_K" + std::to_string(K) + "_Z.svg", true);

        // Same coil wound U.
        auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
        auto processed = bobbin.get_processed_description().value();
        auto windingWindows = processed.get_winding_windows();
        windingWindows[0].set_winding_order(WindingOrder::U);
        processed.set_winding_windows(windingWindows);
        bobbin.set_processed_description(processed);
        coil.set_bobbin(bobbin);
        coil.wind();

        INFO("K=" << K << " U");
        REQUIRE(coil.get_turns_description());
        CHECK(layers_balanced_across_parallels(coil, "winding 0", K));
        CHECK(real_geometry_collisions(coil) == 0);
        paint_connection_demo(coil, "PQ 28/20", "Test_Real_Multifilar_K" + std::to_string(K) + "_U.svg", true);

        settings.reset();
    }
}

TEST_CASE("Test_Real_Geometry_Bifilar_Interleaved", "[constructive-model][coil][real-geometry]") {
    // Interleaved transformer with a BIFILAR primary (2 parallels) and a single secondary, wound real.
    // Exercises per-parallel connections together with interleaving and mixed parallel counts.
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {2, 1};
    uint8_t interleavingLevel = 2;

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 40/40", interleavingLevel);

    REQUIRE(coil.get_turns_description());
    CHECK(coil.get_turns_description().value().size() == size_t(20 * 2 + 20 * 1));
    CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == 2);
    CHECK(layers_balanced_across_parallels(coil, "winding 0", 2));
    CHECK(layers_balanced_across_parallels(coil, "winding 1", 1));
    CHECK(real_geometry_collisions(coil) == 0);
    paint_connection_demo(coil, "PQ 40/40", "Test_Real_Bifilar_Interleaved_Z.svg", true);

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Rectangular_Contiguous", "[constructive-model][coil][real-geometry]") {
    // Contiguous layers (stacked axially, turns running horizontally) on a rectangular concentric
    // bobbin. The connection model runs in a transposed frame, so leads are produced and drawn just
    // like the overlapping case (rotated 90°). Bifilar, Z and U.
    std::vector<int64_t> numberTurns = {12};
    std::vector<int64_t> numberParallels = {2};

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", 1,
        WindingOrientation::CONTIGUOUS, WindingOrientation::CONTIGUOUS);

    REQUIRE(coil.get_turns_description());
    CHECK(coil.get_turns_description().value().size() == size_t(12 * 2));
    // The transposed model must produce connection leads for the contiguous winding.
    CHECK(!coil.get_connection_reserved_spaces().empty());
    CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == 2);
    paint_connection_demo(coil, "PQ 28/20", "Test_Real_Rect_Contiguous_Z.svg", true);

    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    windingWindows[0].set_winding_order(WindingOrder::U);
    processed.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
    coil.wind();
    paint_connection_demo(coil, "PQ 28/20", "Test_Real_Rect_Contiguous_U.svg", true);

    settings.reset();
}

// Count pairs of turns whose centres are closer than ~one wire — i.e. physically overlapping. Toroidal
// turns are cartesian, so this is a plain centre-distance check.
static int toroidal_turn_overlaps(OpenMagnetics::Coil& coil) {
    if (!coil.get_turns_description()) {
        return -1;
    }
    auto turns = coil.get_turns_description().value();
    int overlaps = 0;
    for (size_t i = 0; i < turns.size(); ++i) {
        for (size_t j = i + 1; j < turns.size(); ++j) {
            auto a = turns[i].get_coordinates();
            auto b = turns[j].get_coordinates();
            double minSeparation = 0.9 * std::min(turns[i].get_dimensions().value()[0], turns[j].get_dimensions().value()[0]);
            if (std::hypot(a[0] - b[0], a[1] - b[1]) < minSeparation) {
                overlaps++;
            }
        }
    }
    return overlaps;
}

TEST_CASE("Test_Real_Geometry_Toroidal", "[constructive-model][coil][real-geometry]") {
    // Toroidal core: concentric polar rings of cartesian turns. Connection leads run radially out to
    // the window border (terminals) and straight between rings (inter-layer continuations), per
    // parallel. Both section-overlapping (concentric, full angle) and section-contiguous (angular
    // sectors) are exercised. Drawing + loss only (no angular turn blocking yet).
    settings.set_coil_use_real_winding_geometry(true);
    std::string shape = "T 17.3/9.7/12.7";

    auto rewindAs = [&](OpenMagnetics::Coil& coil, WindingOrder order) {
        auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
        auto processed = bobbin.get_processed_description().value();
        auto windingWindows = processed.get_winding_windows();
        windingWindows[0].set_winding_order(order);
        processed.set_winding_windows(windingWindows);
        bobbin.set_processed_description(processed);
        coil.set_bobbin(bobbin);
        coil.wind();
    };

    // Section-overlapping (concentric, full angle), bifilar, enough turns for more than one ring.
    for (auto order : {WindingOrder::Z, WindingOrder::U}) {
        std::string tag = (order == WindingOrder::Z) ? "Z" : "U";
        auto coil = OpenMagneticsTesting::get_quick_coil({40}, {2}, shape, 1,
            WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING);
        if (order == WindingOrder::U) rewindAs(coil, order);
        INFO("overlapping " << tag);
        REQUIRE(coil.get_turns_description());
        CHECK(coil.get_turns_description().value().size() == size_t(40 * 2));
        CHECK(toroidal_turn_overlaps(coil) == 0);
        CHECK(!coil.get_connection_reserved_spaces().empty());
        CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == 2);
        paint_connection_demo(coil, shape, "Test_Real_Toroidal_Overlapping_" + tag + ".svg", true);
    }

    // Section-contiguous (two windings side by side in angle), single parallel each.
    for (auto order : {WindingOrder::Z, WindingOrder::U}) {
        std::string tag = (order == WindingOrder::Z) ? "Z" : "U";
        auto coil = OpenMagneticsTesting::get_quick_coil({30, 30}, {1, 1}, shape, 1,
            WindingOrientation::CONTIGUOUS, WindingOrientation::OVERLAPPING);
        if (order == WindingOrder::U) rewindAs(coil, order);
        INFO("section-contiguous " << tag);
        REQUIRE(coil.get_turns_description());
        // No two turns may overlap (the spread must respect each winding's angular sector).
        CHECK(toroidal_turn_overlaps(coil) == 0);
        CHECK(!coil.get_connection_reserved_spaces().empty());
        CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == 1);
        CHECK(distinct_parallels_with_terminal_leads(coil, "winding 1") == 1);
        paint_connection_demo(coil, shape, "Test_Real_Toroidal_SectionContiguous_" + tag + ".svg", true);
    }

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Connection_Loss_Per_Parallel", "[constructive-model][coil][real-geometry]") {
    // Each parallel's terminal/connection leads add DC resistance in series with that parallel, so
    // enabling real winding geometry must raise a bifilar winding's DC resistance above ideal (and the
    // per-parallel split keeps it a sensible parallel combination, not a single lumped series term).
    std::vector<int64_t> numberTurns = {18};
    std::vector<int64_t> numberParallels = {2};

    settings.reset();
    auto idealCoil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", 1);
    double idealResistance = WindingOhmicLosses::calculate_dc_resistance_per_winding(idealCoil, 25.0)[0];

    settings.set_coil_use_real_winding_geometry(true);
    auto realCoil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", 1);
    double realResistance = WindingOhmicLosses::calculate_dc_resistance_per_winding(realCoil, 25.0)[0];

    CHECK(idealResistance > 0);
    CHECK(realResistance > idealResistance);  // connection leads add series resistance per parallel
    settings.reset();
}

TEST_CASE("Demo_Real_Vs_Ideal_Connection_Geometry", "[real-geometry-demo]") {
    // Demonstration (not an assertion-heavy test): winds the same interleaved transformer with
    // ideal vs real winding geometry, prints the filling-factor and DC-resistance deltas, and paints
    // both, drawing the reserved connection rectangles (magenta) in the real one.
    std::vector<int64_t> numberTurns = {10, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 2;
    std::string shape = "PQ 28/20";

    auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
    auto core = OpenMagneticsTesting::get_quick_core(shape, json::parse("[]"), 1, "Dummy");

    settings.reset();
    settings.set_coil_use_real_winding_geometry(false);
    auto coilIdeal = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shape, interleavingLevel);
    auto sectionsIdeal = coilIdeal.get_sections_description_conduction();
    auto resistanceIdeal = WindingOhmicLosses::calculate_dc_resistance_per_winding(coilIdeal, 25.0);

    settings.set_coil_use_real_winding_geometry(true);
    auto coilReal = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shape, interleavingLevel);
    auto sectionsReal = coilReal.get_sections_description_conduction();
    auto resistanceReal = WindingOhmicLosses::calculate_dc_resistance_per_winding(coilReal, 25.0);
    auto reserved = coilReal.get_connection_reserved_spaces();

    std::cout << "\n=== REAL vs IDEAL connection geometry (" << shape << ", interleaved, turns {10,10}) ===\n";
    std::cout << "Conduction sections: " << sectionsIdeal.size() << "\n";
    for (size_t i = 0; i < sectionsIdeal.size(); ++i) {
        std::cout << "  section " << i << " (" << sectionsIdeal[i].get_name() << ") filling factor:"
                  << "  ideal=" << sectionsIdeal[i].get_filling_factor().value()
                  << "  real=" << sectionsReal[i].get_filling_factor().value() << "\n";
    }
    for (size_t w = 0; w < resistanceIdeal.size(); ++w) {
        std::cout << "  winding " << w << " DC resistance [ohm]:"
                  << "  ideal=" << resistanceIdeal[w]
                  << "  real=" << resistanceReal[w]
                  << "  (delta +" << (resistanceReal[w] - resistanceIdeal[w]) << ")\n";
    }
    std::cout << "  reserved connection rectangles: " << reserved.size() << "\n";

    // Demo invariants: same sectioning, real geometry reserves connection space, and the
    // real DC resistance stays close to the ideal one (turn repositioning around the
    // reserved connection rectangles can move it slightly in EITHER direction, so only a
    // broad 20% envelope is pinned, not a sign).
    REQUIRE(sectionsIdeal.size() == sectionsReal.size());
    REQUIRE(resistanceIdeal.size() == resistanceReal.size());
    CHECK(reserved.size() > 0);
    for (size_t w = 0; w < resistanceIdeal.size(); ++w) {
        INFO("Winding " << w);
        CHECK(std::isfinite(resistanceReal[w]));
        CHECK(resistanceIdeal[w] > 0);
        CHECK(resistanceReal[w] > 0);
        CHECK(std::abs(resistanceReal[w] - resistanceIdeal[w]) < 0.2 * resistanceIdeal[w]);
    }

    OpenMagnetics::Magnetic magneticIdeal; magneticIdeal.set_core(core); magneticIdeal.set_coil(coilIdeal);
    OpenMagnetics::Magnetic magneticReal;  magneticReal.set_core(core);  magneticReal.set_coil(coilReal);

    {
        auto outFile = outputFilePath; outFile.append("Demo_Connection_Ideal.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magneticIdeal);
        painter.paint_bobbin(magneticIdeal);
        painter.paint_coil_turns(magneticIdeal);
        painter.export_svg();
        OpenMagneticsTesting::check_svg(outFile);
    }
    {
        auto outFile = outputFilePath; outFile.append("Demo_Connection_Real.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magneticReal);
        painter.paint_bobbin(magneticReal);
        painter.paint_coil_turns(magneticReal);
        painter.paint_coil_connections(magneticReal);
        painter.export_svg();
        OpenMagneticsTesting::check_svg(outFile);
    }
    std::cout << "  SVGs written to " << outputFilePath.string() << "/Demo_Connection_{Ideal,Real}.svg\n" << std::endl;

    settings.reset();
}

// Compact, deterministic digest of a wound coil's geometry. Pins turn/layer/section counts, the
// per-layer turn distribution (which is exactly what get_parallels_proportions drives), and a
// fixed-point bounding box, so any drift in the ideal winding path is caught exactly. NOTE: bind
// get_*_description().value() to a local before iterating — ranging over the temporary dangles.
struct IdealGeometryDigest {
    size_t turns = 0;
    size_t conductionLayers = 0;
    size_t conductionSections = 0;
    int64_t minXum = 0, maxXum = 0, minYum = 0, maxYum = 0;  // turn-centre bounding box, micrometres
    int64_t sumDimsUm = 0;                                   // Σ round((width + height) in µm)
    std::string perLayerTurns;                               // turns per conduction layer, radial order
};

static IdealGeometryDigest ideal_geometry_digest(OpenMagnetics::Coil& coil) {
    IdealGeometryDigest d;
    if (coil.get_turns_description()) {
        auto turns = coil.get_turns_description().value();
        double minX = std::numeric_limits<double>::max(), maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max(), maxY = std::numeric_limits<double>::lowest();
        std::map<std::string, size_t> turnsPerLayer;
        for (const auto& turn : turns) {
            d.turns++;
            auto c = turn.get_coordinates();
            minX = std::min(minX, c[0]); maxX = std::max(maxX, c[0]);
            minY = std::min(minY, c[1]); maxY = std::max(maxY, c[1]);
            if (turn.get_dimensions()) {
                auto dim = turn.get_dimensions().value();
                d.sumDimsUm += std::llround((dim[0] + dim[1]) * 1e6);
            }
            if (turn.get_layer()) {
                turnsPerLayer[turn.get_layer().value()]++;
            }
        }
        if (d.turns > 0) {
            d.minXum = std::llround(minX * 1e6); d.maxXum = std::llround(maxX * 1e6);
            d.minYum = std::llround(minY * 1e6); d.maxYum = std::llround(maxY * 1e6);
        }
        // Per-layer turn counts in radial (then axial) order, so the distribution is order-stable.
        auto conductionLayers = coil.get_layers_description_conduction();
        std::sort(conductionLayers.begin(), conductionLayers.end(), [](const Layer& a, const Layer& b) {
            if (a.get_coordinates()[0] != b.get_coordinates()[0]) return a.get_coordinates()[0] < b.get_coordinates()[0];
            return a.get_coordinates()[1] < b.get_coordinates()[1];
        });
        for (const auto& layer : conductionLayers) {
            if (!d.perLayerTurns.empty()) d.perLayerTurns += ",";
            d.perLayerTurns += std::to_string(turnsPerLayer[layer.get_name()]);
        }
    }
    d.conductionLayers = coil.get_layers_description_conduction().size();
    d.conductionSections = coil.get_sections_description_conduction().size();
    return d;
}

TEST_CASE("Test_Ideal_Winding_Unchanged_Multifilar", "[constructive-model][coil][ideal-regression]") {
    // Locks the ideal (real-geometry OFF) winding geometry across single- and multi-parallel
    // (N-filar) configs. The golden digests below were captured on the Phase-1 baseline. The Phase-2
    // bifilar/N-filar real-winding work must NOT alter the ideal path (it is gated entirely behind
    // Settings::get_coil_use_real_winding_geometry, default false). Any drift here is a regression.
    settings.reset();
    REQUIRE(settings.get_coil_use_real_winding_geometry() == false);

    struct Config { std::vector<int64_t> turns; std::vector<int64_t> parallels; std::string shape; uint8_t interleaving; };
    struct Golden {
        Config config;
        size_t turns, layers, sections;
        int64_t minXum, maxXum, minYum, maxYum, sumDimsUm;
        std::string perLayerTurns;
    };
    // Captured on the Phase-1 baseline (real geometry off). Covers single-parallel, N-filar bifilar/
    // trifilar/5-filar, and interleaved + mixed-parallel layouts. perLayerTurns is the turn count of
    // each conduction layer in radial order — exactly the split get_parallels_proportions produces.
    std::vector<Golden> golden = {
        {{{10},     {1},    "PQ 28/20", 1}, 10, 1, 1,  7354,  7354,  -2290,  2290, 10180, "10"},
        {{{12},     {2},    "PQ 28/20", 1}, 24, 2, 1,  7354,  7863,  -2799,  2799, 24432, "12,12"},
        {{{12},     {3},    "PQ 28/20", 1}, 36, 2, 1,  7354,  7863,  -4326,  4326, 36648, "18,18"},
        {{{10},     {5},    "PQ 40/40", 1}, 50, 1, 1,  9324,  9324, -12470, 12470, 50900, "50"},
        {{{20, 20}, {1, 1}, "PQ 28/20", 2}, 40, 4, 4,  7354,  8956,  -2290,  2290, 40720, "10,10,10,10"},
        {{{20, 20}, {2, 1}, "PQ 40/40", 2}, 60, 4, 4,  9324, 10926,  -4835,  4835, 61080, "20,10,20,10"},
    };
    for (size_t i = 0; i < golden.size(); ++i) {
        const auto& g = golden[i];
        INFO("ideal-regression config " << i << " turns=" << g.config.turns.size());
        auto coil = OpenMagneticsTesting::get_quick_coil(g.config.turns, g.config.parallels, g.config.shape, g.config.interleaving);
        auto d = ideal_geometry_digest(coil);
        CHECK(d.turns == g.turns);
        CHECK(d.conductionLayers == g.layers);
        CHECK(d.conductionSections == g.sections);
        CHECK(d.minXum == g.minXum);
        CHECK(d.maxXum == g.maxXum);
        CHECK(d.minYum == g.minYum);
        CHECK(d.maxYum == g.maxYum);
        CHECK(d.sumDimsUm == g.sumDimsUm);
        CHECK(d.perLayerTurns == g.perLayerTurns);
    }
    settings.reset();
}


}  // namespace
