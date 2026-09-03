#include "RandomUtils.h"
#include <catch2/catch_approx.hpp>
#include <source_location>
#include "support/Settings.h"
#include "support/Painter.h"
#include "constructive_models/Coil.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "processors/Inputs.h"
#include "json.hpp"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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

// ABT #840: MARGIN-FILL ASSERTIONS THAT STATE WHICH CASE THEY TEST.
// The family `marginAfterMarginFill[k] > marginAfterMarginNoFill[k]` guards margin fill, but a
// site whose fixture has NO SLACK (the section already fills its band) passes on float dust --
// one was measured at 8.67e-19 m, an attometre -- and would stay green if margin fill stopped
// working entirely. MKF 416718ac split four such sites by hand; this macro measures every site
// so the split can be made on evidence rather than on a pattern match (the blanket rewrite in
// b54b7d51 is exactly what silently removed a genuine feature test).
//
// It preserves the strict '>' contract and adds two things: the site's measured difference under
// MKF_MARGIN_DIAG, and a floor -- a site claiming fill GROWS the margin must show an
// engineering-scale growth (1 nm), not rounding noise. Sites with no slack assert equality
// instead and do not use this macro.
// ABT #840: the NO-SLACK half of the family. The section already fills its band, so margin fill
// has nothing to absorb and the margin must come back EXACTLY as it was. Measured, not assumed:
// every site converted here was reported by MARGIN_FILL_GROWS at a delta of 4.3e-19 to 8.7e-19 m
// -- attometres, i.e. float dust -- so as a strict '>' it was guarding nothing and would have
// stayed green if margin fill stopped working. As an equality it closes the real hole, which is
// SPURIOUS GROWTH (fill enlarging a margin it must not touch); a shrink was always caught, since
// a >= b is false when a < b.
#define MARGIN_FILL_UNCHANGED(fill, noFill)                                                    \
    do {                                                                                       \
        const double _mfFill = (fill), _mfNoFill = (noFill);                                   \
        if (std::getenv("MKF_MARGIN_DIAG")) {                                                  \
            std::cerr << "[margin-fill] line " << __LINE__ << " (no slack) delta="             \
                      << (_mfFill - _mfNoFill) << "\n";                                        \
        }                                                                                      \
        INFO("margin fill must NOT change the margin at line " << __LINE__ << ": fill "         \
             << _mfFill << " vs no-fill " << _mfNoFill);                                        \
        REQUIRE_THAT(_mfFill, Catch::Matchers::WithinAbs(_mfNoFill, 1e-12));                    \
    } while (0)

#define MARGIN_FILL_GROWS(fill, noFill)                                                        \
    do {                                                                                       \
        const double _mfFill = (fill), _mfNoFill = (noFill);                                   \
        if (std::getenv("MKF_MARGIN_DIAG")) {                                                  \
            std::cerr << "[margin-fill] line " << __LINE__ << " delta="                        \
                      << (_mfFill - _mfNoFill) << "\n";                                        \
        }                                                                                      \
        INFO("margin fill must GROW the margin at line " << __LINE__ << ": fill " << _mfFill    \
             << " vs no-fill " << _mfNoFill << " (delta " << (_mfFill - _mfNoFill) << ")");     \
        REQUIRE(_mfFill > _mfNoFill + 1e-9);                                                   \
    } while (0)

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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_2[1], marginAfterMarginNoFill_2[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_2[0], marginAfterMarginNoFill_2[0]);
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
    // ABT #830/#831/#839 (2026-08-20, deliberate): EQUALITY, not >= and not >. Removing the
    // per-station nanometre rounding in compute_spread_turn_stations (owner-approved; those
    // positions are index-based, so the rounding corrected no accumulated error and only gave
    // back up to a nanometre of the spacing just laid) revealed that margin-fill absorbs
    // NOTHING in this SPREAD fixture — its section already fills its band, so there is no slack
    // to reclaim. The old strict inequality was held up entirely by the rounding artifact
    // (0.002 > 0.002 once removed).
    //
    // Equality is asserted rather than >=, per the ABT #839 hand-off: >= would pass both for the
    // right reason (nothing to absorb) and for a wrong one (fill spuriously ENLARGING a margin it
    // should not touch), so it could not fail if that regressed. The tolerance is float dust, not
    // an engineering band — a real change here is millimetre-scale. The FEATURE (fill genuinely
    // enlarging a margin) is asserted by the _Horizontal_Bottom siblings above, whose sections do
    // not fill their bands and which therefore still assert strict growth.
    REQUIRE_THAT(marginAfterMarginFill[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill[0], 1e-12));
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    REQUIRE_THAT(marginAfterMarginFill_2[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_2[1], 0.0001));
    // ABT #830/#831/#839 (2026-08-20, deliberate): EQUALITY, not >= and not >. Removing the
    // per-station nanometre rounding in compute_spread_turn_stations (owner-approved; those
    // positions are index-based, so the rounding corrected no accumulated error and only gave
    // back up to a nanometre of the spacing just laid) revealed that margin-fill absorbs
    // NOTHING in this SPREAD fixture — its section already fills its band, so there is no slack
    // to reclaim. The old strict inequality was held up entirely by the rounding artifact
    // (0.002 > 0.002 once removed).
    //
    // Equality is asserted rather than >=, per the ABT #839 hand-off: >= would pass both for the
    // right reason (nothing to absorb) and for a wrong one (fill spuriously ENLARGING a margin it
    // should not touch), so it could not fail if that regressed. The tolerance is float dust, not
    // an engineering band — a real change here is millimetre-scale. The FEATURE (fill genuinely
    // enlarging a margin) is asserted by the _Horizontal_Bottom siblings above, whose sections do
    // not fill their bands and which therefore still assert strict growth.
    REQUIRE_THAT(marginAfterMarginFill_0[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_0[0], 1e-12));
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 1e-12));
    REQUIRE_THAT(marginAfterMarginFill_2[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_2[0], 1e-12));
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0]);
    REQUIRE_THAT(marginAfterMarginFill_2[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_2[1], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    // Margin tape can no longer GROW this margin (ABT #579). Filling used to reclaim the dead
    // half-gap that SPREAD parked OUTSIDE the outermost turn at each end of the run: measured on
    // this configuration as 57.78 um, exactly what the margin used to gain (500 -> 557.782 um).
    // Fence-post spread puts that space back BETWEEN the turns, so the section now covers the full
    // span the turns were always allotted and the tape has nothing left to absorb -- the filled
    // margin IS the requested one. Equality within a nanometre; section bounds round to 1 nm.
    REQUIRE_THAT(marginAfterMarginFill[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill[1], 1e-9));
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    // Margin tape can no longer GROW this margin (ABT #579). Filling used to reclaim the dead
    // half-gap that SPREAD parked OUTSIDE the outermost turn at each end of the run: measured on
    // this configuration as 57.78 um, exactly what the margin used to gain (500 -> 557.782 um).
    // Fence-post spread puts that space back BETWEEN the turns, so the section now covers the full
    // span the turns were always allotted and the tape has nothing left to absorb -- the filled
    // margin IS the requested one. Equality within a nanometre; section bounds round to 1 nm.
    REQUIRE_THAT(marginAfterMarginFill_0[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_0[1], 1e-9));
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    // Same fence-post consequence as the _0 margin above (ABT #579): with no dead half-gap left
    // at the ends of the run, filling cannot grow THIS section's own margin either. The
    // cross-section comparisons below (_1 against _0) still hold -- those compare different
    // requested margins, not a margin against itself.
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 1e-9));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    // Margin tape can no longer GROW this margin (ABT #579). Filling used to reclaim the dead
    // half-gap that SPREAD parked OUTSIDE the outermost turn at each end of the run: measured on
    // this configuration as 57.78 um, exactly what the margin used to gain (500 -> 557.782 um).
    // Fence-post spread puts that space back BETWEEN the turns, so the section now covers the full
    // span the turns were always allotted and the tape has nothing left to absorb -- the filled
    // margin IS the requested one. Equality within a nanometre; section bounds round to 1 nm.
    REQUIRE_THAT(marginAfterMarginFill[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill[1], 1e-9));
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    // Margin tape can no longer GROW this margin (ABT #579). Filling used to reclaim the dead
    // half-gap that SPREAD parked OUTSIDE the outermost turn at each end of the run: measured on
    // this configuration as 57.78 um, exactly what the margin used to gain (500 -> 557.782 um).
    // Fence-post spread puts that space back BETWEEN the turns, so the section now covers the full
    // span the turns were always allotted and the tape has nothing left to absorb -- the filled
    // margin IS the requested one. Equality within a nanometre; section bounds round to 1 nm.
    REQUIRE_THAT(marginAfterMarginFill_0[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_0[1], 1e-9));
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    // Same fence-post consequence as the _0 margin above (ABT #579): with no dead half-gap left
    // at the ends of the run, filling cannot grow THIS section's own margin either. The
    // cross-section comparisons below (_1 against _0) still hold -- those compare different
    // requested margins, not a margin against itself.
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 1e-9));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    // Margin tape can no longer GROW this margin (ABT #579). Filling used to reclaim the dead
    // half-gap that SPREAD parked OUTSIDE the outermost turn at each end of the run: measured on
    // this configuration as 57.78 um, exactly what the margin used to gain (500 -> 557.782 um).
    // Fence-post spread puts that space back BETWEEN the turns, so the section now covers the full
    // span the turns were always allotted and the tape has nothing left to absorb -- the filled
    // margin IS the requested one. Equality within a nanometre; section bounds round to 1 nm.
    REQUIRE_THAT(marginAfterMarginFill[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill[1], 1e-9));
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    // Margin tape can no longer GROW this margin (ABT #579). Filling used to reclaim the dead
    // half-gap that SPREAD parked OUTSIDE the outermost turn at each end of the run: measured on
    // this configuration as 57.78 um, exactly what the margin used to gain (500 -> 557.782 um).
    // Fence-post spread puts that space back BETWEEN the turns, so the section now covers the full
    // span the turns were always allotted and the tape has nothing left to absorb -- the filled
    // margin IS the requested one. Equality within a nanometre; section bounds round to 1 nm.
    REQUIRE_THAT(marginAfterMarginFill_0[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_0[1], 1e-9));
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    // Same fence-post consequence as the _0 margin above (ABT #579): with no dead half-gap left
    // at the ends of the run, filling cannot grow THIS section's own margin either. The
    // cross-section comparisons below (_1 against _0) still hold -- those compare different
    // requested margins, not a margin against itself.
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 1e-9));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_GROWS(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[1], marginAfterMarginNoFill[1]);
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_GROWS(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1]);
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 0.0001));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill[0], marginAfterMarginNoFill[0]);
    // Margin tape can no longer GROW this margin (ABT #579). Filling used to reclaim the dead
    // half-gap that SPREAD parked OUTSIDE the outermost turn at each end of the run: measured on
    // this configuration as 57.78 um, exactly what the margin used to gain (500 -> 557.782 um).
    // Fence-post spread puts that space back BETWEEN the turns, so the section now covers the full
    // span the turns were always allotted and the tape has nothing left to absorb -- the filled
    // margin IS the requested one. Equality within a nanometre; section bounds round to 1 nm.
    REQUIRE_THAT(marginAfterMarginFill[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill[1], 1e-9));
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
    // ABT #721: this test pins exact hand-set per-section margins; margin equalization
    // (default-on) would proportionally re-split them — disable to keep the pinned geometry.
    settings.set_coil_equalize_margins(false);

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
    MARGIN_FILL_UNCHANGED(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0]);
    // Margin tape can no longer GROW this margin (ABT #579). Filling used to reclaim the dead
    // half-gap that SPREAD parked OUTSIDE the outermost turn at each end of the run: measured on
    // this configuration as 57.78 um, exactly what the margin used to gain (500 -> 557.782 um).
    // Fence-post spread puts that space back BETWEEN the turns, so the section now covers the full
    // span the turns were always allotted and the tape has nothing left to absorb -- the filled
    // margin IS the requested one. Equality within a nanometre; section bounds round to 1 nm.
    REQUIRE_THAT(marginAfterMarginFill_0[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_0[1], 1e-9));
    REQUIRE_THAT(marginAfterMarginFill_1[0], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[0], 0.0001));
    // Same fence-post consequence as the _0 margin above (ABT #579): with no dead half-gap left
    // at the ends of the run, filling cannot grow THIS section's own margin either. The
    // cross-section comparisons below (_1 against _0) still hold -- those compare different
    // requested margins, not a margin against itself.
    REQUIRE_THAT(marginAfterMarginFill_1[1], Catch::Matchers::WithinAbs(marginAfterMarginNoFill_1[1], 1e-9));
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[0], marginAfterMarginNoFill_0[0]);
    MARGIN_FILL_GROWS(marginAfterMarginFill_1[1], marginAfterMarginNoFill_0[1]);
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

// ABT #720: margins on a MULTI-GROUP toroid. _marginsPerSection is keyed by conduction-section
// ordinal flat across groups; under the old flat interleaved keying (with preload_margins'
// duplicated entries) the second group read the FIRST group's margin rows, so its sections
// carried the wrong margins. One winding per group makes the expected mapping unambiguous.
TEST_CASE("Test_Abt720_MultiGroup_Toroidal_Margins_Keyed_Per_Group", "[constructive-model][coil][round-winding-window][margin][abt720]") {
    settings.reset();
    settings.set_coil_equalize_margins(false);
    std::vector<int64_t> numberTurns = {20, 15};
    std::vector<int64_t> numberParallels = {1, 1};
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "T 20/10/7", 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING);

    // Two explicit groups, one winding each, sharing the toroidal window's angle in halves.
    auto bobbin = coil.resolve_bobbin();
    auto windingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
    std::vector<Group> groups;
    for (size_t windingIndex = 0; windingIndex < 2; ++windingIndex) {
        Group group;
        group.set_name("Group " + std::to_string(windingIndex));
        group.set_winding_window(0);
        group.set_coordinates({windingWindow.get_coordinates().value()[0], windingWindow.get_coordinates().value()[1]});
        group.set_dimensions(std::vector<double>{windingWindow.get_radial_height().value(),
                                                 windingWindow.get_angle().value() / 2});
        group.set_coordinate_system(CoordinateSystem::POLAR);
        group.set_sections_orientation(WindingOrientation::OVERLAPPING);
        group.set_type(WiringTechnology::WOUND);
        PartialWinding partialWinding;
        partialWinding.set_winding(coil.get_functional_description()[windingIndex].get_name());
        partialWinding.set_parallels_proportion(std::vector<double>(1, 1));
        group.set_partial_windings({partialWinding});
        groups.push_back(group);
    }
    coil.set_groups_description(groups);

    const std::vector<double> firstMargin{0.0002, 0.0003};
    const std::vector<double> secondMargin{0.0005, 0.0007};
    coil.reset_margins_per_section();
    coil.preload_margins({firstMargin, secondMargin});
    coil.wind_by_sections();

    REQUIRE(coil.get_sections_description());
    const auto sections = coil.get_sections_description().value();
    size_t checkedSections = 0;
    for (const auto& section : sections) {
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        auto margin = coil.resolve_margin(section);
        const auto& expected = section.get_partial_windings()[0].get_winding() ==
                                       coil.get_functional_description()[0].get_name()
                                   ? firstMargin
                                   : secondMargin;
        INFO(section.get_name() << " margin {" << margin[0] << "," << margin[1] << "}");
        CHECK_THAT(margin[0], Catch::Matchers::WithinRel(expected[0], 1e-9));
        CHECK_THAT(margin[1], Catch::Matchers::WithinRel(expected[1], 1e-9));
        ++checkedSections;
    }
    CHECK(checkedSections == 2);
    settings.reset();
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
    
    // ABT #231 + gap-fill lean (user-mandated, 2026-08-03). History of this contract:
    // the ORIGINAL bug test wanted every outer crossing on ONE radius (maximal fill),
    // reachable only by skewing crossings azimuthally into the gaps — but the then-code
    // skewed NON-MONOTONICALLY (26.1, 41.1, 56.0, 48.5, 63.5 against monotonic inners),
    // crossing consecutive turns' top chords in 3D (the actual #231 defect). The interim
    // re-pin locked outer azimuth == inner azimuth, which killed the crossings but ALSO
    // the fill: a crossing azimuthally aligned with a lower-layer wire perched a full OD
    // on top instead of leaning into the free gap beside it, which no real winding does.
    //
    // The physical contract, asserted here:
    //   (1) every outer crossing lies within the LEAN WINDOW of its turn's azimuth (the
    //       outer leg may lean up to ~one wire OD of arc to rest in a gap);
    //   (1b) outer azimuths stay MONOTONIC in winding order within each layer — the lean
    //       may never reorder crossings (that is what crossed the chords in #231);
    //   (2) no two outer crossings are closer than one wire OD;
    //   (3) compaction genuinely fills: the innermost outer ring is used by more than
    //       one layer, and no crossing sits at or beyond a full OD above the base ring
    //       (full-OD stacking means the lean failed to find the gap that must exist on
    //       this spare-roomed fixture).
    double wireOuterDiameter = coil.get_wires()[0].get_maximum_outer_height();

    std::vector<std::vector<double>> outerCrossings;
    std::vector<double> outerAngles(turns.size()), innerAngles(turns.size());
    double minOuterRadius = std::numeric_limits<double>::max();
    for (size_t i = 0; i < turns.size(); ++i) {
        auto addCoords = turns[i].get_additional_coordinates().value()[0];
        innerAngles[i] = atan2(turns[i].get_coordinates()[1], turns[i].get_coordinates()[0]);
        outerAngles[i] = atan2(addCoords[1], addCoords[0]);
        outerCrossings.push_back({addCoords[0], addCoords[1]});
        minOuterRadius = std::min(minOuterRadius, hypot(addCoords[0], addCoords[1]));
    }
    // (1) near its own azimuth. Two wire ODs, not one (Alf, 2026-08-24): a crossing whose own
    // slot is taken goes to the NEAREST FREE RING SLOT ON EITHER SIDE of the occupant — the WE
    // CMC render request — so a cascade of one occupied slot legitimately displaces it by just
    // over one OD. The binding contracts are (1b) monotonic order, (2) one-OD spacing and
    // (3) everything on (or within one OD of) the base ring; this bound only keeps a crossing
    // from wandering.
    const double leanLimit = 2 * wireOuterDiameter / minOuterRadius + 1e-6;
    for (size_t i = 0; i < turns.size(); ++i) {
        double lean = std::remainder(outerAngles[i] - innerAngles[i], 2 * std::numbers::pi);
        CHECK(std::abs(lean) <= leanLimit);
    }
    // (1b) monotonic per layer, in the layer's own winding direction
    for (auto& [layerName, turnIndices] : turnsByLayer) {
        if (turnIndices.size() < 2) continue;
        double direction = std::remainder(innerAngles[turnIndices[1]] - innerAngles[turnIndices[0]],
                                          2 * std::numbers::pi) < 0 ? -1.0 : 1.0;
        for (size_t k = 0; k + 1 < turnIndices.size(); ++k) {
            double step = std::remainder(outerAngles[turnIndices[k + 1]] - outerAngles[turnIndices[k]],
                                         2 * std::numbers::pi) * direction;
            CHECK(step > 0);
        }
    }
    // (3, second half) on this spare-roomed fixture nothing may sit BEYOND one full OD of the
    // base ring, and the gap-fill must genuinely nest the MAJORITY of the later rings'
    // crossings below the full stack. Crossings azimuthally adjacent to a terminal connection
    // legitimately sit at the full stack: their nested candidates are vetoed so the wire never
    // rests on the connection's radial crossing line.
    size_t beyondBase = 0, fullyStacked = 0;
    for (const auto& crossing : outerCrossings) {
        double crossingRadius = hypot(crossing[0], crossing[1]);
        CHECK(crossingRadius <= minOuterRadius + wireOuterDiameter + 1e-9);
        if (crossingRadius > minOuterRadius + 1e-6) {
            ++beyondBase;
            if (crossingRadius > minOuterRadius + wireOuterDiameter - 1e-6) {
                ++fullyStacked;
            }
        }
    }
    if (beyondBase > 0) {
        CHECK(fullyStacked * 2 < beyondBase + 1);   // majority of raised crossings are NESTED
    }

    // (2) mutual clearance of at least one wire OD
    for (size_t i = 0; i < outerCrossings.size(); ++i) {
        for (size_t j = i + 1; j < outerCrossings.size(); ++j) {
            double separation = hypot(outerCrossings[i][0] - outerCrossings[j][0],
                                      outerCrossings[i][1] - outerCrossings[j][1]);
            CHECK(separation >= wireOuterDiameter - 1e-9);
        }
    }

    // (3) the innermost outer ring is actually used by more than one layer's turns — i.e.
    // crossings compact into the gaps rather than every later ring stacking wholesale.
    std::set<double> uniqueRadii;
    for (const auto& crossing : outerCrossings) {
        uniqueRadii.insert(round(hypot(crossing[0], crossing[1]) * 10000) / 10000);
    }
    double innermostOuterRadius = *uniqueRadii.begin();
    std::set<std::string> layersOnInnermostRing;
    for (size_t i = 0; i < turns.size(); ++i) {
        if (round(hypot(outerCrossings[i][0], outerCrossings[i][1]) * 10000) / 10000
                == innermostOuterRadius) {
            layersOnInnermostRing.insert(turns[i].get_layer().value());
        }
    }
    CHECK(layersOnInnermostRing.size() > 1);

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

TEST_CASE("Test_Toroidal_Rewind_Keeps_Additional_Coordinates", "[constructive-model][coil][round-winding-window][multi-column][bug]") {
    // Winding-studio regression: rewind_layers_and_turns (the custom-rect
    // re-flow) skipped delimit_and_compact_round_window, the only pass that
    // generates the toroidal outer return crossings — every studio edit on a
    // toroid silently dropped all additionalCoordinates.
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_include_additional_coordinates(true);

    std::vector<int64_t> numberTurns = {30, 30};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 1;
    std::string coreShape = "T 25/15/10";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel,
                                                     WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD, CoilAlignment::SPREAD);

    auto woundTurns = coil.get_turns_description().value();
    for (auto& turn : woundTurns) {
        REQUIRE(turn.get_additional_coordinates());
    }
    size_t numberTurnsBefore = woundTurns.size();

    REQUIRE(coil.rewind_layers_and_turns());

    auto rewoundTurns = coil.get_turns_description().value();
    CHECK(rewoundTurns.size() == numberTurnsBefore);
    for (auto& turn : rewoundTurns) {
        REQUIRE(turn.get_coordinate_system());
        CHECK(turn.get_coordinate_system().value() == CoordinateSystem::CARTESIAN);
        REQUIRE(turn.get_additional_coordinates());
        // The outer crossing sits outside the ring: strictly farther from the
        // axis than the in-window position.
        auto coordinates = turn.get_coordinates();
        auto additionalCoordinates = turn.get_additional_coordinates().value();
        REQUIRE(additionalCoordinates.size() >= 1);
        double innerRadius = std::hypot(coordinates[0], coordinates[1]);
        double outerRadius = std::hypot(additionalCoordinates[0][0], additionalCoordinates[0][1]);
        CHECK(outerRadius > innerRadius);
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

// ABT #978 / MAS-RFC 0012: with real winding geometry a planar coil is wound to the PCB reference placement —
// turns left-justified from the column at pcb.designRules.coreToTrack, pitch trackWidth + trackToTrack — and the
// caller's pcb survives the wind on the printed group.
TEST_CASE("Test_Wind_Planar_Real_Winding_Left_Justified_From_Pcb_Rules", "[constructive-model][coil][planar][real-geometry]") {
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(false);
    settings.set_coil_use_real_winding_geometry(true);

    const double trackWidth = 0.001, trackToTrack = 0.0002, coreToTrack = 0.0005;
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(0.01, 0.02);

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(trackWidth);
    wire.set_nominal_value_conducting_height(0.00007);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::PLANAR);

    OpenMagnetics::Coil coil;
    OpenMagnetics::Winding winding;
    winding.set_number_turns(6);
    winding.set_number_parallels(1);
    winding.set_name("Primary");
    winding.set_isolation_side(IsolationSide::PRIMARY);
    winding.set_wire(wire);
    coil.get_mutable_functional_description().push_back(winding);
    coil.set_bobbin(bobbin);

    // caller-provided pcb on a printed group (the group geometry is rebuilt by the wind; the pcb must survive)
    Pcb pcb;
    PcbVias vias; DimensionWithTolerance d; d.set_nominal(0.0004); DimensionWithTolerance dr; dr.set_nominal(0.0003);
    vias.set_diameter(d); vias.set_drill_diameter(dr); pcb.set_vias(vias);
    PcbDesignRules rules; rules.set_track_to_track(trackToTrack); rules.set_core_to_track(coreToTrack); rules.set_via_to_via(0.0003); rules.set_via_to_track(0.0003);
    pcb.set_design_rules(rules);
    PcbOutline outline; outline.set_width(0.06); outline.set_depth(0.04); pcb.set_outline(outline);
    Group group;
    group.set_name("Board"); group.set_type(WiringTechnology::PRINTED); group.set_sections_orientation(WindingOrientation::CONTIGUOUS);
    group.set_partial_windings({}); group.set_dimensions({0.02, 0.01}); group.set_coordinates({0.01, 0});
    group.set_pcb(pcb);
    coil.set_groups_description(std::vector<Group>{group});

    REQUIRE(coil.wind_planar({0, 0}, std::nullopt, {}, {}, defaults.coreToLayerDistance));

    auto groups = coil.get_groups_description().value();
    REQUIRE(groups.size() == 1);
    REQUIRE(groups[0].get_pcb().has_value());
    CHECK(groups[0].get_pcb()->get_design_rules().get_core_to_track() == Catch::Approx(coreToTrack));

    auto layers = coil.get_layers_by_type(ElectricalType::CONDUCTION);
    REQUIRE(layers.size() == 2);
    auto turns = coil.get_turns_description().value();
    REQUIRE(turns.size() == 6);
    // the copper region starts coreToTrack from the column: layer left edge = column edge + coreToTrack
    const double columnEdge = bobbin.get_processed_description()->get_winding_windows()[0].get_coordinates().value()[0]
                            - bobbin.get_processed_description()->get_winding_windows()[0].get_width().value() / 2;
    const double layerLeftEdge = layers[0].get_coordinates()[0] - layers[0].get_dimensions()[0] / 2;
    CHECK(layerLeftEdge == Catch::Approx(columnEdge + coreToTrack).margin(1e-9));
    // first turn hugs that edge; following turns at the PCB pitch
    std::vector<double> layer0Radii;
    for (auto& turn : turns) if (turn.get_layer().value() == layers[0].get_name()) layer0Radii.push_back(turn.get_coordinates()[0]);
    std::sort(layer0Radii.begin(), layer0Radii.end());
    REQUIRE(layer0Radii.size() == 3);
    CHECK(layer0Radii[0] == Catch::Approx(columnEdge + coreToTrack + trackWidth / 2).margin(1e-9));
    CHECK(layer0Radii[1] - layer0Radii[0] == Catch::Approx(trackWidth + trackToTrack).margin(1e-9));
    CHECK(layer0Radii[2] - layer0Radii[1] == Catch::Approx(trackWidth + trackToTrack).margin(1e-9));

    settings.set_coil_use_real_winding_geometry(false);
}

TEST_CASE("Test_Wind_Planar_Real_Winding_Requires_Pcb", "[constructive-model][coil][planar][real-geometry]") {
    settings.set_coil_use_real_winding_geometry(true);
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(0.01, 0.02);
    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.001); wire.set_nominal_value_conducting_height(0.00007);
    wire.set_number_conductors(1); wire.set_material("copper"); wire.set_type(WireType::PLANAR);
    OpenMagnetics::Coil coil;
    OpenMagnetics::Winding winding;
    winding.set_number_turns(2); winding.set_number_parallels(1); winding.set_name("Primary");
    winding.set_isolation_side(IsolationSide::PRIMARY); winding.set_wire(wire);
    coil.get_mutable_functional_description().push_back(winding);
    coil.set_bobbin(bobbin);
    CHECK_THROWS(coil.wind_planar({0}));
    settings.set_coil_use_real_winding_geometry(false);
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

TEST_CASE("Test_Toroidal_Compaction_Syncs_Turn_Rotation_To_Polar_Angle", "[toroidal][coil][compaction]") {
    // ABT #186: angular compaction (delimit_and_compact_round_window) shifts each toroidal turn's
    // polar angle, but historically left turn.rotation at its creation value, so rotation -- the
    // cross-section azimuth read by MagneticField's induced-image rotation and by the painters --
    // went stale. A multi-winding contiguous toroid compacts each section by tens of degrees, so
    // this exercises the shift; guard that rotation == polar angle for every turn afterwards.
    auto angularDiffDeg = [](double a, double b) {
        double d = std::fmod(a - b, 360.0);
        if (d > 180.0) d -= 360.0;
        if (d < -180.0) d += 360.0;
        return std::abs(d);
    };

    std::vector<int64_t> numberTurns = {60, 42};
    std::vector<int64_t> numberParallels = {1, 1};

    // Reference wind WITHOUT compaction: rotation trivially equals the polar angle at creation.
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_delimit_and_compact(false);
    auto coilNoCompact = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "T 20/10/7", 1,
                                                              WindingOrientation::CONTIGUOUS, WindingOrientation::OVERLAPPING,
                                                              CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    coilNoCompact.convert_turns_to_polar_coordinates();
    auto turnsNoCompact = coilNoCompact.get_turns_description().value();
    std::map<std::string, double> uncompactedAngleByName;
    for (const auto& turn : turnsNoCompact) {
        uncompactedAngleByName[turn.get_name()] = turn.get_coordinates()[1];
    }
    settings.reset();

    // Compacted wind: the path that used to leave rotation stale.
    clear_databases();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_delimit_and_compact(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "T 20/10/7", 1,
                                                     WindingOrientation::CONTIGUOUS, WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    coil.convert_turns_to_polar_coordinates();
    auto turns = coil.get_turns_description().value();
    settings.reset();

    REQUIRE(turns.size() == uncompactedAngleByName.size());

    bool rotationMatchesAngle = true;
    size_t shiftedTurns = 0;
    for (const auto& turn : turns) {
        REQUIRE(turn.get_rotation());
        double rotation = turn.get_rotation().value();
        double polarAngle = turn.get_coordinates()[1];
        rotationMatchesAngle &= angularDiffDeg(rotation, polarAngle) < 1e-3;

        // Confirm compaction actually moved this turn's angle relative to the uncompacted wind,
        // so this test genuinely exercises (and would fail on) the stale-rotation bug.
        if (angularDiffDeg(polarAngle, uncompactedAngleByName.at(turn.get_name())) > 1e-3) {
            shiftedTurns++;
        }
    }
    CHECK(rotationMatchesAngle);
    CHECK(shiftedTurns > 0);
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
            // Marker dimensions are {X extent, Y extent} for both layer orientations, the same
            // convention the turns use, so this AABB needs no per-orientation handling.
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

// ABT #492: a Z interleaved inter-section return is manufactured as a DRAGBACK on the core's
// front/back (YZ) face — a radial climb over the intervening sections' build, a near-axial run on
// that face, and a climb back down — consuming NO winding-window space. This aggregates the
// FRONT_YZ segments per (winding, parallel): valid for fixtures whose windings have at most one
// inter-section return each (two sections per winding), which all the fixtures below are.
struct ZDragbackGroup {
    std::string winding;
    int64_t parallel = -1;
    int segments = 0;
    double totalLength = 0;    // sum of each segment's long-axis extent (the FRONT_YZ length rule)
    double bumpRadius = std::numeric_limits<double>::quiet_NaN();  // radial position of the axial run
    double axialRunExtent = 0; // the near-axial run's own axial extent (incl. corner overlaps)
    double wireWidth = 0;      // thickness of the radial climbs (= the run wire's width)
    double wireHeight = 0;     // thickness of the axial run (= the run wire's height)
    // The turn-side end of each radial climb: {radius, axial level}. Two climbs = the return's
    // endpoint turns.
    std::vector<std::pair<double, double>> climbEnds;
};
static std::vector<ZDragbackGroup> front_yz_dragback_groups(OpenMagnetics::Coil& coil) {
    std::map<std::pair<std::string, int64_t>, std::vector<const OpenMagnetics::ConnectionReservedSpace*>> byConductor;
    auto spaces = coil.get_connection_reserved_spaces();
    for (const auto& space : spaces) {
        if (space.plane == RoutePlane::FRONT_YZ) {
            byConductor[{space.winding, space.parallel}].push_back(&space);
        }
    }
    std::vector<ZDragbackGroup> groups;
    for (const auto& [key, segments] : byConductor) {
        ZDragbackGroup group;
        group.winding = key.first;
        group.parallel = key.second;
        double axialRunCenter = std::numeric_limits<double>::quiet_NaN();
        for (const auto* seg : segments) {
            group.segments++;
            group.totalLength += std::max(seg->dimensions[0], seg->dimensions[1]);
            if (seg->dimensions[1] > seg->dimensions[0]) {
                // The near-axial run (thin in x, long in y for overlapping layers).
                group.bumpRadius = seg->coordinates[0];
                group.axialRunExtent = seg->dimensions[1];
                group.wireWidth = seg->dimensions[0];
                axialRunCenter = seg->coordinates[1];
            }
        }
        // Second pass: the climbs' turn-side ends are the ends FARTHER from the bump radius.
        for (const auto* seg : segments) {
            if (seg->dimensions[0] >= seg->dimensions[1]) {
                group.wireHeight = seg->dimensions[1];
                double left = seg->coordinates[0] - seg->dimensions[0] / 2;
                double right = seg->coordinates[0] + seg->dimensions[0] / 2;
                double turnSide = (std::abs(left - group.bumpRadius) > std::abs(right - group.bumpRadius))
                    ? left : right;
                group.climbEnds.push_back({turnSide, seg->coordinates[1]});
            }
        }
        // When the destination layer sits exactly at the bump radius (the usual case: the bump
        // clears the build by half a wire plus insulation, which is precisely where the next
        // section's inner layer starts), the climb-down degenerates and the axial run lands
        // straight on the entry turn — reconstruct that endpoint from the run's far end, stripping
        // the half-wire corner overhangs.
        if (group.segments >= 2 && group.climbEnds.size() == 1 && group.axialRunExtent > 0) {
            double sourceLevel = group.climbEnds[0].second;
            double runLowTurnLevel = axialRunCenter - (group.axialRunExtent - group.wireHeight) / 2;
            double runHighTurnLevel = axialRunCenter + (group.axialRunExtent - group.wireHeight) / 2;
            double destinationLevel = (std::abs(runLowTurnLevel - sourceLevel) > std::abs(runHighTurnLevel - sourceLevel))
                ? runLowTurnLevel : runHighTurnLevel;
            group.climbEnds.push_back({group.bumpRadius, destinationLevel});
        }
        groups.push_back(group);
    }
    return groups;
}

// ABT #492 (b): Z returns must not touch the winding window at all — no in-window squeeze marker
// may name a layer for a non-U continuation (that is what would feed blocking / filling factors).
static int z_return_window_markers(OpenMagnetics::Coil& coil) {
    int markers = 0;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.isTerminal && !space.layer.empty()
            && coil.get_winding_order(space.section) != WindingOrder::U) {
            markers++;
        }
    }
    return markers;
}

// ABT #229: number of pairs of drawn HORIZONTAL lead runs belonging to DIFFERENT conductors
// (winding, parallel) that geometrically overlap — the ticket's exact defect was the K parallels'
// terminal leads all drawn on the SAME edge line (coincident centrelines, which 3D consumers' gates
// throw on). Scope matches the per-edge row allocator's domain: terminal leads (edge-routed runs and
// own-level radial exits) and edge-routed U interleaved continuations (edgeDepth > 0). Vertical
// stubs and U adjacent-layer turnaround links are excluded: those are short link segments at one
// angular position, physically separated azimuthally (out of the 2D plane) in a real multifilar
// winding. 0 means every conductor's horizontal run has its own line.
static int coincident_connection_runs(OpenMagnetics::Coil& coil) {
    auto spaces = coil.get_connection_reserved_spaces();
    std::vector<const OpenMagnetics::ConnectionReservedSpace*> runs;
    for (auto& s : spaces) {
        if (!s.layer.empty()) {
            continue;  // per-layer squeeze markers, not drawn runs
        }
        if (std::abs(s.rotation) > 1e-9) {
            continue;  // diagonal (Z) or polar (toroidal) — not row-allocated
        }
        if (s.dimensions[1] > s.dimensions[0]) {
            continue;  // vertical stub/link — azimuthally separated in 3D, not a row
        }
        if (!s.isTerminal) {
            // ABT #615: inter-section continuation runs SHARE one band by design — in 3D they are
            // tangential chords separated in AZIMUTH, so coincident 2D rectangles are the model,
            // not a defect. Only terminal leads still demand one row per conductor (they fan on
            // one connection plane).
            continue;
        }
        runs.push_back(&s);
    }
    int overlaps = 0;
    for (size_t i = 0; i < runs.size(); ++i) {
        for (size_t j = i + 1; j < runs.size(); ++j) {
            const auto* a = runs[i];
            const auto* b = runs[j];
            if (a->winding == b->winding) {
                // ABT #685 (Alf, superseding the per-conductor rows this helper was written
                // against): a winding's parallels SHARE one terminal row and separate in ANGLE
                // (MVB++'s fan), so coincident 2D rectangles within one winding are the model —
                // the same doctrine the #615 continuation-band exemption above already applies.
                // Distinct WINDINGS still demand distinct rows: nothing separates them in 3D.
                continue;
            }
            double overlapX = (a->dimensions[0] + b->dimensions[0]) / 2 - std::abs(a->coordinates[0] - b->coordinates[0]);
            double overlapY = (a->dimensions[1] + b->dimensions[1]) / 2 - std::abs(a->coordinates[1] - b->coordinates[1]);
            if (overlapX > 1e-6 && overlapY > 1e-6) {
                overlaps++;
                std::cout << "[RUNCOLL] run(w=" << a->winding << " p=" << a->parallel << " c=" << a->coordinates[0]
                          << "," << a->coordinates[1] << ") vs run(w=" << b->winding << " p=" << b->parallel
                          << " c=" << b->coordinates[0] << "," << b->coordinates[1] << ")\n";
            }
        }
    }
    return overlaps;
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

// ABT #430: on the OVERLAPPING path the blocking fixpoint already shrinks a crossed layer by exactly
// the room its leads need ("shrink the layer height by the blocked slots ... leaving room for the
// leads", wind_by_rectangular_layers), so the layer's extent — and the filling factor measured against
// it — ALREADY excludes them. Charging the full lead extent again in apply_connection_reserved_space
// counted the same room twice, which stayed invisible while the leads were thin relative to the layer
// and exploded with fine wire and deep lead stacks (13_current_sense_er95_n87: a correctly-packed
// 0.95-full layer reported at 2.97, so the coil read as not-fitting when it fits).
//
// The invariant asserted here is exact rather than a threshold: an OVERLAPPING layer is one wire WIDE,
// so its filling factor is numberTurns * wire / layerHeight, and blocking guarantees the turns fit in
// the height that remains — hence it cannot exceed 1 unless the leads genuinely could not be given
// room (the shrink is capped at one turn slot minimum, and that residual is charged on purpose). Any
// value above 1 on a coil whose blocking converged is the double charge coming back.
static void check_overlapping_layers_not_double_charged(OpenMagnetics::Coil& coil) {
    // Marker dimensions are {X, Y}; an OVERLAPPING layer's turns stack along Y, so its leads are
    // measured on index 1.
    std::map<std::string, double> reservedPerLayer;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.layer.empty()) {
            reservedPerLayer[space.layer] += space.dimensions[1];
        }
    }
    auto layers = coil.get_layers_description().value();
    int crossed = 0;
    for (const auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION
            || layer.get_orientation() != WindingOrientation::OVERLAPPING
            || !reservedPerLayer.count(layer.get_name())) {
            continue;
        }
        INFO("layer " << layer.get_name() << " leads " << reservedPerLayer.at(layer.get_name())
             << " height " << layer.get_dimensions()[1]
             << " turns " << coil.get_turns_by_layer(layer.get_name()).size());
        CHECK(layer.get_filling_factor().value() <= 1.0);
        crossed++;
    }
    CHECK(crossed > 0);  // the fixture must actually have leads crossing an overlapping layer
}

// ABT #685 (Alf, 2026-08-17): REAL WINDING PLACES ONE CROSSING PER LAYER ENTERED. A wire making
// N turns over L layers crosses the winding-window plane N + L times, not N + 1: each layer's
// helix begins on the plane and ends on it, and the connector joining two layers links the END
// crossing of one to the START crossing of the next AT DIFFERENT RADII -- two distinct wire
// sections that nothing merges. True of a U turnaround and a Z dragback alike, since neither
// encircles the column and so neither consumes a turn.
//
// These helpers state that RULE rather than a number, reading L back from the layout, so a test
// keeps meaning what it says when a design's layer count changes.
static size_t layers_of_conductor(OpenMagnetics::Coil& coil, const std::string& windingName,
                                  int64_t parallel) {
    const auto wound = coil.get_turns_description().value();
    std::set<std::string> layers;
    for (const auto& turn : wound) {
        if (turn.get_winding() == windingName && turn.get_parallel() == parallel && turn.get_layer()) {
            layers.insert(turn.get_layer().value());
        }
    }
    return std::max<size_t>(layers.size(), 1);
}
// Stations a whole winding must show: sum over its parallels of (its turns + the layers it spans).
static size_t expected_stations(OpenMagnetics::Coil& coil, const std::string& windingName,
                                int64_t declaredTurnsPerParallel) {
    const auto windingIndex = coil.get_winding_index_by_name(windingName);
    const int64_t parallels = coil.get_functional_description()[windingIndex].get_number_parallels();
    size_t expected = 0;
    for (int64_t parallel = 0; parallel < parallels; ++parallel) {
        expected += size_t(declaredTurnsPerParallel) + layers_of_conductor(coil, windingName, parallel);
    }
    return expected;
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
        // Real winding: N turns cross the window plane N+1 times, one extra slot per parallel.
        CHECK(coil.get_turns_description().value().size() == expected_stations(coil, "winding 0", 18));
        CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == int(K));
        CHECK(layers_balanced_across_parallels(coil, "winding 0", K));
        CHECK(real_geometry_collisions(coil) == 0);
        // ABT #229: the K parallels' leads must not be drawn on the same line.
        CHECK(coincident_connection_runs(coil) == 0);
        paint_connection_demo(coil, "PQ 28/20", "Test_Real_Multifilar_K" + std::to_string(K) + "_Z.svg", true);
        check_overlapping_layers_not_double_charged(coil);

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
        CHECK(coincident_connection_runs(coil) == 0);
        paint_connection_demo(coil, "PQ 28/20", "Test_Real_Multifilar_K" + std::to_string(K) + "_U.svg", true);
        check_overlapping_layers_not_double_charged(coil);

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
    // Real winding: one extra crossing per parallel ((20+1)*2 + (20+1)*1).
    CHECK(coil.get_turns_description().value().size() ==
          expected_stations(coil, "winding 0", 20) + expected_stations(coil, "winding 1", 20));
    CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == 2);
    CHECK(layers_balanced_across_parallels(coil, "winding 0", 2));
    CHECK(layers_balanced_across_parallels(coil, "winding 1", 1));
    CHECK(real_geometry_collisions(coil) == 0);
    CHECK(coincident_connection_runs(coil) == 0);
    paint_connection_demo(coil, "PQ 40/40", "Test_Real_Bifilar_Interleaved_Z.svg", true);
    check_overlapping_layers_not_double_charged(coil);

    settings.reset();
}

// ABT #492: shared assertions for a Z interleaved inter-section return's YZ-face dragback.
// The return is manufactured as a dragback on the core's front/back face — where there are no
// lateral legs — riding as a local radial bump over the intervening sections' build, so in the XY
// window cross-section it simply is not there. Checks, per dragback: (a) its geometry (endpoints
// land exactly on real turns of its conductor, the bump clears the intervening build by at least
// half the run wire and never overshoots the outer endpoint, the axial run spans the endpoints'
// levels) and (c) its length exceeds the straight diagonal it replaces (the detour is physical and
// the loss model pays for it).
static void check_z_dragback_geometry(OpenMagnetics::Coil& coil, const std::vector<ZDragbackGroup>& groups) {
    auto turns = coil.get_turns_description().value();
    auto conductionLayers = coil.get_layers_description_conduction();
    for (const auto& group : groups) {
        INFO("dragback w=" << group.winding << " p=" << group.parallel << " segments=" << group.segments
             << " bump=" << group.bumpRadius * 1e3 << "mm length=" << group.totalLength * 1e3 << "mm");
        CHECK(group.segments >= 2);
        REQUIRE(std::isfinite(group.bumpRadius));
        REQUIRE(group.climbEnds.size() == 2);
        double x1 = group.climbEnds[0].first;
        double y1 = group.climbEnds[0].second;
        double x2 = group.climbEnds[1].first;
        double y2 = group.climbEnds[1].second;
        // Both endpoints must coincide with actual turns of this conductor — the dragback connects
        // the real exit turn to the real entry turn (this also validates the degenerate-climb
        // reconstruction: the axial run really lands on the entry turn).
        for (const auto& [endRadius, endLevel] : group.climbEnds) {
            bool onTurn = false;
            for (const auto& turn : turns) {
                if (turn.get_winding() == group.winding && turn.get_parallel() == group.parallel
                    && std::abs(turn.get_coordinates()[0] - endRadius) < 1e-6
                    && std::abs(turn.get_coordinates()[1] - endLevel) < 1e-6) {
                    onTurn = true;
                    break;
                }
            }
            CHECK(onTurn);
        }
        // The bump rides over the intervening sections' build: above the outer edge of every
        // conduction layer radially between the endpoints by at least half the run wire (plus any
        // inter-winding insulation, not asserted exactly here), and never beyond the outer
        // endpoint's own radius.
        double radialLow = std::min(x1, x2);
        double radialHigh = std::max(x1, x2);
        double interveningOuter = std::numeric_limits<double>::lowest();
        int intervening = 0;
        for (const auto& layer : conductionLayers) {
            double x = layer.get_coordinates()[0];
            if (x > radialLow + 1e-9 && x < radialHigh - 1e-9) {
                intervening++;
                interveningOuter = std::max(interveningOuter, x + layer.get_dimensions()[0] / 2);
            }
        }
        CHECK(intervening > 0);  // vacuity guard: the return really crosses another section's build
        CHECK(group.bumpRadius >= interveningOuter + group.wireWidth / 2 - 1e-9);
        CHECK(group.bumpRadius <= radialHigh + 1e-9);
        // The axial run spans the endpoints' levels plus the half-wire corner overlaps.
        CHECK_THAT(group.axialRunExtent,
                   Catch::Matchers::WithinAbs(std::abs(y2 - y1) + group.wireHeight, 1e-6));
        // (c) The dragback is strictly longer than the straight in-window diagonal it replaces.
        CHECK(group.totalLength > std::hypot(x2 - x1, y2 - y1));
    }
}

TEST_CASE("Test_Real_Geometry_Z_Interleaved_Return_Dragback", "[constructive-model][coil][real-geometry]") {
    // ABT #492: radially-interleaved sections wound Z. The inter-section return is a DRAGBACK on
    // the core's front/back (YZ) face and consumes NO winding-window space: FRONT_YZ segments are
    // emitted for the loss model and 3D consumers, while the XY layout stays exactly the layout of
    // a coil whose Z returns take nothing from the window (an in-window reservation was measured to
    // run away in the blocking fixpoint — corroborating that the return does not belong there).
    std::vector<int64_t> numberTurns = {40, 40};
    std::vector<int64_t> numberParallels = {1, 1};
    uint8_t interleavingLevel = 2;

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", interleavingLevel);

    REQUIRE(coil.get_turns_description());
    // The fixture must genuinely exercise the class: default Z winding order, both windings
    // interleaved among the layers.
    REQUIRE(coil.get_winding_order(coil.get_sections_description_conduction()[0].get_name()) == WindingOrder::Z);
    std::set<std::string> windingsPresent;
    for (const auto& layer : coil.get_layers_description_conduction()) {
        windingsPresent.insert(layer.get_partial_windings()[0].get_winding());
    }
    REQUIRE(windingsPresent.size() == 2);

    // (a) ABT #615 (supersedes #492): the inter-section return routes IN-WINDOW along the shared
    // top-edge band — no FRONT_YZ face dragbacks exist any more.
    auto groups = front_yz_dragback_groups(coil);
    REQUIRE(groups.empty());

    // (b) The return's corridor now BLOCKS the intervening sections (the missing clearance of ABT
    // #612) — in-window markers exist — but per Alf's ruling it must never charge the return's OWN
    // section or the RECEIVING section, both of which belong to the return's own winding: no
    // non-terminal marker may name a layer owned by the marker's winding.
    CHECK(z_return_window_markers(coil) > 0);
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (space.isTerminal || space.layer.empty()) {
            continue;
        }
        for (const auto& layer : coil.get_layers_description_conduction()) {
            if (layer.get_name() == space.layer) {
                INFO("marker of " << space.winding << " names layer " << space.layer);
                CHECK(layer.get_partial_windings()[0].get_winding() != space.winding);
            }
        }
    }
    // Each winding's single inter-section return draws a run in the shared band: a horizontal
    // non-terminal run with an allocated depth, one per conductor.
    {
        std::set<std::pair<std::string, int64_t>> runConductors;
        for (const auto& space : coil.get_connection_reserved_spaces()) {
            if (!space.isTerminal && space.layer.empty() && space.edgeDepth > 0
                && space.dimensions[0] >= space.dimensions[1]) {
                runConductors.insert({space.winding, space.parallel});
            }
        }
        CHECK(runConductors.size() == 2);
    }

    // (c) The connection resistance includes the dragback: real winding geometry must cost more
    // copper than the ideal wind of the same coil.
    auto resistanceReal = WindingOhmicLosses::calculate_dc_resistance_per_winding(coil, 25.0);
    settings.set_coil_use_real_winding_geometry(false);
    auto idealCoil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", interleavingLevel);
    auto resistanceIdeal = WindingOhmicLosses::calculate_dc_resistance_per_winding(idealCoil, 25.0);
    CHECK(resistanceReal[0] > resistanceIdeal[0]);
    CHECK(resistanceReal[1] > resistanceIdeal[1]);

    settings.set_coil_use_real_winding_geometry(true);
    CHECK(real_geometry_collisions(coil) == 0);
    CHECK(coincident_connection_runs(coil) == 0);
    paint_connection_demo(coil, "PQ 28/20", "Test_Real_Z_Interleaved_Return_Dragback.svg", true);

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Z_Interleaved_Return_Dragback_Bifilar", "[constructive-model][coil][real-geometry]") {
    // ABT #492, same class with a BIFILAR primary: each parallel is its own conductor with its own
    // inter-section return, so each emits its own YZ-face dragback (three in total with the single
    // secondary), each connecting its own parallel's exit/entry turns.
    std::vector<int64_t> numberTurns = {20, 20};
    std::vector<int64_t> numberParallels = {2, 1};
    uint8_t interleavingLevel = 2;

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", interleavingLevel);

    REQUIRE(coil.get_turns_description());
    REQUIRE(coil.get_winding_order(coil.get_sections_description_conduction()[0].get_name()) == WindingOrder::Z);

    // ABT #615: no FRONT_YZ dragbacks; each parallel's return draws its own run in the SHARED
    // in-window band instead (three conductors: primary parallel 0 + parallel 1 + secondary).
    auto groups = front_yz_dragback_groups(coil);
    REQUIRE(groups.empty());
    std::set<std::pair<std::string, int64_t>> conductors;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.isTerminal && space.layer.empty() && space.edgeDepth > 0
            && space.dimensions[0] >= space.dimensions[1]) {
            conductors.insert({space.winding, space.parallel});
        }
    }
    CHECK(conductors.size() == 3);

    CHECK(z_return_window_markers(coil) > 0);
    CHECK(real_geometry_collisions(coil) == 0);
    CHECK(coincident_connection_runs(coil) == 0);
    paint_connection_demo(coil, "PQ 28/20", "Test_Real_Z_Interleaved_Return_Dragback_Bifilar.svg", true);

    settings.reset();
}

// ABT #424 + #427: a connection lead is charged against, and given room out of, the axis the crossed
// layer's turns actually run along. For a CONTIGUOUS layer that is its WIDTH — its height is one wire
// by construction (see wind_by_layers). This pins the AXIS geometrically rather than through the
// filling factor, because once #427 gave contiguous layers turn blocking the crossed layer surrenders
// the room and #430's subtraction cancels the charge, so a wrong axis no longer shows up in the fill
// alone: a crossed layer must be NARROWER than an uncrossed sibling (it gave up width for the leads)
// while being exactly as TALL (its height, the axis that must not be touched, is untouched), and must
// hold correspondingly FEWER turns. Charging the height instead — the mirrored, OVERLAPPING rule —
// made the layer surrender its whole thickness for a single lead and reported ~100x fills (92.5).
static void check_contiguous_lead_reservation(OpenMagnetics::Coil& coil) {
    // Marker dimensions are {X, Y}; a CONTIGUOUS layer's turns run along X, so its leads are measured
    // on index 0 — the mirror of the overlapping helper above, one index apart rather than one frame.
    std::map<std::string, double> reservedPerLayer;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.layer.empty()) {
            reservedPerLayer[space.layer] += space.dimensions[0];
        }
    }
    auto layers = coil.get_layers_description().value();
    std::vector<MAS::Layer> crossed;
    for (const auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION
            || layer.get_orientation() != WindingOrientation::CONTIGUOUS) {
            continue;
        }
        if (reservedPerLayer.count(layer.get_name())) {
            crossed.push_back(layer);
        }
    }
    REQUIRE(!crossed.empty());  // the fixture must exercise a crossed contiguous layer
    // The yardstick is each layer's OWN unblocked geometry, not an untouched sibling. ABT #608 gave
    // the U tangential inter-layer link its own reservation, so in a U winding every layer but the
    // first is crossed and there is no full untouched layer left to compare against — the only
    // uncrossed one is the section's final spillover layer, which is narrow because it holds one
    // turn per parallel, not because nothing crosses it (comparing against THAT is what broke: a
    // 7.03 mm layer measured against a 1.6 mm stub). The two invariants the sibling stood for are
    // absolute anyway, and asserting them directly is stronger: the width came out of the SECTION's
    // width, and the height is exactly one wire, untouched.
    auto sections = coil.get_sections_description().value();
    auto sectionWidthOf = [&](const MAS::Layer& layer) {
        for (const auto& section : sections) {
            if (section.get_name() == layer.get_section().value()) {
                return section.get_dimensions()[0];
            }
        }
        throw std::runtime_error("layer " + layer.get_name() + " names no known section");
    };
    // The coil-wide convention holds for contiguous windings too: a marker that names a layer is an
    // AXIS-ALIGNED rectangle whose dimensions are {X extent, Y extent}, exactly as for overlapping
    // layers and exactly as sections, layers and turns are. The contiguous case is one index apart,
    // not one frame apart — nothing is emitted rotated 90 degrees for a consumer to undo. Only Z
    // diagonals, which are genuinely not axis-aligned, carry an angle, and they name no layer.
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (space.layer.empty()) {
            continue;
        }
        INFO("marker on " << space.layer << " rotation " << space.rotation);
        CHECK(space.rotation == 0.0);
    }
    for (const auto& layer : crossed) {
        size_t windingIndex = coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
        double wireWidth = coil.get_wires()[windingIndex].get_maximum_outer_width();
        size_t markersOnLayer = 0;
        for (const auto& space : coil.get_connection_reserved_spaces()) {
            if (space.layer == layer.get_name()) {
                markersOnLayer++;
            }
        }
        double layerWidth = layer.get_dimensions()[0];
        size_t turnsOnLayer = coil.get_turns_by_layer(layer.get_name()).size();
        INFO("crossed layer " << layer.get_name() << " " << layerWidth << "x" << layer.get_dimensions()[1]
             << " leads " << reservedPerLayer.at(layer.get_name()) << " from " << markersOnLayer
             << " markers, wire width " << wireWidth << ", " << turnsOnLayer << " turns");
        // The leads' room came out of the WIDTH, so the layer is narrower than its section...
        CHECK(layerWidth < sectionWidthOf(layer));
        // ... and NOT out of the height, which stays exactly one wire as it was built.
        CHECK_THAT(layer.get_dimensions()[1],
                   Catch::Matchers::WithinRel(coil.get_wires()[windingIndex].get_maximum_outer_height(),
                                              1e-12));
        // Each lead crossing this layer costs one turn slot measured on the TURN axis, which for a
        // contiguous layer is the wire's OUTER WIDTH. This is the assertion round wire cannot make:
        // with a square marker the width and the height are the same number, so reading the wrong
        // index is invisible. With a non-square wire it is not.
        CHECK_THAT(reservedPerLayer.at(layer.get_name()),
                   Catch::Matchers::WithinRel(double(markersOnLayer) * wireWidth, 1e-9));
        // The turns the layer kept actually fit the width it kept — what the blocking is for.
        CHECK(double(turnsOnLayer) * wireWidth <= layerWidth * (1 + 1e-9));
        // Same invariant the overlapping path satisfies: the layer made room for its leads, so it is
        // not reported as over-full.
        CHECK(layer.get_filling_factor().value() <= 1.0);
    }
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
    // Real winding: one extra crossing per parallel ((12+1)*2).
    CHECK(coil.get_turns_description().value().size() == expected_stations(coil, "winding 0", 12));
    // The transposed model must produce connection leads for the contiguous winding.
    CHECK(!coil.get_connection_reserved_spaces().empty());
    CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == 2);
    paint_connection_demo(coil, "PQ 28/20", "Test_Real_Rect_Contiguous_Z.svg", true);
    check_contiguous_lead_reservation(coil);
    // ABT #427: contiguous layers now shed turn slots for the leads crossing them, so no turn may be
    // left sitting under one — the same guarantee the overlapping path has had since ABT #229.
    CHECK(real_geometry_collisions(coil) == 0);

    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    windingWindows[0].set_winding_order(WindingOrder::U);
    processed.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
    coil.wind();
    paint_connection_demo(coil, "PQ 28/20", "Test_Real_Rect_Contiguous_U.svg", true);
    check_contiguous_lead_reservation(coil);
    CHECK(real_geometry_collisions(coil) == 0);

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Rectangular_Contiguous_Rectangular_Wire", "[constructive-model][coil][real-geometry]") {
    // ABT #427, the case round wire CANNOT test. Every axis decision in the connection model — the
    // turn axis in the blocking and in apply_connection_reserved_space, the layer axis for the lead
    // length in WindingOhmicLosses — reads index 0 or index 1 of a marker's dimensions depending on
    // the layer orientation. Round wire has equal outer width and height, so a crossing marker is
    // SQUARE and both indices hold the same number: swap any of those choices and every round-wire
    // fixture still passes. Rectangular wire makes the two extents differ, so a wrong index changes
    // the result and the assertions below actually bite.
    //
    // It is also the only contiguous coil that reaches wind_by_layers' rectangular branch under
    // blocking, where coil_only_one_turn_per_layer_in_contiguous_rectangular caps a layer at a single
    // turn.
    std::vector<int64_t> numberTurns = {10};
    std::vector<int64_t> numberParallels = {2};

    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(0.00076);
    wire.set_nominal_value_conducting_height(0.00038);
    wire.set_nominal_value_outer_width(0.0008);
    wire.set_nominal_value_outer_height(0.0004);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);

    settings.set_coil_use_real_winding_geometry(true);
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 40/40", 1,
        WindingOrientation::CONTIGUOUS, WindingOrientation::CONTIGUOUS,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED, {wire});

    REQUIRE(coil.get_turns_description());
    // The whole point of this fixture: the wire must NOT be square, or the index choices go untested.
    auto wound = coil.get_wires()[0];
    REQUIRE(std::abs(wound.get_maximum_outer_width() - wound.get_maximum_outer_height()) > 1e-6);

    check_contiguous_lead_reservation(coil);
    CHECK(real_geometry_collisions(coil) == 0);
    paint_connection_demo(coil, "PQ 40/40", "Test_Real_Rect_Contiguous_RectWire_Z.svg", true);

    // Same coil wound U, where the interleaved continuations route along the window edge too.
    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    windingWindows[0].set_winding_order(WindingOrder::U);
    processed.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
    coil.wind();
    REQUIRE(coil.get_turns_description());
    check_contiguous_lead_reservation(coil);
    CHECK(real_geometry_collisions(coil) == 0);
    paint_connection_demo(coil, "PQ 40/40", "Test_Real_Rect_Contiguous_RectWire_U.svg", true);

    settings.reset();
}

// Builds a contiguous, real-geometry coil wound with the given rectangular wire.
static OpenMagnetics::Coil contiguous_rectangular_wire_coil(int64_t turns, int64_t parallels,
                                                            double outerWidth, double outerHeight) {
    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(outerWidth * 0.95);
    wire.set_nominal_value_conducting_height(outerHeight * 0.95);
    wire.set_nominal_value_outer_width(outerWidth);
    wire.set_nominal_value_outer_height(outerHeight);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::RECTANGULAR);
    settings.set_coil_use_real_winding_geometry(true);
    return OpenMagneticsTesting::get_quick_coil({turns}, {parallels}, "PQ 40/40", 1,
        WindingOrientation::CONTIGUOUS, WindingOrientation::CONTIGUOUS,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED, {wire});
}

static std::map<std::string, int> markers_per_layer(OpenMagnetics::Coil& coil) {
    // ABT #685: a winding's parallels share one terminal row, so their squeeze markers are
    // COINCIDENT rectangles — one physical row. Count distinct rectangles, mirroring
    // apply_connection_reserved_space's charge (the slope these tests measure is fill per
    // PHYSICAL crossing, and a duplicate marker is not a crossing).
    std::map<std::string, int> markers;
    std::set<std::tuple<std::string, double, double, double, double>> seen;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.layer.empty()) {
            auto key = std::make_tuple(space.layer,
                                       roundFloat(space.coordinates[0], 9),
                                       roundFloat(space.coordinates[1], 9),
                                       roundFloat(space.dimensions[0], 9),
                                       roundFloat(space.dimensions[1], 9));
            if (seen.insert(key).second) {
                markers[space.layer]++;
            }
        }
    }
    return markers;
}

TEST_CASE("Test_Real_Geometry_Contiguous_Blocked_Room_Is_Measured_In_Turn_Axis_Wire", "[constructive-model][coil][real-geometry]") {
    // ABT #449, gap 2: the slot pitch wind_by_rectangular_layers surrenders room in. For a CONTIGUOUS
    // layer a turn slot is one wire OUTER WIDTH (turns run along the width); using the outer height
    // would be the overlapping rule. That choice was invisible to every earlier fixture — round wire
    // makes the two equal, and with more than one parallel the capacity rounds down to a multiple of
    // the parallel count, which absorbs a one-slot error. A SINGLE parallel and a 2:1 wire leave it
    // nowhere to hide.
    //
    // Both layers here are saturated, so their compacted widths are their allotted widths, and the
    // crossed one is narrower by exactly the room it gave up.
    auto coil = contiguous_rectangular_wire_coil(20, 1, 0.0008, 0.0004);
    REQUIRE(coil.get_turns_description());
    auto wire = coil.get_wires()[0];
    double wireWidth = wire.get_maximum_outer_width();
    double wireHeight = wire.get_maximum_outer_height();
    REQUIRE(std::abs(wireWidth - wireHeight) > 1e-6);  // a square wire could not tell the two apart

    auto markers = markers_per_layer(coil);
    auto layers = coil.get_layers_description().value();
    std::optional<MAS::Layer> crossed;
    std::optional<MAS::Layer> uncrossed;
    for (const auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        REQUIRE(layer.get_orientation() == WindingOrientation::CONTIGUOUS);
        if (markers.count(layer.get_name())) {
            if (!crossed) crossed = layer;
        }
        else if (!uncrossed) {
            uncrossed = layer;
        }
    }
    REQUIRE(crossed);
    REQUIRE(uncrossed);

    int crossings = markers.at(crossed->get_name());
    double roomGivenUp = uncrossed->get_dimensions()[0] - crossed->get_dimensions()[0];
    INFO("crossed " << crossed->get_dimensions()[0] << " uncrossed " << uncrossed->get_dimensions()[0]
         << " gave up " << roomGivenUp << " for " << crossings << " lead(s); wire " << wireWidth
         << " x " << wireHeight);
    CHECK(roomGivenUp > 0);
    // The room given up is a whole number of TURN-AXIS slots. Rather than pin a tolerance on
    // compacted geometry, assert it lands nearer the width-based figure than the height-based one:
    // exact enough to fail the moment the wrong wire dimension is used, with nothing to tune.
    CHECK(std::abs(roomGivenUp - crossings * wireWidth) < std::abs(roomGivenUp - crossings * wireHeight));

    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Contiguous_Oversubscribed_Charge_Is_Measured_In_Turn_Axis_Wire", "[constructive-model][coil][real-geometry]") {
    // ABT #449, gap 1: the turn axis apply_connection_reserved_space charges leads against. While the
    // blocking can give the leads all the room they need, ABT #430's subtraction cancels the charge
    // and the axis makes no difference to any reported number. It only becomes visible once the
    // one-slot-minimum cap binds and a residual survives — a wide wire and enough parallels that the
    // lead stack is deeper than a layer can surrender.
    //
    // Here every layer is the same size and they differ only in how many leads cross them, so the
    // filling factor rises by exactly one turn slot per crossing: d(fill)/d(crossings) must be
    // wireWidth / layerWidth. Using the height would make it wireHeight / layerWidth.
    auto coil = contiguous_rectangular_wire_coil(8, 4, 0.0035, 0.0008);
    REQUIRE(coil.get_turns_description());
    auto wire = coil.get_wires()[0];
    double wireWidth = wire.get_maximum_outer_width();
    REQUIRE(std::abs(wireWidth - wire.get_maximum_outer_height()) > 1e-6);

    auto markers = markers_per_layer(coil);
    auto layers = coil.get_layers_description().value();
    // Crossed layers of identical geometry, keyed by how many leads cross them.
    std::map<int, MAS::Layer> byCrossings;
    double layerWidth = 0;
    for (const auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION || !markers.count(layer.get_name())) {
            continue;
        }
        if (layerWidth == 0) {
            layerWidth = layer.get_dimensions()[0];
        }
        if (std::abs(layer.get_dimensions()[0] - layerWidth) < 1e-12) {
            byCrossings.emplace(markers.at(layer.get_name()), layer);
        }
    }
    REQUIRE(byCrossings.size() >= 2);  // need two crossing counts to measure a slope
    auto fewest = byCrossings.begin();
    auto most = std::prev(byCrossings.end());
    // The cap really binds here: the leads could not all be given room, so the layer is honestly
    // reported as over-subscribed rather than made to fit.
    CHECK(most->second.get_filling_factor().value() > 1.0);

    double fillPerCrossing = (most->second.get_filling_factor().value() - fewest->second.get_filling_factor().value())
                           / double(most->first - fewest->first);
    INFO("layers of width " << layerWidth << ": " << fewest->first << " crossings -> "
         << fewest->second.get_filling_factor().value() << ", " << most->first << " crossings -> "
         << most->second.get_filling_factor().value() << "; per crossing " << fillPerCrossing
         << ", one turn slot is " << wireWidth / layerWidth);
    CHECK_THAT(fillPerCrossing, Catch::Matchers::WithinRel(wireWidth / layerWidth, 1e-9));

    settings.reset();
}

// ABT #492 loss-reader audit: expected centerline connection-copper length per (winding,
// parallel), derived INDEPENDENTLY of the markers from turn positions, layer extents, insulation
// sections, the window box and the documented route model — so the reader can no longer be "right"
// by construction (its old index/shape inference misread U stubs as one wire width and tall-wire
// dragback climbs as the wire height). Rectangular windows only; single winding window (the edge
// row allocator is replayed with window index 0). Mirrors the emitter's spec, not its code.
static std::map<std::pair<std::string, int64_t>, double> expected_connection_length_rectangular(OpenMagnetics::Coil& coil) {
    const bool diag = std::getenv("MKF_TEST_DIAG") != nullptr;
    auto piece = [&](const std::string& what, const std::string& winding, double value) {
        if (diag) fprintf(stderr, "[replay] %s %s += %.4f mm\n", winding.c_str(), what.c_str(), value * 1e3);
        return value;
    };
    std::map<std::pair<std::string, int64_t>, double> expected;
    auto turns = coil.get_turns_description().value();
    auto allLayers = coil.get_layers_description().value();
    auto wires = coil.get_wires();

    std::vector<Layer> conductionLayers;
    for (const auto& layer : allLayers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            conductionLayers.push_back(layer);
        }
    }
    REQUIRE(!conductionLayers.empty());
    bool contiguous = (conductionLayers[0].get_orientation() == WindingOrientation::CONTIGUOUS);
    // Virtual frame: the layer axis is x. For contiguous coils transpose turn/layer positions.
    auto vx = [&](const std::vector<double>& c) { return contiguous ? c[1] : c[0]; };
    auto vy = [&](const std::vector<double>& c) { return contiguous ? c[0] : c[1]; };

    auto windingWindow = coil.resolve_bobbin().get_processed_description().value().get_winding_windows()[0];
    double windowCenterLayerAxis = vx(windingWindow.get_coordinates().value());
    double windowCenterTurnAxis = vy(windingWindow.get_coordinates().value());
    double windowSizeLayerAxis = contiguous ? windingWindow.get_height().value() : windingWindow.get_width().value();
    double windowSizeTurnAxis = contiguous ? windingWindow.get_width().value() : windingWindow.get_height().value();
    double windowOuterX = windowCenterLayerAxis + windowSizeLayerAxis / 2;
    double windowTopY = windowCenterTurnAxis + windowSizeTurnAxis / 2;
    double windowBottomY = windowCenterTurnAxis - windowSizeTurnAxis / 2;

    // Per-edge SPAN-AWARE row allocator replay (single window), mirroring the coil's: rows are
    // shared between runs whose radial spans don't overlap (a terminal lead only occupies
    // [connecting turn, border] -- Alf 2026-08-09: it blocks nothing inward of its turn).
    struct ReplayRow { double height; double depthBefore; std::vector<std::pair<double,double>> spans; };
    std::map<int, std::vector<ReplayRow>> replayRows;
    auto allocateEdgeRow = [&](bool atTop, double wireHeight, double spanLo, double spanHi) {
        auto& rows = replayRows[atTop ? 0 : 1];
        ReplayRow* target = nullptr;
        for (auto& row : rows) {
            if (std::abs(row.height - wireHeight) > 1e-12) {
                continue;
            }
            bool overlaps = false;
            for (const auto& [lo, hi] : row.spans) {
                if (spanLo < hi - 1e-12 && lo < spanHi - 1e-12) {
                    overlaps = true;
                    break;
                }
            }
            if (!overlaps) {
                target = &row;
                break;
            }
        }
        if (target == nullptr) {
            double depthBefore = rows.empty() ? 0.0 : rows.back().depthBefore + rows.back().height;
            rows.push_back({wireHeight, depthBefore, {}});
            target = &rows.back();
        }
        target->spans.push_back({spanLo, spanHi});
        return atTop ? windowTopY - target->depthBefore - wireHeight / 2
                     : windowBottomY + target->depthBefore + wireHeight / 2;
    };
    // ABT #615: ALL inter-section continuations of an edge share ONE band (tangential chords
    // separate in azimuth in 3D); the band claims only its runs' spans, registered on its row.
    // Stage 2: the band's EDGE follows the hop's exit turn (top for hop 1, bottom for hop 2...).
    struct ReplayBand { double y; double height; int row; };
    std::map<int, ReplayBand> continuationBands;
    auto continuationRow = [&](bool atTop, double wireHeight, double spanLo, double spanHi) {
        int edge = atTop ? 0 : 1;
        auto found = continuationBands.find(edge);
        if (found == continuationBands.end() || found->second.height + 1e-12 < wireHeight) {
            double y = allocateEdgeRow(atTop, wireHeight, spanLo, spanHi);
            int row = -1;
            for (size_t r = 0; r < replayRows[edge].size(); ++r) {
                double ry = atTop ? windowTopY - replayRows[edge][r].depthBefore - replayRows[edge][r].height / 2
                                  : windowBottomY + replayRows[edge][r].depthBefore + replayRows[edge][r].height / 2;
                if (std::abs(ry - y) < 1e-12) {
                    row = int(r);
                }
            }
            continuationBands[edge] = {y, wireHeight, row};
        }
        else if (found->second.row >= 0) {
            replayRows[edge][found->second.row].spans.push_back({spanLo, spanHi});
        }
        return continuationBands[edge].y;
    };

    // Electrical order of layers, and endpoint turns per (winding|layer, parallel).
    std::map<std::string, size_t> layerElectricalOrder;
    size_t order = 0;
    std::map<std::pair<std::string, int64_t>, Turn> entranceTurn, exitTurn, firstTurnInLayer, lastTurnInLayer;
    std::map<std::pair<std::string, int64_t>, size_t> stationsInLayer;   // ABT #685
    for (const auto& turn : turns) {
        if (turn.get_layer() && !layerElectricalOrder.count(turn.get_layer().value())) {
            layerElectricalOrder[turn.get_layer().value()] = order++;
        }
        auto windingKey = std::make_pair(turn.get_winding(), turn.get_parallel());
        if (!entranceTurn.count(windingKey)) {
            entranceTurn[windingKey] = turn;
        }
        exitTurn[windingKey] = turn;
        if (turn.get_layer()) {
            auto layerKey = std::make_pair(turn.get_layer().value(), turn.get_parallel());
            if (!firstTurnInLayer.count(layerKey)) {
                firstTurnInLayer[layerKey] = turn;
            }
            lastTurnInLayer[layerKey] = turn;
            // ABT #833: this counter was DECLARED AND NEVER FILLED. Every lookup below therefore
            // read 0, realTurnsIn() returned 0, and the ABT #685 "a layer holding a single turn
            // connects as U" override fired on EVERY layer pair -- which forced every adjacent
            // -layer Z return down the orthogonal-L branch and left this helper's own Z-diagonal
            // branch unreachable. The replay then charged |dx| + |dy| where the winder draws the
            // centre-to-centre dragback hypot(dx, dy), over-counting each dragback by the corner
            // it does not cut (0.4932 mm on the PQ 28/20 fixture -> the ~2.6% shortfall the whole
            // ticket was about). Filling it restores the override to the case it was written for.
            ++stationsInLayer[layerKey];
        }
    }

    // GLOBAL terminal pass mirroring the coil: all windings' leads sorted WIDEST-CROSSING FIRST
    // (then insertion order, then nearest-edge within a group) before any row is allocated, so a
    // narrow late lead stacks deep without charging sections it never reaches.
    {
        struct ReplayLead {
            std::string winding;
            int64_t parallel;
            double turnX, turnY;
            size_t crossedCount;
            double edgeDistance;
            double wireW, wireH;
            size_t order;
        };
        std::vector<ReplayLead> leads;
        for (size_t wi = 0; wi < coil.get_functional_description().size(); ++wi) {
            auto wName = coil.get_functional_description()[wi].get_name();
            double wW = contiguous ? wires[wi].get_maximum_outer_height() : wires[wi].get_maximum_outer_width();
            double wH = contiguous ? wires[wi].get_maximum_outer_width() : wires[wi].get_maximum_outer_height();
            for (int64_t parallel = 0; parallel < int64_t(coil.get_number_parallels(wi)); ++parallel) {
                for (bool entrance : {true, false}) {
                    auto key = std::make_pair(wName, parallel);
                    const auto& source = entrance ? entranceTurn : exitTurn;
                    if (!source.count(key)) {
                        continue;
                    }
                    double tx = vx(source.at(key).get_coordinates());
                    double ty = vy(source.at(key).get_coordinates());
                    if (windowOuterX <= tx) {
                        continue;
                    }
                    size_t crossedCount = 0;
                    for (const auto& layer : conductionLayers) {
                        double layerX = vx(layer.get_coordinates());
                        if (layerX > tx + 1e-9 && layerX < windowOuterX) {
                            ++crossedCount;
                        }
                    }
                    leads.push_back({wName, parallel, tx, ty, crossedCount,
                                     ty >= windowCenterTurnAxis ? windowTopY - ty : ty - windowBottomY,
                                     wW, wH, leads.size()});
                }
            }
        }
        std::stable_sort(leads.begin(), leads.end(), [](const ReplayLead& a, const ReplayLead& b) {
            if (a.crossedCount != b.crossedCount) {
                return a.crossedCount > b.crossedCount;
            }
            if (a.order / 1 != b.order / 1 && a.crossedCount == b.crossedCount) {
                // groups keep first-seen order; within identical crossings the nearest edge first
                // is approximated by edgeDistance (exact for the fixtures, which have one lead per
                // group).
                return a.order < b.order;
            }
            return a.edgeDistance < b.edgeDistance;
        });
        for (const auto& lead : leads) {
            double run = windowOuterX - lead.turnX + lead.wireW;
            if (lead.crossedCount == 0) {
                expected[{lead.winding, lead.parallel}] += piece("radial-exit", lead.winding, run);
                continue;
            }
            double edgeY = allocateEdgeRow(lead.turnY >= windowCenterTurnAxis, lead.wireH,
                                           lead.turnX - lead.wireW / 2, windowOuterX);
            double stub = (std::abs(edgeY - lead.turnY) > lead.wireH / 2)
                              ? std::abs(edgeY - lead.turnY) + lead.wireH / 2 : 0.0;
            expected[{lead.winding, lead.parallel}] += piece("lead-stub", lead.winding, stub)
                                                     + piece("lead-run", lead.winding, run);
        }
    }

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto windingName = coil.get_functional_description()[windingIndex].get_name();
        double wireW = contiguous ? wires[windingIndex].get_maximum_outer_height() : wires[windingIndex].get_maximum_outer_width();
        double wireH = contiguous ? wires[windingIndex].get_maximum_outer_width() : wires[windingIndex].get_maximum_outer_height();
        int64_t numberParallels = int64_t(coil.get_number_parallels(windingIndex));

        // Terminal leads are handled in the GLOBAL widest-first pass below.

        // Inter-layer links, in electrical order.
        std::vector<Layer> windingLayers;
        for (const auto& layer : conductionLayers) {
            if (layer.get_partial_windings()[0].get_winding() == windingName) {
                windingLayers.push_back(layer);
            }
        }
        std::sort(windingLayers.begin(), windingLayers.end(), [&](const Layer& a, const Layer& b) {
            return layerElectricalOrder[a.get_name()] < layerElectricalOrder[b.get_name()];
        });
        for (size_t i = 0; i + 1 < windingLayers.size(); ++i) {
            double radialA = vx(windingLayers[i].get_coordinates());
            double radialB = vx(windingLayers[i + 1].get_coordinates());
            std::vector<const Layer*> intervening;
            double interveningOuter = std::numeric_limits<double>::lowest();
            const Layer* outermostIntervening = nullptr;
            for (const auto& crossed : conductionLayers) {
                double radial = vx(crossed.get_coordinates());
                if (radial > std::min(radialA, radialB) + 1e-12 && radial < std::max(radialA, radialB) - 1e-12) {
                    intervening.push_back(&crossed);
                    double radialExtent = contiguous ? crossed.get_dimensions()[1] : crossed.get_dimensions()[0];
                    if (radial + radialExtent / 2 > interveningOuter) {
                        interveningOuter = radial + radialExtent / 2;
                        outermostIntervening = &crossed;
                    }
                }
            }
            const WindingOrder sectionWindingOrder =
                coil.get_winding_order(windingLayers[i].get_section().value());
            for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
                // ABT #685: a layer holding a SINGLE turn of this parallel connects as U to both
                // its neighbours whatever the section's order -- there is nothing to drag back
                // along. Replayed here because this helper derives the geometry independently;
                // reading the kind off the markers instead would make the check a tautology.
                auto realTurnsIn = [&](const std::string& layerName) {
                    const size_t stations = stationsInLayer.count({layerName, parallel})
                                                ? stationsInLayer.at({layerName, parallel}) : 0;
                    return stations >= 2 ? stations - 1 : stations;
                };
                const WindingOrder windingOrder =
                    (realTurnsIn(windingLayers[i].get_name()) <= 1 ||
                     realTurnsIn(windingLayers[i + 1].get_name()) <= 1)
                        ? WindingOrder::U : sectionWindingOrder;
                auto exitKey = std::make_pair(windingLayers[i].get_name(), parallel);
                auto entryKey = std::make_pair(windingLayers[i + 1].get_name(), parallel);
                if (!lastTurnInLayer.count(exitKey) || !firstTurnInLayer.count(entryKey)) {
                    continue;
                }
                double x1 = vx(lastTurnInLayer.at(exitKey).get_coordinates());
                double y1 = vy(lastTurnInLayer.at(exitKey).get_coordinates());
                double x2 = vx(firstTurnInLayer.at(entryKey).get_coordinates());
                double y2 = vy(firstTurnInLayer.at(entryKey).get_coordinates());
                // ABT #685 FINAL LANDING, replayed (ABT #833): when the destination layer holds
                // exactly ONE turn of this parallel and that turn is the conductor's LAST, the wire
                // is not returning anywhere -- it is finishing. One revolution carries it across the
                // window, so the transition is a spiral turn of full-height pitch, drawn (and
                // charged) as the straight centre-to-centre hop, exactly as a dragback is. Without
                // this case the replay billed it as a U turnaround's orthogonal L, which is both a
                // different path and a different length.
                const bool destHoldsOneTurn =
                    firstTurnInLayer.count(entryKey) && lastTurnInLayer.count(entryKey) &&
                    firstTurnInLayer.at(entryKey).get_name() == lastTurnInLayer.at(entryKey).get_name();
                const auto conductorKey = std::make_pair(windingName, parallel);
                const bool destIsConductorExit =
                    exitTurn.count(conductorKey) &&
                    exitTurn.at(conductorKey).get_name() == firstTurnInLayer.at(entryKey).get_name();
                if (intervening.empty() && destHoldsOneTurn && destIsConductorExit) {
                    expected[{windingName, parallel}] += piece("final-landing", windingName, std::hypot(x2 - x1, y2 - y1));
                }
                else if (windingOrder == WindingOrder::Z && intervening.empty()) {
                    expected[{windingName, parallel}] += piece("diagonal", windingName, std::hypot(x2 - x1, y2 - y1));
                }
                else if (!intervening.empty()) {
                    // ABT #615: any inter-section continuation (U or Z) routes in-window along the
                    // SHARED top-edge band: stubs from both end turns to the band + run across.
                    // (The FRONT_YZ face-dragback model for Z is superseded.)
                    double edgeY = continuationRow(y1 >= windowCenterTurnAxis, wireH,
                                                   std::min(x1, x2) - wireW / 2,
                                                   std::max(x1, x2) + wireW / 2);
                    if (std::abs(edgeY - y1) > wireH / 2) {
                        expected[{windingName, parallel}] += piece("band-stub1", windingName, std::abs(edgeY - y1) + wireH / 2);
                    }
                    expected[{windingName, parallel}] += piece("band-run", windingName, std::abs(x2 - x1) + wireW);
                    if (std::abs(edgeY - y2) > wireH / 2) {
                        expected[{windingName, parallel}] += piece("band-stub2", windingName, std::abs(edgeY - y2) + wireH / 2);
                    }
                }
                else {
                    // U adjacent: orthogonal L (horizontal past the corner, vertical pulled back).
                    // ...EXCEPT the same-section U link, which is horizontal ONLY: it is laid as a
                    // tangential run at constant height, and any residual height difference is the
                    // destination turn's own helix, not drawn copper. Mirrors the emitter's own
                    // needVertical guard -- without it this replay charges a stub nobody draws,
                    // which is exactly what the ABT #685 single-turn-layer rule exposed by turning
                    // Z returns into same-section U links.
                    const bool sameSection =
                        windingLayers[i].get_section() == windingLayers[i + 1].get_section();
                    // ABT #833: this mirrors the emitter's guard, so it must read the SAME order the
                    // emitter reads -- the EFFECTIVE one, after the single-turn-layer downgrade
                    // above. Reading the raw section order here while selecting the branch with the
                    // effective order made the two disagree exactly when the downgrade fired: the
                    // emitter laid a same-section U_TANGENTIAL (one constant-height run, the
                    // residual height being the destination turn's own helix, not drawn copper)
                    // while this replay charged that run PLUS a vertical stub nobody draws.
                    bool needVertical = std::abs(y2 - y1) > wireH / 2
                                        && !(windingOrder == WindingOrder::U && sameSection);
                    expected[{windingName, parallel}] += needVertical
                        ? std::abs(x2 - x1) + wireW / 2 : std::abs(x2 - x1);
                    if (needVertical) {
                        expected[{windingName, parallel}] += std::abs(y2 - y1) - wireH / 2;
                    }
                }
            }
        }
    }
    return expected;
}

// Toroidal counterpart: radial terminal runs capped at the bore wall + inter-ring hops.
static std::map<std::pair<std::string, int64_t>, double> expected_connection_length_toroidal(OpenMagnetics::Coil& coil) {
    std::map<std::pair<std::string, int64_t>, double> expected;
    auto turns = coil.get_turns_description().value();
    auto wires = coil.get_wires();
    double maxTurnRadius = 0;
    std::map<std::string, size_t> ringElectricalOrder;
    size_t order = 0;
    std::map<std::string, std::string> ringWinding;
    std::map<std::pair<std::string, int64_t>, Turn> entranceTurn, exitTurn, firstTurnInRing, lastTurnInRing;
    for (const auto& turn : turns) {
        auto c = turn.get_coordinates();
        maxTurnRadius = std::max(maxTurnRadius, std::hypot(c[0], c[1]));
        if (turn.get_layer() && !ringElectricalOrder.count(turn.get_layer().value())) {
            ringElectricalOrder[turn.get_layer().value()] = order++;
            ringWinding[turn.get_layer().value()] = turn.get_winding();
        }
        auto windingKey = std::make_pair(turn.get_winding(), turn.get_parallel());
        if (!entranceTurn.count(windingKey)) {
            entranceTurn[windingKey] = turn;
        }
        exitTurn[windingKey] = turn;
        if (turn.get_layer()) {
            auto ringKey = std::make_pair(turn.get_layer().value(), turn.get_parallel());
            if (!firstTurnInRing.count(ringKey)) {
                firstTurnInRing[ringKey] = turn;
            }
            lastTurnInRing[ringKey] = turn;
        }
    }
    auto windingWindows = coil.resolve_bobbin().get_processed_description().value().get_winding_windows();
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto windingName = coil.get_functional_description()[windingIndex].get_name();
        double wireW = wires[windingIndex].get_maximum_outer_width();
        double border = maxTurnRadius + 1.5 * wireW;
        if (!windingWindows.empty() && windingWindows[0].get_radial_height()) {
            border = std::min(border, windingWindows[0].get_radial_height().value());
        }
        int64_t numberParallels = int64_t(coil.get_number_parallels(windingIndex));
        for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
            auto key = std::make_pair(windingName, parallel);
            if (entranceTurn.count(key)) {
                auto c = entranceTurn.at(key).get_coordinates();
                double radius = std::hypot(c[0], c[1]);
                if (radius > 1e-9 && border > radius) {
                    expected[key] += border - (radius - wireW / 2);
                }
            }
            if (exitTurn.count(key)) {
                auto c = exitTurn.at(key).get_coordinates();
                double radius = std::hypot(c[0], c[1]);
                if (radius > 1e-9 && border > radius) {
                    expected[key] += border - (radius - wireW / 2);
                }
            }
        }
        // Inter-ring hops in electrical order, per parallel.
        std::vector<std::string> rings;
        for (const auto& [ringName, owner] : ringWinding) {
            if (owner == windingName) {
                rings.push_back(ringName);
            }
        }
        std::sort(rings.begin(), rings.end(), [&](const std::string& a, const std::string& b) {
            return ringElectricalOrder[a] < ringElectricalOrder[b];
        });
        for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
            for (size_t i = 0; i + 1 < rings.size(); ++i) {
                auto exitKey = std::make_pair(rings[i], parallel);
                auto entryKey = std::make_pair(rings[i + 1], parallel);
                if (!lastTurnInRing.count(exitKey) || !firstTurnInRing.count(entryKey)) {
                    continue;
                }
                auto a = lastTurnInRing.at(exitKey).get_coordinates();
                auto b = firstTurnInRing.at(entryKey).get_coordinates();
                double length = std::hypot(b[0] - a[0], b[1] - a[1]);
                if (length > 1e-9) {
                    expected[{windingName, parallel}] += length;
                }
            }
        }
    }
    return expected;
}

// Checks the connection resistance of every (winding, parallel) against perMeter x the
// geometry-derived expected length.
static void check_connection_resistance_matches(OpenMagnetics::Coil& coil,
                                                const std::map<std::pair<std::string, int64_t>, double>& expected) {
    auto wires = coil.get_wires();
    auto connectionResistance = OpenMagnetics::WindingOhmicLosses::calculate_connection_resistance_per_winding_per_parallel(coil, 25.0);
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto windingName = coil.get_functional_description()[windingIndex].get_name();
        double perMeter = OpenMagnetics::WindingOhmicLosses::calculate_dc_resistance_per_meter(wires[windingIndex], 25.0);
        for (size_t parallel = 0; parallel < connectionResistance[windingIndex].size(); ++parallel) {
            auto found = expected.find({windingName, int64_t(parallel)});
            REQUIRE(found != expected.end());
            INFO("winding " << windingName << " parallel " << parallel
                 << " expectedLength=" << found->second * 1e3 << " mm");
            if (std::getenv("MKF_TEST_DIAG")) {
                fprintf(stderr, "[code] %s p%zu length=%.4f mm  (replay %.4f)\n", windingName.c_str(),
                        parallel, connectionResistance[windingIndex][parallel] / perMeter * 1e3,
                        found->second * 1e3);
            }
            CHECK(found->second > 0);  // vacuity: every conductor must route some connection copper
            CHECK_THAT(connectionResistance[windingIndex][parallel],
                       Catch::Matchers::WithinRel(perMeter * found->second, 1e-6));
        }
    }
}

TEST_CASE("Test_Real_Geometry_Connection_Resistance_Is_Centerline_Geometry", "[constructive-model][coil][real-geometry]") {
    // ABT #492 loss-reader audit: connection resistance must equal resistance-per-metre times the
    // TRUE centerline copper length of every routed segment — derived here independently from turn
    // positions, layer builds, bump radii and edge routes — for every marker class and both layer
    // orientations plus toroidal, so any regression back to index/shape inference (which misread U
    // stubs as one wire width and tall-wire dragback climbs as the wire height) is caught.
    // (Supersedes the section-orientation lead-length test: with routedLength carried explicitly
    // by every emitter, there is no orientation-derived axis left to follow.)

    SECTION("overlapping round wire, interleaved Z (terminals + adjacent dragbacks + YZ dragbacks)") {
        settings.set_coil_use_real_winding_geometry(true);
        auto coil = OpenMagneticsTesting::get_quick_coil({40, 40}, {1, 1}, "PQ 28/20", 2);
        REQUIRE(coil.get_turns_description());
        check_connection_resistance_matches(coil, expected_connection_length_rectangular(coil));
        settings.reset();
    }

    SECTION("overlapping TALL rectangular wire (climb run shorter than the wire height)") {
        // The case any shape-based inference misreads: the dragback's radial climbs span the
        // intervening build (~2 wire widths) while the wire's own HEIGHT is much larger, so
        // "the longer dimension" reads the wire height instead of the climb and overcounts.
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.00045);
        wire.set_nominal_value_conducting_height(0.0019);
        wire.set_nominal_value_outer_width(0.0005);
        wire.set_nominal_value_outer_height(0.002);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        settings.set_coil_use_real_winding_geometry(true);
        auto coil = OpenMagneticsTesting::get_quick_coil({4, 4}, {1, 1}, "PQ 28/20", 2,
            WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
            CoilAlignment::CENTERED, CoilAlignment::CENTERED, {wire, wire});
        REQUIRE(coil.get_turns_description());
        // The defect only bites while some segment's true copper is SHORTER than its rectangle's
        // longer side — the exact condition under which the old max() shape rule overcounts. (The
        // rectangle alone cannot even classify such a segment as a climb: a 1.3 mm climb of a
        // 2.0 mm-tall wire is a 1.3 x 2.0 rect whose long side is the wire.)
        // ABT #615 moved the return in-window, so the shape-misread case is now the RUN of a
        // TALL wire: its rectangle is (run x wireHeight) with the wire height the longer side,
        // so "the longer dimension" would overcount exactly as the old climb did.
        bool anyShortRun = false;
        for (const auto& space : coil.get_connection_reserved_spaces()) {
            if (space.routedLength
                && space.routedLength.value() < std::max(space.dimensions[0], space.dimensions[1]) - 1e-9) {
                anyShortRun = true;
            }
        }
        REQUIRE(anyShortRun);
        check_connection_resistance_matches(coil, expected_connection_length_rectangular(coil));
        settings.reset();
    }

    SECTION("contiguous rectangular wire, interleaved Z (the transposed frame)") {
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.00076);
        wire.set_nominal_value_conducting_height(0.00038);
        wire.set_nominal_value_outer_width(0.0008);
        wire.set_nominal_value_outer_height(0.0004);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        settings.set_coil_use_real_winding_geometry(true);
        auto coil = OpenMagneticsTesting::get_quick_coil({10, 10}, {1, 1}, "PQ 40/40", 2,
            WindingOrientation::CONTIGUOUS, WindingOrientation::CONTIGUOUS,
            CoilAlignment::CENTERED, CoilAlignment::CENTERED, {wire, wire});
        REQUIRE(coil.get_turns_description());
        check_connection_resistance_matches(coil, expected_connection_length_rectangular(coil));
        settings.reset();
    }

    SECTION("toroidal (radial leads capped at the bore + inter-ring hops)") {
        settings.reset();
        clear_databases();
        settings.set_use_toroidal_cores(true);
        settings.set_coil_use_real_winding_geometry(true);
        auto coil = OpenMagneticsTesting::get_quick_coil({30, 30}, {1, 1}, "T 40/24/16", 1,
            WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
            CoilAlignment::INNER_OR_TOP, CoilAlignment::INNER_OR_TOP);
        REQUIRE(coil.get_turns_description());
        check_connection_resistance_matches(coil, expected_connection_length_toroidal(coil));
        settings.reset();
    }
}

TEST_CASE("Test_Real_Geometry_Connection_Skin_Losses", "[constructive-model][coil][real-geometry]") {
    // Option A (ABT #492): connection copper (terminal leads, layer-to-layer links, YZ-face
    // dragbacks) gets the isolated-conductor skin correction per harmonic — the same per-meter
    // machinery, above-DC convention and routed lengths as the DC stage, at each parallel's DC
    // current divider — reported on the per-WINDING loss element (connections are not turns, and
    // per-turn consumers map elements by turn name). Proximity stays zero for connections BY
    // DESIGN: those runs are perpendicular to the winding and outside the window field.
    double temperature = 25;

    auto runStages = [&](OpenMagnetics::Coil& coil, double frequency, double& connectionSkinTotal,
                         double& turnSkinTotal, double& totalBeforeSkin, double& totalAfterSkin) {
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
            frequency, 100e-6, temperature, WaveformLabel::TRIANGULAR, 2, 0.5, 0, {1});
        auto operatingPoint = inputs.get_operating_point(0);
        auto ohmicOutput = OpenMagnetics::WindingOhmicLosses::calculate_ohmic_losses(coil, operatingPoint, temperature);
        totalBeforeSkin = ohmicOutput.get_winding_losses();
        auto skinOutput = OpenMagnetics::WindingSkinEffectLosses::calculate_skin_effect_losses(coil, temperature, ohmicOutput);
        totalAfterSkin = skinOutput.get_winding_losses();
        connectionSkinTotal = 0;
        auto lossesPerWinding = skinOutput.get_winding_losses_per_winding().value();
        for (auto& element : lossesPerWinding) {
            REQUIRE(element.get_skin_effect_losses());
            auto connectionElement = element.get_skin_effect_losses().value();
            CHECK(connectionElement.get_method_used() == "ConnectionSkin");
            for (auto loss : connectionElement.get_losses_per_harmonic()) {
                connectionSkinTotal += loss;
            }
        }
        turnSkinTotal = 0;
        auto lossesPerTurn = skinOutput.get_winding_losses_per_turn().value();
        for (auto& element : lossesPerTurn) {
            REQUIRE(element.get_skin_effect_losses());
            auto turnElement = element.get_skin_effect_losses().value();
            for (auto loss : turnElement.get_losses_per_harmonic()) {
                turnSkinTotal += loss;
            }
        }
    };

    // Independent value: geometry-derived connection lengths (turn positions, bump radii, edge
    // routes — NOT the markers) x the per-meter skin machinery with this coil's dividers (single
    // parallel -> exactly 1), same harmonic pruning as the stage.
    auto expectedConnectionSkin = [&](OpenMagnetics::Coil& coil, double frequency) {
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
            frequency, 100e-6, temperature, WaveformLabel::TRIANGULAR, 2, 0.5, 0, {1});
        auto prunedOperatingPoint = OpenMagnetics::Inputs::prune_harmonics(
            inputs.get_operating_point(0), Defaults().harmonicAmplitudeThreshold);
        auto expectedLengths = expected_connection_length_rectangular(coil);
        double expected = 0;
        for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
            auto windingName = coil.get_functional_description()[windingIndex].get_name();
            auto current = prunedOperatingPoint.get_excitations_per_winding()[windingIndex].get_current().value();
            double perMeter = OpenMagnetics::WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(
                coil.resolve_wire(windingIndex), current, temperature, 1).first;
            expected += perMeter * expectedLengths.at({windingName, 0});
        }
        return expected;
    };

    SECTION("overlapping round wire: value from geometry-derived lengths, assembly, monotonicity") {
        settings.set_coil_use_real_winding_geometry(true);
        auto coil = OpenMagneticsTesting::get_quick_coil({40, 40}, {1, 1}, "PQ 28/20", 2);
        REQUIRE(coil.get_turns_description());
        double connectionSkin100k, turnSkin, totalBefore, totalAfter;
        runStages(coil, 100000, connectionSkin100k, turnSkin, totalBefore, totalAfter);
        CHECK(connectionSkin100k > 0);  // wiring guard: the term cannot silently drop out
        CHECK_THAT(connectionSkin100k, Catch::Matchers::WithinRel(expectedConnectionSkin(coil, 100000), 1e-6));  // 1e-6 absorbs the emitters' nanometre roundFloat on marker lengths
        // Assembly invariant: total = ohmic (turns + connections) + skin(turns) + skin(connections).
        CHECK_THAT(totalAfter, Catch::Matchers::WithinRel(totalBefore + turnSkin + connectionSkin100k, 1e-9));
        // Monotonic in frequency.
        double connectionSkin400k, turnSkin400, totalBefore400, totalAfter400;
        runStages(coil, 400000, connectionSkin400k, turnSkin400, totalBefore400, totalAfter400);
        CHECK(connectionSkin400k > connectionSkin100k);
        settings.reset();
    }

    SECTION("tall rectangular wire: the dragback-heavy fixture dispatches the rectangular model") {
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.00045);
        wire.set_nominal_value_conducting_height(0.0019);
        wire.set_nominal_value_outer_width(0.0005);
        wire.set_nominal_value_outer_height(0.002);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        settings.set_coil_use_real_winding_geometry(true);
        auto coil = OpenMagneticsTesting::get_quick_coil({4, 4}, {1, 1}, "PQ 28/20", 2,
            WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
            CoilAlignment::CENTERED, CoilAlignment::CENTERED, {wire, wire});
        REQUIRE(coil.get_turns_description());
        double connectionSkin, turnSkin, totalBefore, totalAfter;
        runStages(coil, 100000, connectionSkin, turnSkin, totalBefore, totalAfter);
        CHECK(connectionSkin > 0);
        CHECK_THAT(connectionSkin, Catch::Matchers::WithinRel(expectedConnectionSkin(coil, 100000), 1e-6));  // 1e-6 absorbs the emitters' nanometre roundFloat on marker lengths
        settings.reset();
    }

    SECTION("ideal mode: zero connection skin, per-winding elements untouched") {
        settings.reset();
        auto coil = OpenMagneticsTesting::get_quick_coil({40, 40}, {1, 1}, "PQ 28/20", 2);
        REQUIRE(coil.get_turns_description());
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
            100000, 100e-6, temperature, WaveformLabel::TRIANGULAR, 2, 0.5, 0, {1});
        auto ohmicOutput = OpenMagnetics::WindingOhmicLosses::calculate_ohmic_losses(coil, inputs.get_operating_point(0), temperature);
        auto skinOutput = OpenMagnetics::WindingSkinEffectLosses::calculate_skin_effect_losses(coil, temperature, ohmicOutput);
        auto lossesPerWinding = skinOutput.get_winding_losses_per_winding().value();
        for (auto& element : lossesPerWinding) {
            CHECK(!element.get_skin_effect_losses());
        }
        settings.reset();
    }
}

TEST_CASE("Test_Abt967_Autocomplete_Accepts_A_Database_Foil", "[constructive-model][coil][foil][regression]") {
    // ABT #967: a database FOIL ("Foil 0.2") states only its thickness (conductingWidth); its
    // height is the section's, cut in at wind() time. The ABT #823 pre-wind outer-size derivation
    // asked such a wire for its outer HEIGHT and dereferenced an empty optional -- every design
    // naming a database foil died in autocomplete with INVALID_WIRE_DATA. Autocomplete must derive
    // only what the wire states, and leave the rest to the wind.
    settings.reset();
    clear_databases();
    auto mas = OpenMagneticsTesting::mas_loader(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'))
                                                + "/../MAS/examples/04_forward_xfmr_e3216_n87.json");
    auto magneticIn = mas.get_magnetic();
    auto& windings = magneticIn.get_mutable_coil().get_mutable_functional_description();
    REQUIRE(windings.size() >= 2);
    windings[1].set_wire("Foil 0.2");
    // A fresh wind, so autocomplete has to complete the wire itself.
    magneticIn.get_mutable_coil().set_turns_description(std::nullopt);
    magneticIn.get_mutable_coil().set_layers_description(std::nullopt);
    magneticIn.get_mutable_coil().set_sections_description(std::nullopt);
    auto magnetic = OpenMagnetics::magnetic_autocomplete(magneticIn);
    auto wire = magnetic.get_mutable_coil().resolve_wire(1);
    REQUIRE(wire.get_type() == WireType::FOIL);
    REQUIRE(wire.get_outer_width());
    REQUIRE_THAT(resolve_dimensional_values(wire.get_outer_width().value()),
                 Catch::Matchers::WithinRel(resolve_dimensional_values(wire.get_conducting_width().value()), 1e-9));
    settings.reset();
}

TEST_CASE("Test_Real_Geometry_Planar_Throws", "[constructive-model][coil][real-geometry]") {
    // ABT #492 owner ruling: planar wires are PCBs — the real-winding connection model (leads,
    // markers, blocking, YZ-face dragbacks, connection losses) is for WOUND magnetics only, and
    // real winding for planar has not been started. Enabling the setting on a planar construction
    // must THROW, loudly, at the first point the machinery would engage — no via model, no
    // fallback, no silent skip. (Production planar flows are unaffected: the setting defaults
    // false.)
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto mas = OpenMagneticsTesting::mas_loader(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'))
                                                + "/../MAS/examples/09_planar_xfmr_er2510_3c94.json");
    // The planar example ships fully wound, so autocomplete does not re-wind it; the gate must fire
    // at the first real-winding machinery a planar coil can reach from there: a re-wind, and the
    // connection-resistance path the loss chain uses (the two entries that consult the setting).
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    REQUIRE_THROWS_WITH(magnetic.get_mutable_coil().wind(),
                        Catch::Matchers::ContainsSubstring("not implemented for planar"));
    REQUIRE_THROWS_WITH(OpenMagnetics::WindingOhmicLosses::calculate_connection_resistance_per_winding_per_parallel(
                            magnetic.get_coil(), 25.0),
                        Catch::Matchers::ContainsSubstring("not implemented for planar"));
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
                std::cout << "[TOVER] " << turns[i].get_name() << " (layer " << turns[i].get_layer().value_or("?")
                          << " az=" << std::atan2(a[1], a[0]) * 180.0 / std::numbers::pi << " r=" << std::hypot(a[0], a[1])
                          << ") vs " << turns[j].get_name() << " (layer " << turns[j].get_layer().value_or("?")
                          << " az=" << std::atan2(b[1], b[0]) * 180.0 / std::numbers::pi << " r=" << std::hypot(b[0], b[1]) << ")\n";
            }
        }
    }
    return overlaps;
}

// ABT #187: number of turns sitting inside a connection lead's angular corridor on their own ring.
// A radial terminal lead crossing a ring emits a marker (layer = ring name, rotation = azimuth);
// the ring's turns must clear the corridor: angular distance >= marker half-angle + turn half-angle.
static int toroidal_corridor_intrusions(OpenMagnetics::Coil& coil) {
    if (!coil.get_turns_description()) {
        return -1;
    }
    auto turns = coil.get_turns_description().value();
    auto spaces = coil.get_connection_reserved_spaces();
    auto wires = coil.get_wires();
    std::map<std::string, std::pair<double, size_t>> radiusAccumulator;
    for (auto& turn : turns) {
        if (turn.get_layer()) {
            auto& acc = radiusAccumulator[turn.get_layer().value()];
            acc.first += std::hypot(turn.get_coordinates()[0], turn.get_coordinates()[1]);
            acc.second++;
        }
    }
    auto angularDistance = [](double a, double b) {
        return std::abs(std::fmod(a - b + 540.0, 360.0) - 180.0);
    };
    int intrusions = 0;
    for (auto& space : spaces) {
        if (space.layer.empty() || !radiusAccumulator.count(space.layer)) {
            continue;
        }
        double ringRadius = radiusAccumulator[space.layer].first / double(radiusAccumulator[space.layer].second);
        double markerHalf = OpenMagnetics::wound_distance_to_angle(space.dimensions[1], ringRadius) / 2;
        for (auto& turn : turns) {
            if (!turn.get_layer() || turn.get_layer().value() != space.layer) {
                continue;
            }
            size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
            double turnHalf = OpenMagnetics::wound_distance_to_angle(wires[windingIndex].get_maximum_outer_height(), ringRadius) / 2;
            double turnAngle = std::atan2(turn.get_coordinates()[1], turn.get_coordinates()[0]) * 180.0 / std::numbers::pi;
            if (angularDistance(turnAngle, space.rotation) < markerHalf + turnHalf - 0.01) {
                intrusions++;
                std::cout << "[TORCOLL] turn " << turn.get_name() << " angle=" << turnAngle
                          << " inside corridor of lead(w=" << space.winding << " p=" << space.parallel
                          << ") at " << space.rotation << " on " << space.layer << "\n";
            }
        }
    }
    return intrusions;
}

// ABT #187: number of crossing markers (spaces that name a ring) — proves the radial leads actually
// declared their ring crossings for the corridor machinery.
static int toroidal_crossing_markers(OpenMagnetics::Coil& coil) {
    int markers = 0;
    for (auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.layer.empty()) {
            markers++;
        }
    }
    return markers;
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
        // ABT #833: a TOROID gets NO closing station, so expected_stations() -- which encodes the
        // concentric rule "turns + layers spanned" -- must NOT be used here. Alf, ABT #685
        // 2026-08-18: "make sure that the N_layer + 1 circles only affect concentric cores, not
        // toroidals"; a toroid with 9 turns has 9 inner crossings and 8 outer ones, the first and
        // last inner crossing carrying the terminals. The extra concentric station exists because a
        // concentric layer's wire must come back to its starting azimuth to close the layer; a
        // toroid's turn closes through the bore and needs no such station. So the count is exactly
        // turns x parallels. (The previous expectation, and the stale "(40+1)*2" comment it
        // replaced, both predate that ruling and asserted phantom stations the winder is right not
        // to place.)
        const int64_t toroidalParallels =
            coil.get_functional_description()[coil.get_winding_index_by_name("winding 0")].get_number_parallels();
        CHECK(coil.get_turns_description().value().size() == size_t(40 * toroidalParallels));
        CHECK(toroidal_turn_overlaps(coil) == 0);
        CHECK(!coil.get_connection_reserved_spaces().empty());
        CHECK(distinct_parallels_with_terminal_leads(coil, "winding 0") == 2);
        // ABT #187: the exit leads (innermost ring) cross the outer ring(s) radially — those
        // crossings must be declared as markers, and no turn may sit inside a lead's corridor.
        CHECK(toroidal_crossing_markers(coil) > 0);
        CHECK(toroidal_corridor_intrusions(coil) == 0);
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
        CHECK(toroidal_corridor_intrusions(coil) == 0);
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

TEST_CASE("Test_Single_Layer_Winding_Emits_Terminal_Leads",
          "[constructive-model][coil][real-geometry]") {
    // A single-layer winding has no inter-layer links, but its entrance/exit TERMINAL
    // leads exist all the same (drawn + connection loss). Regression for the early
    // return that skipped ALL reserved spaces when fewer than two conduction layers.
    std::vector<int64_t> numberTurns = {8};
    std::vector<int64_t> numberParallels = {1};
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels,
                                                     "PQ 28/20", 1);
    REQUIRE(coil.get_layers_description_conduction().size() == 1);

    auto spaces = coil.get_connection_reserved_spaces();
    size_t drawnTerminals = 0;
    size_t links = 0;
    for (const auto& s : spaces) {
        if (!s.layer.empty()) continue;
        if (s.isTerminal) ++drawnTerminals; else ++links;
    }
    // One entrance + one exit lead (each may be a single own-level rect or a stub+run L).
    REQUIRE(drawnTerminals >= 2);
    REQUIRE(links == 0);
}

TEST_CASE("Test_Centered_Single_Turn_Toroidal_Emits_Outer_Crossing",
          "[constructive-model][coil][toroid][single-turn]") {
    // Toroid whose single turn's wire OD exceeds the winding-window radial height ->
    // wind() takes the build_centered_single_turn_toroidal() special path (skips the
    // sections/layers fit pipeline). That path must still emit the outer XY-plane
    // crossing (additionalCoordinates), with the same polar-mirror convention as
    // wind_toroidal_additional_turns; without it downstream consumers (Painter, 3D
    // builders) cannot know where the wire wraps the ring.
    std::string coilString = R"({"bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.0,"columnWidth":0.002625,"coordinates":[0.0,0.0,0.0],"wallThickness":0.0,"windingWindows":[{"angle":360.0,"coordinates":[0.0074,0.0,0.0],"radialHeight":0.0074,"sectionsOrientation":"overlapping","shape":"round"}]}},"functionalDescription":[{"isolationSide":"primary","name":"primary","numberParallels":1,"numberTurns":1,"wire":{"coating":{"grade":1,"type":"enamelled"},"conductingDiameter":{"nominal":0.0095},"material":"copper","name":"Round 9.50 - Custom","numberConductors":1,"outerDiameter":{"nominal":0.010},"type":"round"}}]})";

    CoilWindingConfig config;
    config.coilJsonStr = coilString;
    config.pattern = {0};
    config.repetitions = 1;
    auto coil = prepare_and_wind_coil(config);

    auto turnsOpt = coil.get_turns_description();
    REQUIRE(turnsOpt.has_value());
    REQUIRE(turnsOpt->size() == 1);
    const auto& turn = (*turnsOpt)[0];

    // Inner crossing: geometric centre of the hole, converted to cartesian.
    REQUIRE(turn.get_coordinate_system() == CoordinateSystem::CARTESIAN);
    REQUIRE_THAT(turn.get_coordinates()[0], Catch::Matchers::WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(turn.get_coordinates()[1], Catch::Matchers::WithinAbs(0.0, 1e-9));

    // Outer crossing: polar mirror {-2*columnWidth - radialHeight, 0 deg} ->
    // cartesian radius 2*radialHeight + 2*columnWidth at angle 0.
    auto addOpt = turn.get_additional_coordinates();
    REQUIRE(addOpt.has_value());
    REQUIRE(addOpt->size() == 1);
    const auto& outer = (*addOpt)[0];
    REQUIRE(outer.size() >= 2);
    const double expectedRadius = 2 * 0.0074 + 2 * 0.002625;
    REQUIRE_THAT(outer[0], Catch::Matchers::WithinAbs(expectedRadius, 1e-9));
    REQUIRE_THAT(outer[1], Catch::Matchers::WithinAbs(0.0, 1e-9));
}



// ABT #278: the mid-loop measurement gap in the rectangular blocking loop (wind() skipping
// wind_by_turns when the grown layout transiently stopped fitting) silently produced a turnless
// coil for a fitting design. Guard the contract: a design that winds ideally must also wind with
// real geometry when its blocking fixpoint fits, and real winding adds one crossing per conductor.
TEST_CASE("Test_Real_Geometry_Wind_Survives_Transient_Unfit", "[constructive-model][coil][real-geometry]") {
    namespace fs = std::filesystem;
    auto file = fs::path{std::source_location::current().file_name()}.parent_path().append("..").append("MAS").append("examples").append("13_current_sense_er95_n87.json");
    settings.reset();
    std::ifstream f(file);
    std::string data((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    json masJson = json::parse(data);
    OpenMagnetics::compat::migrate_pre_1_0(masJson);
    auto magneticIn = OpenMagnetics::Magnetic(masJson["magnetic"]);
    settings.set_coil_use_real_winding_geometry(true);
    auto magnetic = OpenMagnetics::magnetic_autocomplete(magneticIn);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    // 1-turn primary + 80-turn secondary + one real-winding crossing per conductor.
    // The secondary was 100 turns until ABT #673: 100 x 0.1125 mm OD never fit the ER 9.5/5
    // window (17 per layer where 16 go), so the example was corrected to the 80 turns it can
    // actually hold. The count here follows the design, not the other way round.
    // ABT #685: one crossing per LAYER entered, so the count follows the layout's layers.
    CHECK(coil.get_turns_description().value().size() ==
          expected_stations(coil, coil.get_functional_description()[0].get_name(), 1) +
          expected_stations(coil, coil.get_functional_description()[1].get_name(), 80));
    settings.reset();
}

}  // namespace


// ABT #230: a toroidal terminal lead's reserved rect ran from the connecting turn's crossing
// CENTRELINE radially out to maxTurnRadius + 1.5*wireOuterWidth. For a wall-adjacent outer ring
// that far end lands INSIDE the core annulus (measured 12.96 mm against a 12.0 mm bore on
// T 40/24/16 with ~0.959 mm OD wire), and MVB++ replays these rects verbatim as 3D lead routes,
// so the overrun would place copper inside the core. The near edge was also the centreline, not
// the wire envelope, so it disagreed with the concentric convention which spans
// [turnX - w/2, borderX + w/2]. Guard both ends.
TEST_CASE("Test_Toroidal_Terminal_Lead_Rect_Stays_Inside_Bore",
          "[constructive-model][coil][toroid][real-geometry]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);

    std::vector<int64_t> numberTurns = {12};
    std::vector<int64_t> numberParallels = {1};
    std::vector<OpenMagnetics::Wire> wires;
    wires.push_back(OpenMagnetics::find_wire_by_name("Round 0.90 - Grade 1"));
    auto coil = OpenMagnetics::Coil::create_quick_coil("T 40/24/16", numberTurns, numberParallels, wires,
                                                       WindingOrientation::OVERLAPPING,
                                                       WindingOrientation::OVERLAPPING,
                                                       CoilAlignment::CENTERED, CoilAlignment::CENTERED);

    auto windingWindows = coil.resolve_bobbin().get_processed_description().value().get_winding_windows();
    REQUIRE(!windingWindows.empty());
    REQUIRE(windingWindows[0].get_radial_height().has_value());
    double boreRadius = windingWindows[0].get_radial_height().value();
    double wireOuterWidth = coil.get_wires()[0].get_maximum_outer_width();

    auto spaces = coil.get_connection_reserved_spaces();
    size_t terminalLeadsChecked = 0;
    for (const auto& space : spaces) {
        // The drawn radial lead runs (not the per-ring squeeze markers, which carry a layer).
        if (!space.isTerminal || !space.layer.empty()) {
            continue;
        }
        REQUIRE(space.coordinates.size() >= 2);
        REQUIRE(space.dimensions.size() >= 2);
        double centreRadius = std::hypot(space.coordinates[0], space.coordinates[1]);
        double farEdge = centreRadius + space.dimensions[0] / 2;
        // (1) the rect must not reach past the bore wall into the core annulus.
        CHECK(farEdge <= boreRadius + 1e-9);
        // (2) the near edge is the wire envelope: it must sit at least half a wire inward of the
        // outermost turn it can start from, i.e. the rect is never a bare centreline-to-border span.
        double nearEdge = centreRadius - space.dimensions[0] / 2;
        CHECK(nearEdge < boreRadius - wireOuterWidth / 2 + 1e-9);
        terminalLeadsChecked++;
    }
    REQUIRE(terminalLeadsChecked > 0);
    settings.reset();
}


// ABT #240: a terminal lead crossing a FOREIGN winding's layer used to end up exactly tangent to
// that winding's extreme turn (measured separation 7.6e-13 um on this fixture) — the reserved band
// was one wire deep and blocking freed exactly one slot, so the next turn began where the lead
// ended. Same-winding turns touching is MKF's packing convention and stays legal; two DIFFERENT
// windings touching is not, and the separation must be the mechanical insulation the coil already
// builds between them.
TEST_CASE("Test_Terminal_Lead_Clears_Foreign_Winding_By_Insulation",
          "[constructive-model][coil][real-geometry]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    namespace fs = std::filesystem;
    auto path = fs::path{std::source_location::current().file_name()}.parent_path()
                    .append("..").append("MAS").append("examples")
                    .append("16_coupled_inductor_e2513_dmr95.json");
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto coil = magnetic.get_coil();
    REQUIRE(coil.get_turns_description().has_value());

    auto turns = coil.get_turns_description().value();
    auto wires = coil.get_wires();
    auto spaces = coil.get_connection_reserved_spaces();

    // The insulation the coil placed between the two windings — the clearance contract.
    double insulationThickness = 0;
    for (const auto& insulationSection : coil.get_sections_by_type(ElectricalType::INSULATION)) {
        insulationThickness = std::max(insulationThickness,
                                       coil.get_insulation_section_thickness(insulationSection.get_name()));
    }
    REQUIRE(insulationThickness > 0);

    size_t pairsChecked = 0;
    for (const auto& space : spaces) {
        if (!space.isTerminal || !space.layer.empty()) {
            continue;  // drawn radial runs only
        }
        for (const auto& turn : turns) {
            if (turn.get_winding() == space.winding) {
                continue;  // same net: flush packing is legal
            }
            size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
            double turnWidth = wires[windingIndex].get_maximum_outer_width();
            double turnHeight = wires[windingIndex].get_maximum_outer_height();
            double gapX = std::abs(turn.get_coordinates()[0] - space.coordinates[0])
                          - (space.dimensions[0] + turnWidth) / 2;
            double gapY = std::abs(turn.get_coordinates()[1] - space.coordinates[1])
                          - (space.dimensions[1] + turnHeight) / 2;
            // Clear on at least one axis, by at least the inter-winding insulation.
            CHECK(std::max(gapX, gapY) >= insulationThickness - 1e-12);
            pairsChecked++;
        }
    }
    REQUIRE(pairsChecked > 0);
    settings.reset();
}


// ABT #231: the companion to Test_Additiona_Turns_Bug, which only exercises the case where the
// outer face has SPARE room (T 20/10/7, 60t of 0.509 mm — 62.8 mm of outer circumference for
// 30.5 mm of wire). This fixture fills it: T 40/24/16 has 125.7 mm of outer circumference and
// 65 turns of 2.00 mm OD need 130 mm, so the first outer ring physically cannot hold them all
// and some crossings MUST stack outward. Both branches of the placement rule are then covered.
TEST_CASE("Test_Toroidal_Outer_Crossings_Stack_When_Outer_Face_Is_Full",
          "[constructive-model][coil][toroid][round-winding-window]") {
    clear_databases();
    settings.reset();
    settings.set_use_toroidal_cores(true);
    settings.set_coil_include_additional_coordinates(true);

    std::vector<int64_t> numberTurns = {65};
    std::vector<int64_t> numberParallels = {1};
    std::vector<OpenMagnetics::Wire> wires;
    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_diameter(0.0018);
    wire.set_nominal_value_outer_diameter(0.002);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::ROUND);
    wires.push_back(wire);

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "T 40/24/16", 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::SPREAD, CoilAlignment::SPREAD, wires);

    auto turns = coil.get_turns_description().value();
    REQUIRE(turns.size() == 65);
    double wireOuterDiameter = coil.get_wires()[0].get_maximum_outer_height();

    std::vector<std::vector<double>> outerCrossings;
    double minOuterRadiusFull = std::numeric_limits<double>::max();
    for (const auto& turn : turns) {
        REQUIRE(turn.get_additional_coordinates());
        auto addCoords = turn.get_additional_coordinates().value()[0];
        outerCrossings.push_back({addCoords[0], addCoords[1]});
        minOuterRadiusFull = std::min(minOuterRadiusFull, hypot(addCoords[0], addCoords[1]));
    }
    // Within TWO wire ODs of arc of its own azimuth (Alf, 2026-08-24: a crossing whose own
    // ring slot is taken goes to the NEAREST FREE RING SLOT ON EITHER SIDE of the occupant --
    // own slot + the occupant's possible one-OD lean = two ODs). Beyond that the crossing
    // STACKS at its own azimuth instead of travelling the rim, which the stacking assertions
    // below enforce on this over-subscribed face.
    const double leanLimitFull = 2 * wireOuterDiameter / minOuterRadiusFull + 1e-6;
    for (const auto& turn : turns) {
        auto addCoords = turn.get_additional_coordinates().value()[0];
        double innerAngle = atan2(turn.get_coordinates()[1], turn.get_coordinates()[0]);
        double outerAngle = atan2(addCoords[1], addCoords[0]);
        double lean = std::remainder(outerAngle - innerAngle, 2 * std::numbers::pi);
        CHECK(std::abs(lean) <= leanLimitFull);
    }

    // No two outer crossings closer than one wire OD, even though the face is over-subscribed.
    for (size_t i = 0; i < outerCrossings.size(); ++i) {
        for (size_t j = i + 1; j < outerCrossings.size(); ++j) {
            double separation = hypot(outerCrossings[i][0] - outerCrossings[j][0],
                                      outerCrossings[i][1] - outerCrossings[j][1]);
            CHECK(separation >= wireOuterDiameter - 1e-9);
        }
    }

    // The face is genuinely over-subscribed, so more than one outer ring must be in use.
    std::set<double> uniqueRadii;
    for (const auto& crossing : outerCrossings) {
        uniqueRadii.insert(round(hypot(crossing[0], crossing[1]) * 10000) / 10000);
    }
    CHECK(uniqueRadii.size() > 1);
    settings.reset();
}


// ABT #352 follow-up: a windingOrder set on the CORE's winding window (the MAS core schema
// carries it) used to be silently dropped twice — Core::process_data rebuilt the windows from
// the piece geometry, and Bobbin::create_quick_bobbin did not copy it onto the autocompleted
// bobbin — so the coil always fell back to the Z (dragback) default. Guard the whole chain on
// the WE-TI drum reconstruction: the order must survive re-processing and reach the wound coil.
TEST_CASE("Test_Core_Window_Winding_Order_Reaches_Autocompleted_Bobbin",
          "[constructive-model][coil][drum]") {
    settings.reset();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(),
                                                         "we_ti_7447720470_reconstructed.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    auto masJson = json::parse(file);

    OpenMagnetics::Core core(masJson["magnetic"]["core"]);
    core.process_data();

    // Simulate a MAS file whose CORE window carries the order (the fixture keeps its own copy
    // on the bobbin; here the coil gets a quick bobbin instead, so the core is the only source).
    auto processedDescription = core.get_processed_description().value();
    processedDescription.get_mutable_winding_windows()[0].set_winding_order(MAS::WindingOrder::U);
    core.set_processed_description(processedDescription);

    // (1) The order survives a re-run of process_data, which regenerates the winding windows.
    core.process_data();
    auto reprocessedWindows = core.get_processed_description().value().get_winding_windows();
    REQUIRE(!reprocessedWindows.empty());
    REQUIRE(reprocessedWindows[0].get_winding_order().has_value());
    CHECK(reprocessedWindows[0].get_winding_order().value() == MAS::WindingOrder::U);

    // (2) create_quick_bobbin carries it onto the autocompleted bobbin, where
    //     Coil::get_winding_order picks it up for every section of the wind.
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    auto coilJson = masJson["magnetic"]["coil"];
    coilJson["bobbin"] = "Dummy";
    magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
    auto completed = OpenMagnetics::magnetic_autocomplete(magnetic);

    auto bobbin = completed.get_mutable_coil().resolve_bobbin();
    auto bobbinWindows = bobbin.get_processed_description().value().get_winding_windows();
    REQUIRE(!bobbinWindows.empty());
    REQUIRE(bobbinWindows[0].get_winding_order().has_value());
    CHECK(bobbinWindows[0].get_winding_order().value() == MAS::WindingOrder::U);

    auto sectionsDescription = completed.get_coil().get_sections_description();
    REQUIRE(sectionsDescription.has_value());
    auto sections = sectionsDescription.value();
    REQUIRE(!sections.empty());
    CHECK(completed.get_coil().get_winding_order(sections[0].get_name()) == MAS::WindingOrder::U);
    settings.reset();
}


// ABT #357: molded body — REAL-WIRE turn placement inside the coil cavity. The cavity is the
// winding window (pot-core letters: F post OD, E cavity OD, D cavity height), so the standard
// concentric winder must land every turn's centre inside the annulus with at least a wire
// radius of clearance to the post, the cavity wall, and both plates. Guards the reconstruction
// contract: cavity dims + N x wire must be mutually consistent for a faithful molded record.
// ABT #930: a wind that produces no turns must say WHY. Reported from the field as "Turns not
// created" on a 2.5 x 2.0 x 1.0 mm moulded body (183 WE-PMI records blocked), diagnosed as a
// broken winder. It was not: the reporter's F = 0.50 / E = 0.70 leaves a (E-F)/2 = 0.10 mm radial
// annulus, and the printed conductor is 0.36 mm wide, so MKF was correctly refusing an impossible
// geometry — it just never said so, and wind()'s bool is discarded by magnetic_autocomplete. Pin
// both halves: the refusal still happens, and it now names the conductor, the window and the gap.
TEST_CASE("Test_Molded_Small_Body_Refusal_Names_The_Reason", "[constructive-model][coil][molded][abt930]") {
    settings.reset();
    auto build = [](double postF, double cavityE, int64_t turns, double wireWidth) {
        json shapeJson = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
            {"aliases", json::array()}, {"name", "PMI-like 2520"},
            {"dimensions", {
                {"A", {{"nominal", 0.0025}}}, {"B", {{"nominal", 0.0010}}}, {"C", {{"nominal", 0.0020}}},
                {"D", {{"nominal", 0.0007}}}, {"E", {{"nominal", cavityE}}}, {"F", {{"nominal", postF}}}}}
        };
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", "closedShape"}, {"material", "Kool Mµ 26"}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        OpenMagnetics::Core core(coreJson);
        core.process_data();
        core.process_gap();
        json wireJson = {
            {"type", "rectangular"}, {"material", "copper"}, {"name", "printed"},
            {"conductingWidth", {{"nominal", wireWidth}}}, {"conductingHeight", {{"nominal", 0.000025}}},
            {"outerWidth", {{"nominal", wireWidth}}}, {"outerHeight", {{"nominal", 0.000025}}},
            {"numberConductors", 1}};
        json coilJson;
        coilJson["bobbin"] = "Dummy";
        coilJson["functionalDescription"] = json::array();
        json winding;
        winding["name"] = "winding 0";
        winding["numberTurns"] = turns;
        winding["numberParallels"] = 1;
        winding["isolationSide"] = "primary";
        winding["wire"] = wireJson;
        coilJson["functionalDescription"].push_back(winding);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
        return OpenMagnetics::magnetic_autocomplete(magnetic).get_coil();
    };

    SECTION("a conductor wider than the annulus is refused, and the refusal names the numbers") {
        // 0.36 mm conductor into (0.70 - 0.50)/2 = 0.10 mm of annulus.
        auto coil = build(0.00050, 0.00070, 5, 0.00036);
        REQUIRE_FALSE(coil.get_turns_description().has_value());
        auto reason = coil.get_last_fit_failure();
        UNSCOPED_INFO("reason: " << reason);
        REQUIRE_FALSE(reason.empty());
        CHECK(reason.find("winding 0") != std::string::npos);
        CHECK(reason.find("0.360") != std::string::npos);   // the conductor
        CHECK(reason.find("0.100") != std::string::npos);   // the window it will not fit
    }

    SECTION("turn count is not the lever: one turn fails identically") {
        auto coil = build(0.00050, 0.00070, 1, 0.00036);
        REQUIRE_FALSE(coil.get_turns_description().has_value());
        CHECK_FALSE(coil.get_last_fit_failure().empty());
    }

    SECTION("widen the annulus and the same part winds") {
        // (0.80 - 0.30)/2 = 0.25 mm of annulus for a 0.10 mm conductor.
        auto coil = build(0.00030, 0.00080, 5, 0.00010);
        REQUIRE(coil.get_turns_description().has_value());
        CHECK(coil.get_turns_description()->size() == 5);
        CHECK(coil.get_last_fit_failure().empty());
    }
    settings.reset();
}

TEST_CASE("Test_Molded_Cavity_Turn_Placement", "[constructive-model][coil][molded]") {
    settings.reset();
    json shapeJson = {
        {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
        {"aliases", json::array()}, {"name", "MAPI-like 4020"},
        {"dimensions", {
            {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
            {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}
    };
    json coreJson;
    coreJson["functionalDescription"] = {
        {"type", "closedShape"}, {"material", "Kool Mµ 26"}, {"shape", shapeJson},
        {"gapping", json::array()}, {"numberStacks", 1}};
    OpenMagnetics::Core core(coreJson);
    core.process_data();
    core.process_gap();

    json coilJson;
    coilJson["bobbin"] = "Dummy";
    coilJson["functionalDescription"] = json::array();
    json winding;
    winding["name"] = "winding 0";
    winding["numberTurns"] = 8;
    winding["numberParallels"] = 1;
    winding["isolationSide"] = "primary";
    winding["wire"] = "Round 0.1 - Grade 1";
    coilJson["functionalDescription"].push_back(winding);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
    auto completed = OpenMagnetics::magnetic_autocomplete(magnetic);
    auto coil = completed.get_coil();
    REQUIRE(coil.get_turns_description().has_value());
    auto turns = coil.get_turns_description().value();
    REQUIRE(turns.size() == 8);

    double wireOuterRadius = coil.get_wires()[0].get_maximum_outer_width() / 2;
    double postRadius = 0.0012 / 2;
    double cavityRadius = 0.0030 / 2;
    double cavityHalfHeight = 0.0014 / 2;
    for (auto& turn : turns) {
        double radialPosition = turn.get_coordinates()[0];
        double heightPosition = turn.get_coordinates()[1];
        CHECK(radialPosition >= postRadius + wireOuterRadius * 0.99);
        CHECK(radialPosition <= cavityRadius - wireOuterRadius * 0.99);
        CHECK(std::abs(heightPosition) <= cavityHalfHeight - wireOuterRadius * 0.99);
    }
    settings.reset();
}

// ABT #374: a bore-capacity claim of the form (N_turns + N_leads) * wireOD <= pi * ID says these
// toroidal fixtures cannot hold their own leads. That inequality is the capacity of a SINGLE ring
// at the inner surface, and it is not what MKF winds: the toroidal winder fills concentric rings
// inward, so a winding "1.59x over" on one ring is comfortable on three. This pins what MKF
// actually produces for the densest of those fixtures, so the distinction stays visible: every
// ring must be non-overlapping AND leave azimuthal slack for a lead to pass.
TEST_CASE("Test_Toroidal_Dense_Winding_Fills_Concentric_Rings_With_Slack", "[coil][toroidal]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    // T 40/24/16, 60 turns of 2 mm wire: 124 mm of wire against a 75.4 mm bore circumference.
    auto coil = OpenMagneticsTesting::get_quick_coil({60}, {1}, "T 40/24/16", 1,
                                                     MAS::WindingOrientation::OVERLAPPING,
                                                     MAS::WindingOrientation::OVERLAPPING,
                                                     MAS::CoilAlignment::CENTERED,
                                                     MAS::CoilAlignment::CENTERED,
                                                     {OpenMagnetics::find_wire_by_name("Round 2.00 - Grade 1")});
    REQUIRE(coil.get_turns_description());
    auto turns = coil.get_turns_description().value();
    double wireOuterDiameter = OpenMagnetics::find_wire_by_name("Round 2.00 - Grade 1").get_maximum_outer_width();

    std::map<int, std::vector<std::vector<double>>> turnsPerRing;
    for (auto turn : turns) {
        auto coordinates = turn.get_coordinates();
        turnsPerRing[static_cast<int>(std::round(hypot(coordinates[0], coordinates[1]) * 1e4))].push_back(coordinates);
    }
    UNSCOPED_INFO(turns.size() << " crossings spread over " << turnsPerRing.size() << " rings");
    CHECK(turnsPerRing.size() > 1);  // a single ring is exactly the geometry that cannot fit

    for (auto& [ringKey, ringTurns] : turnsPerRing) {
        double ringRadius = ringKey / 1e4;
        double closestApproach = std::numeric_limits<double>::max();
        for (size_t i = 0; i < ringTurns.size(); ++i) {
            for (size_t j = i + 1; j < ringTurns.size(); ++j) {
                closestApproach = std::min(closestApproach,
                    hypot(ringTurns[i][0] - ringTurns[j][0], ringTurns[i][1] - ringTurns[j][1]));
            }
        }
        double ringCircumference = 2 * std::numbers::pi * ringRadius;
        double occupied = ringTurns.size() * wireOuterDiameter;
        UNSCOPED_INFO("ring r=" << ringRadius * 1000 << " mm holds " << ringTurns.size() << " turns, "
                      << occupied * 1000 << " of " << ringCircumference * 1000 << " mm, closest approach "
                      << closestApproach * 1000 << " mm vs OD " << wireOuterDiameter * 1000 << " mm");
        // No turn may sit closer to its neighbour than one wire diameter...
        CHECK(closestApproach >= wireOuterDiameter - 1e-9);
        // ...and no ring may be packed solid, or a lead would have no azimuth to pass through.
        CHECK(occupied < ringCircumference);
    }
    settings.reset();
}

// Pins the two SPREAD placement rules the turn placer must keep (ABT #578 / ABT #579):
//   #579 fence-post -- the outermost turns' SURFACES sit on the layer edges, so a SPREAD layer
//        leaves no dead margin at the flanges (it used to park half a gap outside each end).
//   #578 bundles    -- when parallels are wound together N-filar (WIND_BY_CONSECUTIVE_PARALLELS)
//        the members of one bundle TOUCH, and only the gaps BETWEEN bundles carry the slack.
// Sweeps 1..4 parallels so both the degenerate (bundle == 1 turn) and the N-filar cases are
// covered, and requires that the N-filar configurations really did produce bundles, so the
// bundling assertions cannot go quietly vacuous if the winding-style heuristic changes.
TEST_CASE("Test_Spread_Is_Fence_Post_And_Keeps_Bundles_Together", "[constructive-model][coil][spread]") {
    settings.reset();
    size_t layersWithBundlesChecked = 0;
    for (int64_t parallels : {1, 2, 3, 4}) {
        auto coil = OpenMagneticsTesting::get_quick_coil(std::vector<int64_t>{12},
                                                         std::vector<int64_t>{parallels},
                                                         "PQ 28/20", 1,
                                                         WindingOrientation::OVERLAPPING,
                                                         WindingOrientation::OVERLAPPING,
                                                         CoilAlignment::SPREAD,
                                                         CoilAlignment::CENTERED);
        REQUIRE(coil.get_turns_description());
        auto turns = coil.get_turns_description().value();
        auto layers = coil.get_layers_description().value();
        for (auto& layer : layers) {
            if (layer.get_type() != ElectricalType::CONDUCTION) {
                continue;
            }
            std::vector<double> centres;
            for (auto& turn : turns) {
                if (turn.get_layer() && turn.get_layer().value() == layer.get_name()) {
                    centres.push_back(turn.get_coordinates()[1]);
                }
            }
            if (centres.size() < 2) {
                continue;
            }
            std::sort(centres.begin(), centres.end());
            auto windingIndex = coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
            double wireOuterHeight = coil.get_wires()[windingIndex].get_maximum_outer_height();
            double layerLow = layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2;
            double layerHigh = layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2;
            INFO("parallels=" << parallels << " layer=" << layer.get_name() << " turns=" << centres.size());

            // ABT #579: outermost turn surfaces flush with the layer edges (1 nm coordinate rounding).
            CHECK_THAT(centres.front() - wireOuterHeight / 2, Catch::Matchers::WithinAbs(layerLow, 1e-9));
            CHECK_THAT(centres.back() + wireOuterHeight / 2, Catch::Matchers::WithinAbs(layerHigh, 1e-9));

            // ABT #578: inside a bundle the turns touch; between bundles every gap is the same.
            int64_t bundleSize = 1;
            if (layer.get_winding_style() &&
                layer.get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS) {
                for (auto proportion : layer.get_partial_windings()[0].get_parallels_proportion()) {
                    if (proportion > 0) {
                        bundleSize++;
                    }
                }
                bundleSize--;
            }
            if (bundleSize > 1) {
                layersWithBundlesChecked++;
            }
            double interBundleGap = -1;
            for (size_t k = 1; k < centres.size(); ++k) {
                double gap = centres[k] - centres[k - 1] - wireOuterHeight;
                if (int64_t(k) % bundleSize != 0) {
                    CHECK_THAT(gap, Catch::Matchers::WithinAbs(0.0, 1e-9));  // inside a bundle: touching
                }
                else {
                    if (interBundleGap < 0) {
                        interBundleGap = gap;
                    }
                    CHECK(gap >= -1e-9);
                    // Turn centres are rounded to 1 nm (roundFloat(..., 9)), so a gap measured
                    // between two of them can differ from its neighbour by up to ~2 nm purely from
                    // that grid. Compare at the quantization, not below it.
                    CHECK_THAT(gap, Catch::Matchers::WithinAbs(interBundleGap, 5e-9));
                }
            }
        }
    }
    CHECK(layersWithBundlesChecked > 0);
    settings.reset();
}

// The five ABT #577 examples whose drawn connection routes used to interpenetrate winding copper
// must keep POSITIVE bare-copper clearance under the SPREAD layout (ABT #578/#579 + continuous
// blocked depths). Replicates the MVB++ gate's 2D test: distance from every drawn TERMINAL
// segment's copper CENTERLINE to every turn centre, minus the sum of bare radii. Z-order
// continuations are exempt (the classic in-plane crossover, same convention as
// compute_connection_blocked_slots_per_layer), and junction contact — a segment ending ON the
// turn it serves — is excluded by endpoint proximity. The authoritative gate is MVB++'s 3D
// realization; this pins the 2D geometry that gate measures, so a regression cannot ship unseen.
// Baselines on the fixed layout (2026-08-07): 06 +0.034, 13 +0.0125, 14 +0.055, 23 +0.034,
// 24 +0.057 mm — and the 38.6 mm full-window exit vertical of 14_dab is gone.
// ABT #676: margin tape is reserved space in the GEOMETRY, not just in the drawing. Real-winding
// blocking re-aligns the turns of every blocked layer, and it used to spread them against the raw
// winding window — so the copper (and the section and terminal leads that follow it) walked into
// the margin band while the painter kept drawing the margin. Checked with the flag off and on:
// the ideal wind always respected the margin, only the blocked one did not.
// ABT #676, second half: the margin has to survive a RE-wind. It is persisted on the sections,
// but the winder only ever read it from the transient _marginsPerSection — so every consumer that
// re-winds an already-wound coil (MVB++'s internal autocomplete behind the 3D view and the STEP
// export, PyOpenMagnetics, anything round-tripping a MAS) wound the copper straight over the tape.
// That is why the 3D ignored a margin the 2D drew correctly.
// ABT #682: U order and margin tape each worked alone and failed together. With a margin the
// window band shrinks, the U landing reservation at the top no longer leaves room, and the
// fence-post spread centred the copper — eating into the reservation at BOTH ends, which put the
// layer's last turn exactly on the entrance lead's row at the bottom. The 2D drew it; the 3D
// conductor builder refused it as a collision, so the user got a core and bobbin with no copper.
// Every conduction turn must stay clear of the lead row its own layer reserved.
// ABT #683 (Alf): "in U windings the second layer starts just after the first layer, that is the
// point of U winding, so we just need to connect to it horizontally... when using real winding all
// layers MUST use SPREAD turns, so that the first turn of the second layer is at the same height as
// the last turn of the previous layer". The U landing used to reserve one wire OD PAST the arrival,
// so the landing layer hung a full turn below the turn it connects to and the link was a diagonal
// drop rather than a horizontal step — and the reservation was measured from the raw window, which
// counted the section's margin a second time on top of that.
TEST_CASE("Test_Abt683_U_Layers_Connect_Level", "[constructive-model][coil][real-geometry][abt683]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_allow_margin_tape(true);
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt683_u_order_margin_e16.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto& coil = magnetic.get_mutable_coil();
    coil.wind();

    REQUIRE(coil.get_layers_description());
    REQUIRE(coil.get_turns_description());
    auto layers = coil.get_layers_description().value();
    auto turns = coil.get_turns_description().value();

    // Layers in ELECTRICAL order, each with its turns in wound order.
    std::vector<std::vector<const Turn*>> turnsPerLayer;
    std::vector<std::string> layerOrder;
    for (const auto& turn : turns) {
        if (!turn.get_layer()) {
            continue;
        }
        auto found = std::find(layerOrder.begin(), layerOrder.end(), turn.get_layer().value());
        if (found == layerOrder.end()) {
            layerOrder.push_back(turn.get_layer().value());
            turnsPerLayer.push_back({});
            found = layerOrder.end() - 1;
        }
        turnsPerLayer[size_t(found - layerOrder.begin())].push_back(&turn);
    }
    REQUIRE(turnsPerLayer.size() >= 2);
    REQUIRE(coil.get_winding_order(layers[0].get_section().value()) == WindingOrder::U);

    double wireOuterHeight = coil.get_wires()[0].get_maximum_outer_height();
    for (size_t layerIndex = 0; layerIndex + 1 < turnsPerLayer.size(); ++layerIndex) {
        double arrival = turnsPerLayer[layerIndex].back()->get_coordinates()[1];
        double landing = turnsPerLayer[layerIndex + 1].front()->get_coordinates()[1];
        INFO("layer " << layerIndex << " leaves at " << arrival * 1000
             << " mm and layer " << layerIndex + 1 << " starts at " << landing * 1000 << " mm");
        // Horizontal step: the two turns the link joins sit at the SAME height, so the connection
        // is a radial hop between neighbours. A tolerance well under one wire keeps this honest
        // without pinning the exact spread arithmetic.
        CHECK(std::abs(landing - arrival) < wireOuterHeight / 4);
    }

    // ...and that only holds because every layer is SPREAD: a packed layer would end short of the
    // band edge and the next one could not meet it.
    for (const auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        INFO(layer.get_name());
        CHECK(layer.get_turns_alignment() == CoilAlignment::SPREAD);
    }
    settings.reset();
}

// ABT #684 (Alf, from the 2D view): "the connection is going through the margin in the bottom".
// Margin tape is reserved for TAPE — and that applies to the leads, not just the turns. The
// terminal edge rows were stacked from the WINDOW edge, so with a margin the lead ran inside the
// band the tape occupies. The turns were already correct (ABT #676), which is what made the
// picture so clear: copper stopped at the tape and the magenta lead crossed it.
TEST_CASE("Test_Abt684_Terminal_Leads_Stay_Out_Of_The_Margin", "[constructive-model][coil][real-geometry][abt684]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_allow_margin_tape(true);
    // Asymmetric on purpose: a single margin would not catch a row measured from the wrong edge.
    const double topMargin = 0.003;
    const double bottomMargin = 0.002;

    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt646_e16_litz_2layer_leadcollision.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    coil.set_functional_description(sourceCoil.get_functional_description());
    coil.preload_margins({{topMargin, bottomMargin}});
    coil.wind();

    auto bobbinOut = coil.resolve_bobbin();
    auto processedDescription = bobbinOut.get_processed_description().value();
    auto windingWindow = processedDescription.get_winding_windows()[0];
    double windowTop = windingWindow.get_coordinates().value()[1] + windingWindow.get_height().value() / 2;
    double windowBottom = windingWindow.get_coordinates().value()[1] - windingWindow.get_height().value() / 2;
    const double usableTop = windowTop - topMargin;
    const double usableBottom = windowBottom + bottomMargin;

    size_t terminals = 0;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.isTerminal) {
            continue;
        }
        ++terminals;
        double low = space.coordinates[1] - space.dimensions[1] / 2;
        double high = space.coordinates[1] + space.dimensions[1] / 2;
        INFO("terminal of " << space.winding << " at y [" << low << "," << high << "] against usable ["
             << usableBottom << "," << usableTop << "]");
        CHECK(high <= usableTop + 1e-9);
        CHECK(low >= usableBottom - 1e-9);
    }
    REQUIRE(terminals > 0);

    // And the turns still respect it (ABT #676) — the lead inset must not have moved the copper.
    auto turns = coil.get_turns_description().value();
    auto wires = coil.get_wires();
    for (const auto& turn : turns) {
        size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        double half = wires[windingIndex].get_maximum_outer_height() / 2;
        INFO(turn.get_name());
        CHECK(turn.get_coordinates()[1] + half <= usableTop + 1e-9);
        CHECK(turn.get_coordinates()[1] - half >= usableBottom - 1e-9);
    }
    settings.reset();
}

// ABT #726: the same #684 contract on the CONTIGUOUS (transposed) path. Contiguous layers run
// their turns along x, so margin[0] ("top or left") owns the LEFT window edge and margin[1] the
// RIGHT — align_blocked_layer_turns already insets with that convention. The reserved-space row
// allocator mapped margin[0] to the virtual-frame top unconditionally, which on this path is the
// real RIGHT: with asymmetric tape the lead rows were inset on the tape-FREE side and stacked
// straight through the tape on the other.
TEST_CASE("Test_Abt726_Contiguous_Terminal_Leads_Stay_Out_Of_The_Margin", "[constructive-model][coil][real-geometry][abt726]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_allow_margin_tape(true);
    // Asymmetric on purpose: a symmetric margin cannot catch a row measured from the wrong edge.
    // Contiguous margins live on the WIDTH axis of the window.
    const double leftMargin = 0.0006;   // margin[0]: "top or left" — LEFT for contiguous layers
    const double rightMargin = 0.0003;  // margin[1]: "bottom or right" — RIGHT

    // Same fixture as Test_Real_Geometry_Rectangular_Contiguous: 12 turns x 2 parallels on
    // PQ 28/20, sections and layers both contiguous. Rebuilt UNWOUND (the quick-coil constructor
    // winds, which sizes _marginsPerSection and would push the preloaded margins past index 0).
    std::vector<int64_t> numberTurns = {12};
    std::vector<int64_t> numberParallels = {2};
    auto sourceCoil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", 1,
        WindingOrientation::CONTIGUOUS, WindingOrientation::CONTIGUOUS);
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    coil.set_functional_description(sourceCoil.get_functional_description());
    coil.set_winding_orientation(WindingOrientation::CONTIGUOUS);
    coil.set_layers_orientation(WindingOrientation::CONTIGUOUS);
    coil.preload_margins({{leftMargin, rightMargin}});
    coil.wind();
    REQUIRE(coil.get_turns_description());

    // The contract only means something if the wind actually produced contiguous conduction layers.
    bool anyContiguousConductionLayer = false;
    const auto layersForOrientationCheck = coil.get_layers_description().value();
    for (const auto& layer : layersForOrientationCheck) {
        if (layer.get_type() == ElectricalType::CONDUCTION
            && layer.get_orientation() == WindingOrientation::CONTIGUOUS) {
            anyContiguousConductionLayer = true;
        }
    }
    REQUIRE(anyContiguousConductionLayer);

    auto bobbinOut = coil.resolve_bobbin();
    auto processedDescription = bobbinOut.get_processed_description().value();
    auto windingWindow = processedDescription.get_winding_windows()[0];
    double windowRight = windingWindow.get_coordinates().value()[0] + windingWindow.get_width().value() / 2;
    double windowLeft = windingWindow.get_coordinates().value()[0] - windingWindow.get_width().value() / 2;
    const double usableRight = windowRight - rightMargin;
    const double usableLeft = windowLeft + leftMargin;

    size_t terminals = 0;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.isTerminal) {
            continue;
        }
        ++terminals;
        double low = space.coordinates[0] - space.dimensions[0] / 2;
        double high = space.coordinates[0] + space.dimensions[0] / 2;
        INFO("terminal of " << space.winding << " at x [" << low << "," << high << "] against usable ["
             << usableLeft << "," << usableRight << "]");
        CHECK(high <= usableRight + 1e-9);
        CHECK(low >= usableLeft - 1e-9);
    }
    REQUIRE(terminals > 0);
    settings.reset();
}

// get_interleaving_level used to return the repetitions of the LAST wind (_currentRepetitions),
// so before any wind it read 1 no matter what was requested, and after set_interleaving_level
// it kept returning the stale wound value — MagneticAdviser's level-1 retry and CoilAdviser's
// WireAdviser call both consumed the wrong number. The wound value is now get_current_repetitions.
TEST_CASE("Test_Get_Interleaving_Level_Reflects_Setter", "[constructive-model][coil]") {
    OpenMagnetics::Coil coil;
    coil.set_interleaving_level(3);
    CHECK(coil.get_interleaving_level() == 3);
    CHECK(coil.get_current_repetitions() == 1);
}

// The per-section set_layers_orientation overload (exposed in PyOM and the web winding studio)
// wrote _layersOrientationPerSection, but nothing ever read the map: sections always got the
// global orientation, making the override a silent no-op.
TEST_CASE("Test_Per_Section_Layers_Orientation_Override_Is_Honored", "[constructive-model][coil]") {
    settings.reset();
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt646_e16_litz_2layer_leadcollision.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    coil.set_functional_description(sourceCoil.get_functional_description());
    coil.wind();

    std::string targetName;
    WindingOrientation baseOrientation = WindingOrientation::OVERLAPPING;
    auto sectionsBefore = coil.get_sections_description().value();
    for (const auto& section : sectionsBefore) {
        if (section.get_type() == ElectricalType::CONDUCTION) {
            targetName = section.get_name();
            baseOrientation = section.get_layers_orientation();
            break;
        }
    }
    REQUIRE(!targetName.empty());
    auto flippedOrientation = baseOrientation == WindingOrientation::OVERLAPPING ? WindingOrientation::CONTIGUOUS
                                                                                 : WindingOrientation::OVERLAPPING;
    coil.set_layers_orientation(flippedOrientation, targetName);
    coil.wind();

    bool foundTarget = false;
    auto sectionsAfter = coil.get_sections_description().value();
    for (const auto& section : sectionsAfter) {
        if (section.get_name() == targetName) {
            foundTarget = true;
            CHECK(section.get_layers_orientation() == flippedOrientation);
        }
        else if (section.get_type() == ElectricalType::CONDUCTION) {
            // Sections without an override keep the global orientation.
            CHECK(section.get_layers_orientation() == baseOrientation);
        }
    }
    REQUIRE(foundTarget);
    settings.reset();
}

// ABT #674: real winding winds one extra STATION per parallel — a wire making N turns crosses the
// window plane N+1 times, and MVB++ builds the copper as the wraps BETWEEN stations (drop the
// station and it builds one wrap fewer: 82 primitives against 87). The station is therefore not a
// turn and carries no wrap, so it must contribute no LENGTH: every consumer that sums turn lengths
// — WindingOhmicLosses computes DC resistance exactly that way — counted one turn too many, ~5% on
// a 20-turn winding.
TEST_CASE("Test_Abt674_Crossing_Station_Adds_No_Length", "[constructive-model][coil][real-geometry][abt674]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto examplesDir = std::filesystem::path{__FILE__}.parent_path().append("..").append("MAS").append("examples");
    auto mas = OpenMagneticsTesting::mas_loader((examplesDir / "17_cllc_xfmr_e5528_3c92a.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto& coil = magnetic.get_mutable_coil();

    REQUIRE(coil.get_turns_description());
    auto turns = coil.get_turns_description().value();
    for (const auto& winding : coil.get_functional_description()) {
        int64_t declared = winding.get_number_turns() * winding.get_number_parallels();
        std::vector<const Turn*> windingTurns;
        for (const auto& turn : turns) {
            if (turn.get_winding() == winding.get_name()) {
                windingTurns.push_back(&turn);
            }
        }
        INFO(winding.get_name() << ": " << windingTurns.size() << " stations for " << declared << " turns");
        // One station per LAYER entered, per parallel (ABT #685).
        CHECK(windingTurns.size() ==
              expected_stations(coil, winding.get_name(), winding.get_number_turns()));

        // ...and it carries no length, so the sum is the conductor's real length: N wraps.
        size_t zeroLength = 0;
        double total = 0;
        for (const auto* turn : windingTurns) {
            total += turn->get_length();
            if (turn->get_length() == 0) {
                ++zeroLength;
            }
        }
        // ...one per (layer, parallel): every layer's OPENING crossing carries no wrap, its
        // copper being the connector, whose length the connection markers charge explicitly.
        size_t expectedZero = 0;
        for (int64_t parallel = 0; parallel < winding.get_number_parallels(); ++parallel) {
            expectedZero += layers_of_conductor(coil, winding.get_name(), parallel);
        }
        CHECK(zeroLength == expectedZero);
        double meanWrap = total / double(declared);
        INFO("total " << total << " m over " << declared << " wraps, mean " << meanWrap);
        CHECK(total > 0);
        // The sum must be N wraps, not N+1: no station's length may be counted twice over.
        CHECK(total <= meanWrap * double(declared) + 1e-9);
    }
    settings.reset();
}

TEST_CASE("Test_Abt685_Parallels_Share_One_Terminal_Row", "[constructive-model][coil][real-geometry][abt685]") {
    // ABT #685 (Alf, on an 8t x 2p E16): "parallels are still using one OD of height block in the
    // terminal connection per parallel, when it should be one, and spread in angle" — and, after
    // the first attempt was reverted, "I insist ... they can start in the angle their terminal is".
    // A winding's parallels leave TOGETHER: ONE row of wire height for the whole bundle, the leads
    // separated azimuthally by MVB++'s fan, and the parallel whose first turn sits higher climbing
    // a vertical stub at its own angle. Stacking a row per parallel spent a whole wire of window
    // height on every extra conductor.
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt646_e16_litz_2layer_leadcollision.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    auto windings = sourceCoil.get_functional_description();
    windings[0].set_number_turns(8);
    windings[0].set_number_parallels(2);
    coil.set_functional_description(windings);
    REQUIRE(coil.wind());

    // The RUNS (radial: wider than tall) are the rows; the tall, narrow markers are the vertical
    // stubs that carry a parallel from the shared row up to its own first turn — one per parallel
    // except the one whose turn sits at the row, exactly as the design intends.
    std::map<double, std::set<int64_t>> parallelsPerRow;
    size_t verticalStubs = 0;
    for (const auto& reservedSpace : coil.get_connection_reserved_spaces()) {
        if (!reservedSpace.isTerminal || !reservedSpace.layer.empty()) {
            continue;
        }
        if (reservedSpace.dimensions[0] <= reservedSpace.dimensions[1]) {
            ++verticalStubs;
            continue;
        }
        parallelsPerRow[roundFloat(reservedSpace.coordinates[1], 6)].insert(reservedSpace.parallel);
    }
    REQUIRE(!parallelsPerRow.empty());
    // The ENTRANCE edge is the one the winding starts from, and it is the expensive one: those
    // leads run inward THROUGH the outer layers, so every extra row costs a turn slot in each of
    // them. (The exits leave from the outermost layer and cross nothing, so their runs are free
    // and stay at their own turns' heights.) On the entrance edge the bundle shares ONE row.
    const auto turns = coil.get_turns_description().value();
    REQUIRE(!turns.empty());
    const bool entranceAtTop = turns.front().get_coordinates()[1] > 0;
    size_t entranceRows = 0;
    for (const auto& [rowY, parallelsOnRow] : parallelsPerRow) {
        if ((rowY > 0) != entranceAtTop) {
            continue;
        }
        ++entranceRows;
        INFO("entrance row at y = " << rowY * 1000 << " mm holds " << parallelsOnRow.size()
                                    << " of the winding's 2 parallels");
        CHECK(parallelsOnRow.size() == 2);
    }
    CHECK(entranceRows == 1);
    INFO(verticalStubs << " vertical stubs carry parallels from their shared row to their turns");
    CHECK(verticalStubs > 0);
    settings.reset();
}

TEST_CASE("Test_Abt685_U_Parallels_Keep_Their_Lane", "[constructive-model][coil][real-geometry][abt685]") {
    // The same 8t x 2p design in U order: the parallels must keep their spatial order in every
    // layer (wires wound side by side turn around together, they do not swap places), and each
    // must land LEVEL with its OWN last turn. Landing level used to be single-conductor-only
    // exactly because the bundle flipped over at the turnaround and the links crossed.
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt646_e16_litz_2layer_leadcollision.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    auto bobbin = sourceCoil.resolve_bobbin();
    auto processedDescription = bobbin.get_processed_description().value();
    auto windingWindows = processedDescription.get_winding_windows();
    windingWindows[0].set_winding_order(MAS::WindingOrder::U);
    processedDescription.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processedDescription);

    OpenMagnetics::Coil coil;
    coil.set_bobbin(bobbin);
    auto windings = sourceCoil.get_functional_description();
    windings[0].set_number_turns(8);
    windings[0].set_number_parallels(2);
    coil.set_functional_description(windings);
    REQUIRE(coil.wind());

    const auto turns = coil.get_turns_description().value();
    std::map<std::string, std::vector<size_t>> turnsPerLayer;
    std::vector<std::string> layerOrder;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        const auto layerName = turns[turnIndex].get_layer().value();
        if (!turnsPerLayer.count(layerName)) {
            layerOrder.push_back(layerName);
        }
        turnsPerLayer[layerName].push_back(turnIndex);
    }
    REQUIRE(layerOrder.size() >= 2);

    // Lane check: in every layer, the topmost turn belongs to the same parallel.
    std::optional<int64_t> topParallel;
    for (const auto& layerName : layerOrder) {
        auto indices = turnsPerLayer[layerName];
        std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
            return turns[a].get_coordinates()[1] > turns[b].get_coordinates()[1];
        });
        if (indices.size() < 2) {
            continue;
        }
        if (!topParallel) {
            topParallel = turns[indices.front()].get_parallel();
            continue;
        }
        INFO("layer " << layerName << " is topped by parallel "
                      << turns[indices.front()].get_parallel());
        CHECK(turns[indices.front()].get_parallel() == topParallel.value());
    }

    // Landing lane and pitch, per parallel: each parallel lands EXACTLY its layer's own pitch below
    // its OWN last turn — K wire ODs for K parallels, because that is what every revolution of a
    // K-filar layer advances. Landing any less puts the landing ramp out of phase with the layer it
    // joins and the sibling, ramping the full pitch, crosses it. What this pins is that the wire
    // lands in its own lane at the layer's own gradient, not on the station its sibling departs
    // from (the flipped bundle) nor half way down it (the one-OD landing).
    const double wireHeight = coil.get_wires()[0].get_maximum_outer_height();
    for (int64_t parallelIndex = 0; parallelIndex < 2; ++parallelIndex) {
        std::optional<double> lastInFirstLayer;
        std::optional<double> firstInSecondLayer;
        for (const auto& turn : turns) {
            if (turn.get_parallel() != parallelIndex) {
                continue;
            }
            if (turn.get_layer().value() == layerOrder[0]) {
                lastInFirstLayer = turn.get_coordinates()[1];
            }
            else if (turn.get_layer().value() == layerOrder[1] && !firstInSecondLayer) {
                firstInSecondLayer = turn.get_coordinates()[1];
            }
        }
        REQUIRE(lastInFirstLayer);
        REQUIRE(firstInSecondLayer);
        INFO("parallel " << parallelIndex << " leaves layer 0 at " << lastInFirstLayer.value() * 1000
                         << " mm and lands at " << firstInSecondLayer.value() * 1000 << " mm");
        // Alf 2026-08-24 (spread-all, abt631 review): "the U connection between layers must be
        // horizontal and the layers must be on the same height" — the landing is LEVEL with the
        // departure, the link a pure radial step, and no landing band is reserved. This
        // supersedes the ABT #608/#685 K-OD descent this pin used to hold (still reachable via
        // MKF_NO_SPREAD_ALL, where the old expectation would apply). The lane guarantee itself —
        // each parallel lands on its OWN turn, not its sibling's — is untouched and pinned by
        // the parallel-identity checks above.
        (void)wireHeight;
        CHECK(std::abs(firstInSecondLayer.value() - lastInFirstLayer.value()) ==
              Catch::Approx(0.0).margin(1e-6));
    }
    settings.reset();
}

// ABT #728: the crossing station must be identified by its EXPLICIT wind-order marker — the
// turn named "turn 0" of each (winding, parallel) — not by first-seen vector order, which
// zeroes a mid-winding wrap the moment any post-wind pass reorders the description.
TEST_CASE("Test_Abt728_Crossing_Station_Is_The_Turn0_Marker", "[constructive-model][coil][real-geometry][abt728]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto examplesDir = std::filesystem::path{__FILE__}.parent_path().append("..").append("MAS").append("examples");
    auto mas = OpenMagneticsTesting::mas_loader((examplesDir / "17_cllc_xfmr_e5528_3c92a.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto& coil = magnetic.get_mutable_coil();

    REQUIRE(coil.get_turns_description());
    const auto turns = coil.get_turns_description().value();
    // ABT #685: a length-free STATION is the opening crossing of each (layer, parallel) — the
    // wire enters every layer on the plane, and the copper before that point is the connector,
    // charged on the connection markers instead.
    std::set<std::pair<std::string, int64_t>> openedLayers;
    for (const auto& turn : turns) {
        const bool isStation =
            turn.get_layer()
                ? openedLayers.insert({turn.get_layer().value(), turn.get_parallel()}).second
                : turn.get_name() == turn.get_winding() + " parallel "
                                         + std::to_string(turn.get_parallel()) + " turn 0";
        INFO(turn.get_name() << " length " << turn.get_length());
        if (isStation) {
            CHECK(turn.get_length() == 0);
        }
        else {
            CHECK(turn.get_length() > 0);
        }
    }
    settings.reset();
}

// ABT #728: the station-zeroing is gated on the turns being THIS wind's output, not on the
// wind FITTING — with windEvenIfNotFit the caller consumes the not-fitting layout and its
// stations must still carry no length (a naive success-gate would re-inflate DC resistance
// by 1/N on exactly those flows). And a wind that fails without producing turns must leave
// nothing behind for the destructor to corrupt.
TEST_CASE("Test_Abt728_Not_Fitting_Wind_Still_Zeroes_Its_Own_Stations", "[constructive-model][coil][real-geometry][abt728]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_try_rewind(false);
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt646_e16_litz_2layer_leadcollision.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    auto functionalDescription = sourceCoil.get_functional_description();
    // Blow the design far past the window so the ideal wind cannot fit.
    for (auto& winding : functionalDescription) {
        winding.set_number_turns(winding.get_number_turns() * 8);
    }
    coil.set_functional_description(functionalDescription);
    bool wound = coil.wind();
    CHECK_FALSE(wound);
    REQUIRE(coil.get_turns_description());
    const auto turns = coil.get_turns_description().value();
    size_t stations = 0;
    std::set<std::pair<std::string, int64_t>> openedLayers;
    for (const auto& turn : turns) {
        const bool isStation =
            turn.get_layer()
                ? openedLayers.insert({turn.get_layer().value(), turn.get_parallel()}).second
                : turn.get_name() == turn.get_winding() + " parallel "
                                         + std::to_string(turn.get_parallel()) + " turn 0";
        INFO(turn.get_name() << " length " << turn.get_length());
        if (isStation) {
            ++stations;
            CHECK(turn.get_length() == 0);
        }
        else {
            CHECK(turn.get_length() > 0);
        }
    }
    // ABT #685: one station per (layer, parallel) -- the opening crossing of every layer the
    // conductor enters, not one per conductor.
    size_t expectedStations = 0;
    for (const auto& winding : coil.get_functional_description()) {
        for (int64_t parallel = 0; parallel < winding.get_number_parallels(); ++parallel) {
            expectedStations += layers_of_conductor(coil, winding.get_name(), parallel);
        }
    }
    CHECK(stations == expectedStations);
    settings.reset();
}

// ABT #725: the per-section layers-orientation override is honored where the section record
// is written, but section sizing and insulation still read the global orientation — so an
// overridden section COULD get layers stacked along an axis its envelope was not sized for.
// This pins the geometric consequence, not just the recorded enum: every conduction layer of
// the overridden section must stay inside the section rectangle (and the window).
TEST_CASE("Test_Abt725_Per_Section_Orientation_Override_Layers_Stay_Inside_Envelope", "[constructive-model][coil][abt725]") {
    settings.reset();
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt646_e16_litz_2layer_leadcollision.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    coil.set_functional_description(sourceCoil.get_functional_description());
    coil.wind();

    std::string targetName;
    WindingOrientation baseOrientation = WindingOrientation::OVERLAPPING;
    auto sectionsBefore = coil.get_sections_description().value();
    for (const auto& section : sectionsBefore) {
        if (section.get_type() == ElectricalType::CONDUCTION) {
            targetName = section.get_name();
            baseOrientation = section.get_layers_orientation();
            break;
        }
    }
    REQUIRE(!targetName.empty());
    auto flippedOrientation = baseOrientation == WindingOrientation::OVERLAPPING ? WindingOrientation::CONTIGUOUS
                                                                                 : WindingOrientation::OVERLAPPING;
    coil.set_layers_orientation(flippedOrientation, targetName);
    coil.wind();

    auto sectionsAfter = coil.get_sections_description().value();
    auto layersAfter = coil.get_layers_description().value();
    const double tol = 1e-9;
    bool checkedAnyLayer = false;
    for (const auto& section : sectionsAfter) {
        if (section.get_name() != targetName) {
            continue;
        }
        const double sx0 = section.get_coordinates()[0] - section.get_dimensions()[0] / 2;
        const double sx1 = section.get_coordinates()[0] + section.get_dimensions()[0] / 2;
        const double sy0 = section.get_coordinates()[1] - section.get_dimensions()[1] / 2;
        const double sy1 = section.get_coordinates()[1] + section.get_dimensions()[1] / 2;
        for (const auto& layer : layersAfter) {
            if (!layer.get_section() || layer.get_section().value() != targetName ||
                layer.get_type() != ElectricalType::CONDUCTION) {
                continue;
            }
            checkedAnyLayer = true;
            INFO(layer.get_name() << " in " << targetName << " (flipped to "
                 << (flippedOrientation == WindingOrientation::CONTIGUOUS ? "CONTIGUOUS" : "OVERLAPPING") << ")");
            CHECK(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 >= sx0 - tol);
            CHECK(layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2 <= sx1 + tol);
            CHECK(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 >= sy0 - tol);
            CHECK(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 <= sy1 + tol);
        }
    }
    REQUIRE(checkedAnyLayer);
    settings.reset();
}

// ABT #725 (remainder): the interlayer-insulation template stores its thickness in the dims slot
// named by the GLOBAL orientation at set_interlayer_insulation time. The winder's insertion
// branches on the section's (possibly overridden) orientation — reading the mismatched slot took
// the full window dimension as "thickness", inserting a window-sized insulation sheet inside the
// overridden section. The sheet must carry the template thickness and the section's orientation,
// and the section's layers must stay inside the window.
TEST_CASE("Test_Abt725_Interlayer_Insulation_Follows_Section_Orientation_Override", "[constructive-model][coil][abt725]") {
    settings.reset();
    const double interlayerThickness = 0.0001;
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt646_e16_litz_2layer_leadcollision.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    coil.set_functional_description(sourceCoil.get_functional_description());
    coil.set_interlayer_insulation(interlayerThickness, std::nullopt, std::nullopt, false);
    coil.wind();

    std::string targetName;
    WindingOrientation baseOrientation = WindingOrientation::OVERLAPPING;
    auto sectionsBefore = coil.get_sections_description().value();
    for (const auto& section : sectionsBefore) {
        if (section.get_type() == ElectricalType::CONDUCTION) {
            targetName = section.get_name();
            baseOrientation = section.get_layers_orientation();
            break;
        }
    }
    REQUIRE(!targetName.empty());
    auto flippedOrientation = baseOrientation == WindingOrientation::OVERLAPPING ? WindingOrientation::CONTIGUOUS
                                                                                 : WindingOrientation::OVERLAPPING;
    coil.set_layers_orientation(flippedOrientation, targetName);
    coil.wind();

    auto bobbinOut = coil.resolve_bobbin();
    auto windingWindow = bobbinOut.get_processed_description().value().get_winding_windows()[0];
    const double windowLeft = windingWindow.get_coordinates().value()[0] - windingWindow.get_width().value() / 2;
    const double windowRight = windingWindow.get_coordinates().value()[0] + windingWindow.get_width().value() / 2;
    const double windowBottom = windingWindow.get_coordinates().value()[1] - windingWindow.get_height().value() / 2;
    const double windowTop = windingWindow.get_coordinates().value()[1] + windingWindow.get_height().value() / 2;
    const double tol = 1e-9;

    const auto layersAfter = coil.get_layers_description().value();
    size_t insulationInTarget = 0;
    for (const auto& layer : layersAfter) {
        if (!layer.get_section() || layer.get_section().value() != targetName) {
            continue;
        }
        INFO(layer.get_name() << " " << layer.get_dimensions()[0] << "x" << layer.get_dimensions()[1]);
        if (layer.get_type() == ElectricalType::INSULATION) {
            ++insulationInTarget;
            CHECK(layer.get_orientation() == flippedOrientation);
            double thickness = flippedOrientation == WindingOrientation::CONTIGUOUS ? layer.get_dimensions()[1]
                                                                                    : layer.get_dimensions()[0];
            CHECK_THAT(thickness, Catch::Matchers::WithinRel(interlayerThickness, 1e-9));
        }
        // Conduction and insulation alike: nothing in the overridden section may leave the window.
        CHECK(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 >= windowLeft - tol);
        CHECK(layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2 <= windowRight + tol);
        CHECK(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 >= windowBottom - tol);
        CHECK(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 <= windowTop + tol);
    }
    REQUIRE(insulationInTarget > 0);
    settings.reset();
}

// ABT #721: margin equalization on RECTANGULAR windows (owner ruling 2026-08-14).
// Rectangular semantics: adjacent sections carry tape on the SAME two window edges, so
// equalization pairs same-index margins (left[0] with right[0], left[1] with right[1]),
// conserving each edge's creepage-path total while splitting it proportionally to each
// section's allotted space. Two windings -> one pair: clean assertions.
TEST_CASE("Test_Abt721_Equalize_Margins_Rectangular_Proportional_Split", "[constructive-model][coil][margin][abt721]") {
    settings.reset();   // coilEqualizeMargins defaults to true — the feature under test
    std::vector<int64_t> numberTurns = {34, 10};   // deliberately asymmetric section spaces
    std::vector<int64_t> numberParallels = {1, 1};
    const double margin = 0.001;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    REQUIRE(coil.get_sections_description_conduction().size() == 2);

    // Hand section 0 a margin on both edges; the re-wind's equalization shares each
    // edge's total with section 1 proportionally to the sections' widths (OVERLAPPING:
    // margins live on the height, so widths ARE the allotted spaces).
    coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});

    auto conductionSections = coil.get_sections_description_conduction();
    auto margin0 = OpenMagnetics::Coil::resolve_margin(conductionSections[0]);
    auto margin1 = OpenMagnetics::Coil::resolve_margin(conductionSections[1]);

    // The split basis is each section's ALLOTTED space, which for equal wires is the
    // winding's share of turn area — i.e. its turn count. (Post-wind section widths
    // drift slightly from the allotment via insulation subtraction, so they are NOT the
    // basis; 1e-5 absorbs the wire-area rounding in the allotment itself.)
    double expectedShare0 = margin * 34.0 / (34.0 + 10.0);
    for (size_t edge : {size_t(0), size_t(1)}) {
        INFO("edge " << edge);
        // Conservation: each edge's creepage-path total is exactly what was handed in.
        CHECK_THAT(margin0[edge] + margin1[edge], Catch::Matchers::WithinAbs(margin, 1e-9));
        // Proportionality: split by allotted space, bigger section carries more tape.
        CHECK_THAT(margin0[edge], Catch::Matchers::WithinAbs(expectedShare0, 1e-5));
        CHECK(margin0[edge] > margin1[edge]);
        CHECK(margin1[edge] > 0);
    }
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

// ABT #613: on a LOW-PROFILE window (EL 25/4.3: 2.80 mm, ~6.5 slots of Round 0.4) the U order
// used to COLLAPSE the section after the first intervening crossing — Primary section 1 kept
// ~2 of 6.5 slots (height 0.86-1.09 mm) and sprawled into 8-13 single/two-turn layers
// (overlapping factor up to 2.06) while the Z twin wound clean. Fixed by the ABT #616 fixpoint
// honesty work and finished by ABT #780's shared terminal rows (one row per winding side
// instead of one per parallel, so the lead cascade stops eating the shallow window). This pins
// the ticket's own acceptance shape: every conduction section keeps a healthy share of the
// window and no section sprawls into degenerate layers.
TEST_CASE("Test_Abt613_U_LowProfile_Sections_Do_Not_Collapse", "[constructive-model][coil][real-geometry][abt613]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    // The fixture carries only 'magnetic' (no inputs), so build the magnetic directly.
    std::ifstream fixtureFile((dataDir / "abt613_psps_el2543_U.json").string());
    json fixtureJson = json::parse(fixtureFile);
    auto magnetic = OpenMagnetics::magnetic_autocomplete(OpenMagnetics::Magnetic(fixtureJson["magnetic"]));
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_turns_description());

    auto bobbin = coil.resolve_bobbin();
    const double windowHeight =
        bobbin.get_processed_description().value().get_winding_windows()[0].get_height().value();
    const auto sections = coil.get_sections_description_conduction();
    REQUIRE(sections.size() == 4);   // PSPS: pattern {0,1} x 2

    const auto layers = coil.get_layers_description().value();
    const auto turns = coil.get_turns_description().value();
    for (const auto& section : sections) {
        size_t conductionLayers = 0;
        size_t minimumTurnsOnALayer = std::numeric_limits<size_t>::max();
        for (const auto& layer : layers) {
            if (layer.get_type() != ElectricalType::CONDUCTION || !layer.get_section() ||
                layer.get_section().value() != section.get_name()) {
                continue;
            }
            ++conductionLayers;
            size_t turnsOnLayer = 0;
            for (const auto& turn : turns) {
                if (turn.get_layer() && turn.get_layer().value() == layer.get_name()) {
                    ++turnsOnLayer;
                }
            }
            minimumTurnsOnALayer = std::min(minimumTurnsOnALayer, turnsOnLayer);
        }
        INFO(section.get_name() << " height " << section.get_dimensions()[1] * 1e3 << " mm of "
             << windowHeight * 1e3 << " mm window, " << conductionLayers
             << " conduction layers, thinnest layer " << minimumTurnsOnALayer << " turns");
        // The ticket's collapse: ~2 of 6.5 slots (< 40% of the window). Healthy sections keep
        // a real share — the pre-fix failure sat at 31-39%, the fixed layouts at 54%+.
        CHECK(section.get_dimensions()[1] > windowHeight * 0.5);
        // ...and do not sprawl: the failing state was 8-13 layers with SINGLE-turn layers
        // (a 2-turn remainder layer on a 12-turn section is a legitimate tail).
        REQUIRE(conductionLayers > 0);
        // ABT #685 raised this from 4: every layer now also holds its own opening crossing,
        // so a section of the same turns legitimately needs one more layer. The point of the
        // check is that it does not SPRAWL (the ticket's failure was 8-13), not the exact bound.
        CHECK(conductionLayers <= 5);
        CHECK(minimumTurnsOnALayer >= 2);
    }
    settings.reset();
}

// ABT #611: OpenMagnetics::Winding shadows MAS::CoilFunctionalDescription::wire (the derived
// member holds the richer OpenMagnetics::Wire) and OpenMagnetics::Coil shadows MAS::Coil::bobbin
// the same way. Serializing through the BASE — the natural way to emit schema-clean MAS json —
// used to see only the defaults: wire 'Dummy' (a file that re-winds with a 12.5 um dummy wire)
// and an empty bobbin. The setters now keep the base members in sync.
TEST_CASE("Test_Abt611_Base_Serialization_Sees_Resolved_Wire_And_Bobbin", "[constructive-model][coil][abt611]") {
    settings.reset();
    auto coil = OpenMagneticsTesting::get_quick_coil({10}, {1}, "PQ 28/20", 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    REQUIRE(coil.get_turns_description());

    const auto windings = coil.get_functional_description();
    REQUIRE(!windings.empty());
    json windingJson;
    MAS::to_json(windingJson, static_cast<const MAS::CoilFunctionalDescription&>(windings[0]));
    INFO("base-serialized wire: " << windingJson.at("wire").dump().substr(0, 120));
    if (windingJson.at("wire").is_string()) {
        CHECK(windingJson.at("wire").get<std::string>() != "Dummy");
    }
    else {
        // The resolved Wire object made it through the base.
        CHECK(windingJson.at("wire").is_object());
        CHECK(windingJson.at("wire").contains("type"));
    }

    json coilJson;
    MAS::to_json(coilJson, static_cast<const MAS::Coil&>(coil));
    INFO("base-serialized bobbin: " << coilJson.at("bobbin").dump().substr(0, 120));
    if (coilJson.at("bobbin").is_string()) {
        CHECK(coilJson.at("bobbin").get<std::string>() != "Dummy");
    }
    else {
        CHECK(coilJson.at("bobbin").is_object());
        CHECK(coilJson.at("bobbin").contains("processedDescription"));
    }
    settings.reset();
}

// ABT #724 (owner ruling, Alf 2026-08-15): margins recovered from a previous wind FOLLOW THE
// WINDING. A consumer that arrives with only the wound description (the MVB++/PyOM round-trip)
// and re-winds with a DIFFERENT layout — here repetitions 1 -> 2 — must see every conduction
// section carry ITS OWN winding's margin, not whatever sat at that ordinal in the old layout
// (positionally, the new layout's second half got {0,0} and P1/S1 wound into the creepage band).
TEST_CASE("Test_Abt724_Recovered_Margins_Follow_The_Winding", "[constructive-model][coil][margin][abt724]") {
    settings.reset();
    settings.set_coil_equalize_margins(false);
    settings.set_coil_allow_margin_tape(true);
    std::vector<int64_t> numberTurns = {34, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    const std::vector<double> primaryMargin{0.0004, 0.0003};
    const std::vector<double> secondaryMargin{0.0002, 0.0001};

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    REQUIRE(coil.get_sections_description_conduction().size() == 2);
    coil.add_margin_to_section_by_index(0, primaryMargin);
    coil.add_margin_to_section_by_index(1, secondaryMargin);

    // The round-trip: a fresh Coil handed only the wound description.
    OpenMagnetics::Coil recoveredCoil;
    recoveredCoil.set_bobbin(coil.resolve_bobbin());
    recoveredCoil.set_functional_description(coil.get_functional_description());
    recoveredCoil.set_sections_description(coil.get_sections_description());
    REQUIRE(recoveredCoil.wind({0.5, 0.5}, {0, 1}, 2));

    const auto conductionSections = recoveredCoil.get_sections_description_conduction();
    REQUIRE(conductionSections.size() == 4);   // interleaved: P0 S0 P1 S1
    const std::string primaryName = coil.get_functional_description()[0].get_name();
    for (const auto& section : conductionSections) {
        auto margin = OpenMagnetics::Coil::resolve_margin(section);
        const auto& expected = section.get_partial_windings()[0].get_winding() == primaryName
                                   ? primaryMargin
                                   : secondaryMargin;
        INFO(section.get_name() << " margin {" << margin[0] << "," << margin[1] << "}");
        CHECK_THAT(margin[0], Catch::Matchers::WithinRel(expected[0], 1e-9));
        CHECK_THAT(margin[1], Catch::Matchers::WithinRel(expected[1], 1e-9));
    }
    settings.reset();
}

// ABT #721: a rectangular window does NOT wrap — the last and first sections meet the
// bobbin walls, not each other. Handing a margin only to the LAST of three sections must
// redistribute it with its real neighbour (the middle section) and leave the FIRST
// section's margins exactly zero; the pre-#721 wrap pairing would have leaked a share of
// the last section's margin onto the first.
TEST_CASE("Test_Abt721_Equalize_Margins_Rectangular_No_Wraparound", "[constructive-model][coil][margin][abt721]") {
    settings.reset();
    std::vector<int64_t> numberTurns = {34, 25, 10};
    std::vector<int64_t> numberParallels = {1, 1, 1};
    const double margin = 0.001;

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, "PQ 28/20", 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    REQUIRE(coil.get_sections_description_conduction().size() == 3);

    coil.add_margin_to_section_by_index(2, std::vector<double>{margin, margin});

    auto conductionSections = coil.get_sections_description_conduction();
    auto margin0 = OpenMagnetics::Coil::resolve_margin(conductionSections[0]);
    auto margin1 = OpenMagnetics::Coil::resolve_margin(conductionSections[1]);
    auto margin2 = OpenMagnetics::Coil::resolve_margin(conductionSections[2]);

    // Allotted-space basis, as in the proportional-split test: same wire, so the
    // middle/last pair splits by their turn counts 25:10.
    double expectedShare1 = margin * 25.0 / (25.0 + 10.0);
    for (size_t edge : {size_t(0), size_t(1)}) {
        INFO("edge " << edge);
        // No wrap: the first section shares no boundary with the last.
        CHECK(margin0[edge] == 0.0);
        // The real neighbour pair conserves and splits the handed-in margin.
        CHECK_THAT(margin1[edge] + margin2[edge], Catch::Matchers::WithinAbs(margin, 1e-9));
        CHECK_THAT(margin1[edge], Catch::Matchers::WithinAbs(expectedShare1, 1e-5));
    }
    OpenMagneticsTesting::check_turns_description(coil);
    settings.reset();
}

TEST_CASE("Test_Abt682_U_Order_With_Margin_Keeps_The_Lead_Row_Clear", "[constructive-model][coil][real-geometry][abt682]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_allow_margin_tape(true);
    const double topMargin = 0.003;

    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    auto mas = OpenMagneticsTesting::mas_loader((dataDir / "abt646_e16_litz_2layer_leadcollision.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    auto bobbin = sourceCoil.resolve_bobbin();
    auto processedDescription = bobbin.get_processed_description().value();
    auto windingWindows = processedDescription.get_winding_windows();
    windingWindows[0].set_winding_order(MAS::WindingOrder::U);
    processedDescription.set_winding_windows(windingWindows);
    bobbin.set_processed_description(processedDescription);

    OpenMagnetics::Coil coil;
    coil.set_bobbin(bobbin);
    coil.set_functional_description(sourceCoil.get_functional_description());
    coil.preload_margins({{topMargin, 0.0}});
    coil.wind();

    REQUIRE(coil.get_turns_description());
    auto turns = coil.get_turns_description().value();
    auto wires = coil.get_wires();
    auto bobbinOut = coil.resolve_bobbin();
    auto processedOut = bobbinOut.get_processed_description().value();
    auto windowOut = processedOut.get_winding_windows()[0];
    double windowBottom = windowOut.get_coordinates().value()[1] - windowOut.get_height().value() / 2;
    double wireHeight = wires[0].get_maximum_outer_height();
    // The entrance lead runs along the bottom edge and CROSSES every layer outboard of the one it
    // attaches to, so those layers must leave that row free. The attach layer (the innermost, which
    // owns the lead's own turn) is exempt: a lead ends ON its turn by construction.
    std::string innermostLayer = turns[0].get_layer().value();
    double worstOverlap = 0;
    std::string worstDescription;
    for (const auto& turn : turns) {
        if (!turn.get_layer() || turn.get_layer().value() == innermostLayer) {
            continue;
        }
        double turnBottom = turn.get_coordinates()[1] - wireHeight / 2;
        double intrusion = (windowBottom + wireHeight) - turnBottom;
        if (intrusion > worstOverlap) {
            worstOverlap = intrusion;
            worstDescription = turn.get_name() + " at y=" + std::to_string(turn.get_coordinates()[1]);
        }
    }
    INFO("deepest intrusion into the bottom lead row " << worstOverlap * 1000 << " mm (" << worstDescription << ")");
    CHECK(worstOverlap < 1e-9);
    settings.reset();
}

TEST_CASE("Test_Abt676_Margin_Survives_A_Rewind", "[constructive-model][coil][real-geometry][abt676]") {
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_allow_margin_tape(true);
    const double topMargin = 0.003;

    auto examplesDir = std::filesystem::path{__FILE__}.parent_path().append("..").append("MAS").append("examples");
    auto mas = OpenMagneticsTesting::mas_loader((examplesDir / "01_simple_inductor_etd34_n87.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto sourceCoil = magnetic.get_mutable_coil();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    coil.set_functional_description(sourceCoil.get_functional_description());
    coil.preload_margins({{topMargin, 0.0}});
    coil.wind();

    auto bobbin = coil.resolve_bobbin();
    auto processedDescription = bobbin.get_processed_description().value();
    auto windingWindow = processedDescription.get_winding_windows()[0];
    double marginInnerFace = windingWindow.get_coordinates().value()[1]
                           + windingWindow.get_height().value() / 2 - topMargin;

    // Round-trip through JSON exactly as a consumer receives it, then wind again WITHOUT
    // preloading anything: the coil itself must be enough.
    json coilJson;
    to_json(coilJson, coil);
    OpenMagnetics::Coil reloaded(coilJson, false);
    reloaded.wind();

    REQUIRE(reloaded.get_sections_description());
    auto sections = reloaded.get_sections_description().value();
    CHECK_THAT(reloaded.resolve_margin(sections[0])[0], Catch::Matchers::WithinAbs(topMargin, 1e-12));

    REQUIRE(reloaded.get_turns_description());
    auto turns = reloaded.get_turns_description().value();
    double highestTurn = std::numeric_limits<double>::lowest();
    for (const auto& turn : turns) {
        highestTurn = std::max(highestTurn, turn.get_coordinates()[1] + turn.get_dimensions().value()[1] / 2);
    }
    INFO("highest turn after the re-wind " << highestTurn << " against " << marginInnerFace);
    CHECK(highestTurn <= marginInnerFace + 1e-9);
    settings.reset();
}

TEST_CASE("Test_Abt676_Real_Winding_Keeps_Turns_Out_Of_The_Margin", "[constructive-model][coil][real-geometry][abt676]") {
    const bool realWinding = GENERATE(false, true);
    settings.reset();
    settings.set_coil_use_real_winding_geometry(realWinding);
    settings.set_coil_allow_margin_tape(true);
    const double topMargin = 0.003;

    auto examplesDir = std::filesystem::path{__FILE__}.parent_path().append("..").append("MAS").append("examples");
    auto mas = OpenMagneticsTesting::mas_loader((examplesDir / "01_simple_inductor_etd34_n87.json").string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    // Margins reach the winder preloaded on a fresh coil, exactly as the web path builds it.
    auto sourceCoil = magnetic.get_mutable_coil();
    OpenMagnetics::Coil coil;
    coil.set_bobbin(sourceCoil.resolve_bobbin());
    coil.set_functional_description(sourceCoil.get_functional_description());
    coil.preload_margins({{topMargin, 0.0}});
    coil.wind();

    auto bobbin = coil.resolve_bobbin();
    auto processedDescription = bobbin.get_processed_description().value();
    auto windingWindow = processedDescription.get_winding_windows()[0];
    double windowTop = windingWindow.get_coordinates().value()[1] + windingWindow.get_height().value() / 2;
    double marginInnerFace = windowTop - topMargin;

    REQUIRE(coil.get_turns_description());
    auto turns = coil.get_turns_description().value();
    double highestTurn = std::numeric_limits<double>::lowest();
    for (const auto& turn : turns) {
        highestTurn = std::max(highestTurn, turn.get_coordinates()[1] + turn.get_dimensions().value()[1] / 2);
    }
    INFO("realWinding=" << realWinding << " highest turn " << highestTurn
         << " against the margin's inner face " << marginInnerFace);
    CHECK(highestTurn <= marginInnerFace + 1e-9);

    // The layers the turns live in must clear it too, or the next re-wind walks back in.
    REQUIRE(coil.get_layers_description());
    auto layers = coil.get_layers_description().value();
    for (const auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        INFO("layer " << layer.get_name());
        CHECK(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 <= marginInnerFace + 1e-9);
    }
    settings.reset();
}

TEST_CASE("Test_Abt577_Terminal_Routes_Clear_Winding_Copper", "[constructive-model][coil][real-geometry][abt577]") {
    std::vector<std::string> exampleFiles = {
        "06_llc_xfmr_eq4128_3c97.json",
        "13_current_sense_er95_n87.json",
        "14_dab_xfmr_pm8770_n97.json",
        "23_interleaved_llc_pq3530_n97.json",
        "24_margin_interleaved_flyback_pq3230_3c94.json",
    };
    auto examplesDir = std::filesystem::path{__FILE__}.parent_path().append("..").append("MAS").append("examples");
    for (auto& exampleFile : exampleFiles) {
        settings.reset();
        settings.set_coil_use_real_winding_geometry(true);
        auto mas = OpenMagneticsTesting::mas_loader((examplesDir / exampleFile).string());
        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        magnetic.get_mutable_coil().wind();
        auto& coil = magnetic.get_mutable_coil();
        INFO(exampleFile);
        REQUIRE(coil.get_turns_description());
        auto turns = coil.get_turns_description().value();
        auto wires = coil.get_wires();

        double worstClearance = std::numeric_limits<double>::max();
        std::string worstDescription;
        for (const auto& space : coil.get_connection_reserved_spaces()) {
            if (!space.layer.empty() || space.plane == OpenMagnetics::RoutePlane::FRONT_YZ) {
                continue;  // squeeze bookkeeping / out-of-plane dragbacks: not drawn in XY
            }
            // Z-order continuations are the classic in-plane crossover — they ride OVER the
            // adjacent turns at one azimuth by convention, and MKF's own blocking exempts them
            // (compute_connection_blocked_slots_per_layer skips non-terminal non-U markers). All
            // five #577 violations were TERMINAL segments; scope this sweep the same way.
            if (!space.isTerminal && coil.get_winding_order(space.section) != WindingOrder::U) {
                continue;
            }
            auto runWindingIndex = coil.get_winding_index_by_name(space.winding);
            double runBareRadius = std::min(wires[runWindingIndex].get_maximum_conducting_width(),
                                            wires[runWindingIndex].get_maximum_conducting_height()) / 2;
            double runOuterRadius = std::min(space.dimensions[0], space.dimensions[1]) / 2;
            // The rectangle is the copper's OUTER envelope; MVB++ realizes its CENTERLINE. Model the
            // copper as the long-axis centerline segment, endpoints inset by the outer radius, rotated
            // by the drawn rotation about the centre. A junction marker (square ~OD x OD on the turn)
            // degenerates to a point at the turn — excluded by the endpoint rule below.
            bool longAxisIsX = space.dimensions[0] >= space.dimensions[1];
            double halfRun = (longAxisIsX ? space.dimensions[0] : space.dimensions[1]) / 2 - runOuterRadius;
            if (halfRun < 0) {
                halfRun = 0;
            }
            double axisX = longAxisIsX ? 1.0 : 0.0;
            double axisY = longAxisIsX ? 0.0 : 1.0;
            if (std::abs(space.rotation) > 1e-9) {
                double rotationRadians = space.rotation * std::numbers::pi / 180.0;
                double rotatedX = axisX * std::cos(rotationRadians) - axisY * std::sin(rotationRadians);
                double rotatedY = axisX * std::sin(rotationRadians) + axisY * std::cos(rotationRadians);
                axisX = rotatedX;
                axisY = rotatedY;
            }
            double endAX = space.coordinates[0] - axisX * halfRun;
            double endAY = space.coordinates[1] - axisY * halfRun;
            double endBX = space.coordinates[0] + axisX * halfRun;
            double endBY = space.coordinates[1] + axisY * halfRun;
            for (auto& turn : turns) {
                // ABT #685 (Alf): a winding's parallels share ONE terminal row, their leads
                // separated in ANGLE — so a lead's 2D rectangle passing its OWN winding's sibling
                // turns is the design (the vertical stub climbs at its own azimuth, and helical
                // turns clear it there). The #577 contract this sweep enforces is about crossing
                // OTHER windings' copper, which no azimuth can fix in-plane.
                if (turn.get_winding() == space.winding) {
                    continue;
                }
                double turnX = turn.get_coordinates()[0];
                double turnY = turn.get_coordinates()[1];
                auto turnWindingIndex = coil.get_winding_index_by_name(turn.get_winding());
                double turnOuterRadius = std::min(turn.get_dimensions().value()[0], turn.get_dimensions().value()[1]) / 2;
                double turnBareRadius = std::min(wires[turnWindingIndex].get_maximum_conducting_width(),
                                                 wires[turnWindingIndex].get_maximum_conducting_height()) / 2;
                // Junction exclusion: the segment ends ON the turn it serves (marker covers the wire
                // surface at the turn), so a turn within outer-radii contact of an ENDPOINT is the
                // intended connection, not a violation.
                double contactDistance = runOuterRadius + turnOuterRadius + 1e-6;
                if (std::hypot(turnX - endAX, turnY - endAY) <= contactDistance ||
                    std::hypot(turnX - endBX, turnY - endBY) <= contactDistance) {
                    continue;
                }
                // Point-to-segment distance, then bare-copper clearance as the ticket measures it.
                double segmentDX = endBX - endAX;
                double segmentDY = endBY - endAY;
                double segmentLengthSquared = segmentDX * segmentDX + segmentDY * segmentDY;
                double projection = segmentLengthSquared > 0
                    ? std::clamp(((turnX - endAX) * segmentDX + (turnY - endAY) * segmentDY) / segmentLengthSquared, 0.0, 1.0)
                    : 0.0;
                double nearestX = endAX + projection * segmentDX;
                double nearestY = endAY + projection * segmentDY;
                double clearance = std::hypot(turnX - nearestX, turnY - nearestY) - runBareRadius - turnBareRadius;
                if (clearance < worstClearance) {
                    worstClearance = clearance;
                    worstDescription = space.winding + " p" + std::to_string(space.parallel) +
                                       (space.isTerminal ? " TERM" : " trans") + " vs " + turn.get_name();
                }
            }
        }
        // Through-layer verticals (the 14_dab class): any drawn rectangle overlapping a conduction
        // layer's radial band while spanning more than two of that layer's pitches vertically.
        // The 14_dab failure class: a terminal segment drawn down a conduction layer's own radial
        // band for a large number of that layer's pitches (the ticket's exit vertical spanned the
        // full 43.7 mm window — 45 pitches — 6.5 um from its own helix). Track the worst such span;
        // the legitimate residue is the short climb stubs at the window edge (<= 3.5 pitches).
        double worstThroughLayerSpanPitches = 0;
        auto layersForProbe = coil.get_layers_description().value();
        for (const auto& space : coil.get_connection_reserved_spaces()) {
            if (!space.layer.empty() || space.plane == OpenMagnetics::RoutePlane::FRONT_YZ) {
                continue;
            }
            if (!space.isTerminal && coil.get_winding_order(space.section) != WindingOrder::U) {
                continue;  // Z crossover convention, as above
            }
            for (auto& layer : layersForProbe) {
                if (layer.get_type() != ElectricalType::CONDUCTION) {
                    continue;
                }
                auto layerWindingIndex = coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
                double layerPitch = wires[layerWindingIndex].get_maximum_outer_height();
                bool overlapsRadially = std::abs(space.coordinates[0] - layer.get_coordinates()[0]) <
                                        (space.dimensions[0] + wires[layerWindingIndex].get_maximum_outer_width()) / 2;
                if (overlapsRadially && layerPitch > 0) {
                    worstThroughLayerSpanPitches = std::max(worstThroughLayerSpanPitches,
                                                            space.dimensions[1] / layerPitch);
                }
            }
        }
        INFO(exampleFile << ": worst bare clearance " << worstClearance * 1000 << " mm (" << worstDescription
                         << "), worst through-layer span " << worstThroughLayerSpanPitches << " pitches");
        // The #577 contract: every drawn terminal route keeps its bare copper clear of every turn...
        CHECK(worstClearance > 0);
        // ...and no terminal segment runs down a layer's radial band for more than a handful of its
        // pitches. The bound tolerates the edge climb stubs (3.5 pitches today) with headroom, while
        // the 14_dab regression it guards against sat at 45.
        CHECK(worstThroughLayerSpanPitches < 8);
    }
    settings.reset();
}

// Visual + numeric preview harness for the SPREAD placer rules: ABT #578 (N-filar bundles stay
// TOUCHING, slack goes in the gaps BETWEEN bundles) and ABT #579 (fence-post: the outermost
// turns' surfaces sit on the layer edges, no dead margin at the flanges). Paints every case with
// real winding geometry OFF and ON so the two layouts can be compared, and prints the turn
// stations of the tallest layer so the spacing can be checked by eye against the numbers.
// Not an assertion test -- it exists to produce the SVGs and the station tables.
TEST_CASE("Test_Spread_Preview_Svg", "[constructive-model][coil][svg-preview]") {
    struct PreviewCase { std::string file; std::string label; };
    std::vector<PreviewCase> previewCases = {
        {"15_gan_inductor_e138_3f4.json", "15_gan"},  // ABT #579 repro: 6 turns, 1 parallel
        {"06_llc_xfmr_eq4128_3c97.json", "06_llc"},   // ABT #577/#578: 3 and 3 parallels
        {"14_dab_xfmr_pm8770_n97.json", "14_dab"},    // ABT #577/#578: 4 and 4 parallels
        {"23_interleaved_llc_pq3530_n97.json", "23_llc"},  // ABT #577: 2 and 4 parallels
    };
    auto examplesDir = std::filesystem::path{__FILE__}.parent_path().append("..").append("MAS").append("examples");

    for (auto& previewCase : previewCases) {
        for (bool realWinding : {false, true}) {
            settings.reset();
            settings.set_coil_use_real_winding_geometry(realWinding);
            auto mas = OpenMagneticsTesting::mas_loader((examplesDir / previewCase.file).string());
            auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
            // Several examples ship fully wound, so autocomplete alone would not re-run the placer
            // and the real-winding flag would never engage. Force the wind so both variants are real.
            magnetic.get_mutable_coil().wind();

            std::string variant = realWinding ? "real" : "ideal";
            auto outFile = outputFilePath;
            outFile.append("spread_" + previewCase.label + "_" + variant + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            // Terminal leads (magenta) and inter-layer/section transitions (blue). Drawn from the
            // coil's connection reserved spaces, which only exist under real winding geometry — the
            // count is reported below so an empty ideal drawing is visibly a real absence of
            // connections rather than a forgotten paint call.
            painter.paint_coil_connections(magnetic);
            painter.export_svg();
            std::cout << previewCase.label << " [" << variant << "] connection reserved spaces: "
                      << magnetic.get_mutable_coil().get_connection_reserved_spaces().size() << std::endl;

            auto coil = magnetic.get_coil();
            if (!coil.get_turns_description()) {
                std::cout << previewCase.label << " [" << variant << "] produced NO turns" << std::endl;
                continue;
            }
            auto turns = coil.get_turns_description().value();
            // Report the conduction layer holding the most turns -- the one where the spacing rule
            // is most visible.
            std::map<std::string, std::vector<double>> centresPerLayer;
            for (auto& turn : turns) {
                if (turn.get_layer()) {
                    centresPerLayer[turn.get_layer().value()].push_back(turn.get_coordinates()[1]);
                }
            }
            std::string tallestLayer;
            size_t tallestCount = 0;
            for (auto& [layerName, centres] : centresPerLayer) {
                if (centres.size() > tallestCount) {
                    tallestCount = centres.size();
                    tallestLayer = layerName;
                }
            }
            if (tallestLayer.empty()) {
                continue;
            }
            auto centres = centresPerLayer[tallestLayer];
            std::sort(centres.begin(), centres.end());
            double wireOuterHeight = 0;
            auto layersForPreview = coil.get_layers_description().value();
            for (auto& layer : layersForPreview) {
                if (layer.get_name() == tallestLayer) {
                    auto windingIndex = coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
                    wireOuterHeight = coil.get_wires()[windingIndex].get_maximum_outer_height();
                    // Measure the dead margin against the WINDING WINDOW, not the layer:
                    // delimit_and_compact shrinks each layer to its turns' own bounding box, so a
                    // layer-relative margin is 0 by construction and hides the flange gap entirely.
                    auto bobbin = coil.resolve_bobbin();
                    auto windowDimensions = bobbin.get_winding_window_dimensions(0);
                    auto windowCoordinates = bobbin.get_winding_window_coordinates(0);
                    double windowLow = windowCoordinates[1] - windowDimensions[1] / 2;
                    double windowHigh = windowCoordinates[1] + windowDimensions[1] / 2;
                    std::cout << previewCase.label << " [" << variant << "] layer '" << tallestLayer
                              << "' h=" << layer.get_dimensions()[1] * 1000 << " mm centred at "
                              << layer.get_coordinates()[1] * 1000 << " mm, OD " << wireOuterHeight * 1000
                              << " mm, " << centres.size() << " turns; window h="
                              << windowDimensions[1] * 1000 << " mm" << std::endl;
                    std::cout << "    dead margin to FLANGE: low "
                              << ((centres.front() - wireOuterHeight / 2) - windowLow) * 1000 << " mm, high "
                              << (windowHigh - (centres.back() + wireOuterHeight / 2)) * 1000 << " mm" << std::endl;
                }
            }
            // Surface-to-surface gap between consecutive turns: 0 means touching (one bundle).
            std::cout << "    gaps (mm):";
            for (size_t k = 1; k < centres.size(); ++k) {
                std::cout << " " << (centres[k] - centres[k - 1] - wireOuterHeight) * 1000;
            }
            std::cout << std::endl;
        }
    }
    settings.reset();
}

// ABT #620: a MAS file that carries only functionalDescription + bobbin (no
// sectionsDescription — the common case for a design saved before winding, e.g. the
// shipped "22_margin_tape_forward" / "24_margin_interleaved_flyback" examples) gets
// rewound by magnetic_autocomplete for painter/3D/simulation. Coil::wind() only applies
// the design's declared insulation standard (calculate_insulation) when _inputs is set
// on the coil; without it, it silently falls back to calculate_mechanical_insulation()
// — a single bare mechanical layer, ZERO margin — even when the design declares
// reinforced insulation. mas_autocomplete has the Inputs available (mas.get_inputs())
// but was not threading it into magnetic_autocomplete before this fix.
TEST_CASE("Test_ABT620_MasAutocomplete_Applies_Declared_Insulation_Margin", "[constructive-model][coil][bug]") {
    settings.reset();
    clear_databases();
    namespace fs = std::filesystem;
    auto examplesDir = fs::path{std::source_location::current().file_name()}.parent_path()
                           .append("..").append("MAS").append("examples");
    auto mas = OpenMagneticsTesting::mas_loader((examplesDir / "22_margin_tape_forward_e4218_3c95.json").string());

    REQUIRE(mas.get_inputs().get_design_requirements().get_insulation());
    REQUIRE(mas.get_inputs().get_design_requirements().get_insulation()->get_insulation_type() == IsolationClass::REINFORCED);
    // The fixture itself must carry no sectionsDescription — otherwise wind() is never
    // reached here and this test would not exercise the reported path at all.
    REQUIRE_FALSE(mas.get_magnetic().get_coil().get_sections_description());

    auto autocompleted = OpenMagnetics::mas_autocomplete(mas, false);
    auto& coil = autocompleted.get_mutable_magnetic().get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    REQUIRE(coil.get_sections_description());

    auto sections = coil.get_sections_description().value();
    bool foundNonzeroMargin = false;
    for (auto& section : sections) {
        if (!section.get_margin()) {
            continue;
        }
        auto margin = std::get<std::vector<double>>(section.get_margin().value());
        if (margin[0] > 1e-6 || margin[1] > 1e-6) {
            foundNonzeroMargin = true;
        }
    }
    // The reported symptom was margin = [0, 0] on every conduction section (single 25 um
    // mechanical layer, no standard applied). A reinforced-insulation design must place a
    // real margin somewhere in the coil once its own declared requirements are honoured.
    REQUIRE(foundNonzeroMargin);
    settings.reset();
}

// ABT #623: split off from ABT #620 (ask 2). 22_margin_tape_forward and
// 24_margin_interleaved_flyback, once their declared insulation reaches the winder
// (ABT #620), produce a layout calculate_filling_factor() reports as fitting
// (windingFits=1) -- yet Coil::wind()'s own are_sections_and_layers_fitting() verdict
// disagreed. Root-caused and fixed under ABT #624 (commit 84752f26): a dangling
// reference in the ABT #616 turn-envelope check silently no-opped it, and packed
// layers could move copper past the section rect wind_by_sections() gave them without
// anything updating that rect, so the section's own filling factor read >1 even when
// the copper fit. Both checks must agree once a design's own summary metric says it fits.
TEST_CASE("Test_ABT623_FitVerdict_Agrees_With_FillingFactor", "[constructive-model][coil][bug]") {
    for (std::string fixture : {"22_margin_tape_forward_e4218_3c95.json", "24_margin_interleaved_flyback_pq3230_3c94.json"}) {
        settings.reset();
        clear_databases();
        namespace fs = std::filesystem;
        auto examplesDir = fs::path{std::source_location::current().file_name()}.parent_path()
                               .append("..").append("MAS").append("examples");
        auto mas = OpenMagneticsTesting::mas_loader((examplesDir / fixture).string());
        auto autocompleted = OpenMagnetics::mas_autocomplete(mas, false);
        auto& coil = autocompleted.get_mutable_magnetic().get_mutable_coil();
        REQUIRE(coil.get_turns_description());

        auto fillingFactor = coil.calculate_filling_factor();
        bool fits = coil.are_sections_and_layers_fitting();
        // The reported symptom: calculate_filling_factor() said windingFits=1 while
        // are_sections_and_layers_fitting() (wind()'s own verdict) said false.
        REQUIRE(fillingFactor.windingFits);
        REQUIRE(fits);
    }
    settings.reset();
}

TEST_CASE("Test_Failed_Wind_Does_Not_Poison_The_Coil", "[coil][wind][smoke-test]") {
    // ABT #850: one FAILED sectioned wind used to poison the object permanently — the
    // ABT #676 margin recovery loaded the infeasible margins into _marginsPerSection,
    // which survived clearing the descriptions, so every later wind re-applied them and
    // failed too. A wind that fails and produces no turns must leave the coil as it
    // found it; a margin-free re-wind after clearing the sections must then succeed.
    //
    // Fixture: WE choke 744821240 as enriched by asgard — a two-winding toroid whose
    // sheet-specified 3 mm spacer fits the real cased part but not the modelled bare
    // bore, the exact catalogue case that exposed the poisoning.
    auto dataDir = std::filesystem::path{__FILE__}.parent_path().append("testData");
    std::ifstream jsonFile((dataDir / "abt850_infeasible_margins_toroid.json").string());
    json magneticJson;
    jsonFile >> magneticJson;
    OpenMagnetics::Magnetic magnetic(magneticJson);
    auto coil = magnetic.get_coil();
    REQUIRE(!coil.get_turns_description());

    bool sectionedWind = coil.wind();
    CHECK(!sectionedWind);
    CHECK(!coil.get_turns_description());

    // The poisoning regression: clearing the infeasible sections and re-winding fresh
    // must succeed. Before the ABT #850 fix this second wind failed as well, because
    // the recovered margins survived in transient members the clear cannot reach.
    coil.set_sections_description(std::nullopt);
    coil.set_layers_description(std::nullopt);
    bool freshWind = coil.wind();
    CHECK(freshWind);
    REQUIRE(coil.get_turns_description());
    CHECK(coil.get_turns_description()->size() > 0);
}
